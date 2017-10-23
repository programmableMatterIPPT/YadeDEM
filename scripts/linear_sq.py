# linear actuator of wall type, square-based

from yade import pack, plot

# 5 materials - 4 for spheres, 1 for walls
# damping coefficients
cn0=50000 
cs0=50000
# cohesion coefficients
nch=1e15
sch=1e15

O.materials.append(PmCohFrictMat(young=1e11, poisson=0.25, density=100, frictionAngle=radians(0), cn=cn0, cs=cs0, normalCohesion=nch, shearCohesion=sch, normalCohesion1=nch, shearCohesion1=sch, num=0, momentRotationLaw=True, label='bluesphere'))

O.materials.append(PmCohFrictMat(young=1e11, poisson=0.25, density=100, frictionAngle=radians(0), cn=cn0, cs=cs0, normalCohesion=nch, shearCohesion=sch, num=1, momentRotationLaw=False, label='redsphere'))

O.materials.append(PmCohFrictMat(young=1e11, poisson=0.25, density=100, frictionAngle=radians(0), cn=cn0, cs=cs0, normalCohesion2=nch, shearCohesion2=sch, normalCohesion3=nch, shearCohesion3=sch, num=2, momentRotationLaw=True, label='greensphere'))

O.materials.append(PmCohFrictMat(young=1e11, poisson=0.25, density=100, frictionAngle=radians(0), cn=cn0, cs=cs0, normalCohesion2=nch, shearCohesion2=sch, num=3, momentRotationLaw=False, label='yellowsphere'))

O.materials.append(FrictMat(young=1e12, poisson=.25, frictionAngle=0, density=1000, label='frictionlessWalls'))


import math
s2=math.sqrt(2)
s3=math.sqrt(3)

Nz=11 # height
Nx=5 # number of concentric walls
r=0.5 # sphere radius
fr=1.001 # interaction distance enlargement factor
trq=20000 # actuation torque on sphere
istep = 1000000
iMax = 1000*istep

Vz=Vector3(0,0,2)*r
Vx=Vector3(2,0,0)*r
Vy=Vector3(0,2,0)*r
V1=Vector3(2,2,0)*r
V2=Vector3(-2,2,0)*r

Mask = [1,3,7,15]

trqId = []
trqV = []

def bounds(top,bottom):
	for i in range(-Nx+1,Nx):
		for j in range(-Nx+1,Nx):
			O.bodies.append([utils.sphere(i*Vx+j*Vy+Nz*Vz,r,material=top[0],color=top[1],mask=-1)])
			O.bodies.append([utils.sphere(i*Vx+j*Vy-Vz,r,material=bottom[0],color=bottom[1],mask=-1)])
#----------------------------------------------------------------------
# elementary cell - with protrusions
#[0 material, 1 color, 2 setMinIndex, 3 step of index in column, 4 torque direction, 5 blockedDOFs ('xyzXYZ')]

#red=['redsphere',Vector3(1,0,0),True,0,1,'']
#green=['greensphere',Vector3(0,1,0),False,-1,0,'']
#blue=['bluesphere',Vector3(0,0,1),False,1,0,'']
#yellow=['yellowsphere',Vector3(1,1,0),True,0,-1,'']

#sph = [[[green,green],[green,green],[blue,blue],[blue,blue]],
#       [[green,green],[red,red],[blue,blue],[yellow,yellow]]]
#----------------------------------------------------------------------
# elementary cell - no protrusions, cohesion outwards
#[0 material, 1 color, 2 setMinIndex, 3 step of index in column, 4 torque direction, 5 blockedDOFs ('xyzXYZ')]

#red=['redsphere',Vector3(1,0,0),True,0,1,'']
#green=['greensphere',Vector3(0,1,0),False,-1,0,'']
#blue=['bluesphere',Vector3(0,0,1),False,1,0,'']
#yellow=['yellowsphere',Vector3(1,1,0),True,0,-1,'']

#sph = [[[green,red],[blue,yellow]],
#       [[green,green],[blue,blue]]]
		
#bounds(green,blue)
#----------------------------------------------------------------------
# elementary cell - no protrusions, cohesion inwards
#[0 material, 1 color, 2 setMinIndex, 3 step of index in column, 4 torque direction, 5 blockedDOFs ('xyzXYZ')]

red=['redsphere',Vector3(1,0,0),True,0,-1,'']
green=['greensphere',Vector3(0,1,0),False,-1,0,'']
blue=['bluesphere',Vector3(0,0,1),False,1,0,'']
yellow=['yellowsphere',Vector3(1,1,0),True,0,1,'']

sph = [[[blue,yellow],[green,red]],
       [[blue,blue],[green,green]]]
		
bounds(green,blue)
#----------------------------------------------------------------------
# add cells to the scene
for k in range(Nz):
	O.bodies.append([utils.sphere(k*Vz,r,material=green[0],color=green[1],mask=-1)])
	O.bodies[-1].state.index = 3*Nx*istep+green[3]*k+1
	O.bodies[-1].state.indexn = 3*Nx*istep+green[3]*k
	O.bodies[-1].state.blockedDOFs = green[5]
	k0 = k%len(sph)
	for i in range(1,Nx):
		i0 = (i-1)%len(sph[0])
		for j in range(0,2*i):
			j0 = j%len(sph[0][0])
			sph0 = sph[k0][i0][j0]
			if sph0 == 0: continue
#			ind = 3*(Nx+i)*istep+sph0[3]*k+1 #when cohesion outwards
			ind = 3*(Nx-i)*istep+sph0[3]*k+1 #when cohesion inwards
#			if j == 0: ind = ind+istep #increase indexes at corners
			Cen = [i*V1-j*Vx,i*V2-j*Vy,-i*V1+j*Vx,-i*V2+j*Vy]
			if sph0[2]:
#				if i==Nx: continue   #remove outer active spheres
				trqId.extend(range(O.bodies[-1].id+1,O.bodies[-1].id+5))
				if j == 0: trqV.extend([sph0[4]*trq*Vector3(-c[1],c[0],0)/c.norm() for c in Cen])
				else: trqV.extend([sph0[4]*trq*c/c.norm() for c in [-Vx,-Vy,Vx,Vy]])
			NewSpheres = [utils.sphere(Cen[l]+k*Vz,r,material=sph0[0],color=sph0[1],mask=Mask[l]) for l in range(4)]
			for b in NewSpheres:
				b.state.setMinIndex = sph0[2]
				b.state.index = ind
				b.state.indexn = ind - 1
				b.state.blockedDOFs = sph0[5]
			O.bodies.append(NewSpheres)
			

O.bodies.append(utils.wall(min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)]),axis=2,sense=1,material='frictionlessWalls'))
#O.bodies.append(utils.wall(max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)]),axis=2,sense=-1,material='frictionlessWalls'))
#plate=O.bodies[-1]
#plate.state.vel=(0,0,-0.5)

#define engines:
O.engines=[
	ForceResetter(),
	PyRunner(command='moment()',iterPeriod=1),
	InsertionSortCollider([Bo1_Sphere_Aabb(aabbEnlargeFactor=fr),Bo1_Wall_Aabb()]),
	PmControlApplier(),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom6D(interactionDetectionFactor=fr),Ig2_Wall_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_FrictPhys(),PmIp2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionNow=True,setCohesionOnNewContacts=True)],
		[Law2_ScGeom_FrictPhys_CundallStrack(),PmLaw2_ScGeom6D_CohFrictPhys_CohesionMoment(useIncrementalForm=False,always_use_moment_law=True)]
	),
	NewtonIntegrator(damping=0,gravity=(0,0,-9.81)),
#	PyRunner(command='addPlotData()',iterPeriod=10,label='plotter'),
]

# function for plotting the vertical force acting on the upper wall, used in PyRunner
#def addPlotData():
#   plot.addData(Fz=O.forces.f(plate.id)[2],w=plate.state.pos[2]-plate.state.refPos[2],unbalanced=utils.unbalancedForce(),i=O.iter)

Gxy=Matrix3(
  1,  0,  0,
  0,  1,  0,
  0,  0,  0)

diffr = r*(2-(s2+s3)/2)

def moment():
	toRemove=[]
	for i in range(len(trqId)):
	      if (Gxy*O.bodies[trqId[i]].state.pos-Gxy*O.bodies[trqId[i]].state.refPos).norm()<diffr:
		      O.forces.addT(trqId[i],trqV[i])
	      else:
		      toRemove.append(i)
		      O.bodies[trqId[i]].state.setMinIndex=False
		      O.bodies[trqId[i]].state.minIndex=iMax
		      print("Switch off. ID=",trqId[i],", old minIndex=",O.bodies[trqId[i]].state.minIndex,", new minIndex=",iMax)
	toRemove1=sorted(toRemove,reverse=True) 
	for i in range(len(toRemove1)):
		del trqId[toRemove1[i]]
		del trqV[toRemove1[i]]


#set time step and run simulation:
O.dt=0.5*utils.PWaveTimeStep()

Gl1_Sphere.stripes=True
Gl1_Sphere.quality=1
from yade import qt
qt.View()
qtr = qt.Renderer()
qtr.lightPos=[-1,-1,1]
#qtr.light2Color=[1,1,1]
#qtr.light2Pos=[-130,-75,30]
qtr.mask=14

#O.run(1000,True)
