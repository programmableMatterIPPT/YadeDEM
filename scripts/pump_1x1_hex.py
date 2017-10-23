# linear actuator of pump type, hexagon-based
from yade import pack, plot

# 6 materials - 5 for spheres, 1 for walls
# damping coefficients
cn0=50000 
cs0=50000
# cohesion coefficients
nch=1e15
sch=1e15

O.materials.append(PmCohFrictMat(young=1e11, poisson=0.25, density=100, frictionAngle=radians(0), cn=cn0, cs=cs0, normalCohesion=nch, shearCohesion=sch, normalCohesion1=nch, shearCohesion1=sch, normalCohesion4=nch, shearCohesion4=sch, num=0, momentRotationLaw=True, label='bluesphere'))
O.materials.append(PmCohFrictMat(young=1e11, poisson=0.25, density=100, frictionAngle=radians(0), cn=cn0, cs=cs0, normalCohesion=nch, shearCohesion=sch, num=1, momentRotationLaw=False, label='redsphere'))
O.materials.append(PmCohFrictMat(young=1e11, poisson=0.25, density=100, frictionAngle=radians(0), cn=cn0, cs=cs0, normalCohesion2=nch, shearCohesion2=sch, normalCohesion3=nch, shearCohesion3=sch, normalCohesion4=nch, shearCohesion4=sch, num=2, momentRotationLaw=True, label='greensphere'))
O.materials.append(PmCohFrictMat(young=1e11, poisson=0.25, density=100, frictionAngle=radians(0), cn=cn0, cs=cs0, normalCohesion2=nch, shearCohesion2=sch, num=3, momentRotationLaw=False, label='yellowsphere'))
O.materials.append(PmCohFrictMat(young=1e11, poisson=0.25, density=1000, frictionAngle=radians(0), cn=cn0, cs=cs0, normalCohesion=nch, shearCohesion=sch, normalCohesion2=nch, shearCohesion2=sch, normalCohesion4=nch, shearCohesion4=sch, num=4, momentRotationLaw=True, label='greysphere'))
O.materials.append(FrictMat(young=1e12, poisson=.25, frictionAngle=0, density=1000, label='frictionlessWalls'))

s2=sqrt(2)
s3=sqrt(3)

Nz=10 # height of the actuator
Nx=9 # number of cells in x-direction
Ny=6 # number of cells in y-direction
r=0.5 # sphere radius
fr=1.001 # interaction distance enlargement factor
trq=26000 # actuation torque on sphere
iactive = 2*Nx*Ny*Nz
igrey = iactive+10
iMax = igrey+10

Vz=Vector3(0,0,2)*r
Vx=Vector3(s3,0,0)*r
Vy=Vector3(0,2,0)*r
dV=Vector3(0,-1,0)*r

# colors
#cgrey1=Vector3(.5,.5,.5); cgrey2=Vector3(.5,.5,.5); cblue=Vector3(0,0,1); cred=Vector3(1,0,0); cyellow=Vector3(1,1,0); cgreen=Vector3(0,1,0);
cgrey1=Vector3(0,0,1); cgrey2=Vector3(.5,.5,.5); cblue=Vector3(.5,.5,.5); cred=Vector3(1,0,0); cyellow=Vector3(1,1,0); cgreen=Vector3(.5,.5,.5);

#[0 material, 1 color, 2 setMinIndex, 3 step of index in column, 4 torque direction, 5 blockedDOFs ('xyzXYZ')]
blue=['bluesphere',cblue,False,1,0,'']
red1=['redsphere',cred,True,0,Vector3(0,-1,0),'']
red2=['redsphere',cred,True,0,Vector3(0,1,0),'']
green=['greensphere',cgreen,False,1,0,'']
yellow1=['yellowsphere',cyellow,True,0,Vector3(0,1,0),'']
yellow2=['yellowsphere',cyellow,True,0,Vector3(0,-1,0),'']
grey1=['greysphere',cgrey1,False,1,0,'']
grey2=['greysphere',cgrey2,False,1,0,'']

trqId = []
trqV = []

# elementary cell of the actuator
sph = [[[blue,blue],[grey2,red1],[green,green],[grey2,red2]],[[blue,blue],[grey2,yellow1],[green,green],[grey2,yellow2]]]

#---------------------------------------------------------------		
# add cells to the scene
for k in range(Nz):
	k0 = k%len(sph)
	for i in range(Nx):
		i0 = i%len(sph[0])
		for j in range(Ny+i%2):
			j0 = j%len(sph[0][0])
			sph0 = sph[k0][i0][j0]
			ind = k+2
			minind = 0
			if sph0[2]:
				if k > Nz/2:
					sph0 = grey1
					ind = igrey
					minind = igrey
				else:
					ind = iactive
					trqId.extend([O.bodies[-1].id+1])
					trqV.extend([sph0[4]*trq])
			NewSphere = utils.sphere(i*Vx+j*Vy+(i%2)*dV+k*Vz,r,material=sph0[0],color=sph0[1])
			NewSphere.state.setMinIndex = sph0[2]
			NewSphere.state.minIndex = minind
			NewSphere.state.index = ind
			NewSphere.state.indexn = ind-1
			NewSphere.state.blockedDOFs = sph0[5]
			O.bodies.append([NewSphere])
			

def bound(top):
	for i in range(Nx):
		for j in range(Ny+i%2):
			O.bodies.append([utils.sphere(i*Vx+j*Vy+(i%2)*dV+(Nz)*Vz,r,material=top[0],color=top[1])])
			O.bodies[-1].state.index = igrey
			O.bodies[-1].state.indexn = igrey-1
			O.bodies[-1].state.minIndex = igrey
		
bound(grey1)
#---------------------------------------------------------------
for b in O.bodies:
	if isinstance(b.shape,Sphere):
		if b.state.pos[1] < r:
			b.mask = 1
		else: b.mask = 3
#---------------------------------------------------------------

O.bodies.append(utils.wall(min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)]),axis=2,sense=1,material='frictionlessWalls',mask=-1))
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

def moment():
	toRemove=[]
	for i in range(len(trqId)):
	      diffr = (Gxy*O.bodies[trqId[i]].state.pos-Gxy*O.bodies[trqId[i]].state.refPos).norm()
	      if diffr < r/2:
		      O.forces.addT(trqId[i],trqV[i])
	      else:
		      toRemove.append(i)
		      O.bodies[trqId[i]].state.setMinIndex=False
		      print("Switch off. ID=",trqId[i],", old minIndex=",O.bodies[trqId[i]].state.minIndex,", diffr=",diffr)
		      O.bodies[trqId[i]].state.minIndex=iMax
	toRemove1=sorted(toRemove,reverse=True) 
	for i in range(len(toRemove1)):
		del trqId[toRemove1[i]]
		del trqV[toRemove1[i]]


#set time step and run simulation:
O.dt=0.5*utils.PWaveTimeStep()

Gl1_Sphere.stripes=True
Gl1_Sphere.quality=1
from yade import qt
qtv = qt.View()
qtv.viewDir = Vector3(0,1,0)
qtv.upVector = Vector3(0,0,1)
qtv.showEntireScene()
qtr = qt.Renderer()
qtr.lightPos=[-1,-1,1]
#qtr.light2Color=[1,1,1]
#qtr.light2Pos=[-130,-75,30]
qtr.mask=2
#qtr.mask=3

#O.run(1000,True)
