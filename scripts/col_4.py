# linear actuator of column type, with 4 active modules around a column
from yade import pack, plot

# 4 materials - 3 for spheres, 1 for walls
# damping coefficients
cn0=50000 
cs0=50000
# cohesion coefficients
nch=1e15
sch=1e15

O.materials.append(PmCohFrictMat(young=1e11, poisson=0.25, density=100, frictionAngle=radians(0), cn=cn0, cs=cs0, normalCohesion=nch, shearCohesion=sch, normalCohesion1=nch, shearCohesion1=sch,  num=0, momentRotationLaw=True, label='blue'))
O.materials.append(PmCohFrictMat(young=1e11, poisson=0.25, density=100, frictionAngle=radians(0), cn=cn0, cs=cs0, normalCohesion=nch, shearCohesion=sch, num=1, momentRotationLaw=False, label='red'))
O.materials.append(PmCohFrictMat(young=1e11, poisson=0.25, density=1000, frictionAngle=radians(0), cn=cn0, cs=cs0, normalCohesion2=nch, shearCohesion2=sch, num=2, momentRotationLaw=True, label='grey'))
O.materials.append(FrictMat(young=1e12, poisson=.25, frictionAngle=0, density=1000, label='frictionlessWalls'))

s2=sqrt(2)
s3=sqrt(3)

Nx=3 # number of cells in x-direction
Ny=3 # number of cells in y-direction
Nz=5 # number of cells in z-direction
r=0.5 # sphere radius
fr=1.001 # interaction distance enlargement factor
trq=26000 # actuation torque on sphere
iMax = 10*Nx*Ny*Nz

# "crystallographic" vectors for the actuator
Vx=Vector3(4,2,0)
Vy=Vector3(-2,4,0)
Vz=Vector3(0,0,4)

cblue = Vector3(0,0,1)
cred = Vector3(1,0,0)
cgrey = Vector3(.5,.5,.5)

#---------------------------------------------------------------
trqId = []
trqV = []

# definition of a cell (repetitive pattern of spheres in the structure) composed of 10 spheres
Mat = ['blue','red','red','red','red','blue','grey','grey','grey','grey'] # materials
Col = [cblue,cred,cred,cred,cred,cblue,cgrey,cgrey,cgrey,cgrey] # colors
Tor = [0,Vector3(0,-1,0),Vector3(1,0,0),Vector3(0,1,0),Vector3(-1,0,0),0,0,0,0,0] # torque directions
Pos = [Vector3(0,0,0),Vector3(2,0,0),Vector3(0,2,0),Vector3(-2,0,0),Vector3(0,-2,0),Vector3(0,0,2),Vector3(2,0,2),Vector3(0,2,2),Vector3(-2,0,2),Vector3(0,-2,2)] # positions in the cell
Inc = [0,iMax,iMax,iMax,iMax,1,iMax,iMax,iMax,iMax] # index increments
Smi = [False,True,True,True,True,False,False,False,False,False] # setMinIndex

def addsphere(pos,mat,col,smi,ind):
	NewSphere = utils.sphere(pos,r,material=mat,color=col)
	NewSphere.state.setMinIndex = smi
	NewSphere.state.minIndex = 0
	NewSphere.state.index = ind
	NewSphere.state.indexn = ind-1
	O.bodies.append([NewSphere])

def addcell(i,j,k):
	for l in range(len(Mat)):
		addsphere((i*Vx+j*Vy+k*Vz+Pos[l])*r,Mat[l],Col[l],Smi[l],2*k+Inc[l]+1)
		if Smi[l]:
			trqId.extend([O.bodies[-1].id])
			trqV.extend([Tor[l]*trq])

def addbottom():
	for i in range(Nx):
		for j in range(Ny):
			for l in range(5): addsphere((i*Vx+j*Vy-Vz/2+Pos[l])*r,'blue',cblue,False,iMax)

def addsideline(v):
	for k in range(2*Nz): addsphere((v+k*Vz/2)*r,'grey',cgrey,False,iMax)

def addsides():
	for i in range(Nx): 
		addsideline(i*Vx+Vector3(2,-2,0)) # bottom
		addsideline((Ny-1)*Vy+i*Vx+Vector3(-2,2,0)) # top
	for j in range(Ny): 
		addsideline(j*Vy+Vector3(-2,-2,0)) # left
		addsideline((Nx-1)*Vx+j*Vy+Vector3(2,2,0)) # right
	
#---------------------------------------------------------------		
# add cells to the scene 
for k in range(Nz):
	for j in range(Ny):
		for i in range(Nx): addcell(i,j,k)
			
addbottom() # add bottom
addsides() # add side bracings
O.bodies[0].state.blockedDOFs = 'xyzXYZ'

#O.bodies.append(utils.wall(min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)]),axis=2,sense=1,material='frictionlessWalls',mask=-1))
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
qtv.viewDir = Vector3(0,-1,0)
qtv.upVector = Vector3(0,0,1)
qtv.showEntireScene()
qtr = qt.Renderer()
qtr.lightPos=[-1,-1,1]
#qtr.light2Color=[1,1,1]
#qtr.light2Pos=[-130,-75,30]
#qtr.mask=2

#O.run(1000,True)
