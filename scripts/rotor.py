# linear actuator of rotor type
from yade import pack, plot

Nz=8 # height of the actuator
Nc=8 # number of concentric rings
r = 0.5 # sphere radius
fr=1.0001 # interaction distance enlargement factor
trq=20000 # torque on sphere
trqId = []
trqV = []

# damping coefficients
cn0=100000 
cs0=100000

# cohesion coefficients
nch=1e15
sch=1e15

# 5 materials - 4 for spheres, 1 for walls

O.materials.append(PmCohFrictMat(young=1e11, poisson=0.25, density=100, frictionAngle=radians(0), cn=cn0, cs=cs0, normalCohesion=nch, shearCohesion=sch, normalCohesion1=nch, shearCohesion1=sch, num=0, momentRotationLaw=True, label='bluesphere'))

O.materials.append(PmCohFrictMat(young=1e11, poisson=0.25, density=100, frictionAngle=radians(0), cn=cn0, cs=cs0, normalCohesion=nch, shearCohesion=sch, num=1, momentRotationLaw=False, label='redsphere'))

O.materials.append(PmCohFrictMat(young=1e11, poisson=0.25, density=100, frictionAngle=radians(0), cn=cn0, cs=cs0, normalCohesion2=nch, shearCohesion2=sch, normalCohesion3=nch, shearCohesion3=sch, num=2, momentRotationLaw=True, label='greensphere'))

O.materials.append(PmCohFrictMat(young=1e11, poisson=0.25, density=100, frictionAngle=radians(0), cn=cn0, cs=cs0, normalCohesion2=nch, shearCohesion2=sch, num=3, momentRotationLaw=False, label='yellowsphere'))

O.materials.append(FrictMat(young=1e12, poisson=.25, frictionAngle=0, density=1000, label='frictionlessWalls'))


#[0 material, 1 color, 2 setMinIndex, 3 step of index along perimeter, 4 torque direction, 5 blockedDOFs ('xyzXYZ')]
blue=['bluesphere',Vector3(0,0,1),False,1,0,'']
red=['redsphere',Vector3(1,0,0),True,0,-1,'']
green=['greensphere',Vector3(0,1,0),False,-1,0,'']
yellow=['yellowsphere',Vector3(1,1,0),True,0,1,'']

sph = [green,green,red,blue,blue,yellow]

#Mask = [1,3,7,15,31,63]


def cylinder(no,R,nz,dz,sp,shiftind,invind):
	fi = 2*pi/no
	if invind: irange = range(no-1,-1,-1)
	else: irange = range(no)
        for i in range(no):
		x = R*cos(i*fi)
		y = R*sin(i*fi)
		for j in range(nz):
			O.bodies.append(utils.sphere(Vector3(x,y,2*j*r+dz),r,material=sp[0],color=sp[1]))
			O.bodies[-1].state.index = shiftind+irange[i]+1
			O.bodies[-1].state.indexn = shiftind+(no+irange[i]-1)%no + 1


def acylinder(n0,n1,R0,R1,dirmot,sp):
	df0 = dirmot*2*pi/n0
	df1 = dirmot*2*pi/n1
	f0 = df0
	f1 = -df1
	R2 = R0 + r*sqrt(3)
        for i in range(n1):
		f1 = f1 + df1
		V0 = Vector3(R2*cos(f0),R2*sin(f0),0)
		V1 = Vector3(R1*cos(f1),R1*sin(f1),0)
		dV = V1 - V0
		D = dV.norm()
		if D <= r*(2+sqrt(3)):
			d = (r*r+D*D)/(2*D)
			h = sqrt(4*r*r-d*d)
			vn = dV/D
			x = V0[0] + vn[0]*d + vn[1]*h*dirmot
			y = V0[1] + vn[1]*d - vn[0]*h*dirmot
			if (Vector3(x,y,0)-Vector3(R1*cos(f1+df1),R1*sin(f1+df1),0)).norm() >= r*sqrt(3):
				f0 = f0 + df0
				for j in range(Nz-1):
					O.bodies.append(utils.sphere(Vector3(x,y,2*j*r+r),r,material=sp[0],color=sp[1]))
					O.bodies[-1].state.setMinIndex = True
					trqId.extend([O.bodies[-1].id])
					trqV.extend([Vector3(0,0,trq*sp[4])])
		

Rot=Matrix3(
  cos(2*pi/3),  sin(2*pi/3),  0,
  -sin(2*pi/3),  cos(2*pi/3),  0,
  0,  0,  1
)

			
def upbracing(R0,R1,sp):
	nbra = ceil((R1-R0)/(2*r)) - 1
	cosa = 1 + (R1-R0)/(4*r) - (nbra+1)/2
	sina = sqrt(1-cosa*cosa)
	vcur = Vector3(R0+2*r*cosa,0,2*r*(Nz+sina))
	i = 0
	while i < nbra:
		O.bodies.append([utils.sphere(M*vcur,r,material=sp[0],color=sp[1]) for M in [1,Rot,Rot*Rot]])
		vcur = vcur + Vector3(2*r,0,0)
		i = i + 1


def downbracing(R0,R1,sp):
	nbra = ceil((R1-R0)/(2*r)) - 1
	cosa = 1 + (R1-R0)/(4*r) - (nbra+1)/2
	sina = sqrt(1-cosa*cosa)
	vcur = Vector3(R0+2*r*cosa,0,-2*r*(1+sina))
	i = 0
	while i < nbra:
		O.bodies.append([utils.sphere(M*vcur,r,material=sp[0],color=sp[1]) for M in [1,Rot,Rot*Rot]])
		vcur = vcur + Vector3(2*r,0,0)
		i = i + 1
		
# add central axis to the scene
for i in range(Nz+1):
	O.bodies.append(utils.sphere(Vector3(0,0,2*i*r),r,material='bluesphere',color=(0,0,1),fixed=False))


# add concentric cylinders and bracings to the scene
cylinder(6,2*r,Nz+1,-2*r,green,0,False)
cylinder(3,r*(2+sqrt(3)),Nz-1,r,green,0,False)
acylinder(3,18,2*r,r/sin(pi/18),1,red)
cylinder(18,r/sin(pi/18),Nz+1,0,blue,1000000,False)
cylinder(9,r/sin(pi/18)+r*sqrt(3),Nz-1,r,blue,0,False)
cylinder(30,r/sin(pi/30),Nz+1,-2*r,green,1000000,True)
acylinder(9,30,r/sin(pi/18),r/sin(pi/30),-1,yellow)
upbracing(0,r/sin(pi/18),blue)
downbracing(2*r,r/sin(pi/30),green)


#O.bodies.append(utils.wall(min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)]),axis=2,sense=1,material='bluesphere'))

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
	NewtonIntegrator(damping=0,gravity=(0,0,-9.81*0)),
#	PyRunner(command='addPlotData()',iterPeriod=10,label='plotter'),
]

def moment():
	for i in range(len(trqId)):
		O.forces.addT(trqId[i],trqV[i])


#set time step and run simulation:
O.dt=0.5*utils.PWaveTimeStep()



