import sys  #for interaction with OS (work with shell input variables in the following program)
from dolfin import *  #computational C++ backend of FEniCS
import numpy  
from ufl_legacy import indices  #for expressing variational problems
import json  #for working with files (save and load a list in a dir)
import os.path  #for searching in OS

#As we have more complex mathematical operations, we may need to use uflacs compiler
parameters["form_compiler"]["representation"] = "uflacs"  
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

#Set Solver Parameters
krylov_params = parameters["krylov_solver"]  #for the following solver
krylov_params["relative_tolerance"] = 1E-9
krylov_params["absolute_tolerance"] = 1E-9
krylov_params["maximum_iterations"] = 150000
krylov_params["monitor_convergence"] = False      
krylov_params["error_on_nonconvergence"] = True

newton_solver_parameters = {"nonlinear_solver": "newton",
	"newton_solver":{"linear_solver": "cg", #"preconditioner":"icc", 
						"convergence_criterion": "residual",
						"relative_tolerance": 1E-6,
						"absolute_tolerance": 1E-6, 
						"report":True,
						"krylov_solver": krylov_params, 
						"relaxation_parameter":1.0, 
						"maximum_iterations":30, 
						"error_on_nonconvergence": True
					} 
	}
	
newton_solver_parameters2 = {"nonlinear_solver": "newton",
	"newton_solver":{"linear_solver": "mumps",  
			"convergence_criterion": "incremental",
			"relative_tolerance": 1E-5,
			"absolute_tolerance": 1E-5, 
			"report":True,
			"krylov_solver": krylov_params, 
			"relaxation_parameter":1.0, 
			"maximum_iterations":30, 
			"error_on_nonconvergence": True
					} 
	}
	
snes_solver_parameters = {"nonlinear_solver": "snes",
	"snes_solver": { "line_search": "basic", 
			"linear_solver": "gmres",
			"preconditioner": "sor",
			"maximum_iterations": 100,
			"relative_tolerance": 1E-4,
			"absolute_tolerance": 1E-5,
			"report": True,
			"error_on_nonconvergence": True,
					}
	}
	
solver_parameters = newton_solver_parameters   #for phi
solver_parameters2 = newton_solver_parameters  #for psi if we wanted to change the solver
set_log_level(50) #we can choose 50, 40, 30, 20, 16, 13, or 10 for different monitoring

data = sys.argv[1]  #mesh file name in xml format
mesh = Mesh(data+".xml")  #import mesh to FEniCS
cells_mesh = MeshFunction('size_t', mesh, data+'_physical_region.xml')  #import subdomains for our two separate materials
# facets_mesh = MeshFunction('size_t', mesh, 'bigCylinderComp_facet_region.xml')  #we'll define faces in code

#from msh file
'''
------------------------------
Dim / ID / Subdomain name
------------------------------
3 2 "mat"  air     
3 1 "inc"  Cubic UC           
------------------------------
'''

Dim = mesh.topology().dim() #find out the dimension is 2D or 3D
inc, mat=1, 2  #according to msh file

CubeLength=float(sys.argv[2])
xmin, xmax = -CubeLength/2, CubeLength/2
ymin, ymax = -CubeLength/2, CubeLength/2
zmin, zmax = -CubeLength/2, CubeLength/2

null = 1E-10  #for geometry and stiffness of air
TOL = null

class PeriodicBoundary(SubDomain):
     #"target domains" G 
	def inside(self, x, on_boundary):
		return bool( on_boundary and
		( 
			(   near(x[0], xmin, TOL) and not ( near(x[1], ymax, TOL) or near(x[2], zmax, TOL) ) ) 
		or	(   near(x[1], ymin, TOL) and not ( near(x[0], xmax, TOL) or near(x[2], zmax, TOL) ) )
		or 	(   near(x[2], zmin, TOL) and not ( near(x[0], xmax, TOL) or near(x[1], ymax, TOL) ) )
		) )
		
	def map(self, x, y):	 #all of vertexes
		if near(x[0], xmin, TOL) and near(x[1], ymin, TOL) and near(x[2], zmax, TOL):
			y[0] = x[0]
			y[1] = x[1]
			y[2] = x[2] - (zmax-zmin)
		elif near(x[0], xmin, TOL) and near(x[1], ymax, TOL) and near(x[2], zmin, TOL):
			y[0] = x[0]
			y[1] = x[1] - (ymax-ymin)
			y[2] = x[2]
		elif near(x[0], xmin, TOL) and near(x[1], ymax, TOL) and near(x[2], zmax, TOL):
			y[0] = x[0]
			y[1] = x[1] - (ymax-ymin)
			y[2] = x[2] - (zmax-zmin)
		elif near(x[0], xmax, TOL) and near(x[1], ymin, TOL) and near(x[2], zmin, TOL):
			y[0] = x[0] - (xmax-xmin)
			y[1] = x[1]
			y[2] = x[2]
		elif near(x[0], xmax, TOL) and near(x[1], ymin, TOL) and near(x[2], zmax, TOL):
			y[0] = x[0] - (xmax-xmin)
			y[1] = x[1]
			y[2] = x[2] - (zmax-zmin)
		elif near(x[0], xmax, TOL) and near(x[1], ymax, TOL) and near(x[2], zmin, TOL):
			y[0] = x[0] - (xmax-xmin)
			y[1] = x[1] - (ymax-ymin)
			y[2] = x[2]
		elif near(x[0], xmax, TOL) and near(x[1], ymax, TOL) and near(x[2], zmax, TOL):
			y[0] = x[0] - (xmax-xmin)
			y[1] = x[1] - (ymax-ymin)
			y[2] = x[2] - (zmax-zmin)
		elif near(x[0], xmin, TOL) and near(x[1], ymin, TOL) and near(x[2], zmin, TOL): #
 			y[0] = x[0]
 			y[1] = x[1]
 			y[2] = x[2]
	################edge
		elif near(x[0], xmax, TOL) and near(x[2], zmin, TOL):
			y[0] = x[0] - (xmax-xmin)
			y[1] = x[1]
			y[2] = x[2]
		elif near(x[0], xmax, TOL) and near(x[2], zmax, TOL):
			y[0] = x[0] - (xmax-xmin)
			y[1] = x[1]
			y[2] = x[2] - (zmax-zmin)
		elif near(x[0], xmin, TOL) and near(x[2], zmax, TOL):
			y[0] = x[0]
			y[1] = x[1]
			y[2] = x[2] - (zmax-zmin)
		elif near(x[1], ymin, TOL) and near(x[2], zmax, TOL): #error? xmin no deferrence
			y[0] = x[0]
			y[1] = x[1]
			y[2] = x[2] - (zmax-zmin)
		elif near(x[1], ymax, TOL) and near(x[2], zmax, TOL):
			y[0] = x[0]
			y[1] = x[1] - (ymax-ymin)
			y[2] = x[2] - (zmax-zmin)
		elif near(x[1], ymax, TOL) and near(x[2], zmin, TOL):
			y[0] = x[0]
			y[1] = x[1] - (ymax-ymin)
			y[2] = x[2]
		elif near(x[0], xmax, TOL) and near(x[1], ymin, TOL):
			y[0] = x[0] - (xmax-xmin)
			y[1] = x[1]
			y[2] = x[2]
		elif near(x[0], xmax, TOL) and near(x[1], ymax, TOL):
			y[0] = x[0] - (xmax-xmin)
			y[1] = x[1] - (ymax-ymin)
			y[2] = x[2]
		elif near(x[0], xmin, TOL) and near(x[1], ymax, TOL):  #error of last authors
			y[0] = x[0]
			y[1] = x[1] - (ymax-ymin)
			y[2] = x[2]
		################surface
		elif near(x[0], xmax, TOL):
			y[0] = x[0] - (xmax-xmin)
			y[1] = x[1]
			y[2] = x[2]
		elif near(x[1], ymax, TOL):
			y[0] = x[0]
			y[1] = x[1] - (ymax-ymin)
			y[2] = x[2]
		#elif near(x[2], zmax, TOL):
		else:
			y[0] = x[0]
			y[1] = x[1]
			y[2] = x[2] - (zmax-zmin)
 	
unk_vector = VectorElement("Lagrange", mesh.ufl_cell(), 1)   
lag_vector = VectorElement("R", mesh.ufl_cell(), 0)  #space of Real numbers is for Lagrange multipliers 
Space = FunctionSpace(mesh, unk_vector*lag_vector, constrained_domain=PeriodicBoundary(TOL))
# VV = FunctionSpace(mesh, unk_vector) #for plot
# dA = Measure('ds', domain=mesh, subdomain_data=facets_mesh, metadata={'quadrature_degree': 2, "quadrature_scheme": "uflacs"})
dV = Measure('dx', domain=mesh, subdomain_data=cells_mesh, metadata={'quadrature_degree': 2, "quadrature_scheme": "uflacs"})


dunk = TrialFunction(Space)    #trials for both unknowns phi and psi with lagrange
del_unk = TestFunction(Space)
del_p, del_l = split(del_unk)  #tests for both unknowns phi and psi(del_p) and for lagrange(del_l)

unk00 = Function(Space)
phi00, lag00 = split(unk00)
unk01 = Function(Space)
phi01, lag01 = split(unk01)
unk02 = Function(Space)
phi02, lag02 = split(unk02)
#unk10 = unk01  #because of symmetry
unk11 = Function(Space)
phi11, lag11 = split(unk11)
unk12 = Function(Space)
phi12, lag12 = split(unk12)
#unk20 = unk02
#unk21 = unk12
unk22 = Function(Space)
phi22, lag22 = split(unk22)


unk000 = Function(Space)
psi000, lag000 = split(unk000)
unk001 = Function(Space)
psi001, lag001 = split(unk001)
unk002 = Function(Space)
psi002, lag002 = split(unk002)
unk010 = Function(Space)
psi010, lag010 = split(unk010)
unk011 = Function(Space)
psi011, lag011 = split(unk011)
unk012 = Function(Space)
psi012, lag012 = split(unk012)
unk020 = Function(Space)
psi020, lag020 = split(unk020)
unk021 = Function(Space)
psi021, lag021 = split(unk021)
unk022 = Function(Space)
psi022, lag022 = split(unk022)

#unk100 = unk010
#unk101 = unk011
#unk102 = unk012
unk110 = Function(Space)
psi110, lag110 = split(unk110)
unk111 = Function(Space)
psi111, lag111 = split(unk111)
unk112 = Function(Space)
psi112, lag112 = split(unk112)
unk120 = Function(Space)
psi120, lag120 = split(unk120)
unk121 = Function(Space)
psi121, lag121 = split(unk121)
unk122 = Function(Space)
psi122, lag122 = split(unk122)

#unk200 = unk020
#unk201 = unk021
#unk202 = unk022
#unk210 = unk120
#unk211 = unk121
#unk212 = unk122
unk220 = Function(Space)
psi220, lag220 = split(unk220)
unk221 = Function(Space)
psi221, lag221 = split(unk221)
unk222 = Function(Space)
psi222, lag222 = split(unk222)

#units: mm, Mg = 1000 kg=ton, s, MPa, mJ, K
i,j,k,l = indices(4)
delta = Identity(Dim)

E_mat = Constant(null) #in MPa
nu_mat = Constant(0.0)
rho_mat = Constant(0.0) #in Mg/mm3 = 1E12 kg/m3

#for Ti6Al4V ASM material data
E_inc = Constant(float(sys.argv[6])) #in MPa 
nu_inc=Constant(float(sys.argv[5]))
rho_inc = Constant(4430E-12)  #in Mg/mm3 = 1E12 kg/m3 #its value doesn't matter

#convert E,nu to la,mu
la_mat = E_mat*nu_mat/(1.+nu_mat)/(1.-2.*nu_mat)
mu_mat = E_mat/2./(1.+nu_mat)
la_inc = E_inc*nu_inc/(1.+nu_inc)/(1.-2.*nu_inc)
mu_inc = E_inc/2./(1.+nu_inc)

#isotropic in microscale
Cmat_m = as_tensor(la_mat*delta[i,j]*delta[k,l] + mu_mat*delta[i,k]*delta[j,l] + mu_mat*delta[i,l]*delta[j,k], (i,j,k,l))
Cinc_m = as_tensor(la_inc*delta[i,j]*delta[k,l] + mu_inc*delta[i,k]*delta[j,l] + mu_inc*delta[i,l]*delta[j,k], (i,j,k,l))

dofs = Space.dim()
Vol = assemble(1.*dV)
x = SpatialCoordinate(mesh) #for M and I later

L00 = as_tensor( delta[i,0]*delta[j,0] + phi00[i].dx(j) , (i,j))
L11 = as_tensor( delta[i,1]*delta[j,1] + phi11[i].dx(j) , (i,j))
L22 = as_tensor( delta[i,2]*delta[j,2] + phi22[i].dx(j) , (i,j))
L12 = as_tensor( delta[i,1]*delta[j,2] + phi12[i].dx(j) , (i,j))
L21 = L12
L02 = as_tensor( delta[i,0]*delta[j,2] + phi02[i].dx(j) , (i,j))
L20 = L02
L01 = as_tensor( delta[i,0]*delta[j,1] + phi01[i].dx(j) , (i,j))
L10 = L01


#Form_phiab = Cmat_m[i,j,k,l]*Lab[k,l]*del_p[i].dx(j)*dV(mat) + Cinc_m[i,j,k,l]*Lab[k,l]*del_p[i].dx(j)*dV(inc)

Form_phi00 = Cmat_m[i,j,k,l]*L00[k,l]*del_p[i].dx(j)*dV(mat) + Cinc_m[i,j,k,l]*L00[k,l]*del_p[i].dx(j)*dV(inc) + (phi00[i]*del_l[i] + del_p[i]*lag00[i])*dV
Form_phi11 = Cmat_m[i,j,k,l]*L11[k,l]*del_p[i].dx(j)*dV(mat) + Cinc_m[i,j,k,l]*L11[k,l]*del_p[i].dx(j)*dV(inc) + (phi11[i]*del_l[i] + del_p[i]*lag11[i])*dV
Form_phi22 = Cmat_m[i,j,k,l]*L22[k,l]*del_p[i].dx(j)*dV(mat) + Cinc_m[i,j,k,l]*L22[k,l]*del_p[i].dx(j)*dV(inc) + (phi22[i]*del_l[i] + del_p[i]*lag22[i])*dV
Form_phi12 = Cmat_m[i,j,k,l]*L12[k,l]*del_p[i].dx(j)*dV(mat) + Cinc_m[i,j,k,l]*L12[k,l]*del_p[i].dx(j)*dV(inc) + (phi12[i]*del_l[i] + del_p[i]*lag12[i])*dV
Form_phi02 = Cmat_m[i,j,k,l]*L02[k,l]*del_p[i].dx(j)*dV(mat) + Cinc_m[i,j,k,l]*L02[k,l]*del_p[i].dx(j)*dV(inc) + (phi02[i]*del_l[i] + del_p[i]*lag02[i])*dV
Form_phi01 = Cmat_m[i,j,k,l]*L01[k,l]*del_p[i].dx(j)*dV(mat) + Cinc_m[i,j,k,l]*L01[k,l]*del_p[i].dx(j)*dV(inc) + (phi01[i]*del_l[i] + del_p[i]*lag01[i])*dV

if True: 
# 	elapsed = int(comp_time.time() - start_time)
# 	e_h, e_m, e_s = int(elapsed/3600), int(elapsed % 3600 / 60), int( ( elapsed % 3600 ) % 60 )
# 	print('%.0f DOFs in %.0f h %.0f min %.0f s' % (dofs, e_h, e_m, e_s)  )
    print('%.0f DOFs' % (dofs))
    
print('Solving phi 00')

Gain = derivative(Form_phi00, unk00, dunk)  
problem = NonlinearVariationalProblem(Form_phi00, unk00, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters)
(iternr, converged) = solver.solve()

# file_phi00 = File('phi00_No_B.pvd')  
# (phi00, lag) = split(unk00)
# phi00 = project (phi00, VV, solver_type='bicgstab', preconditioner_type='hypre_amg')
# file_phi00 << phi00

print('Solving phi 11')

Gain = derivative(Form_phi11, unk11, dunk)
problem = NonlinearVariationalProblem(Form_phi11, unk11, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters)
(iternr, converged) = solver.solve()

# file_Psi11 = File('Psi000_No_B.pvd')
# (PPsi11, lag) = split(unk11)
# Psi11 = project (PPsi11, VV, solver_type='bicgstab', preconditioner_type='hypre_amg')
# file_Psi11 << Psi11

print('Solving phi 22')

Gain = derivative(Form_phi22, unk22, dunk)
problem = NonlinearVariationalProblem(Form_phi22, unk22, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters)
(iternr, converged) = solver.solve()

print('Solving phi 12')

Gain = derivative(Form_phi12, unk12, dunk)
problem = NonlinearVariationalProblem(Form_phi12, unk12, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters)
(iternr, converged) = solver.solve()

print('Solving phi 02')

Gain = derivative(Form_phi02, unk02, dunk)
problem = NonlinearVariationalProblem(Form_phi02, unk02, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters)
(iternr, converged) = solver.solve()

print('Solving phi 01')

Gain = derivative(Form_phi01, unk01, dunk)
problem = NonlinearVariationalProblem(Form_phi01, unk01, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters)
(iternr, converged) = solver.solve()

print('Determining parameters of C') #21 #just intervene symmetry of C (major symmetry)

#Cabcd = 1./Vol * assemble( Cmat_m[i,j,k,l]*Lab[i,j]*Lcd[k,l]*dV(mat) )
C0000 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*L00[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*L00[k,l]*dV(inc) )
C0011 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*L11[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*L11[k,l]*dV(inc) )
C0022 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*L22[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*L22[k,l]*dV(inc) )
C0012 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*L12[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*L12[k,l]*dV(inc) )
C0002 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*L02[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*L02[k,l]*dV(inc) )
C0001 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*L01[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*L01[k,l]*dV(inc) )

C1100 = C0011
C1111 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*L11[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*L11[k,l]*dV(inc) )
C1122 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*L22[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*L22[k,l]*dV(inc) )
C1112 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*L12[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*L12[k,l]*dV(inc) )
C1102 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*L02[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*L02[k,l]*dV(inc) )
C1101 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*L01[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*L01[k,l]*dV(inc) )

C2200 = C0022
C2211 = C1122
C2222 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*L22[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*L22[k,l]*dV(inc) )
C2212 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*L12[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*L12[k,l]*dV(inc) )
C2202 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*L02[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*L02[k,l]*dV(inc) )
C2201 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*L01[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*L01[k,l]*dV(inc) )

C1200 = C0012
C1211 = C1112
C1222 = C2212
C1212 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*L12[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*L12[k,l]*dV(inc) )
C1202 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*L02[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*L02[k,l]*dV(inc) )
C1201 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*L01[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*L01[k,l]*dV(inc) )

C0200 = C0002
C0211 = C1102
C0222 = C2202
C0212 = C1202
C0202 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*L02[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*L02[k,l]*dV(inc) )
C0201 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*L01[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*L01[k,l]*dV(inc) )

C0100 = C0001
C0111 = C1101
C0122 = C2201
C0112 = C1201
C0102 = C0201
C0101 = 1./Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*L01[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*L01[k,l]*dV(inc) )

Sf = 8 #for round
divd = 1
if True: 
	with open(data+'.tex', 'w') as file:
		file.write(r'\documentclass[a4paper,10pt,abstracton,openbib,final]{scrartcl}'+'\n')
		file.write(r'\usepackage{amsmath, amsthm, amsxtra, amsfonts, amssymb, mathtools}'+'\n')
		file.write(r'\setcounter{MaxMatrixCols}{20}'+'\n')
		file.write(r'\newcommand{\begeq}{\begin{equation}\begin{gathered}}'+'\n')
		file.write(r'\newcommand{\eqend}{\end{gathered}\end{equation}}'+'\n')
		file.write(r'\newcommand{\mi}{^\text{m}}'+'\n')
		file.write(r'\newcommand{\ma}{^\text{M}}'+'\n')
		file.write(r'\begin{document}'+'\n')


	with open(data+'.tex', 'a') as file:
		file.write(r'\begeq')
		file.write('\n')
		file.write(r'C\ma_{AB} = \begin{pmatrix}')
		file.write('\n')
		file.write(str(round(C0000/divd,Sf))+' & '+str(round(C0011/divd,Sf))+' & '+str(round(C0022/divd,Sf))+' & '+str(round(C0012/divd,Sf))+' & '+str(round(C0002/divd,Sf))+' & '+str(round(C0001/divd,Sf))+r' \\ ')
		file.write('\n')
		file.write(str(round(C1100/divd,Sf))+' & '+str(round(C1111/divd,Sf))+' & '+str(round(C1122/divd,Sf))+' & '+str(round(C1112/divd,Sf))+' & '+str(round(C1102/divd,Sf))+' & '+str(round(C1101/divd,Sf))+r' \\ ') 
		file.write('\n')
		file.write(str(round(C2200/divd,Sf))+' & '+str(round(C2211/divd,Sf))+' & '+str(round(C2222/divd,Sf))+' & '+str(round(C2212/divd,Sf))+' & '+str(round(C2202/divd,Sf))+' & '+str(round(C2201/divd,Sf))+r' \\ ') 
		file.write('\n')
		file.write(str(round(C1200/divd,Sf))+' & '+str(round(C1211/divd,Sf))+' & '+str(round(C1222/divd,Sf))+' & '+str(round(C1212/divd,Sf))+' & '+str(round(C1202/divd,Sf))+' & '+str(round(C1201/divd,Sf))+r' \\ ')
		file.write('\n')
		file.write(str(round(C0200/divd,Sf))+' & '+str(round(C0211/divd,Sf))+' & '+str(round(C0222/divd,Sf))+' & '+str(round(C0212/divd,Sf))+' & '+str(round(C0202/divd,Sf))+' & '+str(round(C0201/divd,Sf))+r' \\ ')
		file.write('\n')
		file.write(str(round(C0100/divd,Sf))+' & '+str(round(C0111/divd,Sf))+' & '+str(round(C0122/divd,Sf))+' & '+str(round(C0112/divd,Sf))+' & '+str(round(C0102/divd,Sf))+' & '+str(round(C0101/divd,Sf))+r' \\ ')
		file.write('\n')
		file.write(r'\end{pmatrix} \text{\,MPa} \ ,')
		file.write('\n')
		file.write(r'\eqend')
		file.write('\n')



def VoigtToTensor(A):
	A11, A12, A13, A14, A15, A16 = A[0,0], A[0,1], A[0,2], A[0,3], A[0,4], A[0,5]
	A22, A23, A24, A25, A26 = A[1,1], A[1,2], A[1,3], A[1,4], A[1,5]
	A33, A34, A35, A36 = A[2,2], A[2,3], A[2,4], A[2,5]
	A44, A45, A46 = A[3,3], A[3,4], A[3,5]
	A55, A56 = A[4,4], A[4,5]
	A66 = A[5,5]
	A21, A31, A41, A51, A61 = A12, A13, A14, A15, A16
	A32, A42, A52, A62 = A23, A24, A25, A26
	A43, A53, A63 = A34, A35, A36
	A54, A64 = A45, A46
	A65 = A56
	return as_tensor([ \
	[ \
	[ [A11,A16,A15], [A16,A12,A14], [A15,A14,A13]] , \
	[ [A61,A66,A65], [A66,A62,A64], [A65,A64,A63]] , \
	[ [A51,A56,A55], [A56,A52,A54], [A55,A54,A53]] \
	] , [ \
	[ [A61,A66,A65], [A66,A62,A64], [A65,A64,A63]] , \
	[ [A21,A26,A25], [A26,A22,A24], [A25,A24,A23]] , \
	[ [A41,A46,A45], [A46,A42,A44], [A45,A44,A43]] \
	] , [ \
	[ [A51,A56,A55], [A56,A52,A54], [A55,A54,A53]] , \
	[ [A41,A46,A45], [A46,A42,A44], [A45,A44,A43]] , \
	[ [A31,A36,A35], [A36,A32,A34], [A35,A34,A33]] ] \
	])

C_Voigt = numpy.array([ \
[C0000, C0011, C0022, C0012, C0002, C0001],\
[C1100, C1111, C1122, C1112, C1102, C1101],\
[C2200, C2211, C2222, C2212, C2202, C2201],\
[C1200, C1211, C1222, C1212, C1202, C1201],\
[C0200, C0211, C0222, C0212, C0202, C0201],\
[C0100, C0111, C0122, C0112, C0102, C0101]  ])

C_M = VoigtToTensor(C_Voigt)

N000 = as_tensor( phi00[i]*delta[j,0] + psi000[i].dx(j) , (i,j))
N001 = as_tensor( phi00[i]*delta[j,1] + psi001[i].dx(j) , (i,j))
N002 = as_tensor( phi00[i]*delta[j,2] + psi002[i].dx(j) , (i,j))
N010 = as_tensor( phi01[i]*delta[j,0] + psi010[i].dx(j) , (i,j))
N011 = as_tensor( phi01[i]*delta[j,1] + psi011[i].dx(j) , (i,j))
N012 = as_tensor( phi01[i]*delta[j,2] + psi012[i].dx(j) , (i,j))
N020 = as_tensor( phi02[i]*delta[j,0] + psi020[i].dx(j) , (i,j))
N021 = as_tensor( phi02[i]*delta[j,1] + psi021[i].dx(j) , (i,j))
N022 = as_tensor( phi02[i]*delta[j,2] + psi022[i].dx(j) , (i,j))

N100 = N010 
N101 = N011 
N102 = N012 

N110 = as_tensor( phi11[i]*delta[j,0] + psi110[i].dx(j) , (i,j))
N111 = as_tensor( phi11[i]*delta[j,1] + psi111[i].dx(j) , (i,j))
N112 = as_tensor( phi11[i]*delta[j,2] + psi112[i].dx(j) , (i,j))

N120 = as_tensor( phi12[i]*delta[j,0] + psi120[i].dx(j) , (i,j))
N121 = as_tensor( phi12[i]*delta[j,1] + psi121[i].dx(j) , (i,j))
N122 = as_tensor( phi12[i]*delta[j,2] + psi122[i].dx(j) , (i,j))

N200 = N020 
N201 = N021 
N202 = N022 

N210 = N120 
N211 = N121 
N212 = N122 

N220 = as_tensor( phi22[i]*delta[j,0] + psi220[i].dx(j) , (i,j))
N221 = as_tensor( phi22[i]*delta[j,1] + psi221[i].dx(j) , (i,j))
N222 = as_tensor( phi22[i]*delta[j,2] + psi222[i].dx(j) , (i,j))

rho_M = assemble(rho_mat*dV(mat) + rho_inc*dV(inc))
rho_M = rho_M/assemble(1*dV)

#Form_psiabc = ( Cmat_m[i,j,k,l]*(psiabc[k].dx(l)+phiab[k]*delta[l,c])*del_p[i].dx(j) - (phiab[k].dx(l)+delta[k,a]*delta[l,b])*Cmat_m[i,c,k,l]*del_p[i]  )*dV(mat) + ( Cinc_m[i,j,k,l]*(psiabc[k].dx(l)+phiab[k]*delta[l,c])*del_p[i].dx(j) - (phiab[k].dx(l)+delta[k,a]*delta[l,b])*Cinc_m[i,c,k,l]*del_p[i] )*dV(inc) + C_M[a,b,i,c]*del_p[i]*dV

Form_psi000 = (- rho_mat/rho_M *C_M[0,0,0,i]*del_p[i] +  L00[k,l]*Cmat_m[k,l,0,i]*del_p[i] - Cmat_m[i,j,k,l]*N000[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M * C_M[0,0,0,i]*del_p[i]  + L00[k,l]*Cinc_m[k,l,0,i]*del_p[i] - Cinc_m[i,j,k,l]*N000[k,l]*del_p[i].dx(j) )*dV(inc) + (psi000[i]*del_l[i] + del_p[i]*lag000[i])*dV

Form_psi001 = (- rho_mat/rho_M *C_M[0,0,1,i]*del_p[i] +  L00[k,l]*Cmat_m[k,l,1,i]*del_p[i] - Cmat_m[i,j,k,l]*N001[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M *C_M[0,0,1,i]*del_p[i] + L00[k,l]*Cinc_m[k,l,1,i]*del_p[i] - Cinc_m[i,j,k,l]*N001[k,l]*del_p[i].dx(j) )*dV(inc) + (psi001[i]*del_l[i] + del_p[i]*lag001[i])*dV  

Form_psi002 = (- rho_mat/rho_M *C_M[0,0,2,i]*del_p[i] +  L00[k,l]*Cmat_m[k,l,2,i]*del_p[i] - Cmat_m[i,j,k,l]*N002[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M *C_M[0,0,2,i]*del_p[i] + L00[k,l]*Cinc_m[k,l,2,i]*del_p[i] - Cinc_m[i,j,k,l]*N002[k,l]*del_p[i].dx(j) )*dV(inc) + (psi002[i]*del_l[i] + del_p[i]*lag002[i])*dV

Form_psi010 = (- rho_mat/rho_M *C_M[0,1,0,i]*del_p[i] +  L01[k,l]*Cmat_m[k,l,0,i]*del_p[i] - Cmat_m[i,j,k,l]*N010[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M * C_M[0,1,0,i]*del_p[i]  + L01[k,l]*Cinc_m[k,l,0,i]*del_p[i] - Cinc_m[i,j,k,l]*N010[k,l]*del_p[i].dx(j) )*dV(inc) + (psi010[i]*del_l[i] + del_p[i]*lag010[i])*dV

Form_psi011 = (- rho_mat/rho_M *C_M[0,1,1,i]*del_p[i] +  L01[k,l]*Cmat_m[k,l,1,i]*del_p[i] - Cmat_m[i,j,k,l]*N011[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M * C_M[0,1,1,i]*del_p[i]  + L01[k,l]*Cinc_m[k,l,1,i]*del_p[i] - Cinc_m[i,j,k,l]*N011[k,l]*del_p[i].dx(j) )*dV(inc) + (psi011[i]*del_l[i] + del_p[i]*lag011[i])*dV

Form_psi012 = (- rho_mat/rho_M *C_M[0,1,2,i]*del_p[i] +  L01[k,l]*Cmat_m[k,l,2,i]*del_p[i] - Cmat_m[i,j,k,l]*N012[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M * C_M[0,1,2,i]*del_p[i]  + L01[k,l]*Cinc_m[k,l,2,i]*del_p[i] - Cinc_m[i,j,k,l]*N012[k,l]*del_p[i].dx(j) )*dV(inc) + (psi012[i]*del_l[i] + del_p[i]*lag012[i])*dV

Form_psi020 = (- rho_mat/rho_M *C_M[0,2,0,i]*del_p[i] +  L02[k,l]*Cmat_m[k,l,0,i]*del_p[i] - Cmat_m[i,j,k,l]*N020[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M * C_M[0,2,0,i]*del_p[i]  + L02[k,l]*Cinc_m[k,l,0,i]*del_p[i] - Cinc_m[i,j,k,l]*N020[k,l]*del_p[i].dx(j) )*dV(inc) + (psi020[i]*del_l[i] + del_p[i]*lag020[i])*dV

Form_psi021 = (- rho_mat/rho_M *C_M[0,2,1,i]*del_p[i] +  L02[k,l]*Cmat_m[k,l,1,i]*del_p[i] - Cmat_m[i,j,k,l]*N021[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M * C_M[0,2,1,i]*del_p[i]  + L02[k,l]*Cinc_m[k,l,1,i]*del_p[i] - Cinc_m[i,j,k,l]*N021[k,l]*del_p[i].dx(j) )*dV(inc) + (psi021[i]*del_l[i] + del_p[i]*lag021[i])*dV

Form_psi022 = (- rho_mat/rho_M *C_M[0,2,2,i]*del_p[i] +  L02[k,l]*Cmat_m[k,l,2,i]*del_p[i] - Cmat_m[i,j,k,l]*N022[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M * C_M[0,2,2,i]*del_p[i]  + L02[k,l]*Cinc_m[k,l,2,i]*del_p[i] - Cinc_m[i,j,k,l]*N022[k,l]*del_p[i].dx(j) )*dV(inc) + (psi022[i]*del_l[i] + del_p[i]*lag022[i])*dV

Form_psi100 = Form_psi010 
Form_psi101 = Form_psi011 
Form_psi102 = Form_psi012 

Form_psi110 = (- rho_mat/rho_M *C_M[1,1,0,i]*del_p[i] +  L11[k,l]*Cmat_m[k,l,0,i]*del_p[i] - Cmat_m[i,j,k,l]*N110[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M * C_M[1,1,0,i]*del_p[i]  + L11[k,l]*Cinc_m[k,l,0,i]*del_p[i] - Cinc_m[i,j,k,l]*N110[k,l]*del_p[i].dx(j) )*dV(inc) + (psi110[i]*del_l[i] + del_p[i]*lag110[i])*dV

Form_psi111 = (- rho_mat/rho_M *C_M[1,1,1,i]*del_p[i] +  L11[k,l]*Cmat_m[k,l,1,i]*del_p[i] - Cmat_m[i,j,k,l]*N111[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M * C_M[1,1,1,i]*del_p[i]  + L11[k,l]*Cinc_m[k,l,1,i]*del_p[i] - Cinc_m[i,j,k,l]*N111[k,l]*del_p[i].dx(j) )*dV(inc) + (psi111[i]*del_l[i] + del_p[i]*lag111[i])*dV

Form_psi112 = (- rho_mat/rho_M *C_M[1,1,2,i]*del_p[i] +  L11[k,l]*Cmat_m[k,l,2,i]*del_p[i] - Cmat_m[i,j,k,l]*N112[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M * C_M[1,1,2,i]*del_p[i]  + L11[k,l]*Cinc_m[k,l,2,i]*del_p[i] - Cinc_m[i,j,k,l]*N112[k,l]*del_p[i].dx(j) )*dV(inc) + (psi112[i]*del_l[i] + del_p[i]*lag112[i])*dV

Form_psi120 = (- rho_mat/rho_M *C_M[1,2,0,i]*del_p[i] +  L12[k,l]*Cmat_m[k,l,0,i]*del_p[i] - Cmat_m[i,j,k,l]*N120[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M * C_M[1,2,0,i]*del_p[i]  + L12[k,l]*Cinc_m[k,l,0,i]*del_p[i] - Cinc_m[i,j,k,l]*N120[k,l]*del_p[i].dx(j) )*dV(inc) + (psi120[i]*del_l[i] + del_p[i]*lag120[i])*dV

Form_psi121 = (- rho_mat/rho_M *C_M[1,2,1,i]*del_p[i] +  L12[k,l]*Cmat_m[k,l,1,i]*del_p[i] - Cmat_m[i,j,k,l]*N121[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M * C_M[1,2,1,i]*del_p[i]  + L12[k,l]*Cinc_m[k,l,1,i]*del_p[i] - Cinc_m[i,j,k,l]*N121[k,l]*del_p[i].dx(j) )*dV(inc) + (psi121[i]*del_l[i] + del_p[i]*lag121[i])*dV

Form_psi122 = (- rho_mat/rho_M *C_M[1,2,2,i]*del_p[i] +  L12[k,l]*Cmat_m[k,l,2,i]*del_p[i] - Cmat_m[i,j,k,l]*N122[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M * C_M[1,2,2,i]*del_p[i]  + L12[k,l]*Cinc_m[k,l,2,i]*del_p[i] - Cinc_m[i,j,k,l]*N122[k,l]*del_p[i].dx(j) )*dV(inc) + (psi122[i]*del_l[i] + del_p[i]*lag122[i])*dV

Form_psi200 = Form_psi020 
Form_psi201 = Form_psi021 
Form_psi202 = Form_psi022 
Form_psi210 = Form_psi120 
Form_psi211 = Form_psi121 
Form_psi212 = Form_psi122 

Form_psi220 = (- rho_mat/rho_M *C_M[2,2,0,i]*del_p[i] +  L22[k,l]*Cmat_m[k,l,0,i]*del_p[i] - Cmat_m[i,j,k,l]*N220[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M * C_M[2,2,0,i]*del_p[i]  + L22[k,l]*Cinc_m[k,l,0,i]*del_p[i] - Cinc_m[i,j,k,l]*N220[k,l]*del_p[i].dx(j) )*dV(inc) + (psi220[i]*del_l[i] + del_p[i]*lag220[i])*dV

Form_psi221 = (- rho_mat/rho_M *C_M[2,2,1,i]*del_p[i] +  L22[k,l]*Cmat_m[k,l,1,i]*del_p[i] - Cmat_m[i,j,k,l]*N221[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M * C_M[2,2,1,i]*del_p[i]  + L22[k,l]*Cinc_m[k,l,1,i]*del_p[i] - Cinc_m[i,j,k,l]*N221[k,l]*del_p[i].dx(j) )*dV(inc) + (psi221[i]*del_l[i] + del_p[i]*lag221[i])*dV

Form_psi222 = (- rho_mat/rho_M *C_M[2,2,2,i]*del_p[i] +  L22[k,l]*Cmat_m[k,l,2,i]*del_p[i] - Cmat_m[i,j,k,l]*N222[k,l]*del_p[i].dx(j) )*dV(mat) + ( - rho_inc/rho_M * C_M[2,2,2,i]*del_p[i]  + L22[k,l]*Cinc_m[k,l,2,i]*del_p[i] - Cinc_m[i,j,k,l]*N222[k,l]*del_p[i].dx(j) )*dV(inc) + (psi222[i]*del_l[i] + del_p[i]*lag222[i])*dV


print('Solving psi 00i')

Gain = derivative(Form_psi000, unk000, dunk)
problem = NonlinearVariationalProblem(Form_psi000, unk000, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

# file_Psi000 = File('Psi000_No_B.pvd')  
# (PPsi000, lag) = split(unk000)
# Psi000 = project (PPsi000, VV, solver_type='bicgstab', preconditioner_type='hypre_amg')
# file_Psi000 << Psi000


Gain = derivative(Form_psi001, unk001, dunk)
problem = NonlinearVariationalProblem(Form_psi001, unk001, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

Gain = derivative(Form_psi002, unk002, dunk)
problem = NonlinearVariationalProblem(Form_psi002, unk002, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

# if processID == 0: 
# 	elapsed = int(comp_time.time() - start_time)
# 	e_h, e_m, e_s = int(elapsed/3600), int(elapsed % 3600 / 60), int( ( elapsed % 3600 ) % 60 )
# 	print('%.0f DOFs in %.0f h %.0f min %.0f s' % (dofs, e_h, e_m, e_s)  )

print('Solving psi 01i')

Gain = derivative(Form_psi010, unk010, dunk)
problem = NonlinearVariationalProblem(Form_psi010, unk010, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

Gain = derivative(Form_psi011, unk011, dunk)
problem = NonlinearVariationalProblem(Form_psi011, unk011, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

Gain = derivative(Form_psi012, unk012, dunk)
problem = NonlinearVariationalProblem(Form_psi012, unk012, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

# if processID == 0: 
# 	elapsed = int(comp_time.time() - start_time)
# 	e_h, e_m, e_s = int(elapsed/3600), int(elapsed % 3600 / 60), int( ( elapsed % 3600 ) % 60 )
# 	print('%.0f DOFs in %.0f h %.0f min %.0f s' % (dofs, e_h, e_m, e_s)  )

print('Solving psi 02i')

Gain = derivative(Form_psi020, unk020, dunk)
problem = NonlinearVariationalProblem(Form_psi020, unk020, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

Gain = derivative(Form_psi021, unk021, dunk)
problem = NonlinearVariationalProblem(Form_psi021, unk021, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

Gain = derivative(Form_psi022, unk022, dunk)
problem = NonlinearVariationalProblem(Form_psi022, unk022, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

# if processID == 0: 
# 	elapsed = int(comp_time.time() - start_time)
# 	e_h, e_m, e_s = int(elapsed/3600), int(elapsed % 3600 / 60), int( ( elapsed % 3600 ) % 60 )
# 	print('%.0f DOFs in %.0f h %.0f min %.0f s' % (dofs, e_h, e_m, e_s)  )

# print('Skipping, psi 10i = psi 01i ')

print('Solving psi 11i')

Gain = derivative(Form_psi110, unk110, dunk)
problem = NonlinearVariationalProblem(Form_psi110, unk110, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

Gain = derivative(Form_psi111, unk111, dunk)
problem = NonlinearVariationalProblem(Form_psi111, unk111, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

Gain = derivative(Form_psi112, unk112, dunk)
problem = NonlinearVariationalProblem(Form_psi112, unk112, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

# if processID == 0: 
# 	elapsed = int(comp_time.time() - start_time)
# 	e_h, e_m, e_s = int(elapsed/3600), int(elapsed % 3600 / 60), int( ( elapsed % 3600 ) % 60 )
# 	print('%.0f DOFs in %.0f h %.0f min %.0f s' % (dofs, e_h, e_m, e_s)  )

print('Solving psi 12i')

Gain = derivative(Form_psi120, unk120, dunk)
problem = NonlinearVariationalProblem(Form_psi120, unk120, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

Gain = derivative(Form_psi121, unk121, dunk)
problem = NonlinearVariationalProblem(Form_psi121, unk121, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

Gain = derivative(Form_psi122, unk122, dunk)
problem = NonlinearVariationalProblem(Form_psi122, unk122, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

# if processID == 0: 
# 	elapsed = int(comp_time.time() - start_time)
# 	e_h, e_m, e_s = int(elapsed/3600), int(elapsed % 3600 / 60), int( ( elapsed % 3600 ) % 60 )
# 	print('%.0f DOFs in %.0f h %.0f min %.0f s' % (dofs, e_h, e_m, e_s)  )

# print('Skipping, psi 20i = psi 02i')

# print('Skipping, psi 21i = psi 12i')

print('Solving psi 22i')

Gain = derivative(Form_psi220, unk220, dunk)
problem = NonlinearVariationalProblem(Form_psi220, unk220, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

Gain = derivative(Form_psi221, unk221, dunk)
problem = NonlinearVariationalProblem(Form_psi221, unk221, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

Gain = derivative(Form_psi222, unk222, dunk)
problem = NonlinearVariationalProblem(Form_psi222, unk222, bcs=[], J=Gain)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters2)
(iternr, converged) = solver.solve()

   #                    unk11.vector().max(),unk22.vector().max(),unk12.vector().max()
   # ,unk02.vector().max(),unk01.vector().max(), unk000.vector().max(), unk001.vector().max()
   # , unk002.vector().max(), unk010.vector().max(), unk011.vector().max(), unk012.vector().max()
   # , unk020.vector().max(), unk021.vector().max(), unk022.vector().max(), unk110.vector().max()
   # , unk111.vector().max(), unk112.vector().max(), unk120.vector().max(), unk121.vector().max()
   # , unk122.vector().max(), unk220.vector().max(), unk221.vector().max(), unk222.vector().m
   
# if processID == 0: 
# 	elapsed = int(comp_time.time() - start_time)
# 	e_h, e_m, e_s = int(elapsed/3600), int(elapsed % 3600 / 60), int( ( elapsed % 3600 ) % 60 )
# 	print('%.0f DOFs in %.0f h %.0f min %.0f s' % (dofs, e_h, e_m, e_s)  )

# CurrentMaxValuesList=[max(abs(unk00.vector().max()),abs(unk00.vector().min()))
#                       , max(abs(unk11.vector().max()),abs(unk11.vector().min()))
#                       , max(abs(unk22.vector().max()),abs(unk22.vector().min()))
#                       , max(abs(unk12.vector().max()),abs(unk12.vector().min()))
#                       , max(abs(unk02.vector().max()),abs(unk02.vector().min()))
#                       , max(abs(unk01.vector().max()),abs(unk01.vector().min()))
#                       , max(abs(unk000.vector().max()),abs(unk000.vector().min()))
#                       , max(abs(unk001.vector().max()),abs(unk001.vector().min()))
#                       , max(abs(unk002.vector().max()),abs(unk002.vector().min()))
#                       , max(abs(unk010.vector().max()),abs(unk010.vector().min()))
#                       , max(abs(unk011.vector().max()),abs(unk011.vector().min()))
#                       , max(abs(unk012.vector().max()),abs(unk012.vector().min()))
#                       , max(abs(unk020.vector().max()),abs(unk020.vector().min()))
#                       , max(abs(unk021.vector().max()),abs(unk021.vector().min()))
#                       , max(abs(unk022.vector().max()),abs(unk022.vector().min()))
#                       , max(abs(unk110.vector().max()),abs(unk110.vector().min()))
#                       , max(abs(unk111.vector().max()),abs(unk111.vector().min()))
#                       , max(abs(unk112.vector().max()),abs(unk112.vector().min()))
#                       , max(abs(unk120.vector().max()),abs(unk120.vector().min()))
#                       , max(abs(unk121.vector().max()),abs(unk121.vector().min()))
#                       , max(abs(unk122.vector().max()),abs(unk122.vector().min()))
#                       , max(abs(unk220.vector().max()),abs(unk220.vector().min()))
#                       , max(abs(unk221.vector().max()),abs(unk221.vector().min()))
#                       , max(abs(unk222.vector().max()),abs(unk222.vector().min()))]

# if os.path.isfile("ForCheckingConvergence.txt"): #not in first time
#     with open("ForCheckingConvergence.txt", "r") as fp:         
#         LastMaxValuesList = json.load(fp)
        
# with open("ForCheckingConvergence.txt", "w") as fp:  
#     json.dump(CurrentMaxValuesList, fp)
    
# if 'LastMaxValuesList' in globals(): #not in first time
#     for indexes in range(0,24):  
#         print('for', indexes+1, 'last is ', LastMaxValuesList[indexes], ' and current is ',CurrentMaxValuesList[indexes])
#         print(indexes+1 ,' = ', abs((1-(CurrentMaxValuesList[indexes]/LastMaxValuesList[indexes]))*100))
#     for indexes in range(0,24):  
#         if abs((1-(CurrentMaxValuesList[indexes]/LastMaxValuesList[indexes]))*100) >= 10:
#             print(indexes+1 , 'is not good')
#             conv = 0
#             sys.exit(conv)
# else:   #for first time
#     conv = 0
#     sys.exit(conv)
   
# conv = 1
  
M000 = as_tensor( L00[i,j]*x[0] + N000[i,j] , (i,j))
M001 = as_tensor( L00[i,j]*x[1] + N001[i,j] , (i,j))
M002 = as_tensor( L00[i,j]*x[2] + N002[i,j] , (i,j))
M010 = as_tensor( L01[i,j]*x[0] + N010[i,j] , (i,j))
M011 = as_tensor( L01[i,j]*x[1] + N011[i,j] , (i,j))
M012 = as_tensor( L01[i,j]*x[2] + N012[i,j] , (i,j))
M020 = as_tensor( L02[i,j]*x[0] + N020[i,j] , (i,j))
M021 = as_tensor( L02[i,j]*x[1] + N021[i,j] , (i,j))
M022 = as_tensor( L02[i,j]*x[2] + N022[i,j] , (i,j))

M100 = M010
M101 = M011
M102 = M012
M110 = as_tensor( L11[i,j]*x[0] + N110[i,j] , (i,j))
M111 = as_tensor( L11[i,j]*x[1] + N111[i,j] , (i,j))
M112 = as_tensor( L11[i,j]*x[2] + N112[i,j] , (i,j))
M120 = as_tensor( L12[i,j]*x[0] + N120[i,j] , (i,j))
M121 = as_tensor( L12[i,j]*x[1] + N121[i,j] , (i,j))
M122 = as_tensor( L12[i,j]*x[2] + N122[i,j] , (i,j))

M200 = M020
M201 = M021
M202 = M022
M210 = M120
M211 = M121
M212 = M122
M220 = as_tensor( L22[i,j]*x[0] + N220[i,j] , (i,j))
M221 = as_tensor( L22[i,j]*x[1] + N221[i,j] , (i,j))
M222 = as_tensor( L22[i,j]*x[2] + N222[i,j] , (i,j))

epsilon = 1.0 #homothetic ratio, corresponding to the real size of the substructure

print('Determining parameters of G')

#Gabcde = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*Lab[i,j]*Mcde[k,l]*dV(mat) )
G00000 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M000[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M000[k,l]*dV(inc))
G11000 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M000[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M000[k,l]*dV(inc))
G22000 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M000[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M000[k,l]*dV(inc))
G12000 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M000[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M000[k,l]*dV(inc))
G02000 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M000[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M000[k,l]*dV(inc))
G01000 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M000[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M000[k,l]*dV(inc))

G00110 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M110[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M110[k,l]*dV(inc))
G11110 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M110[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M110[k,l]*dV(inc))
G22110 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M110[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M110[k,l]*dV(inc))
G12110 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M110[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M110[k,l]*dV(inc))
G02110 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M110[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M110[k,l]*dV(inc))
G01110 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M110[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M110[k,l]*dV(inc))

G00220 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M220[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M220[k,l]*dV(inc))
G11220 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M220[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M220[k,l]*dV(inc))
G22220 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M220[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M220[k,l]*dV(inc))
G12220 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M220[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M220[k,l]*dV(inc))
G02220 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M220[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M220[k,l]*dV(inc))
G01220 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M220[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M220[k,l]*dV(inc))

G00120 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M120[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M120[k,l]*dV(inc))
G11120 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M120[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M120[k,l]*dV(inc))
G22120 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M120[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M120[k,l]*dV(inc))
G12120 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M120[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M120[k,l]*dV(inc))
G02120 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M120[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M120[k,l]*dV(inc))
G01120 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M120[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M120[k,l]*dV(inc))

G00020 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M020[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M020[k,l]*dV(inc))
G11020 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M020[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M020[k,l]*dV(inc))
G22020 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M020[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M020[k,l]*dV(inc))
G12020 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M020[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M020[k,l]*dV(inc))
G02020 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M020[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M020[k,l]*dV(inc))
G01020 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M020[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M020[k,l]*dV(inc))

G00010 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M010[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M010[k,l]*dV(inc))
G11010 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M010[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M010[k,l]*dV(inc))
G22010 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M010[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M010[k,l]*dV(inc))
G12010 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M010[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M010[k,l]*dV(inc))
G02010 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M010[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M010[k,l]*dV(inc))
G01010 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M010[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M010[k,l]*dV(inc))

G00001 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M001[k,l]*dV(inc))
G11001 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M001[k,l]*dV(inc))
G22001 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M001[k,l]*dV(inc))
G12001 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M001[k,l]*dV(inc))
G02001 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M001[k,l]*dV(inc))
G01001 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M001[k,l]*dV(inc))

G00111 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M111[k,l]*dV(inc))
G11111 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M111[k,l]*dV(inc))
G22111 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M111[k,l]*dV(inc))
G12111 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M111[k,l]*dV(inc))
G02111 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M111[k,l]*dV(inc))
G01111 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M111[k,l]*dV(inc))

G00221 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M221[k,l]*dV(inc))
G11221 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M221[k,l]*dV(inc))
G22221 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M221[k,l]*dV(inc))
G12221 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M221[k,l]*dV(inc))
G02221 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M221[k,l]*dV(inc))
G01221 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M221[k,l]*dV(inc))

G00121 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M121[k,l]*dV(inc))
G11121 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M121[k,l]*dV(inc))
G22121 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M121[k,l]*dV(inc))
G12121 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M121[k,l]*dV(inc))
G02121 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M121[k,l]*dV(inc))
G01121 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M121[k,l]*dV(inc))

G00021 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M021[k,l]*dV(inc))
G11021 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M021[k,l]*dV(inc))
G22021 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M021[k,l]*dV(inc))
G12021 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M021[k,l]*dV(inc))
G02021 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M021[k,l]*dV(inc))
G01021 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M021[k,l]*dV(inc))

G00011 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M011[k,l]*dV(inc))
G11011 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M011[k,l]*dV(inc))
G22011 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M011[k,l]*dV(inc))
G12011 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M011[k,l]*dV(inc))
G02011 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M011[k,l]*dV(inc))
G01011 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M011[k,l]*dV(inc))


G00002 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M002[k,l]*dV(inc))
G11002 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M002[k,l]*dV(inc))
G22002 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M002[k,l]*dV(inc))
G12002 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M002[k,l]*dV(inc))
G02002 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M002[k,l]*dV(inc))
G01002 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M002[k,l]*dV(inc))

G00112 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M112[k,l]*dV(inc))
G11112 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M112[k,l]*dV(inc))
G22112 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M112[k,l]*dV(inc))
G12112 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M112[k,l]*dV(inc))
G02112 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M112[k,l]*dV(inc))
G01112 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M112[k,l]*dV(inc))

G00222 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M222[k,l]*dV(inc))
G11222 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M222[k,l]*dV(inc))
G22222 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M222[k,l]*dV(inc))
G12222 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M222[k,l]*dV(inc))
G02222 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M222[k,l]*dV(inc))
G01222 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M222[k,l]*dV(inc))

G00122 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M122[k,l]*dV(inc))
G11122 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M122[k,l]*dV(inc))
G22122 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M122[k,l]*dV(inc))
G12122 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M122[k,l]*dV(inc))
G02122 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M122[k,l]*dV(inc))
G01122 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M122[k,l]*dV(inc))

G00022 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M022[k,l]*dV(inc))
G11022 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M022[k,l]*dV(inc))
G22022 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M022[k,l]*dV(inc))
G12022 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M022[k,l]*dV(inc))
G02022 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M022[k,l]*dV(inc))
G01022 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M022[k,l]*dV(inc))

G00012 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L00[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L00[i,j]*M012[k,l]*dV(inc))
G11012 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L11[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L11[i,j]*M012[k,l]*dV(inc))
G22012 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L22[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L22[i,j]*M012[k,l]*dV(inc))
G12012 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L12[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L12[i,j]*M012[k,l]*dV(inc))
G02012 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L02[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L02[i,j]*M012[k,l]*dV(inc))
G01012 = epsilon/Vol * assemble( Cmat_m[i,j,k,l]*L01[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*L01[i,j]*M012[k,l]*dV(inc))

# if processID == 0: 
# 	elapsed = int(comp_time.time() - start_time)
# 	e_h, e_m, e_s = int(elapsed/3600), int(elapsed % 3600 / 60), int( ( elapsed % 3600 ) % 60 )
# 	print('%.0f DOFs in %.0f h %.0f min %.0f s' % (dofs, e_h, e_m, e_s)  )

if True: 
	with open(data+'.tex', 'a') as file:
		file.write(r'\begeq')
		file.write('\n')
		file.write(r'G\ma_{A\alpha} =')
		file.write('\n')
		file.write(r'\resizebox{.85\textwidth}{!}{$\displaystyle')
		file.write('\n')
		file.write(r'\begin{pmatrix}')
		file.write('\n')
		file.write(str(round(G00000,1))+' & '+str(round(G00110,1))+' & '+str(round(G00220,1))+' & '+str(round(G00120,1))+' & '+str(round(G00020,1))+' & '+str(round(G00010,1))+' & ')
		file.write('\n')
		file.write(str(round(G00001,1))+' & '+str(round(G00111,1))+' & '+str(round(G00221,1))+' & '+str(round(G00121,1))+' & '+str(round(G00021,1))+' & '+str(round(G00011,1))+' & ')
		file.write('\n')
		file.write(str(round(G00002,1))+' & '+str(round(G00112,1))+' & '+str(round(G00222,1))+' & '+str(round(G00122,1))+' & '+str(round(G00022,1))+' & '+str(round(G00012,1))+' & ')
		file.write('\n')
		file.write(r'\\')
		file.write(str(round(G11000,1))+' & '+str(round(G11110,1))+' & '+str(round(G11220,1))+' & '+str(round(G11120,1))+' & '+str(round(G11020,1))+' & '+str(round(G11010,1))+' & ')
		file.write('\n')
		file.write(str(round(G11001,1))+' & '+str(round(G11111,1))+' & '+str(round(G11221,1))+' & '+str(round(G11121,1))+' & '+str(round(G11021,1))+' & '+str(round(G11011,1))+' & ')
		file.write('\n')
		file.write(str(round(G11002,1))+' & '+str(round(G11112,1))+' & '+str(round(G11222,1))+' & '+str(round(G11122,1))+' & '+str(round(G11022,1))+' & '+str(round(G11012,1))+' & ')
		file.write('\n')
		file.write(r'\\')
		file.write(str(round(G22000,1))+' & '+str(round(G22110,1))+' & '+str(round(G22220,1))+' & '+str(round(G22120,1))+' & '+str(round(G22020,1))+' & '+str(round(G22010,1))+' & ')
		file.write('\n')
		file.write(str(round(G22001,1))+' & '+str(round(G22111,1))+' & '+str(round(G22221,1))+' & '+str(round(G22121,1))+' & '+str(round(G22021,1))+' & '+str(round(G22011,1))+' & ')
		file.write('\n')
		file.write(str(round(G22002,1))+' & '+str(round(G22112,1))+' & '+str(round(G22222,1))+' & '+str(round(G22122,1))+' & '+str(round(G22022,1))+' & '+str(round(G22012,1))+' & ')
		file.write('\n')
		file.write(r'\\')
		file.write(str(round(G12000,1))+' & '+str(round(G12110,1))+' & '+str(round(G12220,1))+' & '+str(round(G12120,1))+' & '+str(round(G12020,1))+' & '+str(round(G12010,1))+' & ')
		file.write('\n')
		file.write(str(round(G12001,1))+' & '+str(round(G12111,1))+' & '+str(round(G12221,1))+' & '+str(round(G12121,1))+' & '+str(round(G12021,1))+' & '+str(round(G12011,1))+' & ')
		file.write('\n')
		file.write(str(round(G12002,1))+' & '+str(round(G12112,1))+' & '+str(round(G12222,1))+' & '+str(round(G12122,1))+' & '+str(round(G12022,1))+' & '+str(round(G12012,1))+' & ')
		file.write('\n')
		file.write(r'\\')
		file.write(str(round(G02000,1))+' & '+str(round(G02110,1))+' & '+str(round(G02220,1))+' & '+str(round(G02120,1))+' & '+str(round(G02020,1))+' & '+str(round(G02010,1))+' & ')
		file.write('\n')
		file.write(str(round(G02001,1))+' & '+str(round(G02111,1))+' & '+str(round(G02221,1))+' & '+str(round(G02121,1))+' & '+str(round(G02021,1))+' & '+str(round(G02011,1))+' & ')
		file.write('\n')
		file.write(str(round(G02002,1))+' & '+str(round(G02112,1))+' & '+str(round(G02222,1))+' & '+str(round(G02122,1))+' & '+str(round(G02022,1))+' & '+str(round(G02012,1))+' & ')
		file.write('\n')
		file.write(r'\\')
		file.write(str(round(G01000,1))+' & '+str(round(G01110,1))+' & '+str(round(G01220,1))+' & '+str(round(G01120,1))+' & '+str(round(G01020,1))+' & '+str(round(G01010,1))+' & ')
		file.write('\n')
		file.write(str(round(G01001,1))+' & '+str(round(G01111,1))+' & '+str(round(G01221,1))+' & '+str(round(G01121,1))+' & '+str(round(G01021,1))+' & '+str(round(G01011,1))+' & ')
		file.write('\n')
		file.write(str(round(G01002,1))+' & '+str(round(G01112,1))+' & '+str(round(G01222,1))+' & '+str(round(G01122,1))+' & '+str(round(G01022,1))+' & '+str(round(G01012,1))+' & ')
		file.write('\n')
		file.write(r'\\')
		file.write('\n')
		file.write(r'\end{pmatrix} $} \text{\,N/mm} \ ,')
		file.write('\n')
		file.write(r'\eqend')





#Icf = assemble(x[c]*x[f]*dV)
I00 = epsilon**2 / Vol * assemble(x[0]*x[0]*dV)
I01 = epsilon**2 / Vol * assemble(x[0]*x[1]*dV)
I02 = epsilon**2 / Vol * assemble(x[0]*x[2]*dV)

I10 = I01
I11 = epsilon**2 / Vol * assemble(x[1]*x[1]*dV)
I12 = epsilon**2 / Vol * assemble(x[1]*x[2]*dV)

I20 = I02
I21 = I12
I22 = epsilon**2 / Vol * assemble(x[2]*x[2]*dV)

print('Determining parameters of D')

#Dabcdef = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*Mabc[i,j]*Mdef[k,l]*dV(mat) ) - C_M[a,b,d,e]*Icf ) 
#c=0
D000000 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M000[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M000[k,l]*dV(inc) ) - C_M[0,0,0,0]*I00
D000110 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M110[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M110[k,l]*dV(inc) ) - C_M[0,0,1,1]*I00
D000220 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M220[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M220[k,l]*dV(inc) ) - C_M[0,0,2,2]*I00
D000120 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M120[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M120[k,l]*dV(inc) ) - C_M[0,0,1,2]*I00
D000020 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M020[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M020[k,l]*dV(inc) ) - C_M[0,0,0,2]*I00
D000010 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M010[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M010[k,l]*dV(inc) ) - C_M[0,0,0,1]*I00

D000001 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M001[k,l]*dV(inc) ) - C_M[0,0,0,0]*I01
D000111 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M111[k,l]*dV(inc) ) - C_M[0,0,1,1]*I01
D000221 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M221[k,l]*dV(inc) ) - C_M[0,0,2,2]*I01
D000121 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M121[k,l]*dV(inc) ) - C_M[0,0,1,2]*I01
D000021 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M021[k,l]*dV(inc) ) - C_M[0,0,0,2]*I01
D000011 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M011[k,l]*dV(inc) ) - C_M[0,0,0,1]*I01

D000002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M002[k,l]*dV(inc) ) - C_M[0,0,0,0]*I02
D000112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M112[k,l]*dV(inc) ) - C_M[0,0,1,1]*I02
D000222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M222[k,l]*dV(inc) ) - C_M[0,0,2,2]*I02
D000122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M122[k,l]*dV(inc) ) - C_M[0,0,1,2]*I02
D000022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M022[k,l]*dV(inc) ) - C_M[0,0,0,2]*I02
D000012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M000[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M000[i,j]*M012[k,l]*dV(inc) ) - C_M[0,0,0,1]*I02

D110000 = D000110 #epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M000[k,l]*dV(mat) ) - C_M[1,1,0,0]*I00
D110110 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M110[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M110[k,l]*dV(inc) ) - C_M[1,1,1,1]*I00
D110220 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M220[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M220[k,l]*dV(inc)) - C_M[1,1,2,2]*I00
D110120 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M120[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M120[k,l]*dV(inc)) - C_M[1,1,1,2]*I00
D110020 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M020[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M020[k,l]*dV(inc)) - C_M[1,1,0,2]*I00
D110010 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M010[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M010[k,l]*dV(inc)) - C_M[1,1,0,1]*I00

D110001 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M001[k,l]*dV(inc)) - C_M[1,1,0,0]*I01
D110111 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M111[k,l]*dV(inc)) - C_M[1,1,1,1]*I01
D110221 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M221[k,l]*dV(inc)) - C_M[1,1,2,2]*I01
D110121 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M121[k,l]*dV(inc)) - C_M[1,1,1,2]*I01
D110021 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M021[k,l]*dV(inc)) - C_M[1,1,0,2]*I01
D110011 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M011[k,l]*dV(inc)) - C_M[1,1,0,1]*I01

D110002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M002[k,l]*dV(inc)) - C_M[1,1,0,0]*I02
D110112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M112[k,l]*dV(inc)) - C_M[1,1,1,1]*I02
D110222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M222[k,l]*dV(inc)) - C_M[1,1,2,2]*I02
D110122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M122[k,l]*dV(inc)) - C_M[1,1,1,2]*I02
D110022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M022[k,l]*dV(inc)) - C_M[1,1,0,2]*I02
D110012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M110[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M110[i,j]*M012[k,l]*dV(inc)) - C_M[1,1,0,1]*I02


D220000 = D000220 
D220110 = D110220 
D220220 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M220[i,j]*M220[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M220[i,j]*M220[k,l]*dV(inc)) - C_M[2,2,2,2]*I00
D220120 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M220[i,j]*M120[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M220[i,j]*M120[k,l]*dV(inc)) - C_M[2,2,1,2]*I00
D220020 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M220[i,j]*M020[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M220[i,j]*M020[k,l]*dV(inc)) - C_M[2,2,0,2]*I00
D220010 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M220[i,j]*M010[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M220[i,j]*M010[k,l]*dV(inc)) - C_M[2,2,0,1]*I00

D220001 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M220[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M220[i,j]*M001[k,l]*dV(inc)) - C_M[2,2,0,0]*I01
D220111 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M220[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M220[i,j]*M111[k,l]*dV(inc)) - C_M[2,2,1,1]*I01
D220221 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M220[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M220[i,j]*M221[k,l]*dV(inc)) - C_M[2,2,2,2]*I01
D220121 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M220[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M220[i,j]*M121[k,l]*dV(inc)) - C_M[2,2,1,2]*I01
D220021 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M220[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M220[i,j]*M021[k,l]*dV(inc)) - C_M[2,2,0,2]*I01
D220011 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M220[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M220[i,j]*M011[k,l]*dV(inc)) - C_M[2,2,0,1]*I01

D220002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M220[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M220[i,j]*M002[k,l]*dV(inc)) - C_M[2,2,0,0]*I02
D220112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M220[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M220[i,j]*M112[k,l]*dV(inc)) - C_M[2,2,1,1]*I02
D220222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M220[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M220[i,j]*M222[k,l]*dV(inc)) - C_M[2,2,2,2]*I02
D220122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M220[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M220[i,j]*M122[k,l]*dV(inc)) - C_M[2,2,1,2]*I02
D220022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M220[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M220[i,j]*M022[k,l]*dV(inc)) - C_M[2,2,0,2]*I02
D220012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M220[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M220[i,j]*M012[k,l]*dV(inc)) - C_M[2,2,0,1]*I02


D120000 = D000120 
D120110 = D110120 
D120220 = D220120 
D120120 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M120[i,j]*M120[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M120[i,j]*M120[k,l]*dV(inc)) - C_M[1,2,1,2]*I00
D120020 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M120[i,j]*M020[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M120[i,j]*M020[k,l]*dV(inc)) - C_M[1,2,0,2]*I00
D120010 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M120[i,j]*M010[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M120[i,j]*M010[k,l]*dV(inc)) - C_M[1,2,0,1]*I00

D120001 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M120[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M120[i,j]*M001[k,l]*dV(inc)) - C_M[1,2,0,0]*I01
D120111 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M120[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M120[i,j]*M111[k,l]*dV(inc)) - C_M[1,2,1,1]*I01
D120221 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M120[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M120[i,j]*M221[k,l]*dV(inc)) - C_M[1,2,2,2]*I01
D120121 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M120[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M120[i,j]*M121[k,l]*dV(inc)) - C_M[1,2,1,2]*I01
D120021 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M120[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M120[i,j]*M021[k,l]*dV(inc)) - C_M[1,2,0,2]*I01
D120011 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M120[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M120[i,j]*M011[k,l]*dV(inc)) - C_M[1,2,0,1]*I01

D120002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M120[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M120[i,j]*M002[k,l]*dV(inc)) - C_M[1,2,0,0]*I02
D120112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M120[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M120[i,j]*M112[k,l]*dV(inc)) - C_M[1,2,1,1]*I02
D120222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M120[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M120[i,j]*M222[k,l]*dV(inc)) - C_M[1,2,2,2]*I02
D120122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M120[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M120[i,j]*M122[k,l]*dV(inc)) - C_M[1,2,1,2]*I02
D120022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M120[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M120[i,j]*M022[k,l]*dV(inc)) - C_M[1,2,0,2]*I02
D120012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M120[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M120[i,j]*M012[k,l]*dV(inc)) - C_M[1,2,0,1]*I02


D020000 = D000020 
D020110 = D110020 
D020220 = D220020 
D020120 = D120020 
D020020 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M020[i,j]*M020[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M020[i,j]*M020[k,l]*dV(inc)) - C_M[0,2,0,2]*I00
D020010 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M020[i,j]*M010[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M020[i,j]*M010[k,l]*dV(inc)) - C_M[0,2,0,1]*I00

D020001 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M020[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M020[i,j]*M001[k,l]*dV(inc)) - C_M[0,2,0,0]*I01
D020111 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M020[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M020[i,j]*M111[k,l]*dV(inc)) - C_M[0,2,1,1]*I01
D020221 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M020[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M020[i,j]*M221[k,l]*dV(inc)) - C_M[0,2,2,2]*I01
D020121 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M020[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M020[i,j]*M121[k,l]*dV(inc)) - C_M[0,2,1,2]*I01
D020021 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M020[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M020[i,j]*M021[k,l]*dV(inc)) - C_M[0,2,0,2]*I01
D020011 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M020[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M020[i,j]*M011[k,l]*dV(inc)) - C_M[0,2,0,1]*I01

D020002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M020[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M020[i,j]*M002[k,l]*dV(inc)) - C_M[0,2,0,0]*I02
D020112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M020[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M020[i,j]*M112[k,l]*dV(inc)) - C_M[0,2,1,1]*I02
D020222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M020[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M020[i,j]*M222[k,l]*dV(inc)) - C_M[0,2,2,2]*I02
D020122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M020[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M020[i,j]*M122[k,l]*dV(inc)) - C_M[0,2,1,2]*I02
D020022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M020[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M020[i,j]*M022[k,l]*dV(inc)) - C_M[0,2,0,2]*I02
D020012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M020[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M020[i,j]*M012[k,l]*dV(inc)) - C_M[0,2,0,1]*I02


D010000 = D000010 
D010110 = D110010 
D010220 = D220010 
D010120 = D120010 
D010020 = D020010 
D010010 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M010[i,j]*M010[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M010[i,j]*M010[k,l]*dV(inc)) - C_M[0,1,0,1]*I00

D010001 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M010[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M010[i,j]*M001[k,l]*dV(inc)) - C_M[0,1,0,0]*I01
D010111 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M010[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M010[i,j]*M111[k,l]*dV(inc)) - C_M[0,1,1,1]*I01
D010221 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M010[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M010[i,j]*M221[k,l]*dV(inc)) - C_M[0,1,2,2]*I01
D010121 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M010[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M010[i,j]*M121[k,l]*dV(inc)) - C_M[0,1,1,2]*I01
D010021 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M010[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M010[i,j]*M021[k,l]*dV(inc)) - C_M[0,1,0,2]*I01
D010011 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M010[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M010[i,j]*M011[k,l]*dV(inc)) - C_M[0,1,0,1]*I01

D010002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M010[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M010[i,j]*M002[k,l]*dV(inc)) - C_M[0,1,0,0]*I02
D010112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M010[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M010[i,j]*M112[k,l]*dV(inc)) - C_M[0,1,1,1]*I02
D010222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M010[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M010[i,j]*M222[k,l]*dV(inc)) - C_M[0,1,2,2]*I02
D010122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M010[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M010[i,j]*M122[k,l]*dV(inc)) - C_M[0,1,1,2]*I02
D010022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M010[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M010[i,j]*M022[k,l]*dV(inc)) - C_M[0,1,0,2]*I02
D010012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M010[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M010[i,j]*M012[k,l]*dV(inc)) - C_M[0,1,0,1]*I02


#c=1
D001000 = D000001 
D001110 = D110001 
D001220 = D220001 
D001120 = D120001 
D001020 = D020001 
D001010 = D010001 
D001001 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M001[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M001[i,j]*M001[k,l]*dV(inc)) - C_M[0,0,0,0]*I11
D001111 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M001[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M001[i,j]*M111[k,l]*dV(inc)) - C_M[0,0,1,1]*I11
D001221 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M001[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M001[i,j]*M221[k,l]*dV(inc)) - C_M[0,0,2,2]*I11
D001121 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M001[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M001[i,j]*M121[k,l]*dV(inc)) - C_M[0,0,1,2]*I11
D001021 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M001[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M001[i,j]*M021[k,l]*dV(inc)) - C_M[0,0,0,2]*I11
D001011 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M001[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M001[i,j]*M011[k,l]*dV(inc)) - C_M[0,0,0,1]*I11
D001002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M001[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M001[i,j]*M002[k,l]*dV(inc)) - C_M[0,0,0,0]*I12
D001112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M001[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M001[i,j]*M112[k,l]*dV(inc)) - C_M[0,0,1,1]*I12
D001222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M001[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M001[i,j]*M222[k,l]*dV(inc)) - C_M[0,0,2,2]*I12
D001122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M001[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M001[i,j]*M122[k,l]*dV(inc)) - C_M[0,0,1,2]*I12
D001022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M001[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M001[i,j]*M022[k,l]*dV(inc)) - C_M[0,0,0,2]*I12
D001012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M001[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M001[i,j]*M012[k,l]*dV(inc)) - C_M[0,0,0,1]*I12
D111000 = D000111
D111110 = D110111 
D111220 = D220111 
D111120 = D120111 
D111020 = D020111 
D111010 = D010111 
D111001 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M111[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M111[i,j]*M001[k,l]*dV(inc)) - C_M[1,1,0,0]*I11
D111111 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M111[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M111[i,j]*M111[k,l]*dV(inc)) - C_M[1,1,1,1]*I11
D111221 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M111[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M111[i,j]*M221[k,l]*dV(inc)) - C_M[1,1,2,2]*I11
D111121 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M111[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M111[i,j]*M121[k,l]*dV(inc)) - C_M[1,1,1,2]*I11
D111021 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M111[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M111[i,j]*M021[k,l]*dV(inc)) - C_M[1,1,0,2]*I11
D111011 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M111[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M111[i,j]*M011[k,l]*dV(inc)) - C_M[1,1,0,1]*I11
D111002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M111[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M111[i,j]*M002[k,l]*dV(inc)) - C_M[1,1,0,0]*I12
D111112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M111[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M111[i,j]*M112[k,l]*dV(inc)) - C_M[1,1,1,1]*I12
D111222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M111[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M111[i,j]*M222[k,l]*dV(inc)) - C_M[1,1,2,2]*I12
D111122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M111[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M111[i,j]*M122[k,l]*dV(inc)) - C_M[1,1,1,2]*I12
D111022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M111[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M111[i,j]*M022[k,l]*dV(inc)) - C_M[1,1,0,2]*I12
D111012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M111[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M111[i,j]*M012[k,l]*dV(inc)) - C_M[1,1,0,1]*I12
D221000 = D000221
D221110 = D110221
D221220 = D220221 
D221120 = D120221 
D221020 = D020221 
D221010 = D010221 
D221001 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M221[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M221[i,j]*M001[k,l]*dV(inc)) - C_M[2,2,0,0]*I11
D221111 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M221[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M221[i,j]*M111[k,l]*dV(inc)) - C_M[2,2,1,1]*I11
D221221 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M221[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M221[i,j]*M221[k,l]*dV(inc)) - C_M[2,2,2,2]*I11
D221121 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M221[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M221[i,j]*M121[k,l]*dV(inc)) - C_M[2,2,1,2]*I11
D221021 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M221[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M221[i,j]*M021[k,l]*dV(inc)) - C_M[2,2,0,2]*I11
D221011 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M221[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M221[i,j]*M011[k,l]*dV(inc)) - C_M[2,2,0,1]*I11
D221002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M221[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M221[i,j]*M002[k,l]*dV(inc)) - C_M[2,2,0,0]*I12
D221112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M221[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M221[i,j]*M112[k,l]*dV(inc)) - C_M[2,2,1,1]*I12
D221222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M221[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M221[i,j]*M222[k,l]*dV(inc)) - C_M[2,2,2,2]*I12
D221122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M221[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M221[i,j]*M122[k,l]*dV(inc)) - C_M[2,2,1,2]*I12
D221022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M221[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M221[i,j]*M022[k,l]*dV(inc)) - C_M[2,2,0,2]*I12
D221012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M221[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M221[i,j]*M012[k,l]*dV(inc)) - C_M[2,2,0,1]*I12
D121000 = D000121
D121110 = D110121
D121220 = D220121
D121120 = D120121 
D121020 = D020121 
D121010 = D010121 
D121001 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M121[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M121[i,j]*M001[k,l]*dV(inc)) - C_M[1,2,0,0]*I11
D121111 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M121[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M121[i,j]*M111[k,l]*dV(inc)) - C_M[1,2,1,1]*I11
D121221 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M121[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M121[i,j]*M221[k,l]*dV(inc)) - C_M[1,2,2,2]*I11
D121121 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M121[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M121[i,j]*M121[k,l]*dV(inc)) - C_M[1,2,1,2]*I11
D121021 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M121[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M121[i,j]*M021[k,l]*dV(inc)) - C_M[1,2,0,2]*I11
D121011 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M121[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M121[i,j]*M011[k,l]*dV(inc)) - C_M[1,2,0,1]*I11
D121002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M121[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M121[i,j]*M002[k,l]*dV(inc)) - C_M[1,2,0,0]*I12
D121112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M121[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M121[i,j]*M112[k,l]*dV(inc)) - C_M[1,2,1,1]*I12
D121222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M121[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M121[i,j]*M222[k,l]*dV(inc)) - C_M[1,2,2,2]*I12
D121122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M121[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M121[i,j]*M122[k,l]*dV(inc)) - C_M[1,2,1,2]*I12
D121022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M121[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M121[i,j]*M022[k,l]*dV(inc)) - C_M[1,2,0,2]*I12
D121012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M121[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M121[i,j]*M012[k,l]*dV(inc)) - C_M[1,2,0,1]*I12
D021000 = D000021
D021110 = D110021
D021220 = D220021
D021120 = D120021
D021020 = D020021 
D021010 = D010021 
D021001 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M021[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M021[i,j]*M001[k,l]*dV(inc)) - C_M[0,2,0,0]*I11
D021111 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M021[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M021[i,j]*M111[k,l]*dV(inc)) - C_M[0,2,1,1]*I11
D021221 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M021[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M021[i,j]*M221[k,l]*dV(inc)) - C_M[0,2,2,2]*I11
D021121 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M021[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M021[i,j]*M121[k,l]*dV(inc)) - C_M[0,2,1,2]*I11
D021021 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M021[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M021[i,j]*M021[k,l]*dV(inc)) - C_M[0,2,0,2]*I11
D021011 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M021[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M021[i,j]*M011[k,l]*dV(inc)) - C_M[0,2,0,1]*I11
D021002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M021[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M021[i,j]*M002[k,l]*dV(inc)) - C_M[0,2,0,0]*I12
D021112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M021[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M021[i,j]*M112[k,l]*dV(inc)) - C_M[0,2,1,1]*I12
D021222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M021[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M021[i,j]*M222[k,l]*dV(inc)) - C_M[0,2,2,2]*I12
D021122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M021[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M021[i,j]*M122[k,l]*dV(inc)) - C_M[0,2,1,2]*I12
D021022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M021[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M021[i,j]*M022[k,l]*dV(inc)) - C_M[0,2,0,2]*I12
D021012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M021[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M021[i,j]*M012[k,l]*dV(inc)) - C_M[0,2,0,1]*I12
D011000 = D000011
D011110 = D110011
D011220 = D220011
D011120 = D120011
D011020 = D020011
D011010 = D010011 
D011001 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M011[i,j]*M001[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M011[i,j]*M001[k,l]*dV(inc)) - C_M[0,1,0,0]*I11
D011111 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M011[i,j]*M111[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M011[i,j]*M111[k,l]*dV(inc)) - C_M[0,1,1,1]*I11
D011221 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M011[i,j]*M221[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M011[i,j]*M221[k,l]*dV(inc)) - C_M[0,1,2,2]*I11
D011121 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M011[i,j]*M121[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M011[i,j]*M121[k,l]*dV(inc)) - C_M[0,1,1,2]*I11
D011021 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M011[i,j]*M021[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M011[i,j]*M021[k,l]*dV(inc)) - C_M[0,1,0,2]*I11
D011011 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M011[i,j]*M011[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M011[i,j]*M011[k,l]*dV(inc)) - C_M[0,1,0,1]*I11
D011002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M011[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M011[i,j]*M002[k,l]*dV(inc)) - C_M[0,1,0,0]*I12
D011112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M011[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M011[i,j]*M112[k,l]*dV(inc)) - C_M[0,1,1,1]*I12
D011222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M011[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M011[i,j]*M222[k,l]*dV(inc)) - C_M[0,1,2,2]*I12
D011122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M011[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M011[i,j]*M122[k,l]*dV(inc)) - C_M[0,1,1,2]*I12
D011022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M011[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M011[i,j]*M022[k,l]*dV(inc)) - C_M[0,1,0,2]*I12
D011012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M011[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M011[i,j]*M012[k,l]*dV(inc)) - C_M[0,1,0,1]*I12


#c=2
D002000 = D000002 
D002110 = D110002 
D002220 = D220002 
D002120 = D120002 
D002020 = D020002 
D002010 = D010002 
D002001 = D001002 
D002111 = D111002 
D002221 = D221002 
D002121 = D121002 
D002021 = D021002 
D002011 = D011002 
D002002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M002[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M002[i,j]*M002[k,l]*dV(inc)) - C_M[0,0,0,0]*I22
D002112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M002[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M002[i,j]*M112[k,l]*dV(inc)) - C_M[0,0,1,1]*I22
D002222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M002[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M002[i,j]*M222[k,l]*dV(inc)) - C_M[0,0,2,2]*I22
D002122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M002[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M002[i,j]*M122[k,l]*dV(inc)) - C_M[0,0,1,2]*I22
D002022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M002[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M002[i,j]*M022[k,l]*dV(inc)) - C_M[0,0,0,2]*I22
D002012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M002[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M002[i,j]*M012[k,l]*dV(inc)) - C_M[0,0,0,1]*I22
D112000 = D000112
D112110 = D110112 
D112220 = D220112 
D112120 = D120112 
D112020 = D020112 
D112010 = D010112 
D112001 = D001112 
D112111 = D111112 
D112221 = D221112 
D112121 = D121112 
D112021 = D021112 
D112011 = D011112 
D112002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M112[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M112[i,j]*M002[k,l]*dV(inc)) - C_M[1,1,0,0]*I22
D112112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M112[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M112[i,j]*M112[k,l]*dV(inc)) - C_M[1,1,1,1]*I22
D112222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M112[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M112[i,j]*M222[k,l]*dV(inc)) - C_M[1,1,2,2]*I22
D112122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M112[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M112[i,j]*M122[k,l]*dV(inc)) - C_M[1,1,1,2]*I22
D112022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M112[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M112[i,j]*M022[k,l]*dV(inc)) - C_M[1,1,0,2]*I22
D112012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M112[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M112[i,j]*M012[k,l]*dV(inc)) - C_M[1,1,0,1]*I22
D222000 = D000222
D222110 = D110222
D222220 = D220222 
D222120 = D120222 
D222020 = D020222 
D222010 = D010222 
D222001 = D001222 
D222111 = D111222 
D222221 = D221222 
D222121 = D121222 
D222021 = D021222 
D222011 = D011222 
D222002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M222[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M222[i,j]*M002[k,l]*dV(inc)) - C_M[2,2,0,0]*I22
D222112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M222[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M222[i,j]*M112[k,l]*dV(inc)) - C_M[2,2,1,1]*I22
D222222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M222[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M222[i,j]*M222[k,l]*dV(inc)) - C_M[2,2,2,2]*I22
D222122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M222[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M222[i,j]*M122[k,l]*dV(inc)) - C_M[2,2,1,2]*I22
D222022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M222[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M222[i,j]*M022[k,l]*dV(inc)) - C_M[2,2,0,2]*I22
D222012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M222[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M222[i,j]*M012[k,l]*dV(inc)) - C_M[2,2,0,1]*I22
D122000 = D000122
D122110 = D110122
D122220 = D220122
D122120 = D120122
D122020 = D020122
D122010 = D010122
D122001 = D001122 
D122111 = D111122 
D122221 = D221122 
D122121 = D121122 
D122021 = D021122 
D122011 = D011122 
D122002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M122[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M122[i,j]*M002[k,l]*dV(inc)) - C_M[1,2,0,0]*I22
D122112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M122[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M122[i,j]*M112[k,l]*dV(inc)) - C_M[1,2,1,1]*I22
D122222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M122[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M122[i,j]*M222[k,l]*dV(inc)) - C_M[1,2,2,2]*I22
D122122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M122[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M122[i,j]*M122[k,l]*dV(inc)) - C_M[1,2,1,2]*I22
D122022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M122[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M122[i,j]*M022[k,l]*dV(inc)) - C_M[1,2,0,2]*I22
D122012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M122[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M122[i,j]*M012[k,l]*dV(inc)) - C_M[1,2,0,1]*I22
D022000 = D000022
D022110 = D110022
D022220 = D220022
D022120 = D120022
D022020 = D020022 
D022010 = D010022 
D022001 = D001022 
D022111 = D111022 
D022221 = D221022 
D022121 = D121022 
D022021 = D021022 
D022011 = D011022 
D022002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M022[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M022[i,j]*M002[k,l]*dV(inc)) - C_M[0,2,0,0]*I22
D022112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M022[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M022[i,j]*M112[k,l]*dV(inc)) - C_M[0,2,1,1]*I22
D022222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M022[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M022[i,j]*M222[k,l]*dV(inc)) - C_M[0,2,2,2]*I22
D022122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M022[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M022[i,j]*M122[k,l]*dV(inc)) - C_M[0,2,1,2]*I22
D022022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M022[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M022[i,j]*M022[k,l]*dV(inc)) - C_M[0,2,0,2]*I22
D022012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M022[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M022[i,j]*M012[k,l]*dV(inc)) - C_M[0,2,0,1]*I22
D012000 = D000012
D012110 = D110012
D012220 = D220012
D012120 = D120012
D012020 = D020012
D012010 = D010012 
D012001 = D001012 
D012111 = D111012 
D012221 = D221012 
D012121 = D121012 
D012021 = D021012 
D012011 = D011012 
D012002 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M012[i,j]*M002[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M012[i,j]*M002[k,l]*dV(inc)) - C_M[0,1,0,0]*I22
D012112 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M012[i,j]*M112[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M012[i,j]*M112[k,l]*dV(inc)) - C_M[0,1,1,1]*I22
D012222 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M012[i,j]*M222[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M012[i,j]*M222[k,l]*dV(inc)) - C_M[0,1,2,2]*I22
D012122 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M012[i,j]*M122[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M012[i,j]*M122[k,l]*dV(inc)) - C_M[0,1,1,2]*I22
D012022 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M012[i,j]*M022[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M012[i,j]*M022[k,l]*dV(inc)) - C_M[0,1,0,2]*I22
D012012 = epsilon**2 / Vol * assemble( Cmat_m[i,j,k,l]*M012[i,j]*M012[k,l]*dV(mat) + Cinc_m[i,j,k,l]*M012[i,j]*M012[k,l]*dV(inc)) - C_M[0,1,0,1]*I22

# if processID == 0: 
# 	elapsed = int(comp_time.time() - start_time)
# 	e_h, e_m, e_s = int(elapsed/3600), int(elapsed % 3600 / 60), int( ( elapsed % 3600 ) % 60 )
# 	print('%.0f DOFs in %.0f h %.0f min %.0f s' % (dofs, e_h, e_m, e_s)  )


divd = 1
if True:  
	with open(data+'.tex', 'a') as file:
		file.write(r'\begeq')
		file.write('\n')
		file.write(r'D\ma_{\alpha\beta} =')
		file.write('\n')
		file.write(r'\resizebox{.85\textwidth}{!}{$\displaystyle')
		file.write('\n')
		file.write(r'\begin{pmatrix}')
		file.write('\n')
		file.write(str(round(D000000/divd,Sf))+' & '+str(round(D000110/divd,Sf))+' & '+str(round(D000011/divd,Sf))+' & '+str(round(D000220/divd,Sf))+' & '+str(round(D000022/divd,Sf))+' & '+str(round(D000111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D000001/divd,Sf))+' & '+str(round(D000010/divd,Sf))+' & '+str(round(D000221/divd,Sf))+' & '+str(round(D000122/divd,Sf))+' & '+str(round(D000222/divd,Sf))+' & '+str(round(D000002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D000020/divd,Sf))+' & '+str(round(D000112/divd,Sf))+' & '+str(round(D000121/divd,Sf))+' & '+str(round(D000120/divd,Sf))+' & '+str(round(D000021/divd,Sf))+' & '+str(round(D000012/divd,Sf))+r' \\')

		file.write('\n')
		file.write(str(round(D110000/divd,Sf))+' & '+str(round(D110110/divd,Sf))+' & '+str(round(D110011/divd,Sf))+' & '+str(round(D110220/divd,Sf))+' & '+str(round(D110022/divd,Sf))+' & '+str(round(D110111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D110001/divd,Sf))+' & '+str(round(D110010/divd,Sf))+' & '+str(round(D110221/divd,Sf))+' & '+str(round(D110122/divd,Sf))+' & '+str(round(D110222/divd,Sf))+' & '+str(round(D110002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D110020/divd,Sf))+' & '+str(round(D110112/divd,Sf))+' & '+str(round(D110121/divd,Sf))+' & '+str(round(D110120/divd,Sf))+' & '+str(round(D110021/divd,Sf))+' & '+str(round(D110012/divd,Sf))+r' \\')

		file.write('\n')
		file.write(str(round(D011000/divd,Sf))+' & '+str(round(D011110/divd,Sf))+' & '+str(round(D011011/divd,Sf))+' & '+str(round(D011220/divd,Sf))+' & '+str(round(D011022/divd,Sf))+' & '+str(round(D011111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D011001/divd,Sf))+' & '+str(round(D011010/divd,Sf))+' & '+str(round(D011221/divd,Sf))+' & '+str(round(D011122/divd,Sf))+' & '+str(round(D011222/divd,Sf))+' & '+str(round(D011002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D011020/divd,Sf))+' & '+str(round(D011112/divd,Sf))+' & '+str(round(D011121/divd,Sf))+' & '+str(round(D011120/divd,Sf))+' & '+str(round(D011021/divd,Sf))+' & '+str(round(D011012/divd,Sf))+r' \\')

		file.write('\n')
		file.write(str(round(D220000/divd,Sf))+' & '+str(round(D220110/divd,Sf))+' & '+str(round(D220011/divd,Sf))+' & '+str(round(D220220/divd,Sf))+' & '+str(round(D220022/divd,Sf))+' & '+str(round(D220111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D220001/divd,Sf))+' & '+str(round(D220010/divd,Sf))+' & '+str(round(D220221/divd,Sf))+' & '+str(round(D220122/divd,Sf))+' & '+str(round(D220222/divd,Sf))+' & '+str(round(D220002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D220020/divd,Sf))+' & '+str(round(D220112/divd,Sf))+' & '+str(round(D220121/divd,Sf))+' & '+str(round(D220120/divd,Sf))+' & '+str(round(D220021/divd,Sf))+' & '+str(round(D220012/divd,Sf))+r' \\')

		file.write('\n')
		file.write(str(round(D022000/divd,Sf))+' & '+str(round(D022110/divd,Sf))+' & '+str(round(D022011/divd,Sf))+' & '+str(round(D022220/divd,Sf))+' & '+str(round(D022022/divd,Sf))+' & '+str(round(D022111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D022001/divd,Sf))+' & '+str(round(D022010/divd,Sf))+' & '+str(round(D022221/divd,Sf))+' & '+str(round(D022122/divd,Sf))+' & '+str(round(D022222/divd,Sf))+' & '+str(round(D022002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D022020/divd,Sf))+' & '+str(round(D022112/divd,Sf))+' & '+str(round(D022121/divd,Sf))+' & '+str(round(D022120/divd,Sf))+' & '+str(round(D022021/divd,Sf))+' & '+str(round(D022012/divd,Sf))+r' \\')

		file.write('\n')
		file.write(str(round(D111000/divd,Sf))+' & '+str(round(D111110/divd,Sf))+' & '+str(round(D111011/divd,Sf))+' & '+str(round(D111220/divd,Sf))+' & '+str(round(D111022/divd,Sf))+' & '+str(round(D111111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D111001/divd,Sf))+' & '+str(round(D111010/divd,Sf))+' & '+str(round(D111221/divd,Sf))+' & '+str(round(D111122/divd,Sf))+' & '+str(round(D111222/divd,Sf))+' & '+str(round(D111002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D111020/divd,Sf))+' & '+str(round(D111112/divd,Sf))+' & '+str(round(D111121/divd,Sf))+' & '+str(round(D111120/divd,Sf))+' & '+str(round(D111021/divd,Sf))+' & '+str(round(D111012/divd,Sf)))
		file.write('\n')
		file.write(r' \\')

		file.write('\n')
		file.write(str(round(D001000/divd,Sf))+' & '+str(round(D001110/divd,Sf))+' & '+str(round(D001011/divd,Sf))+' & '+str(round(D001220/divd,Sf))+' & '+str(round(D001022/divd,Sf))+' & '+str(round(D001111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D001001/divd,Sf))+' & '+str(round(D001010/divd,Sf))+' & '+str(round(D001221/divd,Sf))+' & '+str(round(D001122/divd,Sf))+' & '+str(round(D001222/divd,Sf))+' & '+str(round(D001002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D001020/divd,Sf))+' & '+str(round(D001112/divd,Sf))+' & '+str(round(D001121/divd,Sf))+' & '+str(round(D001120/divd,Sf))+' & '+str(round(D001021/divd,Sf))+' & '+str(round(D001012/divd,Sf))+r' \\')

		file.write('\n')
		file.write(str(round(D010000/divd,Sf))+' & '+str(round(D010110/divd,Sf))+' & '+str(round(D010011/divd,Sf))+' & '+str(round(D010220/divd,Sf))+' & '+str(round(D010022/divd,Sf))+' & '+str(round(D010111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D010001/divd,Sf))+' & '+str(round(D010010/divd,Sf))+' & '+str(round(D010221/divd,Sf))+' & '+str(round(D010122/divd,Sf))+' & '+str(round(D010222/divd,Sf))+' & '+str(round(D010002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D010020/divd,Sf))+' & '+str(round(D010112/divd,Sf))+' & '+str(round(D010121/divd,Sf))+' & '+str(round(D010120/divd,Sf))+' & '+str(round(D010021/divd,Sf))+' & '+str(round(D010012/divd,Sf))+r' \\')

		file.write('\n')
		file.write(str(round(D221000/divd,Sf))+' & '+str(round(D221110/divd,Sf))+' & '+str(round(D221011/divd,Sf))+' & '+str(round(D221220/divd,Sf))+' & '+str(round(D221022/divd,Sf))+' & '+str(round(D221111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D221001/divd,Sf))+' & '+str(round(D221010/divd,Sf))+' & '+str(round(D221221/divd,Sf))+' & '+str(round(D221122/divd,Sf))+' & '+str(round(D221222/divd,Sf))+' & '+str(round(D221002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D221020/divd,Sf))+' & '+str(round(D221112/divd,Sf))+' & '+str(round(D221121/divd,Sf))+' & '+str(round(D221120/divd,Sf))+' & '+str(round(D221021/divd,Sf))+' & '+str(round(D221012/divd,Sf))+r' \\')

		file.write('\n')
		file.write(str(round(D122000/divd,Sf))+' & '+str(round(D122110/divd,Sf))+' & '+str(round(D122011/divd,Sf))+' & '+str(round(D122220/divd,Sf))+' & '+str(round(D122022/divd,Sf))+' & '+str(round(D122111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D122001/divd,Sf))+' & '+str(round(D122010/divd,Sf))+' & '+str(round(D122221/divd,Sf))+' & '+str(round(D122122/divd,Sf))+' & '+str(round(D122222/divd,Sf))+' & '+str(round(D122002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D122020/divd,Sf))+' & '+str(round(D122112/divd,Sf))+' & '+str(round(D122121/divd,Sf))+' & '+str(round(D122120/divd,Sf))+' & '+str(round(D122021/divd,Sf))+' & '+str(round(D122012/divd,Sf))+r' \\')

		file.write('\n')
		file.write(str(round(D222000/divd,Sf))+' & '+str(round(D222110/divd,Sf))+' & '+str(round(D222011/divd,Sf))+' & '+str(round(D222220/divd,Sf))+' & '+str(round(D222022/divd,Sf))+' & '+str(round(D222111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D222001/divd,Sf))+' & '+str(round(D222010/divd,Sf))+' & '+str(round(D222221/divd,Sf))+' & '+str(round(D222122/divd,Sf))+' & '+str(round(D222222/divd,Sf))+' & '+str(round(D222002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D222020/divd,Sf))+' & '+str(round(D222112/divd,Sf))+' & '+str(round(D222121/divd,Sf))+' & '+str(round(D222120/divd,Sf))+' & '+str(round(D222021/divd,Sf))+' & '+str(round(D222012/divd,Sf))+r' \\')

		file.write('\n')
		file.write(str(round(D002000/divd,Sf))+' & '+str(round(D002110/divd,Sf))+' & '+str(round(D002011/divd,Sf))+' & '+str(round(D002220/divd,Sf))+' & '+str(round(D002022/divd,Sf))+' & '+str(round(D002111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D002001/divd,Sf))+' & '+str(round(D002010/divd,Sf))+' & '+str(round(D002221/divd,Sf))+' & '+str(round(D002122/divd,Sf))+' & '+str(round(D002222/divd,Sf))+' & '+str(round(D002002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D002020/divd,Sf))+' & '+str(round(D002112/divd,Sf))+' & '+str(round(D002121/divd,Sf))+' & '+str(round(D002120/divd,Sf))+' & '+str(round(D002021/divd,Sf))+' & '+str(round(D002012/divd,Sf))+r' \\')

		file.write('\n')
		file.write(str(round(D020000/divd,Sf))+' & '+str(round(D020110/divd,Sf))+' & '+str(round(D020011/divd,Sf))+' & '+str(round(D020220/divd,Sf))+' & '+str(round(D020022/divd,Sf))+' & '+str(round(D020111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D020001/divd,Sf))+' & '+str(round(D020010/divd,Sf))+' & '+str(round(D020221/divd,Sf))+' & '+str(round(D020122/divd,Sf))+' & '+str(round(D020222/divd,Sf))+' & '+str(round(D020002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D020020/divd,Sf))+' & '+str(round(D020112/divd,Sf))+' & '+str(round(D020121/divd,Sf))+' & '+str(round(D020120/divd,Sf))+' & '+str(round(D020021/divd,Sf))+' & '+str(round(D020012/divd,Sf))+r' \\')

		file.write('\n')
		file.write(str(round(D112000/divd,Sf))+' & '+str(round(D112110/divd,Sf))+' & '+str(round(D112011/divd,Sf))+' & '+str(round(D112220/divd,Sf))+' & '+str(round(D112022/divd,Sf))+' & '+str(round(D112111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D112001/divd,Sf))+' & '+str(round(D112010/divd,Sf))+' & '+str(round(D112221/divd,Sf))+' & '+str(round(D112122/divd,Sf))+' & '+str(round(D112222/divd,Sf))+' & '+str(round(D112002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D112020/divd,Sf))+' & '+str(round(D112112/divd,Sf))+' & '+str(round(D112121/divd,Sf))+' & '+str(round(D112120/divd,Sf))+' & '+str(round(D112021/divd,Sf))+' & '+str(round(D112012/divd,Sf))+r' \\')

		file.write('\n')
		file.write(str(round(D121000/divd,Sf))+' & '+str(round(D121110/divd,Sf))+' & '+str(round(D121011/divd,Sf))+' & '+str(round(D121220/divd,Sf))+' & '+str(round(D121022/divd,Sf))+' & '+str(round(D121111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D121001/divd,Sf))+' & '+str(round(D121010/divd,Sf))+' & '+str(round(D121221/divd,Sf))+' & '+str(round(D121122/divd,Sf))+' & '+str(round(D121222/divd,Sf))+' & '+str(round(D121002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D121020/divd,Sf))+' & '+str(round(D121112/divd,Sf))+' & '+str(round(D121121/divd,Sf))+' & '+str(round(D121120/divd,Sf))+' & '+str(round(D121021/divd,Sf))+' & '+str(round(D121012/divd,Sf))+r' \\')

		file.write('\n')
		file.write(str(round(D120000/divd,Sf))+' & '+str(round(D120110/divd,Sf))+' & '+str(round(D120011/divd,Sf))+' & '+str(round(D120220/divd,Sf))+' & '+str(round(D120022/divd,Sf))+' & '+str(round(D120111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D120001/divd,Sf))+' & '+str(round(D120010/divd,Sf))+' & '+str(round(D120221/divd,Sf))+' & '+str(round(D120122/divd,Sf))+' & '+str(round(D120222/divd,Sf))+' & '+str(round(D120002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D120020/divd,Sf))+' & '+str(round(D120112/divd,Sf))+' & '+str(round(D120121/divd,Sf))+' & '+str(round(D120120/divd,Sf))+' & '+str(round(D120021/divd,Sf))+' & '+str(round(D120012/divd,Sf))+r' \\')

		file.write('\n')
		file.write(str(round(D021000/divd,Sf))+' & '+str(round(D021110/divd,Sf))+' & '+str(round(D021011/divd,Sf))+' & '+str(round(D021220/divd,Sf))+' & '+str(round(D021022/divd,Sf))+' & '+str(round(D021111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D021001/divd,Sf))+' & '+str(round(D021010/divd,Sf))+' & '+str(round(D021221/divd,Sf))+' & '+str(round(D021122/divd,Sf))+' & '+str(round(D021222/divd,Sf))+' & '+str(round(D021002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D021020/divd,Sf))+' & '+str(round(D021112/divd,Sf))+' & '+str(round(D021121/divd,Sf))+' & '+str(round(D021120/divd,Sf))+' & '+str(round(D021021/divd,Sf))+' & '+str(round(D021012/divd,Sf))+r' \\')

		file.write('\n')
		file.write(str(round(D012000/divd,Sf))+' & '+str(round(D012110/divd,Sf))+' & '+str(round(D012011/divd,Sf))+' & '+str(round(D012220/divd,Sf))+' & '+str(round(D012022/divd,Sf))+' & '+str(round(D012111/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D012001/divd,Sf))+' & '+str(round(D012010/divd,Sf))+' & '+str(round(D012221/divd,Sf))+' & '+str(round(D012122/divd,Sf))+' & '+str(round(D012222/divd,Sf))+' & '+str(round(D012002/divd,Sf))+' & ')
		file.write('\n')
		file.write(str(round(D012020/divd,Sf))+' & '+str(round(D012112/divd,Sf))+' & '+str(round(D012121/divd,Sf))+' & '+str(round(D012120/divd,Sf))+' & '+str(round(D012021/divd,Sf))+' & '+str(round(D012012/divd,Sf)))

		file.write('\n')
		file.write('\n')
		file.write(r'\end{pmatrix} $} \text{\,N} \ ,')
		file.write('\n')
		file.write(r'\eqend')



	with open(data+'.tex', 'a') as file:
		file.write(r'\begeq')
		file.write('\n')
		file.write(r'D\ma_{\alpha\beta} =')
		file.write('\n')
		file.write(r'\resizebox{.85\textwidth}{!}{$\displaystyle')
		file.write('\n')
		file.write(r'\begin{pmatrix}')
		file.write('\n')
		file.write(str(round(D000000,1))+' & '+str(round(D000110,1))+' & '+str(round(D000220,1))+' & '+str(round(D000120,1))+' & '+str(round(D000020,1))+' & '+str(round(D000010,1))+' & ')
		file.write('\n')
		file.write(str(round(D000001,1))+' & '+str(round(D000111,1))+' & '+str(round(D000221,1))+' & '+str(round(D000121,1))+' & '+str(round(D000021,1))+' & '+str(round(D000011,1))+' & ')
		file.write('\n')
		file.write(str(round(D000002,1))+' & '+str(round(D000112,1))+' & '+str(round(D000222,1))+' & '+str(round(D000122,1))+' & '+str(round(D000022,1))+' & '+str(round(D000012,1))+r' \\')
		file.write('\n')
		file.write(str(round(D110000,1))+' & '+str(round(D110110,1))+' & '+str(round(D110220,1))+' & '+str(round(D110120,1))+' & '+str(round(D110020,1))+' & '+str(round(D110010,1))+' & ')
		file.write('\n')
		file.write(str(round(D110001,1))+' & '+str(round(D110111,1))+' & '+str(round(D110221,1))+' & '+str(round(D110121,1))+' & '+str(round(D110021,1))+' & '+str(round(D110011,1))+' & ')
		file.write('\n')
		file.write(str(round(D110002,1))+' & '+str(round(D110112,1))+' & '+str(round(D110222,1))+' & '+str(round(D110122,1))+' & '+str(round(D110022,1))+' & '+str(round(D110012,1))+r' \\')
		file.write('\n')
		file.write(str(round(D220000,1))+' & '+str(round(D220110,1))+' & '+str(round(D220220,1))+' & '+str(round(D220120,1))+' & '+str(round(D220020,1))+' & '+str(round(D220010,1))+' & ')
		file.write('\n')
		file.write(str(round(D220001,1))+' & '+str(round(D220111,1))+' & '+str(round(D220221,1))+' & '+str(round(D220121,1))+' & '+str(round(D220021,1))+' & '+str(round(D220011,1))+' & ')
		file.write('\n')
		file.write(str(round(D220002,1))+' & '+str(round(D220112,1))+' & '+str(round(D220222,1))+' & '+str(round(D220122,1))+' & '+str(round(D220022,1))+' & '+str(round(D220012,1))+r' \\')
		file.write('\n')
		file.write(str(round(D120000,1))+' & '+str(round(D120110,1))+' & '+str(round(D120220,1))+' & '+str(round(D120120,1))+' & '+str(round(D120020,1))+' & '+str(round(D120010,1))+' & ')
		file.write('\n')
		file.write(str(round(D120001,1))+' & '+str(round(D120111,1))+' & '+str(round(D120221,1))+' & '+str(round(D120121,1))+' & '+str(round(D120021,1))+' & '+str(round(D120011,1))+' & ')
		file.write('\n')
		file.write(str(round(D120002,1))+' & '+str(round(D120112,1))+' & '+str(round(D120222,1))+' & '+str(round(D120122,1))+' & '+str(round(D120022,1))+' & '+str(round(D120012,1))+r' \\')
		file.write('\n')
		file.write(str(round(D020000,1))+' & '+str(round(D020110,1))+' & '+str(round(D020220,1))+' & '+str(round(D020120,1))+' & '+str(round(D020020,1))+' & '+str(round(D020010,1))+' & ')
		file.write('\n')
		file.write(str(round(D020001,1))+' & '+str(round(D020111,1))+' & '+str(round(D020221,1))+' & '+str(round(D020121,1))+' & '+str(round(D020021,1))+' & '+str(round(D020011,1))+' & ')
		file.write('\n')
		file.write(str(round(D020002,1))+' & '+str(round(D020112,1))+' & '+str(round(D020222,1))+' & '+str(round(D020122,1))+' & '+str(round(D020022,1))+' & '+str(round(D020012,1))+r' \\')
		file.write('\n')
		file.write(str(round(D010000,1))+' & '+str(round(D010110,1))+' & '+str(round(D010220,1))+' & '+str(round(D010120,1))+' & '+str(round(D010020,1))+' & '+str(round(D010010,1))+' & ')
		file.write('\n')
		file.write(str(round(D010001,1))+' & '+str(round(D010111,1))+' & '+str(round(D010221,1))+' & '+str(round(D010121,1))+' & '+str(round(D010021,1))+' & '+str(round(D010011,1))+' & ')
		file.write('\n')
		file.write(str(round(D010002,1))+' & '+str(round(D010112,1))+' & '+str(round(D010222,1))+' & '+str(round(D010122,1))+' & '+str(round(D010022,1))+' & '+str(round(D010012,1)))
		file.write('\n')
		file.write(r' \\')
		file.write('\n')
		file.write(str(round(D001000,1))+' & '+str(round(D001110,1))+' & '+str(round(D001220,1))+' & '+str(round(D001120,1))+' & '+str(round(D001020,1))+' & '+str(round(D001010,1))+' & ')
		file.write('\n')
		file.write(str(round(D001001,1))+' & '+str(round(D001111,1))+' & '+str(round(D001221,1))+' & '+str(round(D001121,1))+' & '+str(round(D001021,1))+' & '+str(round(D001011,1))+' & ')
		file.write('\n')
		file.write(str(round(D001002,1))+' & '+str(round(D001112,1))+' & '+str(round(D001222,1))+' & '+str(round(D001122,1))+' & '+str(round(D001022,1))+' & '+str(round(D001012,1))+r' \\')
		file.write('\n')
		file.write(str(round(D111000,1))+' & '+str(round(D111110,1))+' & '+str(round(D111220,1))+' & '+str(round(D111120,1))+' & '+str(round(D111020,1))+' & '+str(round(D111010,1))+' & ')
		file.write('\n')
		file.write(str(round(D111001,1))+' & '+str(round(D111111,1))+' & '+str(round(D111221,1))+' & '+str(round(D111121,1))+' & '+str(round(D111021,1))+' & '+str(round(D111011,1))+' & ')
		file.write('\n')
		file.write(str(round(D111002,1))+' & '+str(round(D111112,1))+' & '+str(round(D111222,1))+' & '+str(round(D111122,1))+' & '+str(round(D111022,1))+' & '+str(round(D111012,1))+r' \\')
		file.write('\n')
		file.write(str(round(D221000,1))+' & '+str(round(D221110,1))+' & '+str(round(D221220,1))+' & '+str(round(D221120,1))+' & '+str(round(D221020,1))+' & '+str(round(D221010,1))+' & ')
		file.write('\n')
		file.write(str(round(D221001,1))+' & '+str(round(D221111,1))+' & '+str(round(D221221,1))+' & '+str(round(D221121,1))+' & '+str(round(D221021,1))+' & '+str(round(D221011,1))+' & ')
		file.write('\n')
		file.write(str(round(D221002,1))+' & '+str(round(D221112,1))+' & '+str(round(D221222,1))+' & '+str(round(D221122,1))+' & '+str(round(D221022,1))+' & '+str(round(D221012,1))+r' \\')
		file.write('\n')
		file.write(str(round(D121000,1))+' & '+str(round(D121110,1))+' & '+str(round(D121220,1))+' & '+str(round(D121120,1))+' & '+str(round(D121020,1))+' & '+str(round(D121010,1))+' & ')
		file.write('\n')
		file.write(str(round(D121001,1))+' & '+str(round(D121111,1))+' & '+str(round(D121221,1))+' & '+str(round(D121121,1))+' & '+str(round(D121021,1))+' & '+str(round(D121011,1))+' & ')
		file.write('\n')
		file.write(str(round(D121002,1))+' & '+str(round(D121112,1))+' & '+str(round(D121222,1))+' & '+str(round(D121122,1))+' & '+str(round(D121022,1))+' & '+str(round(D121012,1))+r' \\')
		file.write('\n')
		file.write(str(round(D021000,1))+' & '+str(round(D021110,1))+' & '+str(round(D021220,1))+' & '+str(round(D021120,1))+' & '+str(round(D021020,1))+' & '+str(round(D021010,1))+' & ')
		file.write('\n')
		file.write(str(round(D021001,1))+' & '+str(round(D021111,1))+' & '+str(round(D021221,1))+' & '+str(round(D021121,1))+' & '+str(round(D021021,1))+' & '+str(round(D021011,1))+' & ')
		file.write('\n')
		file.write(str(round(D021002,1))+' & '+str(round(D021112,1))+' & '+str(round(D021222,1))+' & '+str(round(D021122,1))+' & '+str(round(D021022,1))+' & '+str(round(D021012,1))+r' \\')
		file.write('\n')
		file.write(str(round(D011000,1))+' & '+str(round(D011110,1))+' & '+str(round(D011220,1))+' & '+str(round(D011120,1))+' & '+str(round(D011020,1))+' & '+str(round(D011010,1))+' & ')
		file.write('\n')
		file.write(str(round(D011001,1))+' & '+str(round(D011111,1))+' & '+str(round(D011221,1))+' & '+str(round(D011121,1))+' & '+str(round(D011021,1))+' & '+str(round(D011011,1))+' & ')
		file.write('\n')
		file.write(str(round(D011002,1))+' & '+str(round(D011112,1))+' & '+str(round(D011222,1))+' & '+str(round(D011122,1))+' & '+str(round(D011022,1))+' & '+str(round(D011012,1)))
		file.write('\n')
		file.write(r' \\')
		file.write('\n')
		file.write(str(round(D002000,1))+' & '+str(round(D002110,1))+' & '+str(round(D002220,1))+' & '+str(round(D002120,1))+' & '+str(round(D002020,1))+' & '+str(round(D002010,1))+' & ')
		file.write('\n')
		file.write(str(round(D002001,1))+' & '+str(round(D002111,1))+' & '+str(round(D002221,1))+' & '+str(round(D002121,1))+' & '+str(round(D002021,1))+' & '+str(round(D002011,1))+' & ')
		file.write('\n')
		file.write(str(round(D002002,1))+' & '+str(round(D002112,1))+' & '+str(round(D002222,1))+' & '+str(round(D002122,1))+' & '+str(round(D002022,1))+' & '+str(round(D002012,1))+r' \\')
		file.write('\n')
		file.write(str(round(D112000,1))+' & '+str(round(D112110,1))+' & '+str(round(D112220,1))+' & '+str(round(D112120,1))+' & '+str(round(D112020,1))+' & '+str(round(D112010,1))+' & ')
		file.write('\n')
		file.write(str(round(D112001,1))+' & '+str(round(D112111,1))+' & '+str(round(D112221,1))+' & '+str(round(D112121,1))+' & '+str(round(D112021,1))+' & '+str(round(D112011,1))+' & ')
		file.write('\n')
		file.write(str(round(D112002,1))+' & '+str(round(D112112,1))+' & '+str(round(D112222,1))+' & '+str(round(D112122,1))+' & '+str(round(D112022,1))+' & '+str(round(D112012,1))+r' \\')
		file.write('\n')
		file.write(str(round(D222000,1))+' & '+str(round(D222110,1))+' & '+str(round(D222220,1))+' & '+str(round(D222120,1))+' & '+str(round(D222020,1))+' & '+str(round(D222010,1))+' & ')
		file.write('\n')
		file.write(str(round(D222001,1))+' & '+str(round(D222111,1))+' & '+str(round(D222221,1))+' & '+str(round(D222121,1))+' & '+str(round(D222021,1))+' & '+str(round(D222011,1))+' & ')
		file.write('\n')
		file.write(str(round(D222002,1))+' & '+str(round(D222112,1))+' & '+str(round(D222222,1))+' & '+str(round(D222122,1))+' & '+str(round(D222022,1))+' & '+str(round(D222012,1))+r' \\')
		file.write('\n')
		file.write(str(round(D122000,1))+' & '+str(round(D122110,1))+' & '+str(round(D122220,1))+' & '+str(round(D122120,1))+' & '+str(round(D122020,1))+' & '+str(round(D122010,1))+' & ')
		file.write('\n')
		file.write(str(round(D122001,1))+' & '+str(round(D122111,1))+' & '+str(round(D122221,1))+' & '+str(round(D122121,1))+' & '+str(round(D122021,1))+' & '+str(round(D122011,1))+' & ')
		file.write('\n')
		file.write(str(round(D122002,1))+' & '+str(round(D122112,1))+' & '+str(round(D122222,1))+' & '+str(round(D122122,1))+' & '+str(round(D122022,1))+' & '+str(round(D122012,1))+r' \\')
		file.write('\n')
		file.write(str(round(D022000,1))+' & '+str(round(D022110,1))+' & '+str(round(D022220,1))+' & '+str(round(D022120,1))+' & '+str(round(D022020,1))+' & '+str(round(D022010,1))+' & ')
		file.write('\n')
		file.write(str(round(D022001,1))+' & '+str(round(D022111,1))+' & '+str(round(D022221,1))+' & '+str(round(D022121,1))+' & '+str(round(D022021,1))+' & '+str(round(D022011,1))+' & ')
		file.write('\n')
		file.write(str(round(D022002,1))+' & '+str(round(D022112,1))+' & '+str(round(D022222,1))+' & '+str(round(D022122,1))+' & '+str(round(D022022,1))+' & '+str(round(D022012,1))+r' \\')
		file.write('\n')
		file.write(str(round(D012000,1))+' & '+str(round(D012110,1))+' & '+str(round(D012220,1))+' & '+str(round(D012120,1))+' & '+str(round(D012020,1))+' & '+str(round(D012010,1))+' & ')
		file.write('\n')
		file.write(str(round(D012001,1))+' & '+str(round(D012111,1))+' & '+str(round(D012221,1))+' & '+str(round(D012121,1))+' & '+str(round(D012021,1))+' & '+str(round(D012011,1))+' & ')
		file.write('\n')
		file.write(str(round(D012002,1))+' & '+str(round(D012112,1))+' & '+str(round(D012222,1))+' & '+str(round(D012122,1))+' & '+str(round(D012022,1))+' & '+str(round(D012012,1)))
		file.write('\n')
		file.write(r'\end{pmatrix} $} \text{\,N} \ ,')
		file.write('\n')
		file.write(r'\eqend')



	with open(data+'.tex', 'a') as file:
		file.write(r'\end{document}')
		file.write('\n')
        


# a1=100*abs(float((1-C0000/C1111)))   #c1
# a2=100*abs(float((1-C0011/C0022)))   #c2
# a3=100*abs(float((1-C0202/C0101)))   #c3
# b=100*abs(float((1-D000000/D111111))) #d1
# c=100*abs(float((1-D000110/D111001))) #d2
# d=100*abs(float((1-D000011/D222020))) #d3
# e=100*abs(float((1-D001001/D110110))) #d4##
# f=100*abs(float((1-D221122/D001010))) #d5
# g=100*abs(float((1-D110220/D001221))) #d6
# h=100*abs(float((1-D110022/D001122))) #d7
# i=100*abs(float((1-D022022/D121121))) #d8##
# j=100*abs(float((1-D011022/D020121))) #d9
# k=100*abs(float((1-D120120/D012012))) #d10
# l=100*abs(float((1-D120012/D021012))) #d11

# print(a1,a2,a3,b,c,d,e,f,g,h,i,j,k,l)
#print(C0000,C2222)
# conv=0
with open(sys.argv[3]+'FinalValuesConvList.tex', 'a') as file: #for monitoring
    file.write(str(dofs)+','+str(round(C0000,Sf))+','+str(round(C0011,Sf))+','+str(round(C1212,Sf))
                +','+str(round(D000000,Sf))+','+str(round(D000110,Sf))+','+str(round(D000011,Sf))
                +','+str(round(D110110,Sf))+','+str(round(D110011,Sf))+','+str(round(D110220,Sf))
                +','+str(round(D110022,Sf))+','+str(round(D011011,Sf))+','+str(round(D011022,Sf))
                +','+str(round(D120120,Sf))+','+str(round(D120021,Sf))+','+sys.argv[1]+'\n')

CurrentMaxValuesList=[C0000, C0011, C1212, float(D000000), float(D000110), float(D000011)
                      , float(D110110), float(D110011), float(D110220), float(D110022), float(D011011)
                      , float(D011022), float(D120120), float(D120021)]
    
if os.path.isfile("ForCheckingConvergence.txt"): #not in first time
    with open("ForCheckingConvergence.txt", "r") as fp:         
        LastMaxValuesList = json.load(fp)
        
with open("ForCheckingConvergence.txt", "w") as fp:  
    json.dump(CurrentMaxValuesList, fp)
    
if 'LastMaxValuesList' in globals(): #not in first time
    for indexes in range(0,14):  
        print('for', indexes+1, 'last is ', LastMaxValuesList[indexes], ' and current is ',CurrentMaxValuesList[indexes])
        print(indexes+1 ,' = ', abs((CurrentMaxValuesList[indexes]-LastMaxValuesList[indexes])/CurrentMaxValuesList[indexes])*100)
    for indexes in range(0,14):  
        if abs((CurrentMaxValuesList[indexes]-LastMaxValuesList[indexes])/CurrentMaxValuesList[indexes])*100 >= float(sys.argv[7]):
            print(indexes+1 , 'is not good')
            conv = 0
            sys.exit(conv)
else:   #for first time
    conv = 0
    sys.exit(conv)
   
conv = 1
with open('CsvFinalValues.tex', 'a') as file:
    file.write(str(round(C0000,Sf))+','+str(round(C0011,Sf))+','+str(round(C1212,Sf))
                +','+str(round(D000000,Sf))+','+str(round(D000110,Sf))+','+str(round(D000011,Sf))
                +','+str(round(D110110,Sf))+','+str(round(D110011,Sf))+','+str(round(D110220,Sf))
                +','+str(round(D110022,Sf))+','+str(round(D011011,Sf))+','+str(round(D011022,Sf))
                +','+str(round(D120120,Sf))+','+str(round(D120021,Sf))
                +','+str(sys.argv[4])+','+str(float(sys.argv[5]))+'\n')
sys.exit(conv)   #exit code

#if a1<1 and a2<1 and a3<1 and b<2 and c<2 and d<2 and e<2 and f<2 and g<2 and h<2 and i<2 and j<2 and k<2 and l<2:
# if a1<0.1 and a2<0.1 and a3<0.1 and b<0.1 and c<0.1 and d<0.1 and e<0.1 and f<0.1 and g<0.1 and h<0.1 and i<0.1 and j<0.1 and k<0.1 and l<0.1:
# #if True: 
#     conv=1
#     with open(sys.argv[3]+'FinalValues.tex', 'w') as file:
#         file.write(str(round(C0000,Sf))+','+str(round(C0011,Sf))+','+str(round(C1212,Sf))
#                    +','+str(round(D000000,Sf))+','+str(round(D000110,Sf))+','+str(round(D000011,Sf))
#                    +','+str(round(D110110,Sf))+','+str(round(D110011,Sf))+','+str(round(D110220,Sf))
#                    +','+str(round(D110022,Sf))+','+str(round(D011011,Sf))+','+str(round(D011022,Sf))
#                    +','+str(round(D120120,Sf))+','+str(round(D120021,Sf))+'\nDOF = '+str(dofs))
#     with open('CsvFinalValues.tex', 'a') as file:
#         file.write(str(round(C0000,Sf))+','+str(round(C0011,Sf))+','+str(round(C1212,Sf))
#                    +','+str(round(D000000,Sf))+','+str(round(D000110,Sf))+','+str(round(D000011,Sf))
#                    +','+str(round(D110110,Sf))+','+str(round(D110011,Sf))+','+str(round(D110220,Sf))
#                    +','+str(round(D110022,Sf))+','+str(round(D011011,Sf))+','+str(round(D011022,Sf))
#                    +','+str(round(D120120,Sf))+','+str(round(D120021,Sf))
#                    +','+str(sys.argv[4])+','+str(float(sys.argv[5]))+'\n')
          
# else:
#     conv=0
    # with open(sys.argv[1]+'FinalValues.tex', 'w') as file:
    #     file.write(str(round(C0000,1))+','+str(round(C0011,1))+','+str(round(C1212,1))
    #                +','+str(round(D000000,1))+','+str(round(D000110,1))+','+str(round(D000011,1))
    #                +','+str(round(D110110,1))+','+str(round(D110011,1))+','+str(round(D110220,1))
    #                +','+str(round(D110022,1))+','+str(round(D011011,1))+','+str(round(D011022,1))
    #                +','+str(round(D120120,1))+','+str(round(D120021,1))+'\nDOF = '+str(dofs))


# sys.exit(conv)   #exit code

   
#why C?
# D_voigt = numpy.array ( [\
#     [ C0000,     C0011,      C0022,     C0012,     C0002,      C0001, 0., 0., 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., ],\
#     [ C1100,     C1111,      C1122,     C1112,     C1102,      C1101, 0., 0., 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., ],\
#     [ C2200,     C2211,      C2222,     C2212,     C2202,      C2201, 0., 0., 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., ],\
#     [ C1200,     C1211,      C1222,     C1212,     C1202,      C1201, 0., 0., 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., ],\
#     [ C0200,     C0211,      C0222,     C0212,     C0202,      C0201, 0., 0., 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., ],\
#     [ C0100,     C0111,      C0122,     C0112,     C0102,      C0101, 0., 0., 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., ],\
#     [ 0.,0.,0.,0.,0.,0.,D000000, D000110, D000011, D000220, D000022, D000111, D000001, D000010, D000221, D000122, D000222, D000002, D000020, D000112, D000121, D000120, D000021, D000012 ],\
#     [ 0.,0.,0.,0.,0.,0.,D110000, D110110, D110011, D110220, D110022, D110111, D110001, D110010, D110221, D110122, D110222, D110002, D110020, D110112, D110121, D110120, D110021, D110012 ],\
#     [ 0.,0.,0.,0.,0.,0.,D011000, D011110, D011011, D011220, D011022, D011111, D011001, D011010, D011221, D011122, D011222, D011002, D011020, D011112, D011121, D011120, D011021, D011012 ],\
#     [ 0.,0.,0.,0.,0.,0.,D220000, D220110, D220011, D220220, D220022, D220111, D220001, D220010, D220221, D220122, D220222, D220002, D220020, D220112, D220121, D220120, D220021, D220012 ],\
#     [ 0.,0.,0.,0.,0.,0.,D022000, D022110, D022011, D022220, D022022, D022111, D022001, D022010, D022221, D022122, D022222, D022002, D022020, D022112, D022121, D022120, D022021, D022012 ],\
#     [ 0.,0.,0.,0.,0.,0.,D111000, D111110, D111011, D111220, D111022, D111111, D111001, D111010, D111221, D111122, D111222, D111002, D111020, D111112, D111121, D111120, D111021, D111012 ],\
#     [ 0.,0.,0.,0.,0.,0.,D001000, D001110, D001011, D001220, D001022, D001111, D001001, D001010, D001221, D001122, D001222, D001002, D001020, D001112, D001121, D001120, D001021, D001012 ],\
#     [ 0.,0.,0.,0.,0.,0.,D010000, D010110, D010011, D010220, D010022, D010111, D010001, D010010, D010221, D010122, D010222, D010002, D010020, D010112, D010121, D010120, D010021, D010012 ],\
#     [ 0.,0.,0.,0.,0.,0.,D221000, D221110, D221011, D221220, D221022, D221111, D221001, D221010, D221221, D221122, D221222, D221002, D221020, D221112, D221121, D221120, D221021, D221012 ],\
#     [ 0.,0.,0.,0.,0.,0.,D122000, D122110, D122011, D122220, D122022, D122111, D122001, D122010, D122221, D122122, D122222, D122002, D122020, D122112, D122121, D122120, D122021, D122012 ],\
#     [ 0.,0.,0.,0.,0.,0.,D222000, D222110, D222011, D222220, D222022, D222111, D222001, D222010, D222221, D222122, D222222, D222002, D222020, D222112, D222121, D222120, D222021, D222012 ],\
#     [ 0.,0.,0.,0.,0.,0.,D002000, D002110, D002011, D002220, D002022, D002111, D002001, D002010, D002221, D002122, D002222, D002002, D002020, D002112, D002121, D002120, D002021, D002012 ],\
#     [ 0.,0.,0.,0.,0.,0.,D020000, D020110, D020011, D020220, D020022, D020111, D020001, D020010, D020221, D020122, D020222, D020002, D020020, D020112, D020121, D020120, D020021, D020012 ],\
#     [ 0.,0.,0.,0.,0.,0.,D112000, D112110, D112011, D112220, D112022, D112111, D112001, D112010, D112221, D112122, D112222, D112002, D112020, D112112, D112121, D112120, D112021, D112012 ],\
#     [ 0.,0.,0.,0.,0.,0.,D121000, D121110, D121011, D121220, D121022, D121111, D121001, D121010, D121221, D121122, D121222, D121002, D121020, D121112, D121121, D121120, D121021, D121012 ],\
#     [ 0.,0.,0.,0.,0.,0.,D120000, D120110, D120011, D120220, D120022, D120111, D120001, D120010, D120221, D120122, D120222, D120002, D120020, D120112, D120121, D120120, D120021, D120012 ],\
#     [ 0.,0.,0.,0.,0.,0.,D021000, D021110, D021011, D021220, D021022, D021111, D021001, D021010, D021221, D021122, D021222, D021002, D021020, D021112, D021121, D021120, D021021, D021012 ],\
#     [0.,0.,0.,0.,0.,0., D012000, D012110, D012011, D012220, D012022, D012111, D012001, D012010, D012221, D012122, D012222, D012002, D012020, D012112, D012121, D012120, D012021, D012012 ],\
#  ], dtype='float')

# def isPSD(A, tol=1e-8):
#   E = numpy.linalg.eigvalsh(A)
#   return numpy.all(E > -tol)
# a = isPSD(D_voigt,  1e-8)
# print(a)


