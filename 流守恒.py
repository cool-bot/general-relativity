#sample code: calculation of connection, Riemann tensor and Gaussian curvature (if dim=2)
#----------by Zhiqi Huang--for the course ''General Relativity''-------------

import sympy as sym
#optmize printing
sym.init_printing()

#dimension of the space
dim = 4

#define u[0],u[1],... in a single line
u = sym.symarray('u',dim)

A = sym.symbols('A', cls = sym.Function)
B = sym.symbols('B', cls = sym.Function)
P = sym.symbols('P', cls = sym.Function)
rho = sym.symbols('rho', cls = sym.Function)
t, r, theta, phi = sym.symbols(r't, r, theta, phi')

#----define covariant metric
r = sym.symbols('r')
gdown = sym.diag(-sym.exp(2*A(u[1])), sym.exp(2*B(u[1])), u[1]**2, u[1]**2*(sym.sin(u[2]))**2)
EM_TENSOR_up = sym.diag(rho(u[1])*sym.exp(-2*A(u[1])), P(u[1])*sym.exp(-2*B(u[1])), P(u[1])*u[1]**-2, P(u[1])*u[1]**-2*(sym.sin(u[2]))**-2)
#gdown = sym.diag( 1, sym.sin(u[0])**2 )
#gdown = sym.diag(1, sym.exp(2*u[0]))
#gdown = sym.diag(1/(1+u[0]**2+u[1]**2),1/(1+u[0]**2+u[1]**2))
#gdown = sym.diag(1/(1+u[0]**2+u[1]**2)**2,1/(1+u[0]**2+u[1]**2)**2)
#gdown = sym.Matrix([[1/(u[0]**2+u[1]**2+1), 1/(u[0]**2+u[1]**2+1)/2],[1/(u[0]**2+u[1]**2+1)/2, 1/(u[0]**2+u[1]**2+1)]])

#contravariant metric
gup = gdown ** -1

#determinant of the covariant metric
detg = gdown.det()

#\Gamma_{ijk}
def connection_down(i, j, k):
    return (sym.diff(gdown[j, i], u[k]) + sym.diff(gdown[i, k], u[j]) - sym.diff(gdown[j, k], u[i]))/2

#\Gamma^i_{\ jk}
def connection_up(i, j, k):
    gam = 0
    for l in range(dim):
        gam += connection_down(l, j, k) * gup[l, i]
    return sym.simplify(gam)

#R_{ijkl}
def Riemann_tensor_down(i, j, k, l):
    R = sym.diff(connection_down(i, j, k), u[l]) - sym.diff(connection_down(i, j, l), u[k])
    for m in range(dim):
        R += connection_down(m, i, k) * connection_up(m, j, l) - connection_down(m, i, l) * connection_up(m, j, k)
    return sym.simplify(R)

def Riemann_tensor_up(i, j, k, l):
    R = sym.diff(connection_up(i, j, k), u[l]) - sym.diff(connection_up(i, j, l), u[k])
    for m in range(dim):
        R += connection_up(i, m, l) * connection_up(m, j, k) - connection_up(i, m, k) * connection_up(m, j, l)
    return sym.simplify(R)

def J_flow(i):
    J1 = J2 = 0
    for m in range(dim):
        J1 += (sym.diff(EM_TENSOR_up[m, i], u[m]))   
    for a in range(dim):
        for b in range(dim):
            J2 +=  EM_TENSOR_up[a, b] * connection_up(i,a,b) + EM_TENSOR_up[i,a] * connection_up(b,a,b)
    return sym.simplify(J1+J2)


#R_{ij}
def Ricci_tensor_down(i, j):
    R = 0
    for k in range(dim):
        for l in range(k+1):
            if(gup[k, l] != 0):
                if(k==l):
                    R += Riemann_tensor_down(k, i, j, l)*gup[k,l]
                else:
                    R += 2*Riemann_tensor_down(k, i, j, l)*gup[k,l]
    return sym.simplify(R)

#R
def Ricci_scalar():
    R = 0
    for k in range(dim):
        for l in range(k+1):
            if(gup[k, l] != 0):
                if(k==l):
                    R += Ricci_tensor_down(k, l)*gup[k,l]
                else:
                    R += 2*Ricci_tensor_down(k, l)*gup[k,l]    
    return sym.simplify(R)

##this function only works for dim=2
def Gaussian_curvature():
    return sym.simplify(Riemann_tensor_down(0, 1, 1, 0)/detg)


print("---------J_flow--------")
for i in range(dim):
    J = J_flow(i)
    if(J != 0):
       print('J_%d = '%(i)+str(J.subs(u[0], t).subs(u[1], r).subs(u[2], theta).subs(u[3], phi)))
       print(sym.latex(J.subs(u[0], t).subs(u[1], r).subs(u[2], theta).subs(u[3], phi)))


       
