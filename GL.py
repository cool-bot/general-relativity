#sample code: calculation of connection, Riemann tensor and Gaussian curvature (if dim=2)
#----------by Zhiqi Huang--for the course ''General Relativity''-------------

import sympy as sym
#optmize printing
sym.init_printing()

#dimension of the space
dim = 4

#define u[0],u[1],... in a single line
u = sym.symarray('u',dim)

f = sym.symbols('f', cls = sym.Function)
g = sym.symbols('g', cls = sym.Function)
t, r, theta, phi = sym.symbols(r't, r, theta, phi')

#----define covariant metric
r = sym.symbols('r')
gdown = sym.diag(-f(u[1]), 1/g(u[1]), u[1]**2, u[1]**2*(sym.sin(u[2]))**2)
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



print("---------Riemann_tensor_up------")
for i in range(dim):
    for j in range(dim):
        for k in range(dim):
            for l in range(dim):
                R = Riemann_tensor_up(i, j, k, l)
                if(R != 0 ):
                    print('R%d_%d%d%d = '%(i, j, k, l)+str(R.subs(u[0], t).subs(u[1], r).subs(u[2], theta).subs(u[3], phi)))



'''
print("---------Gamma down--------")
for i in range(dim):
    for j in range(dim):
        for k in range(dim):
            gam = connection_down(i, k, j)
            if(gam != 0):
                print(i, k, j, gam)

print("---------Gamma up--------")
for i in range(dim):
    for j in range(dim):
        for k in range(dim):
            gam = connection_up(i, k, j)
            if(gam != 0):
                print(i, k, j, gam)


print("---------Ricci tensor down------")
for i in range(dim):
    for j in range(dim):
        R = Ricci_tensor_down(i, j)
        if(R != 0):
            print(i, j, R)
print("---------Ricci scalar------")
print(Ricci_scalar())
'''

if(dim==2):
    print("--------Gaussian Curvature-------------")
    K= Gaussian_curvature()
    print(K)
    print("--------Gaussian Curvature at (0, 0)----")    
    print(K.subs(u[0],0).subs(u[1],0))
