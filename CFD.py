import numpy as np
import matplotlib.pyplot as plt
x=201
y=11
p=(x-1)*(y-1)
x1=np.linspace(0,20,x)
y1=np.linspace(0,1,y)
dx=.1
dy=.1
pr=7
k=0
dt=.0001
ra=3000
old_w=np.zeros((x,y))
old_temp=np.zeros((x,y))
old_s=np.zeros((x,y))
u=np.zeros((x,y))
v=np.zeros((x,y))
a=np.zeros((x,y))
b=np.zeros((x,y))
c1=np.zeros((x,y))
c=dt/dx
eps=10**-7

o=1
it=0
rms=1
W=[]
S=[]
U=[]
V=[]
T=[]

initial_w = np.zeros((x, y))
initial_s = np.zeros((x, y))
initial_temp = np.zeros((x, y))

W.append(old_w)
T.append(old_temp)
U.append(u)
V.append(v)
S.append(old_s)
iter=0
while it>-1:

    while iter>-1:
        # BOUNDARY CONDITIONS

        j = 0
        for i in range(0, x):
            initial_s[i, j] = k
            initial_w[i, j] = 2 * (initial_s[i, j] - initial_s[i, j + 1]) / dy ** 2
            initial_temp[i, j] = 1

        j = y - 1
        for i in range(0, x):
            initial_s[i, j] = k
            initial_w[i, j] = 2 * (initial_s[i, j] - initial_s[i, j - 1]) / dy ** 2
            initial_temp[i, j] = 0

        i = 0
        for j in range(1, y - 1):
            initial_s[i, j] = k
            initial_w[i, j] = 2 * (initial_s[i, j] - initial_s[i + 1, j]) / dx ** 2
            initial_temp[i, j] = initial_temp[i + 1, j]

        i = x - 1
        for j in range(1, y - 1):
            initial_s[i, j] = k
            initial_w[i, j] = 2 * (initial_s[i, j] - initial_s[i - 1, j]) / dx ** 2
            initial_temp[i, j] =    initial_temp[i - 1, j]

        err=0

        for i in range(1,x-1):
            for j in range(1,y-1):
                a[i,j-1]=-c*(v[i,j-1]/2+pr/dy)
                a[i-1,j]=-c*(u[i-1,j]/2+pr/dx)
                a[i,j]=1+2*pr*dt*(1/dx**2+1/dy**2)
                a[i+1,j]=c*(u[i+1,j]/2-pr/dx)
                a[i,j+1]=c*(v[i,j+1]/2-pr/dy)

                b_w=old_w[i,j]+ra*pr*c/2*(initial_temp[i+1,j]-initial_temp[i-1,j])
                au_w=a[i,j-1]*initial_w[i,j-1]+a[i-1,j]*initial_w[i-1,j]+a[i,j]*initial_w[i,j]+a[i+1,j]*initial_w[i+1,j]+a[i,j+1]*initial_w[i,j+1]
                cor_w=b_w-au_w
                res=cor_w
                err = err + res + res
                initial_w[i,j]=initial_w[i,j]+o*cor_w/a[i,j]




                b[i,j-1]=1/dy**2
                b[i-1,j]=1/dx**2
                b[i,j]=-2*(1/dx**2+1/dy**2)
                b[i+1,j]=1/dx**2
                b[i,j+1]=1/dy**2

                b_s=-initial_w[i,j]
                au_s=b[i,j-1]*initial_s[i,j-1]+b[i-1,j]*initial_s[i-1,j]+b[i,j]*initial_s[i,j]+b[i+1,j]*initial_s[i+1,j]+b[i,j+1]*initial_s[i,j+1]
                cor_s=b_s-au_s
                res=cor_s
                # err = err + res + res
                initial_s[i,j]=initial_s[i,j]+o*cor_s/b[i,j]

                u[i,j]=(initial_s[i,j+1]-initial_s[i,j-1])/(2*dy)
                v[i,j]=(initial_s[i-1,j]-initial_s[i+1,j])/(2*dx)


                c1[i, j - 1] = -c * (v[i, j - 1] / 2 + 1 / dy)
                c1[i - 1, j] = -c * (u[i - 1, j] / 2 + 1 / dx)
                c1[i, j] = 1 + 2 *  dt * (1 / dx ** 2 + 1 / dy ** 2)
                c1[i + 1, j] = c * (u[i + 1, j] / 2 - 1 / dx)
                c1[i, j + 1] = c * (v[i, j + 1] / 2 - 1 / dy)

                b_t=old_temp[i,j]
                au_t=c1[i,j-1]*initial_temp[i,j-1]+c1[i-1,j]*initial_temp[i-1,j]+c1[i,j]*initial_temp[i,j]+c1[i+1,j]*initial_temp[i+1,j]+c1[i,j+1]*initial_temp[i,j+1]
                cor_t=b_t-au_t
                initial_temp[i, j] = initial_temp[i, j] + o * cor_t / c1[i, j]
                res=cor_t
                # err=err+res+res


        sq=err/p
        rms=sq

        if rms<eps:
            old_temp=initial_temp.copy()
            old_w=initial_w.copy()
            break

    T.append(old_temp)
    S.append(initial_s)
    U.append(u)
    V.append(v)
    W.append(old_w)
    err1=0
    for i in range(0,x):
        for j in range(0,y):
            err1=err1+(T[it+1][i,j]-T[it][i,j])**2
    rms1=(err1/(x*y))**.5
    print('rms1 = ',rms1)
    if rms1<eps:
        break
    it=it+1
    print('it = ',it)




print(T[it][100][3],'t')
TEMP=np.transpose(T[it-1])
sci=np.transpose(S[it-1])
omega=np.transpose(W[it-1])
t1=np.zeros(201)
u1=np.zeros(201)
for i in range(0,x):
    t1[i]=T[it][i][5]
    u1[i]=U[it][i][5]
print(t1)

fig=plt.figure()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Temperature plot')
plt.contourf(x1,y1,TEMP)
cp=plt.contourf(x1,y1,TEMP)
fig.colorbar(cp)
plt.show()

fig=plt.figure()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Stream function  plot')
plt.contourf(x1,y1,sci)
cp=plt.contourf(x1,y1,sci)
fig.colorbar(cp)
plt.show()

fig=plt.figure()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Vorticity plot')
plt.contourf(x1,y1,omega)
cp=plt.contourf(x1,y1,omega)
fig.colorbar(cp)
plt.show()

plt.xlabel('y')
plt.ylabel('U*')
plt.title('Velocity profile \n X=10L ')
plt.plot(y1,U[it][100][:])
plt.show()


plt.xlabel('y')
plt.ylabel('Temperature')
plt.title('Temperature profile \n X=10L ')
plt.plot(y1,T[it][100][:])
plt.show()

plt.xlabel('X')
plt.ylabel('TEMPERATURE')
plt.title('TEMPERATURE PROFILE \n Y=L/2')
plt.plot(x1,t1)
plt.show()

plt.xlabel('X')
plt.ylabel('VELOCITY')
plt.title('VELOCITY PROFILE \n Y=L/2')
plt.plot(x1,u1)
plt.show()












