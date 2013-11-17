#nx = 10
#ny = 7
#X,Y = mgrid[0:1:nx*1j,0:1:ny*1j]
#Z = randn(nx,ny)

#X = Eg
#Y = Sg
Z = Hg.astype('float')
#contourf(X,Y,Z)
#show()

newX, newY = mgrid[X.min():X.max():0.01,Y.min():Y.max():0.01]

spl = I.bisplrep(X,Y,Z, s=0.01)
newZ = I.bisplev(newX[:,0], newY[0,:], spl)

contourf(newX,newY,newZ)
show()
