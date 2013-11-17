import scipy.optimize as so

#x = array([ 0.45,  0.5 ,  0.55,  0.6 ,  0.65,  0.7 ,  0.75,  0.8 ,  0.85,
#        0.9 ,  0.95,  1.  ,  1.05,  1.1 ])
#y = array([136, 309, 490, 758, 926, 885, 810, 660, 473, 364, 243, 225, 203,
#    152])

gauss = lambda x, A, sig, m: A*exp(-(x-m)**2/(2*sig**2))

residuals = lambda p,x,y, model: y - model(x,*p)
center_of_mass = lambda x,y: sum(x*y)/sum(y)

def gaussian_fit(x,y):
    A0 = max(x)
    m0 = center_of_mass(x,y)
    sig0 = (x.max()-x.min())/4.
    p0 = array([A0, sig0, m0])
    p,_ = so.leastsq(residuals,p0, args=(x,y,gauss))
    return p


def LAR_residuals(p, x,y,model, alpha=0.1):
    return sqrt(abs(y - model(x,*p)))

def bisq_residuals(p, x,y,model, alpha=5):
    return bisquare_weights((y - model(x,*p))/y)

def bisquare_weights(u):
    w = zeros(u.size)
    w[abs(u) < 1] = (1-u[abs(u) < 1]**2)
    return w

def bisq_gaussian_fit(x,y):
    A0 = max(x)
    m0 = center_of_mass(x,y)
    sig0 = (x.max()-x.min())/4.
    p0 = array([A0, sig0, m0])
    p,_ = so.leastsq(bisq_residuals,p0, args=(x,y,gauss))
    return p


# fname = '2011-06-15/DV_8x1_Dicr_xFRET5b_100pM_ND0_0.dat'
eff0 = d.eff[0]
dir='figure/fitting/'

fit_fun = bisq_gaussian_fit
#fit_fun = gaussian_fit

figure()
step = 0.025
h,x,_=hist(eff0,bins=arange(-0.2,1.2,step))
th=0.45; xp=x[x>th][:-1]+step/2; yp=h[x[:-1]>th]
plot(xp,yp,'-o')
p = fit_fun(xp,yp)
print p
xx = arange(0,1.2,0.005)
plot(xx,gauss(xx,p[0],p[1],p[2]), lw=2)
title('Gaussian fitting of x > %1.2f' % th)
#savefig(dir+'gaussian_fitting_th.png')

figure()
h,x,_=hist(eff0,bins=arange(-0.2,1.2,step))
th=0.45; xi=x[x>th][:-1]+step/2; yi=h[x[:-1]>th]
yth = 0.6; 
n=find(yi>yth*yi.max())[0]
yp =yi[n:]; xp=xi[n:]
#yp =yi[yi>0.5*yi.max()]; xp=xi[yi>0.5*yi.max()]
plot(xp,yp,'-o')
p = fit_fun(xp,yp)
print p
xx = arange(0,1.2,0.005)
plot(xx,gauss(xx,p[0],p[1],p[2]), lw=2)
title('Gaussian fitting of y > %1.2f' % yth)
#savefig(dir+'gaussian_fitting_yth.png')

#figure()
#deg=5
#h,x,_=hist(eff0,bins=arange(-0.2,1.2,step))
#th=0.45; xi=x[x>th][:-1]+step/2; yi=h[x[:-1]>th]
#yth = 0.5; yp =yi[yi>0.5*yi.max()]; xp=xi[yi>0.5*yi.max()]
#plot(xp,yp,'-o')
#p = polyfit(xp,yp,deg)
#print p
#xx = arange(th,1.2,0.005)
#y1,y2 = ylim()
#plot(xx,polyval(p,xx), lw=2)
#ylim(y1,y2)
#title('Polynomial fitting of %d degree' % deg)
#savefig(dir+'polynomial_fitting_th.png')

#from saviztky_golay import savitzky_golay
#figure()
#step = 0.025
#h,x,_=hist(eff0,bins=arange(-0.2,1.2,step))
#th=0.45; xp=x[x>th][:-1]+step/2; yp=h[x[:-1]>th]
#plot(xp,yp,'-o')
#ys = savitzky_golay(yp)
#plot(xp,ys,'-', lw=2)
#title('Saviztky_golay smoothing')
#savefig(dir+'savitzky_golay_smoothing.png')

