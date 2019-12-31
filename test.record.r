library(fda)
rangeval=c(0,1)
nbasis=5
# basisobj=create.monomial.basis(rangeval,nbasis)
#order 4
norder=4
#degree 3
#knots 9 (0,0,0,0,0.5,1,1,1,1) boundary point #knots=#order (typically)
breaks=c(0,0.5,1)
#breakpoints 3, 1 interior breakpoint
#single nodes per breakpoints
# basis functions = order + #interior knots 5
basisobj=create.bspline.basis(rangeval,nbasis,norder,breaks)
splinebasis=create.bspline.basis(c(0,10),13)
plot(splinebasis,xlab='t',ylab='Bspline basis functions B(t)',las=1,lwd=2)
# to make use of nth order derivatives, the order of spline basis need at least to be n+2
# evaluate the basis
eval.basis(1:0.1:10,splinebasis)
eval.basis(1:0.1:10,splinebasis,1)
##function creating
tempfd=fd(coefmat,daybasis65)
thatvec=eval.fd(tvec,thawfd)
D2thatvec=eval.fd(tvec,thawfd,2)
# linear differential operator
harmaccelLfd=vec2Lfd(c(0,omegaˆ2,0),c(0,365))
D2tempfd=deriv.fd(tempfd,2)
heightcoef=lsfit(basismat,heightmat,intercept=FALSE)$coef
# another way by smooth
heightList=smooth.basis(age,heightmat,heightbasis12)
##smoothing splines
#Rmatrix
Rmat=eval.penalty(tempbasis,harmaccelLfd)
#
norder=6
nbasis=length(age)+norder-2
heightbasis=create.bspline.basis(c(1,18),nbasis,norder,age)
heightfdPar=fdPar(heightbasis,4,0.01)
heightfd=smooth.basis(age,heightmat,heightfdPar)$fd
#constraint: positive
VanPrecPos=smooth.pos(day.5,VanPrec,WfdParobj)
#constraint: monotone
result=smooth.monotone(day,tib,WfdPar)
#constraint: probability density function
densityList=density.fd(RegPrec,WfdPar)
#statistics
meanlogprec=mean(logprec.fd)
stddevlogprec=std.fd(logprec.fd)
logprecvar.bifd=var.fd(logprec.fd)#eval.bifd
