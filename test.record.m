% add the directory /Users/mikeaalv/Downloads/fdaM/ into the path
basisobj=create_bspline_basis(rangeval,nbasis,norder,breaks);
basisobj=create_monomial_basis(rangeval,nbasis);
basismatrix=eval_basis(tvec,mybasis)
Dbasismatrix=eval_basis(tvec,mybasis,1)
% SSE can be used to estimate coef as in linear model
heightcoef=basismat\heightmat
% another way by smooth
fdobj,df,gcv,coef,SSE,penmat,y2cMap]=smooth_basis(age, heightmat,heightbasis12);
% smoothing spline GCV
loglam= -6:0.25:0;
gcvsave=zeros(length(loglam),1);
dfsave=gcvsave;
for i=1:length(loglam)
   lambdai=10ˆloglam(i);
   hgtfdPari=fdPar(heightbasis,4,lambdai);
   [hgtfdi,dfi,gcvi]=smooth_basis(age,hgtfmat,hgtfdPari);
   gcvsave(i)=sum(gcvi);
   dfsave(i)=dfi;
end
