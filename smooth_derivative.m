function [returnstr]=smooth_derivative(Ymat,xvec,loglambda_vec,nDer,plotflag)
% a wrap up for smoothing spline in Functional Data Analysis (FDA)
% spline on both the original function and its derivatives are supported
% B-spline will be used, and lambda value will be searched based on GCV (generalized cross-validation)
% Argument:
%%          Ymat: 2d numeric matrix. The value to be smoothed. Each column is different curve and each row is the different index in the curve (such as time), represented by xvec. If there is missing value, deal with it before using this function. Must be provided. Ymat can also be 3d, index * replicate * variables
%%          xvec: 1d numeric array corresponding to the rows in Ymat. Ex: time sequence. Must be provided. Carefully choose the unit, so that the subinterval is not too far from 1.
%%          loglambda_vec: 1d numeric array. the grid for log10 lambda. Default -6:1:2. The user is highly recommended to try different values.
%%          nDer: int. the highest derivative to be analzed in following workflow. Default 0.
%%          plotflag: bool. If true plot the GCV curve. If false do not plot
% Return:   returnstr: struct.
% CITATION:
%         Ramsay, James & Hooker, Giles & Graves, Spencer. (2009). Functional data analysis with R and MATLAB. 10.1007/978-0-387-98185-7.
% More information:
%     http://www.psych.mcgill.ca/misc/fda/
%     https://cran.r-project.org/web/packages/fda/index.html
%     http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/
%
% YUE WU 12312019

if ~exist('Ymat','var')
  error('Ymat is needed');
end
if ~exist('xvec','var')
  error('xvec is needed');
end
if ~exist('loglambda_vec','var')
  loglambda_vec=-6:1:2;
end
if ~exist('nDer','var')
  nDer=0;
end
if ~exist('plotflag','var')
  plotflag=false;
end
% smoothing penalty is based on second derivative
nsmooth=nDer+2;
% to make use of nth order derivatives, the order of spline basis need at least to be n+2
norder=nsmooth+2;
% basis functions = order + #interior knots
nbasis=length(xvec)+norder-2;
bbasis=create_bspline_basis([min(xvec) max(xvec)],nbasis,norder,xvec);
nlambda=length(loglambda_vec);
gcvsave=zeros(nlambda,1);
dfsave=gcvsave;
spfdcell=cell(nlambda,1);
for i=1:nlambda
   lambda=10^loglambda_vec(i);
   fdPari=fdPar(bbasis,int2Lfd(nsmooth),lambda);
   [fdres,dfi,gcvi]=smooth_basis(xvec,Ymat,fdPari);
   gcvsave(i)=sum(gcvi(:));
   dfsave(i)=dfi;
   spfdcell{i}=fdres;
end
if plotflag
  h=figure()
    plot(loglambda_vec,gcvsave,'k-o');
    ylabel('GCV');
    xlabel('loglambda');
end

[gcv_select,minind]=min(gcvsave);
returnstr=struct();
returnstr.dfsave=dfsave;
returnstr.gcvsave=gcvsave;
returnstr.gcv_select=gcv_select;
returnstr.df_selec=dfsave(minind);
returnstr.spfd_selec=spfdcell{minind};
returnstr.loglambda=loglambda_vec(minind);
end
