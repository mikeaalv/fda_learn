% example 1
%% test on the data set Registered Handwriting Data
%% make sure the fda matlab toolbox is downloaded and added to the path

%%TEST CODE FROM FDA matlab examples (http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/)
% close all;
% clear all;
% load fda;
% nbasis92=92;
% norder92=6;
% fdabasis92=create_bspline_basis(fdarange,nbasis92,norder92);
% fdafd92=smooth_basis(fdatime,fdaarray,fdabasis92);
% %  Add suitable names for the dimensions of the data.
% fdafd92_fdnames{1}='Milliseconds';
% fdafd92_fdnames{2}='Replications';
% fdafd92_fdnames{3}='Metres';
% fdafd92=putnames(fdafd92,fdafd92_fdnames);
% lambda=0;
% fdPar92=fdPar(fdabasis92,4,lambda);
% nharm=3;
% fdapcastr=pca_fd(fdafd92,nharm,fdPar92);
% %  plot unrotated harmonics
% subplot(1,1,1)
% plot_pca_fd(fdapcastr)
% %  Varimax rotation
% fdarotpcastr=varmx_pca(fdapcastr);
% %  plot rotated harmonics
% subplot(1,1,1)
% plot_pca_fd(fdarotpcastr)
% %  plot log eigenvalues
% fdaharmeigval=fdapcastr.values;
% plot(1:19,fdaharmeigval(1:19),'-o')
% xlabel('Eigenvalue Number')
% ylabel('Eigenvalue')
% %  plot factor scores
% harmscr=fdarotpcastr.harmscr;
% plot(harmscr(:,1),harmscr(:,2),'o')
% xlabel('pc1')
% ylabel('pc2')
% testrec=struct();
% testrec.fdapcastr=fdapcastr;
% testrec.fdarotpcastr=fdarotpcastr;
% clearvars -except testrec;

%%New code using functions (2d data)
load fda;
nbasis92=92;
norder92=6;
fdabasis92=create_bspline_basis(fdarange,nbasis92,norder92);
fdafd92=smooth_basis(fdatime,fdaarray(:,:,1),fdabasis92);%fdaarray(:,:,1)
%  Add suitable names for the dimensions of the data.
fdafd92_fdnames{1}='Milliseconds';
fdafd92_fdnames{2}='Replications';
fdafd92_fdnames{3}='Metres';
fdafd92=putnames(fdafd92,fdafd92_fdnames);
lambda=0;
fdPar92=fdPar(fdabasis92,4,lambda);
nrep=size(fdaarray,2);
colsele=jet(2);
colvec=zeros([nrep,3]);
colvec(1:10,:)=repmat(colsele(1,:),[10,1]);
colvec(11:20,:)=repmat(colsele(2,:),[10,1]);
% res=fda_pca(fdafd92,fdPar92,3);
res=fda_pca(fdafd92,fdPar92,5,colvec);
harmscr=res.fdarotpcastr.harmscr;
% isequaln(res.fdapcastr,testrec.fdapcastr)
% isequaln(res.fdarotpcastr,testrec.fdarotpcastr)

% h2=figure();
%   scatter(harmscr(:,1,2),harmscr(:,2,2),[],colvec,'filled');
%   xlabel('pc1');
%   ylabel('pc2');
%%one diemension
fdameanfd=mean(fdafd92);
fdameanmat=squeeze(eval_fd(fdatime,fdameanfd));
harmfd=res.fdarotpcastr.harmfd;
harmmat=eval_fd(fdatime,harmfd);
nfine=201;
fdapointtime=linspace(0,2300,nfine)';
fdameanpoint=squeeze(eval_fd(fdapointtime,fdameanfd));
harmpointmat=squeeze(eval_fd(fdapointtime,harmfd));
fac=0.1;
harmplusmat=zeros(nfine,3,1);
harmminsmat=zeros(nfine,3,1);
for j=1:3
    harmplusmat(:,j,:)=fdameanpoint+fac.*squeeze(harmpointmat(:,j,:));
    harmminsmat(:,j,:)= fdameanpoint-fac.*squeeze(harmpointmat(:,j,:));
end
for j=1:2
  fig=figure();
  phdl=plot(fdatime,fdameanmat,'b-');
  set(phdl,'LineWidth',2)
  hold on
  plot(fdapointtime,harmplusmat(:,j,1),'b--');
  plot(fdapointtime,harmminsmat(:,j,1),'b--');
  hold off
  axis('square')
end

%%two dimension for Handwriting data
load fda;
nbasis92=92;
norder92=6;
fdabasis92=create_bspline_basis(fdarange,nbasis92,norder92);
fdafd92=smooth_basis(fdatime,fdaarray,fdabasis92);%fdaarray(:,:,1)
%  Add suitable names for the dimensions of the data.
fdafd92_fdnames{1}='Milliseconds';
fdafd92_fdnames{2}='Replications';
fdafd92_fdnames{3}='Metres';
fdafd92=putnames(fdafd92,fdafd92_fdnames);
lambda=0;
fdPar92=fdPar(fdabasis92,4,lambda);
nrep=size(fdaarray,2);
colsele=jet(2);
colvec=zeros([nrep,3]);
colvec(1:10,:)=repmat(colsele(1,:),[10,1]);
colvec(11:20,:)=repmat(colsele(2,:),[10,1]);
% res=fda_pca(fdafd92,fdPar92,3);
res=fda_pca(fdafd92,fdPar92,5,colvec);

fdameanfd=mean(fdafd92);
fdameanmat=squeeze(eval_fd(fdatime,fdameanfd));
harmfd=res.fdarotpcastr.harmfd;
harmmat=eval_fd(fdatime,harmfd);
nfine=201;
fdapointtime=linspace(0,2300,nfine)';
fdameanpoint=squeeze(eval_fd(fdapointtime,fdameanfd));
harmpointmat=squeeze(eval_fd(fdapointtime,harmfd));
fac=0.1;
harmplusmat=zeros(nfine,3,2);
harmminsmat=zeros(nfine,3,2);
for j=1:3
    harmplusmat(:,j,:)=fdameanpoint+fac.*squeeze(harmpointmat(:,j,:));
    harmminsmat(:,j,:)=fdameanpoint-fac.*squeeze(harmpointmat(:,j,:));
end

for j=1:2
    subplot(1,2,j)
    phdl=plot(fdameanmat(:,1),fdameanmat(:,2),'b-');
    set(phdl,'LineWidth',2)
    hold on
    plot(harmplusmat(:,j,1),harmplusmat(:,j,2),'b--');
    plot(harmminsmat(:,j,1),harmminsmat(:,j,2),'b--');
    hold off
    axis('square')
end

%example 2
%% test on ridge tracking data set
%% smooth the curve first and then PCA (not smoothed)
close all;
clear all;
load('/Users/mikeaalv/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/fda.test/tracing.newmeth.experiment.manual.mat')
nsample=6;
isample=1;
ntime=52;
maxridges=60;
nDer=1;%maximal derivative to consider
xvec=1:ntime;
ymat=NaN(ntime,nsample,maxridges);
for isample=1:nsample
  dataset=Sample(isample).ridges;
  nridges=size(dataset,2)-1;
  for i=1:nridges
    locstr=Sample(isample).ridges(i+1).result;
    timevec=locstr.rowind;
    intensity=locstr.intensity;
    indvec=matchPPMs(timevec',xvec);
    locintvec=NaN(1,ntime);
    locintvec(indvec)=intensity;
    locintvec(find(isnan(locintvec)))=min(intensity);
    ymat(:,isample,i)=locintvec;
  end
end
naid=find(isnan(ymat));
ymat(naid)=min(ymat(:));
%% add small noise for NA ridges in some samples
snoise_sd=min(abs(ymat(:)))*0.01;
ymat(naid)=ymat(naid)+normrnd(0,snoise_sd,1,length(naid))';
%% center and scale for each feature
for i=1:nsample
  for j=1:maxridges
    shiftv=ymat(:,i,j)-mean(ymat(:,i,j));
    ymat(:,i,j)=shiftv./std(shiftv);
  end
end
loglambda_vec= -4:0.25:4;
res=smooth_derivative(ymat,xvec,loglambda_vec,nDer)
%best log lambda 2
% rangcheck=randsample(nridges,5)%1:5;%randsample(nridges,10);%35:40%1:10%14:20%
fdres=res.spfd_selec;
% xfine=linspace(min(xvec),max(xvec),201)';
% intmatfine=eval_fd(xfine,fdres(rangcheck));
% phdl=plot(xfine,intmatfine,'-');
% set(phdl,'LineWidth',2);
% hold on;
% plot(xvec,ymat(:,rangcheck),'o')
% hold off;
% xlabel('\fontsize{19} time')
% ylabel('\fontsize{19} intensity')

% velfmatfine=eval_fd(xfine,fdres(rangcheck), 1);
% figure();
% phdl=plot(xfine,velfmatfine,'-');
% set(phdl,'LineWidth',2)
% xlabel('\fontsize{19} time')
% ylabel('\fontsize{19} intensity Velocity')
nsmooth=nDer+2;
norder=nsmooth+2;
nbasis=length(xvec)+norder-2;
bbasis=create_bspline_basis([min(xvec) max(xvec)],nbasis,norder,xvec);
% fdafd92=smooth_basis(fdatime,fdaarray,fdabasis92);
% lambdapc=0.01;
lambdapc=0.0;
fdParpc=fdPar(bbasis,nsmooth,lambdapc);
colsele=jet(2);
colvec=zeros([nsample,3]);
colvec(1:3,:)=repmat(colsele(1,:),[3,1]);
colvec(4:6,:)=repmat(colsele(2,:),[3,1]);
res=fda_pca(fdres,fdParpc,3,colvec);
% subplot(1,1,1)
% plot_pca_fd(res.fdarotpcastr)
harmscr=res.fdarotpcastr.harmscr;
h2=figure();
  scatter(harmscr(:,2),harmscr(:,3),[],colvec,'filled');
  xlabel('pc2');
  ylabel('pc3');


%% the PCA is done on each sample and different ridges is treated as different samples
close all;
clear all;
load('/Users/mikeaalv/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/fda.test/tracing.newmeth.experiment.manual.mat')
isample=1;%%just use sample 1
ntime=52;
nDer=1;%maximal derivative to consider
xvec=1:ntime;
dataset=Sample(isample).ridges;
nridges=size(dataset,2)-1;
ymat=NaN(ntime,nridges,1);

for i=1:nridges
  locstr=Sample(isample).ridges(i+1).result;
  timevec=locstr.rowind;
  intensity=locstr.intensity;
  indvec=matchPPMs(timevec',xvec);
  locintvec=NaN(1,ntime);
  locintvec(indvec)=intensity;
  locintvec(find(isnan(locintvec)))=min(intensity);
  ymat(:,i,1)=locintvec;
end
naid=find(isnan(ymat));
ymat(naid)=min(ymat(:));
%% add small noise for NA ridges in some samples
snoise_sd=min(abs(ymat(:)))*0.01;
ymat(naid)=ymat(naid)+normrnd(0,snoise_sd,1,length(naid))';
%% center and scale for each feature
for j=1:nridges
  shiftv=ymat(:,j,1)-mean(ymat(:,j,1));
  ymat(:,j,1)=shiftv./std(shiftv);
end
loglambda_vec= -4:0.25:4;
res=smooth_derivative(ymat,xvec,loglambda_vec,nDer);
%best log lambda 2
fdres=res.spfd_selec;
% rangcheck=randsample(nridges,5)%1:5;%randsample(nridges,10);%35:40%1:10%14:20%
% xfine=linspace(min(xvec),max(xvec),201)';
% intmatfine=eval_fd(xfine,fdres(rangcheck));
% phdl=plot(xfine,intmatfine,'-');
% set(phdl,'LineWidth',2);
% hold on;
% plot(xvec,ymat(:,rangcheck),'o')
% hold off;
% xlabel('\fontsize{19} time')
% ylabel('\fontsize{19} intensity')

% velfmatfine=eval_fd(xfine,fdres(rangcheck), 1);
% figure();
% phdl=plot(xfine,velfmatfine,'-');
% set(phdl,'LineWidth',2)
% xlabel('\fontsize{19} time')
% ylabel('\fontsize{19} intensity Velocity')
nsmooth=nDer+2;
norder=nsmooth+2;
nbasis=length(xvec)+norder-2;
bbasis=create_bspline_basis([min(xvec) max(xvec)],nbasis,norder,xvec);
% fdafd92=smooth_basis(fdatime,fdaarray,fdabasis92);
% lambdapc=0.01;
lambdapc=0.0;
fdParpc=fdPar(bbasis,nsmooth,lambdapc);
% colsele=jet(2);
% colvec=zeros([nsample,3]);
% colvec(1:3,:)=repmat(colsele(1,:),[3,1]);
% colvec(4:6,:)=repmat(colsele(2,:),[3,1]);
res=fda_pca(fdres,fdParpc,5);
% subplot(1,1,1)
% plot_pca_fd(res.fdapcastr)
harmscr=res.fdarotpcastr.harmscr;
h2=figure();
  scatter(harmscr(:,2),harmscr(:,3),[],'k','filled');
  xlabel('pc2');
  ylabel('pc3');
