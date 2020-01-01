%example 1
%% test on the data set Berkeley Growth Data
%% make sure the fda matlab toolbox is downloaded and added to the path

%%TEST CODE FROM FDA matlab examples (http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/)
% close all;
% clear all;
% load growth;
% rng      = [1,18];
% knots    = age;
% norder   = 6;
% nbasis   = length(knots) + norder - 2;
% hgtbasis = create_bspline_basis(rng, nbasis, norder, knots);
% Lfdobj   = int2Lfd(4);
% lambda   = 10^(-3.75);
% hgtfdPar = fdPar(hgtbasis, Lfdobj, lambda);
% %  smooth the data
% hgtffd = smooth_basis(age, hgtfmat, hgtfdPar);
% agefine = linspace(1,18,101)';
% hgtfmatfine = eval_fd(agefine, hgtffd(1:10));
% phdl = plot(agefine, hgtfmatfine, '-');
% set(phdl, 'LineWidth', 2)
% hold on
% plot(age, hgtfmat(:,1:10), 'o')
% hold off
% xlabel('\fontsize{19} Age')
% ylabel('\fontsize{19} Height (cm)')
% axis([1,18,60,200])
% testrec=struct();
% testrec.hgtffd=hgtffd;
% testrec.hgtfmatfine=hgtfmatfine;
% clearvars -except testrec;

%%New code using functions
load growth;
loglambda_vec= -6:0.25:0;
res=smooth_derivative(hgtfmat,age,loglambda_vec,2);
fdres=res.spfd_selec;
agefine=linspace(min(age),max(age),101)';
hgtfmatfine=eval_fd(agefine,fdres(1:10));
phdl=plot(agefine,hgtfmatfine,'-');
set(phdl,'LineWidth',2);
hold on;
plot(age,hgtfmat(:,1:10),'o')
hold off;
xlabel('\fontsize{19} Age')
ylabel('\fontsize{19} Height (cm)')
axis([min(age),max(age),60,200])
% isequaln(hgtfmatfine,testrec.hgtfmatfine)
% isequaln(fdres,testrec.hgtffd)

%  Velocity
velfmatfine=eval_fd(agefine,fdres(1:10), 1);
phdl=plot(agefine,velfmatfine,'-');
set(phdl,'LineWidth',2)
xlabel('\fontsize{19} Age')
ylabel('\fontsize{19} Height Velocity (cm/yr)')
axis([min(age),max(age),0,20])

%  Acceleration
accfmatfine=eval_fd(agefine,fdres(1:10),2);
phdl=plot(agefine,accfmatfine,'-');
set(phdl,'LineWidth',2)
xlabel('\fontsize{19} Age')
ylabel('\fontsize{19} Height Acceleration (cm/yr/yr)')
axis([min(age),max(age),-20,10])

%example 2
%% test on ridge tracking data set
load('/Users/mikeaalv/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/fda.test/tracing.newmeth.experiment.manual.mat')
isample=1;
ntime=52;
nridges=60;
xvec=1:ntime;
ymat=NaN(ntime,nridges);
for i=1:nridges
  locstr=Sample(isample).ridges(i+1).result;
  timevec=locstr.rowind;
  intensity=locstr.intensity;
  indvec=matchPPMs(timevec',xvec);
  locintvec=NaN(1,ntime);
  locintvec(indvec)=intensity;
  locintvec(find(isnan(locintvec)))=min(intensity);
  ymat(:,i)=locintvec;
end
loglambda_vec= -4:0.25:4;
res=smooth_derivative(ymat,xvec,loglambda_vec,1)
%best log lambda 2
rangcheck=randsample(nridges,5)%1:5;%randsample(nridges,10);%35:40%1:10%14:20%
fdres=res.spfd_selec;
xfine=linspace(min(xvec),max(xvec),201)';
intmatfine=eval_fd(xfine,fdres(rangcheck));
phdl=plot(xfine,intmatfine,'-');
set(phdl,'LineWidth',2);
hold on;
plot(xvec,ymat(:,rangcheck),'o')
hold off;
xlabel('\fontsize{19} time')
ylabel('\fontsize{19} intensity')

velfmatfine=eval_fd(xfine,fdres(rangcheck), 1);
figure();
phdl=plot(xfine,velfmatfine,'-');
set(phdl,'LineWidth',2)
xlabel('\fontsize{19} time')
ylabel('\fontsize{19} intensity Velocity')
