function [returnstr]=fda_pca(linefda,pcfda,npc,colvec)
% a wrap up for PCA in Functional Data Analysis (FDA)
% Functional pca will be done based on both the fda of original cuve and of the principal elements (PCs), including lambda.
% varmx will be used to rotate PCs to get simplified representation. Rotated PCs will be returned.
% Both score and eigen value (explained variance) will be plot.
% The user can use the following to plot the PC curve
%%              subplot(1,1,1)
%%              plot_pca_fd(fdarotpcastr)
%% Be careful when Ymat is multi-dimension, in which case the score plot is just the preojecting of first feature(dimension) on PC1 and PC2.
%
% Argument:         linefda: struct. The fda structure for original smoothed line. Must be provided.
%         pcfda: struct. The fda structure for PCs. Must be provided. lambda is set here and the user might try multiple choices.
%         npc: numeric. the number of pc to get/retrieve. default 2.
%         colvec: vector. the color vector for score plot. optional
%
% Return: return the rotated functional pca result.
% CITATION:
%         Ramsay, James & Hooker, Giles & Graves, Spencer. (2009). Functional data analysis with R and MATLAB. 10.1007/978-0-387-98185-7.
% More information:
%     http://www.psych.mcgill.ca/misc/fda/
%     https://cran.r-project.org/web/packages/fda/index.html
%     http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/
%
% CAUTION:
% the function seems to have bugs for multi-dim FDA PCA.
% YUE WU 01012020


if ~exist('linefda','var')
  error('linefda is needed');
end
if ~exist('pcfda','var')
  error('pcfda is needed');
end
if ~exist('npc','var')
  npc=2;
end
if ~exist('colvec','var')
  colvec='k';
end
fdapcastr=pca_fd(linefda,npc,pcfda);
%  plot unrotated harmonics
% subplot(1,1,1)
% plot_pca_fd(fdapcastr)
%  Varimax rotation
fdarotpcastr=varmx_pca(fdapcastr);
%  plot rotated harmonics
% subplot(1,1,1)
% plot_pca_fd(fdarotpcastr)
%  plot log eigenvalues
fdaharmeigval=fdapcastr.varprop;
nshow=npc;
h1=figure();
  plot(1:nshow,fdaharmeigval(1:nshow),'-o');
  xlabel('Eigenvalue Number');
  ylabel('Eigenvalue');
%  plot factor scores
harmscr=fdarotpcastr.harmscr;
h2=figure();
  scatter(harmscr(:,1),harmscr(:,2),[],colvec,'filled');
  xlabel('pc1');
  ylabel('pc2');
returnstr=struct();
returnstr.fdapcastr=fdapcastr;
returnstr.fdarotpcastr=fdarotpcastr;
