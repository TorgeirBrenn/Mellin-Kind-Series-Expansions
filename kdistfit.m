%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB code: kdistfit.m
% Last update: 31/05/12
%
% DESCRIPTION:
% [L,M,m] = kdistfit(x) produces estimates of the shape parameters L and M 
% and the location parameter m of the K distribution for the data set x.
%
% OUTPUT
% L : First shape parameter of the Fisher distribution
% M : Second shape parameter of the Fisher distribution
% m : Location parameter of the Fisher distribution
%
% INPUT
% x : Data set
%
% AUTHOR : Stian Normann Anfinsen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L,M,m] = kdistfit(x)

  x = x(:);
  y = log(x);

  % Compute log-cumulants
  k1 = mean(y);
  k2 = var(y);
  k3 = mean((y - k1).^3);

  % Compute initial parameter estimates
  init = fzero(@(x) (psi(1,x*x) - 0.5*k2),1);
  init = init * init;

  % Solve matrix log-cumulant equations
  % NB! Initialisations of L and M must be distinct
  % to ensure numerical stability
  [M,L] = mexKdistGradOpt(1.1*init,0.9*init,k2,k3);
  m = mean(x);
  %m = exp(k1 - psi(L) - psi(M) + log(L*M))
  

return