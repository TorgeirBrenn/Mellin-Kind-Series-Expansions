%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB code: kdistfitexact.m
% Last update: 20/03/17
%
% DESCRIPTION:
% [L,M,m] = kdistfitexact(x) produces estimates of the shape parameters L and M 
% and the location parameter m of the K distribution for the log-cumulants k.
%
% OUTPUT
% L : First shape parameter of the Fisher distribution
% M : Second shape parameter of the Fisher distribution
% m : Location parameter of the Fisher distribution
%
% INPUT
% x : Data set
%
% AUTHOR : Stian Normann Anfinsen, small changes by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L,M,m] = kdistfitexact(k)

  % Compute initial parameter estimates
  init = fzero(@(x) (psi(1,x*x) - 0.5*k(2)),1);
  init = init * init;

  % Solve matrix log-cumulant equations
  % NB! Initialisations of L and M must be distinct
  % to ensure numerical stability
  [M,L] = mexKdistGradOpt(1.1*init,0.9*init,k(2),k(3));
  m = exp(k(1) - psi(L) - psi(M) + log(L*M));

return