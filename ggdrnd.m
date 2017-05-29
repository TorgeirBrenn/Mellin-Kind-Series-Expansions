%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ggdrnd : Generalized gamma distribution random number generator
%
%r = ggdrnd(N, nu, kappa, kappa) returns N random GGD distributed samples. 
%See Li et al. 2011 for details about the distribution.
%
%INPUT
%N : Number of samples in output, positive integer.
%nu    : Power parameter of the GGD distribution, non-zero double.
%k     : Shape parameter of the GGD distribution, positive double.
%sigma : Scale parameter of the GGD distribution, positive double.
%
%OUTPUT
%r  : Random samples, Nx1 vector.
%
%Last update: 2017-03-21
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Alternative notation: kappa->a, sigma->1/b, nu->d
function r=ggdrnd(N, nu, k, sigma)
        
       r = gamrnd(k, 1, [N 1]).^(1/nu)*sigma;       

end