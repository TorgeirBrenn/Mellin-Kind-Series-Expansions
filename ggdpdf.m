%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ggdpdf : generalized gamma distribution probability density function
%
%f = ggdpdf(x, nu, kappa, sigma) returns the GGD-distribution PDF at the 
%values in x. See Stacy et al. 1962 for details.
%
%INPUT
%x     : Values at which the PDF should be evaluated, Nx1 vector.
%nu    : Power parameter of the GGD distribution, non-zero double.
%kappa : Shape parameter of the GGD distribution, positive double.
%sigma : Scale parameter of the GGD distribution, positive double.
%
%Last update: 2017-03-21
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Alternative notation: kappa->a, sigma->1/b, nu->d
function f=ggdpdf(x, nu, kappa, sigma)
    
    
    lf = log(abs(nu)/sigma)-gammaln(kappa) + (kappa*nu-1)*log(x/sigma)-(x/sigma).^nu;
    f = exp(lf);
    
    %f = abs(nu)/(sigma*gamma(kappa))*(x/sigma).^(kappa*nu-1).*exp(-(x/sigma).^nu);
end