%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%molcggdfit : fits a generalized gamma PDF to data using the method of 
%log-cumulants.
%
%
%g = molcggdfit(x,d) evaluates the generalized gamma distribution at the 
%points in x, using the parameter estimates found with the method of log
%cumulants. The method is found in [Li et al., 2011]
%
%INPUT
%x : Points at the which the log-normal PDF will be evaluated, vector.
%d : Data, vector.
%
%OUTPUT
%g : The fitted generalized gamma PDF evaluated at points in x.
%
%Last update: 2017-03-09
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = molcggdfit(x,d)
    
    k = emplc(d, 3); %Empirical log-cumulants of order 1 through 3
    
    a0 = 8*k(3)^2;
    a1 = 4*(3*k(3)^2-2*k(2)^3);
    a2 = 2*(3*k(3)^2-8*k(2)^3);
    a3 = k(3)^2-8*k(2)^3;
    
    p = (3*a0*a2-a1^2)/(3*a0^2);
    q = (2*a1^3-9*a0*a1*a2+27*a0^2*a3)/(27*a0^3);
    
    %Original version, prone to occasional failures
    %kappa = -a1/(3*a0)+nthroot(-q/2+sqrt((q/2)^2+(p/3)^3),3)...
    %    +nthroot(-q/2-sqrt((q/2)^2+(p/3)^3),3);
    
    kappa = -a1/(3*a0);
    if ((q/2)^2+(p/3)^3)>0
        kappa = kappa + nthroot(-q/2+sqrt((q/2)^2+(p/3)^3),3)...
            +nthroot(-q/2-sqrt((q/2)^2+(p/3)^3),3);
    else
        kappa = kappa + 2*nthroot(-q/2,3);
        disp 'Warning: GGD was about fail.'
    end
    
    
    nu = sign(-k(3))*sqrt(psi(1,kappa)/k(2));
    
    sigma = exp(k(1)-(psi(0,kappa)-log(kappa))/nu);
    
    
    
    
    %g = abs(nu)*kappa^kappa/(sigma*gamma(kappa))*(x/sigma).^(kappa*nu-1)...
    %    .*exp(-kappa*(x/sigma).^nu);
    
    %Trick
    C = kappa*log(kappa)-gammaln(kappa);
    lg = log(abs(nu)/sigma) + C + (kappa*nu-1)*log(x/sigma)-kappa*(x/sigma).^nu;
    g = exp(lg);
    
    
    
   
    %COMMENT 2017-02-22
    %This function oftengave Inf or NaN values due to the estimated kappa
    %being too high, which gave kappa^kappa and/or gamma(kappa) as Inf, 
    %resulting in Inf or NaN output. The "trick" is to do all the work in
    %the logarithmic domain and then transform back. The problem arose
    %around kappa=170.
    
end