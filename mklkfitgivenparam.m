%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mklkfitgivenparam : fits a series expansion of the log-normal PDF.
%
%f = mklkfitgivenparam(x, k, N, mu, sigma) evaluates the fitted PDF at the 
%points in x. The output is the Mellin kind log-normal kernel 
%(Gram-Charlier) series expansion of the log-normal PDF based on 
%log-cumulants.
%
%INPUT
%x      : Points at the which the fitted PDF will be evaluated, vector.
%k      : Log-cumulants, vector.
%N      : The highest order of the log-cumulants used in the expansion, N<=18.
%mu     : Kernel log-mean parameter
%logvar : Kernel log-variance parameter
%
%OUTPUT
%f : The fitted PDF evaluated at the points in x.
%
%Last update: 2017-04-19
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = mklkfitgivenparam(x, k, N, mu, logvar)
    
    if N > 18
        disp 'Warning: N>18 not supported. N=18 will be used.'
        N = 18;
    end
    
    sigma = sqrt(logvar);
    
    y = (log(x)-mu)/sigma; %Standardizing the argument
    fln = pdf(makedist('Lognormal', 'mu', mu, 'sigma', sigma), x);
    f = fln; %Just the kernel

    
    for n = 1:N
        f = f + fln.*hermitefast(n,y)/(factorial(n)*sigma^n)...
            *Bpoly([k(1)-mu;k(2)-sigma^2;k(3:n)]);
    end
    
    %hermitefast is an optimised version of hermiteHprob, but is limited to
    %input n<=12. 
end