%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mklkfitexact : fits a series expansion of the log-normal PDF.
%
%f = mklkfitexact(x, k, N) evaluates the fitted PDF at the points in x.
%The output is the novel Mellin kind log-normal kernel (Gram-Charlier)
%series expansion of the log-normal PDF based on log-cumulants.
%
%INPUT
%x : Points at the which the fitted PDF will be evaluated, vector.
%k : Log-cumulants, vector.
%N : The highest order of the log-cumulants used in the expansion, N<=18.
%
%OUTPUT
%f : The fitted PDF evaluated at the points in x.
%
%Last update: 2017-03-20
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = mklkfitexact(x, k, N)
    
    if N > 18
        disp 'Warning: N>18 not supported. N=18 will be used.'
        N = 18;
    end
    
    sigma = sqrt(k(2));
    
    y = (log(x)-k(1))/sigma; %Standardizing the argument
    fln = pdf(makedist('Lognormal', 'mu', k(1), 'sigma', sigma), x);
    f = fln; %Just the kernel
    
    if N < 3
        %disp(['Warning: For N<3, this method is simply the log-normal ' ...
        %    'distribution.'])
        return
    end

    
    for n = 3:N
        f = f + fln.*hermitefast(n,y)/(factorial(n)*sigma^n)...
            *Bpoly([0;0;k(3:n)]);
    end
    
    %hermitefast is an optimised version of hermiteHprob, but is limited to
    %input n<=12. 
end