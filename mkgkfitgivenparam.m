%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mkgkfitgivenparam : fits the Mellin kind Gamma kernel series.
%
%g = mkgkfitgivenparam(x, k, N, L, mu) evaluates the fitted PDF at the 
%points in x. The output is the Mellin kind Gamma kernel expansion of 
%the Gamma PDF  based on the log-cumulants k, and N is the highest order 
%log-cumulant used. The kernel parameters are given as inputs.
%
%INPUT
%x  : Points at the which the fitted PDF will be evaluated, vector.
%k  : Log-cumulants, vector.
%N  : The highest order of the log-cumulants used in the expansion, N<=8.
%L  : Kernel shape parameter
%mu : Kernel location parameter
%
%OUTPUT
%f : The fitted PDF evaluated at the points in x.
%
%Last update: 2017-04-18
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = mkgkfitgivenparam(x, k, N, L, mu)
    
    if N > 8
        disp 'Warning: N>8 not supported. N=8 will be used.'
        N = 8;
    end
    
    b = L/mu;
    a = L;
    
    g = b^a*x.^(a-1).*exp(-b*x)/gamma(a); %The kernel
    f = g;
    dk = zeros(N, 1); %Log-cumulant differences     
    
    if N>=1
        dk(1) = k(1) - psi(0,a) + log(b);
        f = f + g.*Mpoly(1, a, b*x)*Bpoly(dk(1));
    end
    
    for n = 2:N
        %Computing the next log-cumulant difference
        dk(n) = k(n) - psi(n-1,a);
        
        %Adding the term to the MKGK series
        f = f + g.*Mpoly(n,a,b*x)*Bpoly(dk(1:n))/factorial(n);
    end

end