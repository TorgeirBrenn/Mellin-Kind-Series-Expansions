%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mkbkfitexactgivenparam : fits the Mellin kind beta prime kernel series.
%
%g = mkgkfitgivenparam(x, k, N, a1, a2, b) evaluates the fitted PDF at the 
%points in x. The output is the Mellin kind beta prime kernel expansion of 
%the beta prime PDF based on the log-cumulants k, and N is the order of the
%highest order log-cumulant used.
%
%INPUT
%x  : Points at the which the fitted PDF will be evaluated, vector.
%k  : Log-cumualnts, vector.
%N  : The highest order of the log-cumulants used in the expansion, N<=8.
%a1 : First kernel shape parameter
%a2 : Second kernel shape parameter
%b  : Kernel scale parameter
%
%OUTPUT
%f : The fitted PDF evaluated at the points in x.
%
%Last update: 2017-04-19
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = mkbkfitgivenparam(x, k, N, a1, a2, b)
    
    if N > 8
        disp 'Warning: N>8 not supported. N=8 will be used.'
        N = 8;
    end
    
    
    %Creating the kernel
    kernel = betaprimepdf(x, a1, a2, b); %The kernel
    f = kernel;
    y = b*x./(1+b*x);
    
    %Log-cumulant differences
    dk = zeros(N, 1);

    dk(1) = k(1) - psi(0,a1) + psi(0,a2) + log(b);
    f = f + kernel.*Mprimepoly(1, a1, a2, y)*Bpoly(dk(1));
    
    for n = 2:N
        %Computing the next log-cumulant difference
        dk(n) = k(n) - psi(n-1,a1) - (-1)^n*psi(n-1,a2);
  
        %Adding the term to the MKBK series
        f = f + kernel.*Mprimepoly(n, a1, a2, y)*Bpoly(dk(1:n))/factorial(n);
    end

end