%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mkgkfitexact : fits the Mellin kind Gamma kernel series.
%
%g = mkgkfitexact(x, k, N) evaluates the fitted PDF at the points in x.
%The output is the Mellin kind Gamma kernel expansion of the Gamma PDF 
%based on the log-cumulants k, and N is the order of the highest order 
%log-cumulant used.
%
%INPUT
%x : Points at the which the fitted PDF will be evaluated, vector.
%k : Log-cumualnts, vector.
%N : The highest order of the log-cumulants used in the expansion, N<=8.
%
%OUTPUT
%f : The fitted PDF evaluated at the points in x.
%
%Last update: 2017-03-20
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = mkgkfitexact(x, k, N)
    
    if N > 8
        disp 'Warning: N>8 not supported. N=8 will be used.'
        N = 8;
    end

    
    %The kernel is found using the log-cumulants
    L = fsolve(@(L) psi(1,L)-k(2),4, ...
        optimoptions('fsolve', 'Display', 'off'));
    mu = L*exp(k(1)-psi(0,L));
    
    b = L/mu;
    a = L;
    
    g = b^a*x.^(a-1).*exp(-b*x)/gamma(a); %The kernel
    f = g;
    dk = zeros(N, 1); %Log-cumulant differences
    %dk(1) = dk(2) = 0
    
    if N < 3
        %disp(['Warning: For N<3, this method is simply the gamma' ...
        %    'distribution.'])
        return
    end        
    
    for n = 3:N
        %Computing the next log-cumulant difference
        dk(n) = k(n) - psi(n-1,a);
        
        %Adding the term to the MKGK series
        f = f + g.*Mpoly(n,a,b*x)*Bpoly(dk(1:n))/factorial(n);
    end

end