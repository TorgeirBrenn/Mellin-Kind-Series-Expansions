%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mkgkfit : fits the Mellin kind Gamma kernel series to data.
%
%g = mkgkfit(x, d, N) evaluates the fitted PDF at the points in x.
%The output is the Mellin kind Gamma kernel expansion of the Gamma PDF 
%based on log-cumulants. The log-cumulants are based on the data d, and N
%is the order of the highest order log-cumulant used.
%
%INPUT
%x : Points at the which the fitted PDF will be evaluated, vector.
%d : Data, vector.
%N : The highest order of the log-cumulants used in the expansion, N<=8.
%L : Optional shape parameter input, if it is considered known.
%
%OUTPUT
%f : The fitted PDF evaluated at the points in x.
%
%Last update: 2017-01-06
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = mkgkfit(x, d, N)
    
    if N > 8
        disp 'Warning: N>8 not supported. N=8 will be used.'
        N = 8;
    end
    
    if N < 3
        %disp(['Warning: For N<3, this method is functionally ' ...
        %    'equivalent to molcgamfit().'])
        f = molcgamfit(x,d);
        return
    end

    k = emplc(d, N); %Empirical log-cumulants;
    
    L = fsolve(@(L) psi(1,L)-k(2), 4, ...
        optimoptions('fsolve', 'Display', 'off'));
    mu = L*exp(k(1)-psi(0,L));
    %  mu = mean(d);
    
    b = L/mu;
    a = L;
    
    g = b^a*x.^(a-1).*exp(-b*x)/gamma(a); %The kernel
    
    dk = zeros(N, 1); %Log-cumulant differences
    
    %Set f equal to the kernel
    f = g;
    
    
    for n = 3:N
        %Computing the next log-cumulant difference
        dk(n) = k(n) - psi(n-1,a);
        
        %Adding the term to the MKGK series
        f = f + g.*Mpoly(n,a,b*x)*Bpoly(dk(1:n))/factorial(n);
    end

end