%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%partialmkgkfit : fits the Mellin kind Gamma kernel PDF to data.
%
%g = partialmkgkfit(x, d, N) evaluates the fitted PDF at the points in x.
%The output is the Mellin kind Gamma kernel expansion of the Gamma PDF 
%based on log-cumulants. The log-cumulants are based on the data d, and N
%is the order of the highest order log-cumulant used. The partial MKGK
%scales each term with (n+1)! instead of n!, resulting in less propagation.
%
%INPUT
%x : Points at the which the fitted PDF will be evaluated, vector.
%d : Data, vector.
%N : The highest order of the log-cumulants used in the expansion, N<=6.
%
%OUTPUT
%f : The fitted PDF evaluated at the points in x.
%
%Last update: 2017-01-12
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = partialmkgkfit(x, d, N)
    
    if N > 6
        disp 'Warning: N>6 not supported. N=6 will be used.'
        N = 6;
    end
    
    if N < 3
        disp(['Warning: For N<3, this method is functionally ' ...
            'equivalent to molcgamfit().'])
        f = molcgamfit(x,d);
        return
    end

    k = emplc(d, N); %Empirical log-cumulants;
    
    %The kernel is found using the empirical log-cumulants
    L = fsolve(@(L) psi(1,L)-k(2),4, ...
        optimoptions('fsolve', 'Display', 'off'));
    mu = L*exp(k(1)-psi(0,L));
    
    b = L/mu;
    a = L-1;
    
    g = b^(a+1)*x.^a.*exp(-b*x)/gamma(a+1); %The kernel
    
    dk = zeros(N, 1); %Log-cumulant differences
    
    %dk(1) = dk(2) = 0
    
    %For efficiency, the first term is outside the for-loop.
    dk(3) = k(3) - psi(2,a+1);
    f = g.*(1+Mpoly(3,a,b*x)/24*Bpoly(dk(1:3))); %For efficiency.
    
    for n = 4:N
        %Computing the next log-cumulant difference
        dk(n) = k(n) - psi(n-1,a+1);
        
        %Adding the term to the MKGK series
        f = f + g.*Mpoly(n,a,b*x)*Bpoly(dk(1:n))/factorial(n+1);
    end

end