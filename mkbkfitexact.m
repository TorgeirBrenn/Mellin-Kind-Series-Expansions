%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mkbkfitexact : fits the Mellin kind beta prime kernel series.
%
%g = mkgkfitexact(x, k, N) evaluates the fitted PDF at the points in x.
%The output is the Mellin kind beta prime kernel expansion of the beta 
%prime PDF based on the log-cumulants k, and N is the order of the highest 
%order log-cumulant used.
%
%INPUT
%x : Points at the which the fitted PDF will be evaluated, vector.
%k : Log-cumualnts, vector.
%N : The highest order of the log-cumulants used in the expansion, N<=8.
%
%OUTPUT
%f : The fitted PDF evaluated at the points in x.
%
%Last update: 2017-04-08
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = mkbkfitexact(x, k, N)
    
    if N > 8
        disp 'Warning: N>8 not supported. N=8 will be used.'
        N = 8;
    end
    
    %The kernel parameters are found using the empirical log-cumulants
    %syms s1 s2
    %[a1, a2] = solve([psi(1,s1)+psi(1,s2)==k(2), ...
    %    psi(2,s1)-psi(2,s2)==k(3)], [s1, s2]);
    
    [a1,a2] = mexBPdistGradOpt(3,5,k(2),k(3));
    b = exp(psi(0,a1)-psi(0,a2)-k(1));
    
    
    %Creating the kernel
    kernel = betaprimepdf(x, a1, a2, b); %The kernel
    f = kernel;
    y = b*x./(1+b*x);
    
    %Log-cumulant differences
    dk = zeros(N, 1); %dk(1) = dk(2) = dk(3) = 0
          
%     %Mprime polynomials
%     %This is not tested so there might errors, Mprimepoly() is faster
%     anyway.
%     syms u
%     Mps = (a1+a2)*u/(1+u)-a1; %Symbolic Mprime, initated
%     Mps1 = Mps;
%     Mprime = zeros(numel(x), N);
%     Mprime(:, 1) = subs(Mps, y);
%     for n = 2:N
%         Mps = Mps*Mps1 - u*diff(Mps(u),u);
%         Mprime(:, n) = subs(collect(simplify(Mps(u/(1-u)))), y);
%     end
    
    for n = 4:N
        %Computing the next log-cumulant difference
        dk(n) = k(n) - psi(n-1,a1) - (-1)^n*psi(n-1,a2);
  
        %Adding the term to the MKBK series
        %f = f + kernel.*Mprime(:, n)*Bpoly(dk(1:n))/factorial(n);
        f = f + kernel.*Mprimepoly(n, a1, a2, y)*Bpoly(dk(1:n))/factorial(n);
    end

end