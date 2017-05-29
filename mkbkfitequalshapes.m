%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mkbkfitequalshapes : fits the Mellin kind beta prime kernel series with
%equal shapes.
%
%g = mkgkfitequalshapes(x, d, N) evaluates the fitted PDF at the points in
%x. The output is the Mellin kind beta prime kernel expansion of the beta
%prime PDF (with a_1=a_2=L) based on the data d, and N is the order of the
%highest order log-cumulant used.
%
%Note that this can easily be optimized by adapting Mprimepoly() to the
%fact that a_1=a_2=L.
%
%INPUT
%x          : Points at the which the fitted PDF will be evaluated, vector.
%d          : Data.
%N          : The highest order of the log-cumulants used in the expansion.
%UnitScale  : Logical, 0: MoLC value of b used, 1: b = 1. Default = 0
%OddLCDZero : Logical, 0: Corrects for all log-cumulant differences
%                      1: Assumes \Delta\kappa_n = 0 for n odd.
%                      Default: 0
%
%OUTPUT
%f : The fitted PDF evaluated at the points in x.
%
%Last update: 2017-05-18
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = mkbkfitequalshapes(x, d, N, UnitScale, OddLCDZero)

if N > 8
    disp 'Warning: N>8 not supported. N=8 will be used.'
    N = 8;
end

if nargin < 5
    OddLCDZero = 0;
end

if nargin < 4
    UnitScale = 0;
end


k = emplc(d, max(N,3)); %Empirical log-cumulants;

%The kernel parameters are found using the empirical log-cumulants
L = fsolve(@(L) psi(1,L)-k(2)/2, 4, ...
    optimoptions('fsolve', 'Display', 'off'));
if UnitScale
    b = 1;
else
    b = exp(-k(1));
end


%Creating the kernel
kernel = betaprimepdf(x, L, L, b); %The kernel
f = kernel;
y = b*x./(1+b*x);

%Log-cumulant differences
dk = zeros(N, 1); %dk(1) = dk(2) = 0

switch OddLCDZero
    case 0
        for n = 3:N
            %Computing the next log-cumulant difference
            dk(n) = k(n) - psi(n-1,L) - (-1)^n*psi(n-1,L);
            
            f = f + kernel.*Mprimepoly(n, L, L, y)*Bpoly(dk(1:n))/factorial(n);
        end
    case 1
        for n = 4:2:N
            %Computing the next log-cumulant difference
            dk(n) = k(n) - psi(n-1,L) - (-1)^n*psi(n-1,L);
            f = f + kernel.*Mprimepoly(n, L, L, y)*Bpoly(dk(1:n))/factorial(n);
        end
end

end