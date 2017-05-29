%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%targetpdf : produces the target PDF for analysis in my Master's thesis.
%
%[x, f, kappa] = mkgkfit(dist, param, N) creates a PDF f which follows the
%distribution 'dist' with the parameter values specified in 'param'. f is
%evaluated at the values of x, which are (almost) 10,000 equally spaced
%values between the points where f first reaches 0.001 of its maximum, and
%where it first goes below 0.001 of its maximum.
%
%The optional input N indicates that the first N theoretical log-cumulants
%are returned in kappa.
%
%If dist = 'kernel', then the second input ('param') is taken as a data
%vector, from which a kernel density estimate should be produced.
%
%INPUT
%dist  : Distribution type, supported: 'k', 'gamma', 'igam', 'g0', 'ggd'
%param : Parameter struct
%N     : Optional input which indicates that we want the first N
%        theoretical log-cumulants returned in kappa.
%
%OUTPUT
%x     : Points at the which the fitted PDF is evaluated, vector.
%f     : The target PDF evaluated at the points in x.
%kappa : The first N theoretical log-cumulants. If no N was specified in
%        the input, only the first log-cumulant is returned.
%
%Last update: 2017-05-22
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x, f, kappa] = targetpdf(dist, param, N)

numx  = 1e4;  %Number of x values
relth = 1e-3; %Relative threshold which limits the domain

dx = 0.1;
x = dx:dx:1000; x=x';

if nargin == 2
    N = 1;
end
kappa = zeros(N,1); %Exact log-cumulants

switch dist
    case "k"
        if (param.M<0)
            disp 'Warning: M<0 in K distribution, using -M.'
            param.M = -param.M;
        end
        f = kpdf(x, param.m, param.M, param.L);
        threshold = max(f)*relth;
        if (f(numel(x))>threshold)
            disp 'Warning: Increase the domain (larger max. x value).'
            %Warning that a significant part of the PDF might be omitted.
        end
        x = x(f>threshold);
        x = linspace(x(1)-dx, x(end)+dx, numx); x=x';
        dx = x(2)-x(1);
        f = kpdf(x, param.m, param.M, param.L);
        x = x(f>threshold);
        x = linspace(x(1)-dx, x(end)+dx, numx); x=x';
        f = kpdf(x, param.m, param.M, param.L);
        
        %Log-cumulants
        kappa(1) = log(param.m/(param.L*param.M)) ...
            + psi(param.L) + psi(param.M);
        for n = 2:N
            kappa(n) = psi(n-1,param.L) + psi(n-1,param.M);
        end
    case "gamma"
        f = gampdf(x, param.L, param.m/param.L);
        threshold = max(f)*relth;
        if (f(numel(x))>threshold)
            disp 'Warning: Increase the domain (larger max. x value).'
            %Warning that a significant part of the PDF might be omitted.
        end
        x = x(f>threshold);
        x = linspace(x(1)-dx, x(end)+dx, numx); x=x';
        dx = x(2)-x(1);
        f = gampdf(x, param.L, param.m/param.L);
        x = x(f>threshold);
        x = linspace(x(1)-dx, x(end)+dx, numx); x=x';
        f = gampdf(x, param.L, param.m/param.L);
        
        %Log-cumulants
        kappa(1) = log(param.m/param.L) + psi(param.L);
        for n = 2:N
            kappa(n) = psi(n-1,param.L);
        end
    case "g0"
        if (param.M>0)
            disp 'Warning: M>0 in G0 distribution, using -M.'
            param.M = -param.M;
        end
        %G0 distribution, M/mu->0, M<0
        f = g0pdf(x, param.g, param.M, param.L);
        threshold = max(f)*relth;
        if (f(numel(x))>threshold)
            disp 'Warning: Increase the domain (larger max. x value).'
            %Warning that a significant part of the PDF might be omitted.
        end
        x = x(f>threshold);
        x = linspace(x(1)-dx, x(end)+dx, numx); x=x';
        dx = x(2)-x(1);
        f = g0pdf(x, param.g, param.M, param.L);
        x = x(f>threshold);
        x = linspace(x(1)-dx, x(end)+dx, numx); x=x';
        f = g0pdf(x, param.g, param.M, param.L);
        
        %Log-cumulants
        kappa(1) = log(param.g/param.L) + psi(param.L) - psi(-param.M);
        for n = 2:N
            kappa(n) = psi(n-1,param.L) + (-1)^n*psi(n-1,-param.M);
        end
    case "ggd"
        f = ggdpdf(x, param.nu, param.L, param.m/param.L);
        threshold = max(f)*relth;
        if (f(numel(x))>threshold)
            disp 'Warning: Increase the domain (larger max. x value).'
            %Warning that a significant part of the PDF might be omitted.
        end
        x = x(f>threshold);
        x = linspace(x(1)-dx, x(end)+dx, numx); x=x';
        dx = x(2)-x(1);
        f = ggdpdf(x, param.nu, param.L, param.m/param.L);
        x = x(f>threshold);
        x = linspace(x(1)-dx, x(end)+dx, numx); x=x';
        f = ggdpdf(x, param.nu, param.L, param.m/param.L);
        
        %Log-cumulants
        kappa(1) = log(param.m/param.L) + psi(param.L)/param.nu;
        for n = 2:N
            kappa(n) = psi(n-1,param.L)/param.nu^n;
        end
    case "igam"
        f = igampdf(x, param.L, param.m/param.L);
        threshold = max(f)*relth;
        if (f(numel(x))>threshold)
            disp 'Warning: Increase the domain (larger max. x value).'
            %Warning that a significant part of the PDF might be omitted.
        end
        x = x(f>threshold);
        x = linspace(x(1)-dx, x(end)+dx, numx); x=x';
        dx = x(2)-x(1);
        f = igampdf(x, param.L, param.m/param.L);
        x = x(f>threshold);
        x = linspace(x(1)-dx, x(end)+dx, numx); x=x';
        f = igampdf(x, param.L, param.m/param.L);
        
        %Log-cumulants
        kappa(1) = log(param.m/param.L) - psi(param.L);
        for n = 2:N
            kappa(n) = (-1)^n*psi(n-1,param.L);
        end
    case "kernel" %Special case
        dat = param;
        dx = max(dat)/1e4;
        x  = dx:dx:max(dat);
        fPD = fitdist(dat, 'Kernel', 'kernel', 'epanechnikov', ...
            'support', 'positive');
        f = pdf(fPD, x);
        
        %Finding the part of the PDF which is above the threshold
        [~, fmode] = max(f); %Returns the index of max(f)
        threshold = f(fmode)*relth;
        i = fmode; j = fmode;
        fstart = 0; fstop = 0;
        while fstart == 0
            if f(i)<threshold || i == 1
                fstart = i;
            end
            i = i - 1;
        end
        while fstop == 0
            if f(j) < threshold || i == numx
                fstop = j;
            end
            j = j + 1;
        end
        x = linspace(x(fstart)-dx, x(fstop)+dx, numx); x=x';
        f = pdf(fPD, x);
        
        %Just need a final check that f is never 0
        x = x(f>0); 
        f = f(f>0);
    otherwise
        error 'Distribution not supported.'
end

%Limit domain
threshold = max(f)*relth;
if (f(numel(x))>threshold)
    disp 'Warning: Consider increasing the domain (larger max. x value).'
    %Warning that a significant part of the PDF might be omitted.
end
x  = x(f>threshold);
f  = f(f>threshold);

end