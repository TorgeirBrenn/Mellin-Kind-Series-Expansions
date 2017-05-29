%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%emplc : epirical log-cumulants
%
%
%lc = emplc(d, N) returns the N first log-cumulants of the data d. This is
%done by find the log-moments (trivial) and using the combinatoric
%relationship between log-moments and log-cumulants, stated in e.g.
%Anfinsen and Eltoft, 2011.
%
%INPUT
%d : Data, vector.
%N : Number of log-cumulants in output, positive integer.
%
%OUTPUT
%lc : Log-cumulants of d, Nx1 vector.
%
%Last update: 2017-01-04
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lc = emplc(d, N)
    lm = zeros(N,1); %Log-moments
    lc = zeros(N,1); %Log-cumulants
    
    ld = log(d); %Compute this once, for efficiency.
    
    for n = 1:N
        lm(n) = mean(ld.^n);
        lc(n) = lm(n);
        for i=1:(n-1)
            lc(n) = lc(n) - nchoosek(n-1,i-1)*lc(i)*lm(n-i);
        end
    end
end