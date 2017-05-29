%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%molckfit : fits a K PDF to data using the method of log-cumulants.
%
%
%k = molckfit(x,d) evaluates the K distribution at the points
%in x, using the parameter estimates found with the method of log
%cumulants. This function is a shell which relies on functions written by
%Stian Normann Anfinsen, which he has optimised by writing some of the code
%in C.
%
%INPUT
%x : Points at the which the K PDF will be evaluated, vector.
%d : Data, vector.
%
%OUTPUT
%k : The fitted K PDF evaluated at points in x.
%
%Last update: 2017-01-05
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function k = molckfit(x,d)
    
%     Alternative approach, very slow!
%     kappa = emplc(d, 3);
%     syms s1 s2
%     [M,L] = solve([psi(1,s1)+psi(1,s2)==kappa(2), ...
%         psi(2,s1)+psi(2,s2)==kappa(3)], [s1, s2]);
%     L = double(L);
%     M = double(M);
%     m = exp(kappa(1) - psi(L) - psi(M) + log(L*M));
    
    [L,M,m] = kdistfit(d);
    
    k = kdistpdf(x,L,M,m);
    
end