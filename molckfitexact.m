%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%molckfitexact : fits a K PDF to data using the method of log-cumulants.
%
%k = molckfitexact(x,d) evaluates the K distribution at the points
%in x, using the parameter estimates found with the method of log
%cumulants. This function is a shell which relies on functions written by
%Stian Normann Anfinsen, which he has optimised by writing some of the code
%in C.
%
%INPUT
%x : Points at the which the K PDF will be evaluated, vector.
%k : Log-cumulants, vector.
%
%OUTPUT
%f : The fitted K PDF evaluated at points in x.
%
%Last update: 2017-03-20
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = molckfitexact(x,k)

    [L,M,m] = kdistfitexact(k);
    f = kdistpdf(x,L,M,m);
    
end