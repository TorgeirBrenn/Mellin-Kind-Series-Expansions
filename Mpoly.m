%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Mn = Mpoly(n, a, x) returns the n-th M polynomial with
% coefficient a evaluated at the points in the vector x.
%
% In order to scale the polynomials, use b*x instead of x.
%
% This function is made for the FYS-3740 Project Paper in Applied Physics
% and Mathematics.
%
% Made by Torgeir Brenn, november 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mn = Mpoly(n, a, x)
    
    Mn = zeros(size(x));

    for k=0:n
        Mn = Mn + stirling2(n,k)*(-1)^k*factorial(k)*laguerrefast(k,a-1,x);
    end
    %laguerrefast is an optimised version of laguerreL, but is limited to
    %input k<=6. It should easy to make a version which is slightly slower
    %than the hard coded laguerrefast, but much faster than laguerreL and
    %valid for all k.
end