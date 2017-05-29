%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kldist : returns the Kullback Leibler distance between two functions
%
%dkl = kldist(f, g, dx) returns the Kullback Leibler (KL) distance between
%f(x) and g(x), where x has grid resolution dx. f and g must be vectors of
%equal length where each element f(n) corresponds to the same value of x as
%g(n). The KL divergence is non-negative and definite, but not symmetric 
%and the triangle inequality does not hold. The KL distance is the mean of 
%the KL divergences each way and is thus also symmetric.
%
%INPUT
%f  : Input function f(x) (e.g. true function).
%g  : Input function g(x) (e.g. est. function).
%dx : Grid resolution of x. 
%
%OUTPUT
%dkl : The KL distance between f and g.
%
%Last update: 2017-01-12
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dkl = kldist(f, g, dx)
    
    %f = normalizepdf(f,dx);
    g = normalizepdf(g,dx);
    dkl = dx/2*(sum((f-g).*log(f./g)));

end