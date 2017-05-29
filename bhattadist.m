%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bhattadist : returns the Bhattacharyya distance between two functions
%
%db = bhattadist(f, g, dx) returns the Bhattacharyya distance between f(x) 
%and g(x), where x has grid resolution dx. f and g must be vectors of equal
%length where each element f(n) corresponds to the same value of x as g(n).
%The Bhattacharyya distance is non-negative, symmetric and definite, but
%the triangle inequality does not hold. For the Hellinger distance
%sqrt(1-e^-db), the triangle inequality does hold.
%
%INPUT
%f  : Input function f(x) (e.g. true function).
%g  : Input function g(x) (e.g. est. function).
%dx : Grid resolution of x. 
%
%OUTPUT
%db : The Bhattacharyya distance between f and g.
%
%Last update: 2017-01-12
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function db = bhattadist(f, g, dx)
    
    BC = dx * sum(sqrt(f.*normalizepdf(g,dx))); %Bhattacharyya coefficient
    
    %Hellinger distance = sqrt(1-BC);
    
    db = -log(BC);

end