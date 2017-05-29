%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function s = stirling2(n,k) returns the Stirling number of the second
% kind for the non-negative integers n,k.
%
% This function is made for the FYS-3740 Project Paper in Applied Physics
% and Mathematics.
%
% Made by Torgeir Brenn, november 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = stirling2(n,k)
    
    %Using the explicit definition
    s = 0;
    
    for i=0:k
        s = s + (-1)^(k-i)*nchoosek(k,i)*i^n;
    end
    
    s = s/factorial(k);

end