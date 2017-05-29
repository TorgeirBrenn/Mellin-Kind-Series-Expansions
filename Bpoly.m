%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function B = Bpoly(x) returns the complete exponential Bell polynomial of
% degree N evaluated at the values in x, where x is a Nx1 vector.
%
% This function is made for the FYS-3740 Project Paper in Applied Physics
% and Mathematics.
%
% Made by Torgeir Brenn, november 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = Bpoly(x)
    N = numel(x);
    
    %The Bell polynomials can be expressed as the determinant of a matrix.
    
    A = zeros(N);
    A = A + diag(repmat(-1, 1, N-1),-1);
    A = A + x(1)*eye(N);
    
    for i = 1:N %Row index
        for j = (i+1):N %Col index
            A(i,j) = nchoosek(N-i, j-i)*x(j-i+1);
        end
    end
    
    B = det(A);
    
end