%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION C = HermOutProd(v) returns the Hermitian outer product of the
%input vector v of size Nx1.
%
%Made by Torgeir Brenn, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function C = hermoutprod(v)
    s = size(v);
    if s(2)~=1
        disp 'Input must be (column) vector of size Nx1'
        return
    end
    
    C = v*conj(transpose(v));
    %C = v*v' also would work since the '-operator is the conjugate
    %transpose.
end