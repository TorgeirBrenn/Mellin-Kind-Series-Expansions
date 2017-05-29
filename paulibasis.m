%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION sp = paulibasis(sl) returns the Pauli basis target vector given
%the target vector in lexicographic basis. sl is a 3x1 vector.
%
%%%     %%%     %%%     %%%     %%%
%Made by Torgeir Brenn, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sp = paulibasis(sl)
    s = size(sl);

    if (s(1) ~= 3) || (s(2) ~= 1)
        disp 'Input must be a 3x1 (column) vector.'
        return
    end
    P = 1/sqrt(2)*[1 0 1; 1 0 -1; 0 sqrt(2) 0];
    sp = P*sl;
end
