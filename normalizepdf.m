%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normalizepdf : normalizes input to a PDF (positive everywhere, sums to 1)
%
%nf = normalizepdf(f, dx) returns the normalized probability density 
%version of function f. That is, f can have negative values and sum to
%values other than 1. nf is the normalized version of f in the sense that
%it sums to 1 and is positive everywhere. Note that nf=0 is not permitted
%here in order to facilitate use of the Kullback-Leibler distance measure.
%
%INPUT
%f  : Input function f(x) .
%dx : Grid resolution of x. 
%
%OUTPUT
%nf : Normalized PDF version of f(x).
%
%Last update: 2017-02-10
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nf = normalizepdf(f,dx)
    nf = f;
    nf(nf<=0) = eps;
    
    nf = nf/(sum(nf)*dx);
end