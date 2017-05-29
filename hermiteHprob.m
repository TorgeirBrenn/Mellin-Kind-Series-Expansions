%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hermiteHprob : Returns the probabilist's Hermite polynomials.
%
%h = hermiteHprobs(N,x) returns the N-th probabilists Hermite polynomial
%evaluated at the values specified in x. hermiteHprob calls the built-in
%function hermiteH and converts that result to the probabilists
%
%INPUT
%N : Order of the Hermite polynomial.
%x : Datapoints on which the Hermite polynomial should be evaluated.
%
%OUTPUT
%h : The functional value of the Hermite polynomial.
%
%Last update: 2017-02-13
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h=hermiteHprob(N,x)
    h = 2^(-N/2)*hermiteH(N,x/sqrt(2));
end