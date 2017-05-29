%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hermitefast : faster (non-symbolic) version of hermiteHprob with same
%inputs, but limited to N <= 18.
%
%h = hermitefast(N,x) mimics the function of hermiteHprob(), but is limited
%to N<=18. The 18 first generalised polynomials are explicitly written in
%this function.
%
%INPUT
%N : Order of the Hermite polynomial.
%x : Datapoints on which the Hermite polynomial should be evaluated.
%
%OUTPUT
%h : The functional value of the Hermite polynomial.
%
%Last update: 2017-02-22
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function h = hermitefast(N,x)
    switch N
        case 0
            h = ones(size(x));
        case 1
            h = x;
        case 2
            h = x.^2 - 1;
        case 3
            h = x.^3 - 3*x;
        case 4
            h = x.^4 - 6*x.^2 + 3;
        case 5
            h = x.^5 - 10*x.^3 + 15*x;
        case 6
            h = x.^6 - 15*x.^4 + 45*x.^2 - 15;
        case 7
            h = x.^7 - 21*x.^5 + 105*x.^3 - 105*x;
        case 8
            h = x.^8 - 28*x.^6 + 210*x.^4 - 420*x.^2 + 105;
        case 9
            h = x.^9 - 36*x.^7 + 378*x.^5 - 1260*x.^3 + 945*x;
        case 10
            h = x.^10 - 45*x.^8 + 630*x.^6 - 3150*x.^4 ...
                + 4725*x.^2 - 945;
        case 11
            h = x.^11 - 55*x.^9 + 990*x.^7 - 6930*x.^5 ...
            + 17325*x.^3 - 10395*x;
        case 12
            h = x.^12 - 66*x.^10 + 1485*x.^8 - 13860*x.^6 ...
                + 51975*x.^4 - 62370*x.^2 + 10395;
        case 13
            h = x.^13 - 78*x.^11 + 2145*x.^9 - 25740*x.^7 + 135135*x.^5 ...
                - 270270*x.^3 + 135135*x;
        case 14
            h = x.^14 - 91*x.^12 + 3003*x.^10 - 45045*x.^8 ...
                + 315315*x.^6 - 945945*x.^4 + 945945*x.^2 - 135135;
        case 15
            h = x.^15 - 105*x.^13 + 4095*x.^11 - 75075*x.^9 ...
                + 675675*x.^7 - 2837835*x.^5 + 4729725*x.^3 - 2027025*x;
        case 16
            h = x.^16 - 120*x.^14 + 5460*x.^12 - 120120*x.^10 ...
                + 1351350*x.^8 - 7567560*x.^6 + 18918900*x.^4 ...
                - 16216200*x.^2 + 2027025; 
        case 17
            h = x.^17 - 136*x.^15 + 7140*x.^13 - 185640*x.^11 ...
                + 2552550*x.^9 - 18378360*x.^7 + 64324260*x.^5 ...
                - 91891800*x.^3 + 34459425*x;
        case 18
            h = x.^18 - 153*x.^16 + 9180*x.^14 - 278460*x.^12 ...
                + 4594590*x.^10 - 41351310*x.^8 + 192972780*x.^6 ...
                - 413513100*x.^4 + 310134825*x.^2 - 34459425;
        otherwise
            disp 'ERROR: N must take an integer value in {0,...,18}.'
    end
end