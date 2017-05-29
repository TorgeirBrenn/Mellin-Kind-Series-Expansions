%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%g0pdf : g0-distribution probability density function
%
%f = g0pdf(x, g, M, L) returns the G0-distribution PDF at the values in x.
%See Frery et al. 1997 for details.
%
%INPUT
%x : Values at which the PDF should be evaluated, Nx1 vector.
%g : Gamma parameter of the G0 distribution, positive double.
%M : M parameter of G0 distribution, negative double.
%L : L parameter of the G0 distribution (equiv. number of looks), positive
%    double.
%
%Last update: 2017-02-20
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f=g0pdf(x, g, M, L)

    f = L^L*gamma(L-M)*x.^(L-1)/(g^M*gamma(L)*gamma(-M)).*(g+L*x).^(M-L);

end