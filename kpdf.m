%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kpdf : K-distribution probability density function
%
%f = kpdf(x, mu, M, L) returns the K-distribution PDF at the values in x.
%
%INPUT
%x  : Values at which the PDF should be evaluated, Nx1 vector.
%mu : Location parameter of K distribution, positive double.
%M  : First shape (texture) parameter of K distribution, positive double.
%L  : Second shape parameter of K distribution (equiv. number of looks),
%     positive double.
%
%Last update: 2016-08-30
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f=kpdf(x, mu, M, L)
       %Declare brevity variables
       a = M-L;
       b = L+M-1;
       c = L*M/mu;
       
       %Getting the results of the modified Bessel function
       K = besselk(a,2*sqrt(c*x));
       f = 2*c^((b+1)/2)*x.^((b-1)/2)/(gamma(L)*gamma(M)).*K;    
end