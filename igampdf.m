%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%igampdf : generalized gamma distribution probability density function
%
%f = igampdf(x, a, b) returns the inverse gammadistribution PDF at the 
%values in x. See Nicolas, 2002 for details.
%
%INPUT
%x : Values at which the PDF should be evaluated, Nx1 vector.
%L : Shape parameter of the inverse gamma distribution, positive dougle.
%b : Scale parameter of the inverse gamma distribution, positive double.
%
%
%Last update: 2017-03-21
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f=igampdf(x, L, b)

    a = L-1; %Old parameter choices
    f = b^(a+1)/gamma(a+1)*x.^(-a-2).*exp(-b./x);

end