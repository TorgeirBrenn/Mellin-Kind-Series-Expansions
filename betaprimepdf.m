%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%betaprimepdf : beta prime distribution probability density function
%
%f = betaprimepdf(x, a1, a2, b) returns the beta prime distribution PDF at 
%the values in x. 
%
%INPUT
%x  : Values at which the PDF should be evaluated, Nx1 vector.
%a1 : a1 shape parameter of the beta prime distribution, positive double.
%a2 : a2 shape parameter of the beta prime distribution, positive double.
%b  : b scale parameter of the beta prime distribution, positive double.
%
%Last update: 2017-04-13
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f=betaprimepdf(x, a1, a2, b)
    
    if isinf(gamma(a1+a2))
        f = (a1-1)*log(b*x)-(a1+a2)*log(1+b*x);
        f = exp(f);
        %f = (b*x).^(a1-1)./(1+b*x).^(a1+a2);
        f = normalizepdf(f,x(2)-x(1)); %Assumes that the values in x are monotonically rising with identical intervals.
    else
        f = b/beta(a1,a2)*(b*x).^(a1-1)./(1+b*x).^(a1+a2);
    end

    %Below is the logarithmic evaluation, which doesn't help much as
    %gamma(a1+a2) is the major bottleneck which can be infinite.
    %f = log(b)+log(gamma(a1+a2))-log(gamma(a1))-log(gamma(a2))+(a1-1)*log(b*x)-(a1+a2)*log(1+b*x);
    %f = exp(f);
    
     
end