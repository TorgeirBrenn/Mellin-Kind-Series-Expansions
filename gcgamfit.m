%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%gcgamfit : fits a Gamma PDF to data using the Edgeworth expansion
%with the gamma distribution and Laguerre polynomials.
%
%
%g = gcgamfit(x,d, N) evaluates the Gamma distribution at the points
%in x, using parameters based on the method of moments, corrected
%with a Gram-Charlier expansion using the N first (classical) empirical 
%moments.
%
%INPUT
%x : Points at the which the Gamma PDF will be evaluated, vector.
%d : Data, vector.
%N : Order of the highest moment which will be used in the expansion. N up
%    to 8 supported.
%
%OUTPUT
%g : The fitted Gamma PDF expansion evaluated at points in x.
%
%Last update: 2017-02-21
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = gcgamfit(x, d, N)
    if N < 3
        %disp 'Warning: N<3 is just the gamma distribution fitted with the MoM.'
        N = 2; %Just to get the required moments for the MoM. 
    end
    
    if N > 8
        disp 'Warning: N>8 not supported. N=8 will be used.'
        N = 8;
    end
     
    %Moments
    m = zeros(N,1);
    for n = 1:N
        m(n) = mean(d.^n);
    end
    
    %Kernel parameters via MoM
    a = m(1)^2/(m(2)-m(1)^2)-1;
    b = (a+1)/m(1);
    
    %Kernel
    g = gampdf(x, a+1, 1/b);
    
    %A trick
    mb = m;
    for n = 1:N
        mb(n) = mb(n)*b^n;
        mb(n:N) = mb(n:N)/(a+n);
        %This greatly simplifies the coefficients.
    end
    
    %Coefficients, computed using laguerreL(), substituting x^n with m(n).
    %This is hard coded to save time.
    c = zeros(N,1); %c0 = 1, c1 = c2 = 0 always
    
    if N>=3
        c(3) = -mb(3) + 1;
    end
    if N>=4
        c(4) = mb(4)-4*mb(3)+3;
    end
    if N>=5
        c(5) = -mb(5)+5*mb(4)-10*mb(3)+6;
    end
    if N>=6
        c(6) = mb(6)-6*mb(5)+15*mb(4)-20*mb(3)+10;
    end
    if N>=7
        c(7) = -mb(7)+7*mb(6)-21*mb(5)+35*mb(4)-35*mb(3)+15;
    end
    if N>=8
        c(8) = mb(8)-8*mb(7)+28*mb(6)-56*mb(5)+70*mb(4)-56*mb(3)+21;
    end
    
    %Applying the Gram-Charlier expansion
    for n = 3:N
        g = g + c(n)*g.*laguerrefast(n, a, b*x);
    end
    
    

end