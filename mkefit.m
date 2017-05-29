%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mkefit : fits a series expansion of the log-normal PDF to data.
%
%f = mkefit(x, d, N) evaluates the fitted PDF at the points in x.
%The output is the series expansion of the log-normal PDF based on log-
%cumulants, introduced by Pastor (2014, 2016). The log-cumulants are based
%on the data d, and the highest order log-cumulant used is N.
%
%INPUT
%x : Points at the which the fitted PDF will be evaluated, vector.
%d : Data, vector.
%N : The highest order of the log-cumulants used in the expansion, N<=6.
%
%OUTPUT
%f : The fitted PDF evaluated at the points in x.
%
%Last update: 2017-02-22
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = mkefit(x, d, N)
    
    if N > 8
        disp 'Warning: N>8 not supported. N=8 will be used.'
        N = 8;
    end
    
    if N < 3
        %disp(['Warning: For N<3, this method is functionally ' ...
        %    'equivalent to molclognormfit().'])
        f = molclognormfit(x,d);
        return
    end

    k = emplc(d, N); %Empirical log-cumulants;
    sigma = sqrt(k(2));
    
    %Making coefficients which will simplify and speed up computations
    a = zeros(N-2,1);
    for i = 1:(N-2)
        a(i) = k(i+2)/((i+1)*(i+2)*sigma^(i+2));
    end
    
    
    y = (log(x)-k(1))/sigma; %Standardizing the argument
    fln = pdf(makedist('Lognormal', 'mu', k(1), 'sigma', sigma), x);
    
    %B_1
    if N >= 3
        f = fln + fln.*hermitefast(3,y)*a(1);
    end
    
    %B_2
    if N >= 4
        f = f + fln.*(hermitefast(6,y)*a(1)^2 + hermitefast(4,y)*a(2))/2;
    end
    
    %B_3
    if N >= 5
        f = f + fln.*(hermitefast(9,y)*a(1)^3 ...
            + 3*hermitefast(7,y)*a(1)*a(2) + hermitefast(5,y)*a(3))/6;
    end
    
    %B_4
    if N >= 6
        f = f + fln.*(hermitefast(12,y)*a(1)^4 ...
            + hermitefast(10,y)*a(1)^2*a(2)*6 ...
            + hermitefast(8,y)*(4*a(1)*a(3)+3*a(2)^2)... 
            + hermitefast(6,y)*a(4))/24;
    end
    
    %B_5
    if N >= 7
        f = f + fln.*(hermitefast(15,y)*a(1)^5 ...
            + hermitefast(13,y)*a(1)^3*a(2)*10 ...
            + hermitefast(11,y)*(15*a(2)^2*a(1)+10*a(3)*a(1)^2) ...
            + hermitefast(9,y)*(10*a(3)*a(2)+5*a(4)*a(1)) ...
            + hermitefast(7,y)*a(5))/120;
    end
    
    %B_6
    if N >= 8
        f = f + fln.*(hermitefast(18,y)*a(1)^6 ...
            + hermitefast(16,y)*15*a(1)^4*a(2) ...
            + hermitefast(14,y)*(20*a(3)*a(1)^3+45*a(2)^2*a(1)^2) ...
            + hermitefast(12,y)*(15*a(2)^3+60*a(1)*a(2)*a(3)+15*a(4)*a(1)^2) ...
            + hermitefast(10,y)*(10*a(3)^2+15*a(4)*a(2)+6*a(5)*a(1)) ...
            + hermitefast(8,y)*a(6))/720;
    end
    
    %hermitefast is an optimised version of hermiteHprob, but is limited to
    %input n<=12. 
end