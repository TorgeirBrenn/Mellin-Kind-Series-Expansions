%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generatedata : generates data from one of the supported distributions
%
%data = generatedata(n, dist, param) generates n pseudorandom data points
%drawn from the 'dist' distribution with parameters contained in the struct
%'param'.
%
%INPUT
%n     : Number of data points desired.
%dist  : Distribution type desired.
%param : Parameter struct.
%
%OUTPUT
%data : The pseudorandom data, nx1 vector.
%
%Last update: 2017-05-05
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = generatedata(n, dist, param)

switch dist
    case "k"
        %K distribution, g->0, M>0
        if param.M<0
            disp 'Warning: M<0 in the K distribution, using -M instead.'
            param.M = -param.M;
        end
        data = krnd(n, param.m, param.M, param.L);
    case "gamma"
        %Gamma distribution
        data = gamrnd(param.L, param.m/param.L, n, 1);
    case "g0"
        %G0 distribution, M/mu->0, M<0
        if param.M>0
            disp 'Warning: M>0 in the G0 distribution, using -M instead.'
            param.M = -param.M;
        end
        
        data = g0rnd(n, param.g, param.M, param.L);
    case "ggd"
        %Generalized gamma distribution
        data = ggdrnd(n, param.nu, param.L, param.m/param.L);
    case "igam"
        %Inverse gamma distribution
        data = igamrnd(n, param.L, param.m/param.L);
end
end