%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pdfname :
%
%
%INPUT
%dist    : Distribution, string.
%param   : Parameter struct.
%compact : Logical, 1: More compact version, 0: default.
%
%OUTPUT
%s : String: PDF(x;param1, param2,...) prepared for LaTeX formatting.
%
%Last update: 2017-05-18
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = pdfname(dist, param, compact)
if nargin == 2
    compact = 0;
end
if compact == 0
    switch dist
        case "gamma"
            s = ['$\gamma(x;L=',num2str(param.L), ',m=', ...
                num2str(param.m), ')$'];
        case "igam"
            s = ['$\gamma^{-1}(x;L=', num2str(param.L), ',m=', ...
                num2str(param.m), ')$'];
        case "k"
            s = ['$K(x;L=', num2str(param.L), ',m=', num2str(param.m), ...
                ',M=', num2str(param.M), ')$'];
        case "ggd"
            s = ['G$\Gamma$D$(x;L=', num2str(param.L), ',m=', ...
                num2str(param.m), ',\nu=', num2str(param.nu), ')$'];
        case "g0"
            s = ['$G^0(L=', num2str(param.L), ',g=', num2str(param.g), ...
                ',M=', num2str(param.M), ')$'];
        otherwise
            disp 'Warning: function pdfname called with invalid dist'
    end
elseif compact == 1
    switch dist
        case "gamma"
            s = ['$\gamma(x;\!L{=}',num2str(param.L), ',\!m{=}', ...
                num2str(param.m), ')$'];
        case "igam"
            s = ['$\gamma^{-1}(x;\!L{=}', num2str(param.L), ',\!m{=}', ...
                num2str(param.m), ')$'];
        case "k"
            s = ['$K(x;\!L{=}', num2str(param.L), ',\!m{=}', ...
                num2str(param.m), ',\!M{=}', num2str(param.M), ')$'];
        case "ggd"
            s = ['G$\Gamma$D$(x;\!L{=}', num2str(param.L), ',\!m{=}', ...
                num2str(param.m), ',\!\nu{=}', num2str(param.nu), ')$'];
        case "g0"
            s = ['$G^0(x;\!L{=}', num2str(param.L), ',\!g{=}', ...
                num2str(param.g), ',\!M{=}{-}', num2str(-param.M), ')$'];
        otherwise
            disp 'Warning: function pdfname called with invalid dist'
    end
end
end