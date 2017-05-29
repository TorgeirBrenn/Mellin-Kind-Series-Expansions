%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mlgamfit : fits a Gamma PDF to data using MLE of the parameters
%
%
%g = mlgamfit(x,d) evaluates the Gamma distribution at the points in x, 
%using the parameters that are the MLEs of the data in d.
%
%INPUT
%x : Points at the which the Gamma PDF will be evaluated, vector.
%d : Data, vector.
%
%OUTPUT
%g : The fitted Gamma PDF evaluated at points in x.
%
%Last update: 2017-01-04
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = mlgamfit(x,d)

    parmhat = gamfit(d); %Built-in function
    g = gampdf(x, parmhat(1), parmhat(2));

end