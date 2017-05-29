%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%igamrnd : Inverse gamma distribution random number generator
%
%r = igamrnd(N, a, b) returns N random inverse gamma distributed samples.
%See Nicolas, 2002 for details about the distribution.
%
%INPUT
%N : Number of samples in output, positive integer.
%L : Shape parameter of the inverse gamma distribution, positive double.
%b : Scale parameter of the inverse gamma distribution, positive double.
%
%OUTPUT
%r  : Random samples, Nx1 vector.
%
%Last update: 2017-03-21
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r=igamrnd(N, L, b)

    a = L-1; %Old parameter choices
    r = 1./gamrnd(a+1, 1/b, [N 1]);

end