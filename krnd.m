%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%krnd : K-distribution random number generator
%
%r = krnd(N, mu, M, L) returns N random K-distributed samples.
%
%INPUT
%N  : Number of samples in output, positive integer.
%mu : Location parameter of K distribution, positive double.
%M  : First shape (texture) parameter of K distribution, positive double.
%L  : Second shape parameter of K distribution (equiv. number of looks),
%     positive double.
%
%OUTPUT
%r  : Random samples, Nx1 vector.
%
%Last update: 2016-09-07
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function r=krnd(N, mu, M, L)

       r = gamrnd(L, 1/L, [N 1]).*gamrnd(M, mu/M, [N 1]);       

end