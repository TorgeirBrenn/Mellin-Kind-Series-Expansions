%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%g0rnd : G0-distribution random number generator
%
%r = g0rnd(N, g, M, L) returns N random G0-distributed samples. See Frery 
%et al., 1997 for details about the distribution.
%
%INPUT
%N : Number of samples in output, positive integer.
%g : Gamma parameter of G0 distribution, positive double.
%M : M parameter of G0 distribution, negative double.
%L : L parameter of the G0 distribution (equiv. number of looks), positive
%    double.
%
%OUTPUT
%r  : Random samples, Nx1 vector.
%
%Last update: 2017-02-20
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function r=g0rnd(N, g, M, L)
        
       r = gamrnd(L, 1/L, [N 1])./gamrnd(-M, 1/g, [N 1]);       

end