%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Timage : displays the Pauli basis RCS
%
%im = Timage(T) takes the coherence matrices of SAR data and displays the
%diagonals as blue, red and green, respectively. This amounts to displaying
%the Pauli basis RCS. T can be single or multi look. Displaying uses the
%default 2% ENVI stretch, on a logarithmic scale.
%
%INPUT
%T : SAR coherency matrices, NxMx3x3
%
%OUTPUT
%ih : Image handle
%
%Last update: 2017-05-13
%Made by Torgeir Brenn, heavily inspired by equivalent function made by
%Stian N. Anfinsen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ih = Timage(T)
    s = size(T);
    if (s(3) ~= 3) || (s(4) ~= 3)
        disp 'Input must be a MxNx3x3 matrix.'
        return
    end

    B = 4*pi*real(T(:,:,1,1)); %Rough surface    -> Blue 
    R = 4*pi*real(T(:,:,2,2)); %Dihedral         -> Red
    G = 4*pi*real(T(:,:,3,3)); %Rotated dihedral -> Green
    
    %The real() function is used to get rid of 0.000*i-terms. The diagonals
    %of a Hermitian outer product are strictly real. The factor of 4*pi is
    %to convert it to RCS values, but is strictly not needed because the
    %image is scaled later.
    
    R(R<eps) = eps; %Avoid taking log10(0)
    R = 10*log10(R);%Convert to dB.
    
    G(G<eps) = eps;
    G = 10*log10(G);
    
    B(B<eps) = eps;
    B = 10*log10(B);

    %Scaling the bands to [0,1]
    
    eps2 = 10*log10(eps);
    Rmin = min(min(R(R>eps2)));
    Rmax = max(max((R)));
    Gmin = min(min(G(G>eps2)));
    Gmax = max(max((G)));
    Bmin = min(min(B(B>eps2)));
    Bmax = max(max(B));

    R = (R - Rmin) ./ (Rmax - Rmin);
    G = (G - Gmin) ./ (Gmax - Gmin);
    B = (B - Bmin) ./ (Bmax - Bmin);
    
    figure('color', [1 1 1])
    RGB = zeros(s(1),s(2),3);
    RGB(:,:,1) = R;
    RGB(:,:,2) = G;
    RGB(:,:,3) = B;
    
    ih = imagees(RGB); %Custom function, performs default ENVI 2% linear stretch

end