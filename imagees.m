%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%imagees : displays an image using the 2% ENVI stretch
%
%im = imagees(I) displays an image using the 2% ENVI stretch, by simply
%taking an input image, stretching it so that the top 2% and bottom 2% are
%clipped at the maximum and minimum intensities.
%
%INPUT
%I : NxMx3 image
%
%OUTPUT
%ih : Image handle
%
%Last update: 2017-05-13
%Made by Torgeir Brenn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ih = imagees(I)
    s = size(I);
    Y = reshape(I, numel(I), 1); %Turn into vector before sorting
    Y = sort(Y);
    low = Y(round(numel(Y)*0.02)); %The 2%-threshhold
    high= Y(round(numel(Y)*0.98)); %The 98%-threshhold
    
    I(I<low) = low; %Perform threshholds
    I(I>high) = high; %Perform threshholds
    ih = image([1, s(1)], [1, s(2)], I, 'CDataMapping', 'scaled'); %This stretches the image to use the full range
end