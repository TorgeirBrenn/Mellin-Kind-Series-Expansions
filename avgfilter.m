%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION A = avgfilter(I, s) applies a local averaging filter of size s to
%the input image I and returns the filtered image in A. The filter uses the
%arithmetic average. I is an MxN image and s is a 2x1 vector of sizes. NOT
%tested for dim(I)>2, but it might work...
%
%This function is used in the second assignment in FYS 3023.
%%%     %%%     %%%     %%%     %%%
%Made by Candidate 6 - 2016-10-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = avgfilter(I,s)

    A = imfilter(I, ones(s)/numel(s), 'replicate', 'same');
    
    %ones(s)/numel(s) is simply the filter, with all elements equal and
    %summing to 1
    
    %'replicate' ensures that border values are replicated outside I in 
    %order to perform the filtering on the edges of the Image
    
    %Same ensures that the size(A) = size(I)
end