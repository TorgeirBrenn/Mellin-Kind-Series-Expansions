%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script was made to test the performance of the Mellin kind series
% expansions on real data. It is separated into small pieces to allow for
% flexible testing, and also to produce the T matrices.
%
% For compiling .c-files:
% mex -g ...
% /Users/torgeirbrenn/Documents/Skole/Masteroppgave/C/mexKdistGradOpt.c ...
% -I/usr/local/include -lgsl -lgslcblas -lm
%
% Made by Torgeir Brenn, 2017. For the custom functions, the author is
% specified within those functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(groot, 'defaultFigureColor', [1 1 1]);   %White figure background
set(groot,'defaulttextinterpreter','latex'); %Figure text
set(0,'defaultAxesFontSize',16);
set(groot,'defaultFigurePosition', [250 0 560 420])
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/functions'


%% Loading the raw data
%Channels: sHH, sHV, sVH, and sVV
load '/Users/torgeirbrenn/Documents/Skole/data/sanfransisco.mat'


%% Creating the single look coherency matrix
T=zeros(2800, 2800, 3, 3);
sCR = 1/sqrt(2)*(sHV + sVH);
for i = 1:2800
    for j = 1:2800
        v = [sHH(i, j); sCR(i, j); sVV(i,j)];
        v = PauliBasis(v); %The difference from part 6b)
        T(i,j,:,:) = HermOutProd(v);
    end
end

%% Creating multi looked coherency matrix
L = 16; %Must be a squared interger in this case, e.g. 4, 16, 25 etc.
TL16 = zeros(size(T));
for i = 1:3
    for j = 1:3
        TL16(:,:,i,j) = avgfilter(T(:,:,i,j),sqrt(L)); %Window size 5x5
    end
end

%% Saving the coherency matrices we made
save('/Users/torgeirbrenn/Documents/Skole/data/sanfransiscoT.mat', ...
    'T', 'TL4', 'TL16')

%% Reseting the workspace and closing all figures
close all
clear variables
load '/Users/torgeirbrenn/Documents/Skole/data/sanfransiscoT.mat'

%% Water Region of Interest (RoI)
%[70 - 450, 130 - 850]
left = 70;
right = 450;
top = 130;
bottom = 850;

%% Alternative RoI, mixed area
left = 1000;
right = 1500;
top = 2000;
bottom = 2500;

%% Creating the RoI
%Region of Interest
RoI = zeros(2800, 2800);
RoI(top:bottom, left:right) = 1;

%% Displaying the entire image
ih = Timage(T);

%% Possible adjustments to images
%ih.Parent.XLabel.String = '$x$'; 
%ih.Parent.YLabel.String = '$y$'; 
ih.Parent.XTickLabel = cell(0); %Removes the labels
ih.Parent.YTickLabel = cell(0);
ih.Parent.XTick = zeros(0); %Removes the ticks
ih.Parent.YTick = zeros(0);

%% Displaying a frame over the RoI
hold on;
plot(left:right, top*ones(size(left:right)), 'r', 'LineWidth', 2)
plot(left:right, bottom*ones(size(left:right)), 'r', 'LineWidth', 2)
plot(left*ones(size(top:bottom)), top:bottom, 'r', 'LineWidth', 2)
plot(right*ones(size(top:bottom)), top:bottom, 'r', 'LineWidth', 2)

%% Displaying only the part of the image within the RoI
ih = Timage(T(top:bottom, left:right, :,:));

%% Extracting the data
tmp = real(TL16(:,:,1,1)); %The blue channel, discard the zero complex parts
data = tmp(RoI==1); %All blue intensities within WaterROI, Nx1

%% Making a histogram
figure; hold on;
histogram(data, 'Normalization', 'pdf')

%% Creating the x grid
dx = max(data)/1e4;
x  = dx:dx:max(data);

%% Making a KDE
fPD = fitdist(data, 'Kernel', 'kernel', 'epanechnikov', ...
    'support', 'positive');
f = pdf(fPD, x);
plot(x, f, 'k', 'LineWidth', 1)

%% Fitting a K distributions to the image
[L, M, m] = kdistfit(data);
k = kpdf(x, m, M, L);
hold on;
plot(x, k, 'r', 'LineWidth', 1)

%% Fitting the MKGK series expansion
mkgk = mkgkfit(x, data, 2);
plot(x, mkgk, 'g' , 'LineWidth', 1)
