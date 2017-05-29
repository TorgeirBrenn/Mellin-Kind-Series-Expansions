%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script was made to assess the performance of the Mellin kind series
% expansions on real data.
%
% For compiling .c-files:
% mex -g ...
% /Users/torgeirbrenn/Documents/Skole/Masteroppgave/C/mexKdistGradOpt.c ...
% -I/usr/local/include -lgsl -lgslcblas -lm
%
% Made by Torgeir Brenn, 2017. For the custom functions, the author is
% specified within those functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear variables
load '/Users/torgeirbrenn/Documents/Skole/data/sanfransiscoT.mat'

set(groot, 'defaultFigureColor', [1 1 1]);   %White figure background
set(groot,'defaulttextinterpreter','latex'); %Figure text
set(0,'defaultAxesFontSize',18);
set(groot,'defaultFigurePosition', [0 0 2800 2800])
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/functions'

Nd = [1e2 1e3 5e3]; %Number of data points
Ni = 1e2; %Number of iterations
Nt  = 4;   %Highest order log-cumulant corrected for

saveImages   = 0; %1 => stored images will be overwritten
saveTextFile = 1;

CurrentRoI = 'Water'; %Supported: 'Water', 'Park', 'Urban', 'Test'

switch CurrentRoI %RoI = Region of Interest
    case 'Water'
        % Water, [70 - 450, 130 - 850]
        left = 70;
        right = 450;
        top = 130;
        bottom = 850;
        ch = 1; %Blue channel
    case 'Park'
        % Park, [990 - 1106, 1910 - 2040]
        left = 970;
        right = 1106;
        top = 1900;
        bottom = 2030;
        ch = 3; %Green channel
    case 'Urban'
        left = 1580;
        right = 1783;
        top = 1400;
        bottom = 1803;
        ch = 2; %Red channel
    case 'Test'
        % Testing
        left = 1280;
        right = 1500;
        top = 2000;
        bottom = 2500;
        ch = 1;
end

%%% CREATING THE FIGURES %%%

%Making the actual RoI
RoI = zeros(2800, 2800);
RoI(top:bottom, left:right) = 1;

%Displaying the entire image
ihi = Timage(T); hold on
ihi.Parent.XTickLabel = cell(0); %Removes the labels
ihi.Parent.YTickLabel = cell(0);
ihi.Parent.XTick = zeros(0); %Removes the ticks
ihi.Parent.YTick = zeros(0);
ihi.Parent.PlotBoxAspectRatio = [1 1 1]; %Avoids stretched pixels

%Making a border around the RoI
plot(left:right,    top*ones(size(left:right)), 'r', 'LineWidth', 3)
plot(left:right, bottom*ones(size(left:right)), 'r', 'LineWidth', 3)
plot(left*ones(size(top:bottom)),   top:bottom, 'r', 'LineWidth', 3)
plot(right*ones(size(top:bottom)),  top:bottom, 'r', 'LineWidth', 3)

%Making lines to the RoI plot (which will be on the right)
plot([right, 2800], [bottom, 2800], 'r--', 'LineWidth', 3)
plot([right, 2800], [top,       2], 'r--', 'LineWidth', 3)

%Saving the entire image including the RoI frame
if saveImages == 1
    tmp = getframe;
    imwrite(tmp.cdata, ['/Users/torgeirbrenn/Documents/', ...
        ['Skole/Masteroppgave/Innlevering/figures/SARSingleImage', ...
        CurrentRoI, '.png']])
end
%Making the RoI image, resized to match height of the entire image
tmp = getframe;
s = size(tmp.cdata);
figure;
ihr = imshow(imresize(imcrop(ihi.CData, ...
    [left top right-left+1 bottom-top+1]), [s(1) NaN])); hold on
s = size(ihr.CData);

%RoI image border
plot(2:s(2),        2*ones(s(2)-1,1), 'r--', 'LineWidth', 3)
plot(2:s(2), (s(1)-1)*ones(s(2)-1,1), 'r--', 'LineWidth', 3)
plot(       2*ones(s(1)-1,1), 2:s(1), 'r--', 'LineWidth', 3)
plot((s(2)-1)*ones(s(1)-1,1), 2:s(1), 'r--', 'LineWidth', 3)
%truesize;

%Saving the RoI
if saveImages == 1
    tmp = getframe;
    imwrite(tmp.cdata, ['/Users/torgeirbrenn/Documents/', ...
        ['Skole/Masteroppgave/Innlevering/figures/SARSingleRoI', ...
        CurrentRoI, '.png']])
end

%Making a new image which is a combination of the two above
tmp1 = getframe(ihi.Parent.Parent.CurrentAxes);
tmp2 = getframe(ihr.Parent.Parent.CurrentAxes);
figure;
imshow([imresize(tmp1.cdata, [s(1) NaN]), tmp2.cdata]); hold on;

%Saving the combined image
if saveImages == 1
    tmp = getframe;
    imwrite(tmp.cdata, ['/Users/torgeirbrenn/Documents/', ...
        ['Skole/Masteroppgave/Innlevering/figures/SARSingleCombined', ...
        CurrentRoI, '.png']])
end

%Extracting the data
tmp = real(T(:,:,ch,ch)); %Discard the zero complex parts
dataL1 = tmp(RoI==1); %All intensities within the RoI, Nx1

tmp = real(TL4(:,:,ch,ch)); %Discard the zero complex parts
dataL4 = tmp(RoI==1); %All intensities within the RoI, Nx1

tmp = real(TL16(:,:,ch,ch)); %Discard the zero complex parts
dataL16 = tmp(RoI==1); %All intensities within the RoI, Nx1

% Making histograms for L=1,4,16
figure('Position', [0 500 1200 400]); hold on;

%For the subplot positions
lt = 0.07; %left
bt = 0.12; %bottom
wt = 0.275; %width
ht = 0.81; %height
hg = wt+0.05; %horizontal gap
vg = ht+0.07; %vertical gap

subplot('Position', [lt+0*hg, bt+0*vg, wt, ht])
auxfunction1(dataL1);
title 'Histogram and Estimates, $L = 1$'
subplot('Position', [lt+1*hg, bt+0*vg, wt, ht])
auxfunction1(dataL4);
title 'Histogram and Estimates, $L = 4$'
subplot('Position', [lt+2*hg, bt+0*vg, wt, ht])
auxfunction1(dataL16);
title 'Histogram and Estimates, $L = 16$'

%Saving the histograms and PDF estimates
if saveImages == 1
    tmp = getframe(gcf);
    imwrite(tmp.cdata, ['/Users/torgeirbrenn/Documents/', ...
        ['Skole/Masteroppgave/Innlevering/figures/SARSingleHistogram', ...
        CurrentRoI, '.png']])
end

%%% CREATING THE DATA TABLE %%%

%Making KDEs as the target PDFs, disregarding values where f(x)=0
%L = 1
[xL1, fL1] = targetpdf('kernel', dataL1);

%L = 4
[xL4, fL4] = targetpdf('kernel', dataL4);

%L = 16
[xL16, fL16] = targetpdf('kernel', dataL16);

%Setting up the container for the results
results = struct(...
    'L1',  zeros(8,3),...
    'L4',  zeros(8,3),...
    'L16', zeros(8,3)...
    );

for nd=1:numel(Nd)
    results.L1(:, nd)  = auxfunction2(dataL1, xL1, fL1, ...
        Nd(nd), Ni, Nt);
    results.L4(:, nd)  = auxfunction2(dataL4, xL4, fL4, ...
        Nd(nd), Ni, Nt);
    results.L16(:, nd) = auxfunction2(dataL16, xL16, fL16, ...
        Nd(nd), Ni, Nt);
end

if saveTextFile == 1
    %Printing the results in the table format, to a text file
    printArray = [results.L1, results.L4, results.L16];
    formatSpec = ['%s & $%3.2e}$ & $%3.2e}$ & $%3.2e}$ & $%3.2e}$ & ', ...
        '$%3.2e}$ & $%3.2e}$ & $%3.2e}$ & $%3.2e}$ & $%3.2e}$', ...
        '\\\\ \\hline \n'];
    
    %fileID = fopen(['../results/SARSingleImage', CurrentRoI, '.txt'],'w');
    fileID = fopen('../results/SARSingleImageWaterBhatta.txt','w');
    
    names = {'MKGK kernel'; 'MKGK series'; 'MKLK kernel'; ...
        'MKLK series, $N=3$'; 'MKLK series'; 'MKE series'; ...
        'MKBK kernel'; 'MKBK series'};
    
    for m = 1:8
        fprintf(fileID, formatSpec, [string(names(m)), printArray(m,:)]);
    end
    
    %Instructions:
    %Replace 'e-0' with '{\cdot}10^{-', and if necessary
    %replace 'e+0' with '{\cdot}10^{' and pad with \hspace{5pt}
end

%%% AUXILLARY FUNCTIONS %%
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/'
%auxfunction1: plots histogram and the MK series expansions PDF estimates
%of a dataset
function auxfunction1(dat)

%x grid
dx = max(dat)/1e4;
x  = dx:dx:max(dat);

Nt = 4;
PDFest = zeros(numel(x), 7); %PDF estimates

PDFest(:, 1) = mkgkfit(x, dat, 2);
PDFest(:, 2) = mkgkfit(x, dat, Nt);
PDFest(:, 3) = mklkfit(x, dat, 2);
PDFest(:, 4) = mklkfit(x, dat, Nt);
PDFest(:, 5) =  mkefit(x, dat, Nt);
PDFest(:, 6) = mkbkfit(x, dat, 2);
PDFest(:, 7) = mkbkfit(x, dat, Nt);

names = {'MKGK kernel'; 'MKGK series'; 'MKLK kernel'; ...
    'MKGK series'; 'MKE series'; 'MKBK kernel'; 'MKBK series'};
style = {'-'; '--'; '-'; '--'; '-.'; '-'; '--'};
color = {[0 0.4470 0.7410]; [0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; ...
    [0.8500 0.3250 0.0980]; [0.8500 0.3250 0.0980]; ...
    [0.4660 0.6740 0.1880]; [0.4660 0.6740 0.1880]};

%Plotting
hold on; grid on; box on;
histogram(dat, 'Normalization', 'pdf', 'EdgeAlpha', 0, ...
    'FaceColor', [0.4940 0.1840 0.5560], 'FaceAlpha', 0.4)
for m = 1:7
    plot(x, PDFest(:, m), 'LineStyle', style{m}, ...
        'Color', color{m}, 'LineWidth', 1.2)
end
lh = legend(['Histogram', names']);
set(lh,'Interpreter','latex')

%Cosmetics
ax = gca;
ax.XLim(1) = 0;
ax.XLim(2) = x(end);
ax.YLim(1) = 0;
xlabel $x$
ylabel $\hat{f}(x)$

end
%auxfunction2: all methods applied to 'Nd' data points drawn from 'dat',
%returns distance to target PDF 'tgt', on grid 'x', mean of 'Ni' iterations
function res = auxfunction2(dat, x, tgt, Nd, Ni, Nt)

%Setting up a container for the results
distances = zeros(Ni,8);

dx = x(2)-x(1); %Need this to standardize the distance

wb = waitbar(0, 'Computing');
for ni = 1:Ni
    %Draw Nd random sample without replacement
    datasubset = datasample(dat, Nd, 'Replace', false);
    
    %Computing the distances
    distances(ni,1) = bhattadist(tgt, mkgkfit(x, datasubset, 2), dx);
    distances(ni,2) = bhattadist(tgt, mkgkfit(x, datasubset, Nt), dx);
    distances(ni,3) = bhattadist(tgt, mklkfit(x, datasubset, 2), dx);
    distances(ni,4) = bhattadist(tgt, mklkfit(x, datasubset, 3), dx);
    distances(ni,5) = bhattadist(tgt, mklkfit(x, datasubset, Nt), dx);
    distances(ni,6) = bhattadist(tgt,  mkefit(x, datasubset, Nt), dx);
    distances(ni,7) = bhattadist(tgt, mkbkfit(x, datasubset, 3), dx);
    distances(ni,8) = bhattadist(tgt, mkbkfit(x, datasubset, Nt), dx);
    waitbar(ni/Ni, wb)
end
close(wb)
res = mean(distances,1);
end