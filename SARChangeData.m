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
set(0,'defaultAxesFontSize',16);
set(groot,'defaultFigurePosition', [0 0 420 560])
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/functions'

saveImages   = 0; %1 => stored images will be overwritten
saveTextFile = 0;

Ni = 1e2; %Number of iterations
Nt = [2 4 6 8]; %Numbers of terms
Nd = struct(...
    'Water', [1e2 1e3 1e4], ...
    'Park',  [1e2 5e2 1e3], ...
    'Urban', [1e2 1e3 5e3], ...
    'Test',  [1e2 1e3 1e4] ...
    ); %Number of data points


%Creating the RoIs
regions = struct(...
    'Water', zeros(2800, 2800), ...
    'Park',  zeros(2800, 2800), ...
    'Urban', zeros(2800, 2800), ...
    'Test',  zeros(2800, 2800) ...
    );
regions.Water( 130:850,    70:450) = 1; %[top:bottom, left:right];
regions.Park( 1900:2030,  970:1106) = 1; %[top:bottom, left:right];
regions.Urban(1400:1803, 1580:1783) = 1; %[top:bottom, left:right];
regions.Water( 130:850,    70:450) = 1; %[top:bottom, left:right];

channels = struct(...
    'Water', 1, ... %Blue
    'Park',  3, ... %Green
    'Urban', 2, ... %Red
    'Test',  3 ...
    );

%Separating out the relevant data field
data = struct(...
    'WL4', [], 'WL16', [],...
    'PL4', [], 'PL16', [],...
    'UL4', [], 'UL16', []...
    );

data.WL4  = auxfunction1(TL4,  regions.Water, channels.Water);
data.WL16 = auxfunction1(TL16, regions.Water, channels.Water);
data.PL4  = auxfunction1(TL4,  regions.Park,  channels.Park);
data.PL16 = auxfunction1(TL16, regions.Park,  channels.Park);
data.UL4  = auxfunction1(TL4,  regions.Urban, channels.Urban);
data.UL16 = auxfunction1(TL16, regions.Urban, channels.Urban);

%Making histograms for L=4,16
figure('Position', [0 500 1000 600]); hold on;

%For the subplot positions
lt = 0.07; %left
bt = 0.06; %bottom
wt = 0.27; %width
ht = 0.40; %height
hg = wt+0.05; %horizontal gap
vg = ht+0.09; %vertical gap

%Making the histogram figure
subplot('Position', [lt+0*hg, bt+1*vg, wt, ht])
auxfunction2(data.WL4, Nt, 8);
title '"Water" Region, $L = 4$'

subplot('Position', [lt+0*hg, bt+0*vg, wt, ht])
auxfunction2(data.WL16, Nt, 5);
title '"Water" Region, $L = 16$'

subplot('Position', [lt+1*hg, bt+1*vg, wt, ht])
auxfunction2(data.PL4, Nt, 8);
title '"Park" Region, $L = 4$'

subplot('Position', [lt+1*hg, bt+0*vg, wt, ht])
auxfunction2(data.PL16, Nt, 6);
title '"Park" Region, $L = 16$'

subplot('Position', [lt+2*hg, bt+1*vg, wt, ht])
auxfunction2(data.UL4, Nt, 7);
title '"Urban" Region, $L = 4$'

subplot('Position', [lt+2*hg, bt+0*vg, wt, ht])
auxfunction2(data.UL16, Nt, 7);
title '"Urban" Region, $L = 16$'

[~,supxlabel] = suplabel('$x$','x');
[~,supylabel] = suplabel('$\hat{f}(x)$','y');
supxlabel.Position(2) = supxlabel.Position(2) + 0.05;
supylabel.Position(1) = supylabel.Position(1) + 0.03;

%Saving the histograms and PDF estimates
if saveImages == 1
    tmp = getframe(gcf);
    imwrite(tmp.cdata, ['/Users/torgeirbrenn/Documents/', ...
        'Skole/Masteroppgave/Innlevering/figures/SARChangeHistograms.png'])
end

%Making distance plots
figure('Position', [0 500 1000 600]); hold on;

subplot('Position', [lt+0*hg, bt+1*vg, wt, ht])
auxfunction3(data.WL4, Ni, Nd.Water, Nt);
title '"Water" Region, $L = 4$'

subplot('Position', [lt+0*hg, bt+0*vg, wt, ht])
auxfunction3(data.WL16, Ni, Nd.Water, Nt);
title '"Water" Region, $L = 16$'

subplot('Position', [lt+1*hg, bt+1*vg, wt, ht])
auxfunction3(data.PL4, Ni, Nd.Park, Nt);
title '"Park" Region, $L = 4$'

subplot('Position', [lt+1*hg, bt+0*vg, wt, ht])
auxfunction3(data.PL16, Ni, Nd.Park, Nt);
title '"Park" Region, $L = 16$'

subplot('Position', [lt+2*hg, bt+1*vg, wt, ht])
auxfunction3(data.UL4, Ni, Nd.Urban, Nt);
title '"Urban" Region, $L = 4$'

subplot('Position', [lt+2*hg, bt+0*vg, wt, ht])
auxfunction3(data.UL16, Ni, Nd.Urban, Nt);
title '"Urban" Region, $L = 16$'

[~,supxlabel] = suplabel('$N$','x');
[~,supylabel] = suplabel('Distance $d_{\rm KL}$','y');
supxlabel.Position(2) = supxlabel.Position(2) + 0.05;
supylabel.Position(1) = supylabel.Position(1) + 0.03;

%Saving the figure
if saveImages == 1
    tmp = getframe(gcf);
    imwrite(tmp.cdata, ['/Users/torgeirbrenn/Documents/', ...
        'Skole/Masteroppgave/Innlevering/figures/SARChangeDistances.png'])
end

%%% CREATING A DATA TABLE, NOT USED IN THE THESIS %%%
if saveTextFile == 1
    %Setting up the containers
    results = struct(...
        'WL4',  zeros(numel(Nt),numel(Nd.Water)), ... %Water, L=4
        'WL16', zeros(numel(Nt),numel(Nd.Water)), ... %Water, L=16
        'PL4',  zeros(numel(Nt),numel(Nd.Park)), ... %Park, L=4 etc.
        'PL16', zeros(numel(Nt),numel(Nd.Park)), ...
        'UL4',  zeros(numel(Nt),numel(Nd.Urban)), ...
        'UL16', zeros(numel(Nt),numel(Nd.Urban)) ...
        );
    
    x = struct(...
        'WL4', [], 'WL16', [],...
        'PL4', [], 'PL16', [],...
        'UL4', [], 'UL16', []...
        );
    f = x;  %I.e. the same empty structure
    
    %Creating the target PDfs
    [x.WL4,  f.WL4]  = targetpdf('kernel', data.WL4);
    [x.WL16, f.WL16] = targetpdf('kernel', data.WL16);
    [x.PL4,  f.PL4]  = targetpdf('kernel', data.PL4);
    [x.PL16, f.PL16] = targetpdf('kernel', data.PL16);
    [x.UL4,  f.UL4]  = targetpdf('kernel', data.UL4);
    [x.UL16, f.UL16] = targetpdf('kernel', data.UL16);
    
    %Computing the results
    for nd = 1:numel(Nd.Water)
        results.WL4(:, nd)  = auxfunction4(data.WL4, x.WL4, f.WL4, ...
            Nd.Water(nd), Ni, Nt);
        results.WL16(:, nd) = auxfunction4(data.WL16,x.WL16,f.WL16,...
            Nd.Water(nd), Ni, Nt);
    end
    
    for nd = 1:numel(Nd.Park)
        results.PL4(:, nd)  = auxfunction4(data.PL4, x.PL4, f.PL4, ...
            Nd.Park(nd), Ni, Nt);
        results.PL16(:, nd) = auxfunction4(data.PL16,x.PL16,f.PL16,...
            Nd.Park(nd), Ni, Nt);
    end
    
    for nd = 1:numel(Nd.Urban)
        results.UL4(:, nd)  = auxfunction4(data.UL4, x.UL4, f.UL4, ...
            Nd.Urban(nd), Ni, Nt);
        results.UL16(:, nd) = auxfunction4(data.UL16,x.UL16,f.UL16,...
            Nd.Urban(nd), Ni, Nt);
    end
    numcols = numel(Nd.Water) + numel(Nd.Park) + numel(Nd.Urban);
    
    %Printing the results in the table format, to a text file
    printArrayL4  = [results.WL4,  results.PL4,  results.UL4];
    printArrayL16 = [results.WL16, results.PL16, results.UL16];
    formatSpec = ['$N=%i', '$', repmat(' & $%3.2e}$', 1, numcols), ...
        '\\\\ \\hline \n'];
    
    fileIDL4  = fopen('../results/SARChangeDataL4.txt', 'w');
    fileIDL16 = fopen('../results/SARChangeDataL16.txt','w');
    %fileIDL4 = fopen('../results/SARChangeDataTEST.txt','w');
    
    for nt = 1:numel(Nt)
        fprintf(fileIDL4,  formatSpec, [Nt(nt), printArrayL4(nt,:)]);
        fprintf(fileIDL16, formatSpec, [Nt(nt), printArrayL16(nt,:)]);
    end
    
    %Instructions:
    %Replace 'e-0' with '{\cdot}10^{-', and if necessary
    %replace 'e+0' with '{\cdot}10^{' and pad with \hspace{5pt}
end

%%% AUXILLARY FUNCTIONS %%
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/'
%auxfunction1: creates change data from region 'RoI' of data 'T', channel
%'ch'
function qdat = auxfunction1(T, RoI, ch)
tmp = real(T(:,:,ch,ch)); %Discard the zero complex parts
dat = tmp(RoI==1); %All intensities within the RoI, Nx1

%Splitting dataset into two
r = round(rand(size(dat)));
d1 = dat(r==0);
d2 = dat(r==1);

%Creating change data
m = min(numel(d1), numel(d2));
%Randomly permutate the first m data points, creating quotients.
qdat = d1(randperm(m))./d2(randperm(m));
end

%auxfunction2: plots the histogram and MKBK series expansion with
%N=4,6,8 for coherency data T, region RoI and channel ch. Truncates x at
%trunc.
function auxfunction2(dat, Nt, trunc)
%x grid
if nargin < 3
    trunc = max(dat);
end
dx = trunc/1e4;
x = dx:dx:trunc;

PDFest = zeros(numel(x), numel(Nt)); %PDF estimates
for nt = 1:numel(Nt)
    PDFest(:, nt) = mkbkfitequalshapes(x, dat, Nt(nt), 1, 1);
end

%Linestyle
style = {'-', '--', '-', '--'};

%Legend
lgd = strings(numel(Nt)+1, 1);
lgd(1) = 'Histogram';
lgd(2) = 'Kernel';
for nt = 2:numel(Nt)
    lgd(nt+1) = ['$N=', num2str(Nt(nt)), '$'];
end

%Color
firstcolor = [0 0 1]; %Blue
lastcolor  = [1 0 0]; %Red
color = zeros(3, numel(Nt));
color(1,:) = linspace(firstcolor(1), lastcolor(1), numel(Nt));
color(2,:) = linspace(firstcolor(2), lastcolor(2), numel(Nt));
color(3,:) = linspace(firstcolor(3), lastcolor(3), numel(Nt));

%Plotting
hold on; grid on; box on
hh = histogram(dat, 'Normalization', 'pdf', 'EdgeAlpha', 0, ...
    'FaceColor', [0.494 0.184 0.556], 'FaceAlpha', 0.4);
hh.FaceColor = [0.4660 0.6740 0.1880]; %Green color alternative
hh.FaceAlpha = 0.5;
hh.BinWidth = trunc/100;
for nt = 1:numel(Nt)
    plot(x, PDFest(:, nt), style{nt}, 'Color', color(:,nt), ...
        'LineWidth', 1.2)
end

lh = legend(lgd);
set(lh,'Interpreter','latex');

%Cosmetics
ax = gca;
ax.XLim(1) = 0;
ax.XLim(2) = x(end);
ax.XLim(2) = trunc;
ax.YLim(1) = 0;

end

%auxfunction3: plots mean of KL distance of Ni iterations of the MKBK
%series with Nt(1:end) terms, using Nd(1:end) data points drawn from dat.
function auxfunction3(dat, Ni, Nd, Nt)

%Generating the target PDF using kernel density estimation
[x, f]  = targetpdf('kernel', dat);

results = zeros(numel(Nt), numel(Nd)); %results(:,1) empty by design

%Computing the results
for nd = 1:numel(Nd)
    results(:, nd) = auxfunction4(dat, x, f, Nd(nd), Ni, Nt);
end

%Linestyle
switch numel(Nd)
    case 3
        style = {':+', '--p', '-*'};
    case 4
        style = {':+', '-.d', '--p', '-*'};
    case 5
        style = {':+', ':s', '-.d', '--p', '-*'};
end

%Plotting
hold on; 
for nd = 1:numel(Nd)
    plot(Nt, results(:,nd), style{nd}, ...
        'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.2)
end

%Legend
lgd = strings(numel(Nd), 1);
for nd = 1:numel(Nd)
    lgd(nd) = num2bank(Nd(nd));
end
lh = legend(lgd, 'Location', 'northwest');
set(lh,'Interpreter','latex');
lh.Title.String = 'Data points';
lh.Title.FontWeight = 'normal';


%Cosmetics
axis tight; grid on; box on;
ax = gca; ax.XTick = Nt; ax.XLim = [Nt(1), Nt(end)];
auxfunction5();
end

%auxfunction4: the mkbk series with Nt = 2, 4, 6, 8 applied to 'Nd' data
%points drawn from 'dat', returns distance to target PDF 'tgt', on grid
%'x', mean of 'Ni' iterations
function res = auxfunction4(dat, x, tgt, Nd, Ni, Nt)

%Setting up a container for the results
distances = zeros(Ni,numel(Nt));

dx = x(2)-x(1); %Need this to standardize the distance

wb = waitbar(0, 'Computing');
for ni = 1:Ni
    %Draw Nd random sample without replacement
    datasubset = datasample(dat, Nd, 'Replace', false);
    
    %Computing the distances
    for nt = 1:numel(Nt)
        distances(ni,nt) = kldist(tgt, ...
            mkbkfitequalshapes(x, datasubset, Nt(nt), 1, 1), dx);
    end
    waitbar(ni/Ni, wb)
end
close(wb)
res = mean(distances,1);
end

%auxfunction5: Controlling the appearance of the logarithmic Y axis
function auxfunction5(ax, ticks)
if nargin < 2
    ticks = -10:1:2; %I.e. more than enough
end

if nargin == 0
    ax = gca;
end
ax.YScale = 'log';
ax.YTick = 10.^(ticks);
ax.YTickLabel = cell(numel(ticks),1);
for i = 1:numel(ticks)
    ax.YTickLabel{i} = ['10^{', num2str(ticks(i)), '}'];
end
end