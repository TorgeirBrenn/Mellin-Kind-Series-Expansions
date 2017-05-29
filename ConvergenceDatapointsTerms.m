%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we analyze how the MK series expansions converge as the number of
% correcting terms and data points are increased.
%
% Made by Torgeir Brenn, 2017. For the custom functions, the author is
% specified within those functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear variables
set(groot, 'defaultFigureColor', [1 1 1]);   %White figure background
set(groot,'defaulttextinterpreter','latex'); %Figure text
set(0,'defaultAxesFontSize',16);
set(groot,'defaultFigurePosition', [0 500 1500 400]) %Default: 560 420
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/functions'

%Declaring variables
Ni = 1e1;   %No. of iterations
Nd = [1e2 1e2 1e2];   %No. of data points
Nt = 8;     %Highest order (log-)moment/cumulant used. N<=8

parameters = struct(...
        'm', 10, ...
        'g', 10, ...
        'L', 16, ...
        'M', -10, ...
        'nu', 0.5 ...
);
distribution = 'g0'; %Supported: 'k', 'gamma', 'g0', 'ggd', 'igam' 

%Spacing for the figure
lt = 0.05; %left
bt = 0.12; %bottom
wt = 0.21; %width
ht = 0.75; %height
hg = wt+0.03; %horizontal gap


figure;

subplot('Position', [lt+0*hg, bt, wt, ht])
auxfunction1(distribution, parameters, 'mkgk', Ni, Nd, Nt)
title('MKGK Series', 'FontSize', 18)

subplot('Position', [lt+1*hg, bt, wt, ht])
auxfunction1(distribution, parameters, 'mklk', Ni, Nd, Nt)
title('MKLK Series', 'FontSize', 18)

subplot('Position', [lt+2*hg, bt, wt, ht])
auxfunction1(distribution, parameters, 'mke', Ni, Nd, Nt)
title('MKE Series', 'FontSize', 18)

subplot('Position', [lt+3*hg, bt, wt, ht])
auxfunction1(distribution, parameters, 'mkbk', Ni, Nd, Nt)
title('MKBK Series', 'FontSize', 18)

%Title for the entire figure, based on the MathWorks Central function
%suplabel
titletext = pdfname(distribution, parameters, 0);
[~,suptitle] = suplabel(titletext, 't');
suptitle.Position(2) = suptitle.Position(2) + 0.05;
suptitle.FontSize = 20;

[~,supxlabel] = suplabel('Highest order log-cumulant corrected for','x');
[~,supylabel] = suplabel('Distance $d_{\rm KL}$','y');
supxlabel.Position(2) = supxlabel.Position(2) + 0.04;
supxlabel.FontSize = 18;
supylabel.Position(1) = supylabel.Position(1) + 0.035;
supylabel.FontSize = 18;

%NOTE: There are some problems with the suplabel font sizes, use the code
%below to set it in retrospect
% fh = gcf;
% fh.Children(9).YLabel.FontSize = 18;
% fh.Children(10).XLabel.FontSize = 18;
% fh.Children(11).Title.FontSize = 20;

%%% AUXILLARY FUNCTIONS %%
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/'
function auxfunction1(dist, param, MKseries, Ni, Nd, Nt)
[x, f] = targetpdf(dist, param);
dx = x(2) - x(1);

%Waitbar
wb = waitbar(0,'Computing');

%Preparing the holders of the estimates and measures
PDFest = zeros(size(f)); %PDF estimate, temporary holder
R    = zeros(3, Ni, numel(Nd), Nt); %Distance, divergence measures

for nd = 1:numel(Nd)
    waitbar(0, wb, ['Computing ', num2str(nd), ' of ', num2str(numel(Nd))]);
    for ni = 1:Ni
        %Generating data
        data = generatedata(Nd(nd), dist, param);
        
        for nt = 2:Nt
            switch MKseries
                case "mkgk" 
                    PDFest = mkgkfit(x, data, nt);
                case "mklk"
                    PDFest = mklkfit(x, data, nt);
                case "mke"
                    PDFest = mkefit(x, data, nt);
                case "mkbk"
                    PDFest = mkbkfit(x, data, nt);
            end
           
            R(2, ni, nd, nt) = kldist(f, PDFest, dx);
        end
        waitbar(ni/Ni, wb);
    end
end

close(wb);
disp 'Computation complete'


%Color
switch MKseries
    case "mkgk"
        color = [0 0.4470 0.7410];
    case "mklk"
        color = [0.8500 0.3250 0.0980];
    case "mke"
        color = [0 0 0];
    case "mkbk"
        color = [0.4660 0.6740 0.1880];
end

%Line styles
switch numel(Nd)
    case 3
        style = {':+', '--p', '-*'};
    case 4
        style = {':+', '-.d', '--p', '-*'};
    case 5
        style = {':+', ':s', '-.d', '--p', '-*'};
end

%Legend
lgd = strings(numel(Nd), 1);
for nd = 1:numel(Nd)
    lgd(nd) = num2bank(Nd(nd));
end

%To save time
results = squeeze(mean(real(R(2,:,:,:)), 2, 'omitnan'));

%Plotting
for nd = 1:numel(Nd)
    plot(2:Nt, results(nd,2:Nt), ...
        style{nd}, 'Color', color, 'LineWidth', 1.2);
    hold on;
end
%title('Kullback-Leibler', 'FontSize', 18)
%xlabel 'Highest order (log-)cumulant corrected for'
%ylabel 'Distance'
lh = legend (lgd, 'Location', 'northwest');
lh.Title.String = 'Data points';
lh.Title.FontWeight = 'normal';
set(lh,'Interpreter','latex');
%axis([2 Nt 0 0.01]); grid on;
axis tight; grid on;
ax = gca; ax.XTick = 2:Nt; ax.YLim(1) = 0;
ax.YScale = 'log';
end