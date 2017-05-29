%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basically the same as "ConvergenceTerms.m" but here the approximated
% distributions are known
%
% Made by Torgeir Brenn, 2017. For the custom functions, the author is
% specified within those functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear variables
set(groot, 'defaultFigureColor', [1 1 1]);   %White figure background
set(groot,'defaulttextinterpreter','latex'); %Figure text
set(0,'defaultAxesFontSize',16);
set(groot,'defaultFigurePosition', [0 500 1200 420])
%set(groot,'defaultFigurePosition', [0 500 560 840])
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/functions'

saveImages = 1;

%Declaring variables
Nt = 8;     %Highest order (log-)moment/cumulant used. N<=8
parameters = struct(...
    'm', 10, ...
    'g', 10, ...
    'L', 16, ...
    'M', -10, ...
    'nu', 0.5 ...
    );
distribution = 'g0'; %Currently supported: 'k', 'gamma', 'g0', 'ggd', 'igam'

%Waitbar
wb = waitbar(0,'Computing');

[x, f, kappa] = targetpdf(distribution, parameters, Nt(end));
dx = x(2)-x(1);

%Preparing the holders of the approximations
PDFapprox = zeros(numel(x), 4, Nt); %PDF estimates
results = struct(...
    'Bhatta', zeros(4, Nt), ...
    'KL',     zeros(4, Nt), ...
    'Max',    zeros(4, Nt)...
    );


for nt=2:Nt
    
    PDFapprox = mkgkfitexact(x, kappa, nt);
    results.Bhatta(1, nt) = bhattadist(f, PDFapprox, dx);
    results.KL(1, nt) = kldist(f, PDFapprox, dx);
    results.Max(1, nt) = max(abs(f-PDFapprox));
    
    PDFapprox = mklkfitexact(x, kappa, nt);
    results.Bhatta(2, nt) = bhattadist(f, PDFapprox, dx);
    results.KL(2, nt) = kldist(f, PDFapprox, dx);
    results.Max(2, nt) = max(abs(f-PDFapprox));
    
    PDFapprox = mkefitexact(x, kappa, nt);
    results.Bhatta(3, nt) = bhattadist(f, PDFapprox, dx);
    results.KL(3, nt) = kldist(f, PDFapprox, dx);
    results.Max(3, nt) = max(abs(f-PDFapprox));
    
    PDFapprox = mkbkfitexact(x, kappa, nt);
    results.Bhatta(4, nt) = bhattadist(f, PDFapprox, dx);
    results.KL(4, nt) = kldist(f, PDFapprox, dx);
    results.Max(4, nt) = max(abs(f-PDFapprox));
    
    waitbar(nt/Nt, wb);
end

close(wb);
disp 'Computation complete'

%The struct "methods" defines the style.
methods = struct(...
    'name', {'MKGK'; 'MKLK';'MKE';'MKBK'}, ...
    'style', {'--+'; '-.+'; '--d'; '--o'}, ...
    'color', {[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; 'k';...
    [0.4660 0.6740 0.1880]} ...
    );


%Deciding which methods to plot
switch distribution
    case "gamma"
        p = [2 3];
    case "igam"
        p = 1:3;
    case "k"
        p = 1:4;
    case "ggd"
        p = 1:4;
    case "g0"
        p = 1:3;
end

%Subplot positioning
lt = 0.06; %left
bt = 0.12; %bottom
wt = 0.27; %width
ht = 0.76; %height
hg = wt+0.04; %horizontal gap

%Left: Bhattacharyya distance
subplot('Position', [lt+0*hg bt wt ht]); hold on; grid on;
for i = p
    plot(2:Nt, results.Bhatta(i,2:Nt), ...
        methods(i).style, 'Color', methods(i).color, 'LineWidth', 1.2);
end
title('Bhattacharyya', 'FontSize', 18)

lh = legend (methods(p).name, 'Location', 'northeast');
set(lh,'Interpreter','latex');
auxfunction1();


%Center: Kullback-Leibler distance
subplot('Position', [lt+1*hg bt wt ht]); hold on; grid on;
for i = p
    plot(2:Nt, results.KL(i,2:Nt), ...
        methods(i).style, 'Color', methods(i).color, 'LineWidth', 1.2);
end
title('Kullback-Leibler', 'FontSize', 18)
%lh = legend (methods(p).name, 'Location', 'northeast');
set(lh,'Interpreter','latex');
auxfunction1();

%Right: Maximum distance
subplot('Position', [lt+2*hg bt wt ht]); hold on; grid on;
for i = p
    plot(2:Nt, results.Max(i,2:Nt), ...
        methods(i).style, 'Color', methods(i).color, 'LineWidth', 1.2);
end
title('Maximum', 'FontSize', 18)
%lh = legend (methods(p).name, 'Location', 'northeast');
set(lh,'Interpreter','latex');
auxfunction1();

%Title for the entire figure, based on the MathWorks Central function
%suplabel
titletext = pdfname(distribution, parameters, 0);
[~,suptitle] = suplabel(titletext, 't');
suptitle.Position(2) = suptitle.Position(2) + 0.05;
suptitle.FontSize = 22;

[~,supylabel] = suplabel('Distance','y');
supylabel.Position(1) = supylabel.Position(1) + 0.03;
supylabel.FontSize = 18;
[~,supxlabel] = suplabel('Highest order log-cumulant corrected for','x');
supxlabel.Position(2) = supxlabel.Position(2) + 0.04;
supxlabel.FontSize = 18;

%Saving the image
if saveImages == 1
    tmp = getframe(gcf);
    imwrite(tmp.cdata, ['/Users/torgeirbrenn/Documents/', ...
        'Skole/Masteroppgave/Innlevering/figures/',...
        'ConvergenceApprox', distribution, '.png'])
end

%%% AUXILLARY FUNCTIONS %%
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/'

%auxfunction1(): Controlling the appearance of the logarithmic Y axis
function auxfunction1(ax, ticks)
if nargin < 2
    ticks = -10:1:2; %I.e. more than enough
end

if nargin == 0
    ax = gca;
end
box on;

ax.XLim = [ax.Children(1).XData(1) ax.Children(1).XData(end)];
%ax.XTicks = ax.Children(1).XData; %2:Nt

ax.YScale = 'log';
ax.YTick = 10.^(ticks);
ax.YTickLabel = cell(numel(ticks),1);
for i = 1:numel(ticks)
    ax.YTickLabel{i} = ['10^{', num2str(ticks(i)), '}'];
end
ax.YLabel.Position(1) = ax.YLabel.Position(1)-0.06;
end