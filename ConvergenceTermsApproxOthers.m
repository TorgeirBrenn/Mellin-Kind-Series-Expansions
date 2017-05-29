%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots the cases that were not made with ConvergengeTermsApprox.m
%
% Made by Torgeir Brenn, 2017. For the custom functions, the author is
% specified within those functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear variables
set(groot, 'defaultFigureColor', [1 1 1]);   %White figure background
set(groot,'defaulttextinterpreter','latex'); %Figure text
set(0,'defaultAxesFontSize',12);
set(groot,'defaultFigurePosition', [0 500 1000 1500])
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/functions'

saveImages = 0;

%Waitbar
wb = waitbar(0,'Computing');

dist =  cell(16,1);
param = cell(16,1);
pos =  zeros(16,4);

dist{1}  = "gamma"; param{1}  = struct('L', 4, 'm', 1);
dist{2}  = "gamma"; param{2}  = struct('L', 4, 'm', 10);
dist{3}  = "gamma"; param{3}  = struct('L', 16, 'm', 1);
dist{4}  = "gamma"; param{4}  = struct('L', 16, 'm', 10);
dist{5}  = "igam";  param{5}  = struct('L', 4, 'm', 1);
dist{6}  = "igam";  param{6}  = struct('L', 4, 'm', 10);
dist{7}  = "igam";  param{7}  = struct('L', 16, 'm', 1);
dist{8}  = "ggd";   param{8}  = struct('L', 4, 'm', 10, 'nu', 2);
dist{9}  = "ggd";   param{9}  = struct('L', 16, 'm', 10, 'nu', 0.5);
dist{10} = "ggd";   param{10} = struct('L', 16, 'm', 10, 'nu', 2);
dist{11} = "k";     param{11} = struct('L', 4, 'm', 10, 'M', 1);
dist{12} = "k";     param{12} = struct('L', 4, 'm', 10, 'M', 10);
dist{13} = "k";     param{13} = struct('L', 16, 'm', 10, 'M', 1);
dist{14} = "g0";    param{14} = struct('L', 4, 'g', 1, 'M', -10);
dist{15} = "g0";    param{15} = struct('L', 4, 'g', 10, 'M', -10);
dist{16} = "g0";    param{16} = struct('L', 16, 'g', 1, 'M', -10);

lt = 0.05; %left
bt = 0.06; %bottom
wt = 0.21; %width
ht = 0.19; %height
hg = wt+0.03; %horizontal gap
vg = ht+0.05; %vertical gap
pos(1,:)   = [lt+0*hg, bt+3*vg, wt, ht];
pos(2,:)   = [lt+1*hg, bt+3*vg, wt, ht];
pos(3,:)   = [lt+2*hg, bt+3*vg, wt, ht];
pos(4,:)   = [lt+3*hg, bt+3*vg, wt, ht];
pos(5,:)   = [lt+0*hg, bt+2*vg, wt, ht];
pos(6,:)   = [lt+1*hg, bt+2*vg, wt, ht];
pos(7,:)   = [lt+2*hg, bt+2*vg, wt, ht];
pos(8,:)   = [lt+3*hg, bt+2*vg, wt, ht];
pos(9,:)   = [lt+0*hg, bt+1*vg, wt, ht];
pos(10,:)  = [lt+1*hg, bt+1*vg, wt, ht];
pos(11,:)  = [lt+2*hg, bt+1*vg, wt, ht];
pos(12,:)  = [lt+3*hg, bt+1*vg, wt, ht];
pos(13,:)  = [lt+0*hg, bt+0*vg, wt, ht];
pos(14,:)  = [lt+1*hg, bt+0*vg, wt, ht];
pos(15,:)  = [lt+2*hg, bt+0*vg, wt, ht];
pos(16,:)  = [lt+3*hg, bt+0*vg, wt, ht];

%For testing
% dist{1}  = "igam"; param{1}  = [16,10];
% dist{2}  = "ggd"; param{2}  = [4,10,0.5];
% dist{3}  = "k"; param{3}  = [16,10,10];
% dist{4}  = "g0"; param{4}  = [16,10,-10];

for i = 1:16
    subplot('Position', pos(i,:));
    auxfunction1(dist{i}, param{i}, 8);
    waitbar(i/16, wb);
end
close(wb);
disp 'Computation complete'

%Title for the entire figure, based on the MathWorks Central function
%suplabel
%[~,suptitle]  = suplabel('Convergence','t');
[~,supxlabel] = suplabel('Highest order log-cumulant corrected for','x');
[~,supylabel] = suplabel('Distance $d_{\rm KL}$','y');
supxlabel.Position(2) = supxlabel.Position(2) + 0.04;
supxlabel.FontSize = 22;
supylabel.Position(1) = supylabel.Position(1) + 0.04;
supylabel.FontSize = 22;

%Saving the image
if saveImages == 1
    tmp = getframe(gcf);
    imwrite(tmp.cdata, ['/Users/torgeirbrenn/Documents/', ...
        'Skole/Masteroppgave/Innlevering/figures/',...
        'ConvergenceApproxOthers.png'])
end

%%% AUXILLARY FUNCTIONS %%
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/'
%auxfunction1 : performs the actual plotting
function auxfunction1(dist, param, Nt)

%True distribution
[x, f, kappa] = targetpdf(dist, param, Nt);
dx = x(2)-x(1);
%Preparing the holders of the approximations
PDFapprox = zeros(numel(x), 4, Nt); %PDF estimates

%Determining which methods to plot
switch dist
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

%Note: PDFests(:, :, 1) is left empty by design, as it would be
%identical to PDFests(:, :, 2). This is done for the sake of
%clarity in terms of consistent notation. Change if necessary.

for nt=2:Nt
    PDFapprox(:, 1, nt) = mkgkfitexact(x, kappa, nt);
    PDFapprox(:, 2, nt) = mklkfitexact(x, kappa, nt);
    PDFapprox(:, 3, nt) = mkefitexact(x, kappa, nt);
    PDFapprox(:, 4, nt) = mkbkfitexact(x, kappa, nt);
end

%The struct "methods" defines the style.
methods = struct(...
    'name', {'MKGK'; 'MKLK';'MKE';'MKBK'}, ...
    'style', {'--d'; '--+'; '-.+'; '--o'}, ...
    'color', {[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; ...
    'k'; [0.4660 0.6740 0.1880]} ...
    );

results = zeros(4, Nt);

%Computing the distances
for m = p
    for nt = 2:Nt
        results(m,nt) = kldist(f, PDFapprox(:,m,nt), dx);
    end
end

%Potting
hold on;
for i = p
    plot(2:Nt, results(i,2:Nt), ...
        methods(i).style, 'Color', methods(i).color, 'LineWidth', 1.2);
end
%lh = legend (methods(p).name, 'Location', 'northeast');
%set(lh,'Interpreter','latex')
axis tight; grid on; box on
auxfunction2();

titletext = pdfname(dist, param, 1);
title(titletext)

end

%auxfunction2: Controlling the appearance of the logarithmic Y axis
function auxfunction2(ax, ticks)
if nargin < 2
    ticks = -10:1:2; %I.e. more than enough
end

if nargin == 0
    ax = gca;
end
ax.XLim = [ax.Children(1).XData(1) ax.Children(1).XData(end)];
ax.XAxis.TickValues = 2:8; %Limited by 'XLim', i.e. works for Nt<8.
ax.YScale = 'log';
ax.YTickLabelRotation = 90;
ax.TitleFontSizeMultiplier = 1.2;
ax.YTick = 10.^(ticks);
ax.YTickLabel = cell(numel(ticks),1);
for i = 1:numel(ticks)
    ax.YTickLabel{i} = ['10^{', num2str(ticks(i)), '}'];
end
end