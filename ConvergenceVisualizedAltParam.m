%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we analyze how non-trivial choices of the kernel parameters affect
% the mellin kind series expansions.
%
% Made by Torgeir Brenn, 2017. For the custom functions, the author is
% specified within those functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear variables
set(groot, 'defaultFigureColor', [1 1 1]);   %White figure background
set(groot,'defaulttextinterpreter','latex'); %Figure text
set(0,'defaultAxesFontSize',18);
set(groot,'defaultFigurePosition', [0 500 1500 700])
%set(groot,'defaultFigurePosition', [0 500 560 420])
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/functions'

saveImages = 1;

Nt = 8;

distribution = ["k", "k", "k", "k"];

parameters = struct(... %Add g, nu if necessary!
        'm',  { 10  10  10  10}, ...
        'L',  {  4   4   4   4}, ...
        'M',  {  1   5  10  20} ...
);

%Spacing for the figure
lt = 0.06; %left
bt = 0.1; %bottom
wt = 0.185; %width
ht = 0.75; %height
hg = wt+0.065; %horizontal gap

figure;
for i = 1:4
    auxfunction1(distribution(i), parameters(i), ...
        [lt+(i-1)*hg, bt, wt, ht], Nt);
end

%Placing the legend along the bottom of the figure
fh = gcf;
lh = fh.Children(2);
lh.Visible = 'on';
lh.Orientation = 'horizontal';
lh.Position = [0 0 1 0.035];

[~,supxlabel] = suplabel('$x$','x');
supxlabel.Position(2) = supxlabel.Position(2) + 0.045;
supxlabel.FontSize = 18;

%Saving the figure
if saveImages == 1
    tmp = getframe(gcf);
    imwrite(tmp.cdata, ['/Users/torgeirbrenn/Documents/', ...
        'Skole/Masteroppgave/Innlevering/figures/',...
        'AlternativeParametersConvergenceVisualized.png'])
end

%%% AUXILIARY FUNCTIONS %%
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/'
%auxfunction1: Plots MKGK series convergence to 'dist' with 'param', PDF
%with approximations (top), approximation errors (bottom)
function auxfunction1(dist, param, pos, Nt)
%Target distribution
[x, f, kappa] = targetpdf(dist, param, Nt);

%Preparing the holders of the approximations
PDFapprox = zeros(numel(x), Nt+1); %PDF estimates

for nt=0:Nt
    PDFapprox(:, nt+1)= mkgkfitgivenparam(x, kappa, nt, param.L, param.m);
end

%Tailored kernel
tk = mkgkfitexact(x, kappa, 2);

%The line below presents the tailored value (MoLC estimate) of L
%fsolve(@(L) psi(1,L)-kappa(2),4, optimoptions('fsolve', 'Display', 'off'))

%Style
style = {'-', '--', '-', '--', '-', '--', '-', '--', '-', '--'};

%Legend
lgd = strings(Nt+3,1);
if dist == 'k'
    lgd(1) = 'Target $K$';
end
lgd(2) = 'Tailored Kernel';
for nt = 0:Nt
    lgd(nt+3) = ['$N=', num2str(nt), '$'];
end

firstcolor = [0 0 1]; %Blue
lastcolor  = [1 0 0]; %Red
color = zeros(3, 9);
color(1,:) = linspace(firstcolor(1), lastcolor(1), 9);
color(2,:) = linspace(firstcolor(2), lastcolor(2), 9);
color(3,:) = linspace(firstcolor(3), lastcolor(3), 9);

pos = pos .* [1 1 1 0.5]; %Halve the heigth, to get two plots

%All methods plotted
subplot('Position', pos + [0 0.49 0 0]); hold on;
plot(x, f, 'k', 'LineWidth', 0.8);
plot(x, tk, 'g', 'LineWidth', 0.8);
for nt = 0:Nt
    plot(x, PDFapprox(:,nt+1), style{nt+1}, 'Color', color(:,nt+1),...
        'LineWidth', 0.8);
end
ylabel '$f(x),\,\hat{f}(x)$'
lh = legend(lgd); 
set(lh,'Interpreter','latex')
lh.Visible = 'off'; %I.e. we return an invisible legend.
axis tight; grid on; box on;

title(pdfname(dist, param), 'FontSize', 18);

%Approximation errors
subplot('Position', pos); hold on;
plot(x, tk-f, 'g', 'LineWidth', 0.8);
for nt = 0:Nt
    plot(x, PDFapprox(:,nt+1)-f, style{nt+1}, 'Color', color(:,nt+1),...
        'LineWidth', 0.8);
end
plot(x, zeros(size(x)), 'k', 'LineWidth', 0.8);
title('Approximation Error', 'FontSize', 18)
ylabel '$\hat{f}(x)-f(x)$'
axis tight; grid on; box on;
%lh = legend(lgd(2:(N+2)));
%set(lh,'Interpreter','latex')
ax = gca;
ax.YAxis.Limits = [-max(abs(ax.YAxis.Limits)), max(abs(ax.YAxis.Limits))];
end