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
set(0,'defaultAxesFontSize',16);
set(groot,'defaultFigurePosition', [0 500 2*560 800])
%set(groot,'defaultFigurePosition', [0 500 560 420])
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/functions'

saveImages = 0;

Nt = [6 8];

distribution = ["k", "k"];

parameters = struct(... %Add g, nu if necessary!
        'm',  { 10  10}, ...
        'L',  { 16  16}, ...
        'M',  {  1  10} ...
);

%Spacing for the figure
lt = 0.07; %left
bt = 0.06; %bottom
wt = 0.43; %width
ht = 0.81; %height
hg = wt+0.06; %horizontal gap

figure;
auxfunction1(distribution(1), parameters(1),  [lt+0*hg, bt, wt, ht], Nt(1))
auxfunction1(distribution(2), parameters(2),  [lt+1*hg, bt, wt, ht], Nt(2))

%Saving the figure
if saveImages == 1
    tmp = getframe(gcf);
    imwrite(tmp.cdata, ['/Users/torgeirbrenn/Documents/', ...
        'Skole/Masteroppgave/Innlevering/figures/',...
        'ConvergenceVisualized.png'])
end

%%% AUXILIARY FUNCTIONS %%
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/'
%auxfunction1: Plots MKGK series convergence to 'dist' with 'param', PDF
%with approximations (top), approximation errors (bottom)
function auxfunction1(dist, param, pos, Nt)
%Target distribution
[x, f, kappa] = targetpdf(dist, param, Nt);

%Preparing the holders of the approximations
PDFapprox = zeros(numel(x), Nt); %PDF estimates

%Note: PDapprox(:, 1) is left empty by design, as it would be
%identical to PDFests(:, :, 2). This is done for the sake of
%clarity in terms of consistent notation. Change if necessary.

for nt=2:Nt
    PDFapprox(:, nt) = mkgkfitexact(x, kappa, nt);
end

style = {'--', '-', '--', '-', '--', '-', '--', '-'};

lgd = strings(Nt,1);
if dist == 'k'
    lgd(1) = 'Target $K$';
end
for nt = 2:Nt
    lgd(nt) = ['$N=', num2str(nt), '$'];
end

firstcolor = [0 0 1]; %Blue
lastcolor = [1 0 0];  %Red
color = zeros(3, 8);
color(1,2:8) = linspace(firstcolor(1), lastcolor(1), 7);
color(2,2:8) = linspace(firstcolor(2), lastcolor(2), 7);
color(3,2:8) = linspace(firstcolor(3), lastcolor(3), 7);

pos(4) = pos(4)/2; %Halve the height, to get two plots

%All methods plotted
subplot('Position', pos + [0 0.49 0 0]); hold on;
plot(x, f, 'k', 'LineWidth', 0.8);
for nt = 2:Nt
    plot(x, PDFapprox(:,nt), style{nt}, 'Color', color(:,nt),...
        'LineWidth', 0.8);
end
axis tight; grid on; box on;
xlabel '$x$'
ylabel '$f(x),\,\hat{f}(x)$'
lh = legend(lgd);
set(lh,'Interpreter','latex')

title(pdfname(dist, param), 'FontSize', 18)


subplot('Position', pos); hold on;
for nt = 2:Nt
    plot(x, PDFapprox(:,nt)-f, style{nt}, 'Color', color(:,nt),...
        'LineWidth', 0.8);
end
plot(x, zeros(size(x)), 'k', 'LineWidth', 0.8);
title('Approximation Error', 'FontSize', 18)
xlabel '$x$'
ylabel '$f(x)-\hat{f}(x)$'
axis tight; grid on; box on;
lh = legend(lgd(2:Nt));
set(lh,'Interpreter','latex')
ax = gca;
ax.YAxis.Limits = [-max(abs(ax.YAxis.Limits)), max(abs(ax.YAxis.Limits))];
end