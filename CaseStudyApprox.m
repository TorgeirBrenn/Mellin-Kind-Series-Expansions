%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basically the same as "ComparingMethods.m" but here the approximated
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
set(groot,'defaultFigurePosition', [250 0 1000 1000])
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/functions'

saveImages    = 0;
writeTextFile = 0; %1: Will overwrite the given text file with new results.

%Declaring variables
Nt = 4;     %Highest order (log-)moment/cumulant used. N<=8
parameters = struct(...
    'm', 10, ...
    'g', 2, ...
    'L', 4, ...
    'M', 10, ...
    'nu', 0.5 ...
    );
distribution = 'ggd'; %Supported: 'k', 'gamma', 'g0', 'ggd', 'igam'

%True distribution
[x, f, kappa] = targetpdf(distribution, parameters, Nt);
dx = x(2)-x(1);

%Preparing the holders of the estimates and measures
PDFapprox = zeros(numel(x), 9); %PDF estimates

PDFapprox(:, 1) = molckfitexact(x, kappa);
PDFapprox(:, 2) = molcggdfitexact(x, kappa);
PDFapprox(:, 3) = mklkfitexact(x, kappa, 2);
PDFapprox(:, 4) = mklkfitexact(x, kappa, Nt);
PDFapprox(:, 5) = mkefitexact(x, kappa, Nt);
PDFapprox(:, 6) = mkgkfitexact(x, kappa, 2);
PDFapprox(:, 7) = mkgkfitexact(x, kappa, Nt);
PDFapprox(:, 8) = mkbkfitexact(x, kappa, 2);
PDFapprox(:, 9) = mkbkfitexact(x, kappa, Nt);

disp 'Computation complete'

%Name, style, color
methods = struct(...
    'name', {'MoLC K'; 'MoLC G{\Gamma}D'; ...
    'MKLK kernel'; 'MKLK series'; 'MKE series'; 'MKGK kernel'; ...
    'MKGK series'; 'MKBK kernel'; 'MKBK series'}, ...
    'style', {''; ''; ''; '--'; '-.'; ''; '--'; ''; '--'}, ...
    'color', {[0.4940 0.1840 0.5560]; [0.9290 0.6940 0.1250]; ...
    [0.8500 0.3250 0.0980]; [0.8500 0.3250 0.0980]; ...
    [0.8500 0.3250 0.0980]; [0 0.4470 0.7410]; [0 0.4470 0.7410]; ...
    [0.4660 0.6740 0.1880]; [0.4660 0.6740 0.1880]} ...
    );
%Style guide:
%blue = gamma, MK, red = log-normal, cyan = gamma, classical
%solid = kernels and 3-param. molc, dashed = series, dashed-dot = MKE

%target PDF legend text
switch distribution
    case "k"
        targetlgd = "Target $K$";
    case "gamma"
        targetlgd = "Target $\gamma$";
    case "g0"
        targetlgd = "Target $G^0$";
    case "ggd"
        targetlgd = "Target G$\Gamma$D";
    case "igam"
        targetlgd = "Target $\gamma^{-1}$";
end

%Deciding which methods to plot
%p = [1 2 3 4 5 6 7 8 9];    %All methods
p = [6 7 3 4 5 8 9];        %Thesis

%For the positions of the plot
lt = 0.08; %left
bt = 0.07; %bottom
wt = 0.4; %width
ht = 0.4; %height
hg = wt+0.09; %horizontal gap
vg = ht+0.10; %vertical gap

%Upper left: The methods, normal grid
subplot('Position', [lt bt+vg wt ht]); hold on; grid on; box on;
plot(x, f, 'k', 'LineWidth', 1);
for i = p
    plot(x, PDFapprox(:,i), methods(i).style, ...
        'Color', methods(i).color, 'LineWidth', 0.8);
end
title('PDF, True and Approximated', 'FontSize', 18)
xlabel '$x$'
ylabel '$f(x),\,\hat{f}(x)$'
lh = legend ([targetlgd, methods(p).name], 'Location', 'northeast');
set(lh,'Interpreter','latex')
%axis([x(1) x(end) 0 0.7]);
axis tight

%Upper right: The methods, x logarithmic scale
subplot('Position', [lt+hg bt+vg wt ht]); hold on; grid on; box on;
plot(x, f, 'k', 'LineWidth', 1);
for i = p
    plot(x, PDFapprox(:,i), methods(i).style, ...
        'Color', methods(i).color, 'LineWidth', 0.8);
end
title('PDF, True and Approximated, Logarithmic Scale', 'FontSize', 18)
xlabel '$x$, logarithmic scale'
ylabel '$f(x),\,\hat{f}(x)$'
%lh = legend ([targetlgd, methods(p).name], 'Location', 'northwest');
%set(lh,'Interpreter','latex')
axis tight
ax = gca;
auxfunction1(ax, -1:2);

%Lower left: Absolute errors, normal grid
subplot('Position', [lt bt wt ht]); hold on; grid on; box on;
plot([x(1) x(end)], [0 0], 'k', 'LineWidth', 1);
for i = p
    plot(x, (PDFapprox(:,i)-f)./f, methods(i).style, ...
        'Color', methods(i).color, 'LineWidth', 0.8);
end
title('Approximation Errors (Relative)', 'FontSize', 18)
xlabel '$x$'
ylabel '$\frac{\hat{f}(x)-f(x)}{f(x)}$'
%lh = legend ([methods(p).ll], methods(p).name, 'Location', 'southeast');
%set(lh,'Interpreter','latex');
axis([x(1) x(end) -0.31 0.31])
%axis tight

%Lower right: Absolute errors, x logarithmic scale
subplot('Position', [lt+hg bt wt ht]); hold on; grid on; box on;
plot([x(1) x(end)], [0 0], 'k', 'LineWidth', 1);
for i = p
    plot(x, PDFapprox(:,i)-f, methods(i).style, ...
        'Color', methods(i).color, 'LineWidth', 0.8);
end
title('Approximation Errors (Absolute), Logarithmic Scale', 'FontSize', 18)
xlabel '$x$, logarithmic scale'
ylabel '$\hat{f}(x)-f(x)$'
%lh = legend ([methods(p).lr], methods(p).name, 'Location', 'northwest');
%set(lh,'Interpreter','latex')

ax = gca;
auxfunction1(ax, -1:2);
ax.YLim = [-0.0055 0.0055];
ax.YTick = [-0.005 0.005];
ax.YTickLabel = {'-5\times10^{-3}';'5\times10^{-3}'};

%The line below moves the y-label to a nice spot, if necessary.
ax.YAxis.Label.Position = [x(1)/2 0 -1];

%A container for the numerical results
results = struct('Bhatta', cell(9,1), 'KL', cell(9,1), 'max', cell(9,1));

%Computing the distances
for m = 1:9
    results(m).Bhatta = bhattadist(f, PDFapprox(:,m), dx);%Bhattacharyya
    results(m).KL = kldist(f, PDFapprox(:,m), dx);        %Kullback Leibler
    results(m).max = max(abs(f-PDFapprox(:,m)));          %Maximum error
end

%Printing the results in a table
T = table([results.Bhatta]', [results.KL]', [results.max]', 'RowNames', ...
    {methods.name}', 'VariableNames', {'Bhattacharyya', ...
    'KullbackLeibler', 'Maximum'});
format shorte
disp(T);
disp '-----------------------------------------------------------------'

if writeTextFile == 1
    % Writing results to a text file
    fileID = fopen('../results/ComparingMethods.txt','w');
    %'w': discard, 'a': append
    %fprintf(fileID, ['\n \n', 'Maximum error', '\n']);
    fprintf(fileID, ['\t \t', 'Bhattacharyya distance', '\t \t', ...
        'Kullback-Leibler distance', '\t', 'Maximum error']);
    
    methods(2).name = 'MoLC GGD'; %\Gamma gives an error
    for m = 1:9
        fprintf(fileID, ['\n', methods(m).name,  '\t']);
        if (m == 1 || m == 4 || m == 5 || m == 7 || m == 9)
            fprintf(fileID, '\t');
        end
        fprintf(fileID, '%.3e \t \t \t', T{m,:});
        %.2e for two digit precision, .3e for three digits etc.
    end
    
    %Metadata
    fprintf(fileID,['\n \n', 'Highest order (log-)cumulant used: ', ...
        num2str(Nt), '\n',...
        'Data type: ', distribution, '\n', ...
        'Mu = ', num2str(mu), '\n',...
        'g = ', num2str(g), '\n',...
        'L = ', num2str(L), '\n',...
        'M = ', num2str(M), '\n',...
        'nu = ', num2str(nu), '\n',...
        'dx = ', num2str(dx), '\n',...
        'x in [', num2str(x(1)), ',', num2str(x(end)), '] \n',...
        '------ END OF RUN ------\n \n']);
    fclose(fileID);
end

%Saving the figure
if saveImages == 1
    tmp = getframe(gcf);
    imwrite(tmp.cdata, ['/Users/torgeirbrenn/Documents/', ...
        'Skole/Masteroppgave/Innlevering/figures/',...
        'CaseStudyApprox', distribution, '.png'])
end

%%% AUXILLARY FUNCTIONS %%
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/'
%auxfunction1(): Controlling the appearance of the logarithmic X axis
function auxfunction1(ax, ticks)
if nargin < 2
    ticks = -10:0.5:10; %I.e. more than enough
end

if nargin == 0
    ax = gca;
end
ax.XLim = [ax.Children(1).XData(1) ax.Children(1).XData(end)];
ax.XScale = 'log';
ax.XTick = 10.^(ticks);
ax.XTickLabel = cell(numel(ticks),1);
for i = 1:numel(ticks)
    ax.XTickLabel{i} = ['10^{', num2str(ticks(i)), '}'];
end
end