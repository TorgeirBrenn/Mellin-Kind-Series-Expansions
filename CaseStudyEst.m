%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script was made to synthesize non-negative data and estimate the
% (underlying) PDF using different methods. The focus is on a broad
% comparison of the methods.
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
set(groot, 'defaultFigureColor', [1 1 1]);   %White figure background
set(groot,'defaulttextinterpreter','latex'); %Figure text
set(0,'defaultAxesFontSize',16);
set(groot,'defaultFigurePosition', [250 0 1000 1000])
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/functions'

saveImages    = 1;
writeTextFile = 0; %1: Will overwrite the given text file with new results.

%Declaring variables
Nt = 4;     %Highest order (log-)moment/cumulant used. N<=8
Ni = 1e0;   %Number of iterations
Nd = 1e3;   %Number of data points
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
PDFests = zeros(numel(x), 11, Ni); %PDF estimates
results = zeros(3, 11, Ni);        %Distance, divergence measures

%Waitbar
wb = waitbar(0,'Computing');

profile on
for ni = 1:Ni
    %Generating data
    data = generatedata(Nd, distribution, parameters);
    
    %%% METHOD 1 - ML Gamma distribution %%%
    % Gamma distribution fitted using MLEs of parameters. Choi & Wette
    %(1969) discuss this procedure.
    PDFests(:, 1, ni) = mlgamfit(x, data);
    
    %%% METHOD 2 - Laguerre, Gamma Gram-Charlier %%%
    % As mentioned in Kendall's (p. 236), the Gram-Charlier/Edgeworth
    % expansion can be done with the Gamma distribution using the Laguerre
    % polynomials. In this method, the classical moments are used.
    PDFests(:, 2, ni) = gcgamfit(x, data, Nt);
    
    %%% METHOD 3 - MoLC K distribution %%%
    % Estimating the parameters in the K distribution. Analytically and
    % computationally demanding. Optimised by programming parts of the code
    % in C (Author: SNA).
    PDFests(:, 3, ni) = molckfit(x, data);
    
    %%% METHODS 4 - MoLC Generalized Gamma Distribution %%%
    %Fitted generalized Gamma distribution, using MoLC. See Li et al., 2011
    PDFests(:, 4, ni) = molcggdfit(x, data);
    
    %%% METHOD 5 - MoLC Log-normal distribution %%%
    %Estimating the parameters in the log-normal distribution, based on the
    %log-cumulants.
    PDFests(:, 5, ni) = molclognormfit(x, data);
    
    %%% METHOD 6 - MKLK series %%%
    %MKLK Gram-Charlier.
    PDFests(:, 6, ni) = mklkfit(x, data, Nt);
    
    %%% METHOD 7 - MKE series %%%
    %Pastor's method (2014, 2016) which is a series expansion around a
    %fitted log-normal kernel, using the log-cumulants.
    PDFests(:, 7, ni) = mkefit(x, data, Nt);
    
    %%% METHOD 8 - MoLC Gamma distribution %%%
    %Estimating the parameters in the Gamma distribution, based on the
    %log-cumulants.
    PDFests(:, 8, ni) = molcgamfit(x, data);
    
    %%% METHOD 9 - MKGK series %%%
    %Novel series expansion using the Gamma kernel and Laguerre polynomials
    %with the log-cumulants.
    PDFests(:, 9, ni) = mkgkfit(x, data, Nt);
    
    %%% METHOD 10 - MoLC Beta prime distribution %%%
    %Estimating the parameters in the beta prime distribution using the
    %MoLC.
    PDFests(:, 10, ni) = mkbkfit(x, data, 2);
    
    %%% METHOD 11 - MKBK series %%%
    %Novel series expansion using the beta prime kernel, the Mprime
    %polynomials, and the log-cumulants.
    PDFests(:, 11, ni) = mkbkfit(x, data, Nt);
    
    waitbar(ni/Ni, wb);
end
profile off
close(wb);
disp 'Computation complete'

%The struct "methods" defines the style.
methods = struct(...
    'name', {'ML Gamma'; 'Gamma GC'; 'MoLC K'; 'MoLC G{\Gamma}D'; ...
    'MKLK kernel'; 'MKLK series'; 'MKE series'; 'MKGK kernel'; ...
    'MKGK series'; 'MKBK kernel'; 'MKBK series'}, ...
    'style', {''; '--'; ''; ''; ''; '--'; '-.'; ''; '--'; ''; '--'}, ...
    'color', {[0.9290 0.6940 0.1250]; [0.9290 0.6940 0.1250]; ...
    [0.4940 0.1840 0.5560]; [0.6350 0.0780 0.1840]; ...
    [0.8500 0.3250 0.0980]; [0.8500 0.3250 0.0980]; ...
    [0.8500 0.3250 0.0980]; [0 0.4470 0.7410]; [0 0.4470 0.7410];...
    [0.4660 0.6740 0.1880]; [0.4660 0.6740 0.1880]}...
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
%p = [1 2 3 4 5 6 7 8 9]; %All methods
p = [8 9 5 6 7 10 11];    %Thesis

%For the positions of the plot
lt = 0.08; %left
bt = 0.07; %bottom
wt = 0.4; %width
ht = 0.4; %height
hg = wt+0.09; %horizontal gap
vg = ht+0.10; %vertical gap

%Upper left: The methods, normal grid
subplot('Position', [lt bt+vg wt ht]);
plot(x, f, 'k', 'LineWidth', 1); hold on; grid on;
for i = p
    plot(x, PDFests(:,i,Ni), methods(i).style, ...
        'Color', methods(i).color, 'LineWidth', 0.8);
end
title('PDF, True and Estimated', 'FontSize', 18)
xlabel '$x$'
ylabel '$f(x),\,\hat{f}(x)$'
lh = legend ([targetlgd, methods(p).name], 'Location', 'northeast');
set(lh,'Interpreter','latex')
%axis([x(1) x(end) 0 0.7]);
axis tight

%Upper right: The methods, x logarithmic scale
subplot('Position', [lt+hg bt+vg wt ht]);
plot(x, f, 'k', 'LineWidth', 1); hold on; grid on;
for i = p
    plot(x, PDFests(:,i,Ni), methods(i).style, ...
        'Color', methods(i).color, 'LineWidth', 0.8);
end
title('PDF, True and Estimated, Logarithmic Scale', 'FontSize', 18)
xlabel '$x$, logarithmic scale'
ylabel '$f(x),\,\hat{f}(x)$'
%lh = legend([targetlgd, methods(p).name], 'Location', 'northwest');
%set(lh,'Interpreter','latex')
ax = gca;
axis tight;
auxfunction1(ax, -1:2);

%Lower left: Absolute errors, normal grid
subplot('Position', [lt bt wt ht]);
plot([x(1) x(end)], [0 0], 'k', 'LineWidth', 1); hold on; grid on;
for i = p
    plot(x, (PDFests(:,i,Ni)-f)./f, methods(i).style, ...
        'Color', methods(i).color, 'LineWidth', 0.8);
end
title('Estimation Errors (Relative)', 'FontSize', 18)
xlabel '$x$'
ylabel '$\frac{\hat{f}(x)-f(x)}{f(x)}$'
%lh = legend ([methods(p).ll], methods(p).name, 'Location', 'southeast');
%set(lh,'Interpreter','latex')
axis([x(1) x(end) -0.41 0.41])

%Lower right: Absolute errors, x logarithmic scale
subplot('Position', [lt+hg bt wt ht]);
plot([x(1) x(end)], [0 0], 'k', 'LineWidth', 1); hold on; grid on;
for i = p
    plot(x, PDFests(:,i,Ni)-f, methods(i).style, ...
        'Color', methods(i).color, 'LineWidth', 0.8);
end
title('Estimation Errors (Absolute), Logarithmic Scale', 'FontSize', 18)
xlabel '$x$, logarithmic scale'
ylabel '$\hat{f}(x)-f(x)$'
%lh = legend ([methods(p).lr], methods(p).name, 'Location', 'northwest');
%set(lh,'Interpreter','latex')
axis([0 x(end) -0.007 0.007])
ax = gca;
auxfunction1(ax,-1:2);
ax.YTick = [-0.01 -0.005 0.005 0.01];
ax.YTickLabel = {'-1\times10^{-2}';'-5\times10^{-3}';...
    '5\times10^{-3}';'1\times10^{-2}'};
%ax.YAxis.Exponent = 0;

%The line below moves the y-label to a nice spot, if necessary.
ax.YAxis.Label.Position = [x(1)*0.8 0 -1];

% Computing the distances
for m = 1:11
    for ni = 1:Ni
        results(1,m,ni) = bhattadist(f, PDFests(:,m,ni), dx);    %Bhattacharyya
        results(2,m,ni) = kldist(f, PDFests(:,m,ni), dx);        %Kullback Leibler
        results(3,m,ni) = max(abs(f-PDFests(:,m,ni)));           %Maximum error
    end
    %     methods(m).Bhatta = mean(real(R(1,m,:)));
    %     methods(m).KL = mean(real(R(2,m,:)));
    %     methods(m).max = mean(R(3,m,:));
    
    %If NaN's are a problem (usually MoLC K)
    methods(m).Bhatta = mean(real(results(1,m,:)), 'omitnan');
    methods(m).KL = mean(real(results(2,m,:)), 'omitnan');
    methods(m).max = mean(results(3,m,:), 'omitnan');
end

%Printing the results in a table
T = table([methods.Bhatta]', [methods.KL]', [methods.max]', 'RowNames', ...
    {methods.name}', 'VariableNames', {'Bhattacharyya', ...
    'KullbackLeibler', 'Maximum'});
format shorte
disp(T);
disp '-----------------------------------------------------------------'

% Writing results to a text file
if writeTextFile == 1
    fileID = fopen('../results/ComparingMethods.txt','w');
    %'w': discard, 'a': append
    %fprintf(fileID, ['\n \n', 'Maximum error', '\n']);
    fprintf(fileID, ['\t \t', 'Bhattacharyya distance', '\t \t', ...
        'Kullback-Leibler distance', '\t', 'Maximum error']);
    
    methods(4).name = 'MoLC GGD'; %\Gamma gives an error
    for m = 1:11
        fprintf(fileID, ['\n', methods(m).name,  '\t']);
        if (m == 3 || m == 6 || m == 7 || m == 9 || m == 11)
            fprintf(fileID, '\t');
        end
        fprintf(fileID, '%.3e \t \t \t', T{m,:});
        %.2e for two digit precision, .3e for three digits etc.
    end
    
    %Metadata
    fprintf(fileID,['\n \n', 'Number of iterations: ', num2str(Ni), '\n',...
        'Number of data points: ', num2str(Nd), '\n',...
        'Highest order (log-)cumulant used: ', num2str(Nt), '\n',...
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
        'CaseStudyEst', distribution, '.png'])
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