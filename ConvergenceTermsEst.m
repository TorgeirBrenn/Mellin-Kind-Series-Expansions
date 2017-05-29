%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script was made to synthesize non-negative data and estimate the
% (underlying) PDF using series expansion methods. The focus is on
% evaluating the convergence as the number of terms increase.
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
%set(groot,'defaultFigurePosition', [0 500 2000 400])
set(groot,'defaultFigurePosition', [0 500 560 840])
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/functions'

writeTextFile = 1;

%Declaring variables
Ni = 1e1;   %No. of iterations
Nd = 1e2;   %No. of data points
Nt = 8;     %Highest order (log-)moment/cumulant used. N<=8

distribution = 'k'; %Currently supported: 'k', 'gamma', 'g0', 'ggd', 'igam'

parameters = struct(...
    'm', 1, ...
    'g', 10, ...
    'L', 16, ...
    'M', 10, ...
    'nu', 2 ...
    );

[x, f] = targetpdf(distribution, parameters);


%Waitbar
wb = waitbar(0,'Computing');

%Preparing the holders of the estimates and measures
PDFests = zeros(numel(x), 5, Nt, Ni); %PDF estimates
R       = zeros(3, 5, Nt, Ni);        %Distance, divergence measures
%R will play a temporary role, with its means presented in 'results'

profile on
for ni = 1:Ni
    %Generating data
    data = generatedata(Nd, distribution, parameters);
    
    for nt = 2:Nt
        %Note: PDFests(:, :, 1, :) is left empty by design, as it would be
        %identical to PDFests(:, :, 2, :). This is done for the sake of
        %clarity in terms of consistent notation. Change if necessary.
        
        PDFests(:, 1, nt, ni) = gcgamfit(x, data, nt);
        PDFests(:, 2, nt, ni) = mklkfit(x, data, nt);
        PDFests(:, 3, nt, ni) = mkefit(x, data, nt);
        PDFests(:, 4, nt, ni) = mkgkfit(x, data, nt);
        PDFests(:, 5, nt, ni) = mkgkfitgivenparam(x, data, nt, ...
            parameters.L, parameters.m);
        waitbar(ni/Ni, wb);
    end
end
profile off
close(wb);
disp 'Computation complete'

% Plotting and printing results

%The struct "methods" defines the style.
methods = struct(...
    'name', {'Gamma GC'; 'MKLK series'; 'MKE series'; 'MKGK series'; ...
    'MKGK alt.'}, ...
    'style', {'--s'; '--+'; '-.+'; '--d'; '--*'}, ...
    'color', {[0.9290 0.6940 0.1250]; [0.8500 0.3250 0.0980]; ...
    'k'; [0 0.4470 0.7410]; [0.3010 0.7450 0.9330]} ...
    );

%Computing the distances
results = struct(...
    'Bhatta', cell(5,1), 'KL', cell(5,1), 'max', cell(5,1) ...
    );

for m = 1:5
    results(m).Bhatta = zeros(Nt,1);    %Bhattacharyya distance
    results(m).KL = zeros(Nt,1);        %Kullback-Leibler distance
    results(m).max = zeros(Nt,1);       %Maximum error
    for nt = 2:Nt
        for ni = 1:Ni
            R(1,m,nt,ni) = bhattadist(f, PDFests(:,m,nt,ni), x(2)-x(1));
            R(2,m,nt,ni) = kldist(f, PDFests(:,m,nt,ni), x(2)-x(1));
            R(3,m,nt,ni) = max(abs(f-PDFests(:,m,nt,ni)));
        end
        results(m).Bhatta(nt) = mean(real(R(1,m,nt,:)));
        results(m).KL(nt) = mean(real(R(2,m,nt,:)));
        results(m).max(nt) = mean(R(3,m,nt,:));
    end
end

%Deciding which methods to plot
p = [1 2 3 4 5];              %All methods
%p = [2 3 4];                 %MK series expansions
%p = [4 5 2 3];               %The order of the thesis

%Top: Bhattacharyya distance
subplot(2,1,1); hold on;
for i = p
    plot(2:Nt, results(i).Bhatta(2:Nt), ...
        methods(i).style, 'Color', methods(i).color, 'LineWidth', 1.2);
    
end
title 'Bhattacharyya'
lh = legend (methods(p).name, 'Location', 'northwest');
set(lh,'Interpreter','latex')
%axis([2 Nt 0 0.01]); grid on;
axis tight; grid on; box on;
auxfunction1();

%Center: Kullback-Leibler distance
%subplot('Position', [0.37 0.1 0.29 0.82]);
subplot(2,1,2); hold on;
for i = p
    plot(2:Nt, results(i).KL(2:Nt), ...
        methods(i).style, 'Color', methods(i).color, 'LineWidth', 1.2);
    
end
title 'Kullback-Leibler'
%lh = legend (methods(p).name, 'Location', 'northeast');
%set(lh,'Interpreter','latex');
axis tight; grid on; box on;
auxfunction1();
%
% %Right: Maximum distance
%subplot(); hold on; %NEED TO SPECIFY THE SUBPLOT BEFORE USE
% for i = p
%     plot(2:Nt, results(i).max(2:Nt), ...
%         methods(i).style, 'Color', methods(i).color, 'LineWidth', 1.2);
% end
% title 'Maximum'
% lh = legend (methods(p).name, 'Location', 'northeast');
% set(lh,'Interpreter','latex');
% axis tight; grid on; box on;
% auxfunction1();



%Printing the results in a table

rn = cell(Nt,1); %Row names
for nt = 1:Nt
    rn{nt} = ['N = ', num2str(nt)];
end
disp '                            Bhattacharyya distance'
Tb = table(results(1).Bhatta(2:Nt), results(2).Bhatta(2:Nt), ...
    results(3).Bhatta(2:Nt), results(4).Bhatta(2:Nt), ...
    results(5).Bhatta(2:Nt), 'VariableNames',  ...
    matlab.lang.makeValidName({methods.name}), 'RowNames', rn(2:Nt));
disp(Tb); disp ' ';

disp '                           Kullback-Leibler distance'
Tkl = table(results(1).KL(2:Nt), results(2).KL(2:Nt), ...
    results(3).KL(2:Nt), results(4).KL(2:Nt), results(5).KL(2:Nt), ...
    'VariableNames',  ...
    matlab.lang.makeValidName({methods.name}), 'RowNames', rn(2:Nt));
disp(Tkl); disp ' ';

disp '                                 Maximum error'
Tme = table(results(1).max(2:Nt), results(2).max(2:Nt), ...
    results(3).max(2:Nt), results(4).max(2:Nt), results(5).max(2:Nt), ...
    'VariableNames',  ...
    matlab.lang.makeValidName({methods.name}), 'RowNames', rn(2:Nt));
disp(Tme);
disp '-----------------------------------------------------------------'

% Writing results to a text file
if writeTextFile == 1
    fileID = fopen('../results/ConvergenceTerms.txt','w');
    %'w': discard, 'a': append
    
    %Bhattacharyya distance
    fprintf(fileID, 'Bhattacharyya distance \n');
    fprintf(fileID, ['\t \t' methods(1).name, '\t', methods(2).name, '\t \t',...
        methods(3).name, '\t \t', methods(4).name, '\t \t', methods(5).name]);
    for nt = 2:Nt
        fprintf(fileID, ['\n', rn{nt},  '\t \t']);
        fprintf(fileID, '%.3e \t', Tb{nt-1,:});
        %.2e for two digit precision, .3e for three digits etc.
    end
    
    %Kullback-Leibler distance
    fprintf(fileID, ['\n \n', 'Kullback-Leibler distance', '\n']);
    fprintf(fileID, ['\t \t' methods(1).name, '\t', methods(2).name, '\t \t',...
        methods(3).name, '\t \t', methods(4).name, '\t \t', methods(5).name]);
    for nt = 2:Nt
        fprintf(fileID, ['\n', rn{nt},  '\t \t']);
        fprintf(fileID, '%.3e \t', Tkl{nt-1,:});
    end
    
    %Maximum error
    fprintf(fileID, ['\n \n', 'Maximum error', '\n']);
    fprintf(fileID, ['\t \t' methods(1).name, '\t', methods(2).name, '\t \t',...
        methods(3).name, '\t \t', methods(4).name, '\t \t', methods(5).name]);
    for nt = 2:Nt
        fprintf(fileID, ['\n', rn{nt},  '\t \t']);
        fprintf(fileID, '%.3e \t', Tme{nt-1,:});
    end
    
    %Metadata
    fprintf(fileID,['\n \n', 'Number of iterations: ', num2str(Ni), '\n',...
        'Number of data points: ', num2str(Nd), '\n',...
        'Data type: ', distribution, '\n', ...
        'm = ', num2str(parameters.m), '\n',...
        'g = ', num2str(parameters.g), '\n',...
        'L = ', num2str(parameters.L), '\n',...
        'M = ', num2str(parameters.M), '\n',...
        'nu = ', num2str(parameters.nu), '\n',...
        'dx = ', num2str(x(2)-x(1)), '\n',...
        'x in [', num2str(x(1)), ',', num2str(x(numel(x))), '] \n',...
        '------ END OF RUN ------\n \n']);
    fclose(fileID);
end

%%% AUXILLARY FUNCTIONS %%
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/'
%auxfunction1: Controlling the appearance of the logarithmic Y axis
function auxfunction1(ax, ticks)
if nargin < 2
    ticks = -10:1:2; %I.e. more than enough
end

if nargin == 0
    ax = gca;
end

xlabel 'Highest order (log-)cumulant corrected for'
ylabel 'Distance'

ax.YScale = 'log';
ax.YTick = 10.^(ticks);
ax.YTickLabel = cell(numel(ticks),1);
for i = 1:numel(ticks)
    ax.YTickLabel{i} = ['10^{', num2str(ticks(i)), '}'];
end
end