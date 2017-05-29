%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script was made to synthesize non-negative data and estimate the
% (underlying) PDF using series expansion methods. The focus is on
% evaluating the convergence as the number of datapoints increase.
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
set(0,'defaultAxesFontSize',14);
%set(groot,'defaultFigurePosition', [0 500 2000 400])
%set(groot,'defaultFigurePosition', [0 500 560 420])
set(groot,'defaultFigurePosition', [0 500 560 840])
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/functions'

%saveImages    = 1; %Not yet specified a file name (not used in the thesis)
writeTextFile = 0;

%Declaring variables
Ni = 1e1;   %No. of iterations
Nd = [1e2 5e2 1e3];   %No. of data points, vector
Nt = 4; %Highest order (log-)moment/cumulant used. N<=8
parameters = struct(...
    'm', 10, ...
    'g', 2, ...
    'L', 16, ...
    'M', -10, ...
    'nu', 0.5 ...
    );
distribution = 'ggd'; %Supported: 'k', 'gamma', 'g0', 'ggd', 'igam'

%Target distribution
[x, f] = targetpdf(distribution, parameters, Nt(end));
dx = x(2)-x(1);

%Waitbar
wb = waitbar(0,'Preparing...');

%Preparing the holders of the estimates and measures
PDFests = zeros(numel(x), 4, numel(Nd), Ni); %PDF estimates
results = zeros(3, 4, numel(Nd), Ni);        %Distance, divergence measures

profile on
for nd = 1:numel(Nd)
    waitbar(0, wb, ['Computing ',num2str(nd),' of ',num2str(numel(Nd))]);
    for ni = 1:Ni
        %Generating data
        data = generatedata(Nd(nd), distribution, parameters);
        
        PDFests(:, 1, nd, ni) = mlgamfit(x, data);
        PDFests(:, 2, nd, ni) = gcgamfit(x, data, Nt);
        PDFests(:, 3, nd, ni) = molckfit(x, data);
        PDFests(:, 4, nd, ni) = molcggdfit(x, data);
        PDFests(:, 5, nd, ni) = molclognormfit(x, data);
        PDFests(:, 6, nd, ni) = mklkfit(x, data, Nt);
        PDFests(:, 7, nd, ni) = mkefit(x, data, Nt);
        PDFests(:, 8, nd, ni) = molcgamfit(x, data);
        PDFests(:, 9, nd, ni) = mkgkfit(x, data, Nt);
        
        waitbar(ni/Ni, wb);
    end
end
profile off
close(wb);
disp 'Computation complete'

% Plotting and printing results
close all

%The struct "methods" defines the style.
methods = struct(...
    'name', {'ML Gamma'; 'Gamma GC'; 'MoLC K'; 'MoLC G{\Gamma}D'; ...
    'MKLK kernel'; 'MKLK series'; 'MKE series'; 'MKGK kernel'; ...
    'MKGK series'}, ...
    'style', {'-s';'--s';'-x';'-x';'-+';'--+';'-.+';'-d';'--d'}, ...
    'color', {[0.9290 0.6940 0.1250]; [0.9290 0.6940 0.1250]; ...
    [0.4940 0.1840 0.5560]; [0.4660 0.6740 0.1880]; ...
    [0.8500 0.3250 0.0980]; [0.8500 0.3250 0.0980]; ...
    [0.8500 0.3250 0.0980]; [0 0.4470 0.7410]; [0 0.4470 0.7410]} ...
    );

%Computing the distances
for m = 1:9
    for nd = 1:numel(Nd)
        for ni = 1:Ni
            results(1,m,nd,ni) = bhattadist(f, PDFests(:,m,nd,ni), dx);
            results(2,m,nd,ni) = kldist(f, PDFests(:,m,nd,ni), dx);
            results(3,m,nd,ni) = max(abs(f-PDFests(:,m,nd,ni)));
        end
        
        %If NaN's are a problem (usually MoLC K)
        methods(m).Bhatta(nd) = mean(real(results(1,m,nd,:)), 'omitnan');
        methods(m).KL(nd) = mean(real(results(2,m,nd,:)), 'omitnan');
        methods(m).max(nd) = mean(results(3,m,nd,:), 'omitnan');
    end
end

%

%Deciding which methods to plot
p = [1 2 4 5 6 7 8 9]; %All methods
%p = [6 7 9];             %MK series expansions excl. kernels
%p = [5 6 7 8 9];         %MK series expansions incl. kernels
%p = [2 6 7 9];           %All series expansions
%p = [3 4 7];             %MoLC K, GGD, MKE
%p = [1 3 4 5 6 7 8]; %Excluding GC Gamma, MKGK


%Left: Bhattacharyya distance
%subplot('Position', [0.035 0.1 0.29 0.82]);
subplot(2,1,1); hold on;
for i = p
    plot(Nd, methods(i).Bhatta(1:numel(Nd)), ...
        methods(i).style, 'Color', methods(i).color, 'LineWidth', 1.2);
end
title('Bhattacharyya', 'FontSize', 18)
xlabel 'Number of datapoints'
ylabel 'Distance'
legend (methods(p).name, 'Location', 'southwest')
axis([Nd(1) Nd(end) 1e-4 7e-3]); grid on;
%axis tight; grid on;
ax = gca; ax.XTick = Nd; %ax.YLim(1) = 0;
ax.XScale = 'log'; ax.YScale = 'log';

%Center: Kullback-Leibler distance
%subplot('Position', [0.37 0.1 0.29 0.82]);
subplot(2,1,2); hold on;
for i = p
    plot(Nd, methods(i).KL(1:numel(Nd)), ...
        methods(i).style, 'Color', methods(i).color, 'LineWidth', 1.2);
end
title('Kullback-Leibler', 'FontSize', 18)
xlabel 'Number of datapoints'
ylabel 'Distance'
%legend (methods(p).name, 'Location', 'northeast')
axis([Nd(1) Nd(end) 1e-4 5e-2]); grid on;
%axis tight; grid on;
ax = gca; ax.XTick = Nd; %ax.YLim(1) = 0;
ax.XScale = 'log'; ax.YScale = 'log';
%
% %Right: Maximum distance
% subplot('Position', [0.7 0.1 0.29 0.82]); hold on;
% for i = p
%     plot(Nd, methods(i).max(1:numel(Nd)), ...
%         methods(i).style, 'Color', methods(i).color, 'LineWidth', 0.6);
% end
% title('Maximum', 'FontSize', 18)
% xlabel 'Number of datapoints'
% ylabel 'Distance'
% legend (methods(p).name, 'Location', 'northeast')
% axis tight; grid on;
% ax = gca; ax.XTick = Nd; ax.YLim(1) = 0;
% ax.XScale = 'log'; %ax.YScale = 'log';



% Printing the results in a table

rn = cell(1,numel(Nd)); %Row names
for nd = 1:numel(Nd)
    rn{nd} = ['Nd = ', num2str(Nd(nd))];
end
disp 'Bhattacharyya distance'
Tb = table(methods.Bhatta, 'VariableNames',  ...
    matlab.lang.makeValidName({methods.name}), 'RowNames', rn);
disp(Tb); disp ' ';

disp 'Kullback-Leibler distance'
Tkl = table(methods.KL, 'VariableNames',  ...
    matlab.lang.makeValidName({methods.name}), 'RowNames', rn);
disp(Tkl); disp ' ';

disp 'Maximum error'
Tme = table(methods.max, 'VariableNames',  ...
    matlab.lang.makeValidName({methods.name}), 'RowNames', rn);
disp(Tme);
disp '-----------------------------------------------------------------'

if writeTextFile == 1
    % Writing results to a text file
    fileID = fopen('../results/ConvergenceDatapoints.txt','w');
    %'w': discard, 'a': append
    methods(4).name = 'MoLC GGD'; %\Gamma gives an error
    
    %Bhattacharyya distance
    fprintf(fileID, 'Bhattacharyya distance \n \n \t \t');
    fprintf(fileID, 'N = %d \t', Nd);
    for m = 1:9
        fprintf(fileID, ['\n', methods(m).name,  '\t']);
        if (m == 3 || m == 6 || m == 7 || m == 9)
            fprintf(fileID, '\t');
        end
        fprintf(fileID, '%.3e \t', Tb{:,m});
        %.2e for two digit precision, .3e for three digits etc.
    end
    
    
    %Kullback-Leibler distance
    fprintf(fileID, ['\n \n', 'Kullback-Leibler distance \n \n \t \t']);
    fprintf(fileID, 'N = %d \t', Nd);
    for m = 1:9
        fprintf(fileID, ['\n', methods(m).name,  '\t']);
        if (m == 3 || m == 6 || m == 7 || m == 9)
            fprintf(fileID, '\t');
        end
        fprintf(fileID, '%.3e \t', Tkl{:,m});
        %.2e for two digit precision, .3e for three digits etc.
    end
    
    %Maximum error
    fprintf(fileID, ['\n \n', 'Maximum error \n \n \t \t']);
    fprintf(fileID, 'N = %d \t', Nd);
    for m = 1:9
        fprintf(fileID, ['\n', methods(m).name,  '\t']);
        if (m == 3 || m == 6 || m == 7 || m == 9)
            fprintf(fileID, '\t');
        end
        fprintf(fileID, '%.3e \t', Tme{:,m});
        %.2e for two digit precision, .3e for three digits etc.
    end
    
    %Metadata
    fprintf(fileID,['\n \n', 'Number of iterations: ', num2str(Ni), ...
        '\n', 'Highest order (log-)cumulant used: ', num2str(Nt), '\n',...
        'Data type: ', distribution, '\n', ...
        'Mu = ', num2str(mu), '\n',...
        'g = ', num2str(g), '\n',...
        'L = ', num2str(L), '\n',...
        'M = ', num2str(M), '\n',...
        'nu = ', num2str(nu), '\n',...
        'dx = ', num2str(dx), '\n',...
        'x in [', num2str(x(1)), ',', num2str(x(numel(x))), '] \n',...
        '------ END OF RUN ------\n \n']);
    fclose(fileID);
end