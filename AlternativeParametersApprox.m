%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we analyze how non-trivial choices of the kernel parameters affect
% the mellin kind series expansions.
%
% As this script is very long, it would be nice if the actual plotting was
% moved to an auxillary function, if possible.
%
% Made by Torgeir Brenn, 2017. For the custom functions, the author is
% specified within those functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear variables
set(groot, 'defaultFigureColor', [1 1 1]);   %White figure background
set(groot,'defaulttextinterpreter','latex'); %Figure text
set(0,'defaultAxesFontSize',13);
set(groot,'defaultFigurePosition', [0 500 1600 800])
%set(groot,'defaultFigurePosition', [0 500 560 420])
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/functions'

saveImages = 1;

%Declaring variables
Nt = [0 2 4 6];     %Highest order (log-)moment/cumulant used. 
parameters = struct(...
        'm', 10, ...
        'g', 10, ...
        'L', 16, ...
        'M', -10, ...
        'nu', 0.5 ...
);
distribution = 'g0'; %Supported: 'k', 'gamma', 'g0', 'ggd', 'igam'

global ParamVals %Global since we need this in a function
ParamVals = [0.5 0.65 0.8 0.9 1 1.1 1.25 1.5 2]; %Scaled with true values
NoParamVals = numel(ParamVals); %Number of parameter values

[x, f, kappa] = targetpdf(distribution, parameters, Nt(end));
dx = x(2)-x(1);

%Cosmetics, legends
methods = struct(...
    'name', {'MKLK'; 'MKE';'MKGK';'MKBK'}, ...
    'color', {[0.8500 0.3250 0.0980]; 'k'; [0 0.4470 0.7410]; ...
              [0.4660 0.6740 0.1880]} ...
);
style = {'-o', '--^', '-.p', ':*'};

lgdtext = strings(numel(Nt),1);
lgdtext(1) = 'Kernel';
for nt = 2:numel(Nt)
    lgdtext(nt) = ['$N=', num2str(Nt(nt)), '$'];
end

%For the positions of the plot
lt = 0.05; %left
bt = 0.05; %bottom
wt = 0.21; %width
ht = 0.26; %height
hg = wt+0.03; %horizontal gap
vg = ht+0.05; %vertical gap

%In case we want to omit the MKBK series
if distribution == "gamma" || distribution ==  "igam"
    lt = 0.05; %left
    bt = 0.05; %bottom
    wt = 0.29; %width
    ht = 0.26; %height
    hg = wt+0.03; %horizontal gap
    vg = ht+0.05; %vertical gap
end

%MKGK
syms s1
Lc = solve(psi(1,s1)-kappa(2), s1); Lc = double(Lc);
mc = Lc*exp(kappa(1)-psi(0,Lc));
bc = Lc/mc;
Lvec = ParamVals*Lc;
bvec = ParamVals*bc;
mvec = ParamVals*mc;

MKGK = struct(...
    'L', zeros(NoParamVals,numel(Nt)), ...
    'b', zeros(NoParamVals,numel(Nt)), ...
    'm', zeros(NoParamVals,numel(Nt))...
    );

for l = 1:NoParamVals
    for nt = 1:numel(Nt)
    MKGK(l,nt).L = kldist(f, ...
        mkgkfitgivenparam(x, kappa, Nt(nt), Lvec(l), mc), dx);
    
    MKGK(l,nt).b = kldist(f, ...
        mkgkfitgivenparam(x, kappa, Nt(nt), Lc, Lc/bvec(l)), dx);
    
    MKGK(l,nt).m = kldist(f, ...
        mkgkfitgivenparam(x, kappa, Nt(nt), Lc, mvec(l)), dx);
    end
end

%Plotting, L
figure;
subplot('Position', [lt+0*hg, bt+2*vg, wt, ht]); hold on;
for nt = 1:numel(Nt)
    plot(ParamVals, [MKGK(:,nt).L], style{nt}, 'Color', methods(3).color, 'LineWidth', 0.8);
end
auxfunction1();
xlabel(['$a(=L)$, relative to the tailored value ', num2str(Lc)]);
title 'MKGK series'
lh = legend(lgdtext, 'Location', 'southeast');
set(lh,'Interpreter','latex')

%b
subplot('Position', [lt+0*hg, bt+1*vg, wt, ht]); hold on;
for nt = 1:numel(Nt)
    plot(ParamVals, [MKGK(:,nt).b], style{nt}, 'Color', methods(3).color, 'LineWidth', 0.8); 
end
auxfunction1();
xlabel(['$b$, relative to the tailored value ', num2str(Lc/mc)]);
lh = legend('Kernel', ['$N=', num2str(Nt(2)), '$'], ...
    ['$N=', num2str(Nt(3)), '$'], ['$N=', num2str(Nt(4)), '$'], ...
    'Location', 'southeast');
set(lh,'Interpreter','latex')

%m
subplot('Position', [lt+0*hg, bt+0*vg, wt, ht]); hold on;
for nt = 1:numel(Nt)
    plot(ParamVals, [MKGK(:,nt).m], style{nt}, 'Color', methods(3).color, 'LineWidth', 0.8);
end
auxfunction1();
xlabel(['$m$, relative to the tailored value ', num2str(mc)]);
lh = legend(lgdtext, 'Location', 'southeast');
set(lh,'Interpreter','latex')

%MKLK
muc = kappa(1);
logvarc = kappa(2);
muvec = ParamVals*muc;
logvarvec = ParamVals*logvarc;

MKLK = struct(...
    'mu',     zeros(NoParamVals,numel(Nt)), ...
    'logvar', zeros(NoParamVals,numel(Nt)) ...
);

for l = 1:NoParamVals
    for nt = 1:numel(Nt)
        MKLK(l,nt).mu = kldist(f, ...
            mklkfitgivenparam(x, kappa, Nt(nt), muvec(l), logvarc), dx);
        
        MKLK(l,nt).logvar = kldist(f, ...
            mklkfitgivenparam(x, kappa, Nt(nt), muc, logvarvec(l)), dx);
    end
end

%Plotting, mu
subplot('Position', [lt+1*hg, bt+2*vg, wt, ht]); hold on;
for nt = 1:numel(Nt)
    plot(ParamVals, [MKLK(:,nt).mu], style{nt}, 'Color', methods(1).color, 'LineWidth', 0.8);
end
auxfunction1();
xlabel(['$\mu$, relative to the tailored value ', num2str(muc)]);
title 'MKLK series'
lh = legend(lgdtext, 'Location', 'southeast');
set(lh,'Interpreter','latex')

%logvar
subplot('Position', [lt+1*hg, bt+1*vg, wt, ht]); hold on;
for nt = 1:numel(Nt)
    plot(ParamVals, [MKLK(:,nt).logvar], ...
        style{nt}, 'Color', methods(1).color, 'LineWidth', 0.8); 
end
auxfunction1();
xlabel(['$\sigma^2$, relative to the tailored value ', num2str(logvarc)]);
lh = legend(lgdtext, 'Location', 'southeast');
set(lh,'Interpreter','latex')

%MKE
muc = kappa(1);
logvarc = kappa(2);
muvec = ParamVals*muc;
logvarvec = ParamVals*logvarc;

MKE = struct(...
    'mu',     zeros(NoParamVals,numel(Nt)), ...
    'logvar', zeros(NoParamVals,numel(Nt)) ...
    );
    

for l = 1:NoParamVals
    for nt = [1 3 4] %Since Nt = 2 is the same as the kernel
    MKE(l,nt).mu = kldist(f, ...
        mkefitexact(x, [muvec(l); logvarc; kappa(3:Nt(nt))], Nt(nt)), dx);
    
    MKE(l,nt).logvar = kldist(f, ...
        mkefitexact(x, [muc; logvarvec(l); kappa(3:Nt(nt))], Nt(nt)), dx);
    end
end

%Plotting, mu
subplot('Position', [lt+2*hg, bt+2*vg, wt, ht]); hold on;
for nt = [1 3 4]
    plot(ParamVals, [MKE(:,nt).mu], style{nt}, 'Color', methods(2).color, 'LineWidth', 0.8); 
end
auxfunction1();
xlabel(['$\mu$, relative to the tailored value ', num2str(muc)]);
title 'MKE series'
lh = legend('Kernel', 'Discard O$(r^{-3/2})$', ...
    'Discard O$(r^{-5/2})$', 'Location', 'southeast');
set(lh,'Interpreter','latex')

%logvar
subplot('Position', [lt+2*hg, bt+1*vg, wt, ht]); hold on;
for nt = [1 3 4]
    plot(ParamVals, [MKE(:,nt).logvar], ...
        style{nt}, 'Color', methods(2).color, 'LineWidth', 0.8);
end
auxfunction1();
xlabel(['$\sigma^2$, relative to the tailored value ', num2str(logvarc)]);
lh = legend('Kernel', 'Discard O$(r^{-3/2})$', ...
    'Discard O$(r^{-5/2})$', 'Location', 'southeast');
set(lh,'Interpreter','latex')


%MKBK
if distribution == "ggd" || distribution == "k" || distribution == "g0"
    [a1c,a2c] = mexBPdistGradOpt(3,5,kappa(2),kappa(3));
    
    %Alternative, slower but more precise (use for the G^0 distribution!)
    if distribution == "g0"
        syms s1 s2
        [a1c, a2c] = solve([psi(1,s1)+psi(1,s2)==kappa(2), ...
            psi(2,s1)-psi(2,s2)==kappa(3)], [s1, s2]);
        a1c = double(a1c);
        a2c = double(a2c);
    end
    
    bc = exp(psi(0,a1c)-psi(0,a2c)-kappa(1));
    a1vec = ParamVals*a1c;
    a2vec = ParamVals*a2c;
    bvec = ParamVals*bc;
    
    MKBK = struct(...
        'a1', zeros(NoParamVals,numel(Nt)), ...
        'a2', zeros(NoParamVals,numel(Nt)), ...
        'b',  zeros(NoParamVals,numel(Nt)) ...
        );
    
    for l = 1:NoParamVals
        for nt = 1:numel(Nt)
        MKBK(l,nt).a1 = kldist(f, ...
            mkbkfitgivenparam(x, kappa, Nt(nt), a1vec(l), a2c, bc), dx);
        
        MKBK(l,nt).a2 = kldist(f, ...
            mkbkfitgivenparam(x, kappa, Nt(nt), a1c, a2vec(l), bc), dx);
        
        MKBK(l,nt).b = kldist(f, ...
            mkbkfitgivenparam(x, kappa, Nt(nt), a1c, a2c, bvec(l)), dx);
        end
    end
    
    %Plotting, a1
    subplot('Position', [lt+3*hg, bt+2*vg, wt, ht]); hold on;
    for nt = 1:numel(Nt)
        plot(ParamVals, [MKBK(:,nt).a1], ...
            style{nt}, 'Color', methods(4).color, 'LineWidth', 0.8);
    end
    auxfunction1();
    xlabel(['$a_1$, relative to the tailored value ', num2str(a1c)]);
    title 'MKBK series'
    lh = legend(lgdtext, 'Location', 'southeast');
    set(lh,'Interpreter','latex')
    
    %a2
    subplot('Position', [lt+3*hg, bt+1*vg, wt, ht]); hold on;
    for nt = 1:numel(Nt)
        plot(ParamVals, [MKBK(:,nt).a2], ...
            style{1}, 'Color', methods(4).color, 'LineWidth', 0.8);
    end
    auxfunction1();
    xlabel(['$a_2$ relative to the tailored value ', num2str(a2c)]);
    lh = legend(lgdtext, 'Location', 'southeast');
    set(lh,'Interpreter','latex')
    
    %b
    subplot('Position', [lt+3*hg, bt+0*vg, wt, ht]); hold on;
    for nt = 1:numel(Nt)
        plot(ParamVals, [MKBK(:,nt).b], ...
            style{nt}, 'Color', methods(4).color, 'LineWidth', 0.8); 
    end
    auxfunction1();
    xlabel(['$b$, relative to the tailored value ', num2str(bc)]);
    lh = legend(lgdtext, 'Location', 'southeast');
    set(lh,'Interpreter','latex')
    
end

%Title for the entire figure, based on the MathWorks Central function
%suplabel
titletext = pdfname(distribution, parameters, 0);
[~,suptitle] = suplabel(titletext, 't');
%suptitle.Position(2) = suptitle.Position(2) - 0.05;
suptitle.FontSize = 22;
[~,supylabel] = suplabel('Distance $d_{\rm KL}$','y');
supylabel.Position(1) = supylabel.Position(1) + 0.03;
supylabel.FontSize = 22;
% Use "ax.MinorGridAlpha = 0;" to hide minor grids as needed.

%Saving the figure
if saveImages == 1
    tmp = getframe(gcf);
    imwrite(tmp.cdata, ['/Users/torgeirbrenn/Documents/', ...
        'Skole/Masteroppgave/Innlevering/figures/',...
        'AlternativeParametersApprox', distribution, '.png'])
end

%%% AUXILLARY FUNCTIONS %%
cd '/Users/torgeirbrenn/Documents/Skole/Masteroppgave/MATLAB/'
%auxfunction1(): The plotting commands which are shared by all plots
function auxfunction1()
global ParamVals
axis tight; grid on; box on;
ax = gca; ax.XTick = ParamVals; ax.XAxis.TickLabelRotation=0; 
ax.YScale = 'log';
end