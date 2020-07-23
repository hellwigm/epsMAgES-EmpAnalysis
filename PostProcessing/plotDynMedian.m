clear all, clc

%% build median dynamics
%
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

ESvariant{1}        = 'epsMAgES';
ESvariant{2}        = 'epsMAES';
ESvariant{3}        = 'epsMAgESwor';
ESvariant{4}        = 'epsMAgESnl';
ESvariant{5}        = 'epsSAgES';
ESvariant{6}        = 'lexMAgES';
ESvariant{7}        = 'lexMAES';

ESvariantL{1}        = '$\epsilon$MAg-ES';
ESvariantL{2}        = '$\epsilon$MA-ES';
ESvariantL{3}        = '$\epsilon$MAg-ES w/o';
ESvariantL{4}        = '$\epsilon$MAg-ES nl';
ESvariantL{5}        = '$\epsilon$SAg-ES';
ESvariantL{6}        = '$lex$MAg-ES';
ESvariantL{7}        = '$lex$MA-ES';

D                   = 100;
FuncNum             = 21;
if FuncNum/10<1
    owntitle = ['Problem C0' num2str(FuncNum) ' in $N=' num2str(D) '$'];
else
    owntitle = ['Problem C' num2str(FuncNum) ' in $N=' num2str(D) '$'];
end

numberOfStrategies  = 7;
numberOfRuns        = 25;


for j=1:numberOfStrategies
    k           = FuncNum;
    strategy    = ESvariant{j};
    foldername  = ['CEC17_RS_' strategy];
    load([foldername '/' strategy '_RS_on_fun' num2str(k) 'D' num2str(D) '.mat'])
    
    numberOfRuns = input.runs;
    
    for i=1:numberOfRuns
        LengthOfDyn(i)=length(Dyn{i}.gen);
    end
    
    minLength = min(LengthOfDyn);
    
    
    fit = FitT(:,11);
    cov = FitT(:,12); 
    
    fr = sum(cov==0)/length(cov);
    
    ranking =lex_sort(fit,cov);
    
    medRank = ranking((length(ranking)+1)/2);
    
%     fit(medRank)
%     cov(medRank)
    
    MedianE = Dyn{medRank}.fev;
    MedianF = Dyn{medRank}.fit;
    MedianV = Dyn{medRank}.conv;
    MedianS = Dyn{medRank}.sigma;
    
    MeanV = zeros(minLength,1);
    for i=1:numberOfRuns
%         MeanE =  MeanE + Dyn{i}.fev(1:minLength)';
        MeanV =  MeanV + Dyn{i}.conv(1:minLength)';
    end
    MeanV = MeanV./numberOfRuns;

    A{j} = [MedianE', MedianF', MedianV', MedianS']; %, MeanR, MeanN, MeanI
    
    a = min(find(MeanV==0));
    if ~isempty(a)
        minVind(j) = a;
    else
        minVind(j) = -1;
    end
    
    B{j} = [fit, cov, fit+cov];
    C{j} = fr;
end



%% plot preparations
DispN = 'DisplayName';
LineW = 'LineWidth';
LineC = 'Color';
col = [[0    0.4510    0.7412]; 
        [1 0 0]; 
        [ 0.4667    1.0000         0]; 
        [0.9294    0.6902    0.1294]; 
        [0.7176    0.2745    1.0000];
        [0.0588    1   1];
        [0.4314    0.1333    0.1882]];
    
Mlist =['o';
    'd';
    's';
    '>';
    '<';
    'v';
    '^']
    

    %% plot aggregated fitness and constraint violation dynamics   
%
figure('visible', 'on','position',[0, 0, 600, 600]);
set(gcf,'color','w')
for j= numberOfStrategies:-1:1
    strategy    = ESvariant{j};
    strategyL   = ESvariantL{j};
    cc          = col(j,:);
%     eval(['p' num2str(j) '=loglog(A{' num2str(j) '}(:,1),A{' num2str(j) '}(:,2)+A{' num2str(j) '}(:,3),LineC,cc,DispN,strategyL,LineW,1.5)']);
    eval(['p' num2str(j) '=loglog(A{' num2str(j) '}(:,1),A{' num2str(j) '}(:,2)+A{' num2str(j) '}(:,3),LineC,cc,LineW,1.5)']);
    hold on
    
    LH(j) = plot(nan, nan,['-' Mlist(j)],LineC,cc,'MarkerFaceColor',cc,'MarkerSize',8);
    L{j} = strategyL;

    
    if minVind(j) ~=-1
        plot(A{j}(minVind(j),1),A{j}(minVind(j),2)+A{j}(minVind(j),3),Mlist(j),'MarkerEdgeColor','k','MarkerFaceColor',cc,'MarkerSize',12,'HandleVisibility','off')
    end
end
 
title(owntitle)
xlabel('number of function evaluations')
ylabel('fitness $f(y)$ + violation $\nu(y)$')
% axis equal 
xlim([4*D 2*10^4*D])
% ylim([10^-30 10^5])
% ylim([-1 5])
grid on
grid minor
grid minor
lgd=legend(LH, L,'Location','SouthWest','FontSize',14);
legend boxoff
lgd.Position(4)=lgd.Position(4)*1.2;
set(gca,'FontSize',18)
if D==10
    xticks([ 1*10^1*D 1*10^2*D 1*10^3*D 1*10^4*D])
    xticklabels({'$10^2$','$10^3$','$10^4$','$10^5$'})
else
    xticks([ 1*10^1*D 1*10^2*D 1*10^3*D 1*10^4*D])
    xticklabels({'$10^3$','$10^4$','$10^5$','$10^6$'})
end

 

print(gcf,['MedianFitconv-C' num2str(FuncNum) 'N' num2str(D) '.pdf'],'-r300','-dpdf')
 

%% plot mutation strength dynamics
%
figure('visible', 'on','position',[0, 0, 600, 600]);
set(gcf,'color','w')
for j= numberOfStrategies:-1:1
    strategy    = ESvariant{j};
    strategyL   = ESvariantL{j};
    cc          = col(j,:);
    eval(['p' num2str(j) '=loglog(A{' num2str(j) '}(:,1),A{' num2str(j) '}(:,4),LineC,cc,LineW,1.5)']);
    hold on
    LH(j) = plot(nan, nan,['-' Mlist(j)],LineC,cc,'MarkerFaceColor',cc,'MarkerSize',8);
    L{j} = strategyL;
    if minVind(j) ~=-1
        plot(A{j}(minVind(j),1),A{j}(minVind(j),4),Mlist(j),'MarkerEdgeColor','k','MarkerFaceColor',cc,'MarkerSize',12,'HandleVisibility','off')
    end
end
title(owntitle)
xlabel('number of function evaluations')
ylabel('mutation strength $\sigma$')
% axis equal
xlim([4*D 2*10^4*D])
% ylim([10^-8 10^8])
grid on
grid minor
grid minor
lgd2=legend(LH, L,'Location','SouthWest','FontSize',14);
legend boxoff
lgd2.Position(4)=lgd2.Position(4)*1.2;
set(gca,'FontSize',18)
xticks([ 1*10^1*D 1*10^2*D 1*10^3*D 1*10^4*D])
xticklabels({'$10^2$','$10^3$','$10^4$','$10^5$'})
 print(gcf,['MedianSigma-C' num2str(FuncNum) 'N' num2str(D) '.pdf'],'-r300','-dpdf')

% 
% 
% 

% figure('visible', 'on','position',[0, 0, 1200, 600]);
% set(gcf,'color','w')
% b=boxplot([B{1}(:,3) B{2}(:,3) B{3}(:,3) B{4}(:,3) B{5}(:,3) B{6}(:,3) B{7}(:,3)],...%'notch','on',...
%         'labels',{ESvariantL{1},ESvariantL{2},ESvariantL{3},ESvariantL{4},ESvariantL{5},ESvariantL{6},ESvariantL{7}})
% set(b,'LineWidth',2)    
% % ylim([0 10^4])
% bp = gca;
% bp.XAxis.TickLabelInterpreter = 'latex';
% hold on
% yyaxis right
% plot([1:numberOfStrategies],[C{1} C{2} C{3} C{4} C{5} C{6} C{7}],'p','Color',[0 0.4470 0.7410],'MarkerSize',12,'LineWidth',1,'MarkerFaceColor','[0 0.4470 0.7410]')
% ylim([0 1])
% ylabel('feasibilty rate $FR$')
% yyaxis left
% % set(gca, 'YScale', 'log')
% grid on
% XTickMode = 'manual'
% XTickLabelMode = 'manual' 
% title(owntitle)
% xlabel('algorithm variants')
% ylabel('fitness $f(y)$ + violation $\nu(y)$')
% set(gca,'FontSize',14)
% set(gcf,'PaperOrientation','landscape')
% print(gcf,['BoxPlotAlg-C' num2str(FuncNum) 'N' num2str(D) '.pdf'],'-r300','-dpdf','-fillpage')
% 
% 

