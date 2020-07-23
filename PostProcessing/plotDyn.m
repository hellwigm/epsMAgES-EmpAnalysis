clear all, clc

%% build mean dynamics
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

D                   = 10;
FuncNum             = 1;
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
    
%     if input.runs < 25
%         break
%     end
    
    for i=1:numberOfRuns
        LengthOfDyn(i)=length(Dyn{i}.gen);
    end
    
    minLength = min(LengthOfDyn);
    
    MeanE = zeros(minLength,1);
    MeanF = zeros(minLength,1);
    MeanV = zeros(minLength,1);
    MeanS = zeros(minLength,1);
    for i=1:numberOfRuns
        MeanE =  MeanE + Dyn{i}.fev(1:minLength)';
        MeanF =  MeanF + Dyn{i}.fit(1:minLength)';
        MeanV =  MeanV + Dyn{i}.conv(1:minLength)';
        MeanS =  MeanS + Dyn{i}.sigma(1:minLength)';
    end
    MeanE = MeanE./numberOfRuns;
    MeanF = MeanF./numberOfRuns;
    MeanV = MeanV./numberOfRuns;
    MeanS = MeanS./numberOfRuns;

    A{j} = [MeanE, MeanF, MeanV, MeanS]; %, MeanR, MeanN, MeanI
    a = min(find(MeanV==0));
    if ~isempty(a)
        minVind(j) = a;
    else
        minVind(j) = -1;
    end
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

    %% plot aggregated fitness and constraint violation dynamics   
%
figure('visible', 'on','position',[0, 0, 600, 600]);
set(gcf,'color','w')
for j= numberOfStrategies:-1:1
    strategy    = ESvariant{j};
    strategyL   = ESvariantL{j};
    cc          = col(j,:);
    eval(['p' num2str(j) '=loglog(A{' num2str(j) '}(:,1),A{' num2str(j) '}(:,2)+A{' num2str(j) '}(:,3),LineC,cc,DispN,strategyL,LineW,1.5)']);
    hold on
    if minVind(j) ~=-1
        plot(A{j}(minVind(j),1),A{j}(minVind(j),2)+A{j}(minVind(j),3),['o'],'MarkerEdgeColor','k','MarkerFaceColor',cc,'MarkerSize',12,'HandleVisibility','off')
    end
end
title(owntitle)
xlabel('number of function evaluations')
ylabel('fitness $f(y)$ + violation $\nu(y)$')
% axis equal 
xlim([4*D 2*10^4*D])
% ylim([10^1 10^7])
% ylim([-1 5])
grid on
grid minor
grid minor
lgd=legend('show','Location','SouthWest','FontSize',14);
legend boxoff
lgd.Position(4)=lgd.Position(4)*1.2;
set(gca,'FontSize',14)
if D==10
    xticks([ 1*10^1*D 1*10^2*D 1*10^3*D 1*10^4*D])
    xticklabels({'$10^2$','$10^3$','$10^4$','$10^5$'})
else
    xticks([ 1*10^1*D 1*10^2*D 1*10^3*D 1*10^4*D])
    xticklabels({'$10^3$','$10^4$','$10^5$','$10^6$'})
end
% print(gcf,['fitconv-C' num2str(FuncNum) 'N' num2str(D) '.pdf'],'-r300','-dpdf')

%% plot mutation strength dynamics
%
figure('visible', 'on','position',[0, 0, 600, 600]);
set(gcf,'color','w')
for j= numberOfStrategies:-1:1
    strategy    = ESvariant{j};
    strategyL   = ESvariantL{j};
    cc          = col(j,:);
    eval(['p' num2str(j) '=loglog(A{' num2str(j) '}(:,1),A{' num2str(j) '}(:,4),LineC,cc,DispN,strategyL,LineW,1.5)']);
    hold on
    if minVind(j) ~=-1
        plot(A{j}(minVind(j),1),A{j}(minVind(j),4),['o'],'MarkerEdgeColor','k','MarkerFaceColor',cc,'MarkerSize',12,'HandleVisibility','off')
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
lgd2=legend('show','Location','SouthWest','FontSize',14);
legend boxoff
lgd2.Position(4)=lgd2.Position(4)*1.2;
set(gca,'FontSize',14)
xticks([ 1*10^1*D 1*10^2*D 1*10^3*D 1*10^4*D])
xticklabels({'$10^2$','$10^3$','$10^4$','$10^5$'})
% print(gcf,['sigma-C' num2str(FuncNum) 'N' num2str(D) '.pdf'],'-r300','-dpdf')




