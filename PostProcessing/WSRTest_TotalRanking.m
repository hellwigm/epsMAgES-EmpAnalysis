%% statistical testing for significant differences from the original epsMAg-ES
%
clear all, clc

%% set baseline algorithm 
%
ESvariant{1}        = 'epsMAgES';

%% set competitor algorithms
%
ESvariant{2}        = 'epsMAES';
ESvariant{3}        = 'epsMAgESwo';
ESvariant{4}        = 'epsMAgESnl';
ESvariant{5}        = 'epsSAgES';
ESvariant{6}        = 'lexMAgES';
ESvariant{7}        = 'lexMAES';

% dimension
D                   = 10;
owntitle = ['WilcoxonSignedRankTestN' num2str(D)];

numberOfStrategies  = 7;
numberOfFuns        = 28;

fevalsUntilGlobalBest=zeros(numberOfFuns,numberOfStrategies);

%% gather problems specific statistics
%
for j=1:numberOfStrategies
    if j==1
        strategy    = ESvariant{j};
        foldername  = ['CEC17_RS_' strategy];
        for k=1:numberOfFuns
            load([foldername '/' strategy '_RS_on_fun' num2str(k) 'D' num2str(D) '.mat'],'FitT','GB','input','Dyn')
            BL{k}=[FitT(:,11),FitT(:,12)];
            meanGB = 0;
            for i=1:input.runs
                meanGB=meanGB+GB{i}.evals;
            end
            fevalsUntilGlobalBest(k,j)=meanGB/input.runs;
        end
    else
        strategy    = ESvariant{j};
        foldername  = ['CEC17_RS_' strategy];
        for k=1:numberOfFuns
            load([foldername '/' strategy '_RS_on_fun' num2str(k) 'D' num2str(D) '.mat'],'FitT','GB','input','Dyn')
            eval(['S' num2str(j) '{' num2str(k) '}=[FitT(:,11),FitT(:,12)];']);
             meanGB = 0;
            for i=1:input.runs
                meanGB=meanGB+GB{i}.evals;
            end
            fevalsUntilGlobalBest(k,j)=meanGB/input.runs;
        end
    end
end

%% test for statistical differences in the final algorithm realizations over all 28 functions
%
%  using the Wilcoxon Signed Rank Test at a 5% significance level 
 
FSF=[1];%1:numberOfFuns

meanRankBvsS = zeros(numberOfFuns,numberOfStrategies-1);
medianRankBvsS = zeros(numberOfFuns,numberOfStrategies-1);

for kk=1:numberOfFuns
    k = FSF(kk);
    brank      = lex_sort(BL{k}(:,1),BL{k}(:,2));
    B          = BL{k};
    sBL        = BL{k}(brank,:);
    len        = size(sBL,1);
    medBL      = sBL((len+1)/2,:);
    frBL       = sum(sBL(:,2)==0)/len;
    meaBL      = mean(sBL,1);
    
    for j=2:numberOfStrategies
        eval(['srank = lex_sort(S' num2str(j) '{' num2str(k) '}(:,1),'...
                'S' num2str(j) '{' num2str(k) '}(:,2));']);
        eval(['S = S' num2str(j) '{' num2str(k) '};']);
        eval(['sS = S' num2str(j) '{' num2str(k) '}(srank,:);']);
        len       = size(sS,1);
        medS      = sS((len+1)/2,:);
        frS       = sum(sS(:,2)==0)/len;
        meaS      = mean(sS,1);
   
        %% mean ranking according to CEC + significance test
        %
        % set t1 and t2 to zero to omit the Wilcoxon Signed Rank test and
        % to presume with the CEC2017 recommended ranking
        t1 = 1;  % t1 = 0;
        t2 = 1;  % t2 = 0;

        if t1 == 1
            [p,h1] = signrank(B(1:len,1),S(1:len,1));          
        else
            h1=1; % for t1=0 significant differences in the fitness values are always assumed
        end
        
        if t2 == 1
            [p,h2] = signrank(B(1:len,2),S(1:len,2));
        else 
            h2=1; % for t2=0 significant differences in the constraint violation values are always assumed
        end
        
        
        if frS == frBL
            if frS == 1
                if h1 == 0
                    meanRankBvsS(kk,j-1) = 0;
                else
                    if meaBL(1) < meaS(1)  
                        meanRankBvsS(kk,j-1) = 1;
                    elseif meaS(1) < meaBL(1)  
                        meanRankBvsS(kk,j-1) = -1;
                    else 
                        meanRankBvsS(kk,j-1) = 0;
                    end
                end
            elseif frS < 1
                if h2 == 0
                    if h1 == 0
                           meanRankBvsS(kk,j-1) = 0;
                    else
                        if meaBL(1) < meaS(1)
                            meanRankBvsS(kk,j-1) = 1;
                        elseif meaS(1) < meaBL(1)
                            meanRankBvsS(kk,j-1) = -1;
                        else 
                            meanRankBvsS(kk,j-1) = 0;
                        end
                    end
                else
                    if meaBL(2) < meaS(2)  
                        meanRankBvsS(kk,j-1) = 1;
                    elseif meaS(2) < meaBL(2)  
                        meanRankBvsS(kk,j-1) = -1;
                    else 
                        if h1 == 0
                            meanRankBvsS(kk,j-1) = 0;
                        else
                            if meaBL(1) < meaS(1)   
                                meanRankBvsS(kk,j-1) = 1;
                            elseif meaS(1) < meaBL(1)   
                                meanRankBvsS(kk,j-1) = -1;
                            else 
                                meanRankBvsS(kk,j-1) = 0;
                            end
                        end
                    end
                end
            end
        elseif frS < frBL
            meanRankBvsS(kk,j-1) = 1;
        else
            meanRankBvsS(kk,j-1) = -1;
        end
        
        %% median based ranking according to CEC + significance test
        %
        % set t1 and t2 to zero to omit the Wilcoxon Signed Rank test and
        % to presume with the CEC2017 recommended ranking
        t1 = 1;  % t1 = 0;
        t2 = 1;  % t2 = 0;
        
        if t1 == 1
            [p,h1] = signrank(B(1:len,1),S(1:len,1));
        else
            h1=1;
        end
        
        if t2 == 1
            [p,h2] = signrank(B(1:len,2),S(1:len,2));
        else 
            h2=1;
        end
        
        
        if h2 == 0
            if h1 == 0
                medianRankBvsS(kk,j-1) = 0;
            else
                if medBL(1) < medS(1)
                        medianRankBvsS(kk,j-1) = 1;
                elseif medS(1) < medBL(1)
                        medianRankBvsS(kk,j-1) = -1;
                else 
                        medianRankBvsS(kk,j-1) = 0;
                end
            end
        else
            if medBL(2) < medS(2)
                medianRankBvsS(kk,j-1) = 1;
            elseif medS(2) < medBL(2)
                medianRankBvsS(kk,j-1) = -1;
            else
                if h1 == 0
                    medianRankBvsS(kk,j-1) = 0;
                else
                    if medBL(1) < medS(1)
                            medianRankBvsS(kk,j-1) = 1;
                    elseif medS(1) < medBL(1)
                            medianRankBvsS(kk,j-1) = -1;
                    else 
                            medianRankBvsS(kk,j-1) = 0;
                    end
                end
            end
        end
        
    end
end

% compute final (total ranks) by aggregating median and mean rankings
totalRankBvsS = medianRankBvsS+meanRankBvsS;
totalWTL = [sum(totalRankBvsS>0);sum(totalRankBvsS==0);sum(totalRankBvsS<0)]; % counts

% provide the relative running times in terms of functin evaluations
% consumed until best observed candidate solution was found
for i = 2:7
    for j=1:numberOfFuns
        if totalRankBvsS(j,i-1)==0
            TotalList(j,i-1) = string(num2str(fevalsUntilGlobalBest(j,i)./fevalsUntilGlobalBest(j,1),'%1.4f'));
        elseif totalRankBvsS(j,i-1)>0
            TotalList(j,i-1) = string(['+']);
        elseif totalRankBvsS(j,i-1)<0
            TotalList(j,i-1) = string(['-']);
        else
            disp('some unexpected issue! --> check code') % usually unnecessary
        end
    end
end
    
filePh = fopen(['totalRanks' num2str(D) '.txt'],'w');
fprintf(filePh,'%s, %s, %s, %s, %s, %s\n',TotalList');
fclose(filePh);   