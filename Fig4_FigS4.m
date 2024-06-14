load('functional_correlation_pre_stim.mat')

pcStim = nan(81,76);
pcpreStim = nan(81,76);

for tt = 1:81
    [tt 81]
      
   
  
    for ss = 1:2
        
        %load in fc matrix
        if ss == 1
            cc=pre_stim_CC{tt,ss};
        elseif ss == 2
            cc=pre_stim_CC{tt,ss};
           
        end
        
        %community_louvain hyper parameter
gamma = 1.3;

%remove all self correlations
        cc = cc - diag(diag(cc));

        %PC is a regionwise measure
        pcAll = nan(76,100);
        
        for nn = 1:100
            
        [Ci,~] = community_louvain(cc,gamma,[],'negative_asym');
        [Ppos,~]=participation_coef_sign(cc,Ci);
        
        pcAll(:,nn) = (Ppos);
        end

         PcROI = mean(pcAll,2);

         if ss == 1
             pcpreStim(tt,:) = PcROI;
         elseif ss == 2
             pcStim(tt,:) = PcROI;


         end

    end
    
end

%% Plot boxplots of Fig. 4H

figure
pc3HzStim = pcStim(hz3ID,:) - mean(pcpreStim(hz3ID,:),2);
pc3HzStim = mean(pc3HzStim,2);

pc5HzStim = pcStim(hz5ID,:) - mean(pcpreStim(hz5ID,:),2);
pc5HzStim = mean(pc5HzStim,2);

pc15HzStim = pcStim(hz15ID,:) - mean(pcpreStim(hz15ID,:),2);
pc15HzStim = mean(pc15HzStim,2);

pcShamStim = pcStim(shamID,:) - mean(pcpreStim(shamID,:),2);
pcShamStim = mean(pcShamStim,2);

%

dat = [pcShamStim(:); pc3HzStim(:); pc5HzStim(:); pc15HzStim(:)];

gID = [ones(numel(pcShamStim),1); 2.*ones(numel(pc3HzStim),1); 3.*ones(numel(pc5HzStim),1); 4.*ones(numel(pc15HzStim),1)]; 

boxchart(gID,dat)

%boxchart(gID,dat,'GroupByColor',gID)
hold on

swarmchart(gID,dat,20,'filled')

yline(0)
%ylim([-0.2 0.25])

handle = gca;
handle.TickDir = 'out'; 
handle.TickLength = [.01 .01];
handle.Box = 'off';
handle.Color = 'none';
handle.FontSize = 6.24;

pbaspect([2.4 1 1])

%% Statistical analysis for Fig. 4
[p,tbl,stats] = anova1([dat], gID,'off');
[c,~,~]  = multcompare(stats);

%% Supplementary analysis Fig. S4

%Cortical heirarchy score
heirID = [32	7	3	25	26	4	30	22	6	31	5	8	24	2	38	29	23	20	37	34	27	17	21	28	19	1	35	18	33	16	36	13	15	10	14	12	9	11];
TTweight = [1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39];

ROI_names = {'SSp-m'	'SSp-bfd'	'SSp-n'	'SSp-ul'	'SSp-ll'	'SSs'	'SSp-tr'	'AUDd'	'VISp'	'MOp'	'VISal'	'AUDp'	'PTLp'	'VISl'	'RSPv'	'VISpm'	'VISam'	'RSPd'	'MOs'	'AUDv'	'VISpl'	'VISC'	'RSPagl'	'GU'	'ACAd'	'ORBl'	'FRP'	'AId'	'PERI'	'AIp'	'TEa'	'ECT'	'ACAv'	'ORBvl'	'ORBm'	'ILA'	'AIv'	'PL'};


pc3HzStim = pcStim(hz3ID,:);
pc3HzStim = mean(pc3HzStim,1);

pc5HzStim = pcStim(hz5ID,:);
pc5HzStim = mean(pc5HzStim,1);

pc15HzStim = pcStim(hz15ID,:);
pc15HzStim = mean(pc15HzStim,1);


PCdat_3vs15 = pc3HzStim-pc15HzStim;
PCdat_3vs5 = pc3HzStim-pc5HzStim;

subplot(2,3,[1 2])

[r,p] = corr(fliplr(TTweight)',PCdat_3vs15(heirID)','Type','Spearman');

scatter(fliplr(1:numel(ROI_names)),PCdat_3vs15(heirID),20,PCdat_3vs15(heirID),"filled")

set(gca,'XTick',1:numel(ROI_names),'XTickLabel',ROI_names)
xtickangle(90)

colorbar

title(['3Hz - 15Hz r = ' num2str(round(r,3)) ' p =  ' num2str(p)] )

lsline
yline(0)
ylabel('\Delta PC')
xlabel('primary somatosensory ---> transmodal')
tplt

%fixing ticks out
handle = gca;
handle.TickDir = 'out'; 
handle.TickLength = [.01 .01];
handle.Box = 'off';
handle.Color = 'none';
handle.FontSize = 5.5;

pbaspect([2.4 1 1])

h = get(gcf, 'Position');
h(3) = h(3)*1.275;
h(4) = h(4)*1.275;
set(gcf, 'Position', h)


subplot(2,3,[4 5])

[r,p] = corr(fliplr(TTweight)',PCdat_3vs5(heirID)','Type','Spearman');

scatter(fliplr(1:numel(ROI_names)),PCdat_3vs5(heirID),20,PCdat_3vs5(heirID),"filled")

set(gca,'XTick',1:numel(ROI_names),'XTickLabel',ROI_names)
xtickangle(90)

colorbar

title(['3Hz - 5Hz r = ' num2str(round(r,3)) ' p =  ' num2str(p)] )

lsline
yline(0)
ylabel('\Delta PC')
xlabel('primary somatosensory ---> transmodal')
tplt

%fixing ticks out
handle = gca;
handle.TickDir = 'out'; 
handle.TickLength = [.01 .01];
handle.Box = 'off';
handle.Color = 'none';
handle.FontSize = 5.5;

pbaspect([2.4 1 1])

h = get(gcf, 'Position');
h(3) = h(3)*1.275;
h(4) = h(4)*1.275;
set(gcf, 'Position', h)


