%% load, calculate, combine

clear all

load("SynSimResults3.mat")
Ts=0.01;
StFreq=20;
DSample=1/Ts/StFreq;
ParamNum=3;
SynNum=5;
ValNum=100;
DataRange=[50 500];
DSampDataRange=[50 500]/DSample;
DataInd=DataRange(1):DataRange(2);
DSampDataInd=DSampDataRange(1):DSampDataRange(2);

TrackTar=0.5;
TrackVec=ones(ValNum,length(DataInd))*TrackTar;
DSamptrackVec=ones(ValNum,length(DSampDataInd))*TrackTar;

FileLabels={'syn_3filter2_PhandNoReps2','syn_3filter2_SHandNoReps2',...
    'syn_3filter2_OccNoReps2','syn_3filter2_mWaveDelNoReps2','syn_3filter3_occrate7'};
SampTypeLabels=["PChangeMat", "VarMat", "AvgMinMaxMat"];

SynLabels=["PHand", "SHand","Occ","DmWave","OccRate"];


for iSyn=1:SynS.SynNum
    
    FileLbl=FileLabels{iSyn};
    SynLabel=SynLabels{iSyn};
    load(FileLbl)

    for iRep=1:1

        RepLabel=sprintf('Rep_%d',iRep);

        for iFilt=1:3
            FiltLabel=sprintf('Filt_%d',iFilt);
            
            TrackingError= S.(RepLabel).(FiltLabel).TrackingError ; % rows are parameters values
            TargetVec=ones(ValNum,size(TrackingError,2))*TrackTar;
            TrackingVec=TargetVec-TrackingError;
            
            EstEffort= S.(RepLabel).(FiltLabel).EstEffort;
            TrueEffort=S.(RepLabel).(FiltLabel).TrueEffort;
            
            T.(SynLabel).RSquared.(FiltLabel).EstEffort=EstEffort;
            T.(SynLabel).RSquared.(FiltLabel).TrueEffort=TrueEffort;
            
            T.(SynLabel).RSquared.(FiltLabel).TrackingError=TrackingError;
            T.(SynLabel).RSquared.(FiltLabel).TrackingVec=TrackingVec;
            T.(SynLabel).RSquared.(FiltLabel).TargetVec=TargetVec;

            SynS.(SynLabel).TrackingError.wc=wc;


            for iVal=1:NumPoints
                ValLabel=sprintf('Val_%d',iVal);
                
                mdl = fitlm(EstEffort(iVal,DataInd),TrueEffort(iVal,DataInd));
                r_sqr=mdl.Rsquared.Ordinary;
%                 coef_p_val=mdl.Coefficients.pValue;
                AnovaTable=anova(mdl,'summary');
%                 F=table2array(AnovaTable(2,4));
                p_val=table2array(AnovaTable(2,5));
                
                T.(SynLabel).RSquared.Effdata.r_sqr(iFilt,iVal)=r_sqr;
                T.(SynLabel).RSquared.data(iFilt,iVal)=r_sqr;
                T.(SynLabel).RSquared.wc=SynS.(SynLabel).TrackingError.wc;
                
                T.(SynLabel).RSquared.Effdata.p_val(iFilt,iVal)=p_val;
%                 T.(SynLabel).RSquared.data.AnovaTable(iFilt,iVal)=AnovaTable;
                
                T.(SynLabel).RSquared.Effmdl.(FiltLabel).(ValLabel)=mdl; 
                T.(SynLabel).RSquared.Effatable.(FiltLabel).(ValLabel)=AnovaTable;   
                
                %
                mdl = fitlm(TargetVec(iVal,DataInd),TrackingVec(iVal,DataInd));
                r_sqr=mdl.Rsquared.Ordinary;
                AnovaTable=anova(mdl,'summary');
%                 p_val=table2array(AnovaTable(2,5));
                
                T.(SynLabel).RSquared.Trackdata.r_sqr(iFilt,iVal)=r_sqr;
                T.(SynLabel).RSquared.Trackdata.p_val(iFilt,iVal)=p_val;
                
                T.(SynLabel).RSquared.Trackmdl.(FiltLabel).(ValLabel)=mdl; 
                T.(SynLabel).RSquared.Trackatable.(FiltLabel).(ValLabel)=AnovaTable;   

                %

                DSampEstEffort(iVal,:)=downsample(EstEffort(iVal,:),DSample);
                DSampTrueEffort(iVal,:)=downsample(TrueEffort(iVal,:),DSample);
                
                mdl = fitlm(DSampEstEffort(iVal,DSampDataInd),DSampTrueEffort(iVal,DSampDataInd));
                r_sqr=mdl.Rsquared.Ordinary;
                coef_p_val=mdl.Coefficients.pValue;
                AnovaTable=anova(mdl,'summary');
                F=table2array(AnovaTable(2,4));
                p_val=table2array(AnovaTable(2,5));

                T.(SynLabel).RSquared.DSamp.Effdata.r_sqr(iFilt,iVal)=r_sqr;
                T.(SynLabel).RSquared.DSamp.Effdata.p_val(iFilt,iVal)=p_val;
%                 T.(SynLabel).RSquared.data.AnovaTable(iFilt,iVal)=AnovaTable;
                
                T.(SynLabel).RSquared.DSamp.Effmdl.(FiltLabel).(ValLabel)=mdl; 
                T.(SynLabel).RSquared.DSamp.Effatable.(FiltLabel).(ValLabel)=AnovaTable;
                
                DSampTrackingVec(iVal,:)=downsample(TrackingVec(iVal,:),DSample);
                DSampTargetVec(iVal,:)=downsample(TargetVec(iVal,:),DSample);
                
                mdl = fitlm(DSampTrackingVec(iVal,DSampDataInd),DSampTargetVec(iVal,DSampDataInd));
                r_sqr=mdl.Rsquared.Ordinary;
                coef_p_val=mdl.Coefficients.pValue;
                AnovaTable=anova(mdl,'summary');
                F=table2array(AnovaTable(2,4));
                p_val=table2array(AnovaTable(2,5));
                
                T.(SynLabel).RSquared.DSamp.Trackdata.r_sqr(iFilt,iVal)=r_sqr;
                T.(SynLabel).RSquared.DSamp.Trackdata.p_val(iFilt,iVal)=p_val;
%                 T.(SynLabel).RSquared.data.AnovaTable(iFilt,iVal)=AnovaTable;
                
                T.(SynLabel).RSquared.DSamp.Trackmdl.(FiltLabel).(ValLabel)=mdl; 
                T.(SynLabel).RSquared.DSamp.Trackatable.(FiltLabel).(ValLabel)=AnovaTable;
                
                T.(SynLabel).RSquared.(FiltLabel).EstEffort=EstEffort;
                T.(SynLabel).RSquared.(FiltLabel).TrueEffort=TrueEffort;
                T.(SynLabel).RSquared.(FiltLabel).DSamp.DSampEstEffort=DSampEstEffort;
                T.(SynLabel).RSquared.(FiltLabel).DSamp.DSampTrueEffort=DSampTrueEffort;
                T.(SynLabel).RSquared.(FiltLabel).DSamp.DSampTrackingVec=DSampTrackingVec;
                T.(SynLabel).RSquared.(FiltLabel).DSamp.DSampTargetVec=DSampTargetVec;
                
            end
        end
    end
    
    SynS.(SynLabel).RSquared=T.(SynLabel).RSquared;
    clear T
end


SynS.Ts=Ts;
SynS.DSample=DSample;
SynS.StFreq=StFreq;
SynS.DataRange=DataRange;
SynS.DataInd=DataInd;
SynS.TrackTar=TrackTar;
SynS.SynLabels=SynLabels;
 

save('SynSimResults3_rsqr','SynS')

%% plot rsquared
clear all
load('SynSimResults2')
%%
% BaseLine=50;
iRep=1;
RepLabel=sprintf('Rep_%d',iRep);
fit_n=4;

for iSyn=1:5
    SynLabel=SynS.SynLabels(iSyn);

    r_sqr=SynS.(SynLabel).RSquared.Effdata.r_sqr;
    p_val=SynS.(SynLabel).RSquared.Effdata.p_val;

    if iSyn==5
        for FilterNum=1:3
            Fit_r_sqr(FilterNum,:)=nPolyFit(1:length(r_sqr),r_sqr(FilterNum,:),fit_n);
        end
        r_sqr=Fit_r_sqr;
    end

    figure (1)
    subplot(2,5,iSyn)
    plot(r_sqr')
    xlabel(SynS.SynLabels(iSyn))
    ylabel('R-Squared')
    
    subplot(2,5,iSyn+5)
    plot(p_val')
    title(SynS.SynLabels(iSyn))
    xlabel(SynS.SynLabels(iSyn))
    ylabel('p-value')
    
    hold 
    plot([1 100],[0.05 0.05],'k')
    ylim([0 0.1])
    
%     r_sqr=SynS.(SynLabel).RSquared.Trackdata.r_sqr;
%     p_val=SynS.(SynLabel).RSquared.Trackdata.p_val;
% 
%     
%     figure (2)
%     subplot(2,5,iSyn)
%     plot(r_sqr')
%     xlabel(SynS.SynLabels(iSyn))
%     ylabel('R-Squared')
%     
%     subplot(2,5,iSyn+5)
%     plot(p_val')
%     
%     hold 
%     plot([1 100],[0.05 0.05],'k')
%     xlabel(SynS.SynLabels(iSyn))
%     ylabel('p-value')
    
        
%     ylim([0 0.1])
%     legend(["Comb", "GS", "Blanking" ])
    r_sqr=SynS.(SynLabel).RSquared.DSamp.Effdata.r_sqr;
    p_val=SynS.(SynLabel).RSquared.DSamp.Effdata.p_val;
    if iSyn==5
        for FilterNum=1:3
            Fit_r_sqr(FilterNum,:)=nPolyFit(1:length(r_sqr),r_sqr(FilterNum,:),fit_n);
        end
        r_sqr=Fit_r_sqr;
    end
    figure (3)
    subplot(2,5,iSyn)
    plot(r_sqr')
    xlabel(SynS.SynLabels(iSyn))
    ylabel('R-Squared')
    subplot(2,5,iSyn+5)
    plot(p_val')
    hold 
    plot([1 100],[0.05 0.05],'k')
    ylim([0 0.1])
    xlabel(SynS.SynLabels(iSyn))
    ylabel('p-value')
    
% 
%     r_sqr=SynS.(SynLabel).RSquared.DSamp.Trackdata.r_sqr;
%     p_val=SynS.(SynLabel).RSquared.DSamp.Trackdata.p_val;
% 
%     figure (4)
%     subplot(2,5,iSyn)
%     plot(r_sqr')
%     subplot(2,5,iSyn+5)
%     plot(p_val')
%     hold 
%     plot([1 100],[0.05 0.05],'k')
%     ylim([0 0.1])
%     legend(["Comb", "GS", "Blanking" ])


end
% figure(1)
% subplot(2,5,1)
% legend(["Comb", "GS", "Blanking" ])
% figure(3)
% subplot(2,5,1)
% legend(["Comb", "GS", "Blanking" ])
% figure(2)
% subplot(2,5,1)
% legend(["Comb", "GS", "Blanking" ])
% figure(4)
% subplot(2,5,1)
% legend(["Comb", "GS", "Blanking" ])

%%
clear EstEffort TrueEffort DSampEstEffort DSampTrueEffort

iVal=50;
for iSyn=1:5
    SynLabel=SynS.SynLabels(iSyn);

    for iFilt=1:3
        FiltLabel=sprintf('Filt_%d',iFilt);

        EstEffort(iFilt,:)=SynS.(SynLabel).RSquared.(FiltLabel).EstEffort(iVal,:);
        TrueEffort(iFilt,:)=SynS.(SynLabel).RSquared.(FiltLabel).TrueEffort(iVal,:);
        DSampEstEffort(iFilt,:)=SynS.(SynLabel).RSquared.(FiltLabel).DSamp.DSampEstEffort(iVal,:);
        DSampTrueEffort(iFilt,:)=SynS.(SynLabel).RSquared.(FiltLabel).DSamp.DSampTrueEffort(iVal,:);
        
        for iVal=50:50
            ValLabel=sprintf('Val_%d',iVal);
            
            DSamp_mdl=SynS.(SynLabel).RSquared.DSampMdl.(FiltLabel).(ValLabel);
            mdl=SynS.(SynLabel).RSquared.mdl.(FiltLabel).(ValLabel);
            
            figure (1)
            subplot(2,5,iSyn)
            plot(mdl)
            subplot(2,5,iSyn+5)
            plot(DSamp_mdl)

            
        end

    end
    
    
        figure (2)
        subplot(2,5,iSyn)
        plot(TrueEffort',EstEffort')
        subplot(2,5,iSyn+5)
        plot(DSampTrueEffort',DSampEstEffort')
        
        
        
        
end


%%
BaseLine=50;
Target=0.5;
fit_n=8;
ParamNumber=3;
ParamLabels={'TrackingError','TrackingSNR','EstEffSNR','EffortCorr'};
ParamTitles={'TrackingError','TrackingSNR',...
    'Estimated Effort SNR','Correlation: True vs Estimated Effort'};
xlabels={'Volitional Paretic Hand Opening Constant' };
ylabels={'Tracking Error (a.u.)','Tracking SNR',...
    'Estimtad Effort SNR','Pearson Correlation'};
FilterLabels={'Comb','GS','Blanking'};
LineLabels={'-','--','-.'};

for iParam=1:ParamNumber

    ParamLabel=ParamLabels{iParam};
    
    if iParam==4
        fit_n=0;
    end
    
    for FilterNum=1:3
        y(FilterNum,:)=O.(ParamLabel).data(FilterNum,:);
        x=O.wc.^-1;
        PercentChange(FilterNum,:)=-100+100*abs(y(FilterNum,:))/abs(y(1,BaseLine));
        FittedVal(FilterNum,:)=nPolyFit(x,y(FilterNum,:),fit_n);
        FittedPercent(FilterNum,:)=nPolyFit(x,PercentChange(FilterNum,:),fit_n);
    end
    
    figure(1)
    subplot(ParamNumber,2,iParam*2-1)
%     yyaxis left
    for FilterNum=1:3
%         p=semilogx(x,y(FilterNum,:),LineLabels{FilterNum},'LineWidth',1,'Color','k');
%        p.Color(4) = 0.5;
       pright(iParam,FilterNum)=semilogx(x,FittedVal(FilterNum,:),LineLabels{FilterNum},'LineWidth',2,'Color','k');
        hold on

    end
    
    a = get(gca,'Children');
    y1data = get(a, 'YData');
    y1min=min( [min(y1data{1}) min(y1data{2}) min(y1data{3}) ]);
    y1max=max( [max(y1data{1}) max(y1data{2}) max(y1data{3}) ]);
    scaledif=abs(y1max-y1min);
    ylim([y1min-scaledif*0.2 y1max+scaledif*0.2]);
    ylabel(ylabels{iParam});
    ax = gca;
    ax.YColor = 'k';
    
%     figure(1)
%     yyaxis right
%     nmin=1;
%     nmax=3;
% 
%     for FilterNum=nmin:nmax
%     %     semilogx(x,PercentChange(FilterNum,:),'-','LineWidth',2,'Color','k')
%         pleft(iParam,FilterNum)=semilogx(x,FittedPercent(FilterNum,:),LineLabels{FilterNum},'LineWidth',2,'Color','k');
%     end
%     
%     b = get(gca,'Children');
%     y2data = get(b, 'YData');
% 
%     y2min=min( [min(y2data{1}) min(y2data{2}) min(y2data{3}) ]);
%     y2max=max( [max(y2data{1}) max(y2data{2}) max(y2data{3}) ]);
%     scaledif=abs(y2max-y2min);
%     ylim([y2min-scaledif*0.2 y2max+scaledif*0.2]);
%     title(ParamTitles{iParam});
%     ax = gca;
%     ax.YColor = 'k';
%     
    if iParam==1
        legend(FilterLabels,'AutoUpdate','off');
    end
    
    for FilterNum=1:3
        plot(x(BaseLine),FittedVal(FilterNum,BaseLine),'*','Color','k')
    end
    
end



subplot(ParamNumber,2,5)
% yyaxis right
% lbl=ylabel('Change from Nominal Value (% change)');
% lbl.Position(2) = 50; 
xlabel('Volitional Paretic Hand Opening Constant')

% delete(pleft(4,1))