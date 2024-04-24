%% Combine Files
clear all
FileLabels={'work_paretic_noseed2', 'work_stim_noseed2',...
    'work_occ_noseed2','work_mwave_noseed3','syn_3filter3_occrate7'};
SynLabels=["PHand", "SHand","Occ","DmWave","OccRate"];
ParamLabels=["TrackingError","TrackingSNR","EstEffSNR","RSquared"];
FilterLabels=["Comb","GS","Blanking"];
SampTypeLabels=["PChangeMat", "VarMat", "AvgMinMaxMat"];

ParamNum=3;
SynS=struct;
SynNum=5;
iRep=1;
RepLabel=sprintf('Rep_%d',iRep);
for iSyn=1:SynNum
    
    FileLbl=FileLabels{iSyn};
    SynLabel=SynLabels{iSyn};

    load(FileLbl)
    clear M

                
    if iSyn==5 

        M.TrackingError.data=S.(RepLabel).TrackingError;
        M.StdTe=S.(RepLabel).StdTe;
        M.EstEffort=S.(RepLabel).EstEffort;
        M.TrueEffort=S.(RepLabel).TrueEffort;
        M.NonPareticAngle=S.(RepLabel).NonPareticAngle;
        M.StdEE=S.(RepLabel).StdEE;
        M.EffortCorr.data=S.(RepLabel).EffortCorr;
        M.TrackingSNR.data=(S.(RepLabel).StdTe/Target).^-1;
        M.EstEffSNR.data=(S.(RepLabel).EstEffort./S.(RepLabel).StdEE);
        M.wc=wc;
    
    else
        M=struct;
        M.TrackingError.data=TrackingError;
        M.StdTe=StdTe;
        M.EstEffort=EstEffort;
        M.TrueEffort=TrueEffort;
        M.NonPareticAngle=NonPareticAngle;
        M.StdEE=StdEE;
        M.EffortCorr.data=EffortCorr;
        M.TrackingSNR.data=(StdTe/Target).^-1;
        M.EstEffSNR.data=(EstEffort./StdEE);
        M.wc=wc;
    end
    

    
    if iSyn==4
        for iParam=1:ParamNum

            ParamLabel=ParamLabels{iParam};

            for FilterNum=1:3   

                y(FilterNum,:)=M.(ParamLabel).data(FilterNum,:);
                newy(FilterNum,:)= y(FilterNum,:)/y(FilterNum,50)*y(FilterNum,10);
                M.(ParamLabel).data(FilterNum,:)=newy(FilterNum,:);
                
            end
        end
    end

    for iParam=1:ParamNum

        ParamLabel=ParamLabels{iParam};

        for iFilt=1:3   
            
            y(iFilt,:)=M.(ParamLabel).data(iFilt,:);
            x=M.wc;
            SynS.(SynLabel).(ParamLabel).data(iFilt,:)=y(iFilt,:);
            SynS.(SynLabel).(ParamLabel).wc=x;

        end
    end

end


SynS.FileLabels=FileLabels;
SynS.SynLabels=SynLabels;
SynS.ParamLabels=ParamLabels;
SynS.ParamNum=ParamNum;
SynS.SynNum=SynNum;
SynS.FilterLabels=FilterLabels;
SynS.SampTypeLabels=SampTypeLabels;


% save('SynSimResults3','SynS')

%%Calc rsqr

% clear all
% load("SynSimResults3.mat")


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
 
%
% save('SynSimResults3_rsqr','SynS')

%%SynStats- calculate the metrics

% clear all
% load('SynSimResults3_rsqr');
ParamLabels=["TrackingError","TrackingSNR","EstEffSNR","RSquared"];
SynS.ParamLabels=ParamLabels;

fit_n=8;
BaseLineVal=[ 50 50 50 50 50];
SynS.BaseLineVal=BaseLineVal;
SynS.fit_n=fit_n;
SynNum=SynS.SynNum;
FilterLabels=SynS.FilterLabels;
PercentVar=[ 20 30 40 35 35]/100;
SynS.PercentVar=PercentVar;
GenSampleSize=10;
SynS.GenSampleSize=GenSampleSize;
ParamNum=4;
SynS.ParamNum=ParamNum;


for iSyn=1:SynNum
    
    SynLabel=SynS.SynLabels{iSyn};

    for iParam=1:ParamNum

        ParamLabel=ParamLabels(iParam);

        for iFilt=1:3

            % calculate fitting
            y(iFilt,:)=SynS.(SynLabel).(ParamLabel).data(iFilt,:);
            x=SynS.(SynLabel).(ParamLabel).wc;
            PercentChange(iFilt,:)=-100+100*abs(y(iFilt,:))/abs(y(1,BaseLineVal(iSyn)));
            FittedVal(iFilt,:)=nPolyFit(x,y(iFilt,:),fit_n);
            FittedPercent(iFilt,:)=nPolyFit(x,PercentChange(iFilt,:),fit_n);
            
            Pmin=FittedVal(iFilt,10)/FittedVal(iFilt,BaseLineVal(iSyn))*100-100;
            Pmax=FittedVal(iFilt,90)/FittedVal(iFilt,BaseLineVal(iSyn))*100-100;
            
            % calculate the change metrics
            m=(abs([Pmin-Pmax]));
            v=randn(1,GenSampleSize)*m*PercentVar(iParam)+m;
            PChangeMat(iParam).(SynLabel).PercentChange(iFilt,:)=m;
            PChangeMat(iParam).(SynLabel).GeneratedValues(iFilt,:)=v;
            
            m=var(abs(FittedPercent(iFilt,:)));
            v=randn(1,GenSampleSize)*m*PercentVar(iParam)+m;
            VarMat(iParam).(SynLabel).PercentChangeVar(iFilt,:)=m;
            VarMat(iParam).(SynLabel).GeneratedValues(iFilt,:)=v;
            
            m=abs(min(FittedPercent(iFilt,:))-max(FittedPercent(iFilt,:)));
            v=randn(1,GenSampleSize)*m*PercentVar(iParam)+m;
            AvgMinMaxMat(iParam).(SynLabel).AvgMinMax(iFilt,:)=m;
            AvgMinMaxMat(iParam).(SynLabel).GeneratedValues(iFilt,:)=v;
            
            %Baseline values
            m=FittedVal(iFilt,BaseLineVal(iSyn));
            v=randn(1,GenSampleSize)*m*PercentVar(iParam)+m;
            BaseValMat(iParam).(SynLabel).BaseLineValue(iFilt,:)=m; 
            BaseValMat(iParam).(SynLabel).GeneratedValues(iFilt,:)=v;           
            
            % saving to struct
%             SynS.(SynLabel).(ParamLabel).PercentChange(iFilt,:)=PercentChange(iFilt,:);
%             SynS.(SynLabel).(ParamLabel).FittedVal(iFilt,:)=FittedVal(iFilt,:);
%             SynS.(SynLabel).(ParamLabel).FittedPercent(iFilt,:)=FittedPercent(iFilt,:);
%             SynS.(SynLabel).(ParamLabel).PMinMax=[ Pmin Pmax];

        end
    end
end
SynS.PChangeMat=PChangeMat; 
SynS.VarMat=VarMat; 
SynS.AvgMinMaxMat=AvgMinMaxMat; 
SynS.BaseValMat=BaseValMat; 

% range of change
% load('SynStatsVec3')

ParamLabels=["TrackingError","TrackingSNR","EstEffSNR","RSquared"];
SynS.ParamLabels=ParamLabels;

fit_n=8;
BaseLineVal=[ 50 50 50 50 50];
SynS.BaseLineVal=BaseLineVal;
SynS.fit_n=fit_n;
SynNum=SynS.SynNum;
FilterLabels=SynS.FilterLabels;
PercentVar=[ 20 30 40 35 35]/100;
SynS.PercentVar=PercentVar;
GenSampleSize=10;
SynS.GenSampleSize=GenSampleSize;
ParamNum=4;
SynS.ParamNum=ParamNum;

for iSyn=1:SynNum
    
    SynLabel=SynS.SynLabels{iSyn};

    for iParam=1:ParamNum

        ParamLabel=ParamLabels(iParam);

        for iFilt=1:3

            % calculate fitting
            y(iFilt,:)=SynS.(SynLabel).(ParamLabel).data(iFilt,:);
            x=SynS.(SynLabel).(ParamLabel).wc;
            PercentChange(iFilt,:)=-100+100*abs(y(iFilt,:))/abs(y(1,BaseLineVal(iSyn)));
            FittedVal(iFilt,:)=nPolyFit(x,y(iFilt,:),fit_n);
            FittedPercent(iFilt,:)=nPolyFit(x,PercentChange(iFilt,:),fit_n);
            
            Pmin=FittedVal(iFilt,10)/FittedVal(iFilt,BaseLineVal(iSyn))*100-100;
            Pmax=FittedVal(iFilt,90)/FittedVal(iFilt,BaseLineVal(iSyn))*100-100;
            
            % calculate the change metrics and generate
            m=abs(Pmin-BaseLineVal(iSyn)) + abs(Pmax-BaseLineVal(iSyn));
            v=randn(1,GenSampleSize)*m*PercentVar(iParam)+m;
            RangeChangeMat(iParam).(SynLabel).RangeChange(iFilt,:)=m;
            RangeChangeMat(iParam).(SynLabel).GeneratedValues(iFilt,:)=v;  

        end
    end
end
SynS.RangeChangeMat=RangeChangeMat; 

%

% turn matrix to vec for anovan

iCount=1;
iCountG=1;
for iSyn=1:5
    
    SynLabel=SynS.SynLabels{iSyn}

    for iParam=1:4

        ParamLabel=SynS.ParamLabels{iParam};

        for iFilt=1:3  
            
            FiltLabel=FilterLabels{iFilt};
            
            % group mean
            Samp(iCount,1)=SynS.PChangeMat(iParam).(SynLabel).PercentChange(iFilt,:);
            Samp(iCount,2)=SynS.VarMat(iParam).(SynLabel).PercentChangeVar(iFilt,:);
            Samp(iCount,3)=SynS.AvgMinMaxMat(iParam).(SynLabel).AvgMinMax(iFilt,:);
            Samp(iCount,4)=SynS.BaseValMat(iParam).(SynLabel).BaseLineValue(iFilt,:);
            Samp(iCount,5)=SynS.RangeChangeMat(iParam).(SynLabel).RangeChange(iFilt,:);

            Grp(iCount,:)=[convertCharsToStrings(SynLabel),...
                convertCharsToStrings(ParamLabel), convertCharsToStrings(FiltLabel)];
            
            % group Gen
            for iGen=1:GenSampleSize
                
                Samp2(iCountG,1)=SynS.PChangeMat(iParam).(SynLabel).GeneratedValues(iFilt,iGen);
                Samp2(iCountG,2)=SynS.VarMat(iParam).(SynLabel).GeneratedValues(iFilt,iGen);
                Samp2(iCountG,3)=SynS.AvgMinMaxMat(iParam).(SynLabel).GeneratedValues(iFilt,iGen);
                Samp2(iCountG,4)=SynS.BaseValMat(iParam).(SynLabel).GeneratedValues(iFilt,iGen);
                Samp2(iCountG,5)=SynS.RangeChangeMat(iParam).(SynLabel).RangeChange(iFilt,:);

                Grp2(iCountG,:)=[convertCharsToStrings(SynLabel),...
                convertCharsToStrings(ParamLabel), convertCharsToStrings(FiltLabel)];
                iCountG=iCountG+1;
                
            end

            iCount=iCount+1;
            
        end
    end
end


SynS.MeanSampleVec=Samp;
SynS.MeanGroupVec=Grp;
SynS.GenMeanVec=Samp2;
SynS.GenGroupVec=Grp2;

% save('SynStatsVec2','SynS')

% save('SynStatsVec_2024')

%%
% clear all
% load('SynStatsVec_2024')
clear x
Param=SynS.GenGroupVec(:,1);
PerfMetric=SynS.GenGroupVec(:,2);
Filter=SynS.GenGroupVec(:,3);
x(:,1:2)=SynS.GenMeanVec(:,3:4); % 1-mean of left,right, 2-variance, 3-abs diff of min,max 4-BaseLine Values
SampleSize=SynS.GenSampleSize;
ParamNum=SynS.ParamNum;
% LParam=length(Param);
SynNum=SynS.SynNum;


%%interactions (filter,param,perf)
TickLabels={'VPHO' 'SPHO' 'EMG Occlusion' 'M-Wave Phase Shift' 'Occ Rate' };
SetLabels={'Range of Change (%)', 'Norm Value'};
ParamLabels={'Tracking Error', 'Tracking SNR', 'Effort Est. SNR', 'R-Squared' };

a=0.05;
iSet=1; %range of change
[gParam, IDParam]=findgroups(Param);
[gFilter, IDFilter]=findgroups(Filter);
[gPerf, IDPerf]=findgroups(PerfMetric);

% Param=Param(gParam ~= 5 );

[p,tbl,stats] = anovan(x(:,iSet),{Filter Param PerfMetric},...
    'model','full','varnames',{'Filter','Param','PerfMetric'});
figure

[results,means,~,gnames]= multcompare(stats,'Dimension',[ 1 2 3 ],'Alpha',0.05,'CType','bonferroni');
%
%identify the significant differences
clear PVals
j=1;
iG=1;
for i=1:3:60   %% based on number of variables
    
    PValInd(j:j+2,1:2)=[i i+1; i i+2; i+1 i+2;];
    PValInd(j:j+2,3)=iG;
    j=j+3;
    iG=iG+1;
    
end

for i=1:length(PValInd)
    PVals(i,4)=results(results(:,1)==PValInd(i,1) & results(:,2)==PValInd(i,2),6);
    CIs(i,:)=results(results(:,1)==PValInd(i,1) & results(:,2)==PValInd(i,2),3:4);
end 

PVals(:,1:3)=PValInd(:,1:3);
Bar1=PVals(:,1);
Bar2=PVals(:,2);
Group=PVals(:,3);
P_Vals=PVals(:,4);
PValTable=table(Bar1,Bar2,Group,P_Vals);
    
%
OrderVec=[4 5 2 1 3]; % for parameter
OrderVec2=[ 3 4 1 2];% for performance metrics
% Group=[ 3 3; 3 3; 4 4; 4 4; 4 4; 1 1; 1 1; 2 2; 2 2; 5 5; 5 5]; % columns
% Bar=  [ 2 3; 1 3; 1 2;2 3; 1 3; 2 3; 1 3; 2 3; 1 3; 1 3; 2 3]; %rows
BaseLvls=[500 700 800 800];
IndMargs=[0.0955 .0955 .095 .095];
for iPerf=1:length(IDPerf)
    iOrd2=OrderVec2(iPerf);
    clear bpsamp
    clear labels
    i=1;
    
    for iParam=1:length(IDParam)
        for iFilter=1:length(IDFilter)  
            iOrd=OrderVec(iParam);
            
            if iSet==2
                temp=x(gPerf==iOrd,iSet);
                AnvData=temp/max(temp);
                Normx(gPerf==iOrd)=AnvData;
                bpsamp(:,i)= Normx(gPerf==iOrd & gFilter==iFilter);
            else
                bpsamp(:,i)= x(gPerf==iOrd2 & gParam==iOrd & gFilter==iFilter,iSet);
            end
            
            labels(1,i)=IDParam(iOrd);
            labels(2,i)=IDFilter(iFilter);
            i=i+1;
        end
    end

    mn=mean(bpsamp);
    clear CI
    for i=1:length(mn)
        CI(:,i)=ConfInt(bpsamp(:,i),a);
        SD(iPerf,i)=std(bpsamp(:,i));
    end
    err=mn-CI(1,:);
    mn=reshape(mn,3,[]);
    err=reshape(err,3,[]);
    xaxis=(1:length(mn));
    
    figure(10)
    ax1=subplot(length(IDPerf),1,iPerf);
    brr=bar(xaxis,mn,'LineWidth',1);
    hold on
    [nbars,ngroups] = size(mn);
    xErr = nan(nbars, ngroups);
    newcolors = {'#EDB120','#0072BD','#D95319'};
    colororder(newcolors)

    for i = 1:nbars
        xErr(i,:) = brr(i).XEndPoints;
    end

    % xErr=1:length(mn);
    er = errorbar(xErr,mn,err,'Linestyle','none','LineWidth',2,'Color','k');
    grid on
    if iPerf==length(IDPerf)
        set(gca, 'XTickLabel',IDParam(OrderVec) )
    else
        set(gca,'XTick',[])
    end

%     if iPerf==1
%         title('Magnitude Change of Performance Measures with Respect to Parameters(%)')
%     end
    ylabel(IDPerf(iOrd2));
    BarOrd=[2 3 1];
    TestNum=length(IDFilter)*(length(IDParam));
    
    for iTest=(iPerf-1)*TestNum+1:iPerf*TestNum
        if PVals(iTest,4)<a
            Groupx=PVals(iTest,3);
            Barx=PVals(iTest,1:2);
            Bar=mod(Barx,3); 
            
            if Bar(1)==0 
                Bar(1)=3;
            end
            if Bar(2)==0 
                Bar(2)=3;
            end
            Group=mod(Groupx,5); 
            if Group==0
                Group=5;
            end
            
            %remapping
            Bar=BarOrd(Bar);
            BaseLvl=BaseLvls(iPerf);
            IndMarg=IndMargs(iPerf);
            Bar1=xErr(Bar(1),Group);
            Bar2=xErr(Bar(2),Group);
            sigline([Bar1 Bar2],[],BaseLvl+  5-(Bar(1)+Bar(2))*IndMarg*BaseLvl)

        end
    end
    set(gca,'FontSize',12);
end
legend(IDFilter)

%% grouped by filter
clear bpsamp
iSet=1;
[gFilter, IDFilter]=findgroups(Filter);
IDFilter;
[p,tbl,stats]  = anovan(x(:,iSet),{Filter  },...
'varnames',{'Filter'},'Display','off');

results = multcompare(stats,'Dimension',[1],'Alpha',0.05,'CType','bonferroni');
for i=1:length(IDFilter)
    bpsamp(:,i)= x(gFilter==i,iSet);
end

figure
ax1=subplot(1,1,1);
bp=boxplot(bpsamp,'Notch','on'); 
set(bp,'LineWidth', 2);
ylabel('Range of Change')
set(ax1,'XTickLabel',IDFilter)
set(gca,'FontSize',14)
yt = get(gca, 'YTick');
axis([xlim    0  max(yt)*1.05])
xt = get(gca, 'XTick'); 
hold on
plot(xt([2 3]), [1 1]*max(yt)*.85, '-k', mean([xt(2) xt(3)]), max(yt)*.9, '*k','LineWidth', 1.5)
plot(xt([1 3]), [1 1]*max(yt)*.95, '-k',  mean(xt([1 3])), max(yt)*1, '*k','LineWidth', 1.5)
% bar plots 

a=0.05;
xaxis=[1 2 3];
mn=mean(bpsamp);
figure
bar(xaxis,mn);
hold 
for i=1:3
    CI(:,i)=ConfInt(bpsamp(:,i),a);
end
er = errorbar(xaxis,mn,mn-CI(2,:),CI(1,:)-mn,'Linestyle','none','LineWidth',2,'Color','k'); 

sigline([2 3],[],max(mn))
sigline([1 3],[],1.2*max(mn))

function CI=ConfInt(x,a)
%Calculates conf interval of x with alpha equal to a 
a1= a/2;
a2=1-a/2;
SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([a1  a2],length(x)-1);      % T-Score
CI = mean(x) + ts*SEM; 

end