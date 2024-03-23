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
    
    FileLbl=FileLabels{iSyn}
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

%% save after combining
save('SynSimResults','SynS')


%% SynStats- calculate the metrics

clear all
load('SynSimResults3_rsqr');
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
    
    SynLabel=SynS.SynLabels{iSyn}

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

%% range of change
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

%%

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

save('SynStatsVec_2024')

%% group by metrics 
clear all
load('SynStatsVec_2024')
clear x
Param=SynS.GenGroupVec(:,1);
PerfMetric=SynS.GenGroupVec(:,2);
Filter=SynS.GenGroupVec(:,3);
x(:,1:2)=SynS.GenMeanVec(:,3:4); % 1-mean of left,right, 2-variance, 3-abs diff of min,max 4-BaseLine Values
SampleSize=SynS.GenSampleSize;
ParamNum=SynS.ParamNum;
LParam=length(Param);
SynNum=SynS.SynNum;

%%
iSet=2; % 2 for perf metric values 


DimNames={'Param','Filter'};
[gParam, IDParam]=findgroups(Param);
[gFilter, IDFilter]=findgroups(Filter);
[gPerf, IDPerf]=findgroups(PerfMetric);
IDPerf


for iPerf=1:4
    
    IndPerf= gPerf==iPerf;
    [p(iSet,iPerf),tbl,stats]  = anovan(x(IndPerf,iSet),{Filter(IndPerf)  },...
        'varnames',{'Filter'},'Display','off');
    
    results = multcompare(stats,'Dimension',[1],'Alpha',0.05,'CType','bonferroni','Display','off');
    multi_p_new(:,1:2)=results(:,1:2); 
    multi_p_new(:,3+iPerf)=results(:,6); 

%     Fstat(1,iPerf*3-2:iPerf*3)=[ tbl(2,3) tbl(4,3) tbl(2,6)];
%     Mstat(1,iPerf*3-2:iPerf*3)=[mean( NormGbyParamPerf(1:10,iD)) mean( NormGbyParamPerf(11:20,iD)) mean( NormGbyParamPerf(21:30,iD))];
%     
    f1=figure(100);
    xmargin=0.15;
    ymargin=0.15;
    distx=1-2*xmargin;
    disty=1-2*ymargin;
    xmrg=0.01;
    ymrg=0.07;
    NumRow=4;
    NumCol=2;
    iRow=1;
    iCol=iPerf;
    PlotLengthY=disty/NumCol;
    PlotLengthX=distx/NumRow;
    ax1=subplot('Position',[(xmargin+(xmrg+PlotLengthX)*(iCol-1)) (ymargin+(ymrg+PlotLengthY)*(NumCol-iRow))...
    PlotLengthX PlotLengthY],'Parent',f1);

%    subplot(5,5,iParam*5+iSyn-5)
%     ax1=subplot('Position',[(xmargin+(mrg+PlotLengthX)*(iPerf-1)) (ymargin+(mrg+PlotLengthY)*(2-1+1))...
%         PlotLengthX PlotLengthY],'Parent',f1);
%     ax1=subplot(2,length(gIDPerf),iPerf,'Parent',f1);

    Perf=x(IndPerf,iSet);
    gPerfFilter=gFilter(IndPerf);
    
    bpsamp(:,1)=x(gFilter==1 & gPerf==iPerf,iSet);
    bpsamp(:,2)=x(gFilter==2 & gPerf==iPerf,iSet);
    bpsamp(:,3)=x(gFilter==3 & gPerf==iPerf,iSet);
    
    Normbpsamp=abs(bpsamp./max(max(bpsamp)));
    bp=boxplot(Normbpsamp,'symbol', ''); set(bp,'Parent',ax1);

%     legend(ax1,'hide')
    set(bp,'LineWidth', 2);
    
    set(ax1,'XTickLabel',IDFilter)
%     set(ax,'YAxis','Normalized Perf. Metric')
    title(ax1,IDPerf(iPerf))
    
    yt = get(gca, 'YTick');
    axis([xlim    0  max(yt)*1.4])
    xt = get(gca, 'XTick');
    hold on
    plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k', mean([xt(2) xt(3)]), max(yt)*1.15, '*k','LineWidth', 1.5)
%     text(xt(2), max(yt)*1.16,sprintf('  p=%.2d',multi_p(2,iD)));
    plot(xt([1 3]), [1 1]*max(yt)*1.2, '-k',  mean(xt([1 3])), max(yt)*1.25, '*k','LineWidth', 1.5)
%     text(mean(xt([1 2])), max(yt)*1.26,sprintf('  p=%.2d',multi_p(3,iD)));

    if iPerf ~= 1
    set(gca,'YTick',[])
    end
    
    if iPerf==1
        ylabel(ax1,'Norm. Value of Metric')
        yt = get(gca, 'YTick');
        h=text(-1,max(yt)/2,'(b)','FontSize', 14);
        set(h,'Rotation',90);
    else 
        set (gca,'YTick', []);
    end
    
    set(gca,'FontSize',14)

    

end


%%
iSet =1;
[gParam, IDParam]=findgroups(Param);
[gFilter, IDFilter]=findgroups(Filter);
[gPerf, IDPerf]=findgroups(PerfMetric);

ID=[length(IDParam) length( IDFilter)];
% Tittels=["EESNR ", "R-Squared","Tracking Error","Tracking SNR"];


for iPerf=1:4
    
    IndPerf= gPerf==iPerf;
    [p(iSet,iPerf),tbl,stats]  = anovan(x(IndPerf,iSet),{Filter(IndPerf)  },...
        'varnames',{'Filter'});
    figure
    results = multcompare(stats,'Dimension',[1],'Alpha',0.05,'CType','bonferroni');
    multi_p_new(:,1:2)=results(:,1:2); 
    multi_p_new(:,3+iPerf)=results(:,6); 

%     Fstat(1,iPerf*3-2:iPerf*3)=[ tbl(2,3) tbl(4,3) tbl(2,6)];
%     Mstat(1,iPerf*3-2:iPerf*3)=[mean( NormGbyParamPerf(1:10,iD)) mean( NormGbyParamPerf(11:20,iD)) mean( NormGbyParamPerf(21:30,iD))];
%     
    f1=figure(100);
    xmargin=0.15;
    ymargin=0.15;
    distx=1-2*xmargin;
    disty=1-2*ymargin;
    xmrg=0.01;
    ymrg=0.07;
    NumRow=4;
    NumCol=2;
    iRow=2;
    iCol=iPerf;
    PlotLengthY=disty/NumCol;
    PlotLengthX=distx/NumRow;
    ax1=subplot('Position',[(xmargin+(xmrg+PlotLengthX)*(iCol-1)) (ymargin+(ymrg+PlotLengthY)*(NumCol-iRow))...
    PlotLengthX PlotLengthY],'Parent',f1);

%    subplot(5,5,iParam*5+iSyn-5)
%     ax1=subplot('Position',[(xmargin+(mrg+PlotLengthX)*(iPerf-1)) (ymargin+(mrg+PlotLengthY)*(2-1+1))...
%         PlotLengthX PlotLengthY],'Parent',f1);
%     ax1=subplot(2,length(gIDPerf),iPerf,'Parent',f1);

    Perf=x(IndPerf,iSet);
    gPerfFilter=gFilter(IndPerf);
    
    bpsamp(:,1)=x(gFilter==1 & gPerf==iPerf,iSet);
    bpsamp(:,2)=x(gFilter==2 & gPerf==iPerf,iSet);
    bpsamp(:,3)=x(gFilter==3 & gPerf==iPerf,iSet);
    
    bp=boxplot(bpsamp,'symbol', ''); set(bp,'Parent',ax1);

%     legend(ax1,'hide')
    set(bp,'LineWidth', 2);
    
    set(ax1,'XTickLabel',IDFilter)
%     set(ax,'YAxis','Normalized Perf. Metric')
    title(ax1,IDPerf(iPerf))
    
    yt = get(gca, 'YTick');
    axis([xlim    0  max(yt)*1.4])
    xt = get(gca, 'XTick');
    hold on
    plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k', mean([xt(2) xt(3)]), max(yt)*1.15, '*k','LineWidth', 1.5)
%     text(xt(2), max(yt)*1.16,sprintf('  p=%.2d',multi_p(2,iD)));
    plot(xt([1 3]), [1 1]*max(yt)*1.2, '-k',  mean(xt([1 3])), max(yt)*1.25, '*k','LineWidth', 1.5)
%     text(mean(xt([1 2])), max(yt)*1.26,sprintf('  p=%.2d',multi_p(3,iD)));

    if iPerf ~= 1
    set(gca,'YTick',[])
    end
    
    if iPerf==1
        ylabel(ax1,'Percent Change of Metric')
        yt = get(gca, 'YTick');
        h=text(-1,max(yt)/2,'(a)','FontSize', 14);
        set(h,'Rotation',90);
    else 
        set (gca,'YTick', []);
    end
    
    set(gca,'FontSize',14)
    
end

%% Interactions 

iSet =2;
[gParam, IDParam]=findgroups(Param);
[gFilter, IDFilter]=findgroups(Filter);
[gPerf, IDPerf]=findgroups(PerfMetric);

ID=[length(IDParam) length( IDFilter)];
Tittels=["EESNR ", "R-Squared","Tracking Error","Tracking SNR"];

for iPerf=1:1

    IndPerf= gPerf==iPerf;
    [~,~,stats]  = anovan(x(IndPerf,iSet),{ Param(IndPerf)   Filter(IndPerf) },...
        'model','interaction', 'varnames',{'Param','Filter'});
    figure
    results = multcompare(stats,'Dimension',[ 2 1],'Alpha',0.05,'CType','bonferroni');
%     multi_p_new(:,1:2)=results(:,1:2); 
%     multi_p_new(:,3+iPerf)=results(:,6); 

    bplot

end

%%
iSet=1; % 1 for percent change

[~,~,stats] = anovan(x(:,iSet),{ Param PerfMetric Filter},...
    'model','interaction','varnames',{'Param','PerfMetric','Filter'});
figure
    results = multcompare(stats,'Dimension',[3 2 1],'Alpha',0.05,'CType','bonferroni');


% Test the two-factor interactions. This time specify the variable names.
% [~,~,stats] = anovan(x(:,iSet),{Param PerfMetric Filter },'model','interaction','varnames',{'Param','PerfMetric','Filter'});
% figure
% results = multcompare(stats,'Dimension',[  1 3]);


%%
[GbyParam,gID]=GroupBy(x(:,iSet),Param);
% [p,tbl,stats] = anovan(GbyParam(:,1),{Filter(1:90) PerfMetric(1:90)});

[GbyParamFilt,gIDFilt]=GroupBy(GbyParam(:,1),Filter(1:LParam/SynNum));
[GbyParamPerf,gIDPerf]=GroupBy(GbyParam(:,1),PerfMetric(1:LParam/SynNum));

NormGbyParamPerf=GbyParamPerf./max(GbyParamPerf);  % use for normalized 

%%
clear multi_p
close all
DataSet= [ 3 4 1 2 ]; % change this to get the order in columns of the figure
Titels=["EESNR ", "R-Squared","Tracking Error","Tracking SNR"];
% Titels=gIDPerf;
for iPerf=1:length(gIDPerf)
    iD=DataSet(iPerf);

    [GbyParamPerfFilt,gIDPerfFilt]=GroupBy(NormGbyParamPerf(:,iD),Filter(1:LParam/SynNum/ParamNum));
    [p(iSet,iD),tbl,stats] = anovan(NormGbyParamPerf(:,iD),{Filter(1:LParam/SynNum/ParamNum)  },'display','off');

    %   figure
    results = multcompare(stats,'Dimension',[ 1 ],'Display','off');
    multi_p(:,1:2)=results(:,1:2); % 1- comb 2-GS 3-Blanking
    multi_p(:,3+iPerf)=results(:,6); 
    Fstat(1,iPerf*3-2:iPerf*3)=[ tbl(2,3) tbl(4,3) tbl(2,6)];
    Mstat(1,iPerf*3-2:iPerf*3)=[mean( NormGbyParamPerf(1:10,iD)) mean( NormGbyParamPerf(11:20,iD)) mean( NormGbyParamPerf(21:30,iD))];
    f1=figure(100);
    xmargin=0.15;
    ymargin=0.15;
    distx=1-2*xmargin;
    disty=1-2*ymargin;
    xmrg=0.01;
    ymrg=0.07;
    NumRow=4;
    NumCol=2;
    iRow=1;
    iCol=iPerf;
    PlotLengthY=disty/NumCol;
    PlotLengthX=distx/NumRow;
    ax1=subplot('Position',[(xmargin+(xmrg+PlotLengthX)*(iCol-1)) (ymargin+(ymrg+PlotLengthY)*(NumCol-iRow))...
    PlotLengthX PlotLengthY],'Parent',f1);

%    subplot(5,5,iParam*5+iSyn-5)
%     ax1=subplot('Position',[(xmargin+(mrg+PlotLengthX)*(iPerf-1)) (ymargin+(mrg+PlotLengthY)*(2-1+1))...
%         PlotLengthX PlotLengthY],'Parent',f1);
%     ax1=subplot(2,length(gIDPerf),iPerf,'Parent',f1);
    bp=boxplot(GbyParamPerfFilt,'Notch','off','symbol', ''); set(bp,'Parent',ax1);

%     legend(ax1,'hide')
    set(bp,'LineWidth', 2);

    
    set(ax1,'XTickLabel',gIDPerfFilt)
%     set(ax,'YAxis','Normalized Perf. Metric')
    title(ax1,Titels(iD))
    
    yt = get(gca, 'YTick');
    axis([xlim    0  max(yt)*1.4])
    xt = get(gca, 'XTick');
    hold on
    plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k', mean([xt(2) xt(3)]), max(yt)*1.15, '*k','LineWidth', 1.5)
%     text(xt(2), max(yt)*1.16,sprintf('  p=%.2d',multi_p(2,iD)));
    plot(xt([1 3]), [1 1]*max(yt)*1.2, '-k',  mean(xt([1 3])), max(yt)*1.25, '*k','LineWidth', 1.5)
%     text(mean(xt([1 2])), max(yt)*1.26,sprintf('  p=%.2d',multi_p(3,iD)));

    if iPerf ~= 1
    set(gca,'YTick',[])
    end
    
    if iPerf==1
        ylabel(ax1,'Norm. Value of Metric')
        yt = get(gca, 'YTick');
        h=text(-1,max(yt)/2,'(a)','FontSize', 14);
        set(h,'Rotation',90);
    else 
        set (gca,'YTick', []);
    end
    
    set(gca,'FontSize',14)

end
%
% boxplot with groups 
% Ord= [ 3 4 1 2 ]; % change this to get the order in columns of the figure
% Titels=gIDPerf;
% for iPerf=1:length(gIDPerf)
% %     [d(:,),gIDPerfFilt]=GroupBy(GbyParamPerf(:,iPerf),Filter(1:LParam/SynNum/ParamNum));
% end
% data = {rand(100,2), rand(100,2)+.2, rand(100,2)-.2}; 
% boxplotGroup(data, 'PrimaryLabels', {'a' 'b' 'c'}, ...
%   'SecondaryLabels',{'Group1', 'Group2'}, 'InterGroupSpace', 2)


iSet=1;
% [GbyParam,gID]=GroupBy(x(:,1),Param);
[GbyPerf,gIDPerf]=GroupBy(x(:,iSet),PerfMetric(1:LParam/SynNum));

Titles=gIDPerf;
for iPerf=1:length(gIDPerf)
    iD=DataSet(iPerf);
    
    [GbyPerfFilt,gIDPerfFilt]=GroupBy(GbyPerf(:,iD),Filter(1:LParam/SynNum/ParamNum));
    [p(iSet,iD),tbl,stats] = anovan(GbyPerf(:,iD),{Filter(1:LParam/SynNum/ParamNum)  },'Display','off');
%     figure
    results = multcompare(stats,'Dimension',[ 1 ],'Display','off');
    multi_p(4:6,1:2)=results(:,1:2); % 1- comb 2-GS 3-Blanking
    multi_p(4:6,3+iPerf)=results(:,6); 
    Fstat(2,iPerf*3-2:iPerf*3)=[ tbl(2,3) tbl(4,3) tbl(2,6)];
    Mstat(2,iPerf*3-2:iPerf*3)=[mean( GbyPerf(1:10,iD)) mean( GbyPerf(11:20,iD)) mean( GbyPerf(21:30,iD))];

    f1=figure(100);
    
%     xmargin=0.1;
%     ymargin=0.2;
%     distx=1-2*xmargin;
%     disty=1-2*ymargin;
%     mrg=0.01;
%     NumRow=4;
%     NumCol=2;
    iRow=2;
    iCol=iPerf;
    PlotLengthY=disty/(NumCol);
    PlotLengthX=distx/NumRow;
    ax2=subplot('Position',[(xmargin+(xmrg+PlotLengthX)*(iCol-1)) (ymargin+(ymrg+PlotLengthY)*(NumCol-iRow))...
    PlotLengthX PlotLengthY],'Parent',f1);

%   ax2=subplot(2,4,iPerf+4,'Parent',f1);
    bp2=boxplot(GbyPerfFilt,'Notch','off','symbol', ''); set(bp2,'Parent',ax2);
    set(bp2,'LineWidth', 2);

    set(gca,'XTickLabel',gIDPerfFilt)
%     title(Titels(iD))
    yt = get(gca, 'YTick');
    axis([xlim 0 ceil(max(yt)*1.4)])
    xt = get(gca, 'XTick');
    
    if iD == 3
        
        hold on
        plot(xt([2 1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([2 1])), max(yt)*1.15, '*k','LineWidth', 1.5)
%         text(xt(1), max(yt)*1.16,sprintf('  p=%.2d',multi_p(1,iD)));
        plot(xt([2 3]), [1 1]*max(yt)*1.2, '-k',  mean(xt([2 3])), max(yt)*1.25, '*k','LineWidth', 1.5)
%         text(xt(2), max(yt)*1.26,sprintf('  p=%.2d',multi_p(3,iD)));

    else
        
        hold on
        plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k',  mean(xt([2 3])), max(yt)*1.15, '*k','LineWidth', 1.5)
%         text(xt(2), max(yt)*1.16,sprintf('  p=%.2d',multi_p(3,iD)));
        plot(xt([1 3]), [1 1]*max(yt)*1.2, '-k',  mean(xt([1 3])), max(yt)*1.25, '*k','LineWidth', 1.5)
%         text(mean(xt([1 2])), max(yt)*1.26,sprintf('  p=%.2d',multi_p(2,iD)));

    end
    
    if iPerf == 4|| iPerf == 3
        
        hold on
        plot(xt([1 2]), [1 1]*max(yt)*1.3, '-k',  mean(xt([1 2])), max(yt)*1.35, '*k','LineWidth', 1.5)
%         text(mean(xt([1 ])), max(yt)*1.36,sprintf('  p=%.2d',multi_p(1,iD)));
    end
    
        
    if iPerf==1
        ylabel(ax2,'Percent Change of Metric')
        yt = get(gca, 'YTick');
        h=text(-1,max(yt)/2,'(b)','FontSize', 14);
        set(h,'Rotation',90);
    else 
        set (gca,'YTick', []);
    end

    set(gca,'FontSize',14)

end


%% group by each syn or par --- Dont use

clear all
load('SynStatsVec2')

Param=SynS.GenGroupVec(:,1);
PerfMetric=SynS.GenGroupVec(:,2);
Filter=SynS.GenGroupVec(:,3);
x(:,1:2)=SynS.GenMeanVec(:,3:4); % 1-mean of left,right, 2-variance, 3-abs diff of min,max 4-BaseLine Values
SampleSize=SynS.GenSampleSize;
ParamNum=SynS.ParamNum;
LParam=length(Param);
SynNum=SynS.SynNum;
FilterNum=3;
%
% [p,tbl,stats] = anovan(x,{g1 g3,g2});
%Test the two-factor interactions. This time specify the variable names.
% [~,~,stats] = anovan(y,{g1 g2 g3 },'model','interaction','varnames',{'g1','g2','g3'})
% figure
% results = multcompare(stats,'Dimension',[  2 3]);
%%
iSet=2;

[GbyParam,gIDParam]=GroupBy(x(:,iSet),Param);
% [p,tbl,stats] = anovan(GbyParam(:,1),{Filter(1:90) PerfMetric(1:90)});


[GbyParamFilt,gIDFilt]=GroupBy(GbyParam(:,1),Filter(1:LParam/SynNum));
NormGbyParam=GbyParam./max(GbyParam);  % use for normalized 

%%
DataSet= [ 3 4 1 2 5]; % change this to get the order in columns of the figure
Titels=gIDParam;
for iParam=1:length(gIDParam)
    iD=DataSet(iParam);

    [GbyParamFilt,gIDFilt]=GroupBy(NormGbyParam(:,iD),Filter(1:LParam/SynNum));
    [p(iSet,iD),tbl,stats] = anovan(NormGbyParam(:,iD),{Filter(1:LParam/SynNum)  },'display','off');
  
    results = multcompare(stats,'Dimension',[ 1 ],'Display','off');
    multi_p(4:6,1:2)=results(:,1:2); % 1- comb 2-GS 3-Blanking
    multi_p(4:6,3+iParam)=results(:,6); 
    Fstat(2,iParam*3-2:iParam*3)=[ tbl(2,3) tbl(4,3) tbl(2,6)];
    Mstat(2,iParam*3-2:iParam*3)=[mean( NormGbyParam(1:10,iD)) mean( NormGbyParam(11:20,iD)) mean( NormGbyParam(21:30,iD))];

    f1=figure(100); 
    xmargin=0.1;
    ymargin=0.1;
    distx=1-2*xmargin;
    disty=1-2*ymargin;
    xmrg=0.01;
    ymrg=0.05;
    NumRow=5;
    NumCol=2;
    iRow=1;
    iCol=iParam;
    PlotLengthY=disty/NumCol;
    PlotLengthX=distx/NumRow;
    ax1=subplot('Position',[(xmargin+(xmrg+PlotLengthX)*(iCol-1)) (ymargin+(ymrg+PlotLengthY)*(NumCol-iRow))...
    PlotLengthX PlotLengthY],'Parent',f1);

%    subplot(5,5,iParam*5+iSyn-5)
%     ax1=subplot('Position',[(xmargin+(mrg+PlotLengthX)*(iPerf-1)) (ymargin+(mrg+PlotLengthY)*(2-1+1))...
%         PlotLengthX PlotLengthY],'Parent',f1);
%     ax1=subplot(2,length(gIDPerf),iPerf,'Parent',f1);
    bp=boxplot(GbyParamFilt,'Notch','on'); set(bp,'Parent',ax1);
%     legend(ax1,'hide')
    set(bp,'LineWidth', 2);

    
    set(ax1,'XTickLabel',gIDFilt)
%     set(ax,'YAxis','Normalized Perf. Metric')
    title(ax1,Titels(iD))
    
    yt = get(gca, 'YTick');
    axis([xlim    0  max(yt)*1.4])
    xt = get(gca, 'XTick');
    hold on
    plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k',  mean(xt([2 3])), max(yt)*1.15, '*k', '*k','LineWidth', 2)
    plot(xt([1 3]), [1 1]*max(yt)*1.2, '-k',  mean(xt([1 3])), max(yt)*1.25, '*k', '*k','LineWidth', 2)
    if iParam ~= 1
    set(gca,'YTick',[])
    end
    
    if iParam==1
        ylabel(ax1,'Norm. Value of Metric')
    else 
        set (gca,'YTick', []);
    end
    set(gca,'FontSize',12)

end
%
% boxplot with groups 
% Ord= [ 3 4 1 2 ]; % change this to get the order in columns of the figure
% Titels=gIDPerf;
% for iPerf=1:length(gIDPerf)
% %     [d(:,),gIDPerfFilt]=GroupBy(GbyParamPerf(:,iPerf),Filter(1:LParam/SynNum/ParamNum));
% end
% data = {rand(100,2), rand(100,2)+.2, rand(100,2)-.2}; 
% boxplotGroup(data, 'PrimaryLabels', {'a' 'b' 'c'}, ...
%   'SecondaryLabels',{'Group1', 'Group2'}, 'InterGroupSpace', 2)

%%
clear multi_p
DataSet= [ 3 4 2 1 5];  % change this for different order of plots
iSet=1;
[GbyParam,gID]=GroupBy(x(:,1),Param);
lblplot={'(a)', '(b)','(c)','(d)','(e)',};
Titles=gID;
for iParam=1:length(gID)
    iD=DataSet(iParam);
    
    [GbyParamFilt,gIDPerfFilt]=GroupBy(GbyParam(:,iD),Filter(1:LParam/SynNum));
    [p(iSet,iD),tbl,stats] = anovan(GbyParam(:,iD),{Filter(1:LParam/SynNum)  },'Display','off');
%     results = multcompare(stats,'Dimension',[ 1 ],'Display','off');
%     multi_p(:,iD)=results(:,6);
%     figure
    results = multcompare(stats,'Dimension',[ 1 ],'Display','off');
    multi_p(1:3,1:2)=results(:,1:2); % 1- comb 2-GS 3-Blanking
    multi_p(1:3,3+iParam)=results(:,6); 
    Fstat(1,iParam*3-2:iParam*3)=[ tbl(2,3) tbl(4,3) tbl(2,6)];
    Mstat(1,iParam*3-2:iParam*3)=[mean( GbyParam(1:10,iD)) mean( GbyParam(11:20,iD)) mean( GbyParam(21:30,iD))];

    
    f1=figure(100);
    xmargin=0.1;
    ymargin=0.2;
    distx=1-2*xmargin;
    disty=1-2*ymargin;
    mrg=0.01;
    xmrg=0.01;
    ymrg=0.05;
    NumRow=5;
    NumCol=2;
    iRow=2;
    iCol=iParam;
    PlotLengthY=disty/(NumCol);
    PlotLengthX=distx/NumRow;
    ax2=subplot('Position',[(xmargin+(xmrg+PlotLengthX)*(iCol-1)) (ymargin+(ymrg+PlotLengthY)*(NumCol-iRow))...
    PlotLengthX PlotLengthY],'Parent',f1);

%   ax2=subplot(2,4,iPerf+4,'Parent',f1);
    bp2=boxplot(GbyParamFilt,'Notch','off','symbol', ''); set(bp2,'Parent',ax2);
    text(1.8,-130,lblplot(iParam),'FontSize', 14);
    ylim([0 400]);
    
    set(bp2,'LineWidth', 2);

    set(gca,'XTickLabel',gIDPerfFilt)
    title(Titles(iD))
    yt = get(gca, 'YTick');
    axis([xlim 0 ceil(max(yt)*1.4)])
    xt = get(gca, 'XTick');
        if iParam==1
        ylabel(ax2,'Percent Change of Metric')
    else 
        set (gca,'YTick', []);
    end
        set(gca,'FontSize',12)

    if iParam == 1
        
        hold on
        plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k',  mean(xt([2 3])), max(yt)*1.15, '*k','LineWidth', 1.5)
        plot(xt([1 3]), [1 1]*max(yt)*1.2, '-k',  mean(xt([1 3])), max(yt)*1.25, '*k','LineWidth', 1.5)

%         text(mean(xt([2])), max(yt)*1.16,sprintf('  p=%.2d',multi_p(3,iD)));

    elseif iParam==2
        
        hold on
        plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k',  mean(xt([2 3])), max(yt)*1.15, '*k','LineWidth', 1.5)
        plot(xt([1 3]), [1 1]*max(yt)*1.2, '-k',  mean(xt([1 3])), max(yt)*1.25, '*k','LineWidth', 1.5)

%         plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k','LineWidth', 1.5)
%         text(mean(xt([1])), max(yt)*1.16,sprintf('  p=%.2d',multi_p(1,iD)));

    elseif iParam==3
        hold on


        plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k',  mean(xt([2 3])), max(yt)*1.15, '*k','LineWidth', 1.5)
%         text(mean(xt([2])), max(yt)*1.16,sprintf('  p=%.2d',multi_p(3,iD)));
        plot(xt([1 3]), [1 1]*max(yt)*1.2, '-k',  mean(xt([1 3])), max(yt)*1.25, '*k','LineWidth', 1.5)
%         text(mean(xt([1 2])), max(yt)*1.26,sprintf('  p=%.2d',multi_p(2,iD)));
    elseif iParam==5
        hold on
        plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k',  mean(xt([2 3])), max(yt)*1.15, '*k','LineWidth', 1.5)
%         text(mean(xt([2])), max(yt)*1.16,sprintf('  p=%.2d',multi_p(3,iD)));
        plot(xt([1 3]), [1 1]*max(yt)*1.2, '-k',  mean(xt([1 3])), max(yt)*1.25, '*k','LineWidth', 1.5)
%         text(mean(xt([1 2])), max(yt)*1.26,sprintf('  p=%.2d',multi_p(2,iD)));

    else 
        
        hold on
        plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k',  mean(xt([2 3])), max(yt)*1.15, '*k','LineWidth', 1.5)
%         text(mean(xt([2])), max(yt)*1.16,sprintf('  p=%.2d',multi_p(3,iD)));
        plot(xt([1 3]), [1 1]*max(yt)*1.2, '-k',  mean(xt([1 3])), max(yt)*1.25, '*k','LineWidth', 1.5)
%         text(mean(xt([1 2])), max(yt)*1.26,sprintf('  p=%.2d',multi_p(2,iD)));
        plot(xt([1 2]), [1 1]*max(yt)*1.3, '-k',  mean(xt([1 2])), max(yt)*1.35, '*k','LineWidth', 1.5)
%         text(mean(xt([1])), max(yt)*1.36,sprintf('  p=%.2d',multi_p(1,iD)));
    end
    
        


end




%%
% clear all
% load('SynStatsVec')

Param=SynS.GenGroupVec(:,1);
PerfMetric=SynS.GenGroupVec(:,2);
Filter=SynS.GenGroupVec(:,3);
x(:,1:2)=SynS.GenMeanVec(:,3:4); % 1-mean of left,right, 2-variance, 3-abs diff of min,max 4-BaseLine Values
iSet=2;
iParam=1;
[GbyParam,gIDParam1]=GroupBy(x(:,iSet),Param);
[GbyParamPerf,gIDPerf]=GroupBy(GbyParam(:,iParam),PerfMetric(1:LParam/SynNum));
NormGbyParamPerf=GbyParamPerf./max(GbyParamPerf);

[NormGbyParam,g]=BreakGroup(NormGbyParamPerf,gIDPerf);
[NormGbyParamFilt,gIDFilt]=GroupBy(NormGbyParam,Filter(1:LParam/SynNum));

figure(101)
subplot(1,2,1)
boxplot(NormGbyParamFilt,'Notch','on')
set(gca,'XTickLabel',gIDFilt)
title(sprintf('Norm. %s ',gIDParam1(iParam)));
[p(1),~,stats] = anovan(NormGbyParam,{Filter(1:LParam/SynNum)  });
%%
iSet=1;
[GbyParam,gIDParam2]=GroupBy(x(:,iSet),Param);

[GbyParamFilt,gIDFilt]=GroupBy(GbyParam(:,iParam),Filter(1:LParam/SynNum));
% [GbyParamPerf,gIDPerf]=GroupBy(GbyParam(:,1),PerfMetric(1:LParam/SynNum));

figure(101)
subplot(1,2,2)
boxplot(GbyParamFilt,'Notch','on')
set(gca,'XTickLabel',gIDFilt)
title('Percent Change in Perf. Metrics (%s)',gIDParam2(iParam));
[p(2),~,stats] = anovan(GbyParam(:,iParam),{Filter(1:LParam/SynNum)});

%%






function [GroupedVal,gID]=GroupBy (x,g)

    [GroupNum gID]=findgroups(g);
    for i=1:max(GroupNum)
        GroupedVal(:,i)=x(GroupNum==i);
    end
end

function [y,g]=BreakGroup (x,gID)

    y=reshape(x,[],1);

    for i=1:size(x,2)
        gg(1:size(x,1),i)=i;
    end
    gg=reshape(gg,[],1);
    for i=1:length(gg)
        for j=1:length(gID)
            if gg(i)==j
                g(i)=gID(j);
            end
        end  
    end
    g=g';
end
