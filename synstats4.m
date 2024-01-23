clear all
load('SynStatsVec2')

Param=SynS.GenGroupVec(:,1);
PerfMetric=SynS.GenGroupVec(:,2);
Filter=SynS.GenGroupVec(:,3);
x(:,1:2)=SynS.GenMeanVec(:,3:4); % 1-mean of left,right, 2-variance, 3-abs diff of min,max 4-BaseLine Values
SampleSize=SynS.GenSampleSize;
ParamNum=SynS.ParamNum;
% LParam=length(Param);
SynNum=SynS.SynNum;
%% grouped by filter
iSet=1;
[gFilter, IDFilter]=findgroups(Filter);
IDFilter
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
%%
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
%% interactions (filter,param,perf)
TickLabels={'VPHO' 'SPHO' 'EMG Occlusion' 'M-Wave Phase Shift' 'Stim. Freq.' };
SetLabels={'Range of Change (%)', 'Norm Value'};
ParamLabels={'Tracking Error', 'Tracking SNR', 'Effort Est. SNR', 'R-Squared' };

a=0.05;
iSet=1; %range of change
[gParam, IDParam]=findgroups(Param);
[gFilter, IDFilter]=findgroups(Filter);
[gPerf, IDPerf]=findgroups(PerfMetric);

Param=Param(gParam ~= 5 );

[p,tbl,stats] = anovan(x(gParam ~= 5 ,iSet),{Filter(gParam ~= 5 ) Param(gParam ~= 5 ) PerfMetric(gParam ~= 5 )},...
    'model','full','varnames',{'Filter','Param','PerfMetric'});
figure
results = multcompare(stats,'Dimension',[ 1 2 3 ],'Alpha',0.05,'CType','bonferroni');
%%
%identify the significant differences
clear PVals
j=1;
iG=1;
for i=1:3:48   %% based on number of variables
    
    PValInd(j:j+2,1:2)=[i i+1; i i+2; i+1 i+2;];
    PValInd(j:j+2,3)=iG;
    j=j+3;
    iG=iG+1;
    
end

for i=1:length(PValInd)
    PVals(i,4)=results(results(:,1)==PValInd(i,1) & results(:,2)==PValInd(i,2),6);
end

PVals(:,1:3)=PValInd(:,1:3);
Bar1=PVals(:,1);
Bar2=PVals(:,2);
Group=PVals(:,3);
P_Vals=PVals(:,4);
PValTable=table(Bar1,Bar2,Group,P_Vals);
    
%%
OrderVec=[3 4 2 1  ]; % for parameter
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
    
    for iParam=1:length(IDParam) -1
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
    brr=bar(xaxis,mn,'LineWidth',2);
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
        set(gca, 'XTickLabel',TickLabels(1:5) )
    else
        set(gca,'XTick',[])
    end

    if iPerf==1
        title('Magnitude Change of Performance Measures with Respect to Parameters(%)')
    end
    ylabel(ParamLabels(iPerf));
    BarOrd=[2 3 1];
    TestNum=length(IDFilter)*(length(IDParam)-1);
    
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
            Group=mod(Groupx,4); 
            if Group==0
                Group=4;
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

%% interactions (filter,param)
iSet=1;
[gParam, IDParam]=findgroups(Param);
[gFilter, IDFilter]=findgroups(Filter);
[gPerf, IDPerf]=findgroups(PerfMetric);
[p,tbl,stats] = anovan(x(:,iSet),{Filter PerfMetric Param},...
    'model','interaction','varnames',{'Filter','PerfMetric','Param'});
figure
results = multcompare(stats,'Dimension',[ 1 3 ],'Alpha',0.05,'CType','bonferroni');

clear bpsamp
clear labels

i=1;
for iParam=1:length(IDParam) 
    for iFilter=1:length(IDFilter)  
        bpsamp(:,i)= x(gParam==iParam & gFilter==iFilter,iSet);
        labels(1,i)=IDParam(iParam);
        labels(2,i)=IDFilter(iFilter);
        i=i+1;
    end
end

IDFilter
IDParam
figure
ax1=subplot(2,1,1);
bp=boxplot(bpsamp,'Notch','on'); 

yt = get(gca, 'YTick');
set(ax1,'XTickLabel',labels(2,:))

for iText=1:length(IDParam)
    text(iText*3-1,-max(yt)/4,IDParam(iText),'FontSize', 10,'HorizontalAlignment', 'center');
    annotation('line',[.13+0.155*(iText-1) .13+0.155*(iText-1)],[.92 .5])
end
annotation('line',[.13+0.155*(length(IDParam)) .13+0.155*(length(IDParam))],[.92 .5])
yt = get(gca, 'YTick');
axis([xlim    0  max(yt)*1.15])
xt = get(gca, 'XTick');
hold on
plot(xt([2 3]), [1 1]*max(yt)*.85, '-k', mean([xt(2) xt(3)]), max(yt)*.9, '*k','LineWidth', 1.5)
plot(xt([1 3]), [1 1]*max(yt)*.95, '-k',  mean(xt([1 3])), max(yt)*1, '*k','LineWidth', 1.5)
plot(xt([1 2]), [1 1]*max(yt)*1.05, '-k',  mean(xt([1 2])), max(yt)*1.1, '*k','LineWidth', 1.5)

plot(xt([5 6]), [1 1]*max(yt)*.85, '-k', mean(xt([5 6])), max(yt)*.9, '*k','LineWidth', 1.5)
plot(xt([4 6]), [1 1]*max(yt)*.95, '-k', mean(xt([4 6])), max(yt)*1, '*k','LineWidth', 1.5)

plot(xt([8 9]), [1 1]*max(yt)*.85, '-k', mean(xt([8 9])), max(yt)*.9, '*k','LineWidth', 1.5)
plot(xt([7 9]), [1 1]*max(yt)*.95, '-k', mean(xt([7 9])), max(yt)*1, '*k','LineWidth', 1.5)

plot(xt([11 12]), [1 1]*max(yt)*.85, '-k', mean(xt([11 12])), max(yt)*.9, '*k','LineWidth', 1.5)
plot(xt([10 12]), [1 1]*max(yt)*.95, '-k', mean(xt([10 12])), max(yt)*1, '*k','LineWidth', 1.5)

plot(xt([14 15]), [1 1]*max(yt)*.85, '-k', mean(xt([14 15])), max(yt)*.9, '*k','LineWidth', 1.5)
plot(xt([13 15]), [1 1]*max(yt)*.95, '-k', mean(xt([13 15])), max(yt)*1, '*k','LineWidth', 1.5)
%%
TickLabels={'VPHO' 'SPHO'  'EMG Occlusion' 'M-Wave Delay' 'Stim. Freq.' };
OrderVec=[3 4 2 1 5 ];
a=0.05;
clear bpsamp
clear labels

i=1;
for iParam=1:length(IDParam) 
    for iFilter=1:length(IDFilter)  
        iOrd=OrderVec(iParam);
        bpsamp(:,i)= x(gParam==iOrd & gFilter==iFilter,iSet);
        labels(1,i)=IDParam(iOrd);
        labels(2,i)=IDFilter(iFilter);
        i=i+1;
    end
end


% ReOrdbpsamp(1:5,:)=bpsamp(OrderVec);
mn=mean(bpsamp);
clear CI
for i=1:length(mn)
    CI(:,i)=ConfInt(bpsamp(:,i),a);
end
err=mn-CI(1,:);
mn=reshape(mn,3,[]);
err=reshape(err,3,[]);
xaxis=(1:length(mn));

figure(10)
subplot(5,1,5)
brr=bar(xaxis,mn,'LineWidth',2);
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
set(gca, 'XTickLabel',TickLabels(1:5) )
ylabel('All Above Merged');

% Group=[ 1 1; 1 1; 1 1;2 2; 2 2; 3 3; 3 3; 4 4; 4 4; 5 5; 5 5]; % columns
Group=[ 3 3; 3 3; 4 4; 4 4; 4 4; 1 1; 1 1; 2 2; 2 2; 5 5; 5 5]; % columns
Bar=  [ 2 3; 1 3; 1 2;2 3; 1 3; 2 3; 1 3; 2 3; 1 3; 1 3; 2 3]; %rows
[IndNum,~]=size(Group);
BaseLvl=420;
IndMarg=0.075;
for iInd=1:IndNum
    
    Bar1=xErr(Bar(iInd,1),Group(iInd,1));
    Bar2=xErr(Bar(iInd,2),Group(iInd,2));
    sigline([Bar1 Bar2],[],BaseLvl+  5-(Bar(iInd,1)+Bar(iInd,2))*IndMarg*BaseLvl)
    
end


set(gca,'FontSize',12);
legend(IDFilter)

%% interactions (filter,perf)
iSet=1;
[gParam, IDParam]=findgroups(Param);
[gFilter, IDFilter]=findgroups(Filter);
[gPerf, IDPerf]=findgroups(PerfMetric);
[p,tbl,stats] = anovan(x(:,iSet),{  Param  Filter PerfMetric },...
    'model','interaction','varnames',{'Param','Filter','PerfMetric'},'Display','off');
figure
results = multcompare(stats,'Dimension',[ 2 3 ],'Alpha',0.05,'CType','bonferroni');

clear bpsamp
clear labels
i=1;
for iPerf=1:length(IDPerf) 
    for iFilter=1:length(IDFilter)  
        
        bpsamp(:,i)= x(gPerf==iPerf & gFilter==iFilter,iSet);
        labels(1,i)=IDPerf(iPerf);
        labels(2,i)=IDFilter(iFilter);

        i=i+1;
    end
end

IDFilter
IDPerf
figure
ax1=subplot(2,1,1);
bp=boxplot(bpsamp,'Notch','on'); 

yt = get(gca, 'YTick');
set(ax1,'XTickLabel',labels(2,:))

initDist=.13;
incDist=0.1935;
for iText=1:length(IDPerf)
    
    text(iText*3-1,-max(yt)/4,IDPerf(iText),'FontSize', 10,'HorizontalAlignment', 'center');
    annotation('line',[initDist+incDist*(iText-1) initDist+incDist*(iText-1)],[.92 .5])

end
annotation('line',[initDist+incDist*(length(IDPerf)) initDist+incDist*(length(IDPerf))],[.92 .5])
yt = get(gca, 'YTick');
axis([xlim    0  max(yt)*1.15])
xt = get(gca, 'XTick');
hold on
plot(xt([2 3]), [1 1]*max(yt)*.85, '-k', mean([xt(2) xt(3)]), max(yt)*.9, '*k','LineWidth', 1.5)
plot(xt([1 3]), [1 1]*max(yt)*.95, '-k',  mean(xt([1 3])), max(yt)*1, '*k','LineWidth', 1.5)
plot(xt([1 2]), [1 1]*max(yt)*1.05, '-k',  mean(xt([1 2])), max(yt)*1.1, '*k','LineWidth', 1.5)

plot(xt([5 6]), [1 1]*max(yt)*.85, '-k', mean(xt([5 6])), max(yt)*.9, '*k','LineWidth', 1.5)
plot(xt([4 6]), [1 1]*max(yt)*.95, '-k', mean(xt([4 6])), max(yt)*1, '*k','LineWidth', 1.5)

plot(xt([8 9]), [1 1]*max(yt)*.85, '-k', mean(xt([8 9])), max(yt)*.9, '*k','LineWidth', 1.5)
plot(xt([7 9]), [1 1]*max(yt)*.95, '-k', mean(xt([7 9])), max(yt)*1, '*k','LineWidth', 1.5)

plot(xt([11 12]), [1 1]*max(yt)*.85, '-k', mean(xt([11 12])), max(yt)*.9, '*k','LineWidth', 1.5)
plot(xt([10 12]), [1 1]*max(yt)*.95, '-k', mean(xt([10 12])), max(yt)*1, '*k','LineWidth', 1.5)
%%
TickLabels={'Tracking Error' 'Tracking SNR' 'EESNR' 'EEA'};
OrderVec=[3 4 1 2 ];
a=0.05;
clear bpsamp
clear labels
i=1;
for iPerf=1:length(IDPerf) 
    for iFilter=1:length(IDFilter)  
        iOrd=OrderVec(iPerf);
        bpsamp(:,i)= x(gPerf==iOrd & gFilter==iFilter,iSet);
        labels(1,i)=IDPerf(iOrd);
        labels(2,i)=IDFilter(iFilter);
        i=i+1;
    end
end

mn=mean(bpsamp);
clear CI
for i=1:length(mn)
    CI(:,i)=ConfInt(bpsamp(:,i),a);
end
err=mn-CI(1,:);
mn=reshape(mn,3,[]);
err=reshape(err,3,[]);

xaxis=(1:length(mn));

figure(1)
subplot(2,1,2)
brr=bar(xaxis,mn,'LineWidth',2);
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
set(gca, 'XTickLabel',TickLabels )
ylabel({'(b)'; 'Range of Change (%)'})

% Group=[ 1 1; 1 1; 2 2; 2 2; 4 4; 4 4]; %columns
Group=[ 3 3; 3 3; 4 4; 4 4; 2 2; 2 2]; %columns
Bar=  [ 2 3; 1 3; 2 3; 1 3; 2 3; 1 3]; %rows
[IndNum,~]=size(Group);
BaseLvl=420;
IndMarg=0.065;

for iInd=1:IndNum
    
    Bar1=xErr(Bar(iInd,1),Group(iInd,1));
    Bar2=xErr(Bar(iInd,2),Group(iInd,2));
    sigline([Bar1 Bar2],[],BaseLvl + 5-(Bar(iInd,1)+Bar(iInd,2))*IndMarg*BaseLvl)
    
end
set(gca,'FontSize',12);
legend(IDFilter)

%% Metric values (perf, filter)
TickLabels={'Tracking Error' 'Tracking SNR' 'Est. Effort SNR' 'Effort Est. Acc.'};
OrderVec=[3 4 1 2 ];
iSet=2;
[gParam, IDParam]=findgroups(Param);
[gFilter, IDFilter]=findgroups(Filter);
[gPerf, IDPerf]=findgroups(PerfMetric);
clear results p tbl
%%
for iPerf=1:4
    
    temp=x(gPerf==iPerf,iSet);
    AnvData=temp/max(temp);
    AnvFilt=Filter(gPerf==iPerf);
    figure
    [p(iSet,iPerf),tbl,stats] = anovan(AnvData,{AnvFilt },...
    'varnames',{'Filter'},'Display','off');
    results((iPerf-1)*3+1:(iPerf-1)*3+3,:) = multcompare(stats,'Dimension',[1],'Alpha',0.05,'CType','bonferroni');
    
end
%%
clear bpsamp
clear labels
i=1;
for iPerf=1:length(IDPerf)
    for iFilter=1:length(IDFilter)  
        iOrd=OrderVec(iPerf);
        temp=x(gPerf==iOrd,iSet);
        AnvData=temp/max(temp);
        Normx(gPerf==iOrd)=AnvData;
        bpsamp(:,i)= Normx(gPerf==iOrd & gFilter==iFilter);
        labels(1:2,i)=[IDPerf(iOrd) IDFilter(iFilter)];
        labels(2,i)=IDFilter(iFilter);
        i=i+1;
    end
end
%%
a=0.05;
mn=mean(bpsamp);
clear CI
for i=1:length(mn)
    CI(:,i)=ConfInt(bpsamp(:,i),a);
end
err=mn-CI(1,:);
mn=reshape(mn,3,[]);
err=reshape(err,3,[]);

xaxis=(1:length(mn));

figure(1)
subplot(2,1,1)
brr=bar(xaxis,mn,'LineWidth',2);
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
yt = get(gca, 'YTick');

grid on

initDist=.13;
incDist=0.1935;
for iText=1:length(IDPerf)+1
    
%     text(iText*3-1,-max(yt)/4,IDPerf(iText),'FontSize', 10,'HorizontalAlignment', 'center');
    annotation('line',[initDist+incDist*(iText-1) initDist+incDist*(iText-1)],[0.95 .55],'LineWidth',1)

end
set(gca, 'XTickLabel',TickLabels)
% ylabel({'(a)'; 'Normalized Performance '})
ylabel({'Normalized Performance '})


Group=[ 1 1; 1 1; 1 1; 2 2 ;2 2; 2 2; 3 3; 3 3 ;4 4; 4 4]; %columns
Bar=  [ 2 3; 1 3; 1 2; 1 2 ;2 3; 1 3; 1 3; 2 3 ;2 3; 1 3]; %rows
[IndNum,~]=size(Group);
BaseLvl=1.30;
IndMarg=0.09;

for iInd=1:IndNum
    
    Bar1=xErr(Bar(iInd,1),Group(iInd,1));
    Bar2=xErr(Bar(iInd,2),Group(iInd,2));
    sigline([Bar1 Bar2],[],BaseLvl+  .1-(Bar(iInd,1)+Bar(iInd,2))*IndMarg*BaseLvl)
    
end

set(gca,'FontSize',12);
legend(IDFilter)
%%
model_series = [10 40 50 60; 20 50 60 70; 30 60 80 90]; 
model_error = [1 4 8 6; 2 5 9 12; 3 6 10 13]; 
b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',model_series,model_error,'k','linestyle','none');
hold off
%%
x = 1:13;
data = [37.6 24.5 14.6 18.1 19.5 8.1 28.5 7.9 3.3 4.1 7.9 1.9 4.3]';
errhigh = [2.1 4.4 0.4 3.3 2.5 0.4 1.6 0.8 0.6 0.8 2.2 0.9 1.5];
errlow  = [4.4 2.4 2.3 0.5 1.6 1.5 4.5 1.5 0.4 1.2 1.3 0.8 1.9];
bar(x,data)                

hold on

er = errorbar(x,data,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

%% Bar plots

a=0.05;
x = pi/10:2*pi;
y = sin(x);
e = std(y)*ones(size(x));
bar(x,y)
hold on
errorbar(x,y,e,'rx')

%%

rng default  % For reproducibility
x1 = normrnd(5,1,100,1);
x2 = normrnd(6,1,100,1);
figure
boxplot([x1,x2],'Notch','on','Labels',{'mu = 5','mu = 6'})
title('Compare Random Data from Different Distributions')
cust_colr = [0, 0.5, 1
    0.60156, 0.80078, 0.19531
    0.5, 0, 0];
h = findobj(gca,'Tag','Box');
 for j=1:length(h)
     patch(get(h(j),'XData'),get(h(j),'YData'),cust_colr(j,:));
 end
 
 %%
 
 % Generate random data
 X = rand(10);
 % Create a new figure and draw a box plot
 figure;
 boxplot(X,'Notch', 'on')
 % Overlay the mean as green diamonds
 hold on
 plot(mean(X), 'dg')
 hold off
%%

y = [0.1423 0.3203; 0.1232 0.1325; 0.1297 0.1302];
hBar = bar(y)
Labels = {'first', 'second', 'third'};
set(gca, 'XTick', 1:6, 'XTickLabel', Labels);
colormap(gray)
ctr2 = bsxfun(@plus, hBar(2).XData, [hBar(2).XOffset]');
hold on
plot(ctr2(1:2), [1 1]*y(1,2)*1.1, '-k', 'LineWidth',2)
plot(mean(ctr2(1:2)), y(1,2)*1.15, '*k')
hold off
%%
X = rand(10);
% Create a new figure and draw a box plot
figure;
% boxplot(X,'Notch', 'on')
  bar(rand(5,1))
sigline([2,4],[],0.9)



function CI=ConfInt(x,a)
%Calculates conf interval of x with alpha equal to a 
a1= a/2;
a2=1-a/2;
SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([a1  a2],length(x)-1);      % T-Score
CI = mean(x) + ts*SEM; 

end

function y=gety(h)
    %Returns the largest single value of ydata in a given graphic handle
    %h= figure,axes,line. Note that y=h if h is not a graphics
    if isgraphics(h) 
        switch(get(h,'type'))
            case {'line','hggroup','patch'},
                y=max(get(h,'ydata'));
                return;
            otherwise
                ys=[];
                hs=get(h,'children');
                for n=1:length(hs)
                    ys=[ys,gety(hs(n))];
                end
                y=max(ys(:));
        end
    else
        y=h;
    end
end

