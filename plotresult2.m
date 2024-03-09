clear all
SimTime=10;
NSFsim_syn_plots
NumPoints=100;
iwc=50; %index for baseline value
FilterNum=1; %2 for GS
str1=sprintf('%d',FilterNum);
set_param('NSFsim_syn_plots/Effort Estimator/FilterNum','Value',str1)

low=1;
high=4;
wc=logspace(log10(low),log10(high),NumPoints);
str1=sprintf('%0.2f',wc(iwc));
set_param('NSFsim_syn_plots/GainP','Gain',str1)

low=0.3/2;
high=0.3*2;
wc=logspace(log10(low),log10(high),NumPoints);
str1=sprintf('[%0.3f 1]',wc(iwc));
set_param('NSFsim_syn_plots/Paretic Hand/Stim Hand Angle/Transfer Fcn','Denominator',str1)

low=20/4;
high=20*4;
wc=logspace(log10(low),log10(high),NumPoints);
str1=sprintf('%0.2f',wc(iwc));
set_param('NSFsim_syn_plots/Paretic Hand/Occlusion/Noise','Value',str1)

low=0.01;
high=.2;
wc=logspace(log10(low),log10(high),NumPoints);
str1=sprintf('%0.3f',wc(iwc));
set_param('NSFsim_syn_plots/Paretic Hand/mWave/GainVar','Gain',str1)

low=10;
high=30;
stimfreqs=logspace(log10(low),log10(high),NumPoints);
tstep=0.0001;
wcc=1/tstep./stimfreqs;
a=round(wcc)*tstep;
% wc=1./a;
wc=round(stimfreqs,4);
mdlWks = get_param('NSFsim_syn_plots','ModelWorkspace');
assignin(mdlWks,'stimfreq',wc(iwc));

sim('NSFsim_syn_plots',SimTime);

tstep=ans.SimulationMetadata.ModelInfo.SolverInfo.FixedStepSize*100;

%% 
UnFiltEMG=ans.plots(:,11);
mWaves=ans.plots(:,12);
TimeVec=ans.plots(:,9);

TimeRange=[6 6.29];
Ind=TimeVec>=TimeRange(1) & TimeVec<=TimeRange(2);

figure
subplot(2,1,1)
plot(TimeVec(Ind),UnFiltEMG(Ind),'LineWidth',2,'Color','k');
% legend( {'Target Level', 'SPHO','VPHO','Total Hand Opening'})
% xlabel('Time (s)')
title('Generated EMG (vEMG+M-Waves)')
ylabel('Voltage (a.u)')
set(gca,'FontSize',12)
xlabel('Time (s)')

% set(gca,'XTick',[])
grid
subplot(2,1,2)
plot(TimeVec(Ind),mWaves(Ind),'LineWidth',2,'Color','k');
xlabel('Time (s)')
ylabel('Voltage (a.u)')
title('Generated M-Waves')
% legend( {'Unfiltered EMG', 'M-Waves'})
set(gca,'FontSize',12)
grid
%%
seed=1001;
mdlWks = get_param('NSFsim_syn_plots','ModelWorkspace');
assignin(mdlWks,'rngseed',seed);

FilterNum=1; %2 for GS
str1=sprintf('%d',FilterNum);
set_param('NSFsim_syn_plots/Effort Estimator/FilterNum','Value',str1)

val=1/0.7;
str1=sprintf('%0.2f',val);
set_param('NSFsim_syn_plots/GainP','Gain',str1)

val=0.2;
str1=sprintf('[%0.3f 1]',val);
set_param('NSFsim_syn_plots/Paretic Hand/Stim Hand Angle/Transfer Fcn','Denominator',str1)

sim('NSFsim_syn_plots',SimTime);

Target=ans.plots(:,1);
StimHand=ans.plots(:,2);
VoliHand=ans.plots(:,3);
TotalHand=ans.plots(:,4);
TrueEff=ans.plots(:,10);
CombEff=ans.plots(:,6);
GSEff=ans.plots(:,7);
BlankingEff=ans.plots(:,8);
StimInt=ans.plots(:,13);
TimeVec=ans.plots(:,9);
TimeRange=[0.5 5];
Ind=TimeVec>=TimeRange(1) & TimeVec<=TimeRange(2);

figure(22)
subplot(3,1,1)
plot(TimeVec(Ind), TotalHand(Ind),'LineWidth',2,'DisplayName','Participant A');
hold on


subplot(3,1,2)
plot(TimeVec(Ind),StimInt(Ind),'LineWidth',2,'DisplayName','Participant A');
hold on
%-----------------------------------
seed=1002;
mdlWks = get_param('NSFsim_syn_plots','ModelWorkspace');
assignin(mdlWks,'rngseed',seed);

val=1/0.3;
str1=sprintf('%0.2f',val);
set_param('NSFsim_syn_plots/GainP','Gain',str1)

val=0.5;
str1=sprintf('[%0.3f 1]',val);
set_param('NSFsim_syn_plots/Paretic Hand/Stim Hand Angle/Transfer Fcn','Denominator',str1)

sim('NSFsim_syn_plots',SimTime);

Target=ans.plots(:,1);
StimHand=ans.plots(:,2);
VoliHand=ans.plots(:,3);
TotalHand=ans.plots(:,4);
TrueEff=ans.plots(:,10);
CombEff=ans.plots(:,6);
GSEff=ans.plots(:,7);
BlankingEff=ans.plots(:,8);
StimInt=ans.plots(:,13);
TimeVec=ans.plots(:,9);

TimeRange=[0.5 5];
Ind=TimeVec>=TimeRange(1) & TimeVec<=TimeRange(2);

figure(22)
subplot(3,1,1)
plot(TimeVec(Ind), TotalHand(Ind),'LineWidth',2,'DisplayName','Participant B');
subplot(3,1,2)
plot(TimeVec(Ind),StimInt(Ind),'LineWidth',2,'DisplayName','Participant B');

%-----------------------------------
seed=1003;
mdlWks = get_param('NSFsim_syn_plots','ModelWorkspace');
assignin(mdlWks,'rngseed',seed);

val=1/0.7;
str1=sprintf('%0.2f',val);
set_param('NSFsim_syn_plots/GainP','Gain',str1)

val=0.5;
str1=sprintf('[%0.3f 1]',val);
set_param('NSFsim_syn_plots/Paretic Hand/Stim Hand Angle/Transfer Fcn','Denominator',str1)

sim('NSFsim_syn_plots',SimTime);

Target=ans.plots(:,1);
StimHand=ans.plots(:,2);
VoliHand=ans.plots(:,3);
TotalHand=ans.plots(:,4);
TrueEff=ans.plots(:,10);
CombEff=ans.plots(:,6);
GSEff=ans.plots(:,7);
BlankingEff=ans.plots(:,8);
StimInt=ans.plots(:,13);
TimeVec=ans.plots(:,9);

TimeRange=[0.5 5];
Ind=TimeVec>=TimeRange(1) & TimeVec<=TimeRange(2);

figure(22)
subplot(3,1,1)
plot(TimeVec(Ind), TotalHand(Ind),'LineWidth',2,'DisplayName','Participant C');
subplot(3,1,2)
plot(TimeVec(Ind),StimInt(Ind),'LineWidth',2,'DisplayName','Participant C');
%-----------------------------------

seed=1003;
mdlWks = get_param('NSFsim_syn_plots','ModelWorkspace');
assignin(mdlWks,'rngseed',seed);


iwc=50; %index for baseline value
FilterNum=5; %2 for GS
str1=sprintf('%d',FilterNum);
mWavePercent=0.8;
set_param('NSFsim_syn_plots/Effort Estimator/FilterNum','Value',str1)
str1=sprintf('%f',mWavePercent);
set_param('NSFsim_syn_plots/Effort Estimator/HypFilter/Constant','Value',str1);

low=1;
high=4;
wc=logspace(log10(low),log10(high),NumPoints);
str1=sprintf('%0.2f',wc(iwc));
set_param('NSFsim_syn_plots/GainP','Gain',str1)


low=0.3/2;
high=0.3*2;
wc=logspace(log10(low),log10(high),NumPoints);
str1=sprintf('[%0.3f 1]',wc(iwc));
set_param('NSFsim_syn_plots/Paretic Hand/Stim Hand Angle/Transfer Fcn','Denominator',str1)

seed=1004;
mdlWks = get_param('NSFsim_syn_plots','ModelWorkspace');
assignin(mdlWks,'rngseed',seed);

sim('NSFsim_syn_plots',SimTime);

Target=ans.plots(:,1);
StimHand=ans.plots(:,2);
VoliHand=ans.plots(:,3);
TotalHand=ans.plots(:,4);
TrueEff=ans.plots(:,10);
CombEff=ans.plots(:,6);
GSEff=ans.plots(:,7);
BlankingEff=ans.plots(:,8);
StimInt=ans.plots(:,13);
TimeVec=ans.plots(:,9);

TimeRange=[0.5 5];
Ind=TimeVec>=TimeRange(1) & TimeVec<=TimeRange(2);

figure(22)
subplot(3,1,1)
plot(TimeVec(Ind), TotalHand(Ind),'LineWidth',2,'DisplayName','Case A');
legend()
subplot(3,1,2)
plot(TimeVec(Ind),StimInt(Ind),'LineWidth',2,'DisplayName','Case A');
legend()
%-----------------------------------
iwc=50; %index for baseline value
FilterNum=5; %2 for GS
str1=sprintf('%d',FilterNum);
mWavePercent=0.2;
set_param('NSFsim_syn_plots/Effort Estimator/FilterNum','Value',str1)
str1=sprintf('%f',mWavePercent);
set_param('NSFsim_syn_plots/Effort Estimator/HypFilter/Constant','Value',str1);

seed=1005;
mdlWks = get_param('NSFsim_syn_plots','ModelWorkspace');
assignin(mdlWks,'rngseed',seed);

sim('NSFsim_syn_plots',SimTime);

Target=ans.plots(:,1);
StimHand=ans.plots(:,2);
VoliHand=ans.plots(:,3);
TotalHand=ans.plots(:,4);
TrueEff=ans.plots(:,10);
CombEff=ans.plots(:,6);
GSEff=ans.plots(:,7);
BlankingEff=ans.plots(:,8);
StimInt=ans.plots(:,13);
TimeVec=ans.plots(:,9);

TimeRange=[0.5 5];
Ind=TimeVec>=TimeRange(1) & TimeVec<=TimeRange(2);

figure(22)
subplot(3,1,1)
plot(TimeVec(Ind), TotalHand(Ind),'LineWidth',2,'DisplayName','Case B');
xlabel('Time (s)')
ylabel('Hand Opening Angle (a.u)') 
set(gca,'FontSize',12)
yt = get(gca, 'YTick');
h=text(0.1,max(yt)/2,'(a)','FontSize', 12);
set(h,'Rotation',90);
grid on

subplot(3,1,2)
plot(TimeVec(Ind),StimInt(Ind),'LineWidth',2,'DisplayName','Case B');
legend()
xlabel('Time (s)')
ylabel('Stimulation Intensity (a.u)')
set(gca,'FontSize',12)
yt = get(gca, 'YTick');
h=text(0.1,max(yt)/2,'(b)','FontSize', 12);
set(h,'Rotation',90);
grid on


subplot(3,1,3)
plot(TimeVec(Ind),0.8*CombEff(Ind),'LineWidth',2,'DisplayName','Comb');
hold on
plot(TimeVec(Ind),GSEff(Ind),'LineWidth',2,'DisplayName','GS');
plot(TimeVec(Ind),BlankingEff(Ind),'LineWidth',2,'DisplayName','Blanking');
plot(TimeVec(Ind),TrueEff(Ind),'LineWidth',2,'DisplayName','True Effort');
legend()
xlabel('Time (s)')
ylabel('Estimated Effort (a.u)')
set(gca,'FontSize',12)
yt = get(gca, 'YTick');
h=text(0.1,max(yt)/2,'(c)','FontSize', 12);
set(h,'Rotation',90);
grid on

%-----------------------------

seed=1010;
mdlWks = get_param('NSFsim_syn_plots','ModelWorkspace');
assignin(mdlWks,'rngseed',seed);

val=1/0.6;
str1=sprintf('%0.2f',val);
set_param('NSFsim_syn_plots/GainP','Gain',str1)

val=0.3;
str1=sprintf('[%0.3f 1]',val);
set_param('NSFsim_syn_plots/Paretic Hand/Stim Hand Angle/Transfer Fcn','Denominator',str1)

sim('NSFsim_syn_plots',SimTime);

Target=ans.plots(:,1);
StimHand=ans.plots(:,2);
VoliHand=ans.plots(:,3);
TotalHand=ans.plots(:,4);
TrueEff=ans.plots(:,10);
CombEff=ans.plots(:,6);
GSEff=ans.plots(:,7);
BlankingEff=ans.plots(:,8);
StimInt=ans.plots(:,13);
TimeVec=ans.plots(:,9);

TimeRange=[0.5 5];
Ind=TimeVec>=TimeRange(1) & TimeVec<=TimeRange(2);

figure(22)
subplot(3,1,1)
plot(TimeVec(Ind), TotalHand(Ind),'LineWidth',2,'DisplayName','Case C');
plot(TimeVec(Ind), Target(Ind),'LineWidth',2,'DisplayName','Target Level');

subplot(3,1,2)
plot(TimeVec(Ind),StimInt(Ind),'LineWidth',2,'DisplayName','Case C');
%-----------------------------------



