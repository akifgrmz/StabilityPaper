% new sysn analysis jan20.2024

%1. Preliminaries
clear lgTau TrackingError EstEffort TrueEffort NonPareticAngle PareticAngle StdTe StdEE EffortCorr
SimTime=10;
TimeRange=[2 SimTime];
TimeRange2=[1 SimTime];
NumPoints=100;
Target=0.5;
TrackingError=zeros(1,NumPoints);
EstEffort=zeros(1,NumPoints);
TrueEffort=zeros(1,NumPoints);
NonPareticAngle=zeros(1,NumPoints);
PareticAngle=zeros(1,NumPoints);
StdTe=zeros(1,NumPoints);
StdEE=zeros(1,NumPoints);
EffortCorr=zeros(1,NumPoints);

mdlWks = get_param('NSFsim_syn_2023','ModelWorkspace');
assignin(mdlWks,'rngseed',1000);


%% Define baseline parameters 


stim_hand_tau=0.3;
GainNP=5;
Occ_Var=20; % percent variance 
target=0.5;
stim_gain=5;

str1=sprintf('[%0.3f 1]',stim_hand_tau);
set_param('NSFsim_syn_2023/Paretic Hand/Stim Hand Angle/Transfer Fcn','Denominator',str1)
str1=sprintf('%0.2f',GainNP);
set_param('NSFsim_syn_2023/GainNP','Gain',str1)
str1=sprintf('%0.2f',Occ_Var);
set_param('NSFsim_syn_2023/Paretic Hand/Occlusion/Noise','Value',str1)

str1=sprintf('%d',target);
set_param('NSFsim_syn_2023/DC','Value',str1)

str1=sprintf('%d',stim_gain);
set_param('NSFsim_syn_2023/Gain1','Gain',str1)

GainP=2;
str1=sprintf('%0.2f',GainP);
set_param('NSFsim_syn_2023/GainP','Gain',str1)

RepNum=1;
rng('default');
a=1000;
b=100;
rngseed=round(a+(b-a)*randn(1,RepNum));
mdlWks = get_param('NSFsim_syn_2023','ModelWorkspace');

%% Define changing parameters and run the sim

BaseLineVar=2;
low=BaseLineVar;
high=BaseLineVar*6;
Vars=logspace(log10(low),log10(high),NumPoints);
FiltOrder=[2 1 3 4];
FiltNum=1;
for iRep=1:RepNum
    
    RepLabel=sprintf('Rep_%d',iRep);
    
    for iFilt=1:FiltNum

    FiltLabel=sprintf('Filt_%d',iFilt);
    str1=sprintf('%d',FiltOrder(iFilt));
    set_param('NSFsim_syn_2023/Effort Estimator/FilterNum','Value',str1)

        for iVar=1:length(Vars)

            str1=sprintf('%0.2f',Vars(iVar));
            set_param('NSFsim_syn_2023/Gain1','Gain',str1)

            sim('NSFsim_syn_2023',SimTime);

            tstep=ans.SimulationMetadata.ModelInfo.SolverInfo.FixedStepSize*100;
            S.(RepLabel).TrackingError(iFilt,iVar)=mean(abs(ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),1)));
            S.(RepLabel).StdTe(iFilt,iVar)=std(ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),1)-Target);
            S.(RepLabel).EstEffort(iFilt,iVar)=mean(ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),2));
            S.(RepLabel).TrueEffort(iFilt,iVar)=mean(ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),3));
            S.(RepLabel).NonPareticAngle(iFilt,iVar)=mean(ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),4));
            S.(RepLabel).PareticAngle(iFilt,iVar)=mean(ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),5));
            S.(RepLabel).StdEE(iFilt,iVar)=std(ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),2)...
                -ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),3));
            
            S.(RepLabel).EffortCorr(iFilt,iVar)=corr2(ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),2),...
                ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),3));
            
            S.(RepLabel).(FiltLabel).TrackingError(iVar,:)=ans.Outcome((TimeRange2(1)/tstep:TimeRange2(2)/tstep),1);
            S.(RepLabel).(FiltLabel).EstEffort(iVar,:)=ans.Outcome((TimeRange2(1)/tstep:TimeRange2(2)/tstep),2);
            S.(RepLabel).(FiltLabel).TrueEffort(iVar,:)=ans.Outcome((TimeRange2(1)/tstep:TimeRange2(2)/tstep),3);
        end
    end
end

str1=sprintf('%0.2f',BaseLineVar);
set_param('NSFsim_syn_2023/Gain1','Gain',str1)

%%
save('syn2023_stimgain_test')
%% Frequency response form input output data

NumFilter=3;
NumVar=4;
tstep=0.0001;
fs=1/(tstep*100);
NyqFreq=fs/2;
ChirpFreq=10;
Amp=0.5;
SimTime=10;
x_sampled=zeros(NumFilter,SimTime/tstep/100+1);
y_sampled=zeros(NumFilter,SimTime/tstep/100+1);

    %% StimGain
    
BaseLineVar=2;
low=BaseLineVar;
high=BaseLineVar*6;
Vars=logspace(log10(low),log10(high),NumPoints);
NumPoints=100;


ParLabels(1)="StimGain";
ParValues=logspace(log10(low),log10(high),NumPoints);
AdressString(1)="SysIden_2023/Gain1";
BlockString(1)="Gain";
BaseLineVarString(1)=num2str(BaseLineVar);


for i=1:length(ParValues)
    StimHandValueString(i)=num2str(ParValues(i));
end

O.(ParLabels(1)).ParValues=ParValues;
O.(ParLabels(1)).BaseLineVar=BaseLineVar;
O.(ParLabels(1)).AdressString=AdressString(1);
O.(ParLabels(1)).BaseLineVarString=BaseLineVarString(1);
O.(ParLabels(1)).BlockString=BlockString(1);
O.(ParLabels(1)).StimHandValueString=StimHandValueString;

O.ParLabels=ParLabels;


%%
for iPar=1:1
    
    ParLabel=ParLabels(iPar);
    ParValues=O.(ParLabel).StimHandValueString;

for iFilter=1:1
    
    FiltLabel=sprintf('Filt%d',iFilter);
    str1=sprintf('%d',iFilter);
    set_param('SysIden_NSFsim/Effort Estimator/FilterNum','Value',str1)
    
    for iRepeat=1:1
            
        RepeatLabel=sprintf('Rep_%d',iRepeat);
        assignin(mdlWks,'rngseed',rngseed(iRepeat));
        O.(ParLabel).(FiltLabel).(RepeatLabel).RNGSeed=rngseed(iRepeat);
        
        for iValue=1:NumPoints

            set_param(O.(ParLabel).AdressString,O.(ParLabel).BlockString,O.(ParLabel).ParValues);
            
            sim('SysIden_NSFsim',SimTime);

            O.(ParLabel).(FiltLabel).(RepeatLabel).x_sampled(iValue,:)=ans.Outcome(:,1);
            O.(ParLabel).(FiltLabel).(RepeatLabel).y_sampled(iValue,:)=ans.Outcome(:,5);
        
        end
    end
end
end

