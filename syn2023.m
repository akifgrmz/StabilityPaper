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

mdlWks = get_param('NSFsim_syn','ModelWorkspace');
assignin(mdlWks,'rngseed',1000);

RepNum=10;
rng('default');
a=1000;
b=100;
rngseed=round(a+(b-a)*randn(1,10));
mdlWks = get_param('NSFsim_syn','ModelWorkspace');


% Define baseline parameters 
low=1;
high=4;
wc=logspace(log10(low),log10(high),NumPoints);

str1=sprintf('[%0.3f 1]',0.3);
set_param('NSFsim_syn/Paretic Hand/Stim Hand Angle/Transfer Fcn','Denominator',str1)
str1=sprintf('%0.2f',10);
set_param('NSFsim_syn/GainNP','Gain',str1)
str1=sprintf('%0.2f',20);
set_param('NSFsim_syn/Paretic Hand/Occlusion/Noise','Value',str1)



%% Define changing parameters and run the sim
%% pary hand
% 2. Variable: paretic hand opening time constant
   % vary Paretic hand

low=1;
high=4;
wc=logspace(log10(low),log10(high),NumPoints);


str1=sprintf('[%0.3f 1]',0.3);
set_param('NSFsim_syn/Paretic Hand/Stim Hand Angle/Transfer Fcn','Denominator',str1)
str1=sprintf('%0.2f',10);
set_param('NSFsim_syn/GainNP','Gain',str1)
str1=sprintf('%0.2f',20);
set_param('NSFsim_syn/Paretic Hand/Occlusion/Noise','Value',str1)
 
% low=0.0016;
% high=0.016;
% str1=sprintf('[%0.4f 1]',low);
% str2=sprintf('[%0.4f 1]',high);
% str3=sprintf('[%0.4f 0]',high);
% set_param('NSFsim_syn/Paretic Hand/vEMG/Transfer Fcn1','Denominator',str1)
% set_param('NSFsim_syn/Paretic Hand/vEMG/Transfer Fcn2','Denominator',str2)
% set_param('NSFsim_syn/Paretic Hand/vEMG/Transfer Fcn2','Numerator',str3)

for iRep=1:RepNum
    
    RepLabel=sprintf('Rep_%d',iRep);
%     assignin(mdlWks,'rngseed',rngseed(iRep));
    
    for FilterNum=1:3

    str1=sprintf('%d',FilterNum);
    FiltLabel=sprintf('Filt_%d',FilterNum);

    set_param('NSFsim_syn/Effort Estimator/FilterNum','Value',str1)

        for iwc=1:length(wc)

            str1=sprintf('%0.2f',wc(iwc));
            set_param('NSFsim_syn/GainP','Gain',str1)

            sim('NSFsim_syn',SimTime);

            tstep=ans.SimulationMetadata.ModelInfo.SolverInfo.FixedStepSize*100;
            S.(RepLabel).TrackingError(FilterNum,iwc)=mean(abs(ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),1)));
            S.(RepLabel).StdTe(FilterNum,iwc)=std(ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),1)-Target);
            S.(RepLabel).EstEffort(FilterNum,iwc)=mean(ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),2));
            S.(RepLabel).TrueEffort(FilterNum,iwc)=mean(ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),3));
            S.(RepLabel).NonPareticAngle(FilterNum,iwc)=mean(ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),4));
            S.(RepLabel).PareticAngle(FilterNum,iwc)=mean(ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),5));
            S.(RepLabel).StdEE(FilterNum,iwc)=std(ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),2)...
                -ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),3));

            S.(RepLabel).EffortCorr(FilterNum,iwc)=corr2(ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),2),...
                ans.Outcome((TimeRange(1)/tstep:TimeRange(2)/tstep),3));
            
            S.(RepLabel).(FiltLabel).TrackingError(iwc,:)=ans.Outcome((TimeRange2(1)/tstep:TimeRange2(2)/tstep),1);
            S.(RepLabel).(FiltLabel).EstEffort(iwc,:)=ans.Outcome((TimeRange2(1)/tstep:TimeRange2(2)/tstep),2);
            S.(RepLabel).(FiltLabel).TrueEffort(iwc,:)=ans.Outcome((TimeRange2(1)/tstep:TimeRange2(2)/tstep),3);
            

        end

    end
end

str1=sprintf('%0.2f',2);
set_param('NSFsim_syn/GainP','Gain',str1)

save('syn_3filter2_PhandAllReps2')

