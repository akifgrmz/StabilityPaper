%% system inde. for NSF project with filters below:
%1-Comb
%2-GS
%3-Banking
%Input:tracking error
%output: Hand angle

% clear lgTau TrackingError EstEffort TrueEffort NonPareticAngle PareticAngle StdTe StdEE EffortCorr
SimTime=10;
TimeRange=[2 SimTime];
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
mdlWks = get_param('SysIden_2023_ver1','ModelWorkspace');
assignin(mdlWks,'rngseed',1000);


NumFilter=3;
NumVar=4;
tstep=0.0001;
fs=1/(tstep*100);
NyqFreq=fs/2;
ChirpFreq=5;
Amp=0.5;
SimTime=10;
x_sampled=zeros(NumFilter,SimTime/tstep/100+1);
y_sampled=zeros(NumFilter,SimTime/tstep/100+1);
RepeatNum=10;
% O=struct;
rngseed=1000:1000+RepeatNum;
a=0.1;

O.ParamLabels=["Target","Paresis"];
O.AdressString=["SysIden_2023_ver1/DC","SysIden_2023_ver1/Saturation"];
O.BlockParam=["Value","UpperLimit"];

low=[0.5 0.2];
base=[0.5 1];
high=[1 0.8];
    
for iParam=1:length(O.ParamLabels)
    ParamLabel=O.ParamLabels(iParam);
    
    O.(ParamLabel).ParamValues=logspace(log10(low(iParam)),log10(high(iParam)),NumPoints);
    O.(ParamLabel).AdressString=O.AdressString(iParam);
    O.(ParamLabel).BaseValue=base(iParam);
%     O.(ParamLabel).BaseString=num2str(base(iParam));
    O.(ParamLabel).BlockParam=O.BlockParam(iParam);
    
end

% %%
%     % Paresis 
% low=0.2;
% base=1;
% high=0.8;
% Values=logspace(log10(low),log10(high),NumPoints);
% AdressString="SysIden_2023/Saturation";
% BlockString="Upper Limit";
% NominalString=num2str(base);
% ParamLabel="Paresis";
% O.(ParamLabel).ParamValues=Values;
% O.(ParamLabel).AdressString=AdressString;
% O.(ParamLabel).BaseValue=base;
% 
% 
% % for i=1:length(StimHandValues)
% % %     str1=sprintf('[%0.4f 1]',NonParHandValues(i));
% %     StimHandValueString{i}=sprintf('[%0.8f 1]',StimHandValues(i));
% % end
% % 
% % ValueStringMat=[ParHandValueString;StimHandValueString;KValuesString;mWaveVarString;StimFreqVarString;TargetVarString];



Param1Label=O.ParamLabels(1);
Param1Address= O.(Param1Label).AdressString;
for iVal1=1:length(O.(Param1Label).ParamValues)
    Val1=O.(Param1Label).ParamValues(iVal1);
    Block1Param=O.(Param1Label).BlockParam;
    
    set_param(Param1Address,Block1Param,num2str(Val1)) 

    Param2Label=O.ParamLabels(2);
    Param2Address= O.(Param2Label).AdressString;
    for iVal2=1:length(O.(Param2Label).ParamValues)
        Val2=O.(Param2Label).ParamValues(iVal2);
        Block2Param=O.(Param2Label).BlockParam;
        
        set_param(Param2Address,Block2Param,num2str(Val2)) 

        for iFilter=1:NumFilter

            FiltLabel=sprintf('Filt%d',iFilter);
            str1=sprintf('%d',iFilter);
            set_param('SysIden_2023_ver1/Effort Estimator/FilterNum','Value',str1)

            for iRepeat=1:10

                RepeatLabel=sprintf('Rep_%d',iRepeat);
                assignin(mdlWks,'rngseed',rngseed(iRepeat));
                O.Sims.(FiltLabel).(RepeatLabel).RNGSeed=rngseed(iRepeat);

                sim('SysIden_2023_ver1',SimTime);

                O.Sims.(FiltLabel).(RepeatLabel).x_sampled(iVal1,iVal2,:)=ans.Outcome(:,1);
                O.Sims.(FiltLabel).(RepeatLabel).y_sampled(iVal1,iVal2,:)=ans.Outcome(:,5);

            end
        end
    end
end
filename=sprintf('SysIden_ver1');
save(filename);
%%
set_param(AdressString{iPar},BlockString{iPar},NominalString{iPar}) 
assignin(mdlWks,'rngseed',1000);
filename=sprintf('worksim2_ipar%d',iPar);
save(filename);
 
