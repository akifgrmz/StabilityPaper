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
SysIden_orig_2023
mdlWks = get_param('SysIden_orig_2023','ModelWorkspace');
assignin(mdlWks,'rngseed',1000);


FilterNum=2;
tstep=0.0001;
fs=1/(tstep*100);
NyqFreq=fs/2;
ChirpFreq=10;
% Amp=0.5;
SimTime=10;
x_sampled=zeros(FilterNum,SimTime/tstep/100+1);
y_sampled=zeros(FilterNum,SimTime/tstep/100+1);
RepNum=1;
% O=struct;
rngseed=1000:1000+RepNum;
a=0.1;

O.ParamLabels=["Target","Paresis"];
O.AdressString=["SysIden_orig_2023/DC","SysIden_orig_2023/Saturation"];
O.BlockParam=["Value","UpperLimit"];

low=[0.5 0.2];
base=[0.5 1];
high=[1 0.8];
    NumPoints=[100 10];
for iParam=1:length(O.ParamLabels)
    ParamLabel=O.ParamLabels(iParam);
    
    O.(ParamLabel).ParamValues=logspace(log10(low(iParam)),log10(high(iParam)),NumPoints(iParam));
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
Param2Label=O.ParamLabels(2);

RepNum=1;
Val1Num=length(O.(Param1Label).ParamValues);
Val2Num=length(O.(Param2Label).ParamValues);
NumofRun=RepNum*FilterNum*Val1Num*Val2Num;
Val2Num=length(O.(Param2Label).ParamValues);
Count=1;

Param1Address= O.(Param1Label).AdressString;
for iVal1=1:Val1Num
    Val1=O.(Param1Label).ParamValues(iVal1);
    Block1Param=O.(Param1Label).BlockParam;
    
    set_param(Param1Address,Block1Param,num2str(Val1)) 

    Param2Address= O.(Param2Label).AdressString;
    for iVal2=1:Val2Num
        Val2=O.(Param2Label).ParamValues(iVal2);
        Block2Param=O.(Param2Label).BlockParam;
        
        set_param(Param2Address,Block2Param,num2str(Val2)) 

        for iFilter=1:FilterNum

            FiltLabel=sprintf('Filt_%d',iFilter);
            str1=sprintf('%d',iFilter+1);
            set_param('SysIden_orig_2023/Effort Estimator/FilterNum','Value',str1)

            for iRepeat=1:RepNum

                RepeatLabel=sprintf('Rep_%d',iRepeat);
                assignin(mdlWks,'rngseed',rngseed(iRepeat));
                O.Sims.(FiltLabel).(RepeatLabel).RNGSeed=rngseed(iRepeat);

                sim('SysIden_orig_2023',SimTime);

                O.Sims.(FiltLabel).(RepeatLabel).x_sampled(iVal1,iVal2,:)=ans.Outcome(:,1);
                O.Sims.(FiltLabel).(RepeatLabel).y_sampled(iVal1,iVal2,:)=ans.Outcome(:,5);
                Count=Count+1;
                sprintf('%.2f%% Complete',Count/NumofRun*100)
            end
        end
    end
end

filename=sprintf('SysIden_ver2');
save(filename);

%% Freq response
ParNum=1;
load('SysIden_ver1');
ConfIntLabel={ 'Mean','Lb','Ub'};
FiltTitle={'Comb Filter', 'GS Filter' 'Blanking Method'};
NominalInd=50;  %assumption



% Reordering
      
for iFilt=1:FilterNum

    FiltLabel=sprintf('Filt_%d',iFilt);
%         FiltLabel2=sprintf('Filt_%d',iFilt);

    for iRep=1:RepNum

        RepeatLabel=sprintf('Rep_%d',iRep);

        for iVal1=1:Val1Num
            Val1Label=sprintf('Val_%d',iVal1);

            for iVal2=1:Val2Num

                Val2Label=sprintf('Val_%d',iVal2);

                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).x_sampled=squeeze(O.Sims.(FiltLabel).(RepeatLabel).x_sampled(iVal1,iVal2,:));
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).y_sampled=squeeze(O.Sims.(FiltLabel).(RepeatLabel).y_sampled(iVal1,iVal2,:));
            end
        end
    end
end



Ts=(1/fs);
tol = 1e-4;
IoDelay=0.1;
DelayNumSamp=0.1/Ts;
NumFFT=(SimTime/Ts)+1-DelayNumSamp;  % removing the delay 100 ms from the response
f_full=linspace(-1/2/Ts,1/2/Ts,NumFFT);
NumUniquePts = ceil(length(f_full)/2);
f_onesided=f_full(NumUniquePts:end);
TimeVec=0:Ts:SimTime;
TimeIndice=DelayNumSamp+1:length(TimeVec); %this is for ignoring the first 100ms due to the fixed delay

for iFilter=1:FilterNum

    FiltLabel=sprintf('Filt_%d',iFilter);

    for iVal1=1:Val1Num
        Val1Label=sprintf('Val_%d',iVal1);
        
        for iVal2=1:Val2Num
            Val2Label=sprintf('Val_%d',iVal2);
            
            for iRep=1:RepNum
                RepeatLabel=sprintf('Rep_%d',iRep);

                x_v=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).x_sampled;
                y_v=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).y_sampled;

                [FreqRespWelch, f_welch]=tfestimate(x_v,y_v,200,[100],[f_onesided],fs);
        %             [FreqRespWelch, f_welch]=tfestimate(x_v,y_v,200,[100],[f_onesided],fs/DownSampleRate);

                FT_L=fft(y_v,NumFFT)./fft(x_v,NumFFT);
                FT_L(abs(FT_L) < tol) = 0;
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).FT=FT_L;
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).FreqRespOneSided=FT_L(1:NumUniquePts);
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).MagAllRepsOneSided=abs(FT_L(1:NumUniquePts));
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).PhaseAllRepsOneSided=angle(FT_L(1:NumUniquePts));

                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).FreqRespWelch=FreqRespWelch;
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).FreqRespWelchOneSided=FreqRespWelch(1:NumUniquePts);
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).MagWelchAllRepsOneSided=abs(FreqRespWelch(1:NumUniquePts));
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).PhaseWelchAllRepsOneSided=angle(FreqRespWelch(1:NumUniquePts));

            end

            a=0.05;
            for i=1:fs*SimTime

                for iRep=1:RepNum

                    RepeatLabel=sprintf('Rep_%d',iRep);

                    x_v(iRep)=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).x_sampled(i);
                    y_v(iRep)=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).y_sampled(i);
                end

                CI(:,i)=ConfInt(x_v,a);
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgTimeIn(1,i)=mean(x_v);
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgTimeIn(2:3,i)=CI(:,i);

                CI(:,i)=ConfInt(y_v,a);
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgTimeOut(1,i)=mean(y_v);
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgTimeOut(2:3,i)=CI(:,i);

             end

            for i=1:NumFFT

                for iRep=1:RepNum

                    RepeatLabel=sprintf('Rep_%d',iRep);
                    x(iRep)= O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).FT(i);

                end

        %         x=O.(ParLabel).(FiltLabel).(ValLabel).FT(:,i);

                x(abs(x) < tol) = 0; % removing insignificant valuese for better accu in phase

                CI(:,i)=ConfInt(x,a);
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgFreqResp(1,i)=mean(x);
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgFreqResp(2:3,i)=CI(:,i);

                CIMag(:,i)=ConfInt(abs(x),a);
                CIAngle(:,i)=ConfInt(angle(x),a);

                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgMag(1,i)=mean(abs(x));
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgMag(2:3,i)=(CIMag(:,i));
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgPhase(1,i)=mean(angle(x));
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgPhase(2:3,i)=CIAngle(:,i);

            end

            O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgMagOneSided= O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgMag(:,1:NumUniquePts);
            O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgPhaseOneSided= O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgPhase(:,1:NumUniquePts);
            O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgFreqRespOneSided= O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgFreqResp(:,1:NumUniquePts);

            O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).f_full=f_full;
            O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).NumUniquePts=NumUniquePts;
            O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).f_onesided=f_onesided;
        
        end
    end
end

%% plotting

    
RepeatLabel=sprintf('Rep_%d',1);

for iFilt=1:2

    FiltLabel=sprintf('Filt_%d',iFilt);

    for iVal1=20:20

        Val1Label=sprintf('Val_%d',iVal1);
        
        for iVal2=7:7

            Val2Label=sprintf('Val_%d',iVal2);

            MagL=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgMagOneSided(:,2:end);
            PhaseL=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgPhaseOneSided(:,2:end);
            f_onesided=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).f_onesided(2:end);
            MagWelchAllRepsOneSided=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).MagWelchAllRepsOneSided(2:end);

    %             PhaseAllRepsOneSided=O.(ParLabel).(FiltLabel).(ValLabel).PhaseAllRepsOneSided(1,:);
    %             MagAllRepsOneSided=O.(ParLabel).(FiltLabel).(ValLabel).MagAllRepsOneSided(1,:);

            plotrange=1:NumUniquePts;
            figure(iVal1)
            subplot(3,1,iFilt)
            semilogx(f_onesided,20*log10((MagWelchAllRepsOneSided)),'DisplayName','Welch')
            hold 
            semilogx(f_onesided,20*log10((MagL)),'DisplayName','FFT')

            xlabel(' Freq. (Hz) ')
            ylabel('Mag. (dB)')
            title(FiltTitle{iFilt})

            ylim([-50 10 ])
%             xlim([0 10])
            legend
%             legend(ConfIntLabel)
%             subplot(2,1,2)
%             semilogx(f_onesided,PhaseL)
%             xlabel(' Freq. (Hz) ')
%             ylabel('Phase (rad)')
%             xlim([0 NyqFreq])
        end
    end
end



%% Estimate TF
global Estnp Estnz FreqResp FreqVec Ts StabThres minfThres IterLim

FreqRange=1.5; % Hz
FreqIndice=2:NumUniquePts*FreqRange/NyqFreq;
ConfIntLabel={ 'FitAvg','FitLb','FitUb'};
lgd={'Avg Fitted', 'Avg ','Lb Fitted','Lb', 'Ub Fitted', 'Ub '};
ColorLabel={'b','r','g'};
FiltTitle={'Comb Filter', 'GS Filter' 'Blanking Method'};
Ts=1/fs;

%run for the nominal values 

for iFilt=1:FilterNum
    
    FiltLabel=sprintf('Filt_%d',iFilt);
    
    for iVal1=1:Val1Num
        
        Val1Label=sprintf('Val_%d',iVal1);
        for iVal2=1:Val2Num
        
            Val2Label=sprintf('Val_%d',iVal2);

            FreqVec=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).f_onesided(FreqIndice);

            for iRep=1:RepNum % fitting for all the repeats 

                RepeatLabel=sprintf('Rep_%d',iRep);

                Phase=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).PhaseAllRepsOneSided(FreqIndice);
                Mag=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).MagAllRepsOneSided(FreqIndice);
                FreqResp=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).FreqRespOneSided(FreqIndice);

                Estnp=1;
                Estnz=0;

                minfThres=2;
                StabThres=1;
                IterLim=20;

                lb= [0     .5];
                ub= [0.1    .6];
                x0=0.01*ones(1,(Estnp+Estnz+1));
                [FittedFreqResp SysTf]= sysiden_basic(Mag,x0,lb,ub,FreqVec,Estnp,Estnz,Ts);
                SysFreq = frd(FreqResp,FreqVec,'FrequencyUnit','Hz','Ts',Ts);

                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).nz = Estnz;
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).lb = lb;
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).ub = ub;
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).SysFreq=SysFreq;
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).FreqVec=FreqVec;

                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).SysTf=SysTf;
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).FittedFreqResp=FittedFreqResp';

            end


            for iConf=1:3  % fitting for Conf intervals.

    %             RepeatLabel=sprintf('Rep_%d',iRep);
                  ConfLabel=sprintf('ConfIntFit_%d',iConf);

    %             Mag=O.(ParLabel).(FiltLabel).(ValLabel).MagOneSided(i,FreqIndice);
    %             Phase=O.(ParLabel).(FiltLabel).(ValLabel).PhaseOneSided(i,FreqIndice);
    %             FreqResp=O.(ParLabel).(FiltLabel).(ValLabel).AvgFreqRespOneSided(i,FreqIndice);

                Phase=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgPhaseOneSided(iConf,FreqIndice);
                Mag=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgMagOneSided(iConf,FreqIndice);
                FreqResp=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).AvgFreqRespOneSided(iConf,FreqIndice);

                Estnp=1;
                Estnz=0;

                minfThres=2;
                StabThres=1;
                IterLim=20;


                lb= [0     .5];
                ub= [0.1    .6];
                x0=0.01*ones(1,(Estnp+Estnz+1));
    %             x0=Detx0(Estnp+Estnz+1,lb,ub);
                [FittedFreqResp SysTf]= sysiden_basic(Mag,x0,lb,ub,FreqVec,Estnp,Estnz,Ts);
                SysFreq = frd(FreqResp,FreqVec,'FrequencyUnit','Hz','Ts',Ts);
    %             SysTf=tfest(SysFreq,1,0,'Ts',Ts);

                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).(ConfLabel).np = Estnp;
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).(ConfLabel).Mag = Mag;
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).(ConfLabel).Phase = Phase;
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).(ConfLabel).FreqResp = FreqResp;

                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).(ConfLabel).nz = Estnz;
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).(ConfLabel).lb = lb;
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).(ConfLabel).ub = ub;
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).(ConfLabel).SysFreq=SysFreq;
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).FreqVec=FreqVec;

                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).(ConfLabel).SysTf=SysTf;
                O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).(ConfLabel).FittedFreqResp=FittedFreqResp';
            end
        end
    end
end

S.(ParLabel)=O.(ParLabel);



save('sysiden2023_target_ana')

%% plotting

for iPar=1:1
    
    ParLabel=VarLabels{iPar};
    ParLabel=strcat(ParLabel,sprintf('_%d',fs/DownSampleRate));

for iFilt=1:3
    
    FiltLabel=sprintf('Filt_%d',iFilt);
    
    for iVal=50:50
        
        ValLabel=sprintf('Val_%d',iVal);
        
        for iConf=1:3
            
            ConfLabel=sprintf('ConfIntFit_%d',iConf);

            FittedMag=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).(ConfLabel).FittedFreqResp;
            Mag=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).(ConfLabel).Mag;
            
            SysTf=O.(FiltLabel).(RepeatLabel).(Val1Label).(Val2Label).(ConfLabel).SysTf;

            figure(iVal)
            subplot(3,1,iFilt)
            semilogx(FreqVec,20*log10(abs(FittedMag)),"LineWidth",2,"Color",ColorLabel{iConf});
            hold on
            semilogx(FreqVec,20*log10(abs(Mag)),'--',"Color",ColorLabel{iConf});

            legend(lgd)
            title(FiltTitle(iFilt))
            xlabel(' Freq. (Hz) ')
            ylabel('Mag.')
            ylim([-40 40 ])
                      
            figure(100+iVal)
            [re,im,wout,sdre,sdim] = nyquist(SysTf);
            
            re=squeeze(re);
            im=squeeze(im);
            
            subplot(3,1,iFilt)
            plot(re,im)
            hold on 

        end
            
%             figure
%             MagL=O.(ParLabel).(FiltLabel).(ValLabel).MagAllRepsOneSided(1,:);
%             PhaseL=O.(ParLabel).(FiltLabel).(ValLabel).PhaseAllRepsOneSided(1,:);
%             plotrange=1:NumUniquePts;
%             figure(iFilt)
%             subplot(2,1,1)
%             semilogx(f_onesided,20*log10((MagL)))
% %             hold 
% %             semilogx(f_onesided,zeros(1,length(f_onesided)),'k')
%             xlabel(' Freq. (Hz) ')
%             ylabel('Mag. (dB)')
% %             ylim([-60 40 ])
% %             xlim([0 NyqFreq])
%             title(FiltTitle{iFilt})cl
%             legend(ConfIntLabel)
%             subplot(2,1,2)
%             semilogx(f_onesided,PhaseL)
%             xlabel(' Freq. (Hz) ')
%             ylabel('Phase (rad)')
%             xlim([0 NyqFreq])
%             
%             figure
%             nyquist(SysTf)
%             rlocus(SysTf)
%             bode(SysTf)
    end
end
end


function [Wtest Wtf]=sysiden_basic(FreqResp,x0,lb,ub,FreqVector,np,nz,Ts)
% system indentification using fmincon\
global minfThres IterLim FreqVec
minf=10000000;
nIter=0;
A=[];
b=[];
Aeq=[];
beq=[];


while minf>minfThres && nIter<IterLim
    
%     options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
    
    [x minf outt s] = fmincon(@objective_fit2,x0,A,b,Aeq,beq,lb,ub,@stability_z);
    nIter=nIter+1;

    MinFs(nIter,:)=minf;
    xs(nIter,:)=x;
    x0=0.01*ones(1,length(x));
    
end

[minFval I]=min(MinFs);
minFval
x=xs(I,:)
nIter
% [Wtf num den]=ConstructTF(x,nz);
[Wtf num den]=SpecialTF(x,nz);
zrs=roots(num);
pls=roots(den);
num 
den
Wtest=freqresp(Wtf,FreqVec,'Hz');
Wtest=squeeze(Wtest)';

end


function obj=objective_fit2(x)

global FreqVec FreqResp Estnz

% Stf=ConstructTF(x,Estnz);
Stf=SpecialTF(x,Estnz);


S=freqresp(Stf,FreqVec,'Hz');
S=squeeze(S);

obj1=sum((abs(FreqResp)-abs(S')).^2);
% obj2=sum(((angle(FreqResp)-angle(S'))).^2);
obj=obj1;
end

function [TFun,num,den]=SpecialTF(x,z)
global Ts
num=[ x(1)   ] ;
den=[ 1 -x(2)-1 x(2) ];
z;
TFun=tf(num,den,Ts);

end


function [c,ceq]=stability_z(x)
% constraints for roots(den)<1 and roots (num)<1
global Estnz Estnp StabThres
% Estnz=0;

c=zeros(1,length(x));
Estnp=length(x)-Estnz-1;
% [~,num,den]=ConstructTF(x,Estnz);
[~,num,den]=SpecialTF(x,Estnz);

zrs=roots(num);
pls=roots(den);
c(1:length(zrs))=abs(zrs)-StabThres;
c(length(zrs)+1:length(zrs)+length(pls))=abs(pls)-StabThres;
ceq=[];

end
function CI=ConfInt(x,a)
%Calculates conf interval of x with alpha equal to a 
a1= a/2;
a2=1-a/2;
SEM = std(x)/sqrt(length(x));         % Standard Error
ts = tinv([a1  a2],length(x)-1);      % T-Score
CI = mean(x) + ts*SEM; 

end