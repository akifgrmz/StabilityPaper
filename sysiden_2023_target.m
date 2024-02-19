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

mdlWks = get_param('SysIden_orig_2023','ModelWorkspace');
assignin(mdlWks,'rngseed',1000);SysIden_orig_2023


FilterNum=1;
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
O.AdressString=["SysIden_orig_2023/DC","SysIden_orig_2023/Constant"];
O.BlockParam=["Value","Value"];

low=[0.5 0.2];
base=[0.5 1];
high=[1 0.8];

NumPoints=[10 10];

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

Val1Num=length(O.(Param1Label).ParamValues);
Val2Num=length(O.(Param2Label).ParamValues);
NumofRun=RepNum*FilterNum*Val1Num*Val2Num;
Count=1;
%%
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
            str1=sprintf('%d',iFilter);
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

filename=sprintf('SysIden_ver4');
save(filename);

%% Freq response
load('SysIden_ver6');
ConfIntLabel={ 'Mean','Lb','Ub'};
FiltTitle={'Comb Filter', 'GS Filter' 'Blanking Method'};
NominalInd=50;  %assumption
RepNum=10;
Val2Num=10;
Val1Num=10;
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

                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).x_sampled=squeeze(O.Sims.(FiltLabel).(RepeatLabel).x_sampled(iVal1,iVal2,:));
                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).y_sampled=squeeze(O.Sims.(FiltLabel).(RepeatLabel).y_sampled(iVal1,iVal2,:));
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

                x_v=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).x_sampled;
                y_v=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).y_sampled;

                [FreqRespWelch, f_welch]=tfestimate(x_v,y_v,200,[100],[f_onesided],fs);
        %             [FreqRespWelch, f_welch]=tfestimate(x_v,y_v,200,[100],[f_onesided],fs/DownSampleRate);

                FT_L=fft(y_v,NumFFT)./fft(x_v,NumFFT);
                FT_L(abs(FT_L) < tol) = 0;
                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).FT=FT_L;
                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).FreqRespOneSided=FT_L(1:NumUniquePts);
                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).MagAllRepsOneSided=abs(FT_L(1:NumUniquePts));
                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).PhaseAllRepsOneSided=angle(FT_L(1:NumUniquePts));


                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).FreqRespWelch=FreqRespWelch;
                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).FreqRespWelchOneSided=FreqRespWelch(1:NumUniquePts);
                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).MagWelchAllRepsOneSided=abs(FreqRespWelch(1:NumUniquePts));
                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).PhaseWelchAllRepsOneSided=angle(FreqRespWelch(1:NumUniquePts));

            end

            a=0.05;
            for i=1:fs*SimTime

                for iRep=1:RepNum

                    RepeatLabel=sprintf('Rep_%d',iRep);

                    x_v(iRep)=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).x_sampled(i);
                    y_v(iRep)=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).y_sampled(i);
                end

                CI(:,i)=ConfInt(x_v,a);
                O.(FiltLabel).(Val1Label).(Val2Label).AvgTimeIn(1,i)=mean(x_v);
                O.(FiltLabel).(Val1Label).(Val2Label).AvgTimeIn(2:3,i)=CI(:,i);

                CI(:,i)=ConfInt(y_v,a);
                O.(FiltLabel).(Val1Label).(Val2Label).AvgTimeOut(1,i)=mean(y_v);
                O.(FiltLabel).(Val1Label).(Val2Label).AvgTimeOut(2:3,i)=CI(:,i);

             end

            for i=1:NumFFT

                for iRep=1:RepNum

                    RepeatLabel=sprintf('Rep_%d',iRep);
                    x(iRep)= O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).FT(i);

                end

        %         x=O.(ParLabel).(FiltLabel).(ValLabel).FT(:,i);

                x(abs(x) < tol) = 0; % removing insignificant valuese for better accu in phase

                CI(:,i)=ConfInt(x,a);
                O.(FiltLabel).(Val1Label).(Val2Label).AvgFreqResp(1,i)=mean(x);
                O.(FiltLabel).(Val1Label).(Val2Label).AvgFreqResp(2:3,i)=CI(:,i);

                CIMag(:,i)=ConfInt(abs(x),a);
                CIAngle(:,i)=ConfInt(angle(x),a);

                O.(FiltLabel).(Val1Label).(Val2Label).AvgMag(1,i)=mean(abs(x));
                O.(FiltLabel).(Val1Label).(Val2Label).AvgMag(2:3,i)=(CIMag(:,i));
                O.(FiltLabel).(Val1Label).(Val2Label).AvgPhase(1,i)=mean(angle(x));
                O.(FiltLabel).(Val1Label).(Val2Label).AvgPhase(2:3,i)=CIAngle(:,i);

            end

            O.(FiltLabel).(Val1Label).(Val2Label).AvgMagOneSided= O.(FiltLabel).(Val1Label).(Val2Label).AvgMag(:,1:NumUniquePts);
            O.(FiltLabel).(Val1Label).(Val2Label).AvgPhaseOneSided= O.(FiltLabel).(Val1Label).(Val2Label).AvgPhase(:,1:NumUniquePts);
            O.(FiltLabel).(Val1Label).(Val2Label).AvgFreqRespOneSided= O.(FiltLabel).(Val1Label).(Val2Label).AvgFreqResp(:,1:NumUniquePts);

            O.(FiltLabel).(Val1Label).(Val2Label).f_full=f_full;
            O.(FiltLabel).(Val1Label).(Val2Label).NumUniquePts=NumUniquePts;
            O.(FiltLabel).(Val1Label).(Val2Label).f_onesided=f_onesided;
        
        end
    end
end

%% plotting

    
RepeatLabel=sprintf('Rep_%d',10);

for iFilt=1:1
    FiltLabel=sprintf('Filt_%d',iFilt);

    for iVal1=1:1

        Val1Label=sprintf('Val_%d',iVal1);
        
        for iVal2=1:1

            Val2Label=sprintf('Val_%d',iVal2);

            MagL=O.(FiltLabel).(Val1Label).(Val2Label).AvgMagOneSided(:,2:end);
%             AvgMag=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).AvgMag(1,:);

            PhaseL=O.(FiltLabel).(Val1Label).(Val2Label).AvgPhaseOneSided(:,2:end);
%             AvgPhase=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).AvgPhase(1,:);

            f_onesided=O.(FiltLabel).(Val1Label).(Val2Label).f_onesided(2:end);
            MagWelchAllRepsOneSided=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).MagWelchAllRepsOneSided(2:end);

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

FreqRange=9; % Hz
FreqIndice=2:NumUniquePts*FreqRange/NyqFreq;
ConfIntLabel={ 'FitAvg','FitLb','FitUb'};
lgd={'Avg Fitted', 'Avg ','Lb Fitted','Lb', 'Ub Fitted', 'Ub '};
ColorLabel={'b','r','g'};
FiltTitle={'Comb Filter', 'GS Filter' 'Blanking Method'};
Ts=1/fs;
% FilterNum=1;
% Val2Num=1;
%run for the nominal values 
% RepeatLabel=sprintf('Rep_%d',10);
Val1Num=10;
Val2Num=1;

FilterNum=1;
for iFilt=1:FilterNum
    
    FiltLabel=sprintf('Filt_%d',iFilt);
    
    for iVal1=1:Val1Num
        iVal1;
        Val1Label=sprintf('Val_%d',iVal1);
        for iVal2=1:Val2Num
        
            Val2Label=sprintf('Val_%d',iVal2);

            FreqVec=O.(FiltLabel).(Val1Label).(Val2Label).f_onesided(FreqIndice);

            for iRep=1:1 % fitting for all the repeats 

                RepeatLabel=sprintf('Rep_%d',iRep);

                Phase=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).PhaseAllRepsOneSided(FreqIndice);
                Mag=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).MagAllRepsOneSided(FreqIndice);
                FreqResp=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).FreqRespOneSided(FreqIndice);
                
                Estnp=3;
                Estnz=2;

                minfThres=1000;
                StabThres=1;
                IterLim=20;

                lb= [0 -2  -1 0 -10 4];
                ub= [1  0   0 3  -8 5];
                x0=0.01*ones(1,(Estnp+Estnz+1));
                [FittedFreqResp SysTf]= sysiden_basic2(Mag,x0,lb,ub,FreqVec,Estnp,Estnz,Ts);
                SysFreq = frd(FreqResp,FreqVec,'FrequencyUnit','Hz','Ts',Ts);
% 
%                 Estnp=1;
%                 Estnz=0;
% 
%                 minfThres=2;
%                 StabThres=1;
%                 IterLim=20;
% 
%                 lb= [-10    10];
%                 ub= [-10    10];
%                 x0=0.01*ones(1,(Estnp+Estnz+1));
%                 [FittedFreqResp SysTf]= sysiden_basic(Mag,x0,lb,ub,FreqVec,Estnp,Estnz,Ts);
%                 SysFreq = frd(FreqResp,FreqVec,'FrequencyUnit','Hz','Ts',Ts);

                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).nz = Estnz;
                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).lb = lb;
                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).ub = ub;
                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).SysFreq=SysFreq;
                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).FreqVec=FreqVec;

                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).SysTf=SysTf;
                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).FittedFreqResp=FittedFreqResp';
                
                MagWelch=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).MagWelchAllRepsOneSided(FreqIndice);
                
                [WelchFittedFreqResp WelchSysTf]= sysiden_basic2(MagWelch,x0,lb,ub,FreqVec,Estnp,Estnz,Ts);
                WelchSysFreq = frd(FreqResp,FreqVec,'FrequencyUnit','Hz','Ts',Ts);
                
                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).WelchSysFreq=WelchSysFreq;
                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).WelchSysTf=WelchSysTf;
                O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).WelchFittedFreqResp=WelchFittedFreqResp';
% 
            end


            for iConf=1:1  % fitting for Conf intervals.

    %             RepeatLabel=sprintf('Rep_%d',iRep);
                  ConfLabel=sprintf('ConfIntFit_%d',iConf);

    %             Mag=O.(ParLabel).(FiltLabel).(ValLabel).MagOneSided(i,FreqIndice);
    %             Phase=O.(ParLabel).(FiltLabel).(ValLabel).PhaseOneSided(i,FreqIndice);
    %             FreqResp=O.(ParLabel).(FiltLabel).(ValLabel).AvgFreqRespOneSided(i,FreqIndice);

                Phase=O.(FiltLabel).(Val1Label).(Val2Label).AvgPhaseOneSided(iConf,FreqIndice);
                Mag=O.(FiltLabel).(Val1Label).(Val2Label).AvgMagOneSided(iConf,FreqIndice);
                FreqResp=O.(FiltLabel).(Val1Label).(Val2Label).AvgFreqRespOneSided(iConf,FreqIndice);

                Estnp=3;
                Estnz=2;

                minfThres=1000;
                StabThres=1;
                IterLim=20;

                lb= [0.1 0.1 0.1 -1 -3 1 ];
                ub= [0.2  .2  .2  0 -2 2 ];
                x0=0.01*ones(1,(Estnp+Estnz+1));
    %             x0=Detx0(Estnp+Estnz+1,lb,ub);
                [FittedFreqResp SysTf]= sysiden_basic(Mag,x0,lb,ub,FreqVec,Estnp,Estnz,Ts);
                SysFreq = frd(FreqResp,FreqVec,'FrequencyUnit','Hz','Ts',Ts);
    %             SysTf=tfest(SysFreq,1,0,'Ts',Ts);

                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).np = Estnp;
                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).Mag = Mag;
                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).Phase = Phase;
                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).FreqResp = FreqResp;

                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).nz = Estnz;
                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).lb = lb;
                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).ub = ub;
                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).SysFreq=SysFreq;
                O.(FiltLabel).(Val1Label).(Val2Label).FreqVec=FreqVec;

                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).SysTf=SysTf;
                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).FittedFreqResp=FittedFreqResp';
            end
        end
    end
end


S=O;


%%
save('sysiden2023_target_ana4')

%% plotting
close all
RepeatLabel=sprintf('Rep_%d',1);
FiltLabels=["Comb" "GS", "Blanking"];
for iFilt=1:1
   
    FiltLabel=sprintf('Filt_%d',iFilt);
    FiltLabel2=FiltLabels(iFilt);

    for iVal1=1:10
        
        Val1Label=sprintf('Val_%d',iVal1);
        rndfixed=randi([90 110],1)/100;

        for iVal2=Val2Num
            
%         
            Val2Label=sprintf('Val_%d',iVal2);
            
            MagL=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).MagAllRepsOneSided;
            FittedFreqResp=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).FittedFreqResp;
            WelchFittedFreqResp=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).WelchFittedFreqResp;
            MagWelch=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).MagWelchAllRepsOneSided(FreqIndice);

            FreqVec=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).FreqVec;

            PhaseL=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).PhaseAllRepsOneSided;

            plotrange=1:NumUniquePts;
            figure(iFilt)
            subplot(2,1,1)
            semilogx(f_onesided,20*log10((MagL(2:end))))
            hold on
            semilogx(FreqVec,20*log10(FittedFreqResp))
%             hold 
%             semilogx(f_onesided,zeros(1,length(f_onesided)),'k')
            xlabel(' Freq. (Hz) ')
            ylabel('Mag. (dB)')
            ylim([-60 40 ])
            xlim([0 NyqFreq])
            title(FiltTitle{iFilt})
            legend(ConfIntLabel)
            subplot(2,1,2)
            semilogx(FreqVec,20*log10((MagWelch)))
            hold on
            semilogx(FreqVec,20*log10(WelchFittedFreqResp))
%             hold 
%             semilogx(f_onesided,zeros(1,length(f_onesided)),'k')
            xlabel(' Freq. (Hz) ')
            ylabel('Mag. (dB)')
            ylim([-60 40 ])
            xlim([0 NyqFreq])

            for iConf=1:1

%                 ConfLabel=sprintf('ConfIntFit_%d',iConf);
% 
%                 FittedMag=O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).FittedFreqResp*rndfixed;
%                 Mag=O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).Mag'*rndfixed;
%                 SysTf=O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).SysTf*rndfixed;
%                 
%                 WelchFittedFreqResp=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).WelchFittedFreqResp;
%                 MagWelchAllRepsOneSided=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).MagWelchAllRepsOneSided(FreqIndice);
% 
%                 figure(1)
%                 ColorLabel=lines(10);
%                 subplot(2,1,1)
%                 semilogx(FreqVec,(abs(FittedMag)),"LineWidth",2,"DisplayName",Val1Label);
%                 hold on
%                 semilogx(FreqVec,(abs(Mag)),'--',"DisplayName",Val1Label);
% 
% %                 title(FiltTitle(1))
%                 xlabel(' Freq. (Hz) ')
%                 ylabel('Mag.')
%                 
%                 subplot(2,1,2)
%                 semilogx(FreqVec,(abs(WelchFittedFreqResp)),"LineWidth",2,"DisplayName",Val1Label);
%                 hold on
%                 semilogx(FreqVec,(abs(MagWelchAllRepsOneSided)),'--',"DisplayName",Val1Label);
% 
% %                 title(FiltTitle(1))
%                 xlabel(' Freq. (Hz) ')
%                 ylabel('Mag.')
% %                 ylim([-40 40 ])

%                 figure(100+iVal2)
%                 [re,im,wout,sdre,sdim] = nyquist(SysTf);
% 
%                 re=squeeze(re);
%                 im=squeeze(im);
% 
%                 subplot(3,1,iFilt)
%                 plot(re,im)
%                 hold on 
            end
        end
    end
end
legend

%% Nyquist Analysis

global Wdiff FreqVec np nz 

i=1;
c_range=1:1:88;
WMargin=[1 1 2.2];
WCoef= [10 50 2];
nps=[3 3 3];
nzs=[1 1 1];
TitleLabels={'Mean Freq Response', 'Lower Bound Response', 'Upper Bound Resoponse'};
NyquistLabel={'Nyquist with Comb Filter', 'Nyquist with GS Filter', 'Nyquist with Blanking Method'};
FiltTitle={'Comb Filter', 'GS Filter' 'Blanking Method'};
ylimValues=[ 0.5 0.2 0.2 ];
LineLabel={'--','-','--o'};

for iFilt=1:FilterNum

    FiltLabel=sprintf('Filt_%d',iFilt);
    TitleLabel=sprintf('%s',FiltTitle{iFilt});

    for iVal1=1:Val1Num
        Val1Label=sprintf('Val_%d',iVal1);

        for iVal2=1:Val2Num
            Val2Label=sprintf('Val_%d',iVal2);
            
%             for iRep=1:RepNum
%                 RepeatLabel=sprintf('Rep_%d',iRep);
% 
% 
% %                 L0= O.(ParLabel).(FiltLabel).(ValLabel).(ConfIntLabel{i}).FittedFreqResp;
% %                 L=O.(ParLabel).(FiltLabel).(ValLabel).MagAllRepsOneSided(1,FreqIndice);
% %                 L0Tf=O.(ParLabel).(FiltLabel).(ValLabel).(ConfIntLabel{i}).SysTf;
% 
%                 L=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).MagAllRepsOneSided(FreqIndice);
%                 L0Tf=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).SysTf;
%                 L0=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).FittedFreqResp;
%                 Wdiff=abs(abs(L)-abs(L0))';
% 
%         %     Find W using this block
% 
%             %     x0=randn(1,nps(iFilt)+nzs(iFilt));
%                 x0=10*ones(1,nps(iFilt)+nzs(iFilt));
%                 np=nps(iFilt);
%                 nz=nzs(iFilt);
%                 [Wtest, Wtf]=fminsysiden(Wdiff,x0,FreqVec,nps(iFilt),nzs(iFilt),Ts);
%                    
%                 [re,im,wout,~,~]=nyquist(L0Tf,FreqVec*2*pi);
%                 re=squeeze(re);
%                 im=squeeze(im);
%                 nyq=re+im*1i;
% 
%                 ConfIntDist=abs(nyq+1);
%                 [MinConfIntMar,I(iVal1)]=min(ConfIntDist);
%                 W=Wtest;
%                 r = abs(W);
%                 WDist=abs(nyq'+1-r);
%                 [WMar,I]=min(WDist);  
%                 I
%                 
%                 O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).Wdiff=Wdiff;
%                 O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).Wfit=Wtest;
%                 O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).Wtf=Wtf;
%                 O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).r=r;
%                 O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).WDist=WDist;
%                 O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).WMar=WMar;
%                 
%                 O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).ConfIntDist=ConfIntDist;
%                 O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).MinConfIntMar=MinConfIntMar;
%             end
           
            for iConf=1:1
                
                ConfLabel=sprintf('ConfIntFit_%d',iConf);

%                 L0= O.(ParLabel).(FiltLabel).(ValLabel).(ConfIntLabel{i}).FittedFreqResp;
%                 L=O.(ParLabel).(FiltLabel).(ValLabel).MagAllRepsOneSided(1,FreqIndice);
%                 L0Tf=O.(ParLabel).(FiltLabel).(ValLabel).(ConfIntLabel{i}).SysTf;

                L=O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).Mag;
                L0Tf=O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).SysTf;
                L0=O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).FittedFreqResp;
                Wdiff=abs(abs(L)-abs(L0'));

        %     Find W using this block

            %     x0=randn(1,nps(iFilt)+nzs(iFilt));
                x0=10*ones(1,nps(iFilt)+nzs(iFilt));
                np=nps(iFilt);
                nz=nzs(iFilt);
                [Wtest, Wtf]=fminsysiden(Wdiff,x0,FreqVec,nps(iFilt),nzs(iFilt),Ts);
                   
                [re,im,wout,~,~]=nyquist(L0Tf,FreqVec*2*pi);
                re=squeeze(re);
                im=squeeze(im);
                nyq=re+im*1i;

                ConfIntDist=abs(nyq+1);
                [MinConfIntMar,I(iVal1)]=min(ConfIntDist);
                W=Wtest; 
                r = abs(W);
                WDist=abs(nyq'+1-r);
                [WMar,I]=min(WDist);  
                I
                
                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).Wdiff=Wdiff;
                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).Wfit=Wtest;
                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).Wtf=Wtf;
                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).r=r;
                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).WDist=WDist;
                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).WMar=WMar;

                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).ConfIntDist=ConfIntDist;
                O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).MinConfIntMar=MinConfIntMar;

%                 xlabel('Freq. (Hz')
%                 ylabel('Mag')
%                 legend({'L','L0'});
%                 subplot(6,iFilt,2)
%                 semilogx(FreqVec,abs(Wtest))
%                 hold on
%                 semilogx(FreqVec, abs(Wdiff))
% 
%                 legend({'Wfit','L-L0'});
%                 title(FiltTitle(iFilt))
%                 xlabel('Freq. (Hz')
%                 ylabel('Mag')


% %           
    %         pause
    % 
    %         calculate the margins using this block
    %         Wmin=(L-L0);
    %         MagWmin=abs(Wmin);
    %         PhaseWmin=angle(Wmin);
    %         Wmult=1:2/length(FreqVec):3-2/length(FreqVec);
    %         SysW = frd(WCoef(iFilt)*Wmult.*WMargin(iFilt).*MagWmin.*exp(1j*PhaseWmin/10),FreqVec,'FrequencyUnit','Hz','Ts',Ts);

    %     np=5;
    %     nz=5;
    %     iodelay = 0.1/Ts;
    % %    opt = tfestOptions('EnforceStability', true,'SearchMethod','gna');    
    %     Wtf = tfest(SysW,np,nz,iodelay,'Ts',Ts);
    %     
    %     WFreq=freqresp(Wtf,2*pi*FreqVec);
    %     WFreq=squeeze(WFreq);
    %     W=WFreq;

    %     figure(iFilt)
    %     semilogx(FreqVec,Wmult.*MagWmin)
    %     hold on
    %     semilogx(FreqVec, abs(WFreq))
    %     semilogx(FreqVec, MagWmin)
    %     hold off

    
    
    
%         figure(1)
%         subplot(2,1,1)
%         semilogx(FreqVec,(abs(L)),'LineWidth',1,'Color','k');
%         hold on
%         xlabel('Freq Hz')
%         ylabel('Mag (dB)')
%         title('Frequency Response with Different Filters')
%         semilogx(FreqVec,(abs(L0)),LineLabel{iFilt},'LineWidth',2,'Color','k');
%         lgd={'L (GS)', 'L_0 (GS)','L (Comb)', 'L_0 (Comb)','L (Blanking)', 'L_0 (Blanking)'};
%         legend(lgd);
% 
%         subplot(2,1,2)
%         semilogx(FreqVec,(abs(W)),LineLabel{iFilt},'Color','k','LineWidth',2);
%         hold on
%         semilogx(FreqVec,(abs(Wdiff)),LineLabel{iFilt},'Color','k','LineWidth',2);
% 
%         ylabel('|W = L-L_0|  (dB)')
%     %     ylim([-40 20])
%         xlim([0 9])
%         xlabel('Frequency (Hz)')
%         lgd={'|W| (GS)','|W| (Comb)','|W| (Blanking)'};
%         legend(lgd);
%         title('|W| for each Frequency Response ')    
% 
%         figure(2)
%         subplot(3,1,iFilt)
%         [re,im,wout,~,~]=nyquist(L0Tf,FreqVec*2*pi);
%         re=squeeze(re);
%         im=squeeze(im);
%         plot(re,-im,'--','Color','k','Linewidth',2)
%         hold on
%         viscircles(cc,rr,'LineStyle','-','Color','k','Linewidth',.1)
%         plot(100,100,'LineStyle','-','Color','k','Linewidth',.1)
%     %     plot(cc(:,1),cc(:,2),'*')
%         viscircles([re(I) -im(I)],r(I),'LineStyle','-','Color','k','Linewidth',2)
%         plot(100,100,'LineStyle','-','Color','k','Linewidth',2)
%         plot(-1,0,'+','LineWidth',2,'Color', 'k')
%         title(NyquistLabel{iFilt})
%         xlim([-1.01 2]) 
%         ylim([-ylimValues(iFilt) ylimValues(iFilt) ]);
%         xlabel('Real Axis')
%         ylabel('Im Axis')
%         grid
%         subplot(3,1,1)
%         legend('Nyquist Plot','Uncertainty Circles','Critical Circle','Critical Point','Location','West')

% %         [re,im,wout,~,~]=nyquist(L0Tf,FreqVec*2*pi);
% %         re=squeeze(re);
% %         im=squeeze(im);

    %     subplot(3,1,iFilt)

    %    plot(re,im)
    %    hold
    %    plot([re(I) -1], [im(I) 0],'r')

    
            end
        end
    end
end

% S.(ParLabel)=O.(ParLabel);
%%
save('SysIden_nyq_ver3')

%% plotting
iVal1=9;
iVal2=1;

iFilt=1;
Val1Label=sprintf('Val_%d',iVal1);
Val2Label=sprintf('Val_%d',iVal2);
FiltLabel=sprintf('Filt_%d',iFilt);

for iConf=1:1
    
    ConfLabel=sprintf('ConfIntFit_%d',iConf);

    figure(iVal1)
    subplot(2,3,iFilt)
    
    L= O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).Mag;
    L0= O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).FittedFreqResp;

    semilogx(FreqVec,abs(L))
    hold on
    semilogx(FreqVec, abs(L0))
    title(FiltTitle(iFilt))

    subplot(2,3,iFilt+3)
    semilogx(FreqVec,abs(Wtest))
    hold on
    semilogx(FreqVec, abs(Wdiff))

    legend({'Wfit','L-L0'});
    title(FiltTitle(iFilt))
    xlabel('Freq. (Hz')
    ylabel('Mag')


end

%% PLotting target_vs_paresis
FilterNum=1;
RepeatLabel=sprintf('Rep_%d',10);
Val1Range=[1 10];
Val2Range=[1 10];
clear SysDen
for iFilt=1:FilterNum
    FiltLabel=sprintf('Filt_%d',iFilt);

    for iVal1=Val1Range(1):Val1Range(2)
        Val1Label=sprintf('Val_%d',iVal1);

        for iVal2=Val2Range(1):Val2Range(2)
            Val2Label=sprintf('Val_%d',iVal2);
            
%             RepMar(iVal2,iVal1)=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).MinConfIntMar;
%             WMar(iVal2,iVal1)=O.(FiltLabel).(Val1Label).(Val2Label).(RepeatLabel).WMar;

            for iConf=1:1
                
                ConfLabel=sprintf('ConfIntFit_%d',iConf);
                
                ConfMar(iVal1,iVal2)=O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).MinConfIntMar;
                WMar(iVal1,iVal2)=O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).WMar;
                SysNum(iVal1,iVal2)=O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).SysTf.Numerator;
                SysDen(iVal1,iVal2)=O.(FiltLabel).(Val1Label).(Val2Label).(ConfLabel).SysTf.Denominator;
                    
            end
        end
    end
    
    for iVal=1:length(ConfMar(:,1))
        [den1(iVal,:)]=SysDen{iVal};
    
    end
    figure(1000) 
    subplot(3,1,1)
    semilogx(1:length(ConfMar(:,1)),ConfMar(1,:))
    hold
    subplot(3,1,2)
%     semilogx(1:length(ConfMar(:,1)),WMar(1,:))
    semilogx(1:length(ConfMar(:,1)),den1(:,2))

    hold
    subplot(3,1,3)

    semilogx(1:length(ConfMar(:,1)),den1(:,3))
    hold
    
%     ParLabelTitle=VarLabels{iPar};
%     xlabel(ParLabelTitle)
%     ylabel(FiltLabel)
%     
%     subplot(2,1,1)
%     semilogx(1:iVal1Num,RepMar)
%     

qy=flip(O.('Paresis').ParamValues);
qx=flip(O.('Target').ParamValues);


    figure
    mesh(qx,qy,ConfMar)
    xlabel('Target')
    ylabel('Paresis Level')
    zlabel('Stability Margin')
    
        figure
    mesh(qx,qy,WMar)
    xlabel('Target')
    ylabel('Paresis Level')
    zlabel('Stability Margin')

end
              
%% 2d black figure




%%


function [Wtest Wtf]=fminsysiden(Wdiff,x0,FreqVec,np,nz,Ts)
% system indentification for with minimized upper bound for the freq
% respons Wdiff
% Wdiff=W;
% nz=2;
% np=3;
A=[];
b=[];
Aeq=[];
beq=[];

lb=-100*ones(np+nz,1);
ub=100*ones(np+nz,1);

[x minf outt s] = fmincon(@objective_f,x0,A,b,Aeq,beq,lb,ub,@const_f);

x 
minf

Wtf=ConstructTF(x,nz);

Wtest=freqresp(Wtf,FreqVec,'Hz');
Wtest=squeeze(Wtest)';
 
% figure(100)
% semilogx(FreqVec,abs(Wtest))
% hold
% semilogx(FreqVec,abs(Wdiff))

% 
% NumUniquePts = ceil(NumFFT/2);
% f_onesided=f_vec(1:NumUniquePts);
% L0_FreqResp=freqresp(L0,f_onesided,'Hz');
% L0_FreqResp=squeeze(L0_FreqResp);
% 
% L_onesided=L(1:NumUniquePts);
% % Ldiff=L0-L;
% 
% figure
% subplot(2,1,1)
% plot(t,y_vec)
% hold on
% plot(t,x_vec)
% plot(t,y);
% hold off
% subplot(2,1,2)
% semilogx(f_onesided,20*log10(abs(L0_FreqResp)))
% hold 
% semilogx(f_onesided,20*log10(abs(L_onesided)));
 
end
function [TFun,num,den]=ConstructTF(x,z)
global Ts
num=1;
den=1;

p=length(x)-z-1;
for iz=1:z
    
    num(iz+1)=x(iz+1);
 
end

for ip=1:p    
    
    den(ip+1)=x(1+z+ip);

end

num=num*x(1);
TFun=tf(num,den,Ts);

end

function obj=objective_f(x)
global FreqVec Wdiff nz Ts

% n=length(Wdiff);
Stf=ConstructTF(x,nz);

S=freqresp(Stf,FreqVec,'Hz');
S=squeeze(S);

obj=sum((abs(Wdiff)-abs(S')).^2);

end

function [c,ceq] = const_f(x)

    global Wdiff FreqVec nz 

    % Construct TF

    xtf=ConstructTF(x,nz);

    FreqTf=freqresp(xtf,FreqVec,'Hz');
    FreqTf=squeeze(FreqTf)';

    MinVal=max(abs(Wdiff)-abs(FreqTf));  % mignt change to max over sum

    for i=1:length(x)
        
        c(i)=MinVal;

    end

    % Nonlinear equality constraints
    ceq = [];

end

function [Wtest Wtf]=sysiden_basic(FreqResp,x0,lb,ub,FreqVector,np,nz,Ts)
% system indentification using fmincon\
global minfThres IterLim FreqVec
minf=100000;
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
[Wtf num den]=ConstructTF(x,nz);
% [Wtf num den]=SpecialTF(x,nz);
zrs=roots(num);
pls=roots(den);
num 
den
Wtest=freqresp(Wtf,FreqVec,'Hz');
Wtest=squeeze(Wtest)';

end
function [Wtest Wtf]=sysiden_basic2(FreqResp,x0,lb,ub,FreqVector,np,nz,Ts)
% system indentification using fmincon\
global minfThres IterLim FreqVec
minf=100000;
nIter=0;
A=[];
b=[];
Aeq=[];
beq=[];


while minf>minfThres && nIter<IterLim
    
%     options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
    
    [x minf outt s] = fmincon(@objective_fit3,x0,A,b,Aeq,beq,lb,ub,@stability_z);
    nIter=nIter+1;

    MinFs(nIter,:)=minf;
    xs(nIter,:)=x;
    x0=0.01*ones(1,length(x));
    
end

[minFval I]=min(MinFs);
minFval
x=xs(I,:)
nIter
[Wtf num den]=ConstructTF(x,nz);
% [Wtf num den]=SpecialTF(x,nz);
zrs=roots(num);
pls=roots(den);
num 
den
Wtest=freqresp(Wtf,FreqVec,'Hz');
Wtest=squeeze(Wtest)';

end

function obj=objective_fit2(x)

global FreqVec FreqResp Estnz

Stf=ConstructTF(x,Estnz);
% Stf=SpecialTF(x,Estnz);


S=freqresp(Stf,FreqVec,'Hz');
S=squeeze(S);

obj1=sum((abs(FreqResp)-abs(S')).^2);
% obj2=sum(((angle(FreqResp)-angle(S'))).^2);
obj=obj1;
end
function obj=objective_fit3(x)

global FreqVec FreqResp Estnz

Stf=ConstructTF(x,Estnz);
% Stf=SpecialTF(x,Estnz);


S=freqresp(Stf,FreqVec,'Hz');
S=squeeze(S);

obj1=sum((abs(FreqResp)-abs(S)).^2);
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