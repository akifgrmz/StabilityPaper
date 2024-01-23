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
mdlWks = get_param('SysIden_NSFsim','ModelWorkspace');
assignin(mdlWks,'rngseed',1000);

    %paretic hand
low=1;
high=4;
ParHandValues=logspace(log10(low),log10(high),NumPoints);
AdressString{1}='SysIden_NSFsim/GainP';
BlockString{1}='Gain';
NominalString{1}='2';
VarLabels{1}='Phand';
O.(VarLabels{1}).ParValues=ParHandValues;
O.(VarLabels{1}).AdressString=AdressString{1};
O.(VarLabels{1}).NomValue=2;

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
RepeatNum=10;
% O=struct;
rngseed=1000:1000+RepeatNum;
a=0.1;


for i=1:length(ParHandValues)
    ParHandValueString{i}=sprintf('%0.8f',ParHandValues(i));
end
    %stim hand
low=0.3/2;
high=0.3*2;
StimHandValues=logspace(log10(low),log10(high),NumPoints);
AdressString{2}='SysIden_NSFsim/Paretic Hand/Stim Hand Angle/Transfer Fcn';
BlockString{2}='Denominator';
NominalString{2}='[0.3 1]';
VarLabels{2}='NPhand';
O.(VarLabels{2}).ParValues=StimHandValues;
O.(VarLabels{2}).AdressString=AdressString{2};
O.(VarLabels{2}).NomValue=0.3;


for i=1:length(StimHandValues)
%     str1=sprintf('[%0.4f 1]',NonParHandValues(i));
    StimHandValueString{i}=sprintf('[%0.8f 1]',StimHandValues(i));
end

    %occlusion K
low=20/4;
high=20*4;
OccKValues=logspace(log10(low),log10(high),NumPoints);

AdressString{3}='SysIden_NSFsim/Paretic Hand/Occlusion/Noise';
BlockString{3}='Value';
NominalString{3}='20';
VarLabels{3}='OccK';
O.(VarLabels{3}).ParValues=OccKValues;
O.(VarLabels{3}).AdressString=AdressString{3};
O.(VarLabels{3}).NomValue=20;

for i=1:length(OccKValues)
    KValuesString{i}=sprintf('%0.8f',OccKValues(i));
end

% m-Wave
low=0.01;
high=0.2;
mWaveVar=logspace(log10(low),log10(high),NumPoints);
AdressString{4}='SysIden_NSFsim/Paretic Hand/mWave/GainVar';
BlockString{4}='Gain';
NominalString{4}=num2str(low);
VarLabels{4}='mWaveVar';
O.(VarLabels{4}).ParValues=mWaveVar;
O.(VarLabels{4}).AdressString=AdressString{4};
O.(VarLabels{4}).NomValue=low;

for i=1:length(mWaveVar)
    mWaveVarString{i}=sprintf('%0.8f',mWaveVar(i));
end

% stim freq
low=10;
high=40;
StimFreqVar=logspace(log10(low),log10(high),NumPoints);
AdressString{5}='SysIden_NSFsim/Paretic Hand/mWave/GainVar';
BlockString{5}='Gain';
NominalString{5}='20';
VarLabels{5}='StimFreq';
O.(VarLabels{5}).ParValues=StimFreqVar;
O.(VarLabels{5}).AdressString=AdressString{5};
O.(VarLabels{5}).NomValue=20;

for i=1:length(StimFreqVar)
    StimFreqVarString{i}=sprintf('%0.8f',StimFreqVar(i));
end

ValueStringMat=[ParHandValueString;StimHandValueString;KValuesString;mWaveVarString;StimFreqVarString];

for iPar=5:5
    
    ParLabel=VarLabels{iPar};
    ParValues=ValueStringMat(iPar,:);

for iFilter=1:NumFilter
    
    FiltLabel=sprintf('Filt%d',iFilter);
    str1=sprintf('%d',iFilter);
    set_param('SysIden_NSFsim/Effort Estimator/FilterNum','Value',str1)
    
    for iRepeat=1:10
            
        RepeatLabel=sprintf('Rep_%d',iRepeat);
        assignin(mdlWks,'rngseed',rngseed(iRepeat));
        O.(ParLabel).(FiltLabel).(RepeatLabel).RNGSeed=rngseed(iRepeat);
        
        for iValue=1:NumPoints
            
%             sprintf('%d% completed',100*(iValue+(iRepeat-1)*iValue+iValue*iRepeat*(iFilter-1))/(RepeatNum*NumPoints*NumFilter))
            
            if iPar==5
                assignin(mdlWks,'stimfreq',StimFreqVar(iValue));
            else
                set_param(AdressString{iPar},BlockString{iPar},ParValues{iValue});
            end
            
            sim('SysIden_NSFsim',SimTime);

            O.(ParLabel).(FiltLabel).(RepeatLabel).x_sampled(iValue,:)=ans.Outcome(:,1);
            O.(ParLabel).(FiltLabel).(RepeatLabel).y_sampled(iValue,:)=ans.Outcome(:,5);
        
        end
    end
end

set_param(AdressString{iPar},BlockString{iPar},NominalString{iPar}) 
assignin(mdlWks,'rngseed',1000);
filename=sprintf('worksim2_ipar%d',iPar);
save(filename);
 
end


%% Freq response
% clear all
ParNum=5;
ParFileName=sprintf('worksim2_ipar%d',ParNum);
load(ParFileName);
ConfIntLabel={ 'Mean','Lb','Ub'};
FiltTitle={'Comb Filter', 'GS Filter' 'Blanking Method'};
NominalInd=50;  %assumption
ParVec=[ParNum ParNum];
% Reordering

for iPar=ParVec(1):ParVec(2)
    
        ParLabel=VarLabels{iPar};
        
    for iFilt=1:3
        
        FiltLabel=sprintf('Filt%d',iFilt);
        FiltLabel2=sprintf('Filt_%d',iFilt);
%         0.(FiltLabel2).NumUniquePts=NumUniquePts;

        for iRep=1:10
            
            RepeatLabel=sprintf('Rep_%d',iRep);
            
            for iValue=1:NumPoints
            
                ValLabel=sprintf('Val_%d',iValue);

                O.(ParLabel).(FiltLabel2).(ValLabel).(RepeatLabel).x_sampled=O.(ParLabel).(FiltLabel).(RepeatLabel).x_sampled(iValue,:);
                O.(ParLabel).(FiltLabel2).(ValLabel).(RepeatLabel).y_sampled=O.(ParLabel).(FiltLabel).(RepeatLabel).y_sampled(iValue,:);
    %             clear O.(ParLabel).(FiltLabel).(RepeatLabel).x_sampled(iValue,:)
    %             clear O.(ParLabel).(FiltLabel).(RepeatLabel).y_sampled(iValue,:)

            end
        end
    end
end

% downsampling 
DownSampleRate=1;

for iPar=ParVec(1):ParVec(2)
    
        ParLabel=VarLabels{iPar};
        ParLabel2=strcat(ParLabel,sprintf('_%d',fs/DownSampleRate));
        
        O.(ParLabel2).ParValues=O.(ParLabel).ParValues;
        O.(ParLabel2).NumValue=O.(ParLabel).NomValue;
        O.(ParLabel2).AddressString=O.(ParLabel).AdressString;

    for iFilt=1:3
        
        FiltLabel=sprintf('Filt_%d',iFilt);
%         0.(FiltLabel2).NumUniquePts=NumUniquePts;

        for iRep=1:10
            
            RepeatLabel=sprintf('Rep_%d',iRep);
            
            for iValue=1:NumPoints
            
            ValLabel=sprintf('Val_%d',iValue);
            O.(ParLabel2).(FiltLabel).(ValLabel).(RepeatLabel).x_sampled=...
                downsample(O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).x_sampled,DownSampleRate);
            
            O.(ParLabel2).(FiltLabel).(ValLabel).(RepeatLabel).y_sampled=...
                downsample(O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).y_sampled,DownSampleRate);

            end
        end
    end
end

O.DownSampled.DownSampleRate=DownSampleRate;

Ts=(1/fs)*DownSampleRate;
tol = 1e-5;
IoDelay=0.1;
DelayNumSamp=0.1/Ts;
NumFFT=(SimTime/Ts)+1-DelayNumSamp;  % removing the delay 100 ms from the response
f_full=linspace(-1/2/Ts,1/2/Ts,NumFFT);
NumUniquePts = ceil(length(f_full)/2);
f_onesided=f_full(NumUniquePts:end);
TimeVec=0:Ts:SimTime;
TimeIndice=DelayNumSamp+1:length(TimeVec); %this is for ignoring the first 100ms due to the fixed delay


for iPar=ParVec(1):ParVec(2)
    
    ParLabel=VarLabels{iPar};
    ParLabel=strcat(ParLabel,sprintf('_%d',fs/DownSampleRate));

    for iFilter=1:3

        FiltLabel=sprintf('Filt_%d',iFilter);

        for iVal=1:100

            ValLabel=sprintf('Val_%d',iVal);

                %           O.(FiltLabel).FT(iRepeat,:)=zeros(1,NumFFT);
                %           x_v=O.(FiltLabel).x_shf(iRepeat,:);
                %           y_v=O.(FiltLabel).y_shf(iRepeat,:);
                
            for iRep=1:10

            RepeatLabel=sprintf('Rep_%d',iRep);

            x_v=O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).x_sampled;
            y_v=O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).y_sampled;
            
            [FreqRespWelch, f_welch]=tfestimate(x_v,y_v,200,[100],[f_onesided],fs/DownSampleRate);
%             [FreqRespWelch, f_welch]=tfestimate(x_v,y_v,200,[100],[f_onesided],fs/DownSampleRate);
            
            FT_L=fft(y_v,NumFFT)./fft(x_v,NumFFT);
            FT_L(abs(FT_L) < tol) = 0;
            O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).FT=FT_L;
            O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).FreqRespOneSided=FT_L(1:NumUniquePts);
            O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).MagAllRepsOneSided=abs(FT_L(1:NumUniquePts));
            O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).PhaseAllRepsOneSided=angle(FT_L(1:NumUniquePts));
            
            O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).FreqRespWelch=FreqRespWelch;
            O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).FreqRespWelchOneSided=FreqRespWelch(1:NumUniquePts);
            O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).MagWelchAllRepsOneSided=abs(FreqRespWelch(1:NumUniquePts));
            O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).PhaseWelchAllRepsOneSided=angle(FreqRespWelch(1:NumUniquePts));
            
            end

        a=0.05;
        for i=1:fs*SimTime/DownSampleRate

            for iRep=1:10

                RepeatLabel=sprintf('Rep_%d',iRep);

                x_v(iRep)=O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).x_sampled(i);
                y_v(iRep)=O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).y_sampled(i);
            end

            CI(:,i)=ConfInt(x_v,a);
            O.(ParLabel).(FiltLabel).(ValLabel).AvgTimeIn(1,i)=mean(x_v);
            O.(ParLabel).(FiltLabel).(ValLabel).AvgTimeIn(2:3,i)=CI(:,i);

            CI(:,i)=ConfInt(y_v,a);
            O.(ParLabel).(FiltLabel).(ValLabel).AvgTimeOut(1,i)=mean(y_v);
            O.(ParLabel).(FiltLabel).(ValLabel).AvgTimeOut(2:3,i)=CI(:,i);

         end

        for i=1:NumFFT

            for iRep=1:10

                RepeatLabel=sprintf('Rep_%d',iRep);

                x(iRep)=O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).FT(i);
                
            end

    %         x=O.(ParLabel).(FiltLabel).(ValLabel).FT(:,i);

            x(abs(x) < tol) = 0; % removing insignificant valuese for better accu in phase

            CI(:,i)=ConfInt(x,a);
            O.(ParLabel).(FiltLabel).(ValLabel).AvgFreqResp(1,i)=mean(x);
            O.(ParLabel).(FiltLabel).(ValLabel).AvgFreqResp(2:3,i)=CI(:,i);

            CIMag(:,i)=ConfInt(abs(x),a);
            CIAngle(:,i)=ConfInt(angle(x),a);

            O.(ParLabel).(FiltLabel).(ValLabel).AvgMag(1,i)=mean(abs(x));
            O.(ParLabel).(FiltLabel).(ValLabel).AvgMag(2:3,i)=(CIMag(:,i));
            O.(ParLabel).(FiltLabel).(ValLabel).AvgPhase(1,i)=mean(angle(x));
            O.(ParLabel).(FiltLabel).(ValLabel).AvgPhase(2:3,i)=CIAngle(:,i);

        end

        O.(ParLabel).(FiltLabel).(ValLabel).AvgMagOneSided= O.(ParLabel).(FiltLabel).(ValLabel).AvgMag(:,1:NumUniquePts);
        O.(ParLabel).(FiltLabel).(ValLabel).AvgPhaseOneSided= O.(ParLabel).(FiltLabel).(ValLabel).AvgPhase(:,1:NumUniquePts);
        O.(ParLabel).(FiltLabel).(ValLabel).AvgFreqRespOneSided=O.(ParLabel).(FiltLabel).(ValLabel).AvgFreqResp(:,1:NumUniquePts);

        O.(ParLabel).(FiltLabel).(ValLabel).f_full=f_full;
        O.(ParLabel).(FiltLabel).(ValLabel).NumUniquePts=NumUniquePts;
        O.(ParLabel).(FiltLabel).(ValLabel).f_onesided=f_onesided;

        end
    end
end

%% plotting

for iPar=ParVec(1):ParVec(2)
    
    ParLabel=VarLabels{iPar};
    ParLabel=strcat(ParLabel,sprintf('_%d',fs/DownSampleRate));

    for iFilt=1:3

        FiltLabel=sprintf('Filt_%d',iFilt);
        
        for iVal=20:20

            ValLabel=sprintf('Val_%d',iVal);

            MagL=O.(ParLabel).(FiltLabel).(ValLabel).AvgMagOneSided(:,2:end);
            PhaseL=O.(ParLabel).(FiltLabel).(ValLabel).AvgPhaseOneSided(:,2:end);
            f_onesided=O.(ParLabel).(FiltLabel).(ValLabel).f_onesided(2:end);
            MagWelchAllRepsOneSided=O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).MagWelchAllRepsOneSided(2:end);

%             PhaseAllRepsOneSided=O.(ParLabel).(FiltLabel).(ValLabel).PhaseAllRepsOneSided(1,:);
%             MagAllRepsOneSided=O.(ParLabel).(FiltLabel).(ValLabel).MagAllRepsOneSided(1,:);
            
            plotrange=1:NumUniquePts;
            figure(iVal)
            subplot(3,1,iFilt)
            semilogx(f_onesided,20*log10((MagWelchAllRepsOneSided)))
            hold 
            semilogx(f_onesided,20*log10((MagL)))

            xlabel(' Freq. (Hz) ')
            ylabel('Mag. (dB)')
            title(FiltTitle{iFilt})

            ylim([-50 10 ])
            xlim([0 10])

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
FreqIndice=2:NumUniquePts*FreqRange/NyqFreq*DownSampleRate;
ConfIntLabel={ 'FitAvg','FitLb','FitUb'};
lgd={'Avg Fitted', 'Avg ','Lb Fitted','Lb', 'Ub Fitted', 'Ub '};
ColorLabel={'b','r','g'};
FiltTitle={'Comb Filter', 'GS Filter' 'Blanking Method'};
Ts=1/fs*DownSampleRate;

%run for the nominal values 

for iPar=ParVec(1):ParVec(2)
    
    ParLabel=VarLabels{iPar};
    ParLabel=strcat(ParLabel,sprintf('_%d',fs/DownSampleRate));

for iFilt=1:3
    
    FiltLabel=sprintf('Filt_%d',iFilt);
    
    for iVal=1:100
        
        ValLabel=sprintf('Val_%d',iVal);
        FreqVec=O.(ParLabel).(FiltLabel).(ValLabel).f_onesided(FreqIndice);
        
        for iRep=1:10 % fitting for all the repeats 
            
            RepeatLabel=sprintf('Rep_%d',iRep);
 
            Phase=O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).PhaseAllRepsOneSided(FreqIndice);
            Mag=O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).MagAllRepsOneSided(FreqIndice);
            FreqResp=O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).FreqRespOneSided(FreqIndice);

            Estnp=1;
            Estnz=0;

            minfThres=2;
            StabThres=1;
            IterLim=20;

            lb= [0     .5];
            ub= [0.1    .6];c
            x0=0.01*ones(1,(Estnp+Estnz+1));
            [FittedFreqResp SysTf]= sysiden_basic(Mag,x0,lb,ub,FreqVec,Estnp,Estnz,Ts);
            SysFreq = frd(FreqResp,FreqVec,'FrequencyUnit','Hz','Ts',Ts);

            O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).nz = Estnz;
            O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).lb = lb;
            O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).ub = ub;
            O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).SysFreq=SysFreq;
            O.(ParLabel).(FiltLabel).(ValLabel).FreqVec=FreqVec;


            O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).SysTf=SysTf;
            O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).FittedFreqResp=FittedFreqResp';
            
        end
        
        
        for iConf=1:3  % fitting for Conf intervals.
            
%             RepeatLabel=sprintf('Rep_%d',iRep);
              ConfLabel=sprintf('ConfIntFit_%d',iConf);
    
%             Mag=O.(ParLabel).(FiltLabel).(ValLabel).MagOneSided(i,FreqIndice);
%             Phase=O.(ParLabel).(FiltLabel).(ValLabel).PhaseOneSided(i,FreqIndice);
%             FreqResp=O.(ParLabel).(FiltLabel).(ValLabel).AvgFreqRespOneSided(i,FreqIndice);

            Phase=O.(ParLabel).(FiltLabel).(ValLabel).AvgPhaseOneSided(iConf,FreqIndice);
            Mag=O.(ParLabel).(FiltLabel).(ValLabel).AvgMagOneSided(iConf,FreqIndice);
            FreqResp=O.(ParLabel).(FiltLabel).(ValLabel).AvgFreqRespOneSided(iConf,FreqIndice);

%             SysFreq = frd(Mag.*exp(1j*Phase/pi/1),FreqVec,'FrequencyUnit','Hz','Ts',Ts);
%             SysFreq = frd(FreqResp,FreqVec,'FrequencyUnit','Hz','Ts',Ts);
%             SysTime = frd(FreqResp,fs);     

%             init_sys = idtf([1 0],[1 1]);
%             
%             init_sys.Structure.Numerator.Minimum = eps;
%             init_sys.Structure.Denominator.Minimum = eps;
%             init_sys.Structure.Denominator.Free(end) = true;
%             init_sys.Ts=Ts;
%             iodelay = 0.1/Ts;
%             Estnp=nps(iFilt);
%             Estnz=nzs(iFilt);
%             opt = tfestOptions;
%             opt.SearchMethod='lm';
%             opt.Advanced.StabilityThreshold.z=0.99;
%             opt.EnforceStability=1;
%             SysTf = tfest(SysFreq,init_sys,opt);
%             SysTf = tfest(SysFreq,Estnp,Estnz,'Ts',Ts,opt);
% 
%             poles=roots(SysTf.Denominator)
%             zeros=roots(SysTf.Numerator)

%             nps=[3 1 1];
%             nzs=[1 1 1];
            Estnp=1;
            Estnz=0;

            minfThres=2;
            StabThres=1;
            IterLim=20;

%             lb=-10*ones((Estnp+Estnz+1),1);
%             ub=10*ones((Estnp+Estnz+1),1);
            lb= [0     .5];
            ub= [0.1    .6];
            x0=0.01*ones(1,(Estnp+Estnz+1));
%             x0=Detx0(Estnp+Estnz+1,lb,ub);
            [FittedFreqResp SysTf]= sysiden_basic(Mag,x0,lb,ub,FreqVec,Estnp,Estnz,Ts);
            SysFreq = frd(FreqResp,FreqVec,'FrequencyUnit','Hz','Ts',Ts);
%             SysTf=tfest(SysFreq,1,0,'Ts',Ts);

            O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).np = Estnp;
            O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).Mag = Mag;
            O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).Phase = Phase;
            O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).FreqResp = FreqResp;

            O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).nz = Estnz;
            O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).lb = lb;
            O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).ub = ub;
            O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).SysFreq=SysFreq;
            O.(ParLabel).(FiltLabel).(ValLabel).FreqVec=FreqVec;
                
%             FittedFreqResp=freqresp(SysTf,2*pi*FreqVec);
%             FittedFreqResp=squeeze(FittedFreqResp);
            O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).SysTf=SysTf;
            O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).FittedFreqResp=FittedFreqResp';

%             subplot(2,1,2)
%             semilogx(FreqVec,angle(FittedFreqResp),"LineWidth",2,"Color",ColorLabel{i});
%             hold on
%             semilogx(FreqVec,Phase,'--',"Color",ColorLabel{i});
%             xlabel('Freq. (Hz) ')
%             ylabel('Phase (rad)')

%             figure(iFilt+3)
%             bode(SysTf)
%             legend(ConfIntLabel)
%             title(FiltTitle(iFilt))
%             hold on
%             figure(iFilt+6)
%             nyquist(SysTf)
%             hold on

% %             figure(iVal+iFilt)
% %             subplot(2,1,1)
% %             semilogx(FreqVec,20*log10(abs(FittedFreqResp)),"LineWidth",2,"Color",ColorLabel{iConf});
% %             title(FiltTitle(iFilt))
% %             hold on
% %             semilogx(FreqVec,20*log10(Mag),'--',"Color",ColorLabel{iConf});
% %             legend(lgd)
% %             xlabel(' Freq. (Hz) ')
% %             ylabel('Mag. (dB)')
% %             ylim([-40 40 ])
% %             
% %             subplot(2,1,2)
% %             semilogx(FreqVec,20*log10(abs(FittedFreqResp')-abs(Mag)),"LineWidth",2,"Color",ColorLabel{iConf});
% %             hold on
% %             ylabel('|FittedMag|-|Mag| (dB)')


%             figure(iFilt+6)
%             rlocus(SysTf)
%             legend(Label)
%             title(FiltTitle(iFilt))
%             hold on
% 
%             [mag,phase,wout,sdmag,sdphase] = bode(SysTf,2*pi*FreqVec);
%             mag=squeeze(mag);
%             semilogx(wout/pi/2,mag)
    
        end
    end
end
S.(ParLabel)=O.(ParLabel);

end

save('work2022_par5')

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

            FittedMag=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).FittedFreqResp;
            Mag=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).Mag;
            
            SysTf=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).SysTf;

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


for iPar=ParVec(1):ParVec(2)
    
    ParLabel=VarLabels{iPar};
    ParLabel=strcat(ParLabel,sprintf('_%d',fs/DownSampleRate));

    for iFilt=1:3
        
        FiltLabel=sprintf('Filt_%d',iFilt);
        TitleLabel=sprintf('%s',FiltTitle{iFilt});

        for iVal=1:100

            ValLabel=sprintf('Val_%d',iVal);
            
            for iRep=1:10
                RepeatLabel=sprintf('Rep_%d',iRep);


%                 L0= O.(ParLabel).(FiltLabel).(ValLabel).(ConfIntLabel{i}).FittedFreqResp;
%                 L=O.(ParLabel).(FiltLabel).(ValLabel).MagAllRepsOneSided(1,FreqIndice);
%                 L0Tf=O.(ParLabel).(FiltLabel).(ValLabel).(ConfIntLabel{i}).SysTf;

                L=O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).MagAllRepsOneSided(FreqIndice);
                L0Tf=O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).SysTf;
                L0=O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).FittedFreqResp;
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
                [MinConfIntMar,I(iVal)]=min(ConfIntDist);
                W=Wtest;
                r = abs(W);
                WDist=abs(nyq'+1-r);
                [WMar,I]=min(WDist);  
                I
                
                O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).Wdiff=Wdiff;
                O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).Wfit=Wtest;
                O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).Wtf=Wtf;
                O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).r=r;
                O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).WDist=WDist;
                O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).WMar=WMar;
                
                O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).ConfIntDist=ConfIntDist;
                O.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).MinConfIntMar=MinConfIntMar;
            end
           
            for iConf=1:3
                
                ConfLabel=sprintf('ConfIntFit_%d',iConf);

%                 L0= O.(ParLabel).(FiltLabel).(ValLabel).(ConfIntLabel{i}).FittedFreqResp;
%                 L=O.(ParLabel).(FiltLabel).(ValLabel).MagAllRepsOneSided(1,FreqIndice);
%                 L0Tf=O.(ParLabel).(FiltLabel).(ValLabel).(ConfIntLabel{i}).SysTf;

                 L=O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).Mag;
                L0Tf=O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).SysTf;
                L0=O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).FittedFreqResp;
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
                [MinConfIntMar,I(iVal)]=min(ConfIntDist);
                W=Wtest;
                r = abs(W);
                WDist=abs(nyq'+1-r);
                [WMar,I]=min(WDist);  
                I
                
                O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).Wdiff=Wdiff;
                O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).Wfit=Wtest;
                O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).Wtf=Wtf;
                O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).r=r;
                O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).WDist=WDist;
                O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).WMar=WMar;
                
                O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).ConfIntDist=ConfIntDist;
                O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).MinConfIntMar=MinConfIntMar;

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

S.(ParLabel)=O.(ParLabel);

save('work2022Nyq_ipar5')
%%


%% plotting
iVal=80;
iFilt=1;
iPar=5;
ParLabel=VarLabels{iPar};
ParLabel=strcat(ParLabel,sprintf('_%d',fs/DownSampleRate));
ValLabel=sprintf('Val_%d',iVal);
FiltLabel=sprintf('Filt_%d',iFilt);

for iConf=1:1
    
    ConfLabel=sprintf('ConfIntFit_%d',iConf);

    figure(iVal)
    subplot(2,3,iFilt)
    
    L=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).Mag;
    L0=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).FittedFreqResp;

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



%%

for iPar=4:4
    
    ParLabel=VarLabels{iPar};
    ParLabel=strcat(ParLabel,sprintf('_%d',fs/DownSampleRate));
    ParValues=S.(ParLabel).ParValues;

    for iFilt=1:3

        FiltLabel=sprintf('Filt_%d',iFilt);

        for iVal=1:100

            ValLabel=sprintf('Val_%d',iVal);
            
            for iRep=1:10
                RepeatLabel=sprintf('Rep_%d',iRep);

                RepMar(iRep,iVal)=S.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).MinConfIntMar;
                WMar(iRep,iVal)=S.(ParLabel).(FiltLabel).(ValLabel).(RepeatLabel).WMar;
                
            end

            for iConf=1:3
                
                ConfLabel=sprintf('ConfIntFit_%d',iConf);
%                 ConfMar(iConf,iVal)=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).MinConfIntMar;
%                 WMar(iConf,iVal)=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).WMar;
%                 SysNum(iConf,iVal)=O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).SysTf.Numerator;
%                 SysDen(iConf,iVal)=O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).SysTf.Denominator;

            end
        end
        
        figure(iPar+1000) 
        subplot(3,1,iFilt)
        semilogx(ParValues,RepMar)
        ParLabelTitle=VarLabels{iPar};

        xlabel(ParLabelTitle)
        ylabel(FiltLabel)
        hold on
%         semilogx(ParValues,WMar)
        
    end
end
                

% x_vec=x_sampled(iFilter,:);
% y_vec=y_sampled(iFilter,:);
% 
% f_full=-fs/2:1/SimTime:fs/2;
% FreqRes=1000;
% f_vec=-ChirpFreq:2*ChirpFreq/FreqRes:ChirpFreq;
% 
% % [Sm fm]=tfestimate(x_vec,y_vec,500,[],[f_vec],fs);
% 
% L=fft(y_vec)./fft(x_vec);
% NumUniquePts = ceil(length(L)/2);
% f_onesided=f_full(NumUniquePts:end);
% L_onesided=L(1:NumUniquePts);
% 
% figure
% semilogx(f_onesided,abs(L_onesided))
% 
% % frd(L_onesided,f_onesided)
%%

function CI=ConfInt(x,a)
%Calculates conf interval of x with alpha equal to a 
a1= a/2;
a2=1-a/2;
SEM = std(x)/sqrt(length(x));         % Standard Error
ts = tinv([a1  a2],length(x)-1);      % T-Score
CI = mean(x) + ts*SEM; 

end

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


function obj=objective_f(x)
global FreqVec Wdiff nz Ts

% n=length(Wdiff);
Stf=ConstructTF(x,nz);

S=freqresp(Stf,FreqVec,'Hz');
S=squeeze(S)';

obj=sum(((abs(Wdiff)-abs(S))).^2);

end

function obj=objective_f2(x)
global FreqVec Wdiff nz Ts

nzterm=nz+1;
Stf=tf(x(1)*[1 x(2:nzterm)],[1 x(nzterm+1:end)]);


S=freqresp(Stf,FreqVec,'Hz');
S=squeeze(S)';

obj=sum(((abs(Wdiff)-abs(S))).^2);

end



% function TransFunc=ConstructTF(x,np,Ts)
% %construct a tf with given poles and zeros
% % zeros comes first in x then poles 
% num=1;
% den=1;
% 
% for iNum=1:length(x)-np
%     num=[num x(iNum)];
% end
% 
% for iDen=1:np
%     den=[den x(length(x)-np+iDen)];
% end
% 
% TransFunc=tf(num,den,Ts);
% 
% end


%%

function obj=objective_t(x)
global t x_vec y_vec
s=tf('s');
delay=exp(-0.1*s);
tfunc=tf(x(1)*[1 x(2) ],[1 x(3) x(4)])*delay;
y=lsim(tfunc,x_vec,t);
obj=sum((y_vec-y').^2);
end


function [c,ceq] = con_t(x)

% Nonlinear inequality constraints
% global S0tf

c(1:2)=real(roots([1 x(3) x(4) ]));
% c(3)=x(3)^2-4*x(4);
% c(4)=-abs(sqrt(x(2)^2-4*x(3)))+2;
%  c(4)=-x(4);
% Nonlinear equality constraints
ceq = [];


end

% function [Wtest Wtf]=sysiden(FreqResp,x0,lb,ub,FreqVec,np,nz,Ts)
% % system indentification using fmincon\
% 
% A=[];
% b=[];
% Aeq=[];
% beq=[];
% 
% [x minf outt s] = fmincon(@objective_fit,x0,A,b,Aeq,beq,lb,ub,@const_stab);
% 
% NumZeros=length(x)/2-np;
% zrs=ConvComplex(x(1:2*NumZeros));
% pls=ConvComplex(x(2*np+1:end));
% zrs
% pls
% minf
% 
% num=poly(zrs)
% den=poly(pls)
% Wtf=tf(num, den,Ts);
% Wtest=freqresp(Wtf,FreqVec,'Hz');
% Wtest=squeeze(Wtest)';
% 
% end

function obj=objective_fit(x)

global FreqVec FreqResp Estnz Ts

ZeroValues=x(1:2*Estnz);
zrs=ConvComplex(ZeroValues);
pls=ConvComplex(x(2*Estnz+1:end));

num=poly(zrs);
den=poly(pls);

Stf=tf(num,den,Ts);

S=freqresp(Stf,FreqVec,'Hz');
S=squeeze(S);

obj=sum(((abs(FreqResp)-abs(S'))).^2);

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
function [TFun,num,den]=SpecialTF(x,z)
global Ts
num=[ x(1)   ] ;
den=[ 1 -x(2)-1 x(2) ];
z;
TFun=tf(num,den,Ts);

end

function [c,ceq] = const_stab(x)

    global Estnp
%     den=1;
%     num=1;
%     Numzeros=length(x)-Estnp;
%     c=zeros(1,Numzeros+Estnp);
%     den=[den x(Numzeros+1:end)];
%     num=[num x(1:Numzeros)];
%     zers=roots(num);
%     poles=roots(den);
%     c(Numzeros+1:end)=abs(poles)-0.8;  % unit circle poles
%     c(1:Numzeros)=abs(zers)-1;      % unit circle zeros

%     pz=ConvComplex(x);
%     
%     c1=abs(pz)-1;
% %     c2=
%     c=[ c1 c1];

    % Nonlinear equality constraints
    ceq = [];
    c=[];
end

function pz_vector=ConvComplex(x)

% complex conjugate poles

npp=length(x);
a=zeros(1,npp);
b=zeros(1,npp);

if mod(npp,2)==0
    
    for i=1:2:length(x)
        
        a(i)=x(i);
        a(i+1)=x(i);
        b(i)=x(i+1);
        b(i+1)=-x(i+1);
        
    end

else
    
    for i=1:2:length(x)-1
        
        a(i)=x(i);
        a(i+1)=x(i);
        b(i)=x(i+1);
        b(i+1)=-x(i+1);
    
    end
    
    a(end)=x(end);
    b(end)=0;
   
end

ab=[ a; b];

pz_vector=complex(ab(1,:),ab(2,:));

end


function x0=Detx0(n,lb,ub)

for i=1:n
    
    
x0(i)=ub(i)+(lb(i)-ub(i))*rand(1,1);

end

end
