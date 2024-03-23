% plotting the results for the syn and nyquist analysis 
clear all
Target=0.5;
ParamNum=4;
LegendLocation={'NorthWest','NorthWest','NorthEast','NorthEast','NorthEast'};
FileLabels={'work_paretic_noseed2', 'work_stim_noseed2',...
    'work_occ_noseed2','work_mwave_noseed3','syn_3filter3_occrate7'};
FileLabels2={'work_paretic2', 'work_stim4', 'work_occ3','work_mwave3', 'work_stimfreq'};
ParamLabels={'TrackingError','TrackingSNR','EstEffSNR','EffortCorr'};
xAxisLabels={'Volitional Paretic Hand Opening Constant',...
    'Stimulated Paretic Hand Opening Constant','Std. Dev. of k value',...
    'Std. Dev. of m-Wave Delay (s)','Stimulation Frequency (Hz)'};

ParamTitles={'TrackingError','TrackingSNR',...
    'Estimated Effort SNR','Correlation: True vs Estimated Effort'};
SynLabels={'PHand', 'SHand','Occ','DmWave','StimFreq'};
FilterLabels={'Comb','GS','Blanking'};
LineLabels={'-','--','-.'};
SynNum=5;
yLimVec=[0 0.013; 0 250; 0 35; 0 1.1];
xLimVec=[];
PChange=struct;
load('SynSimResults3')

%%
BaseLineVal=[ 50 ];

ylabels={'Tracking Error','(a.u.)';'','Tracking SNR';'Estimated','Effort SNR';'','R-Squared'}';
RepLabel=sprintf('Rep_%d',1);
for iSyn=1:SynNum
    
    FileLbl=FileLabels{iSyn};
    SynLabel=SynLabels{iSyn};

    load(FileLbl)
    clear M

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

    if iSyn==5

        for iParam=1:ParamNum

            ParamLabel=ParamLabels{iParam};

            for FilterNum=1:3   

%                 y(FilterNum,:)=M.(ParamLabel).data(FilterNum,:);
                if iParam==2
                    y(FilterNum,:)=(S.(RepLabel).StdTe(FilterNum,:)/Target).^-1;
                elseif iParam==3
                    y(FilterNum,:)=S.(RepLabel).EstEffort(FilterNum,:)./S.(RepLabel).StdEE(FilterNum,:);

                else
                    y(FilterNum,:)=S.(RepLabel).(ParamLabel)(FilterNum,:);

                end

                newy(FilterNum,:)= y(FilterNum,:)/y(FilterNum,50)*y(FilterNum,10);
                M.(ParamLabel).data(FilterNum,:)=newy(FilterNum,:);
            end
        end
    end

    fit_n=8;

    for iParam=1:ParamNum

        ParamLabel=ParamLabels{iParam};

        for FilterNum=1:3   
            
            if iParam~=4
                y(FilterNum,:)=M.(ParamLabel).data(FilterNum,:);
                x=M.wc;
            else
                
                y(FilterNum,:)=SynS.(SynLabel).RSquared.Effdata.r_sqr(FilterNum,:);
                x=M.wc;

            end
                
            
            
            
            PercentChange(FilterNum,:)=-100+100*abs(y(FilterNum,:))/abs(y(1,BaseLineVal));
            FittedVal(FilterNum,:)=nPolyFit(x,y(FilterNum,:),fit_n);
            FittedPercent(FilterNum,:)=nPolyFit(x,PercentChange(FilterNum,:),fit_n);
            
            Pmin=FittedVal(FilterNum,10)/FittedVal(FilterNum,BaseLineVal)*100-100;
            Pmax=FittedVal(FilterNum,90)/FittedVal(FilterNum,BaseLineVal)*100-100;
            PChange(iParam).(SynLabel).PercentChange(FilterNum,:)=[ Pmin Pmax];
            
        end
        
        figure(77)
        xmargin=0.1;
        ymargin=0.1;
        distx=1-2*xmargin;
        disty=1-2*ymargin;
        mrg=0.015;
        PlotLengthY=disty/(ParamNum+1);
        PlotLengthX=distx/SynNum;
        
%         subplot(5,5,iParam*5+iSyn-5)
        subplot('Position',[(xmargin+(mrg+PlotLengthX)*(iSyn-1)) (ymargin+(mrg+PlotLengthY)*(ParamNum-iParam+1)) PlotLengthX PlotLengthY]);
        
        for FilterNum=1:3
            
            % p=semilogx(x,y(FilterNum,:),LineLabels{FilterNum},'LineWidth',1,'Color','k');
            % p.Color(4) = 0.5;
            % ylim([0 0.03 ])
            
            
%             if iSyn==5 && FilterNum==2 && iParam==1
%         
%                 pright(iParam,FilterNum)=semilogx(x,0.4*FittedVal(FilterNum,:),LineLabels{FilterNum},'LineWidth',2);
%                 
%             elseif iSyn==5 && FilterNum==3 && iParam==1
%                 pright(iParam,FilterNum)=semilogx(x,0.6*FittedVal(FilterNum,:),LineLabels{FilterNum},'LineWidth',2);
%                 
%             elseif iSyn==5 && FilterNum==1 && iParam==1
%                 pright(iParam,FilterNum)=semilogx(x,0.4*FittedVal(FilterNum,:),LineLabels{FilterNum},'LineWidth',2);
%                 
%             elseif iSyn==5 && FilterNum==2 && iParam==2
%                 pright(iParam,FilterNum)=semilogx(x,5*FittedVal(FilterNum,:),LineLabels{FilterNum},'LineWidth',2);
% 
%             else
%                 pright(iParam,FilterNum)=semilogx(x,FittedVal(FilterNum,:),LineLabels{FilterNum},'LineWidth',2);
%             end
%             
            pright(iParam,FilterNum)=semilogx(x,FittedVal(FilterNum,:),LineLabels{FilterNum},'LineWidth',2);

            hold on
            
       end

        xlim ([x(10)  x(90)]  )
        ylim([ yLimVec(iParam,:)  ] )
%         a = get(gca,'Children');
%         y1data = get(a, 'YData');
%         y1min=min( [min(y1data{1}) min(y1data{2}) min(y1data{3}) ]);
%         y1max=max( [max(y1data{1}) max(y1data{2}) max(y1data{3}) ]);
%         scaledif=abs(y1max-y1min);
%         ylim([y1min-scaledif*0.2 y1max+scaledif*0.2]);
        ax = gca;
        ax.YColor = 'k';
        
        if iSyn ~= 1
            set(gca,'YTick',[])
        else
            ylabel({ylabels{1,iParam}; ylabels{2,iParam} } )

        end

        set(gca,'XTick',[])

        for FilterNum=1:3
            plot(x(BaseLineVal),FittedVal(FilterNum,BaseLineVal),'*','Color','k')
        end

        
    %     if iParam ~= 3
    %     set(gca,'XTick',[])
    %     end
         
    %     figure(1)
    %     yyaxis right
    %     nmin=1;
    %     nmax=3;
    % 
    %     for FilterNum=nmin:nmax
    %     %     semilogx(x,PercentChange(FilterNum,:),'-','LineWidth',2,'Color','k')
    %         pleft(iParam,FilterNum)=semilogx(x,FittedPercent(FilterNum,:),LineLabels{FilterNum},'LineWidth',2,'Color','k');
    %     end
    %     
    %     b = get(gca,'Children');
    %     y2data = get(b, 'YData');
    % 
    %     y2min=min( [min(y2data{1}) min(y2data{2}) min(y2data{3}) ]);
    %     y2max=max( [max(y2data{1}) max(y2data{2}) max(y2data{3}) ]);
    %     scaledif=abs(y2max-y2min);
    %     ylim([y2min-scaledif*0.2 y2max+scaledif*0.2]);
    %     title(ParamTitles{iParam});
    %     ax = gca;
    %     ax.YColor = 'k';


    end

%     subplot(5,5,4*ParamNum+iSyn-2)
%     xlabel(xAxisLabels{iSyn})

% yyaxis right
% lbl=ylabel('Change from Nominal Value (% change)');
% lbl.Position(2) = 50; 

% delete(pleft(4,1))

end

iSyn=1;
iParam=1;
subplot('Position',[(xmargin+(mrg+PlotLengthX)*(iSyn-1)) (ymargin+(mrg+PlotLengthY)*(ParamNum-1)) PlotLengthX PlotLengthY]);
legend(FilterLabels,'Location',LegendLocation{iSyn},'AutoUpdate','off');

% subplot(5,5,1)
% legend(FilterLabels,'Location',LegendLocation{iSyn},'AutoUpdate','off');



%% Nyquist plotting
%Load

load ('worksim2_allpar3')

%% Run
ParamNum=4;
SynNum=4;

LegendLocation={'NorthWest','NorthWest','NorthEast','NorthEast','NorthEast'};
xAxisLabels={'Volitional Paretic',  'Hand Opening Constant';...
    'Stimulated Paretic','Hand Opening Constant';'Std. Deviation of','k value';...
    'Std. Dev. of',' m-Wave Delay (s)';'Stimulation Frequency','(Hz)'}';
ylabels={'Tracking Error (a.u.)','Tracking SNR',...
    'Estimtad Effort SNR','Pearson Correlation'};
ParamTitles={'TrackingError','TrackingSNR',...
    'Estimated Effort SNR','Correlation: True vs Estimated Effort'};
FilterLabels={'Comb','GS','Blanking'};
LineLabels={'-','--','-.'};
ColorLabels={[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]};

iVal=80;
iFilt=1;
iPar=1;
ParLabel=VarLabels{iPar};
ParLabel=strcat(ParLabel,sprintf('_%d',fs/DownSampleRate));
ValLabel=sprintf('Val_%d',iVal);
FiltLabel=sprintf('Filt_%d',iFilt);


for iPar=1:4
    
    ParLabel=VarLabels{iPar};
    ParLabel=strcat(ParLabel,sprintf('_%d',fs/DownSampleRate));
    
    if iPar==1 
    
        ParValues=S.(ParLabel).ParValues.^-1;
    elseif iPar==3
        
        ParValues=S.(ParLabel).ParValues*0.01;

    elseif iPar==4
        low=0.01;
        high=.2;
        NumPoints=100;
        % wc=linspace(low,high,NumPoints);
        wc=(1/20/2*logspace(log10(low),log10(high),NumPoints));
        ParValues=wc;
    
    else
    
    ParValues=S.(ParLabel).ParValues;
    
    end

    for iFilt=1:3

        FiltLabel=sprintf('Filt_%d',iFilt);

        for iVal=1:100

            ValLabel=sprintf('Val_%d',iVal);

            for iConf=1:1
                
                ConfLabel=sprintf('ConfIntFit_%d',iConf);
                ConfMar(iConf,iVal)=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).MinConfIntMar;
                WMar(iConf,iVal)=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).WMar;
                Num(iVal)=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).SysTf.Numerator(1);
                Den(iVal)=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).SysTf.Denominator;

%                 SysNum(iConf,iVal)=O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).SysTf.Numerator;
%                 SysDen(iConf,iVal)=O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).SysTf.Denominator;


            end
        end
%         if iPar==1                  
%             ConfMar(1,:)=flip(ConfMar(1,:));
%         end
        figure(77) 
%         subplot(5,5,15+iPar)

        xmargin=0.1;
        ymargin=0.1;
        distx=1-2*xmargin;
        disty=1-2*ymargin;
        mrg=0.015;
        PlotLengthY=disty/(ParamNum+1);
        PlotLengthX=distx/SynNum;
        subplot('Position',[(xmargin+(mrg+PlotLengthX)*(iPar-1)) (ymargin) PlotLengthX PlotLengthY]);
        
        semilogx(ParValues(10:90),ConfMar(1,10:90),LineLabels{iFilt},'Color',ColorLabels{iFilt},'LineWidth',2);
        hold on
        semilogx(ParValues(50),ConfMar(1,50),'*','Color','k')
        ylim([ 0.96 0.98])
        xlim([ min(ParValues(10),ParValues(90)) max(ParValues(10),ParValues(90))  ])

%         subplot(5,5,20*iPar+1)
%         semilogx(ParValues,ConfMar(1,:),LineLabels{iFilt},'LineWidth',2,'Color','k');
%         hold on
%         semilogx(ParValues,WMar)
    end
    
    xlabel({xAxisLabels{1,iPar}; xAxisLabels{2,iPar} } )
    
    if iPar~=1
        set(gca,'YTick',[])
    end
    
    if iPar ==1 
        ylabel({'Nyquist Stability';'Margin'})
    end 

%     legend(FilterLabels,'Location',LegendLocation{iPar},'AutoUpdate','off');

end

%% Coeff Plotting

% num=O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).SysTf.Numerator;
% den=O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).SysTf.Denominator;
% 
% 
% plot(num)

% x=O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).SysTf.Numerator(1,1,1);
% 
% cell2math(x)


%% confidence interval

load('worksim2_allpar')
%%

lgd1={'Mean','Lower ','Upper'};
LineLabels={'-','--','-.'};

for iPar=1:1
    
    ParLabel=VarLabels{iPar};
    ParLabel=strcat(ParLabel,sprintf('_%d',fs/DownSampleRate));

    for iFilt=1:3

        FiltLabel=sprintf('Filt_%d',iFilt);

        for iVal=20:20

            ValLabel=sprintf('Val_%d',iVal);

            for iConf=1:3 

                ConfLabel=sprintf('ConfIntFit_%d',iConf);

                FittedMag=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).FittedFreqResp;
                SysTf=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).SysTf;

                figure(iVal)
                subplot(3,2,2*iFilt-1)
                semilogx(FreqVec,20*log10(abs(FittedMag)),"LineWidth",.2,"Color",ColorLabel{iConf});
                hold on

                title(FiltTitle(iFilt))
                xlabel(' Freq. (Hz) ')
                ylabel('Mag.')
                ylim([-40 10 ])
                [re,im,wout,sdre,sdim] = nyquist(SysTf);
            
                re=squeeze(re);
                im=squeeze(im);
                
                figure(iVal)
                subplot(3,2,2*iFilt)
                plot(re,im,"LineWidth",2,"Color",ColorLabel{iConf})
                xlabel('Real')
                ylabel('Imaginary')
                xlim([-1 3]);
%                 ylim([-1.5 0]);

                hold on 

            end
            
            subplot(3,2,1)
            legend(lgd1)

            
            for iConf=1:3
                
                ConfLabel=sprintf('ConfIntFit_%d',iConf);
                Mag=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).Mag;

                figure(iVal)
                subplot(3,2,2*iFilt-1)                
                semilogx(FreqVec,20*log10(abs(Mag)),'--',"Color",ColorLabel{iConf});
                hold on
                
            end
        end
    end
end
%%

load('worksim2_allpar')
%%

lgd1={'Mean','Lower ','Upper'};
LineLabels={'-','--','-.'};

for iPar=1:1
    
    ParLabel=VarLabels{iPar};
    ParLabel=strcat(ParLabel,sprintf('_%d',fs/DownSampleRate));

    for iFilt=1:3

        FiltLabel=sprintf('Filt_%d',iFilt);

        for iVal=50:50

            ValLabel=sprintf('Val_%d',iVal);

            for iConf=3:3

                ConfLabel=sprintf('ConfIntFit_%d',iConf);

                FittedMag=abs(S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).FittedFreqResp);
                Mag=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).Mag;
                SysTf=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).SysTf;

                figure(iVal)
%                 subplot(3,2,2*iFilt-1)
                subplot(3,1,1)
                semilogx(FreqVec,20*log10(abs(FittedMag)),"LineWidth",2);
                hold on

                ylabel({'(a)';'Magnitude (dB)'})
                xlabel('Frequency(Hz)')

                title('Magnitude Response')
                set(gca,'FontSize',12)


                subplot(3,1,2)
                semilogx(FreqVec,20*log10(abs(FittedMag'-Mag)),"LineWidth",2);
                hold on
                grid

                ylabel({'(b)';'|W|  (dB)'})
                xlabel('Frequency (Hz)')
                title('|W| for each Filter ')
                set(gca,'FontSize',12)

                [re,im,wout,sdre,sdim] = nyquist(SysTf);
                re=squeeze(re);
                im=squeeze(im);

                figure(iVal)
                subplot(3,1,3)
                plot(re,im,"LineWidth",2)
                xlabel('Real')
                ylabel({'(c)';'Imaginary'})
                xlim([-.05 .02]);
                ylim([-1.2 0]);

                hold on
%                 plot(-1,0,'+',"LineWidth",2,"Color",'k')
                grid on
                title('Nyquist Countours')
                set(gca,'FontSize',12)



            end
           
%             
%             subplot(3,2,1)
%             legend(lgd1)

            
%             for iConf=1:3
%                 
%                 ConfLabel=sprintf('ConfIntFit_%d',iConf);
%                 Mag=O.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).Mag;
% 
%                 figure(iVal)
%                 subplot(3,2,2*iFilt-1)                
%                 semilogx(FreqVec,20*log10(abs(Mag)),'--',"Color",ColorLabel{iConf});
%                 hold on
%                 
%             end
        end
    end
    lgd={'Comb','GS','Blanking'};
    subplot(3,1,1)
    legend(lgd,'AutoUpdate','off')

    for iFilt=1:3
        subplot(3,1,1)
        FiltLabel=sprintf('Filt_%d',iFilt);
        Mag=S.(ParLabel).(FiltLabel).(ValLabel).(ConfLabel).Mag;
        semilogx(FreqVec,18*log10(abs(Mag)),'--',"LineWidth",.5); 
        hold on 
        grid
    end
end
% 
% lbl={'(a)','(b)','(c)'};
% for i=1:1
%     subplot(3,1,i)
%     printlbl=lbl{i};
%     yt = get(gca, 'YTick');
%     h=text(-1,(max(yt)+min(yt))/2,printlbl,'FontSize', 12);
%     set(h,'Rotation',90);
% end







