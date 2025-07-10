%Master_SwingTask_ScaledFF
% Master code to pass parameters to odeSwingTask_ScaledFF and plot feedforward
% response times. 

clear all;close all;clc
%%
Exp=[-3,log10(0.005),-2,-1,0,1,2,3,log10(5000),4];ind0=find(Exp==0);
M=10.^Exp;

run_opt=0;% set to 0 to simulate initial guess/1 to optimize tswitch
parms.tend=2;% simulation max time. usualy ode event stops sim
parms.plotfig=0;% to switch on and off figure plotting within the odeSwingTask function
parms.tdec=1e-4;% decimation in data output
parms.dp=1000; % number of data points in output vectors

parms.IangleD=-15.03;% degrees
parms.XangleD=(-parms.IangleD);% target angle in degrees
optimizerMethod = 'fminsearch';
limbnam='Forelimb';

%=====================================
% Initial guess

load('Data_SwingTaskFF.mat','OPvals');%

% 
TswitchI=OPvals.Table(11,:)./1000;% 1ms after sensorimotor delays end
clear OPvals

%%

tic
for i=1:length(M)
    disp(['Mass: ' num2str(M(i))]);
    
     [OP,tnew,Angle,AngleV,uMusc,Ttot]=odeSwingTask_ScaledFF(M(i),TswitchI(i),parms,run_opt);

    
OPvals.Table(:,i)=OP;
OPvals.Time(i,:)=tnew;
OPvals.Angle(i,:)=Angle;
OPvals.AngleV(i,:)=AngleV;
OPvals.uMusc(i,:)=uMusc;
OPvals.Ttot(i,:)=Ttot;

clear OP tnew Angle AngleV uMusc Ttot
    
end

AA.Tablehead={'Mass (kg)';'Mlimb (kg)';' MOI (kg.m^2)';'Lcom (m)';'Tmusc (Nm)';'SMdelay (ms)';'Target (deg)';' ';'Final angle (deg)';'Final angvel (deg/s)';'Tswitch (ms)';'Tend (ms)'};

AA.Table=OPvals.Table;

AAtable=struct2table(AA);

%% Saving data
%{
t=datetime;
coderutime=toc;
notes={'Swing Task-Feedforward response times with freefall during delay';
   'Dataset: Data_SwingTaskFF.mat' ';
    'Master code: ';
   'singlemass code:';
    '';
    'Set Tmusc=0 to check natural time periods'};
save('STFFdataV');
%}


%% Fitting data and graphing
%
MOIdat=OPvals.Table(3,:);
Tmuscdat=OPvals.Table(5,:);
RTopt=OPvals.Table(12,:);
Tswitchopt=OPvals.Table(11,:);
Rat_tresp_tswitch=AA.Table(12,:)./AA.Table(11,:);% is tswitch half of tresp?



[p,S] = polyfit(log10(M),log10(RTopt),1);
Exponent.RT=p(1);
Coeff.RT=10^p(2);

[p,S] = polyfit(log10(M),log10(MOIdat),1);
Exponent.MOIdat=p(1);
Coeff.MOIdat=10^p(2);

[p,S] = polyfit(log10(M),log10(Tmuscdat),1);
Exponent.Tmusc=p(1);
Coeff.Tmusc=10^p(2);

[p,S] = polyfit(log10(M),log10(Tswitchopt),1);
Exponent.Tswitch=p(1);
Coeff.Tswitch=10^p(2);

% disp(['Kp: Exponent=' num2str(Exponent.KP) ' & ' 'Coefficient=' num2str(Coeff.KP) ])
% disp(['Kd: Exponent=' num2str(Exponent.KD) ' & ' 'Coefficient=' num2str(Coeff.KD) ])
% disp(['Ki: Exponent=' num2str(Exponent.KI) ' & ' 'Coefficient=' num2str(Coeff.KI) ])

disp(['Response time: Exponent=' num2str(Exponent.RT) ' & ' 'Coefficient=' num2str(Coeff.RT) ])
disp(['Tswitch: Exponent=' num2str(Exponent.Tswitch) ' & ' 'Coefficient=' num2str(Coeff.Tswitch) ])

disp(['MOIdat: Exponent=' num2str(Exponent.MOIdat) ' & ' 'Coefficient=' num2str(Coeff.MOIdat) ])
disp(['Tmuscdat: Exponent=' num2str(Exponent.Tmusc) ' & ' 'Coefficient=' num2str(Coeff.Tmusc) ])

%--------------------------------------------------------------------------
% Output and figures
close all;

 %-------------------------------------------------------------------------   
    
    nam=['Angle vs time-' limbnam];
    figure('name',nam)
    hold on;
    for n=1:length(M)
        plot(OPvals.Time(n,:),radtodeg(OPvals.Angle(n,:)))
    end
    grid on;
    xlabel('time(sec)')
    ylabel('angle(deg)')
    legend('1 gm','5 gm','10 gms','100 gms','1 kg','10 kg','100 kg','1 ton','5 ton','10 tons')
    %axis([0 max(OP.Tdelay) radtodeg(min(min(OP.angle))) 0])
    title(nam);
 %-------------------------------------------------------------------------
     nam=['Angular Velocity vs time-' limbnam];
    figure('name',nam)
    hold on;
    for n=1:length(M)
        plot(OPvals.Time(n,:),radtodeg(OPvals.AngleV(n,:)))
    end
    grid on;
    xlabel('time(sec)')
    ylabel('angular velocity (deg/s)')
    legend('1 gm','5 gm','10 gms','100 gms','1 kg','10 kg','100 kg','1 ton','5 ton','10 tons')
    axis([0 max(max(OPvals.Time)) radtodeg(min(min(OPvals.AngleV))) radtodeg(max(max(OPvals.AngleV)))])
    title(nam);   
 %-------------------------------------------------------------------------
 
  %-------------------------------------------------------------------------
     nam=['Total torque vs time-' limbnam];
    figure('name',nam)
    hold on;
    for n=1:length(M)
        plot(OPvals.Time(n,:),OPvals.Ttot(n,:))
    end
    grid on;
    xlabel('time(sec)')
    ylabel('Total Torque (Nm)')
    legend('1 gm','5 gm','10 gms','100 gms','1 kg','10 kg','100 kg','1 ton','5 ton','10 tons')
    axis([0 max(max(OPvals.Time)) (min(min(OPvals.Ttot))) (max(max(OPvals.Ttot)))])
    title(nam);   
 %-------------------------------------------------------------------------
 
 
  
nam='Log log-Response time';
figure('name',nam);
plot(log10(M),log10(RTopt),'mo-')
grid on;
ylabel('Settling time (ms)')
xlabel('log10(Mass) (kg)')
ylabel('Settling time (s)')
grid on;
%}