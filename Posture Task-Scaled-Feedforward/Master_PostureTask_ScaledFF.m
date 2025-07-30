%Master_PostureTask_ScaledFF
% Master code to pass parameters to odePostureTask_ScaledFF and plot feedforward
% response times. 
clear all;close all;clc
%%
Exp=[-3,log10(0.005),-2,-1,0,1,2,3,log10(5000),4];ind0=find(Exp==0);
M=10.^Exp;

run_opt=0;
parms.tend=2;% simulation max time. usualy ode event stops sim
parms.plotfig=0;% to switch on and off figure plotting within the SimPendFF function
parms.tdec=1e-4;% decimation in data output
parms.dp=1000; % number of data points in output vectors


parms.Fr=-0.21;
optimizerMethod = 'fminsearch';
limbnam='Forelimb';

%=====================================
% %Initial guess-method1
load('Data_PostureTaskFF.mat','OPvals');nTs=12;%

TswitchI=OPvals.Table(nTs,:)./1000;% 1ms after sensorimotor delays end
clear OPvals

%%

tic
for i=1:length(M)
    disp(['Mass: ' num2str(M(i))]);
    
    
     [OP,tnew,Angle,AngleV,uMusc,Ttot]=odePostureTask_ScaledFF(M(i),TswitchI(i),parms,run_opt);


    
OPvals.Table(:,i)=OP;
OPvals.Time(i,:)=0:tnew(end)/1000:tnew(end);
OPvals.Angle(i,:)=interp1(tnew,Angle,OPvals.Time(i,:),'spline');
OPvals.AngleV(i,:)=interp1(tnew,AngleV,OPvals.Time(i,:),'spline');
OPvals.uMusc(i,:)=interp1(tnew,uMusc,OPvals.Time(i,:),'linear');
OPvals.Ttot(i,:)=interp1(tnew,Ttot,OPvals.Time(i,:),'linear');

clear OP tnew Angle AngleV uMusc Ttot
    
end

AA.Tablehead={'Mass (kg)';'Height (m)';' MOI (kg.m^2)';'Tmusc (Nm)';'SMdelay (ms)';'Inertial delay (ms)';'FFRT est (ms)';'Froude no';' ';'Final angle (deg)';'Final angvel (deg/s)';'Tswitch (ms)';'Tend (ms)'};

AA.Table=OPvals.Table;

AAtable=struct2table(AA);

%% Saving data
%{
t=datetime;
coderutime=toc;
notes={'Posture Task-Feedforward response times with freefall during delay';
   '';
    'Master code:';
   'singlemass code:';
    'Best dataset:';
    ''};
save('');
%}


%% Fitting data and graphing
%
MOIdat=OPvals.Table(3,:);
Ldat=OPvals.Table(2,:);
Tmuscdat=OPvals.Table(4,:);
%RTopt=OPvals.Table(11,:);
Tswitchopt=OPvals.Table(12,:);
RTopt=OPvals.Table(13,:);




[p,S] = polyfit(log10(M),log10(RTopt),1);
Exponent.RT=p(1);
Coeff.RT=10^p(2);

[p,S] = polyfit(log10(M),log10(MOIdat),1);
Exponent.MOIdat=p(1);
Coeff.MOIdat=10^p(2);

[p,S] = polyfit(log10(M),log10(Ldat),1);
Exponent.L=p(1);
Coeff.L=10^p(2);

[p,S] = polyfit(log10(M),log10(Tmuscdat),1);
Exponent.Tmusc=p(1);
Coeff.Tmusc=10^p(2);

[p,S] = polyfit(log10(M),log10(Tswitchopt),1);
Exponent.Tswitch=p(1);
Coeff.Tswitch=10^p(2);


disp(['Response time: Exponent=' num2str(Exponent.RT) ' & ' 'Coefficient=' num2str(Coeff.RT) ])
disp(['Tswitch: Exponent=' num2str(Exponent.Tswitch) ' & ' 'Coefficient=' num2str(Coeff.Tswitch) ])
disp(['Average Leg length: Exponent=' num2str(Exponent.L) ' & ' 'Coefficient=' num2str(Coeff.L) ])

disp(['MOIdat: Exponent=' num2str(Exponent.MOIdat) ' & ' 'Coefficient=' num2str(Coeff.MOIdat) ])
disp(['Tmuscdat: Exponent=' num2str(Exponent.Tmusc) ' & ' 'Coefficient=' num2str(Coeff.Tmusc) ])

%--------------------------------------------------------------------------
% Output and figures
close all;
newcolors=[240 240 220
   240 240 200
    255 220 170
    254 227 145
    254 196 79
   251 154 41
    236 112 20
    204 76 2
    153 52 4
    181 18 27];
newcolors2=newcolors./255;
 %-------------------------------------------------------------------------   
    
    nam=['Angle vs time-' limbnam];
    figure('name',nam)
    hold on;
    for n=1:length(M)
        plot(OPvals.Time(n,:),radtodeg(OPvals.Angle(n,:)),'Color',newcolors2(n,:))
    end
    %grid on;
    xlabel('time(sec)')
    ylabel('angle(deg)')
    legend('1 gm','5 gm','10 gms','100 gms','1 kg','10 kg','100 kg','1 ton','5 ton','10 tons')
    %axis([0 max(OP.Tdelay) radtodeg(min(min(OP.angle))) 0])
    title(nam);
    ylim([-8 0])
 %-------------------------------------------------------------------------
     nam=['Angular Velocity vs time-' limbnam];
    figure('name',nam)
    hold on;
    for n=1:length(M)
        plot(OPvals.Time(n,:),radtodeg(OPvals.AngleV(n,:)))
    end
    %grid on;
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
 
 
FFline=Coeff.RT.*M.^Exponent.RT;

  nam='Log log-Response time';
figure('name',nam);
hold on;
plot(log10(M),log10(RTopt),'mo')
plot(log10(M),log10(FFline),'k-')

grid on;
ylabel('Settling time (ms)')
xlabel('log10(Mass) (kg)')
ylabel('Settling time (s)')
grid on;
legend('vals','fit')
%}