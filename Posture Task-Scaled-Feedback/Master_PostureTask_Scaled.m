%Master_PostureTaks_Scaled
% Code to simulate Posture Task perturbation responses for animal sizes from
% 1 gram to 10 tons. 

clear all;close all;clc

%% Test values
% % % if KpIC and KdIC are [], then initial guess is KpIC=parms.I,KdIC=2*sqrt(parm.I*KdIC);
% %KpIC=double.empty(length(M),0);
% %KdIC=double.empty(length(M),0);
% KpIC=[];
% KdIC=[];

%fminsearch optimization from initial values 
load('Data_PostureTask.mat','OPvals');%


KpIC=OPvals.Table(7,:);
KdIC=OPvals.Table(8,:);
iMOI=OPvals.Table(2,:);
iSMdelay=OPvals.Table(4,:)./1000;
iTresp=OPvals.Table(9,:)./1000;
parms.Kpmaxvec=4*iMOI.*(0.647./iSMdelay).^2;% increasing max gains for posture task
parms.Kdmaxvec=4.*sqrt(parms.Kpmaxvec.*iMOI);
%========================================================
clear OPvals AA;

%=====================================================

%% parameters
Exp=[-3,log10(0.005),-2,-1,0,1,2,3,log10(5000),4];ind0=find(Exp==0);
M=10.^Exp;
Fr=-0.21;% Froude number of perturbation
%Fr=-0.3;% Froude number of perturbation

run_opt=0;% set to 1 to run an optimziation from the above initial guesses, if set to 0, the model will be simulated with the initial guess and output response
%run_single=~(run_opt);% if optimziation is on, no single run. 
plotfig=0;% to switch on and off figure plotting within the ddeBlock function, for each individual mass
parms.tendvec=20*iSMdelay;% runtime of 20 Td.
%parms.tdec=1e-4;% output timestep size, only used for final time and value output, does not affect optimizations
parms.tdec2=2000;% number of output datapoints between 0 and parms.tend,linspace
parms.r=0;% reference input;final angle for the pendulum is 0 degrees and 0 velocity 
%parms.STpc=0.001;% ratio for settling time range of final. 0.1% is 0.001.Default is 0.02
parms.STpc=0.02;% ratio for settling time range of final. 0.1% is 0.001.Default is 0.02
parms.deadz=1;% turn on/off torques for the initial delay period.


optimizerMethod = 'fminsearch';
%optimizerMethod = 'fmincon';
%%
tic
for i=1:length(M)
    disp(['Mass: ' num2str(M(i))]);
    
    if isempty(KpIC) && isempty(KdIC)
    Kp_init=[];
    Kd_init=[];
    else
    Kp_init=KpIC(i);
    Kd_init=KdIC(i);
    end
    parms.tend=parms.tendvec(i);
    parms.Kpmax=parms.Kpmaxvec(i);
    parms.Kdmax=parms.Kdmaxvec(i);

             [OP,tnew,Angle,AngleV,uMusc,Ttot]=ddePostureTask_Scaled(Fr,M(i),Kp_init,Kd_init,run_opt,parms,optimizerMethod,plotfig);   

    
    OPvals.Table(:,i)=OP;
    OPvals.Time(i,:)=tnew;
    OPvals.Angle(i,:)=Angle;
    OPvals.AngleV(i,:)=AngleV;
    OPvals.uMusc(i,:)=uMusc;
    OPvals.Ttot(i,:)=Ttot;
    
end
toc

clear OP tnew Angle AngleV uMusc Ttot

OPvals.Tablehead={'Mass (Kg)';'MOI';'Tmusc max (Nm)';'SM delay (s)';' Inertial delay 0.21 Fr (s)';' ';'Kp';'Kd';'Settling time (s)';'%Overshoot';'Error'};

AA.Tablehead=OPvals.Tablehead;
AA.Table=OPvals.Table;
AAtable=struct2table(AA);



%% Data graphing and processing
%{
clear all;close all;clc
load('Data_PostureTask','OPvals','M');
%}
%% Fitting data and graphing
%{1
MOIdat=OPvals.Table(2,:);
Tmuscdat=OPvals.Table(3,:);
KPopt=OPvals.Table(7,:);
KDopt=OPvals.Table(8,:);
STopt=OPvals.Table(9,:);

[p,S] = polyfit(log10(M),log10(MOIdat),1);
Exponent.MOIdat=p(1);
Coeff.MOIdat=10^p(2);

[p,S] = polyfit(log10(M),log10(Tmuscdat),1);
Exponent.Tsat=p(1);
Coeff.Tsat=10^p(2);

[p,S] = polyfit(log10(M),log10(KPopt),1);
Exponent.KP=p(1);
Coeff.KP=10^p(2);

[p,S] = polyfit(log10(M),log10(KDopt),1);
Exponent.KD=p(1);
Coeff.KD=10^p(2);

[p,S] = polyfit(log10(M),log10(STopt),1);
Exponent.ST=p(1);
Coeff.ST=10^p(2);

disp(['Kp: Exponent=' num2str(Exponent.KP) ' & ' 'Coefficient=' num2str(Coeff.KP) ])
disp(['Kd: Exponent=' num2str(Exponent.KD) ' & ' 'Coefficient=' num2str(Coeff.KD) ])
disp(['Settling time: Exponent=' num2str(Exponent.ST) ' & ' 'Coefficient=' num2str(Coeff.ST) ])
disp(['MOIdat: Exponent=' num2str(Exponent.MOIdat) ' & ' 'Coefficient=' num2str(Coeff.MOIdat) ])
disp(['Tmuscdat: Exponent=' num2str(Exponent.Tsat) ' & ' 'Coefficient=' num2str(Coeff.Tsat) ])
%--------------------------------------------------------------------------
close all;
% Output and figures
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
%colororder(newcolors)
nam='Angle vs time';
figure('name',nam);

hold on;
for jj=1:size(OPvals.Time,1)
plot(OPvals.Time(jj,:),rad2deg(OPvals.Angle(jj,:)),'LineWidth',2,'Color',newcolors2(jj,:))
end
xlabel('time (s)')
ylabel('Angle (deg)')
%grid on;
title(nam)
xlim([0 3.5])
%newcolors = {'181,18,27','#F80','#FF0','#0B0','#00F','#50F','#A0F'};


nam='Angular V vs time';
figure('name',nam);
hold on;
for jj=1:size(OPvals.Time,1)
plot(OPvals.Time(jj,:),rad2deg(OPvals.AngleV(jj,:)),'LineWidth',2)
end
xlabel('time (s)')
ylabel('Angular Velocity (deg/sec)')
grid on;
title(nam)


nam='Torque (saturated) vs time';
figure('name',nam);
hold on;
for jj=1:size(OPvals.Time,1)
plot(OPvals.Time(jj,:),rad2deg(OPvals.Ttot(jj,:)),'LineWidth',2)
end
xlabel('time (s)')
ylabel('Torque (Nm)')
grid on;
title(nam)



nam='log-log plots of Kp, Kd and Settling time';
figure('name',nam);
subplot(311)
plot(log10(M),log10(KPopt),'ko')
ylabel('log Kp')
grid on;
title(nam)
subplot(312)
plot(log10(M),log10(KDopt),'ko')
ylabel('log Kd')
grid on;
subplot(313)
plot(log10(M),log10(STopt),'ko')
grid on;
xlabel('log10(Mass)')
ylabel('log Settling time (s)')
grid on;

%}
%--------------------------------------------------------------------------
% Max torque vs isometric max

PeakTrat(1,:)=OPvals.Table(3,:);% Maximum torque produced by isometric contraction at peak activation .saturation limit
PeakTrat(2,:)=max(OPvals.uMusc,[],2)';%  Peak value from the torque with saturation vs time profile. 
PeakTrat(3,:)=PeakTrat(2,:)./PeakTrat(1,:);
nSat=find(PeakTrat(3,:)>1);% indices where peak torque >saturation
PeakTrat(4,:)=PeakTrat(3,:);
PeakTrat(4,nSat)=1;

figure('name','Max Torque/Saturation limit');
hold on;
plot(log10(M),PeakTrat(4,:),'b-','LineWidth',2)
plot(log10(M),PeakTrat(3,:),'k--')
ylabel('Fraction of Max Torque')
xlabel('log10(Mass) (kg)')
grid on;
%--------------------------------------------------------------------------
% Kd vs theoretical Kd=2*sqrt(Kp*MOI)

MOI=OPvals.Table(2,:);
Kdcrit=2*sqrt(KPopt.*MOI);

Kd_Ratio=KDopt./Kdcrit;


figure('name','Kd_OPT / 2*sqrt(Kp*MOI)');
plot(log10(M),Kd_Ratio)

ylabel('Kdopt/Kdcritical ratio')
xlabel('log10(Mass)')
%% David normalizations.


MOI=OPvals.Table(2,:);
Tiso=OPvals.Table(3,:);
SMdelay=OPvals.Table(4,:)./1000;
Tresp=OPvals.Table(9,:)./1000;
OSpc=OPvals.Table(10,:);
AngVel1=-OPvals.AngleV(:,1)';

Kp_normF=MOI./SMdelay.^2;
Kd_normF=MOI./SMdelay;
%Tiso_normF=(MOI.*parms.r)./(SMdelay.^2);
Tiso_normF=(MOI.*AngVel1.*SMdelay)./(SMdelay.^2);% changed for posture task



Kp_norm=KPopt./Kp_normF;
Kd_norm=KDopt./Kd_normF;
Tiso_norm=Tiso./Tiso_normF;
Tresp_norm=Tresp./SMdelay;

Tresp_pred_norm=ones(size(Tiso_norm)).*7.3851;%from TisofitPT_v3
Tresp_pred_norm(Tiso_norm<0.80)=1.606.*(Tiso_norm(Tiso_norm<0.80).^-1.102)+5.174;% from TisofitPT_v3 power2 fit
Tresp_pred=Tresp_pred_norm.*SMdelay;
ErrorT=(Tresp-Tresp_pred)./Tresp*100;
ErrorTmean=mean(ErrorT);

NormVals.Tablehead={'Mass (Kg)';'Tiso_norm';'Kp_norm';'Kd_norm';'Tresp_norm';'%OS';'Error(%)'};
NormVals.Table=[M;Tiso_norm;Kp_norm;Kd_norm;Tresp_norm;OSpc;ErrorT];

ABTable=struct2table(NormVals);


nam='Normalized results';
figure('name',nam)
subplot(2,2,1)
plot(log10(M),Kp_norm)
grid on;
ylabel('Kp_norm')
title(nam)
axis([-3 4 0 0.5])

subplot(2,2,2)
hold on;
plot(log10(M),Kd_norm)
grid on;
ylabel('Kd_norm')
axis([-3 4 0 1.5])


subplot(2,2,3)
hold on;
plot(log10(M),Tiso_norm)
grid on;
ylabel('Tiso_norm')
xlabel('log10(Mass)')
axis([-3 4 0 10])

subplot(2,2,4)
hold on;
plot(log10(M),Tresp_norm)
grid on;
ylabel('Tresp_norm')
xlabel('log10(Mass)')
axis([-3 4 0 15])
%--------------------------------------------------------------------------

nam='Predicted Tresp vs actual Error-Posture task';
figure('name',nam)
hold on;
plot(log10(M),Tresp,'k.-')
plot(log10(M),Tresp_pred,'r.-')
%grid on;
ylabel('Tresp (s)')
xlabel('Log Mass (kg)')
title(nam)
legend('Actual','Normalized prediction')
%axis([-3 4 0 0.5])

%}
%% Saving data
%{
t=datetime;
coderutime=toc;
notes={'From PrevDataset';
    'Change1: ';
    'Master code: Master_PostureTask_Scaled';
   'singlemass code:ddePostureTask_Scaled';
      'Change2: ';
      'Optimizer: fminsearch';
    'Change3: '};
save('Data_PostureTask')
%}
