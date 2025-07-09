%Master_SwingTask_Scaled
% Code to simulate Swing Task perturbation responses for animal sizes from
% 1 gram to 10 tons. 

clear all;close all;clc
%% 
Exp=[-3,log10(0.005),-2,-1,0,1,2,3,log10(5000),4];ind0=find(Exp==0);
M=10.^Exp;

%==========================================================================
% % if KpIC and KdIC are [], then initial guess is KpIC=parms.I,KdIC=2*sqrt(parm.I*KdIC);
% KpIC=double.empty(length(M),0);
% KdIC=double.empty(length(M),0);
% KpIC=[];
% KdIC=[];
% KiIC=[];
%==========================================================================

%fminsearch optimization from initial values 
load('Data_SwingTask','OPvals');

KpIC=OPvals.Table(7,:);
KdIC=OPvals.Table(8,:);
iMOI=OPvals.Table(2,:);
iSMdelay=OPvals.Table(4,:)./1000;
parms.Kpmaxvec=iMOI.*(0.647./iSMdelay).^2;
parms.Kdmaxvec=4.*sqrt(parms.Kpmaxvec.*iMOI);

% 
% KpIC=OPvals.Table(7,:);
% KdIC=OPvals.Table(8,:);

% KpIC=[];
% KdIC=[];


clear OPvals;


% KpIC=KpIC+0.01*KpIC;
% KdIC=KpIC-0.01*KdIC;
% KiIC=KiIC-0.01*KiIC;

% % KiIC=KpIC;
% % KiIC(1)=1e-5;


%% 
run_opt=1;% set to 1 to run an optimziation from the above initial guesses, if set to 0, the model will be simulated with the initial guess and output response
%run_single=~(run_opt);% if optimziation is on, no single run. 
plotfig=0;% to switch on and off figure plotting within the ddeBlock function
%parms.tend=4;% End time of simulations
parms.tendvec=20*iSMdelay;% runtime of 20 Td.
%parms.tdec=1e-4;% output timestep size, only used for final time and value output, does not affect optimizations 
parms.tdec2=2000;% number of output datapoints between 0 and parms.tend,linspace
%parms.IC=[degtorad(-20),0];% initial angle in degrees (anticlock +ve). Final angle will be -(parms.init). 
%parms.IC=[degtorad(-15.13),0];% initial angle in degrees (anticlock +ve). Final angle will be -(parms.init).
parms.IC=[deg2rad(-15.03),0];% initial angle in degrees (anticlock +ve). Final angle will be -(parms.init).
parms.r=-(parms.IC(1));% reference value
parms.MovR=parms.r-parms.IC(1,1);
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
    
  
[OP,tnew,Angle,AngleV,uMusc,Ttot]=ddeSwingTask_Scaled(M(i),Kp_init,Kd_init,run_opt,parms,optimizerMethod,plotfig);

% uMusc: muscle torque before saturation 
% Ttot: total torque from acceleration*inertia from solver output. Matches uSat+uGrav. considers saturation and gravity


OPvals.Table(:,i)=OP;
OPvals.Time(i,:)=tnew;
OPvals.Angle(i,:)=Angle;
OPvals.AngleV(i,:)=AngleV;
OPvals.uMusc(i,:)=uMusc;
OPvals.Ttot(i,:)=Ttot;
end

clear OP tnew Angle AngleV uMusc Ttot

%OPvals.Tablehead={'Mass (Kg)';'MOI';'Tmusc max (Nm)';'SM delay (ms)';' Inertial delay 40deg (ms)';' ';'Kp';'Kd';'Ki';'Settling time (ms)';'Overshoot (%)';'Time constant'};
OPvals.Tablehead={'Mass (Kg)';'MOI';'Tmusc max (Nm)';'SM delay (ms)';' Inertial delay (ms)';' ';'Kp';'Kd';'Settling time (ms)';'Overshoot (%)';'Error'};


AA.Tablehead=OPvals.Tablehead;
AA.Table=OPvals.Table;
AAtable=struct2table(AA);

%% Clear and reload
%{
clear all;close all;clc
load('Data_SwingTask','OPvals','M','parms');
%}
%% Fitting data and graphing
%
close all;
MOIdat=OPvals.Table(2,:);
Tmuscdat=OPvals.Table(3,:);
KPopt=OPvals.Table(7,:);
KDopt=OPvals.Table(8,:);
STopt=OPvals.Table(9,:);


[p,S] = polyfit(log10(M),log10(KPopt),1);
Exponent.KP=p(1);
Coeff.KP=10^p(2);

[p,S] = polyfit(log10(M),log10(KDopt),1);
Exponent.KD=p(1);
Coeff.KD=10^p(2);


[p,S] = polyfit(log10(M),log10(STopt),1);
Exponent.ST=p(1);
Coeff.ST=10^p(2);

[p,S] = polyfit(log10(M),log10(MOIdat),1);
Exponent.MOIdat=p(1);
Coeff.MOIdat=10^p(2);

[p,S] = polyfit(log10(M),log10(Tmuscdat),1);
Exponent.Tsat=p(1);
Coeff.Tsat=10^p(2);

disp(['Kp: Exponent=' num2str(Exponent.KP) ' & ' 'Coefficient=' num2str(Coeff.KP) ])
disp(['Kd: Exponent=' num2str(Exponent.KD) ' & ' 'Coefficient=' num2str(Coeff.KD) ])

disp(['Settling time: Exponent=' num2str(Exponent.ST) ' & ' 'Coefficient=' num2str(Coeff.ST) ])
disp(['MOIdat: Exponent=' num2str(Exponent.MOIdat) ' & ' 'Coefficient=' num2str(Coeff.MOIdat) ])
disp(['Tmuscdat: Exponent=' num2str(Exponent.Tsat) ' & ' 'Coefficient=' num2str(Coeff.Tsat) ])

%--------------------------------------------------------------------------
% Output and figures
close all;
nam='Angle vs time';
figure('name',nam);
hold on;
for jj=1:size(OPvals.Time,1)
plot(OPvals.Time(jj,:),radtodeg(OPvals.Angle(jj,:)),'LineWidth',2)
end
%plot(OPvals.Time(1,:),0.63*radtodeg(parms.r)*ones(size(OPvals.Time(1,:))),'k--')
xlabel('time (s)')
ylabel('Angle (deg)')
%grid on;
title(nam)
xlim([0 3])

nam='Angular V vs time';
figure('name',nam);
hold on;
for jj=1:size(OPvals.Time,1)
plot(OPvals.Time(jj,:),radtodeg(OPvals.AngleV(jj,:)),'LineWidth',2)
end
xlabel('time (s)')
ylabel('Angular Velocity (deg/sec)')
grid on;
title(nam)



nam='Total torque vs time';
figure('name',nam);
hold on;
for jj=1:size(OPvals.Time,1)
plot(OPvals.Time(jj,:),radtodeg(OPvals.Ttot(jj,:)),'LineWidth',2)
end
xlabel('time (s)')
ylabel('Torque (Nm)')
grid on;
title(nam)


nam='Log log plots of Kp, Kd and Settling time';
figure('name',nam);
subplot(311)
plot(log10(M),log10(KPopt),'ro-')
ylabel('Kp')
grid on;
title(nam)
subplot(312)
plot(log10(M),log10(KDopt),'go-')
ylabel('Kd')
grid on;
subplot(313)
plot(log10(M),log10(STopt),'mo-')
grid on;
ylabel('Settling time (s)')
xlabel('log10(Mass)')
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
%--------------------------------------------------------------------------
%% Saving data
%{
t=datetime;
coderutime=toc;
notes={'From PrevDataset';
    'Change1: ';
    'Master code: Master_SwingTask_Scaled';
   'singlemass code:ddeSwingTask_Scaled';
      'Change2: ';
      'Optimizer: fminsearch';
    'Change3: '};
save('Data_SwingTaskv2');
%}
