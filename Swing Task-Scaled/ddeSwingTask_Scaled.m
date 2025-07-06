function [OP,tnew,Angle,AngleV,uMusc,Ttot]=ddeSwingTask_Scaled(M,KpIC,KdIC,run_opt,parms,optimizerMethod,plotfig)
% mod of optimizeControlBlockWDelay
%March 7 2019: paramterizing using swing task values. 
% This script runs code to simulate a sliding block under PD control with
% delay and actuator limits. It uses fmincon to search for optimal
% controller gains.
%v3: tightening stepinfo conditions. sim time same for all sizes. 
%original code, no gravity
%v1quer: trying to query dde
    % % ddeBlockSTnoGv01quer variables: 
    % % sol is states in solver timesteps
    % % solnew is states in equal spaced timesteps, tdec spacing, tnew is time vector
    % % simOnce outputs Angle, AngleV and Tsat, created from Solnew. 
    % % Ttot is torque=acceleration * Inertia, based on solver output, so saturation is considered
    % % uKp, uKd, uMusc, uSat are torque components recreated from solnew Angle, AngleV data. 
    % % uMusc does not consider saturation, uSat considers it if saturation threshold is exceeded. 
    % %v01 quer does not have Torque as an output, so outputting uMusc in its place
%v05reset: Intead of integral action, using Feedforward rest term
%ddeBlockSTGv06reset now has fmincon based on run_SayedNEWmainST_v1,
    %-David’s settling time & overshoot calc methods
    %-E for fminsearch = stepStats.SettlingTime + (stepStats.Overshoot*1e6);
    %-runtime is 20 Td
    %-optimzier tolerances match RelTol, now changing to AbsTol
%% To run as script, comment out function definition on top and uncomment below. Comment out end at line 193
%{
clear all;close all;clc
run_opt=0;

plotfig=1;% to switch on and off figure plotting
% parms.tend=4;
% parms.tdec=1e-4;% timestep size, only used for final time and value output
parms.IC=[degtorad(-15.03),0];% initial angle in degrees (anticlock +ve). Final angle will be -(parms.init). 
parms.r=-(parms.IC(1));% reference value
parms.MovR=parms.r-parms.IC(1,1);
%parms.STpc=0.001;% ratio for settling time range of final. 0.1% is 0.001.Default is 0.02
parms.STpc=0.02;% ratio for settling time range of final. 0.1% is 0.001.Default is 0.02
parms.deadz=1;% turn on/off torques for the initial delay period.

optimizerMethod = 'fminsearch';
%optimizerMethod = 'fmincon';

% set KpIC and KdIC=[] to use Kp=MOI and Kd=2sqrt(MOI*Kp);
% M=0.001;KpIC=4.44e-06;KdIC=1.26e-07;
%M=0.005;KpIC=[4.03115611637303e-05];KdIC=[1.57501817394001e-06];% 339.8M^0.20 vals
%M=0.005;KpIC=[2.973327081761118e-05];KdIC=[1.219997030281229e-06];% inteprolated from fit, because of kink at 0.005

% M=0.01;KpIC=9.43e-05; KdIC=4.35e-06;
% M=0.1;KpIC=0.00262906817028930; KdIC=0.000189239699910916;% Optimum E=0.368367321071723
% M=0.1;KpIC=0.002090; KdIC=0.000155;
%M=0.1;KpIC=0.00201093164329174; KdIC=0.000150338498085238;
%M=1;KpIC=0.042790691;KdIC=0.005188444;
% M=10;KpIC=0.910770283;KdIC=0.179092085;
%M=1;KpIC=19.42040341;KdIC=6.185140301;
% M=1000;KpIC=428.962983;KdIC=215.7495279;
%M=5000;KpIC=3677.32826;KdIC=2567.015794;
%M=5000;KpIC=3677.32826;KdIC=10567.015794;
%M=10000;KpIC=[];KdIC=[];
%M=10000;KpIC=9268.007157;KdIC=7470.401284;
%M=10000;KpIC=[6412.04632593989];KdIC=[6163.18944703624];

M=10000;
%load('DDEdataSTGv5','OPvals');
%load('DDEdataST1_SMb21v5','OPvals');
%load('DDEdataST1_SMb21v5 STpc02','OPvals');
%load('DDEdataSTnoGv8_STpc02','OPvals');
%load('DDEdataSTnoGv8_STpc02','OPvals');
%load('DDEdataSTGv19resetB','OPvals');
%load('DDEdataSTGv20reset','OPvals');
load('DDEdataSTGv24resetB','OPvals');


indM=find(abs(OPvals.Table(1,:)-M)<1e-10);
%indM=indM;
KpIC=OPvals.Table(7,indM);
KdIC=OPvals.Table(8,indM);
iMOI=OPvals.Table(2,indM);
iSMdelay=OPvals.Table(4,indM)./1000;
parms.Kpmax=iMOI*(0.647/iSMdelay)^2;
parms.Kdmax=4*sqrt(parms.Kpmax*iMOI);
parms.tend=20*iSMdelay;% runtime of 20 Td.
parms.tdec2=2000;% number of output datapoints between 0 and parms.tend,linspace
clear OPvals
%========================================================

%KpIC=[];KdIC=[];
%Kp=Kp+0.1*Kp;
%Kp=Kp-0.1*Kp;
%Kd=Kd+0.1*Kd;
%Kd=Kd-0.1*Kd;
%}
run_single=~(run_opt);
%% Swing task parameters
%{1
limbnam='Forelimb';
% Limb inertial properties, all final values in SI units, from fitlm fit to
% Kilbourne raw data
LlimbA=0.161;LlimbB=0.384;% Mammalian fore limb
Llimb=LlimbA*M.^LlimbB;
MlimbA=0.058;MlimbB=1.001;
Mlimb=MlimbA*M.^MlimbB;% 
MOIA=0.000252413780732547;MOIB=1.748;
MOI=MOIA*M.^MOIB;% 
LcomA=0.056;LcomB=0.359;
Lcom=LcomA*M.^LcomB;% 

% Muscle scaling. Alexander et. al.  (1981), Table 1 
% Assuming muscle force (100 N)=muscle mass (kg)/ fiber length (m)
MmuscA=6.2/1000;MmuscB=1.11;
Mmusc=(MmuscA*M.^MmuscB);% mass of muscle (kg). Triceps All non hopper
FLmuscA=18.7/1000;FLmuscB=0.33;
FLmusc=(FLmuscA*M.^FLmuscB);% fiber length muscle (m). Triceps All non hoppers
rho=1060;%density of muscle in kg/m^3. (Mendeys and keys 1960, Alexander 1977)
CSmusc=Mmusc./rho./FLmusc.*(100^2);% Cross sectional area of muscle in cm^2
Fmusc=CSmusc.*20;% assuming 20N/cm^2 muscle isometric stress. varies from 20-30 N/cm^2
% REfs: Allometry of muscle,tendon and elastic energy storage capacity in
% mammals Shadwick and Pollock 1994
% Close R 1972 Dynamic Properties of mammalian skeletal muscles 
AmuscA=8.7/1000;AmuscB=0.41;
Amusc=(AmuscA*M.^AmuscB);% Moment arm value from Triceps, all non hoppers Table II
Tmusc=Fmusc.*Amusc;% Torque by the muscle (Nm)
%}
%% Sensorimotor delay and its components. 
%{1
SenseD=(0.6/1000);
NCD=(5.3/1000)*M.^0.30;
SynapticD=(0.7/1000);
NMJD=(0.9/1000);
EMD=(4.3/1000)*M.^0.21;
FGD=(17.6/1000)*M.^0.20;
SMdelay=31/1000*M.^0.21;% sensorimotor delays in MS
%SMdelay=1e-20;
% Ref fig 6 proposal
Delay1=SynapticD+(NCD/2)+NMJD+EMD+FGD;
Delay2=SenseD+(NCD/2);

load('SimpendT1-ForeLimb_copy','InitVal','PowerLaw');%
indM=find(abs((-InitVal)-rad2deg(parms.r))<1e-10);
IDA=PowerLaw(indM,1);
IDB=PowerLaw(indM,2);
clear InitVal PowerLaw
Inertialdelay=IDA*M.^IDB;%for a movement from -ref to +ref
%Inertialdelay=35.5/1000*M.^0.2769;%for a movement from -20 to +20 degrees
%}

%% initial values
parms.Mass=M;parms.I=MOI;parms.Td=SMdelay;parms.umax=Tmusc;
parms.Mlimb=Mlimb;parms.Lcom=Lcom;


% Trial values
% if KpIC and KdIC are [], then initial guess is KpIC=parms.I,KdIC=2*sqrt(parm.I*KdIC);
if isempty(KpIC) && isempty(KdIC)
    parms.Kp = 1.*parms.I; parms.Kd = 2*sqrt(parms.Kp.*parms.I);%% Initial guess
else
    parms.Kp=KpIC;parms.Kd=KdIC;
end



% solver settings
parms.natPeriod = 2*pi*sqrt(parms.I/parms.Kp);
%parms.tspan = [0 5*parms.natPeriod]; % simulation start and stop time
parms.tspanhist = [0 parms.Td];
parms.tspan = [parms.Td parms.tend]; % 5 seconds is sufficient for the 10,000 kg animal
parms.solverOptions = ddeset('RelTol',1e-6, 'AbsTol',1e-9);

%% Generating history

%  optode45 = odeset('RelTol',1e-9,'AbsTol',1e-12);
%  HISTsol= ode45(@(t,y) pendEOM(t,y,parms),parms.tspanhist,parms.IC,optode45);
%  HISTsol.history=[];
%  HISTsol.discont=[];

HISTparms=parms;

HISTparms.Kp=0;
HISTparms.Kd=0;

HISTsol = dde23(@HistddeBlock,[],parms.IC,parms.tspanhist,parms.solverOptions,HISTparms);
%HISTsol = dde23(@ddeBlock,eps,parms.IC,parms.tspanhist,parms.solverOptions,HISTparms);
%HISTsol = dde23(@ddeBlock,0,parms.IC,parms.tspanhist,parms.solverOptions,parms);
%% optimize


switch optimizerMethod
    case 'fmincon'
        optimizerOptions = optimoptions('fmincon');
        optimizerOptions.Display = 'iter';
        optimizerOptions.MaxIterations = 200;
        optimizerOptions.MaxFunEvals = 1000;
        optimizerOptions.OptimalityTolerance = 1e-9;
        optimizerOptions.ConstraintTolerance = 1e-9;
        optimizerOptions.StepTolerance = 1e-9;
        objfun=@(x0) ddeObjective_fmincon(x0,parms,HISTsol);
        constrfun=@(x0) ddeConstraint_fmincon(x0,parms,HISTsol);
       optimizerOptions.Algorithm = 'interior-point';
%          optimizerOptions.Algorithm = 'trust-region-reflective';
%          optimizerOptions.Algorithm = 'sqp';
%       optimizerOptions.Algorithm = 'sqp-legacy';
%optimizerOptions.Algorithm = 'active-set';
        
        %A =[-1 0 0; 0 -1]; % consider only positive values
        %b = [0;0];
        A=[];
        b=[];
        Aeq = [];
        beq = [];
        %lb = ones(size(x0)).*1e-10; 
        lb=[1e-10,1e-10];
        ub = [parms.Kpmax,parms.Kdmax];

    case 'fminsearch'
        optimizerOptions = optimset('fminsearch');
        optimizerOptions.Display = 'iter';
        %optimizerOptions.Display = 'none';

        optimizerOptions.TolX = 1e-9;
        optimizerOptions.TolFun = 1e-9;
        optimizerOptions.MaxIterations = 1000;
        optimizerOptions.MaxFunEvals = 2500;
        objfun=@(x0)ddeObjective_fminsearch(x0,parms,HISTsol);
end

%%
% x0 = [parms.Kp]; % initital guess at the derivative gain
x0 = [parms.Kp;parms.Kd]; % initital guess at the derivative gain

% TDs = linspace(0.001,0.1,20);

if run_opt==1

    switch optimizerMethod
        case 'fmincon'
 
            % [x,fval, exitflag,output]= fmincon(objfun,x0,A,b,Aeq,beq,lb,ub,constrfun,optimizerOptions);
            [x,fval, exitflag,output]= fmincon(objfun,x0,A,b,Aeq,beq,lb,ub,constrfun,optimizerOptions);
            
        case 'fminsearch'

            [x,fval, exitflag,output]= fminsearch(objfun,x0,optimizerOptions);

    end
    
    
    %% test out optimized result
    parms.Kp = x(1);
    parms.Kd = x(2);
%     % parms.Kd = 2*sqrt(parms.Kp.*parms.I);
%     parms.natPeriod = 2*pi*sqrt(parms.I/parms.Kp);
%     parms.tspan = [0 5*parms.natPeriod]; % simulation start and stop time
    

    % parms.Kd = x(2);
[E,stepStats,tnew,Angle,AngleV,Ttot,sol,solnew] = simOnce(parms,plotfig,HISTsol);

    
OP=[M;MOI;Tmusc;SMdelay*1000;Inertialdelay*1000;0 ;parms.Kp;parms.Kd;stepStats.SettlingTime*1000;stepStats.Overshoot;E];



end% run_opt

%%
if run_single==1
  
[E,stepStats,tnew,Angle,AngleV,Ttot,sol,solnew] = simOnce(parms,plotfig,HISTsol);



OP=[M;MOI;Tmusc;SMdelay*1000;Inertialdelay*1000;0 ;parms.Kp;parms.Kd;stepStats.SettlingTime*1000;stepStats.Overshoot;E];
end


AA.Tablehead={'Mass';'MOI';'Tmusc';'SMdelay (ms)';'Inertialdelay (ms)';'';'Kp';'Kd';'Settling Time (ms)';'Overshoot(%)';'Error'};
AA.vals=OP;
AATable=struct2table(AA);


%% NOrmalized vals
Kp_normF=MOI./(SMdelay.^2);
Kd_normF=MOI./SMdelay;
Tiso_normF=(MOI.*parms.MovR)./(SMdelay.^2);

Kp_norm=parms.Kp/Kp_normF;
Kd_norm=parms.Kd/Kd_normF;
Tiso_norm=Tmusc./Tiso_normF;
Tresp_norm=stepStats.SettlingTime/SMdelay;

NormVals.Tablehead={'Mass (Kg)';'Tiso_norm';'Kp_norm';'Kd_norm';'Tresp_norm';'%OS'};
NormVals.Table=[M;Tiso_norm;Kp_norm;Kd_norm;Tresp_norm;stepStats.Overshoot];

ABTable=struct2table(NormVals);

%% Querying dde, trying to get proper torque profiles

%{
% is sol and solnew the same?
nam='sol vs solnew: Angle and Angular Velocity';
figure('name',nam)
hold on;
plot(tnew,solnew.s,'.')
plot(sol.x,sol.y,'*')
xlabel('Time (s)')
ylabel('Angle & AngVel')
legend('Angle solnew','Angvel solnew','Angle sol','AngVel sol')
grid on;

nam='sol vs solnew: Angular velocity and Acceleration';
figure('name',nam)
hold on;
plot(tnew,solnew.dts,'-')
plot(sol.x,sol.yp,'*')
xlabel('Time (s)')
ylabel('Angle & AngVel')
legend('AngVel solnew','Accel solnew','AngVel sol','Accel sol')
grid on;

%}

% recreating torque profile from solnew info, splitting Kp and Kd  contributions
%indN_Delay=round(SMdelay/parms.tdec);% How many indices long is the delay
indN_Delay=find(tnew>=SMdelay,1);
indN=length(tnew);
indN_F=indN-indN_Delay;
uKp=zeros(1,indN);
uKd=zeros(1,indN);
uFF=zeros(1,indN);

uKp(1:indN_Delay)=0;
uKd(1:indN_Delay)=0;
uFF(1:indN_Delay)=0;


uKp(indN_Delay+1:end)=parms.Kp.*(parms.r-solnew.s(1,1:indN_F));% proportional contribution
uKd(indN_Delay+1:end)=parms.Kd.*(-solnew.s(2,1:indN_F));% Derivative contribution
uFF(indN_Delay+1:end)=-(parms.Mlimb*-9.8066*parms.Lcom*sin(parms.r)).*ones(size(uFF(indN_Delay+1:end)));

uMusc=uKp+uKd+uFF;% Muscle torque before saturation removed
uSat=uMusc;
uSat(uMusc>=parms.umax) = parms.umax;
uSat(uMusc<=-parms.umax) = -parms.umax;
uGrav=parms.Mlimb*-9.8066*parms.Lcom*sin(solnew.s(1,:));
uTot=uSat+uGrav;



Ttot_error1=uTot-Ttot;
Ttot_error2=sum(abs(Ttot_error1));

%{1
if plotfig==1
    %{1
% figure('name','Components of total applied torque')
%     subplot(2,2,1)
%     hold on;
%     plot(tnew,uKp,'c-','LineWidth',1);
%     plot(tnew,uKd,'m-','LineWidth',1);
%     plot(tnew,uFF,'Color',[1 0.5 0.75],'LineWidth',1);
%     plot(tnew,uMusc,'b-','LineWidth',2);
%     xlabel('time (s)')
%     ylabel('Torque (Nm)')
%     legend('uKp','uKd','ureset','uMusc')
%     grid on;
%     title('Components of uMusc')
%     
%     subplot(2,2,2)
%     hold on;
%     plot(tnew,uMusc,'b-','LineWidth',2);
%     plot(tnew,uSat,'y--','LineWidth',2);
%     xlabel('time (s)')
%     ylabel('Torque (Nm)')
%     legend('uMusc','uSat')
%     grid on;
%     title('Is muscle torque saturating?')
%     
%     subplot(2,2,3)
%     hold on;
%     plot(tnew,uSat,'y--','LineWidth',2);
%     plot(tnew,uGrav,'Color',[1 0.5 0],'LineWidth',2);
%     plot(tnew,Ttot,'r-','LineWidth',2)
%     plot(tnew,uTot,'g--','LineWidth',2);
% 
%     xlabel('time (s)')
%     ylabel('Torque (Nm)')
%     legend('uSat','uGrav','Actual Torque','Recreated torque')
%     grid on;
%     title('Actual torque vs recreated torque')
%     
%     subplot(2,2,4)
%     hold on;
%     plot(tnew,Ttot_error1,'k.-')
%     xlabel('Time (s)')
%     ylabel('Torque error (Nm)')
%     grid on;
%     title('Error b/w actual and recreated torque')
    
%-------------------------------------------------------------------------
%close all

%{1
nam='Components of total applied torque-FBC paper';
figure('name',nam)

x0=10;
y0=0;
width=9;
height=8;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

    hold on;
    plot(tnew,uKp,'c-','LineWidth',1);
    plot(tnew,uKd,'m-','LineWidth',1);
    plot(tnew,uFF,'-','Color',[1 0.5 0.75],'LineWidth',1);

%    plot(tnew,uSat,'r-','LineWidth',2);
    plot(tnew,uMusc,'b-','LineWidth',1);
    plot(tnew,uGrav,'Color',[1 0.5 0],'LineWidth',2);
    plot(tnew,Ttot,'k-','LineWidth',2)
%         yl=ylim;
%     plot([parms.Td parms.Td],yl,'k-')
%     plot([2*parms.Td 2*parms.Td],yl,'k-')
%     plot([stepStats.SettlingTime stepStats.SettlingTime],yl,'k-')
%     grid on;
    xlabel('time (s)')
    ylabel('Torque (Nm)')
   % legend('uKp','uKd','uFF','uSat','uMusc','uGrav','Total')
    lgd=legend('uKp','uKd','uSteadystate','uMusc','uGrav','Total');
    lgd.NumColumns = 3;

    %grid on;
  %  title(nam)

%}

%{
nam='Components of total applied torque';
figure('name',nam)
    hold on;
    plot(tnew,uKp,'c--','LineWidth',1);
    plot(tnew,uKd,'m--','LineWidth',1);
    plot(tnew,uFF,'--','Color',[1 0.5 0.75],'LineWidth',1);

    plot(tnew,uSat,'r-','LineWidth',2);
    plot(tnew,uMusc,'b--','LineWidth',1);
    plot(tnew,uGrav,'Color',[1 0.5 0],'LineWidth',2);
    plot(tnew,Ttot,'k-','LineWidth',2)
    xlabel('time (s)')
    ylabel('Torque (Nm)')
    legend('uKp','uKd','uFF','uSat','uMusc','uGrav','Total')
    grid on;
    title(nam)

%}

%{
nam='Feedforward reset term';
figure('name',nam)
subplot(2,1,1)
hold on;
plot(tnew,uFF,'--','Color',[1 0.5 1],'LineWidth',1);
yl=ylim;
plot([parms.Td parms.Td],yl,'c-')
plot([2*parms.Td 2*parms.Td],yl,'m-')
xlabel('time (s)')
ylabel('Torque (Nm)')
grid on;
title('uFF')

subplot(2,1,2)
hold on;
plot(tnew,uFF,'k-','LineWidth',2)
yl=ylim;
plot([parms.Td parms.Td],yl,'c-')
plot([2*parms.Td 2*parms.Td],yl,'m-')
xlabel('time (s)')
ylabel('Feedforward reset torque')
grid on;
title('Feedforward reset torque'')

%}

%------------------------------------------------------------------------- 
%}
end  % plotfig
%% Skogestad best tuning rules gain estimates
% k2prime=1/MOI;% plant gain is inverse of MOI
% Kp_cascade=(1/k2prime)*(1/(4*(stepStats.tauC+SMdelay)^2));
% Kd_cascade=4*(stepStats.tauC+SMdelay);
% Ki_cascade=4*(stepStats.tauC+SMdelay);
% 
% Kp_ideal=Kp_cascade*(1+(Kd_cascade/Ki_cascade));
% Ki_ideal=Ki_cascade*(1+(Kd_cascade/Ki_cascade));
% Kd_ideal=Kd_cascade/(1+(Kd_cascade/Ki_cascade));
% 
% Kp_parallel=Kp_cascade*(1+(Kd_cascade/Ki_cascade));
% Ki_parallel=Kp_cascade/Ki_cascade;
% Kd_parallel=Kp_cascade*Kd_cascade;
% 
% Gainratios=[parms.Kp/Kp_parallel parms.Kd/Ki_parallel parms.Kd/Kd_parallel];


% COMMENT THIS END OUT IF RUNNING AS SCRIPT
end % function ddeBlockv03 

%%
function [E,stepStats,tnew,Angle,AngleV,Ttot,sol,solnew] = simOnce(parms,plotfig,HISTsol)
% this code simulates the response to the step input just a single time.
sol = dde23(@ddeBlock,parms.Td,HISTsol,parms.tspan,parms.solverOptions,parms);

% creating data with equally spaced time for export to Master
%tnew = [0:parms.tdec:parms.tend];
tnew=linspace(0,parms.tend,parms.tdec2);% different runtimes for each mass, 20Td
[solnew.s,solnew.dts] = deval(sol,tnew);% re evaluating sol to get state& state derivatives at regular time intervals
Angle=solnew.s(1,:);
AngleV=solnew.s(2,:);
%idx_Tc = find(Angle>=0.63*parms.r, 1, 'first');% finding the time constant
idx_Tc = find(Angle>=(parms.IC(1)+0.63*parms.MovR), 1, 'first');% finding the time constant
tauC=tnew(idx_Tc);%added to StepStats structure in stepinfo
%Tsat=solnew.dts(2,:).*parms.I;
Ttot=solnew.dts(2,:).*parms.I;% total applied torque, considering gravity,saturation etc.


% these torque calcs are wrong, torque should be calculated on delayed
% state values, anyway sol output acceleration so no need 4 this
%{ 
Torque=parms.Kp.*(parms.r-Angle) + parms.Kd.*(-AngleV);
Torque(tnew<parms.Td) = 0;
Tsat=Torque;
Tsat(Torque>=parms.umax) = parms.umax;
Tsat(Torque<=-parms.umax) = -parms.umax;

% u = parms.r-(parms.Kp).*sol.y(1,:)-(parms.Kd).*sol.y(2,:);
u = parms.Kp.*(parms.r-sol.y(1,:)) + parms.Kd.*(-sol.y(2,:));
%size(u)
u2=u;
u2(sol.x<parms.Td) = 0;
%size(u2)
usat=u2;
usat(u2>=parms.umax) = parms.umax;
usat(u2<=-parms.umax) = -parms.umax;
%}
%--------------------------------------------------------------------------

%stepStats = stepinfo(sol.y(1,:),sol.x,parms.r,'SettlingTimeThreshold',parms.STpc);% parms.r is the final value, settling time range set to 0.1% of final value

%Calculate settling time David method
        %isValid = ~any(abs(theta.Data)>theta_r*(1+maxOvershoot));
        % Compute settle time (within 2% of target):
        % Compute the last index in which the signal was above the band
        indexIsAbove = find(sol.y(1,:)>(parms.r+parms.STpc*parms.MovR),1,'last');
        if isempty(indexIsAbove)
            indexIsAbove = 1;
        end
        % Compute the last index in which the signal was below the band
        indexIsBelow = find(sol.y(1,:)<(parms.r-parms.STpc*parms.MovR),1,'last');
        if isempty(indexIsBelow)
            indexIsBelow = 1;
        end
        % Compute the last index in which the signal was outside the band
        indexIsOutside = max(indexIsAbove,indexIsBelow);
        % Check if the trial is valid and this index is not the end of the time
        % series (which also implies an invalid trial)
        if indexIsOutside<size(sol.y(1,:),2)% && isValid
           stepStats.SettlingTime= sol.x(indexIsOutside+1);
        else
            stepStats.SettlingTime = NaN;
        end
         % finding overshoot
        Peakval=max(sol.y(1,:));
        if Peakval>parms.r
            stepStats.Overshoot=((Peakval-parms.r)/parms.MovR)*100;
        else
            stepStats.Overshoot=0;
        end % finding overshoot



stepStats.tauC=tauC;
%--------------------------------------------------------------------------

E = ddeObjective_fminsearch([parms.Kp,parms.Kd],parms,HISTsol);

if plotfig==1
    %{
    nam=['Torque,AngVel,Ang of swing task. M=' num2str(parms.Mass) ' kg' '-solver timesteps'];
    figure('name',nam);
    subplot(311); hold on
    hold on;
    plot(sol.x,sol.yp(2,:)*parms.I,'k-','LineWidth',2);
    %     plot([min(sol.x) max(sol.x)],[parms.umax parms.umax],'b-','LineWidth',2)
    %     plot([min(sol.x) max(sol.x)],[-parms.umax -parms.umax],'b-','LineWidth',2)
    xlabel('time (s)');
    ylabel('Total Torque(Nm)');
    %    legend('Torque','Limit +','Limit -')
    grid on;
    yl=ylim;
    plot([parms.Td parms.Td],yl,'c-')
    plot([2*parms.Td 2*parms.Td],yl,'m-')
    plot([stepStats.SettlingTime stepStats.SettlingTime],yl)
    grid on;
    %title(nam)
    
    
    subplot(312); hold on
    plot(sol.x,rad2deg(sol.y(2,:)),'k-','LineWidth',2);
    xlabel('time (s)');
    ylabel('AngVel (deg/sec)');
    grid on;
    yl=ylim;
    plot([parms.Td parms.Td],yl,'c-')
    plot([2*parms.Td 2*parms.Td],yl,'m-')
    plot([stepStats.SettlingTime stepStats.SettlingTime],yl)
    
    subplot(313); hold on
    plot(sol.x,rad2deg(sol.y(1,:)),'k-','LineWidth',2);
    xlabel('time (s)');
    ylabel('Angle (deg)');
    grid on;
    yl=ylim;
    plot([parms.Td parms.Td],yl,'c-')
    plot([2*parms.Td 2*parms.Td],yl,'m-')
    plot([stepStats.SettlingTime stepStats.SettlingTime],yl)
    legend('Angle','Tdelay','2Tdelay','Settling time')
    
    %}
    % =========================================================================
    %{1
    nam=['FBC paper figure 4: Torque,AngVel,Ang of swing task. M=' num2str(parms.Mass) ' kg' '-solver timesteps'];
    figure('name',nam);
    
    x0=10;
    y0=0;
    width=6.5;
    height=9;
    set(gcf,'units','centimeters','position',[x0,y0,width,height])
    Ttotal=sol.yp(2,:)*parms.I;
    AngVeldat=rad2deg(sol.y(2,:));
    Angdat=rad2deg(sol.y(1,:));
    subplot(311); hold on
    hold on;
    plot(sol.x.*1000,Ttotal,'k-','LineWidth',2);
    %     plot([min(sol.x) max(sol.x)],[parms.umax parms.umax],'b-','LineWidth',2)
    %     plot([min(sol.x) max(sol.x)],[-parms.umax -parms.umax],'b-','LineWidth',2)
    %   xlabel('time (s)');
    %   ylabel('Total Torque(Nm)');
    %    legend('Torque','Limit +','Limit -')
    %   grid on;
    yl=ylim;
    plot([parms.Td*1000 parms.Td*1000],yl,'c-')
    %     plot([2*parms.Td 2*parms.Td],yl,'m-')
    plot([stepStats.SettlingTime*1000 stepStats.SettlingTime*1000],yl)
    %title(nam)
    axis ([0 parms.tend*1000 min(Ttotal)+0.1*min(Ttotal) max(Ttotal)+0.1*max(Ttotal)])
    
    subplot(312); hold on
    plot(sol.x.*1000,AngVeldat,'k-','LineWidth',2);
    %    xlabel('time (s)');
    %    ylabel('AngVel (deg/sec)');
    %    grid on;
    yl=ylim;
    plot([parms.Td*1000 parms.Td*1000],yl,'c-')
    %     plot([2*parms.Td 2*parms.Td],yl,'m-')
    plot([stepStats.SettlingTime*1000 stepStats.SettlingTime*1000],yl)
    axis ([0 parms.tend*1000 min(AngVeldat)+0.1*min(AngVeldat) max(AngVeldat)+0.1*max(AngVeldat)])
    
    subplot(313); hold on
    plot(sol.x.*1000,Angdat,'k-','LineWidth',2);
    xlabel('time (ms)');
    %    ylabel('Angle (deg)');
    %    grid on;
    (parms.r+parms.STpc*parms.MovR)
    plot([0 parms.tend*1000], rad2deg([ (parms.r+parms.STpc*parms.MovR)  (parms.r+parms.STpc*parms.MovR)]),'r--')
    plot([0 parms.tend*1000], rad2deg([ (parms.r-parms.STpc*parms.MovR)  (parms.r-parms.STpc*parms.MovR)]),'r--')
    yl=ylim;
    plot([parms.Td*1000 parms.Td*1000],yl,'c-')
    %    plot([2*parms.Td 2*parms.Td],yl,'m-')
    plot([stepStats.SettlingTime*1000 stepStats.SettlingTime*1000],yl)
    %  axis ([0 parms.tend*1000 min(Angdat)+0.1*min(Angdat) max(Angdat)+0.1*max(Angdat)])
    axis ([0 parms.tend*1000 -20 +20])
    
    %    legend('Angle','Tdelay','2Tdelay','Settling time')
    %}
    
end % plotting

end %SimOnce

%%
function E = ddeObjective_fminsearch(x0,parms,HISTsol)
% objective is to minimize the settling time
if x0(1)<=0 || x0(2)<=0
E=1e10;
else
parms.Kp = x0(1);
parms.Kd = x0(2);
% % parms.Kd = 2*sqrt(parms.Kp.*parms.I);
% parms.natPeriod = 2*pi*sqrt(parms.I/parms.Kp);
% parms.tspan = [0 5*parms.natPeriod] ;% simulation start and stop time
sol = dde23(@ddeBlock,[parms.Td],HISTsol,parms.tspan,parms.solverOptions,parms);

%--------------------------------------------------------------------------
%stepStats = stepinfo(sol.y(1,:),sol.x,parms.r,'SettlingTimeThreshold',parms.STpc);% parms.r is the final value, settling time range set to 2% of final value


%Calculate settling time David method
        %isValid = ~any(abs(theta.Data)>theta_r*(1+maxOvershoot));
        % Compute settle time (within 2% of target):
        % Compute the last index in which the signal was above the band
        indexIsAbove = find(sol.y(1,:)>(parms.r+parms.STpc*parms.MovR),1,'last');
        if isempty(indexIsAbove)
            indexIsAbove = 1;
        end
        % Compute the last index in which the signal was below the band
        indexIsBelow = find(sol.y(1,:)<(parms.r-parms.STpc*parms.MovR),1,'last');
        if isempty(indexIsBelow)
            indexIsBelow = 1;
        end
        % Compute the last index in which the signal was outside the band
        indexIsOutside = max(indexIsAbove,indexIsBelow);
        % Check if the trial is valid and this index is not the end of the time
        % series (which also implies an invalid trial)
        if indexIsOutside<size(sol.y(1,:),2)% && isValid
           stepStats.SettlingTime= sol.x(indexIsOutside+1);
        else
            stepStats.SettlingTime = NaN;
        end
         % finding overshoot
        Peakval=max(sol.y(1,:));
        if Peakval>parms.r
            stepStats.Overshoot=((Peakval-parms.r)/parms.MovR)*100;
        else
            stepStats.Overshoot=0;
        end % finding overshoot
%--------------------------------------------------------------------------
E = stepStats.SettlingTime + (stepStats.Overshoot*1e6);

end % if x0<0
%pause
end
%%
function E = ddeObjective_fmincon(x0,parms,HISTsol)
% objective is to minimize the settling time

parms.Kp = x0(1);
parms.Kd = x0(2);
% % parms.Kd = 2*sqrt(parms.Kp.*parms.I);
% parms.natPeriod = 2*pi*sqrt(parms.I/parms.Kp);
% parms.tspan = [0 5*parms.natPeriod] ;% simulation start and stop time
sol = dde23(@ddeBlock,[parms.Td],HISTsol,parms.tspan,parms.solverOptions,parms);

%--------------------------------------------------------------------------
%stepStats = stepinfo(sol.y(1,:),sol.x,parms.r,'SettlingTimeThreshold',parms.STpc);% parms.r is the final value, settling time range set to 2% of final value


%Calculate settling time David method
        %isValid = ~any(abs(theta.Data)>theta_r*(1+maxOvershoot));
        % Compute settle time (within 2% of target):
        % Compute the last index in which the signal was above the band
        indexIsAbove = find(sol.y(1,:)>(parms.r+parms.STpc*parms.MovR),1,'last');
        if isempty(indexIsAbove)
            indexIsAbove = 1;
        end
        % Compute the last index in which the signal was below the band
        indexIsBelow = find(sol.y(1,:)<(parms.r-parms.STpc*parms.MovR),1,'last');
        if isempty(indexIsBelow)
            indexIsBelow = 1;
        end
        % Compute the last index in which the signal was outside the band
        indexIsOutside = max(indexIsAbove,indexIsBelow);
        % Check if the trial is valid and this index is not the end of the time
        % series (which also implies an invalid trial)
        if indexIsOutside<size(sol.y(1,:),2)% && isValid
           stepStats.SettlingTime= sol.x(indexIsOutside+1);
        else
            stepStats.SettlingTime = NaN;
        end
%          % finding overshoot
%         Peakval=max(sol.y(1,:));
%         if Peakval>parms.r
%             stepStats.Overshoot=((Peakval-parms.r)/parms.MovR)*100;
%         else
%             stepStats.Overshoot=0;
%         end % finding overshoot
%--------------------------------------------------------------------------
E = stepStats.SettlingTime;


end
%%
%{1
function [cineq,ceq] = ddeConstraint_fmincon(x0,parms,HISTsol)
parms.Kp = x0(1);
parms.Kd = x0(2);
% parms.Kd = 2*sqrt(parms.Kp.*parms.I);

sol = dde23(@ddeBlock,[parms.Td],HISTsol,parms.tspan,parms.solverOptions,parms);

%--------------------------------------------------------------------------
%stepStats = stepinfo(sol.y(1,:),sol.x,parms.r,'SettlingTimeThreshold',parms.STpc);% parms.r is the final value, settling time range set to 2% of final value


         % finding overshoot
        Peakval=max(sol.y(1,:));
        if Peakval>parms.r
            stepStats.Overshoot=((Peakval-parms.r)/parms.MovR)*100;
        else
            stepStats.Overshoot=0;
        end % finding overshoot
%--------------------------------------------------------------------------



cineq = stepStats.Overshoot;
ceq = [];
end %ddeConstraint_fmincon
%}

%%
function dxdt = ddeBlock(t,x,Z,parms)
% this is our delay differential equation
Kp = parms.Kp;
Kd = parms.Kd;
I = parms.I;
Td = parms.Td;
r = parms.r;
Mlimb=parms.Mlimb;
Lcom=parms.Lcom;

umax = parms.umax;
umin = -umax;

% here is our delayed states
x1lag = Z(1);
x2lag = Z(2);

% if sim time is less than time delay, make the reference still equal to
% zero. This feels like a kluge, but I couldn't determine a better way to
% do it. 
if parms.deadz==1% turn off torque during the initial delay period
if t<=parms.Td, r = 0; end % note there is never a need to remove the reference. 
end

% feedforward (reset) term to counter final gravity and prevent steady state error
%uFF=0;
uFF=-(Mlimb*-9.8066*Lcom*sin(r));

% controller output force
% u = r-(Kp).*x1lag-(Kd).*x2lag;
u = Kp*(r-x1lag) + Kd*(-x2lag)+uFF;
% limit the max and min force. This appears to generate simulations that don't make
% sense. The trajectories look flattened and blocky. I am not sure why. 
if u >= umax, u = umax; elseif u <= umin, u = umin; end

if parms.deadz==1% turn off torque during the initial delay period
if t<=parms.Td, u = 0; end
end


gravT=Mlimb*-9.8066*Lcom*sin(x(1));

% state derivatives. 
dxdt(1,1) = x(2);
dxdt(2,1) = (u+gravT)/I;

%disp(['t=' num2str(t) ', x1lag=' num2str(x1lag) ', x2lag=' num2str(x2lag) ', U=' num2str(u)])

end
%%
function dxdt = HistddeBlock(t,x,Z,parms)
% this is our delay differential equation
Kp = parms.Kp;
Kd = parms.Kd;
I = parms.I;
Td = parms.Td;
r = parms.r;
Mlimb=parms.Mlimb;
Lcom=parms.Lcom;


gravT=Mlimb*-9.8066*Lcom*sin(x(1));

% state derivatives. 
dxdt(1,1) = x(2);
dxdt(2,1) = (gravT)/I;
end