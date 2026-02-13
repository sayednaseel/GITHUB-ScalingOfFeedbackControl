function [OP,tnew,Angle,AngleV,uMusc,Ttot]=ddePostureTask_Scaled(Fr,M,KpIC,KdIC,run_opt,parms,optimizerMethod,plotfig)
% %Function to take initial guess for controller gains, simulate single run,
%or optimize controller gains. Optimizes for fastest settling time (2%
%thresholds) on angular velocity, with 0 overshoot on angle curve.  
%% To run as script, comment out function definition on top and uncomment below. Comment out end at line 219
%{
clear all;close all;clc
run_opt=0;
plotfig=1;% to switch on and off figure plotting
%parms.STpc=0.001;% ratio for settling time range of final. 0.1% is 0.001.Default is 0.02
parms.STpc=0.02;% ratio for settling time range of final. 0.1% is 0.001.Default is 0.02
optimizerMethod = 'fminsearch';
%optimizerMethod = 'fmincon';
parms.deadz=1;% turn on/off torques for the initial delay period.
%mode=0;Iangle=-5;modnam='Angle Perturbation';% set mode to 0 for initial position/0 velocity  
mode=1;Fr=-0.21;modnam='Froude Perturbation';%set mode to 1 for initial velocity/0 position

% Trial values
% if KpIC and KdIC are [], then initial guess is KpIC=parms.I,KdIC=2*sqrt(parm.I*KdIC);
%M=0.001;KpIC=[];KdIC=[];

M=10000;

load('Data_PostureTask','OPvals');%

indM=find(abs(OPvals.Table(1,:)-M)<1e-10);
indM=indM;
KpIC=OPvals.Table(7,indM);
KdIC=OPvals.Table(8,indM);
iMOI=OPvals.Table(2,indM);
iSMdelay=OPvals.Table(4,indM)./1000;
parms.Kpmax=4*iMOI*(0.647/iSMdelay)^2;
parms.Kdmax=4*sqrt(parms.Kpmax*iMOI);
parms.tend=20*iSMdelay;% runtime of 20 Td.
parms.tdec2=2000;% number of output datapoints between 0 and parms.tend,linspace
%clear OPvals
%}
run_single=~(run_opt);
%% Posture task parameters
%{1
% Limb inertial properties from fitlm Type 1 regression, all final values in SI units
L_HlimbA=0.163;L_HlimbB=0.357;% Mammalian fore limb
L_Hlimb=L_HlimbA*M.^L_HlimbB;
L_FlimbA=0.161;L_FlimbB=0.384;% Mammalian fore limb
L_Flimb=L_FlimbA*M.^L_FlimbB;
%Llimb=L_Hlimb;% old method
L=(L_Hlimb+L_Flimb)./2;
MOI=M*L^2;
% Muscle scaling. Alexander et. al.  (1981), Table 1 All non-hoppers, ankle extensors
% Assuming muscle force (100 N)=muscle mass (kg)/ fiber length (m)
Mmusc=(5.1*M.^0.97)/1000;% mass of muscle (kg), using ankle extensors
FLmusc=(10.6*M.^0.14)/1000;% fiber length muscle (m)
rho=1060;%density of muscle in kg/m^3 (Mendez and keys 1960)
CSmusc=Mmusc./rho./FLmusc.*(100^2);% Cross sectional area of muscle in cm^2
Fmusc=CSmusc.*20;% assuming 20N/cm^2 muscle isometric stress. varies from 20-30 N/cm^2
% REfs: Allometry of muscle,tendon and elastic energy storage capacity in
% mammals Shadwick and Pollock 1994
% Close R 1972 Dynamic Properties of mammalian skeletal muscles 
Amusc=(9.4*M.^0.38)/1000;% Moment arm value from Table II ankle extensors-all non hoppers
Tmusc=4*Fmusc.*Amusc;% Torque by the 4 leg muscles (Nm)

mode=1;
IC=[0;0];% Initial conditions
if mode==0 % sim recovery from angle perturbation
    IC(1) = deg2rad(Iangle);
    IC(2) = 0;
elseif mode==1 % sim recovery from angular velocity perturbation
    IC(1) = 0;
    g=-9.8066;
    Vpert=Fr.*sqrt(L.*-g);
    Wpert=Vpert./L; 
    IC(2)=Wpert;
    parms.Wpert=Wpert;
end
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
%Inertialdelay=30.79/1000*M.^0.35;%for a movement from 0.21 Fr number

% %Fr=-0.21;
load('Data_Inertialdelay_PostureTask','InitVal','PowerLaw');%
indFr=find(abs((InitVal)-Fr)<1e-10);
IDA=PowerLaw(indFr,1);
IDB=PowerLaw(indFr,2);
clear InitVal PowerLaw;
Inertialdelay=IDA*M.^IDB;%for a posture correction perturbation response

%}

%% initial values

parms.Mass=M;parms.I=MOI;parms.Td=SMdelay;parms.umax=Tmusc;parms.IC=IC;parms.r=0;parms.Length=L;

% Trial values
% if KpIC and KdIC are [], then initial guess is KpIC=parms.I,KdIC=2*sqrt(parm.I*KdIC);
if isempty(KpIC) && isempty(KdIC)
    parms.Kp = 1.*parms.I; parms.Kd = 2*sqrt(parms.Kp.*parms.I);%% Initial guess
else
    parms.Kp=KpIC;parms.Kd=KdIC;
end


% solver settings
parms.natPeriod = 2*pi*sqrt(parms.I/parms.Kp);
parms.tspanhist = [0 parms.Td];
%parms.tspan = [0 5*parms.natPeriod]; % simulation start and stop time
parms.tspan = [parms.Td parms.tend]; % 5 seconds is sufficient for the 10,000 kg animal
%parms.solverOptions = ddeset('AbsTol',1e-6, 'RelTol',1e-9);
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
[E,stepStats,tnew,Angle,AngleV,Ttot,sol,solnew] = simOnce(parms,plotfig,optimizerMethod,HISTsol);
    
OP=[M;MOI;Tmusc;SMdelay*1000;Inertialdelay*1000;0 ;parms.Kp;parms.Kd;stepStats.SettlingTime*1000;stepStats.Overshoot;E];

end% run_opt

%%
if run_single==1
  
[E,stepStats,tnew,Angle,AngleV,Ttot,sol,solnew] = simOnce(parms,plotfig,optimizerMethod,HISTsol);
    

OP=[M;MOI;Tmusc;SMdelay*1000;Inertialdelay*1000;0 ;parms.Kp;parms.Kd;stepStats.SettlingTime*1000;stepStats.Overshoot;E];

end

AA.Tablehead={'Mass';'MOI';'Tmusc';'SMdelay (ms)';'Inertialdelay (ms)';'';'Kp';'Kd';'Settling Time (ms)';'%Overshoot';'Error'};
AA.vals=OP;
AATable=struct2table(AA);

%% NOrmalized vals
Kp_normF=MOI./SMdelay.^2;
Kd_normF=MOI./SMdelay;
Tiso_normF=(MOI.*abs(Wpert).*SMdelay)./(SMdelay.^2);% changed for posture task

Kp_norm=parms.Kp/Kp_normF;
Kd_norm=parms.Kd/Kd_normF;
Tresp_norm=stepStats.SettlingTime/SMdelay;
Tiso_norm=Tmusc./Tiso_normF;


NormVals.Tablehead={'Mass (Kg)';'Tiso_norm';'Kp_norm';'Kd_norm';'Tresp_norm';'%OS'};
NormVals.Table=[M;Tiso_norm;Kp_norm;Kd_norm;Tresp_norm;stepStats.Overshoot];

ABTable=struct2table(NormVals);
%% Querying dde, trying to get proper torque profiles



% recreating torque profile from solnew info, splitting Kp and Kd  contributions
%indN_Delay=round(SMdelay/parms.tdec);% How many indices long is the delay
indN_Delay=find(tnew>=SMdelay,1);
indN=length(tnew);
indN_F=indN-indN_Delay;
uKp=zeros(1,indN);
uKd=zeros(1,indN);

uKp(1:indN_Delay)=0;
uKd(1:indN_Delay)=0;
uKp(indN_Delay+1:end)=parms.Kp.*(parms.r-solnew.s(1,1:indN_F));% proportional contribution
uKd(indN_Delay+1:end)=parms.Kd.*(-solnew.s(2,1:indN_F));% Derivative contribution

uMusc=uKp+uKd;% Muscle torque before saturation removed
uSat=uMusc;
uSat(uMusc>=parms.umax) = parms.umax;
uSat(uMusc<=-parms.umax) = -parms.umax;
uGrav=parms.Mass*9.8066*parms.Length*sin(solnew.s(1,:));
uTot=uSat+uGrav;



Ttot_error1=uTot-Ttot;
Ttot_error2=sum(abs(Ttot_error1));

%{1
if plotfig==1
%-------------------------------------------------------------------------
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
    
 %   plot(tnew,uSat,'r-','LineWidth',2);
    plot(tnew,uMusc,'b-','LineWidth',1);
    plot(tnew,uGrav,'Color',[1 0.5 0],'LineWidth',2);
    plot(tnew,Ttot,'k-','LineWidth',2)
    xlabel('time (s)')
    ylabel('Torque (Nm)')
    %legend('uKp','uKd','uSat','uMusc','uGrav','Total')
    
    lgd=legend('uKp','uKd','uMusc','uGrav','Total');
    lgd.NumColumns = 3;
 %   grid on;
  %  title(nam)

%}

%-------------------------------------------------------------------------    
end  % plotfig
%}

% COMMENT THIS END OUT IF RUNNING AS SCRIPT
end % function ddeBlockv03
%%
function [E,stepStats,tnew,Angle,AngleV,Ttot,sol,solnew] = simOnce(parms,plotfig,optimizerMethod,HISTsol)
% this code simulates the response to the step input just a single time.
sol = dde23(@ddeBlock,parms.Td,HISTsol,parms.tspan,parms.solverOptions,parms);


%tnew = [0:parms.tdec:parms.tend];
tnew=linspace(0,parms.tend,parms.tdec2);% different runtimes for each mass, 20Td
[solnew.s,solnew.dts]= deval(sol,tnew);% re evaluating sol to get state at regular time intervals
Angle=solnew.s(1,:);
AngleV=solnew.s(2,:);
Ttot=solnew.dts(2,:).*parms.I;% total applied torque, considering gravity,saturation etc.


%--------------------------------------------------------------------------
%stepStats = stepinfo(sol.y(2,:),sol.x,parms.r,'SettlingTimeThreshold',parms.STpc);% parms.r is the final value, settling time range set to 2% of final value


%Calculate settling time David method
        %isValid = ~any(abs(theta.Data)>theta_r*(1+maxOvershoot));
        % Compute settle time (within 2% of target):
        % Compute the last index in which the signal was above the band
        indexIsAbove = find(sol.y(2,:)>(0+abs(parms.Wpert)*parms.STpc),1,'last');
        if isempty(indexIsAbove)
            indexIsAbove = 1;
        end
        % Compute the last index in which the signal was below the band
        indexIsBelow = find(sol.y(2,:)<(0-abs(parms.Wpert)*parms.STpc),1,'last');
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
            stepStats.Overshoot=(Peakval/deg2rad(360))*100;
        else
            stepStats.Overshoot=0;
        end % finding overshoot
        
%         indexIsAbove
%         indexIsBelow
%         stepStats.SettlingTime
%--------------------------------------------------------------------------
E = ddeObjective_fminsearch([parms.Kp,parms.Kd],parms,HISTsol);

    

if plotfig==1
    nam=['Torque,AngVel,Ang of swing task. M=' num2str(parms.Mass) ' kg'  '-solvertimesteps'];
    figure('name',nam);

    subplot(311); hold on
    hold on;
    plot(sol.x,sol.yp(2,:)*parms.I,'b.-');
    plot([min(sol.x) max(sol.x)],[parms.umax parms.umax],'k-','LineWidth',2)
   plot([min(sol.x) max(sol.x)],[-parms.umax -parms.umax],'k-','LineWidth',2)
    xlabel('time (s)');
    ylabel('Controller Effort or Torque (Nm)');
 legend('Torque','Limit +','Limit -')
    grid on;
    yl=ylim;
    plot([parms.Td parms.Td],yl,'c-')
    plot([2*parms.Td 2*parms.Td],yl,'m-')
    plot([stepStats.SettlingTime stepStats.SettlingTime],yl)
    grid on;     

    subplot(312); hold on
    plot(sol.x,rad2deg(sol.y(2,:)),'b.-');
    xlabel('time (s)');
    ylabel('Velocity (deg/sec)');
    grid on;
    yl=ylim;
    plot([parms.Td parms.Td],yl,'c-')
    plot([2*parms.Td 2*parms.Td],yl,'m-')
    plot([stepStats.SettlingTime stepStats.SettlingTime],yl)

    subplot(313); hold on
    plot(sol.x,rad2deg(sol.y(1,:)),'b.-');
    xlabel('time (s)');
    ylabel('Position (deg)');
  %  grid on;
    title(nam)
    yl=ylim;
    plot([parms.Td parms.Td],yl,'c-')
    plot([2*parms.Td 2*parms.Td],yl,'m-')
    plot([stepStats.SettlingTime stepStats.SettlingTime],yl)
    legend('Angle','Tdelay','2Tdelay','Settling time')
    %-----------------------
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
plot([0 parms.tend*1000], [(0+abs(rad2deg(parms.Wpert))*parms.STpc) (0+abs(rad2deg(parms.Wpert))*parms.STpc)],'r--')
plot([0 parms.tend*1000], [(0-abs(rad2deg(parms.Wpert))*parms.STpc) (0-abs(rad2deg(parms.Wpert))*parms.STpc)],'r--')
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
    yl=ylim;
    plot([parms.Td*1000 parms.Td*1000],yl,'c-')
%    plot([2*parms.Td 2*parms.Td],yl,'m-')
    plot([stepStats.SettlingTime*1000 stepStats.SettlingTime*1000],yl)
   axis ([0 parms.tend*1000 min(Angdat)+0.1*min(Angdat) max(Angdat)+0.1*max(Angdat)])
 %axis ([0 parms.tend*1000 -20 +20])

%    legend('Angle','Tdelay','2Tdelay','Settling time')
%}
    
end % plotting

end

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
%stepStats = stepinfo(sol.y(2,:),sol.x,parms.r,'SettlingTimeThreshold',parms.STpc);% parms.r is the final value, settling time range set to 2% of final value


%Calculate settling time David method
        %isValid = ~any(abs(theta.Data)>theta_r*(1+maxOvershoot));
        % Compute settle time (within 2% of target):
        % Compute the last index in which the signal was above the band
        indexIsAbove = find(sol.y(2,:)>(0+abs(parms.Wpert)*parms.STpc),1,'last');
        if isempty(indexIsAbove)
            indexIsAbove = 1;
        end
        % Compute the last index in which the signal was below the band
        indexIsBelow = find(sol.y(2,:)<(0-abs(parms.Wpert)*parms.STpc),1,'last');
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
            stepStats.Overshoot=(Peakval/deg2rad(360))*100;
        else
            stepStats.Overshoot=0;
        end % finding overshoot
        
%--------------------------------------------------------------------------
E = stepStats.SettlingTime + (stepStats.Overshoot*1e6);


end
%pause
end

%%
function E = ddeObjective_fmincon(x0,parms,HISTsol)
% objective is to minimize the settling time
% if x0(1)<=0 || x0(2)<=0
% E=1e10;
%else
parms.Kp = x0(1);
parms.Kd = x0(2);
% % parms.Kd = 2*sqrt(parms.Kp.*parms.I);
% parms.natPeriod = 2*pi*sqrt(parms.I/parms.Kp);
% parms.tspan = [0 5*parms.natPeriod] ;% simulation start and stop time
sol = dde23(@ddeBlock,[parms.Td],HISTsol,parms.tspan,parms.solverOptions,parms);
%--------------------------------------------------------------------------
%stepStats = stepinfo(sol.y(1,:),sol.x,parms.r,'SettlingTimeThreshold',parms.STpc);% parms.r is the final value, settling time range set to 0.1% of final value

%Calculate settling time David method
        %isValid = ~any(abs(theta.Data)>theta_r*(1+maxOvershoot));
        % Compute settle time (within 2% of target):
        % Compute the last index in which the signal was above the band
        indexIsAbove = find(sol.y(2,:)>(0+abs(parms.Wpert)*parms.STpc),1,'last');
        if isempty(indexIsAbove)
            indexIsAbove = 1;
        end
        % Compute the last index in which the signal was below the band
        indexIsBelow = find(sol.y(2,:)<(0-abs(parms.Wpert)*parms.STpc),1,'last');
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
%--------------------------------------------------------------------------
E = stepStats.SettlingTime;
%end  % if x0<0 
%pause
end
%%
function [cineq,ceq] = ddeConstraint_fmincon(x0,parms,HISTsol)
% constraint is to make overshoot zero
parms.Kp = x0(1);
parms.Kd = x0(2);
% parms.Kd = 2*sqrt(parms.Kp.*parms.I);

sol = dde23(@ddeBlock,[parms.Td],HISTsol,parms.tspan,parms.solverOptions,parms);
%--------------------------------------------------------------------------
%stepStats = stepinfo(sol.y(2,:),sol.x,parms.r,'SettlingTimeThreshold',parms.STpc);% parms.r is the final value, settling time range set to 2% of final value

         % finding overshoot
        Peakval=max(sol.y(1,:));
        if Peakval>parms.r
            stepStats.Overshoot=(Peakval/deg2rad(360))*100;
        else
            stepStats.Overshoot=0;
        end % finding overshoot
%--------------------------------------------------------------------------
cineq = stepStats.Overshoot;
ceq = [];
end
%%
function dxdt = ddeBlock(t,x,Z,parms)
% this is our delay differential equation
Kp = parms.Kp;
Kd = parms.Kd;
I = parms.I;
Td = parms.Td;
r = parms.r;% reference position is 0,0
umax = parms.umax;
umin = -umax;

% here is our delayed states
x1lag = Z(1);
x2lag = Z(2);

% if sim time is less than time delay, make the reference still equal to
% zero. 
if parms.deadz==1% turn off torque during the initial delay period
if t<=parms.Td, r = 0; end
end
%Note: since tspan=[td tend], this is not needed in STG,PTG


% controller output force
% u = r-(Kp).*x1lag-(Kd).*x2lag;
u = Kp.*(r-x1lag) + Kd.*(-x2lag);
gravT=parms.Mass*9.8066*parms.Length*sin(x(1));

% limit the max and min force. 
if u >= umax, u = umax; elseif u <= umin, u = umin; end

if parms.deadz==1% turn off torque during the initial delay period
if t<=parms.Td, u = 0; end
end


% state derivatives. 
dxdt(1,1) = x(2);
dxdt(2,1) = (u+gravT)/I;

end

%%
function dxdt = HistddeBlock(t,x,Z,parms)
% this is our delay differential equation
Kp = parms.Kp;
Kd = parms.Kd;
I = parms.I;

gravT=parms.Mass*9.8066*parms.Length*sin(x(1));
% state derivatives. 
dxdt(1,1) = x(2);
dxdt(2,1) = gravT/I;

end
