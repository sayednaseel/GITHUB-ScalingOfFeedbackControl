function [OPtable,OP]=run_SwingTask_Norm(KpIC,KdIC,parms,Stoptime,hws,run_opt,mdl,plotfig,optimizerMethod)
%run_SwingTask_Norm
% Code to simulate normalized PD control model and return settling times,
% overshoot and angle, angular velocity and torque profiles. Optimizes
% controller gains Kp and Kd to minimize settling time with 0 overshoot. 


%% To run as script, comment out function definition on top and uncomment below. Comment out end at line XX

%{
clear all;close all;clc;
mdl = 'SwingTaskNorm';% 
load_system(mdl)
solver_variable=1; % 1 for variable step, 2 for fixed step, modify settings within init_SwingTaskCL
init_Norm(mdl,solver_variable);% applies solver settings to mdl
hws = get_param(mdl,'modelworkspace');%handle to model workspace
hws.clear;
Stoptime=20;% simulation runtime
run_opt=1;% to run cost landscape mapping
plotfig=1;% to plot figures
optimizerMethod = 'fminsearch';
%optimizerMethod = 'fsolve';
%optimizerMethod = 'fmincon';

% Define parameters
parms.I        = 1;     % normalized
parms.t_delay  = 1;     % normalized
parms.theta_r  = 1;     % normalized
parms.theta_0  = 0;     % doesn't change result, only theta_r - theta_0 matters.
parms.dtheta_0 = 0;     % ignore for now
%parms.dtheta_0 = -0.21*sqrt(9.8066*1);% Initial angular velocity of 0.21 froude number
parms.tau_iso  = 0.063;  % Torque limit

parms.d        = 0;     % ignore for now
parms.g        = 0;     % ignore for now
parms.m        = 1;     % doesn't matter for g = 0
parms.l        = 1;     % doesn't matter for g = 0
parms.inverse  = +1;    % doesn't matter for g = 0
% Define borders: Exclude sims which exceed this limit
parms.maxOvershoot = 1e-6;%Same as solver reltol, see init code
parms.band = 0.02;% settling time band ()
parms.Kpmax=4*parms.I*(0.647/parms.t_delay)^2;
parms.Kdmax=4*sqrt(parms.Kpmax*parms.I);


% % Optimal gains:
load('TisofitST_v13','OPvals');
KpICVec=OPvals.Table(:,5);KdICVec=OPvals.Table(:,6);Ki=0;%2st column is fastest settling time, limiting overshoot to maxOvershoot
TisoVec_old=OPvals.Table(:,2);
clear OPvals;

%When seeding future sims, I just need 1 set of gains for the first Tiso
parms.tau_iso1=0.064;
%parms.tau_iso1=parms.tau_iso;
% Gains loaded from dataset
indTiso=find(abs(TisoVec_old-parms.tau_iso1)<1e-10);
%indM=indM;
KpIC=KpICVec(indTiso);
KdIC=KdICVec(indTiso);
%--------------------------------------------------------------------------
%}
%%
hws.assignin('I',parms.I);% 
hws.assignin('t_delay',parms.t_delay);% 
hws.assignin('theta_r',parms.theta_r);% 
hws.assignin('theta_0',parms.theta_0);% 
hws.assignin('dtheta_0',parms.dtheta_0);% 
hws.assignin('tau_iso',parms.tau_iso);% 
hws.assignin('d',parms.d);% 
hws.assignin('g',parms.g);% 
hws.assignin('m',parms.m);% 
hws.assignin('l',parms.l);% 
hws.assignin('inverse',parms.inverse);% 
hws.assignin('maxOvershoot',parms.maxOvershoot);% 
hws.assignin('band',parms.band);% 

Gains=[KpIC,KdIC];Ki=0;
%% Optimizer settings

switch optimizerMethod
    
    case 'fmincon'
        optimizerOptions = optimoptions('fmincon');
        optimizerOptions.Display = 'iter';
        optimizerOptions.MaxIterations = 200;
        optimizerOptions.MaxFunEvals = 1000;
        optimizerOptions.OptimalityTolerance = 1e-12;
        optimizerOptions.ConstraintTolerance = 1e-12;
        optimizerOptions.StepTolerance = 1e-12;
        objfun=@(Gains) Objective_fmincon(Gains,mdl,hws,parms);
        constrfun=@(Gains) Constraint_fmincon(Gains,mdl,hws,parms);
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
        optimizerOptions.TolX = 1e-12;
        optimizerOptions.TolFun = 1e-12;
        optimizerOptions.MaxIterations = 200;
        optimizerOptions.MaxFunEvals = 1000;
        objfun=@(Gains) Objective_fminsearch(Gains,mdl,hws,parms);

        %options = optimset('Display','iter','MaxFunEvals',1000,'MaxIter',200,'TolFun',1e-9,'TolX',1e-9);

end % switch

%% Run simulations


%Gains=Kp_init;

if run_opt==1
    warning('off','all')
    
    switch optimizerMethod
        
        case 'fminsearch'
            tic
            [GainsOPT,Error, exitflag,output]= fminsearch(objfun,Gains,optimizerOptions);
            runtime = toc;
            
        case 'fmincon'
            [GainsOPT,Error, exitflag,output]= fmincon(objfun,Gains,A,b,Aeq,beq,lb,ub,constrfun,optimizerOptions);
    end
    
else
    GainsOPT=Gains;
end


% singlerun
warning('off','all')
parms.mode=1;
[Error,OP]=Objective_fminsearch(GainsOPT,mdl,hws,parms);


%% Figures

if plotfig==1
     
     nam=['Step response-Swing task'];
     figure('name',nam);
     %clf
     hold on
     plot(OP.Time,OP.Angle1,'LineWidth',2);
     plot([0,Stoptime],[parms.theta_r,parms.theta_r])
     plot([0,Stoptime],[parms.theta_r*(1+parms.band),parms.theta_r*(1+parms.band)])
     plot([0,Stoptime],[parms.theta_r*(1-parms.band),parms.theta_r*(1-parms.band)])
     plot([OP.tSettleD OP.tSettleD],ylim,'k.-','LineWidth',2)
     grid on;
   legend('Step response','Reference','Settling band +','Settling band -','Settling time David')
     xlabel('Time')
     ylabel('Angle')

     %--------------------------------------------------------------------------
     nam=['Torque,AngVel,Ang of Swing task'];
     figure('name',nam);
     
     subplot(311); hold on
     hold on;
     plot(OP.Time,OP.Tsat,'LineWidth',2);
     ylabel('Torque (ND)');
     %    legend('Torque','Limit +','Limit -')
     grid on;
     yl=ylim;
     plot([parms.t_delay parms.t_delay],yl,'c-')
     plot([2*parms.t_delay 2*parms.t_delay],yl,'m-')
     plot([OP.tSettleD OP.tSettleD],yl)
     grid on;
     title(nam)

     
     subplot(312); hold on
     plot(OP.Time,OP.AngVel,'LineWidth',2);
     %xlabel('Time (Td)');
     ylabel('AngVel (ND)');
     grid on;
     yl=ylim;
     plot([parms.t_delay parms.t_delay],yl,'c-')
     plot([2*parms.t_delay 2*parms.t_delay],yl,'m-')
     plot([OP.tSettleD OP.tSettleD],yl)
     
     subplot(313); hold on
     plot(OP.Time,OP.Angle1,'LineWidth',2);
     xlabel('Time (Td)');
     ylabel('Angle (ND)');
     grid on;
     yl=ylim;
     plot([parms.t_delay parms.t_delay],yl,'c-')
     plot([2*parms.t_delay 2*parms.t_delay],yl,'m-')
     plot([OP.tSettleD OP.tSettleD],yl)
     legend('Angle','Tdelay','2Tdelay','Settling time')
 end % plotting

%% Output

KpOPT=GainsOPT(1);KdOPT=GainsOPT(2);
OPtable=[parms.t_delay parms.tau_iso parms.band parms.maxOvershoot  KpOPT KdOPT OP.tSettleD OP.Overshootval OP.Peakval];


AA.Tablehead={'Kp';'Ki';'Kd';'Settling Time-David';'Overshoot-David';'Settling Time-stepinfo';'Overshoot-stepinfo'};
AA.vals=[KpOPT;Ki;KdOPT;OP.tSettleD;OP.Overshootval;OP.stepStats.SettlingTime;OP.stepStats.Overshoot];
AATable=struct2table(AA);


%%
% Saving data
%{
t=datetime;
notes={'code: run_SwingTask_Norm';
       'slx: SwingTaskNorm';
       'maxOvershoot = 1e-6';
       'band = 0.02';
       'fminsearch opt';
       'Error=tSettleD+(Overshootval*1e6);';
       ''};
save('.mat');
%}
%%

%IF RUNNING AS SCRIPT, COMMENT OUT THIS END
%end % function run_SwingTask_Norm

%% 
function [Error,OP]=Objective_fminsearch(Gains,mdl,hws,parms)


Kp=Gains(1);Kd=Gains(2);Ki=0;
hws.assignin('Kp',Kp);% 
hws.assignin('Kd',Kd);% 
hws.assignin('Ki',Ki);% 

sim(mdl)


%Calculate settling time David method
        %isValid = ~any(abs(theta.Data)>theta_r*(1+maxOvershoot));
        % Compute settle time (within 2% of target):
        % Compute the last index in which the signal was above the band
        indexIsAbove = find(theta.Data>parms.theta_r*(1+parms.band),1,'last');
        if isempty(indexIsAbove)
            indexIsAbove = 1;
        end
        % Compute the last index in which the signal was below the band
        indexIsBelow = find(theta.Data<parms.theta_r*(1-parms.band),1,'last');
        if isempty(indexIsBelow)
            indexIsBelow = 1;
        end
        % Compute the last index in which the signal was outside the band
        indexIsOutside = max(indexIsAbove,indexIsBelow);
        % Check if the trial is valid and this index is not the end of the time
        % series (which also implies an invalid trial)
        if indexIsOutside<size(theta.Data,1)% && isValid
            tSettleD = theta.time(indexIsOutside+1);
        else
            tSettleD = NaN;
        end
         % finding overshoot
        Peakval=max(theta.Data);
        if Peakval>parms.theta_r
            Overshootval=((Peakval-parms.theta_r)/parms.theta_r)*100;
        else
            Overshootval=0;
        end % finding overshoot

stepStats= stepinfo(theta.Data,theta.time,parms.theta_r,'SettlingTimeThreshold',parms.band);% considering the entire angle curve
%--------------------------------------------------------------------------
Error=tSettleD+(Overshootval*1e6);
 
%-------------------------------------------------------------------------- 
OP.Time=theta.Time;
OP.Angle1=theta.Data;
OP.AngVel=thetadot.Data;
OP.Tsat=Tsat.Data;

OP.tSettleD=tSettleD;
OP.Overshootval=Overshootval;
OP.stepStats=stepStats;
OP.Peakval=Peakval;

end % Objective_fminsearch
%%
function [Error]=Objective_fmincon(Gains,mdl,hws,parms)

Kp=Gains(1);Kd=Gains(2);Ki=0;
hws.assignin('Kp',Kp);% 
hws.assignin('Kd',Kd);% 
hws.assignin('Ki',Ki);% 

sim(mdl)


%Calculate settling time David method
        %isValid = ~any(abs(theta.Data)>theta_r*(1+maxOvershoot));
        % Compute settle time (within 2% of target):
        % Compute the last index in which the signal was above the band
        indexIsAbove = find(theta.Data>parms.theta_r*(1+parms.band),1,'last');
        if isempty(indexIsAbove)
            indexIsAbove = 1;
        end
        % Compute the last index in which the signal was below the band
        indexIsBelow = find(theta.Data<parms.theta_r*(1-parms.band),1,'last');
        if isempty(indexIsBelow)
            indexIsBelow = 1;
        end
        % Compute the last index in which the signal was outside the band
        indexIsOutside = max(indexIsAbove,indexIsBelow);
        % Check if the trial is valid and this index is not the end of the time
        % series (which also implies an invalid trial)
        if indexIsOutside<size(theta.Data,1)% && isValid
            tSettleD = theta.time(indexIsOutside+1);
        else
            tSettleD = NaN;
        end
%          % finding overshoot
%         Peakval=max(theta.Data);
%         if Peakval>parms.theta_r
%             Overshootval=((Peakval-parms.theta_r)/parms.theta_r)*100;
%         else
%             Overshootval=0;
%         end % finding overshoot

%stepStats= stepinfo(theta.Data,theta.time,parms.theta_r,'SettlingTimeThreshold',parms.band);% considering the entire angle curve
%--------------------------------------------------------------------------
Error=tSettleD;%+(Overshootval*1e6);
%-------------------------------------------------------------------------- 
end % Objective_fmincon

%%
function [cineq,ceq]=Constraint_fmincon(Gains,mdl,hws,parms)


Kp=Gains(1);Kd=Gains(2);Ki=0;
hws.assignin('Kp',Kp);% 
hws.assignin('Kd',Kd);% 
hws.assignin('Ki',Ki);% 

sim(mdl)



         % finding overshoot
        Peakval=max(theta.Data);
        if Peakval>parms.theta_r
            Overshootval=((Peakval-parms.theta_r)/parms.theta_r)*100;
        else
            Overshootval=0;
        end % finding overshoot

%stepStats= stepinfo(theta.Data,theta.time,parms.theta_r,'SettlingTimeThreshold',parms.band);% considering the entire angle curve
%--------------------------------------------------------------------------
% cineq = [];
% ceq = Overshootval;

cineq = Overshootval;
ceq = [];
%-------------------------------------------------------------------------- 
end % Constraint_fmincon