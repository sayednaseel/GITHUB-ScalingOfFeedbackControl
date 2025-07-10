function [OP,tnew,Angle,AngleV,uMusc,Ttot]=odeSwingTask_ScaledFF(M,TswitchI,parms,run_opt)
% Function to simulate swing task under feedforward control for a single
% mass. Run initial guess, or optimzie switch time to change direction of
% torque under bang-bang control

%% Test case : Comment out function definition above and uncomment below section. Aso comment out the "end" at line ~156
%{
clear all;close all;clc;
run_opt=0;
parms.tend=1;% simulation max time. usualy ode event stops sim
parms.dp=1000; % number of data points in output vectors
parms.plotfig=1;% to switch on and off figure plotting within the SimPendFF function
parms.tdec=1e-3;% decimation in data output
parms.IangleD=-15.03;% degrees
parms.XangleD=(-parms.IangleD);% target angle in degrees
optimizerMethod = 'fminsearch';

M=10000;XX=1;Tunit='(ms)';
load('Data_SwingTaskFF.mat','OPvals');


indM=find(abs(OPvals.Table(1,:)-M)<1e-10);
indM=indM;
TswitchI=OPvals.Table(11,indM)/1000-2e-7;
clear OPvals
%}
%% Scaling of input parameters
%{1
MuscRatio=1;
XX=1000;Tunit='(ms)';

limbnam='Forelimb';
% Limb inertial properties, all final values in SI units, from fitting to
% Kilborne 2013 raw data
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
Tmusc=Tmusc*MuscRatio;
%Tmusc=0;% no muscle torque, natural time period
SMdelayA=31/1000;
SMdelayB=0.21;
SMdelay=SMdelayA*M^SMdelayB;% More 2018 PRS B second sub
%SMdelay=0;% Setting delay before torque starts to 0. 

yi(1) = deg2rad(parms.IangleD);
yi(2) = 0;
Xangle=deg2rad(parms.XangleD);
%}


%%

if run_opt==1
        %optimizerOptions = optimset('fminsearch');
        optimizerOptions.Display = 'iter';
        optimizerOptions.TolX = 1e-10;
        optimizerOptions.TolFun = 1e-10;
        optimizerOptions.MaxIterations = 1000;
        optimizerOptions.MaxFunEvals = 2500;

            objfun=@(Tswitch) SimPendFFObjective(Tswitch, yi,Tmusc,Lcom,Mlimb,MOI,SMdelay,Xangle,parms);
            [TswitchF,~, exitflag,output]= fminsearch(objfun,TswitchI,optimizerOptions);

            
                        sol = pendSim(yi,Tmusc,Lcom,Mlimb,MOI,SMdelay,TswitchF,Xangle,parms);

elseif run_opt==0
    
    TswitchF=TswitchI;
                            sol = pendSim(yi,Tmusc,Lcom,Mlimb,MOI,SMdelay,TswitchF,Xangle,parms);
end

%% Outputs

t=sol.x;
y=sol.y';
tnew = [0:(t(end)/parms.dp):t(end)];% Since end times are different, need different decimation to get same number of data points

[solnew.s,solnew.dts] = deval(sol,tnew);% re evaluating sol to get state& state derivatives at regular time intervals
Angle=solnew.s(1,:);
AngleV=solnew.s(2,:);
Ttot=solnew.dts(2,:).*MOI;
Tgrav=Mlimb*-9.8066*Lcom.*sin(Angle);
uMusc=Ttot-Tgrav;

OP=[M;Mlimb;MOI;Lcom;Tmusc;SMdelay*1000;parms.XangleD;0;rad2deg(y(end,1));rad2deg(y(end,2));TswitchF*1000;t(end)*1000];


AA.Tablehead={'Mass (kg)';'Mlimb (kg)';' MOI (kg.m^2)';'Lcom (m)';'Tmusc (Nm)';'SMdelay (ms)';'Target (deg)';' ';'Final angle (deg)';'Final angvel (deg/s)';'Tswitch (ms)';'Tend (ms)'};
AA.vals=OP;
AATable=struct2table(AA);

%disp(['End states: ' num2str(radtodeg(y(end,1))) ' degrees. '  num2str(radtodeg(y(end,2))) ' degrees/sec']);

%% Graphing
% plot results



if parms.plotfig==1

    nam=['Torque,AngVel,Ang of swing task. M=' num2str(M) ' kg'];
    figure('name',nam);

    subplot(311); hold on
    hold on;
    plot(tnew,Ttot,'b-');%note ode45 sol does not have sol.yp, i.e. acceleration data. so had to use tnew
    plot(tnew,Tgrav,'g-');
    plot(tnew,uMusc,'k-');
    xlabel('time (s)');
    ylabel('Controller Effort or Torque (Nm)');
    legend('Ttotal','Tgrav','uMusc')
    grid on;
    
    subplot(312); hold on
    plot(tnew,radtodeg(AngleV));
    xlabel('time (s)');
    ylabel('Velocity (deg/sec)');
    grid on;
    yl=ylim;

    
    subplot(313); hold on
    plot(tnew,radtodeg(Angle));
    xlabel('time (s)');
    ylabel('Position (deg)');
    grid on;
    title(nam)
    yl=ylim;

end % plotting


% TO RUN AS SCRIPT, COMMENT OUT THIS END DEF
end% function definition

%% In line function definitions



function dy = pendEOM(t,y,Tmusc,Lcom,Mlimb,MOI,SMdelay,Tswitch,Xangle)
% pendulum equations of motion.

if t<SMdelay
    u=0;
elseif t<Tswitch
    u=Tmusc;
elseif t>Tswitch
    u=-Tmusc;
end

%disp(u)
dy(1,1) = y(2);
dy(2,1) = (u+Mlimb*-9.8066*Lcom*sin(y(1)))/(MOI);
end

function [position,isterminal,direction] = pendEvents(t,y,Xangle)
% the event function. We could halt the sim when either angle or velocity
% is zero, but the optimization works much better when we set the event to
% zero velocity and the objective function to zero angle.
%position = y(1)-Xangle; % The value that we want to be zero
position = y(2); % The value that we want to be zero

isterminal = 1;  % Halt integration 
direction = -1;   % zero detected when event function is decresaing
end

function sol = pendSim(yi,Tmusc,Lcom,Mlimb,MOI,SMdelay,Tswitch,Xangle,parms)
% Once we have an optimal Tswitch, redo the simulation one last time to
% inspect the results. 

% first simulate from t=0 to t=Tswitch. 
stopFunc=@(t,y) pendEvents(t,y,Xangle);

options = odeset('RelTol',1e-6,'AbsTol',1e-9,'Events',stopFunc);

sol= ode45(@(t,y) pendEOM(t,y,Tmusc,Lcom,Mlimb,MOI,SMdelay,Tswitch,Xangle), [0:parms.tdec:parms.tend], yi,options);


end


function E = SimPendFFObjective(Tswitch, yi,Tmusc,Lcom,Mlimb,MOI,SMdelay,Xangle,parms)

sol = pendSim(yi,Tmusc,Lcom,Mlimb,MOI,SMdelay,Tswitch,Xangle,parms);

t=sol.x;
y=sol.y';

E=abs(rad2deg(Xangle-y(end,1)))*1000;% Error is final angle difference in millidegrees. 


end


