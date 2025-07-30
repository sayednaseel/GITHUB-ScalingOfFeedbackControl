function [OP,tnew,Angle,AngleV,uMusc,Ttot]=odePostureTask_ScaledFF(M,TswitchI,parms,run_opt)
% Function to simulate posture task under feedforward control for a single
% mass. Run initial guess, or optimzie switch time to change direction of
% torque under bang-bang control
global g
graph=0;
%% Test case : Comment out function definition above and uncomment below section. Aso comment out the "end" at line XXX
%{
clear all;close all;clc;
graph=1;

run_opt=0;
parms.tend=2;% simulation max time. usualy ode event stops sim
parms.dp=1000; % number of data points in output vectors
parms.plotfig=1;% to switch on and off figure plotting within the SimPendFF function
parms.tdec=1e-3;% decimation in data output
optimizerMethod = 'fminsearch';

%mode=0;Iangle=-5;modnam='Angle Perturbation';% set mode to 0 for initial position/0 velocity  
mode=1;parms.Fr=-0.21;modnam='Froude Perturbation';%set mode to 1 for initial velocity/0 position


M=10000;
load('Data_PostureTaskFF.mat','OPvals');


indM=find(abs(OPvals.Table(1,:)-M)<1e-10);
indM=indM;
TswitchI=OPvals.Table(12,indM)/1000;
clear OPvals

%}
%% Scaling of input parameters
% Limb inertial properties from fitlm Type 1 regression, all final values in SI units
L_HlimbA=0.163;L_HlimbB=0.357;% Mammalian fore limb
L_Hlimb=L_HlimbA*M.^L_HlimbB;
L_FlimbA=0.161;L_FlimbB=0.384;% Mammalian fore limb
L_Flimb=L_FlimbA*M.^L_FlimbB;
L=(L_Hlimb+L_Flimb)./2;
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
B=4*Fmusc.*Amusc;% Torque by the 4 leg muscles (Nm)
SMdelayA=31/1000;
SMdelayB=0.21;
SMdelay=SMdelayA*M^SMdelayB;% More 2018 PRS B second sub
g=-9.8066;

mode=1;
if mode==0 % sim recovery from angle perturbation
    yi(1) = deg2rad(parms.Iangle);
    yi(2) = 0;
elseif mode==1 % sim recovery from angular velocity perturbation
    yi(1) = 0;
    Vpert=parms.Fr*sqrt(L.*-g);
    Wpert=Vpert./L; 
    yi(2)=Wpert;
end
%}
%%
load( 'Inertialdelay_PT','InitVal','PowerLaw');

% Initial condition of movement for feedforward inertial delay scaling values
% -0.21 Froude number is the perturbation size in the posture task for which the
% inertial delay for a 1 kg animal equals 31 ms , the sensorimotor delay
% for a 1 kg animal. 
CrossPT=0.21;
ind=find(abs(-InitVal-CrossPT)<1e-5);
IDA=PowerLaw(ind,1)*1000;
IDB=PowerLaw(ind,2);
ID=IDA*M.^IDB;% ID from Inertialdelay_PT, feedforward bang bang control. in ms

%%
if run_opt==1
        %optimizerOptions = optimset('fminsearch');
        optimizerOptions.Display = 'iter';
        optimizerOptions.TolX = 1e-15;
        optimizerOptions.TolFun = 1e-15;
        optimizerOptions.MaxIterations = 1000;
        optimizerOptions.MaxFunEvals = 2500;

         objfun=@(Tswitch) pendObjective(Tswitch,yi,B,L,M,SMdelay);
           [TswitchF,~, exitflag,output]= fminsearch(objfun,TswitchI,optimizerOptions);
         %[TswitchF,fval, exitflag,output]= fsolve(objfun,TswitchI,optimizerOptions);

            
                        [t,y] = pendSim(TswitchF,yi,B,L,M,SMdelay);

elseif run_opt==0
    
    TswitchF=TswitchI;
                           [t,y] = pendSim(TswitchF,yi,B,L,M,SMdelay);
end
%% outputs
tnew=t;
Angle=y(:,1);
AngleV=y(:,2);
% creating torque vector
indTsmdelay=find(tnew>SMdelay,1,'first');
indTswitch=find(tnew>TswitchF,1,'first');
uMusc=ones(length(t),1);
uMusc(1:indTswitch)=B;
uMusc(1:indTsmdelay)=0;

uMusc(indTswitch:end)=-B;
Tgrav=M*-g*L.*sin(y(:,1));
Ttot=uMusc+Tgrav;

% Ttot=solnew.dts(2,:).*MOI;
% Tgrav=Mlimb*-9.8066*Lcom.*sin(Angle);
% uMusc=Ttot-Tgrav;

OP=[M;L;(M.*L.^2);B;SMdelay*1000;ID;SMdelay*1000+ID;parms.Fr;0;rad2deg(y(end,1));rad2deg(y(end,2));TswitchF*1000;t(end)*1000];


AA.Tablehead={'Mass (kg)';'Height (m)';' MOI (kg.m^2)';'Tmusc (Nm)';'SMdelay (ms)';'Inertial delay (ms)';'FFRT est (ms)';'Froude no';' ';'Final angle (deg)';'Final angvel (deg/s)';'Tswitch (ms)';'Tend (ms)'};
AA.vals=OP;
AATable=struct2table(AA);

%% plot results
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

%%
% if running as script, comment out end below
end


%% In line function definitions

function f = pendObjective(Tswitch,yi,B,L,M,SMdelay)
% this is the objective function that we are trying to minimize.
[t,y] = pendSim(Tswitch,yi,B,L,M,SMdelay);

% minimize this objective function. Note that if the event crossing
% occurred, the second term in this objective function will be nearly zero.
f = abs(y(end,1))*1000;

end

function dy = pendEOM(t,y,B,L,M,SMdelay)
% pendulum equations of motion.
global g


if t<SMdelay
    B=0;
end

dy(1,1) = y(2);
dy(2,1) = (B+M*9.8066*L*sin(y(1)))/(M*L^2);
end

function [position,isterminal,direction] = pendEvents(t,y)
% the event function. We could halt the sim when either angle or velocity
% is zero, but the optimization works much better when we set the event to
% zero velocity and the objective function to zero angle.
position = y(2); % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = -1;   % The zero can be while decreasing
end

function [t,y] = pendSim(Tswitch,yi,B,L,M,SMdelay)
% Once we have an optimal Tswitch, redo the simulation one last time to
% inspect the results. 

% first simulate from t=0 to t=Tswitch. 
options = odeset('RelTol',1e-9,'AbsTol',1e-12);

[t,y] = ode45(@(t,y) pendEOM(t,y,B,L,M,SMdelay), [0 Tswitch], yi,options);

% use the end state and time as the initial state and time for the new sim
yinew = y(end,:);
tspannew = [t(end) 20];
B = -B; % switch the torque direction.

% this sim stops the sim when the angular velocity is zero, and not at a
% specific time.
options = odeset('RelTol',1e-9,'AbsTol',1e-12,'Events',@pendEvents);
[tnew,ynew,te,ye] = ode45(@(t,y) pendEOM(t,y,B,L,M,SMdelay), tspannew, yinew,options);

% final time and states are the two sim results joined together.
t = [t(1:end-1);tnew];
y = [y(1:end-1,:);ynew];

end

