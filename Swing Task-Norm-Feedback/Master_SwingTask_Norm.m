%Master_SwingTask_Norm
% Code to evaluate the relationship between normalized torque limits
% (tau_iso) and response time (tresp) for the normalized swing task model.
% For each tau_iso, run_SwingTask_Norm optimizes controller gains to find
% the fastest settling time without overshoot. Tries fitting different
% curves to the tau_iso vs tresp relationship.
clear all;close all;clc
%%

mdl = 'SwingTaskNorm';% delays moved to feedback pathways


load_system(mdl)
solver_variable=1; % 1 for variable step, 2 for fixed step, modify settings within init_SwingTaskCL
init_Norm(mdl,solver_variable);% applies solver settings to mdl
hws = get_param(mdl,'modelworkspace');%handle to model workspace
hws.clear;
Stoptime=20;% simulation runtime
run_opt=1;% to run optimization
plotfig=0;% to plot figures within Run code.
optimizerMethod = 'fminsearch';
%optimizerMethod = 'fsolve';
%optimizerMethod = 'fmincon';

%% Define parameters
parms.I        = 1;     % normalized
parms.t_delay  = 1;     % normalized
parms.theta_r  = 1;     % normalized
parms.theta_0  = 0;     % doesn't change result, only theta_r - theta_0 matters.
parms.dtheta_0 = 0;     % ignore for now
%parms.dtheta_0 = -0.21*sqrt(9.8066*1);% Initial angular velocity of 0.21 froude number
%parms.tau_iso  = 1e12;  % ignore for now
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
load('TisofitST_v13','OPvals_new');OPvals.Table=OPvals_new.Table;

KpICVec=OPvals.Table(:,5);KdICVec=OPvals.Table(:,6);Ki=0;%2st column is fastest settling time, limiting overshoot to maxOvershoot
TisoVec_old=OPvals.Table(:,2);
clear OPvals;

%%

tnew=0:1e-1:Stoptime;
TisoVec=0.200:-0.001:0.001;% vector of torque limits

tic
for i=1:length(TisoVec)
disp(['The Torque limit is ' num2str(TisoVec(i))]);
parms.tau_iso=TisoVec(i);

% Gains loaded from dataset
indTiso=find(abs(TisoVec_old-parms.tau_iso)<1e-10);
%indM=indM;
KpIC=KpICVec(indTiso);
KdIC=KdICVec(indTiso);


[OPtable1,OP1]=run_SwingTask_Norm(KpIC,KdIC,parms,Stoptime,hws,run_opt,mdl,plotfig,optimizerMethod);

OPvals.Table(i,:)=OPtable1;
OPvals.Time(i,:)=tnew;
OPvals.Angle(i,:)=interp1(OP1.Time,OP1.Angle1,tnew);
OPvals.AngVel(i,:)=interp1(OP1.Time,OP1.AngVel,tnew);
OPvals.Tsat(i,:)=interp1(OP1.Time,OP1.Tsat,tnew);
end 

% Note response time is settling time+SMdelay, since simulation starts when
% reference is 1, ignoring the initial delay period. this prevents the
% impulse due to differentiating a step input
AAOUT=table(OPvals.Table(:,1),OPvals.Table(:,2),OPvals.Table(:,3),OPvals.Table(:,4),OPvals.Table(:,5),OPvals.Table(:,6),OPvals.Table(:,7),OPvals.Table(:,8),OPvals.Table(:,9),OPvals.Table(:,1)+OPvals.Table(:,7));
AAOUT.Properties.VariableNames={'SMdelay','Tiso' ,'band' ,'maxOvershoot', 'KpOPT', 'KdOPT', 'Tsettling' ,'Overshoot (%)' ,'Peak','Tresp'};

%% Output and figures
close all
nam='Angle vs time-normalized model';
figure('name',nam);
hold on;
plot(OPvals.Time(1,:),(OPvals.Angle),'LineWidth',2)
xlabel('time (Td)')
ylabel('Angle')
grid on;
title(nam)
colormap(copper)
%legend(num2str(TisoVec'))

Curve1=0.0790;% highlight 1 curve
indTiso2=find(abs(OPvals.Table(:,2)-Curve1)<1e-10);
% %indM=indM;
plot(OPvals.Time(indTiso2,:),(OPvals.Angle(indTiso2,:)),'k-','LineWidth',2)
grid on;
%title(['Tisofit_v1b:' num2str(Curve1)])


%--------------------------------------------------------------------------
nam='Angular V vs time-normalized model';
figure('name',nam);
hold on;
plot(OPvals.Time(1,:),(OPvals.AngVel),'LineWidth',2)
xlabel('Time (Td)')
ylabel('Angular Velocity ')
grid on;
title(nam)
%legend(num2str(TisoVec'))

%indTiso2=find(abs(OPvals.Table(:,2)-Curve1)<1e-10);
% %indM=indM;
plot(OPvals.Time(indTiso2,:),(OPvals.AngVel(indTiso2,:)),'k-','LineWidth',2)
grid on;
%title(['Tisofit_v1b:' num2str(Curve1)])

%--------------------------------------------------------------------------


nam='Torque (saturated) vs time-normalized model';
figure('name',nam);
hold on;
plot(OPvals.Time(1,:),OPvals.Tsat,'LineWidth',2)
xlabel('time (Td)')
ylabel('Torque')
grid on;
title(nam)
%legend(num2str(TisoVec'))

%indTiso2=find(abs(OPvals.Table(:,2)-Curve1)<1e-10);
% %indM=indM;
plot(OPvals.Time(indTiso2,:),(OPvals.Tsat(indTiso2,:)),'k-','LineWidth',2)
grid on;
%title(['Tisofit_v1b:' num2str(Curve1)])

%--------------------------------------------------------------------------
nam='Torque (saturated) vs time-normalized model-3 regions';
figure('name',nam);
hold on;
%plot(OPvals.Time(1,:),OPvals.Tsat,'LineWidth',2)
xlabel('time (N.D.)')
ylabel('Torque (N.D.)')
grid on;
title(nam)
%legend(num2str(TisoVec'))

TisoAval=0.20;
TisoBval=0.10;
TisoCval=0.05;

indTisoA=find(abs(OPvals.Table(:,2)-TisoAval)<1e-10);
indTisoB=find(abs(OPvals.Table(:,2)-TisoBval)<1e-10);
indTisoC=find(abs(OPvals.Table(:,2)-TisoCval)<1e-10);

% %indM=indM;
plot(OPvals.Time(indTisoA,:),(OPvals.Tsat(indTisoA,:)),'k-','LineWidth',2)
plot(OPvals.Time(indTisoB,:),(OPvals.Tsat(indTisoB,:)),'r-','LineWidth',2)
plot(OPvals.Time(indTisoC,:),(OPvals.Tsat(indTisoC,:)),'b-','LineWidth',2)
grid on;
legend(['Tau iso:' num2str(TisoAval)],['Tau iso:' num2str(TisoBval)],['Tau iso:' num2str(TisoCval)])
%title(['Tisofit_v1b:' num2str(Curve1)])
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Expanding dataset to show delay-limited range
Data1=OPvals.Table(1,:);
OPvals.TableB=repmat(Data1,5,1);
OPvals.TableB(1:5,2)=[0.21:0.01:0.25];
OPvals.TableC=[OPvals.TableB;OPvals.Table];
%--------------------------------------------------------------------------
%Tnosat=0.4;
Line1=0.1610;% When Tiso limits clips top part of Torque curve
Line2=0.0790;% When Tiso limits clip both top and bottom parts of torque curve.
clip=10;% REmoving the last few Tiso values which had Nan settling times
nam='Normalized model response for increasing Tau iso-Step response';
figure('name',nam);
subplot(2,1,1)
hold on;
% plot(OPvals.Table(:,1),OPvals.Table(:,2),'r.-')
plot(OPvals.TableC(1:end-clip,2),OPvals.TableC(1:end-clip,7),'r-')
%yl=ylim;
%plot([Tnosat Tnosat],yl)
ylabel('Settling time')
%legend()
title(nam)
%grid on;
%axis([min(OPvals.TableC(:,1)) max(OPvals.TableC(:,1)) 0 max(OPvals.TableC(:,2))+0.1*max(OPvals.TableC(:,2))])
yl=ylim;
% plot([0.1610 0.1610],yl,'b-')
% plot([0.0710 0.0710],yl,'b-')
plot([Line1 Line1],yl,'b-')
plot([Line2 Line2],yl,'b-')


subplot(2,1,2)
hold on;
%plot(OPvals.Table(1:end-0,2),OPvals.Table(1:end-0,8),'k.-')
plot(OPvals.TableC(1:end-clip,2),OPvals.TableC(1:end-clip,8),'k-')
ylabel('Overshoot')
%grid on;
%yl=ylim;
%plot([Tnosat Tnosat],yl)
axis([min(OPvals.TableC(:,2)) max(OPvals.TableC(:,2)) -0.5 1])

yl=ylim;
plot([Line1 Line1],yl,'b-')
plot([Line2 Line2],yl,'b-')
 xlabel('Tau iso')

%--------------------------------------------------------------------------

nam='Normalized model response for increasing Tau iso-Gains';
figure('name',nam);
subplot(2,1,1)
hold on;
%plot(OPvals.Table(1:end-0,2),OPvals.Table(1:end-0,5),'r.-')
plot(OPvals.TableC(1:end-clip,2),OPvals.TableC(1:end-clip,5),'r-')
ylabel('Kp')
title(nam)
%grid on;
axis([min(OPvals.TableC(:,2)) max(OPvals.TableC(:,2)) 0.15 max(OPvals.TableC(:,5))+0.1*max(OPvals.TableC(:,5))])
yl=ylim;
plot([Line1 Line1],yl,'b-')
plot([Line2 Line2],yl,'b-')
% plot([0.0620 0.0620],yl,'b-')
% plot([0.0500 0.0500],yl,'b-')

subplot(2,1,2)
hold on;
%plot(OPvals.Table(1:end-0,2),OPvals.Table(1:end-0,6),'k.-')
plot(OPvals.TableC(1:end-clip,2),OPvals.TableC(1:end-clip,6),'k-')
ylabel('Kd')
xlabel('Tau iso')
%grid on;
axis([min(OPvals.TableC(:,2)) max(OPvals.TableC(:,2)) 0 max(OPvals.TableC(:,6))+0.1*max(OPvals.TableC(:,6))])
yl=ylim;
plot([Line1 Line1],yl,'b-')
plot([Line2 Line2],yl,'b-')
% plot([0.0620 0.0620],yl,'b-')
% plot([0.0500 0.0500],yl,'b-')

%% fitting curve
%{
%run after saving dataset, in case Nans crash fit
% find and ignore Nans
Bstart=40;Bend=122;Cstart=130;Cend=191;
Vals_TisoB=OPvals.Table(Bstart:Bend,2);
Vals_TrespB=OPvals.Table(Bstart:Bend,7);
Vals_TisoC=OPvals.Table(Cstart:Cend,2);
Vals_TrespC=OPvals.Table(Cstart:Cend,7);

%[fitobject,gof,output]=fit(Vals_Tiso,Vals_Tresp,'poly1');
%[fitobject,gof,output]=fit(Vals_Tiso,Vals_Tresp,'poly2');
%[fitobject,gof,output]=fit(Vals_Tiso,Vals_Tresp,'poly3');
%[fitobject,gof,output] = fit(Vals_Tiso,Vals_Tresp,'exp1');
%[fitobject,gof,output] = fit(Vals_Tiso,Vals_Tresp,'power1');
[fitobjectB,gofB,outputB] = fit(Vals_TisoB,Vals_TrespB,'exp2');
[fitobjectC,gofC,outputC] = fit(Vals_TisoC,Vals_TrespC,'exp2');


close all;
nam='Fit to Tresp vs Tiso curve';
figure('name',nam);
hold on;
plot(fitobjectB,Vals_TisoB,Vals_TrespB)
plot(fitobjectC,Vals_TisoC,Vals_TrespC)
%plot(fitobjectD,Vals_TisoD,Vals_TrespD)
xlabel('Tiso');
ylabel('Tresp');
grid on;

clc
fitobjectB
gofB

fitobjectC
gofC

%}
%-------------------------------------------------------------------------
%{
%Dstart=40;Dend=191;% combine B &C regions

% Vals_TisoD=OPvals.Table(Dstart:Dend,2);
% Vals_TrespD=OPvals.Table(Dstart:Dend,7);

%[fitobjectD,gofD,outputD] = fit(Vals_TisoD,Vals_TrespD,'power2');

close all;
nam='Fit to Tresp vs Tiso curve';
figure('name',nam);
hold on;
plot(fitobjectD,Vals_TisoD,Vals_TrespD)
xlabel('Tiso');
ylabel('Tresp');
grid on;

clc
% fitobjectD
% gofD
%}

%% Saving data
%{
t=datetime;
coderuntime=toc;

notes={'ST norm model-Tiso fit';
       'slx: SwingTaskNorm';
        'Master code: Master_SwingTask_Norm';                                              
    'singlerun code:run_SwingTask_Norm';
       'maxOvershoot = 1e-6';
       'band = 0.02';
       'tnew=0:1e-1:Stoptime. This reduces output dataset size from 500 Mb to 3Mb, does not affect Tresp and %OS';
    ' ';
    'TisoVec=0.20:-0.001:0.001, below 0.009 is Nan initial  guess settling times';
    'Optimizer: none';
    'Fitting curve to normalized Tresp vs Tiso curve';
    'Solver Tols: ';
    ''};
save('TisofitST_v13.mat', '-v7.3');
%}
