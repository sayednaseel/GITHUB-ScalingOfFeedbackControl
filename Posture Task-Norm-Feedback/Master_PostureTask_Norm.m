%Master_PostureTask_Norm
% Code to evaluate the relationship between normalized torque limits
% (tau_iso) and response time (tresp) for the normalized posture task model.
% For each tau_iso, run_PostureTask_Norm optimizes controller gains to find
% the fastest settling time (on the angular velocity curve) without overshoot (on the angle curve). Tries fitting different
% curves to the tau_iso vs tresp relationship.

clear all;close all;clc
%%
mdl = 'PostureTaskNorm';

load_system(mdl)
solver_variable=2; % 1 for variable step, 2 for fixed step, modify settings within init_SwingTaskCL
init_Norm(mdl,solver_variable);% applies solver settings to mdl
hws = get_param(mdl,'modelworkspace');%handle to model workspace
hws.clear;
Stoptime=20;% simulation runtime
run_opt=0;% to run optimization
plotfig=0;% to plot figures within Run code.
optimizerMethod = 'fminsearch';
%optimizerMethod = 'fsolve';
%optimizerMethod = 'fmincon';
parms.StepRespM=2;% Method to find Tsettl: 1 sets bands based on Angle movement range. 2 sets bands on AngVel. 

%% Define parameters
parms.I        = 1;     % normalized
parms.t_delay  = 1;     % normalized
parms.theta_r  = 1;     % normalized
parms.theta_0  = 1;     % doesn't change result, only theta_r - theta_0 matters.
parms.dtheta_0 = -1;     % Posture task initial AngVel
%parms.tau_iso  = 1e12;  % ignore for now
parms.d        = 0;     % ignore for now
parms.g        = 0;     % ignore for now
parms.m        = 1;     % doesn't matter for g = 0
parms.l        = 1;     % doesn't matter for g = 0
parms.inverse  = +1;    % doesn't matter for g = 0
% Define borders: Exclude sims which exceed this limit
parms.maxOvershoot = 1e-6;%Same as solver reltol, see init code
parms.band = 0.02;% settling time band ()
parms.Kpmin=0.05;% below this, there is a region with bad sims with low settling times
parms.Kdmin=0.3;
parms.Kpmax=4*parms.I*(0.647/parms.t_delay)^2;
parms.Kdmax=4*sqrt(parms.Kpmax*parms.I);




%% Optimal gains:

load('TisofitPT_v3small','OPvals');
KpICVec=OPvals.Table(:,5);KdICVec=OPvals.Table(:,6);Ki=0;%2st column is fastest settling time, limiting overshoot to maxOvershoot
TisoVec_old=OPvals.Table(:,2);
clear OPvals;

%%
tnew=0:1e-4:Stoptime;
TisoVec=0.14:0.01:1.00;% vector of torque limits

tic
for i=1:length(TisoVec)
disp(['The Torque limit is ' num2str(TisoVec(i))]);
parms.tau_iso=TisoVec(i);

% Gains loaded from dataset
indTiso=find(abs(TisoVec_old-parms.tau_iso)<1e-10);
KpIC=KpICVec(indTiso);
KdIC=KdICVec(indTiso);


[OPtable1,OP1]=run_PostureTask_Norm(KpIC,KdIC,parms,Stoptime,hws,run_opt,mdl,plotfig,optimizerMethod);

OPvals.Table(i,:)=OPtable1;
OPvals.Time(i,:)=tnew;
OPvals.Angle(i,:)=interp1(OP1.Time,OP1.Angle1,tnew);
OPvals.AngVel(i,:)=interp1(OP1.Time,OP1.AngVel,tnew);
OPvals.Tsat(i,:)=interp1(OP1.Time,OP1.Tsat,tnew);
end 

% Note response time is settling time+SMdelay, since simulation starts when
% reference is 1, ignoring the initial delay period. this prevents the
% impulse due to differentiating a step input
AAOUT=table(OPvals.Table(:,1),OPvals.Table(:,2),OPvals.Table(:,3),OPvals.Table(:,4),OPvals.Table(:,5),OPvals.Table(:,6),OPvals.Table(:,7),OPvals.Table(:,8),OPvals.Table(:,9));
AAOUT.Properties.VariableNames={'SMdelay','Tiso' ,'band' ,'maxOvershoot', 'KpOPT', 'KdOPT', 'Tsettling' ,'Overshoot (%)' ,'Peak'};



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

Curve1=0.24;% highlight 1 curve
indTiso2=find(abs(OPvals.Table(:,2)-Curve1)<1e-10);
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
TisoAval=0.90;
TisoBval=0.50;
TisoCval=0.20;

indTisoA=find(abs(OPvals.Table(:,2)-TisoAval)<1e-10);
indTisoB=find(abs(OPvals.Table(:,2)-TisoBval)<1e-10);
indTisoC=find(abs(OPvals.Table(:,2)-TisoCval)<1e-10);

% %indM=indM;
plot(OPvals.Time(indTisoA,:),(OPvals.Tsat(indTisoA,:)),'k-','LineWidth',2)
plot(OPvals.Time(indTisoB,:),(OPvals.Tsat(indTisoB,:)),'r-','LineWidth',2)
plot(OPvals.Time(indTisoC,:),(OPvals.Tsat(indTisoC,:)),'b-','LineWidth',2)
grid on;
legend(['Tau iso:' num2str(TisoAval)],['Tau iso:' num2str(TisoBval)],['Tau iso:' num2str(TisoCval)])
%title(['TisofitPT_v3:' num2str(Curve1)])



% %Expanding dataset to show delay-limited range
% Data1=OPvals.Table(1,:);
% OPvals.TableB=repmat(Data1,20,1);
% OPvals.TableB(1:20,2)=[1.01:0.01:1.20];
% OPvals.TableC=[OPvals.TableB;OPvals.Table];
OPvals.TableC=OPvals.Table;
%--------------------------------------------------------------------------
%Tnosat=0.4;
Line1=0.80;% When Tiso limits clips top part of Torque curve
Line2=0.24;% When Tiso limits clip both top and bottom parts of torque curve.
clip=0;% TisofitPTv3 only goes to 0.14 to avoid Nans. 
nam='Normalized model response for increasing Tau iso-Step response';
figure('name',nam);
subplot(2,1,1)
hold on;
plot(OPvals.TableC(1:end-clip,2),OPvals.TableC(1:end-clip,7),'r-')
%yl=ylim;
%plot([Tnosat Tnosat],yl)
ylabel('Settling time')
%legend()
title(nam)
%grid on;
%axis([min(OPvals.Table(:,1)) max(OPvals.Table(:,1)) 0 max(OPvals.Table(:,2))+0.1*max(OPvals.Table(:,2))])
yl=ylim;
plot([Line1 Line1],yl,'b-')
plot([Line2 Line2],yl,'b-')


subplot(2,1,2)
hold on;
plot(OPvals.TableC(1:end-clip,2),OPvals.TableC(1:end-clip,8),'k-')
ylabel('Overshoot')
%grid on;
%axis([0 max(OPvals.TableC(:,2)) -0.5 1])
%yl=ylim;
%plot([Tnosat Tnosat],yl)
yl=ylim;
plot([Line1 Line1],yl,'b-')
plot([Line2 Line2],yl,'b-')

 xlabel('Tau iso')


%--------------------------------------------------------------------------

nam='Normalized model response for increasing Tau iso-Gains';
figure('name',nam);
subplot(2,1,1)
hold on;
plot(OPvals.TableC(1:end-clip,2),OPvals.TableC(1:end-clip,5),'r-')
ylabel('Kp')
title(nam)
%grid on;
axis([0 max(OPvals.TableC(:,2)) 0 max(OPvals.TableC(:,5))+0.1*max(OPvals.TableC(:,5))])
yl=ylim;
plot([Line1 Line1],yl,'b-')
plot([Line2 Line2],yl,'b-')


subplot(2,1,2)
hold on;
plot(OPvals.TableC(1:end-clip,2),OPvals.TableC(1:end-clip,6),'k-')
ylabel('Kd')
xlabel('Tau iso')
%grid on;
axis([0 max(OPvals.TableC(:,2)) 0 max(OPvals.TableC(:,6))+0.1*max(OPvals.TableC(:,6))])
yl=ylim;
plot([Line1 Line1],yl,'b-')
plot([Line2 Line2],yl,'b-')

%% Saving data
%{
t=datetime;
coderuntime=toc;

notes={'PT norm model-Tiso fit';
       'slx: PostureTaskNorm';
       'maxOvershoot = 1e-6';
       'band = 0.02';
    ' ';
    'TisoVec=0.14:0.01:1.00';
    'Master code: Master_';                                              
    'singlerun code:run_';
    'Optimizer: fminsearch';
    'Solver tols: fixed timestep 1e-4';
    'parms.Kpmax=4*parms.I*(0.647/parms.t_delay)^2';
    '';
    '';
    '';
    ''};
save('TisofitPT_v3.mat', '-v7.3');
%}
%% fitting curve
%{
%run after saving dataset, in case Nans crash fit
%TisofitPT_v3:'TisoVec=1:-0.01:0.14. fit indices 22 to end
% find and ignore Nans
Vals_Tiso=OPvals.Table(22:end,2);
Vals_Tresp=OPvals.Table(22:end,7);

%[fitobject,gof,output]=fit(Vals_Tiso,Vals_Tresp,'poly1');
%[fitobject,gof,output]=fit(Vals_Tiso,Vals_Tresp,'poly2');
%[fitobject,gof,output]=fit(Vals_Tiso,Vals_Tresp,'poly3');
%[fitobject,gof,output] = fit(Vals_Tiso,Vals_Tresp,'exp1');
%[fitobject,gof,output] = fit(Vals_Tiso,Vals_Tresp,'power1');
[fitobject,gof,output] = fit(Vals_Tiso,Vals_Tresp,'power2');
%[fitobject,gof,output] = fit(Vals_Tiso,Vals_Tresp,'exp2');


close all;
nam='Fit to Tresp vs Tiso curve';
figure('name',nam);
plot(fitobject,Vals_Tiso,Vals_Tresp)
xlabel('Tiso');
ylabel('Tresp');
grid on;

clc
fitobject
gof
%}
%% checking predictions
%{
DAT=load('','OPvals');%
DAT2=load('','NormVals');%

iMass=DAT.OPvals.Table(1,:);
iMOI=DAT.OPvals.Table(2,:);
iTiso=DAT.OPvals.Table(3,:);
iTiso_norm=DAT2.NormVals.Table(2,:);
iTiso_norm2=iTiso_norm;
TmaxInf=OPvals.Table(1,5);%the maximum Torque when no torque limits
iTiso_norm2(iTiso_norm>=TmaxInf)=TmaxInf;
iSMdelay=DAT.OPvals.Table(4,:)./1000;
iThetaR=deg2rad(30);
iTresp_gold=DAT.OPvals.Table(9,:)./1000;
clear DAT DAT2

% Term1=(iTiso.^f.b.*iSMdelay.^(2.*f.b).*iMOI.^(-f.b)*iThetaR.^(-f.b));
% iTresp_fit(2,:)=iSMdelay.*(f.a.*Term1+f.c);

iTresp_fit(2,:)=iSMdelay.*((f.a.*iTiso_norm2.^f.b)+f.c);


iTresp_fit(1,:)=iMass;
iTresp_fit(3,:)=iTresp_gold;
iTresp_fit(4,:)=((iTresp_fit(2,:)-iTresp_fit(3,:))./iTresp_fit(3,:))*100;

AB.Tablehead={'Mass (Kg)';'Predicted Tresp (s)';'Actual Tresp (s)';'Error (%)'};



AB.Table=iTresp_fit;
ABOUT=struct2table(AB);
%}

