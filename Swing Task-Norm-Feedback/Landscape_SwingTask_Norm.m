%Landscape_SwingTask_Norm
%Code to simulate normalized Swing Task model for input controller gains. To generate landscape of
%overshoot and settling time for brute force search through combinations of
%Kp and Kd. 

clear all;close all;clc;
%%
%{1
mdl = 'SwingTaskNorm';

load_system(mdl)
solver_variable=1; % 1 for variable step, 2 for fixed step, modify settings within init_SwingTaskCL
init_Norm(mdl,solver_variable);% applies solver settings to mdl
hws = get_param(mdl,'modelworkspace');%handle to model workspace
hws.clear;
run_Elandscape=1;% to run cost landscape mapping
plotfig=1;% to plot figures
%}
%% Define parameters
Stoptime=20;% simulation runtime
I        = 1;     % normalized
theta_r  = 1;     % normalized
t_delay  = 1;     % normalized
theta_0  = 0;     % doesn't change result, only theta_r - theta_0 matters.
dtheta_0 = 0;     % ignore for now
%dtheta_0 = -0.21*sqrt(9.8066*1);% Initial angular velocity of 0.21 froude number
tau_iso  = 1e12; 


d        = 0;     % ignore for now
g        = 0;     % ignore for now
m        = 1;     % doesn't matter for g = 0
l        = 1;     % doesn't matter for g = 0
inverse  = +1;    % doesn't matter for g = 0
% Define borders: Exclude sims which exceed this limit
maxOvershoot = 1e-6;% Same as solver reltol, see init code
band = 0.02;% settling time band ()


% % For comparision, the paper gives:
% Kp = 0.647^2*I./(t_delay^2);
% Ki = 0;
% Kd = 2*sqrt(I*Kp);

%% Optimal gains:
load('optGains_SwingTask_Norm.mat','AA','AATable');
Kp=AA.vals(1,1);Kd=AA.vals(3,1);Ki=0;
clear AA AATable;
Kdcrit=2*sqrt(Kp*I);
Kdrat=Kd/Kdcrit;
%% Simulate once

sim(mdl)


%Calculate settling time David method
        %isValid = ~any(abs(theta.Data)>theta_r*(1+maxOvershoot));
        %isValid = ~any(abs(theta.Data)>theta_r*(1+maxOvershoot));

        % Compute settle time (within defined band):
        % Compute the last index in which the signal was above the band
        indexIsAbove = find(theta.Data>theta_r*(1+band),1,'last');
        if isempty(indexIsAbove)
            indexIsAbove = 1;
        end
        % Compute the last index in which the signal was below the band
        indexIsBelow = find(theta.Data<theta_r*(1-band),1,'last');
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
        if Peakval>theta_r
            Overshootval=((Peakval-theta_r)/theta_r)*100;
        else
            Overshootval=0;
        end % finding overshoot
        
 stepStats= stepinfo(theta.Data,theta.time,theta_r,'SettlingTimeThreshold',band);% considering the entire angle curve

%--------------------------------------------------------------------------
 if plotfig==1
     
     nam=['Step response-Swing task'];
     figure('name',nam);
     %clf
     hold on
     plot(theta.Time,theta.Data,'LineWidth',2);
     plot([0,Stoptime],[theta_r,theta_r])
     %plot([0,20],[theta_r*(1+maxOvershoot),theta_r*(1+maxOvershoot)])
     plot([0,Stoptime],[theta_r*(1+band),theta_r*(1+band)])
     plot([0,Stoptime],[theta_r*(1-band),theta_r*(1-band)])
     plot([tSettleD tSettleD],ylim,'k.-','LineWidth',2)
     grid on;
     legend('Step response','Reference','Settling band +','Settling band -','Settling time David')
     xlabel('Time (td)')
     ylabel('Angle (ND)')

%--------------------------------------------------------------------------
     nam=['Torque,AngVel,Ang of Swing task'];
     figure('name',nam);
     
     subplot(311); hold on
     hold on;
     plot(Tsat.Time,Tsat.Data,'LineWidth',2);
     ylabel('Torque (ND)');
     yl=ylim;
     plot([t_delay t_delay],yl,'c-')
     plot([2*t_delay 2*t_delay],yl,'m-')
     plot([tSettleD tSettleD],yl)
     %grid on;
     title(nam)

     subplot(312); hold on
     plot(thetadot.Time,thetadot.Data,'LineWidth',2);
     %xlabel('Time (Td)');
     ylabel('AngVel (ND)');
     %grid on;
     yl=ylim;
     plot([t_delay t_delay],yl,'c-')
     plot([2*t_delay 2*t_delay],yl,'m-')
     plot([tSettleD tSettleD],yl)
     
     subplot(313); hold on
     plot(theta.Time,theta.Data,'LineWidth',2);
     plot([0,Stoptime],[theta_r*(1+band),theta_r*(1+band)],'k--')
     plot([0,Stoptime],[theta_r*(1-band),theta_r*(1-band)],'k--')
     xlabel('Time (Td)');
     ylabel('Angle (ND)');
    % grid on;
     yl=ylim;
     plot([t_delay t_delay],yl,'c-')
     plot([2*t_delay 2*t_delay],yl,'m-')
     plot([tSettleD tSettleD],yl)
     legend('Angle','','','Tdelay','2Tdelay','Settling time')
     axis([0 20 0 1.2])
 end % plotting


AA.Tablehead={'Kp';'Ki';'Kd';'Settling Time-David';'Overshoot-David';'Settling Time-stepinfo';'Overshoot-stepinfo'};
AA.vals=[Kp;Ki;Kd;tSettleD;Overshootval;stepStats.SettlingTime;stepStats.Overshoot];
AATable=struct2table(AA);


%% Run a parameter study:
tic
if run_Elandscape==1
%{1
% Mesh gains:
n = 300;
% Overview for Ki=0;
%[KPs,KDs] = meshgrid(linspace(0,0.3,n),linspace(0.3,1.3,n)); 
[KPs,KDs] = meshgrid(linspace(0,0.5,n),linspace(0.3,1.5,n)); 

% Zoom in for Ki = 0:
%[KPs,KDs] = meshgrid(linspace(0.17,0.19,n),linspace(0.64,0.66,n));
OP_TSettle = zeros(size(KPs)); 
OP_Peak=zeros(size(KPs));
OP_Overshoot=zeros(size(KPs));
tic
for i = 1:size(KPs,1)
    disp(i);
    for j = 1:size(KPs,2)
        Kp = KPs(i,j);
        Kd = KDs(i,j);
        % run simulation
        sim(mdl)
        % isValid = ~any(abs(theta.Data)>theta_r*(1+maxOvershoot));
        % Compute settle time (within 2% of target):
        % Compute the last index in which the signal was above the band
        indexIsAbove = find(theta.Data>theta_r*(1+band),1,'last');
        if isempty(indexIsAbove)
            indexIsAbove = 1;
        end
        % Compute the last index in which the signal was below the band
        indexIsBelow = find(theta.Data<theta_r*(1-band),1,'last');
        if isempty(indexIsBelow)
            indexIsBelow = 1;
        end
        % Compute the last index in which the signal was outside the band
        indexIsOutside = max(indexIsAbove,indexIsBelow);
        % Check if the trial is valid and this index is not the end of the time
        % series (which also implies an invalid trial)
        if indexIsOutside<size(theta.Data,1)% && isValid
            tSettle = theta.time(indexIsOutside+1);
        else
            tSettle = NaN;
        end
             
        OP_TSettle(i,j) = tSettle;
        
        % finding overshoot
        Peakval=max(theta.Data);
        OP_Peak(i,j)=Peakval;
        if Peakval>theta_r
            Overshootval=((Peakval-theta_r)/theta_r)*100;
        else
            Overshootval=0;
        end % finding overshoot
        
        OP_Overshoot(i,j)=Overshootval;

    end% KDs
end % KPs
runtime=toc;
%}

%% Processing info

ElandMin=min(min(OP_TSettle));
[Ind_Kpmin,Ind_Kdmin]=find(OP_TSettle==ElandMin);


% Finding sims which stay within overshoot limits
OP_OvershootLIM=OP_Overshoot>maxOvershoot;
OP_TSettleLIM=OP_TSettle;
OP_TSettleLIM(OP_OvershootLIM)=nan;% taking only sims which stay below overshoot limit

ElandMinLIM=min(min(OP_TSettleLIM));
[Ind_KpminLIM,Ind_KdminLIM]=find(OP_TSettleLIM==ElandMinLIM);


AC.Tablehead={'Kp_min';'Kd_min';'Settling Time-David';'Overshoot-David'};
AC.vals(:,1)=[KPs(Ind_Kpmin,Ind_Kdmin);KDs(Ind_Kpmin,Ind_Kdmin);OP_TSettle(Ind_Kpmin,Ind_Kdmin);OP_Overshoot(Ind_Kpmin,Ind_Kdmin)];
AC.vals(:,2)=[KPs(Ind_KpminLIM,Ind_KdminLIM);KDs(Ind_KpminLIM,Ind_KdminLIM);OP_TSettle(Ind_KpminLIM,Ind_KdminLIM);OP_Overshoot(Ind_KpminLIM,Ind_KdminLIM)];
ACTable=struct2table(AC);

%% plotting settling time and overshoot brute force search landscapes
close all;

if plotfig==1
    figure('name','Settling Time')
    hold on;
    surf(KPs,KDs,OP_TSettle,'EdgeColor','flat')
    plot3(KPs(Ind_Kpmin,Ind_Kdmin),KDs(Ind_Kpmin,Ind_Kdmin),OP_TSettle(Ind_Kpmin,Ind_Kdmin),'ko','MarkerSize',10,'MarkerSize',10,'MarkerFaceColor','k')
    plot3(KPs(Ind_KpminLIM,Ind_KdminLIM),KDs(Ind_KpminLIM,Ind_KdminLIM),OP_TSettle(Ind_KpminLIM,Ind_KdminLIM),'ko','MarkerSize',10,'MarkerSize',10,'MarkerFaceColor','r')
    xlabel('\it K_p','Fontname','Times New Roman','Fontsize',25)
    ylabel('\it K_d','Fontname','Times New Roman','Fontsize',25)
    zlabel('\it t_{settling}','Fontname','Times New Roman','Fontsize',25)
    %legend('','min(Tsettle),with OS','min(Tsettle), constrained OS')
    %title(['Torque limit:' num2str(tau_iso) ])
    
    figure('name','Peak')
    hold on;
    surf(KPs,KDs,OP_Peak,'EdgeColor','flat')
    plot3(KPs(Ind_Kpmin,Ind_Kdmin),KDs(Ind_Kpmin,Ind_Kdmin),OP_Peak(Ind_Kpmin,Ind_Kdmin),'ko','MarkerSize',10,'MarkerSize',10,'MarkerFaceColor','k')
    plot3(KPs(Ind_KpminLIM,Ind_KdminLIM),KDs(Ind_KpminLIM,Ind_KdminLIM),OP_Peak(Ind_KpminLIM,Ind_KdminLIM),'ko','MarkerSize',10,'MarkerSize',10,'MarkerFaceColor','r')
    xlabel('Kp')
    ylabel('Kd')
    zlabel('Peak')
    title(['Torque limit:' num2str(tau_iso) ])

    
    figure('name','Overshoot')
    hold on;
    surf(KPs,KDs,OP_Overshoot,'EdgeColor','flat')
    plot3(KPs(Ind_Kpmin,Ind_Kdmin),KDs(Ind_Kpmin,Ind_Kdmin),OP_Overshoot(Ind_Kpmin,Ind_Kdmin),'ko','MarkerSize',10,'MarkerSize',10,'MarkerFaceColor','k')
    plot3(KPs(Ind_KpminLIM,Ind_KdminLIM),KDs(Ind_KpminLIM,Ind_KdminLIM),OP_Overshoot(Ind_KpminLIM,Ind_KdminLIM),'ko','MarkerSize',10,'MarkerSize',10,'MarkerFaceColor','r')
    xlabel('\it K_p','Fontname','Times New Roman','Fontsize',25)
    ylabel('\it K_d','Fontname','Times New Roman','Fontsize',25)
    zlabel('\it Overshoot','Fontname','Times New Roman','Fontsize',25)
    %title(['Torque limit:' num2str(tau_iso) ])

end
%% Saving data
%{
t=datetime;
notes={'code: Landscape_SwingTask_Norm';
       'slx: SwingTaskNorm';
       'maxOvershoot = 1e-6';
       'band = 0.02, Tiso=Inf';
       'n=200';
       '[KPs,KDs] = meshgrid(linspace(0,0.5,n),linspace(0.3,1.5,n))';
       '';
       ''};
save('Landscape_SwingTask_Norm.mat');
%}
%%
end % Elandscape
