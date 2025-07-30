function []=init_Norm(mdl,solver_variable)

% Code to set solver parameters and settings of the Simulink model
% Take care of the error messages 

%%
load_system(mdl)


% % here's a way to get solver setting names.
% configSet = getActiveConfigSet(mdl);
% configSetNames = get_param(configSet, 'ObjectParameters');


if solver_variable==1
%     TolRel='1e-3';TolAbs='1e-6';
 TolRel='1e-6';TolAbs='1e-9';% consistent with dde23 mass codes
%  TolRel='1e-18';TolAbs='1e-18';


    set_param(mdl,'SolverType','Variable-step')
    %set_param(mdl,'Solver','ode15s')
    set_param(mdl,'Solver','ode45')

    set_param(mdl,'RelTol',TolRel)
    set_param(mdl,'AbsTol',TolAbs)
    set_param(mdl,'MaxStep','1e-3')
    set_param(mdl,'MaxConsecutiveMinStep','10000');% To prevent sims stuck in optimization, limiting min step and consecutive min steps.
    set_param(mdl,'MinStep','1e-21');%

elseif solver_variable==2
    % fixed timestep solver
    
    set_param(mdl,'SolverType','Fixed-step')
    set_param(mdl,'Solver','ode4')
%    set_param(mdl,'FixedStep','1e-3')
    set_param(mdl,'FixedStep','1e-4')
    
end

set_param(mdl,'ZeroCrossControl','UseLocalSettings');%'UseLocalSettings' | 'EnableAll' | 'DisableAll'
set_param(mdl,'ZeroCrossAlgorithm','Adaptive');% 'Nonadaptive' | 'Adaptive'
set_param(mdl,'IgnoredZcDiagnostic','none');% When using adaptive
set_param(mdl,'DiscreteInheritContinuousMsg','none');%timestep counter warning off
set_param(mdl,'MaskedZcDiagnostic','none');% Masked zero crossing off
set_param(mdl,'ReturnWorkspaceOutputs','off');
set_param(mdl,'SimulationMode','normal')%normal, rapid, accelerator
set_param(mdl,'AccelVerboseBuild','on'); % see build process
set_param(mdl,'AlgebraicLoopMsg','warning'); % Setting algebraic loop warning off
set_param(mdl,'AlgebraicLoopSolver','LineSearch');%'TrustRegion'   'LineSearch'

set_param(mdl,'LoadInitialState','off');% If set on, run code cannot get xInitial
save_system(mdl)
