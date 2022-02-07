%% Clear Workspace
close all;
clear;
clc

%% Load SBML Model
Model = sbmlimport('SBML_Model.xml');

%% Controller Parameters
% % Positive
% mu = 15;
% theta = 0.1;
% eta = 0.1;
% k = 1;
% delta_c = 0;
% Negative
mu = 90;
theta = 0.1;
eta = 0.0001;
k = 1;
delta_c = 0;

%% Add Controller Compartment
Controller = addcompartment(Model, 'Controller');
% Controller Species
addspecies(Controller, 'Z_1', 0);
addspecies(Controller, 'Z_2', 0);
% Reference Reaction
Reaction = addreaction(Model, 'null -> Controller.Z_1');
KineticLaw = addkineticlaw (Reaction, 'Unknown');
addparameter(KineticLaw, 'mu', mu);
Reaction.ReactionRate = 'mu';
% Sensing Reaction
Reaction = addreaction(Model, 'null -> Controller.Z_2');
KineticLaw = addkineticlaw (Reaction, 'Unknown');
addparameter(KineticLaw, 'theta', theta);
Reaction.ReactionRate = 'theta * Plasma.LDLC';
% Sequestration Reaction
Reaction = addreaction(Model, 'Controller.Z_1 + Controller.Z_2 -> null');
KineticLaw = addkineticlaw (Reaction, 'Unknown');
addparameter(KineticLaw, 'eta', eta);
Reaction.ReactionRate = 'eta * Controller.Z_1 * Controller.Z_2';
% % Positive Actuation Reaction
% Reaction = addreaction(Model, 'null -> Intestine.IC');
% KineticLaw = addkineticlaw (Reaction, 'Unknown');
% addparameter(KineticLaw, 'k', k);
% Reaction.ReactionRate = 'k * Controller.Z_1';
% Negative Actuation Reaction
Reaction = addreaction(Model, 'Controller.Z_2 + Intestine.IC -> Controller.Z_2');
KineticLaw = addkineticlaw (Reaction, 'Unknown');
addparameter(KineticLaw, 'k', k);
Reaction.ReactionRate = 'k * Controller.Z_2 * Intestine.IC';
% % Dilution Reactions
% Reaction = addreaction(Model, 'Controller.Z_1 -> null');
% KineticLaw = addkineticlaw (Reaction, 'Unknown');
% addparameter(KineticLaw, 'delta_c', delta_c);
% Reaction.ReactionRate = 'delta_c * Controller.Z_1';
% Reaction = addreaction(Model, 'Controller.Z_2 -> null');
% KineticLaw = addkineticlaw (Reaction, 'Unknown');
% addparameter(KineticLaw, 'delta_c', delta_c);
% Reaction.ReactionRate = 'delta_c * Controller.Z_2';

%% Get Equations
Equations = getequations(Model);

%% Simulate Model
csObj = getconfigset(Model,'active');
    csObj.SolverType = 'ode15s';
    csObj.Stoptime = 800;
SimData = sbiosimulate(Model);

%% Plot
Index = [2, 28, 35, 36];
figure();
Handle_Axis1 = subplot(2,2,1);
    hold(Handle_Axis1, 'on');
    Handle_Axis1.XLabel.String = 't';
    Handle_Axis1.YLabel.String = SimData.DataNames{Index(1)};
plot(Handle_Axis1, SimData.Time, SimData.Data(:,Index(1)));
Handle_Axis2 = subplot(2,2,2);
    hold(Handle_Axis2, 'on');
    Handle_Axis2.XLabel.String = 't';
    Handle_Axis2.YLabel.String = SimData.DataNames{Index(2)};
plot(Handle_Axis2, SimData.Time, SimData.Data(:,Index(2)));
Handle_Axis3 = subplot(2,2,3);
    hold(Handle_Axis3, 'on');
    Handle_Axis3.XLabel.String = 't';
    Handle_Axis3.YLabel.String = SimData.DataNames{Index(3)};
plot(Handle_Axis3, SimData.Time, SimData.Data(:,Index(3)));
Handle_Axis4 = subplot(2,2,4);
    hold(Handle_Axis4, 'on');
    Handle_Axis4.XLabel.String = 't';
    Handle_Axis4.YLabel.String = SimData.DataNames{Index(4)};
plot(Handle_Axis4, SimData.Time, SimData.Data(:,Index(4)));