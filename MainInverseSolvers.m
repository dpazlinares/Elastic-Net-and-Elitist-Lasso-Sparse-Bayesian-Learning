%% Run Inverse Solvers
close all;
clear all;
%% Loading Simulation Substrate or Real Data
load('VisioMotor.mat')
%% Pick Inverse Solution Method
IS = 1; %User defined, pick numbers from 1-2
%% Pick Time/Frequency components
figure;plot(G')
components = [1,2]; %User defined, pick numbers from 1-201
%% Run Unidimensional Inverse Solvers 
if IS == 1
%% enet-ssbl
sol_enetssbl = InverseSolver_enet_ssbl(J_sim,V_sim,V0_sim,LeadFields,vertices,faces,components);
save Solutions_enetssbl.mat sol_enetssbl -v7.3
elseif IS == 2
%% elasso-ssbl
sol_elassossbl = InverseSolver_elasso_ssbl(J_sim,V_sim,V0_sim,LeadFields,vertices,faces,components);
save Solutions_elassossbl.mat sol_elassossbl -v7.3
end