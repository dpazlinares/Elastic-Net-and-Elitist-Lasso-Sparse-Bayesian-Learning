%% Pick Lead Fields for analysis
load ('E:\SSBL\Codes\vector_regression\LeadFields.mat');
subject    = [1,2]; %User defined, pick numbers from 1-9
LeadFields = LeadFields(1,subject);
%% Surface
load('ind-wMC0000005_mid_surface_total_81920-red-6000.mat');
load('faces-red-6000.mat');
%% Computing Simulation Dimensions
Nv      = size(vertices,1);    % Number of vertices
Nnoise  = 500;                 % Number of noisy vertices
Nsubj   = size(LeadFields,2);  % Number of Subjects
Nt      = 201;                 % Number of Time points
t       = 1:Nt;
%%
% The Simulated Sources Configurations are built by marking the coordenates in the
% 3d visualization using function Plot_sources_Haufe:
%     Plot_sources_Haufe(Jtest,vertices,faces,'simple')
% Later the coordinates are find in the vertices vector:
% Ej. coordinate = (19.53,72.34,32.39), then:
%     ans = find(abs(vertices(:,2)-19.53)<0.01)
%% Simulated Sources Locations at Cortical Surface Mesh
%%  Definition of Cortical points
OL = 1824; OR = 4001;
TL = 185;  TR = 5894;
PL = 628;  PR = 5248;
ML = 1077; MR = 4092;
FL = 2251; FR = 4234;
%% Generate Combinations Mask
comb = [0 0 0 0 0];
comb = cat(1,comb,unique(perms([1 0 0 0 0]),'rows'));
comb = cat(1,comb,unique(perms([1 1 0 0 0]),'rows'));
comb = cat(1,comb,unique(perms([1 1 1 0 0]),'rows'));
comb = cat(1,comb,unique(perms([1 1 1 1 0]),'rows'));
comb = cat(1,comb,[1 1 1 1 1]);
comb = comb';
%% Generate all possible interhemisferic combinations
for cont1 = 1:size(comb,2)
    point_sim{1,cont1} = [OL;TL;PL;ML;FL].*(1 - comb(:,cont1)) + [OR;TR;PR;MR;FR].*comb(:,cont1);
end
%% Geodesic Size of the Patches E0 CNEURO E-3 Eduardo
d = [20E0 10E0 10E0 10E0 15E0];
%% Amplitude of Activations
m = [1 1 0 1 0];
%% Pick configuration of sources for analysis
for cont1 = 1:size(comb,2)
    %% 1-Occipital Lobe 2-Temporal Lobe 3-Parietal Lobe 4-Motor Cortex 5-Frontal Lobe
    J = zeros(Nv,1);
    for cont2 = 1:length(m)
        Source           = point_sim{1,cont1}(cont2);
        [ind,findex]   = surfpatch(Source,vertices,faces,d(cont2));
        J(ind)         = 1;
    end
    figure
    Plot_sources_Haufe(J(:,1),vertices,faces,'simple');
    caxis([-1 1]);
end
config     = [1,32]; %User defined, pick numbers from 1-32 
close all
Nsim       = length(config);
V0_sim     = cell(1,Nsim);
Svv0_sim   = cell(1,Nsim);
V_sim      = cell(1,Nsim);
Svv_sim    = cell(1,Nsim);
J_sim      = cell(1,Nsim);
%% Simulating Time Series
G(1,:) = m(1)*cos((2*pi/75)*(t-1));
G(2,:) = m(2)*cos((2*pi/25)*(t-1));
G(3,:) = m(3)*cos((2*pi/50)*(t-1));
G(4,:) = m(4)*cos((2*pi/125)*(t-1));
G(5,:) = m(5)*cos((2*pi/225)*(t-1));
%% Pick sample number
Nsamp   = 3;