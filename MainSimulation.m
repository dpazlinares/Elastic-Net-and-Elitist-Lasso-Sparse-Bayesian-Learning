close all;
clear all;
%% Simulation initialization
Initialize_simulation
%% Run Simulation for specific configurations of sources and subjects
cont11 = 1;
for cont1 = config
    %% 1-Occipital Lobe 2-Temporal Lobe 3-Parietal Lobe 4-Motor Cortex 5-Frontal Lobe
    for cont2 = 1:length(m)
        Source                  = point_sim{1,cont1}(cont2);
        [index{cont2},findex]   = surfpatch(Source,vertices,faces,d(cont2));
        extensions(cont2,1)     = length(index{cont2});
    end
    index_sim{cont11} = index;
    %% Setting correlation structure between patches
    options.extensions      = extensions;
    options.connections     = [1 2;1 4;3 4;3 5;1 5];
    [S,Data,X,Index_cell]   = CSGGM(Nsamp,sum(extensions),options);
    pointer = 1;
    for cont2 = 1:length(m)
        M(:,cont2,:)            = zeros(Nv,1,Nsamp);
        M(index{cont2},cont2,:) = real(Data(pointer:pointer+extensions(cont2)-1,:));
        pointer                 = pointer + extensions(cont2); 
    end
    %% Simulating Sensors Time Series Data
    V0     = cell(1,Nsubj);
    for cont4 = 1:Nsubj
        K  = LeadFields{cont4};
        for cont3 = 1:Nsamp   
            V0{cont4}(:,:,cont3)  = K*M(:,:,cont3)*G;
        end
    end
    %% Computing Sources Variances
    J = 0;
    for cont3 = 1:Nsamp
    J  = J + abs(M(:,:,cont3)*G);
    end
    figure
    Plot_sources_Haufe(J(:,1),vertices,faces,'simple')
    caxis([-1 1])
    th = 0.001; J(abs(J)<th*(max(abs(J(:))))) = 0;
    %% Simulating Biological Noise
    rs0                = zeros(Nv,Nt*Nsamp*Nsubj);
    for cont4 = 1:Nsubj
    index_noise        = randperm(Nv);
    index_noise        = index_noise(1:Nnoise);
    rs0(index_noise,(cont4-1)*Nt*Nsamp+1:cont4*Nt*Nsamp) = restingstatenoise(200,10,2,Nt*Nsamp,Nnoise);
    end
    rs0                = reshape(rs0,Nv,Nt,Nsamp,Nsubj);
    %% Projecting Biological Noise to Sensors
    noisesources            = cell(1,Nsubj);
    for cont4 = 1:Nsubj
        K                   = LeadFields{cont4};
        Ne                  = size(K,1);
        rs_tmp              = zeros(Ne,Nt,Nsamp);
        for cont3 = 1:Nsamp
            rs_tmp(:,:,cont3) = K*squeeze(rs0(:,:,cont3,cont4));
        end
        noisesources{cont4} = sum(V0{cont4}(:).^2)^(1/2)*rs_tmp/sum(rs_tmp(:).^2)^(1/2);
    end
    %% Simulating Sensors Noise
    noisesensors            = cell(1,Nsubj);
    for cont4 = 1:Nsubj
        K                   = LeadFields{cont4};
        Ne                  = size(K,1);
        rs_tmp              = randn(Ne,Nt,Nsamp);
        noisesensors{cont4} = sum(V0{cont4}(:).^2)^(1/2)*rs_tmp/sum(rs_tmp(:).^2)^(1/2);
    end
    %% Simulating Data by K*J + SensorsNoise + BiologicalNoise
    V = cell(1,Nsubj);
    for cont4 = 1:Nsubj
    V{cont4} = V0{cont4} + 0.25*noisesources{cont4} + 0.25*noisesensors{cont4};
    end
    %% Computing Covariance across Data samples
    Svv  = cell(1,Nsubj);
    Svv0 = cell(1,Nsubj);
    for cont4 = 1:Nsubj
        for cont2 = 1:Nt
            Svv0{cont4}(:,:,cont2) = squeeze(V0{cont4}(:,cont2,:))*squeeze(V0{cont4}(:,cont2,:))'/Nsamp;
            Svv{cont4}(:,:,cont2)  = squeeze(V{cont4}(:,cont2,:))*squeeze(V{cont4}(:,cont2,:))'/Nsamp;
        end
    end
    %%
    J_sim{:,cont11}    = J;
    V0_sim{:,cont11}   = V0;
    V_sim{:,cont11}    = V;
    Svv0_sim{:,cont11} = Svv0;
    Svv_sim{:,cont11}  = Svv;
    cont11             = cont11 + 1;
end
save VisioMotor.mat G J_sim index_sim V0_sim Svv0_sim V_sim Svv_sim Nsamp vertices faces LeadFields