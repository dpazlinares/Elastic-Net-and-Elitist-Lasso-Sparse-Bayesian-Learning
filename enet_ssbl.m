function [miu,alpha1,alpha2]=enet_ssbl(V,K,vertices,faces)
% Elastic Net_Sparse Bayesian Learning
%% Loading Dimensions
Ne           = size(K,1);                                             % Number of electrodes
Nv           = size(K,2);                                             % Number of nodes
Ine          = spdiags(ones(Ne,1),0,Ne,Ne);                           % Identity matrix with Ne size
sigma        = ones(Nv,1);                                            % Prior Variances
sigma_post   = ones(Nv,1);                                            % Posterior Variances
%% User defined parameters
maxiter1     = 15;                                                    % Number of Iterations of outer cycle
maxiter11    = 15;                                                    % Number of Iterations of inner cycle
ealpha       = 1E1;                                                   % Rate parameter of the gamma pdf of alpha (higher -> smoother)
ek           = 1E1;                                                   % Rate parameter of the gamma pdf of k (higher -> smoother)
beta         = 1;                                                     % Data variance
%% Initialization of parameters
Kt           = transpose(K); 
scale_K      = sqrt(trace(K*Kt)/Ne);                                  % Lead Field Scale Norm inf
K            = K/scale_K;                                             % Lead Field scaling
Kt           = Kt/scale_K;
KtK          = Kt*K;
%% Calibration Inverse Solution
sigmaKt      = spdiags(sigma,0,Nv,Nv)*Kt;
sigma_post1  = sigmaKt/(K*sigmaKt+beta*Ine);
sigma_post2  = K*spdiags(sigma,0,Nv,Nv);  
for jj=1:Nv
    sigma_post(jj)  = sigma(jj)-sigma_post1(jj,:)*sigma_post2(:,jj);
end
miu_cal      = (1/beta).*(sigmaKt - sigma_post1*(sigma_post2*Kt))*V;
%% Data Scaling
scale_V      = sqrt(sum(abs(miu_cal).^2)/Nv)/sqrt(max(sigma_post));                               
V            = V./scale_V;                                      
%% Initialization of Hyperparameters
alpha1       = 1E0;                                                     % Hyperparameter of the L2 norm
k            = 1E0;                                                     % Hyperparameter of the Truncate Gamma pdf
%% Main Cycle
for cont1 = 1:maxiter1
    for cont11 =1:maxiter11
        %% Update Posterior Mean and Covariance matrix
        sigmaKt     = spdiags(sigma,0,Nv,Nv)*Kt;
        sigma_post1 = sigmaKt/(K*sigmaKt+beta*Ine);
        sigma_post2 = K*spdiags(sigma,0,Nv,Nv);
        % Only save the diagonals of the Posterior Covariance
        for jj=1:Nv
            sigma_post(jj)  = sigma(jj)-sigma_post1(jj,:)*sigma_post2(:,jj);
        end
        miu         = (1/beta).*(sigmaKt - sigma_post1*(sigma_post2*Kt))*V;
        %% Update Gammas
        index_h     = find((abs(miu).^2 + sigma_post)<0);
        h           = sqrt((ones(Nv,1)./4).^2+alpha1.*k.*(abs(miu).^2+sigma_post))-ones(Nv,1)./4;
        h(index_h)  = 0;
        gamma       = k + h;
        sigma_bar   = h./gamma;
        sigma       = (1/(2*alpha1))*sigma_bar;
    end
    %% Update alpha_1
    index_alpha = find(sigma_bar>0);
    alpha1      = (length(index_alpha)/2 + Nv)/(sum((abs(miu(index_alpha)).^2 + sigma_post(index_alpha))./(sigma_bar(index_alpha))) + ealpha*Nv);
    %% Update alpha_2
    f_aux       = @(k_aux) ek + sum(ones(Nv,1)./(1-sigma_bar))/Nv - (1/2)/k_aux - trascendent_term(k_aux);
    k           = fzero(f_aux,[0.000001 70000]);
    alpha2      = (4*alpha1*k)^(1/2);
    sigma       = (1/(2*alpha1))*sigma_bar;
    hyp_convergence(cont1,:) = [alpha1 k];
end
miu = (miu*scale_V/scale_K);
end