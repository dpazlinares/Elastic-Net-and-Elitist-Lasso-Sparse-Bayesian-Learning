function [miu,alpha] = elasso_ssbl(V,K,vertices,faces)
% Elitist Lasso_Sparse Bayesian Learning
%% Loading Dimensions
Ne           = size(K,1);                                             % Number of electrodes
Nv           = size(K,2);                                             % Number of nodes
Nt           = size(V,2);                                             % Number of Time points or Frequency bins
Ine          = spdiags(ones(Ne,1),0,Ne,Ne);                           % Identity matrix with Ne size
sigma        = ones(Nv,Nt);                                           % Prior Variances
sigma_post   = ones(Nv,Nt);                                           % Posterior Variances
miu          = ones(Nv,Nt);                                           % Inverse Solution
miu_cal      = ones(Nv,Nt);                                           % Calibration Inverse Solution
%% User defined parameters
maxiter1     = 15;                                                    % Number of Iterations of outer cycle
maxiter11    = 15;                                                    % Number of Iterations of inner cycle
ealpha       = 1E4;                                                   % Rate parameter of the gamma pdf of alpha (higher -> smoother)
delta_0      = 1E-2;                                                  % Minimum value of delta's
beta         = 1;                                                     % Data variance
w            = ones(Nv,1);                                            % Vector of parameters' weights
W            = repmat(transpose(w),length(w),1)-diag(w);              % Parameters to Deltas' Transformation matrix
w_array      = repmat(w,1,Nt);                                        % Repmat transformation matrix along time
%% Initialization of parameters
Kt           = transpose(K); 
scale_K      = sqrt(trace(K*Kt)/Ne);                                  % Lead Field Scale Norm inf
K            = K/scale_K;                                             % Lead Field scaling
Kt           = Kt/scale_K;
KtK          = Kt*K;
%% Calibration Inverse Solution
for time = 1:Nt
    sigmaKt      = spdiags(sigma(:,time),0,Nv,Nv)*Kt;
    sigma_post1  = sigmaKt/(K*sigmaKt+beta*Ine);
    sigma_post2  = K*spdiags(sigma(:,time),0,Nv,Nv);
    for jj=1:Nv
        sigma_post(jj,time)  = sigma(jj,time)-sigma_post1(jj,:)*sigma_post2(:,jj);
    end
    miu_cal(:,time) = (1/beta).*(sigmaKt - sigma_post1*(sigma_post2*Kt))*V(:,time);
end
%% Data Scaling
scale_V      = sqrt(sum(abs(miu_cal(:)).^2)/(Nv*Nt))/sqrt(max(mean(sigma_post,2)));                               
V            = V./scale_V;                                      
%% Initialization of Hyperparameters
alpha        = 1E-8;                                                     % Elitist LASSO norm hyperparameter
%% Main Cycle
for cont1 = 1:maxiter1
    for cont11 = 1:maxiter11
        for time = 1:Nt
            %% Update Posterior Mean and Covariance matrix
            sigmaKt     = spdiags(sigma(:,time),0,Nv,Nv)*Kt;
            sigma_post1 = sigmaKt/(K*sigmaKt+beta*Ine);
            sigma_post2 = K*spdiags(sigma(:,time),0,Nv,Nv);
            % Only save the diagonals of the Posterior Covariance
        for jj=1:Nv
            sigma_post(jj,time)  = sigma(jj,time)-sigma_post1(jj,:)*sigma_post2(:,jj);
        end
        miu(:,time)      = (1/beta).*(sigmaKt - sigma_post1*(sigma_post2*Kt))*V(:,time);
        end
        %% Update Gammas
        index_h          = find((abs(miu).^2+sigma_post)<0);
        delta            = W*abs(miu)+delta_0.*ones(Nv,Nt);
        h                = sqrt((ones(Nv,Nt)./4).^2+alpha^2.*(w_array.^2).*(abs(miu).^2+sigma_post).*delta.^2)-ones(Nv,Nt)./4;
        h(index_h)       = 0;
        gamma            = alpha.*delta.^2+h;
        sigma            = 1/(2*alpha)-delta.^2./(2*gamma);
        sigma_bar        = 2*alpha*sigma;
    end
    %% Update alpha
    miu_aux             = miu;
    index_alpha         = find(sigma_bar>0.0005*max(max(sigma_bar)));
    mixed_norm_array    = abs(miu_aux);
    mixed_norm          = sum(sum(mixed_norm_array).^2);
    c1                  =(sum((abs(miu_aux(index_alpha)).^2+sigma_post(index_alpha))./(sigma_bar(index_alpha)))+sum(sum(delta.^2./(1-sigma_bar)))+mixed_norm)/(Nv*Nt);
    f_aux               = @(alpha_aux) c1+ealpha-1/alpha_aux-(1/2)*(1+length(index_alpha)/(Nv*Nt))*(1/alpha_aux)-(1/(Nv*Nt))*sum(sum((trascendent_term(alpha_aux*delta(:).^2))'.*(delta(:)).^2));
    alpha               = fzero(f_aux,[10^(-8) 1000]);
    hyp_convergence(cont1) = alpha;
end
miu = (miu*scale_V/scale_K);
end