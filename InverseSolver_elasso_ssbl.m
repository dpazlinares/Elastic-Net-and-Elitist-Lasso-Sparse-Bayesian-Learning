function [sol_elassossbl] = InverseSolver_elasso_ssbl(J_sim,V_sim,V0_sim,LeadFields,vertices,faces,components)
Nsim         = size(V_sim,2);
sol_elassossbl = cell(2,Nsim);
for cont1 = 1:Nsim
    J       = J_sim{cont1};
    V_subj  = V_sim{cont1};
    V0_subj = V0_sim{cont1};
    Nsubj   = size(V_subj,2);
    for cont4 = 1:Nsubj
        V      = V_subj{cont4}; 
        V0     = V0_subj{cont4};
        Nt     = size(V,2);
        Nsamp  = size(V,3);
        cont22 = 1;
        for cont2 = components
            for cont3 = 1:Nsamp
                [iv,alpha] = elasso_ssbl(squeeze(V(:,cont2,cont3)),squeeze(LeadFields{cont4}),vertices,faces);
                sol(:,cont22,cont3,cont4)       = iv;
                alpha_est(cont22,cont3,cont4)  = alpha;
            end
            cont22 = cont22 + 1;
        end
    end
    %% Pick Time/Frequency Component and Subject to Plot 
    comp = 1;
    subj = 2;
    test = mean(sol.^2,3).^(1/2);
    figure
    subplot(1,2,1)
    Plot_sources_Haufe(J(:,comp)/max(abs(J(:,comp))),vertices,faces,'simple')
    caxis([-1 1])
    subplot(1,2,2)
    Plot_sources_Haufe(squeeze(test(:,comp,1,subj))/max(abs(squeeze(test(:,comp,1,subj)))),vertices,faces,'simple')
    caxis([-1 1])
    %%
    sol_elassossbl{1,cont1} = sol;
    sol_elassossbl{2,cont1} = alpha_est;
end