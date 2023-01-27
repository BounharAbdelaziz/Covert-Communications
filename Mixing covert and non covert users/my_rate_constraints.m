function[c,ceq] = my_rate_constraints(probas_and_eps, sk_budget, W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, T_cardinality, X2_cardinality, optimization_P_X2, DEBUG_covert)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extracting probabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % probas_and_eps is in the form of [P_X2_mid_T(0 | 0), P_X2_mid_T(0 | 1), P_X2_mid_T(1 | 0), P_X2_mid_T(1 | 1), P_T(0), P_T(1), Epsilon_T(0), Epsilon_T(1)]
    % so we extract and construct P_T, P_X2_mid_T, Epsilon_T from it to
    % compute the secret key rate and compare it to the budget.
    
    P_X2_mid_T = zeros(X2_cardinality, T_cardinality);
    
    % Don't optimize P_X2_mid_T and set P_{X_2|T}(1|t)=1. (line 1)
    if (optimization_P_X2 == 1)        
        for t=1:T_cardinality
            P_X2_mid_T(2,t) = 1; % index starts with 1 so x2=1 is in index 2
        end
        index = 1; % to use it latter for P_T and Epsilon_T
    % Don't optimize P_X2_mid_T and set P_{X_2|T}(1|t)=0, i.e. set P_{X_2|T}(0|t)=1. (line 2)
    elseif (optimization_P_X2 == 2)
        for t=1:T_cardinality
            P_X2_mid_T(1,t) = 1; % index starts with 1 so x2=0 is in index 1
        end
        index = 1; % to use it latter for P_T and Epsilon_T
    % Optimize P_X2_mid_T
    elseif (optimization_P_X2 == 3) 
        index = 1;
        for x2=1:X2_cardinality
            for t=1:T_cardinality
                P_X2_mid_T(x2,t) = probas_and_eps(index);
                index = index + 1;
            end
        end
    end

    % P_T and Epsilon_T are vectors of size (T_cardinality, 1)
    P_T = zeros(T_cardinality, 1);
    for t=1:T_cardinality
        P_T(t) = probas_and_eps(index);
        index = index + 1;
    end

    Epsilon_T = zeros(T_cardinality, 1);
    for t=1:T_cardinality
        Epsilon_T(t) = probas_and_eps(index);
        index = index + 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Secret key budget constraint (rk \leq sk_buget) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rk = CovertCommunication.covert_sk_rate(P_T, P_X2_mid_T, Epsilon_T, W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, DEBUG_covert); 
    c_sk_budget = rk - sk_budget; % rk < sk_budget => rk - sk_budget < 0 (c <0 for fmincon)
    c = c_sk_budget;
%     c = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Probabilities constraint (should sum to 1)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (optimization_P_X2 == 3)
        % \sum_{x2} P_{X2 \mid T} (x2 \mid T=t) = 1
        for t=1:T_cardinality
            ceq(t) = sum(P_X2_mid_T(:,t)) - 1; 
        end
        % \sum_{t} P_{T} (T=t) = 1
        ceq(t+1) = sum(P_T(:,1)) - 1; 
    else
        % only \sum_{t} P_{T} (T=t) = 1 as P_X2_mid_T is not optimized in
        % the other cases
        ceq(1) = sum(P_T(:,1)) - 1;
    end    

%     disp('P_X2_mid_T')
%     disp(P_X2_mid_T)
%     disp('P_T')
%     disp(P_T)
%     disp('Epsilon_T')
%     disp(Epsilon_T)
%     disp('ceq')
%     disp(ceq)
%     disp('-------------------------------')
    
end
