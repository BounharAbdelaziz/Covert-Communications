function [r1, r2, rk] = compute_rates_for_fmincon(probas_and_eps, W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, T_cardinality, X2_cardinality, Y_cardinality, X1_X2_cardinality, sk_budget, DEBUG_covert)
    
    % For card(T)=card(X2)=2, probas_and_eps is in the form of [P_X2_mid_T(0 | 0), P_X2_mid_T(0 | 1), P_X2_mid_T(1 | 0), P_X2_mid_T(1 | 1), P_T(0), P_T(1), Epsilon_T(0), Epsilon_T(1)]
    % so we extract and construct P_T, P_X2_mid_T, Epsilon_T from it.
    P_X2_mid_T = zeros(X2_cardinality, T_cardinality);
    index = 1;
    for x2=1:X2_cardinality
        for t=1:T_cardinality
            P_X2_mid_T(x2,t) = probas_and_eps(index);
            index = index + 1;
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

%     P_X2_mid_T = zeros(2, 1);
%     % First line is for X2=0
%     P_X2_mid_T(1,1) = probas_and_eps(1);
%     % Second line is for X2= 1
%     P_X2_mid_T(2,1) = probas_and_eps(2);
%     
%     % P_T and Epsilon_T are vectors of size (T_cardinality, 1)
%     P_T = [probas_and_eps(3)];
%     Epsilon_T = [probas_and_eps(4)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Marginal distributions P_X1 and P_X2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    P_X2 = zeros(1,X2_cardinality);
    for x2=1:X2_cardinality
        avg = 0;
        for t=1:T_cardinality
            avg = avg + P_X2_mid_T(x2,t) * P_T(t);
        end
        P_X2(x2) = avg;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Marginal distributions W_Y_X1_0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    W_Y_X1_0 = zeros(1, Y_cardinality);
    for y=1:Y_cardinality
        W_y_X1_0 = 0;
        for x2=1:X2_cardinality
            W_y_X1_0 = W_y_X1_0 + P_X2(x2) * W_Y_X1_0_X2(x2,y);
        end
        W_Y_X1_0(y) = W_y_X1_0;
    end    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute rates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rk = CovertCommunication.covert_sk_rate(P_T, P_X2_mid_T, Epsilon_T, W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, DEBUG_covert);    
    if (rk <= sk_budget)
        r1 = CovertCommunication.covert_message_rate(P_T, P_X2_mid_T, Epsilon_T, W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, Y_cardinality, X1_X2_cardinality, DEBUG_covert);  
%         r2 = CovertCommunication.conditional_MI(P_T, P_X2_mid_T, W_Y_X1_0_X2, X2_cardinality, Y_cardinality);
        r2 = CovertCommunication.non_covert_rate(P_T, P_X2_mid_T, W_Y_X1_0, W_Y_X1_0_X2, X2_cardinality, DEBUG_covert);
    else
        r1 =0;
        r2= 0;
    end
end