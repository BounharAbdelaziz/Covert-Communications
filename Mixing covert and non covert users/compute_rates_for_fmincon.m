function [r1, r2, rk] = compute_rates_for_fmincon(P_T, P_X2_mid_T, Epsilon_T, W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, Y_cardinality, X1_X2_cardinality, sk_budget, DEBUG_covert)
    
    % multiply with -1 for fmincon to the maximization
    rk = -1 * CovertCommunication.covert_sk_rate(P_T, P_X2_mid_T, Epsilon_T, W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, DEBUG_covert);    
    if (rk < sk_budget)
        r1 = -1 * CovertCommunication.covert_message_rate(P_T, P_X2_mid_T, Epsilon_T, W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, Y_cardinality, X1_X2_cardinality, DEBUG_covert);  
        r2 = -1 * CovertCommunication.conditional_MI(P_T, P_X2_mid_T, W_Y_X1_0_X2, X2_cardinality, Y_cardinality);
%         r2 = CovertCommunication.non_covert_rate(P_T, P_X2_mid_T, W_Y_X1_0, W_Y_X1_0_X2, X2_cardinality, DEBUG);
    else
        % we cannot achieve positive rates
        r1 = 0;
        r2 = 0;
        rk = 0;
    end
end