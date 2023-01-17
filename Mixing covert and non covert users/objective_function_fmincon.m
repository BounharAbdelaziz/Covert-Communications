function [loss] = objective_function_fmincon(probas_and_eps, W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, T_cardinality, X2_cardinality, Y_cardinality, X1_X2_cardinality, mu_1, mu_2, DEBUG_covert)

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

    % multiply with -1 for fmincon to the maximization and add mu_1 and
    % mu_2 to allow getting all points in the boudary.

    r1 = -1 * mu_1 * CovertCommunication.covert_message_rate(P_T, P_X2_mid_T, Epsilon_T, W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, Y_cardinality, X1_X2_cardinality, DEBUG_covert);  
    r2 = -1 * mu_2 * CovertCommunication.conditional_MI(P_T, P_X2_mid_T, W_Y_X1_0_X2, X2_cardinality, Y_cardinality);

    loss = r1 + r2;

end