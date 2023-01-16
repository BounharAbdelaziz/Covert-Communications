function covert_point = compute_covert_point(P_T, P_X2_mid_T, Epsilon_T, W_Z_X1_X2, X2_cardinality, Y_cardinality, DEBUG)
    % Input:
    %   P_T             : a probabilities vector for the time sharing rv
    %   P_X2_given_T    : conditional probability vector P_{X_2 \mid T}
    %   epsilons        : vector of \epsilon_{t} values 
    %   W_Y_X1_X2       : channel law of the legitimate user W_{Y \mid X_1, X_2}
    %   W_Z_X1_X2       : channel law of the adversary W_{Z \mid X_1, X_2}
    % Output:
    %   The secret key rate

    % We first compute the relative entropy for all X_2 for the legitimate user --------------------------------------------------------
    relative_entropy_vect   = zeros(X2_cardinality,1);
    
    % Eve's channel
    W_Z_X1_0_X2 = zeros(X2_cardinality,Y_cardinality);
    W_Z_X1_1_X2 = zeros(X2_cardinality,Y_cardinality);
    
    W_Z_X1_0_X2(1,:) = W_Z_X1_X2(1,:); % W_Z_X1_0_X2_0
    W_Z_X1_0_X2(2,:) = W_Z_X1_X2(2,:); % W_Z_X1_0_X2_1
    W_Z_X1_1_X2(1,:) = W_Z_X1_X2(3,:); % W_Z_X1_1_X2_0
    W_Z_X1_1_X2(2,:) = W_Z_X1_X2(4,:); % W_Z_X1_1_X2_1

    for x2=1:X2_cardinality                
        % Eve
        W_Z_X1_0_X2_x2 = W_Z_X1_0_X2(x2,:);
        W_Z_X1_1_X2_x2 = W_Z_X1_1_X2(x2,:);
        relative_entropy_eve = InformationTheory.relative_entropy(W_Z_X1_1_X2_x2, W_Z_X1_0_X2_x2);
        relative_entropy_vect(x2) = relative_entropy_eve;
    end

    % Then, we compute the chi_2 distance for all X_2 for Eve's channel --------------------------------------------------------
    chi_2_vect = zeros(X2_cardinality,1);
    
    for x2=1:X2_cardinality
        W_Z_X1_0_X2_x2 = W_Z_X1_0_X2(x2,:);
        W_Z_X1_1_X2_x2 = W_Z_X1_1_X2(x2,:);
        chi_2 = InformationTheory.chi_2_distance(W_Z_X1_1_X2_x2, W_Z_X1_0_X2_x2);
        chi_2_vect(x2) = chi_2;
    end
    if DEBUG
        disp('chi_2_vect is: ');
        disp(chi_2_vect)
    end
    % Finaly, we compute the average of the nominator and denominator w.r.t P_T and P_{X_2 \mid T} --------------------------------------------------------
    avg_nominator = 0;
    avg_denominator = 0;
    for t=1:length(P_T)
        % computing the averages of the relative entropy and chi_2 for a fixed t
        average_relative_entropy_fixed_t = dot(P_X2_mid_T(:,t), relative_entropy_vect);
        average_chi_2_fixed_t = dot(P_X2_mid_T(:,t), chi_2_vect);
        
        % taking the average over T
        avg_nominator = avg_nominator + P_T(t)*Epsilon_T(t)*average_relative_entropy_fixed_t;
        avg_denominator = avg_denominator + P_T(t)*Epsilon_T(t)*Epsilon_T(t)*average_chi_2_fixed_t;
    end
    
    % Final expression of the rate
    ratio = avg_nominator/sqrt(avg_denominator);

    covert_point = sqrt(2)*ratio;
end
