function [r1_vect, r2_vect, rk_vect, higher_than_sk_budgets, cpt, cpt_higher_than_max_P_X1_1, W_Y_X1_X2, W_Z_X1_X2] = run_simulation( ...
                                            uniform_P_T,  W_Y_X1_X2, W_Z_X1_X2, fixed_channel_laws, ...
                                            T_cardinality, X1_cardinality, X2_cardinality, Y_cardinality, X1_X2_cardinality, ...
                                            sk_budget, swap_bob_eve_channels, max_epsilon_t, optimize_epsilons_T, ...
                                            max_P_X1_1, compute_marginal_PY, ...
                                            DEBUG, DEBUG_covert, DEBUG_covert_theorem_contraints, tolerance, seed,...
                                            N_epochs, draw_covert_point)

    % seed for reproducibility
    rng(seed);
    
    % to count the number of times we didn't verify the conditions of the theorem inside the loop.
    cpt = 0;
    % to count the number of times P_{X1}(1) > max_P_X1_1
    cpt_higher_than_max_P_X1_1      = 0;

    % for fixed epsilon_T
    if (ismembertol(optimize_epsilons_T, 0, tolerance))
        disp('[INFO] Running experiment with fixed Epsilon_T')
        Epsilon_T           = max_epsilon_t*rand(T_cardinality,1); %[0.3 0.1 0.2 0.4];
    else
        disp('[INFO] Running experiment with different Epsilon_T')
    end

    if (fixed_channel_laws)
        disp('[INFO] Running experiment with fixed channel laws')
        if (swap_bob_eve_channels)
            disp('[INFO] Swaping Eve and Bob channels.');
            tmp = W_Y_X1_X2;
            W_Y_X1_X2 = W_Z_X1_X2;
            W_Z_X1_X2 = tmp;
        end
        [W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2] = CovertCommunication.extract_laws(W_Y_X1_X2, W_Z_X1_X2, X2_cardinality,Y_cardinality);
    else
        disp('[INFO] Running experiment with random channel laws')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fix a random channel law %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % MAC Channel matrix : We consider the channel Y = (2*X1 + X2 + w) % 4 where all variables are binary => \mathcal{Y} = \{0,1,2,3}
        % for each Y and fixed (x1, x2) we have a vector of probabilities [y=0, y=1, y=2, y=3]

        W_Y_X1_X2 = transpose(InformationTheory.generate_probability_vector(Y_cardinality, X1_X2_cardinality,1,0,1));
        W_Z_X1_X2 = transpose(InformationTheory.generate_probability_vector(Y_cardinality, X1_X2_cardinality,1,0,1));
        [W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2] = CovertCommunication.extract_laws(W_Y_X1_X2, W_Z_X1_X2, X2_cardinality,Y_cardinality);
        [verified_conditions, absolute_continuity_bob, absolute_continuity_eve, different_output_distributions_eve] = CovertCommunication.check_theorem_conditions(W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, DEBUG_covert_theorem_contraints);
          
        % loop while untill constraint on absolute continuity and difference are met 
        while (verified_conditions < 1)
            % to get matrices in the size we want, we do a transpose as
            % generate_probability_vector doesn't work if a < b.
            W_Y_X1_X2 = transpose(InformationTheory.generate_probability_vector(Y_cardinality, X1_X2_cardinality,1,0,1));
            W_Z_X1_X2 = transpose(InformationTheory.generate_probability_vector(Y_cardinality, X1_X2_cardinality,1,0,1));
            [W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2] = CovertCommunication.extract_laws(W_Y_X1_X2, W_Z_X1_X2, X2_cardinality,Y_cardinality);
            [verified_conditions, absolute_continuity_bob, absolute_continuity_eve, different_output_distributions_eve] = CovertCommunication.check_theorem_conditions(W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, DEBUG_covert_theorem_contraints);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rate Region simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % we add 3 because we add manually missing points:
    %     (r1=0,r2=0) (r1=0,max r2) (max r1, r2=0)
    N = N_epochs+3;
    % Covert user rates
    r1_vect = zeros(N,1);
    
    % Non covert user rates
    r2_vect = zeros(N,1);
%     r2_vect2 = zeros(N,1);

    % Secret key rates
    rk_vect = zeros(N,1);
    higher_than_sk_budgets = zeros(N_epochs,1);
    
    % covert region
    if draw_covert_point
        covert_points = zeros(N_epochs,1);
    end

    % if not fixed and we need to swap    
    if (~fixed_channel_laws)
        if (swap_bob_eve_channels)
            disp('[INFO] Swaping Eve and Bob channels.');
            tmp = W_Y_X1_X2;
            W_Y_X1_X2 = W_Z_X1_X2;
            W_Z_X1_X2 = tmp;
        end
    end

    disp('[INFO] Using the following channel matrix for Bob')
    disp(W_Y_X1_X2)
    
    disp('[INFO] Using the following channel matrix for Eve')
    disp(W_Z_X1_X2)

    % iterate over all input distributions.
    disp('[INFO] Starting simulation now...');
    for epoch = progress(1:N_epochs)

        % Generate new P_T
        if (uniform_P_T)
            P_T                 = (1/T_cardinality)*ones(T_cardinality,1);
        else
            P_T                 = InformationTheory.generate_probability_vector(T_cardinality,1,1,0,1);
        end
    
        if (ismembertol(optimize_epsilons_T, 1, tolerance))
        % The epsilon's that normalizes the probability of sending one for the covert user (we choose them summing to one for now but it's not necessary)
            Epsilon_T           = max_epsilon_t*rand(T_cardinality,1);
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Conditional input disbutions P_X1_mid_T and P_X2_mid_T %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        P_X1_mid_T = InformationTheory.generate_probability_vector(X1_cardinality, T_cardinality,1,0,1);

        if (X2_cardinality > 1)
            P_X2_mid_T          = InformationTheory.generate_probability_vector(X2_cardinality, T_cardinality,1,0,1);
        else
            % As it's deterministic
            P_X2_mid_T          = ones(X2_cardinality, T_cardinality);
        end
%         disp('P_X2_mid_T');
%         disp(P_X2_mid_T);
%         disp('size(P_X2_mid_T)');
%         disp(size(P_X2_mid_T));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Verify if conditions are met with these choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        [verified_conditions, absolute_continuity_bob, absolute_continuity_eve, different_output_distributions_eve] = CovertCommunication.check_theorem_conditions(W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, DEBUG_covert_theorem_contraints);
    
        if (verified_conditions < 1)
            cpt = cpt+1;

            if (DEBUG)
                disp('[ERROR-Constraint] Absolute continuity is not met! Check bellow which one you need to fix.');
                disp(['[INFO] absolute_continuity_bob : ', num2str(absolute_continuity_bob)]);
                disp(['[INFO] absolute_continuity_eve : ', num2str(absolute_continuity_eve)]);
                disp(['[INFO] different_output_distributions_eve : ', num2str(different_output_distributions_eve)]);
            end
            continue % go to next iteration
        end   
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Marginal distributions P_X1 and P_X2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        P_X1 = zeros(1,X1_cardinality);
        
        for x1=1:X1_cardinality
            avg = 0;
            for t=1:T_cardinality
                avg = avg + P_X1_mid_T(x1,t) * P_T(t);
            end
            P_X1(x1) = avg;
        end
    
        if (P_X1(1) > max_P_X1_1)
            cpt_higher_than_max_P_X1_1 = cpt_higher_than_max_P_X1_1 +1;
            continue
        end
    
    
        P_X2 = zeros(1,X2_cardinality);
        for x2=1:X2_cardinality
            avg = 0;
            for t=1:T_cardinality
                avg = avg + P_X2_mid_T(x2,t) * P_T(t);
            end
            P_X2(x2) = avg;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Marginal distributions P_Y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        if (compute_marginal_PY)
            P_X1_X2 = [P_X1(1)*P_X2(1) P_X1(1)*P_X2(2) P_X1(2)*P_X2(1) P_X1(2)*P_X2(2)]; % [P1(1)*P2(1) P1(0)*P2(1) P1(1)*P2(0) P1(1)*P2(1)]
       
            P_Y = zeros(1,Y_cardinality);
            for i=1:length(P_Y)
                P_Y(i) = dot(P_X1_X2, W_Y_X1_X2(:,i));
            end
        end

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Marginal distributions W_Y_X1_0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        W_Y_X1_0 = zeros(1, Y_cardinality);
        for y=1:Y_cardinality
            W_y_X1_0 = 0;
            for x2=1:X2_cardinality
                W_y_X1_0 = W_y_X1_0 + P_X2(x2) * W_Y_X1_0_X2(x2,y);
            end
            W_Y_X1_0(y) = W_y_X1_0;
        end        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Computing the rates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % start with secret key rate for efficiency of the running time
        rk = CovertCommunication.covert_sk_rate(P_T, P_X2_mid_T, Epsilon_T, W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, DEBUG_covert);
        
        if (rk < sk_budget)
            rk_vect(epoch) = rk;

            r1 = CovertCommunication.covert_message_rate(P_T, P_X2_mid_T, Epsilon_T, W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, Y_cardinality, X1_X2_cardinality, DEBUG_covert);
            r1_vect(epoch) = r1;
            
            r2 = CovertCommunication.conditional_MI(P_T, P_X2_mid_T, W_Y_X1_0_X2, X2_cardinality, Y_cardinality);
%             r2 = CovertCommunication.non_covert_rate(P_T, P_X2_mid_T, W_Y_X1_0, W_Y_X1_0_X2, X2_cardinality, DEBUG);
            r2_vect(epoch) = r2;
%             r2_vect2(epoch) = r22;

        else
            higher_than_sk_budgets(epoch) = rk; % just for analysis
            continue
        end
        
        
        if draw_covert_point
            covert_point = CovertCommunication.compute_covert_point(P_T, P_X2_mid_T, Epsilon_T, W_Z_X1_X2, X2_cardinality, Y_cardinality, DEBUG);
            covert_points(epoch) = covert_point;
        end
    end

    min_r2 = min(r2_vect);
    max_r2 = max(r2_vect);
    min_r1 = min(r1_vect);
    max_r1 = max(r1_vect);

    % add (r1=0,r2=0) 
    r1_vect(N_epochs+1)=0;
    r2_vect(N_epochs+1)=0;
    % add (r1=0,max r2) 
    r1_vect(N_epochs+2)=0;
    r2_vect(N_epochs+2)=max_r2;
    % add  (max r1, r2=0)
    r1_vect(N_epochs+3)=max_r1;
    r2_vect(N_epochs+3)=0;

    disp('--------------------------------------')    
    disp(['The minimum covert rate r1 is: ', num2str(min_r1)]);
    disp(['The maximum covert rate r1 is: ', num2str(max_r1)]);    
    disp('--------------------------------------')    
    disp(['The minimum non covert rate r2 is: ', num2str(min_r2)]);
    disp(['The maximum non covert rate r2 is: ', num2str(max_r2)]); 
    disp('--------------------------------------')
    disp(['[Exp-1] The number of time P_{X1}(1) >', num2str(max_P_X1_1), ' is: ', num2str(cpt_higher_than_max_P_X1_1), '. Which corresponds to a percentage of: ', num2str(100*cpt_higher_than_max_P_X1_1/N_epochs), '%']);  
    disp('--------------------------------------')
    disp(['[Exp-1] The number of time if did not verfiy the theorem conditions is: ', num2str(cpt), '. Which corresponds to a percentage of: ', num2str(100*cpt/N_epochs), '%']);  
    disp('--------------------------------------')

    
end