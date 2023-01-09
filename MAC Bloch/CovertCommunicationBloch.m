classdef CovertCommunicationBloch

    methods (Static)
        
        function covert_points = compute_covert_point(W_Z_X1_X2, rhos_star, n_users, Y_cardinality)
            covert_points = zeros(1,2);
            
            Q1 = W_Z_X1_X2(3,:);
            Q2 = W_Z_X1_X2(2,:);
            Q0 = W_Z_X1_X2(1,:);

            Q = zeros(n_users, Y_cardinality);
            Q(1,:) = Q1;
            Q(2,:) = Q2;

            relative_entropy_Q1_Q0 = InformationTheory.relative_entropy(Q1, Q0);
            relative_entropy_Q2_Q0 = InformationTheory.relative_entropy(Q2, Q0);
            
            chi_2_distance = CovertCommunicationBloch.chi_2_distance_bloch(Q, Q0, rhos_star);
             
            ratio = sqrt(2)/chi_2_distance;
            covert_point_1 = ratio * rhos_star(1) * relative_entropy_Q1_Q0;
            covert_point_2 = ratio * rhos_star(2) * relative_entropy_Q2_Q0;

            covert_points(1) = covert_point_1;
            covert_points(2) = covert_point_2;
        end

        function chi = chi_2_distance_bloch(Q, Q0, rhos)
            % Input:
            %   Qk: a probabilities matrix of size (2,Z_cardinality)
            %   Q0: a probabilities vector
            % Output:
            %   The Chi 2 distance \Chi_{2}(Qk||Q0)
            assert(numel(Q(1,:)) == numel(Q0));
            
            chi = 0;
            for z=1:length(Q(1,:))
                % computing the average nominator
                difference_power_2_avg_rho = 0;
                for k=1:length(rhos)
                    difference_power_2_avg_rho = difference_power_2_avg_rho + rhos(k)*(Q(k,z)-Q0(z));
                end
                difference_power_2 = difference_power_2_avg_rho.^2;

                % denominator
                Q0_z = Q0(z);

                % ratio
                ratio = difference_power_2/Q0_z;

                % sum over z
                chi = chi + ratio;
            end
                        
        end % end chi_2_distance_bloch

        function without_secret_key = without_secret_key_condition(P0, P1, P2, Q0, Q1, Q2)
            % Output:
            %   Checks if D(PK \mid P0) >= D(Qk \mid Q0)

            relative_entropy_P1_P0 = InformationTheory.relative_entropy(P1, P0);
            relative_entropy_P2_P0 = InformationTheory.relative_entropy(P2, P0);
            relative_entropy_Q1_P0 = InformationTheory.relative_entropy(Q1, Q0);
            relative_entropy_Q2_P0 = InformationTheory.relative_entropy(Q2, Q0);
            
            without_secret_key_user_1 = relative_entropy_P1_P0 >= relative_entropy_Q1_P0;
            without_secret_key_user_2 = relative_entropy_P2_P0 >= relative_entropy_Q2_P0;

            if (without_secret_key_user_1 && without_secret_key_user_2)
                without_secret_key = 1;
            else
                without_secret_key = 0;
            end

        end %end check without secret key condition

        function bool = check_absolute_continuity(Qk, Q0, DEBUG)
            % Output:
            %   Checks if Qk is absolutely continuous w.r.t Q0 (Qk << Q0)
            %   that is if for a fixed x2, for all y, Qk = 0 if and only if Q0 = 0
            tolerance = 1e-6;
            bool=1;

            for z=1:length(Qk)
                if ( ismembertol(Qk(z), 0, tolerance) &&  Q0(z) > 0)
                    % found one example, we can stop
                    bool = 0;
                    if DEBUG
                        disp('[ERROR-Constraint] Qk is not absolutely continuous with Q0!');
                    end
                    return
                end
            end
        end % end check_absolute_continuity
        
        function bool = check_theorem_conditions_bloch(W_Y_X1_X2, W_Z_X1_X2, DEBUG)
            tolerance = 1e-8;
            bool = 0;
            
            P0 = W_Y_X1_X2(1,:); %W_Y_X1_0_X2_0 when both users are off
            Q0 = W_Z_X1_X2(1,:); %W_Z_X1_0_X2_0 when both users are off

            P1 = W_Y_X1_X2(3,:); %W_Y_X1_1_X2_0 when user 1 is communicating while user 2 is off
            P2 = W_Y_X1_X2(2,:); %W_Y_X1_0_X2_1 when user 2 is communicating while user 1 is off

            Q1 = W_Z_X1_X2(3,:); %W_Z_X1_1_X2_0 when user 1 is communicating while user 2 is off
            Q2 = W_Z_X1_X2(2,:); %W_Z_X1_0_X2_1 when user 2 is communicating while user 1 is off
            %Q12 = W_Z_X1_X2(4,:); %W_Z_X1_1_X2_1 when both users are communicating

            % checking conditions
            absolute_continuity_eve_Q1_Q0 = CovertCommunicationBloch.check_absolute_continuity(Q1, Q0, DEBUG);
            absolute_continuity_eve_Q2_Q0 = CovertCommunicationBloch.check_absolute_continuity(Q2, Q0, DEBUG);
            %absolute_continuity_eve_Q12_Q0 = check_absolute_continuity_bloch(Q12, Q0, DEBUG);
            total = 2;
            
            absolute_continuity_eve = absolute_continuity_eve_Q1_Q0 + absolute_continuity_eve_Q2_Q0 ; %+ absolute_continuity_eve_Q12_Q0;
            without_secret_key = CovertCommunicationBloch.without_secret_key_condition(P0, P1, P2, Q0, Q1, Q2);
            
            if (ismembertol(absolute_continuity_eve, total, tolerance) && ismembertol(without_secret_key, 1, tolerance))
                bool = 1;
            end
            if DEBUG
                disp(['[INFO] absolute_continuity_eve_Q1_Q0 : ', num2str(absolute_continuity_eve_Q1_Q0)]);
                disp(['[INFO] absolute_continuity_eve_Q2_Q0 : ', num2str(absolute_continuity_eve_Q2_Q0)]);
%                 disp(['[INFO] absolute_continuity_eve_Q12_Q0 : ', num2str(absolute_continuity_eve_Q12_Q0)]);

                disp(['[INFO] without_secret_key : ', num2str(without_secret_key)]);
                disp(['[INFO] bool : ', num2str(bool)]);
            end
            
        end % end check_theorem_conditions_bloch

        function users_rates = covert_message_rate(rhos, W_Y_X1_X2, W_Z_X1_X2, DEBUG)
            % Output:
            %   The rate of the covert users
            n_users = 2;
            Y_cardinality = 4;
            users_rates = zeros(n_users,1); % 2 users

            P0 = W_Y_X1_X2(1,:); %W_Y_X1_0_X2_0 when both users are off
            Q0 = W_Z_X1_X2(1,:); %W_Z_X1_0_X2_0 when both users are off

            P1 = W_Y_X1_X2(3,:); %W_Y_X1_1_X2_0 when user 1 is communicating while user 2 is off
            P2 = W_Y_X1_X2(2,:); %W_Y_X1_0_X2_1 when user 2 is communicating while user 1 is off

            Q1 = W_Z_X1_X2(3,:); %W_Z_X1_1_X2_0 when user 1 is communicating while user 2 is off
            Q2 = W_Z_X1_X2(2,:); %W_Z_X1_0_X2_1 when user 2 is communicating while user 1 is off
            
            % Grouping them in a matrice
            P = zeros(n_users, Y_cardinality);
            P(1,:) = P1;
            P(2,:) = P2;

            Q = zeros(n_users, Y_cardinality);
            Q(1,:) = Q1;
            Q(2,:) = Q2;

            chi_rho = CovertCommunicationBloch.chi_2_distance_bloch(Q, Q0, rhos);
            
            for k=1:n_users                
                relative_entropy = InformationTheory.relative_entropy(P(k,:), P0);
                ratio = sqrt(2)/sqrt(chi_rho);
                users_rates(k) = ratio*rhos(k)*relative_entropy;
            end

            if DEBUG
                disp('users_rates');
                disp(users_rates);
            end
        end % end covert_message_rate

    end % end methods
end