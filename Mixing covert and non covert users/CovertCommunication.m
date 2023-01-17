classdef CovertCommunication

    methods (Static)

        function DEBUG = set_debug(debug)
            CovertCommunication.DEBUG =debug;
        end

        function [W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2] = extract_laws(W_Y_X1_X2, W_Z_X1_X2, X2_cardinality,Y_cardinality)
            
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extracting distributions we need (W_Y_X1_0_X2 and W_Y_X1_1_X2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
            % Init with zeros
            W_Y_X1_0_X2 = zeros(X2_cardinality,Y_cardinality);
            W_Y_X1_1_X2 = zeros(X2_cardinality,Y_cardinality);
            W_Z_X1_0_X2 = zeros(X2_cardinality,Y_cardinality);
            W_Z_X1_1_X2 = zeros(X2_cardinality,Y_cardinality);
    
            % the channel matrix has 2*X2_cardinality lines from which 
            % the first X2_cardinality lines are for X1=0
            % followed by X2_cardinality lines for X1=1
    
            % For X1=0
            for x2=1:X2_cardinality
                % Bob
                W_Y_X1_0_X2(x2,:) = W_Y_X1_X2(x2,:);  % W_Y_X1_0_X2_x2
                % Eve
                W_Z_X1_0_X2(x2,:) = W_Z_X1_X2(x2,:); % W_Z_X1_0_X2_x2
            end
                            
            % For X1=1
            x2 = X2_cardinality+1;
            for i=1:X2_cardinality
                % Bob
                W_Y_X1_1_X2(i,:) = W_Y_X1_X2(x2,:); % W_Y_X1_1_X2_x2
                % Eve
                W_Z_X1_1_X2(i,:) = W_Z_X1_X2(x2,:); % W_Z_X1_1_X2_x2
                x2 = x2 + 1;
            end
        end
           
        
        function bool = check_absolute_continuity(W_Y_X1_1_X2, W_Y_X1_0_X2, X2_cardinality, DEBUG)
            % Input:
            %   W_Y_X1_1_X2       : channel law of the legitimate user W_{Y \mid X_1, X_2}
            %   W_Y_X1_0_X2       : channel law of the adversary W_{Z \mid X_1, X_2}
            % Output:
            %   Checks if W_Y_X1_1_X2 is absolutely continuous w.r.t W_Y_X1_0_X2 (W_Y_X1_1_X2 << W_Y_X1_0_X2)
            %   that is if for a fixed x2, for all y, W_Y_X1_1_X2(y \mid X1=1, X2=x2) = 0 if and only if W_Y_X1_1_X2(y \mid X1=1, X2=x2)
            tolerance = 1e-6;
            bool=1;
            for x2=1:X2_cardinality
                
                W_Y_X1_0_X2_x2 = W_Y_X1_0_X2(x2,:);
                W_Y_X1_1_X2_x2 = W_Y_X1_1_X2(x2,:);

                for y=1:length(W_Y_X1_1_X2_x2)
                    if ( ismembertol(W_Y_X1_1_X2_x2(y), 0, tolerance) &&  W_Y_X1_0_X2_x2(y) > 0)
                        % found one example, we can stop
                        bool = 0;
                        if DEBUG
                            disp('[ERROR-Constraint] W_Y/Z_X1_1_X2 is not absolutely continuous with W_Y/Z_X1_0_X2!');
                        end
                        return
                    end
                end
            end
        end

        function bool = check_different_distributions(W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, DEBUG)
            tolerance = 1e-6;
            bool=1;
            for x2=1:X2_cardinality

                W_Z_X1_0_X2_x2 = W_Z_X1_1_X2(x2,:);
                W_Z_X1_1_X2_x2 = W_Z_X1_0_X2(x2,:);

                for z=1:length(W_Z_X1_1_X2_x2)
                    if ( ismembertol(W_Z_X1_0_X2_x2(z), W_Z_X1_1_X2_x2(z), tolerance) )
                        % found one example, we can stop
                        bool = 0;
                        if DEBUG
                            disp('[ERROR-Constraint] W_Z_X1_1_X2 is not different from W_Z_X1_0_X2!');
                        end
                        return
                    end
                end
            end
        end % different distributions

        function without_secret_key = without_secret_key_condition(sk_budget, P_T, Epsilon_T, P_X2_mid_T, W_Y_X1_0_X2, W_Y_X1_1_X2, W_Z_X1_0_X2, W_Z_X1_1_X2, X2_cardinality, DEBUG) % what we need to do instead(P_T, Epsilon_T, P_X2_mid_T, W_Y_X1_0_X2, W_Y_X1_1_X2, W_Z_X1_0_X2, W_Z_X1_1_X2, X2_cardinality)
            % Input:
            %   W_Y_X1_0/1_X2       : channel law of the legitimate user W_{Y \mid X_1=0/1, X_2}
            %   W_Z_X1_0/1_X2       : channel law of the adversary W_{Z \mid X_1=0/1, X_2}
            %   P_X2_mid_T      : conditional probability vector P_{X_2 \mid T}
            % Output:
            %   Checks if D(W_Y_X1_1_X2 \mid W_Y_X1_0_X2) >= D(W_Z_X1_1_X2 \mid W_Z_X1_0_X2) on average over P_T \epsilon_Tcheck_theorem_conditions

            relative_entropy_bob_vect = zeros(X2_cardinality,1);
            relative_entropy_eve_vect = zeros(X2_cardinality,1);
            
            for x2=1:X2_cardinality

                W_Y_X1_0_X2_x2 = W_Y_X1_0_X2(x2,:);
                W_Y_X1_1_X2_x2 = W_Y_X1_1_X2(x2,:);

                W_Z_X1_0_X2_x2 = W_Z_X1_0_X2(x2,:);
                W_Z_X1_1_X2_x2 = W_Z_X1_1_X2(x2,:);

                relative_entropy_bob = InformationTheory.relative_entropy(W_Y_X1_1_X2_x2, W_Y_X1_0_X2_x2);
                relative_entropy_eve = InformationTheory.relative_entropy(W_Z_X1_1_X2_x2, W_Z_X1_0_X2_x2);

                relative_entropy_bob_vect(x2) = relative_entropy_bob;
                relative_entropy_eve_vect(x2) = relative_entropy_eve;

            end

            difference_relative_entropy_vect = relative_entropy_eve_vect-relative_entropy_bob_vect;
            if DEBUG
                disp('difference_relative_entropy_vect')
                disp(difference_relative_entropy_vect)
            end

            % this is what I think should be done but we varry P_X2_mid_T
            % and P_T so I don't know how we can check this condition in
            % our setup. Bloch didn't have this averaging. I am taking a
            % simpler thing, to ensure 100% it's feasible, both should be
            % positive so that the average is two. When one of them is
            % negative, we cannot conclude. So basically, we don't take the
            % risk and try again untill we find a channel law where both
            % are positive

            avg = 0;
            for t=1:length(P_T)
                % computing the averages of the relative entropy and chi_2 for a fixed t
                average_relative_entropy_fixed_t = dot(P_X2_mid_T(:,t), difference_relative_entropy_vect);
                
                % taking the average over T
                avg = avg + P_T(t)*Epsilon_T(t)*average_relative_entropy_fixed_t;
            end

            if (avg <= sk_budget)
                without_secret_key = 1;
            else
                without_secret_key = 0;
            end            
        end %end check without secret key condition
        
        function [bool, absolute_continuity_bob, absolute_continuity_eve, different_output_distributions_eve] = check_theorem_conditions(W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, DEBUG)
            tolerance = 1e-8;
            bool = 0;

            % checking conditions
            absolute_continuity_bob = CovertCommunication.check_absolute_continuity(W_Y_X1_1_X2, W_Y_X1_0_X2, X2_cardinality, DEBUG);
            absolute_continuity_eve = CovertCommunication.check_absolute_continuity(W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, DEBUG);
            different_output_distributions_eve = CovertCommunication.check_different_distributions(W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, DEBUG);

%             if (swap_channels_for_need_sk)
%                 without_secret_key_condition = 1;
%             else
%                 without_secret_key_condition = CovertCommunication.without_secret_key_condition(sk_budget, P_T, Epsilon_T, P_X2_mid_T, W_Y_X1_0_X2, W_Y_X1_1_X2, W_Z_X1_0_X2, W_Z_X1_1_X2, X2_cardinality, DEBUG);
% && ismembertol(without_secret_key_condition, 1, tolerance))
%             end
            if (ismembertol(absolute_continuity_bob, 1, tolerance) && ismembertol(absolute_continuity_eve, 1, tolerance) && ismembertol(different_output_distributions_eve, 1, tolerance)) 
                bool = 1;
            end
            if DEBUG
                disp(['[INFO] absolute_continuity_bob : ', num2str(absolute_continuity_bob)]);
                disp(['[INFO] absolute_continuity_eve : ', num2str(absolute_continuity_eve)]);
                disp(['[INFO] different_output_distributions_eve : ', num2str(different_output_distributions_eve)]);
%                 disp(['[INFO] without_secret_key_condition : ', num2str(without_secret_key_condition)]);
                disp(['[INFO] bool : ', num2str(bool)]);
            end
            
        end
        
        function r1 = covert_message_rate(P_T, P_X2_mid_T, Epsilon_T, W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, Y_cardinality, X1_X2_cardinality, DEBUG)
            % Input:
            %   P_T             : a probabilities vector for the time sharing rv
            %   P_X2_mid_T      : conditional probability vector P_{X_2 \mid T}
            %   Epsilon_T       : vector of \epsilon_{t} values 
            %   W_Y_X1_X2       : channel law of the legitimate user W_{Y \mid X_1, X_2}
            %   W_Z_X1_X2       : channel law of the adversary W_{Z \mid X_1, X_2}
            % Output:
            %   The rate of the covert user U1
            
            % We first compute the relative entropy for all X_2 for the legitimate user --------------------------------------------------------
            relative_entropy_vect = zeros(X2_cardinality,1);                
            
            for x2=1:X2_cardinality
                W_Y_X1_0_X2_x2 = W_Y_X1_0_X2(x2,:);
                W_Y_X1_1_X2_x2 = W_Y_X1_1_X2(x2,:);
                relative_entropy = InformationTheory.relative_entropy(W_Y_X1_1_X2_x2, W_Y_X1_0_X2_x2);
                relative_entropy_vect(x2) = relative_entropy;
            end

            % Then, we compute the chi_2 distance for all X_2 for Eve's channel --------------------------------------------------------
            chi_2_vect = zeros(X2_cardinality,1);     

            for x2=1:X2_cardinality
                W_Z_X1_0_X2_x2 = W_Z_X1_0_X2(x2,:);
                W_Z_X1_1_X2_x2 = W_Z_X1_1_X2(x2,:);
                chi_2 = InformationTheory.chi_k_distance(W_Z_X1_1_X2_x2, W_Z_X1_0_X2_x2, 2);
                chi_2_vect(x2) = chi_2;
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
            r1 = sqrt(2)*ratio;

        end % end covert_message_rate

        function rk = covert_sk_rate(P_T, P_X2_mid_T, Epsilon_T, W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, DEBUG)
            % Input:
            %   P_T             : a probabilities vector for the time sharing rv
            %   P_X2_given_T    : conditional probability vector P_{X_2 \mid T}
            %   epsilons        : vector of \epsilon_{t} values 
            %   W_Y_X1_X2       : channel law of the legitimate user W_{Y \mid X_1, X_2}
            %   W_Z_X1_X2       : channel law of the adversary W_{Z \mid X_1, X_2}
            % Output:
            %   The secret key rate
            for x2=1:X2_cardinality
                difference_Ds(x2)=InformationTheory.relative_entropy(W_Z_X1_1_X2(x2,:), W_Z_X1_0_X2(x2,:))-InformationTheory.relative_entropy(W_Y_X1_1_X2(x2,:), W_Y_X1_0_X2(x2,:));
                chi_2 (x2)= InformationTheory.chi_k_distance(W_Z_X1_1_X2(x2,:), W_Z_X1_0_X2(x2,:), 2);
            end
            for t=1:length(P_T)
                numerator_T(t)= Epsilon_T(t)*difference_Ds*P_X2_mid_T(:,t) ;
                denominator_T(t)= chi_2* P_X2_mid_T(:,t)*(Epsilon_T(t))^2;    
            end
            
          numerator=numerator_T*P_T;
          denominator=denominator_T*P_T;
          rk=sqrt(2)*numerator./(sqrt(denominator));

%             % We first compute the relative entropy for all X_2 for the legitimate user --------------------------------------------------------
%             relative_entropy_vect_eve   = zeros(X2_cardinality,1);
%             relative_entropy_vect_bob   = zeros(X2_cardinality,1);
%             
%             for x2=1:X2_cardinality
%                 % Bob
%                 W_Y_X1_0_X2_x2 = W_Y_X1_0_X2(x2,:);
%                 W_Y_X1_1_X2_x2 = W_Y_X1_1_X2(x2,:);
%                 relative_entropy_bob = InformationTheory.relative_entropy(W_Y_X1_1_X2_x2, W_Y_X1_0_X2_x2);
%                 relative_entropy_vect_bob(x2) = relative_entropy_bob;
%                 
%                 % Eve
%                 W_Z_X1_0_X2_x2 = W_Z_X1_0_X2(x2,:);
%                 W_Z_X1_1_X2_x2 = W_Z_X1_1_X2(x2,:);
%                 relative_entropy_eve = InformationTheory.relative_entropy(W_Z_X1_1_X2_x2, W_Z_X1_0_X2_x2);
%                 relative_entropy_vect_eve(x2) = relative_entropy_eve;
%             end
%             
%             % take the difference in vectors
%             relative_entropy_vect = relative_entropy_vect_eve - relative_entropy_vect_bob;
% 
%             % Then, we compute the chi_2 distance for all X_2 for Eve's channel --------------------------------------------------------
%             chi_2_vect = zeros(X2_cardinality,1);
%             
%             for x2=1:X2_cardinality
%                 W_Z_X1_0_X2_x2 = W_Z_X1_0_X2(x2,:);
%                 W_Z_X1_1_X2_x2 = W_Z_X1_1_X2(x2,:);
%                 chi_2 = InformationTheory.chi_2_distance(W_Z_X1_1_X2_x2, W_Z_X1_0_X2_x2);
%                 chi_2_vect(x2) = chi_2;
%             end
% 
%             % Finaly, we compute the average of the nominator and denominator w.r.t P_T and P_{X_2 \mid T} --------------------------------------------------------
%             avg_nominator = 0;
%             avg_denominator = 0;
%             for t=1:length(P_T)
%                 % computing the averages of the relative entropy and chi_2 for a fixed t
%                 average_relative_entropy_fixed_t = dot(P_X2_mid_T(:,t), relative_entropy_vect);
%                 average_chi_2_fixed_t = dot(P_X2_mid_T(:,t), chi_2_vect);
%                 
%                 % taking the average over T
%                 avg_nominator = avg_nominator + P_T(t)*Epsilon_T(t)*average_relative_entropy_fixed_t;
%                 avg_denominator = avg_denominator + P_T(t)*Epsilon_T(t)*Epsilon_T(t)*average_chi_2_fixed_t;
%             end
%             
%             % Final expression of the rate
%             ratio = avg_nominator/sqrt(avg_denominator);
% 
%             rk = sqrt(2)*ratio;
% 
%             if (rk < 0)
%                 % we don't need a secret key if the rate is negative
%                 % (because Bob's channel is better than Eve's channel in such case)
%                 rk = 0; 
%             end

        end % end covert_sk_rate        
        
         function I = conditional_MI(P_T, P_X2_mid_T, W_Y_X1_0_X2, X2_cardinality, Y_cardinality)
            % Input:
            %   P_T             : a probabilities vector for the time sharing rv
            %   P_X2_mid_T    : conditional probability vector P_{X_2 \mid T} 
            %   W_Y_X1_0_X2     : channel law 
            % Output: mutual information I(X2;Y \mid X1=0, T)
            %   The rate of the non-covert user U2
            % R_2 = \sum_{(t,x_2)} p(t)p(x_2 \mid t) \mathbb{D}(W_{Y \mid X_2}(\cdot \mid X_2=x_2) \| \sum_x2' P_X2(x2') W_{Y \mid X2}(\cdot \mid X_2=x2'))
            
            
            for t=1:length(P_T)
                PY(t,:)= (P_X2_mid_T(:,t)')*W_Y_X1_0_X2;
                H_Y_X_0_T(t)=-sum( PY(t,:).*log2(PY(t,:)));
                for x2=1:X2_cardinality
                    H_Y_X1_0_X2_T(x2,t)=-sum( W_Y_X1_0_X2(x2,:).*log2(W_Y_X1_0_X2(x2,:)));
                end
                avg_H_Y_X1_0_X2_T(t)= (P_X2_mid_T(:,t)')* H_Y_X1_0_X2_T(:,t);
                I_T(t)=H_Y_X_0_T(t) - avg_H_Y_X1_0_X2_T(t);
            end
              
              I=I_T*P_T;
              
         end

        function r2 = non_covert_rate(P_T, P_X2_mid_T, W_Y_X1_0, W_Y_X1_0_X2, X2_cardinality, DEBUG)
            % Input:
            %   P_T             : a probabilities vector for the time sharing rv
            %   P_X2_mid_T    : conditional probability vector P_{X_2 \mid T} 
            %   P_Y             : a probabilities vector for the output Y rv
            %   W_Y_X1_0_X2     : channel law of the legitimate user W_{Y \mid X_1=0, X_2}
            %   W_Y_X1_0_X2     : channel law of the legitimate user W_{Y \mid X_1=0, X_2}
            % Output:
            %   The rate of the non-covert user U2
            % R_2 = \sum_{(t,x_2)} p(t)p(x_2 \mid t) \mathbb{D}(W_{Y \mid X_1, X_2}(\cdot \mid X_1=0, X_2=x_2) \| W_{Y \mid X_1}(\cdot \mid X_1=0))
                       
            relative_entropy_vect = zeros(X2_cardinality,1);

            for x2=1:X2_cardinality
                W_Y_X1_0_X2_x2 = W_Y_X1_0_X2(x2,:);
                % We use abs()<1e-15 to avoid problems with precision, it can 
                % happen that they are the same but due to precision, they 
                % create a problem and we get a negative relative entropy. 
                % To reproduce, use X2_card = 1 and T_card = 1 
                % vs X2_card = 1 and T_card= 10.
                relative_entropy = InformationTheory.relative_entropy(W_Y_X1_0_X2_x2, W_Y_X1_0); %fix(InformationTheory.relative_entropy(W_Y_X1_0_X2_x2, W_Y_X1_0));
                if (abs(relative_entropy) < 1e-15 && relative_entropy~=0)
                    disp(['[INFO] W_Y_X1_0_X2_x2 and W_Y_X1_0 are not the same in precision leading to a relative_entropy of', num2str(relative_entropy)])
                    relative_entropy = 0;
                end
                relative_entropy_vect(x2) = relative_entropy;
            end

            r2 =0;
            for t=1:length(P_T)
                average_relative_entropy_fixed_t = dot(P_X2_mid_T(:,t), relative_entropy_vect);
                r2 = r2 + P_T(t)*average_relative_entropy_fixed_t;
            end
        end % end non covert rate
        
        % computes the rate as I(X2;Y \mid T) which is asymptotically I(X2;Y \mid X1=0, T)
        function r2 = non_covert_rate_off(P_T, P_X2_mid_T, P_X2, W_Y_X2, Y_cardinality, DEBUG)
            % Input:
            %   P_T             : a probabilities vector for the time sharing rv
            %   P_X2_mid_T    : conditional probability vector P_{X_2 \mid T} 
            %   P_Y             : a probabilities vector for the output Y rv
            %   W_Y_X1_0_X2     : channel law of the legitimate user W_{Y \mid X_1=0, X_2}
            %   W_Y_X1_0_X2     : channel law of the legitimate user W_{Y \mid X_1=0, X_2}
            % Output:
            %   The rate of the non-covert user U2
            % R_2 = \sum_{(t,x_2)} p(t)p(x_2 \mid t) \mathbb{D}(W_{Y \mid X_2}(\cdot \mid X_2=x_2) \| \sum_x2' P_X2(x2') W_{Y \mid X2}(\cdot \mid X_2=x2'))
            
            r2 =0;
            for t=1:length(P_T)
                for x2=1:length(P_X2)
                    for y=1:Y_cardinality
                        sum_y_x2 = 0;
                        for x22=1:length(P_X2)
                            sum_y_x2 = sum_y_x2 + P_X2(x22)*W_Y_X2(x22,y);
                        end
                        r2 = r2 +P_T(t)*P_X2_mid_T(x2,t)*W_Y_X2(x2,y)*log2(W_Y_X2(x2,y)/sum_y_x2);
                    end
                end
            end
        end % end non covert rate    

    end % end methods
end