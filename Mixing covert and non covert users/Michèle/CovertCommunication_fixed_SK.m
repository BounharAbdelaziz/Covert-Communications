classdef CovertCommunication_fixed_SK

    methods (Static)

        function DEBUG = set_debug(debug)
            CovertCommunication.DEBUG =debug;
        end
        
           
        function covert_point = compute_covert_point(P_T, P_X2_mid_T, Epsilon_T, W_Z_X1_X2, X2_cardinality, Y_cardinality, DEBUG)
            % Input:
            %   P_T             : a probabilities vector for the time sharing rv
            %   P_X2_given_T    : conditional probability vector P_{X_2 \mid T}
            %   epsilons        : vector of \epsilon_{t} values 
            %   W_Y_X1_X2       : channel law of the legitimate user W_{Y \mid X_1, X_2}: first X2_cardinality rows contain X_1=0 values and the latter rows the X_1=1 values
            %   W_Z_X1_X2       : channel law of the adversary W_{Z \mid X_1, X_2}:first X2_cardinality rows contain X_1=0 values and the latter rows the X_1=1 values
            

            % We first compute the relative entropy for all X_2 for the legitimate user --------------------------------------------------------
            relative_entropy_vect   = zeros(X2_cardinality,1);
            
         
            % Then, we compute the chi_2 distance for all X_2 for Eve's channel --------------------------------------------------------
            chi_2_vect = zeros(X2_cardinality,1);
            
            for x2=1:X2_cardinality
                W_Z_X1_0_X2_x2 = W_Z_X1_0_X2(x2,:);
                W_Z_X1_1_X2_x2 = W_Z_X1_1_X2(x2+X2_cardinality,:);
                chi_2 = InformationTheory.chi_k_distance(W_Z_X1_1_X2_x2, W_Z_X1_0_X2_x2, 2);
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



   
        function [bool, absolute_continuity_bob, absolute_continuity_eve, different_output_distributions_eve] = check_theorem_conditions(W_Y_X1_X2, W_Z_X1_X2, P_T, Epsilon_T, P_X2_mid_T, X2_cardinality, Y_cardinality, X1_X2_cardinality, swap_channels_for_need_sk, DEBUG)
            tolerance = 1e-8;
            bool = 0;
        
             % the channel matrices have 2*X2_cardinality lines from which 
                % the first X2_cardinality lines are for X1=0
                % followed by X2_cardinality lines for X1=1

            absolute_continuity_bob= (sum(sum((W_Y_X1_1_X2([1:X2_cardinality],:) >0) | (W_Y_X1_1_X2([1+X2_cardinality:2*X2_cardinality],:) ==0)))>0);
            absolute_continuity_eve= (sum(sum((W_Z_X1_1_X2([1:X2_cardinality],:) >0) | (W_Z_X1_1_X2([1+X2_cardinality:2*X2_cardinality],:) ==0)))>0);
            different_output_distributions_eve = sum(sum(( W_Z_X1_1_X2([1:X2_cardinality],:)  - W_Z_X1_1_X2([1+X2_cardinality:2*X2_cardinality] ))==0));  
         

            
               
            


            if DEBUG
                disp(['[INFO] absolute_continuity_bob : ', num2str(absolute_continuity_bob)]);
                disp(['[INFO] absolute_continuity_eve : ', num2str(absolute_continuity_eve)]);
                disp(['[INFO] different_output_distributions_eve : ', num2str(different_output_distributions_eve)]);
                disp(['[INFO] without_secret_key_condition : ', num2str(without_secret_key_condition)]);
                disp(['[INFO] bool : ', num2str(bool)]);
            end
            
        end
        
        
        function r1 = covert_message_rate(P_T, P_X2_mid_T, Epsilon_T, W_Y_X1_X2, W_Z_X1_X2, X2_cardinality, Y_cardinality, X1_X2_cardinality, DEBUG)
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
            W_Y_X1_0_X2 = zeros(X2_cardinality,Y_cardinality);
            W_Y_X1_1_X2 = zeros(X2_cardinality,Y_cardinality);

            if (X2_cardinality > 0)
                % the channel matrix has 2*X2_cardinality lines from which 
                % the first X2_cardinality lines are for X1=0
                % followed by X2_cardinality lines for X1=1

                % For X1=0
                for x2=1:X2_cardinality
                    W_Y_X1_0_X2(x2,:) = W_Y_X1_X2(x2,:);  % W_Y_X1_0_X2_x2
                end
                                
                % For X1=1
                x2 = X2_cardinality+1;
                for i=1:X2_cardinality
                    W_Y_X1_1_X2(i,:) = W_Y_X1_X2(x2,:); % W_Y_X1_1_X2_x2
                    x2 = x2 + 1;
                end
                
            else               
                W_Y_X1_0_X2(1,:) = W_Y_X1_X2(1,:); % W_Y_X1_0_X2_0
                W_Y_X1_0_X2(2,:) = W_Y_X1_X2(2,:); % W_Y_X1_0_X2_1
    
                W_Y_X1_1_X2(1,:) = W_Y_X1_X2(3,:); % W_Y_X1_1_X2_0
                W_Y_X1_1_X2(2,:) = W_Y_X1_X2(4,:); % W_Y_X1_1_X2_1
            end

            for x2=1:X2_cardinality
                W_Y_X1_0_X2_x2 = W_Y_X1_0_X2(x2,:);
                W_Y_X1_1_X2_x2 = W_Y_X1_1_X2(x2,:);
                relative_entropy = InformationTheory.relative_entropy(W_Y_X1_1_X2_x2, W_Y_X1_0_X2_x2);
                relative_entropy_vect(x2) = relative_entropy;
            end
            if DEBUG
                disp('')
                disp('relative_entropy_vect is: ');
                disp(relative_entropy_vect)
            end
            % Then, we compute the chi_2 distance for all X_2 for Eve's channel --------------------------------------------------------
            chi_2_vect = zeros(X2_cardinality,1);
            
            W_Z_X1_0_X2 = zeros(X2_cardinality,Y_cardinality);
            W_Z_X1_1_X2 = zeros(X2_cardinality,Y_cardinality);

            if (X2_cardinality > 0)
                % the channel matrix has 2*X2_cardinality lines from which 
                % the first X2_cardinality lines are for X1=0
                % followed by X2_cardinality lines for X1=1

                % For X1=0
                for x2=1:X2_cardinality
                    W_Z_X1_0_X2(x2,:) = W_Z_X1_X2(x2,:);  % W_Z_X1_0_X2_x2
                end
                                
                % For X1=1
                x2 = X2_cardinality+1;
                for i=1:X2_cardinality
                    W_Z_X1_1_X2(i,:) = W_Z_X1_X2(x2,:); % W_Z_X1_1_X2_x2
                    x2 = x2 + 1;
                end                
            else            
                W_Z_X1_0_X2(1,:) = W_Z_X1_X2(1,:); % W_Z_X1_0_X2_0
                W_Z_X1_0_X2(2,:) = W_Z_X1_X2(2,:); % W_Z_X1_0_X2_1
    
                W_Z_X1_1_X2(1,:) = W_Z_X1_X2(3,:); % W_Z_X1_1_X2_0
                W_Z_X1_1_X2(2,:) = W_Z_X1_X2(4,:); % W_Z_X1_1_X2_1
            end

            for x2=1:X2_cardinality
                W_Z_X1_0_X2_x2 = W_Z_X1_0_X2(x2,:);
                W_Z_X1_1_X2_x2 = W_Z_X1_1_X2(x2,:);
                chi_2 = InformationTheory.chi_k_distance(W_Z_X1_1_X2_x2, W_Z_X1_0_X2_x2, 2);
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
            r1 = sqrt(2)*ratio;

        end % end covert_message_rate

        function rk = covert_sk_rate(P_T, P_X2_mid_T, Epsilon_T, W_Y_X1_X2, W_Z_X1_X2, X2_cardinality, Y_cardinality, X1_X2_cardinality, DEBUG)
            % Input:
            %   P_T             : a probabilities vector for the time sharing rv
            %   P_X2_given_T    : conditional probability vector P_{X_2 \mid T}
            %   epsilons        : vector of \epsilon_{t} values 
            %   W_Y_X1_X2       : channel law of the legitimate user W_{Y \mid X_1, X_2}
            %   W_Z_X1_X2       : channel law of the adversary W_{Z \mid X_1, X_2}
            % Output:
            %   The secret key rate

            % We first compute the relative entropy for all X_2 for the legitimate user --------------------------------------------------------
            relative_entropy_vect       = zeros(X2_cardinality,1);
            relative_entropy_vect_eve   = zeros(X2_cardinality,1);
            relative_entropy_vect_bob   = zeros(X2_cardinality,1);

            % Bob's channel
            W_Y_X1_0_X2 = zeros(X2_cardinality,Y_cardinality);
            W_Y_X1_1_X2 = zeros(X2_cardinality,Y_cardinality);
            % Eve's channel
            W_Z_X1_0_X2 = zeros(X2_cardinality,Y_cardinality);
            W_Z_X1_1_X2 = zeros(X2_cardinality,Y_cardinality);

            if (X2_cardinality > 0)
                % the channel matrix has 2*X2_cardinality lines from which 
                % the first X2_cardinality lines are for X1=0
                % followed by X2_cardinality lines for X1=1

                % For X1=0
                for x2=1:X2_cardinality
                    W_Y_X1_0_X2(x2,:) = W_Y_X1_X2(x2,:);  % W_Y_X1_0_X2_x2
                    W_Z_X1_0_X2(x2,:) = W_Z_X1_X2(x2,:); % W_Z_X1_0_X2_x2
                end
                                
                % For X1=1
                x2 = X2_cardinality+1;
                for i=1:X2_cardinality
                    W_Y_X1_1_X2(i,:) = W_Y_X1_X2(x2,:); % W_Y_X1_1_X2_x2
                    W_Z_X1_1_X2(i,:) = W_Z_X1_X2(x2,:); % W_Z_X1_1_X2_x2
                    x2 = x2 + 1;
                end

            else
                W_Y_X1_0_X2(1,:) = W_Y_X1_X2(1,:); % W_Y_X1_0_X2_0
                W_Y_X1_0_X2(2,:) = W_Y_X1_X2(2,:); % W_Y_X1_0_X2_1
                W_Y_X1_1_X2(1,:) = W_Y_X1_X2(3,:); % W_Y_X1_1_X2_0
                W_Y_X1_1_X2(2,:) = W_Y_X1_X2(4,:); % W_Y_X1_1_X2_1
                
                W_Z_X1_0_X2(1,:) = W_Z_X1_X2(1,:); % W_Z_X1_0_X2_0
                W_Z_X1_0_X2(2,:) = W_Z_X1_X2(2,:); % W_Z_X1_0_X2_1
                W_Z_X1_1_X2(1,:) = W_Z_X1_X2(3,:); % W_Z_X1_1_X2_0
                W_Z_X1_1_X2(2,:) = W_Z_X1_X2(4,:); % W_Z_X1_1_X2_1
            end

            for x2=1:X2_cardinality
                % Bob
                W_Y_X1_0_X2_x2 = W_Y_X1_0_X2(x2,:);
                W_Y_X1_1_X2_x2 = W_Y_X1_1_X2(x2,:);
                relative_entropy_bob = InformationTheory.relative_entropy(W_Y_X1_1_X2_x2, W_Y_X1_0_X2_x2);
                relative_entropy_vect_bob(x2) = relative_entropy_bob;
                
                % Eve
                W_Z_X1_0_X2_x2 = W_Z_X1_0_X2(x2,:);
                W_Z_X1_1_X2_x2 = W_Z_X1_1_X2(x2,:);
                relative_entropy_eve = InformationTheory.relative_entropy(W_Z_X1_1_X2_x2, W_Z_X1_0_X2_x2);
                relative_entropy_vect_eve(x2) = relative_entropy_eve;
            end
            
            % take the difference in vectors
            relative_entropy_vect = relative_entropy_vect_eve - relative_entropy_vect_bob;
            if DEBUG
                disp('')
                disp('relative_entropy_vect is: ');
                disp(relative_entropy_vect)
            end
            % Then, we compute the chi_2 distance for all X_2 for Eve's channel --------------------------------------------------------
            chi_2_vect = zeros(X2_cardinality,1);
            
            for x2=1:X2_cardinality
                W_Z_X1_0_X2_x2 = W_Z_X1_0_X2(x2,:);
                W_Z_X1_1_X2_x2 = W_Z_X1_1_X2(x2,:);
                chi_2 = InformationTheory.chi_k_distance(W_Z_X1_1_X2_x2, W_Z_X1_0_X2_x2, 2);
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

            rk = sqrt(2)*ratio;

            if (rk < 0)
                % we don't need a secret key if the rate is negative
                % (because Bob's channel is better than Eve's channel in such case)
                rk = 0; 
            end

        end % end covert_sk_rate        

       
        
        % computes   I(X2;Y \mid X1=0, T)
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
                end;
                avg_H_Y_X1_0_X2_T(t)= (P_X2_mid_T(:,t)')* H_Y_X1_0_X2_T(:,t);
                I_T(t)=H_Y_X_0_T(t) - avg_H_Y_X1_0_X2_T(t);
              end;
              
              I=I_T*P_T;
              
        end
    end
end
    