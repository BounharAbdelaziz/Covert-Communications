classdef CovertCommunication

    methods (Static)

        function DEBUG = set_debug(debug)
            CovertCommunication.DEBUG =debug;
        end
        
        function [bool, absolute_continuity_bob, absolute_continuity_eve, different_output_distributions_eve, without_secret_key_condition] = check_theorem_conditions(W_Y_X1_X2, W_Z_X1_X2, P_T, Epsilon_T, P_X2_mid_T, X2_cardinality, Y_cardinality, DEBUG)
            tolerance = 1e-8;
            bool = 0;
            % extracting laws
            W_Y_X1_0_X2 = zeros(2,Y_cardinality);
            W_Y_X1_1_X2 = zeros(2,Y_cardinality);
            W_Z_X1_0_X2 = zeros(2,Y_cardinality);
            W_Z_X1_1_X2 = zeros(2,Y_cardinality);
            
            disp('W_Y_X1_0_X2')
            disp(W_Y_X1_0_X2)

            disp(' W_Y_X1_0_X2(1,:)')
            disp( W_Y_X1_0_X2(1,:))

            disp(' W_Y_X1_X2(1,:)')
            disp( W_Y_X1_X2(1,:))
            

            W_Y_X1_0_X2(1,:) = W_Y_X1_X2(1,:); % W_Y_X1_0_X2_0
            W_Y_X1_0_X2(2,:) = W_Y_X1_X2(2,:); % W_Y_X1_0_X2_1

            W_Y_X1_1_X2(1,:) = W_Y_X1_X2(3,:); % W_Y_X1_1_X2_0
            W_Y_X1_1_X2(2,:) = W_Y_X1_X2(4,:); % W_Y_X1_1_X2_1

            W_Z_X1_0_X2(1,:) = W_Z_X1_X2(1,:); % W_Z_X1_0_X2_0
            W_Z_X1_0_X2(2,:) = W_Z_X1_X2(2,:); % W_Z_X1_0_X2_1

            W_Z_X1_1_X2(1,:) = W_Z_X1_X2(3,:); % W_Z_X1_1_X2_0
            W_Z_X1_1_X2(2,:) = W_Z_X1_X2(4,:); % W_Z_X1_1_X2_1

            % checking conditions
            absolute_continuity_bob = InformationTheory.check_absolute_continuity(W_Y_X1_1_X2, W_Y_X1_0_X2, X2_cardinality, DEBUG);
            absolute_continuity_eve = InformationTheory.check_absolute_continuity(W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, DEBUG);
            different_output_distributions_eve = InformationTheory.check_different_distributions(W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, DEBUG);
            without_secret_key_condition = InformationTheory.without_secret_key_condition(W_Y_X1_0_X2, W_Y_X1_1_X2, W_Z_X1_0_X2, W_Z_X1_1_X2, X2_cardinality);
            
            if (ismembertol(absolute_continuity_bob, 1, tolerance) && ismembertol(absolute_continuity_eve, 1, tolerance) && ismembertol(different_output_distributions_eve, 1, tolerance) && ismembertol(without_secret_key_condition, 1, tolerance))
                bool = 1;
            end
            if DEBUG
                disp(['[INFO] absolute_continuity_bob : ', num2str(absolute_continuity_bob)]);
                disp(['[INFO] absolute_continuity_eve : ', num2str(absolute_continuity_eve)]);
                disp(['[INFO] different_output_distributions_eve : ', num2str(different_output_distributions_eve)]);
                disp(['[INFO] without_secret_key_condition : ', num2str(without_secret_key_condition)]);
                disp(['[INFO] bool : ', num2str(bool)]);
            end
            
        end
        function r1 = covert_message_rate(P_T, P_X2_mid_T, Epsilon_T, W_Y_X1_X2, W_Z_X1_X2, X2_cardinality, Y_cardinality, DEBUG)
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

            W_Y_X1_0_X2(1,:) = W_Y_X1_X2(1,:); % W_Y_X1_0_X2_0
            W_Y_X1_0_X2(2,:) = W_Y_X1_X2(2,:); % W_Y_X1_0_X2_1

            W_Y_X1_1_X2(1,:) = W_Y_X1_X2(3,:); % W_Y_X1_1_X2_0
            W_Y_X1_1_X2(2,:) = W_Y_X1_X2(4,:); % W_Y_X1_1_X2_1

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
            
            W_Z_X1_0_X2(1,:) = W_Z_X1_X2(1,:); % W_Z_X1_0_X2_0
            W_Z_X1_0_X2(2,:) = W_Z_X1_X2(2,:); % W_Z_X1_0_X2_1

            W_Z_X1_1_X2(1,:) = W_Z_X1_X2(3,:); % W_Z_X1_1_X2_0
            W_Z_X1_1_X2(2,:) = W_Z_X1_X2(4,:); % W_Z_X1_1_X2_1

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

        function rk = covert_sk_rate(P_T, P_X2_given_T, epsilons, W_Y_X1_X2, W_Z_X1_X2, DEBUG)
            % Input:
            %   P_T             : a probabilities vector for the time sharing rv
            %   P_X2_given_T    : conditional probability vector P_{X_2 \mid T}
            %   epsilons        : vector of \epsilon_{t} values 
            %   W_Y_X1_X2       : channel law of the legitimate user W_{Y \mid X_1, X_2}
            %   W_Z_X1_X2       : channel law of the adversary W_{Z \mid X_1, X_2}
            % Output:
            %   The secret key rate
            rk = sqrt(2);
        end % end covert_sk_rate

        

        function r2 = non_covert_rate(P_T, P_X2_mid_T, W_Y_X1_0, W_Y_X1_0_X2, DEBUG)
            % Input:
            %   P_T             : a probabilities vector for the time sharing rv
            %   P_X2_mid_T    : conditional probability vector P_{X_2 \mid T} 
            %   P_Y             : a probabilities vector for the output Y rv
            %   W_Y_X1_0_X2     : channel law of the legitimate user W_{Y \mid X_1=0, X_2}
            %   W_Y_X1_0_X2     : channel law of the legitimate user W_{Y \mid X_1=0, X_2}
            % Output:
            %   The rate of the non-covert user U2
            % R_2 = \sum_{(t,x_2)} p(t)p(x_2 \mid t) \mathbb{D}(W_{Y \mid X_1, X_2}(\cdot \mid X_1=0, X_2=x_2) \| W_{Y \mid X_1}(\cdot \mid X_1=0))
            
            relative_entropy_vect = zeros(2,1);
            for x2=1:2
                W_Y_X1_0_X2_x2 = W_Y_X1_0_X2(x2,:);
                relative_entropy = InformationTheory.relative_entropy(W_Y_X1_0_X2_x2, W_Y_X1_0);
                relative_entropy_vect(x2) = relative_entropy;
            end

            r2 =0;
            for t=1:length(P_T)
                average_relative_entropy_fixed_t = dot(P_X2_mid_T(:,t), relative_entropy_vect);
                %disp(P_X2_mid_T(:,t));
                r2 = r2 + P_T(t)*average_relative_entropy_fixed_t;
            end

%             r2 =0;
%             for t=1:length(P_T)
%                 W_Y_X1_0_X2_x2 = W_Y_X1_0_X2(t,:);
%                 relative_entropy = InformationTheory.relative_entropy(W_Y_X1_0_X2_x2, W_Y_X1_0);
%                 r2 = r2 + P_T(t)*P_X2_mid_T(t,:)*relative_entropy;
%             end

%             r2 =0;
%             for t=1:length(P_T)
%                 for x2=1:length(P_X2_mid_T)
%                     W_Y_X1_0_X2_x2 = W_Y_X1_0_X2(x2,:);
%                     
%                     disp('-----------------------------');
%                     disp('W_Y_X1_0_X2_x2 is: ');
%                     disp(W_Y_X1_0_X2_x2);
%                     disp('W_Y_X1_0 is: ');
%                     disp(W_Y_X1_0);
%                     disp('P_T is: ');
%                     disp(P_T);
%                     disp('P_X2_mid_T is: ');
%                     disp(P_X2_mid_T);
%                     relative_entropy = InformationTheory.relative_entropy(W_Y_X1_0_X2_x2, W_Y_X1_0);
%                     disp(['The relative entropy is: ', num2str(relative_entropy)]);
%                     r2 = r2 + P_T(t)*P_X2_mid_T(:,x2)*relative_entropy;
%                 end
%             end
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

        function r2 = old_non_covert_rate(P_T, P_X2, P_Y, W_Y_X1_0, W_Y_X1_0_X2)
            % Input:
            %   P_T             : a probabilities vector for the time sharing rv
            %   P_X2            : a probabilities vector for the input X_2 rv
            %   P_Y             : a probabilities vector for the output Y rv
            %   W_Y_X1_0_X2     : channel law of the legitimate user W_{Y \mid X_1=0, X_2}
            % Output:
            %   The rate of the non-covert user U2
            % R_2 \leq I(X_2;Y \mid X_1=0, T) = \sum_{t} p(t) \left[ - \sum_{x_2} p(x_2) \log(p(x_2 \mid t)) 
            %                                  + \sum_{(x_2,y)} p(x_2)p(y) \log(p(x_2 \mid t)) \right ]
            % By continuity it's also = I(X_2;Y \mid T). We use the later
            % as it's simpler.

            sum_one = -dot(P_Y, log2(W_Y_X1_0));
            %sum_two_y = dot(P_Y,log2(W_Y_X1_0_X2)) ;
            %sum_two = dot(P_X2, sum_two_y);

            sum_two = 0;
            for i=1:length(P_X2)
                for j=1:length(P_Y)
                     sum_two =  sum_two + P_X2(i)*P_Y(j)*log2(W_Y_X1_0_X2(i*j));
                end
            end
            disp(sum_one);
            disp(sum_two);
            
            r2 = dot(P_T, sum_one+sum_two);
        end % end non covert rate

    end % end methods
end