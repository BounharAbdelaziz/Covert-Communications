classdef InformationTheory
    methods (Static)

        function DEBUG = set_debug(debug)
            DEBUG =debug;
        end

        function bool = valid_distribution(P)
            bool = assert(sum(P) == 1);
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

        function without_secret_key = without_secret_key_condition(W_Y_X1_0_X2, W_Y_X1_1_X2, W_Z_X1_0_X2, W_Z_X1_1_X2, X2_cardinality) % what we need to do instead(P_T, Epsilon_T, P_X2_mid_T, W_Y_X1_0_X2, W_Y_X1_1_X2, W_Z_X1_0_X2, W_Z_X1_1_X2, X2_cardinality)
            % Input:
            %   W_Y_X1_0/1_X2       : channel law of the legitimate user W_{Y \mid X_1=0/1, X_2}
            %   W_Z_X1_0/1_X2       : channel law of the adversary W_{Z \mid X_1=0/1, X_2}
            %   P_X2_mid_T      : conditional probability vector P_{X_2 \mid T}
            % Output:
            %   Checks if D(W_Y_X1_1_X2 \mid W_Y_X1_0_X2) >= D(W_Z_X1_1_X2 \mid W_Z_X1_0_X2) on average over P_T \epsilon_T (for now
            %   checks only if the two divergences (for each x2) are positive as we cannot do the average)

            relative_entropy_bob_vect = zeros(2,1);
            relative_entropy_eve_vect = zeros(2,1);
            
            for x2=1:X2_cardinality

                W_Y_X1_0_X2_x2 = W_Y_X1_0_X2(x2,:);
                W_Y_X1_1_X2_x2 = W_Y_X1_1_X2(x2,:);
                W_Z_X1_0_X2_x2 = W_Z_X1_0_X2(x2,:);
                W_Z_X1_1_X2_x2 = W_Z_X1_1_X2(x2,:);

                relative_entropy_bob = InformationTheory.relative_entropy(W_Y_X1_1_X2_x2, W_Y_X1_0_X2_x2);
                relative_entropy_eve = InformationTheory.relative_entropy(W_Z_X1_0_X2_x2, W_Z_X1_1_X2_x2);

                relative_entropy_bob_vect(x2) = relative_entropy_bob;
                relative_entropy_eve_vect(x2) = relative_entropy_eve;

            end

            difference_relative_entropy_vect = relative_entropy_eve_vect-relative_entropy_bob_vect;
            disp('difference_relative_entropy_vect')
            disp(difference_relative_entropy_vect)
            sum_signs = sign(difference_relative_entropy_vect(1)) + sign(difference_relative_entropy_vect(2)); 
            % if both are positive the sum is 2
            if (ismembertol(sum_signs, 2, 1e-8)) 
                without_secret_key = 1;
            else
                without_secret_key = 0;
            end
            
            % this is what I think should be done but we varry P_X2_mid_T
            % and P_T so I don't know how we can check this condition in
            % our setup. Bloch didn't have this averaging. I am taking a
            % simpler thing, to ensure 100% it's feasible, both should be
            % positive so that the average is two. When one of them is
            % negative, we cannot conclude. So basically, we don't take the
            % risk and try again untill we find a channel law where both
            % are positive

% % %             avg = 0;
% % %             for t=1:length(P_T)
% % %                 % computing the averages of the relative entropy and chi_2 for a fixed t
% % %                 average_relative_entropy_fixed_t = dot(P_X2_mid_T(:,t), difference_relative_entropy_vect);
% % %                 
% % %                 % taking the average over T
% % %                 avg = avg + P_T(t)*Epsilon_T(t)*average_relative_entropy_fixed_t;
% % %             end
% % % 
% % %             if (avg < 0)
% % %                 without_secret_key = 1;
% % %             else
% % %                 without_secret_key = 0;
% % %             end

            
        end %end check without secret key condition
        function H = entropy(P)
            % Input:
            %   P: a probabilities vector
            % Output:
            %   The entropy H(X), X~P
            H = -dot(P,log2(P));
        end % end entropy
        
        function d = relative_entropy(P, Q)
            % Input:
            %   P: a probabilities vector
            %   Q: a probabilities vector
            % Output:
            %   The relative entropy (KL div) D(P||Q)
            assert(numel(P) == numel(Q));
            
            d = 0;
            for i=1:length(P)
%                 if (Q(i) == 0)
%                     Q(i) = 1e-3; % for numerical stability
%                 end
                log_ratio = log2(P(i)/Q(i));
%                 if (log_ratio == -inf)
%                     log_ratio =0;
%                 end
                d = d + P(i)*log_ratio;
%                 disp('-----------------------------');
%                 disp('d is')
%                 disp(d);
%                 disp('P(i) is')
%                 disp(P(i));
%                 disp('Q(i) is')
%                 disp(Q(i));
%                 disp('log2(P(i)/Q(i)) is')
%                 disp(log2(P(i)/Q(i)));
            end
%             disp('end d is')
%             disp(d);
                        
        end % end relative_entropy

        function chi = chi_k_distance(P, Q, k)
            % Input:
            %   P: a probabilities vector
            %   Q: a probabilities vector
            %   k: k>=1 is the order of the distance. When k=2, it's the Chi square test
            % Output:
            %   The Chi k distance \Chi_{k}(P||Q)
            assert(numel(P) == numel(Q));
            
            chi = 0;
            for i=1:length(P)
                if (P(i) == 0)
                    P(i) = 1e-3; % for numerical stability
                end
                difference_power_k = (P(i)-Q(i)).^k;
                P_i_power_k_1 = (P(i)).^(k-1);
                ratio = difference_power_k/P_i_power_k_1;
                chi = chi + ratio;
%                 disp('-----------------------------');
%                 disp('chi is')
%                 disp(chi);
%                 disp('P(i) is')
%                 disp(P(i));
%                 disp('Q(i) is')
%                 disp(Q(i));
%                 disp('ratio is')
%                 disp(ratio);
            end
%             disp('end chi is')
%             disp(chi);
                        
        end % end relative_entropy

        

        function d = old_relative_entropy(P, Q)
            % Input:
            %   P: a probabilities vector
            %   Q: a probabilities vector
            % Output:
            %   The relative entropy (KL div) D(P||Q)
            assert(numel(P) == numel(Q));

            d= dot(P, log2(P)-log2(Q));
            entropy_P = InformationTheory.entropy(P);
            sum_p_log_q = dot(P, log2(Q));
            d2 = -entropy_P - sum_p_log_q;
            disp('d2 is')
            disp(d2);
                        
        end % end relative_entropy

        function CH = conditional_entropy(x)

            CH=x;
        end

        function expec = expectation(P, X)
            % Input:
            %   P: a probabilities vector
            %   X: realizations of X~P
            % Output:
            %   The expectation of X, X~P

            assert(numel(P) == numel(X));

            n=length(P);

            expec=0;

            for i=1:n
                expec = expec + P(i) * X(i);
            end
        end % end expectation

        function W = find_mac_channel_law_from_p2p_channel_law(P_X, W_Y_X)

            % get input distribution pobabilities
            P0 = P_X(1);
            P1 = P_X(2);
            % get individual probabilities from channel matrix
            W_0_0 = W_Y_X(1,1);
            W_0_1 = W_Y_X(1,2);
            W_1_0 = W_Y_X(2,1);
            W_1_1 = W_Y_X(2,2);

            syms W_0_0_0 W_0_1_0 W_0_0_1 W_0_1_1 W_1_0_0 W_1_1_0 W_1_0_1 W_1_1_1
            
            eqn1 = P0*W_0_0_0 + P1*W_0_1_0 == W_0_0;
            eqn2 = P0*W_0_0_1 + P1*W_0_1_1 == W_0_1;
            eqn3 = P0*W_1_0_0 + P1*W_1_1_0 == W_1_0;
            eqn4 = P0*W_1_0_1 + P1*W_1_1_1 == W_1_1;

            eqn5 = W_0_0_0 + W_1_0_0 ==1;
            eqn6 = W_0_1_0 + W_1_1_0 ==1;
            eqn7 = W_0_0_1 + W_1_0_1 ==1;
            eqn8 = W_0_1_1 + W_1_1_1 ==1;

            % transform the equations to matrix
            [A,B] = equationsToMatrix( ...
                    [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8], ...
                    [W_0_0_0, W_0_1_0, W_0_0_1, W_0_1_1, W_1_0_0, W_1_1_0, W_1_0_1, W_1_1_1] ...
                   );
            disp(A);
            disp(B);

            % solve the linear system to get the individual probabilities
            W = linsolve(A,B);
        end

        function W_Y_X1 = find_marginal_channel_law_from_mac_channel_law(P_X2, W_Y_X1_X2)
            W_0_0_0 = W_Y_X1_X2(1);
            W_0_1_0 = W_Y_X1_X2(2); 
            W_0_0_1 = W_Y_X1_X2(3); 
            W_0_1_1 = W_Y_X1_X2(4);
            W_1_0_0 = W_Y_X1_X2(5);
            W_1_1_0 = W_Y_X1_X2(6);
            W_1_0_1 = W_Y_X1_X2(7);
            W_1_1_1 = W_Y_X1_X2(8);

            W_Y_0_X1_0 = dot(P_X2, [W_0_0_0, W_0_1_0]);
            W_Y_0_X1_1 = dot(P_X2, [W_0_0_1, W_0_1_1]);
            W_Y_1_X1_0 = dot(P_X2, [W_1_0_0, W_1_1_0]);
            W_Y_1_X1_1 = dot(P_X2, [W_1_0_1, W_1_1_1]);

            W_Y_X1      = zeros(2,2);
            W_Y_X1(1,1) = W_Y_0_X1_0;
            W_Y_X1(1,2) = W_Y_0_X1_1;
            W_Y_X1(2,1) = W_Y_1_X1_0;
            W_Y_X1(2,2) = W_Y_1_X1_1;
            
        end

        function P_Y = find_marginal_output_law_from_p2p_channel_law(P_X1, W_Y_X1)
            
            W_Y_0_X1_0 = W_Y_X1(1,1);
            W_Y_0_X1_1 = W_Y_X1(1,2);
            W_Y_1_X1_0 = W_Y_X1(2,1);
            W_Y_1_X1_1 = W_Y_X1(2,2);

            P_Y_0= dot(P_X1, [W_Y_0_X1_0, W_Y_0_X1_1]);
            P_Y_1 = dot(P_X1, [W_Y_1_X1_0, W_Y_1_X1_1]);

            P_Y = [P_Y_0, P_Y_1];
            
        end

        function [x,v] = generate_probability_vector(n,m,s,a,b)
        % adapted from https://fr.mathworks.com/matlabcentral/fileexchange/9700-random-vectors-with-fixed-sum
        % [x,v] = randfixedsum(n,m,s,a,b)
        %
        %   This generates an n by m array x, each of whose m columns
        % contains n random values lying in the interval [a,b], but
        % subject to the condition that their sum be equal to s.  The
        % scalar value s must accordingly satisfy n*a <= s <= n*b.  The
        % distribution of values is uniform in the sense that it has the
        % conditional probability distribution of a uniform distribution
        % over the whole n-cube, given that the sum of the x's is s.
        %
        %   The scalar v, if requested, returns with the total
        % n-1 dimensional volume (content) of the subset satisfying
        % this condition.  Consequently if v, considered as a function
        % of s and divided by sqrt(n), is integrated with respect to s
        % from s = a to s = b, the result would necessarily be the
        % n-dimensional volume of the whole cube, namely (b-a)^n.
        %
        %   This algorithm does no "rejecting" on the sets of x's it
        % obtains.  It is designed to generate only those that satisfy all
        % the above conditions and to do so with a uniform distribution.
        % It accomplishes this by decomposing the space of all possible x
        % sets (columns) into n-1 dimensional simplexes.  (Line segments,
        % triangles, and tetrahedra, are one-, two-, and three-dimensional
        % examples of simplexes, respectively.)  It makes use of three
        % different sets of 'rand' variables, one to locate values
        % uniformly within each type of simplex, another to randomly
        % select representatives of each different type of simplex in
        % proportion to their volume, and a third to perform random
        % permutations to provide an even distribution of simplex choices
        % among like types.  For example, with n equal to 3 and s set at,
        % say, 40% of the way from a towards b, there will be 2 different
        % types of simplex, in this case triangles, each with its own
        % area, and 6 different versions of each from permutations, for
        % a total of 12 triangles, and these all fit together to form a
        % particular planar non-regular hexagon in 3 dimensions, with v
        % returned set equal to the hexagon's area.
        %
        % Roger Stafford - Jan. 19, 2006
        % Check the arguments.
        if (m~=round(m))|(n~=round(n))|(m<0)|(n<1)
         error('n must be a whole number and m a non-negative integer.')
        elseif (s<n*a)|(s>n*b)|(a>=b)
         error('Inequalities n*a <= s <= n*b and a < b must hold.')
        end
        % Rescale to a unit cube: 0 <= x(i) <= 1
        s = (s-n*a)/(b-a);
        % Construct the transition probability table, t.
        % t(i,j) will be utilized only in the region where j <= i + 1.
        k = max(min(floor(s),n-1),0); % Must have 0 <= k <= n-1
        s = max(min(s,k+1),k); % Must have k <= s <= k+1
        s1 = s - [k:-1:k-n+1]; % s1 & s2 will never be negative
        s2 = [k+n:-1:k+1] - s;
        w = zeros(n,n+1); w(1,2) = realmax; % Scale for full 'double' range
        t = zeros(n-1,n);
        tiny = 2^(-1074); % The smallest positive matlab 'double' no.
        for i = 2:n
         tmp1 = w(i-1,2:i+1).*s1(1:i)/i;
         tmp2 = w(i-1,1:i).*s2(n-i+1:n)/i;
         w(i,2:i+1) = tmp1 + tmp2;
         tmp3 = w(i,2:i+1) + tiny; % In case tmp1 & tmp2 are both 0,
         tmp4 = (s2(n-i+1:n) > s1(1:i)); % then t is 0 on left & 1 on right
         t(i-1,1:i) = (tmp2./tmp3).*tmp4 + (1-tmp1./tmp3).*(~tmp4);
        end
        % Derive the polytope volume v from the appropriate
        % element in the bottom row of w.
        v = n^(3/2)*(w(n,k+2)/realmax)*(b-a)^(n-1);
        % Now compute the matrix x.
        x = zeros(n,m);
        if m == 0, return, end % If m is zero, quit with x = []
        rt = rand(n-1,m); % For random selection of simplex type
        rs = rand(n-1,m); % For random location within a simplex
        s = repmat(s,1,m);
        j = repmat(k+1,1,m); % For indexing in the t table
        sm = zeros(1,m); pr = ones(1,m); % Start with sum zero & product 1
        for i = n-1:-1:1  % Work backwards in the t table
         e = (rt(n-i,:)<=t(i,j)); % Use rt to choose a transition
         sx = rs(n-i,:).^(1/i); % Use rs to compute next simplex coord.
         sm = sm + (1-sx).*pr.*s/(i+1); % Update sum
         pr = sx.*pr; % Update product
         x(n-i,:) = sm + pr.*e; % Calculate x using simplex coords.
         s = s - e; j = j - e; % Transition adjustment
        end
        x(n,:) = sm + pr.*s; % Compute the last x
        % Randomly permute the order in the columns of x and rescale.
        rp = rand(n,m); % Use rp to carry out a matrix 'randperm'
        [ig,p] = sort(rp); % The values placed in ig are ignored
        x = (b-a)*x(p+repmat([0:n:n*(m-1)],n,1))+a; % Permute & rescale x
        end


        function A = generate_random_channel_matrix(m, n, smallest_poba_value)
            % adapted from https://fr.mathworks.com/matlabcentral/answers/476398-how-to-generate-matrices-that-satisfies-constraints-on-sum-of-row-elements-and-sum-of-column-element

            rs = ones(m,1); % target row sum
            cs = ones(n,1); % target column sum

            if abs(sum(rs)-sum(cs)) > 1e-10*abs(sum(rs))
                error('sum(rs) must be equal to sum(cs)'); 
            end
            % This part need to be done once if rs and cs are unchanged
            I = zeros(m,n);
            I(:) = 1:m*n;
            Aeq = accumarray([repmat((1:m)',n,1) I(:); repelem(m+(1:n)',m,1) I(:)],1);
            beq = [rs(:); cs(:)];
            lb = ones(m*n,1)*smallest_poba_value;
            ub = inf(m*n,1);
            B = zeros(m*n);
            for k=1:m*n
                f = accumarray(k, 1, [m*n 1]);
                B(:,k) = linprog(f, [], [], Aeq, beq, lb, ub);   
            end
            % This part generate new random A 
            x=-log(rand(m*n,1)); x=x./sum(x);
            % x = randfixedsum(m*n,1,1,0,1); % Fex
            A = reshape(B*x,[m n]);
            % Check
%             sum(A,2) % ~ rs
%             sum(A,1) % ~ cs
        end

    end % end methods
end