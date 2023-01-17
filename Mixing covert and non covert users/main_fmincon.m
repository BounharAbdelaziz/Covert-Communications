%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rate Region simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Steps to follow:
% - define a simple channel
% - verify the constaints (absolute continuity)
% - compute the rate
% - do a loop on P_T, P_{X2 \mid T}, P_{X1 \mid T} and \epsilon_T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
seed    = 1998;

% seed for reproducibility
rng(seed);

% simulation parameters
step_size = 0.5;
lb_mu_1 = 0;
ub_mu_1 = 1;
lb_mu_2 = 0;
ub_mu_2 = 1;
total_runs = ((ub_mu_1-lb_mu_1)/step_size)*((ub_mu_2-lb_mu_2)/step_size);

% \matchal{T} is of cardinality \leq 4
T_cardinalities               = [2,2,2,2];
X2_cardinalities              = [2,2,2,2];
sk_budgets                    = [0.2,0.4,0.5,1];
X1_cardinalities              = 2*ones(length(T_cardinalities)); % always 2.

% ploting parameters
plot_3d                       = 1;
draw_convhull                 = 1;

% Prints if debug mode
DEBUG                           = 0;
DEBUG_covert                    = 0;
DEBUG_covert_theorem_contraints = 0;

% make sure we don't have a missing value
assert(length(T_cardinalities) == length(X2_cardinalities))
assert(length(T_cardinalities) == length(sk_budgets))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rate Region simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% non linear constraints
nonlcon = @(probas_and_eps) my_rate_constraints(probas_and_eps);

% optimization options
options = optimoptions(@fmincon,'StepTolerance',1e-15,'FunctionTolerance',1e-15,'OptimalityTolerance',1e-15,'MaxFunctionEvaluations',1e+10, 'MaxIterations', 1e3, 'Display','off');

% initial point (r1,r2,rk) (should it be the distributions instead?)
probas_and_eps_guess = rand(8,1);% [P_X2_mid_T(0 | 0), P_X2_mid_T(0 | 1), P_X2_mid_T(1 | 0), P_X2_mid_T(1 | 1), P_T(0), P_T(1), Epsilon_T(0), Epsilon_T(1)]

% The optimization is over probabilities so they must be in [0,1]
lb = zeros(length(probas_and_eps_guess),1);
ub = ones(length(probas_and_eps_guess),1);

% we add 3 because we add manually missing points:
%     (r1=0,r2=0) (r1=0,max r2) (max r1, r2=0)
N = total_runs +3;

% The rates are stored in a matrix, each row is for one experiment
r1_vects = zeros(N, length(T_cardinalities)); % length(T_cardinalities) is the number of experiments
r2_vects = zeros(N, length(T_cardinalities));
rk_vects = zeros(N, length(T_cardinalities));

% For each experiment from the list of comparative experiments.
for experiment=1:length(T_cardinalities)

    % for now only one experiment, generalization to be done after.
    T_cardinality           = T_cardinalities(experiment);
    X2_cardinality          = X2_cardinalities(experiment);
    sk_budget               = sk_budgets(experiment);
    X1_cardinality          = X1_cardinalities(experiment);    
    X1_X2_cardinality       = X1_cardinality*X2_cardinality; % cartesian product
    Y_cardinality           = X1_X2_cardinality;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fix the channel law %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % For this channel, the diverence in divergence is always positive.
    W_Y_X1_X2 = [0.2, 0.3, 0.2, 0.3; 0.2, 0.2, 0.3, 0.3;
                0.23,0.26, 0.22, 0.29; 0.23,0.26, 0.22, 0.29];   %% first X2_cardinality rows for x1=0 and latter for  x1=1 
    W_Z_X1_X2 = [0.3, 0.2, 0.1, 0.4; 0.3, 0.2, 0.15, 0.35;
                0.43,0.05, 0.33, 0.19; 0.23,0.16, 0.42, 0.19];  %% first X2_cardinality rows for x1=0 and latter for  x1=1 
    
    % Extract specific channel laws
    [W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2]    = CovertCommunication.extract_laws(W_Y_X1_X2, W_Z_X1_X2, X2_cardinality,Y_cardinality);


    % The loss is \mu_1*r_1 + \mu_2*r_2. We vary \mu_1 and \mu_2 in order
    % to get all points in the boundary. If \mu_1 = 0, we get the maximum
    % rate for r_2 and vice versa. Otherwise, we get the points in between.
    
    index_rates_experiment = 1; % where to store the rates we compute

    for mu_1=progress(lb_mu_1:step_size:ub_mu_1)

        for mu_2=lb_mu_2:step_size:ub_mu_2
            
            % optimization function
            function_to_maximize = @(probas_and_eps) objective_function_fmincon(probas_and_eps(1), probas_and_eps(2), probas_and_eps(3), probas_and_eps(4), probas_and_eps(5), probas_and_eps(6), probas_and_eps(7), probas_and_eps(8), ...
                                                                   W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, Y_cardinality, X1_X2_cardinality, sk_budget, mu_1, mu_2, DEBUG_covert);
            
            % constraints
            nonlcon = @(probas_and_eps) my_rate_constraints(probas_and_eps, sk_budget, W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, DEBUG_covert);

            % choose a new random point
            probas_and_eps_guess = zeros(length(probas_and_eps_guess), 1);
            
            % optimize probas_and_eps_guess
            [b_min,fval] = fmincon(function_to_maximize,probas_and_eps_guess,[],[],[],[],lb,ub,nonlcon,options);
            
%             disp(b_min)
            % compute the rates with the optimal probas_and_eps_guess for
            % the fixed mu_1 and mu_2.

            [r1, r2, rk] = compute_rates_for_fmincon(b_min(1), b_min(2), b_min(3), b_min(4), b_min(5), b_min(6), b_min(7), b_min(8), ...
                                                        W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, Y_cardinality, X1_X2_cardinality, sk_budget, DEBUG_covert);

            r1_vects(index_rates_experiment, experiment) = r1;
            r2_vects(index_rates_experiment, experiment) = r2;
            rk_vects(index_rates_experiment, experiment) = rk;
%             disp(r1);
%             disp(r2);
%             disp(rk);
            index_rates_experiment = index_rates_experiment +1;
            [c,ceq] = my_rate_constraints(b_min, sk_budget, W_Y_X1_1_X2, W_Y_X1_0_X2, W_Z_X1_1_X2, W_Z_X1_0_X2, X2_cardinality, DEBUG_covert);
%             disp(c);
%             disp(ceq);
%             disp('--------------------------')
        end
    end
    min_r2 = min(r2_vects(:, experiment));
    max_r2 = max(r2_vects(:, experiment));
    min_r1 = min(r1_vects(:, experiment));
    max_r1 = max(r1_vects(:, experiment));
    
    % add (r1=0,r2=0) 
    r1_vects(total_runs+1, experiment)=0;
    r2_vects(total_runs+1, experiment)=0;
    % add (r1=0,max r2) 
    r1_vects(total_runs+2, experiment)=0;
    r2_vects(total_runs+2, experiment)=max_r2;
    % add  (max r1, r2=0)
    r1_vects(total_runs+3, experiment)=max_r1;
    r2_vects(total_runs+3, experiment)=0;
end


% rates=-1*fval;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Legends %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To know which size for legends_convhull cell.
sum_non_zero_all = 0;
for experiment=1:length(T_cardinalities)
    check_sum = r1_vects(N, experiment) + r2_vects(N-1, experiment); % max_r1 + max_r2
    if (check_sum >0)
        sum_non_zero_all = sum_non_zero_all + 1;
    end
end
% for legends
legends = cell(length(T_cardinalities), 1);
legends_convhull = cell(sum_non_zero_all, 1);

for experiment=1:length(T_cardinalities)
    legends(experiment) = {['X2 cardinality: ', num2str(X2_cardinalities(experiment)), ' // T cardinality: ', num2str(T_cardinalities(experiment)), ' // Secret-Key budget: ',  num2str(sk_budgets(experiment))]};
end

% we need to go through all of them but the indexes should be \leq
% sum_non_zero_all which is < length(T_cardinalities).
index_ = 1;
for experiment=1:length(T_cardinalities)
    check_sum = r1_vects(N, experiment) + r2_vects(N-1, experiment); % max_r1 + max_r2
    if (check_sum >0)
        legends_convhull(index_) = {['X2 cardinality: ', num2str(X2_cardinalities(experiment)), ' // T cardinality: ', num2str(T_cardinalities(experiment)), ' // Secret-Key budget: ',  num2str(sk_budgets(experiment))]};
        index_ = index_ + 1;
    end
end

title_plot = ['Rate region [seed=', num2str(seed), ']'];
x_title = 'Covert user square-root rate';
y_title = 'Non-Covert user rate';
z_title = 'Secret-Key square-root rate';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3D plot of the region (secret key rate, covert rate, non-covert rate)
if (plot_3d)
    figure 
    ax1 = nexttile;
    for experiment=1:length(T_cardinalities)
        scatter3(ax1, r1_vects(:, experiment), r2_vects(:, experiment), rk_vects(:, experiment));
        hold on
    end
    title(ax1,title_plot);
    xlabel(ax1,x_title);
    ylabel(ax1,y_title);
    zlabel(ax1,z_title);
    legend(legends);
end

% covert rate vs non-covert rate
ax2 = nexttile;
for experiment=1:length(T_cardinalities)
    scatter(ax2, r1_vects(:, experiment), r2_vects(:, experiment));
    hold on
end
title(ax2, title_plot);
xlabel(ax2,x_title);
ylabel(ax2,y_title);
legend(legends);


% secret key rate vs covert rate
ax3 = nexttile;
for experiment=1:length(T_cardinalities)
    scatter(ax3, rk_vects(:, experiment), r1_vects(:, experiment));
    hold on
end
title(ax3, title_plot);
xlabel(ax3,z_title);
ylabel(ax3,x_title);
legend(legends);

% secret key rate vs non-covert rate
ax4 = nexttile;
for experiment=1:length(T_cardinalities)
    scatter(ax4, rk_vects(:, experiment), r2_vects(:, experiment));
    hold on
end
title(ax4, title_plot);
xlabel(ax4,z_title);
ylabel(ax4,y_title);
legend(legends);

% convhull of the rate region (covert and non-covert rates)
if draw_convhull
    figure 
    ax5 = nexttile;
    hold on
    for experiment=1:length(T_cardinalities)
        % if all are zeros then the max is zero, so the convexhull doesn't
        % exist
        check_sum = r1_vects(N, experiment) + r2_vects(N-1, experiment); % max_r1 + max_r2
        if (check_sum >0)
            P_i = [r1_vects(:, experiment) , r2_vects(:, experiment)];
            [k_i,av_i] = convhull(P_i);
            plot(ax5, P_i(k_i,1),P_i(k_i,2))
        else
            disp(['[INFO] Experiment ', num2str(experiment), ' has only one point (0,0,0). No convex hull can be obtained!']);
            continue
        end
        hold on
    end
    title(ax5, title_plot);
    xlabel(ax5,x_title);
    ylabel(ax5,y_title);
    legend(legends_convhull);
end