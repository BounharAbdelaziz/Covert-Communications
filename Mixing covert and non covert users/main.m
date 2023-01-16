%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rate Region simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Steps to follow:
% - define a simple channel
% - verify the constaints (absolute continuity)
% - compute the rate
% - do a loop on P_T, P_{X2 \mid T}, P_{X1 \mid T} and \epsilon_T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
% seed for reproducibility
seed                            = 1998;
rng(seed);

% maximum value of P_{X1}(1)
max_P_X1_1                      = 1;

% error in the if conditions (1==1.00)
tolerance                       = 1e-16;

% Prints if debug mode
DEBUG                           = 0;
DEBUG_covert                    = 0;
DEBUG_covert_theorem_contraints = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rate Region simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulations parameters
N_epochs                        = 500;
draw_covert_point               = 0; % covert region (can be ignored)
draw_convhull                   = 1;
plot_3d                         = 1;

% \matchal{T} is of cardinality \leq 4
T_cardinalities                 = [4000];
X2_cardinalities                = [2];
sk_budgets                      = [1];
X1_cardinalities                = 2*ones(length(T_cardinalities)); % always 2.
max_epsilon_t_s                 = ones(length(T_cardinalities));
optimize_epsilons_T_s           = ones(length(T_cardinalities));
fixed_channel_laws_vect         = ones(length(T_cardinalities)); % or zeros
uniform_P_T_s                   = zeros(length(T_cardinalities));
compute_marginal_PY             = 0;

%(need to be done on the same seed in order to compare) This will swap the channel of eve and bob, so that the difference in relative entropies becoms negative.
swap_bob_eve_channels           = 0;

% make sure we don't have a missing value
assert(length(T_cardinalities) == length(X2_cardinalities))
assert(length(T_cardinalities) == length(sk_budgets))
assert(length(T_cardinalities) == length(max_epsilon_t_s))
assert(length(T_cardinalities) == length(optimize_epsilons_T_s))
assert(length(T_cardinalities) == length(fixed_channel_laws_vect))

% we add 3 because we add manually missing points:
%     (r1=0,r2=0) (r1=0,max r2) (max r1, r2=0)
N = N_epochs + 3;
% The rates are stored in a matrix, each row is for one experiment
r1_vects = zeros(N, length(T_cardinalities));  % length(T_cardinalities) is the number of experiments
r2_vects = zeros(N, length(T_cardinalities));
rk_vects = zeros(N, length(T_cardinalities));
higher_than_sk_budgets = zeros(N_epochs, length(T_cardinalities));

% counters
cpt = zeros(length(T_cardinalities));
cpt_higher_than_max_P_X1_1 = zeros(length(T_cardinalities));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Running all experiments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dummy variables as we generate a channel law in the first experiment and
% force the same to be used in the second one. The channel laws
% will be generated in run_simulation(). We do them first so that if it's
% not fixed, then the new values outputed by the first experiment will be 
% used for the other ones.

W_Y_X1_X2 = 0;
W_Z_X1_X2 = 0;
for experiment=1:length(T_cardinalities)

    T_cardinality                   = T_cardinalities(experiment);
    uniform_P_T                     = uniform_P_T_s(experiment);
    % channel parameters
    X1_cardinality                  = X1_cardinalities(experiment);
    X2_cardinality                  = X2_cardinalities(experiment);
    X1_X2_cardinality               = X1_cardinality*X2_cardinality; % cartesian product
    Y_cardinality                   = X1_X2_cardinality;
    % simulations parameters
    max_epsilon_t                   = max_epsilon_t_s(experiment);
    optimize_epsilons_T             = optimize_epsilons_T_s(experiment);
    fixed_channel_laws              = fixed_channel_laws_vect(experiment);
    sk_budget                       = sk_budgets(experiment);
    max_P_X1_1                      = 1;
    % channel law initialization
    if (fixed_channel_laws)
        % Fixing the channel law instead of generating it
%         W_Y_X1_X2 = [0.2, 0.3, 0.2, 0.3; 0.1, 0.2, 0.3, 0.4;
%                     0.23,0.46, 0.12, 0.19; 0.33,0.26, 0.22, 0.19];   %% first X2_cardinality rows for x1=0 and latter for  x1=1 
%         W_Z_X1_X2 = [0.3, 0.2, 0.1, 0.4; 0.3, 0.2, 0.15, 0.35;
%                     0.33,0.15, 0.23, 0.29; 0.23,0.26, 0.22, 0.29];  %% first X2_cardinality rows for x1=0 and latter for  x1=1 
        % For this channel, the diverence in divergence is always positive.
        W_Y_X1_X2 = [0.2, 0.3, 0.2, 0.3; 0.2, 0.2, 0.3, 0.3;
                    0.23,0.26, 0.22, 0.29; 0.23,0.26, 0.22, 0.29];   %% first X2_cardinality rows for x1=0 and latter for  x1=1 
        W_Z_X1_X2 = [0.3, 0.2, 0.1, 0.4; 0.3, 0.2, 0.15, 0.35;
                    0.43,0.05, 0.33, 0.19; 0.23,0.16, 0.42, 0.19];  %% first X2_cardinality rows for x1=0 and latter for  x1=1 
    end
    
    % running the experiment
    [r1_vect_i, r2_vect_i, rk_vect_i, higher_than_sk_budgets_i, cpt_i, cpt_higher_than_max_P_X1_1_i, W_Y_X1_X2, W_Z_X1_X2] = run_simulation( ...
                                            uniform_P_T, W_Y_X1_X2, W_Z_X1_X2, fixed_channel_laws,...
                                            T_cardinality, X1_cardinality, X2_cardinality, Y_cardinality, X1_X2_cardinality, ...
                                            sk_budget, swap_bob_eve_channels, max_epsilon_t, optimize_epsilons_T, ...
                                            max_P_X1_1, compute_marginal_PY, ...
                                            DEBUG, DEBUG_covert, DEBUG_covert_theorem_contraints, tolerance, seed, ...
                                            N_epochs, draw_covert_point);
    
    % storing the rates
    r1_vects(:, experiment) = r1_vect_i;
    r2_vects(:, experiment) = r2_vect_i;
    rk_vects(:, experiment) = rk_vect_i;
    higher_than_sk_budgets(:, experiment) = higher_than_sk_budgets_i;

    min_rk_higher_than_sk_budgets = min(higher_than_sk_budgets(:, experiment));
    max_rk_higher_than_sk_budgets = max(higher_than_sk_budgets(:, experiment));
    disp('--------------------------------------')    
    disp(['The secret key rate budget for covert communication is: ', num2str(sk_budget)]);
    disp(['The minimum secret key rate required for covert communication is: ', num2str(min_rk_higher_than_sk_budgets)]);
    disp(['The maximum secret key rate required for covert communication is: ', num2str(max_rk_higher_than_sk_budgets)]);
    disp('--------------------------------------')  
    
    % storing the counters
    cpt(experiment) = cpt_i;
    cpt_higher_than_max_P_X1_1(experiment) = cpt_higher_than_max_P_X1_1_i;
end

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

% if (swap_bob_eve_channels)
%     need_sk = 'when we need a Secret-Key';
% else
%     need_sk = 'when we do not need a Secret-Key';
% end

% title_plot = ['Rate region ', need_sk, ' [seed=', num2str(seed), ']'];
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