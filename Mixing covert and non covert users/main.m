%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rate Region simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Steps to follow:
% - define a simple channel
% - verify the constaints (absolute continuity)
% - compute the rate
% - do a loop on P_T, P_{X2 \mid T}, P_{X1 \mid T} and \epsilon_T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
% seed for reproducibility
seed                            = 10;
rng(seed);

sk_budget_1 = 0;
sk_budget_2 = 5;

% simulations parameters
N_epochs                        = 50;

% covert region
draw_covert_point               = 0;
draw_convhull                   = 0;
plot_3d                         = 1;

% maximum value of P_{X1}(1)
max_P_X1_1                      = 1;

% error in the if conditions (1==1.00)
tolerance                       = 1e-16;

% Prints if debug mode
DEBUG                           = 0;
DEBUG_covert                    = 0;
DEBUG_covert_theorem_contraints = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rate Region simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for simulation
% \matchal{T} is of cardinality 4
T_cardinality                   = 4;
uniform_P_T                     = 0;
%(need to be done on the same seed in order to compare) This will swap the channel of eve and bob, so that the difference in relative entropies becoms negative.
swap_channels_for_need_sk       = 0; 

% channel parameters
X1_cardinality                  = 2;
X2_cardinality                  = 2;
X1_X2_cardinality               = X1_cardinality*X2_cardinality; % cartesian product
Y_cardinality                   = X1_X2_cardinality;
compute_marginal_PY             = 0;

% simulations parameters
max_epsilon_t                   = 1;
optimize_epsilons_T             = 1;
fixed_channel_laws              = 0;

% dummy variables as we generate a channel law in the first experiment and
% force the same to be used in the second one.
W_Y_X1_X2 = 0;
W_Z_X1_X2 = 0;

[r1_vect_1, r2_vect_1, rk_vect_1, cpt_1, cpt_higher_than_max_P_X1_1_1, W_Y_X1_X2, W_Z_X1_X2] = run_simulation( ...
                                            uniform_P_T, W_Y_X1_X2, W_Z_X1_X2, fixed_channel_laws,...
                                            T_cardinality, X1_cardinality, X2_cardinality, Y_cardinality, X1_X2_cardinality, ...
                                            swap_channels_for_need_sk, sk_budget_1, max_epsilon_t, optimize_epsilons_T, ...
                                            max_P_X1_1, compute_marginal_PY, ...
                                            DEBUG, DEBUG_covert, DEBUG_covert_theorem_contraints, tolerance, seed, ...
                                            N_epochs, draw_covert_point);

% for the comparative simulation
T_cardinality_2                   = 4;
uniform_P_T_2                     = 0;
%(need to be done on the same seed in order to compare) This will swap the channel of eve and bob, so that the difference in relative entropies becoms negative.
swap_channels_for_need_sk_2       = 0;

% channel parameters
X1_cardinality_2                  = 2;
X2_cardinality_2                  = 2;
X1_X2_cardinality_2               = X1_cardinality_2*X2_cardinality_2; % cartesian product
Y_cardinality_2                   = X1_X2_cardinality_2;
compute_marginal_PY_2             = 0;

% simulations parameters
max_epsilon_t_2                   = 1;
optimize_epsilons_T_2             = 1;
fixed_channel_laws                = 1;

% if the cardinalities are not the same, we cannot enforce them to have the
% same channel laws.
if (~ismembertol(sk_budget_1, sk_budget_2, tolerance))
    fixed_channel_laws = 0;
end
if (~ismembertol(X2_cardinality_2, X2_cardinality, tolerance))
    fixed_channel_laws = 0;
end
[r1_vect_2, r2_vect_2, rk_vect_2, cpt_2, cpt_higher_than_max_P_X1_1_2, W_Y_X1_X2, W_Z_X1_X2] = run_simulation( ...
                                            uniform_P_T, W_Y_X1_X2, W_Z_X1_X2, fixed_channel_laws, ...
                                            T_cardinality_2, X1_cardinality_2, X2_cardinality_2, Y_cardinality_2, X1_X2_cardinality_2, ...
                                            swap_channels_for_need_sk_2, sk_budget_2, max_epsilon_t_2, optimize_epsilons_T_2, ...
                                            max_P_X1_1, compute_marginal_PY_2, ...
                                            DEBUG, DEBUG_covert, DEBUG_covert_theorem_contraints, tolerance, seed,...
                                            N_epochs, draw_covert_point);

% for legends
% legends = { ['X2 cardinality: ', num2str(X2_cardinality)], ['X2 cardinality: ', num2str(X2_cardinality_2)]};
legends =   {   ['X2 cardinality: ', num2str(X2_cardinality), ' T cardinality: ', num2str(T_cardinality), ' Secret-Key budget: ',  num2str(sk_budget_1)], 
                ['X2 cardinality: ', num2str(X2_cardinality_2), ' T cardinality: ', num2str(T_cardinality_2), ' Secret-Key budget: ',  num2str(sk_budget_2)]};

if (swap_channels_for_need_sk || swap_channels_for_need_sk_2)
    need_sk = 'when we need a sk';
else
    need_sk = 'when we do not need a sk';
end

title_plot = ['Rate region ', need_sk, ' [seed=', num2str(seed), ']'];
x_title = 'Covert user square-root rate';
y_title = 'Non-Covert user rate';
z_title = 'Secret-Key square-root rate';

if length(r2_vect_1) > 1

    min_r1_1 = min(r1_vect_1);
    max_r1_1 = max(r1_vect_1);
    disp('--------------------------------------')
    disp(['The minimum covert rate r1 is: ', num2str(min_r1_1)]);
    disp(['The maximum covert rate r1 is: ', num2str(max_r1_1)]);
    disp('--------------------------------------')
    min_r2_1 = min(r2_vect_1);
    max_r2_1 = max(r2_vect_1);
    disp(['The minimum non covert rate r2 is: ', num2str(min_r2_1)]);
    disp(['The maximum non covert rate r2 is: ', num2str(max_r2_1)]); 
    disp('--------------------------------------')
    disp(['[Exp-1] The number of time P_{X1}(1) >', num2str(max_P_X1_1), ' is: ', num2str(cpt_higher_than_max_P_X1_1_1), '. Which corresponds to a percentage of: ', num2str(100*cpt_higher_than_max_P_X1_1_1/N_epochs), '%']);  
    disp(['[Exp-2] The number of time P_{X1}(1) >', num2str(max_P_X1_1), ' is: ', num2str(cpt_higher_than_max_P_X1_1_2), '. Which corresponds to a percentage of: ', num2str(100*cpt_higher_than_max_P_X1_1_2/N_epochs), '%']);  
    disp('--------------------------------------')
    disp(['[Exp-1] The number of time if did not verfiy the theorem conditions is: ', num2str(cpt_1), '. Which corresponds to a percentage of: ', num2str(100*cpt_1/N_epochs), '%']);  
    disp(['[Exp-2] The number of time if did not verfiy the theorem conditions is: ', num2str(cpt_2), '. Which corresponds to a percentage of: ', num2str(100*cpt_2/N_epochs), '%']);  
    disp('--------------------------------------')

    

    figure 
    ax2 = nexttile;
    scatter(ax2, r1_vect_1, r2_vect_1,'black');
    hold on
    scatter(ax2, r1_vect_2, r2_vect_2,'red');
    title(ax2, title_plot);
    xlabel(ax2,x_title);
    ylabel(ax2,y_title);
    legend(legends);
    
    

    ax4 = nexttile;
    scatter(ax4, rk_vect_1, r1_vect_1, 'black');
    hold on
    scatter(ax4, rk_vect_2, r1_vect_2, 'red'); 
    title(ax4, title_plot);
    xlabel(ax4,z_title);
    ylabel(ax4,x_title);
    legend(legends);

    ax5 = nexttile;
    scatter(ax5, rk_vect_1, r2_vect_1, 'black');
    hold on
    scatter(ax5, rk_vect_2, r2_vect_2, 'red');  
    title(ax5, title_plot);
    xlabel(ax5,z_title);
    ylabel(ax5,y_title);
    legend(legends);

    
    
    if (plot_3d)
        figure 
        ax3 = nexttile;
        scatter3(ax3, r1_vect_1, r2_vect_1, rk_vect_1, 'red');
        title(ax3,title_plot);
        xlabel(ax3,x_title);
        ylabel(ax3,y_title);
        zlabel(ax3,z_title);
        legend(legends(1));

        figure 
        ax6 = nexttile;
        scatter3(ax6, r1_vect_2, r2_vect_2, rk_vect_2, 'green');
        title(ax6,title_plot);
        xlabel(ax6,x_title);
        ylabel(ax6,y_title);
        zlabel(ax6,z_title);
        legend(legends(2));
    end

    if draw_convhull
        P = [r1_vect_1 r2_vect_1];
        [k,av] = convhull(P);
        
        figure 
        ax1 = nexttile;
        scatter(ax1, P(:,1),P(:,2),'')
        hold on
        plot(ax1, P(k,1),P(k,2))
        title(ax1,title_plot);
        xlabel(ax1,x_title);
        ylabel(ax1,y_title);
    end
end