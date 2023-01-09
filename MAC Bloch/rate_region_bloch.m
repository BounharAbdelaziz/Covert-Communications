% seed for reproducibility
% seed                            = 1998;
% rng(seed);

% for simulation
n_users                         = 2;
X1_cardinality                  = 2;
X2_cardinality                  = 2;
X1_X2_cardinality               = 4;
Y_cardinality                   = 4;
N_epochs                        = 5000;
optimize_rhos                   = 1;
generate_random_channel_laws    = 0;
% maximum value of P_{X1}(1)
max_P_X1_1                      = 1;
cpt_higher_than_max_P_X1_1      = 0;
% noise probabilities
pw                              = 0.1;
pw_eve                          = 0.15;

% error in the if conditions (1==1.00)
tolerance                       = 1e-6;

% Prints if debug mode
DEBUG                           = 0;
DEBUG_covert                    = 1;
DEBUG_covert_theorem_contraints = 1;

% for fixed epsilon_T
if (ismembertol(optimize_rhos, 0, tolerance))
    disp('[INFO] Running experiment with fixed rhos')
    rhos           = InformationTheory.generate_probability_vector(n_users,1,1,0,1);
else
    disp('[INFO] Running experiment with different rhos')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fix the channel law %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MAC Channel matrix : We consider the channel Y = (2*X1 + X2 + w) % 4 where all variables are binary => \mathcal{Y} = \{0,1,2,3}
% for each Y and fixed (x1, x2) we have a vector of probabilities [y=0, y=1, y=2, y=3]

if (ismembertol(generate_random_channel_laws, 1, tolerance))
    
    W_Y_X1_X2 = transpose(InformationTheory.generate_probability_vector(X1_X2_cardinality,Y_cardinality,1,0,1)); % InformationTheory.generate_random_channel_matrix(m, n, smallest_poba_value_bob);
    W_Z_X1_X2 = transpose(InformationTheory.generate_probability_vector(X1_X2_cardinality,Y_cardinality,1,0,1)); % InformationTheory.generate_random_channel_matrix(m, n, smallest_poba_value_eve);
    verified_conditions = CovertCommunicationBloch.check_theorem_conditions_bloch(W_Y_X1_X2, W_Z_X1_X2, DEBUG_covert_theorem_contraints);
    
    % loop while untill constraint on absolute continuity and difference are met  
    while (verified_conditions < 1)
        W_Y_X1_X2 = transpose(InformationTheory.generate_probability_vector(X1_X2_cardinality,Y_cardinality,1,0,1)); %generate_random_channel_matrix
        W_Z_X1_X2 = transpose(InformationTheory.generate_probability_vector(X1_X2_cardinality,Y_cardinality,1,0,1)); %generate_random_channel_matrix
        verified_conditions = CovertCommunicationBloch.check_theorem_conditions_bloch(W_Y_X1_X2, W_Z_X1_X2, DEBUG_covert_theorem_contraints);
    end
else
    % Binary noise rv follows a bernoulli distribution    
    P_W                 = [1-pw pw]; % [P_W(0) P_W(1)], low probability of noise    
    P_W_eve             = [1-pw_eve pw_eve]; % [P_W_eve(0) P_W_eve(1)], higher probability of noise for the adversary channel

    % Legitimate user channel
    eps=pw;
    W_Y_X1_0_X2_0 = [1-pw-eps/2 pw-eps/2 eps/2 eps/2]; % [y=0, y=1, y=2, y=3]
    W_Y_X1_0_X2_1 = [eps/4 1-pw-eps/2 pw-eps/2 3*eps/4]; % [y=0, y=1, y=2, y=3]
    W_Y_X1_1_X2_0 = [eps/2 eps/2 1-pw-eps/2 pw-eps/2]; % [y=0, y=1, y=2, y=3]
    W_Y_X1_1_X2_1 = [pw-eps/2 3*eps/4 eps/4 1-pw-eps/2]; % [y=0, y=1, y=2, y=3]
    
    W_Y_X1_X2 = zeros(4,4);
    W_Y_X1_X2(1,:) = W_Y_X1_0_X2_0;
    W_Y_X1_X2(2,:) = W_Y_X1_0_X2_1;
    W_Y_X1_X2(3,:) = W_Y_X1_1_X2_0;
    W_Y_X1_X2(4,:) = W_Y_X1_1_X2_1;
    
    % Eve's channel    
    eps=pw_eve;
    W_Z_X1_0_X2_0 = [1-pw_eve-eps/2     pw_eve-eps/2        eps/4               3*eps/4]; % [y=0, y=1, y=2, y=3]
    W_Z_X1_0_X2_1 = [eps/4              1-pw_eve-eps/2      pw_eve-eps/2        3*eps/4]; % [y=0, y=1, y=2, y=3]
    W_Z_X1_1_X2_0 = [3*eps/4            eps/4               1-pw_eve-eps/2      pw_eve-eps/2]; % [y=0, y=1, y=2, y=3]
    W_Z_X1_1_X2_1 = [pw_eve-eps/2       3*eps/4             eps/4               1-pw_eve-eps/2]; % [y=0, y=1, y=2, y=3]
    
    W_Z_X1_X2 = zeros(4,4);
    W_Z_X1_X2(1,:) = W_Z_X1_0_X2_0;
    W_Z_X1_X2(2,:) = W_Z_X1_0_X2_1;
    W_Z_X1_X2(3,:) = W_Z_X1_1_X2_0;
    W_Z_X1_X2(4,:) = W_Z_X1_1_X2_1;

    verified_conditions = CovertCommunicationBloch.check_theorem_conditions_bloch(W_Y_X1_X2, W_Z_X1_X2, DEBUG_covert_theorem_contraints);
    if (verified_conditions < 1)
        disp('[ERROR-Constraint] Absolute continuity is not met! Check bellow which one you need to fix.');
        disp(['[INFO] absolute_continuity_bob : ', num2str(absolute_continuity_bob)]);
        disp(['[INFO] absolute_continuity_eve : ', num2str(absolute_continuity_eve)]);
        disp(['[INFO] different_output_distributions_eve : ', num2str(different_output_distributions_eve)]);
        disp(['[INFO] without_secret_key_condition : ', num2str(without_secret_key_condition)]);
        return
    end

    % Verify if all distributions are valid
    sum_pw = sum(P_W);
    assert(ismembertol(sum_pw, 1, tolerance), 'P_W should sum to 1.');
    % ismembertol(a, b, tolerance) compute the difference |a-b| and check if it's \leq to the tolerance. It's used here to avoid problem 1.000 == 1
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rate Region simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Covert user rates
users_rates = zeros(N_epochs,n_users);

% saving the rho's as we want the find the \rho^*
rhos_saved = zeros(N_epochs,n_users);

% iterate over all input distributions.
for epoch = progress(1:N_epochs)
    
    if (ismembertol(optimize_rhos, 1, tolerance))
    % The epsilon's that normalizes the probability of sending one for the covert user (we choose them summing to one for now but it's not necessary)
        rhos           = InformationTheory.generate_probability_vector(n_users,1,1,0,1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rate Region simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if DEBUG
        disp('W_Y_X1_X2 is: ');
        disp(W_Y_X1_X2);
    end

    users_rates_epoch = CovertCommunicationBloch.covert_message_rate(rhos, W_Y_X1_X2, W_Z_X1_X2, DEBUG);
    users_rates(epoch,:) = transpose(users_rates_epoch);

    if DEBUG
        disp('For P_T: ');
        disp(P_T);
        disp(['The non covert rate r2 is: ', num2str(r2)]);
        disp(['The covert rate r1 is: ', num2str(r1)]);
    end

    rhos_saved(epoch,1)= rhos(1);
    rhos_saved(epoch,2)= rhos(2);
end

if length(users_rates(:,1)) > 1
    r1_vect = users_rates(:,1);
    r2_vect = users_rates(:,2);
    
    % we will take the minimum difference, as it reprents the \rho^*
    difference_rates = abs(r2_vect-r1_vect);
    [tmp, indice_star] = min(difference_rates);
    rho_star_1 = rhos_saved(indice_star, 1);
    rho_star_2 = rhos_saved(indice_star, 2);
    r1_star = r1_vect(indice_star);
    r2_star = r2_vect(indice_star);

    % region of covertness point
    rhos_star = zeros(1,2);
    rhos_star(1) = rho_star_1;
    rhos_star(2) = rho_star_2;
    covert_rhos_starpoints(2) = rho_star_2;

    covert_points = CovertCommunicationBloch.compute_covert_point(W_Z_X1_X2, rhos_star, n_users, Y_cardinality);
    covert_points_1 = covert_points(1);
    covert_points_2 = covert_points(2);

    min_r1 = min(r1_vect);
    max_r1 = max(r1_vect);

    min_r2 = min(r2_vect);
    max_r2 = max(r2_vect);

    disp('--------------------------------------')
    disp(['The minimum covert rate r1 is: ', num2str(min_r1)]);
    disp(['The maximum covert rate r1 is: ', num2str(max_r1)]);
    disp('--------------------------------------')
    disp(['The minimum covert rate r2 is: ', num2str(min_r2)]);
    disp(['The maximum covert rate r2 is: ', num2str(max_r2)]);
    disp('--------------------------------------')
    disp(['The best \rho for user 1 is: ', num2str(rho_star_1), 'the corresponding rate for user 1 is: ',num2str(r1_star)]);
    disp(['The best \rho for user 2 is: ', num2str(rho_star_2), 'the corresponding rate for user 2 is: ',num2str(r2_star)]);
    disp('--------------------------------------')
    disp(['The number of time P_{X1}(1) >', num2str(max_P_X1_1), ' is: ', num2str(cpt_higher_than_max_P_X1_1), '. Which corresponds to a percentage of: ', num2str(100*cpt_higher_than_max_P_X1_1/N_epochs), '%']);  
    disp('--------------------------------------')


%     P = [r1_vect r2_vect];
%     [k,av] = convhull(P);

%     figure 
%     ax1 = nexttile;
%     scatter(ax1, P(:,1),P(:,2),'')
%     hold on
%     plot(ax1, P(k,1),P(k,2))
%     title(ax1,'MAC Rate region - Bloch');
%     xlabel(ax1,'Covert user 1');
%     ylabel(ax1,'Covert user 2');

    figure 
    ax2 = nexttile;
    scatter(ax2, r1_vect, r2_vect,'black');
    hold on;

    % Ploting the reliability constraint
    plot(r1_star, r2_star, 'o','MarkerSize',12);

    x1 = linspace(0,r1_star,100);
    y1=r2_star;
    
    x2 = r1_star;
    y2 = linspace(0,r2_star,100);

    plot(x1,y1*ones(size(x1)), '--g')
    plot(x2*ones(size(y2)),y2, '--g')

    % Ploting the covert constraint
    plot(covert_points_1, covert_points_2, '*','MarkerSize',12);
%     yline(covert_points_1);
%     xline(covert_points_2)

    x1 = linspace(covert_points_1,max_r1,100);
    y1=covert_points_2;
    
    x2 = covert_points_1;
    y2 = linspace(max_r2, covert_points_2,100);

    plot(x1,y1*ones(size(x1)), '--r')
    plot(x2*ones(size(y2)),y2, '--r')

    legend( 'Boundary of covert capicity region', ...
            'maximum achievable covert throughput', ...
            'Reliability constraint on (\rho_1^*, \rho_2^*)', ...
            'Reliability constraint on (\rho_1^*, \rho_2^*)', ...
            'Covert constraint on (\rho_1^*, \rho_2^*)')
    
    title(ax2,'MAC Rate region - Bloch');
    xlabel(ax2,'Covert user 1');
    ylabel(ax2,'Covert user 2');

end