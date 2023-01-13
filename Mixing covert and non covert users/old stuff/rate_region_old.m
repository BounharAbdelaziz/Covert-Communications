
% \matchal{T} is of cardinality 4, choose it uniform for simplicity
P_T         = [0.25 0.01 0.49 0.25]; 
% Binary alphabet for the covert user
P_X1        = [0.95 0.05]; % [P_X1(0) P_X1(1)]
% Binary alphabet for the non covert user
P_X2        = [0.95 0.05]; % [P_X2(0) P_X2(1)]
% Binary output alphabet
%P_Y         = [1 0]; % [P_Y(0) P_Y(1)]

% BSC channel with cross probability epsilon
epsilon = 0.05;
W_Y_X2    = zeros(2,2);
W_Y_X2(1,1) = 1-epsilon;
W_Y_X2(1,2) = epsilon;
W_Y_X2(2,1) = epsilon;
W_Y_X2(2,2) = 1-epsilon;

Q           = [0.5 0.25 0.125 0.125];

% sanity check of the implemented functions so far
h = InformationTheory.entropy(P_T);
T = [1 1 1 1];
e = InformationTheory.expectation(P_T, T);

d1= InformationTheory.relative_entropy(P_T, Q);

% MAC Channel matrix
W_Y_X1_X2 = InformationTheory.find_mac_channel_law_from_p2p_channel_law(P_X1, W_Y_X2);
W_0_0_0 = double(W_Y_X1_X2(1));
W_0_1_0 = double(W_Y_X1_X2(2)); 
W_0_0_1 = double(W_Y_X1_X2(3)); 
W_0_1_1 = double(W_Y_X1_X2(4));
W_1_0_0 = double(W_Y_X1_X2(5));
W_1_1_0 = double(W_Y_X1_X2(6));
W_1_0_1 = double(W_Y_X1_X2(7));
W_1_1_1 = double(W_Y_X1_X2(8));

W_Y_X1 = InformationTheory.find_marginal_channel_law_from_mac_channel_law(P_X2, W_Y_X1_X2);

disp(['The entropy of T~P_T is: ', num2str(h)]);
disp(['The expected value of T: ', num2str(e)]);
disp(['The relative entropy d1 is: ', num2str(d1)]);
disp('------------------')
disp('W_Y_X1_X2 is: ');
disp(W_Y_X1_X2);
disp('------------------')
disp('W_Y_X1 is: ');
disp(W_Y_X1);
disp('------------------')
P_Y = InformationTheory.find_marginal_output_law_from_p2p_channel_law(P_X1, W_Y_X1);% [P_Y(0) P_Y(1)]
disp('P_Y is: ');
disp(P_Y); 
disp('------------------')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rate Region simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W_Y_X1_0 = [W_Y_X1(1,1), W_Y_X1(2,1)];
W_Y_X1_0_X2 = [W_0_0_0, W_0_0_1, W_1_0_0, W_1_0_1];
W_Y_X1_1_X2 = [W_0_1_0, W_0_1_1, W_1_1_0, W_1_1_1];
disp('W_Y_X1_0 is: ');
disp(W_Y_X1_0);
disp('------------------')
disp('W_Y_X1_0_X2 is: ');
disp(W_Y_X1_0_X2);
disp('------------------')

diff_log_P_Q = log2(W_Y_X1_1_X2)-log2(W_Y_X1_0_X2);

d1= InformationTheory.relative_entropy(W_Y_X1_1_X2, W_Y_X1_0_X2);
disp(d1)
disp(['The relative entropy d1 is: ', num2str(d1)]);

%r2 = CovertCommunication.non_covert_rate(P_T, P_X2, P_Y, W_Y_X1_0, W_Y_X1_0_X2);
%disp(['The non covert rate r2 is: ', num2str(r2)]);

% Fix the blocklength n
n = 1e6;

% define a simple channel

% covert constraint


% verify the constaints (absolute continuity)


% compute the rate


% do a loop on the covert constraint