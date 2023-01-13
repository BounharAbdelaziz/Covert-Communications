% Legitimate user channel
%     W_Y_X1_0_X2_0 = [0.1 0.7 0.1 0.1]; % [y=0, y=1, y=2, y=3]
%     W_Y_X1_0_X2_1 = [0.3 0.2 0.4 0.1]; % [y=0, y=1, y=2, y=3]
%     W_Y_X1_1_X2_0 = [0.15 0.05 0.2 0.6]; % [y=0, y=1, y=2, y=3]
%     W_Y_X1_1_X2_1 = [0.45 0.05 0.3 0.2]; % [y=0, y=1, y=2, y=3]

%     W_Y_X1_0_X2_0 = [0.2 0.5 0.2 0.1]; % [y=0, y=1, y=2, y=3]
%     W_Y_X1_0_X2_1 = [0.3 0.3 0.2 0.2]; % [y=0, y=1, y=2, y=3]
%     W_Y_X1_1_X2_0 = [0.15 0.05 0.3 0.5]; % [y=0, y=1, y=2, y=3]
%     W_Y_X1_1_X2_1 = [0.45 0.15 0.13 0.27]; % [y=0, y=1, y=2, y=3]

% if W_Y_X1_0_X2 = W_Y_X1_1_X2 then the covert rate should be zero
%     W_Y_X1_0_X2_0 = [0.15 0.05 0.2 0.6]; % [y=0, y=1, y=2, y=3]
%     W_Y_X1_0_X2_1 = [0.45 0.05 0.3 0.2]; % [y=0, y=1, y=2, y=3]
%     W_Y_X1_1_X2_0 = [0.15 0.05 0.2 0.6]; % [y=0, y=1, y=2, y=3]
%     W_Y_X1_1_X2_1 = [0.45 0.05 0.3 0.2]; % [y=0, y=1, y=2, y=3]

% Pentagonal channel
W_Y_X1_0_X2_0 = [1-pw pw 0 0]; % [y=0, y=1, y=2, y=3]
W_Y_X1_0_X2_1 = [0 1-pw pw 0]; % [y=0, y=1, y=2, y=3]
W_Y_X1_1_X2_0 = [0 0 1-pw pw]; % [y=0, y=1, y=2, y=3]
W_Y_X1_1_X2_1 = [pw 0 0 1-pw]; % [y=0, y=1, y=2, y=3]
%     
%     eps=pw/80;
%     W_Y_X1_0_X2_0 = [1-pw-eps/2 pw-eps/2 eps/2 eps/2]; % [y=0, y=1, y=2, y=3]
%     W_Y_X1_0_X2_1 = [eps/4 1-pw-eps/2 pw-eps/2 3*eps/4]; % [y=0, y=1, y=2, y=3]
%     W_Y_X1_1_X2_0 = [eps/2 eps/2 1-pw-eps/2 pw-eps/2]; % [y=0, y=1, y=2, y=3]
%     W_Y_X1_1_X2_1 = [pw-eps/2 3*eps/4 eps/4 1-pw-eps/2]; % [y=0, y=1, y=2, y=3]

% intervert role of X1 and X2 to see if the rates will be swiped
% W_Y_X1_0_X2_0 = [1-pw pw 0 0]; % [y=0, y=1, y=2, y=3]
% W_Y_X1_0_X2_1 = [0 0 1-pw pw]; % [y=0, y=1, y=2, y=3]
% W_Y_X1_1_X2_0 = [0 1-pw pw 0]; % [y=0, y=1, y=2, y=3]
% W_Y_X1_1_X2_1 = [pw 0 0 1-pw]; % [y=0, y=1, y=2, y=3]

W_Y_X1_X2 = zeros(4,4);
W_Y_X1_X2(1,:) = W_Y_X1_0_X2_0;
W_Y_X1_X2(2,:) = W_Y_X1_0_X2_1;
W_Y_X1_X2(3,:) = W_Y_X1_1_X2_0;
W_Y_X1_X2(4,:) = W_Y_X1_1_X2_1;

% Eve's channel
%     W_Z_X1_0_X2_0 = [0.2 0.6 0.1 0.1]; % [y=0, y=1, y=2, y=3]
%     W_Z_X1_0_X2_1 = [0.1 0.3 0.5 0.1]; % [y=0, y=1, y=2, y=3]
%     W_Z_X1_1_X2_0 = [0.05 0.05 0.35 0.6]; % [y=0, y=1, y=2, y=3]
%     W_Z_X1_1_X2_1 = [0.65 0.05 0.05 0.25]; % [y=0, y=1, y=2, y=3]
%     
% intervert role of X1 and X2 to see if the rates will be swiped
% W_Z_X1_0_X2_0 = [1-pw_eve pw_eve 0 0]; % [z=0, z=1, z=2, z=3]
% W_Z_X1_0_X2_1 = [0 0 1-pw_eve pw_eve]; % [z=0, z=1, z=2, z=3]
% W_Z_X1_1_X2_0 = [0 1-pw_eve pw_eve 0]; % [z=0, z=1, z=2, z=3]
% W_Z_X1_1_X2_1 = [pw_eve 0 0 1-pw_eve]; % [z=0, z=1, z=2, z=3]

% pentagonal channel
%     W_Z_X1_0_X2_0 = [1-pw_eve pw_eve 0 0]; % [z=0, z=1, z=2, z=3]
%     W_Z_X1_0_X2_1 = [0 1-pw_eve pw_eve 0]; % [z=0, z=1, z=2, z=3]
%     W_Z_X1_1_X2_0 = [0 0 1-pw_eve pw_eve]; % [z=0, z=1, z=2, z=3]
%     W_Z_X1_1_X2_1 = [pw_eve 0 0 1-pw_eve]; % [z=0, z=1, z=2, z=3]

eps=pw_eve/20;
W_Z_X1_0_X2_0 = [1-pw_eve-eps/2 pw_eve-eps/2 eps/2 eps/2]; % [y=0, y=1, y=2, y=3]
W_Z_X1_0_X2_1 = [eps/4 1-pw_eve-eps/2 pw_eve-eps/2 3*eps/4]; % [y=0, y=1, y=2, y=3]
W_Z_X1_1_X2_0 = [eps/2 eps/2 1-pw_eve-eps/2 pw_eve-eps/2]; % [y=0, y=1, y=2, y=3]
W_Z_X1_1_X2_1 = [pw_eve-eps/2 3*eps/4 eps/4 1-pw_eve-eps/2]; % [y=0, y=1, y=2, y=3]

W_Z_X1_X2 = zeros(4,4);
W_Z_X1_X2(1,:) = W_Z_X1_0_X2_0;
W_Z_X1_X2(2,:) = W_Z_X1_0_X2_1;
W_Z_X1_X2(3,:) = W_Z_X1_1_X2_0;
W_Z_X1_X2(4,:) = W_Z_X1_1_X2_1;