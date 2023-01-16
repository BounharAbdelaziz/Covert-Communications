%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sannity checks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q           = [0.5 0.5];
P_T         = [0.1 0.9];

% sanity check of the implemented functions so far
h = InformationTheory.entropy(P_T);
T = [1 1];
e = InformationTheory.expectation(P_T, T); % should be 1

d1= InformationTheory.relative_entropy(P_T, Q); % should be 0.531

Q           = [0.5 0.5];
P_T         = [0.1 0.9];
k=2;
chi = InformationTheory.chi_k_distance(P_T, Q, k); % should be 1.77

disp(['The entropy of T~P_T is: ', num2str(h)]);
disp(['The expected value of T: ', num2str(e)]);
disp(['The relative entropy d1 is: ', num2str(d1)]);
disp(['The chi square is: ', num2str(chi)]);
disp('------------------')
