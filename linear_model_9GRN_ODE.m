function dC = linear_model_9GRN_ODE(t,C)

global theta W L

dC = zeros(9,1);

% ORA59
dC(1) = theta(1)*C(2) + theta(2)*C(5) + theta(3)*C(6) + theta(4)*C(1) + theta(5);

% MYB51
dC(2) = theta(6)*C(2) + theta(7) + theta(8)*W;

% LOL1
dC(3) = theta(9)*C(3) + theta(10) + theta(11)*W;

% AT1G79150
dC(4) = theta(12)*C(8) + theta(13)*C(9) + theta(14)*C(4) + theta(15) + theta(16)*W;

% ANAC055
dC(5) = theta(17)*C(9) + theta(18)*C(5) + theta(19) + theta(20)*W;

% a-ERF-1
dC(6) = theta(21)*C(6) + theta(22) + theta(23)*W;

% ATML1
dC(7) = theta(24)*C(9) + theta(25)*C(7) + theta(26) + W*theta(27) + theta(37)*L;

% CHE
dC(8) = theta(28)*C(3) + theta(29)*C(4) + theta(30)*C(7) + theta(31)*C(8) + theta(32) + theta(38)*L;

% RAP2.6L
dC(9) = theta(33)*C(5) + theta(34)*C(9) + theta(35) + W*theta(36);

