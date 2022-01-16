function dC = s_system_model_ODE(t,C)


global theta W L

dC = zeros(9,1);

% ORA59
dC(1) = theta(1)*C(2)^(theta(2))*C(5)^(theta(3))*C(6)^(theta(4)) + theta(5)*C(1);

% MYB51
dC(2) = theta(6)*W^(theta(7)) + theta(8)*C(2);

% LOL1
dC(3) = theta(9)*W^(theta(10)) + theta(11)*C(3);

% AT1G79150
dC(4) = theta(12)*C(8)^(theta(13))*C(9)^(theta(14))*W^(theta(15)) + theta(16)*C(4);

% ANAC055
dC(5) = theta(17)*C(9)^(theta(18))*W^(theta(19)) + theta(20)*C(5);

% a-ERF-1
dC(6) = theta(21)*W^(theta(22)) + theta(23)*C(6);

% ATML1
dC(7) = theta(24)*C(9)^(theta(25))*W^(theta(26))*L^(theta(27)) + theta(28)*C(7);

% CHE
dC(8) = theta(29)*C(3)^(theta(30))*C(4)^(theta(31))*C(7)^(theta(32))*L^(theta(33)) + theta(34)*C(8);

% RAP2.6L
dC(9) = theta(35)*C(5)^(theta(36))*W^(theta(37)) + theta(38)*C(9);

