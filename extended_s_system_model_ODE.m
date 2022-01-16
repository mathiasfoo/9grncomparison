function dC = extended_s_system_model_ODE(t,C)


global theta W L

dC = zeros(9,1);

% ORA59
dC(1) = real(theta(1)*C(2)^(theta(2))*C(5)^(theta(3))*C(6)^(theta(4)) + theta(5)*C(1) + theta(6));

% MYB51
dC(2) = real(theta(7)*C(2) + theta(8) + theta(9)*W);

% LOL1
dC(3) = real(theta(10)*C(3) + theta(11) + theta(12)*W);

% AT1G79150
dC(4) = real(theta(13)*C(8)^(theta(14))*C(9)^(theta(15)) + theta(16)*C(4) + theta(17) + theta(18)*W);

% ANAC055
dC(5) = real(theta(19)*C(9)^(theta(20)) + theta(21)*C(5) + theta(22) + theta(23)*W);

% a-ERF-1
dC(6) = real(theta(24)*C(6) + theta(25) + theta(26)*W);

% ATML1
dC(7) = real(theta(27)*C(9)^(theta(28)) + theta(29)*C(7) + theta(30) + theta(31)*W + theta(32)*L);

% CHE
dC(8) = real(theta(33)*C(3)^(theta(34))*C(4)^(theta(35))*C(7)^(theta(36)) + theta(37)*C(8) + theta(38) + theta(39)*L);

% RAP2.6L
dC(9) = real(theta(40)*C(5)^(theta(41)) + theta(42)*C(9) + theta(43) + theta(44)*W);

