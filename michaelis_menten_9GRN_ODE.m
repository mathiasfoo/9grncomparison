function dC = michaelis_menten_9GRN_ODE(t,C)


global theta W L

dC = zeros(9,1);

% ORA59
dC(1) = theta(1) - theta(2)*C(1) + (theta(3)*C(2)^2)/(theta(4)^2+C(2)^2) + theta(5)/(theta(6)^2+C(5)^2) + theta(7)*C(6)^2/(theta(8)^2+C(6)^2);

% MYB51
dC(2) = theta(9) - theta(10)*C(2) + (theta(11))/(theta(12)^2+W^2);

% LOL1
dC(3) = theta(13) - theta(14)*C(3) + (theta(15)*W^2)/(theta(16)^2+W^2);

% AT1G79150
dC(4) = theta(17) - theta(18)*C(4) + theta(19)/(theta(20)^2+W^2) + (theta(21)*C(8)^2)/(theta(22)^2+C(8)^2) + theta(23)/(theta(24)^2+C(9)^2);

% ANAC055
dC(5) = theta(25) - theta(26)*C(5) + (theta(27)*W^2)/(theta(28)^2+W^2) + (theta(29))/(theta(30)^2+C(9)^2);

% a-ERF-1
dC(6) = theta(31) - theta(32)*C(6) + (theta(33)*W^2)/(theta(34)^2+W^2);

% ATML1
dC(7) = theta(35) - theta(36)*C(7) + (theta(37)*L^2)/(theta(38)^2+L^2) + theta(39)/(theta(40)^2+W^2) + theta(41)/(theta(42)^2+C(9)^2);

% CHE
dC(8) = theta(43) - theta(44)*C(8) + (theta(45)*L^2)/(theta(46)^2+L^2) + (theta(47)*C(3)^2)/(theta(48)^2+C(3)^2) +(theta(49)*C(4)^2)/(theta(50)^2+C(4)^2) +(theta(51)*C(7)^2)/(theta(52)^2+C(7)^2);

% RAP2.6L
dC(9) = theta(53) - theta(54)*C(9) + (theta(55)*W^2)/(theta(56)^2+W^2) + (theta(57)*C(5)^2)/(theta(58)^2+C(5)^2);
