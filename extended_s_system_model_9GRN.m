clc
clear all

global theta W L
%% Loading time series data

datatype = 'training';

if strcmp(datatype, 'training')
    load 9GRN_trainingdataset.mat
    p01 = p01mean; p02 = p02mean; p03 = p03mean;
    p04 = p04mean; p05 = p05mean; p06 = p06mean;
    p07 = p07mean; p08 = p08mean; p09 = p09mean;
elseif strcmp(datatype, 'validation')
    load 9GRN_validationdataset.mat
    p01 = p01val; p02 = p02val; p03 = p03val;
    p04 = p04val; p05 = p05val; p06 = p06val;
    p07 = p07val; p08 = p08val; p09 = p09val;
else
    error('Please select either training or validation')
end

pALL = [p01; p02; p03;...
    p04; p05; p06;...
    p07; p08; p09];

%% Defining external perturbation

tm = 48:1:72;
E1 = 0.001*ones(1,48);
E2 = (1/24)*tm - 2;
E3 = 1*ones(1,24);
E = [E1 E2 E3];
E = E(1:2:end-1); % Data sampled every 2 hours.

tche = 0:length(E)-1;
Light = sin((2*pi*tche/24)-pi/6) + 1.0001;

%% Estimated Parameters

pr01 = [1.0174,0.64823,-0.24169,0.82716,-1.8866,1.0577];
pr02 = [-1.0694,9.0585,1.9728];
pr03 = [-0.04984,0.49937,-0.1328];
pr04 = [3.2664,0.78375,0.069228,-2.3804,-0.78027,0.19118];
pr05 = [1.0481,1.0618,-1.011,-2.5859,2.3646];
pr06 = [-0.27401,2.5325,1.0798];
pr07 = [29.7754,-0.3429,-0.67803,-6.6253,0.0047053,0.0074363];
pr08 = [1.7845,0.17608,0.067548,0.70536,-1.6321,2.1995,0.70909];
pr09 = [1.6086,0.63347,-0.88802,1.9096,0.25114];

lineartheta = [pr01 pr02 pr03 pr04 pr05 pr06 pr07 pr08 pr09];

theta = lineartheta;

%% Initial condition and solving ODE
C = [p01(1) p02(1) p03(1) p04(1) p05(1) p06(1) p07(1) p08(1) p09(1)]*1;

Cinit = [p01(1) p02(1) p03(1) p04(1) p05(1) p06(1) p07(1) p08(1) p09(1)]*1;

ProteinLevel = [];
for t = 1:length(E)
    tspan = [t t+1];
    W = E(t);
    L = Light(t);
    [T,C] = ode45('extended_s_system_model_ODE',tspan,C(end,:));
    ProteinLevel = [ProteinLevel; C(end,:)];
end
ProteinLevel = [Cinit; ProteinLevel(1:end-1,:)];
%% Computing Error
M = length(p01);
% 01-09
ypredp01 = ProteinLevel(:,1)'; ysimp01 = ypredp01; yrealp01 = p01; ep01 = sum((ysimp01 - yrealp01).^2)/M;
ypredp02 = ProteinLevel(:,2)'; ysimp02 = ypredp02; yrealp02 = p02; ep02 = sum((ysimp02 - yrealp02).^2)/M;
ypredp03 = ProteinLevel(:,3)'; ysimp03 = ypredp03; yrealp03 = p03; ep03 = sum((ysimp03 - yrealp03).^2)/M;
ypredp04 = ProteinLevel(:,4)'; ysimp04 = ypredp04; yrealp04 = p04; ep04 = sum((ysimp04 - yrealp04).^2)/M;
ypredp05 = ProteinLevel(:,5)'; ysimp05 = ypredp05; yrealp05 = p05; ep05 = sum((ysimp05 - yrealp05).^2)/M;
ypredp06 = ProteinLevel(:,6)'; ysimp06 = ypredp06; yrealp06 = p06; ep06 = sum((ysimp06 - yrealp06).^2)/M;
ypredp07 = ProteinLevel(:,7)'; ysimp07 = ypredp07; yrealp07 = p07; ep07 = sum((ysimp07 - yrealp07).^2)/M;
ypredp08 = ProteinLevel(:,8)'; ysimp08 = ypredp08; yrealp08 = p08; ep08 = sum((ysimp08 - yrealp08).^2)/M;
ypredp09 = ProteinLevel(:,9)'; ysimp09 = ypredp09; yrealp09 = p09; ep09 = sum((ysimp09 - yrealp09).^2)/M;

eMax = [(max(p01)),(max(p02)),(max(p03)),(max(p04)),(max(p05)),...
    (max(p06)),(max(p07)),(max(p08)),(max(p09))];

wmse = (ep01./(max(p01)^2)) + (ep02./(max(p02)^2)) + (ep03./(max(p03)^2)) + (ep04./(max(p04)^2)) + (ep05./(max(p05)^2))...
    + (ep06./(max(p06)^2)) + (ep07./(max(p07)^2)) + (ep08./(max(p08)^2)) + (ep09./(max(p09)^2))

eIndv = [(ep01./(max(p01)^2)),(ep02./(max(p02)^2)),(ep03./(max(p03)^2)),(ep04./(max(p04)^2)),(ep05./(max(p05)^2)),...
    (ep06./(max(p06)^2)),(ep07./(max(p07)^2)),(ep08./(max(p08)^2)),(ep09./(max(p09)^2))];

average_e = wmse/9

%% Figure plotting
tp = 0:length(E)-1;

gN = {'ORA59', 'MYB51', 'LOL1', 'AT1G79150', 'ANAC055', 'a-ERF1', 'ATML1', 'CHE', 'RAP2.6L'};

for i = 1:9
    figure(1)
    subplot(3,3,i)
    plot(tp,pALL(i,:),'s-','Color',[0.5,0.5,0.5])
    hold on
    plot(tp,ProteinLevel(:,i),'LineWidth',2)
    
    title(gN{i})
    xticks([0 12 24 36 48])
    xticklabels({'0','24','48','72','96'})
end

