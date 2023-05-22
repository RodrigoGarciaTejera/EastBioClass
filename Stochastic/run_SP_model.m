%%-------------------------------------------------------------------------
% Rodrigo García-Tejera
% Contact info: rodrigo.garcia@ed.ac.uk ; rgarcia@fisica.edu.uy
% Affiliations: 
% Centre for Regenerative Medicine, University of Edinburgh, UK. 
% Facultad de Ciencias, UdelaR, Uruguay.
% Created: May-2023
%%-------------------------------------------------------------------------


%Stochastic simulations of the SP model. The S species represent naïve 
%(unlicensed) states, and P species represent the licensed states.

%SP model: S --> 2S ; S <--> P; 2P --> 0

%All the rates are considered with respect to the proliferation rate. This
%is equivalent to defining the non-dimensional time $\tau=k_1 t$, with
%$k_1$ the proliferation rate. 

% clears workspace
clear

% vector of times when number of cells are going to be assessed
timeVector=0:0.05:40;
  
%% Model SP parameters 
nS=18; %homeostatic number of S cells 
nP=26; %homeostatic number of P cells

k1=1; %birth rate 
k2=5; % S-->P switching rate
k3=nS/nP*(k2-1); %sets k3 so nS and nP are steady states
k4=nS/2/nP^2;  %sets k3 so nS and nP are steady states

Omega=1; %system volume

rateVector=[k1,k2,k3,k4];%Reaction rates
stoichiometricMatrix=[1,-1,1,0;0,1,-1,-2]; %Stoichiometric matrix
propensityFunctions= @(x) [x(1),x(1),x(2),x(2)*(x(2)-1)/Omega]; %Propensity 
%functions. Set to mass-action propensities, change if needed.
 
%creation of RxNet object
SPmodel=RxnNet(stoichiometricMatrix,rateVector,propensityFunctions);

%% SSA (stochastic simulation algorithm) calculations
 
initialNumbers=[nS,nP]; %initial numbers of S and P cells
% 
trajectories=SPmodel.SSA(timeVector,initialNumbers,true); %performs one 
%realization of SSA.
 
trajectoriesS=trajectories(:,1); %trajectories of S species
trajectoriesP=trajectories(:,2); %trajectories of P species


%% Plotting 

%This can be commented after one run
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex'); 

close all

%This figure plots the stochastic trajectories-----------------------------
figure,
subplot(2,1,1)
ax=gca; ax.FontSize=22; ax.LineWidth=1.5;
xlabel('time','FontSize',24,'Interpreter','latex');
ylabel('$n_S$','FontSize',24,'Interpreter','latex');
yline(nS,'-.k','LineWidth',2)
xlim([0,timeVector(end)])
ylim([0,50])


hold on
subplot(2,1,2)
ax=gca; ax.FontSize=22; ax.LineWidth=1.5;
xlabel('time','FontSize',24,'Interpreter','latex');
ylabel('$n_P$','FontSize',24,'Interpreter','latex');
yline(nP,'-.k','LineWidth',2)
xlim([0,timeVector(end)])
ylim([0,50])


hold on

for r=1:numel(timeVector)   
    subplot(2,1,1)
    plot(timeVector(1:r),trajectoriesS(1:r),'Color',[0 0.4470 0.7410],'LineWidth',2)
    subplot(2,1,2)
    plot(timeVector(1:r),trajectoriesP(1:r),'Color',[0.8500 0.3250 0.0980],'LineWidth',2)   
    pause(0.0002)
end

%--------------------------------------------------------------------------

%This figure plots the histograms------------------------------------------
figure
subplot(2,1,1)
histogram(trajectoriesS,(min(trajectoriesS)-0.5):(max(trajectoriesS)+0.5), ...
'Normalization','pdf'),
ax=gca; ax.FontSize=22; ax.LineWidth=1.5;
hold on,
xlabel('number of cells','FontSize',24,'Interpreter','latex');
ylabel('probability','FontSize',24,'Interpreter','latex');
title(['S species (naive stem cells)'],'FontSize',22,'Interpreter','latex');

subplot(2,1,2),histogram(trajectoriesP,(min(trajectoriesP)-0.5): ... 
(max(trajectoriesP)+0.5),'Normalization','pdf'),
ax=gca; ax.FontSize=22; ax.LineWidth=1.5;
hold on,
xlabel('number of cells','FontSize',24,'Interpreter','latex');
ylabel('probability','FontSize',24,'Interpreter','latex');
title(['P species (licensed stem cells)'],'FontSize',22,'Interpreter','latex');
