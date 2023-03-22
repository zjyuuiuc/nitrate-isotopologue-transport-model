%% function for calculating nitrate isotopologue concentrations and for conducting nitrate transport modeling
% Input: 
% [1] data: input data contain water age distribution of tile discharge (pQsoy),
% measured nitrate concentration in tile discharge (NO3), index
% corresponding to dates where tile discharge was sampled for nitrate
% measurements (sindex), diffusive mixing constant of nitrate (rd),
% measured d18O of nitrate (d18O), initial d18O of nitrate produced from
% nitrification (d18Oi), and isotope fractionation factor of
% denitrification (alpha).
% [2] theta: parameters to be calibrated. theta(1) is the equilibrium
% nitrate concentration in the soil (Ceq) and theta(2) is the
% denitrification rate constant (kd)

% output:
% MAPE: the average value of mean absolute percentage error of simulated
% nitrate concentration and d18O
% KGE: the average value of KGE of simulated nitrate concentration and d18O
% NO3sim: simulated nitrate concentration
% d18Osim: simulated d18O of nitrate

function [ MAPE, KGE, NO3sim, d18Osim] = nitratetransport(data, theta)

pQsoy = data.pQsoy; % water age distribution of tile discharge
NO3 = data.NO3; % measured nitrate concentration in tile discharge [mg N/L]
sindex = data.sindex; % index corresponding to dates where tile discharge was sampled for nitrate measurements
rd = data.rd; % diffusive mixing constant of nitrate [d-1]
d18O = data.d18O; % measured delta 18-O of nitrate [per mil]
d18Oi = data.d18Oi; % initial delta 18-O of nitrate produced from nitrification [per mil]
alpha = data.alpha; % isotope fractionation factor of denitrification

n = length(pQsoy);
VSMOW17 = 0.0003799; % 17O/16O of VSMOW
VSMOW18 = 0.0020052; % 18O/16O of VSMOW

Ceq = theta(1); % the equilibrium nitrate concentration in the soil 
O18i = (3*Ceq*(d18Oi/1000+1)*VSMOW18)/(1 + (0.52*d18Oi/1000+1)*VSMOW17 + (d18Oi/1000+1)*VSMOW18); % initial atomic concentration of the 18O-substituted nitrate isotopologue 
C18i = O18i; % initial molecular concentration of the 18O-substituted nitrate isotopologue
O17i = (3*Ceq*(0.52*d18Oi/1000+1)*VSMOW17)/(1 + (0.52*d18Oi/1000+1)*VSMOW17 + (d18Oi/1000+1)*VSMOW18); % initial atomic concentration of the 17O-substituted nitrate isotopologue 
C17i = O17i; % initial molecular concentration of the 17O-substituted nitrate isotopologue
C16i = Ceq - C17i - C18i; % initial molecular concentration of the major nitrate isotopologue

kdnf = theta(2); % first-order rate constant of denitrification [d-1]

C16t = zeros(n,1); % set up a vector for simulated molecular concentration of the major nitrate isotopologue
C18t = zeros(n,1); % set up a vector for simulated molecular concentration of the 18O-substituted nitrate isotopologue
for i = 1:n
    
    cage = 1:length(pQsoy{i,1}); % age span at each time step
    
    C16t(i) = sum(pQsoy{i,1}.*(C16i.*(1-exp(-rd*cage)).*exp(-kdnf*cage))'); % simulated molecular concentration of the major nitrate isotopologue at time t
    C18t(i) = sum(pQsoy{i,1}.*(C18i.*(1-exp(-rd*cage)).*exp(-(kdnf/alpha)*cage))'); % simulated molecular concentration of the 18O-substituted nitrate isotopologue at time t
    
end

C17t = (0.52.*VSMOW17.*C18t + 1.44.*VSMOW17.*VSMOW18.*C16t + 0.96.*VSMOW17.*VSMOW18.*C18t)./(VSMOW18 - 0.96*VSMOW18*VSMOW17); % calculate simulated molecular concentration of the 17O-substituted nitrate isotopologue
NO3sim = (C16t + C17t + C18t).*62; % simulated nitrate concentration
O17sim = C17t; % simulated atomic concentration of 17O
O18sim = C18t; % simulated atomic concentration of 18O
O16sim = NO3sim./62.*3 - O17sim - O18sim; % simulated atomic concentration of 16O
d18Osim = ((O18sim./O16sim)./VSMOW18 - 1).*1000;

NO3sim_sindex = NO3sim(sindex);
d18Osim_sindex = d18Osim(sindex);

F(:,1) = abs((NO3 - NO3sim_sindex)./NO3);
F(:,2) = abs((d18O - d18Osim_sindex)./d18O);
MAPE = mean(mean(F));

kgeC = of_KGE(NO3,NO3sim_sindex);
kgeO = of_KGE(d18O,d18Osim_sindex);
KGE = (kgeC+kgeO)/2;



