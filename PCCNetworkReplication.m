%%% Replication of Credit Network model in dynamic partner selection

clear 

%% Simulation size parameters
T = 1000; %Number of time steps
D = 500; %Number of D firms
U = 250; %Number of U firms
B = 100; %Number of B banks


%% Simulation parameters
phi = 2; %Finance constraints on production
beta = 0.9; %Finance constraints on production
alpha = 0.01; %Interest rate setting
gamma = 0.5; % Intermediate goods requirements
deltad = 0.5; %Labour requirements
deltau = 1; %Labour requirements
wage = 1; %Real wage
lambda = 1; %Partner choice
m = 0.1; %Market participant visibility
umin = 0.5; %Lower bound for stochastic consumer good price variation

%% Agent initialization
%D firms
Ad = ones(T,D); %Initial net worth of D firms set to one
Yd = zeros(T,D); %Scale of production of D firms
Nd = zeros(T,D); %Labour demand for U firms
Qd = zeros(T,D); %Intermediate goods requirements of D firms
Bd = zeros(T,D); %Loan demand of D firms
Ld = zeros(T,D); %Leverage of D firms
Rud = zeros(T,D); %Interest rate on trade credit for D firms
Rbd = zeros(T,D); %Interest rate on loans for D firms
PId = zeros(T,D); %Profit of D firms

%U firms
Au = ones(T,U); %Initial net worth of U firms set to one
Nu = zeros(T,U); %Labour demand for U firms
Qu = zeros(T,U); %Scale of production of U firms
Bu = zeros(T,U); %Loan demand of U firms
Lu = zeros(T,U); %Leverage of U firms
Rbu = zeros(T,U); %Interest rate on loans for U firms
PIu = zeros(T,U); %Profit of U firms

%B banks
Ab = ones(T,B); %Initial net worth of B banks set to one
PIb = zeros(T,B); %Profit of B banks

% Partner networks
BU = zeros(T,U);
BU(1,:) = randi([1 B],1,U); %Network between U firms and B banks (each col one U, number represents index of B bank)
    
BD = zeros(T,D);
BD(1,:) = randi([1 B],1,D); %Network between D firms and B banks (each col one D, number represents index of B bank)

UD = zeros(T,D);
UD(1,:) = randi([1 U],1,D); %Network between U firms and D firms (each col one D, number represents index of U firm)

%% Aggregate variables
%YD = sum(Yd); %Aggregate D firm production
%YU = sum(Qu); %Aggregate U firm production

%% Main programm time step
for s = 1:T
   %% D firms: Total 9 variables
   Yd(s,:) = (Ad(s,:) .^ beta) .* phi;
   Qd(s,:) = Yd(s,:) .* gamma;
   Nd(s,:) = Yd(s,:) .* deltad;
   Bd(s,:) = Nd(s,:) .* wage - Ad(s,:);
   Ld(s,:) = Bd(s,:) ./ Ad(s,:);
   %Rud(s,:) = (network_partners(Au(s,:), UD(s,:), D) .^ (alpha*-1)) .* alpha + (Ld(s,:) .^ alpha) .* alpha;
   %Rbd(s,:) = (network_partners(Ab(s,:), BD(s,:), D) .^ (alpha*-1)) .* alpha + (Ld(s,:) .^ alpha) .* alpha;
   
   %% U firms: Total 7 variables
   Qu(s,:) = network_worth(Yd(s,:), UD(s,:), U, D) .* gamma;
   Nu(s,:) = Qu(s,:).* deltau;
   Bu(s,:) = Nu(s,:) .* wage - Au(s,:);
   Lu(s,:) = Bu(s,:) ./ Au(s,:);
   %Rbu(s,:) = (network_partners(Ab(s,:), BU(s,:), U) .^ (alpha*-1)) .* alpha + (Lu(s,:) .^ alpha) .* alpha;
   
   %% Partner choice for s+1
   %UD(s+1,:) =
   %BU(s+1,:) =
   %BD(s+1,:) =
   
   %% Profit calculations
   %PId(s,:) = u .* Yd(s,:) - (1+Rbd(s,:)) .* Bd(s,:) - (1+Rud(s,:)) .* Qd(s,:)
   %PIu(s,:) = 
   %PIb(s,:) =
   
   %% Net worth s+1 calculation
   %Ad(s+1,:) = Ad(s,:) + PId(s,:);
   %Au(s+1,:) = Au(s,:) + PIu(s,:);
   %Ab(s+1,:) = Ab(s,:) + PIb(s,:);
end

save ABM_Replicator