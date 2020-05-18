%%% Replication of Credit Network model in dynamic partner selection
% The simulation runs under the data saving paradigm of each simulation
% time step represents one matrix row in output, while each agent has the
% same column number over all tables.

% Incorporate Re-Entry of substitutes

% Look into substitute search behaviour of banks and corporations to improve partner selection algorithm

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
u = zeros(T,D); %Initial consumer good prices
Ad = ones(T,D); %Initial net worth of D firms set to one
Yd = zeros(T,D); %Scale of production of D firms
Nd = zeros(T,D); %Labour demand for U firms
Qd = zeros(T,D); %Intermediate goods requirements of D firms
Bd = zeros(T,D); %Loan demand of D firms
Ld = zeros(T,D); %Leverage of D firms
Rud = zeros(T,D); %Interest rate on trade credit for D firms
Rbd = zeros(T,D); %Interest rate on loans for D firms
PId = zeros(T,D); %Profit of D firms
BRd = zeros(T,D); %Bankruptcy indicator of D firms

%U firms
Au = ones(T,U); %Initial net worth of U firms set to one
Nu = zeros(T,U); %Labour demand for U firms
Qu = zeros(T,U); %Scale of production of U firms
Bu = zeros(T,U); %Loan demand of U firms
Lu = zeros(T,U); %Leverage of U firms
Rbu = zeros(T,U); %Interest rate on loans for U firms
PIu = zeros(T,U); %Profit of U firms
BDu = zeros(T,U); %Bad debts of U firms
BRu = zeros(T,U); %Bankruptcy indicator of U firms

%B banks
Ab = ones(T,B); %Initial net worth of B banks set to one
PIb = zeros(T,B); %Profit of B banks
BDb = zeros(T,B); %Bad debt of B banks
BRb = zeros(T,B); %Bankruptcy indicator of B banks

% Partner networks
BU = zeros(T,U);
BU(1,:) = randi([1 B],1,U); %Network between U firms and B banks (each col one U, number represents index of B bank)
    
BD = zeros(T,D);
BD(1,:) = randi([1 B],1,D); %Network between D firms and B banks (each col one D, number represents index of B bank)

UD = zeros(T,D);
UD(1,:) = randi([1 U],1,D); %Network between U firms and D firms (each col one D, number represents index of U firm)


%% Main programm time step

for s = 1:T
   %% Consumer good price setting
   u(s,:) = ((2-umin)-umin) .* rand(D,1) + umin;
   
   %% D firms:
   Yd(s,:) = (Ad(s,:) .^ beta) .* phi;
   Qd(s,:) = Yd(s,:) .* gamma;
   Nd(s,:) = Yd(s,:) .* deltad;
   Bd(s,:) = Nd(s,:) .* wage - Ad(s,:);
   Bd(s,Bd(s,:)<0) = 0; %Seems that most values at T2 are 0 because firms are self-financed
   %for i = 1:D
   %    if Bd(s,i) < 0
   %        Bd(s,i) = 0;
   %    end
   %end
   
   Ld(s,:) = Bd(s,:) ./ Ad(s,:);
   Rud(s,:) = (Au(s,UD(s,:)) .^ (alpha*-1)) .* alpha + (Ld(s,:) .^ alpha) .* alpha;
   Rbd(s,:) = (Ab(s,BD(s,:)) .^ (alpha*-1)) .* alpha + (Ld(s,:) .^ alpha) .* alpha;
   
   PId(s,:) = u(s,:) .* Yd(s,:) - (1+Rbd(s,:)) .* Bd(s,:) - (1+Rud(s,:)) .* Qd(s,:);
   Ad(s+1,:) = Ad(s,:) + PId(s,:);
   BRd(s,Ad(s+1,:)<0) = 1;
   
   
   %% U firms:
   Qu(s,:) = network_worth(Yd(s,:), UD(s,:), U) .* gamma;
   Nu(s,:) = Qu(s,:).* deltau;
   Bu(s,:) = Nu(s,:) .* wage - Au(s,:);
   Bu(s,Bu(s,:)<0) = 0;

   Lu(s,:) = Bu(s,:) ./ Au(s,:);
   Rbu(s,:) = (Ab(s,BU(s,:)) .^ (alpha*-1)) .*alpha + (Lu(s,:) .^ alpha) .* alpha;

   PIu(s,:) = network_worth(((Rud(s,:) + 1) .* Qd(s,:)), UD(s,:), U) - ((Rbu(s,:) + 1) .* Bu(s,:));
   BDu(s,:) = bad_debt(Ad(s,:), Qd(s,:), UD(s,:), U);
   Au(s+1,:) = Au(s,:) + PIu(s,:) - BDu(s,:);
   BRu(s,Au(s+1,:)<0) = 1;
   
   
   %% B banks
   PIb(s,:) = network_worth(((Rbd(s,:) + 1) .* Bd(s,:)), BD(s,:), B) + network_worth(((Rbu(s,:) + 1) .* Bu(s,:)), BU(s,:), B);
   BDb(s,:) = bad_debt(Au(s,:), Bu(s,:), BU(s,:), B) + bad_debt(Ad(s,:), Bd(s,:), BD(s,:), B);
   Ab(s+1,:) = Ab(s,:) + PIb(s,:) - BDb(s,:);
   BRb(s,Ab(s+1,:)<0) = 1;
   
   
   %% Partner choice for s+1
   % Original model replace U partner randomly if it just went bankrupt
   % Does not replace bank partners at all if gone bankrupt
   UD(s+1,:) = partner_choice(Au(s,:), Ld(s,:), m, UD(s,:), lambda); 
   BU(s+1,:) = partner_choice(Ab(s,:), Lu(s,:), m, BU(s,:), lambda);
   BD(s+1,:) = partner_choice(Ab(s,:), Ld(s,:), m, BD(s,:), lambda);
   
   
   %% Replace bankrupt agents by new ones
   %Ad(s+1, BRd(s,:)==1) = 1;
   %UD(s+1, BRd(s,:)==1) = randi([1,U],1,1);
   %BD(s+1, BRd(s,:)==1) = randi([1,B],1,1);

   %Au(s+1, BRu(s,:)==1) = 1;
   %BU(s+1, BRu(s,:)==1) = randi([1,B],1,1);
   
   %Ab(s+1, BRb(s,:)==1) = 1;
   
end

%Cut the last row from networth tables & network documentation
Ad = Ad(1:1000,:);
Au = Au(1:1000,:);
Ab = Ab(1:1000,:);
UD = UD(1:1000,:);
BU = BU(1:1000,:);
BD = BD(1:1000,:);

save ABM_Replicator