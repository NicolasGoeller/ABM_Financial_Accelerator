clear

T = 3;
U = 5;
D = 10;
Au = zeros(T,U);
Au(1,:) = randi([5, 20], 1, U);

UD = zeros(T,D);
UD(1,:) = randi([1 U],1,D);

UD(1,:)
Au(1,:)
network_partners(Au(1,:), UD(1,:), D)


