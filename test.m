clear

T = 100;
U = 5;
D = 10;
Au = ones(T,U);
%Au(1,:) = randi([5, 20], 1, U);

UD = zeros(T,D);
UD(1,:) = randi([1 U],1,D);

UD(1,:)
Au(1,:)
network_partners(Au(1,:), UD(1,:), D)


UD = zeros(T,D);
UD(1,:) = randi([1 U],1,D);

for s = 1:3
    size(UD(s,:))
end

for s = 1:T
    s
end