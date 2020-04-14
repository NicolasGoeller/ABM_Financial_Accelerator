clear

B = 4;
U = 10;

BU = randi([1,B], 1, U);
Au = randi([1, 20], 1, U);
Au([2 5 9]) = -1;
Bu = randi([1,40], 1, U);

bad_debt(Au, Bu, BU, B)