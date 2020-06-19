clear

%D = 10;
%U = 4;
%a = randi([1,10],1,4);
%b = randi([1,10],1,4);
%c = a * b;
%d = a .* b

UD = randi([1,10],1,15);
Ad = randi([1,10],1,15);
s = [1,4,7,12];
Ad(s) = -1;
Bd = randi([1,10],1,15);
r = rand(1,15);
r(s) = 0.01;

BD = bad_debt(Ad, Bd, r, UD, 10)

%x = bad_debt(Ad, Bd, UD, U);
%y = bad_debt2(Ad, Bd, UD, U);

%mc = 1

%file = ['test', num2str(mc), '.mat']