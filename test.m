clear

%D = 10;
%U = 4;
%a = randi([1,10],1,4);
%b = randi([1,10],1,4);
%c = a * b;
%d = a .* b
%UD = randi([1,U],1,10);

%Ad = randi([1,15],1,10);
%s = [1,4,7];
%Ad(s) = -1;

%Bd = randi([1,15],1,10);



%x = bad_debt(Ad, Bd, UD, U);
%y = bad_debt2(Ad, Bd, UD, U);

UD = zeros(5,100);
UD(2,:) = randi([1 50],1,100); %Network between U firms and D firms (each col one D, number represents index of U firm)
