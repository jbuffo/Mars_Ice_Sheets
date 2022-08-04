%%% Jacob Buffo
%%% GTid 903073891

%% Function for Thomas algorithm, input a,b,c and y
%% vectors, output x vector
function x=Thomas_Trid(a,b,c,y)

%% Adding zeros into the b and c vectors
b=[b 0];
c=[0 c];

%% Creating P_1 and Q_1 in P&Q vectors that will be modified
%% in the for loop
n=length(y);
P=[-b(1)/a(1)];
Q=[y(1)/a(1)];

%% Loop over 2 to n-1 to solve for P_i and Q_i
for i=2:n;
    P(end+1)=-b(i)/(a(i)+c(i)*P(i-1));
    Q(end+1)=(y(i)-c(i)*Q(i-1))/(a(i)+c(i)*P(i-1));
end;

x=[Q(n)];
%% Loop over n to solve for x_i vector
for i=n-1:-1:1;
    xnew=P(i)*x(end)+Q(i);
    x(end+1)=xnew;
end;

%% Solution vector x_i
x=fliplr(x)';
    