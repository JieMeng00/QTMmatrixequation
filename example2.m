addpath /Users/littleapple/Documents/MATLAB/cqt-toolbox

n=4;
delta=0.5;
a=rand(n,1);
b=a+delta*rand(n,1);
c=rand(n,1);

M=[a b c];

W=M/(sum(sum(M)));


a = W(:,1);
b = W(:,2);
c = W(:,3);

B=cqt(b(1),b); % construction of matrix B

m = floor(n/2);

a1 = flipud(a(1:m)); % a1 and a2 are vecotrs that consist of  coefficnets of a(z)
a2 = a(m:end);

c1 = flipud(c(1:m+1)); %  c1 and c2 are the vectors that consist of coefficnets of c(z)
c2 = c(m+1:end);


corA = zeros(m-1,1); % correction part of A
corC = zeros(m,1); % correction part of C

for j = 1:m-1
    corA(j) = sum(a1(j+1:end));
end


for j = 1:m
    corC(j) = sum(c1(j+1:end));
end

A = cqt(a1,a2,corA); 
C = cqt(c1,c2,corC);


tic;
[G1,r1] = fixedpoint(A,B,C);
t1=toc;

tic;
[G3,r3] = CR0(A,B,C);
t2=toc;


