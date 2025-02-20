I=cqt(1,1);

% for j=1,2,...,10, construction of the matrices A, B, C
j = 10;

[Am1,A0,A1,B0,B1,gm,gp] = MotyerTaylor(j,0);

if j == 10 % for problem 10, we consider the flipped version
   [Am1, A0, A1, B0, B1] = flipqbd(Am1, A0, A1, B0, B1);
end

A = A1;
B = Am1;
C = A0;
Q = I-C;

%% behaviour of the fixed-point iteration
tic;
[G1,r1] = fixedpoint(A,B,C);
t1=toc;
err1 = G1+A*G1^(-1)*B-Q;
tic;

%% behaviour of the fixed-point iteration
[G2,r3] = CR0(A,B,C);
t3=toc;
err2 = G2+A*G2^(-1)*B-Q;


% information of the Toeplitz part and correction part of the solution 

[G2n,G2p]=symbol(G2);
band=length(G2n)+length(G2p)-1;

Cor=correction(G2);
vs=size(Cor);
rows=vs(1);
columns=vs(2);
r=rank(Cor);
 