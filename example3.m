
lambda=0.2;
mu=0.4;

a1n=[0 lambda];
a1p=0;
A1=cqt(a1n,a1p);

a0n=[0 mu];
a0p=0;

A0=cqt(a0n,a0p,lambda);

am1n=0;
am1p=[0 mu];

Am1=cqt(am1n,am1p,mu);

I=cqt(1,1);

A=A1;
B=Am1;
C=A0;
Q=I-A0;


[W,r] = fixedpoint(A,B,C);
[w1,r1] = fixedpoint2(A,B,C);


G=W^(-1)*Am1;
R=A1*W^(-1);

[rn, rp]=symbol(R);
[gn, gp]=symbol(G);


figure;

x1=1:length(rn);
x=-x1;
semilogy(x, abs(rn),'b','linewidth',1.5);
hold on

semilogy(abs(rp),'b','linewidth',1.5);
hold off
xlabel('$i$','Interpreter','latex','fontsize',24)
ylabel('$\log10(|r_i|)$','Interpreter','latex','fontsize',24)


figure; C= correction(R); mesh(log10(abs(C(1:end,1:end))));


