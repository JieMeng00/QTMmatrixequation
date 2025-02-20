dx=0.1;
fm=exp(-(0:dx:log(1/sqrt(eps))).^2);

F=cqt('hankel',fm)+0.1*cqt(fm,fm);

A=cqt([-2 1], [-2 1]);

dt=0.15;
timesteps=2;

M=dx^2/2*cqt(1,1)-dt*A;

[Mn,Mp]=symbol(M);
Im=cqt(Mn(1),Mn(1));


B=Im-M;

beta=Mn(1);

X0=0*cqt(1,1);

%% computing square root using Binomial iteration
tic;
[Sroot,err1] = sroot(B/beta,X0);
t1=toc;



%NSroot = newtonsroot(M);

%% computing square root using CR

tic;
[CSroot, err2] =crsroot(M);
t2=toc;

 
 
 
err1;
err2;
t1
t2

 



%figure; semilogy(abs(sn));%title('Symbol of the Karcher mean');
%figure; semilogy(abs(sp));%title('Symbol of the Karcher mean');
%figure; C= correction(Sroot); mesh(log10(abs(C(1:end,1:end))));
%title('Correction of the Karcher mean');
