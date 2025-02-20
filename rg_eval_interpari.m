function [gm,gp]=rg_eval_interpari(am1m,am1p,a0m,a0p,a1m,a1p)
%function[gm,gp]=rg_eval_inter(am1m,am1p,a0m,a0p,a1m,a1p)
% a1*g^2+(a0-1)*g+am1=0


if sum(am1m)+sum(am1p)==0
  gm=0;gp=0;
  return
end


  % compute am1(1) and a1(1)
am1=sum(am1m)+sum(am1p)-am1m(1);
a1=sum(a1m)+sum(a1p)-a1m(1);
a0=sum(a0m)+sum(a0p)-a0m(1);
g1=min(1,am1/a1);
 % compute first derivatives
dam1=sum(am1p)-sum(am1m);
da1=sum(a1p)-sum(a1m);
da0=sum(a0p)-sum(a0m);
  % compute second derivatives
ddam1=2*am1m(2);
dda1=2*a1m(2);
dda0=2*a0m(2);
  % compute g'(1)
denom = (1-a0-2*a1*g1);
if denom==0
  disp('Infinite derivative')
return
end 
dg=(dam1+da0*g1+da1*g1^2)/denom;

  % compute g''(1)
ddg=ddam1+dda0*g1+dda1*g1^2+2*a1*dg^2+dg*(4*g1*da1+2*da0);
ddg=ddg/(1-a0-2*a1*g1);
  % lengths of the input
lm1m = length(am1m); lm1p = length(am1p);
l0m = length(a0m); l0p = length(a0p);
l1m = length(a1m); l1p = length(a1p);

  % starting values
n=16; debug=true;
delta=1.d100;deltaold=1.d200;
epsi=1.d-16;
cnt=0;
g=ones(4*n,1);gg=g;
%while abs(delta)>epsi*n^2**ddg  && delta<deltaold % ||2>1

while abs(delta)>epsi*ddg*n  && abs(delta)<abs(deltaold) % ||2>1
    cnt=cnt+1;
  n=2*n; N=2*n; med=n;
  if debug
    fprintf('N=%d, delta/n^2=%d\n',N,delta/n^2);
  end
  am1 = zeros(N,1); a0 = am1; a1 = am1;
  am1(1:lm1p) = am1p; am1(end:-1:end-lm1m+2) = am1m(2:end);
  a0(1:l0p) = a0p; a0(end:-1:end-l0m+2) = a0m(2:end);
  a1(1:l1p) = a1p; a1(end:-1:end-l1m+2) = a1m(2:end);

  fam1 = fft(am1);  a0(1)=a0(1)-1; fa0 = fft(a0);  fa1 = fft(a1);
  fx = fam1; 
  fy = fx;
  fx=myroots(fa1,fa0,fam1);
  g = ifft(fx);
  gm=zeros(n,1);gp=zeros(n+1,1);
  gp = g(1:n+1); gm(1) = g(1);
  gm(2:n) = g(end:-1:n+2);
  wm=[0:n-1].*([0:n-1]+1);wm=wm';
  wp=[0:n].*([0:n]-1);wp=wp';
  deltaold=delta;
  delta=sum(gm.*wm)+sum(gp.*wp);
  delta=ddg-delta;
  M=N/2;
  semilogy(abs(g(1:M/2)-gg(1:M/2)))
  gg=g;
%[abs(delta),epsi*ddg*n,deltaold]  
end
  figure; plot(fx,'.')

%clean
nn=min(find(gm<0));
if isempty(nn)
    nn=length(gm)+1;
end
gm=gm(1:nn-1);
nn=min(find(gp<0));
if isempty(nn)
    nn=length(gp)+1;
end
gp=gp(1:nn-1);

%nn=min(find(rm<0));
%if isempty(nn)
%    nn=length(rm)+1;
%end
%rm=rm(1:nn-1);
%nn=min(find(rp<0));
%if isempty(nn)
%    nn=length(rp)+1;
%end
%rp=rp(1:nn-1);


return

%clean
ck=gm(32:end)<gm(31:end-1);nn=min(find(ck==0));
if isempty(nn)
   nn=length(gm)+1;
end
gm=gm(1:nn-1);
ck=gp(32:end)<gp(31:end-1);nn=min(find(ck==0));
if isempty(nn)
   nn=length(gp)+1;
end
gp=gp(1:nn-1);
return


   for i=30:length(gm)
     if gm(i)>gm(i-1)
        gm=gm(1:i-1);
        break
     end
   end

   for i=30:length(gp)
     if gp(i)>gp(i-1)
        gp=gp(1:i-1);
        break
     end
   end
  if size(gm,1)>1
    gm=gm.';
  end
  if size(gp,1)>1
    gp=gp.';
  end


% residual
if debug
  norm(res)
  norm(ifft(res),inf)
end
  
