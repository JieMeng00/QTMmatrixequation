function [gm,gp]=interpari(am1m,am1p)
%function[gm,gp]=inter(am1m,am1p) 
%[am1m,ap1p]=symbol(A);
%A=(I-B)^2;
%compute b, where b satisfies
% b^2-2b+1=a


if sum(am1m)+sum(am1p)==0
  gm=0;gp=0;
  return
end

lm1m = length(am1m); lm1p = length(am1p);

  % compute a(1)
am1=sum(am1m)+sum(am1p)-am1m(1);
 % compute first derivatives
 
% compute the first derivative at 1, negative part
 sm=0;
 for i=1:lm1m-1
     sm=sm-i*am1m(i+1);
 end
 
 % compute the first derivative at 1, positive part 
  sp=0;
 for i=1:lm1p-1
     sp=sp+i*am1p(i+1);
 end
     
 
dam1=sp+sm;

  % compute second derivatives
  ddsp=0;
  for i=1:lm1p-1
      ddsp=ddsp+(i-1)*i*am1p(i+1);
  end
  
  ddsm=0;
  for i=1:lm1m-1
      ddsm=ddsm+(i+1)*i*am1m(i+1);
  end
  
ddam1=ddsp+ddsm;

  % compute a'(1)
denom = 2*(-sqrt(am1));

if denom==0
  disp('Infinite derivative')
return
end 
dg=dam1/denom;

  % compute g''(1)
ddg=(ddam1-2*dg^2)/denom;
  % lengths of the input
  
lm1m = length(am1m); lm1p = length(am1p);

  % starting values
n=16; debug=true;
delta=1.d100;deltaold=1.d200;
epsi=1.d-14;
cnt=0;
g=ones(4*n,1);gg=g;
%while abs(delta)>epsi*n^2**ddg  && delta<deltaold % ||2>1

while abs(delta)>epsi*ddg*n  && abs(delta)<abs(deltaold) % ||2>1
    cnt=cnt+1;
  n=2*n; N=2*n; med=n;
  if debug
    fprintf('N=%d, delta/n^2=%d\n',N,delta/n^2);
  end
  am1 = zeros(N,1);
  am1(1:lm1p) = am1p; am1(end:-1:end-lm1m+2) = am1m(2:end);
 

  fam1 = fft(am1); 
  fx = fam1; 
  fy = fx;
  fx=roots_b(fam1);
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
  
