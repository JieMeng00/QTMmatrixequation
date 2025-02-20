function r=myroots(a,b,c)
% function r=myroots(a,b,c)
% r(i) is the solution of minimum modulus of the equation
% a(i)x^2+b(i)x+c(i)=0
% a,b,c column vectors 

d = sqrt(b.^2-4*a.*c);
a2 = 2*a;
x1 = (-b+d)./a2;
x2 = (-b-d)./a2;
A = [abs(x1),abs(x2)]';
[val,pos]=max(A);
r = x1.*(pos==1)'+x2.*(pos==2)';
r = c./(a.*r);
idx = find(a==0); k = length(idx);
for i=1:k
    ii = idx(i);
    r(ii)=-c(ii)/b(ii);
end


