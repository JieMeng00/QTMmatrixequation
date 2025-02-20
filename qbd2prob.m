function [H,Bx,By,Bo] = qbd2prob(Am1,A0,A1,B0,B1)
% function [H,Rx,Ry,Ro] = qbd2prob(Am1,A0,A1,B0,B1)
% from qbd representation to probability representation
% Input: the blocks defining the QBD
% Bo B1
% Am1 A0  A1
%     Am1 A0 A1
% Output: the probabilities of transition
% H : probabilities in the inner part
% Bx: probabilities along the x axis
% By: probabilities along the y axis
% Bo: probabilities at the origin

H=zeros(3);
[an,ap] = symbol(A1);
if length(an)==1
  an = [an,0];
end
if length(ap)==1
  ap = [ap,0];
end
H(1,:)=[an(2),an(1),ap(2)];
%
[an,ap] = symbol(A0);
if length(an)==1
  an = [an,0];
end
if length(ap)==1
  ap = [ap,0];
end
H(2,:)=[an(2),an(1),ap(2)];
%
[an,ap] = symbol(Am1);
if length(an)==1
  an = [an,0];
end
if length(ap)==1
  ap = [ap,0];
end
H(3,:)=[an(2),an(1),ap(2)];

%By
c = correction(A1);
if size(c,2)==1
  c = [c 0];
elseif size(c,2)==0
  c=[0 0];
end
c = c(1,:);
By(1,:) = c+H(1,2:3);
%
c = correction(A0);
if size(c,2)==1
  c = [c 0];
elseif size(c,2)==0
  c=[0 0];
end
c = c(1,:);
By(2,:) = c+H(2,2:3);
%
c = correction(Am1);
if size(c,2)==1
  c = [c 0];
elseif size(c,2)==0
  c=[0 0];
end
c = c(1,:);
By(3,:) = c+H(3,2:3);
%
% Bx
[syb0n, syb0p] = symbol(B0);
[syb1n, syb1p] = symbol(B1);
if length(syb0n)==1
  syb0n=[syb0n,0];
end
%
if length(syb0p)==1
  syb0p=[syb0p,0];
end
%
if length(syb1n)==1
  syb1n=[syb1n,0];
end
%
if length(syb1p)==1
  syb1p=[syb1p,0];
end
Bx(2,:) = [syb0n(2), syb0n(1), syb0p(2)];
Bx(1,:) = [syb1n(2), syb1n(1), syb1p(2)];

%Bo
c1 = correction(B1);
c0 = correction(B0);
c1 = c1(1,:); c0 = c0(1,:);
if length(c1)==1
  c1 = [c1,0];
elseif length(c1)==0
  c1=[0, 0];
end
if length(c0)==1
  c0 = [c0,0];
elseif length(c0)==0
  c0=[0, 0];
end

Bo = [c1;c0]+[syb1p(1), syb1p(2);syb0p(1), syb0p(2)];

