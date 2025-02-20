function [Am1,A0,A1,B0,B1] = flipqbd(Am1,A0,A1,B0,B1);
% function [Am1,A0,A1,B0,B1] = flipqbd(Am1,A0,A1,B0,B1);
% compute the QBD blocks of the flipped problem
% obtained by exchanging the x and the y axis

[H,Bx,By,Bo] = qbd2prob(Am1,A0,A1,B0,B1);
[H,Bx,By,Bo] = flipxy(H,Bx,By,Bo);
[Am1,A0,A1,B0,B1] = prob2qbd(H,Bx,By,Bo);
