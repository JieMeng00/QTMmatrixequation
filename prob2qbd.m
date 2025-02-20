function [Am1,A0,A1,B0,B1] = prob2qdb(H, Bx, By, Bo)
% function [Am1,A0,A1,B0,B1] = prob2qbd(H, Bx, By, Bo)
% This function generates the blocks Am1, A0, A1, B0, B1
% of the QBD 
% B0 B1
% Am1 A0  A1
%     Am1 A0 A1
% associated with the probabilities given in
% the matrices  H, Bx, By, Bo, where
% The blocks are QT matrices, The CQT-toolbox is needed 
% H: 3x3 matrix with the probabilities in the inner part
% Bx: 2x3 matrix with the probabilities in the x boundary
% By: 3x2 matrix with the probabilities in the y boundary
% Bo: 2x2 matrix with the probabilities in the origin
% more precisely,

% H: contains the Toeplitz part of 
%      A1 (first row), A0 (second row), Am1 (last row)
% By: contains the boundary conditions (first row) of
%      A1 (first row), A0 (second row), Am1 (last row)
% Bx: contains the Toeplitz part of
%      B1 (first row), B0 (second row)
% Bo: contains the boundary conditions (first row) of
%      B1 (first row), B0 (second row)


Am1 = cqt([H(3,2),H(3,1)], [H(3,2),H(3,3)], [By(3,1)-H(3,2), By(3,2)-H(3,3); 0,0]);
A0  = cqt([H(2,2),H(2,1)], [H(2,2),H(2,3)], [By(2,1)-H(2,2), By(2,2)-H(2,3); 0,0]);
A1  = cqt([H(1,2),H(1,1)], [H(1,2),H(1,3)], [By(1,1)-H(1,2), By(1,2)-H(1,3); 0,0]);

B0 = cqt([Bx(2,2),Bx(2,1)], [Bx(2,2), Bx(2,3)], [Bo(2,1)-Bx(2,2), Bo(2,2)-Bx(2,3);0,0] );
B1 = cqt([Bx(1,2),Bx(1,1)], [Bx(1,2), Bx(1,3)], [Bo(1,1)-Bx(1,2), Bo(1,2)-Bx(1,3);0,0]);
