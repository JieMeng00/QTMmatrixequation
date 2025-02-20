function [H,Bx,By,Bo] = flipxy(H, Bx, By, Bo)
% function [H,Bx,By,Bo] = flip(H, Bx, By, Bo)
% This function exchanges the x and he y axis in a bi-dimensional
% random walk
% H: 3x3 matrix with the probabilities in the inner part
% Bx: 2x3 matrix with the probabilities in the x boundary
% By: 3x2 matrix with the probabilities in the y boundary
% Bo: 2x2 matrix with the probabilities in the origin

H = H(:,end:-1:1);
H = H';
H = H(:,end:-1:1);

S = Bx(:,end:-1:1);
S = S';
S = S(:,end:-1:1);

T = By(:,end:-1:1);
T = T';
T = T(:,end:-1:1);

Bx = T;
By = S;

Bo = Bo(:,end:-1:1);
Bo = Bo';
Bo = Bo(:,end:-1:1);


