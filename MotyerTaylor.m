function [Am1,A0,A1,B0,B1,gm,gp]=MotyerTaylor(j,epsilon)
% function [Am1,A0,A1,B0,B1,gm,gp]=MotyerTaylor(j,epsilon)
% provides the blocks of the j-th QBD problem in the
% Motyer-Taylor tests, where j=1,2,...,10
% perturbed with a relative perturbation at most epsilon
% more precisely, the parametrs mu and lambda are peturbed as
% lambda*(1+epsilon*(1+rand)), mu*(1+epsilon*(1+rand))
% In output Am1,A0,A1 are the matrix coefficients while gm,gp are
% the coefficient vectors of the negative and positive powers of g(z) 
% Notation in the paper and standard notation
% QT1 QT0           B0  B1
% Q2  Q1  Q0     =  Am1 A0  A1
%     Q2  Q1 Q0         Am1  A0 A1

% The function construct matrices H, Bx,By,Bo such that
% H: contains the Toeplitz part of 
%      Q0 (first row), Q1 (second row), Q2 (last row)
% By: contains the boundary conditions (first row) of
%      Q0 (first row), Q1 (second row), Q2 (last row)
% Bx: contains the Toeplitz part of
%      QT0 (first row), QT1 (second row)
% Bo: contains the boundary conditions (first row) of
%      Qt0 (first row), QT1 (second row)
% and then invokes the function qbd_blocks

% Correspondence with Leonardo list
% Leo: 1 2 3 4 5 6 7 8 9 10
% MT:  4 3 2 1 6 5 7 8 10 9

MT=[1 0 1.5 2 1 0;
    1 0 2 1.5 1 0;
    0 1 1.5 2 0 1;
    0 1 2 1.5 0 1;
    1 1 2 2 0.1 0.8;
    1 1 2 2 0.8 0.1;
    1 1 2 2 0.4 0.4;
    1 1 10 10 0.5 0.5;
    1 5 10 15 0.4 0.9;
    5 1 15 10 0.9 0.4];
% Q2=Am1 cqt([(1-q)mu2 0],[(1-q)mu2 qmu2]);
% Q1=A0  cqt([-l1-l2-mu1-mu2 (1-p)mu1],[-l1-l2-mu1-mu2, l1],[mu1 0; 0 0]);
% Q0=A1  cqt([l2 pmu1], [l2 0])
% QT1=B0  cqt([-l1-l2-mu1 (1-p)mu1], [-l1-l2-mu1 l1],[mu1 0;0 0]);
% QT0=B1 = A1;

l1=MT(j,1);
l2=MT(j,2);
mu1=MT(j,3);
mu2=MT(j,4);

mu1=mu1*(1+epsilon*(1+rand));
mu2=mu2*(1+epsilon*(1+rand));
l1=l1*(1+epsilon*(1+rand));
l2=l2*(1+epsilon*(1+rand));

p=MT(j,5);
q=MT(j,6);
H = [p*mu1            l2        0;    ...
       (1-p)*mu1  -l1-l2-mu1-mu2  l1;  ...
        0          (1-q)*mu2    q*mu2];      % Toeplitz part of Q0,Q1,Q2
Bx = [p*mu1            l2        0;   ...
        (1-p)*mu1  -l1-l2-mu1  l1];          % Qtilde0, Qtilde1 Toeplitz part
By = [l2 0; -l1-l2-mu2 l1; (1-q)*mu2 q*mu2 ]; % Boundary vaues in Q0,Q1,Q2
Bo = [l2 0; -l1-l2 l1]; % Boundary values in Qtilde0, Qtilde 1

[Am1, A0, A1, B0, B1] = prob2qbd(H,Bx,By,Bo);


% normalization so that Am1+A0+A1 is stochastic
 alpha = l1+l2+mu1+mu2;
 %nrm=-A0(10,10); nrm = max(abs(A0(1,1)),abs(A0(2,2)));
 Am1=Am1/alpha; A0=A0/alpha+cqt(1,1,0); A1=A1/alpha; 
 B0=B0/alpha +cqt(1,1,0);B1=B1/alpha;

% for problems 2 and 6 consider the flipped version
if j==2 || j==6 
   [Am1, A0, A1, B0, B1] = flipqbd(Am1, A0, A1, B0, B1);
end


[am1m,am1p]=symbol(Am1);
[a1m,a1p]=symbol(A1);
[a0m,a0p]=symbol(A0);
%length(a0m)
if length(a0m)==1
    tmp=a0m;
    a0m=[tmp,0];
end

tic
disp('Computing g(z)...')
[gm,gp]=rg_eval_interpari(am1m,am1p,a0m,a0p,a1m,a1p);
time=toc;
fprintf('cpu time for computing gm,gp: %d seconds\n',time);
%lengthgm=length(gm)
%lengthgp=length(gp)
