function [G,r] = fixedpoint2(A,B,C)
% TO COMPUTE THE SOLUTION OF THE EQUATION X+AX^{-1}B=Q; Q=I-C
% X=I-G;
verb = true; maxit = 100; epsi = 1.e-14; cqtoption('threshold',10^(-15));
I = cqt( 1, 1 );
X = 0 * cqt(1,1);   
err = 1;
r = zeros( maxit, 1 );

for k = 1:maxit
    Xold = X;
       X = C + A * (I-Xold)^(-1) * B;
  errold = err;
     err = norm( -X + A * (I-X)^(-1) * B + C, inf );
     
    if verb
         fprintf( 'step=%d, err=%d\n', k, err ); 
    end
    
     r(k)   = err;
    
     if err < epsi || (err - errold > 0 && k > 1), break; end
   
end
G = X;
if (k == maxit)
fprintf( 'Warning: reached the max number of iterations' );
end