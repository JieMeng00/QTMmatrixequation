function [G,r] = CR0(A,B,C)
% TO COMPUTE THE SOLUTION OF THE EQUATION X+AX^{-1}B=Q;Q=I-C
verb = true; maxit = 100; epsi = 1.e-14; cqtoption('threshold',10^(-15));
I = cqt( 1, 1 );
Q=I-C;
X = Q;  % 初始值X0=I
err = 1;
r = zeros( maxit, 1 );

for k = 1:maxit
    Xold   =  X;
    Aold   =  A;
    Bold   =  B;
    Qold   =  Q;

    Wold   =  Aold * Qold^(-1);
    X      =  Xold - Wold * Bold;
    A      =  Wold * Aold;

     
    errold = err;

    err=norm(A,inf);
    
    if verb
         fprintf( 'step=%d, err=%d\n', k, err ); 
    end
    
 
     r(k)   = err;
    
     if err < epsi || (err - errold > 0 && k > 1), break; end

     B = Bold * Qold^(-1) * Bold;
     
     Q = Qold-Bold * Qold^(-1) * Aold-Wold * Bold;
         
    
    

end
G = X;
if (k == maxit)
fprintf( 'Warning: reached the max number of iterations' );
end