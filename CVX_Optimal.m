function theta_Opt = CVX_Optimal( G, Hd, Dh, rho, N)


R = zeros(N+1 , N+1);
R(1:N , :) = [ (rho^2 * Dh'* (G'* G) * Dh), (rho * Dh'* G'* Hd) ];
R( N+1 , :) = [ (rho * Hd'* G * Dh)  , 0 ];

V = zeros(N+1 , N+1);
V(1:N , 1:N) =  rho*(1-rho)* diag(diag( (rho^2 * Dh'* (G'* G) * Dh) ));

Pi = R + V ;

%%
cvx_begin  sdp quiet
    variable Q(N+1,N+1) complex 
    expression QQ(N+1, N+1)
    QQ = ( Q + Q' )./2;
    maximize( real(trace( Pi * QQ )) )
    subject to
        QQ >= 0;
        vec( diag(QQ) ) == vec( ones(N + 1, 1) );        
cvx_end
[Uq, Dq, ~] = svd(QQ);

A_Opt = 0;
for iter = 1 : 1  
    Rq = sqrt(0.5)*(randn(N+1, 1) + 1i*randn(N+1, 1));
    theta = Uq * sqrt(Dq) * Rq ;
    theta = theta./(theta(N+1));
    theta = theta(1:N)./(abs(theta(1:N)));
    A_ite = norm( G * Dh * theta , 'fro') ;
    if A_ite > A_Opt
       A_Opt = A_ite;
       theta_Opt = theta ;
    end  
end

end







