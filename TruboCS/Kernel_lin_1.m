%--------------------------------------------------------------------------
% First kind of piece-wise linear function. See
% "Near Optimal compressed sensing without priors: parametric SURE approximate
% message passing"

% F:        3 * N
% F_div:    3 * 1
% y:        N * 1
%--------------------------------------------------------------------------
function [F,F_div] = Kernel_lin_1(y, v)


alpha_1 = 2 * sqrt(v);      % see (20) in the above paper
alpha_2 = 4 * sqrt(v);

if(alpha_2 / alpha_1 ~= 2 )
%     fprintf('the denoisers are different from the paper if alpha_2/alpha_1!=2\n');
end

N = length(y);
y = reshape(y, 1, N);

%---------------------------------------------
% The endpoints of the line segments
% See the figure
%---------------------------------------------

index_1 = find( y < - alpha_2);
index_2 = find( y > -alpha_2 & y < -alpha_1 );
index_3 = find( y > -alpha_1 & y < alpha_1 );
index_4 = find( y > alpha_1 & y < alpha_2 );
index_5 = find( y > alpha_2 );


if v==0
    F_div(1) = 0;
    F_div(2) = 0;
    F_div(3) = 1;
    F(1,:) = zeros(1,N);
    F(2,:) = zeros(1,N);
    f3 = zeros(1,N);
    f3(index_1) = y(index_1) + alpha_2;
    f3(index_5) = y(index_5) - alpha_2;
    F(3,:) = f3;
else
% function 1
f1 = zeros(1,N);
f1(index_2) = -1/alpha_1 * y(index_2) - 2;
f1(index_3) = 1/alpha_1 * y(index_3);
f1(index_4) = -1/alpha_1 * y(index_4) + 2;

% function 2
f2 = zeros(1,N);
f2(index_1) = -1;
f2(index_2) = ( y(index_2) + alpha_1 ) / ( alpha_2 - alpha_1 );
f2(index_3) = 0;
f2(index_4) = ( y(index_4) - alpha_1 ) / ( alpha_2 - alpha_1 );
f2(index_5) = 1;

% function 3
f3 = zeros(1,N);
f3(index_1) = y(index_1) + alpha_2;
f3(index_5) = y(index_5) - alpha_2;

F(1,:) = f1;
F(2,:) = f2;
F(3,:) = f3;

%-----------------------------------------------------
% Divergence, namely, the sum of derivatives
%-----------------------------------------------------
    F_div(1) = 1/N * ( -1/alpha_1 * length(index_2) + 1/alpha_1 * length(index_3) -1/alpha_1 * length(index_4) );
    F_div(2) = 1/N * ( 1/( alpha_2 - alpha_1 ) * length(index_2) + 1/( alpha_2 - alpha_1 ) * length(index_4) );
    F_div(3) = 1/N * ( length(index_1) + length(index_5) );
end



