%% Q generator: for denser Q plus and sparse Q minus for the problem

function [Q] = Q_neg_generator(n,x)

%while 1
R = exprnd(x,[n,n,n]);
R(R>1) = 1;
Q = round(R);

% Q1 = squeeze(sum(Q_t,1));
% Q2 = squeeze(sum(Q_t,2));
% Q3 = squeeze(sum(Q_t,2));
% 
% if (all(Q1(:) >0)) && (all(Q2(:) >0)) && (all(Q3(:) >0))
%    Q = Q_t;
%    break;
% end
   
%end