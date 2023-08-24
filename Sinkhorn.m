%% Negative Probabilities in High Dimensions

% make up a simple example with a three dimensional signed tensor, 
% that you can color code in 3-D as a signed distribution,
% with blue negative and red positive. Split it two, 
% the positive Q^+ and Q^-, as in our previous paper. 
% Now these are not transition probabilities,
% they are probabilities but indefinite, 
% so that Q=Q^+-Q^-. Then, compute marginals. 
% Btw, Q^- has to have fewer elements so that the marginals are positive,
% and also, Q^+ is "richer". 
% Then, compute marginals, perturb them a bit, 
% and iterate as we have in the previous paper.

clear all
clc

addpath('function')

%% Generate Q+ and Q-

n = 3;
x1 = 1;
x2 = 0.2;

Q_pos = Q_generator(n,x1);
Q_neg = Q_neg_generator(n,x2);





