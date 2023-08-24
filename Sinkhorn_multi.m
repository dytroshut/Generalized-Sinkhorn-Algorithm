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
clf 

addpath('function')



%% Generate Q+ and Q-

n = 3;


% Q_pos_rd(:,:,1) = [1, 1, 0; 1, 0, 1; 1, 0, 1];
% Q_pos_rd(:,:,2) = [0, 0, 0; 1, 1, 0; 1, 1, 0];
% Q_pos_rd(:,:,3) = [1, 1, 1; 1, 0, 1; 0, 1, 1];
% 
% 
% Q_neg_rd(:,:,1) = [0, 0, 1; 0, 0, 0; 0, 0, 0];
% Q_neg_rd(:,:,2) = [1, 0, 0; 0, 0, 0; 0, 0, 0];
% Q_neg_rd(:,:,3) = [0, 0, 0; 0, 1, 0; 0, 0, 0];

Q_pos_rd(:,:,1) = [1, 1, 0; 1, 0, 1; 1, 0, 1];
Q_pos_rd(:,:,2) = [0, 0, 0; 1, 1, 0; 1, 1, 0];
Q_pos_rd(:,:,3) = [0, 1, 1; 1, 0, 1; 0, 0, 1];


Q_neg_rd(:,:,1) = [0, 0, 1; 0, 0, 0; 0, 0, 0];
Q_neg_rd(:,:,2) = [1, 0, 0; 0, 0, 0; 0, 0, 0];
Q_neg_rd(:,:,3) = [0, 0, 0; 0, 0, 0; 0, 1, 0];


% %%
% n = 10;   % works for the 3-D case
% x1 = 1;  % index for generating posotive elements 
% x2 = 0.2; % index for generating negative elements
% 
% 
% Q_pos_rd = Q_generator(n,x1);
% Q_neg_rd = Q_neg_generator(n,x2);

c_1 = [0.2;0.3;0.5];
c_2 = [0.4;0.4;0.2];
c_3 = [0.1;0.6;0.3];
%%
[num_node,~,~] = size(Q_pos_rd);

bQ = Q_pos_rd + Q_neg_rd;
bQ(bQ>=1) = 1;

% P = A;
%% (need to be sorted out, the Q need to be normalized as a given prior.)
Q = bQ;
Q(Q_neg_rd~=0) = -1;

Q_pos = (bQ + Q)./2;
Q_neg = (bQ - Q)./2;

X_plus = Q_pos;
X_minu = Q_neg;
X = X_plus - X_minu;

%% Three marginals

% x = (0:n-1)'/n-1;
% 
% Gaussian = @(x,t0,sigma)exp( -(x-t0).^2/(2*sigma^2) );
% normalize = @(p)p/sum(p(:));
% sigma = .06;
% 
% c_1 = Gaussian(x, .2, sigma); 
% c_2 = Gaussian(x, .5, sigma);
% c_3 = Gaussian(x, .8, sigma);
% 
% c_1 = normalize(c_1);
% c_2 = normalize(c_2);
% c_3 = normalize(c_3);


%% algorithm initialization

% intial value of variable \mu and \nu, set as all-ones vectors.
alpha_1 = ones(num_node,1);
alpha_2 = ones(num_node,1);
alpha_3 = ones(num_node,1);

% element scaling by exp(-1). % 3D tensor 
Qe_pos = Q_pos.*exp(-1);
Qe_neg = Q_neg.*exp(-1);

% empty vector for iteration tracking.
cost = [];
Err_1 = [];
Err_2 = [];
Err_3 = [];

% iteration
iter = 100;

%% so far so good!! redo the Sinkhorn iteration, with tensor structure (3D array)

%% forward-backward (sinkhorn-like) iteration
for t = 1:iter
%% 1. update of alpha1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute a,b,c for alpha_1  (to be changed)
a_23 = tensorprod(alpha_2,alpha_3,2);
Sp_1 = tensorprod(Qe_pos,a_23,[2,3],[1,2]);
Sm_1 = tensorprod(Qe_neg,1./a_23,[2,3],[1,2]);

% address the positive & negative elements, using location indicator
alpha_neg_loc_1 = (Sm_1~=0);
alpha_pos_loc_1 = ~alpha_neg_loc_1;

% update \nu_i in two cases, i.e.,
alpha_var_1 = zeros(num_node,1);
% 1) when a_{ij} < 0, compute \nu_i using function 'FB_solver'
alpha_var_1(alpha_neg_loc_1) = Gsolver(Sp_1(alpha_neg_loc_1),Sm_1(alpha_neg_loc_1),c_1(alpha_neg_loc_1));
% 2) when a_{ij} > 0, compute \nu_i using sinkhorn iteration
alpha_var_1(alpha_pos_loc_1) = c_1(alpha_pos_loc_1)./Sp_1(alpha_pos_loc_1); %%??
% 3) update \nu from the two cases
alpha_1 = alpha_var_1;


%% 2. update of alpha2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute a,b,c for alpha_2   (to be changed)
a_13 = tensorprod(alpha_1,alpha_3,2);
Sp_2 = tensorprod(Qe_pos,a_13,[1,3],[1,2]);
Sm_2 = tensorprod(Qe_neg,1./a_13,[1,3],[1,2]);

% address the positive & negative elements, using location indicator
alpha_neg_loc_2 = (Sm_2~=0);
alpha_pos_loc_2 = ~alpha_neg_loc_2;

% update \mu_j in two cases, i.e.,
alpha_var_2 = zeros(num_node,1);
% 1) when a_{ij} < 0, compute \mu_j using function 'FB_solver'
alpha_var_2(alpha_neg_loc_2) = Gsolver(Sp_2(alpha_neg_loc_2),Sm_2(alpha_neg_loc_2),c_2(alpha_neg_loc_2));
% 2) when a_{ij} > 0, compute \mu_j using sinkhorn iteration
alpha_var_2(alpha_pos_loc_2) = c_2(alpha_pos_loc_2)./Sp_2(alpha_pos_loc_2);
% 3) update \nu from the two cases
alpha_2 = alpha_var_2;


%% 3. update of alpha3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute a,b,c for alpha_3  (to be changed)
a_12 = tensorprod(alpha_1,alpha_2,2);
Sp_3 = tensorprod(Qe_pos,a_12,[1,2],[1,2]);
Sm_3 = tensorprod(Qe_neg,1./a_12,[1,2],[1,2]);

% address the positive & negative elements, using location indicator
alpha_neg_loc_3 = (Sm_3~=0);
alpha_pos_loc_3 = ~alpha_neg_loc_3;

% update \mu_j in two cases, i.e.,
alpha_var_3 = zeros(num_node,1);
% 1) when a_{ij} < 0, compute \mu_j using function 'FB_solver'
alpha_var_3(alpha_neg_loc_3) = Gsolver(Sp_3(alpha_neg_loc_3),Sm_3(alpha_neg_loc_3),c_3(alpha_neg_loc_3));
% 2) when a_{ij} > 0, compute \mu_j using sinkhorn iteration
alpha_var_3(alpha_pos_loc_3) = c_3(alpha_pos_loc_3)./Sp_3(alpha_pos_loc_3);
% 3) update \nu from the two cases
alpha_3 = alpha_var_3;


%% 4. convergence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update \P according to the element-wise closed-form expression. (to be changed)
% P = diag(exp(-mu))*A_pos*diag(exp(-alpha_1)) + diag(exp(mu))*A_neg*diag(exp(alpha_1));

a12  = tensorprod(alpha_1,alpha_2,2);
a12(isnan(a12)) = 0;
a_123 = tensorprod(a12,alpha_3);
a_123(isnan(a_123)) = 0;
P  = Qe_pos.*a_123 + Qe_neg.*(1./a_123);

% compute and save the value of the objective function.
obj = P.*log(P./bQ);
obj(isnan(obj)) = 0;
obj_value = sum(real(obj),'all');
cost(end+1) = obj_value;


%% 5.Error storage
% compuate and save the violation/error of the two marginals (marginal constraints).
sign_P = X.*P;
Err_1(end+1) = norm(squeeze(sum(sum(permute(sign_P,[2,3,1])))) - c_1)/norm(c_1);
Err_2(end+1) = norm(squeeze(sum(sum(permute(sign_P,[3,1,2])))) - c_2)/norm(c_2);
Err_3(end+1) = norm(squeeze(sum(sum(sign_P)))-c_3)/norm(c_3);


end

squeeze(sum(sum(permute(sign_P,[2,3,1])))) - c_1
squeeze(sum(sum(permute(sign_P,[3,1,2])))) - c_2
squeeze(sum(sum(sign_P)))-c_3


%% Marginal convergence plot

% figure(1)
% subplot(3,1,1);
% plot(log10(Err_1),'LineWidth',2,'Color',"#0072BD"); 
% axis tight; 
% % title('log|\sum P - c_1|');
% 
% subplot(3,1,2);
% plot(log10(Err_2),'LineWidth',2,'Color',"#D95319"); 
% axis tight; 
% % title('log|\sum P  - c_2|');
% 
% subplot(3,1,3);
% plot(log10(Err_3),'LineWidth',2,'Color',"#77AC30"); 
% axis tight; 
% % title('log|\sum P - c_3|');


%%

figure(2)
pos_point = sign_P;
pos_point(pos_point<0) = 0;

neg_point = sign_P;
neg_point(neg_point>0) = 0;


list1 = find (pos_point>0);
[I1,J1,K1] = ind2sub(size(pos_point),list1);

list2 = find (neg_point<0);
[I2,J2,K2] = ind2sub(size(neg_point),list2);

hold on

for i=1:n
     marginal1 = scatter3(i,0.5,0.5,300.*c_1(i),'filled','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75]);
end
% 
for i=1:n
     marginal2 = scatter3(0.5,i,0.5,300.*c_2(i),'filled','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75]);
end
% 
for i=1:n
     marginal3 = scatter3(0.5,0.5,i,300.*c_3(i),'filled','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75]);
end

C = ones(length(I1),1)*[0 0 1];
scatter3sph(I1,J1,K1,'size',(1/3).*sqrt(pos_point(list1)),'color',C,'trans',0.6);

C = [1 0 0];
scatter3sph(I2,J2,K2,'size',(1/3).*sqrt(abs(neg_point(list2))),'color',C,'trans',0.6);


hold off
axis tight
axis equal;
grid on

xticks([1:3])
yticks([1:3])
zticks([1:3])
view([66 16])


%%
figure(3)

Q = Q./sum(Q,"all");
q_1 = squeeze(sum(sum(permute(Q,[2,3,1]))));
q_2 = squeeze(sum(sum(permute(Q,[3,1,2]))));
q_3 = squeeze(sum(sum(Q)));



Qpos_point = Q;
Qpos_point(Qpos_point<0) = 0;

Qneg_point = Q;
Qneg_point(Qneg_point>0) = 0;


list1 = find (Qpos_point>0);
[I1,J1,K1] = ind2sub(size(Qpos_point),list1);

list2 = find (Qneg_point<0);
[I2,J2,K2] = ind2sub(size(Qneg_point),list2);

hold on

for i=1:n
     marginal1 = scatter3(i,0.5,0.5,300.*q_1(i),'filled','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75]);
end
% 
for i=1:n
     marginal2 = scatter3(0.5,i,0.5,300.*q_2(i),'filled','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75]);
end
% 
for i=1:n
     marginal3 = scatter3(0.5,0.5,i,300.*q_3(i),'filled','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75]);
end


C = ones(length(I1),1)*[0 0 1];
scatter3sph(I1,J1,K1,'size',(1/3).*sqrt(Qpos_point(list1)),'color',C,'trans',0.6);

C = [1 0 0];
scatter3sph(I2,J2,K2,'size',(1/3).*sqrt(abs(Qneg_point(list2))),'color',C,'trans',0.6);


hold off
axis tight
grid on
xticks([1:3])
yticks([1:3])
zticks([1:3])
view([66 16])
















