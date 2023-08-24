clear all 
clc
clf 


addpath('function')

format short

% p = [0.2;0.3;0.1;0.4];
% a = [0,-1,1,1; 1,0,1,1; 1,1,0,1; 1,1,1,0];

% p = [0.3;0.3;0.4];
% a = [1,1,1;
%      -1,0,1;
%      1,0,0];

c_1 = [0.1;0.05;0.05;0.15;0.2;0.05;0.03;0.07;0.25;0.05];
c_2 = [0.05;0.1;0.05;0.2;0.07;0.15;0.05;0.25;0.03;0.05];

% Q = [1  0 0 0 -1 0 1 0  1 0;
%      0  0 1 0  0 0 0 0  0 1;
%      0  1 0 0  1 0 1 0  0 0;
%      0  0 0 0  1 0 0 1  0 0;
%     -1  0 1 1  1 -1 0 0 0 0;
%      0  0 0 0  -1 0 1 0 1 0;
%      1  0 1 0  0 1 1 0  0 0;
%      0  0 0 1  0 0 0 0  1 -1;
%      1  0 0 0  0 1 0 1  0 1;
%      0  1 0 0  0 0 0 -1 1 0];

Q = [1  -1 1 0  1 0 -1 0  1 1;
    -1   0 1 0  0 0  1 0  0 1;
     1   1 0 0  1 0  1 1  0 0;
     0   0 0 0  1 0  0 1  0 0;
     1   0 1 1  1 1  0 0 0 0;
     0   0 0 0  1 0  1 0 1 0;
     -1  1 1 0  0 1  1 0  0 0;
     0   0 1 1  0 0  0 0  1 1;
     1   0 0 0  0 1  0 1  0 1;
     1   1 0 0  0 0  0 1 1 0];

X = sign(Q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bQ = abs(Q);
Q_pos = (bQ + Q)./2;
Q_neg = (bQ - Q)./2;


% Initialization
[n,~] = size(Q);
alpha1 = ones(n,1);
alpha2 = ones(n,1);

T = 200;
cost = [];

bQ_pos = Q_pos.*exp(-1);
bQ_neg = Q_neg.*exp(-1);

Err_1 = [];
Err_2 = [];

%% Iteration

for t = 1:T

%%%%%%%%%%%%%%%%%%
% update of alpha1
a_1 = bQ_pos'*alpha2;
b_1 = bQ_neg'*(1./alpha2);

% location indicator
alpha1_neg_loc = (b_1~=0);
alpha1_pos_loc = ~alpha1_neg_loc;
% two cases 
alpha1_var = zeros(n,1);
alpha1_var(alpha1_neg_loc) = Gsolver(a_1(alpha1_neg_loc),b_1(alpha1_neg_loc),c_1(alpha1_neg_loc));
alpha1_var(alpha1_pos_loc) = c_1(alpha1_pos_loc)./a_1(alpha1_pos_loc);
alpha1 = alpha1_var;

% alpha1 = Gsolver(a_1,b_1,c_1);
%%
%%%%%%%%%%%%%%%%%%
% update of alpha2
a_2 = bQ_pos*alpha1;
b_2 = bQ_neg*(1./alpha1);

% location indicator
alpha2_neg_loc = (b_2~=0);
alpha2_pos_loc = ~alpha2_neg_loc;
% two cases 
alpha2_var = zeros(n,1);
alpha2_var(alpha2_neg_loc) = Gsolver(a_2(alpha2_neg_loc),b_2(alpha2_neg_loc),c_2(alpha2_neg_loc));
alpha2_var(alpha2_pos_loc) = c_2(alpha2_pos_loc)./a_2(alpha2_pos_loc);
alpha2 = alpha2_var;

%%%alpha2 = Gsolver(a_2,b_2,c_2);

%%

%%%%%%%%%%%%%%%%%%
% convergence 
bP = diag(alpha1)*bQ_pos*diag(alpha2) + diag(1./alpha1)*bQ_neg*diag(1./alpha2);

obj = bP.*log(bP);
obj(isnan(obj)) = 0;
obj_value = sum(real(obj),'all');
cost(end+1) = obj_value;

P = X.*bP;

% Error storage
Err_1(end+1) = norm(sum(P,1)'-c_2)/norm(c_2);
Err_2(end+1) = norm(sum(P,2)-c_1)/norm(c_1);

end

%% Result checking -- prelimaries
disp('objective value:')
sum(obj_value,'all')

disp(sum(P,1)'-c_2)
disp(sum(P,2)-c_1)

%%
figure(1);
plot(1:T,cost);
axis tight;
xlabel('Iteration') 
ylabel('Objective') 
title('Convergence')

figure(2);
subplot(2,1,1)
plot(1:T,log10(Err_1));axis tight; title('log|| P^T 1 - c_2||');
subplot(2,1,2)
plot(1:T,log(Err_2));axis tight; title('log|| P 1 - c_1||');


%% 
figure(3)

pos_point = P;
pos_point(pos_point<0) = 0;

neg_point = P;
neg_point(neg_point>0) = 0;


list1 = find (pos_point>0);
[I1,J1] = ind2sub(size(pos_point),list1);

list2 = find (neg_point<0);
[I2,J2] = ind2sub(size(neg_point),list2);

c_1 = [0.1;0.05;0.05;0.15;0.2;0.05;0.03;0.07;0.25;0.05];
c_2 = [0.05;0.1;0.05;0.2;0.07;0.15;0.05;0.25;0.03;0.05];


hold on

for i=1:n
     marginal1 = plot(i,-1,'o','MarkerSize',100.*c_1(i),'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75]);
end

for i=1:n
     marginal2 = plot(-1,i,'o','MarkerSize',100.*c_2(i),'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75]);
end


for i=1:length(I1)
     pt_pos = plot(I1(i),J1(i),'o','MarkerSize',...
     100*pos_point(I1(i),J1(i)),'MarkerFaceColor','b');
end

for i=1:length(I2)
     pt_neg = plot(I2(i),J2(i),'o','MarkerSize',...
     100*abs(neg_point(I2(i),J2(i))),'MarkerFaceColor','r');
end

hold off

grid on
xticks([1:10])
yticks([1:10])
xlim([-2,10])
ylim([-2,10])
set(gca,"FontSize",20)
% xlim([-1 10])
% ylim([-1 10])
% axis off


%%
figure(4)

Q = Q./sum(Q,'all');

Qpos_point = Q;

Qpos_point(Qpos_point<0) = 0;

Qneg_point = Q;
Qneg_point(Qneg_point>0) = 0;


list1 = find (Qpos_point>0);
[I1,J1] = ind2sub(size(Qpos_point),list1);

list2 = find (Qneg_point<0);
[I2,J2] = ind2sub(size(Qneg_point),list2);

q_1 = sum(Q,2);
q_2 = sum(Q,1)';


hold on

for i=1:n
    marginal1 = plot(i,-1,'o','MarkerSize',100.*q_1(i),'LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75]);
end

for i=1:n
     marginal2 = plot(-1,i,'o','MarkerSize',100.*q_2(i),'LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75]);
end


for i=1:length(I1)
     pt_pos = plot(I1(i),J1(i),'o','MarkerSize',...
     400*Qpos_point(I1(i),J1(i)),'MarkerFaceColor','b');
end

for i=1:length(I2)
     pt_neg = plot(I2(i),J2(i),'o','MarkerSize',...
     400*abs(Qneg_point(I2(i),J2(i))),'MarkerFaceColor','r');
end

hold off

grid on
xticks([1:10])
yticks([1:10])
xlim([-2,10])
ylim([-2,10])
set(gca,"FontSize",20)






