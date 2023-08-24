clear all 
clc
clf

% Read data
addpath('function')
gene_data = readcell('RPPA.csv');
%[num1,num2] = size(gene_data);
num1 = 150;
num2 = 40;

ccl = gene_data(2:num1,1);  % cancer cell line
mrna = gene_data(1,2:num2); % mrna list

data =  gene_data(2:num1,2:num2);
index = find(strcmp(data,'NA' ));
data(index) = {0};
data_mat = cell2mat(data);

sign_mat = sign(abs(data_mat));
abs_sign = abs(sign_mat);

sign_pos = (abs_sign + sign_mat)./2;
sign_neg = (abs_sign - sign_mat)./2;

prior_pos = sign_pos.*data_mat;
prior_neg = -sign_neg.*data_mat;


%% cancer_cell and mrna distributions generator

nc = num1-1;
nm = num2-1;

cancer_axis = (0:nc-1)'/nc;
mrna_axis = (0:nm-1)'/nm;

Gaussian = @(x,t0,sigma)exp( -(x-t0).^2/(2*sigma^2) );
normalize = @(p)p/sum(p(:));
sigma1 = .3;
sigma2 = .6;

cancer_dist = Gaussian(cancer_axis, .2, sigma1) + Gaussian(cancer_axis, .8, sigma2);
mrna_dist = Gaussian(mrna_axis, .5, sigma2);

c_1 = normalize(cancer_dist);
c_2 = normalize(mrna_dist);


%%

[n1,n2] = size(data_mat);
alpha1 = ones(n1,1);
alpha2 = ones(n2,1);

T = 200;
cost = [];

bQ_pos = prior_pos.*exp(-1);
bQ_neg = prior_neg.*exp(-1);

Err_1 = [];
Err_2 = [];

%% Iteration

for t = 1:T

%%%%%%%%%%%%%%%%%%
% update of alpha1
a_1 = bQ_pos*alpha2;
b_1 = bQ_neg*(1./alpha2);

% location indicator
alpha1_neg_loc = (b_1~=0);
alpha1_pos_loc = ~alpha1_neg_loc;
% two cases 
alpha1_var = zeros(n1,1);
alpha1_var(alpha1_neg_loc) = Gsolver(a_1(alpha1_neg_loc),b_1(alpha1_neg_loc),c_1(alpha1_neg_loc));
alpha1_var(alpha1_pos_loc) = c_1(alpha1_pos_loc)./a_1(alpha1_pos_loc);
alpha1 = alpha1_var;

bP = diag(alpha1)*bQ_pos*diag(alpha2) + diag(1./alpha1)*bQ_neg*diag(1./alpha2);
P = sign_mat.*bP;
Err_1(end+1) = norm(sum(P,1)'-c_2)/norm(c_2);


%%
%%%%%%%%%%%%%%%%%%
% update of alpha2
a_2 = bQ_pos'*alpha1;
b_2 = bQ_neg'*(1./alpha1);


% location indicator
alpha2_neg_loc = (b_2~=0);
alpha2_pos_loc = ~alpha2_neg_loc;
% two cases 
alpha2_var = zeros(n2,1);
alpha2_var(alpha2_neg_loc) = Gsolver(a_2(alpha2_neg_loc),b_2(alpha2_neg_loc),c_2(alpha2_neg_loc));
alpha2_var(alpha2_pos_loc) = c_2(alpha2_pos_loc)./a_2(alpha2_pos_loc);
alpha2 = alpha2_var;


%%
%%%%%%%%%%%%%%%%%%
% convergence 
bP = diag(alpha1)*bQ_pos*diag(alpha2) + diag(1./alpha1)*bQ_neg*diag(1./alpha2);

obj = bP.*log(bP);
obj(isnan(obj)) = 0;
obj_value = sum(real(obj),'all');
cost(end+1) = obj_value;

P = sign_mat.*bP;
Err_2(end+1) = norm(sum(P,2)-c_1)/norm(c_1);

end

%% Result checking -- prelimaries

disp('objective value:')
sum(obj_value,'all')
sum(sum(P,1)'-c_2,'all')
sum(sum(P,2)-c_1,'all')
% disp(sum(P,1)'-c_2);
% disp(sum(P,2)-c_1);


%% Plot a graph 



cancer_graph = (P>0);
g = graph(1:n1+n2,1:n1+n2,1,'omitselfloops');


[sp,tp] = find(P>0);
[sn,tn] = find(P<0);


colormap1 = linspecer(n1,'sequential');
colormap1 = lines(n1);

colormap2 = lines(n2);

gc = subgraph(g,1:n1);
gm = subgraph(g,n1+1:n1+n2);

angc = linspace(0,2*pi-0.06,n1);
rc = 1;
cx = rc * cos(angc);
cy = rc * sin(angc);

angm = linspace(0,2*pi-0.06,n2);
rm = 0.4;
mx = rm * cos(angm);
my = rm * sin(angm);

% Sample

sample_number = 1000;
rand_list1 = randsample(length(sp),sample_number);
rand_list1 = sort(rand_list1);

rand_list2 = randsample(length(sn),sample_number);
rand_list2 = sort(rand_list2);

% rand_list1 = 1:length(sp);
% 
% rand_list2 = 1:length(sn);


figure(1)

hold on 
for ii = 1:length(rand_list1)
i = rand_list1(ii);
coord1 = [cx(sp(i)),mx(tp(i))];
coord2 = [cy(sp(i)),my(tp(i))];
% check1 = cx(sp(i)).^2 + cy(sp(i)).^2
% check2 = mx(tp(i)).^2 + my(tp(i)).^2
plot(coord1,coord2,'Color',[0, 0, 1, 0.16],'LineWidth',1)
end

for jj = 1:length(rand_list2)
j = rand_list2(jj);
coord1 = [cx(sn(j)),mx(tn(j))];
coord2 = [cy(sn(j)),my(tn(j))];
plot(coord1,coord2,'Color',[1, 0, 0, 0.16],'LineWidth',1,'LineStyle','--')
end

gcplot = plot(gc,'XData',cx,'YData',cy);
gcplot.NodeColor = colormap1;
gcplot.MarkerSize = 8.*cancer_dist;
gcplot.NodeFontSize = 2;
gcplot.Marker = 'o';
gcplot.LineWidth = 0.8;
hc = text(cx.*1.05, cy.*1.05, ccl, 'FontSize',8);
set(hc,{'Rotation'},num2cell(angc'*180/pi))


gmplot = plot(gm,'XData',mx,'YData',my);
gmplot.NodeColor = colormap2;
gmplot.MarkerSize = 8.*mrna_dist;
gmplot.NodeFontSize = 2;
gmplot.Marker = 'o';
gmplot.LineWidth = 0.8;
hm = text(mx.*1.05, my.*1.05, mrna, 'FontSize',9);
set(hm,{'Rotation'},num2cell(angm'*180/pi))


set(gca,'XColor', 'none','YColor','none')
% set(gca, 'color', 'none');
hold off