clear all
clc
addpath('function')


%% Read Table
%opts.SelectedVariableNames = ["Gene","Tissue","TPM","pTPM","nTPM"];
T = readtable("rna_gtex.tsv","FileType","text",'Delimiter', '\t');

num_node = 50; 
num = num_node*37;

% num = 407;
% num_node = 11;
% the first 832 genes for the sake of simplicity
T = T(1:num,:); % 30748


%% Tissue Type
tissue_list = T(:,"Tissue");
tissue_name = unique(tissue_list);

gene_list = T(:,"Gene");
gene_name = unique(gene_list);

gene = T(:,"GeneName");
name = unique(gene);


%% tissue cell

for i = 1:37
tissue{i} = T(i:37:num,:);
end

%%
names = tissue{1}.GeneName;

figure(1)
% p = plot(tissue{1}.nTPM);
p = bar(tissue{1}.nTPM);
% p.LineStyle = ":";
% p.Marker = ".";
% p.LineWidth = 3;
set(gca,'xtick',[1:num_node],'xticklabel',names,'Fontsize',14)
ylabel('Transcripts Per Million')

normalize = @(p)p/sum(p(:));


%% Correlation Matrix
% R = corrcoef(A)
% C = cov(A)
% C = corr(A)
num_tissue = 30;

dist1 = normalize(tissue{1}.nTPM);
dist2 = normalize(tissue{1}.pTPM);
dist3 = normalize(tissue{1}.TPM);

for i = 1:num_tissue
Obv_TPM(:,i) = normalize(tissue{i}.TPM);
end

for i = 1:num_tissue
Obv_nTPM(:,i) = normalize(tissue{i}.nTPM);
end

for i = 1:num_tissue
Obv_pTPM(:,i) = normalize(tissue{i}.pTPM);
end

cov = corr(Obv_pTPM',Obv_nTPM','Type','Pearson');
% cov = corr(Obv_pTPM',Obv_nTPM','Type','Kendall');

% figure()
imagesc(cov)  
colormap jet

%%
cov = cov-diag(diag(cov));


%% Adjacency matrix
treshold = 0.2;  % 0.2-50 
cov( abs(cov)<= treshold) = 0;
% cov = cov-diag(diag(cov));

% figure()
imagesc(cov)


signA = sign(cov);
adj = abs(signA);




%% Network Visualization (symmetric adjancenct matrix A is required)
up_A = triu(signA);

g = graph(adj,name,'upper');

[neg_row,neg_col] = find(up_A < 0 );
redEdge = [neg_row, neg_col];
colormap = lines(num_node);
% colormap = linspecer(num_node,'sequential');
figure(2)
gplot = plot(g,'Layout','circle');
gplot.NodeColor = colormap;
gplot.MarkerSize = 10;
gplot.NodeFontSize = 14;
gplot.Marker = 'o';
gplot.LineWidth = 0.1;
highlight(gplot,neg_row,neg_col,'EdgeColor','r','LineWidth',2.5)
labelnode(gplot,[1:num_node],names)


%% Sinkhorn-Initial
test_num = 37;
p = normalize(tissue{test_num}.pTPM);

Adj_pos = (signA>0);
Adj_neg = (signA<0);

Q = abs(cov);

%% main function
iter= 100;

[P,Pi,cost,Err_1,Err_2] = fb(Adj_pos,Adj_neg,Q,p,iter);


%% Result & Constraints check

disp('..............Forward-Backward.............')

disp('..............objective value')
disp(cost(end))

disp('.....Marginal constraint and totol sum')
disp(Pi'*p-p)
sum( Pi,2 )
sum(Pi,'all')
% min(Pi)'
% max(Pi)'

figure(3);
plot(1:iter,cost);
axis tight;
xlabel('Iteration') 
ylabel('Objective') 
title('Convergence')

figure(4);
subplot(2,1,1)
plot(1:iter,log10(Err_1));axis tight; title('log||\Pi^T p - p||');
subplot(2,1,2)
plot(1:iter,log(Err_2));axis tight; title('log||\Pi 1 - 1||');


%% Error plot

 error1 = Pi'*normalize(tissue{test_num}.nTPM) - normalize(tissue{test_num}.pTPM);
% error2 = [Pi'*normalize(tissue{test_num}.nTPM),normalize(tissue{test_num}.nTPM),error1];
% error3 = [normalize(tissue{test_num}.nTPM),error1];
error4 = [Pi'*normalize(tissue{test_num}.nTPM),normalize(tissue{test_num}.nTPM)];
% disp(error2)
figure(5)

%p = bar(error3,'stacked');
b = bar(error4);
b(1).BarWidth = 1.5;
b(2).BarWidth = 1.5;
set(gca,'xtick',[1:num_node],'xticklabel',names,'Fontsize',12)
ylabel('Normalized TPM')

% figure(6)
% imagesc(cov)








