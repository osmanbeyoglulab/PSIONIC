clear;
clc;
close;

addpath(genpath('./L1General'));
load('GYN_ccle.mat'); %%sample dataset

task=size(Y,2);
y = cell(task, 1); %gene expression
x = cell(task, 1); %motif scores

for t = 1:task
    %Log?10-transformed RNA-seq TPM gene expression values were unit-normalized 
    y{t} = normalize_column(Y(:,t)); 
    x{t} = zscore(log10(D+1)); 
end

O = 'R';
rng('default'); 
k = 7;
Maxsteps=100;
lambda=0.001;
mu = 0.001;
O = 'R';
[L,S,W] = GO_MTL(x,y,k,Maxsteps,lambda,mu,O);
What = L*S;

final_performance = zeros(task,1);
for t= 1:task
    final_performance(t) = diag(corr(y{t}, x{t}*What(:,t) , 'type', 'Spearman'));
end
mean(final_performance)

colormap(gray);
imagesc(abs(S));
colorbar;
