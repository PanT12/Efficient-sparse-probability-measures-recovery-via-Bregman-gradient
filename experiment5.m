%% demo of the experiment 5 (hyperspectral unmixing)
% The following codes are provided in https://github.com/ricardoborsoi/MUA_SparseUnmixing/tree/master
% A Fast Multiscale Spatial Regularization for Sparse Hyperspectral Unmixing
%     R.A. Borsoi, T. Imbiriba, J.C.M. Bermudez, C. Richard.
%     IEEE Geoscience and Remote Sensing Letters, 2018.
%% HERE YOU NEED TO LOAD THE CUPRITE IMAGE FILE. YOU CAN DOWNLOAD IT FROM http://www.lx.it.pt/~bioucas/code.htm
addpath("data/Experiment5/")
load("data/Experiment5/cuprite_ref.mat")

Yimage = reshape(x',Lines,Columns,L);

Y = reshape(Yimage, [size(Yimage,1)*size(Yimage,2) L])';


%% build the library
load("data/Experiment5/USGS_1995_Library.mat")
%  order bands by increasing wavelength
[dummy index] = sort(datalib(:,1));
A =  datalib(index,4:end);
names = names(4:end,:);

% order  the columns of A by decreasing angles 
[A, index, angles] = sort_library_by_angle(A);
names = names(index',:);
namesStr = char(names);

A = A(BANDS,:);
%% training
con = sqrt(mean(Y.^2,"all"));

% our proposed algorithm
X_hat = zeros(size(A,2),(Lines*Columns));
lambda = 5;
parfor i = 1:(Lines*Columns)
    out = sparse_l0(A, lambda, Y(:, i), 1, 1e-6, false);
    X_hat(:,i) = out.x;
end
% save("result/Experiment5/X_hat_l0.mat",'X_hat'); 

% SUNSAL algorithm
lambda = 1e-3;

tic
[X_hat_l1] =  sunsal(A,Y,'lambda',lambda,'ADDONE','no','POSITIVITY','yes', ...
                    'TOL',1e-4, 'AL_iters',2000,'verbose','yes');
timeSunsal = toc;
% save("result/Experiment5/X_hat_l1.mat",'X_hat_l1'); 
%% showing result

clear all;clc
load("data/Experiment5/cuprite_ref.mat")
load("result/Experiment5/X_hat_GPG.mat")
load("result/Experiment5/X_hat_l0.mat")
load("result/Experiment5/X_hat_l1.mat")

n = size(X_hat_GPG,1);

X_hat_l1 = X_hat_l1./sum(X_hat_l1);
X_hat_l0_im = reshape(X_hat_l0', Lines,Columns,n); 
X_hat_l1_im = reshape(X_hat_l1', Lines,Columns,n); 
X_hat_GPG_im = reshape(X_hat_GPG', Lines,Columns,n); 

ind = zeros((Lines*Columns),3);
for i = 1:(Lines*Columns)
    ind(i,1) = sum(X_hat_l1(:,i) > 0);
    ind(i,2) = sum(X_hat_l0(:,i) > 0);
    ind(i,3) = sum(X_hat_GPG(:,i) > 0);
end

min(ind)
max(ind)
mean(ind)

% material_idx = [420 336 297];
% material_idx = [285 368 465];
material_idx = [420 465 297];

% range_mat = [0 0.38; ... % Alunite
%              0 0.30; ... % Buddingtonite
%              0 0.65];    % Chalcedony
range_mat = [0 0.38; ... % Alunite
             0 0.50; ... % Kaolinite
             0 0.65];    % Chalcedony
% name_mat = {'Alunite','Buddingtonite','Chalcedony'};
name_mat = {'Alunite','Kaolinite','Chalcedony'};

count = 0;
ii = 1;
figure
[ha, pos] = tight_subplot(3,3,[.025 .015],[.05 .05],[.05 .15]);

for idx=material_idx
    count = count + 1;
    y = name_mat(count);

    axes(ha(ii)); ii = ii+1;
    imagesc(X_hat_l1_im(:,:,idx), range_mat(count,:))
    axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
    if count == 3
        xlabel("SUnSAL", "fontsize", 16)
    end
    ylabel(y, "fontsize", 16)

    axes(ha(ii)); ii = ii+1;
    imagesc(X_hat_l0_im(:,:,idx), range_mat(count,:))
    axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
    if count == 3
        xlabel("Alg.2", "fontsize", 16)
    end

    axes(ha(ii)); ii = ii+1;
    imagesc(X_hat_GPG_im(:,:,idx), range_mat(count,:))
    axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
    if count == 3
        xlabel("GPG", "fontsize", 16)
    end
    originalSize2 = get(gca, 'Position');

    h=colorbar; 
    set(ha(ii-1), 'Position', originalSize2);
    % set(h,'fontsize',5);
    set(h,'fontsize',10);

    colormap jet
    % print(strcat('examples/realImg/estim_abundances_',name_mat{count},'_tght'),'-dpdf')

end
