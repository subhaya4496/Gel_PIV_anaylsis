%% Finding Eigenvectors
clear all
dirIn = 'D:\Active_gel_Project\Data\Gel_23_08_2023\Strain_uf_new_3';
dirOut = 'D:\Active_gel_Project\Data\Gel_23_08_2023\Strain_eigen_max_min_3';
fileExtension = '.mat';
directoryContents = dir([dirIn, filesep, ['*' fileExtension]]);
filenames={};
[filenames{1:length(directoryContents),1}] = deal(directoryContents.name);

if exist(dirOut, 'dir')
    rmdir(dirOut, 's');
end

mkdir(dirOut);

amount = length(filenames);
for file_ord = 1:amount
    loadfile = sprintf('Strain_%g.mat',file_ord);
    lf = fullfile(dirIn, loadfile);
    load(lf, 'x', 'y', 'u_xx', 'u_yy', 'u_xy')
    EigenvectorxT = NaN(height(u_xx), width(u_xx));
    EigenvectoryT = NaN(height(u_xx), width(u_xx));
    Eigenvectorxt = NaN(height(u_xx), width(u_xx));
    Eigenvectoryt = NaN(height(u_xx), width(u_xx));
    EigenvalueT = NaN(height(u_xx), width(u_xx));
    Eigenvaluet = NaN(height(u_xx), width(u_xx));
    for i = 1:size(u_xx,1)
        for j = 1:size(u_xx,2)
            A = [u_xx(i,j) u_xy(i,j); u_xy(i,j) u_yy(i,j)];
%             A = [1/2*(u_xx(i,j)-u_yy(i,j)) u_xy(i,j); u_xy(i,j) 1/2*(u_yy(i,j)-u_xx(i,j))];
            [V,D] = eig(A);
            Eigvec = reshape(V, [4,1]);
            Eigval = reshape(D, [4,1]);
            [T,I] = max(abs(Eigval), [], "all");
            Tmax= Eigval(I);
            dd = Eigval;
            dd = abs(dd);
            dd(dd==0) = inf;
            [t,J] = min(dd, [], "all");
            Tmin = Eigval(J);
            if I==1
                EigenvectorxT(i, j) = Eigvec(I);
                EigenvectoryT(i, j) = Eigvec(I+1);
                Eigenvectorxt(i, j) = Eigvec(I+2);
                Eigenvectoryt(i, j) = Eigvec(I+3);
            elseif I==4
                EigenvectorxT(i, j) = Eigvec(I-1);
                EigenvectoryT(i, j) = Eigvec(I);
                Eigenvectorxt(i, j) = Eigvec(I-3);
                Eigenvectoryt(i, j) = Eigvec(I-2);
            end
            EigenvalueT(i,j) = Tmax;
            Eigenvaluet(i,j) = Tmin;
        end
    end
    save_file_name = sprintf('Strain_%g.mat', file_ord);
    matfile = fullfile(dirOut, save_file_name);
    save(matfile, 'u_xx', 'u_yy', 'u_xy', 'x', 'y', 'EigenvectorxT','EigenvectoryT','Eigenvectorxt','Eigenvectoryt', 'Eigenvaluet','EigenvalueT')
end

%% (Not used - Just to check if the calculations are okay)

f = figure(1);
EigenvalueT(EigenvalueT==0) = NaN;    
contourf(X,Y,EigenvalueT, 50, 'EdgeColor','none')
%     caxis([-0.3 0.3])
    
    hold on
    
    quiver(X,Y,EigenvectorxT,EigenvectoryT,'r','ShowArrowHead','off')
    hold off
    set(gca, 'Color', 'k')
    colorbar
    ylim([-270 270])
    xlim([-270 270])
    axis square
    
%% Mean Eigenvalues Max and Min vs R

clear all
dirIn = 'D:\Active_gel_Project\Data\300frames_strain\Strain_eigen_max_min_3';
fileExtension = '.mat';
directoryContents = dir([dirIn, filesep, ['*' fileExtension]]);
filenames={};
[filenames{1:length(directoryContents),1}] = deal(directoryContents.name);
amount = length(filenames);
Eigen_abs_max = [];
Eigen_abs_min = [];
% filenames = cell2mat(filenames);
for file = 1:amount
    loadfile = sprintf('Strain_%g.mat',file);
    lf = fullfile(dirIn, loadfile);
    load(lf, 'x', 'y', 'EigenvalueT', 'Eigenvaluet')
    R = sqrt(x.^2 + y.^2);
    cos_the = x./R;
    sin_the = y./R;
    numbin = 27;
    Eigenvaluet(Eigenvaluet==0) = NaN;
    EigenvalueT(EigenvalueT==0) = NaN;
    A = 6.4*reshape(R,1,[]);
    B = reshape(EigenvalueT,1,[]);
    C = reshape(Eigenvaluet,1,[]);
    [A,I] = sort(A);
    B = B(I);
    C = C(I);
    Rmax = ceil(max(A(~isnan(B)))/numbin)*numbin;
    binsize = Rmax/numbin;
    rmean = zeros(1,numbin);
    rstd = zeros(1,numbin);
    eigen_max = zeros(1,numbin);
    eigen_min = zeros(1,numbin);
    eigen_max_std = zeros(1,numbin);
    eigen_min_std = zeros(1,numbin);
    i = 1;
    for span_r = 1: numbin  
        A_r = A(A<span_r*binsize);
        A_sr = A(i:length(A_r));
        B_sr = B(i:length(A_r));
        rmean(span_r) = mean(A_sr(~isnan(B_sr)));
        rstd(span_r) = std(A_sr(~isnan(B_sr)));
        eigen_max(span_r) = mean(B(i:length(A_r)),'omitnan');
        eigen_max_std(span_r) = std(B(i:length(A_r)),'omitnan');
        eigen_min(span_r) = mean(C(i:length(A_r)),'omitnan');
        eigen_min_std(span_r) = std(C(i:length(A_r)),'omitnan');
        i = length(A_r) + 1;
    end
    Eigen_abs_max = [Eigen_abs_max, rmean', eigen_max'];
    Eigen_abs_min = [Eigen_abs_min, rmean', eigen_min'];
    R_all(file) = max(A(~isnan(B)));
    rad_file =  0:R_all(file)/1000:R_all(file);
    Rad(file,:) = rad_file;
%     emaxs = csaps(rmean, eigen_max, 0.0000005, rad_file);
%     emins = csaps(rmean, eigen_min, 0.0000005, rad_file);
%     eig_max_spline(file,:) = emaxs;
%     eig_min_spline(file,:) = emins;
%                    
end
%%
dirOut1 = 'D:\Active_gel_Project\Data\300frames_strain\Strain_eigenvalue_tiff_3';
if exist(dirOut1, 'dir')
    rmdir(dirOut1, 's');
end
mkdir(dirOut1)
f = figure('Visible','off');
f.Position = [100 50 750 570];
for Fra = 1:amount
    clf
    plot(Eigen_abs_max(:,2*Fra-1)/R_all(Fra),Eigen_abs_max(:,2*Fra), 'Marker', 'o', 'Linestyle', 'none', 'Color',[0 0.4470 0.7410], 'Linewidth', 2,'DisplayName','|Eigenvalue_{max}|')
    hold on
%     plot(Rad(Fra,:)/R_all(Fra), eig_max_spline(Fra,:),'Color',[0 0.4470 0.7410], 'Linewidth', 2, 'MarkerSize', 9 ,'DisplayName','Max Eigen Spline')
    plot(Eigen_abs_min(:,2*Fra-1)/R_all(Fra),Eigen_abs_min(:,2*Fra), 'Marker', 'diamond', 'Linestyle', 'none', 'Color',[0.8500 0.3250 0.0980], 'Linewidth', 2,'DisplayName','|Eigenvalue_{min}|')
    hold on
%     plot(Rad(Fra,:)/R_all(Fra), eig_min_spline(Fra,:),'Color',[0.8500 0.3250 0.0980], 'Linewidth', 2, 'MarkerSize', 9,'DisplayName','Min Eigen Spline' )
    tt = sprintf('Frame = %g', Fra+199);
    title(tt, 'FontSize',40)
    xlabel('R/R_{disk}')
    ylabel('Strain Eigenvalue')
    ymin = floor(-abs(min(Eigen_abs_max(:,2*Fra))/0.005))*0.005;
    ymax = ceil(max(abs(Eigen_abs_max(:,2*Fra))/0.005))*0.005;
    xlim([0 1])
    ylim([-ymax ymax])
    fontsize(f,20,"points")
    fontname(f, 'Cambria Math')
    legend
    nametiff = sprintf('Eigval_vs_R_%g.tiff', Fra);
    figfile_tiff = fullfile(dirOut1, nametiff);
    saveas(f, figfile_tiff)
end
%% Orientational order from eigen vectors 
clear all
dirIn = 'D:\Active_gel_Project\Data\300frames_strain\Strain_eigen_max_min_3';
fileExtension = '.mat';
directoryContents = dir([dirIn, filesep, ['*' fileExtension]]);
filenames={};
[filenames{1:length(directoryContents),1}] = deal(directoryContents.name);
amount = length(filenames);
dirOut1 = 'D:\Active_gel_Project\Data\300frames_strain\Strain_eigenvector_order';
if exist(dirOut1, 'dir')
    rmdir(dirOut1, 's');
end
mkdir(dirOut1)

Mean_order = [];
Error_order = [];
% filenames = cell2mat(filenames);
for file = 1:amount
    loadfile = sprintf('Strain_%g.mat',file);
    lf = fullfile(dirIn, loadfile);
    load(lf, 'x', 'y', 'EigenvalueT', 'Eigenvaluet','EigenvectorxT', 'EigenvectoryT', 'Eigenvectorxt', 'Eigenvectoryt')
    R = sqrt(x.^2 + y.^2);
    cos_phi = x./R;
    sin_phi = y./R;
    numbin = 25;
%     Eigenvaluet(Eigenvaluet==0) = NaN;
%     EigenvalueT(EigenvalueT==0) = NaN;
    A = 6.4*reshape(R,1,[]);
    B = reshape(EigenvectorxT,1,[]);
    C = reshape(EigenvectoryT,1,[]);
    Bt = reshape(Eigenvectorxt,1,[]);
    Ct = reshape(Eigenvectoryt,1,[]);
    D = reshape(cos_phi,1,[]);
    E = reshape(sin_phi,1,[]);
    ET = reshape(EigenvalueT,1,[]);
    Et = reshape(Eigenvaluet,1,[]);
    [A,I] = sort(A);
    B = B(I);
    C = C(I);
    Bt = Bt(I);
    Ct = Ct(I);
    D = D(I);
    E = E(I);
    ET = ET(I);
    Et = Et(I);
    cos_sq_psi = (C.*D - B.*E).^2;
    Order = (abs(ET).*cos_sq_psi + abs(Et).*(1-cos_sq_psi))./(abs(ET)+abs(Et));
    Rmax = ceil(max(A(~isnan(B)))/numbin)*numbin;
    binsize = Rmax/numbin;
    rmean = zeros(1,numbin);
    rstd = zeros(1,numbin);
    order_mean = zeros(1,numbin);
    order_std = zeros(1,numbin);
    i = 1;
    for span_r = 1: numbin  
        A_r = A(A<span_r*binsize);
        A_sr = A(i:length(A_r));
        B_sr = B(i:length(A_r));
        rmean(span_r) = mean(A_sr(~isnan(B_sr)));
        rstd(span_r) = std(A_sr(~isnan(B_sr)));
        order_mean(span_r) = mean(Order(i:length(A_r)),'omitnan');
        order_std(span_r) = std(Order(i:length(A_r)),'omitnan');
        i = length(A_r) + 1;
    end
    Mean_order = [Mean_order, rmean', order_mean'];
    Error_order = [Error_order, order_std'];
    
    R_all(file) = max(A(~isnan(B)));
    rad_file =  0:R_all(file)/1000:R_all(file);
    Rad(file,:) = rad_file;
end
%%
f = figure('Visible','off');
f.Position = [100 50 750 570];
for Fra = 1:amount
    clf
    plot(Mean_order(:,2*Fra-1)/R_all(Fra),Mean_order(:,2*Fra), 'Marker', 'o', 'Linestyle', 'none', 'Color',[0 0.4470 0.7410], 'Linewidth', 2)       % Error_order, ,'DisplayName','|Eigenvalue_{max}|'
%     hold on
    tt = sprintf('Frame = %g', Fra+199);
    title(tt, 'FontSize',40)
    xlabel('R/R_{disk}')
    ylabel('Orientation Order')
%     ymin = floor(-abs(min(Eigen_abs_max(:,2*Fra))/0.005))*0.005;
%     ymax = ceil(max(abs(Eigen_abs_max(:,2*Fra))/0.005))*0.005;
    xlim([0 1])
    ylim([0 1])
    fontsize(f,20,"points")
    fontname(f, 'Cambria Math')
%     legend
    nametiff = sprintf('Ori_Order_vs_R_%g.tiff', Fra);
    figfile_tiff = fullfile(dirOut1, nametiff);
    saveas(f, figfile_tiff)
end
