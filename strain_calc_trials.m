clear all
dirIn = 'D:\Active_gel_Project\Data\Gel_28_07_2023\Gel';
dirOut = 'D:\Active_gel_Project\Data\Gel_28_07_2023\Strain_uf_stuck';
dirMask = 'D:\Active_gel_Project\Data\Gel_28_07_2023\Mask_file';
fileExtension = '.mat';
directoryContents = dir([dirIn, filesep, ['*' fileExtension]]);
filenames={};
[filenames{1:length(directoryContents),1}] = deal(directoryContents.name);
fileExtension = '.tif';
directoryContents = dir([dirMask, filesep, ['*' fileExtension]]);
masknames={};
[masknames{1:length(directoryContents),1}] = deal(directoryContents.name);

if exist(dirOut, 'dir')
    rmdir(dirOut, 's');
end

mkdir(dirOut);

amount = length(filenames); 
XCM_new = 0;
YCM_new = 0;
for file_ord = 1:amount
    file1 = cell2mat(filenames(file_ord));
    mask1 = cell2mat(masknames(file_ord+103));
    loadfile = fullfile(dirIn, file1);
    loadmask = fullfile(dirMask, mask1);
    mask = imread(loadmask);
    if exist(loadfile,"file")
        load(loadfile, 'x', 'y', 'u_filt', 'v_filt')
    else
        continue
    end
    xx = x(1,:);
    yy = y(:,1);
    small_mask = mask(:,xx);
    small_mask = small_mask(yy,:);
    dx = xx(2) - xx(1);
    dy = yy(2) - yy(1);
    X_prime = x;
    Y_prime = y;
    X_prime(small_mask>0) = NaN;
    Y_prime(small_mask>0) = NaN;
    u_filt(small_mask>0) = NaN;
    v_filt(small_mask>0) = NaN;
    XCM_new = mean(X_prime,'all','omitnan');
    x = x - XCM_new;
    YCM_new = mean(Y_prime,'all','omitnan');
    y = y - YCM_new;
    R = sqrt(x.^2 + y.^2);
    cos_phi = x./R;
    sin_phi = y./R;
    % % add more lines
    u_x = u_filt;
    u_y = v_filt;
    u_xx = NaN(height(x),width(x));
    u_yy = u_xx;
    u_xy = u_xx;
    for i = 2:height(u_x)-1
        for j = 2:width(u_y)-1
            u_xx(i,j) = (u_x(i,j+1)-u_x(i,j-1))/(2*dx);
            u_yy(i,j) = (u_y(i+1,j)-u_y(i-1,j))/(2*dy);
            u_xy(i,j) = 1/2*((u_x(i+1,j)-u_x(i-1,j))/(2*dy)+(u_y(i,j+1)-u_y(i,j-1))/(2*dx));
        end
    end
    % UU = abs(u_xx+u_yy);
    % UU = heaviside(-(u_xx+u_yy)).*u_yy;
    % X_CoC = sum(UU.*x,'all', 'omitnan')/sum(UU, 'all', 'omitnan');
    % Y_CoC = sum(UU.*y,'all', 'omitnan')/sum(UU, 'all', 'omitnan');
    Mean_u_x = mean(u_x,'all','omitnan');
    Mean_u_y = mean(u_y,'all','omitnan');
    a = (x.*u_y - y.*u_x)./(Mean_u_x*u_y - Mean_u_y*u_x);
    ax = x(~isnan(a));
    dia_ax = max(ax,[],'all') - min(ax,[],'all');
    X_int = a*Mean_u_x;
    Y_int = a*Mean_u_y;
    check = sum(u_xx,'all','omitnan');
    if check~=0
        X_int(X_int>abs(dia_ax/2)) =NaN;
        Y_int(X_int>abs(dia_ax/2)) =NaN;
        X_CoC = median(X_int,'all','omitnan');
        Y_CoC = median(Y_int,'all','omitnan');
        uxnew = u_x - Mean_u_x;
        uynew = u_y - Mean_u_y;
        aprime = (x.*uynew - y.*uxnew)./(Mean_u_x*uynew - Mean_u_y*uxnew);
        X_int2 = aprime*Mean_u_x;
        Y_int2 = aprime*Mean_u_y;
        X_int2(X_int2>abs(dia_ax/2)) =NaN;
        Y_int2(X_int>abs(dia_ax/2)) =NaN;
        X_CoC2 = median(X_int2,'all','omitnan');
        Y_CoC2 = median(Y_int2,'all','omitnan');
    
    u_rr = u_xx.*cos_phi.^2 + u_yy.*sin_phi.^2 + 2*u_xy.*sin_phi.*cos_phi;
    u_phiphi = u_xx.*sin_phi.^2 + u_yy.*cos_phi.^2 - 2*u_xy.*sin_phi.*cos_phi;
    u_rphi = -(u_xx - u_yy).*sin_phi.*cos_phi + u_xy.*(cos_phi.^2 + sin_phi.^2);
    save_file_name = sprintf('Strain_%g.mat', file_ord);
    matfile = fullfile(dirOut, save_file_name);
    save(matfile, 'u_xx', 'u_yy', 'u_xy', 'u_rphi', 'u_phiphi', 'u_rr', 'R', 'u_filt','v_filt','x','y', 'u_x', 'u_y','XCM_new','YCM_new','uxnew','uynew','X_int','Y_int','a','X_CoC','Y_CoC','X_int2','Y_int2','aprime','X_CoC2','Y_CoC2')
    end
end

%%
f = figure(2);
contourf(x,-y, u_xx+u_yy, 50, 'EdgeColor','none')
hold on
quiver(x,-y,u_x,-u_y,'r', 'ShowArrowHead','on')
hold off
set(gca, 'Color', 'k')
colorbar
ylim([-270 270])
xlim([-270 270])
axis square
        
%%
dirIn = 'D:\Active_gel_Project\Data\Gel_23_08_2023\Processed_raw';
dirIn2 = 'D:\Active_gel_Project\Data\Gel_23_08_2023\Strain_uf_stuck';
dirOut = 'D:\Active_gel_Project\Data\Gel_23_08_2023\Processed_good';
fileExtension = '.mat';
directoryContents = dir([dirIn2, filesep, ['*' fileExtension]]);
filenames2={};
[filenames2{1:length(directoryContents),1}] = deal(directoryContents.name);
% filenames2 = cell2mat(filenames2);
if exist(dirOut, 'dir')
    rmdir(dirOut, 's');
end
mkdir(dirOut);
fileExtension = '.tif';
directoryContents = dir([dirIn, filesep, ['*' fileExtension]]);
filenames={};
[filenames{1:length(directoryContents),1}] = deal(directoryContents.name);
filenames = cell2mat(filenames);
f=figure('Visible','off');
f.Position = [100 50 640 670];
for file_ord = 1:length(filenames)
    loadfile = fullfile(dirIn, filenames(file_ord,:));
    matfile = sprintf('Strain_%g.mat',file_ord);
    loadfile2 = fullfile(dirIn2,matfile);
    if exist(loadfile2,'file')
        load(loadfile2)
        clf
        imshow(loadfile)
        hold on
        quiver(x+XCM_new,y+YCM_new,uxnew,uynew,'Color','g')
        plot(XCM_new,YCM_new,'o','MarkerSize', 7,'MarkerEdgeColor','r','LineWidth',2)
        % plot(XCM_new+X_CoC,YCM_new+Y_CoC,'o','MarkerSize', 7,'MarkerEdgeColor','w','LineWidth',2)
        hold off
        tt = sprintf('Frame = %g', file_ord+199);
        title(tt, 'FontSize',15)
        savefile = filenames(file_ord,:);
        save_tiff_file = fullfile(dirOut, savefile);
        imagewd = getframe(gcf); 
        imwrite(imagewd.cdata, save_tiff_file);
    end
end
