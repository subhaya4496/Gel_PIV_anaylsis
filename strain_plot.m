clear all
dirIn = 'D:\Active_gel_Project\Data\Gel_28_07_2023\Strain_u_r_u_phi_stuck_2';
% dirOut = 'D:\Active_gel_Project\Data\Gel_12_01_2024\Strain_uxx+uyy_fig';
dirOut2 = 'D:\Active_gel_Project\Data\Gel_28_07_2023\Strain_urr+upp_tiff';
fileExtension = '.mat';
directoryContents = dir([dirIn, filesep, ['*' fileExtension]]);

% if exist(dirOut, 'dir')
%     rmdir(dirOut, 's');
% end
if exist(dirOut2, 'dir')
    rmdir(dirOut2, 's');
end

% mkdir(dirOut)
mkdir(dirOut2)
cd(dirIn)
filenames={};
[filenames{1:length(directoryContents),1}] = deal(directoryContents.name);
    
% filenames = sort_nat(filenames);

amount = length(filenames);
f = figure('visible','off');
f.Position = [100 50 640 470];
for file_ord = 1:amount
    clf
    loadfile = sprintf('Strain_%g.mat',file_ord);
%     loadfile = cell2mat(loadfile);
    lf = fullfile(dirIn, loadfile);
    if exist(lf,'file')
        load(lf)
    else
        continue
    end
    U = u_rr+u_phiphi;
    Umin = min(U,[],'all');
    Umax = max(U,[],'all');
    U(U>0)=U(U>0)/Umax;
    U(U<0)=-U(U<0)/Umin;
    t = U(U~=0);
    if ~isempty(t)
       U(U==0) = NaN;
        contourf(x,-y, U, 50, 'EdgeColor','none')
        hold on
        % quiver(x,-y,+EigenvectorxT,-EigenvectoryT,0.5,'r', 'ShowArrowHead','off')
        % hold on
        % quiver(x,-y,-EigenvectorxT,+EigenvectoryT,0.5,'r', 'ShowArrowHead','off')
        quiver(x,-y,uxnew,-uynew,'r', 'ShowArrowHead','on')
        % plot(0,0,'o','MarkerSize', 7,'MarkerEdgeColor','r','LineWidth',2)
        plot(X_CoC2,-Y_CoC2,'o','MarkerSize', 7,'MarkerEdgeColor','r','LineWidth',2)
        hold off
        set(gca, 'Color', 'k')
        colorbar
        clim([-1 1])
        ylim([-150 150])
        xlim([-150 150])
        tt = sprintf('Frame = %g', file_ord-1);
        title(tt, 'FontSize',32)
        xlabel('x (px)')
        ylabel('y (px)')
        fontsize(f,18,"points")
        fontname(f, 'Cambria Math')
        axis square
        fig_file = loadfile;
        fig_file(end-3:end) = [];
        tiff_file = fig_file;
        fig_file(end+1:end+4) = '.fig';
        tiff_file(end+1:end+5) = '.tiff';
    %     fig_file = cell2mat(fig_file);
%         save_fig_file = fullfile(dirOut, fig_file);
        save_tiff_file = fullfile(dirOut2, tiff_file);
%         saveas(f, save_fig_file)
        saveas(f, save_tiff_file)
    end
end

%%
a = decon_Stream0293_tif_mask;
a = double(a);
x = [];
y = [];
a(a > 0) = -1;
a(a==0) = 1;
for i = 2:width(decon_Stream0293_tif_mask)-1
    for j = 2:height(decon_Stream0293_tif_mask)-1
        if (a(i,j)*(a(i,j+1)+a(i,j-1)+a(i+1,j)+a(i-1,j)) < 4)
            x = [x, i];
            y = [y, j];
        end
    end
end
XY = [x;y];
a(a<0) = NaN;


