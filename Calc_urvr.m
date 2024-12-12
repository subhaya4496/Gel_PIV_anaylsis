%% Edge length Calculation
clear all
dirIn = 'D:\Active_gel_Project\Data\Gel_23_08_2023\Strain_u_r_u_phi_stuck';
fileExtension = '.mat';
directoryContents = dir([dirIn, filesep, ['*' fileExtension]]);
filenames={};
[filenames{1:length(directoryContents),1}] = deal(directoryContents.name);
    
filenames = sort_nat(filenames);
amount = length(filenames);
Rad = zeros(amount,1001);
ur_v_r = [];
uphi_v_r = [];
urstd = [];
upstd = [];
R_all = [];
trough = [];
Frame = [];
MinIdx = [];
for file = 1:amount
    % number = file;
    filen = cell2mat(filenames(file));
    % loadfile = sprintf('Strain_%g_u_r.mat',number);
    loadfile = fullfile(dirIn, filen);
    if exist(loadfile,'file')
    load(loadfile)
    numbin = 25;
    phi = atan2(Yp,Xp);
    Xp(isnan(u_r)) = NaN;
    Yp(isnan(u_r)) = NaN;
    R = sqrt(Xp.^2 + Yp.^2);
    Phi = phi;
    Phi(abs(Phi)>pi/3) = NaN;
    % Phi(Phi>0) = NaN;
    A = 6.4*reshape(R,1,[]);
    B = 6.4*reshape(u_r,1,[]);
    C = 6.4*reshape(u_phi,1,[]);
    D = reshape(Phi,1,[]);
    A = A(~isnan(D));
    B = B(~isnan(D));
    C = C(~isnan(D));

    [A,I] = sort(A);
    B = B(I);
    C = C(I);
    R_bound = R*6.4;
    R_bound(isnan(u_r)) = NaN;
    R_disk = R_bound;
    for j = 2:width(R)-1
        for i = 2:height(R)-1
            if (~isnan(R_disk(i+1,j)))&&(~isnan(R_disk(i-1,j)))&&(~isnan(R_disk(i,j-1)))&&(~isnan(R_disk(i,j+1)))
                R_bound(i,j) = NaN;
            end
        end
    end
    Rmax = ceil(max(R_bound,[],"all")/numbin)*numbin;
    % if 6.4*sqrt(X_CoC^2 + Y_CoC^2)<Rmax/12
        binsize = Rmax/numbin;
        rmean = zeros(1,numbin);
        rstd = zeros(1,numbin);
        u_r_mean = zeros(1,numbin);
        u_phi_mean = zeros(1,numbin);
        u_r_std = zeros(1,numbin);
        u_phi_std = zeros(1,numbin);
        i = 1;
        % check = sum(R,"all",'omitnan');
        % if check>0
            for span_r = 1: numbin  
                A_r = A(A<span_r*binsize);
                A_sr = A(i:length(A_r));
                B_sr = B(i:length(A_r));
                rmean(span_r) = mean(A_sr(~isnan(B_sr)));
                rstd(span_r) = std(A_sr(~isnan(B_sr)));
                u_r_mean(span_r) = mean(B(i:length(A_r)),'omitnan');
                u_r_std(span_r) = std(B(i:length(A_r)),'omitnan');
                u_phi_mean(span_r) = mean(C(i:length(A_r)),'omitnan');
                u_phi_std(span_r) = std(C(i:length(A_r)),'omitnan');
                i = length(A_r) + 1;
            end
            ur_v_r = [ur_v_r, rmean', u_r_mean'];
            urstd = [urstd, u_r_std'];
            uphi_v_r = [uphi_v_r, rmean', u_phi_mean'];
            upstd = [upstd, u_phi_std'];
            if isempty(A)
                R_all(file) = NaN;
            else
                R_all(file)=max(A);
            end
            rad_file =  0:R_all(file)/1000:R_all(file);
            Rad(file,:) = rad_file;
            rmean = rmean(~isnan(rmean));
            u_r_mean = u_r_mean(~isnan(u_r_mean));
            u_phi_mean = u_phi_mean(~isnan(u_phi_mean));
            if length(u_r_mean)>1
                urs = csaps(rmean, u_r_mean, 0.000001, rad_file);
                ups = csaps(rmean, u_phi_mean, 0.000001, rad_file);
                u_r_spline(file,:) = urs;
                u_p_spline(file,:) = ups;
                DataInv = 1.01*max(urs) - urs;
                [Minima,MinIdx] = findpeaks(DataInv);
            end
            if isempty(MinIdx)
                trough(file) = NaN;
            else
                trough(file) = rad_file(urs == urs(max(MinIdx)));
            end
        % end
        fileno = filen;
        fileno(1:7) = [];
        fileno(end-3:end)=[];
        fileno = str2double(fileno);
        Frame = [Frame,fileno];
        MinIdx = [];
    end
end
% % R_all(end) = [];
% % trough(end) = [];
delta =  R_all - trough;
D_over_R = delta./R_all;
% vR_all = gradient(R_all,Frames);
%%

ur_v_r_ave40_R = [];
ur_v_r_ave40_uR = [];
no_of_I_ave = 40;
for k = 1:ceil(amount/no_of_I_ave)
    ur_v_r_ave40_R(:,k) = mean(ur_v_r(:,2*(k-1)*no_of_I_ave+1:2:2*(k)*no_of_I_ave),2,'omitnan');
    ur_v_r_ave40_uR(:,k) = mean(ur_v_r(:,2*(k-1)*no_of_I_ave+2:2:2*(k)*no_of_I_ave),2,'omitnan');
end

%%
% dirOut = 'D:\Active_gel_Project\Data\Gel_23_08_2023\Strain_u_phi_v_r_3';
dirOut1 = 'D:\Active_gel_Project\Data\Gel_23_08_2023\Strain_u_r_v_r_s';
% if exist(dirOut, 'dir')
%     rmdir(dirOut, 's');
% end
% mkdir(dirOut);
if exist(dirOut1, 'dir')
    rmdir(dirOut1, 's');
end
mkdir(dirOut1)
f = figure('visible','off');
f.Position = [100 50 640 470];
for Fra = 1:length(Frame)
    clf
    if sum(~isnan(uphi_v_r(:,2*Fra-1)))>0
        plot(ur_v_r(:,2*Fra-1)/R_all(Fra),ur_v_r(:,2*Fra), 'Marker', 'o', 'Linestyle', 'none', 'Color',[0 0.4470 0.7410], 'Linewidth', 2, 'MarkerSize', 9)
        hold on
        % plot(uphi_v_r(:,2*Fra-1)/R_all(Frame(Fra)),sign(min(u_p_spline(Frame(Fra),:)))*uphi_v_r(:,2*Fra)/min(u_p_spline(Frame(Fra),:)), 'Marker', 'diamond', 'Linestyle', 'none', 'Color',[0.8500 0.3250 0.0980], 'Linewidth', 2, 'MarkerSize', 9)
%         plot(Rad(Fra,:)/R_all(Fra), u_r_spline(Fra,:),'Color',[0 0.4470 0.7410], 'Linewidth', 2)
        tt = sprintf('Frame = %g', Frame(Fra)+199);
        title(tt, 'FontSize',40)
        xlabel('R/R_{disk}')
        ylabel('u_r ({\mu}m)')
        % ymin = ceil(max(abs(ur_v_r(:,2*Fra)))/0.1)*0.1+0.1;
        % ymax = ceil(max(ur_v_r(:,2*Fra))/0.1)*0.1+0.1;
        xlim([0 1])
        % ylim([-1.5 0.5])
        fontsize(f,20,"points")
        fontname(f, 'Cambria Math')
%         name = sprintf('u_phi_vs_R_%g.fig', Fra);
        nametiff = sprintf('u_phi_vs_R_%g.tif', Fra);
%         figfile = fullfile(dirOut, name);
        figfile_tiff = fullfile(dirOut1, nametiff);
%         saveas(f, figfile)
        saveas(f, figfile_tiff)
    end
end
%% uR vs R averaged over frames

dirOut = 'D:\Active_gel_Project\Data\Gel_24_06_2024\Strain_u_r_v_r_10';
dirOut1 = 'D:\Active_gel_Project\Data\Gel_24_06_2024\Strain_u_r_v_r_10_tiff';
if exist(dirOut, 'dir')
    rmdir(dirOut, 's');
end
mkdir(dirOut);
if exist(dirOut1, 'dir')
    rmdir(dirOut1, 's');
end
mkdir(dirOut1)
f = figure('visible','off');
f.Position = [100 50 640 470];
for Fra = 1:width(ur_v_r_ave10_uR)
    clf
%     if sum(~isnan(ur_v_r(:,2*Fra-1)))>0
        plot(ur_v_r_ave10_R(:,Fra)/max(ur_v_r_ave10_R(:,Fra)),ur_v_r_ave10_uR(:,Fra), 'Marker', 'o', 'Linestyle', 'none', 'Color',[0 0.4470 0.7410], 'Linewidth', 2)
        hold on
    %     plot(Rad(Fra,:)/R_all(Fra), u_p_spline(Fra,:),'Color',[0 0.4470 0.7410], 'Linewidth', 2, 'MarkerSize', 9 )
        tt = sprintf('Frame = %g - %g', 10*(Fra-1)+1000, 10*(Fra)+1000);
        title(tt, 'FontSize',40)
        xlabel('R/R_{disk}')
        ylabel('u_{R}')
        ymin = floor(min(ur_v_r_ave10_uR(:,Fra))/0.1)*0.1-0.05;
        ymax = ceil(max(ur_v_r_ave10_uR(:,Fra))/0.1)*0.1+.05;
        xlim([0 1])
        ylim([ymin ymax])
        fontsize(f,20,"points")
        fontname(f, 'Cambria Math')
        name = sprintf('u_R_vs_R_%g.fig', Fra);
        nametiff = sprintf('u_R_vs_R_%g.tif', Fra);
        figfile = fullfile(dirOut, name);
        figfile_tiff = fullfile(dirOut1, nametiff);
        saveas(f, figfile)
        saveas(f, figfile_tiff)
%     end
end

%%
    figure
    for Fra = [1, 51, 101, 151, 201, 251, 301]   %102,152,
        tt = sprintf('Frame = %g', Fra+499);
        plot(ur_v_r(:,2*Fra-1),-ur_v_r(:,2*Fra)*10, 'Marker','o', 'Linewidth', 2,'DisplayName',tt)
        hold on
    end