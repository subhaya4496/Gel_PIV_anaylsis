%% u_RR and u_phiphi vs R

clear all
dirIn = 'D:\Active_gel_Project\Data\300frames_strain\Strain_u_r_u_phi_stuck_3';
fileExtension = '.mat';
directoryContents = dir([dirIn, filesep, ['*' fileExtension]]);
filenames={};
[filenames{1:length(directoryContents),1}] = deal(directoryContents.name);
    
filenames = sort_nat(filenames);
amount = length(filenames);
Rad = zeros(amount,1001);
urr_v_r = [];
upp_v_r = [];
trace_strain = [];
Frame = [];
for file = 1:amount
    number = file;
    % loadfile = sprintf('Strain_%g_u_r.mat',number);
    % loadfile = fullfile(dirIn, loadfile);
    filen = cell2mat(filenames(file));
    % loadfile = sprintf('Strain_%g_u_r.mat',number);
    loadfile = fullfile(dirIn, filen);
    if exist(loadfile,'file')
    load(loadfile)
    numbin = 25;
    phi = atan2(Yp,Xp);
    Phi = phi;
    Phi(abs(Phi)>pi/3) = NaN;
    A = 6.4*reshape(R,1,[]);
    B = reshape(u_rr,1,[]);
    C = reshape(u_phiphi,1,[]);
    D = reshape(u_xx+u_yy,1,[]);
    E = reshape(Phi,1,[]);
    A = A(~isnan(E));
    B = B(~isnan(E));
    C = C(~isnan(E));
    D = D(~isnan(E));
    [A,I] = sort(A);
    B = B(I);
    C = C(I);
    D = D(I);
    R_bound = R*6.4;
    R_bound(isnan(u_rr)) = NaN;
    R_disk = R_bound;
    for j = 2:width(R)-1
        for i = 2:height(R)-1
            if (~isnan(R_disk(i+1,j)))&&(~isnan(R_disk(i-1,j)))&&(~isnan(R_disk(i,j-1)))&&(~isnan(R_disk(i,j+1)))
                R_bound(i,j) = NaN;
            end
        end
    end
    Rmax = ceil(max(R_bound,[],"all")/numbin)*numbin;
    binsize = Rmax/numbin;
    rmean = zeros(1,numbin);
    phimean = zeros(1,numbin);
    rstd = zeros(1,numbin);
    u_rr_mean = zeros(1,numbin);
    u_phiphi_mean = zeros(1,numbin);
    STRAINT = zeros(1,numbin);
    u_r_std = zeros(1,numbin);
    u_phi_std = zeros(1,numbin);
    i = 1;
    % check = sum(R,"all",'omitnan');
        for span_r = 1: numbin  
            A_r = A(A<span_r*binsize);
            A_sr = A(i:length(A_r));
            B_sr = B(i:length(A_r));
%             D_sr = D(i:length(A_r));
            rmean(span_r) = mean(A_sr(~isnan(B_sr)));
            rstd(span_r) = std(A_sr(~isnan(B_sr)));
            u_rr_mean(span_r) = mean(B(i:length(A_r)),'omitnan');
            u_r_std(span_r) = std(B(i:length(A_r)),'omitnan');
            u_phiphi_mean(span_r) = mean(C(i:length(A_r)),'omitnan');
            u_phi_std(span_r) = std(C(i:length(A_r)),'omitnan');
            STRAINT(span_r) = mean(D(i:length(A_r)));
            i = length(A_r) + 1;
        end
%         u_rr = gradient(u_rr_mean,rmean);
%         u_pp = gradient(u_phiphi_mean,phimean)./rmean + u_rr_mean./rmean;
        urr_v_r = [urr_v_r, rmean', u_rr_mean'];
        upp_v_r = [upp_v_r, rmean', u_phiphi_mean'];
        trace_strain = [trace_strain, rmean', STRAINT'];
        R_all(file) = max(A);
        rad_file =  0:R_all(file)/1000:R_all(file);
        Rad(file,:) = rad_file;
        CapiR(file) = max(rmean);
    end  
    Frame = [Frame,file];
end

%% Plotting data
Frames = (1:amount);
dirOut1 = 'D:\Active_gel_Project\Data\300frames_strain\Strain_urr_v_r_tiff';
if exist(dirOut1, 'dir')
    rmdir(dirOut1, 's');
end
mkdir(dirOut1)
f = figure('visible','off');
f.Position = [100 50 640 470];
for Fra = 1:amount
    clf
    if sum(~isnan(urr_v_r(:,2*Fra-1)))>0
        plot(urr_v_r(:,2*Fra-1)/R_all(Fra),urr_v_r(:,2*Fra), 'Marker', '*', 'MarkerSize', 10, 'Linestyle', 'none', 'Color',[0 0.4470 0.7410], 'Linewidth', 2,'DisplayName','\epsilon_{rr}+\epsilon_{\phi\phi}')
        hold on
        plot(upp_v_r(:,2*Fra-1)/R_all(Fra),upp_v_r(:,2*Fra), 'Marker', 'diamond', 'MarkerSize', 10,'Linestyle', 'none', 'Color',[0.8500 0.3250 0.0980], 'Linewidth', 2,'DisplayName','\epsilon_{rr}-\epsilon_{\phi\phi}')
%         plot(trace_strain(:,2*Fra-1)/R_all(Fra),trace_strain(:,2*Fra), 'Marker', 'hexagram', 'MarkerSize', 10,'Linestyle', 'none', 'Color',[0.8500 0.3250 0.0980], 'Linewidth', 2,'DisplayName','Strain trace')
        tt = sprintf('Frame = %g', Fra+100);
        title(tt, 'FontSize',40)
        xlabel('R/R_{disk}')
        ylabel('Strain')
        ymin = floor(min(urr_v_r(:,2*Fra)+upp_v_r(:,2*Fra))/0.001)*0.001-0.001;
        % ymax = ceil(max(uphi_v_r(:,2*Fra))/0.1)*0.1+0.001;
        xlim([0 1])
        ylim([ymin 0.008])
        fontsize(f,20,"points")
        fontname(f, 'Cambria Math')
        legend('Location','southeast')
        nametiff = sprintf('Strain_vs_R_%g.tif', Fra);
        figfile_tiff = fullfile(dirOut1, nametiff);
        saveas(f, figfile_tiff)
    end
end

%% Time averaged strain over a fixed number of Frames (getting data)
No_frames_mean = 10;
u_rr_fra_mean = [];
u_pp_fra_mean = [];
r_fra_mean = [];
for fra = 1:floor(length(Frame)/No_frames_mean)
    u_rr_fra_mean = [u_rr_fra_mean, mean(urr_v_r(:,2*(fra-1)*No_frames_mean+2:2:2*(fra)*No_frames_mean),2,'omitnan')];
    u_pp_fra_mean = [u_pp_fra_mean, mean(upp_v_r(:,2*(fra-1)*No_frames_mean+2:2:2*(fra)*No_frames_mean),2,'omitnan')];
    r_fra_mean = [r_fra_mean, mean(urr_v_r(:,2*(fra-1)*No_frames_mean+1:2:2*(fra)*No_frames_mean),2,'omitnan')];
end

%% Time averaged strain over a fixed number of Frames (plotting data)
Frames = (1:amount);
dirOut1 = 'D:\Active_gel_Project\Data\300frames_strain\Strain_urr_v_r_tiff';
if exist(dirOut1, 'dir')
    rmdir(dirOut1, 's');
end
mkdir(dirOut1)
f = figure('visible','off');
f.Position = [100 50 640 470];
for Fra = 1:amount
    clf
    if sum(~isnan(urr_v_r(:,2*Fra-1)))>0
        plot(urr_v_r(:,2*Fra-1)/R_all(Fra),urr_v_r(:,2*Fra), 'Marker', '*', 'MarkerSize', 10, 'Linestyle', 'none', 'Color',[0 0.4470 0.7410], 'Linewidth', 2,'DisplayName','\epsilon_{rr}+\epsilon_{\phi\phi}')
        hold on
        plot(upp_v_r(:,2*Fra-1)/R_all(Fra),upp_v_r(:,2*Fra), 'Marker', 'diamond', 'MarkerSize', 10,'Linestyle', 'none', 'Color',[0.8500 0.3250 0.0980], 'Linewidth', 2,'DisplayName','\epsilon_{rr}-\epsilon_{\phi\phi}')
%         plot(trace_strain(:,2*Fra-1)/R_all(Fra),trace_strain(:,2*Fra), 'Marker', 'hexagram', 'MarkerSize', 10,'Linestyle', 'none', 'Color',[0.8500 0.3250 0.0980], 'Linewidth', 2,'DisplayName','Strain trace')
        tt = sprintf('Frame = %g', Fra+100);
        title(tt, 'FontSize',40)
        xlabel('R/R_{disk}')
        ylabel('Strain')
        ymin = floor(min(urr_v_r(:,2*Fra)+upp_v_r(:,2*Fra))/0.001)*0.001-0.001;
        % ymax = ceil(max(uphi_v_r(:,2*Fra))/0.1)*0.1+0.001;
        xlim([0 1])
        ylim([ymin 0.008])
        fontsize(f,20,"points")
        fontname(f, 'Cambria Math')
        legend('Location','southeast')
        nametiff = sprintf('Strain_vs_R_%g.tif', Fra);
        figfile_tiff = fullfile(dirOut1, nametiff);
        saveas(f, figfile_tiff)
    end
end

%% 
clear all
dirIn = 'D:\Active_gel_Project\Data\300frames_strain\Strain_u_r_u_phi_new_3\';
fileExtension = '.mat';
directoryContents = dir([dirIn, filesep, ['*' fileExtension]]);
filenames={};
[filenames{1:length(directoryContents),1}] = deal(directoryContents.name);
% dirOut = 'D:\Active_gel_Project\Data\Gel_19_10_2023\Strain_u_phiphi_v_phi\';
UPP_all = [];
P_all = [];
for fileord = 1:length(filenames)
    loadfile = (fullfile(dirIn, filenames(fileord,:)));
    loadfile = cell2mat(loadfile);
    load(loadfile)
    phi = atan2(Yp,Xp);
    Phi = phi;
    Phi(abs(Phi)>pi/3) = NaN;
    Phi(isnan(u_r)) = NaN;
    R = sqrt(Xp.^2 + Yp.^2);
    A = 6.4*reshape(R,1,[]);
    B = 6.4*reshape(u_r,1,[]);
    C = reshape(u_phi,1,[]);
    D = reshape(Phi,1,[]);
    D = mod(D,2*pi);
    [A,I] = sort(A);
    B = B(I);
    C = C(I);
    D = D(I);
    A = A(~isnan(D));
    B = B(~isnan(D));
    C = C(~isnan(D));
    D = D(~isnan(D));
    RR = 0;
    d_bound = max(A) - A;
    U_R = 0;
    D_bound = 0;
    U_PHI = 0;
    PHI = 0;
    d = max(A)/6;
    delta_R = d;
    j = 1;
    for i = 1:length(A)
        if d_bound(i) > d - delta_R && d_bound(i) < d +delta_R
            RR(j) = A(i);
            D_bound(j) = d_bound(i);
            U_R(j) = B(i);
            U_PHI(j) = C(i);
            PHI(j) = D(i);
            j = j+1;
        end
    end
    [PHI,K] = sort(PHI);
    U_PHI = U_PHI(K);
    U_R = U_R(K);
    D_bound = D_bound(K);
    RR = RR(K);
    
    phiii = zeros(1,40);
    U_rad = zeros(1,40);
    U_phi = zeros(1,40);
    j=0;
    for i = 1:length(phiii)
        subphi = PHI(PHI<2*pi/(length(phiii))*i);
        k = length(subphi);
        subphi = subphi(subphi>2*pi/(length(phiii))*(i-1));
        phiii(i) = mean(subphi, 'omitnan');
    %     phiii = phiii + pi/4;
        U_rad(i) = mean(U_R(j+1:k), 'omitnan');
        U_phi(i) = mean(U_PHI(j+1:k), 'omitnan');
        j=k;
    end
    phiii = mod(phiii+pi/3,2*pi);
    [phiii,sort_p] = sort(phiii);
    U_rad = U_rad(sort_p);
    U_phi = U_phi(sort_p);
    Nana = U_phi(~isnan(U_phi));
    if length(Nana)>=2
        upp = gradient(U_phi,phiii)+ U_rad/(max(A)-d);
%         urr = gradient(U_rad,)
        UPP_all = [UPP_all;upp];
        P_all = [P_all;phiii];
    end
end

%%
dirOut2 = 'D:\Active_gel_Project\Data\Gel_19_10_2023\Strain_u_phiphi_v_phi_time_ave_25_tiff';
if exist(dirOut2, 'dir')
    rmdir(dirOut2, 's');
end
mkdir(dirOut2)
averageframe = 25;
P = length(phiii);
for frame = 1:floor(height(UPP_all)/averageframe)
    f = figure(1);
    clf
    f.Position = [100 50 950 470];      
    hold on
    dista = sprintf('R = 11R_{disk}/12');
    P = mean(P_all((frame-1)*averageframe+1:frame*averageframe,:),'omitnan');
    UPP = mean(UPP_all((frame-1)*averageframe+1:frame*averageframe,:),'omitnan');
    plot(P, UPP, 'LineStyle', ':', 'Marker', '*', 'MarkerSize', 9, 'LineWidth',2.5,'DisplayName',dista)   %'Color','[0 0.4470 0.7410]'
    xlim([0 2*pi/3+0.001])
    % ymax = max(abs(upp));
    ylim([-1 0.5])
    title_name = sprintf('Frame = %g to %g', (frame-1)*averageframe+1+499,frame*averageframe+1+499);
    title(title_name)
    ylabel({'u_{\phi\phi}'});
    xlabel({'\phi'});
    set(gca,'FontName','Cambria Math','FontSize',25,'XTick',...
        [0 pi/3 2*pi/3],'XTickLabel',...
        {'-\pi/3', '0','\pi/3'});
    
    saveimg = sprintf('U_phiphi_frame=%g.fig',frame);
    savetiff = sprintf('U_phiphi_frame=%g.tif',frame);
    % savefile = (fullfile(dirOut, saveimg));
    savetotiff = (fullfile(dirOut2, savetiff));
    saveas(f,savetotiff)
    hold off
end

%% Mean u_phiphi vs R (extra - Not used)
clear all
dirIn = 'D:\Active_gel_Project\Data\300frames_strain\Strain_u_r_u_phi_new_3';
% dirOut = 'D:\Active_gel_Project\Data\Gel_23_08_2023\Strain_u_phiphi_v_r';
dirOut1 = 'D:\Active_gel_Project\Data\300frames_strain\Strain_v_r_tiff_3';
fileExtension = '.mat';
directoryContents = dir([dirIn, filesep, ['*' fileExtension]]);
filenames={};
[filenames{1:length(directoryContents),1}] = deal(directoryContents.name);
if exist(dirOut1, 'dir')
    rmdir(dirOut1, 's');
end

mkdir(dirOut1);

amount = length(filenames);
% filenames = cell2mat(filenames);
f = figure('Visible','off');
f.Position = [100 50 740 570];
for file_ord = 1:amount
    Frame = file_ord;
    loadfile = sprintf('Strain_%g_u_r.mat', Frame);
    load(loadfile)
    A = reshape(R,1,[])*6.4;
    B = reshape(u_r,1,[])*6.4;
    C = reshape(u_phi,1,[]);
    phi = atan2(Yp, Xp);
    Phi = phi;
    Phi(abs(Phi)>pi/3) = NaN;
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
    D = reshape(phi,1,[]);
    D = mod(D,2*pi);
    [A,I] = sort(A);
    B = B(I);
    C = C(I);
    D = D(I);
    
    A = A(~isnan(B));
    B = B(~isnan(B));
    C = C(~isnan(B));
    D = D(~isnan(B));
    phiii = 30;
    last_index = ceil(length(A)/(phiii/2))-1;
    Rmax = ceil(max(A(~isnan(B)))/last_index)*last_index;
    R_mean = NaN(1, last_index);
    u_rr_mean = NaN(1, last_index);
    phi_r = NaN(phiii, last_index);
    u_phi = NaN(phiii, last_index);
    A_index = 0;
    for r_index = 1:last_index-1
        A_r = A(A_index+1:A_index+phiii);
        B_r = B(A_index+1:A_index+phiii);
        C_r = C(A_index+1:A_index+phiii);
        D_r = D(A_index+1:A_index+phiii);
        [D_r,J] = sort(D_r);
        A_r = A_r(J);
        B_r = B_r(J);
        C_r = C_r(J);
        R_mean(r_index) = mean(A_r);
        u_rr_mean(r_index) = mean(B_r);
        u_phi(:, r_index) = C_r';
        phi_r(:, r_index) = D_r';
    
        A_index = A_index + phiii/2;
    end
    A_r = A(A_index+1:end);
    B_r = B(A_index+1:end);
    C_r = C(A_index+1:end);
    D_r = D(A_index+1:end);
    [D_r,J] = sort(D_r);
    A_r = A_r(J);
    B_r = B_r(J);
    C_r = C_r(J);
    R_mean(last_index) = mean(A_r);
    u_rr_mean(last_index) = mean(B_r);
    u_phi(1:length(C_r), last_index) = C_r';
    phi_r(1:length(D_r), last_index) = D_r';
    
    phi_prog = phi_r(2:end,:);
    phi_prog(end+1,:) = phi_r(1,:) + 2*pi;
    phi_pre = phi_r(end,:) - 2*pi;
    phi_pre(end+1:end+phiii-1,:) = phi_r(1:end-1,:);
    u_phi_prog = u_phi(2:end,:);
    u_phi_prog(end+1,:) = u_phi(1,:);
    u_phi_pre = u_phi(end,:);
    u_phi_pre(end+1:end+phiii-1,:) = u_phi(1:end-1,:);
    
    delta_u_phi_prog = (u_phi_prog - u_phi);
    delta_phi_prog = (phi_prog - phi_r);
    delta_u_phi_pre = (u_phi - u_phi_pre);
    delta_phi_pre = (phi_r - phi_pre);
    
    u_phiphi1 = (delta_u_phi_prog./delta_phi_prog + delta_u_phi_pre./delta_phi_pre)/2;
    u_phiphi = u_phiphi1./R_mean + u_rr_mean./R_mean;
    
    mean_uphiphi = mean(u_phiphi,'omitnan');
    mean_uphiphi1 = mean(u_phiphi1,'omitnan');
    u_rr = gradient(u_rr_mean, R_mean);
    
    R_spline = (0:1000)/1000*max(R_mean);
    u_phiphi1_spline = csaps(R_mean, mean_uphiphi1, 0.0001, R_spline);
    u_phiphi_spline = csaps(R_mean, mean_uphiphi, 0.0001, R_spline);
    
    clf
    plot(R_mean/Rmax, u_rr,'Marker','pentagram','MarkerSize',9, 'LineStyle','None','LineWidth',2,'DisplayName','u_{RR}')
    hold on
    plot(R_mean/Rmax, u_rr_mean./R_mean,'Marker','diamond','MarkerSize',9, 'LineStyle','None','LineWidth',2,'DisplayName','u_{R}/R')
    plot(R_mean/Rmax, mean_uphiphi,'Marker','hexagram','MarkerSize',9, 'LineStyle','None','LineWidth',2,'DisplayName','u_{\phi\phi}')
    legend
    tt = sprintf('Frame = %g', Frame+199);
    title(tt, 'FontSize',40)
    xlabel('R/R_{disk}')
    ylabel('u_{RR}')
%     ymin = floor(min(uphi_v_r(:,2*Frame))/0.5)*0.5;
%     ymax = ceil(max(uphi_v_r(:,2*Frame))/0.5)*0.5;
    xlim([0 1])
%     ylim([ymin 0.5])
    fontsize(f,20,"points")
    fontname(f, 'Cambria Math')
    nametiff = sprintf('u_phiphi_vs_R_%g.tiff', Frame);
    figfile_tiff = fullfile(dirOut1, nametiff);
    saveas(f, figfile_tiff)
end
