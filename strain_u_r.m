%% u_r and u_phi better
clear all
dirIn = 'D:\Active_gel_Project\Data\Gel_23_08_2023\Strain_uf_stuck';
dirOut = 'D:\Active_gel_Project\Data\Gel_23_08_2023\Strain_u_r_u_phi_stuck';
fileExtension = '.mat';
directoryContents = dir([dirIn, filesep, ['*' fileExtension]]);
filenames={};
[filenames{1:length(directoryContents),1}] = deal(directoryContents.name);
if exist(dirOut, 'dir')
    rmdir(dirOut, 's');
end

mkdir(dirOut);

amount = length(filenames);
% filenames = cell2mat(filenames);
for file_ord = 1:amount
    filen = cell2mat(filenames(file_ord));
    % loadfile = fullfile(dirIn,filen);
    lf = fullfile(dirIn, filen);
    if exist(lf,'file')
        load(lf, 'x', 'y', 'u_x', 'u_y', 'uxnew', 'uynew', 'u_xx', 'u_yy', 'u_xy', 'X_CoC', 'Y_CoC','X_CoC2','Y_CoC2')
    else
        continue
    end
    Xp = x-X_CoC;
    Yp = y-Y_CoC;
    phi = atan2(Yp,Xp);
    R = sqrt(Xp.^2 + Yp.^2);
    R(isnan(u_x)) = NaN;
    % phi(abs(phi)>pi/3)=NaN;
    % R(isnan(phi))=NaN;
    Rmax = max(R,[],'all');
    % R(R<3/5*Rmax)=NaN;
    cos_phi = Xp./R;
    sin_phi = Yp./R;
    u_r = uxnew.*cos_phi + uynew.*sin_phi;
    u_phi = -uxnew.*sin_phi + uynew.*cos_phi;
    u_rr = u_xx.*cos_phi.^2 + u_yy.*sin_phi.^2 + 2*u_xy.*sin_phi.*cos_phi;
    u_phiphi = u_xx.*sin_phi.^2 + u_yy.*cos_phi.^2 - 2*u_xy.*sin_phi.*cos_phi;
    u_rphi = -(u_xx - u_yy).*sin_phi.*cos_phi + u_xy.*(cos_phi.^2 + sin_phi.^2);
    % save_file_name = sprintf('Strain_%g_u_r.mat', file_ord);
    matfile = fullfile(dirOut, filen);
    save(matfile, 'u_r', 'u_phi', 'u_x', 'u_y', 'R','Xp','Yp', 'u_xx', 'u_yy', 'u_xy', 'u_rphi', 'u_phiphi', 'u_rr','x','y','X_CoC2','Y_CoC2','uxnew','uynew')
    
end