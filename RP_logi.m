function RP_logi

%%
%
%
%
%

%% Identify subjects
HomeDir = '/media/HDPC-UT/dMRI_data';

RP1_1 = 'RP1-TT-2013-11-01';
RP1_2 = 'RP1-TT-20150623';

RP2_1 = 'RP2-KI-2013-11-01';
RP2_2 = 'RP2-KI-20150531';
RP2_3 = 'RP2-KI-20160326';

RP8_1 = 'RP8-YT-2014-03-14-dMRI-Anatomy';
RP8_2 = 'RP8-YT-20150919';

JMD1_1 = 'JMD1-MM-20121025-DWI';
JMD1_2 ='JMD1-MM-20150628';

JMD2_1 = 'JMD2-KK-20121025-DWI';
JMD2_2 = 'JMD2-KK-20150623';
JMD2_3 = 'JMD2-KK-2016-02-13';

JMD3_1 = 'JMD3-AK-20121026-DWI';
JMD3_2 = 'JMD3-AK-20140228-dMRI';
% JMD3_3 = 'JMD3-AK-20150726'; % measured once only

JMD4_1 ='JMD4-AM-20121026-DWI';
JMD4_2 ='JMD4-AM-20150531';

LHON2_1 = 'LHON2-SO-20121130-DWI';
LHON2_2 = 'LHON2-SO-20151013';

LHON4_1 = 'LHON4-GK-dMRI-2014-11-25';
LHON4_2 = 'LHON4-GK-20150628';

LHON5_1 = 'LHON5-HS-20121220-DWI';
LHON5_2 = 'LHON5-HS-IDBN-20160516';

LHON7_1 = 'LHON7-TT-2014-12-20';
LHON7_2 = 'LHON7-TT-2nd-20150222';

% cur_Pt ={RP1_1,RP1_2,RP2_1,RP2_2,RP2_3,RP8_1,RP8_2,JMD1_1,JMD1_2,...
%     JMD2_1,JMD2_2,JMD2_3,JMD3_1,JMD3_2,JMD3_3,JMD4_1,JMD4_2,...
%     LHON2_1,LHON2_2,LHON4_1,LHON4_2,LHON5_1,LHON5_2,LHON7_1,LHON7_2};

% JMD3_3 removed
cur_Pt ={RP1_1,RP1_2,RP2_1,RP2_2,RP2_3,RP8_1,RP8_2,JMD1_1,JMD1_2,...
    JMD2_1,JMD2_2,JMD2_3,JMD3_1,JMD3_2,JMD4_1,JMD4_2,...
    LHON2_1,LHON2_2,LHON4_1,LHON4_2,LHON5_1,LHON5_2,LHON7_1,LHON7_2};

%%

parfor ii = 22:length(cur_Pt)
    % cur_dwi1st = dtiLoadDt6(fullfile(HomeDir,cur_Pt{ii},'dwi_1st','dt6.mat'));
    % cur_dwi2nd = dtiLoadDt6(fullfile(HomeDir,cur_Pt{ii},'dwi_2nd','dt6.mat'));
    
    cd(fullfile(HomeDir,cur_Pt{ii},'raw'))
    
    % need a white matter Mask
%     if ~exist('wmMask.nii.gz','file')
        cd ../dwi_1st/bin
        !cp wmMask.nii.gz ../../raw
        cd(fullfile(HomeDir,cur_Pt{ii},'raw'))
%     end
    
%     %% command run Osmosis-dti-cod with wmMask
%     !osmosis-dti-cod.py dwi1st_aligned_trilin.nii.gz dwi1st_aligned_trilin.bvecs dwi1st_aligned_trilin.bvals dwi2nd_aligned_trilin.nii.gz dwi2nd_aligned_trilin.bvecs dwi2nd_aligned_trilin.bvals dti_cod_wmMask.nii.gz --mask_file wmMask.nii.gz
%     
%     %% Have a pair of dMRI scan, command run Osmosis-dti-rsquered with wmMask
%     !osmosis-dti-rsquared.py dwi1st_aligned_trilin.nii.gz dwi1st_aligned_trilin.bvecs dwi1st_aligned_trilin.bvals dwi2nd_aligned_trilin.nii.gz dwi2nd_aligned_trilin.bvecs dwi2nd_aligned_trilin.bvals dti_rsquared_wmMask.nii.gz --mask_file wmMask.nii.gz
%     
%     %% Have a dMRI scan command run Osmosis-dti-rsquared SO
%     !osmosis-dti-rsquared.py dwi1st_aligned_trilin.nii.gz dwi1st_aligned_trilin.bvecs dwi1st_aligned_trilin.bvals dwi1st_aligned_trilin.nii.gz dwi1st_aligned_trilin.bvecs dwi1st_aligned_trilin.bvals dti1st_rsquared_wmMask.nii.gz --mask_file wmMask.nii.gz
    
    %% command run osmosis-dti-rrmse.py with wmMask
    if ~exist('dti_rrmse_wmMask.nii.gz','file')
    !osmosis-dti-rrmse.py dwi1st_aligned_trilin.nii.gz dwi1st_aligned_trilin.bvecs dwi1st_aligned_trilin.bvals dwi2nd_aligned_trilin.nii.gz dwi2nd_aligned_trilin.bvecs dwi2nd_aligned_trilin.bvals dti_rrmse_wmMask.nii.gz --mask_file wmMask.nii.gz
    end
end
return


%% Intraclass Correlation Coefficient

for ii = 21:length(cur_Pt)

    cd(fullfile(HomeDir,cur_Pt{ii},'raw'))
    
    ni1 = niftiRead('dwi1st_aligned_trilin.nii.gz');
    ni2 = niftiRead('dwi2nd_aligned_trilin.nii.gz'); 
    
    ni1 = niftiRead('dwi1st.nii.gz');
    ni2 = niftiRead('dwi2nd.nii.gz'); 
    
    
    A =  horzcat(ni1.data(:),ni2.data(:))
    
    %% ICC
    type = '1-1';
    M = A;
    
    [r, LB, UB, F, df1, df2, p] = ICC(M, type)
    




%% load cod.nii and create fiberROI.nii
for  ii =1:length(cur_Pt)
    %  load coeffcient of determination
    nicodDir = fullfile(HomeDir,cur_Pt{ii},'raw');
    %     cd(nicodDir)
    nicod = niftiRead(fullfile(nicodDir,'dti_cod_wmMask.nii.gz'));
    nirrmse = niftiRead(fullfile(nicodDir,'dti_rrmse_wmMask.nii.gz'));
    
    % load the fiber group and dt6 files
    fgDir  = fullfile(HomeDir,cur_Pt{ii},'/dwi_1st/fibers/conTrack/OR_100K');
    dt  =  dtiLoadDt6(fullfile(HomeDir,cur_Pt{ii},'dwi_1st','dt6.mat'));
    
    % fiber names need to load
    fgN ={'*Rh_NOT_MD4.pdb','*Lh_NOT_MD4.pdb'};%,'ROTD4L4_1206.pdb','LOTD4L4_1206.pdb',...
    %         'ROCF_D4L4.pdb','LOCF_D4L4.pdb'};
    %     for j = 1: length(fgN)
    fgLName = dir(fullfile(fgDir,fgN{1}));
    fgRName = dir(fullfile(fgDir,fgN{2}));
    fgL = fgRead(fullfile(fgDir,fgLName.name));
    fgR =  fgRead(fullfile(fgDir,fgRName.name));
    %%
    % Now let's get all of the coordinates that the fibers go through
    coords = horzcat(fgL.fibers{:}, fgR.fibers{:});
    
    % get the unique coordinates
    coords_unique = unique(floor(coords'),'rows');
    
    % These coordsinates are in ac-pc (millimeter) space. We want to transform
    % them to image indices.
    img_coords = unique(floor(mrAnatXformCoords(inv(dt.xformToAcpc), coords_unique)), 'rows');
    
    %         % Now we can calculate FA
    %         fa = dtiComputeFA(dt.dt6);
    
    % Now lets take these coordinates and turn them into an image. First we
    % will create an image of zeros
    OR_img = zeros(size(nicod.data));
    
    % Convert these coordinates to image indices
    ind = sub2ind(size(nicod.data), img_coords(:,1), img_coords(:,2),img_coords(:,3));
    
    % Now replace every coordinate that has the optic radiations with a 1
    OR_img(ind) = 1;
    
    %         % % Now you have an image. Just for your own interest if you want to make a
    %         % % 3d rendering
    %         isosurface(OR_img,.5);
    
    % For each voxel that does not contain the optic radiations we will zero
    % out its value
    %         nicod.data(~OR_img) = 0;
    nirrmse.data(~OR_img) = 0;
    
    
    %         nicod.data(nicod.data<0)=0;
    %
    %         showMontage(nicod.data,[],'hot')
    
    showMontage( nirrmse.data,[],'hot')
    
    
    m = mean(nirrmse.data(:));
    s = std(nirrmse.data(:));
    
    
    scatter(ones( size(nicod.data(:))),nicod.data(:))
    
    
    
    % For each voxel that does not contain the optic radiations we will zero
    % out its value
    %         fa(~OR_img) = 0;
    
    % Now we want to save this as a nifti image; The easiest way to do this is
    % just to steal all the information from another image. For example the b0
    % image
    niFg = dtiWriteNiftiWrapper(fa, dt.xformToAcpc, sprintf('%s.nii.gz',fg.name));
end









%% load cod.nii and create fiberROI.nii

for  ii =1:length(cur_Pt)
    %  load coeffcient of determination
    nicodDir = fullfile(HomeDir,cur_Pt{ii},'raw');
    cd(nicodDir)
    nicod = niftiRead(fullfile(nicodDir,'dti_cod.nii.gz'));
    
    % load the fiber group and dt6 files
    fgDir  = fullfile(homeDir,subDir{i},'/dwi_2nd/fibers');
    dt  = fullfile(homeDir,subDir{i},'dwi_2nd','dt6.mat');
    dt  = dtiLoadDt6(dt);
    
    %     'OCFV1V2Not3mm_MD4Al.pdb'
    fgN ={'ROR1206_D4L4.pdb','LOR1206_D4L4.pdb','ROTD4L4_1206.pdb','LOTD4L4_1206.pdb',...
        'ROCF_D4L4.pdb','LOCF_D4L4.pdb'};
    for j = 1: length(fgN)
        fg = fgRead(fullfile(fgDir,fgN{j}));
        
        %%
        % Now let's get all of the coordinates that the fibers go through
        coords = horzcat(fg.fibers{:});
        
        % get the unique coordinates
        coords_unique = unique(floor(coords'),'rows');
        
        % These coordsinates are in ac-pc (millimeter) space. We want to transform
        % them to image indices.
        img_coords = unique(floor(mrAnatXformCoords(inv(dt.xformToAcpc), coords_unique)), 'rows');
        
        % Now we can calculate FA
        fa = dtiComputeFA(dt.dt6);
        
        % Now lets take these coordinates and turn them into an image. First we
        % will create an image of zeros
        OR_img = zeros(size(fa));
        % Convert these coordinates to image indices
        ind = sub2ind(size(fa), img_coords(:,1), img_coords(:,2),img_coords(:,3));
        % Now replace every coordinate that has the optic radiations with a 1
        OR_img(ind) = 1;
        
        % % Now you have an image. Just for your own interest if you want to make a
        % % 3d rendering
        isosurface(OR_img,.5);
        
        % For each voxel that does not contain the optic radiations we will zero
        % out its value
        fa(~OR_img) = 0;
        
        % Now we want to save this as a nifti image; The easiest way to do this is
        % just to steal all the information from another image. For example the b0
        % image
        niFg = dtiWriteNiftiWrapper(fa, dt.xformToAcpc, sprintf('%s.nii.gz',fg.name));
        niftiWrite(niFg,niFg.fname)
    end
end


