% hard-coded paths for masks and images
epi_stroma_mask_path = "../../ovarian_cancer_results/epithelium_stroma_masks/TCGA-OY-A56P_33000_45000.png"
nuclei_mask_path = "../../ovarian_cancer_results/nuclei_masks/TCGA-OY-A56P_33000_45000.png"
rbc_path = "../../ovarian_cancer_results/rbc_masks/TCGA-OY-A56P_33000_45000.png"
patch_path = "../../ovarian_cancer_results/patches/TCGA-OY-A56P_33000_45000.png"

%% parameters setting, for different task, an appropriate  winSize/filterScale value is also different
epi_stroma_mask = imread(epi_stroma_mask_path);
nuclei_mask = imread(nuclei_mask_path);
rbc_mask = imread(rbc_path);
current_patch = imread(patch_path);
winSize = 200; %size of windows from which the collagen fiber disorder was measured
filter_scale = 3; %kernel size of the BIF model
orientCooccurScheme = 1; %considering the area of each detected collagen fiber when quatifying the orientation disorder.
featureDescriptor = 6; %the 6th category of the Basic Image Feature corresponds to the dark linear structure.
orientBinInterval = 10; %discritize the continous orientation angle by interval of 10. 
orientNum = 180 / orientBinInterval;

%% extract collagen fiber mask
frag_thresh = filter_scale*15; %remove detected collagen fragments with an area lower than the predefined threshold
[bifs] = compute_bifs(current_patch, filter_scale, .1, 1); %use BIF based model to extract the linear structures which was used as the representative of collagen fibers
collagen_mask = bifs == featureDescriptor;
[height, width] = size(collagen_mask);       
collagen_mask = (collagen_mask & (1 - epi_stroma_mask)); %remove the detected collagen fibers in nonstroma region
collagen_mask = (collagen_mask & (1 - nuclei_mask)); %remove the detected collagen fibers in stroma nuclei and epithelial nuclei
collagen_mask = (collagen_mask & (1 - rbc_mask));
collagen_mask = bwareaopen(collagen_mask, frag_thresh);
imshow(labeloverlay(current_patch, collagen_mask, 'transparency', 0, 'Colormap', [0,0,1])) %overlay the collagen fiber mask on top of the tumor sample.
hold on

%% collagen centroid and orientation information extraction
collogen_props = regionprops('table', collagen_mask, 'Centroid', 'Orientation', 'Area');
colg_center = collogen_props.Centroid;
colg_area = collogen_props.Area;
colg_orient = collogen_props.Orientation;

%% arrows to show the orientation of the collagen fibers
for colg_ind=1:length(colg_area)
    u= colg_area(colg_ind) * cosd(colg_orient(colg_ind));
    v= colg_area(colg_ind) * sind(colg_orient(colg_ind));
    quiver(colg_center(colg_ind, 1), colg_center(colg_ind,2), u, -v, 'color', 'b', 'LineWidth', 1, 'AutoScaleFactor', 0.25, 'MaxHeadSize', 1)
    hold on
    quiver(colg_center(colg_ind,1), colg_center(colg_ind, 2), -u, v, 'color', 'b', 'LineWidth', 1, 'AutoScaleFactor', 0.25, 'MaxHeadSize', 1)
    hold on
end

%% CFOD feature extraction    
%colgOrientBin=fix(colgOrient/orientBinInterval);
%colgOrientBin=colgOrientBin+9;
%win_size_ind=0;  
%stepSize=winSize/2; % cfod features were extracted based on a sliding window fashion, 
%disp('calculating CFOD')
%win_x_ind=0;
%for win_x=1:stepSize:width-winSize+1 % a window slides across the tumor sample, a set of orientation disorder related features were then extracted from the collagen fibers in each window
%    win_x_ind=win_x_ind+1;
%    win_y_ind=0;
%    for win_y=1:stepSize:height-winSize+1
%       win_y_ind=win_y_ind+1;
%       p_orient_occur=[];
%       strRatioMap(win_y_ind,win_x_ind)=length(find(strMask(win_y:win_y+winSize-1,win_x:win_x+winSize-1)))/(winSize^2);
%       inwinColgInd=find(colgCenter(:,1)>=win_x & colgCenter(:,1)<win_x+winSize-1 & colgCenter(:,2)>=win_y & colgCenter(:,2)<win_y+winSize-1); 
%       inwinColgOrient=colgOrientBin(inwinColgInd); 
%       inwinColgArea=colgArea(inwinColgInd);
%       if length(inwinColgOrient)>=5
%        [orientOccurFeats]=disorder_feat_extract(inwinColgOrient,inwinColgArea,orientNum,orientCooccurScheme); % disorder feature extraction    
%         if isfield(orientOccurFeats,'val')
%            cfodMap(win_y_ind,win_x_ind,:)=orientOccurFeats.val; % 13 feature map for the tumor sample
%         end 
%       end
%   end   
%end

%cfodMap2=[];
%for featInd=1:13
%cfodMapMoveAvg=movmean(cfodMap(:,:,featInd),2,1,'omitnan','Endpoint','discard'); % feature was calculated in sliding window version
%cfodMap2(:,:,featInd)=movmean(cfodMapMoveAvg,2,2,'omitnan','Endpoint','discard');
%end
%% after acquisation of the feature map, a set of statistics e.g. mean, std, sknewness could be calculated. 