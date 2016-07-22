% This script illustrates how the PORTS texture computation functions are
% used. 
%
% This example is similar to 'example_script_1.m' except that it
% illustrates how to read in a series of DICOM files to create the initial
% image volume and ROI mask and also illustrates 3D linear interpolation 
% to a cubic voxel grid.
%


%%%%%%%%%%%%
%
% Source code developed by :
% The Imaging Research Laboratory - University of Washington
%
% Copyright 2016 Department of Radiology
% University of Washington
% All Right Reserved
% 
%
%%%%%%%%%%%%


%%%%%%%%%%%%
%
% This software is issued without express warranty, no express guarantee of
% fidelty, and the authors are not responsible for the intended or
% unintended results of usage of this software. Quality verification of
% data obtained using PORTS and results drawn from that data are the sole
% responsibility of the end user.
%
% This software is intended for use in whole, and shall not be altered,
% used in part, or modified without full and proper disclosure by end
% parties. 
%
% All publication that use the PORTS software must cite the version number
% and PORTS website: 
%
% https://nciphub.org/groups/ports
% 
%
%%%%%%%%%%%%

%%%%%%%%%%%%
%
% PET Oncology Radiomics Test Suite (PORTS) version 1.1
% 
% 'example_script_3.m' version 1.0 - 18 April 2016
%
% Programmer: Larry Pierce - University of Washington - lapierce@uw.edu
% 
%
%%%%%%%%%%%%


tic % Start clock
%% Data Paths

% Add the matlab path to the PORTS package:
% PORTS_DIR = 'C:/Users/aprilbal/Documents/matlab/sarcoma';
% addpath(PORTS_DIR);

% Add directory name for mask and image locations:
dcm_dir = 'C:/Users/aprilbal/Documents/matlab/sarcoma';
%dcm_img_dir = 'C:/Users/aprilbal/Documents/matlab/sarcoma';


%% Loop through the patients:
for k = 4
% Add paths to the mask and t1gd folders
addpath(strcat(dcm_dir, sprintf('/Sarc%.3d/mask', k)));    
addpath(strcat(dcm_dir, sprintf('/Sarc%.3d/t1gd', k)));    
    
% Now grab names of DICOMs in each patients folders:
mask_file_names = dir(strcat(dcm_dir, sprintf('/Sarc%.3d/mask/*dcm', k)));
img_file_names = dir(strcat(dcm_dir, sprintf('/Sarc%.3d/t1gd/*dcm', k)));


%% Read DICOM files for the image and mask

% We are assuming the image DICOM files are stored in a cell array named
% 'img_file_names' and the mask file names are in a cell called
% 'mask_file_names'.


% Read in the first DICOM header to determine the size of the image volume
% (any of the DICOM files will have this information and should be the same 
% for image and mask):
hi = dicominfo(sprintf('%s/Sarc%.3d/mask/%s', dcm_dir, k, mask_file_names(1).name));

% Record the number of rows, columns, and slices for the image volume:
num_cols = double(hi.Columns);
num_rows = double(hi.Rows);
% *** headers don't have NumberOfSlices
%num_slcs = double(hi.NumberOfSlices);

num_slcs = size(mask_file_names,1);

% Check that the number of slices in the mask is same as in the image vol
if num_slcs ~= size(img_file_names,1)
    error('Number of Files do not match number of DICOM slices -- 2')
end


% Create a placeholder for the mask and image volumes:
mask_vol = zeros(num_rows,num_cols,num_slcs);
img_vol = zeros(num_rows,num_cols,num_slcs);


% We also need to keep track of the z-values of each slice for the linear
% interpolation:
DICOM_z_values = zeros(num_slcs,1);

% Loop over files:
for n = 1:num_slcs

    % The name of this file:
    %this_mask_file_name = mask_file_names{n};
    %this_img_file_name = img_file_names{n};
    this_mask_file_name = mask_file_names(n).name;
    this_img_file_name = img_file_names(n).name;

    % Get the header info for this file:
    hi_mask = dicominfo(this_mask_file_name);
    hi_img = dicominfo(this_img_file_name);

    %% Commented out this part for now.
%     % Test if the image cosines indicate orthogonal voxel faces:
%     if sum(abs(hi_img.ImageOrientationPatient) ~= [1;0;0;0;1;0])
%         error('Image Cosines are not parallel to spatial coordinates')
%     end

    % Record the z-value of this slice:
    DICOM_z_values(hi_img.InstanceNumber) = hi_img.ImagePositionPatient(3);


    % Read the mask and image file, apply DICOM slope/intercept rescaling and 
    % put it into the proper slice within the volume (note these must be 
    % converted to type 'double') :
%     mask_vol(:,:,hi_mask.InstanceNumber) = double(hi_mask.RescaleSlope)*double(dicomread(this_mask_file_name)) + double(hi_mask.RescaleIntercept);
%     img_vol(:,:,hi_img.InstanceNumber) = double(hi_img.RescaleSlope)*double(dicomread(this_img_file_name)) + double(hi_img.RescaleIntercept);

    %*** TODO: Added temp. for scans without RescaleSlope or RescaleIntercept.
    % Essentially assumes images are already rescaled properly
    mask_vol(:,:,hi_mask.InstanceNumber) = double(dicomread(this_mask_file_name));
    img_vol(:,:,hi_img.InstanceNumber) = double(dicomread(this_img_file_name));

end % Loop over DICOM files in directory

% Make sure the mask volume is binary 0/1:
mask_max = max(mask_vol(:));
mask_vol(mask_vol<mask_max) = 0;
mask_vol(mask_vol==mask_max) = 1;


%% *** Commented out this section for non-interpolation. Change names of

interp_mask_vol = mask_vol;
interp_img_vol = img_vol;


%% Texture Parameters


% Choose the number of distinct graytone values to use for computing 
% texture metrics:
num_img_values = 64;



%% Get the names of all metrics into a cell array:

% This step creates a cell structure with all of the metric names. Calling
% a 'compute_*_metrics.m' function without any argument returns a cell 
% structure with the names of the metrics that function computes.


% Put all of the computed metric names into a single list. There are 42
% texture metrics computed in this script:
metric_names = [compute_histogram_metrics() ; ...
                compute_GTSDM_metrics() ; ...
                compute_NGTDM_metrics() ; ...
                compute_zone_size_metrics() ];


% Number the metric names to make it easy on the user:
temp_metric_names = cell(size(metric_names));
for this_metric = 1:size(metric_names,1)
    temp_metric_names{this_metric} = sprintf('(%d) %s',this_metric,metric_names{this_metric});
end

metric_names = temp_metric_names;


% Clear unused variable:
clear temp_metric_names



%% Overhead Computations for the Masks


% The function 'determine_ROI_3D_connectivity.m' is a pre-processing step
% to speed up computation of the co-occurance and neighborhood-dependence
% matrices. It determines which voxels in the mask are connected and how 
% they are connected. It also determines a bounding box around the mask 
% that can be used to speed up computations.


% Determine connectivity and bounding box of this ROI:
[bounding_box,ROI_conn_3D_6,ROI_conn_3D_26,binary_dir_connectivity] = ...
    determine_ROI_3D_connectivity(mask_vol);


% Take the ROI sub-volume within the bounding box:
mask_vol_subvol = mask_vol(bounding_box(1,1):bounding_box(1,2) , ...
                                  bounding_box(2,1):bounding_box(2,2) , ...
                                  bounding_box(3,1):bounding_box(3,2) );      


% Now take the image sub-volume that corresponds to this mask:
img_vol_subvol = img_vol(bounding_box(1,1):bounding_box(1,2) , ...
                                bounding_box(2,1):bounding_box(2,2) , ...
                                bounding_box(3,1):bounding_box(3,2) );      





%% Texture Metric Computation Loop



% Determine the number of voxels in the ROI. This is used to compute the
% histogram-based probabilities and also in the Size-Zone metrics
% computations:
num_ROI_voxels = length(find(mask_vol_subvol));




%%% Discretize the image volume to the desired number of graytones. The
%%% PORTS functions require the image voxel values to be {1,2,3,...,N}. 

% Find the min and max within only the ROI:
img_min = min(img_vol_subvol(logical(mask_vol_subvol))); 
img_max = max(img_vol_subvol(logical(mask_vol_subvol))); 


% Rescale to image volume to [0,N]:
img_vol_subvol = num_img_values .* (img_vol_subvol - img_min)/(img_max - img_min) ;

% Discretize and add 1 to get values {1,2,...,N+1}:
img_vol_subvol = floor(img_vol_subvol) + 1;

% The max value is currently one higher than it should be (N+1), so put 
% those voxels at the max value:
img_vol_subvol(img_vol_subvol==num_img_values+1) = num_img_values;




%%%%%
%%%%% Histogram-based computations:
%%%%%

% Compute the histogram of the ROI and probability of each voxel value:
vox_val_hist = zeros(num_img_values,1);
for this_vox_value = 1:num_img_values
    vox_val_hist(this_vox_value) = length(find((img_vol_subvol == this_vox_value) & (mask_vol_subvol == 1) ));
end

% Compute the relative probabilities from the histogram:
vox_val_probs = vox_val_hist / num_ROI_voxels;


% Compute the histogram_based metrics:
histogram_metrics = compute_histogram_metrics(vox_val_probs,num_img_values);




%%%%%
%%%%% GTDSM (Co-occurance) Matrix calculations:
%%%%%

% Create the Gray-Tone-Spatial-Dependence-Matrix (GTSDM):
GTSDM = compute_3D_GTSDM(mask_vol_subvol,img_vol_subvol,binary_dir_connectivity,num_img_values);


% It is common to compute the mean value over the 13 directions (for
% distance-1 voxels). The following loop computes the metrics for each
% direction:
GTSDM_metrics = zeros(13,19); % There are 19 metrics output from 'compute_GTSDM_metrics.m'
for this_direction = 1:13
    % Compute the metrics for this combination:
    GTSDM_metrics(this_direction,:) = compute_GTSDM_metrics(GTSDM(:,:,this_direction));
end

% Now take the mean texture metric value over the all directions:
GTSDM_metrics = mean(GTSDM_metrics,1);




%%%%%
%%%%% NTGDM Matrices and metrics
%%%%%

% Create the Neighborhood-Gray-Tone-Difference-Matrix(NGTDM):
[NGTDM,vox_occurances_NGD26] = compute_3D_NGTDM(mask_vol_subvol,img_vol_subvol,binary_dir_connectivity,num_img_values);


% Compute NGTDM metrics:
NGTDM_metrics = compute_NGTDM_metrics(NGTDM,num_img_values,vox_occurances_NGD26);




%%%%%
%%%%% Zone Size matrix and Metrics:
%%%%%

% Create the Zone Size Matrix:
GLZSM = compute_GLZSM(mask_vol_subvol,img_vol_subvol,num_img_values);

% Compute the Zone Size Metrics:
GLSZM_metrics = compute_zone_size_metrics(GLZSM,num_ROI_voxels);



%%% Now concatenate all metrics into a single vector:
all_texture_metrics = [ histogram_metrics(:) ; GTSDM_metrics(:) ; NGTDM_metrics(:) ; GLSZM_metrics(:) ];


end

toc % end clock

