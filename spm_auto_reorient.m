function spm_auto_reorient(p,img_type,smooth_factor)

% FORMAT spm_auto_reorient(p,img_type,smooth_factor)
% Automatically reorient a T1 image (or any other usual image modality) in the MNI space
% This uses a non-linear coregistration on a template,
% although the reorientation only applies a rigid-body transform.
%
% This is useful as group analysis rely on between-subjects coregistration
% which is for most methods, such as in "unified segmentation", sensitive
% to initial conditions (the starting orientation of the image).
%
% It is advised to check (and fix if necessary) manually the result.
%
% - p       : filepath of NIfTI image to reorient (as `ls` returns).
% - img_type    : template image type 'T1group' (default), 'T1canonical', 'T1', 'T2', 'PET', 'EPI',...
%               i.e. any of the templates provided by SPM/canonical. 'T1group' is a custom template we provide,
%               computed as the average T1 from normalized T1 images from 10 subjects with intensity normalization
%               and without skull stripping using CAT12 (rm* file). 'T1canonical' is the T1 brain image of a single subject
%               as provided by SPM. Note that you can adjust the respective template if you want to tweak
%               the target realignment to your liking. Special option: can also use the path to a nifti file
%               (this allows for cross-modality coregistration, no smoothing is applied here then, in this case
%               it is advised to use mode 'mi' only).
% - smooth_factor : smoothing kernel (isotropic). Default: 20. Usually, a big kernel helps in coregistering
%             to template, particularly for brain damaged patients, but might also help with healthy volunteers.
%             However, a too big kernel will also result in suboptimal coregistration. 20 is good for T1.
%
% Returns: nothing, but the input image's headers are modified.
%__________________________________________________________________________

% Code originally written by John Ashburner & Carlton Chu, FIL, UCL, London, UK
% Extended by Stephen Karl Larroque, Coma Science Group & GIGA-Consciousness, University Hospital of Liege, Belgium

% Init input variables
if ~exist('p', 'var')
    p = spm_select(inf, 'image', 'Select nifti files to reorient');
end
p = char(p);
imgcount = size(p,1);

if ~exist('img_type', 'var')
    img_type = 't1group';
end
img_type = char(img_type);
if size(img_type, 1) ~= imgcount
    error('Invalid number of template images, the number template does not match the number of input images (if you want to apply one different template for each input image, else you can just use one template for all reorientations).');
end


if ~exist('smooth_factor', 'var')
    smooth_factor = 20;
end

%% Select template to reorient to
img_type = lower(img_type);
% Manage special cases (ie, template name is too long so we allow a shorthand)
if strcmp(img_type, 't1canonical')
    img_template = fullfile(spmDir, 'canonical', 'single_subj_T1.nii')
elseif strcmp(img_type, 't1group')
    img_template = fullfile(spmDir, 'canonical', 'T1_template_CAT12_rm_withskull.nii');  % you need to add this file into spm/canonical
    if ~exist(img_template, 'file') == 2  % if template cannot be found in spm folder, try to look locally, in same folder as current script
        % Build the path to current script (because pwd is unreliable)
        scriptpath = mfilename('fullpath');
        scriptdir = fileparts(scriptpath); % get the parent directory of the current script
        img_template = fullfile(scriptdir, 'T1_template_CAT12_rm_withskull.nii');  % this file needs to be in the same folder as this script
        if ~exist(img_template, 'file') == 2
            error('Cannot find template t1group, please make sure the nifti file is at the appropriate place (see readme)')
        end %endif
    end %endif
elseif exist(img_type(1,:), 'file') == 2
    % The template is directly a file (or list of files), we use it as a template (can be used to coregister across modalities, eg EPI BOLD on T1)
    % Note that then we will use img_type as a cellarr of filepaths to the templates to use
    img_template = 'file';
else
    img_template = fullfile(spmDir, 'toolbox', 'OldNorm', [upper(img_type) '.nii']);
    % if template cannot be found in spm folder, try to use the provided img_type as a filepath
    if exist(img_template, 'file') ~= 2
        error('Cannot find template %s, please make sure the nifti file can be found either in SPM/toolbox/OldNorm or you can provide a full filepath.', img_template)
    end %endif
end %endif
vg = spm_vol(img_template);  % get template image

% Configure coregistration to template (will be the basis of the reorientation)
flags.sep = 5;  % sampling distance. Reducing this enhances a bit the reorientation but significantly increases processing time.
flags.regtype = 'mni';  % can be 'none', 'rigid', 'subj' or 'mni'. On brain damaged patients, 'mni' seems to give the best results (non-affine transform), but we don't use the scaling factor anyway. See also a comparison in: Liu, Yuan, and Benoit M. Dawant. "Automatic detection of the anterior and posterior commissures on MRI scans using regression forests." 2014 36th Annual International Conference of the IEEE Engineering in Medicine and Biology Society. IEEE, 2014.

% For each input image
for i = 1:imgcount
    % get image, create a temporary file (to avoid permission issues) and smooth to ease coregistration to template
    f = strtrim(p(i,:));
    spm_smooth(f,'temp.nii',[smooth_factor smooth_factor smooth_factor]);
    vf = spm_vol('temp.nii');
    % calculate the reorientation matrix from source image to match template image
    [M, scal] = spm_affreg(vg,vf,flags);
    M3 = M(1:3,1:3);
    [u s v] = svd(M3);
    M3 = u*v';
    M(1:3,1:3) = M3;
    % reload source image to apply the transform on it
    N = nifti(f);
    N.mat = M*N.mat;
    create(N);
end

end

