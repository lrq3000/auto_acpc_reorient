function spm_auto_reorient(p,img_type,p_others,mode,smooth_factor,flags_affine,flags_mi)

% FORMAT spm_auto_reorient(p,img_type,p_others,mode,smooth_factor)
% Automatically reorient a T1 image (or any other usual image modality) in the MNI space
% This uses a non-linear coregistration on a template,
% although the reorientation only applies a rigid-body transform.
% This supports both euclidian (spm_affreg) and joint histogram
% (spm_coreg) methods.
%
% This is useful as group analysis rely on between-subjects coregistration
% which is for most methods, such as in "unified segmentation", sensitive
% to initial conditions (the starting orientation of the image).
%
% If mode 'mi' is selected, the origin will also be changed to match AC.
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
% - p_others    : cell array of chararrays filenames of other images to be reoriented with the same transform as the input imgpath image.
%               The cell array should be of the same length as imgpath (but the chararrays can be of arbitrary size), or can be set empty with [].
%               imgpath_other should only include OTHER files than the one specified in imgpath, as SPM12 has a peculiar way to handle
%               the reorientaton of 4D NIfTI files, since the orientation matrix seem to be stored only once in the headers, reorienting the first volume
%               will also reorient all other volumes. Hence, if you don't have any other NIfTI file (and not volume) to reorient, simply set imgpath_other=[].
% - mode    : coregister using the old 'affine' method, or the new 'mi' Mutual Information method or 'both' (first affine then mi) (default)
% - smooth_factor : smoothing kernel (isotropic) for the affine coregistration. Default: 20. Usually, a big kernel helps in coregistering
%             to template, particularly for brain damaged patients, but might also help with healthy volunteers.
%             However, a too big kernel will also result in suboptimal coregistration. 20 is good for T1.
% - flags_affine: provide your custom flags for the affine coregistration
% - flags_mi    : provide your custom flags for the mutual information coregistration
%
% Returns: nothing, but the input image's headers are modified.
%__________________________________________________________________________
% v1.3.6
% Licensed under GPL (General Public License) v2
% Code originally written by John Ashburner & Carlton Chu, FIL, UCL, London, UK
% Extended by Stephen Karl Larroque, Coma Science Group & GIGA-Consciousness, University Hospital of Liege, Belgium

% Init input variables
if ~exist('p', 'var')
    p = spm_select(inf, 'image', 'Select nifti files to reorient');
end
p = char(p);
imgcount = size(p,1);
    error('Wrong number of template images, does not match the number of source images to reorient!');
end

if ~exist('img_type', 'var')
    img_type = 't1group';
end
img_type = char(img_type);
if size(img_type, 1) ~= imgcount
    error('Invalid number of template images, the number template does not match the number of input images (if you want to apply one different template for each input image, else you can just use one template for all reorientations).');
end

if size(p_others, 1)~= imgcount
    error('Invalid number of other images to reorient.');
end

if ~exist('smooth_factor', 'var')
    smooth_factor = 20;
end

if ~exist('flags_affine', 'var')
    flags_affine = [];
end

if ~exist('flags_mi', 'var')
    flags_mi = [];
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

% AFFINE COREGISTRATION
M_aff_mem = {};
if strcmp(mode,'affine') | strcmp(mode,'both')
    fprintf('Affine reorientation, please wait...\n');
    % Configure coregistration to template (will be the basis of the reorientation)
    if ~isempty(flags_affine)
        flags = flags_affine
    else
        flags.sep = 5;  % sampling distance. Reducing this enhances a bit the reorientation but significantly increases processing time.
        flags.regtype = 'mni';  % can be 'none', 'rigid', 'subj' or 'mni'. On brain damaged patients, 'mni' seems to give the best results (non-affine transform), but we don't use the scaling factor anyway. See also a comparison in: Liu, Yuan, and Benoit M. Dawant. "Automatic detection of the anterior and posterior commissures on MRI scans using regression forests." 2014 36th Annual International Conference of the IEEE Engineering in Medicine and Biology Society. IEEE, 2014.
    end %endif

    % For each input image
    for i = 1:imgcount
        % Load template image
        if strcmp(img_template, 'file')
            Vtemplate = spm_vol(img_type(i,:));
        else
            Vtemplate = spm_vol(img_template);
        end %endif
        % Load source image to reorient to template
        source = strtrim(p(i,:));
        Vsource = spm_vol(source);
        % smooth to ease coregistration to template
        Vsourcesmoothed = zeros(Vsource.dim(1:3));  % prevent spm_smooth() from saving the smoothing back to the file, by creating a temporary variable where to store the smoothed image
        spm_smooth(Vsource,Vsourcesmoothed,[smooth_factor smooth_factor smooth_factor]);  % TODO: should update to spm_smoothkern()? Check spm_realign()
        % put the smoothed data back to the original struct, from https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1808&L=spm&P=R27983&1=spm&9=A&I=-3&J=on&d=No+Match%3BMatch%3BMatches&z=4 and https://github.com/spm/spm12/blob/r7219/spm_spm.m#L522
        Vsource.dat = Vsourcesmoothed;
        Vsource.dt = [spm_type('float64') spm_platform('bigend')];  % necessary to make the data readable in-memory
        Vsource.pinfo = [1 0 0]';  % necessary to make the data readable in-memory
        % Calculate the reorientation matrix from source image to match template image
        [M, scal] = spm_affreg(Vtemplate,Vsource,flags);
        M3 = M(1:3,1:3);
        [u s v] = svd(M3);
        M3 = u*v';
        M(1:3,1:3) = M3;
        % Memorize to apply on other images later
        M_aff_mem{i} = M;
        % Reload source image to apply the transform on it
        N = nifti(source);
        N.mat = M*N.mat;
        % Save the transform into nifti file headers
        create(N);
    end %endfor
end %endif

% MUTUAL INFORMATION COREGISTRATION
M_mi_mem = {};
if strcmp(mode,'mi') | strcmp(mode,'both')
    fprintf('Mutual information reorientation, please wait...\n');
    % Configure coregistration
    if ~isempty(flags_mi)
        flags2 = flags_mi;
    else
        flags2.cost_fun = 'nmi';  % ncc works remarkably well, when it works, else it fails very badly... Also ncc should only be used for within-modality coregistration (TODO: test if for reorientation it works well, even on very damaged/artefacted brains?)
        flags2.tol = [0.02, 0.02, 0.02, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001];  % VERY important to get good results, these are defaults from the GUI
    end
    % For each input image
    for i = 1:imgcount
        % Load template image
        if strcmp(img_template, 'file')
            Vtemplate = spm_vol(img_type(i,:));
        else
            Vtemplate = spm_vol(img_template);
        end %endif
        % Load source image to reorient to template
        % NB: no need for smoothing here since the joint histogram is smoothed
        source = strtrim(p(i,:));
        Vsource = spm_vol(source);
        % Calculate the reorientation matrix from source image to match template image
        M_mi = spm_coreg(Vtemplate,Vsource,flags2);
        % Memorize to apply on other images later
        M_mi_mem{i} = M_mi;
        % Reload source image to apply the transform on it
        N = nifti(source);
        N.mat = spm_matrix(M_mi)\N.mat;
        % Save the transform into nifti file headers
        create(N);
    end %endfor
end %endif

% Apply the reorientation transform onto other images (if specified), without recalculating, so that we keep motion information if any
if ~isempty(p_others)
    fprintf('Applying transform to other images...\n');
    for i = 1:imgcount
        % Load the appropriate transforms
        if strcmp(mode,'affine') | strcmp(mode,'both'), M = M_aff_mem{i}; end
        if strcmp(mode,'mi') | strcmp(mode,'both'), M_mi = M_mi_mem{i}; end
        % For each other image
        for j = 1:size(p_others{i},1);
            % Get file path
            if iscell(p_others{i})
                source_other = strtrim(p_others{i}{j});
            elseif ischar(p_other{i})
                source_other = strtrim(p_others{i}(j,:));
            else
                error('Malformatted p_others, please check the structure used');
            end %endif
            if ~isempty(source_other) && ~strcmp(source,source_other)  % If filepath is empty or same as source functional, just skip
                % Load volume
                N = nifti(source_other);
                if strcmp(mode,'affine') | strcmp(mode,'both')
                    % Apply affine transform
                    N.mat = M*N.mat;
                end %endif
                if strcmp(mode,'mi') | strcmp(mode,'both')
                    % Apply Mutual Information rigid-body transform
                    N.mat = spm_matrix(M_mi)\N.mat;
                end %endif
                % Save the transform into nifti file headers
                create(N);
            end %endif
        end %endfor
    end %endfor
end %endif

fprintf('Automatic reorientation done.\n');

end %endfunction

