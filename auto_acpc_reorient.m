function auto_acpc_reorient(imgpath, img_type, imgpath_other, mode, smooth_factor, flags_affine, flags_mi)
% Cross-platform automatic AC-PC realignment/reorientation and coregistration
% for both healthy volunteers and brain damaged patients using template matching
% using SPM 12.
%
% FORMAT auto_acpc_reorient(imgpath, img_type, imgpath_other, mode, smooth_factor, flags_affine, flags_mi)
%
% imgpath       - filepath or chararray of filepaths of NIfTI images to reorient (as `ls` returns).
%               For 4D nifti files, please select only the first volume (eg, p='bold.nii,1'), and not the others (they will also be reoriented).
% img_type      - template image type 'T1group' (default), 'T1', 'T2', 'PET', 'EPI',...
%               i.e. any of the templates provided by SPM/toolbox/OldNorm. 'T1group' is a custom template we provide,
%               computed as the average T1 from normalized T1 images from 10 subjects with intensity normalization
%               and without skull stripping using CAT12 (rm* file). Note that you can adjust the respective templates
%               if you want to tweak the target realignment to your liking. Special option: can also use the path to a
%               nifti file (this allows for cross-modality coregistration, no smoothing is applied here then, in this case
%               it is advised to use mode 'mi' only).
% imgpath_other - cell array of chararrays filenames of other images to be reoriented with the same transform as the input imgpath image.
%               The cell array should be of the same length as imgpath (but the chararrays can be of arbitrary size), or can be set empty with [].
%               imgpath_other should only include OTHER files than the one specified in imgpath, as SPM12 has a peculiar way to handle
%               the reorientaton of 4D NIfTI files, since the orientation matrix seem to be stored only once in the headers, reorienting the first volume
%               will also reorient all other volumes. Hence, if you don't have any other NIfTI file (and not volume) to reorient, simply set imgpath_other=[].
% mode          - coregister using the old 'affine' euclidian method, or the new 'mi' Mutual Information on Joint Histogram method or 'both' (first affine then mi) (default)
% smooth_factor - smoothing kernel (isotropic) for the affine coregistration. Default: 20. Usually, a big kernel helps in coregistering
%               to template, particularly for brain damaged patients, but might also help with healthy volunteers.
%               However, a too big kernel will also result in suboptimal coregistration. 20 is good for T1.
% flags_affine  - provide your custom flags for the affine coregistration
% flags_mi      - provide your custom flags for the mutual information coregistration
%
% Returns: nothing, but the input image's headers are modified.
% _________________________________________________________________________
%
% This uses a non-linear coregistration on a template,
% although the reorientation only applies a rigid-body transform.
% This supports both euclidian (spm_affreg) and joint histogram
% (spm_coreg) methods, combined or alone.
%
% This is useful as group analysis rely on between-subjects coregistration
% which is for most methods, such as in "unified segmentation", sensitive
% to initial conditions (the starting orientation of the image).
%
% If mode 'mi' is selected, the origin will also be changed to match AC.
%
% It is advised to check (and fix if necessary) manually the result.
%__________________________________________________________________________
% v1.4.6
% License: GPL (General Public License) v3 or later, except otherwise noted in comments around the code the other license pertains to
% Copyright (C) 2008 John Ashburner, FIL, UCL, London, UK (source: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=SPM;d1f675f1.0810 )
% Copyright (C) 2008 Carlton Chu, FIL, UCL, London, UK (source: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=SPM;d1f675f1.0810 )
% Copyright (C) 2019-2020 Stephen Karl Larroque, Coma Science Group & GIGA-Consciousness, University Hospital of Liege, Belgium
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License along
%    with this program; if not, write to the Free Software Foundation, Inc.,
%    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

% Init input variables
if ~exist('imgpath', 'var')
    imgpath = spm_select(inf, 'image', 'Select nifti files to reorient');
end
imgpath = char(imgpath);
imgcount = size(imgpath,1);

if ~exist('img_type', 'var')
    img_type = 't1group';
end
img_type = char(img_type);
if size(img_type, 1) ~= imgcount
    error('Invalid number of template images, the number template does not match the number of input images (if you want to apply one different template for each input image, else you can just use one template for all reorientations).');
end

if ~exist('imgpath_other', 'var')
    imgpath_other = [];
end
if ischar(imgpath_other)
    % imgpath_others should be a cellarray of chararrays, so with one such chararray per input imgpath (so that we can provide a list of images to reorient)
    imgpath_other = {imgpath_other};
end
if ~isempty(imgpath_other) & (size(imgpath_other, 1) ~= imgcount)
    error('Invalid number of other images, does not match the number of source images.');
end

if ~exist('mode', 'var')
    mode = 'both';
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

% Get SPM directory path
spmDir = which('spm');
spmDir = spmDir(1:end-5);
oldNormDir = fullfile(spmDir, 'toolbox/OldNorm');

%% Select template to reorient to
img_type = lower(img_type);
% Manage special cases (ie, template name is too long so we allow a shorthand)
if strcmp(img_type, 't1group')
    img_template = fullfile(spmDir, 'canonical', 'T1_template_CAT12_rm_withskull.nii');  % you need to add this file into spm/canonical
    if exist(img_template, 'file') ~= 2  % if template cannot be found in spm folder, try to look locally, in same folder as current script
        % Build the path to current script (because pwd is unreliable)
        scriptpath = mfilename('fullpath');
        scriptdir = fileparts(scriptpath); % get the parent directory of the current script
        img_template = fullfile(scriptdir, 'T1_template_CAT12_rm_withskull.nii');  % this file needs to be in the same folder as this script
        if exist(img_template, 'file') ~= 2
            error('Cannot find template t1group, please make sure the nifti file is at the appropriate place (see readme)')
        end %endif
    end %endif
elseif exist(img_type(1,:), 'file') == 2
    % The template is directly a file (or list of files), we use it as a template (can be used to coregister across modalities, eg EPI BOLD on T1)
    % Note that then we will use img_type as a cellarr of filepaths to the templates to use
    img_template = 'file';
else
    img_template = fullfile(oldNormDir, [upper(img_type) '.nii']);
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
            Vtemplate = spm_vol(img_type(i, :));
        else
            Vtemplate = spm_vol(img_template);
        end %endif
        % Load source image to reorient to template
        source = strtrim(imgpath(i,:));
        Vsource = spm_vol(source);
        % Smooth to ease coregistration to template
        Vsourcesmoothed = zeros(Vsource.dim(1:3));  % prevent spm_smooth() from saving the smoothing back to the file, by creating a temporary variable where to store the smoothed image
        spm_smooth(Vsource, Vsourcesmoothed, [smooth_factor smooth_factor smooth_factor]);  % TODO: should update to spm_smoothkern()? Check spm_realign()
        % Put the smoothed data back to the original struct, from https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1808&L=spm&P=R27983&1=spm&9=A&I=-3&J=on&d=No+Match%3BMatch%3BMatches&z=4 and https://github.com/spm/spm12/blob/r7219/spm_spm.m#L522
        Vsource.dat = Vsourcesmoothed;
        Vsource.dt = [spm_type('float64') spm_platform('bigend')];  % necessary to make the data readable in-memory
        Vsource.pinfo = [1 0 0]';  % necessary to make the data readable in-memory
        % Calculate the reorientation matrix from source image to match template image
        %%% Content licensed under CC-BY-SA 3.0 by John Ashburner and Carlton Chu : https://en.m.wikibooks.org/wiki/Special:MobileDiff/2764852
        [M, scal] = spm_affreg(Vtemplate, Vsource, flags);
        M3 = M(1:3,1:3);
        [u s v] = svd(M3);
        M3 = u*v';
        M(1:3,1:3) = M3;
        %%% End of CC-BY-SA 3.0 licensed content
        % Memorize to apply on other images later
        M_aff_mem{i} = M;
        % Reload source image to apply the transform on it
        %%% Content licensed under CC-BY-SA 3.0 by John Ashburner and Carlton Chu : https://en.m.wikibooks.org/wiki/Special:MobileDiff/2764852
        N = nifti(source);
        N.mat = M*N.mat;
        % Save the transform into nifti file headers
        create(N);
        %%% End of CC-BY-SA 3.0 licensed content
    end %endfor
end %endif

% JOINT HISTOGRAM (MUTUAL INFORMATION) COREGISTRATION
M_mi_mem = {};
if strcmp(mode, 'mi') | strcmp(mode, 'both')
    fprintf('Mutual information reorientation, please wait...\n');
    % Configure coregistration
    if ~isempty(flags_mi)
        flags2 = flags_mi;
    else
        flags2.cost_fun = 'ecc';  % ncc works remarkably well, when it works, else it fails very badly... Also ncc should only be used for within-modality coregistration (TODO: test if for reorientation it works well, even on very damaged/artefacted brains?)
        flags2.tol = [0.02, 0.02, 0.02, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001];  % VERY important to get good results, these are defaults from the GUI
    end
    %% Treat each image p at a time
    for i = 1:imgcount
        % Load template image
        if strcmp(img_template, 'file')
            Vtemplate = spm_vol(img_type(i,:));
        else
            Vtemplate = spm_vol(img_template);
        end %endif
        % Load source image to reorient to template
        % NB: no need for smoothing here since the joint histogram is smoothed
        source = strtrim(imgpath(i, :));
        Vsource = spm_vol(source);
        % Calculate the reorientation matrix from source image to match template image
        M_mi = spm_coreg(Vtemplate,Vsource,flags2);
        % Memorize to apply on other images later
        M_mi_mem{i} = M_mi;
        % Reload source image to apply the transform on it
        N = nifti(source);
        N.mat = spm_matrix(M_mi) \ N.mat;
        % Save the transform into nifti file headers
        create(N);
    end %endfor
end %endif

% Apply the reorientation transform onto other images (if specified), without recalculating, so that we keep motion information if any
if ~isempty(imgpath_other)
    fprintf('Applying transform to other images...\n');
    for i = 1:imgcount
        % Load the appropriate transforms
        if strcmp(mode, 'affine') | strcmp(mode, 'both'), M = M_aff_mem{i}; end
        if strcmp(mode, 'mi') | strcmp(mode, 'both'), M_mi = M_mi_mem{i}; end
        % For each other image
        for j = 1:size(imgpath_other{i}, 1);
            % Get file path
            if iscell(imgpath_other{i})
                source_other = strtrim(imgpath_other{i}{j});
            elseif ischar(imgpath_other{i})
                source_other = strtrim(imgpath_other{i}(j, :));
            else
                error('Malformatted imgpath_other, please check the structure used (needs to be a cellarray).');
            end %endif
            if ~isempty(source_other) && ~strcmp(source, source_other)  % If filepath is empty or same as source functional, just skip
                % Load volume
                N = nifti(source_other);
                if strcmp(mode,'affine') | strcmp(mode,'both')
                    % Apply affine transform
                    N.mat = M*N.mat;
                end %endif
                if strcmp(mode, 'mi') | strcmp(mode, 'both')
                    % Apply Mutual Information rigid-body transform
                    N.mat = spm_matrix(M_mi) \ N.mat;
                end %endif
                % Save the transform into nifti file headers
                create(N);
            end %endif
        end %endfor
    end %endfor
end %endif

fprintf('Automatic reorientation done.\n');

end %endfunction
