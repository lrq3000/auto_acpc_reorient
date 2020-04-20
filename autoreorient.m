function autoreorient(inputpath, mode, flags_affine)
    input_vol = spm_vol(strtrim(inputpath));
    smoothed_vol = spm_vol(inputpath);
    smoothed_vol2 = spm_smoothto8bit(smoothed_vol, 20);
    % Y = spm_read_vols(smoothed_vol2);  % direct access to the data matrix
    % imagesc(Y(:,:,20))  % show a slice
    template_vol = spm_vol(strtrim('C:\git\auto_acpc_reorient\T1_template_CAT12_rm_withskull.nii'));
    if strcmp(mode, 'affine')
        fprintf('Affine reorientation\n');
        %flags = struct('sep', 5, 'regtype', 'rigid');
        if exist('flags_affine', 'var') & ~isempty(flags_affine)
            flags = flags_affine;
            if isarray(flags_affine.sep)
                % If we are provided a vector of sampling steps, spm_affreg() does not support multiple sampling steps (contrary to spm_coreg()), so we manage that manually
                sampling_steps = flags_affine.sep;
            else
                % Default sampling steps
                sampling_steps = [10 5];
            end
        else
            % Default flags for affine coregistration
            flags_affine = struct('regtype', 'mni');
            % Default sampling steps
            sampling_steps = [10 5];
        end
        
        % Initialize the transform matrix M (using eye() as is done in spm_affreg())
        M = eye(4);
        scal = 1.0;
        % For each sampling step, coregister the smoothed input volume onto the template by an affine transform
        for s = 1:numel(sampling_steps)
            % Load up the sampling step for this iteration
            flags.sep = sampling_steps(s);
            % Coarse/fine affine coregistration
            [M, scal] = spm_affreg(template_vol, smoothed_vol2, flags, M, scal);  % reuse M, the previous coarse transform, so that we can go further (this avoids the need to apply the reorientation transform in-memory)
            %smoothed_vol2.mat = template_vol.mat \ M * smoothed_vol2.mat;  % apply the reorientation transform in-memory
        end %endfor
        % Calculate the transform from the input to the template (I think that spm_affreg is giving us the inverse: from the template to the input)
        %MTransform = template_vol.mat \ M * input_vol.mat;  % from spm_affreg doc
        % TODO: shearing is applied, decompose the transform to avoid shearing (but allow isotropic scaling?)
        M = template_vol.mat \ M;
        % Apply the reorientation and save back the new coordinations in the original input file
        spm_get_space(inputpath, M*input_vol.mat);
    else
        fprintf('Mutual information reorientation\n');
        flags = struct('sep', 5);
        x = spm_coreg(template_vol, smoothed_vol2, flags);
        M = inv(spm_matrix(x));
        % Apply the reorientation and save back the new coordinations
        spm_get_space(inputpath, M*spm_get_space(inputpath));
    end %endfunction
    fprintf('Autoreorientation done!\n');
end %endfunction

% Partially inspired under fair use by: https://github.com/jimmyshen007/NeMo/blob/c7cc775e7fc84f84255b0d3fa474b7f955e817f5/mymfiles/eve_tools/Structural_Connectivity/approx_coreg2MNI.m - licensed under BSD-2-Clause (compatible with MIT license anyway)
% BESTTUTO: simple calculations (like imcalc but more freedom) in SPM: http://imaging.mrc-cbu.cam.ac.uk/imaging/BasicProgrammingWithSpm
% fieldtrip also implements align_ctf2spm(mri, opt, template)
% hMRI also implements an auto-reorient: https://www.researchgate.net/figure/Left-After-installation-Section-31-the-hMRI-toolbox-can-be-started-from-the-SPM-menu_fig1_330530392 - extension of Christophe Phillips code by Evelyne Balteau, direct link: https://github.com/hMRI-group/hMRI-toolbox/blob/master/hmri_autoreorient.m
