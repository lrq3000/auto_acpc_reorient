function autoreorient(inputpath, mode)
    input_vol = spm_vol(strtrim(inputpath));
    smoothed_vol = spm_vol(inputpath);
    smoothed_vol2 = spm_smoothto8bit(smoothed_vol);
    % Y = spm_read_vols(smoothed_vol2);  % direct access to the data matrix
    % imagesc(Y(:,:,20))  % show a slice
    template_vol = spm_vol(strtrim('C:\git\auto_acpc_reorient\T1_template_CAT12_rm_withskull.nii'));
    if strcmp(mode, 'affine')
        fprintf('Affine reorientation\n');
        flags = struct('sep', 5, 'regtype', 'rigid');
        M = eye(4);
        % Coarse affine coregistration
        [M, scal] = spm_affreg(template_vol, smoothed_vol2, flags, M);
        % Fine affine coregistration (sep divided by 2 and using rescaling)
        flags.sep = flags.sep / 2;
        [M, scal] = spm_affreg(template_vol, smoothed_vol2, flags, M, scal);
        MTransform = template_vol.mat \ M * input_vol.mat;  % from spm_affreg doc
        % Apply the reorientation and save back the new coordinations
        spm_get_space(inputpath, MTransform);
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
% hMRI also implements an auto-reorient: https://www.researchgate.net/figure/Left-After-installation-Section-31-the-hMRI-toolbox-can-be-started-from-the-SPM-menu_fig1_330530392
