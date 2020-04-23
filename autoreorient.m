function autoreorient(inputpath, mode, flags_affine, noshearing, isorescale, affdecomposition, precoreg, precoreg_reset_orientation, debug)
% DEVELOPMENT FUNCTION
% This function is kept as a minified version of the core routine to do autoreorientation using SPM12 functions. This is kept for development purposes (to quickly debug out only the core routines), do not use it for production.
% affdecomposition defines the decomposition done to ensure structure is maintained (ie, no reflection nor shearing). Can be: qr (default), svd (deprecated), imatrix or none. Only the qr decomposition ensures the structure is maintained.

    if ~exist('noshearing', 'var')
        noshearing = true;
    end
    if ~exist('isorescale', 'var')
        isorescale = true;
    end
    if ~exist('affdecomposition', 'var')
        affdecomposition = 'qr';
    end
    if ~exist('precoreg', 'var')
        precoreg = true;
    end
    if ~exist('debug', 'var')
        debug = false;
    end
    if ~exist('precoreg_reset_orientation', 'var')
        % before coregistration, reset to scanner orientation?
        % value can be: 'raw' or 'scanner' (equals to 'mat0') or false
        precoreg_reset_orientation = false;
    end

    % Load input and template images
    input_vol = spm_vol(strtrim(inputpath));
    template_vol = spm_vol(strtrim('T1_template_CAT12_rm_withskull.nii'));

    if precoreg_reset_orientation ~= false
        if strcmp(precoreg_reset_orientation, 'raw')
            % raw original voxel space with no orientation
            % get voxel size
            vox = sqrt(sum(input_vol.mat(1:3, 1:3).^2));
            % build resetted matrix basis, based simply on the voxel scaling, but reset the whole orientation
            M = diag([vox 1]);
            % get centroid, we will set the origin on it (ie, the translation part of the M orientation matrix)
            input_centroid = get_centroid(input_vol, true, false, debug);
            M(1:3,4) = (-vox .* input_centroid(1:3))';
            % save back into the file
            spm_get_space(inputpath, M);
            input_vol.mat = M;
        elseif strcmp(precoreg_reset_orientation, 'scanner') || strcmp(precoreg_reset_orientation, 'mat0')
            % original orientation set by the scanner
            % simply reload mat0 and overwrite the current orientation matrix
            input_vol.mat = input_vol.private.mat0;
            spm_get_space(inputpath, input_vol.mat);
        end
    end

    % Smooth input image, helps a lot with the coregistration (which is inherently noisy since there is no perfect match)
    smoothed_vol = spm_vol(inputpath);
    smoothed_vol2 = spm_smoothto8bit(smoothed_vol, 20);
    % Y = spm_read_vols(smoothed_vol2);  % direct access to the data matrix
    % imagesc(Y(:,:,20))  % show a slice

    % Manual/Pre coregistration (without using SPM)
    if precoreg
        fprintf('Pre-coregistration by translation of centroid, please wait...\n');
        % Find the center of mass
        input_centroid = [get_centroid(input_vol, true, false, debug) 1];  % add a unit factor to be able to add the translation of the world-to-voxel mapping in the nifti headers, see: https://www.youtube.com/watch?v=UvevHXITVG4
        template_centroid = [get_centroid(template_vol, true, false, debug) 1];
        input_centroid_voxel = (input_vol.mat * input_centroid')';  % convert from raw/scanner voxel space to current orientation voxel space (ie, this represents the distance in voxels from the current origin to the centroid, in other words this is the centroid position relative to the currently set origin)
        template_centroid_voxel = (template_vol.mat * template_centroid')';
        % Calculate the euclidian distance = translation transform from the input to the template centroid, we also resize to match the input volume's space (which may be bigger or smaller - if bigger, then we need to reduce the number of voxels we move)
        i2t_dist = (template_centroid_voxel(1:3) - input_centroid_voxel(1:3)) .* (template_vol.dim ./ input_vol.dim);
        % Apply the translation on the nifti header voxel-to-world transform
        % This effectively resets the origin onto the centroid (nifti viewers such as MRIcron or SPM will apply the transform in a vol.premul matrix before the rest)
        M = input_vol.mat;
        M(1:3,4) = M(1:3,4) - input_centroid_voxel(1:3)';  % set the origin onto the centroid
        %M(1:3,4) = M(1:3,4) + i2t_dist(1:3)';  % shift some more to match where the origin is in the template compared to its own centroid. Not sure this step is necessary since anyway we certainly won't end up in the AC-PC, but well why not, it may bridge some more the relative distance between the template origin and input volume origin. DEPRECATED: unreliable, it's preferable to stick to the centroid
        % Note: at this point, we should also set the template's origin on its centroid to leave only the rotation and scale to be estimated afterwards. But since we provide a template with the origin on the AC-PC, it's an even better origin, so we skip this step here. If the template did not have the origin on AC-PC, then setting the origin on its centroid would be a good thing to do, see: http://nghiaho.com/?page_id=671.
        % Save into the original file
        spm_get_space(inputpath, M);
        if debug, spm_get_space([inputpath(1:end-4) '-centroid.nii'], M); end;
        %spm_get_space([inputpath(1:end-4) '-centroid.nii'], M);  % debug line, save the same transform on the image where the centroid is made visible
        % Apply the new mapping in the already loaded images including smoothed (this avoids the need to reload the input volume and resmooth it)
        input_vol.mat = M;
        smoothed_vol.mat = M;
        smoothed_vol2.mat = M;
        % SPM provides 3 ways to save nifti files:
        % * spm_get_space(vol, transform) to only modify the world-to-voxel mapping headers (vol.mat)
        % * spm_write_vol(vol, vol_data) where vol_data is a 3D matrix of the content of the nifti, and vol a nifti structure as accepted or created by spm using spm_read_vols().
        % * create(vol) where vol is a spm nifti structure as created by nifti(), and it will be saved in vol.fname path. This last option is the most complete, as it will save everything.
        %spm_write_vol(input_vol, input_vol.mat);
        fprintf('Pre-coregistration done!\n');
    end

    % Affine coregistration (with translation)
    if strcmp(mode, 'affine')
        fprintf('Affine reorientation\n');
        if exist('flags_affine', 'var') && ~isempty(flags_affine)
            if isarray(flags_affine.sep)
                % If we are provided a vector of sampling steps, spm_affreg() does not support multiple sampling steps (contrary to spm_coreg()), so we manage that manually
                sampling_steps = flags_affine.sep;
            else
                % Default sampling steps
                sampling_steps = [10 5];
            end
        else
            % Default flags for affine coregistration
            flags_affine = struct('regtype', 'mni');  % BE CAREFUL: spm_affreg() can introduce reflections, even when using 'rigid' mode! Then decomposition is necessary to ensure left and right sides of the brain are maintained in the correct orientation without mirroring!
            % Default sampling steps
            sampling_steps = [10 5];
        end
        % Load the flags (will change the sep field after each iteration)
        flags = flags_affine;

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
        M = template_vol.mat \ M;  % from spm_affreg() doc, then do M * input_vol.mat to apply the transform (like any column-major transform matrix). Beware, test det(M) to check if the rotation is proper or not (ie, if there is a reflection), it can appear after this step rather than as the result of spm_affreg().

        % Post-processing: constraint the transform to be semi-rigid-body (ie, no shearing, only isotropic scaling allowed, etc)
        % Tuto for debugging: simple calculations (like imcalc but more freedom) in SPM: http://imaging.mrc-cbu.cam.ac.uk/imaging/BasicProgrammingWithSpm
        % Can be used to tag a part of the brain to better highlight what transforms are done (and check there is no reflection):
        % V1=spm_vol('t1.nii');
        % Y1=spm_read_vols(V1);
        % Y1(100:140, 160:200, 91:140) = max(max(max(Y1)));
        % V3 = V1;
        % V3.fname='t12.nii';
        % spm_write_vol(V3,Y1);
        if ~isempty(affdecomposition)
            % Factor out the translation, by simply removing 4th column (4th column being the translation, so it's factored out, as hinted by https://math.stackexchange.com/questions/1120209/decomposition-of-4x4-or-larger-affine-transformation-matrix-to-individual-variab)
            % Note that SPM applies 3D transform matriecs in column-major form (translations are on the right side of the transform matrix M - note that SPM uses the column-major form usually, whereas the MATLAB doc uses the row-major form, for more infos about the difference, see: https://www.youtube.com/watch?v=UvevHXITVG4 )
            % About 3D affine transforms, see also: https://www.youtube.com/watch?v=RqZH-7hlI48
            Mnotrans = M(1:3,1:3);  % remove the translation vector (4th column) from the matrix, that's an easy way to factor out the translation vector (so we only have the other transforms to decompose)
            T = M(:, 4);  % the 3D translation vector is the last column of the M matrix because we use the column-major form
            % Sanity check
            assert(size(Mnotrans, 1) == size(Mnotrans, 2));  % only allow square matrices, not sure if it works on non-square ones
            assert(single(det(Mnotrans)) ~= 0.0);  % determinant is non-zero

            if strcmp(affdecomposition, 'qr')
                % In QR decomposition mode, we will decompose the 3D affine transform matrix M into a rotation matrix Q, a scaling matrix K, and a shearing matrix S.
                % Note that QR decomposition does not prevent improper rotation, if the original transform requires it, there will be a reflection
                fprintf('Using QR decomposition...\n')
                % QR Decomposition to extract the rotation matrix Q and the scaling+shearing matrix R
                % See https://math.stackexchange.com/questions/1120209/decomposition-of-4x4-or-larger-affine-transformation-matrix-to-individual-variab/2353782#2353782
                [Q, R] = qr(Mnotrans);  % pure MATLAB implementations can be found at: http://web.mit.edu/18.06/www/Essays/gramschmidtmat.pdf and https://blogs.mathworks.com/cleve/2016/07/25/compare-gram-schmidt-and-householder-orthogonalization-algorithms/ and https://blogs.mathworks.com/cleve/2016/10/03/householder-reflections-and-the-qr-decomposition/#742ad311-5d1b-4b5d-8768-11799d40f72d
                %[Q, R] = lu(Mnotrans);  % alternative
                if (single(det(Q)) == 1.0)
                    % If Q is a proper rotation matrix (ie, no reflections, as shown by determinant being 1), then all is fine, we can even do a small unit test
                    debug_from_proper_rotation_to_improper_and_back(Q, R);
                elseif (single(det(Q)) == -1.0)
                    % If Q is an improper rotation matrix (ie, with an odd number of reflections), then we need to fix it to a proper rotation matrix and put the reflections in the scaling matrix
                    % We can do so by calculating another QR decomposition on Q, which should separate the reflections out of the proper rotation matrix
                    % Note that this is a necessary step as some QR decomposition algorithms will generate an improper rotation matrix
                    % TODO: would be nice to use our own QR decomposition algo in pure MATLAB with always enforced proper rotation matrix to ensure we never run into this issue
                    %error('QR decomposition led to Q being an improper rotation matrix, cannot continue.');
                    [Q, R] = improper_rotation_to_proper(Q, R);
                end
                % Fix reflections in R (they are stored on the main diagonal as negative scaling factors)
                R = remove_reflections(R);
                % Decompose R into the scaling matrix K and shearing matrix S
                % Compute the scaling matrix K and its inverse, by simply creating K from the main diagonal of R (the rest of the matrix being zero)
                K = eye(size(Mnotrans));
                K(1:size(R,1)+1:end) = diag(R);  % anisotropic rescaling (necessary to compute the accurate inverse and then the shearing matrix S)
                if isorescale
                    % Compute the isotropic scaling matrix K if required
                    Kiso = eye(size(Mnotrans));
                    Kiso(1:size(R,1)+1:end) = mean(abs(diag(R))) .* sign(diag(R));  % compute the average rescaling factor, so we rescale but isotropically (ie, the same rescaling factor in all directions). We first calculate the mean of absolute values, and then we restore the sign. The sign is important because the decomposition may place in the rescaling matrix K negative factors that are necessary to apply reflections.
                    %Kiso(1:size(R,1)+1:end) = mean(abs(diag(R)));
                end
                Kinv = eye(size(Mnotrans));
                Kinv(1:size(R,1)+1:end) = 1 ./ diag(R);  % the inv(K) == the reciprocals (inverse) of the diagonal values in R
                % Compute the shearing matrix S, by removing K from R (so the main diagonal will be zero in S)
                S = Kinv*R;
                % Sanity checks
                % See also what each 3D affine transform look like: https://www.tutorialspoint.com/computer_graphics/3d_transformation.htm and https://stackoverflow.com/questions/5107134/find-the-rotation-and-skew-of-a-matrix-transformation?rq=1
                % Check if there is any reflection, we don't want that
                if (single(det(Q)) == -1.0) || (any(sign(diag(K)) == -1) && are_reflections_imbalanced(K))
                    error('QR decomposition led to an improper rotational matrix Q (rotation+reflection) or scaling matrix K (odd number of negative factors), cannot continue.');
                end
                % Q is a proper rotation matrix if its determinant equals to 1
                assert(single(det(Q)) == 1.0);  % a proper rotation matrix will have a determinant of 1. This means that it's a "pure" rotational matrix, with no other transform (such as reflection/mirroring). See: https://www.tau.ac.il/~tsirel/dump/Static/knowino.org/wiki/Euler%27s_theorem_(rotation).html#Matrix_proof
                % K is a scaling matrix if all non-diagonal values are 0
                assert(all((triu(K, 1)+tril(K,-1)) == 0, [1 2]));
                % S is a shearing matrix if:
                assert(single(det(S)) == 1.0);  % its determinant equals to 1
                assert(all(single(diag(S)) == 1));  % all its main diagonal values equal to 1
                assert(all(single(eig(S)) == 1));  % all eigenvalues equal to 1
                assert((rank(S) == trace(S)) && (trace(S) == size(Mnotrans, 1)));  % its trace equals its rank equals the NxN size of the 3D affine transform matrix
                % Final check of the decomposition, we should find the original 3D affine transform matrix (without the translation of course)
                % DEPRECATED: since we correct reflections in R, which propagates into K, we of course won't necessarily find the original matrix
                %assert(all(single(Mnotrans) == single(Q*K*S), [1 2]));  % convert from doubles to float to clean up the rounding errors during the decomposition

                % Construct the final semi-rigid-body (= rigid-body + isotropic scale, or with anisotropic scale but no shearing) transform
                % Remember that M = (Q*K*S), and with the appending of T in the 4th column
                M2 = Q;
                if ~isempty(isorescale)  % if isorescale == [], then we don't apply the rescaling at all (for debug purposes)
                    if isorescale
                        M2 = M2 * Kiso;
                        disp(Kiso);
                    else
                        M2 = M2 * K;
                    end %endif
                end %endif
                if ~noshearing
                    M2 = M2 * S;
                end %endif
                % Reconstruct a 4x4 affine transform matrix (so that we can put back the translation column vector)
                M3 = eye(size(M));  % the basis needs to be an identity matrix, see: https://www.youtube.com/watch?v=UvevHXITVG4
                M3(1:3,1:3) = M2;  % apply other transforms apart from translation (the precise transforms being defined by user variables)
                M3(:,4) = T;  % Restore the translation
                disp(M);
                disp(M2);
                disp(M3);
                M = M3;
            elseif strcmp(affdecomposition, 'svd')
                % In SVD decomposition mode, we will decompose the affine without translation into 2 rotation matrices + 1 anisotropic rescaling in the middle (with the combination of the rescaling and subsequent rotation to a shearing): https://en.wikipedia.org/wiki/Singular_value_decomposition#/media/File:Singular-Value-Decomposition.svg
                % Pros and cons of using a SVD: unique decomposition (main diagonal is positive), no shearing (containing in the scaling matrix S * last rotational matrix), con is that shearing cannot be separated and hence although it can be removed, the transform will be less faithful to the original non constrained affine transform (but it's ok if regtype = rigid, but not if using mni for example)
                % Note that SVD does NOT prevent reflections, as it can generate improper rotations if required to reproduce the original transform: https://math.stackexchange.com/questions/2331207/relation-between-svd-and-affine-transformations-2d?rq=1
                fprintf('Using SVD decomposition...\n')
                % Calculate the SVD decomposition
                [u, s, v] = svd(Mnotrans);
                % Fix improper rotations in rotation matrices v' or u, by transferring the reflections to the scaling matrix s
                if single(det(u)) == -1
                    [u, s] = improper_rotation_to_proper(u, s);
                end %endif
                if single(det(v')) == -1
                    % Since reflections are non-commutative (because of the order in matrix multiplication), then we can't just transfer the reflections in v' to s, we need to transfer to a new scaling matrix s2
                    % See for more details this excellent document on rotations and reflections in 3D: http://scipp.ucsc.edu/~haber/ph251/rotreflect_17.pdf
                    [v2, s2] = qr(v');
                    [v2, s2] = improper_rotation_to_proper(v2, s2);
                    v = v2';
                else
                    % if v' is already a proper rotation matrix, we simply create a dummy scaling matrix s2 as an identity matrix
                    s2 = eye(size(v'));
                end %endif
                % Remove/balance out the reflections in s
                s = remove_reflections(s);
                s2 = remove_reflections(s2);
                % Isotropic rescaling / no shearing?
                if isorescale || noshearing
                    % Compute the mean of the rescaling factor to do only isotropic rescaling (and not anisotropic)
                    % In SVD mode, the shearing is done in the scaling matrix, so we can't separate both transforms, if we remove one we remove both
                    s = diag(mean(abs(diag(s))) * sign(diag(s)));  % don't forget to restore the sign to balance out the reflections
                    s2 = diag(mean(abs(diag(s2))) * sign(diag(s2)));
                end
                % Sanity check: check if there is any remaining imbalanced reflection/mirroring
                if single(det(u)) == -1 || single(det(v')) == -1 || (any(sign(diag(s)) == -1) && are_reflections_imbalanced(s))
                    error('Reflections detected in the decomposed matrices, cannot continue!');
                end
                % Combine back all transforms
                M2 = u*s*v'*s2;  % apply other transforms apart from translation (the precise transforms being defined by user variables)
                % Reconstruct a 4x4 affine transform matrix (so that we can put back the translation column vector)
                M3 = eye(size(M));  % the basis needs to be an identity matrix, see: https://www.youtube.com/watch?v=UvevHXITVG4
                M3(1:3,1:3) = M2;  % apply other transforms apart from translation (the precise transforms being defined by user variables)
                M3(:,4) = T;  % Restore the translation
                disp(M);
                disp(M2);
                disp(M3);
                M = M3;
            elseif strcmp(affdecomposition, 'imatrix')
                % Decompose and constraint reflections/shearing using spm_imatrix (fixed)
                fprintf('Using spm_imatrix decomposition...\n')

                % Extract the transform parameters as a vector from the 3D affine transform matrix
                iM = spm_imatrix(M);
                assert(all(single(spm_matrix(iM)) == single(M), [1 2]));  % sanity check: check that we can reconstruct the original transform matrix
                % Multiple assign each parameter into its own variable (see help spm_matrix: T for translation, R rotation, Z zoom/scale, S shearing)
                iM2 = mat2cell(iM, 1, ones(1,4)*3);
                [T, R, Z, S] = iM2{:};
                % Fix spm_imatrix inappropriate calculation of the reflections (only the first dimension is ever set to be a reflection, whether or not the reflection happens on this plane or another)
                [~, B] = qr(M(1:3,1:3));
                Z = sign(diag(B))' .* abs(Z);
                % Nullify shearing if option is enabled
                if noshearing
                    S = zeros(1, 3);
                end
                % If there are any reflections in the scaling matrix Z (negative factors), then apply a 180° rotation (pi) instead in the same axis (this is not equivalent but it's better than nothing - the best is to afterward do another coregistration, based on joint histogram for example, to clean up the residual coregistration errors)
                %R = R + (sign(Z) == -1)*pi;  % TODO: buggy
                % Turn negative scaling factors to positive (ie, this removes reflections)
                %Z = Z .* sign(Z);
                Z = diag(remove_reflections(diag(Z)))';
                if isorescale
                    % If isotropic scaling is required by user, compute the mean rescaling and apply the same scaling factor in all dimensions
                    Z = sign(Z) .* mean(abs(Z));
                end
                % Convert the vector back to a 3D affine transform matrix
                % Note that the order of the operations can be changed, see help spm_matrix for more information about default order and how to change
                M = spm_matrix([T R Z S], 'T*R*Z*S');
            end %endif
            % TODO: recenter the origin to middle of the brain?
        end %endif
        % Sanity check: check if there is no reflections nor shearing (or if there are reflections, there is an even number so they cancel out into rotations)
        [T, R, Z, S] = get_imat(M);
        if strcmp(affdecomposition, 'none')
            fprintf('Warning: affdecomposition is set to none, no sanity check will be done to ensure that the object structure is maintained (ie, no reflection nor shearing - meaning that the brain can see its left and right sides be flipped).\n');
        else
            assert(~any(sign(Z) == -1) || ~isodd(sum(sign(Z) == -1)));
        end
        if noshearing
            assert(sum(S) < 0.000001);
        end
        % Apply the reorientation on the input image header and save back the new coordinations in the original input file
        % Some inspiration by BSD-2 licensed: https://github.com/jimmyshen007/NeMo/blob/c7cc775e7fc84f84255b0d3fa474b7f955e817f5/mymfiles/eve_tools/Structural_Connectivity/approx_coreg2MNI.m - licensed under BSD-2-Clause (compatible with MIT license anyway)
        %spm_get_space(inputpath, M*input_vol.mat);
        spm_get_space(inputpath, M*input_vol.mat);
    else
        % Joint Histogram coregistration
        % Keep in mind only rotations can be done, which excludes reflections (that's good), but also rescaling, including isotropic rescaling (that's bad), so we need another step (such as affine reorientation) before to be able to rescale as necessary, and also translate (although this coregistration is invariant to translations, if you try to coregister manually after or do between-subjects coregistration this will be an issue)
        % Note also that spm_coreg() will try to set the input volume's origin to match the template's origin (hence on AC-PC)
        fprintf('Joint Histogram (mutual information) reorientation\n');
        flags = struct('sep', [4, 2]);
        flags.cost_fun = 'ecc';  % ncc works remarkably well, when it works, else it fails very badly... Also ncc should only be used for within-modality coregistration (TODO: test if for reorientation it works well, even on very damaged/artefacted brains?)
        flags.tol = [0.02, 0.02, 0.02, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001];  % VERY important to get good results, these are defaults from the GUI
        % spm_coreg can only find a rigid-body rotation, but pre-translation is not necessary for spm_coreg to work because it's working directly on joint histograms, so empty space is meaningless. Also there is no need to do a QR decomposition here since it's rigid-body.
        % however if later we want to coregister another modality (eg, functional) over the structural, we need to fix the translation too. In this case, spm_affreg is necessary for translation and (isotropic) scaling first on the structural.
        x = spm_coreg(template_vol, smoothed_vol2, flags);
        % Convert the transform parameters into a 3D affine transform matrix (but it's only rigid-body)
        M = spm_matrix(x);
        % Apply the reorientation and save back the new coordinations
        % also apply an inversion to project into the input volume space (instead of template space)
        spm_get_space(inputpath, pinv(M)*spm_get_space(inputpath));  % input_vol.mat == spm_get_space(inputpath)
        % alternative to : M \ input_vol.mat
    end %endif
    fprintf('Autoreorientation done!\n');
end %endfunction

function [T, R, Z, S] = get_imat(M)
% [T, R, Z, S] = get_imat(M)
% Calculate the parameters of an affine transform using spm_imatrix and returns 4 vectors (translation T, rotation R, scale/zoom Z and shearing S)
% Note: this is NOT like a decomposition, because the parameters are dependent on the order in which they are applied (default: inverse of 'T*R*Z*S' in SPM, see help spm_matrix)

    % Extract the transform parameters as a vector from the 3D affine transform matrix
    iM = spm_imatrix(M);
    % Multiple assign each parameter into its own variable (see help spm_matrix: T for translation, R rotation, Z zoom/scale, S shearing)
    iM2 = mat2cell(iM, 1, ones(1,4)*3);
    [T, R, Z, S] = iM2{:};
end %endfunction

function res = isodd(i)
% res = isodd(i)
% From an integer i, returns if odd (1) or even (0)
    res = logical(mod(i, 2));
end %endfunction

function res = are_reflections_imbalanced(K)
% res = are_reflections_imbalanced(K)
% From a scaling matrix K, returns if the reflections are balanced (even number of negative scaling factors on the main diagonal) or not (odd number)
    res = isodd(sum(sign(diag(K)) == -1));
end %endfunction

function R = remove_reflections(R)
% Fix reflections in scale matrix R, generated from the QR decomposition or manually constructed with spm_imatrix (they are stored on the main diagonal as negative scaling factors)
% This is based on Clifford's geometric algebra, which shows that 2 successive reflections in two planes equal one rotation (as it equals to a point-wise reflection): Baylis, William. (2004). Applications of Clifford Algebras in Physics. 10.1007/978-0-8176-8190-6_4. https://www.researchgate.net/figure/Successive-reflections-in-two-planes-is-equivalent-to-a-rotation-by-twice-the-dihedral_fig2_228772873
    if isodd(sum(sign(diag(R)) == -1))  % if the number of reflections is even, it's ok, every 2 reflections is equal to a rotation: https://en.wikipedia.org/wiki/Rotations_and_reflections_in_two_dimensions
        % Else, there's an odd number of reflections, we need to fix it by adding or removing one reflection to make it an even number
        reflections = (sign(diag(R)) == -1);
        % In the following, we assume we want to be in the ac-pc orientation, so that the x axis should be the sagittal axis
        if sum(reflections) == 3
            % If there are 3 reflections (ie, all scaling dimensions are negative and hence reflections), we simply switch off the reflection in the x axis
            R(1,1) = abs(R(1,1));
        elseif sum(reflections) == 1
            % If there is only 1 reflection, ...
            if reflections(1)
                % If the reflection is on the x axis, we simply disable it (as it mostly reverses left-right)
                R(1,1) = abs(R(1,1));
            else
                % If the reflection is on another axis (y or z), then we add a reflection in the x axis to restore the left-right orientation
                R(1,1) = -1.0 .* R(1,1);
            end
        end
    end
end %endfunction

function [Q2, R2] = proper_rotation_to_improper(Q, R)
    D = diag(sign(diag(R)));
    Q2 = Q*D;
    R2 = D*R;
end

function [Q2, R2] = improper_rotation_to_proper(Q, R)
    [~, D] = qr(Q);
    D = diag(sign(diag(D)));
    Q2 = Q*D;
    R2 = D*R;
end

function debug_from_proper_rotation_to_improper_and_back(Q, R)
% Showcases how a proper rotation can be converted to an improper rotation matrix (ie, rotation+odd number of reflections) and back into a proper rotation matrix
    if (single(det(Q)) ~= 1.0)
        error('Cannot run the test since Q is not a proper rotation matrix\n');
    end

    % Get the scaling matrix from R's main diagonal, and then extract the signs of the scaling factors (negative scaling factors are reflections, and only them will modify Q and R)
    % This enforces a positive main diagonal in R, but at the expense of making an improper rotational matrix Q (ie, with reflection), see for example: https://www.mathworks.com/matlabcentral/answers/83798-sign-differences-in-qr-decomposition#answer_93385
    D = diag(sign(diag(R)));
    Q2 = Q*D;
    R2 = D*R;

    % Use spm_imatrix to get the reflections
    % DEPRECATED: can't use spm_imatrix, there's a cheat: if Q is detected as an improper matrix, P(7) (x axis scaling) is set as negative, that's a cheap trick, so we don't know which dimension in fact really has a reflection: if det(R)<0, P(7)=-P(7); end % Fix for -ve determinants
    %Q2full = eye(4,4)
    %Q2full(1:3,1:3) = Q2;
    %iQ2 = spm_imatrixfix(Q2full);
    %Z = iQ2(7:9);
    %D2 = diag(sign(Z));
    %Q3 = Q2*D2;
    %R3 = D2*R2;

    % Correct way to convert an improper rotation matrix back to a proper rotation matrix and transfer the reflections onto the scaling matrix: the R matrix in the QR decomposition of Q will be in fact D. We can then multiply the improper matrix with D to get back a proper rotation matrix.
    %[~, B] = qr(Q * diag([1 1 -1])) % debug: B should have in its main diag the sign provided in the diag in the input
    [~, D2] = qr(Q2);
    D2 = diag(sign(diag(D2)));
    Q3 = Q2*D2;
    R3 = D2*R2;
    % Since we started with a proper rotation matrix, we should get back the same after we correct the improper rotation matrix
    assert(all(Q3 == Q, [1 2]));
    assert(all(R3 == R, [1 2]));
end

function wcentroid1 = get_centroid(svol, binarize, invert, debug)
% wcentroid1 = get_centroid(svol, binarize, invert)
% From an spm volume (use spm_vol()), returns the center of mass (if non binarized, else the centroid/geometric mean)
% Note that if the whole head is acquired (including the neck), then if non binarized, the centroid will have a tendency to lean towards the neck, which we don't necessarily want
% Invert allows to work on images where the brain tissue is in black and the background white

    % Optional arguments
    if ~exist('binarize', 'var')
        binarize = true;
    end
    if ~exist('invert', 'var')
        invert = false;
    end
    if ~exist('debug', 'var')
        debug = false;
    end

    % Access the data in the body of the nifti volume
    input_data = spm_read_vols(svol);

    % Invert image if necessary beforehand (because we expect for the calculations to work that the background to be 0 and only the forefront object to have a higher intensity)
    if invert
        % Calculate the complement of the image (ie, invert the image)
        input_data = input_data + min(input_data, [], 'all');  % first we shift any negative value to be null
        input_data = max(input_data, [], 'all') - input_data;  % then we invert
    end

    % 0- threshold data first, so that the background is not used to weight
    bak_thresh = mean(input_data, 'all'); %median(input_data, 'all');  % use mean instead of median, because we don't need to keep the whole brain, we just want to kick out the background and keep what is most likely the brain
    input_data(input_data <= bak_thresh) = 0;
    if binarize
        input_data(input_data > bak_thresh) = 1;
    end %endif
    % 1- extract the coordinates of each voxel and its associated value, each in a different vector (we can combine into a long-form table later)
    idx1 = find(input_data);  % extract linear/vector style indices
    v1 = input_data(idx1);  % extract values in a vector
    [x1, y1, z1] = ind2sub(size(input_data), find(input_data));  % convert linear indices to 3D indices/coordinates
    % 2- calculate the weighted sum = center of mass in each dimension
    wcentroid1 = (v1'*[x1 y1 z1])./sum(v1);  % weighted centroid
    %wcentroid1 = sum([x1 y1 z1], 1) ./ size([x1 y1 z1], 1);  % non weighted (binarized) centroid. If image is binarized before, the weighted centroid and non-weighted centroid should be the same.

    % Debug test: paint a square around the centroid and save the new nifti image just to see where it is
    if debug
        disp(wcentroid1);
        rwcentroid1 = round(wcentroid1);
        input_data(rwcentroid1(1)-20:rwcentroid1(1)+20, rwcentroid1(2)-20:rwcentroid1(2)+20, rwcentroid1(3)-20:rwcentroid1(3)+20) = 0; % max(input_data, [], 'all', 'omitnan');
        V3 = svol;
        V3.fname = [svol.fname(1:end-4) '-centroid.nii'];
        spm_write_vol(V3, input_data);
    end %endif
end  %endfunction

%%% Additional resources
% infinite dimentional qr decomposition: https://link.springer.com/article/10.1007/s00211-019-01047-5
% history of QR decomposition: https://www.math.unipd.it/~alvise/AN_2016/LETTURE/QR_50_years_later.pdf
% https://www.math.uci.edu/~chenlong/RNLA/LSQRSVD.pdf

% molecules transforms etc: https://www.nagwa.com/en/videos/916195283074/
% improper rotation leads to a different molecule: http://chemistryschool.net/symmetry-operation-proper-and-improper-rotations/
% https://www.youtube.com/watch?v=ygB9QbbaM-I

% reflections do not preserve orientation: https://4.files.edl.io/b164/01/17/19/172157-1eb4742a-f536-4b4b-bcc2-0c13410a9fa7.pdf
% rotations of 180° are equivalent to a reflection through the origin? Not for reflection asymmetrical objects!
% 2 reflections = 1 rotation: https://www.euclideanspace.com/maths/geometry/rotations/theory/as2reflections/index.htm

% Manual rotational and scaling realignment? Without reflections?
% BEST RESOURCE: superimpose3D in Python, but no isometric rescaling nor anti-reflection countermeasures (but can implement them + dicom support, would be a pure python implementation of reorientation): https://github.com/jewettaij/superpose3d
% See also this great tutorial: http://nghiaho.com/?page_id=671
% and this great paper: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.173.2196&rep=rep1&type=pdf
% and this SO question: https://scicomp.stackexchange.com/questions/10584/purely-rotational-least-squares-match
% and An Approximate and Efficient Method for Optimal Rotation Alignment of 3D Models
% BEST: and Least-Squares Rigid Motion Using SVD, by Olga Sorkine-Hornung and Michael Rabinovich, 2017 - in particular the entry: "Orientation rectification"
% SPM's normalize function uses the origin as a starting estimate, according to Chris Rorden: https://github.com/rordenlab/spmScripts/blob/master/nii_setOrigin.m - BTW it reuses the coregistration parameters from K.Nemoto https://web.archive.org/web/20180727093129/http://www.nemotos.net/scripts/acpc_coreg.m

% best results: autoreorient('t1.nii', 'mi', [], false, false, 'none', true, 'scanner', true)
% good result but not exactly on AC-PC: best results: autoreorient('t1.nii', 'mi', [], false, false, 'none', true, 'raw', true)
