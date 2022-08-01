function returncode = autoreorient(varargin)
% autoreorient(inputpath, mode, flags_affine, noshearing, isorescale, affdecomposition, precoreg, precoreg_reset_orientation, just_check, debug)
% 
% autoreorient of structural and BOLD and other MRI modalities for SPM12 and MATLAB (tested on v2018b)
%
% This function supports named arguments, use it like this:
% autoreorient('inputpath', 'path/to/file.nii', 'mode', 'jointhistogram', 'debug', true)
% 
% DEVELOPMENT FUNCTION
% This function is kept as a barebone version of the core routines to do autoreorientation using SPM12 functions. This is kept for development purposes (to quickly debug out only the core routines), do not use it for production.
%
% ## Input variables:
% * inputpath: path to the input nifti file to reorient to template.
% * mode: defines how the registration to template will be done: 'affine' or 'jointhistogram'. Note that jointhistogram uses spm_coreg(), it is the recommended mode as it produces better results generally and produces true rigid-body transforms without reflections nor anisotropic scaling nor shearing, and is equivalent to the 'ecc' or 'mi' mode in the antecedant package spm_auto_reorient_coregister. 'affine' uses spm_affreg() and is not recommended as it can produce reflections (inverted hemispheres) and anisotropic scaling, even with the 'rigid-body' regtype flag.
% * flags_affine: custom flags to pass to the SPM12's spm_affreg() function.
% * noshearing: if true, if shearing is detected after reorientation to template, try to nullify the shearing. Default: true.
% * affdecomposition defines the decomposition done to ensure structure is maintained (ie, no reflection nor shearing). Can be: qr (default), svd (deprecated), imatrix or none. Only the qr decomposition ensures the structure is maintained.
% * precoreg_reset_orientation: before coregistration, reset orientation? Value can be: true (null orientation but keep scale/voxel-size), a vector of 3 values (null orientation and reset scale/voxel-size with the provided vector - this is the only way to reset the voxel-size, as the nifti format does NOT store the scanner's original orientation nor voxel-size/dimension, so the user must provide a vector of voxel-size from eg the MRI machine printout of sequences parameters) or 'mat0' (previous orientation matrix if available) or false (to disable)
% * just_check: if true, the software will only check if the input image has an issue with its structure (ie, reflection or shearing), and will return 0 if no issue, or >0 if there is an issue (1: shearing, 2: anisotropic scaling, 4: reflection(s), greater values than 4 represent the sum of multiple issues found, eg, 5 = shearing + reflections).
%
% License: MIT License, except otherwise noted in comments around the code the other license pertains to.
% Copyright (C) 2020-2022 Stephen Karl Larroque, Coma Science Group & GIGA-Consciousness, University Hospital of Liege, Belgium
% v0.6.18
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, includin without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    % == Load auxiliary functions
    % Addpath of the auxiliary library
    if ~exist('auto_acpc_reorient_aux.m','file')
        %restoredefaultpath;
        addpath(genpath(strcat(cd(fileparts(mfilename('fullpath'))),'/auto_acpc_reorient_aux/')));
    end

    % Importing auxiliary functions
    % this way of importing works with both MatLab and Octave
    spmaux = auto_acpc_reorient_spmaux;  % Auxiliary functions to interface with SPM
    aux = auto_acpc_reorient_aux;  % Our own auxiliary functions

    % == Arguments processing
    % List of possible arguments and their default values
    arguments_defaults = struct( ...
        ... % Mandatory
        'inputpath', 0, ...
        'mode', 0, ...
        ... % Optional
        'flags_affine', [], ...
        'noshearing', true, ...
        'isorescale', true, ...
        'affdecomposition', 'qr', ...
        'precoreg', true, ...
        'precoreg_reset_orientation', false, ...
        'just_check', false, ...
        ...
        ... % Debug stuffs
        'debug', false);

    % Process the arguments
    arguments = aux.getnargs(varargin, arguments_defaults, true);

    % Load variables into local namespace (called workspace in MatLab)
    aux.varspull(arguments);
    
    % Sanity Checks on input variables
    if inputpath == 0 or mode == 0
        error('Missing arguments: inputpath and mode are mandatory!');
    end

    % == Load brain templates
    % Load input and template images
    input_vol = spm_vol(strtrim(inputpath));
    template_vol = spm_vol(strtrim('T1_template_CAT12_rm_withskull.nii'));

    % == Sanity checks on input image
    iM = spm_imatrix(input_vol.mat);
    structure_errors_counter = 0;
    %iM_orig = spm_imatrix(input_vol.private.mat0);  % NOT reliable: svol.private.mat0 is simply the previously saved orientation matrix, not necessarily the original one
    if (sum(abs(iM(10:12))) > 1E-5)
        fprintf('Warning: Shearing detected in input image, this can sometimes be fixed by setting precoreg_reset_orientation = "scanner" or noshearing = true.\n');
        structure_errors_counter = structure_errors_counter + 1;
    end
    if (sum(abs(iM(7:9))) > 1E-5)  % alternative: calculate the QR decomposition
        fprintf('Warning: Anisotropic scaling detected in input image, this can sometimes be fixed by setting precoreg_reset_orientation = "scanner".\n');
        structure_errors_counter = structure_errors_counter + 2;
    end
    [Q, R] = qr(input_vol.mat(1:3,1:3));
    Z = sign(diag(R)) .* abs(diag(R));
    if ((single(det(Q)) == -1.0) || (any(sign(Z) == -1, 'all') && isodd(sum(sign(Z) == -1))))
        fprintf('Warning: Reflections detected in input image (i.e., the orientation is not preserved from the original orientation, so the left side of the brain _may_ be inversed with the right side!), this can be fixed by setting precoreg_reset_orientation = "scanner".\n');  % note that some scanner (such as Siemens) may set a reflection in the initial orientation matrix they set, without impacting the left-right orientation (ie, the reflection is in another dimension)
        structure_errors_counter = structure_errors_counter + 4;
    end
    if just_check == true
        % just checking structural errors, we return a unique error code to specify all errors found
        returncode = structure_errors_counter;
        return;
    end

    % == Reset orientation?
    if precoreg_reset_orientation ~= false
        % Important note: there is NO way to recover the original dimensions and orientation, as the initial orientation matrix or quaternions are not saved but replaced by newer values. So if this was tampered, the user needs to manually specify the original dimensions (eg, by using values such as the voxel size in the printout of the sequence parameters from the MRI machine).
        % For more infos, see:
        % Tuto on nifti orientations systems: http://www.grahamwideman.com/gw/brain/orientation/orientterms.htm
        % http://www.grahamwideman.com/gw/brain/tools/gworc/index.htm
        % Tuto on original rationale behind orientation and specification in nifti format: https://nifti.nimh.nih.gov/dfwg/presentations/nifti-1-rationale
        % Test with MRIcron by unchecking "Reorient image before loading" option in the preferences, this will load the raw image with an identity matrix. If you have an anisotropic image (eg, bold, dti), without reorientation you will see the image will be anisotropic, showing there is no way to infer what was the original voxel-size if the orientation matrix is not used.
        % info = niftiinfo('t1.nii'); info.Transform.T' shows the orientation matrix (alternative to SPM).
        if precoreg_reset_orientation == true
            fprintf('Reset orientation but keep scale/voxel-size\n');
            % raw original voxel space with no orientation
            % get voxel size
            % NOT reliable, this may have been manipulated. Unfortunately, we then have no way to get back the original dimensions...
            vox = sqrt(sum(input_vol.mat(1:3, 1:3).^2));
            % build resetted matrix basis, based simply on the voxel scaling, but reset the whole orientation (ie, everything else is 0 apart from the main diagonal)
            M = diag([vox 1]);
            % calculate centroid, set the origin on it and save back the orientation matrix (including new origin on centroid) in the volume and the file
            M = set_origin_to_centroid_and_save(input_vol, M, debug);
            % update input volume in-memory
            input_vol.mat = M;
        elseif isvector(precoreg_reset_orientation) && (numel(precoreg_reset_orientation) == 3)
            fprintf('Reset orientation and scale/voxel-size to the specified vector:\n');
            disp(precoreg_reset_orientation);
            % Reset with a blank orientation matrix but with the scale factors as provided by user
            M = diag([precoreg_reset_orientation(:); 1]);
            % calculate centroid, set the origin on it and save back the orientation matrix (including new origin on centroid) in the volume and the file
            M = set_origin_to_centroid_and_save(input_vol, M, debug);
            % update input volume in-memory
            input_vol.mat = M;
        elseif strcmp(precoreg_reset_orientation, 'mat0')
            fprintf('Reset orientation to previous orientation\n');
            % reuse previously saved orientation matrix
            % simply reload mat0 and overwrite the current orientation matrix
            input_vol.mat = input_vol.private.mat0;
            spm_get_space(inputpath, input_vol.mat);
        end
        fprintf('Reset orientation done!\n');
    end

    % == Smooth input image, helps a lot with the coregistration (which is inherently noisy since there is no perfect match)
    smoothed_vol = spm_vol(inputpath);
    smoothed_vol2 = spm_smoothto8bit(smoothed_vol, 20);
    % Y = spm_read_vols(smoothed_vol2);  % direct access to the data matrix
    % imagesc(Y(:,:,20))  % show a slice
    keyboard

    % == Manual/Pre coregistration (without using SPM)
    if precoreg
        fprintf('Pre-coregistration by translation of centroid (resetting origin), please wait...\n');
        % Find the center of mass
        input_centroid = [get_centroid(input_vol, true, false, debug) 1];  % add a unit factor to be able to add the translation of the world-to-voxel mapping in the nifti headers, see: https://www.youtube.com/watch?v=UvevHXITVG4
        %template_centroid = [get_centroid(template_vol, true, false, debug) 1];
        %input_centroid_voxel = (input_vol.mat * input_centroid')';  % convert from raw/scanner voxel space to current orientation voxel space (ie, this represents the distance in voxels from the current origin to the centroid, in other words this is the centroid position relative to the currently set origin)
        %template_centroid_voxel = (template_vol.mat * template_centroid')';
        % Calculate the euclidian distance = translation transform from the input to the template centroid, we also resize to match the input volume's space (which may be bigger or smaller - if bigger, then we need to reduce the number of voxels we move)
        %i2t_dist = (template_centroid_voxel(1:3) - input_centroid_voxel(1:3)) .* (template_vol.dim ./ input_vol.dim);

        % Apply the translation on the nifti header voxel-to-world transform
        M = set_origin(input_vol.mat, input_centroid, true);
        %M(1:3,4) = M(1:3,4) + i2t_dist(1:3)';  % shift some more to match where the origin is in the template compared to its own centroid. Not sure this step is necessary since anyway we certainly won't end up in the AC-PC, but well why not, it may bridge some more the relative distance between the template origin and input volume origin. DEPRECATED: unreliable, it's preferable to stick to the centroid
        % Note: at this point, we should also set the template's origin on its centroid to leave only the rotation and scale to be estimated afterwards. But since we provide a template with the origin on the AC-PC, it's an even better origin, so we skip this step here. If the template did not have the origin on AC-PC, then setting the origin on its centroid would be a good thing to do, see: http://nghiaho.com/?page_id=671.
        % Save into the original file
        spm_get_space(inputpath, M);
        if debug, spm_get_space([inputpath(1:end-4) '-centroid.nii'], M); end;  % debug line, save the same transform on the image where the centroid is made visible
        % Apply the new mapping in the already loaded images including smoothed (this avoids the need to reload the input volume and resmooth it)
        input_vol.mat = M;
        smoothed_vol.mat = M;
        smoothed_vol2.mat = M;
        % SPM provides 3 ways to save nifti files:
        % * spm_get_space(vol, transform) to only modify the world-to-voxel mapping headers (vol.mat)
        % * spm_write_vol(vol, vol_data) where vol is a spm nifti structure as created by spm_vol(), and vol_data is a 3D matrix of the content of the nifti, and vol a nifti structure as accepted or created by spm using spm_read_vols().
        % * create(vol) where vol is a spm nifti structure as created by nifti(), and it will be saved in vol.fname path. This last option is the most complete, as it will save everything.
        %spm_write_vol(input_vol, input_vol.mat);
        fprintf('Pre-coregistration done!\n');
    end

    % == Affine coregistration (with translation)
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
                Z = sign(diag(B))' .* abs(diag(Z));  % extract the rescaling factor signs from the QR decomposition - TODO: depending on the QR algorithm, some of the reflections may be in Q, so we should also check Q
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
        spm_get_space(inputpath, M*input_vol.mat);
        % spm_affreg does not set the origin, so we do it manually on the centroid/center of mass (but it's strongly advised to call spm_coreg() after spm_affreg instead, as spm_coreg() will also match and project the template's origin onto the input volume's origin
        M = set_origin_to_centroid_and_save(input_vol, M*input_vol.mat, debug);
        % update input volume in-memory
        input_vol.mat = M;
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
        % TODO: to enhance performance, first do a spm_coreg using a binarized input image (ie, a brain mask), using the mean as is done for the centroids calculation? Just to ensure to have a good starting point for orientation.
    end %endif
    fprintf('Autoreorientation done!\n');
    returncode = 0;
    return;
end %endfunction, end of main routine

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
% The method used here is to decompose a rotation matrix to find in which dimensions are the negative scaling factors and balance the number of negative factors to always be even.
%
% Note that there is a theoretically simpler and equivalent alternative, to simply calculate R = V*diag([1 1 1 det(V*U')])*U';
% References for this simpler method:
% * Sorkine-Hornung, O., & Rabinovich, M. (2017). Least-squares rigid motion using svd. Computing, 1(1). https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
% * Eggert, D. W., Lorusso, A., & Fisher, R. B. (1997). Estimating 3-D rigid body transformations: a comparison of four major algorithms. Machine vision and applications, 9(5-6), 272-290. http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.173.2196&rep=rep1&type=pdf
% * http://nghiaho.com/?page_id=671
%
% Either method is giving equivalent results. Note however that fixing reflections will only lead to an approximate orientation without reflection, it won't be exactly the same as with the reflection.

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
% This is invariant from the current orientation matrix M, as we are working directly on the voxels intensities and coordinates
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

function M = set_origin(M, orig, absolutecoord)
% M = set_origin(M, orig, absolutecoord)
% Sets the origin to the coordinates orig and save into svol file.
% if absolutecoord is true (default), the orig coordinate is taken to be in the raw scanner space (before applying the orientation matrix M or svol.mat), hence some additional calculations are done. If the orig coordinates are already relative to the current origin (and hence orientation matrix M), set absolutecoord to false.
% Return M, so that the volume can be updated in-memory using svol.mat = M; or saved back into a file using spm_get_space(svol.fname, M);

    if ~exist('absolutecoord', 'var')
        absolutecoord = true;
    end

    % make it a column vector in any case
    orig = orig(:);

    % add a unit factor to be able to add the translation of the world-to-voxel mapping in the nifti headers, see: https://www.youtube.com/watch?v=UvevHXITVG4
    if numel(orig) == 3
        orig = [orig; 1];
    end

    % convert from raw/scanner voxel space to current orientation voxel space (ie, this represents the distance in voxels from the current origin to the centroid, in other words this is the centroid position relative to the currently set origin - this also includes rescaling of voxels from the main diagonal values, and rotation and shearing)
    if absolutecoord
        orig = M * orig;  % should be a column vector (in and out)
    end

    % Resets the origin onto the centroid, by translating the previous origin to the new one
    % (nifti viewers such as MRIcron or SPM will apply the transform in a vol.premul matrix before the rest)
    M(1:3,4) = M(1:3,4) - orig(1:3);  % set the origin onto the centroid
end  % endfunction

function M = set_origin_to_centroid_and_save(svol, M, debug);
% set_origin_to_centroid_and_save(svol, M)
% calculate centroid, set the origin on it relatively to the provided orientation matrix M (can be svol.mat), and save back the orientation matrix (including new origin on centroid) in the volume and the file
% Returns nothing (saves to file svol.fname directly) - note you should probably also update svol.mat in-memory using the returned M

    if ~exist('debug', 'var')
        debug = false;
    end

    % calculate centroid, will be the origin (this is invariant to the orientation matrix M)
    orig = get_centroid(svol, true, false, debug);
    % set origin to centroid
    M = set_origin(M, orig);
    % save back to the file and update the volume
    spm_get_space(svol.fname, M);
    if debug, spm_get_space([svol.fname(1:end-4) '-centroid.nii'], M); end;
end  %endfunction

function [best, bestfit, bestfitstd, bestpool, bestpoolfit] = mga_imatrix(genestart, popsize, demesize, fitnessfunc, recombinationrate, mutationrate, tournamentcount, multistart, randomstart, maxdim, randomgen, mutatefunc, truerandom, onlyrigid, debug)
% [best, bestfit, bestfitstd] = mga_imatrix(genestart, popsize, demesize, fitnessfunc, recombinationrate, mutationrate, tournamentcount, multistart, randomstart, maxdim, randomgen, mutatefunc, debug)
% From a spm_imatrix() vector, will find a better/locally optimal translation and orientation given the provided fitnessfunc by using the microbial genetic algorithm (MGA)
% genestart is a vector of a single candidate that can be fed to fitnessfunc(genestart). This will both be used as the reference candidate (so we ensure we don't return a worse solution), and to generate new candidates. genestart should be a 12 elements vector as made by spm_imatrix(svol.mat), but don't forget that spm_coreg.optfun and spm_hist2 will apply it on top of input_vol.mat, such that orientation_that_will_be_tested = pinv(spm_matrix(genestart))*input_vol.mat -- in spm_coreg you can see the full call to spm_hist2 being: VF.mat\spm_matrix(x(:)')*VG.mat -- so by default you want genestart to be [0 0 0 0 0 0] (as is done in spm_coreg, which means we start from input_vol.mat as-is without any modification), NOT spm_imatrix(input_vol.mat) which would be a transform that would be applied ON TOP of input_vol.mat
% maxdim is the maximum dimension to generate rotations when combined with randomstart = true. In general, set maxdim = svol.dim
% popsize is the size of the population. This number should be relative to the number of tournamentcount, a good number is 1/10th of tournamentcount, as to allow all candidates to be explored and recombined by the end of the algorithm.
% demesize is the size of the deme, which is the subpopulation size inside the ring that will be used to find a candidate B to compare candidate A. See Inman Harvey's paper. A good value is demesize = 0.1 to 0.3 * popsize.
% fitnessfunc is the fitness function to use, it should accept a single variable x which will be a vector of the same size as genestart. Tip: try to optimize this function to be super fast, as although this genetic algorithm uses caching to try to save some time, a computation expensive fitnessfunc will limit the number of tournamentcount you can do (and hence the likelihood of converging to a satisficing solution). Specify fitnessfunc like this: fitnessfunc = @(x) jointhistogramfitness(spmaux, template_vol, input_vol, x);
% recombination rate is the crossover rate, the probability that one loser's gene will be overwritten by the winner's.
% mutation rate is the rate of randomly assigning a random value to one loser's gene.
% tournamentcount is the number of rounds to find a winner and a loser candidates and update the loser's genes. A high value increases the likelihood of converging to a satisficing solution.
% multistart is the number of independent rounds to relaunch the whole genetic algorithm, with a brand new population. A high value reduces the variance of the final result. This can be used to test the influence of different parameters values (ie, to hyper-optimize). Multistart allows to reduce variance since we restart from multiple starting points and travel different paths, but this does not help with converging to a better result, it's better to increase tournamentcount than multistartcount to improve performances. Increasing multistart is great to test different parameters (eg, recombinationrate or mutationrate) and see in one result whether the new parameter really improves the performance.
% randomstart defines how the initial population is generated: false will generate a whole population as a copy of the input candidate genestart, true will keep genestart as the reference candidate but all other individuals will be randomly generated according to the randomgen function. If randomstart = true, randomgen needs to be a function accepting popsize as 1st argument and length of genestart as 2nd argument, and which generates popsize individuals to add in the pool. For example, if fitnessfunc = @(x)sum(x), randomgen can be @(popsize,numel_genestart)randi(10, popsize, numel_genestart)
% randomgen is a function that defines how the initial population is randomly generated when randomstart == true. randomgen should expect 2 parameters: popsize and numel(genestart).
% mutatefunc defines how the mutation is done (ie, how random values are calculated to be assigned to loser's genes).
% truerandom defines whether the function is deterministic (false) or non-deterministic (true or a rand('seed') integer)
% onlyrigid - if true (default), only rigid-body transforms will be generated (ie, translation and rotation). If false, scale and shearing will also be explored.
%
% Output: best is a imatrix with all the transform parameters as a 12 items vector (use spm_matrix(best) to transform to an orientation matrix), bestfit fitness score of the best candidate, bestfitstd the standard deviation over all multistart rounds (allows to assess performance variance by evaluating the variability of the current set of parameters)
%
% Example usage: [best, bestfit, bestfitstd] = autogen([1 2 3 4 5 6 7 8], 100, 10, @(x)sum(x), 0.7, 0.25, 1000, 1000, true, @(popsize,numel_genestart)randi(10, popsize, numel_genestart))
% another example with a smaller number of tournament rounds and hence smaller population, with similar performances: [best, bestfit, bestfitstd] = autogen([1 2 3 4 5 6 7 8], 10, 3, @(x)sum(x), 0.7, 0.25, 100, 1000, true, @(popsize,numel_genestart)randi(10, popsize, numel_genestart))
%
% Tips: The 4 most important parameters to hyper-optimize for accuracy performance are: popsize, demesize, recombinationrate and mutationrate. Assess using (maximizing) bestfit and (minimize) bestfitstd with a high number of (>1000) multistart rounds, which would mean that the genetic algorithm with your parameters and fitnessfunc can reach a high score while having low variance. Also, increasing tournamentcount of course allows the algorithm to converge to a better solution, so optimize your fitnessfunc to be super fast (can use caching for example).
%
% Reference: algorithm from: The Microbial Genetic Algorithm, Inman Harvey, 1996
% There exists a single-line version, see "Evolvable Neuronal Paths: A Novel Basis for Information and Search in the Brain", 2011, Fernando C, Vasas V, Szathmáry E, Husbands P. PLoS ONE 6(8): e23534. doi: 10.1371/journal.pone.0023534
%
% For a more general implementation (not tailored specifically for the reorientation task), take a look at aux/autogen.m and aux/autogen_mini.m

    % Default values
    if ~exist('multistart', 'var') || isempty(multistart)
        multistart = 1;
    end
    if ~exist('debug', 'var')
        debug = false;
    end
    if ~exist('randomstart', 'var') || isempty(randomstart)
        randomstart = false;
    end
    if ~exist('mutatefunc', 'var') || isempty(mutatefunc)
        mutatefunc = false;
    end
    if ~exist('randomgen', 'var') || isempty(randomgen)
        randomgen = false;
    end
    if ~exist('truerandom', 'var') || isempty(randomgen)
        truerandom = false;
    end
    if ~exist('onlyrigid', 'var') || isempty(onlyrigid)
        onlyrigid = true;
    end

    % Sanity checks
    if ~isvector(genestart)
        % ensure genestart is a vector
        genestart = spm_imatrix(genestart);
    else
        % else if it's already a vector, convert to a matrix and back to a vector to ensure that all parameters are provided (so that we can simply provide [0] and the rest will be filled)
        genestart = spm_imatrix(spm_matrix(genestart));
    end

    % Init multistart gene pool
    bestpool = repmat(genestart, [multistart, 1]);
    bestpoolfit = zeros(multistart, 1);
    if onlyrigid
        genescount = 6;  % restrict to rotation and translation (first 6 items in spm_matrix())
    else
        genescount = numel(genestart);  % explore all genes (rotation, translation, scale and shearing)
    end

    % Save current seed for the random number generator, because SPM functions can reset it (although it is bad practice, their goal is to transform non-deterministic functions to deterministic like that... But this may explain why the same analysis can give different results when ran on different computers/OSes, as this changes the underlying random number generators)
    % We will do 3 things at once here:
    % 1- save the original rng so we can restore it at the end
    orig_rng = rng();
    % 2- set the rng so this function is deterministic
    % Note: setting the generator type is discouraged but we want that so that results are reproducible in the future, we don't need secure pseudorandomness. For more infos, see: https://www.mathworks.com/help/matlab/math/updating-your-random-number-generator-syntax.html
    if isstruct(truerandom)
        %rng(truerandom);
        rand('seed', truerandom);
    elseif ~truerandom
        %rng(0,'philox');  % philox works better than twister. It's important to have a good rng, else the exploration of new solutions by random mutations and random start population will be less effective
        % DEPRECATED: rng() is advised by mathworks but it's magnitude of times slower. It's not made to be in a loop, but we need to reset the rng inside the loop since any call to the fitnessfunc may reset the rng (eg, if fitnessfunc calls spm functions). So we are back to using rand('seed'), which is the fastest option: https://stackoverflow.com/questions/39251206/speed-up-random-number-generation-in-matlab
        % See also these threads about this performance issue of rng():
        % https://www.mathworks.com/matlabcentral/answers/128411-rng-is-slow-to-control-random-number-generation-can-i-go-back-to-seed-or-state
        % https://www.mathworks.com/matlabcentral/answers/5320-speed-improvement-of-the-random-generator
        % https://www.mathworks.com/matlabcentral/answers/67667-performance-degradation-of-random-number-generation-in-matlab-r2013a
        rand('seed', 0);
    % else: we don't reset the rng, we leave as default and hence we get true pseudorandomness (non-deterministic output since we continue with the previous rng state)
    end
    % 3- save the current rng state, so that we can restore it after each spm call
    cur_rng = rand('seed');  %cur_rng = rng();

    % For each multistart round (where we restart anew from scratch - this is different from tournaments where each tournament round reuses the same population, here we change the whole population)
    for m=1:multistart
        fprintf('Multistart: %i/%i\n', m, multistart);
        % Init population's genes pool vars
        if ~randomstart
            % Initialize the gene pool by simply replicating the provided genestart vector
            genepool = repmat(genestart, [popsize 1]);
            genepoolfit = ones(popsize,1) .* fitnessfunc(genestart);  % cache the fitness scores, this is a trick to speed up calculations if fitnessfunc is computation expensive
            bestfit = NaN;  % it's useless to compute a bestfit from input because anyway we replicated the input so we are guaranteed to only find a better solution
            best = NaN;
        else
            % Else initialize the gene pool by generating a random translation and orientation for each individual, but retain the original scale and shear from the genestart vector
            if randomgen ~= false
                % User-defined function
                genepool = randomgen(popsize, numel(genestart));
            else
                % Default function
                if onlyrigid
                    genepool = [round(rand(popsize, 3).*maxdim) - maxdim, rand(popsize,3).*(2*pi) - pi, repmat(genestart(7:end), [popsize 1])];  % random rotation and translation
                    %genepool = [repmat(genestart(1:3), [popsize 1]), rand(popsize,3).*(2*pi) - pi, repmat(genestart(7:end), [popsize 1])];  % random rotation only
                else
                    genepool = [round(rand(popsize, 3).*maxdim) - maxdim, rand(popsize,3).*(2*pi) - pi, 0.5 + rand(popsize, 6)*1.5];  % random rotation, translation, scale and shearing
                end
            end  % endif
            genepoolfit = NaN(popsize,1);  % cache the fitness scores, this is a trick to speed up calculations if fitnessfunc is computation expensive
            % but keep the first candidate as the initial one, so we ensure that any candidate we choose is not worse than the input
            genepool(1, :) = genestart;
            best = 1;
            bestfit = fitnessfunc(genestart);
            %bestfit = NaN;
            %best = NaN;
        end  % endif
        % Launch the tournament
        for t=1:tournamentcount
            if mod(t, 10) == 0, fprintf('Tournament: %i/%i\n', t, tournamentcount); end;
            % Restore previous state of rng (to avoid SPM meddling with rng - each call to fitnessfunc calls spm and hence meddles with the rng)
            rand('seed', cur_rng);  % rng(cur_rng); -- much slower...
            % Randomly select one individual A
            A = randi(popsize);  % alternative: A = ceil(rand()*popsize);
            % Randomly select another individual B in the deme just after A
            B = mod((A+1+randi(demesize)), popsize)+1;  % alternative: B = mod((A+1+ceil(rand()*demesize)), popsize)+1;
            % Save current rng state (to restore at next tournament)
            cur_rng = rand('seed');  %cur_rng = rng();
            if debug, disp([A B]); end;
            % Compute fitness cost for each candidate
            if ~isnan(genepoolfit(A))
                % If there is a cache for this candidate, use it
                Afit = genepoolfit(A);
            else
                % Else compute the fitness cost
                Afit = fitnessfunc(genepool(A,:));  % memoize fitness to optimize, so that we can reuse directly at the end
            end
            if ~isnan(genepoolfit(B))
                Bfit = genepoolfit(B);
            else
                Bfit = fitnessfunc(genepool(B,:));
            end
            % Find the winner
            if debug, disp(genepool(A,:)); disp(genepool(B,:)); disp(Afit); disp(Bfit); end;
            if (Afit > Bfit)
                winner = A;
                loser = B;
                winnerfit = Afit;
            else
                winner = B;
                loser = A;
                winnerfit = Bfit;
            end  % endif
            % Update fitness cost cache ...
            genepoolfit(loser) = NaN;  % ... by deleting (NaN) the loser's fitness cost, so next time it will be recomputed...
            genepoolfit(winner) = winnerfit;  % ... and by caching the winner's fitness cost.
            % Compare winner with the best fit (memoization)
            if winnerfit >= bestfit || isnan(bestfit)  % note: it's crucial to use >= (and not >) because it can happen that both candidates are equal (they maxxed out), in this case there will still be a loser who will be mutated, which can hence become suboptimal. In that case, if the loser was the best candidate, we need to pass the best label to the winner (despite them being equal), because the winner will not change.
                bestfit = winnerfit;
                best = winner;
            end  % endif
            % Recombine and mutate for each gene
            rand('seed', cur_rng);  %rng(cur_rng);  % Restore previous state of rng (to avoid SPM meddling with rng)
            for i=1:genescount  % should be i=1:genescount but we limit to the rotation parameters in spm_matrix() (ie, no rescaling nor shearing, we only do rotation and translation = rigid-body transform)
                r = rand();
                if r < (recombinationrate+mutationrate)  % optimization, see slide 20 of: https://fr.slideshare.net/lrq3000/pathway-evolution-algorithm-in-netlogo
                    if r < recombinationrate
                        % Recombine/crossover (ie, take the allele from the winner)
                        genepool(loser, i) = genepool(winner, i);
                    else
                        % Mutate
                        if mutatefunc ~= false
                            % User-defined mutation function
                            genepool(loser, i) = mutatefunc(i);
                        else
                            % Default mutation function
                            if 1 <= i && i <= 3
                                % Translation, constrained to image dimension
                                genepool(loser, i) = round(rand()*(maxdim(i).*2) - maxdim(i));
                            elseif 4 <= i && i <= 6
                                % Rotation, constrained to radians
                                genepool(loser, i) = rand()*2*pi - pi;
                            elseif i >= 7
                                % Scaling and shearing, bounded inside 0.5 to 2.0 here (but in reality it's not bounded, but then the exploration space is way to big)
                                genepool(loser, i) = 0.5 + rand()*1.5;
                            end  % endif
                        end
                    end  % endif
                end  % endif
            end  % endfor genes walking
            cur_rng = rand('seed');  %cur_rng = rng(); % Save current rng state (to restore at next tournament)
        end  % endfor tournament rounds

        % Select the best candidate
        % DEPRECATED: manual comparison by iterating and comparing all candidates. This does not use memoization. If fitnessfunc is expensive, this will take a long time to compute
        %best = 1;
        %bestfit = fitnessfunc(genepool(best, :));
        %for i=2:popsize
        %    newfit = fitnessfunc(genepool(i, :));
        %    if newfit > bestfit
        %        best = i;
        %        bestfit = newfit;
        %    end  % endif
        %end  % endif

        % Save best candidate of this multistart run
        bestpool(m,:) = genepool(best, :);
        bestpoolfit(m) = bestfit;
    end  %endfor multistart rounds

    % End of multistart: select the best candidate over all runs
    [~, bestidx] = max(bestpoolfit);
    best = bestpool(bestidx, :);
    bestfit = bestpoolfit(bestidx);
    bestfitstd = std(bestpoolfit);

    % Restore original rng
    rng(orig_rng);
end  % endfunction

function fitness = jointhistogramfitness_core(template_vol_uint8, input_vol_uint8, M, mi_algo, voxelsize, fwhm_hist)
% fitness = jointhistogramfitness_core(template_vol_uint8, input_vol_uint8, M, mi_algo, voxelsize, fwhm)
% Compute a fitness cost from the mutual information of the joint histogram of two images
% by calling spm_coreg.optfun() (private method which is exposed if varargin > 4)
%
% template_vol_uint8 and input_vol_uint8 must be smoothed volumens and converted to uint8 format (ie, with a .uint8 field) beforehand, using spmaux.smoothvol2uint8().
% M is an orientation matrix or a vector as returned by spm_imatrix() to apply on input_vol and ON TOP of input_vol.mat to match template_vol.mat. Indeed, spm_coreg.optfun and spm_hist2 will apply it on top of input_vol.mat, such that orientation_that_will_be_tested = pinv(spm_matrix(genestart))*input_vol.mat -- in spm_coreg you can see the full call to spm_hist2 being: VF.mat\spm_matrix(x(:)')*VG.mat -- so by default you want M to be [0 0 0 0 0 0] (as is done in spm_coreg, which means we start from input_vol.mat as-is without any modification), NOT spm_imatrix(input_vol.mat) which would be a transform that would be applied ON TOP of input_vol.mat
% mi_algo is the mutual information like measure that will be applied to the joint histogram to calculate a score. Can be: 'mi', 'nmi', 'ecc' (default), 'ncc'.
% voxelsize is the sampling density for the joint histogram calculation, see spm_hist2.m. With [1 1 1] (default), this means that approximately each voxel will be sampled. With [4 4 4], 1 every 64 voxels will be sampled, which is less accurate but multiple of magnitude faster to compute.
% fwhm_hist is the smoothing kernel that will be applied on the histogram. By default [7 7].
%
% DEPRECATED: auxsmoothvol2uint8 is spmaux.smoothvol2uint8(), which can be memoized using a previously instanciated memoize(@aux.smoothvol2uint8) and reuse the same memoized function for each call to this function, this significantly speed up the processing here

    if ~exist('M', 'var')
        M = [0 0 0 0 0 0];
    end
    if ~isvector(M)
        iM = spm_imatrix(M);
    else
        iM = M;
    end
    if ~exist('mi_algo', 'var') || isempty(mi_algo)
        mi_algo = 'ncc';  % ncc works WAY better than other measures for this purpose
    end
    if ~exist('voxelsize', 'var') || isempty(voxelsize)
        % Sampling density for the joint histogram calculation, see spm_hist2.m. With [1 1 1], this means that approximately each voxel will be sampled. With [4 4 4], 1 every 64 voxels will be sampled.
        voxelsize = [1 1 1];
        % spm_coreg.optfun() expects a voxelsize, which will be kind of the sampling rate at which we will sample both images coregistration. To make it faster, it's useless to use a finer resolution than the images, so the best is to autodetect the input image's resolution and use that as the basis (we could also use the lowest resolution of the 2 images).
        %iM = spm_imatrix(input_vol.mat);
        %voxelsize = abs(iM(7:9));
        
        %iM2 = spm_imatrix(template_vol.mat);
        %voxelsize2 = abs(iM2(7:9));
        %if sum(abs(voxelsize - voxelsize2)) > 1E-6
        %    error('Input images do not have the same dimensions, please first rescale them to the same dimension before attempting a joint histogram!\n');
        %end
    end
    if ~exist('fwhm_hist', 'var') || isempty(fwhm_hist)
        fwhm_hist = [1 1]; % Default in spm_coreg is [7 7], but it's less effective in practice in our tests
    end

    % Memoize the smoothing & converter to uint8 function so that it will significantly speed up the calculations
    %smoothtouint8 = memoize(@spmaux.smoothvol2uint8);

    % spm_coreg.optfun expects images uint8 formatted (and int8 smoothed)
    %fitness = spm_coreg(iM, auxsmoothvol2uint8(template_vol, fwhm(1)), auxsmoothvol2uint8(input_vol, fwhm(2)), voxelsize, mi_algo, fwhm);
    fitness = spm_coreg(iM, template_vol_uint8, input_vol_uint8, voxelsize, mi_algo, fwhm_hist);  % this is the most time consuming part, and it's already a precompiled function so I don't know what we can do better
    % This will call spm_hist2, which we can call manually like so: spm_hist2(input_vol_uint8.uint8, template_vol_uint8.uint8, input_vol.mat, [1 1 1])
    % TODO: manually compute joint histogram and mutual information in pure matlab using accumarray? https://stackoverflow.com/questions/23691398/mutual-information-and-joint-entropy-of-two-images-matlab/23691992#23691992 - could then use fmincg or fminunc (native in matlab) or spm_powell (as is done in spm_coreg)

    % spm_coreg.optfun always negates the resulting fitness (to minimize), so we need to flip the sign around to maximize
    fitness = -1.0 * fitness;
end  % endfunction

function fitnessfunc = jointhistogramfitness(spmaux, template_vol, input_vol, fwhm)
% fitnessfunc = jointhistogramfitness(aux, template_vol, input_vol, mi_algo, voxelsize, fwhm, fwhm_hist)
% Preprocess template_vol and input_vol to uint8 (and smooth them), and return a function handle to compute the mutual information score over the joint histogram
% The main purpose of this function is to do the uint8 & smoothing preprocessing of the volumes only once, as this is computation expensive and thus slows down the genetic algorithm a lot if we don't memoize or precompute. Memoization did not work as intended, so we precompute here.
%
% fwhm is a vector of 2 elements being the value of the isotropic smoothing kernel for each image: fwhm(1) for template_vol and fwhm(2) for input_vol. Default is [20 20].
% See help jointhistogramfitness_core for the other parameters.

    % Default values
    if ~exist('fwhm', 'var')
        fwhm = [20 20];
    end

    % spm_coreg.optfun expects images uint8 formatted (and int8 smoothed), so we precompute here once and for all since we only expect the orientation matrix to change (and not the volumes)
    input_vol_uint8 = spmaux.smoothvol2uint8(input_vol, fwhm(2));
    template_vol_uint8 = spmaux.smoothvol2uint8(template_vol, fwhm(1));

    % return the jointhistogramfitness_core function properly instanciated with the preprocessed uint8 volumes
    fitnessfunc = @(x, mi_algo, voxelsize, fwhm_hist)jointhistogramfitness_core(template_vol_uint8, input_vol_uint8, x, mi_algo, voxelsize, fwhm_hist);
end

function [svol_down, svol_down_data] = niirescale(svol, scalefactor, savetofile)
% [svol_down, svol_down_data] = niirescale(svol, scalefactor, savetofile)
% Downscale or upscale any nifti volume in-memory
% svol is a nifti volume structure as provided by spm_vol()
% scalefactor is a float number with < 1.0 for downscaling and > 1.0 for upscaling. It can also be a vector of 3 scaling factors for x, y and z (eg, size(template_data) or template_vol.dim)
% The rescaled volume can be saved in a file by setting savetofile = 'filename.nii' or manually by: svol_down.fname = 'newfilename.nii'; spm_write_vol(svol_down, svol_down_data);
% The rescaled volume will have its origin repositionned correctly. Also, the orientation matrix will report the correct rescaling factors, hence the voxel-to-world mapping will be correct in nifti viewers.
% TODO: need to implement interpolation for upscaling. For the moment, this function works great for downscaling (most common use case), but not for upscaling. See interp2.m or better waifu2x.

    % Default values
    if ~exist('savetofile', 'var')
        savetofile = false;
    end  % endif

    % Read the voxels content of the input volume
    svol_data = spm_read_vols(svol);
    % Resize image's data with interpolation
    svol_down_data = imresize3(svol_data, scalefactor, 'method', 'lanczos3', 'Antialiasing', true);  % works only with MATLAB >= v2016 - alternative for older versions: use imresizen.m from: https://www.mathworks.com/matlabcentral/answers/358043-resizing-a-3d-image-without-using-imresize3
    %Y2 = real(ifftn(fftn(Y), input_vol.dim ./ 2.0));  % alternative but huge ringing artifacts, see https://www.mathworks.com/matlabcentral/answers/358043-resizing-a-3d-image-without-using-imresize3 and https://www.researchgate.net/post/How_can_I_avoid_Ringing_artifacts_in_Fourier_transformation_in_SCILAB and https://en.wikipedia.org/wiki/Ringing_artifacts to enhance by filtering
    % Make a copy of the input volume
    svol_down = svol;
    % Scale its dimensions
    svol_down.dim = round(svol.dim .* scalefactor);
    % Translate origin first, by dividing distances
    % after this, the image is already looking good, exactly like the original. But the voxel dimension reported in the orientation matrix is incorrect.
    svol_down.mat = svol.mat;
    svol_down.mat(1:3, 4) = svol_down.mat(1:3, 4) .* scalefactor;
    % Calculate the inverse of the scale factor, because what we want to do is to do the opposite of the real operation we did on the data: if we downscale the data, the voxel dimensions in the orientation matrix must be upscaled to compensate, eg, if we downscale to 0.5, the voxel size will be 2x as big as before. This will also change the rotation and translation values of course so that we compensate for the upscaling in order to maintain the same origin and orientation.
    invscale = (1/scalefactor);
    % Scale its orientation matrix (including the origin/translation again! Hence the translation/origin part of the matrix needs to be scaled in square - eg, if scalefactor == 0.5, then the translation factors need to be .* 0.5^2 == 0.25)
    % This corrects the voxel dimension, so that nifti viewers know that the voxels are rescaled and hence can do a correct voxel-to-world remapping
    svol_down.mat = diag([ones(1,3) .* invscale, 1]) * svol_down.mat;  % set 1 for the last value in the rescaling vector so that we maintain the intercept/constant factor to 1 - Deprecated: not working correctly for setting the origin/translation!
    % Save to a new nifti file if option is enabled
    if savetofile ~= false
        svol_down.fname = savetofile;
        spm_write_vol(svol_down, svol_down_data);
    end
end  % endfunction

function mga_allinone
    %autoreorient('t1.nii', 'affine', [], false, false, 'svd', false, false, true);

    %downscaled image
    [input_vol_down, input_vol_down_data] = niirescale(input_vol, 0.5);
    iM_down = spm_imatrix(input_vol_down.mat);
    template_vol_down = niirescale(template_vol, 0.5);
    fitnessfunc_core = jointhistogramfitness(spmaux, template_vol, input_vol, [1 10]);
    fitnessfunc = @(x)fitnessfunc_core(x, 'ncc', [3 3 3], [1 1]);
    [best, bestfit] = mga_imatrix([0 0 0 0 0 0], 10, 3, fitnessfunc, 0.7, 0.25, 1000, 2, true, input_vol_down.dim, [], [], false)
    newM = (spm_matrix(best)/input_vol_down.mat)*input_vol.mat;
    if (sum(abs(newM - input_vol.mat), 'all') < 1E-6) || (abs(bestfit - fitnessfunc(input_vol_down.mat)) < 1E-6)
        fprintf('Could not find a better orientation!');
    else
        spm_get_space(inputpath, pinv(spm_matrix(best))*input_vol.mat);  % simple rule of 3: project the ratio between the best fit orientation matrix over the downscaled image's orientation matrix, and then apply this on the original fullscale volume, which leads to the correct projection of the best orientation found on the downscaled volume but here in the fullscale space.
        input_vol_down.fname = 't1_down.nii';
        input_vol_down.mat = pinv(spm_matrix(best))*input_vol_down.mat;
        spm_write_vol(input_vol_down, input_vol_down_data);
    end

    % Bonus
    [input_vol_down, input_vol_down_data] = niirescale(input_vol, 0.25);
    iM_down = spm_imatrix(input_vol_down.mat);
    input_vol_down.fname = 't1_down.nii';
    spm_write_vol(input_vol_down, input_vol_down_data);

    %------------
    % full image - works best, no need to rescale
    fitnessfunc_core = jointhistogramfitness(spmaux, template_vol, input_vol, [1 20]);
    fitnessfunc = @(x)fitnessfunc_core(x, 'ncc', [8 8 8], [1 1]);  % working best with ncc
    % TODO: redo another try with [1 1 1] resampling step?
    rng(0,'philox');
    custrng = rng;
    [best, bestfit] = mga_imatrix([0 0 0 0 0 0], 10, 3, fitnessfunc, 0.7, 0.25, 500, 5, true, input_vol.dim, [], [], custrng, false)
    [best, bestfit] = mga_imatrix(best, 10, 3, fitnessfunc, 0.7, 0.25, 500, 5, true, input_vol.dim, [], [], custrng, false)
    fitnessfunc_precise = @(x)fitnessfunc_core(x, 'ncc', [4 4 4], [1 1]);  % working best with ncc
    [best, bestfit] = mga_imatrix(best, 10, 3, fitnessfunc_precise, 0.7, 0.25, 300, 2, true, input_vol.dim, [], [], custrng, false)
    fitnessfunc_precise = @(x)fitnessfunc_core(x, 'ncc', [2 2 2], [1 1]);  % working best with ncc
    [best, bestfit] = mga_imatrix(best, 10, 3, fitnessfunc_precise, 0.7, 0.25, 200, 1, true, input_vol.dim, [], [], custrng, false)
    %fitnessfunc_precise = @(x)fitnessfunc_core(x, 'ncc', [4 4 4], [1 1]);  % working best with ncc
    %[best2, bestfit2] = mga_imatrix(best, 10, 3, fitnessfunc_precise, 0.7, 0.25, 500, 2, true, input_vol.dim, [], [], false)
    %fitnessfunc_precise = @(x)fitnessfunc_core(x, 'ncc', [2 2 2], [1 1]);  % working best with ncc
    %[best3, bestfit3] = mga_imatrix(best2, 10, 3, fitnessfunc_precise, 0.7, 0.25, 200, 2, true, input_vol.dim, [], [], false)
    newM = spm_matrix(best);
    if (sum(abs(newM - input_vol.mat), 'all') < 1E-6) || (abs(bestfit - fitnessfunc(input_vol.mat)) < 1E-6)
        fprintf('Could not find a better orientation!');
    else
        spm_get_space(inputpath, pinv(newM) * input_vol.mat);  % spm_coreg returns an orientation vector that needs to be applied ON TOP of the initial input_vol.mat orientation matrix
    end

    % to redo from previous best
    %[best, bestfit] = mga_imatrix(best, 10, 3, fitnessfunc, 0.7, 0.25, 1000, 1, false, input_vol.dim, [], [], false)

    % just test fitness with current input_vol
    fitnessfunc_core = jointhistogramfitness(spmaux, template_vol, input_vol, [1 20]);
    fitnessfunc = @(x)fitnessfunc_core(x, 'ncc', [8 8 8], [1 1]);
    fitnessfunc([0 0 0 0 0 0])
    input_orig = spm_vol('t1_orig.nii');
    fitnessfunc(spm_imatrix(input_vol.mat/input_orig.mat))  % can then do pinv(input_vol.mat/input_orig.mat) * input_vol.mat
    fitnessfunc(best)

    % ------------------

    % full image v2 - works best, no need to rescale
    fitnessfunc_core = jointhistogramfitness(spmaux, template_vol, input_vol, [1 20]);  % use [1 20] to smooth the input_vol to a kernel of 20 as is the default in auto_reorient by John Ashburner and Carlton Chu, but in our experience [1 1] (no smoothing) works better
    fitnessfunc = @(x)fitnessfunc_core(x, 'ncc', [8 8 8], [1 1]);  % working best with ncc
    % TODO: redo another try with [1 1 1] resampling step?
    [best, bestfit, ~, bestpool, bestpoolfit] = mga_imatrix([0 0 0 0 0 0], 10, 3, fitnessfunc, 0.7, 0.25, 500, 20, true, input_vol.dim, [], [], false, true, false)
    %[best, bestfit] = mga_imatrix(best, 10, 3, fitnessfunc, 0.7, 0.25, 500, 5, true, input_vol.dim, [], [], false, false)
    %[best, bestfit] = mga_imatrix(best, 10, 3, fitnessfunc, 0.7, 0.25, 300, 2, true, input_vol.dim, [], [], false, false)
    % select best amongst all candidates
    fitnessfunc_precise = @(x)fitnessfunc_core(x, 'ncc', [1 1 1], [1 1]);
    bestpoolfit2 = bestpoolfit;
    for i=1:size(bestpool, 1)
        bestpoolfit2(i) = fitnessfunc_precise(bestpool(i,:));
    end
    [~, bestidx] = max(bestpoolfit2);
    best = bestpool(bestidx, :);
    newM = spm_matrix(best);
    if (sum(abs(newM - input_vol.mat), 'all') < 1E-6) || (abs(bestfit - fitnessfunc(input_vol.mat)) < 1E-6)
        fprintf('Could not find a better orientation!');
    %else
        %spm_get_space(inputpath, pinv(newM) * input_vol.mat);  % spm_coreg returns an orientation vector that needs to be applied ON TOP of the initial input_vol.mat orientation matrix
    end

    % to redo from previous best
    %[best, bestfit] = mga_imatrix(best, 10, 3, fitnessfunc, 0.7, 0.25, 1000, 1, false, input_vol.dim, [], [], false)

    % just test fitness with current input_vol
    %fitnessfunc_core = jointhistogramfitness(spmaux, template_vol, input_vol, [1 20]);
    fitnessfunc = @(x)fitnessfunc_core(x, 'ncc', [1 1 1], [1 1]);
    fitnessfunc([0 0 0 0 0 0])
    input_orig = spm_vol('t1_orig.nii');
    fitnessfunc(spm_imatrix(input_vol.mat/input_orig.mat))  % can then do pinv(input_vol.mat/input_orig.mat) * input_vol.mat
    fitnessfunc(best)



    % ------------

    input_orig = spm_vol('t1_orig.nii');
    origbest = spm_imatrix(input_vol.mat/input_orig.mat);
    H = spmaux.jointhist(template_vol, input_vol, [0]);
    H2 = spmaux.jointhist(template_vol, input_vol, origbest);
    H3 = spmaux.jointhist(template_vol, input_vol, best);



    % ----------

    % Compute the joint histogram between two 3D matrices of voxels
    % https://stackoverflow.com/questions/23691398/mutual-information-and-joint-entropy-of-two-images-matlab/23691992#23691992
    template_data = spm_read_vols(template_vol);
    input_data = spm_read_vols(input_vol);
    % Rescale input data to match the template size (it's better to do that before rounding to uint8 for precision)
    input_data_down = imresize3(input_data, size(template_data), 'method', 'lanczos3', 'Antialiasing', true);
    % Intensity normalization via feature scaling by normalizing into a int8 range (0-255), so both images can be matched (else they can have very different values), so we end up with a 256*256 joint histogram - see https://en.wikipedia.org/wiki/Feature_scaling#Rescaling_(min-max_normalization)
    input_data_norm = aux.featurescale_antialias(input_data_down, 0, 255);
    template_data_norm = aux.featurescale_antialias(template_data, 0, 255);
    % Compute the joint histogram
    jhist = accumarray([input_data_norm(:)+1, template_data_norm(:)+1], 1);  % shift values by 1 to account for the fact that MATLAB starts indexing at 1. Note also that both arrays need to be exactly the same size.
    jprob = jhist / numel(input_data_norm);
    % Display the joint histogram as a greyscale chart of the log
    aux.plot_jointhist(jhist);
    % Compute the measure
    weightedcoveragemeasure = @(jhist)log(nnz(jhist)) + log(sum(sum(jhist)));
    weightedcoveragemeasure(jhist);
    % TODO apply orientation matrix before calculating the joint hist

    aux.display_mri3d(aux.apply_transform3d(input_data, origbest, true));

    % TODO: tester les reorientation matrix et plot pour voir si correct, et ensuite plug dans code ci-dessus avant jointhist calculation
    %-------

    spmaux.jointhist_display(template_vol, input_vol, [0], [], struct('cost_fun', 'ncc'))

    spmaux.jointhist_score(H, 'ncc')
    spmaux.jointhist_score(H2, 'ncc')
    spmaux.jointhist_score(H3, 'ncc')

    spmaux.jointhist_display(template_vol, input_vol, x, [], struct('cost_fun', 'ncc'))
    spmaux.jointhist_display(template_vol, input_vol, best, [], struct('cost_fun', 'ncc'))

    log(nnz(H3)) + log(sum(H3, 'all'))

    makehgtform

    sign(data) .* log( abs( data ) ); % https://www.mathworks.com/matlabcentral/answers/253936-logarithmic-range-colorbar-imagesc
    % BEST: cptcmap.m for custom color maps with custom ranges: https://www.mathworks.com/matlabcentral/answers/253936-logarithmic-range-colorbar-imagesc
end


%%% Additional resources
% spm_get_space usage found in this BSD-2 licensed script? Can't remember https://github.com/jimmyshen007/NeMo/blob/c7cc775e7fc84f84255b0d3fa474b7f955e817f5/mymfiles/eve_tools/Structural_Connectivity/approx_coreg2MNI.m

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

% tosee https://github.com/jewettaij/superpose3d

% best results, rotation-only and no rescaling, retaining voxel-to-world mapping: autoreorient('t1.nii', 'mi', [], false, false, 'none', false, 'raw', true)
% works awesomely for hard cases even with a hard sheared input image, but losing voxel-to-world mapping: autoreorient('t1.nii', 'mi', [], false, false, 'none', false, 'raw_scale', true)
