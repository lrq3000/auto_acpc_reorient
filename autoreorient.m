function autoreorient(inputpath, mode, flags_affine, noshearing, isorescale, affdecomposition)
    if ~exist('noshearing', 'var')
        noshearing = true;
    end
    if ~exist('isorescale', 'var')
        isorescale = true;
    end
    if ~exist('affdecomposition', 'var')
        affdecomposition = 'qr';
    end

    input_vol = spm_vol(strtrim(inputpath));
    smoothed_vol = spm_vol(inputpath);
    smoothed_vol2 = spm_smoothto8bit(smoothed_vol, 20);
    % Y = spm_read_vols(smoothed_vol2);  % direct access to the data matrix
    % imagesc(Y(:,:,20))  % show a slice
    template_vol = spm_vol(strtrim('C:\git\auto_acpc_reorient\T1_template_CAT12_rm_withskull.nii'));
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
            assert(size(Mnotrans, 1) == size(Mnotrans, 2))  % only allow square matrices, not sure if it works on non-square ones
            assert(single(det(Mnotrans)) ~= 0.0)  % determinant is non-zero

            if strcmp(affdecomposition, 'qr') && (noshearing || (isempty(isorescale) || isorescale))
                % In QR decomposition mode, we will decompose the 3D affine transform matrix M into a rotation matrix Q, a scaling matrix K, and a shearing matrix S.
                fprintf('Using QR decomposition...\n')
                % QR Decomposition to extract the rotation matrix Q and the scaling+shearing matrix R
                % See https://math.stackexchange.com/questions/1120209/decomposition-of-4x4-or-larger-affine-transformation-matrix-to-individual-variab/2353782#2353782
                [Q, R] = qr(Mnotrans);  % pure MATLAB implementations can be found at: http://web.mit.edu/18.06/www/Essays/gramschmidtmat.pdf and https://blogs.mathworks.com/cleve/2016/07/25/compare-gram-schmidt-and-householder-orthogonalization-algorithms/ and https://blogs.mathworks.com/cleve/2016/10/03/householder-reflections-and-the-qr-decomposition/#742ad311-5d1b-4b5d-8768-11799d40f72d
                %[Q, R] = lu(Mnotrans);  % alternative
                % Enforce positive main diagonal, but at the expense of making an improper rotational matrix Q (ie, with reflection), so it's not satisfactory either... https://www.mathworks.com/matlabcentral/answers/83798-sign-differences-in-qr-decomposition#answer_93385
                %D = diag(sign(diag(R)));
                %Q = Q*D;
                %R = D*R;
                % Decompose R into the scaling matrix K and shearing matrix S
                % Compute the scaling matrix K and its inverse, by simply creating K from the main diagonal of R (the rest of the matrix being zero)
                K = eye(size(Mnotrans));
                K(1:size(R,1)+1:end) = diag(R);  % anisotropic rescaling (necessary to compute the accurate inverse and then the shearing matrix S)
                if isorescale
                    % Compute the isotropic scaling matrix K if required
                    Kiso = eye(size(Mnotrans));
                    %Kiso(1:size(R,1)+1:end) = mean(abs(diag(R))) .* sign(diag(R));  % compute the average rescaling factor, so we rescale but isotropically (ie, the same rescaling factor in all directions). We first calculate the mean of absolute values, and then we restore the sign. The sign is important because the decomposition may place in the rescaling matrix K a mirroring rescale (ie, the rotation Q will rotate in the opposite direction as the one we want, but rescaling in the negative direction -1 will mirror the brain and put it back in place). - DEPRECATED: big issue: if some scaling values have a negative sign, this will produce a mirroring which will inverse the place of gray matter and the rest of the brain! Not good, left side of brain will be on the right side, we certainly don't want that!
                    Kiso(1:size(R,1)+1:end) = mean(abs(diag(R)));
                end
                Kinv = eye(size(Mnotrans));
                Kinv(1:size(R,1)+1:end) = 1 ./ diag(R);  % the inv(K) == the reciprocals (inverse) of the diagonal values in R
                % Compute the shearing matrix S, by removing K from R (so the main diagonal will be zero in S)
                S = Kinv*R;
                % Sanity checks
                % See also what each 3D affine transform look like: https://www.tutorialspoint.com/computer_graphics/3d_transformation.htm and https://stackoverflow.com/questions/5107134/find-the-rotation-and-skew-of-a-matrix-transformation?rq=1
                % Check if there is any reflection, we don't want that
                if (single(det(Q)) == -1.0) || any(sign(diag(K)) == -1)
                % TODO: rotate 180 degrees when reflection is detected and skip error
                    error('QR decomposition led to an improper rotational matrix Q (rotation+reflection) or scaling matrix K (negative factors), cannot continue.');
                end
                % Q is a proper rotation matrix if its determinant equals to 1
                assert(single(det(Q)) == 1.0);  % a proper rotation matrix will have a determinant of 1. This means that it's a "pure" rotational matrix, with no other transform (such as reflection/mirroring). See: https://www.tau.ac.il/~tsirel/dump/Static/knowino.org/wiki/Euler%27s_theorem_(rotation).html#Matrix_proof
                % K is a scaling matrix if all non-diagonal values are 0
                assert(all((triu(K, 1)+tril(K,-1)) == 0, [1 2]));
                % S is a shearing matrix if:
                assert(single(det(S)) == 1.0);  % its determinant equals to 1
                assert(all(diag(S) == 1));  % all its main diagonal values equal to 1
                assert(all(eig(S) == 1));  % all eigenvalues equal to 1
                assert((rank(S) == trace(S)) && (trace(S) == size(Mnotrans, 1)));  % its trace equals its rank equals the NxN size of the 3D affine transform matrix
                % Final check of the decomposition, we should find the original 3D affine transform matrix (without the translation of course)
                assert(all(single(Mnotrans) == single(Q*K*S), [1 2]));  % convert from doubles to float to clean up the rounding errors during the decomposition

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
                fprintf('Using SVD decomposition...\n')
                % Calculate the SVD decomposition
                [u, s, v] = svd(Mnotrans);
                if isorescale || noshearing
                    % Compute the mean of the rescaling factor to do only isotropic rescaling (and not anisotropic)
                    s = diag(mean(abs(diag(s))) * [1 1 1]);
                end
                % Sanity check: prevent any reflection/mirroring
                if (single(det(u)) == -1 || single(det(v')) == -1 || any(sign(diag(s)) == -1))
                    % TODO: rotate 180 degrees when reflection is detected and skip error
                    error('Reflections detected in the decomposed matrices, cannot continue!');
                end
                M2 = u*s*v';  % apply other transforms apart from translation (the precise transforms being defined by user variables)
                % Reconstruct a 4x4 affine transform matrix (so that we can put back the translation column vector)
                M3 = eye(size(M));  % the basis needs to be an identity matrix, see: https://www.youtube.com/watch?v=UvevHXITVG4
                M3(1:3,1:3) = M2;  % apply other transforms apart from translation (the precise transforms being defined by user variables)
                M3(:,4) = T;  % Restore the translation
                disp(M);
                disp(M2);
                disp(M3);
                M = M3;
            elseif strcmp(affdecomposition, 'imatrix')
                % Decompose and constraint reflections/shearing using spm_imatrix
                % DEPRECATED: this does NOT work because the parameters are dependent on the order in which they are applied (default: inverse of 'T*R*Z*S' in SPM, see help spm_matrix), so it's not possible to simply modify some parameters and leave the rest as-is, the whole final transform gets out of control. A decomposition is necessary.

                % Extract the transform parameters as a vector from the 3D affine transform matrix
                iM = spm_imatrix(M);
                % Multiple assign each parameter into its own variable (see help spm_matrix: T for translation, R rotation, Z zoom/scale, S shearing)
                iM2 = mat2cell(iM, 1, ones(1,4)*3);
                [T, R, Z, S] = iM2{:};
                % Nullify shearing if option is enabled
                if noshearing
                    S = zeros(1, 3);
                end
                % If there are any reflections in the scaling matrix Z (negative factors), then apply a 180° rotation (pi) instead in the same axis (this is not equivalent but it's better than nothing - the best is to afterward do another coregistration, based on joint histogram for example, to clean up the residual coregistration errors)
                R = R + (sign(Z) == -1)*pi;  % TODO: buggy
                % Turn negative scaling factors to positive (ie, this removes reflections)
                Z = Z .* sign(Z);
                if isorescale
                    % If isotropic scaling is required by user, compute the mean rescaling and apply the same scaling factor in all dimensions
                    Z = ones(1,3) * mean(abs(Z));
                end
                % Convert the vector back to a 3D affine transform matrix
                % Note that the order of the operations can be changed, see help spm_matrix for more information about default order and how to change
                M = spm_matrix([T R Z S], 'T*R*Z*S');
            end %endif
            % TODO: recenter the origin to middle of the brain?
        end %endif
        % Sanity check: check if there is no reflections nor shearing
        [T, R, Z, S] = get_imat(M);
        assert(~any(sign(Z) == -1));
        if noshearing
            assert(sum(S) < 0.000001);
        end
        % Apply the reorientation on the input image header and save back the new coordinations in the original input file
        % Some inspiration by BSD-2 licensed: https://github.com/jimmyshen007/NeMo/blob/c7cc775e7fc84f84255b0d3fa474b7f955e817f5/mymfiles/eve_tools/Structural_Connectivity/approx_coreg2MNI.m - licensed under BSD-2-Clause (compatible with MIT license anyway)
        %spm_get_space(inputpath, M*input_vol.mat);
        spm_get_space(inputpath, M*input_vol.mat);
    else
        fprintf('Mutual information reorientation\n');
        flags = struct('sep', 5);
        % spm_coreg can only find a rigid-body rotation, hence spm_affreg is necessary for translation and (isotropic) scaling first, and also there is no need to do a QR decomposition here
        x = spm_coreg(template_vol, smoothed_vol2, flags);
        M = inv(spm_matrix(x));
        % Apply the reorientation and save back the new coordinations
        spm_get_space(inputpath, M*spm_get_space(inputpath));
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
