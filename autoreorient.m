function autoreorient(inputpath, mode, flags_affine, noshearing, isorescale, affdecomposition)
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
                    [v, s] = improper_rotation_to_proper(v', s);
                    v = v';
                end %endif
                % Remove/balance out the reflections in s
                s = remove_reflections(s);
                % Isotropic rescaling / no shearing?
                if isorescale || noshearing
                    % Compute the mean of the rescaling factor to do only isotropic rescaling (and not anisotropic)
                    % In SVD mode, the shearing is done in the scaling matrix, so we can't separate both transforms, if we remove one we remove both
                    s = diag(mean(abs(diag(s))) * sign(diag(s)));  % don't forget to restore the sign to balance out the reflections
                end
                % Sanity check: check if there is any remaining imbalanced reflection/mirroring
                if single(det(u)) == -1 || single(det(v')) == -1 || (any(sign(diag(s)) == -1) && are_reflections_imbalanced(s))
                    error('Reflections detected in the decomposed matrices, cannot continue!');
                end
                % Combine back all transforms
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
                % Decompose and constraint reflections/shearing using spm_imatrix (fixed)
                fprintf('Using spm_imatrix decomposition...\n')

                % Extract the transform parameters as a vector from the 3D affine transform matrix
                iM = spm_imatrixfix(M);
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
        fprintf('Mutual information reorientation\n');
        flags = struct('sep', 5);
        % spm_coreg can only find a rigid-body rotation, but pre-translation is not necessary for spm_coreg to work because it's working directly on joint histograms, so empty space is meaningless. Also there is no need to do a QR decomposition here since it's rigid-body.
        % however if later we want to coregister another modality (eg, functional) over the structural, we need to fix the translation too. In this case, spm_affreg is necessary for translation and (isotropic) scaling first on the structural.
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
