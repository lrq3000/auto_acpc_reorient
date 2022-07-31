function funs = auto_acpc_reorient_aux
% funs = auto_acpc_reorient_aux
% Auxiliary functions to plot and such things for acpc_auto_reorient
%
% License: MIT License
% Copyright (C) 2020 Stephen Karl Larroque, Coma Science Group & GIGA-Consciousness, University Hospital of Liege, Belgium
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, includin without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    funs = importFunctions; % For Octave compatibility: we need the first function to have the same name as the filename.
end

% For MatLab compatibility, we use a function to return other functions handlers as properties, this is a workaround since it cannot load multiple functions in one .m file (contrarywise to Octave using source())
% This way, we can just do the following to call any function here from any other script:
% to load the aux lib: addpath(genpath(strcat(cd(fileparts(mfilename('fullpath'))),'/../gbnn-core/'))); aux = gbnn_aux;
% to use a function: aux.func_name(args);
function funs = importFunctions
    funs.display_mri3d = @display_mri3d;
    funs.plot_jointhist = @plot_jointhist;
    funs.matrix2coords = @matrix2coords;
    funs.apply_transform3d = @apply_transform3d;
    % Intensity normalization via feature scaling by normalizing into a int8 range (0-255), so both images can be matched (else they can have very different values), so we end up with a 256*256 joint histogram - see https://en.wikipedia.org/wiki/Feature_scaling#Rescaling_(min-max_normalization)
    funs.featurescale_antialias = @(svol_data, rangemin, rangemax)uint8(round(((svol_data + (rand(size(svol_data)) .* min(svol_data(nnz(svol_data)), [], 'all')) - min(svol_data, [], 'all')) .* (rangemax - rangemin)) ./ (max(svol_data, [], 'all') - min(svol_data, [], 'all'))));  % we add a random value to avoid aliasing artifacts, which is scaled according to the minimal non-zero value found in the input data (so we don't mess up significantly with the data)
    funs.featurescale = @(svol_data, rangemin, rangemax)uint8(round(((svol_data - min(svol_data, [], 'all')) .* (rangemax - rangemin)) ./ (max(svol_data, [], 'all') - min(svol_data, [], 'all'))));
    funs.getnargs=@getnargs; % to process named optional arguments
    funs.varspull = @varspull; % to load arguments into local namespace/workspace
    funs.editarg=@editarg;
    funs.delarg=@delarg;
    funs.addarg=@addarg;
    funs.printcputime=@printcputime;
    funs.printtime=@printtime;
    funs.printeta=@printeta;
end  % endfunction

function fh = display_mri3d(svol_data, M, figtitle)
% fh = display_mri3d(svol)
% Display 3 orthogonal slices of the provided volume data.
% svol_data can either be spm_get_vols(svol) or svol. In the first case, this allows to apply a reorientation matrix beforehand.
% M is a reorientation matrix, or a spm_imatrix
% figtitle is the title of the figure (can be set manually after using the fh handle)
%
% Outputs fh, the handle of the figure

    % Sanity checks
    if isstruct(svol_data)
        % Ensure svol_data is the volume's data, not a volume struct
        svol_data = spm_get_vols(svol_data);
    end
    if exist('M', 'var') && ~isempty(M)
        % If a reorientation matrix M was supplied
        if isvector(M)
            % If it's a spm_imatrix, convert it to an affine transform matrix
            M = spm_matrix(M);
        end
        % Apply the orientation matrix M on svol_data
        svol_data = M*svol_data;
    end  % endif

    % Plot the 3D orthogonal slices
    % From: https://www.mathworks.com/help/images/ref/imresize3.html
    fh = figure;
    sizeD = size(svol_data);
    slice(double(svol_data), sizeD(2)/2, sizeD(1)/2, sizeD(3)/2);
    shading interp, colormap gray;
    if exist('figtitle', 'var'), title(figtitle); end;
end  % endfunction

function fh = plot_jointhist(jhist, precisegrey)
% fh = plot_jointhist(jhist, precisegrey)
% Plot a joint histogram in greyscale
% jhist is a joint histogram computed with the jhist function
% if precisegrey is true, we will use a larger greyscale colormap to represent the histogram with more precision (ie, more shades of grey). False by default.
%
% Outputs fh, the figure handle

    if ~exist('precisegrey', 'var') || isempty(precisegrey)
        precisegrey = false;
    end

    fh = figure;
    if precisegrey
        imagesc(log(jhist));
    else
        imagesc(uint8(round(log(jhist))));
    end
    % Plot as greyscale
    % https://www.mathworks.com/help/matlab/ref/contrast.html
    colormap(contrast(log(jhist)));
end  % endfunction

function [x, y, z, v] = matrix2coords(svol_data)
% [x, y, z] = matrix2coords(svol_data)
% Extract the coordinates of each voxel and its associated value, each in a different vector (we can combine into a long-form table later). This can be used to apply transform on coordinates (or do calculations such a finding the centroid).
% svol_data is the volume's voxels data as returned by spm_get_vols(svol)
% Outputs: x, y and z vectors of coordinates, and v the vector of values

    if isstruct(svol_data)
        % If a volume structure was supplied, extract the voxels content
        svol_data = spm_get_vols(svol_data);
    end

    idx = find(svol_data);  % extract linear/vector style indices
    v = svol_data(idx);  % extract values in a vector
    [x, y, z] = ind2sub(size(svol_data), find(svol_data));  % convert linear indices to 3D indices/coordinates
end  %endfunction

function svol_data_reoriented = apply_transform3d(svol_data, M, puremode)
% svol_data_reoriented = apply_transform3d(svol_data, M, puremode)
% Apply/reorient a 3D matrix of voxels with the provided transform matrix M
% svol_data can either be spm_get_vols(svol) or svol. In the first case, this allows to apply a reorientation matrix beforehand.
% M is a reorientation matrix, or a spm_imatrix.
% puremode set to true will use the pure MATLAB implementation we have written, else if false (default) will use MATLAB's native function imwrap() which should theoretically be faster.

    % Default values
    if ~exist('puremode', 'var') || isempty(puremode)
        puremode = false;
    end
    % Sanity checks
    if isvector(M)
        % If M is a spm_imatrix vector, convert to a transform matrix
        M = spm_matrix(M);
    end
    if isstruct(svol_data)
        % If svol_data is a nifti volume struct, we extract the voxels content
        svol_data = spm_read_vols(svol_data);
    end

    if ~puremode
        % Simpler: use MATLAB native functions
        svol_data_reoriented = imwarp(svol_data, affine3d(M'));  % in MATLAB, the default is for transform matrix to be in row-major format, whereas for SPM it's column-major format, so we need to transpose
    else
        % More complex but in pure MATLAB: calculate by ourselves
        % TODO: upscale before rotation to lose less quality, see: https://www.mathworks.com/matlabcentral/answers/387742-using-imrotate-without-losing-quality
        % see also: https://blogs.mathworks.com/steve/2019/04/09/multiresolution-image-pyramids-and-impyramid-part-2/?s_tid=answers_rc2-3_p6_Topic
        % and: https://blogs.mathworks.com/steve/2007/08/13/image-deblurring-introduction/
        % First transform the matrix from dense format to long format, in other words extract the coordinates in a 4-rows/columns table
        [x, y, z, v] = matrix2coords(svol_data);
        % Concatenate all coordinates axes in one 4-rows table/matrix (1st row is x, 2nd row is y, etc).
        coords = [x';y';z';ones(1,numel(x))];  % also add a 4th row for the constant term 1 (intercept)
        % Apply the transform matrix by a simple matrix multiplication
        reoriented_sub = round(M * coords);  % This is a forward mapping, contrary to most software implementation that use the inverse mapping, see: https://stackoverflow.com/questions/30437467/how-imwarp-transfer-points-in-matlab
        % DEPRECATED: Wrap-around, which ensures we do not exceed the boundaries of the matrix dimensions
        %reoriented_sub = mod(reoriented_sub, [size(svol_data)'; 1]) + 1;  % need to add 1 as matlab indexes from 1, not 0
        % Resize the picture so it stays in bounds and we can properly plot it
        featurescale_perrow = @(svol_data, rangesmin, rangesmax)uint8(round(((svol_data - min(svol_data, [], 2)) .* (rangesmax - rangesmin)) ./ (max(svol_data, [], 2) - min(svol_data, [], 2)))+1);
        reoriented_sub(1:3,:) = featurescale_perrow(reoriented_sub(1:3,:), ones(numel(size(svol_data)), 1), size(svol_data)');
        % Sanity checks
        assert(all(min(reoriented_sub(1:3,:), [], 2)' == 1));  % minimum index is 1
        assert(all(max(reoriented_sub(1:3,:), [], 2)'  <= size(svol_data)));  % max index is not more than original size
        %assert(all(max(reoriented_sub(1:3,:), [], 2)' == size(svol_data)));  % max index is the same as original size (not necessarily the case)
        % Convert the new coordinates (in subscript format) to indices
        %a = reshape(spm_matrix([0])*svoldata(:), size(input_data_down));
        reoriented_idx = sub2ind(size(svol_data), reoriented_sub(1,:), reoriented_sub(2,:), reoriented_sub(3,:));
        % Interpolation
        xo = reoriented_sub(1,:);
        yo = reoriented_sub(2,:);
        zo = reoriented_sub(3,:);
        % prepare a list of all indices
        all_idxs = 1:numel(svol_data);
        % find the indices that are not in reoriented_idx (will need to be interpolated)
        missing_idxs = setdiff(all_idxs, reoriented_idx);
        % convert the indices to subscripts
        [x2, y2, z2] = ind2sub(size(svol_data), missing_idxs);
        % Interpolate missing voxels
        % see also manual interpolation (but slower, in a for loop) in pure matlab: https://www.youtube.com/watch?v=m0a5EDt12cQ
        % TOSEE how to use: https://www.mathworks.com/matlabcentral/answers/410105-interpolate-missing-data-and-make-line-smooth
        %v2 = interp3(xo,yo,zo,v,x2,y2,z2);
        % Create a zero-filled volume
        svol_data_reoriented = zeros(size(svol_data));
        % Assign values from the old volume to the new reoriented indices
        svol_data_reoriented(reoriented_idx) = v;
    end  % endif
end  % endfunction

function temporary_reorient()
    % WORKS with my own implementation of transform! But 1- need interpolation, 2- rescaling need to keep proportions (hence need to rescale according to the biggest axis rescaling (ie, axis with biggest trespassing of dimension), not differently for each!), 3- origin is not set properly (except in SliceBrowser but that may be a lucky guess)
    input_data = spm_read_vols(input_vol);
    input_orig = spm_vol('t1_orig.nii');
    origbest = spm_imatrix(input_vol.mat/input_orig.mat);
    %warpedim = imresize3(imwarp(input_data, affine3d(input_orig.mat')), size(input_data)); % not the same as spm_get_space(inputpath, input_orig.mat)
    warpedim = aux.apply_transform3d(input_data, input_orig.mat, true);  % working!!!
    input_vol2 = input_vol;
    input_vol2.fname = 't1_warp.nii';
    tonly = spm_matrix([0]);
    %tonly(:,4) = input_orig.mat(:,4);
    input_vol2.mat = tonly;
    spm_write_vol(input_vol2, warpedim); % same a slicebrowser, correct orientation and all!
    SliceBrowser(warpedim) % same as writing the volume EXCEPT that the orientation convention is different, so don't compare with MRIcron
    aux.display_mri3d(warpedim)
    % https://www.mathworks.com/matlabcentral/fileexchange/20604-3d-slice-viewer
end

function argStruct = getnargs(varargin, defaults, restrict_flag)
%GETNARGS Converts name/value pairs to a struct (this allows to process named optional arguments).
%
% ARGSTRUCT = GETNARGS(VARARGIN, DEFAULTS, restrict_flag) converts
% name/value pairs to a struct, with defaults.  The function expects an
% even number of arguments in VARARGIN, alternating NAME then VALUE.
% (Each NAME should be a valid variable name and is case sensitive.)
% Also VARARGIN should be a cell, and defaults should be a struct().
% Optionally: you can set restrict_flag to true if you want that only arguments names specified in defaults be allowed. Also, if restrict_flag = 2, arguments that aren't in the defaults will just be ignored.
% After calling this function, you can access your arguments using: argstruct.your_argument_name
%
% Examples:
%
% No defaults
% getnargs( {'foo', 123, 'bar', 'qwerty'} )
%
% With defaults
% getnargs( {'foo', 123, 'bar', 'qwerty'} , ...
%               struct('foo', 987, 'bar', magic(3)) )
%
% See also: inputParser
%
% Authors: Jonas, Richie Cotton and LRQ3000
%

    % Extract the arguments if it's inside a sub-struct (happens on Octave), because anyway it's impossible that the number of argument be 1 (you need at least a couple, thus two)
    if (numel(varargin) == 1)
        varargin = varargin{:};
    end

    % Sanity check: we need a multiple of couples, if we get an odd number of arguments then that's wrong (probably missing a value somewhere)
    nArgs = length(varargin);
    if rem(nArgs, 2) ~= 0
        error('NameValuePairToStruct:NotNameValuePairs', ...
            'Inputs were not name/value pairs');
    end

    % Sanity check: if defaults is not supplied, it's by default an empty struct
    if ~exist('defaults', 'var')
        defaults = struct;
    end
    if ~exist('restrict_flag', 'var')
        restrict_flag = false;
    end

    % Syntactic sugar: if defaults is also a cell instead of a struct, we convert it on-the-fly
    if iscell(defaults)
        defaults = struct(defaults{:});
    end

    optionNames = fieldnames(defaults); % extract all default arguments names (useful for restrict_flag)

    argStruct = defaults; % copy over the defaults: by default, all arguments will have the default value.After we will simply overwrite the defaults with the user specified values.
    for i = 1:2:nArgs % iterate over couples of argument/value
        varname = varargin{i};
        % check that the supplied name is a valid variable identifier (it does not check if the variable is allowed/declared in defaults, just that it's a possible variable name!)
        if ~isvarname(varname)
          error('NameValuePairToStruct:InvalidName', ...
             'A variable name was not valid: %s position %i', varname, i);
        % if options are restricted, check that the argument's name exists in the supplied defaults, else we throw an error. With this we can allow only a restricted range of arguments by specifying in the defaults.
        elseif restrict_flag && ~isempty(defaults) && ~any(strmatch(varname, optionNames))
            if restrict_flag ~= 2 % restrict_flag = 2 means that we just ignore this argument, else we show an error
                error('%s is not a recognized argument name', varname);
            end
        % else alright, we replace the default value for this argument with the user supplied one (or we create the variable if it wasn't in the defaults and there's no restrict_flag)
        else
            argStruct = setfield(argStruct, varname, varargin{i + 1});  %#ok<SFLD>
        end
    end

end

function varspull(s)
% Import variables in a structures into the local namespace/workspace
% eg: s = struct('foo', 1, 'bar', 'qwerty'); varspull(s); disp(foo); disp(bar);
% Will print: 1 and qwerty
%
%
% Author: Jason S
%
    for n = fieldnames(s)'
        name = n{1};
        value = s.(name);
        assignin('caller',name,value);
    end
end

function varargin = delarg(varname, varargin)
% varargin = delarg(varname, varargin)
% Removes an argument from varargin with name varname (varname must be either a string or a cell array of strings)

    % Extract the arguments if it's inside a sub-struct (happens on Octave), because anyway it's impossible that the number of argument be 1 (you need at least a couple, thus two)
    if (numel(varargin) == 1)
        varargin = varargin{:};
    end

    nArgs = length(varargin);
    for i = nArgs-1:-2:1 % iterate over couples of argument/value
        vname = varargin{i};
        if iscell(varname) && any(ismember(varname, vname))
            varargin(i:i+1) = [];
        elseif strcmp(vname, varname)
            varargin(i:i+1) = [];
            break;
        end
    end
end

function varargin = editarg(varname, varvalue, varargin)
% varargin = editarg(varname, varargin)
% Replaces an argument from varargin with name varname (varname must either be a string or a cell array of strings, same for varvalue) with the content varvalue

    % Extract the arguments if it's inside a sub-struct (happens on Octave), because anyway it's impossible that the number of argument be 1 (you need at least a couple, thus two)
    if (numel(varargin) == 1)
        varargin = varargin{:};
    end

    nArgs = length(varargin);
    for i = 1:2:nArgs % iterate over couples of argument/value
        vname = varargin{i};
        if iscell(varname) && any(ismember(varname, vname))
            idx = find(ismember(varname, vname));
            varargin{i+1} = varvalue{idx};
        elseif strcmp(vname, varname)
            varargin{i+1} = varvalue;
            break;
        end
    end
end

function varargin = addarg(varname, varvalue, varargin)
    varargin = {varargin ; varname ; varvalue};
end

function printcputime(perf, sometext)
    if ~exist('sometext', 'var') || isempty(sometext)
        fprintf('Elapsed cpu time is %g seconds.\n', perf);
    else
        fprintf(sometext, perf);
    end
end

function printtime(perf, sometext)
    if ~exist('sometext', 'var') || isempty(sometext)
        fprintf('Elapsed time is %g seconds.\n', perf);
    else
        fprintf(sometext, perf);
    end
end

function printeta(current_it, total_it, starttime, silent)
% Print a progress meter with an estimation of the remaining time required to perform all iterations of the experiment
% ETA is calculated simply by using a linear interpolation (ie, rule of 3).
% The only arguments needed are: the current number of iterations done, the total number of iterations to do to end the experiment, and the timestamp before the experiment started (using cputime()).
    elapsed = cputime() - starttime;
    if current_it == total_it
        eta = 0;
    else
        eta = elapsed / current_it * (total_it-current_it);
    end
    if silent; fprintf('\r'); end;
    fprintf('%i/%i [elapsed: %s, eta: %s, %.2f it/s]', current_it, total_it, sec2hms(elapsed), sec2hms(eta), (current_it / elapsed)); flushout();
    if ~silent; fprintf('\n'); aux.flushout(); end;
end
