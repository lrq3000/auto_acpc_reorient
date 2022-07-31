function funs = auto_acpc_reorient_spmaux
% funs = auto_acpc_reorient_spmaux
% Exposes private functions from SPM that are necessary to call optfun() inside spm_coreg()
%
% The content of this file should be considered the intellectual property of Wellcome Trust Centre for Neuroimaging
% with some slight modifications by Stephen Karl Larroque, Coma Science Group & GIGA-Consciousness, University Hospital of Liege, Belgium
% and licensed under GPLv2+ (GNU Public License v2 or later)

    funs.loaduint8 = @loaduint8;
    funs.smooth_uint8 = @smooth_uint8;
    funs.smoothvol2uint8 = @smoothvol2uint8;
    funs.jointhist = @jointhist;
    funs.jointhist_display = @jointhist_display;
    funs.jointhist_score = @jointhist_score;
end  % endfunction

function V = smooth_uint8(V,fwhm)
% Convolve the volume in memory (fwhm in voxels).
% From spm12/toolbox/FieldMap/FieldMap.m
% Copyright (C) 2006-2016 Wellcome Trust Centre for Neuroimaging
% Authors: Jesper Andersson and Chloe Hutton 
% Licensed under GPLv2+ (GNU Public License v2 or later)
% TODO: can be replaced with spm_smoothto8bit()?

lim = ceil(2*fwhm);
s  = fwhm/sqrt(8*log(2));
x  = [-lim(1):lim(1)]; x = smoothing_kernel(fwhm(1),x); x  = x/sum(x);
y  = [-lim(2):lim(2)]; y = smoothing_kernel(fwhm(2),y); y  = y/sum(y);
z  = [-lim(3):lim(3)]; z = smoothing_kernel(fwhm(3),z); z  = z/sum(z);
i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
spm_conv_vol(V.uint8,V.uint8,x,y,z,-[i j k]);
return;
end  % endfunction

function krn = smoothing_kernel(fwhm,x)
% From spm12/toolbox/FieldMap/FieldMap.m
% Copyright (C) 2006-2016 Wellcome Trust Centre for Neuroimaging
% Authors: Jesper Andersson and Chloe Hutton 
% Licensed under GPLv2+ (GNU Public License v2 or later)
% TODO: can be replaced with spm_smoothkern()?

% Variance from FWHM
s = (fwhm/sqrt(8*log(2)))^2+eps;

% The simple way to do it. Not good for small FWHM
% krn = (1/sqrt(2*pi*s))*exp(-(x.^2)/(2*s));

% For smoothing images, one should really convolve a Gaussian
% with a sinc function.  For smoothing histograms, the
% kernel should be a Gaussian convolved with the histogram
% basis function used. This function returns a Gaussian
% convolved with a triangular (1st degree B-spline) basis
% function.

% Gaussian convolved with 0th degree B-spline
% int(exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t= -0.5..0.5)
% w1  = 1/sqrt(2*s);
% krn = 0.5*(erf(w1*(x+0.5))-erf(w1*(x-0.5)));

% Gaussian convolved with 1st degree B-spline
%  int((1-t)*exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t= 0..1)
% +int((t+1)*exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t=-1..0)
w1  =  0.5*sqrt(2/s);
w2  = -0.5/s;
w3  = sqrt(s/2/pi);
krn = 0.5*(erf(w1*(x+1)).*(x+1) + erf(w1*(x-1)).*(x-1) - 2*erf(w1*x   ).* x)...
      +w3*(exp(w2*(x+1).^2)     + exp(w2*(x-1).^2)     - 2*exp(w2*x.^2));

krn(krn<0) = 0;
return;
end  % endfunction

function udat = loaduint8(V)
% Load data from file indicated by V into an array of unsigned bytes.
% From spm12/toolbox/FieldMap/FieldMap.m
% Copyright (C) 2006-2016 Wellcome Trust Centre for Neuroimaging
% Authors: Jesper Andersson and Chloe Hutton 
% Licensed under GPLv2+ (GNU Public License v2 or later)

if size(V.pinfo,2)==1 && V.pinfo(1) == 2
    mx = 255*V.pinfo(1) + V.pinfo(2);
    mn = V.pinfo(2);
else
    %spm_progress_bar('Init',V.dim(3),...
    %    ['Computing max/min of ' spm_file(V.fname,'filename')],...
    %    'Planes complete');
    mx = -Inf; mn =  Inf;
    for p=1:V.dim(3)
        img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
        mx  = max([max(img(:))+paccuracy(V,p) mx]);
        mn  = min([min(img(:)) mn]);
        %spm_progress_bar('Set',p);
    end
end
%spm_progress_bar('Init',V.dim(3),...
%    ['Loading ' spm_file(V.fname,'filename')],...
%    'Planes loaded');

udat = uint8(0);
udat(V.dim(1),V.dim(2),V.dim(3))=0;
rand('state',100);
for p=1:V.dim(3)
    img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
    acc = paccuracy(V,p);
    if acc==0
        udat(:,:,p) = uint8(round((img-mn)*(255/(mx-mn))));
    else
        % Add random numbers before rounding to reduce aliasing artifact
        r = rand(size(img))*acc;
        udat(:,:,p) = uint8(round((img+r-mn)*(255/(mx-mn))));
    end
    %spm_progress_bar('Set',p);
end
%spm_progress_bar('Clear');
return;
end  % endfunction

function acc = paccuracy(V,p)
% From spm12/toolbox/FieldMap/FieldMap.m
% Copyright (C) 2006-2016 Wellcome Trust Centre for Neuroimaging
% Authors: Jesper Andersson and Chloe Hutton 
% Licensed under GPLv2+ (GNU Public License v2 or later)

if ~spm_type(V.dt(1),'intt')
    acc = 0;
else
    if size(V.pinfo,2)==1
        acc = abs(V.pinfo(1,1));
    else
        acc = abs(V.pinfo(1,p));
    end
end
end  %endfunction

function Vo = smoothvol2uint8(VF, fwhmf)
% Smooth and add required uint8 field for spm_coreg
% From spm12/toolbox/FieldMap/FieldMap.m
% Copyright (C) 2006-2016 Wellcome Trust Centre for Neuroimaging
% Authors: Jesper Andersson and Chloe Hutton 
% Licensed under GPLv2+ (GNU Public License v2 or later)

    if ~isfield(VF,'uint8')
        %VF       = smooth_uint8(VF,fwhmf); % Note side effects
        Vo = spm_smoothto8bit(VF, fwhmf);
        Vo.uint8 = loaduint8(VF);
    end
end  % endfunction

function H = jointhist(template_vol, input_vol, iM, fwhm)
% Compute a joint histogram using spm_hist2
% iM is a spm_imatrix so that pinv(iM)*input_vol.mat will match template_vol

    % Default values
    if ~exist('fwhm') || isempty(fwhm)
        fwhm = [1 1];
    end

    % Convert volumes to uint8
    template_vol_uint8 = smoothvol2uint8(template_vol, fwhm(1));
    input_vol_uint8 = smoothvol2uint8(input_vol, fwhm(2));
    % Compute histogram
    H = spm_hist2(template_vol_uint8.uint8,input_vol_uint8.uint8, input_vol_uint8.mat\spm_matrix(iM(:)') * template_vol_uint8.mat, [1 1 1]);
end  %endfunction

function jointhist_display(template_vol, input_vol, x, fwhm, flags)
% Display the joint histogram of two volumes, with and without a transform x
% From spm12/spm_coreg.m , function display_results()

    % Default values
    if ~exist('fwhm') || isempty(fwhm)
        fwhm = [1 1];
    end
    if ~exist('flags') || isempty(flags)
        flags = struct('cost_fun', 'ncc');
    end

    % Convert input volumes to uint8
    VG = smoothvol2uint8(template_vol, fwhm(1));
    VF = smoothvol2uint8(input_vol, fwhm(2));

    % Launch SPM window (necessary to show figures)
    spm('fmri');

    % Get control of the side window
    fig = spm_figure('FindWin','Graphics');
    if isempty(fig), return; end;
    set(0,'CurrentFigure',fig);
    spm_figure('Clear','Graphics');

    %txt = 'Information Theoretic Coregistration';
    switch lower(flags.cost_fun)
        case 'mi',  txt = 'Mutual Information Coregistration';
        case 'ecc', txt = 'Entropy Correlation Coefficient Registration';
        case 'nmi', txt = 'Normalised Mutual Information Coregistration';
        case 'ncc', txt = 'Normalised Cross Correlation';
        otherwise, error('Invalid cost function specified');
    end

    % Display text
    %--------------------------------------------------------------------------
    ax = axes('Position',[0.1 0.8 0.8 0.15],'Visible','off','Parent',fig);
    text(0.5,0.7, txt,'FontSize',16,...
        'FontWeight','Bold','HorizontalAlignment','center','Parent',ax);

    Q = inv(VF.mat\spm_matrix(x(:)')*VG.mat);
    text(0,0.5, sprintf('X1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(1,:)),'Parent',ax);
    text(0,0.3, sprintf('Y1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(2,:)),'Parent',ax);
    text(0,0.1, sprintf('Z1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(3,:)),'Parent',ax);

    % Display joint histograms
    %--------------------------------------------------------------------------
    ax  = axes('Position',[0.1 0.5 0.35 0.3],'Visible','off','Parent',fig);
    H   = spm_hist2(VG.uint8,VF.uint8,VF.mat\VG.mat,[1 1 1]);
    tmp = log(H+1);
    image(tmp*(64/max(tmp(:))),'Parent',ax');
    set(ax,'DataAspectRatio',[1 1 1],...
        'PlotBoxAspectRatioMode','auto','XDir','normal','YDir','normal',...
        'XTick',[],'YTick',[]);
    title('Original Joint Histogram','Parent',ax);
    xlabel(spm_file(VG.fname,'short22'),'Parent',ax,'Interpreter','none');
    ylabel(spm_file(VF.fname,'short22'),'Parent',ax,'Interpreter','none');

    H   = spm_hist2(VG.uint8,VF.uint8,VF.mat\spm_matrix(x(:)')*VG.mat,[1 1 1]);
    ax  = axes('Position',[0.6 0.5 0.35 0.3],'Visible','off','Parent',fig);
    tmp = log(H+1);
    image(tmp*(64/max(tmp(:))),'Parent',ax');
    set(ax,'DataAspectRatio',[1 1 1],...
        'PlotBoxAspectRatioMode','auto','XDir','normal','YDir','normal',...
        'XTick',[],'YTick',[]);
    title('Final Joint Histogram','Parent',ax);
    xlabel(spm_file(VG.fname,'short22'),'Parent',ax,'Interpreter','none');
    ylabel(spm_file(VF.fname,'short22'),'Parent',ax,'Interpreter','none');

    % Display ortho-views
    %--------------------------------------------------------------------------
    spm_orthviews('Reset');
         spm_orthviews('Image',VG,[0.01 0.01 .48 .49]);
    h2 = spm_orthviews('Image',VF,[0.51 0.01 .48 .49]);
    global st
    st.vols{h2}.premul = inv(spm_matrix(x(:)'));
    spm_orthviews('Space');

    spm_print;
end  %endfunction

function o = jointhist_score(H, cf)
% Compute a scalar score from a joint histogram
% From spm12/spm_coreg.m inside optfun() function

    % Compute cost function from histogram
    H  = H+eps;
    sh = sum(H(:));
    H  = H/sh;
    s1 = sum(H,1);
    s2 = sum(H,2);

    switch lower(cf)
        case 'mi'
            % Mutual Information:
            H   = H.*log2(H./(s2*s1));
            mi  = sum(H(:));
            o   = -mi;
        case 'ecc'
            % Entropy Correlation Coefficient of:
            % Maes, Collignon, Vandermeulen, Marchal & Suetens (1997).
            % "Multimodality image registration by maximisation of mutual
            % information". IEEE Transactions on Medical Imaging 16(2):187-198
            H   = H.*log2(H./(s2*s1));
            mi  = sum(H(:));
            ecc = -2*mi/(sum(s1.*log2(s1))+sum(s2.*log2(s2)));
            o   = -ecc;
        case 'nmi'
            % Normalised Mutual Information of:
            % Studholme,  Hill & Hawkes (1998).
            % "A normalized entropy measure of 3-D medical image alignment".
            % in Proc. Medical Imaging 1998, vol. 3338, San Diego, CA, pp. 132-143.
            nmi = (sum(s1.*log2(s1))+sum(s2.*log2(s2)))/sum(sum(H.*log2(H)));
            o   = -nmi;
        case 'ncc'
            % Normalised Cross Correlation
            i     = 1:size(H,1);
            j     = 1:size(H,2);
            m1    = sum(s2.*i');
            m2    = sum(s1.*j);
            sig1  = sqrt(sum(s2.*(i'-m1).^2));
            sig2  = sqrt(sum(s1.*(j -m2).^2));
            [i,j] = ndgrid(i-m1,j-m2);
            ncc   = sum(sum(H.*i.*j))/(sig1*sig2);
            o     = -ncc;
        otherwise
            error('Invalid cost function specified');
    end
end  %endfunction
