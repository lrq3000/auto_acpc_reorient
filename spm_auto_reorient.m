function auto_reorient(p) 
spmDir=which('spm');
spmDir=spmDir(1:end-5);
tmpl=[spmDir 'canonical\avg152T1.nii'];
vg=spm_vol(tmpl);
flags.regtype='rigid';
%p=spm_select(inf,'image');
for i=1:size(p,1)
    f=strtrim(p(i,:));
    spm_smooth(f,'temp.nii',[12 12 12]);
    vf=spm_vol('temp.nii');
    [M,scal] = spm_affreg(vg,vf,flags);
    M3=M(1:3,1:3);
    [u s v]=svd(M3);
    M3=u*v';
    M(1:3,1:3)=M3;
    N=nifti(f);
    N.mat=M*N.mat;
    create(N);
end

