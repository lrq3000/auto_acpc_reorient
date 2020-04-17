# spm_reorient_auto : Cross-platform automatic AC-PC realignment/reorientation and coregistration robust to brain damage in SPM

[![GitHub release](https://img.shields.io/github/release/lrq3000/spm_reorient_auto.svg)](https://github.com/lrq3000/spm_reorient_auto/releases/)

Cross-platform automatic AC-PC realignment/reorientation and coregistration for both healthy volunteers and brain damaged patients using template matching in SPM 12.

This is a set of routines to perform automatic reorient and coregistration with the toolbox [Statistical Parametric Mapping 12 (SPM12)](https://www.fil.ion.ucl.ac.uk/spm/).

![Automatic coregistration example using spm_coreg_auto.m](img/coreg.png)

## Description

Setting up the AC-PC and reorienting images is a recurrent issue in between-subjects group analyses, since they rely on coregistration methods that, like the "unified segmentation" of SPM12, are for most sensitive to initial conditions (the starting orientation of the image). These scripts are built on top of the original auto reorientation routines by [John Ashburner & Carlton Chu](https://en.wikibooks.org/wiki/SPM/How-to#How_to_automatically_reorient_images?).

The main function, `spm_reorient_auto.m`, automatically (but approximately) calculates a reorientation transform onto a target template in MNI space, in two steps:

1. a non-linear coregistration of the input image onto a target template in MNI space is calculated using `spm_affreg`,
2. then another transform is calculated using Mutual Information on a joint histogram (spm_coreg), and then applies only the rigid-body transform part of both coregistrations to reorient the input image. This allows to mainly set the origin on the AC and correct for head rotation, in order to further proceed with the segmentation/normalisation of the image.

This whole reorientation scheme relies on the "template matching" principle (as in the old normalize function), you therefore need to specify the appropriate template/reference image (we provide one, `t1group`, by default).

In any case, it is advised to check the automatically reoriented images afterwards, and [fix the orientation manually, using SPM -> Display](https://en.wikibooks.org/wiki/SPM/How-to#How_to_manually_change_the_orientation_of_an_image?) if necessary.

Another function, `spm_coreg_auto.m`, expands on the same ideas to allow coregistration between modalities (eg, functional BOLD on structural MPRAGE). It is advised that `spm_reorient_auto()` to be first applied on the structural before applying `spm_coreg_auto()` on the other modality (even if you do manually fix the reorientation, as this ensures that the T1 is somewhat in the MNI space, making it easier for `spm_coreg_auto()` to find the correct translation matrix).

## Install

To install this tool :

* copy `spm_reorient_auto.m` and `spm_coreg_auto.m` in your `spm` folder. This will allow the command `spm_reorient_auto()` and `spm_coreg_auto()` to be called from command-line (if no argument is given, a file selector dialog will open).
* copy `T1_template_CAT12_rm_withskull.nii` to your `spm/canonical` folder. This is a template generated on 10 subjects using CAT12 that were manually reoriented to AC-PC and averaged, this provides better performance for reorientation than the more blurry MNI template. Note that this template is slightly better aligned to the AC-PC plane than the original MNI template, so that there may be a slight rotation bias compared to MNI if you use this custom template (usually it's mostly unnoticeable and this should have no influence if afterwards you do the SPM normalization on MNI on your data).
* Add `spm12` and `spm12\toolbox\OldNorm` to the path in MATLAB.

## Usage

Type `help spm_reorient_auto`, for all the details and various options for reorientation of a T1.

Type `help spm_coreg_auto` for coregistering options.

Both scripts allows to use SPM filedialogs GUI, by simply typing `spm_reorient_auto()` or `spm_coreg_auto` in the MATLAB prompt.

General guideline:

* If you want to reorient isotropic T1, use `spm_reorient_auto`.
* If you want to reorient another modality (usually with less resolution), or an anisotropic T1, or an isotropic T1 but on your own custom T1 template, use `spm_coreg_auto`.

Note that the scripts can't be used from SPM12 GUI nor the BATCH system.

## Guarantee

There is no guarantee that this will work 100% of the times, although it was observed to produce good results with our own data (young and old healthy subjects, AD/PD patients, most of brain damaged patients even with significant movement or metal artefacts).

The best results we got were by doing the following steps:

1. Call `spm_reorient_auto()` on the structural in MNI space, so that it is also matching other SPM templates
2. Manually review and fix the misoriented structural images
3. Coregister the functional to SPM12 EPI template (this allows a correct translation matrix and a good basis for rotation)
4. Coregister the functional onto the structural (this fine-tunes rotation to precisely match the subject's structural)

The last 2 steps can be done by calling `spm_coreg_auto()`, which has optimized default parameters for this task. For indication, on a dataset of 100 subjects with some heavily brain damaged or artifacted, the coregistration success rate was more than 95%.

For a comparison of various methods for AC-PC reorientation, the following article is a good read:

`Liu, Yuan, and Benoit M. Dawant. "Automatic detection of the anterior and posterior commissures on MRI scans using regression forests." 2014 36th Annual International Conference of the IEEE Engineering in Medicine and Biology Society. IEEE, 2014.`

## Citation

Please cite this work as following (a paper is in preparation but not available yet):

> spm_reorient_auto. [https://github.com/lrq3000/spm_reorient_auto](https://github.com/lrq3000/spm_reorient_auto)

## Authors

This software was developed by Stephen Karl Larroque (Coma Science Group, GIGA-Consciousness, University of Liege, Belgium).

The code of spm_reorient_auto.m is a fork from [auto_reorient.m](https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=SPM;d1f675f1.0810) written by John Ashburner (FIL, UCL, London, UK) and Carlton Chu (FIL, UCL, London, UK). Unfortunately, the original code was not released under a version control system, so the history in this repository does not include the history of the original commits. The history is hence a bit muddy, so if there is any inaccuracies, please let us know, we will happily update it.

## Similar projects

We later found that K. Nemoto made a similar enhanced algorithm back in 2017, based on centering to the origin (instead of translating to match a template for us - what we call a pre-coregistration step) before applying the original auto_reorient.m script but tweaking it to apply a coregistration on a DARTEL template (instead of SPM12 ones in the original script, or a custom made template using CAT12 for this one). The resulting script, [acpc_coreg.m](https://web.archive.org/web/20180727093129/http://www.nemotos.net/scripts/acpc_coreg.m), can be found on [nemotos.net](https://www.nemotos.net/?p=1892).

Another earlier fork of auto_reorient.m named [spm_auto_reorient](https://github.com/CyclotronResearchCentre/spm_auto_reorient) was done by Christophe Phillips (Cyclotron Research Centre, University of Liege, Belgium), which adds over the original supports for more templates, option to apply the reorientation matrix to other files, sanitization of input variables, integration in the SPM12 main GUI using cfg files and usage in the BATCH system.

## License

General Public License (GPL) v2.

## Contact

For questions or suggestions, contact Stephen Karl Larroque at: stephen dot larroque at uliege dot be.
