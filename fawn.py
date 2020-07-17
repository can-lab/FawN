"""FawN.

FSL analysis with NiPype.

"""

__author__ = "Florian Krause <f.krause@donders.ru.nl>
__version__ = "0.1.0"
__date__ = "2020-07-17"


from nipype.pipeline import engine as pe
from nipype.algorithms.modelgen import SpecifyModel
from nipype.interfaces import fsl, utility, io
from nipype.interfaces.fsl.base import FSLCommand, FSLCommandInputSpec
from nipype.interfaces.base import traits, TraitedSpec


class _PtoZInputSpec(FSLCommandInputSpec):
    pvalue = traits.Float(
        0.05, argstr='%f', usedefault=True, position=0,
        desc='p-value to compute z-statistic of (default 0.05)')
    two_tailed = traits.Bool(argstr='-2', position=1,
                             desc='use 2-tailed conversion (default 1-tailed)')
    resels = traits.Float(
        argstr='-g %f', position=2,
        desc='number of resels for GRF maximum-height theory')


class _PtoZOutputSpec(TraitedSpec):
    zstat = traits.Float(
        desc='z-statistic corresponding to specified p-value')


class _PtoZ(FSLCommand):
    """Determine the z-value corresponding to an observed p-value."""
    input_spec = _PtoZInputSpec
    output_spec = _PtoZOutputSpec
    _cmd = 'ptoz'

    def aggregate_outputs(self, runtime=None, needed_outputs=None):
        outputs = self._outputs()
        outputs.zstat = float(runtime.stdout.strip())
        return outputs


def create_first_level_workflow(tr,
                                unit="scans",
                                high_pass_filter_cutoff=128,
                                normalized_timecourse_mean=None,
                                bases={'dgamma': {'derivs': False}},
                                f_contrasts=False,
                                mem_gb=8,
                                name='first_level'):

    """Create a first level analysis workflow.

    Parameters
    ----------
    tr : float
        the repetition time
    unit : str, optional
        the unit of the data points (one of "scans" or "secs"; default="scans")
    high_pass_filter_cutoff : int, optional
        the length of the high pass filter in seconds (default=128)
    normalized_timecourse_mean : float, optional
        the mean of the normalized timecourse if applicable (default=None)
    bases : dict, optional
        the basis functions (default={'dgamma': {'derivs': False}})
    f_contrasts : bool, optional
        whether f-contrasts are specified to be estimated (default=False)
    mem_gb : int, optional
        the amount of memory in GB to be allocated to the model fittin nodes
        (default=8)
    name : str
        the name of the workflow (default="first_level")

    Inputs
    ------
    in_files : list
        the functional images of the runs
    in_masks : list
        the functional images of the brain masks for each run
    models : list
        the models for each functional run as dictionaries with the fields:
            conditions       - list of names
            onsets           - list of onsets corresponding to conditions
            durations        - list of durations corresponding to conditions
            regressors_names - list of names
            regressors       - list of values for each regressor
                               (will be automatically z-transformed)
    contrasts : list
        the contrasts to be tested in the format:
            ['name', 'T',['reg1', 'reg2'],[1, -1]]

    Ouputs
    ------
    cope_images : list
        the contrast estimates images
    dof_file : str
        the temporal degrees of freedom images for each contrast
    pe_images : list
        the parameter estimates images
    stats_dirs : list
        the stats directories
    tstat_images : list
        the t-statistic images for each contrast
    varcope_images : list
        the variance estimates images for each contrast
    zstat_images : list
        the z-statistic images for each contrast

    Returns
    -------
    wf : `nipype.Workflow` object

    """

    wf = pe.Workflow(name=name)

    def make_bunch(dicts):
        from scipy.stats import zscore
        from nipype.interfaces.base.support import Bunch
        # Z-transform confound regressors
        for d in dicts:
            d["regressors"] = [zscore(x) for x in d["regressors"]]
        return [Bunch(d) for d in dicts if type(d) is not Bunch]

    def make_maskdata_opstrings(in_masks):
        return ['-mas {0}'.format(m) for m in in_masks]

    def make_medianval_opstrings(in_masks):'copes',
                                                           'fstats',
                                                           'mask',
                                                           'mrefvars',
                                                           'pes',
                                                           'res4d',
                                                           'stats_dir',
                                                           'tdof',
                                                           'tstats',
                                                           'var_copes',
                                                           'weights',
                                                           'zfstats',
                                                           'zstats']
        return ['-k {0} -p 50'.format(m) for m in in_masks]

    def get_intnorm_opstrings(medianvals):
        return ['-mul {:.10f}'.format(10000. / x) for x in medianvals]

    inputspec = pe.Node(utility.IdentityInterface(fields=['in_files',
                                                          'in_masks',
                                                          'models',
                                                          'contrasts']),
                        name='inputspec')

    makefloat = pe.MapNode(fsl.ChangeDataType(output_datatype="float"),
                           iterfield=['in_file'],
                           name='makefloat')

    maskdata_opstrings = pe.Node(utility.Function(
        input_names=["in_masks"],
        output_names=["opstrings"],
        function=make_maskdata_opstrings),
                                 name="maskdata_opstrings")

    maskdata = pe.MapNode(fsl.ImageMaths(out_data_type='float',
                                         suffix='_masked'),
                          iterfield=['in_file', 'op_string'],
                          name='maskdata')

    medianval_opstrings = pe.Node(utility.Function(
        input_names=["in_masks"],
        output_names=["opstrings"],
        function=make_medianval_opstrings),
                                  name="medianval_opstrings")

    medianval = pe.MapNode(fsl.ImageStats(),
                           iterfield=['in_file', 'op_string'],
                           name='medianval')

    intnorm = pe.MapNode(fsl.ImageMaths(suffix='_intnorm'),
                         iterfield=['in_file', 'op_string'],
                         name='intnorm')

    dicts_to_bunches = pe.Node(utility.Function(input_names=["dicts"],
                                                output_names=["bunches"],
                                                function=make_bunch),
                               name="dict_to_bunch")

    l1_spec = pe.Node(SpecifyModel(
        parameter_source='FSL',
        time_repetition=tr,
        input_units=unit,
        high_pass_filter_cutoff=high_pass_filter_cutoff),
                      name='l1_spec')

    l1_model = pe.Node(fsl.Level1Design(interscan_interval=tr,
                                        bases=bases,
                                        model_serial_correlations=True),
                       name='l1_model')

    if not normalized_timecourse_mean:
        threshold = 1000
    else:
        threshold = normalized_timecourse_mean / 10
    l1_gen = pe.MapNode(fsl.FEATModel(),
                        iterfield=['fsf_file', 'ev_files'],
                        name='l1_gen')

    if f_contrasts:
        iterfield = ['design_file', 'in_file', 'tcon_file', 'fcon_file']
    else:
        iterfield = ['design_file', 'in_file', 'tcon_file']
    l1_fit = pe.MapNode(fsl.FILMGLS(threshold=threshold,
                                    results_dir='stats',
                                    smooth_autocorr=True,
                                    mask_size=5),
                        iterfield=iterfield,
                        mem_gb=mem_gb,
                        name='l1_fit')

    results_select = pe.MapNode(io.SelectFiles({'cope': 'cope*.nii.gz',
                                                'pe': 'pe*.nii.gz',
                                                'tstat': 'tstat*.nii.gz',
                                                'varcope': 'varcope*.nii.gz',
                                                'zstat': 'zstat*.nii.gz'},
                                               force_lists=True,
                                               sort_filelist=True),
                                iterfield=["base_directory"],
                                name='results_select')

    outputspec = pe.Node(utility.IdentityInterface(fields=['cope_images',
                                                           'dof_file',
                                                           'pe_images',
                                                           'stats_dirs',
                                                           'tstat_images',
                                                           'varcope_images',
                                                           'zstat_images']),
                         name='outputspec')

    if normalized_timecourse_mean is None:
        wf.connect(inputspec, 'in_files', medianval, 'in_file')
        wf.connect(inputspec, 'in_masks', medianval_opstrings, 'in_masks')
        #wf.connect(inputspec, 'in_files', maskdata, 'in_file')
        #wf.connect(inputspec, 'in_masks', maskdata_opstrings, 'in_masks')
        #wf.connect(maskdata_opstrings, 'opstrings', maskdata, 'op_string')
        #wf.connect(maskdata, 'out_file', medianval, 'in_file')
        #wf.connect(inputspec, 'in_masks', medianval_opstrings, 'in_masks')
        wf.connect(medianval_opstrings, 'opstrings', medianval, 'op_string')
        wf.connect(inputspec, 'in_files', intnorm, 'in_file')
        #wf.connect(maskdata, 'out_file', intnorm, 'in_file')
        wf.connect(medianval, ('out_stat', get_intnorm_opstrings),
                   intnorm, 'op_string')
        wf.connect(intnorm, 'out_file', l1_spec, 'functional_runs')
        wf.connect(intnorm, 'out_file', l1_fit, 'in_file')
    else:
        wf.connect(inputspec, 'in_files', l1_spec, 'functional_runs')
        wf.connect(inputspec, 'in_files', l1_fit, 'in_file')
        #wf.connect(inputspec, 'in_files', maskdata, 'in_file')
        #wf.connect(inputspec, 'in_masks', maskdata_opstrings, 'in_masks')
        #wf.connect(maskdata_opstrings, 'opstrings', maskdata, 'op_string')
        #wf.connect(maskdata, 'out_file', l1_spec, 'functional_runs')
        #wf.connect(maskdata, 'out_file', l1_fit, 'in_file')
    wf.connect(inputspec, 'models' , dicts_to_bunches, 'dicts')
    wf.connect(dicts_to_bunches, 'bunches', l1_spec, 'subject_info')
    wf.connect(inputspec, 'contrasts', l1_model, 'contrasts')
    wf.connect(l1_spec, 'session_info', l1_model, 'session_info')
    wf.connect(l1_model, 'fsf_files', l1_gen, 'fsf_file')
    wf.connect(l1_model, 'ev_files', l1_gen, 'ev_files')
    wf.connect(l1_gen, 'design_file', l1_fit, 'design_file')
    wf.connect(l1_gen, 'con_file', l1_fit, 'tcon_file')
    wf.connect(l1_gen, 'fcon_file', l1_fit, 'fcon_file')
    wf.connect(l1_fit, 'results_dir', results_select, 'base_directory')
    wf.connect(results_select, 'cope', outputspec, 'cope_images')
    wf.connect(results_select, 'pe', outputspec, 'pe_images')
    wf.connect(l1_fit, 'results_dir', outputspec, 'stats_dirs')
    wf.connect(results_select, 'tstat', outputspec, 'tstat_images')
    wf.connect(results_select, 'varcope', outputspec, 'varcope_images')
    wf.connect(results_select, 'zstat', outputspec, 'zstat_images')
    wf.connect(l1_fit, 'dof_file', outputspec, 'dof_file')

    return wf

def _create_fixed_effects_workflow():
    """Combine runs of a session."""

    wf = pe.Workflow(name='fixed_effects')

    copemerge = pe.MapNode(fsl.Merge(dimension='t'),
                           iterfield=['in_files'],
                           name="copemerge")

    varcopemerge = pe.MapNode(fsl.Merge(dimension='t'),
                              iterfield=['in_files'],
                              name="varcopemerge")

    level2model = pe.Node(fsl.L2Model(), name='l2model')

    flameo = pe.MapNode(fsl.FLAMEO(run_mode='fe'),
                        iterfield=['cope_file', 'var_cope_file'],
                        name="flameo")

    def get_dofvolumes(dof_files, cope_files):
        import os
        import nibabel as nb
        import numpy as np
        img = nb.load(cope_files[0])
        if len(img.shape) > 3:
            out_data = np.zeros(img.shape)
        else:
            out_data = np.zeros(list(img.shape) + [1])
        for i in range(out_data.shape[-1]):
            dof = np.loadtxt(dof_files[i])
            out_data[:, :, :, i] = dof
        filename = os.path.join(os.getcwd(), 'dof_file.nii.gz')
        newimg = nb.Nifti1Image(out_data, None, img.header)
        newimg.to_filename(filename)
        return filename

    gendof = pe.Node(
        utility.Function(
            input_names=['dof_files', 'cope_files'],
            output_names=['dof_volume'],
            function=get_dofvolumes),
        name='gendofvolume')

    wf.connect(copemerge, 'merged_file', flameo, 'cope_file')
    wf.connect(varcopemerge, 'merged_file', flameo, 'var_cope_file')
    wf.connect(copemerge, 'merged_file', gendof, 'cope_files')
    wf.connect(level2model, 'design_mat', flameo, 'design_file')
    wf.connect(level2model, 'design_con', flameo, 't_con_file')
    wf.connect(level2model, 'design_grp', flameo, 'cov_split_file')

    return wf

def create_session_level_workflow(tr,
                                  unit="scans",
                                  high_pass_filter_cutoff=128,
                                  normalized_timecourse_mean=None,
                                  bases={'dgamma': {'derivs': False}},
                                  mem_gb=8,
                                  name='session_level'):
    """Create a session level analysis workflow.

    This is a convenience worklfow for doing single run (first level) analyses
    followed by averaging over all runs ("intermediate" level in FSL).

    Parameters
    ----------
    tr : float
        the repetition time
    unit : str, optional
        the unit of the data points (one of "scans" or "secs"; default="scans")
    high_pass_filter_cutoff : int, optional
        the length of the high pass filter in seconds (default=128)
    normalized_timecourse_mean : float, optional
        the mean of the normalized timecourse if applicable (default=None)
    bases : dict, optional
        the basis functions (default={'dgamma': {'derivs': False}})
    mem_gb : int, optional
        the amount of memory in GB to be allocated to the model fittin nodes
        (default=8)
    name : str
        the name of the worklfow (default="session_level")

    Inputs
    ------
    in_files : list
        the functional images of the runs
    in_masks : list
        the functional images of the brain masks for each run
    models : list
        the models for each functional run as dictionaries with the fields:
            conditions       - list of names
            onsets           - list of onsets corresponding to conditions
            durations        - list of durations corresponding to conditions
            regressors_names - list of names
            regressors       - list of values for each regressor
                               (will be automatically z-transformed)
    contrasts : list
        the contrasts to be tested in the format:
            ['name', 'T',['reg1', 'reg2'],[1, -1]]

    Ouputs
    ------
    copes : list
        the contrast estimates images
    pes : list
        the parameter estimates images
    res4d : list
        the model fit residual mean-squared error images
    stats_dir : str
        the stats directory
    tdof : list
        the temporal degrees of freedom images for each contrast
    tstats : list
        the t-statistic images for each contrast
    var_copes : list
        the variance estimates images for each contrast
    weights : list
        the weights images for each contrast
    zstats : list
        the z-statistic images for each contrast

    Returns
    -------
    wf : `nipype.Workflow` object

    """

    #TODO: f-contrasts

    wf = pe.Workflow(name=name)

    def make_add_masks_inputs(mask_files):
        first_file = mask_files[0]
        remaining_files = mask_files[1:]
        op_string = " ".join(["-add %s" for f in mask_files[1:]]) + " -bin"
        return first_file, remaining_files, op_string

    def sort_copes(files):
        numelements = len(files[0])
        outfiles = []
        for i in range(numelements):
            outfiles.insert(i, [])
            for j, elements in enumerate(files):
                outfiles[i].append(elements[i])
        return outfiles

    def num_copes(files):
        return len(files)

    inputspec = pe.Node(utility.IdentityInterface(fields=["in_files",
                                                          "in_masks",
                                                          "models",
                                                          "contrasts"]),
                        name="inputspec")

    first_level = create_first_level_workflow(
        tr,
        unit=unit,
        high_pass_filter_cutoff=high_pass_filter_cutoff,
        normalized_timecourse_mean=normalized_timecourse_mean,
        bases=bases,
        mem_gb=mem_gb)

    add_masks_inputs = pe.Node(utility.Function(
        input_names=["mask_files"],
        output_names=["first_file",
                      "remaining_files",
                      "op_string"],
        function=make_add_masks_inputs),
                               name="add_masks_inputs")

    add_masks = pe.Node(fsl.MultiImageMaths(), name="add_masks")

    fixed_effects = _create_fixed_effects_workflow()

    outputspec = pe.Node(utility.IdentityInterface(fields=['copes',
                                                           'mrefvars',
                                                           'pes',
                                                           'res4d',
                                                           'stats_dir',
                                                           'tdof',
                                                           'tstats',
                                                           'var_copes',
                                                           'weights',
                                                           'zstats']),
                         name='outputspec')

    wf.connect(inputspec, "in_files", first_level, "inputspec.in_files")
    wf.connect(inputspec, "in_masks", first_level, "inputspec.in_masks")
    wf.connect(inputspec, "models", first_level, "inputspec.models")
    wf.connect(inputspec, "contrasts", first_level, "inputspec.contrasts")
    wf.connect(first_level, ("outputspec.cope_images", sort_copes),
               fixed_effects, "copemerge.in_files")
    wf.connect(first_level, ("outputspec.varcope_images", sort_copes),
               fixed_effects, "varcopemerge.in_files")
    wf.connect(first_level, ("outputspec.cope_images", num_copes),
               fixed_effects, "l2model.num_copes")
    wf.connect(first_level, 'outputspec.dof_file',
               fixed_effects, 'gendofvolume.dof_files')
    wf.connect(inputspec, "in_masks", add_masks_inputs, "mask_files")
    wf.connect(add_masks_inputs, "first_file", add_masks, "in_file")
    wf.connect(add_masks_inputs, "remaining_files", add_masks, "operand_files")
    wf.connect(add_masks_inputs, "op_string", add_masks, "op_string")
    wf.connect(add_masks, "out_file", fixed_effects, "flameo.mask_file")
    wf.connect(fixed_effects, 'flameo.copes', outputspec, 'copes')
    wf.connect(fixed_effects, 'flameo.mrefvars', outputspec, 'mrefvars')
    wf.connect(fixed_effects, 'flameo.pes', outputspec, 'pes')
    wf.connect(fixed_effects, 'flameo.res4d', outputspec, 'res4d')
    wf.connect(fixed_effects, 'flameo.stats_dir', outputspec, 'stats_dir')
    wf.connect(fixed_effects, 'flameo.tdof', outputspec, 'tdof')
    wf.connect(fixed_effects, 'flameo.tstats', outputspec, 'tstats')
    wf.connect(fixed_effects, 'flameo.var_copes', outputspec, 'var_copes')
    wf.connect(fixed_effects, 'flameo.weights', outputspec, 'weights')
    wf.connect(fixed_effects, 'flameo.zstats', outputspec, 'zstats')

    return wf

def create_higher_level_workflow(mode="flame1", name="higher_level"):
    """Create a higher level analysis workflow.

    Parameters
    ----------
    mode : str, optional
        the mode to run in (one of "fe", "ols", "flame1", flame12";
        "default="flame1")
    name : str
        the name of the workflow (default="higher_level")

    Inputs
    ------
    in_copes : list
        the contrast estimate images
    var_copes : list
        the variance estimate images for each contrast
    models : list
        the models for each functional run as dictionaries with the fields:
            regressors_names - list of names
            regressors       - list of values for each regressor
    contrasts : list
        the contrasts to be tested in the format:
            ['name', 'T',['reg1', 'reg2'],[1, -1]]

    Ouputs
    ------
    copes : list
        the contrast estimates images
    mask : str
        the functional mask image used for the analysis
    mrefvars : list
        the mean random effect variances images for each contrast
    pes : list
        the parameter estimates images
    res4d : list
        the model fit residual mean-squared error images
    stats_dir : str
        the stats directory
    tdof : list
        the temporal degrees of freedom images for each contrast
    tstats : list
        the t-statistic images for each contrast
    var_copes : list
        the variance estimates images for each contrast
    weights : list
        the weights images for each contrast
    zfstats : list
        the z-statistic images for each f-contrast
    zstats : list
        the z-statistic images for each contrast

    Returns
    -------
    wf : `nipype.Workflow` object

    """

    wf = pe.Workflow(name=name)

    def make_mul_masks_inputs(mask_files):
        first_file = mask_files[0]
        remaining_files = mask_files[1:]
        op_string = " ".join(["-add %s" for f in mask_files[1:]]) + " -bin"
        return first_file, remaining_files, op_string

    def get_regressors(model):
        return {x: model["regressors"][c] \
                for c, x in enumerate(model["regressor_names"])}

    def get_groups(model):
        if "groups" in model:
            return model["groups"]
        else:
            return [1] * len(model["regressors"][0])

    def make_list(item):
        if type(item) == str:
            item = [item]
        return item

    inputspec = pe.Node(utility.IdentityInterface(fields=['in_copes',
                                                          'in_varcopes',
                                                          'model',
                                                          'contrasts']),
                        name='inputspec')

    mask_varcopes = pe.MapNode(fsl.UnaryMaths(operation="bin"),
                               iterfield=['in_file'],
                               name='mask_varcopes')

    mul_masks_inputs = pe.Node(utility.Function(
        input_names=["mask_files"],
        output_names=["first_file",
                      "remaining_files",
                      "op_string"],
        function=make_mul_masks_inputs),
                               name="mul_masks_inputs")

    mul_masks = pe.Node(fsl.MultiImageMaths(), name="mul_masks")

    merge_copes = pe.Node(fsl.Merge(dimension='t'), name='merge_copes')
    merge_varcopes = pe.Node(fsl.Merge(dimension='t'), name='merge_varcopes')

    multreg_model = pe.Node(fsl.MultipleRegressDesign(), name='multreg_model')
    flameo = pe.Node(fsl.FLAMEO(run_mode=mode), name='flameo')

    outputspec = pe.Node(utility.IdentityInterface(fields=['copes',
                                                           'fstats',
                                                           'mask',
                                                           'mrefvars',
                                                           'pes',
                                                           'res4d',
                                                           'stats_dir',
                                                           'tdof',
                                                           'tstats',
                                                           'var_copes',
                                                           'weights',
                                                           'zfstats',
                                                           'zstats']),
                        name='outputspec')

    wf.connect(inputspec, 'in_copes', merge_copes, 'in_files')
    wf.connect(inputspec, 'in_varcopes', merge_varcopes, 'in_files')
    wf.connect(inputspec, 'in_varcopes', mask_varcopes, 'in_file')
    wf.connect(mask_varcopes, 'out_file', mul_masks_inputs, 'mask_files')
    wf.connect(mul_masks_inputs, "first_file", mul_masks, "in_file")
    wf.connect(mul_masks_inputs, "remaining_files", mul_masks, "operand_files")
    wf.connect(mul_masks_inputs, "op_string", mul_masks, "op_string")
    wf.connect(merge_copes, 'merged_file', flameo, 'cope_file')
    wf.connect(merge_varcopes, 'merged_file', flameo, 'var_cope_file')
    wf.connect(mul_masks, 'out_file', flameo, 'mask_file')
    wf.connect(inputspec, ('model', get_regressors),
               multreg_model, 'regressors')
    wf.connect(inputspec, ('model', get_groups), multreg_model, 'groups')
    wf.connect(inputspec, 'contrasts', multreg_model, 'contrasts')
    wf.connect(multreg_model, 'design_mat', flameo, 'design_file')
    wf.connect(multreg_model, 'design_con', flameo, 't_con_file')
    wf.connect(multreg_model, 'design_fts', flameo, 'f_con_file')
    wf.connect(multreg_model, 'design_grp', flameo, 'cov_split_file')
    wf.connect(flameo, ('copes', make_list), outputspec, 'copes')
    wf.connect(flameo, 'fstats', outputspec, 'fstats')
    wf.connect(mul_masks, 'out_file', outputspec, 'mask')
    wf.connect(flameo, 'mrefvars', outputspec, 'mrefvars')
    wf.connect(flameo, ('pes', make_list), outputspec, 'pes')
    wf.connect(flameo, 'res4d', outputspec, 'res4d')
    wf.connect(flameo, 'stats_dir', outputspec, 'stats_dir')
    wf.connect(flameo, ('tdof', make_list), outputspec, 'tdof')
    wf.connect(flameo, ('tstats', make_list), outputspec, 'tstats')
    wf.connect(flameo, ('var_copes', make_list), outputspec, 'var_copes')
    wf.connect(flameo, 'weights', outputspec, 'weights')
    wf.connect(flameo, ('zfstats', make_list), outputspec, 'zfstats')
    wf.connect(flameo, ('zstats', make_list), outputspec, 'zstats')

    return wf

def create_thresholding_workflow(pvalue=0.05, two_tailed=True,
                                 cluster_connectivity=26,
                                 cluster_threshold=3.2,
                                 name="thresholding"):
    """Create a thresholding workflow.

    Thresholds at voxel-level (FWE corrected) and cluster-level.

    Parameters
    ----------
    pvalue : float, optional
        the pvalue to threshold at (default=0.05)
    two-tailed : bool, optional
        whether to do two-tailed tests (default=True)
    cluster_connectivity : int, optional
        the cluster connectivity value (one of 6, 18, 26; default=26)
    cluster_threshold : float, optional
        the initial cluster threshold (if None, use single voxel FEW threshold;
        default=3.2)
    name : str
        the name of the workflow (default="thresholding")

    Inputs
    ------
    in_zstats : list
        the z-statistic images
    in_copes : list
        the estimate contrast images
    in_res4d : list
        the model fit residual mean-squared error images
    in_mask : str
        the mask image used to mask the analysis

    Ouputs
    ------
    fwe_zstats : str
        the FWE corrected voxel-level thresholded z-statistic image
    cluster_zstats : str
        the cluster-level thresholded z-statistic image
    cluster_pos_idx : str
        the positive cluster index image
    cluster_pos_max : str
        the positive cluster local maxima file
    cluster_neg_idx : str
        the negative cluster index image
    cluster_neg_max : str
        the negative negative local maixma file
    cluster_pos_max_thresh : str
        the thresholded positive cluster local maxima file
    cluster_neg_max_thresh : str
        the thresholded negative cluster local maxima file

    Returns
    -------
    wf : `nipype.Workflow` object

    """

    wf = pe.Workflow(name=name)

    def _dof(copes):
        return len(copes) - 1

    def _neg(val):
        return -val

    def _resel_count(volume, resel_size):
        return volume/resel_size

    def _threshold_localmax(localmax_file, fwe_threshold):
        import os
        with open(localmax_file) as f_in:
            with open("localmax_thresholded.txt", 'w') as f_out:
                fsl_cluster = 1000000
                new_cluster = 0
                for line in f_in.readlines():
                    tmp = line.split("\t")
                    if tmp[0].startswith("Cluster"):
                        f_out.write(line)
                    elif abs(float(tmp[1])) > abs(fwe_threshold):
                        if int(tmp[0]) < fsl_cluster:
                            fsl_cluster = int(tmp[0])
                            new_cluster += 1
                        tmp[0] = str(new_cluster)
                        f_out.write("\t".join(tmp))
        return os.path.abspath("localmax_thresholded.txt")

    inputspec = pe.MapNode(utility.IdentityInterface(fields=['in_zstats',
                                                             'in_copes',
                                                             'in_res4d',
                                                             'in_mask']),
                           iterfield=["in_zstats"],
                           name='inputspec')

    merge_copes = pe.Node(fsl.Merge(dimension='t'), name='merge_copes')

    dof = pe.Node(utility.Function(input_names=["copes"],
                                   output_names=["dof"],
                                   function=_dof),
                  name="dof")

    if two_tailed:
        pvalue /= 2

    # FWE
    smoothness = pe.Node(fsl.SmoothEstimate(), name='smoothness')
    resel_count = pe.Node(utility.Function(
        input_names=["volume", "resel_size"],
        output_names=["resel_count"],
        function=_resel_count),
                          name="resel_count")
    fwe_ptoz = pe.Node(_PtoZ(pvalue=pvalue), name='fwe_ptoz')
    fwe_nonsig0 = pe.MapNode(fsl.Threshold(direction='above'),
                             iterfield=["in_file"],
                             name='fwe_nonsig0')
    fwe_nonsig1 = pe.MapNode(fsl.Threshold(direction='below'),
                             iterfield=["in_file"],
                             name='fwe_nonsig1')
    fwe_thresh = pe.MapNode(fsl.BinaryMaths(operation='sub'),
                            iterfield=["in_file",
                                       "operand_file"],
                            name='fwe_thresh')

    # Cluster
    cluster_kwargs = {
        'connectivity': cluster_connectivity,
        'out_threshold_file': True,
        'out_index_file': True,
        'out_localmax_txt_file': True
    }
    cluster_pos = pe.MapNode(fsl.Cluster(**cluster_kwargs),
                             iterfield=["in_file"],
                             name='cluster_pos')
    if cluster_threshold is not None:
        cluster_pos.inputs.threshold = cluster_threshold
    cluster_pos.inputs.pthreshold = pvalue

    cluster_neg = pe.MapNode(fsl.Cluster(**cluster_kwargs),
                             iterfield=["in_file"],
                             name='cluster_neg')
    if cluster_threshold is not None:
        cluster_neg.inputs.threshold = cluster_threshold
    cluster_neg.inputs.pthreshold = pvalue

    zstat_inv = pe.MapNode(fsl.BinaryMaths(operation='mul', operand_value=-1),
                           iterfield=["in_file"],
                           name='zstat_inv')
    cluster_inv = pe.MapNode(
        fsl.BinaryMaths(operation='mul', operand_value=-1),
        iterfield=["in_file"],
        name='cluster_inv')
    cluster_all = pe.MapNode(fsl.BinaryMaths(operation='add'),
                             iterfield=["in_file",
                                        "operand_file"],
                             name='cluster_all')

    threshold_localmax_pos = pe.MapNode(utility.Function(
        input_names=["localmax_file",
                     "fwe_threshold"],
        output_names=["localmax_thresholded"],
        function=_threshold_localmax),
                                        iterfield=["localmax_file"],
                                        name="threshold_localmax_pos")
    threshold_localmax_neg = pe.MapNode(utility.Function(
        input_names=["localmax_file",
                     "fwe_threshold"],
        output_names=["localmax_thresholded"],
        function=_threshold_localmax),
                                        iterfield=["localmax_file"],
                                        name="threshold_localmax_neg")

    outputspec = pe.Node(utility.IdentityInterface(
        fields=['fwe_zstats',
                'cluster_zstats',
                'cluster_pos_idx',
                'cluster_pos_max',
                'cluster_neg_idx',
                'cluster_neg_max',
                'cluster_pos_max_thresh',
                'cluster_neg_max_thresh']),
                         name='outputspec')

    wf.connect(inputspec, 'in_res4d', smoothness, 'residual_fit_file')
    wf.connect(inputspec, 'in_copes', dof, 'copes')
    wf.connect(dof, 'dof', smoothness, 'dof')
    wf.connect(inputspec, 'in_mask', smoothness, 'mask_file')
    wf.connect(inputspec, 'in_zstats', fwe_nonsig0, 'in_file')
    wf.connect(fwe_nonsig0, 'out_file', fwe_nonsig1, "in_file")
    wf.connect(smoothness, 'volume', resel_count, 'volume')
    wf.connect(smoothness, 'resels', resel_count, 'resel_size')
    wf.connect(resel_count, 'resel_count', fwe_ptoz, 'resels')
    wf.connect(fwe_ptoz, 'zstat', fwe_nonsig0, 'thresh')
    wf.connect(fwe_ptoz, ('zstat', _neg), fwe_nonsig1, 'thresh')
    wf.connect(inputspec, 'in_zstats', fwe_thresh, 'in_file')
    wf.connect(fwe_nonsig1, 'out_file', fwe_thresh, 'operand_file')
    wf.connect(inputspec, 'in_zstats', cluster_pos, 'in_file')
    wf.connect(inputspec, 'in_copes', merge_copes, 'in_files')
    wf.connect(merge_copes, 'merged_file', cluster_pos, 'cope_file')
    wf.connect(smoothness, 'volume', cluster_pos, 'volume')
    wf.connect(smoothness, 'dlh', cluster_pos, 'dlh')
    wf.connect(inputspec, 'in_zstats', zstat_inv, 'in_file')
    wf.connect(zstat_inv, 'out_file', cluster_neg, 'in_file')
    wf.connect(cluster_neg, 'threshold_file', cluster_inv, 'in_file')
    wf.connect(merge_copes, 'merged_file', cluster_neg, 'cope_file')
    wf.connect(smoothness, 'volume', cluster_neg, 'volume')
    wf.connect(smoothness, 'dlh', cluster_neg, 'dlh')
    if cluster_threshold is None:
        wf.connect(fwe_ptoz, 'zstat', cluster_pos, 'threshold')
        wf.connect(fwe_ptoz, 'zstat', cluster_neg, 'threshold')
    wf.connect(cluster_pos, 'threshold_file', cluster_all, 'in_file')
    wf.connect(cluster_inv, 'out_file', cluster_all, 'operand_file')
    wf.connect(fwe_thresh, 'out_file', outputspec, 'fwe_zstats')
    wf.connect(cluster_all, 'out_file', outputspec, 'cluster_zstats')
    wf.connect(cluster_pos, 'index_file', outputspec, 'cluster_pos_idx')
    wf.connect(cluster_pos, 'localmax_txt_file', outputspec, 'cluster_pos_max')
    wf.connect(cluster_neg, 'index_file', outputspec, 'cluster_neg_idx')
    wf.connect(cluster_neg, 'localmax_txt_file', outputspec, 'cluster_neg_max')
    wf.connect(cluster_pos, 'localmax_txt_file',
               threshold_localmax_pos, 'localmax_file')
    wf.connect(fwe_ptoz, 'zstat', threshold_localmax_pos, 'fwe_threshold')
    wf.connect(threshold_localmax_pos, 'localmax_thresholded',
               outputspec, 'cluster_pos_max_thresh')
    wf.connect(cluster_neg, 'localmax_txt_file',
               threshold_localmax_neg, 'localmax_file')
    wf.connect(fwe_ptoz, 'zstat', threshold_localmax_neg, 'fwe_threshold')
    wf.connect(threshold_localmax_neg, 'localmax_thresholded',
               outputspec, 'cluster_neg_max_thresh')
    return wf
