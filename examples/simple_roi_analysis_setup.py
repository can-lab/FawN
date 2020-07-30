"""Simple ROI analysis setup.

This is an example of setting up a simple ROI analysis for two regions.
The resulting images (and dummy masks) contain two voxels (each representing
the timecourse of one ROI) and can be used as inputs for first and session
level analyses (with `autoregression_smoothing=None`)

"""


import os
from glob import glob

from nipype.pipeline import engine as pe
from nipype.interfaces import io

import fawn


FUNC_DIR = "/path/to/data"
SEQUENCE_NAME = "MYSEQUENCE"
PREPROC = "5mm100sNone"

wf = pe.Workflow(name="roi_analysis")

# Resample ROI masks into same size as data (must already be in same space!)
resample = fawn.create_resampling_workflow(thresholding=0.5, binarise=True)
resample.inputs.inputspec.in_files = ["ROI1.nii.gz",
                                      "ROI2.nii.gz"]
resample.inputs.inputspec.reference = sorted(glob(os.path.join(
    FUNC_DIR, "*{0}*desc-preproc{1}*_bold.nii.gz".format(SEQUENCE_NAME,
                                                         PREPROC))))[0]

# Extract ROI timecourses
timecourses = fawn.create_timecourse_extraction_workflow(method="mean")
timecourses.inputs.inputspec.in_files = sorted(glob(os.path.join(
    FUNC_DIR, "*{0}*desc-preproc{1}*_bold.nii.gz".format(SEQUENCE_NAME,
                                                         PREPROC))))

# Save results
results = pe.Node(io.DataSink(parameterization=False), name='results')
results.inputs.base_directory = os.path.join(wf.base_dir, "results")

# Connect and run
wf.connect(resample, "outputspec.out_files", timecourses, "inputspec.in_masks")
wf.connect(timecourses, "outputspec.out_files", results, "timecourses.@files")
wf.connect(timecourses, "outputspec.out_masks", results, "timecourses.@masks")
wf.run('MultiProc', plugin_args={'n_procs': 8})
