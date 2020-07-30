"""Simple group analysis.

This is an example of a simple group level analysis including thresholding.

"""

import os
import re
from glob import glob

from nipype.pipeline import engine as pe
from nipype.interfaces import io

import fawn


DATA_DIR = "/path/to/data"
SESSION = "ses-mri01"
N_SUBJECTS = 10
N_CONTRASTS = 3
MODE = "FLAME1"  #"FE"

# Helper function to sort contrasts correctly (releveant for N_CONTRASTS > 9)
def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in _nsre.split(s)]

# Get contrasts from all subjects
copes = [[] for _ in range(N_CONTRASTS)]
varcopes = [[] for _ in range(N_CONTRASTS)]
for subject_dir in sorted(glob(os.path.join(DATA_DIR, "sub-*", SESSION))):
    for c, contrast_dir in enumerate(sorted(
            glob(os.path.join(subject_dir, "session_stats", "contrast*")),
            key=natural_sort_key)):
        copes[c].append(os.path.join(contrast_dir, "stats", "cope1.nii.gz"))
        varcopes[c].append(os.path.join(contrast_dir, "stats",
                                        "varcope1.nii.gz"))

for contrast in range(len(copes)):
    wf = pe.Workflow(name="group")
    # Analyze contrasts of all subjects
    group = fawn.create_higher_level_workflow(mode=mode)
    group.inputs.inputspec.in_copes = copes[contrast]
    group.inputs.inputspec.in_varcopes = varcopes[contrast]
    group.inputs.inputspec.model = {"regressor_names": ["effect"],
                                    "regressors": [[1] * N_SUBJECTS]}
    group.inputs.inputspec.contrasts = [['group mean', 'T',['effect'],[1]]]

    # Threshold results
    thresholding = fawn.create_thresholding_workflow()

    # Save results
    results = pe.Node(io.DataSink(), base_directory = os.path.join(
        fmriprep_dir, "secondlevel", "transfer_rest", mode), name='results')

    # Connect and run
    wf.connect(group, 'outputspec.zstats', thresholding, 'inputspec.in_zstats')
    wf.connect(group, 'outputspec.copes', thresholding, 'inputspec.in_copes')
    wf.connect(group, 'outputspec.res4d', thresholding, 'inputspec.in_res4d')
    wf.connect(group, 'outputspec.mask', thresholding, 'inputspec.in_mask')
    wf.connect(group, 'outputspec.fwe_zstats', results,
               'contrast{0}_fwe'.format(contrast))
    wf.connect(thresholding, 'outputspec.cluster_zstats',
               results, 'contrast{0}_cluster'.format(contrast))
    wf.connect(thresholding, 'outputspec.cluster_pos_idx',
               results, 'contrast{0}_cluster.@pos_idx'.format(contrast))
    wf.connect(thresholding, 'outputspec.cluster_pos_max',
               results, 'contrast{0}_cluster.@pos_max'.format(contrast))
    wf.connect(thresholding, 'outputspec.cluster_pos_max_thresh',
               results, 'contrast{0}_cluster.@pos_max_thresh'.format(contrast))
    wf.connect(thresholding, 'outputspec.cluster_neg_idx',
               results, 'contrast{0}_cluster.@neg_idx'.format(contrast))
    wf.connect(thresholding, 'outputspec.cluster_neg_max',
               results, 'contrast{0}_cluster.@neg_max'.format(contrast))
    wf.connect(thresholding, 'outputspec.cluster_neg_max_thresh',
               results, 'contrast{0}_cluster.@neg_max_thresh'.format(contrast))
    wf.run()
