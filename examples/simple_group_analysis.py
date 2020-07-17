"""Simple group analysis.

This is an example of a simple group level analysis including thresholding.

"""

import os
import re
from glob import glob

from nipype.pipeline import engine as pe
from nipype.interfaces import io

from fawn import workflows as fw


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
    for c, contrast_dir in enumerate(sorted(glob(os.path.join(subject_dir, "session_level_stats", "contrast*")), key=natural_sort_key)):
        copes[c].append(os.path.join(contrast_dir, "stats", "cope1.nii.gz"))
        varcopes[c].append(os.path.join(contrast_dir, "stats", "varcope1.nii.gz"))

for contrast in range(len(copes)):
    # Group level analysis
    wf = pe.Workflow(name="group")
    group_level = fw.create_higher_level_workflow(mode=mode)
    group_level.inputs.inputspec.in_copes = copes[contrast]
    group_level.inputs.inputspec.in_varcopes = varcopes[contrast]
    group_level.inputs.inputspec.model = model = {"regressor_names": ["effect"],
                                                  "regressors": [[1] * N_SUBJECTS]}
    group_level.inputs.inputspec.contrasts = [['group mean', 'T',['effect'],[1]]]

    # Thresholding
    thresholding = fw.create_thresholding_workflow()
    results = pe.Node(io.DataSink(),
                      base_directory = os.path.join(fmriprep_dir, "secondlevel", "transfer_rest", mode),
                      name='results')

    wf.connect(second_level, 'outputspec.zstats', thresholding, 'inputspec.in_zstats')
    wf.connect(second_level, 'outputspec.copes', thresholding, 'inputspec.in_copes')
    wf.connect(second_level, 'outputspec.res4d', thresholding, 'inputspec.in_res4d')
    wf.connect(second_level, 'outputspec.mask', thresholding, 'inputspec.in_mask')
    wf.connect(thresholding, 'outputspec.fwe_zstats', results, 'contrast{0}_fwe'.format(contrast))
    wf.connect(thresholding, 'outputspec.cluster_zstats', results, 'contrast{0}_cluster'.format(contrast))
    wf.connect(thresholding, 'outputspec.cluster_pos_idx', results, 'contrast{0}_cluster.@pos_idx'.format(contrast))
    wf.connect(thresholding, 'outputspec.cluster_pos_max', results, 'contrast{0}_cluster.@pos_max'.format(contrast))
    wf.connect(thresholding, 'outputspec.cluster_pos_max_thresh', results, 'contrast{0}_cluster.@pos_max_thresh'.format(contrast))
    wf.connect(thresholding, 'outputspec.cluster_neg_idx', results, 'contrast{0}_cluster.@neg_idx'.format(contrast))
    wf.connect(thresholding, 'outputspec.cluster_neg_max', results, 'contrast{0}_cluster.@neg_max'.format(contrast))
    wf.connect(thresholding, 'outputspec.cluster_neg_max_thresh', results, 'contrast{0}_cluster.@neg_max_thresh'.format(contrast))

    wf.run()
