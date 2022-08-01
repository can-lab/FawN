"""Simple single subject analysis.

This is an example of a simple single subject analysis with a single session of
multiple (identical) runs.

"""


import os
from glob import glob

import pandas as pd
from nipype.pipeline import engine as pe
from nipype.interfaces import io

import fawn

FUNC_DIR = "/path/to/data"
SEQUENCE_NAME = "MYSEQUENCE"
PREPROC = "5mm100sNone"
TR = 1.0

wf = pe.Workflow(name="single_subject")

# Analyze all runs in the session in one go
session = fawn.create_session_level_workflow(tr=TR)
session.inputs.inputspec.in_files = sorted(glob(os.path.join(
    FUNC_DIR, "*{0}*desc-preproc{1}*_bold.nii.gz".format(SEQUENCE_NAME,
                                                         PREPROC))))
session.inputs.inputspec.in_masks = sorted(glob(os.path.join(
    FUNC_DIR, "*{0}*desc-brain_mask.nii.gz".format(SEQUENCE_NAME))))

# Specify design
session.inputs.inputspec.models = []
motion_files = sorted(glob(os.path.join(
    FUNC_DIR, "*{0}*confounds_regressors.tsv".format(SEUQENCE_NAME))))
event_files = sorted(glob(os.path.join(
    FUNC_DIR, "*{0}*events.tsv".format(SEQUENCE_NAME))))
files = zip(event_files, motion_files)
for f in files:
    events = pd.read_csv(f[0], sep='\t')
    condition_names = list(set(events.trial_type))
    motion = pd.read_csv(f[1], sep='\t')
    motion_names = ['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']
    session.inputs.inputspec.models.append(
        {"conditions": condition_names,
         "onsets": [events[events.trial_type == c]["onset"].values \
                    for c in condition_names],
         "durations": [events[events.trial_type == c]["duration"].values \
                       for c in condition_names],
         "regressor_names": motion_names,
         "regressors": motion[motion_names].fillna(0).values.T})

# Specify contrasts
session.inputs.inputspec.contrasts = [
    ('cond1', 'T', ['cond1', 'cond2'], [1, 0]),
    ('cond2', 'T', ['cond1', 'cond2'], [0, 1]),
    ('cond1 > cond2', 'T', ['cond1', 'cond2'], [1, -1])]

# Save results
results = pe.Node(io.DataSink(regexp_substitutions=[("_flameo", "contrast")]),
                  base_directory=FUNC_DIR, name='results')

# Connect and run
wf.connect(session, 'outputspec.stats_dir', results, 'session_stats')
wf.run('MultiProc', plugin_args={'n_procs': 8})
