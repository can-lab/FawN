from setuptools import setup

def get_version():
    """Get version and version_info from fawn/__meta__.py file."""

    import os
    module_path = os.path.join(os.path.dirname('__file__'), 'fawn',
                               '__meta__.py')

    import importlib.util
    spec = importlib.util.spec_from_file_location('__meta__', module_path)
    meta = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(meta)
    return meta.__version__

__version__ = get_version()

setup(
    name = 'FawN',
    version = __version__,
    packages = ['fawn'],
    package_data = {'fawn': ['examples/*.*']},
    install_requires = ['nipype', 'fslpy'])
