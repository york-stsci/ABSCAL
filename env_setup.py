#! /usr/bin/env python
"""
This file sets up conda shell scripts that will activate and deactivate needed environment
variables as the abscal environment is started/stopped. In particular:

- CRDS_PATH will be set as follows:
    - if there is a CRDS_CACHE_TYPE variable, and CRDS_CACHE_TYPE='local', set to $HOME/crds_cache
    - else (internal) set to /grp/crds/cache
- oref is set to $CRDS_PATH/references/hst/oref

Author
-------
- Brian York

Use
---
This file is intended to be run from the command line::

    abscal_setup [options] table
"""

import json
import os
import platform
import shutil
import subprocess
import sys
from urllib.request import urlopen

try:
    import yaml
except Exception as e:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pyyaml'])
    import yaml

try:
    from ruamel.yaml import YAML
except Exception as e:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'ruamel.yaml'])
    from ruamel.yaml import YAML

def set_var(name, value, setup_file, remove_file):
    """
    For an environment variable, write the necessary text to:
    - Preserve any existing variable by that name
    - Create the variable in activate.sh
    - Remove the variable in deactivate.sh
    - Restore any existing variable that was created
    """
    with open(setup_file, 'a+') as f:
        write_str = 'export {}_CONDA_BACKUP='.format(name)
        write_str += '${'+name+':-}\n'
        write_str += 'export {0}={1}\n'.format(name, value)
        f.write(write_str)

    with open(remove_file, 'a+') as f:
        write_str = 'export {}='.format(name)
        write_str += "${" + "{}_CONDA_BACKUP:-".format(name) + '}\nunset '
        write_str += "{}_CONDA_BACKUP\n".format(name)
        f.write(write_str)
        f.write('if [ -z ${0} ]; then\n    unset {0}\nfi\n'.format(name))


def main(**kwargs):
    with open("env_config.yml") as config_file:
        config = yaml.safe_load(config_file)
    
    conf = config['base']
    
    msg = "Creating {} environment for python {}\n"
    sys.stdout.write(msg.format(conf['name'], conf['python_version']))
    
    sys.stdout.write("Fetching latest tag from github.\n")
    try:
        # Find the URL of the latest STScI release tag
        json_data = json.loads(urlopen(conf['info_url']).read())
        stenv_tag = json_data['tag_name']
        
        if stenv_tag > conf["latest_known_tag"]:
            msg = "Updating latest tag from {} to {}\n"
            sys.stdout.write(msg.format(conf['latest_known_tag'], stenv_tag))
            
            ryaml = YAML()
            with open("env_config.yml") as inf:
                new_config = ryaml.load(inf)
            new_config["latest_known_tag"] = stenv_tag
            with open("env_config.yml", mode="w") as outf:
                ryaml.dump(new_config, outf)
    except Exception as e:
        sys.stderr.write("Retrieval error: {}\n".format(e))
        stenv_tag = conf["latest_known_tag"]
        sys.stdout.write("\tFalling back to tag {}\n".format(stenv_tag))
    
    sys.stdout.write("Latest STENV is {}.\n".format(stenv_tag))
    
    telescope = conf['telescope']
    
    sys.stdout.write("Setting up an environment for {} use".format(telescope))
    
    # Find the platform we're running on. Currently support Linux and MacOS, because that's
    # what ``stenv`` supports
    system = platform.system()
    if system == "Darwin":
        stenv_sys = "macOS"
        arch_sys = "MacOSX"
        arch = platform.processor()
        if arch == "arm":
            stenv_sys += "-ARM64"
        else:
            stenv_sys += "-X64"
    elif system == "Linux":
        stenv_sys = "Linux"
        arch_sys = "Linux"
    else:
        msg = "ERROR: stenv only supports macOS and Linux, not {}\n."
        sys.stderr.write(msg.format(system))
        sys.exit(1)
    
    machine = platform.machine()
    
    sys.stdout.write("Running on a {} {} system.\n".format(stenv_sys, machine))
    
    stenv_url = "{0}/{1}/stenv-{2}-py{3}-{1}.yaml"
    stenv_url = stenv_url.format(conf['stenv_url'], stenv_tag, stenv_sys, 
        conf['python_version'])
    
    # Check if we need to install conda
    if shutil.which("conda") is None:
        conda_url = conf['conda_url'].format(arch_sys, machine)
        conda_path = os.path.expandvars(conf['conda_path'])
        conda_file = "conda_install.sh"
        subprocess.run(["wget", conda_url, "-O", conda_file])
        subprocess.run(["bash", conda_file, "-b", "-p", conda_path])
        conda_prefix = conda_path
    else:
        conda_prefix = os.path.expandvars(os.environ["CONDA_PREFIX"])
        if os.environ.get("CONDA_DEFAULT_ENV", "") != "base":
            msg = "ERROR: This setup script must be run from the base environment.\n"
            sys.stderr.write(msg)
            sys.exit(1)
    conda_cmd = os.path.join(conda_prefix, "bin", "conda-env")
    pip_install_prefix = ['conda', 'run', '-n', conf['name'], 'python', '-m', 'pip', 'install']
    
    sys.stdout.write("Initializing {} from {}.\n".format(conf['name'], stenv_url))
    
    env_path = os.path.join(conda_prefix, "envs", conf['name'])
    if os.path.exists(env_path):
        msg = "ERROR: conda environment {} already exists. You must delete your "
        msg += "existing environment before creating a new one.\n"
        sys.stderr.write(msg)
        sys.exit(1)
    
    # Create basic STENV environment
    if system == "Darwin" and arch == "arm":
        sys.stdout.write("Creating initial empty environment.\n")
        base_conda_cmd = os.path.join(conda_prefix, "bin", "conda")
        subprocess.run([base_conda_cmd, "create", "-n", conf['name'], 'python={}'.format(conf['python_version'])])
        sys.stdout.write("Installing numpy with apple silicon acceleration.\n")
        subprocess.check_call(pip_install_prefix+['numpy'])
        sys.stdout.write("Configuring conda.\n")
        subprocess.run(["conda", "config", "--set", "pip_interop_enabled", "true"])
        sys.stdout.write("Updating environment from URL.\n")
        subprocess.run(["conda", "env", "update", "--name", conf['name'], "--file", stenv_url])
    else:
        sys.stdout.write("Creating initial environment from URL.\n")
        subprocess.run([conda_cmd, "create", "-n", conf['name'], "--file", stenv_url], 
                       env=os.environ.copy())

    # Elaborate with additional configuration
    if "local_env" in conf:
        with open(conf['local_env']) as inf:
            additional_config = yaml.safe_load(inf)
        for item in additional_config.get("dependencies", []):
            if isinstance(item, dict) and "pip" in item:
                for package in item["pip"]:
                    if package == "-e .":
                        sys.stdout.write("Installing local module with pip.\n")
                        subprocess.check_call(pip_install_prefix+['-e', os.getcwd()])
                    else:
                        sys.stdout.write("Installing {} with pip.\n".format(package))
                        subprocess.check_call(pip_install_prefix+[package])
# This is commented out because currently the constructed STENV is at least sometimes
# inconsistent, and so I'm installing pip things manually.
#         msg = "Updating {} environment from local file {}.\n"
#         sys.stdout.write(msg.format(conf['name'], conf['local_env']))
#         subprocess.run([conda_cmd, "update", "--name", conf['name'], "--file", 
#             conf['local_env']], env=os.environ.copy())
    
    activate_path = os.path.join(env_path, "etc", "conda", "activate.d")
    if not os.path.exists(activate_path):
        os.makedirs(activate_path)
    activate_file = os.path.join(activate_path, "{}_activate.sh".format(conf['name']))

    deactivate_path = os.path.join(env_path, "etc", "conda", "deactivate.d")
    if not os.path.exists(deactivate_path):
        os.makedirs(deactivate_path)
    deactivate_file = os.path.join(deactivate_path, "{}_deactivate.sh".format(conf['name']))
    
    # Set up info variables for the environment
    set_var('{}_version'.format(conf['name']), conf['package_version'], activate_file, deactivate_file)
    set_var('{}_date'.format(conf['name']), conf['env_date'], activate_file, deactivate_file)
    set_var('stenv_version', stenv_tag, activate_file, deactivate_file)
        
    if 'env' in config:
        internal_crds = False
        if 'CRDS_PATH' in os.environ and os.environ['CRDS_PATH'][:4] == "/grp":
            internal_crds = True
        
        if 'special' in config['env']:
            if "package" in config['env']['special']:
                # Add the local module
                spkg = os.path.join(env_path, "lib", "python{}".format(conf['python_version']),
                                    "site-packages")
                if os.path.isfile(os.path.join(spkg, "{}.egg-link".format(conf['package_name']))):
                    # Editable installation.
                    with open(os.path.join(spkg, "{}.egg-link".format(conf['package_name']))) as f:
                        var_value = f.readline().strip()
                else:
                    var_value = os.path.join(spkg, conf['package_name'])
                var_name = "{}_PKG".format(conf['package_name'])
                set_var(var_name, var_value, activate_file, deactivate_file)
            if 'CRDS_PATH' in config['env']['special']:
                # This is important because, especially for internal STScI people, the CRDS
                # cache path isn't *quite* the same as it is otherwise.
                if internal_crds:
                    value = '/grp/crds/{}'.format(telescope)
                elif 'CRDS_PATH' not in os.environ:
                    value = '$HOME/crds_cache'
                else:
                    value = os.environ['CRDS_PATH']
                set_var('CRDS_PATH', value, activate_file, deactivate_file)
        
        if 'crds_vars' in config['env']:
            if 'CRDS_PATH' not in os.environ:
                set_var('CRDS_PATH', '$HOME/crds_cache', activate_file, deactivate_file)
            for item in config['env']['crds_vars']:
                if internal_crds:
                    value = '$CRDS_PATH/references/{}'.format(telescope)
                else:
                    value = '$CRDS_PATH/references/{}/{}'.format(telescope, item)
                set_var(item, value, activate_file, deactivate_file)
        
        if 'vars' in config['env']:
            for item in config['env']['vars']:
                for key in item:
                    msg = "Adding environment variable {}={}\n"
                    sys.stdout.write(msg.format(key, item[key]))
                    set_var(key, item[key], activate_file, deactivate_file)


if __name__ == "__main__":
    main()
