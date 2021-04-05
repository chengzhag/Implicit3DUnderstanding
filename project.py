import wandb
import argparse
import glob
from tqdm import tqdm
import shutil
import os
import subprocess
import git
import sys
python_bin = sys.executable

def clean(args):
    api = wandb.Api()

    print("Cleaning project outputs to sync with wandb...")
    print("Retrieving wandb runs...")
    runs = api.runs(f"{args.username}/{args.project}")
    wandb_paths = []
    for run in runs:
        wandb_path = run.config.get('log', {}).get('path', None)
        if wandb_path is None:
            print(f"run '{run.id}' do not have log.path config")
            continue
        if not wandb_path.endswith('/'):
            wandb_path += '/'
        wandb_paths.append(wandb_path)
    print(f"Total: {len(wandb_paths)}")
    print(wandb_paths)

    print("Scanning local runs...")
    local_paths = glob.glob("out/*/*/")
    print(f"Total: {len(local_paths)}")
    print(local_paths)

    to_deletes = set(local_paths) - set(wandb_paths)
    to_deletes = [p for p in to_deletes if 'pretrained' not in p]
    to_deletes.sort()
    for to_delete in to_deletes:
        print(to_delete)
    print(f"{len(to_deletes)} to be deleted")
    input("Press Enter to continue...")

    for to_delete in tqdm(to_deletes):
        tqdm.write(f"deleting {to_delete}...")
        try:
            shutil.rmtree(to_delete)
        except:
            print(f"Failed deleting {to_delete}...")

def clone(args):
    dst = args.dst
    if not dst.endswith('/'):
        dst += '/'
    src = os.path.split(__file__)[0]

    print("Cloning project...")
    repo = git.Repo(src)
    repo.clone(path=dst)
    shutil.rmtree(os.path.join(dst, '.idea/'))

    print("Copying Pycharm project folder...")
    shutil.copytree(os.path.join(src, '.idea/'), os.path.join(dst, '.idea/'))
    remote_mappings = os.path.join(dst, '.idea', 'remote-mappings.xml')
    if os.path.exists(remote_mappings):
        subprocess.check_output(f"sed -i s#{src}#{os.path.abspath(os.path.join(src, dst))}#g {remote_mappings}",
                                shell=True)

    print("Adding soft link for essential files...")
    os.symlink(os.path.join(src, 'out'), os.path.join(dst, 'out'))

    rel_dirs = glob.glob(os.path.join(src, 'data', '*', '*/'))
    rel_dirs = [os.path.relpath(d, src) for d in rel_dirs]
    for rel_dir in rel_dirs:
        dst_folder = os.path.join(dst, rel_dir)
        if not os.path.exists(os.path.join(dst, rel_dir)):
            src_folder = os.path.join(src, rel_dir)
            print(f"linking {rel_dir} to {dst_folder}")
            os.symlink(src_folder, dst_folder)

    print("Building project...")
    cmd = ' && '.join([
        f"{python_bin} {os.path.join(dst, 'project.py')} build"
    ])
    print(f"Run: {cmd}")
    subprocess.check_output(cmd, shell=True)

def build(args):
    root = os.path.split(__file__)[0]

    if args.subwork is None or 'gaps' in args.subwork:
        print("Building gaps...")
        cmd = ' && '.join([
            f"cd {os.path.join(root, 'external')}",
            f"bash build_gaps.sh"
        ])
        print(f"Run: {cmd}")
        subprocess.check_output(cmd, shell=True)

    if args.subwork is None or 'mesh_fusion' in args.subwork:
        print("Building pyfusion...")
        mesh_fusion_path = os.path.join(root, 'external', 'mesh_fusion')
        cmd = ' && '.join([
            f"cd {os.path.join(mesh_fusion_path, 'libfusiongpu')}",
            "mkdir -p build",
            "cd build",
            "cmake ..",
            "make",
            "cd ..",
            f"{python_bin} setup.py build_ext --inplace"
        ])
        print(f"Run: {cmd}")
        subprocess.check_output(cmd, shell=True)

        print("Building pyrender...")
        cmd = ' && '.join([
            f"cd {os.path.join(mesh_fusion_path, 'librender')}",
            f"{python_bin} setup.py build_ext --inplace"
        ])
        print(f"Run: {cmd}")
        subprocess.check_output(cmd, shell=True)

        print("Building PyMCubes...")
        cmd = ' && '.join([
            f"cd {os.path.join(mesh_fusion_path, 'libmcubes')}",
            f"{python_bin} setup.py build_ext --inplace"
        ])
        print(f"Run: {cmd}")
        subprocess.check_output(cmd, shell=True)

    if args.subwork is None or 'ldif2mesh' in args.subwork:
        print("Building ldif2mesh...")
        ldif2mesh_path = os.path.join(root, 'external', 'ldif', 'ldif2mesh')
        cmd = ' && '.join([
            f"cd {ldif2mesh_path}",
            f"bash build.sh"
        ])
        print(f"Run: {cmd}")
        subprocess.check_output(cmd, shell=True)
        subprocess.check_output(f"chmod 744 {os.path.join(ldif2mesh_path, 'ldif2mesh')}", shell=True)

        print(f"Run: {cmd}")
        os.system(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Project utils.')
    parser.add_argument('work', type=str, default='clean',
                        help='type of work to do.')
    parser.add_argument('--subwork', type=str, nargs='+', default=None,
                        help='list of subworks to do.')
    parser.add_argument('--username', type=str, default='pidan1231239',
                        help='wandb username.')
    parser.add_argument('--project', type=str, default='implicit3dunderstanding',
                        help='wandb project.')
    parser.add_argument('--dst', type=str, default='../Implicit3DUnderstanding_clone')
    args = parser.parse_args()

    globals()[args.work](args)
