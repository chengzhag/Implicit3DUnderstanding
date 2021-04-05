# Lint as: python3
"""A shell utility to get the CUDA version."""

import subprocess as sp

# nvcc output for Cuda 11
# nvcc: NVIDIA (R) Cuda compiler driver
# Copyright (c) 2005-2020 NVIDIA Corporation
# Built on Mon_Nov_30_19:08:53_PST_2020
# Cuda compilation tools, release 11.2, V11.2.67
# Build cuda_11.2.r11.2/compiler.29373293_0

# nvcc output for Cuda 10
# output = "nvcc: NVIDIA (R) Cuda compiler driver\n" \
#          "Copyright (c) 2005-2019 NVIDIA Corporation\n" \
#          "Built on Sun_Jul_28_19:07:16_PDT_2019\n" \
#          "Cuda compilation tools, release 10.1, V10.1.243\n"

def get_cuda_version():
  try:
    output = sp.check_output(['nvcc', '-V']).decode('utf-8')
    lines = output.split('\n')
    if output.find('V11') < 0:
      version_str = lines[-2].split(',')[1].split(' ')[-1]
    else:
      version_str = lines[-3].split(',')[1].split(' ')[-1]
      
    major_version = int(version_str.split('.')[0])
    minor_version = int(version_str.split('.')[1])
    return major_version, minor_version
  except Exception as e:
    raise ValueError(f'Failed to get cuda version with error: {e}')

if __name__ == '__main__':
  major, minor = get_cuda_version()
  print(major)
  print(minor)
