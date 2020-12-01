# Copyright 2020 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# Lint as: python3
"""Utilties for working with the GAPS geometric processing library."""

import os
import numpy as np


# LDIF is an internal package, should be imported last.
# pylint: disable=g-bad-import-order
from . import file_util
from .file_util import log
# pylint: enable=g-bad-import-order


def read_pts_file(path):
  """Reads a .pts or a .sdf point samples file."""
  _, ext = os.path.splitext(path)
  assert ext in ['.sdf', '.pts']
  l = 4 if ext == '.sdf' else 6
  with file_util.open_file(path, 'rb') as f:
    points = np.fromfile(f, dtype=np.float32)
  points = np.reshape(points, [-1, l])
  return points
