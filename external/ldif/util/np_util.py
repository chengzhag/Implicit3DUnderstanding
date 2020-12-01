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
"""Numpy utility functions."""

# Could penalize minimum volume for RBFs.
# Could penalize heavily if another quadric has higher weight within this quad's
# radius.
# Plot fraction in which 1 quadrics has 90% or more of the total weight.
# If total weight is below epsilon, don't show (or show which those are). Set
# this using a colab notebook, which will need to work first...
# Replace exp with some better falloff function!

import numpy as np


def make_coordinate_grid_3d(length, height, width, is_screen_space,
                            is_homogeneous):
    """Returns an array containing the coordinate grid values for a volume.

    Outputs a numpy array to avoid adding unnecessary operations to the graph.

    Args:
      length: int containing the length of the output volume.
      height: int containing the height of the output volume.
      width: int containing the width of the output volume.
      is_screen_space: bool. If true, then the coordinates are measured in pixels.
        If false, they are in the range (0-1).
      is_homogeneous: bool. If true, then a 1 is appended to the end of each
        coordinate.

    Returns:
      coords: numpy array of shape [length, height, width, 3] or
        [length, height, width, 4], depending on whether is_homogeneous is true.
        The value at location [i, j, k, :] is the (x,y,z) or (x,y,z,1) coordinate
        value at that location.
    """
    x_coords = np.linspace(0.5, width - 0.5, width)
    y_coords = np.linspace(0.5, height - 0.5, height)
    z_coords = np.linspace(0.5, length - 0.5, length)
    if not is_screen_space:
        x_coords /= width
        y_coords /= height
        z_coords /= length
    grid_x, grid_y, grid_z = np.meshgrid(
        x_coords, y_coords, z_coords, sparse=False, indexing='ij')
    if is_homogeneous:
        homogeneous_coords = np.ones_like(grid_x)
        grid = np.stack([grid_x, grid_y, grid_z, homogeneous_coords], axis=3)
    else:
        grid = np.stack([grid_x, grid_y, grid_z], axis=3)
    # Currently the order is (w, h, l), but we need (l, h, w) for
    # TensorFlow compatibility:
    return np.swapaxes(grid, 0, 2)
    # return grid
