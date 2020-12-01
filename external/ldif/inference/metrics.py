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
"""Computes metrics given predicted and ground truth shape."""

import numpy as np
import scipy

# ldif is an internal package, and should be imported last.
# pylint: disable=g-bad-import-order
# pylint: enable=g-bad-import-order

OCCNET_FSCORE_EPS = 1e-09


def sample_points_and_face_normals(mesh, sample_count):
    points, indices = mesh.sample(sample_count, return_index=True)
    points = points.astype(np.float32)
    normals = mesh.face_normals[indices]
    return points, normals


def pointcloud_neighbor_distances_indices(source_points, target_points):
    target_kdtree = scipy.spatial.cKDTree(target_points)
    distances, indices = target_kdtree.query(source_points, n_jobs=-1)
    return distances, indices


def dot_product(a, b):
    if len(a.shape) != 2:
        raise ValueError('Dot Product with input shape: %s' % repr(a.shape))
    if len(b.shape) != 2:
        raise ValueError('Dot Product with input shape: %s' % repr(b.shape))
    return np.sum(a * b, axis=1)


def normal_consistency(mesh1, mesh2, sample_count=100000, return_points=False):
    """Computes the normal consistency metric between two meshes."""
    points1, normals1 = sample_points_and_face_normals(mesh1, sample_count)
    points2, normals2 = sample_points_and_face_normals(mesh2, sample_count)

    _, indices12 = pointcloud_neighbor_distances_indices(points1, points2)
    _, indices21 = pointcloud_neighbor_distances_indices(points2, points1)

    normals12 = normals2[indices12]
    normals21 = normals1[indices21]

    # We take abs because the OccNet code takes abs...
    nc12 = np.abs(dot_product(normals1, normals12))
    nc21 = np.abs(dot_product(normals2, normals21))
    nc = 0.5 * np.mean(nc12) + 0.5 * np.mean(nc21)
    if return_points:
        return nc, points1, points2
    return nc


def percent_below(dists, thresh):
    return np.mean((dists ** 2 <= thresh).astype(np.float32)) * 100.0


def f_score(a_to_b, b_to_a, thresh):
    precision = percent_below(a_to_b, thresh)
    recall = percent_below(b_to_a, thresh)

    return (2 * precision * recall) / (precision + recall + OCCNET_FSCORE_EPS)


def fscore(mesh1,
           mesh2,
           sample_count=100000,
           tau=1e-04,
           points1=None,
           points2=None):
    """Computes the F-Score at tau between two meshes."""
    points1, points2 = get_points(mesh1, mesh2, points1, points2, sample_count)
    dist12, _ = pointcloud_neighbor_distances_indices(points1, points2)
    dist21, _ = pointcloud_neighbor_distances_indices(points2, points1)
    f_score_tau = f_score(dist12, dist21, tau)
    return f_score_tau


def mesh_chamfer_via_points(mesh1=None,
                            mesh2=None,
                            sample_count=100000,
                            points1=None,
                            points2=None):
    points1, points2 = get_points(mesh1, mesh2, points1, points2, sample_count)
    dist12, _ = pointcloud_neighbor_distances_indices(points1, points2)
    dist21, _ = pointcloud_neighbor_distances_indices(points2, points1)
    chamfer = 1000.0 * (np.mean(dist12 ** 2) + np.mean(dist21 ** 2))
    return chamfer


def get_points(mesh1, mesh2, points1, points2, sample_count):
    if points1 is not None or points2 is not None:
        assert points1 is not None and points2 is not None
    else:
        points1, _ = sample_points_and_face_normals(mesh1, sample_count)
        points2, _ = sample_points_and_face_normals(mesh2, sample_count)
    return points1, points2


def all_mesh_metrics(mesh1, mesh2, sample_count=100000):
    nc, points1, points2 = normal_consistency(
        mesh1, mesh2, sample_count, return_points=True)
    fs_tau = fscore(mesh1, mesh2, sample_count, 1e-04, points1, points2)
    fs_2tau = fscore(mesh1, mesh2, sample_count, 2.0 * 1e-04, points1, points2)
    chamfer = mesh_chamfer_via_points(mesh1, mesh2, sample_count, points1,
                                      points2)
    return {
        'fscore_tau': fs_tau,
        'fscore_2tau': fs_2tau,
        'chamfer': chamfer,
        'normal_consistency': nc
    }
