from skimage import measure
import numpy as np
import torch
from .sdf import create_grid, eval_grid_octree, eval_grid
from skimage import measure
import trimesh


def reconstruction(structured_implicit,
                   resolution, b_min, b_max,
                   use_octree=False, num_samples=10000,
                   marching_cube=True):
    '''
    Reconstruct meshes from sdf predicted by the network.
    :param structured_implicit: a StructuredImplicit object.
    :param resolution: resolution of the grid cell
    :param b_min: bounding box corner [x_min, y_min, z_min]
    :param b_max: bounding box corner [x_max, y_max, z_max]
    :param use_octree: whether to use octree acceleration
    :param num_samples: how many points to query each gpu iteration
    :return: marching cubes results.
    '''
    # First we create a grid by resolution
    # and transforming matrix for grid coordinates to real world xyz
    coords, mat = create_grid(resolution, resolution, resolution,
                              b_min, b_max)

    # Then we define the lambda function for cell evaluation
    def eval_func(points, structured_implicit):
        points = np.expand_dims(points, axis=0)
        samples = torch.from_numpy(points).to(device=structured_implicit.device).float()
        samples = samples.transpose(-1, -2)
        samples = samples.expand(structured_implicit.batch_size, -1, -1)
        pred = structured_implicit.class_at_samples(samples, apply_class_transfer=False)[0][..., 0]
        return pred.detach().cpu().numpy()

    # Then we evaluate the grid
    if use_octree:
        sdf = []
        for s in structured_implicit.unbind():
            sdf.append(eval_grid_octree(coords, lambda p:eval_func(p, s), num_samples=num_samples).squeeze())
        sdf = np.stack(sdf)
    else:
        sdf = eval_grid(coords, lambda p:eval_func(p, structured_implicit),
                        num_samples=num_samples, batch_size=structured_implicit.batch_size)

    # Finally we do marching cubes
    if marching_cube:
        mesh = []
        for s in sdf:
            try:
                verts, faces, _, _ = measure.marching_cubes(s, -0.07)
                # transform verts into world coordinate system
                verts = np.matmul(mat[:3, :3], verts.T) + mat[:3, 3:4]
                verts = verts.T
                mesh.append(trimesh.Trimesh(vertices=verts, faces=faces))
            except (ValueError, RuntimeError) as e:
                print('Failed to extract mesh with error %s. Setting to unit sphere.' % repr(e))
                mesh.append(trimesh.primitives.Sphere(radius=0.5))
        return mesh
    else:
        return sdf, mat


def save_obj_mesh(mesh_path, verts, faces):
    file = open(mesh_path, 'w')

    for v in verts:
        file.write('v %.4f %.4f %.4f\n' % (v[0], v[1], v[2]))
    for f in faces:
        f_plus = f + 1
        file.write('f %d %d %d\n' % (f_plus[0], f_plus[2], f_plus[1]))
    file.close()


def save_obj_mesh_with_color(mesh_path, verts, faces, colors):
    file = open(mesh_path, 'w')

    for idx, v in enumerate(verts):
        c = colors[idx]
        file.write('v %.4f %.4f %.4f %.4f %.4f %.4f\n' % (v[0], v[1], v[2], c[0], c[1], c[2]))
    for f in faces:
        f_plus = f + 1
        file.write('f %d %d %d\n' % (f_plus[0], f_plus[2], f_plus[1]))
    file.close()


def save_obj_mesh_with_uv(mesh_path, verts, faces, uvs):
    file = open(mesh_path, 'w')

    for idx, v in enumerate(verts):
        vt = uvs[idx]
        file.write('v %.4f %.4f %.4f\n' % (v[0], v[1], v[2]))
        file.write('vt %.4f %.4f\n' % (vt[0], vt[1]))

    for f in faces:
        f_plus = f + 1
        file.write('f %d/%d %d/%d %d/%d\n' % (f_plus[0], f_plus[0],
                                              f_plus[2], f_plus[2],
                                              f_plus[1], f_plus[1]))
    file.close()
