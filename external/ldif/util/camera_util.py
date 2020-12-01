import torch


def roll_pitch_yaw_to_rotation_matrices(roll_pitch_yaw):
    # Converts roll-pitch-yaw angles to rotation matrices.
    # (ldif.util.camera_util.roll_pitch_yaw_to_rotation_matrices)
    cosines = torch.cos(roll_pitch_yaw)
    sines = torch.sin(roll_pitch_yaw)
    cx, cy, cz = cosines.unbind(-1)
    sx, sy, sz = sines.unbind(-1)
    rotation = torch.stack(
        [cz * cy, cz * sy * sx - sz * cx, cz * sy * cx + sz * sx,
         sz * cy, sz * sy * sx + cz * cx, sz * sy * cx - cz * sx,
         -sy, cy * sx, cy * cx], -1
    )
    rotation = torch.reshape(rotation, [rotation.shape[0], -1, 3, 3])
    return rotation
