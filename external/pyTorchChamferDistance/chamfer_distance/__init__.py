import os
os.makedirs(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'build')), exist_ok=True)
from .chamfer_distance import ChamferDistance
