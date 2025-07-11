from .naive_physics import compute_forces as naive_compute_forces, update_bodies
from .barnes_hut import compute_forces as barnes_hut_compute_forces

class PhysicsEngine:
    def __init__(self, algorithm="naive", theta=0.5):
        self.algorithm = algorithm
        self.theta = theta

    def compute_forces(self, bodies, G):
        if self.algorithm == "naive":
            return naive_compute_forces(bodies, G)
        elif self.algorithm == "barnes_hut":
            return barnes_hut_compute_forces(bodies, G, self.theta)
        else:
            raise ValueError("Unknown algorithm: " + self.algorithm)

    def update_bodies(self, bodies, forces, dt):
        update_bodies(bodies, forces, dt) 