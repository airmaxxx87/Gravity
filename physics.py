# physics.py

def compute_forces(bodies, G):
    """
    Вычисляет гравитационные силы для каждого тела.
    Возвращает список (fx, fy) для каждого тела.
    """
    forces = [(0.0, 0.0) for _ in bodies]
    for i, body in enumerate(bodies):
        fx, fy = 0.0, 0.0
        for j, other in enumerate(bodies):
            if i == j:
                continue
            dx = other.x - body.x
            dy = other.y - body.y
            dist2 = dx * dx + dy * dy
            if dist2 == 0:
                continue
            force = G * body.m * other.m / dist2
            dist = dist2 ** 0.5
            fx += force * dx / dist
            fy += force * dy / dist
        forces[i] = (fx, fy)
    return forces

def update_bodies(bodies, forces, dt):
    """
    Обновляет состояния тел по рассчитанным силам.
    """
    for body, (fx, fy) in zip(bodies, forces):
        body.update(fx, fy, dt) 