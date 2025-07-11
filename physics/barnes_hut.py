import math

class QuadTreeNode:
    def __init__(self, x_min, x_max, y_min, y_max):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.body = None
        self.children = []
        self.mass = 0.0
        self.com_x = 0.0
        self.com_y = 0.0

    def insert(self, body):
        # Если в узле нет тела, просто вставляем
        if self.body is None and not self.children:
            self.body = body
            self.mass = body.m
            self.com_x = body.x
            self.com_y = body.y
            return

        # Если уже есть тело и нет детей — делим узел
        if not self.children:
            self.subdivide()
            self._insert_body_to_child(self.body)
            self.body = None

        self._insert_body_to_child(body)
        # После вставки обновляем массу и центр масс
        self.update_mass_and_com()

    def _insert_body_to_child(self, body):
        mx = (self.x_min + self.x_max) / 2
        my = (self.y_min + self.y_max) / 2
        # Определяем квадрант
        idx = 0
        if body.x > mx:
            idx += 1
        if body.y > my:
            idx += 2
        self.children[idx].insert(body)

    def subdivide(self):
        mx = (self.x_min + self.x_max) / 2
        my = (self.y_min + self.y_max) / 2
        self.children = [
            QuadTreeNode(self.x_min, mx, self.y_min, my),  # NW
            QuadTreeNode(mx, self.x_max, self.y_min, my),  # NE
            QuadTreeNode(self.x_min, mx, my, self.y_max),  # SW
            QuadTreeNode(mx, self.x_max, my, self.y_max),  # SE
        ]

    def update_mass_and_com(self):
        total_mass = 0.0
        com_x = 0.0
        com_y = 0.0
        for child in self.children:
            total_mass += child.mass
            com_x += child.com_x * child.mass
            com_y += child.com_y * child.mass
        if total_mass > 0:
            com_x /= total_mass
            com_y /= total_mass
        self.mass = total_mass
        self.com_x = com_x
        self.com_y = com_y

def compute_forces(bodies, G, theta=0.5):
    if not bodies:
        return []

    # 1. Определяем границы пространства
    min_x = min(b.x for b in bodies)
    max_x = max(b.x for b in bodies)
    min_y = min(b.y for b in bodies)
    max_y = max(b.y for b in bodies)
    size = max(max_x - min_x, max_y - min_y)
    cx = (min_x + max_x) / 2
    cy = (min_y + max_y) / 2
    half = size / 2 + 1e-5

    root = QuadTreeNode(cx - half, cx + half, cy - half, cy + half)
    for body in bodies:
        root.insert(body)

    # 2. Для каждого тела вычисляем силу через обход дерева
    forces = []
    for body in bodies:
        fx, fy = _compute_force_on_body(body, root, G, theta)
        forces.append((fx, fy))
    return forces

def _compute_force_on_body(body, node, G, theta):
    # Если узел пустой или содержит только текущее тело — силы нет
    if (node.body is body and not node.children) or node.mass == 0:
        return 0.0, 0.0

    dx = node.com_x - body.x
    dy = node.com_y - body.y
    dist2 = dx * dx + dy * dy + 1e-8  # чтобы не делить на 0
    dist = math.sqrt(dist2)
    width = node.x_max - node.x_min

    # Barnes-Hut criterion
    if (not node.children) or (width / dist < theta):
        # Аппроксимация: вся масса в центре масс
        if node.com_x == body.x and node.com_y == body.y:
            return 0.0, 0.0
        force = G * body.m * node.mass / dist2
        fx = force * dx / dist
        fy = force * dy / dist
        return fx, fy
    else:
        # Рекурсивно обходим детей
        fx, fy = 0.0, 0.0
        for child in node.children:
            cfx, cfy = _compute_force_on_body(body, child, G, theta)
            fx += cfx
            fy += cfy
        return fx, fy 