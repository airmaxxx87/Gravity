import tkinter as tk
import math
import copy  # Для глубокого копирования состояний
import random  # Для случайных орбит спутников
import itertools
from physics.physics_interface import PhysicsEngine

# --- Параметры системы ---
# WIDTH, HEIGHT убираем, теперь берём из размеров окна
# CENTER теперь вычисляется динамически
PLANET_COLORS = ['#1e90ff', '#32cd32', '#ff6347', '#ffb300', '#8e44ad', '#00b894', '#e17055', '#0984e3']
SATELLITE_COLOR = '#aaaaaa'
G = 0.2  # гравитационная постоянная (масштабированная)
DT_PHYSICS = 0.5  # фиксированный шаг времени для физики

TRAIL_POINT_STEP = 10  # шаг между точками траектории

class Body:
    """
    Класс для небесного тела: id, Name, масса (m), радиус (r), цвет (Color), координаты (x, y),
    вектор скорости (xv, yv), вектор силы (xf, yf)
    """
    _id_counter = itertools.count()
    def __init__(self, Name, x, y, xv, yv, m, Color):
        self.id = next(Body._id_counter)
        self.Name = Name
        self.m = m
        self.r = 2.5 * (m ** (1/3))
        self.Color = Color
        self.x = x
        self.y = y
        self.xv = xv
        self.yv = yv
        self.xf = 0.0
        self.yf = 0.0
        self.trail = []  # след для визуализации орбиты (мировые координаты)

    def update(self, fx, fy, dt):
        # Обновление скорости и положения
        self.xf = fx
        self.yf = fy
        self.xv += fx / self.m * dt
        self.yv += fy / self.m * dt
        self.x += self.xv * dt
        self.y += self.yv * dt
        self.trail.append((self.x, self.y))
        if len(self.trail) > 200:
            self.trail.pop(0)

class GravitySim:
    def __init__(self, master):
        self.iteration_counter = 0
        self.master = master
        self.master.title('Гравитационный симулятор')
        self.master.state('zoomed')
        self.width = self.master.winfo_screenwidth()
        self.height = self.master.winfo_screenheight()
        self.show_trails = tk.BooleanVar(value=True)  # Новая переменная состояния

        # --- Панель управления ---
        self.top_frame = tk.Frame(master, bg='black')
        self.top_frame.pack(side=tk.TOP, fill=tk.X)
        self.show_vectors_var = tk.BooleanVar(value=True)
        self.cb_vectors = tk.Checkbutton(
            self.top_frame, text='Показывать векторы', variable=self.show_vectors_var,
            command=self.toggle_vectors, fg='white', bg='black', selectcolor='black',
            activebackground='black', activeforeground='white'
        )
        self.cb_vectors.pack(side=tk.LEFT, padx=10, pady=2)

        # --- Новый чекбокс для траекторий ---
        self.cb_trails = tk.Checkbutton(
            self.top_frame, text='Показывать траектории', variable=self.show_trails,
            command=self.toggle_trails, fg='white', bg='black', selectcolor='black',
            activebackground='black', activeforeground='white'
        )
        self.cb_trails.pack(side=tk.LEFT, padx=10, pady=2)
        # Кнопка сброса вида
        self.reset_btn = tk.Button(self.top_frame, text='Сбросить вид', command=self.reset_view, fg='white', bg='black', activebackground='gray', activeforeground='white')
        self.reset_btn.pack(side=tk.LEFT, padx=10, pady=2)
        self.master.bind('r', self.reset_view)
        self.master.bind('R', self.reset_view)
        # --- Canvas ---
        self.canvas = tk.Canvas(master, width=self.width, height=self.height, bg='black', highlightthickness=0)
        self.canvas.pack(fill=tk.BOTH, expand=True)
        # --- Новая система навигации (как в Google Maps) ---
        self.focus_x = 0  # точка фокуса в мировых координатах
        self.focus_y = 0
        self.scale = 1.0
        self.MIN_SCALE = 0.01  # минимальный масштаб (видна вся галактика)
        self.MAX_SCALE = 100.0  # максимальный масштаб (детали планет)
        # --- Перемещение карты ---
        self._drag_active = False
        self._drag_start_screen = None  # экранные координаты точки начала перетаскивания
        self._drag_start_focus = None   # мировые координаты фокуса в начале перетаскивания
        self.canvas.bind('<ButtonPress-1>', self.start_drag)  # левая кнопка
        self.canvas.bind('<B1-Motion>', self.do_drag)
        self.canvas.bind('<ButtonRelease-1>', self.stop_drag)
        # --- Выбор объекта ---
        self.selected_body = None  # выбранный объект
        self.selected_body_id = None  # id выбранного объекта
        # --- Центрация на объекте ---
        self.center_on_selected_var = tk.BooleanVar(value=False)
        self.center_on_selected = False
        self.center_checkbutton_widget = None  # сам Checkbutton
        self.center_checkbutton_id = None  # id виджета на Canvas
        self.canvas.bind('<ButtonPress-3>', self.select_body)  # правая кнопка мыши
        # --- Галочка "Траектории относительно объекта" ---
        self.trails_relative_var = tk.BooleanVar(value=False)
        self.trails_relative = False
        # --- Симуляция ---
        self.bodies = []
        self.create_system()
        self.running = True
        self.sim_speed = 1.0  # множитель скорости отображения
        self.dt_physics = DT_PHYSICS
        self._physics_accum = 0.0  # накопитель для дробных шагов
        self.forces = [(0, 0) for _ in self.bodies]
        self.history = []
        self.history_index = 0
        self.save_state()

        # --- Выбор физического движка ---
        self.physics_engine = PhysicsEngine(algorithm="naive")  # или "barnes_hut"
        # Для теста: переключение алгоритма по клавише "b"
        self.master.bind('b', self.toggle_physics_algorithm)
        self.master.bind('B', self.toggle_physics_algorithm)
        self.canvas.bind('<MouseWheel>', self.on_mousewheel)
        # --- Управление клавиатурой ---
        self.master.bind('<KeyPress-plus>', self.increase_speed)
        self.master.bind('<KeyPress-minus>', self.decrease_speed)
        self.master.bind('<KeyPress-equal>', self.increase_speed)
        self.master.bind('<KeyPress-underscore>', self.decrease_speed)
        self.master.bind('<space>', self.toggle_pause)
        self.master.bind('<Left>', self.step_back)
        self.master.bind('<Right>', self.step_forward)
        self.master.bind('<BackSpace>', self.step_back)
        self.master.bind('<Configure>', self.on_resize)
        # --- Навигация клавиатурой ---
        self.master.bind('<Up>', self.move_up)
        self.master.bind('<Down>', self.move_down)
        self.master.bind('<KeyPress-0>', self.zoom_to_fit)
        self.master.bind('<KeyPress-1>', self.zoom_to_sun)
        self.master.bind('<Escape>', self.clear_selection)
        self.animate()

    def screen_to_world(self, screen_x, screen_y):
        """Экранные координаты → мировые координаты"""
        width = self.canvas.winfo_width()
        height = self.canvas.winfo_height()
        world_x = (screen_x - width/2) / self.scale + self.focus_x
        world_y = (screen_y - height/2) / self.scale + self.focus_y
        return world_x, world_y

    def world_to_screen(self, world_x, world_y):
        """Мировые координаты → экранные координаты"""
        width = self.canvas.winfo_width()
        height = self.canvas.winfo_height()
        screen_x = (world_x - self.focus_x) * self.scale + width/2
        screen_y = (world_y - self.focus_y) * self.scale + height/2
        return screen_x, screen_y

    def start_drag(self, event):
        if event.num == 1:  # только левая кнопка мыши
            self._drag_active = True
            self._drag_start_screen = (event.x, event.y)
            self._drag_start_focus = (self.focus_x, self.focus_y)
            self.canvas.config(cursor='hand2')

    def do_drag(self, event):
        if self._drag_active and self._drag_start_screen and self._drag_start_focus:
            # Вычисляем смещение в экранных координатах
            dx_screen = event.x - self._drag_start_screen[0]
            dy_screen = event.y - self._drag_start_screen[1]
            
            # Преобразуем экранное смещение в мировое смещение
            dx_world = dx_screen / self.scale
            dy_world = dy_screen / self.scale
            
            # Применяем смещение к фокусу
            self.focus_x = self._drag_start_focus[0] - dx_world
            self.focus_y = self._drag_start_focus[1] - dy_world

    def stop_drag(self, event):
        self._drag_active = False
        self._drag_start_screen = None
        self._drag_start_focus = None
        self.canvas.config(cursor='')

    def select_body(self, event):
        """Выбор объекта по клику правой кнопкой мыши"""
        # Сбросить центровку при выборе нового объекта
        self.center_on_selected_var.set(False)
        self.center_on_selected = False
        # Преобразуем экранные координаты в мировые
        world_x, world_y = self.screen_to_world(event.x, event.y)
        
        # Ищем ближайший объект в радиусе клика
        click_radius = 20 / self.scale  # радиус клика в мировых координатах
        closest_body = None
        min_distance = float('inf')
        
        for body in self.bodies:
            distance = math.sqrt((body.x - world_x)**2 + (body.y - world_y)**2)
            if distance < click_radius and distance < min_distance:
                min_distance = distance
                closest_body = body
        
        self.selected_body = closest_body
        self.selected_body_id = closest_body.id if closest_body else None

    def clear_selection(self, event=None):
        """Сброс выбора объекта"""
        self.selected_body = None
        self.selected_body_id = None
        self.center_on_selected_var.set(False)
        self.center_on_selected = False
        self.trails_relative_var.set(False)
        self.trails_relative = False
        self.hide_center_checkbutton()
        if hasattr(self, 'trails_relative_checkbutton_widget') and self.trails_relative_checkbutton_widget is not None:
            self.trails_relative_checkbutton_widget.place_forget()

    def hide_center_checkbutton(self):
        if self.center_checkbutton_id is not None:
            self.canvas.delete(self.center_checkbutton_id)
            self.center_checkbutton_id = None
        if self.center_checkbutton_widget is not None:
            self.center_checkbutton_widget.place_forget()

    def on_resize(self, event):
        self.width = event.width
        self.height = event.height
        self.draw()

    def get_center(self):
        """Возвращает центр экрана в мировых координатах"""
        return self.focus_x, self.focus_y

    def create_system(self):
        """Создает систему небесных тел в мировых координатах"""
        # Солнце
        sun = Body("Солнце", 0, 0, 0, 0, 10000, 'yellow')
        self.bodies.append(sun)
        # Планеты и их спутники
        planet_data = [
            (180, 10, 3),    # 1 — малая, быстрая, 3 спутника
            (260, 30, 2),    # 2 — массивная, 2 спутника
            (350, 18, 1),    # 3 — средняя, 1 спутник
            (470, 40, 3),    # 4 — очень массивная, 3 спутника
            (600, 8, 0),     # 5 — малая, без спутников
            (750, 25, 2),    # 6 — массивная, 2 спутника
            (900, 15, 1),    # 7 — средняя, 1 спутник
            (1100, 22, 4),   # 8 — крупная, 4 спутника
        ]
        for i, (orbit_r, mass, n_sat) in enumerate(planet_data):
            angle = math.radians(360 / len(planet_data) * i)
            px = orbit_r * math.cos(angle)
            py = orbit_r * math.sin(angle)
            v = math.sqrt(G * sun.m / orbit_r)
            vx = -v * math.sin(angle)
            vy = v * math.cos(angle)
            planet_name = f"Планета {i+1}"
            planet = Body(planet_name, px, py, vx, vy, mass, PLANET_COLORS[i % len(PLANET_COLORS)])
            self.bodies.append(planet)
            # Спутники
            if n_sat > 0:
                min_sat_r = random.randint(30, 60)
                sat_r_step = random.randint(10, 40)
                for j in range(n_sat):
                    sat_angle = angle + math.radians(360 / n_sat * j)
                    sat_r = min_sat_r + j * sat_r_step
                    sx = px + sat_r * math.cos(sat_angle)
                    sy = py + sat_r * math.sin(sat_angle)
                    v_sat = math.sqrt(G * planet.m / sat_r)
                    svx = planet.xv - v_sat * math.sin(sat_angle)
                    svy = planet.yv + v_sat * math.cos(sat_angle)
                    sat_name = f"Спутник {i+1}-{j+1}"
                    satellite = Body(sat_name, sx, sy, svx, svy, 1, SATELLITE_COLOR)
                    self.bodies.append(satellite)

    def save_state(self):
        # Сохраняет глубокую копию текущего состояния системы в историю
        state = copy.deepcopy(self.bodies)
        # Если отмотали назад и сделали новый шаг — обрезаем "будущее"
        if self.history_index < len(self.history) - 1:
            self.history = self.history[:self.history_index+1]
        self.history.append(state)
        self.history_index = len(self.history) - 1
        # Ограничим длину истории (например, 2000 шагов)
        if len(self.history) > 2000:
            self.history.pop(0)
            self.history_index -= 1

    def restore_state(self, index):
        # Восстанавливает состояние системы из истории по индексу
        self.bodies = copy.deepcopy(self.history[index])
        # Переназначаем выделенный объект по id
        if self.selected_body_id is not None:
            for body in self.bodies:
                if getattr(body, 'id', None) == self.selected_body_id:
                    self.selected_body = body
                    break
            else:
                self.selected_body = None
                self.selected_body_id = None
        else:
            self.selected_body = None

    def animate(self):
        # Если обе галочки включены — фокус камеры всегда (0, 0)
        if self.center_on_selected and self.selected_body and self.trails_relative:
            self.focus_x = 0
            self.focus_y = 0
        # Если только центровка — фокус на выбранном объекте
        elif self.center_on_selected and self.selected_body:
            self.focus_x = self.selected_body.x
            self.focus_y = self.selected_body.y
        if self.running:
            self._physics_accum += self.sim_speed
            steps = int(self._physics_accum)
            for _ in range(steps):
                self.update_physics()
                self.save_state()
                self._physics_accum -= 1.0
        self.draw()
        self.master.after(16, self.animate)

    def update_physics(self):
        self.forces = self.physics_engine.compute_forces(self.bodies, G)
        self.physics_engine.update_bodies(self.bodies, self.forces, self.dt_physics)
        self.iteration_counter += 1

    def draw(self):
        if not self.selected_body:
            self.hide_center_checkbutton()
            if hasattr(self, 'trails_relative_checkbutton_widget') and self.trails_relative_checkbutton_widget is not None:
                self.trails_relative_checkbutton_widget.place_forget()
        self.canvas.delete('all')
        y0 = self.top_frame.winfo_height() + 10
        line_h = 20
        info_bg_width = 420
        info_bg_height = 4 * line_h + 16
        # Фон под информационным блоком (чёрный с белой рамкой)
        self.canvas.create_rectangle(0, y0-6, info_bg_width, y0 + info_bg_height, fill='black', outline='white', width=2)

        # Масштаб и скорость
        self.canvas.create_text(16, y0+4, anchor='nw', fill='white', font=('Arial', 11, 'bold'),
                               text=f'Масштаб: {self.scale:.2f}x   Скорость: {self.sim_speed:.2f}x')
        # Статус "Пауза" справа, цветом
        if not self.running:
            self.canvas.create_text(info_bg_width-16, y0+4, anchor='ne', fill='yellow', font=('Arial', 11, 'bold'),
                                   text='Пауза')

        # Счетчик итераций
        self.canvas.create_text(16, y0 + line_h + 4, anchor='nw', fill='white', font=('Arial', 10),
                               text=f'Итерация: {self.iteration_counter}')
        # Информация о навигации
        nav_info = f'Фокус: ({self.focus_x:.0f}, {self.focus_y:.0f})'
        self.canvas.create_text(16, y0 + 2*line_h + 4, anchor='nw', fill='white', font=('Arial', 9),
                               text=nav_info)
        # Информация о масштабе векторов
        vector_info = f'Векторы: Скорость (2.0px/ед), Сила (0.1px/ед)'
        self.canvas.create_text(16, y0 + 3*line_h + 4, anchor='nw', fill='white', font=('Arial', 9),
                               text=vector_info)

        # --- Новое: отображение текущего алгоритма ---
        algo_name = "Barnes-Hut (O(N log N))" if self.physics_engine.algorithm == "barnes_hut" else "Наивный (O(N²))"
        self.canvas.create_text(info_bg_width + 20, y0+4, anchor='nw', fill='cyan', font=('Arial', 11, 'bold'),
                               text=f'Физика: {algo_name}')
        # --- Новое: подсказка ---
        self.canvas.create_text(info_bg_width + 20, y0 + line_h + 4, anchor='nw', fill='gray', font=('Arial', 9),
                               text='B — сменить физический алгоритм')
        
        # Следы (точки по мировым координатам, шаг TRAIL_POINT_STEP)
        if self.show_trails.get():  # Показывать только если включено
            rel_body = self.selected_body if self.trails_relative and self.selected_body else None
            for body in self.bodies:
                if len(body.trail) > 2:
                    for i in range(0, len(body.trail), TRAIL_POINT_STEP):
                        if rel_body and len(rel_body.trail) > i:
                            x = body.trail[i][0] - rel_body.trail[i][0]
                            y = body.trail[i][1] - rel_body.trail[i][1]
                        else:
                            x, y = body.trail[i]
                        x, y = self.world_to_screen(x, y)
                        self.canvas.create_oval(x-1, y-1, x+1, y+1, fill=body.Color, outline='')
        # Тела и векторы
        for idx, body in enumerate(self.bodies):
            if self.show_trails.get() and self.trails_relative and self.selected_body:
                rel_body = self.selected_body
                bx = body.x - rel_body.x
                by = body.y - rel_body.y
            else:
                bx, by = body.x, body.y
            x, y = self.world_to_screen(bx, by)
            r = max(2, int(body.r * self.scale))
            x1 = x - r
            y1 = y - r
            x2 = x + r
            y2 = y + r
            
            # Выделяем выбранный объект белой рамкой
            outline_color = 'white' if body == self.selected_body else ''
            outline_width = 3 if body == self.selected_body else 0
            
            self.canvas.create_oval(x1, y1, x2, y2, fill=body.Color, outline=outline_color, width=outline_width)
            # Векторы скорости и силы
            if getattr(self, 'show_vectors', True):
                # Вектор скорости (синий) - длина прямо пропорциональна модулю скорости
                vlen = math.hypot(body.xv, body.yv)
                if vlen > 0:
                    # Базовый масштаб для векторов скорости (пикселей на единицу скорости)
                    base_velocity_scale = 2.0
                    # Длина вектора прямо пропорциональна скорости
                    scale_v = base_velocity_scale * self.scale
                    vx_draw = body.xv * scale_v
                    vy_draw = body.yv * scale_v
                    self.draw_arrow(x, y, x + vx_draw, y + vy_draw, color='#00aaff')
                # Вектор силы (красный) - длина прямо пропорциональна модулю силы
                fx, fy = self.forces[idx]
                flen = math.hypot(fx, fy)
                if flen > 0:
                    # Базовый масштаб для векторов силы (пикселей на единицу силы)
                    base_force_scale = 5.0
                    # Длина вектора прямо пропорциональна силе
                    scale_f = base_force_scale * self.scale
                    fx_draw = fx * scale_f
                    fy_draw = fy * scale_f
                    self.draw_arrow(x, y, x + fx_draw, y + fy_draw, color='#ff3333')
        # Информация о выбранном объекте (правый верхний угол)
        if self.selected_body:
            self.draw_selected_body_info()
        
        # Подсказки по управлению
        controls = 'Управление: Мышь - перетаскивание, Колесо - масштаб, R - сброс, 0 - вся система, 1 - солнце'
        self.canvas.create_text(10, self.height - 20, anchor='sw', fill='white', font=('Arial', 8),
                               text=controls)
        
        # Добавляем подсказку о выборе объекта
        selection_hint = 'Правая кнопка мыши - выбрать объект, Escape - сброс выбора'
        self.canvas.create_text(10, self.height - 40, anchor='sw', fill='white', font=('Arial', 8),
                               text=selection_hint)

    def draw_arrow(self, x1, y1, x2, y2, color='white'):
        # Рисует стрелку от (x1, y1) к (x2, y2)
        self.canvas.create_line(x1, y1, x2, y2, fill=color, width=2, arrow=tk.LAST, arrowshape=(8,10,4))

    def draw_selected_body_info(self):
        """Отображает информацию о выбранном объекте в правом верхнем углу"""
        if not self.selected_body:
            self.hide_center_checkbutton()
            if hasattr(self, 'trails_relative_checkbutton_widget') and self.trails_relative_checkbutton_widget is not None:
                self.trails_relative_checkbutton_widget.place_forget()
            return
            
        # Безопасная проверка: есть ли объект в self.bodies
        if self.selected_body not in self.bodies:
            self.selected_body = None
            self.selected_body_id = None
            return

        body = self.selected_body
        
        # Определяем тип объекта
        if body == self.bodies[0]:  # Солнце всегда первое
            body_type = "Солнце"
        elif body.m > 5:  # Планеты имеют массу больше 5
            body_type = "Планета"
        else:
            body_type = "Спутник"
        
        # Вычисляем физические величины
        velocity = math.sqrt(body.xv**2 + body.yv**2)
        kinetic_energy = 0.5 * body.m * velocity**2
        
        # Находим индекс объекта для получения силы
        body_index = self.bodies.index(body)
        force_x, force_y = self.forces[body_index] if body_index < len(self.forces) else (0, 0)
        force_magnitude = math.sqrt(force_x**2 + force_y**2)
        
        # Создаем информационную панель
        text_lines = [
            f"ID: {body.id}",
            f"=== {body_type} ===",
            f"Позиция: ({body.x:.1f}, {body.y:.1f})",
            f"Скорость: ({body.xv:.2f}, {body.yv:.2f})",
            f"Модуль скорости: {velocity:.2f}",
            f"Масса: {body.m:.1f}",
            f"Радиус: {body.r:.1f}",
            f"Кинетическая энергия: {kinetic_energy:.1f}",
            f"Сила: ({force_x:.1f}, {force_y:.1f})",
            f"Модуль силы: {force_magnitude:.1f}",
            f"Цвет: {body.Color}"
        ]
        
        # Вычисляем размеры панели
        max_width = max(len(line) for line in text_lines) * 8  # примерная ширина символа
        panel_height = len(text_lines) * 20 + 10
        canvas_width = self.canvas.winfo_width()
        panel_x = canvas_width - max_width - 20  # 10px отступ справа
        panel_y = 10
        
        # Рисуем фон
        self.canvas.create_rectangle(
            panel_x, panel_y - 5,
            panel_x + max_width + 10, panel_y + panel_height + 30,
            fill='black', outline='white', width=2
        )
        
        # Рисуем текст
        for i, line in enumerate(text_lines):
            y_pos = panel_y + i * 20
            self.canvas.create_text(
                panel_x + 5, y_pos,
                anchor='nw', fill='white', font=('Arial', 9),
                text=line
            )
        # --- Галочка "Центрация на объекте" ---
        def on_center_check():
            self.center_on_selected = self.center_on_selected_var.get()
            if self.center_on_selected and self.selected_body:
                self.focus_x = self.selected_body.x
                self.focus_y = self.selected_body.y
                self.draw()
        # --- Галочка "Траектории относительно объекта" ---
        def on_trails_relative_check():
            self.trails_relative = self.trails_relative_var.get()
            self.draw()
        # Создаём Checkbutton только один раз
        if self.center_checkbutton_widget is None:
            cb = tk.Checkbutton(self.canvas, text='Центрация на объекте', variable=self.center_on_selected_var,
                                command=on_center_check, fg='white', bg='black', selectcolor='black',
                                activebackground='black', activeforeground='white', font=('Arial', 9))
            self.center_checkbutton_widget = cb
        if not hasattr(self, 'trails_relative_checkbutton_widget') or self.trails_relative_checkbutton_widget is None:
            cb2 = tk.Checkbutton(self.canvas, text='Траектории относительно объекта', variable=self.trails_relative_var,
                                 command=on_trails_relative_check, fg='white', bg='black', selectcolor='black',
                                 activebackground='black', activeforeground='white', font=('Arial', 9))
            self.trails_relative_checkbutton_widget = cb2
        # Размещение чекбоксов друг под другом
        x = panel_x + 5
        y = panel_y + panel_height - 5
        self.center_checkbutton_widget.place(x=x, y=y)
        self.center_on_selected_var.set(self.center_on_selected)
        x2 = x
        y2 = y + 28  # аккуратный вертикальный отступ
        self.trails_relative_checkbutton_widget.place(x=x2, y=y2)
        self.trails_relative_var.set(self.trails_relative)

    # --- Управление ---
    def to_screen(self, pos):
        """Совместимость со старым кодом - использует новую систему координат"""
        x, y = pos
        return self.world_to_screen(x, y)

    def on_mousewheel(self, event):
        """Масштабирование относительно курсора (как в Google Maps)"""
        if getattr(self, '_drag_active', False):
            return  # Запрещаем масштабирование во время drag

        zoom_factor = 1.1 if event.delta > 0 else 1/1.1

        mouse_screen_x = event.x
        mouse_screen_y = event.y

        old_focus_x = self.focus_x
        old_focus_y = self.focus_y
        old_scale = self.scale

        new_scale = old_scale * zoom_factor
        new_scale = max(self.MIN_SCALE, min(new_scale, self.MAX_SCALE))
        if abs(new_scale - old_scale) < 1e-6:
            return

        width = self.canvas.winfo_width()
        height = self.canvas.winfo_height()
        # 1. Мировые координаты точки под курсором до масштабирования
        world_x = (mouse_screen_x - width/2) / old_scale + old_focus_x
        world_y = (mouse_screen_y - height/2) / old_scale + old_focus_y

        # 2. Применяем новый масштаб
        self.scale = new_scale

        # 3. Новый фокус так, чтобы world_x/world_y остались под курсором
        self.focus_x = world_x - (mouse_screen_x - width/2) / new_scale
        self.focus_y = world_y - (mouse_screen_y - height/2) / new_scale

        # Удалена отладочная визуализация точки масштабирования

    def increase_speed(self, event=None):
        if self.sim_speed < 1.0:
            self.sim_speed += 0.1
            if self.sim_speed > 1.0:
                self.sim_speed = 1.0
        elif self.sim_speed < 10.0:
            self.sim_speed += 1.0
            if self.sim_speed > 10.0:
                self.sim_speed = 10.0
        elif self.sim_speed < 100.0:
            self.sim_speed += 10.0
            if self.sim_speed > 100.0:
                self.sim_speed = 100.0

    def decrease_speed(self, event=None):
        if self.sim_speed > 10.0:
            self.sim_speed -= 10.0
            if self.sim_speed < 10.0:
                self.sim_speed = 10.0
        elif self.sim_speed > 1.0:
            self.sim_speed -= 1.0
            if self.sim_speed < 1.0:
                self.sim_speed = 1.0
        elif self.sim_speed > 0.1:
            self.sim_speed -= 0.1
            if self.sim_speed < 0.1:
                self.sim_speed = 0.1

    def toggle_pause(self, event=None):
        self.running = not self.running

    def step_back(self, event=None):
        if self.history_index > 0:
            self.history_index -= 1
            self.restore_state(self.history_index)
            self.running = False  # ставим на паузу при перемотке

    def step_forward(self, event=None):
        if self.history_index < len(self.history) - 1:
            self.history_index += 1
            self.restore_state(self.history_index)
            self.running = False  # ставим на паузу при перемотке

    def toggle_vectors(self):
        self.show_vectors = self.show_vectors_var.get()

    def toggle_trails(self):
        # Можно добавить дополнительные действия при переключении, если нужно
        pass

    def reset_view(self, event=None):
        """Сброс вида к начальному состоянию"""
        self.focus_x = 0
        self.focus_y = 0
        self.scale = 1.0

    def move_up(self, event=None):
        """Перемещение вверх"""
        self.focus_y -= 100 / self.scale

    def move_down(self, event=None):
        """Перемещение вниз"""
        self.focus_y += 100 / self.scale

    def zoom_to_fit(self, event=None):
        """Масштабирование для показа всей системы"""
        if not self.bodies:
            return
            
        # Находим границы системы
        min_x = min(body.x for body in self.bodies)
        max_x = max(body.x for body in self.bodies)
        min_y = min(body.y for body in self.bodies)
        max_y = max(body.y for body in self.bodies)
        
        # Добавляем отступы
        padding = 100
        width = max_x - min_x + 2 * padding
        height = max_y - min_y + 2 * padding
        
        # Вычисляем масштаб
        scale_x = self.width / width
        scale_y = self.height / height
        self.scale = min(scale_x, scale_y, self.MAX_SCALE)
        
        # Центрируем систему
        self.focus_x = (min_x + max_x) / 2
        self.focus_y = (min_y + max_y) / 2

    def zoom_to_sun(self, event=None):
        """Масштабирование к солнцу"""
        if self.bodies:
            sun = self.bodies[0]  # Солнце всегда первое
            self.focus_x = sun.x
            self.focus_y = sun.y
            self.scale = 2.0  # Увеличенный масштаб для солнца

    def toggle_physics_algorithm(self, event=None):
        if self.physics_engine.algorithm == "naive":
            self.physics_engine.algorithm = "barnes_hut"
        else:
            self.physics_engine.algorithm = "naive"
        print(f"Сменён алгоритм физики: {self.physics_engine.algorithm}")

if __name__ == '__main__':
    root = tk.Tk()
    app = GravitySim(root)
    root.mainloop() 