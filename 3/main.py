import math
from collections.abc import Callable
from dataclasses import dataclass


@dataclass
class segment:
    a: float
    b: float


@dataclass
class equation:
    """Клас для нелінійного рівняння f(x) = 0 або x = phi(x)"""
    f: Callable[[float], float]
    df: Callable[[float], float] | None = None
    d2f: Callable[[float], float] | None = None
    phi: Callable[[float], float] | None = None
    dphi: Callable[[float], float] | None = None
    seg: segment | None = None
    eps: float = 1e-4
    x0: float | None = None


@dataclass
class equationSystem:
    """Клас для системи рівнянь F(x) = 0 або x = Phi(x)"""
    F: Callable[[list[float]], list[float]]
    J: Callable[[list[float]], list[list[float]]] | None = None
    Phi: Callable[[list[float]], list[float]] | None = None
    J_Phi: Callable[[list[float]], list[list[float]]] | None = None
    x0: list[float] | None = None
    eps: float = 1e-4


def norm_inf(v: list[float]) -> float:
    """Нескінченна норма вектора: ||v||_inf = max |v_i|"""
    return max(abs(x) for x in v)


def norm_inf_matrix(A: list[list[float]]) -> float:
    """Нескінченна норма матриці: ||A||_inf = max_i sum_j |A_ij|"""
    return max(sum(abs(x) for x in row) for row in A)


def solve_linear_system(A: list[list[float]], b: list[float]) -> list[float]:
    """Метод Гаусса для розв'язання A * delta = b на чистому Python"""
    n = len(b)
    M = [A[i][:] + [b[i]] for i in range(n)]

    for i in range(n):
        max_row = max(range(i, n), key=lambda r: abs(M[r][i]))
        if abs(M[max_row][i]) < 1e-14:
            raise ValueError("Матриця Якобі вироджена (det J(x) = 0)!")
        M[i], M[max_row] = M[max_row], M[i]

        pivot = M[i][i]
        for col in range(i, n + 1):
            M[i][col] /= pivot

        for r in range(n):
            if r != i:
                factor = M[r][i]
                for c in range(i, n + 1):
                    M[r][c] -= factor * M[i][c]

    return [M[i][n] for i in range(n)]


def simple_iter_method(eq: equation | equationSystem):
    if isinstance(eq, equation):
        _simple_iter_equation(eq)
    elif isinstance(eq, equationSystem):
        _simple_iter_system(eq)
    else:
        raise TypeError("Невідомий тип рівняння!")


def _simple_iter_equation(eq: equation):
    print("=" * 80)
    print("МЕТОД ПРОСТОЇ ІТЕРАЦІЇ ДЛЯ ОДНОГО РІВНЯННЯ")
    print("=" * 80)

    if eq.phi is None:
        raise ValueError("Для МПІ необхідно задати функцію phi(x)!")

    # обчислення q = max |phi'(x)| на відрізку [a, b] та перевірка q < 1
    if eq.dphi is not None and eq.seg is not None:
        steps = 500
        step_sz = (eq.seg.b - eq.seg.a) / steps
        grid = [eq.seg.a + i * step_sz for i in range(steps + 1)]
        q = max(abs(eq.dphi(x)) for x in grid)

        print(f"Відрізок ізоляції: [{eq.seg.a}, {eq.seg.b}]")
        print(f"Параметр стиску q = max |phi'(x)|: {q:.6f}")
        if q >= 1.0:
            print(f"УВАГА: q = {q:.4f} >= 1, умова збіжності МПІ не виконується!")
            return
        else:
            print(f"Перевірка q < 1: ВИКОНАНО ({q:.4f} < 1)")
    else:
        q = 0.5

    # початкове наближення
    if eq.x0 is not None:
        x_prev = eq.x0
    elif eq.seg is not None:
        x_prev = (eq.seg.a + eq.seg.b) / 2.0
    else:
        x_prev = 0.0

    print(f"Початкове наближення x_0 = {x_prev}")
    print("-" * 80)
    print(f"{'Крок (k)':<10} | {'x^(k)':<18} | {'|x^(k) - x^(k-1)|':<22} | {'Оцінка похибки':<20}")
    print("-" * 80)

    k = 0
    while True:
        k += 1
        x_curr = eq.phi(x_prev)
        diff = abs(x_curr - x_prev)
        err_est = (q / (1.0 - q)) * diff  # Критерій з Фото 1 (п. 4)

        print(f"{k:<10} | {x_curr:<18.8f} | {diff:<22.8e} | {err_est:<20.8e}")

        if err_est < eq.eps:
            break
        x_prev = x_curr

    residual = abs(eq.f(x_curr))
    print("-" * 80)
    print(f"Наближений розв'язок x* : {x_curr:.8f}")
    print(f"Кількість ітерацій     : {k}")
    print(f"Нев'язка |f(x*)|       : {residual:.8e}")
    print("=" * 80 + "\n")


def _simple_iter_system(eq: equationSystem):
    print("=" * 80)
    print("МЕТОД ПРОСТОЇ ІТЕРАЦІЇ ДЛЯ СИСТЕМИ")
    print("=" * 80)

    if eq.Phi is None or eq.x0 is None:
        raise ValueError("Для системи необхідно задати Phi(x) та x0!")

    # обчислення q = ||Phi'(x0)||_inf в околі початкового наближення та перевірка q < 1
    if eq.J_Phi is not None:
        J_val = eq.J_Phi(eq.x0)
        q = norm_inf_matrix(J_val)
        print(f"Початкове наближення x_0: {eq.x0}")
        print("Матриця Якобі Phi'(x0):")
        for row in J_val:
            print(f"  {[round(val, 6) for val in row]}")
        print(f"Параметр стиску q = ||Phi'(x0)||_inf: {q:.6f}")
        if q >= 1.0:
            print(f"УВАГА: q = {q:.4f} >= 1, умова збіжності не гарантується!")
            return
        else:
            print(f"Перевірка q < 1: ВИКОНАНО ({q:.4f} < 1)")
    else:
        q = 0.5

    x_prev = list(eq.x0)
    k = 0

    print("-" * 80)
    print(f"{'Крок (k)':<10} | {'x^(k)':<32} | {'||x^(k) - x^(k-1)||_inf':<24} | {'Оцінка похибки':<18}")
    print("-" * 80)

    while True:
        k += 1
        x_curr = eq.Phi(x_prev)
        diff = norm_inf([x_curr[i] - x_prev[i] for i in range(len(x_curr))])
        err_est = (q / (1.0 - q)) * diff  # Зупинка за критерієм: q / (1 - q) * ||...|| < eps

        x_str = "[" + ", ".join(f"{val:11.8f}" for val in x_curr) + "]"
        print(f"{k:<10} | {x_str:<32} | {diff:<24.8e} | {err_est:<18.8e}")

        if err_est < eq.eps:
            break
        x_prev = x_curr

    res_vec = eq.F(x_curr)
    residual = norm_inf(res_vec)
    print("-" * 80)
    print(f"Наближений вектор x*   : {[round(val, 8) for val in x_curr]}")
    print(f"Кількість ітерацій     : {k}")
    print(f"Нев'язка ||F(x*)||_inf : {residual:.8e}")
    print("=" * 80 + "\n")


def newton_method(eq: equation | equationSystem):
    if isinstance(eq, equation):
        _newton_equation(eq)
    elif isinstance(eq, equationSystem):
        _newton_system(eq)
    else:
        raise TypeError("Невідомий тип рівняння!")


def _newton_equation(eq: equation):
    print("=" * 80)
    print("МЕТОД НЬЮТОНА ДЛЯ ОДНОГО РІВНЯННЯ")
    print("=" * 80)

    if eq.df is None:
        raise ValueError("Для методу Ньютона потрібна похідна df!")

    # вибір початкового наближення за умовою Фур'є: f(x0) * f''(x0) > 0
    x0 = eq.x0
    if x0 is None and eq.seg is not None and eq.d2f is not None:
        f_a_cond = eq.f(eq.seg.a) * eq.d2f(eq.seg.a)
        f_b_cond = eq.f(eq.seg.b) * eq.d2f(eq.seg.b)
        print("Перевірка умови збіжності f(x0)*f''(x0) > 0:")
        print(f"  f(a)*f''(a) = {f_a_cond:.4f}")
        print(f"  f(b)*f''(b) = {f_b_cond:.4f}")

        if f_a_cond > 0:
            x0 = eq.seg.a
            print(f"Обрано x0 = a = {x0} (f(a)*f''(a) > 0)")
        elif f_b_cond > 0:
            x0 = eq.seg.b
            print(f"Обрано x0 = b = {x0} (f(b)*f''(b) > 0)")
        else:
            x0 = (eq.seg.a + eq.seg.b) / 2.0
            print(f"Обрано x0 = {x0}")
    elif x0 is None and eq.seg is not None:
        x0 = eq.seg.b

    x_prev = x0
    print(f"Початкове наближення: x0 = {x_prev}")
    print("-" * 80)
    print(f"{'Крок (k)':<10} | {'x^(k)':<18} | {'|x^(k) - x^(k-1)|':<22} | {'Оцінка похибки':<20}")
    print("-" * 80)

    k = 0
    while True:
        k += 1
        df_val = eq.df(x_prev)
        if abs(df_val) < 1e-14:
            raise ZeroDivisionError(f"Похідна f'(x) близька до 0 на кроці {k}!")

        # delta = - f(x) / f'(x)
        delta_x = -eq.f(x_prev) / df_val
        x_curr = x_prev + delta_x
        diff = abs(delta_x)
        err_est = diff

        print(f"{k:<10} | {x_curr:<18.8f} | {diff:<22.8e} | {err_est:<20.8e}")

        # Критерій зупинки: ||x^(k) - x^(k-1)|| < eps
        if diff < eq.eps:
            break
        x_prev = x_curr

    residual = abs(eq.f(x_curr))
    print("-" * 80)
    print(f"Наближений розв'язок x* : {x_curr:.8f}")
    print(f"Кількість ітерацій     : {k}")
    print(f"Нев'язка |f(x*)|       : {residual:.8e}")
    print("=" * 80 + "\n")


def _newton_system(eq: equationSystem):
    print("=" * 80)
    print("МЕТОД НЬЮТОНА ДЛЯ СИСТЕМИ")
    print("=" * 80)

    if eq.J is None or eq.x0 is None:
        raise ValueError("Для системи потрібні F, J та x0!")

    x_prev = list(eq.x0)
    print(f"Початкове наближення x0: {x_prev}")
    print("-" * 80)
    print(f"{'Крок (k)':<10} | {'x^(k)':<32} | {'||x^(k) - x^(k-1)||_inf':<24} | {'Оцінка похибки':<18}")
    print("-" * 80)

    k = 0
    while True:
        k += 1
        J_val = eq.J(x_prev)
        F_val = eq.F(x_prev)
        minus_F = [-val for val in F_val]

        # J(x^(k-1)) * delta = -F(x^(k-1)) з перевіркою det J != 0
        delta = solve_linear_system(J_val, minus_F)
        x_curr = [x_prev[i] + delta[i] for i in range(len(x_prev))]

        diff = norm_inf(delta)
        err_est = diff

        x_str = "[" + ", ".join(f"{val:11.8f}" for val in x_curr) + "]"
        print(f"{k:<10} | {x_str:<32} | {diff:<24.8e} | {err_est:<18.8e}")

        if diff < eq.eps:
            break
        x_prev = x_curr

    residual = norm_inf(eq.F(x_curr))
    print("-" * 80)
    print(f"Наближений розв'язок   : {[round(val, 8) for val in x_curr]}")
    print(f"Кількість ітерацій     : {k}")
    print(f"Нев'язка ||F(x*)||_inf : {residual:.8e}")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    # Завдання 1: Метод Ньютона
    # Рівняння: e^x - x - 3 = 0 на відрізку [1; 2]
    # Точність: eps = 10^-6
    print(">>> ВАРІАНТ 15. 1: Метод Ньютона для e^x - x - 3 = 0")
    eq_task1 = equation(
        f=lambda x: math.exp(x) - x - 3.0,
        df=lambda x: math.exp(x) - 1.0,
        d2f=lambda x: math.exp(x),
        seg=segment(a=1.0, b=2.0),
        eps=1e-6
    )
    newton_method(eq_task1)

    # 2: Метод простої ітерації для системи
    # 3x - cos(x + y) = 1      => x = (1 + cos(x + y)) / 3
    # 2y + sin(x - y) = 0.5    => y = (0.5 - sin(x - y)) / 2
    # Початкове наближення: (0.6; -0.1)
    # Точність: eps = 10^-5
    print(">>> ВАРІАНТ 15. 2: Метод простої ітерації для системи")
    eq_task2 = equationSystem(
        F=lambda v: [
            3.0 * v[0] - math.cos(v[0] + v[1]) - 1.0,
            2.0 * v[1] + math.sin(v[0] - v[1]) - 0.5
        ],
        Phi=lambda v: [
            (1.0 + math.cos(v[0] + v[1])) / 3.0,
            (0.5 - math.sin(v[0] - v[1])) / 2.0
        ],
        J_Phi=lambda v: [
            [-math.sin(v[0] + v[1]) / 3.0, -math.sin(v[0] + v[1]) / 3.0],
            [-math.cos(v[0] - v[1]) / 2.0,  math.cos(v[0] - v[1]) / 2.0]
        ],
        x0=[0.6, -0.1],
        eps=1e-5
    )
    simple_iter_method(eq_task2)
