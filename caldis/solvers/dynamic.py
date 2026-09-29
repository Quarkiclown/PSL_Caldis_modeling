import numpy as np

from caldis.solvers.steady import solve_steady
from caldis.solvers.newton import newton
from caldis.experiments.results import Results


def simulate(system, t_end, dt, update_inputs=None, init="steady", tol=1e-8):
    if not system._assembled:
        system.assemble()
    system.check()
    states = [v for v in system.unknowns if v.differential]

    if update_inputs is not None:
        update_inputs(system, 0.0)
    if init == "steady":
        sol0 = solve_steady(system)
        if not sol0.success:
            raise RuntimeError(f"Initialisation stationnaire échouée : {sol0.message}")
    for v in states:
        v.der = 0.0

    system.residuals()
    rec0 = system.collect_outputs()
    hist = {k: [v] for k, v in rec0.items()}
    times = [0.0]
    x_prev = {v.name: v.value for v in states}

    n_steps = int(round(t_end / dt))
    for k in range(1, n_steps + 1):
        t = k * dt
        if update_inputs is not None:
            update_inputs(system, t)

        def G(x):
            system.set_x(x)
            for v in states:
                v.der = (v.value - x_prev[v.name]) / dt
            return system.residuals()

        guess = np.array([v.value for v in system.unknowns], dtype=float)
        sol = newton(G, guess, ftol=tol, xtol=tol)
        if not sol.success:
            raise RuntimeError(f"Pas t={t:.4g}s : {sol.message}")
        system.set_x(sol.x)
        for v in states:
            v.der = (v.value - x_prev[v.name]) / dt

        system.residuals()
        for kk, val in system.collect_outputs().items():
            hist[kk].append(val)
        times.append(t)
        x_prev = {v.name: v.value for v in states}

    return Results(np.array(times), {kk: np.array(vv) for kk, vv in hist.items()})