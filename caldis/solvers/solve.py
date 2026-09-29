from caldis.core.dof import check_dof, check_startpoint
from caldis.solvers.homotopy import solve_homotopy


def solve(system, steps=10, verbose=True, show_homotopy=False,
          diagnostics=True, raise_on_fail=False):
    """Full solve workflow, ordered cheap-structural -> expensive-numerical.
    Each stage prints a one-line OK on success, or detailed diagnostics on failure.
    """
    if not system._assembled:
        system.assemble()

    # -- Stage 1: DOF (structural). Cheapest check; a non-square system can't be
    #    solved, so stop here. Note: n_equations() calls residuals() once, which is
    #    why residuals() MUST be evaluation-safe (return NaN, never raise).
    dof = check_dof(system)
    if not dof.ok:
        if verbose:
            print(f"✗ DOF: {dof}")
            print("  -> the model is not square. Add/remove a constraint or boundary "
                  "spec before solving. (This is a MODELLING issue, not a start issue.)")
        if raise_on_fail:
            raise ValueError(str(dof))
        return None
    if verbose:
        print(f"✓ DOF: square ({dof.n_unknowns} unknowns)")

    # -- Stage 2: start evaluable? (can every residual be computed at the start?)
    #    This catches out-of-domain (p,h), degenerate DTLM (all-equal T), etc.,
    #    and pinpoints the failing residual + port states.
    try:
        check_startpoint(system)
        if verbose:
            print("✓ start point: evaluable")
    except ValueError as e:
        if verbose:
            print("✗ start point NOT evaluable:")
            print("   " + str(e).replace("\n", "\n   "))
            print("  -> add start hints (set_start_port) for the ports shown above. "
                  "Most often: give distinct inlet/outlet temperatures so DTLM is defined.")
        if raise_on_fail:
            raise
        return None

    # -- Stage 3: physical preconditions (warnings only; non-blocking).
    warnings = system.check_start_sanity(verbose=False)
    if verbose and warnings:
        print(f"⚠ {len(warnings)} start-sanity warning(s):")
        for w in warnings:
            print("   - " + w)
        print("  (non-blocking: the solver may still converge; shown in case it fails.)")

    # -- Stage 4: solve (expensive). Homotopy trace muted unless show_homotopy.
    sol = solve_homotopy(system, steps=steps, verbose=show_homotopy)

    if sol.success:
        if verbose:
            fs = getattr(sol, "res_scaled", None)
            acc = f"|F_scaled|={fs:.1e}" if fs is not None else "converged"
            print(f"✓ SOLVED ({acc})")
    else:
        if verbose:
            print(f"✗ FAILED: {sol.message}")
        if diagnostics:
            from caldis.core.diagnostics import convergence_report, jacobian_report
            convergence_report(system, sol)
            if "singul" in sol.message.lower():
                jacobian_report(system)
            if verbose:
                print("  -> re-anchor the variable(s) named above (set_start_port), "
                      "then solve again.")
        if raise_on_fail:
            raise RuntimeError(f"solve failed: {sol.message}")
    
    # ensure the system state always reflects sol.x, so a user's manual
    # inspection after solve() matches what the report showed (success or fail).
    # ... après la résolution et le convergence_report éventuel ...
    system.set_x(sol.x)                      # state reflects the last solution
    for c in system.components:
        _update_outputs_recursive(c)
    return sol


def _update_outputs_recursive(comp):
    """Fill outputs on a component and its children (composites)."""
    try:
        comp.update_outputs()
    except Exception:
        pass                                 # a degenerate stop state may be uncomputable
    for child in comp.children:
        _update_outputs_recursive(child)
