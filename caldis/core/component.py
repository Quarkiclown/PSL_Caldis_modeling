from caldis.core.port import FluidPort
from caldis.core.variable import Variable


class Component:
    """Composant = contributeur de résidus, éventuellement composite (avec enfants).

    Nommage : nom LOCAL + parent -> chemin calculé (comp.enfant.petitenfant.var).
    Un composite EST un Component ordinaire ; le système ne voit qu'une liste plate.
    """

    def __init__(self, name):
        self.local_name = name
        self.parent = None
        self.ports = []
        self._internal = []
        self.children = []
        self.params = {}
        self._cache = {}                # dynamic sim: values cached in residuals(), read step by step
        self._outputs = {}              # steady derived outputs, filled by update_outputs()
        # --- dynamic-simulation output selection (used by collect_outputs / Results) ---
        self.save_ports = True          # cat.1: save all port (p,h,mdot); can be disabled
        self.outputs = []               # cat.3: own variables/derived to save
        self.saved_children = None      # cat.2: None = all children; else list of local names

    # --- chemin / nom --------------------------------------------------------
    @property
    def path(self):
        if self.parent is None:
            return self.local_name
        return f"{self.parent.path}.{self.local_name}"

    @property
    def name(self):                     # rétro-compat : identique au chemin (== local si racine)
        return self.path

    def add_port(self, name, backend=None):
        import keyword
        if not name.isidentifier() or keyword.iskeyword(name):
            raise ValueError(
                f"invalid port name {name!r}: must be a non-reserved Python identifier")
        port = FluidPort(self, name)
        self.ports.append(port)
        port.backend = backend if backend is not None else self._backend_for_port(port)
        return port
   
    def add_variable(self, local_name, start=1.0, differential=False,
                     lower=None, upper=None, scale=1.0):
        var = Variable(self, local_name, start=start, differential=differential,
                       lower=lower, upper=upper, scale=scale)
        self._internal.append(var)
        return var

    def add_child(self, child):
        if child.parent is not None:
            raise ValueError(
                f"{child.local_name!r} a déjà un parent ({child.parent.local_name!r})")
        if child.local_name in [c.local_name for c in self.children]:
            raise ValueError(
                f"{self.local_name!r} a déjà un enfant nommé {child.local_name!r}")
        child.parent = self
        self.children.append(child)
        return child

    # --- assemblage (récursif) ----------------------------------------------
    def all_variables(self):
        vs = list(self._internal)
        for port in self.ports:
            vs += port.variables()
        for child in self.children:
            vs += child.all_variables()
        return vs

    # --- résidus -------------------------------------------------------------
    def residuals(self):
        """Résidus propres + ceux des enfants (récursif). Les FEUILLES surchargent
        généralement residuals() directement ; les COMPOSITES surchargent
        _own_residuals() et laissent cette récursion agir."""
        r = list(self._own_residuals())
        for child in self.children:
            r += list(child.residuals())
        return r

    def _own_residuals(self):
        return []

    # --- sorties -------------------------------------------------------------

    def collect_outputs(self, into):
        """Dynamic sim: record port states and derived outputs at the current step.
        Uses update_outputs() as the single source for derived quantities."""
        if self.save_ports:
            for port in self.ports:
                into[f"{port.path}.p"] = port.p
                into[f"{port.path}.mdot"] = port.mdot
                into[f"{port.path}.h"] = port.h
        self.update_outputs()                     # fill _outputs from current state
        for name, val in self._outputs.items():
            into[f"{self.path}.{name}"] = val
        for v in self._internal:                  # internal variables (M, U, T...)
            into[f"{self.path}.{v.local_name}"] = v.value
        for child in self.children:
            if self.saved_children is None or child.local_name in self.saved_children:
                child.collect_outputs(into)
        return into
    # --- garde-fous (récursif) ----------------------------------------------
    def check(self):
        """Recurse into children (no output-list validation anymore)."""
        for child in self.children:
            child.check()

    def __repr__(self):

        return f"{type(self).__name__}({self.path!r})"


    def _backend_for_port(self, port):
        """Legacy fallback: infer a port's backend from its name prefix
        (r_* -> backend_r, s_* -> backend_s, else backend). Prefer passing the
        backend explicitly to add_port(name, backend=...)."""
        name = port.name
        if name.startswith("r_") and getattr(self, "backend_r", None) is not None:
            return self.backend_r
        if name.startswith("s_") and getattr(self, "backend_s", None) is not None:
            return self.backend_s
        return getattr(self, "backend", None)

    def apply_default_starts(self, T_hint):
        """Level-1 fluid-aware default: in-domain (p, h) for each fluid port."""
        for port in self.ports:
            be = port.backend
            if be is None:
                continue
            p, h = be.default_state(T_hint)
            port._p.start = p            # write on the Variable OBJECT
            port._h.start = h
        for child in self.children:
            child.apply_default_starts(T_hint)

    def start_sanity(self):
        """Return a list of (ok: bool, message: str) checking this component's own
        physical preconditions at the current start point (condition 2: equation
        validity, beyond mere CoolProp evaluability). Default: nothing to check.
        Recurses into children."""
        checks = list(self._own_start_sanity())
        for child in self.children:
            checks += child.start_sanity()
        return checks

    def _own_start_sanity(self):
        return []



    def describe(self, verbose=True):
        """Readable dump after solve: internal variables, port states, derived outputs."""
        info = {"variables": {}, "ports": {}, "outputs": dict(self._outputs)}
        for v in self._internal:
            info["variables"][v.local_name] = v.value
        for port in self.ports:
            info["ports"][port.name] = {"p": port.p, "mdot": port.mdot, "h": port.h}
        if verbose:
            print(f"=== {type(self).__name__}({self.path!r}) ===")
            if info["variables"]:
                print("  variables:")
                for k, val in info["variables"].items():
                    print(f"    {k:14s} = {val:.6g}")
            print("  ports (p [bar], mdot [kg/s], T [°C]):")
            for port in self.ports:
                tinfo = ""
                if port.backend is not None:
                    try: tinfo = f"  T={port.T-273.15:.2f}"
                    except Exception: tinfo = "  T=n/a"
                print(f"    {port.name:10s} p={port.p/1e5:.3f}  mdot={port.mdot:+.4g}"
                      f"  h={port.h/1e3:.1f}kJ/kg{tinfo}")
            if info["outputs"]:
                print("  outputs:")
                for k, val in info["outputs"].items():
                    sval = f"{val:.6g}" if isinstance(val, (int, float)) else str(val)
                    print(f"    {k:14s} = {sval}")
            elif not self.children:
                print("  (no outputs yet -- run solve())")
            for child in self.children:
                child.describe(verbose=True)
        return info



    def structure(self, depth=None, _level=0):
        """Print structure: internal variables, ports, sub-components. Derived output
        names appear only after a solve (they are computed by update_outputs)."""
        pad = "  " * _level
        print(f"{pad}{type(self).__name__} '{self.local_name}'")
        for v in self._internal:
            print(f"{pad}  var    {v.local_name}"
                  f"{'  [state]' if v.differential else ''}")
        for port in self.ports:
            be = port.backend
            info = f"  ({be.fluid})" if (be is not None and hasattr(be, "fluid")) \
                   else "  (NO backend: .T/.x/.s unavailable)"
            print(f"{pad}  port   {port.name}{info}   -> .p .mdot .h")
        if self._outputs:
            print(f"{pad}  outputs {sorted(self._outputs)}")
        if depth is None or _level < depth:
            for child in self.children:
                child.structure(depth=depth, _level=_level + 1)
        elif self.children:
            print(f"{pad}  children (not expanded): {[c.local_name for c in self.children]}")
        return None

    def entropy_generation(self):
        """Sgen [W/K], adiabatic: - sum_ports (mdot * s)."""
        sgen = 0.0
        for port in self.ports:
            be = port.backend
            if be is None:
                continue
            s = be.s(port.p, port.h)
            sgen -= port.mdot * s
        return sgen

    def exergy_balance(self, T0, p0):
        X_destroyed = T0 * self.entropy_generation()
        self.X_supplied = self.X_recovered = self.epsilon_II = self.exergy_check = None
        self.X_destroyed, self.T_ref, self.p_ref = X_destroyed, T0, p0
        return {"X_supplied": None, "X_recovered": None,
                "X_destroyed": X_destroyed, "epsilon_II": None, "check": None}
    

    def _finish_balance(self, T0, p0, X_supplied, X_recovered):
        X_destroyed = T0 * self.entropy_generation()
        eps = X_recovered / X_supplied if abs(X_supplied) > 1e-12 else float("nan")
        check = ((X_supplied - X_recovered) - X_destroyed) / X_destroyed \
                if abs(X_destroyed) > 1e-12 else float("nan")
        self.X_supplied, self.X_recovered, self.X_destroyed = X_supplied, X_recovered, X_destroyed
        self.epsilon_II, self.exergy_check = eps, check
        self.T_ref, self.p_ref = T0, p0          # both reference values recorded
        return {"X_supplied": X_supplied, "X_recovered": X_recovered,
                "X_destroyed": X_destroyed, "epsilon_II": eps, "check": check}


    def update_outputs(self):
        """Override: fill self._outputs = {name: value} from the current state.
        Called by solve() after set_x(sol.x). Default: nothing."""
        pass

    def output_names(self):
        """Derived output names currently exposed (populated after update_outputs)."""
        return list(self._outputs.keys())

    def __getattr__(self, name):
        # Called only when normal lookup fails -> never shadows real attributes.
        # Redirects a bare name to a derived output stored in self._outputs.
        if name.startswith("_"):
            raise AttributeError(name)
        outs = self.__dict__.get("_outputs", {})
        if name in outs:
            return outs[name]
        raise AttributeError(
            f"{type(self).__name__} {getattr(self, 'local_name', '?')!r} has no "
            f"attribute or output {name!r}"
            + (f" (outputs: {sorted(outs)})" if outs else
               " (no outputs yet -- run solve() first)"))

