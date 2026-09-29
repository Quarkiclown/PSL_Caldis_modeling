import math

from caldis.core.component import Component
from caldis.core.connection import Connection

def _dtlm(dT1, dT2, name):
    """Log-mean temperature difference, sign-agnostic.
    dT1 == dT2 (equal, non-zero) is a normal case -> returns that value.
    Only dT1 == dT2 == 0 or a genuine cross (opposite signs) is ill-defined -> NaN."""
    a, b = abs(dT1), abs(dT2)
    if dT1 * dT2 < 0.0:
        return float("nan")
    if a < 1e-12 and b < 1e-12:
        return float("nan")
    if abs(a - b) < 1e-9:
        return max(a, b)
    return (a - b) / math.log(a / b)


import math

from caldis.core.component import Component


def _dtlm(dT1, dT2, name):
    """Log-mean temperature difference, sign-agnostic.
    dT1 == dT2 (equal, non-zero) is a normal case -> returns that value.
    Only dT1 == dT2 == 0 or a genuine cross (opposite signs) is ill-defined -> NaN."""
    a, b = abs(dT1), abs(dT2)
    if dT1 * dT2 < 0.0:
        return float("nan")
    if a < 1e-12 and b < 1e-12:
        return float("nan")
    if abs(a - b) < 1e-9:
        return max(a, b)
    return (a - b) / math.log(a / b)


class SizingHX(Component):
    """Generic single-zone counter-flow heat exchanger for sizing.

    ONE model for every heat exchanger (condenser, evaporator, recuperator,
    gas/gas, liquid/liquid), no 'role'. Single-zone LMTD; approximately valid even
    across a phase change, never fails on sign conventions.

    Channels 'hot'/'cold' are a documentation convention: the maths is sign-agnostic
    (a reversed hookup still computes), the names keep results readable and let
    start-sanity warn about a mis-wire.

    Ports (4): hot_in, hot_out (backend_hot); cold_in, cold_out (backend_cold).

    Parameters:
        h_hot, h_cold          : per-side CONVECTIVE (film) coefficients [W/m2/K].
        area_ratio_cold_to_hot : A_cold / A_hot (>=1 finned on the cold side, =1 plate).

    Overall coefficient (both films in series, wall neglected), referred to A_hot:
        U_overall_hot  = 1 / (1/h_hot + 1/(area_ratio*h_cold))
        A_hot          = UA / U_overall_hot
        A_cold         = area_ratio * A_hot ;  A_total = A_hot + A_cold
        U_overall_cold = U_overall_hot / area_ratio      (same UA, cold-side reference)

    Unknowns (2, solved by the system): UA, A_hot (read as hx.UA, hx.A_hot -> value).
    Derived outputs (after solve): U_overall_hot, U_overall_cold, A_cold, A_total,
        Qdot, DTLM.

    Residuals (7): 1-4 mass/dp per side; 5 coupling Qdot_hot + Qdot_cold = 0;
    6 UA = |Qdot_hot| / LMTD; 7 A_hot = UA / U_overall_hot. Pressures come from
    the system."""

    def __init__(self, name, backend_hot, backend_cold,
                 h_hot, h_cold, area_ratio_cold_to_hot):
        super().__init__(name)
        self.backend_hot = backend_hot
        self.backend_cold = backend_cold
        self.params.update(h_hot=h_hot, h_cold=h_cold,
                           area_ratio_cold_to_hot=area_ratio_cold_to_hot)
        self.hot_in   = self.add_port("hot_in",   backend=backend_hot)
        self.hot_out  = self.add_port("hot_out",  backend=backend_hot)
        self.cold_in  = self.add_port("cold_in",  backend=backend_cold)
        self.cold_out = self.add_port("cold_out", backend=backend_cold)
        self._UA    = self.add_variable("UA",    start=1000.0, lower=0.0, scale=1e3)
        self._A_hot = self.add_variable("A_hot", start=1.0,    lower=0.0, scale=1.0)

    # --- public value access for the unknowns ---
    @property
    def UA(self): return self._UA.value
    @property
    def A_hot(self): return self._A_hot.value

    def _U_overall_hot(self):
        """Overall heat-transfer coefficient referred to A_hot [W/m2/K].
        NOT a film coefficient: it combines both films in series."""
        p = self.params
        return 1.0 / (1.0 / p["h_hot"]
                      + 1.0 / (p["area_ratio_cold_to_hot"] * p["h_cold"]))

    def _compute(self):
        """Shared physics for residuals, outputs and sanity: the 4 port temperatures,
        the two heat rates and the LMTD. No side effect (stores nothing)."""
        bh, bc = self.backend_hot, self.backend_cold
        Th_in  = bh.T(self.hot_in.p,   self.hot_in.h)
        Th_out = bh.T(self.hot_out.p,  self.hot_out.h)
        Tc_in  = bc.T(self.cold_in.p,  self.cold_in.h)
        Tc_out = bc.T(self.cold_out.p, self.cold_out.h)
        Qdot_hot  = self.hot_in.mdot   * (self.hot_in.h  - self.hot_out.h)
        Qdot_cold = self.cold_out.mdot * (self.cold_out.h - self.cold_in.h)
        DTLM = _dtlm(Th_in - Tc_out, Th_out - Tc_in, self.path)
        return Th_in, Th_out, Tc_in, Tc_out, Qdot_hot, Qdot_cold, DTLM

    def residuals(self):
        U_overall_hot = self._U_overall_hot()
        _, _, _, _, Qdot_hot, Qdot_cold, DTLM = self._compute()
        return [
            self.hot_in.mdot + self.hot_out.mdot,       # 1 mass hot
            self.hot_out.p - self.hot_in.p,             # 2 dp hot
            self.cold_in.mdot + self.cold_out.mdot,     # 3 mass cold
            self.cold_out.p - self.cold_in.p,           # 4 dp cold
            Qdot_hot + Qdot_cold,                       # 5 coupling
            self._UA.value - abs(Qdot_hot) / DTLM,      # 6 UA def
            self._A_hot.value - self._UA.value / U_overall_hot,   # 7 A def
        ]

    def residual_names(self):
        return ["mass_hot", "dp_hot", "mass_cold", "dp_cold",
                "coupling_Q", "def_UA", "def_A"]

    def update_outputs(self):
        beta = self.params["area_ratio_cold_to_hot"]
        _, _, _, _, Qdot_hot, _, DTLM = self._compute()
        U_overall_hot = self._U_overall_hot()
        self._outputs = {
            "U_overall_hot": U_overall_hot,
            "U_overall_cold": U_overall_hot / beta,
            "A_cold": beta * self._A_hot.value,
            "A_total": (1.0 + beta) * self._A_hot.value,
            "Qdot": Qdot_hot,
            "DTLM": DTLM,
        }

    @property
    def U_overall_hot(self): return self._outputs["U_overall_hot"]
    @property
    def U_overall_cold(self): return self._outputs["U_overall_cold"]
    @property
    def A_cold(self): return self._outputs["A_cold"]
    @property
    def A_total(self): return self._outputs["A_total"]
    @property
    def Qdot(self): return self._outputs["Qdot"]
    @property
    def DTLM(self): return self._outputs["DTLM"]

    def _own_start_sanity(self):
        Th_in, Th_out, Tc_in, Tc_out, *_ = self._compute()
        dT1 = Th_in - Tc_out
        dT2 = Th_out - Tc_in
        return [
            (Th_in >= Tc_out and Th_out >= Tc_in,
             f"{self.path}: 'hot' channel is not hotter than 'cold' at start "
             f"(hot {Th_in-273.15:.0f}/{Th_out-273.15:.0f}°C vs "
             f"cold {Tc_in-273.15:.0f}/{Tc_out-273.15:.0f}°C) -- reversed hookup?"),
            (dT1 * dT2 > 0.0,
             f"{self.path}: temperature cross at start (dT1={dT1:.1f}, dT2={dT2:.1f} K)"),
            (abs(Th_in - Th_out) > 0.1,
             f"{self.path}: hot inlet/outlet nearly equal -> LMTD/flow degenerate"),
            (abs(Tc_out - Tc_in) > 0.1,
             f"{self.path}: cold inlet/outlet nearly equal -> cold MASS FLOW is "
             f"undetermined (anchor the cold-side inlet/outlet TEMPERATURES)"),
        ]

    def exergy_balance(self, T0, p0):
        Xh = self.hot_in.exergy_flow(T0, p0) - self.hot_out.exergy_flow(T0, p0)
        Xc = self.cold_out.exergy_flow(T0, p0) - self.cold_in.exergy_flow(T0, p0)
        return self._finish_balance(T0, p0, Xh, Xc)
        Xh = self.hot_in.exergy_flow(T0, p0) - self.hot_out.exergy_flow(T0, p0)
        Xc = self.cold_out.exergy_flow(T0, p0) - self.cold_in.exergy_flow(T0, p0)
        return self._finish_balance(T0, p0, Xh, Xc)