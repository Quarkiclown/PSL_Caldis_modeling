from caldis.core.component import Component
from caldis.components.vessels.control_volume import ControlVolume
from caldis.components.vessels.thermal_mass import ThermalMass


class Tank(Component):
    """Composite: fluid volume (ControlVolume) + wall (ThermalMass), coupled by
    Q = UA_wall*(T_wall - T_fluid). A heater Q_heater is applied to the wall.
    The tank's fluid ports ARE the inner volume's ports (aliases). The coupling is
    injected into the children via their Qdot (computed in _own_residuals, re-
    evaluated each iteration -> captured by the Jacobian)."""

    def __init__(self, name, backend, V, p_init, h_init,
                 C_wall, T_wall_init, UA_wall, Q_heater=0.0):
        super().__init__(name)
        self.backend = backend
        self.params.update(UA_wall=UA_wall, Q_heater=Q_heater)
        self.fluid = self.add_child(
            ControlVolume("fluid", backend=backend, V=V,
                          p_init=p_init, h_init=h_init, Qdot=0.0))
        self.wall = self.add_child(
            ThermalMass("wall", C=C_wall, UA=0.0,
                        T_init=T_wall_init, T_env=T_wall_init, Qdot_ext=0.0))
        self.inl = self.fluid.inl        # alias: tank exposes the volume's ports
        self.out = self.fluid.out

    def _own_residuals(self):
        UA = self.params["UA_wall"]
        Qh = self.params["Q_heater"]
        T_f = self.backend.T(self.fluid.p, self.fluid.h)
        T_w = self.wall.T
        Q = UA * (T_w - T_f)                          # wall -> fluid
        self._Q_coupling = Q                          # kept for update_outputs
        self.fluid.params["Qdot"] = Q                 # fluid receives +Q
        self.wall.params["Qdot_ext"] = Qh - Q         # wall: heater - Q released
        return []

    def update_outputs(self):
        self._outputs = {
            "Q": getattr(self, "_Q_coupling", 0.0),
            "T_fluid": self.backend.T(self.fluid.p, self.fluid.h),
            "T_wall": self.wall.T,
        }

    @property
    def Q(self): return self._outputs["Q"]
    @property
    def T_fluid(self): return self._outputs["T_fluid"]
    @property
    def T_wall(self): return self._outputs["T_wall"]