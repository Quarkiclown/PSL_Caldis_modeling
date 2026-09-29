import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

from caldis.core import System
from caldis.core import check_dof
from caldis.components import FluidBoundary
from caldis.components import Tank
from caldis.fluids import CoolPropBackend
from caldis.solvers import simulate

FLUID = "Air"
be = CoolPropBackend(FLUID)
p0, T0 = 1.0e5, 293.15
h0 = PropsSI("H", "P", p0, "T", T0, FLUID)
V = 0.05
C_wall, UA_wall, Q_heater = 800.0, 30.0, 150.0

sys_ = System("tank_test")
feed = FluidBoundary("feed", mdot=0.0, h=h0)       # tank ferme (debit nul)
tank = Tank("tank", backend=be, V=V, p_init=p0, h_init=h0,
            C_wall=C_wall, T_wall_init=T0, UA_wall=UA_wall, Q_heater=Q_heater)
cap = FluidBoundary("cap")
sys_.add(feed, tank, cap)
sys_.connect(feed.out, tank.inl)
sys_.connect(tank.out, cap.inl)
sys_.assemble()

# 1) DOF
print(check_dof(sys_))

# 2) noms hierarchiques
names = [v.name for v in sys_.unknowns]
for exp in ("tank.fluid.M", "tank.fluid.U", "tank.wall.T", "tank.fluid.inl.p"):
    assert exp in names, f"nom hierarchique manquant : {exp}"
print("noms hierarchiques OK :", [n for n in names if n.startswith("tank.")][:6], "...")

# 3) garde-fou : sortie inconnue -> erreur claire
tank.wall.outputs = ["Temperature"]             # faute volontaire (c'est "T")
try:
    sys_.check()
    print("ERREUR : garde-fou inactif")
except ValueError as e:
    print("garde-fou OK ->", e)
tank.wall.outputs = ["T"]

# 4) simulation : tank ferme chauffe -> pas d'etat permanent -> init=None
t, hist = simulate(sys_, t_end=400.0, dt=2.0, update_inputs=None, init=None)
print("\ncles enregistrees :", sorted(hist.keys()))

Tw = hist["tank.wall.T"] - 273.15
Tf = hist["tank.fluid.T"] - 273.15
p = hist["tank.fluid.p"] / 1e5
Q = hist["tank.Q"]
U = hist["tank.fluid.U"]

# conservation : dU_fluide + C_wall*dT_wall = Q_heater * t
dE = (U - U[0]) + C_wall * (hist["tank.wall.T"] - hist["tank.wall.T"][0])
print(f"conservation : ecart max |dE - Q_heater*t| = "
      f"{np.max(np.abs(dE - Q_heater * t)):.3f} J")

fig, ax = plt.subplots(1, 3, figsize=(14, 4))
ax[0].plot(t, Tw, 'r-', label="paroi"); ax[0].plot(t, Tf, 'b-', label="fluide")
ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("T (°C)"); ax[0].legend(); ax[0].grid(True)
ax[0].set_title("paroi mene, fluide suit")
ax[1].plot(t, p, 'g-'); ax[1].set_xlabel("t (s)"); ax[1].set_ylabel("p fluide (bar)")
ax[1].grid(True); ax[1].set_title("chauffage a V constant -> p monte")
ax[2].plot(t, Q, 'm-'); ax[2].set_xlabel("t (s)"); ax[2].set_ylabel("Q paroi->fluide (W)")
ax[2].grid(True); ax[2].set_title("flux de couplage")
fig.suptitle(f"Tank composite {FLUID} — Q_heater={Q_heater} W")
fig.tight_layout(); plt.savefig("test_tank.png", dpi=110); plt.show()