# Caldis — bibliothèque n°1 (petits systèmes : PAC, lignes fluides)

Squelette minimal, **numpy + scipy uniquement**. Ouvert *et* fermé : rien n'y suppose une
boucle. On complète couche par couche.

## Lancer la démo
```bash
cd /mnt/user-data/outputs
python demo_open_system.py
```
Système ouvert `source -> détendeur -> puits` : 12 inconnues, 12 équations, résolu.

## Les couches
```
caldis/
  fluids/       propriétés fluides — état construit depuis (p, h)
    backend.py            interface abstraite (+ dérivées, à venir)
    constant_backend.py   PLACEHOLDER sans dépendance (à remplacer par CoolProp)
    state.py              FluidState(p, h) -> T, rho, x
  core/         abstractions structurantes, sans physique
    variable.py           inconnue : bornes, échelle, départ
    port.py               FluidPort (p, mdot, h) ; SignalPort (stub)
    component.py          base : déclare ports/vars + residuals()
    connection.py         connect() -> équations de nœud
    system.py             assemble x et R(x)
    dof.py                vérificateur de degrés de liberté (diagnostic clair)
  components/   bibliothèque métier
    boundary.py           Source / Sink (=> systèmes ouverts)
    valve.py              détendeur isenthalpique (1er composant réel)
  solvers/      stratégie interchangeable, ignore la physique
    steady.py             scipy.root (à faire évoluer : CasADi, bornage, homotopie)
  experiments/  scénarios = System + solveur + conditions
    scenario.py           assemble + check DOF + résout
```

## Conventions
- `FluidPort.mdot > 0` = masse **entrant** dans le composant.
- États toujours en **(p, h)**, jamais (p, T).

## Feuille de route (incréments, dans l'ordre suggéré)
1. **CoolProp** : `CoolPropBackend` (AbstractState + dérivées) derrière `FluidBackend`.
2. **Compresseur** + **échangeur** : premiers composants diphasiques (évaporateur = cas
   qui casse le plus — jeu d'états (p, h), résidus, discrétisation).
3. **Boucle fermée PAC** : premier système fermé, références internes au lieu de Source/Sink.
4. **CasADi** : résidus symboliques + jacobienne exacte (callbacks propriétés).
5. **Bornage + homotopie** : solveurs robustes vers la branche physique.
6. **Tearing physique** : boucle externe sur (p_bas, p_haut, surchauffe, sous-refroidissement).
7. **Dynamique** : résidu DAE F(x, ẋ, t)=0 + IDAS ; le permanent devient le cas ẋ=0.
```
```
