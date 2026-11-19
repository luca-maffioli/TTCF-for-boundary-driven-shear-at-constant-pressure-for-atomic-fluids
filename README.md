# Mother / child NEMD simulations of a fluid sheared between two walls

> **Reference.** This code was used to produce the results published in
>
> L. Maffioli, E. R. Smith, J. P. Ewen, P. J. Daivis, D. Dini and B. D. Todd,
> *Slip and stress from low shear rate nonequilibrium molecular dynamics: The
> transient-time correlation function technique*,
> J. Chem. Phys. **156**, 184111 (2022).
> [doi:10.1063/5.0088127](https://doi.org/10.1063/5.0088127)
>


```bibtex
@article{Maffioli2022,
  author  = {Maffioli, Luca and Smith, Edward R. and Ewen, James P. and Daivis, Peter J. and Dini, Daniele and Todd, B. D.},
  title   = {Slip and stress from low shear rate nonequilibrium molecular dynamics: The transient-time correlation function technique},
  journal = {The Journal of Chemical Physics},
  volume  = {156},
  number  = {18},
  pages   = {184111},
  year    = {2022},
  doi     = {10.1063/5.0088127}
}
```

This program runs non-equilibrium molecular dynamics (NEMD) simulations of a
three-dimensional fluid confined between two atomistic walls. It is built for
transient-time correlation function (TTCF) calculations, and alternates between two
kinds of trajectory:

- **Mother:** a long equilibrium trajectory. It samples starting states, and every
  `delay` steps it saves a snapshot of the whole system.
- **Child:** started from each snapshot. It applies four phase-space mappings to the
  snapshot and runs a short sheared trajectory from each mapped state, recording:
  - the shear and normal pressure on a set of planes (method of planes);
  - the work done by the walls;
  - the velocity, temperature and density profiles across the channel.

The mother and the child run in the same program. After every snapshot the child runs
to completion, and the mother then resumes from exactly the state it was in. The child
works on its own copies of all data structures, so it never touches the mother's
trajectory.

---

## Contents

1. [Building and running](#building-and-running)
2. [System](#system)
3. [Equations of motion](#equations-of-motion)
4. [Integrators](#integrators)
5. [Mother/child workflow and mappings](#motherchild-workflow-and-mappings)
6. [Measured quantities](#measured-quantities)
7. [Parameters](#parameters)
8. [Data structures](#data-structures)
9. [Function reference](#function-reference)
10. [Input and output files](#input-and-output-files)
11. [Notes and limitations](#notes-and-limitations)

---

## Building and running

```bash
g++ -O3 -march=native motherchild.cpp -o motherchild
./motherchild integrator c_pressure nchild PE index gammadot
```

| Argument     | Meaning |
|--------------|---------|
| `integrator` | 1 = velocity Verlet, any other value = 4th-order Runge–Kutta |
| `c_pressure` | 1 = constant pressure (the upper wall moves in y), any other value = constant volume |
| `nchild`     | Number of snapshots taken by the mother, i.e. number of child runs |
| `PE`         | Imposed pressure (used only when `c_pressure = 1`) |
| `index`      | Line of `seed.txt` that holds the random seed of the mother |
| `gammadot`   | Shear rate of the child runs |

- **Integrator and pressure control:** the mother and the child use the same
  integrator and the same `c_pressure` setting.
- **Shear:** the mother is always at equilibrium; only the children are sheared.
- **Run length:** the run takes `warmup + (nchild − 1)·delay + 1` mother steps, plus
  4·`simtime` child steps per snapshot.

---

## System

```
        y
        ▲   ● ● ● ● ● ● ●     upper wall ("h"): fcc slab, sites move at +wallvel  →
        │   ● ● ● ● ● ● ●     (optional barostat: the sites move rigidly in y)
        ├──────────────────
        │    ○  ○  ○  ○       fluid: fcc start, periodic in x and z
        │  ○  ○  ○  ○  ○
        ├──────────────────
        │   ● ● ● ● ● ● ●     lower wall ("l"): fcc slab, sites move at −wallvel  ←
        │   ● ● ● ● ● ● ●
        └──────────────────▶ x
```

**Interactions.** All pairs interact through a Lennard-Jones potential cut at 2.2 σ:

| Pair         | σ          | ε            | cutoff  |
|--------------|------------|--------------|---------|
| wall–wall    | `sigmaWW`  | `epsWW` = 1  | `rWW`   |
| wall–fluid   | `sigmaWF`  | `epsWF` = 0.4| `rWF`   |
| fluid–fluid  | `sigmaFF`  | `epsFF` = 1  | `rFF`   |

The force is not shifted at the cutoff. The particle type (`ptype`: 0 = wall,
1 = fluid) selects the row and column of the tables `seff`, `eeff` and `reff`.

**Walls.**

- **Structure:** each wall is an fcc slab of `wlayer` layers.
- **Tethers:** every wall particle is tied by a harmonic spring of constant `K` to its
  own lattice site.
- **Shear:** to drive the flow, the lower-wall sites move at −`wallvel` and the
  upper-wall sites at +`wallvel` along x.
- **Interactions:** the two walls are too far apart to interact with each other.

**Periodicity.** The box is periodic in x and z. There are no periodic boundaries in y.

**Particle indices:**

| Range                      | Particles  |
|----------------------------|------------|
| `[0, nwall)`               | lower wall |
| `[nwall, nwall+nfluid)`    | fluid      |
| `[nwall+nfluid, npart)`    | upper wall |

**Size with the default parameters:**

| Quantity           | Value |
|--------------------|-------|
| Box `Lx × Ly × Lz` | 30.78 × 13.42 × 30.78 |
| Fluid particles    | 6632 (density 0.700) |
| Wall particles     | 2 × 1296 |
| Total `npart`      | 9224 |
| Initial fluid layers | y ∈ [−4.70, 4.70] |
| Wall lattice sites | ±[5.70, 6.56] |

**Initial state** (`setLattice`, `randomMomenta`).

- **Fluid:** placed on an fcc lattice of 18 × 6 × 18 cells (7776 sites). The 1144
  sites in the middle of the list are left empty, which leaves `nfluid` particles.
- **Walls:** each fcc slab is placed at (σ<sub>FF</sub> + σ<sub>WW</sub>)/2 from the
  outermost fluid layer.
- **Momenta:** random, with zero total momentum, rescaled to temperature `T` with
  3·npart − 3 degrees of freedom.

---

## Equations of motion

All particles have unit mass. The lower wall (L) uses α<sub>0</sub> and the upper wall
(H) uses α<sub>1</sub>.

```
fluid:   dq/dt = p          dp/dt = F
walls:   dq/dt = p          dp/dt = F − α p                  (velocity Verlet)

thermostat (Nosé–Hoover, one per wall):
         dα/dt = ( Σ_wall p'² − 3·nwall·T ) / Q
```

- **Peculiar momentum.** p' is the wall momentum with the wall's mean momentum
  removed along x.
- **Thermostat mass.** Q = `QVV` with velocity Verlet and Q = `QRK` with RK4.
- **RK4.** The thermostat variables are integrated, but the friction term is not
  applied to the wall momenta (see [Integrators](#integrators)).

**Forces on the wall particles.** These include the springs
−K (q<sub>i</sub> − site<sub>i</sub>), with the x and z separations wrapped
periodically.

**Wall sites.** The sites move rigidly in x at ∓`wallvel`:

- **Mother:** `wallvel = 0` (equilibrium).
- **Child:** `wallvel = gammadot · DL`, where DL is half the distance between the
  innermost site layers of the two walls. This gives a linear streaming profile with
  shear rate `gammadot`.

**Barostat** (`c_pressure = 1`). The upper-wall sites move together in y, as a piston
of mass `M` = 1:

```
d yh/dt  = pyh
d pyh/dt = fyh = < K (y_i − yh_i) >_upper wall  −  FP,      FP = PE·Lx·Lz / nwall
```

With `c_pressure = 0`, fyh = 0 and the walls are fixed in y.

---

## Integrators

### Velocity Verlet (`integrator = 1`)

1. Compute the forces at t.
2. Move the positions to t + dt. The walls include the friction term −αp.
3. Move the wall sites by ±`wallvel`·dt. With the barostat, also move the upper sites
   in y.
4. Advance the momenta by a half step with the forces at t, and α by a half step.
5. Compute the forces at t + dt, then advance α by the second half step.
6. Advance the momenta by the second half step. For the walls, the friction is treated
   implicitly:

   p = (p + dt/2·F) / (1 + dt/2·α)

### Fourth-order Runge–Kutta (any other value)

**What is integrated.** Classical RK4 on positions, momenta, the thermostat variables
and the barostat variables (`yh`, `pyh`).

**Stages.**

- The wall sites move by ±`wallvel`·dt/2 before stage 2 and again before stage 4, so
  they sit at the mid-step position in stages 2–3 and at the end-point in stage 4.
- The thermostat variables are integrated, but their friction term does not enter the
  equations of the wall momenta.

---

## Mother/child workflow and mappings

```
mother(par)
 ├─ setLattice, randomMomenta          initial state
 ├─ warmup steps                       equilibration
 └─ until nchild snapshots are taken:
      one step
      every delay steps:
         printphase(k)  → child<k>.bin
         child(par, k)  → 4 mapped sheared trajectories
         (the mother resumes from its own, unchanged state)

child(par, k)
 └─ for each mapping m = 1..4:
      importinstate(k)                 reload child<k>.bin, set wallvel
      mappings(m)                      apply the mapping
      simtime steps:
         set the plane positions
         MOP_pressure                  pressures and profiles at this step
         computedissipation            wall work at this step
         one integration step
      printonfile                      append the results of mapping m
```

**The four mappings** (`mappings`), applied to the reloaded snapshot:

| m | Mapping                       | Transformation |
|---|-------------------------------|----------------|
| 1 | identity                      | — |
| 2 | time reversal                 | p → −p (particles and barostat) |
| 3 | reflection in x               | x → −x, p<sub>x</sub> → −p<sub>x</sub> (particles and wall sites) |
| 4 | reflection in x + time reversal | x → −x, p<sub>y</sub> → −p<sub>y</sub>, p<sub>z</sub> → −p<sub>z</sub> (particles), p<sub>y</sub> of the barostat |

Each mapping takes an equilibrium state to another equally probable one, so all four
trajectories are valid TTCF starting points. Averaging over the four cancels the
contributions that are odd under the mappings.

The thermostat variables α are not changed by the mappings. That is still a valid
mapping, because the equilibrium distribution does not depend on the sign of p alone.
With α flipped in mappings 2 and 4, those trajectories would be the exact time-reversed
images of mapping 1.

---

## Measured quantities

All of these are recorded by the child, at every step of every mapped trajectory.

### Pressure on planes (method of planes)

**Plane positions.** There are `nplane` = 7 planes, and they move with the walls:

| Plane | Position |
|-------|----------|
| 0 | between the two layers of the lower-wall lattice sites |
| 1 | on the fluid side of the lower wall, σ<sub>WW</sub>/2 above its top site layer |
| 2 … 4 | equally spaced between planes 1 and 5 |
| 5 | on the fluid side of the upper wall, σ<sub>WW</sub>/2 below its bottom site layer |
| 6 | between the two layers of the upper-wall lattice sites |

**Pressure.** For a plane at y<sub>p</sub> with area A = Lx·Lz:

```
P_αy(y_p) =  (1/A)    Σ_{pairs i,j across the plane} F_ij,α · sgn(y_i − y_j)
          +  (1/(A dt)) Σ_{particles crossing during the last step} p'_α p'_y
```

for α = x, y.

- **Configurational term:** sums over the pairs that interact across the plane.
- **Kinetic term:** uses the peculiar momenta p' (streaming velocity of the particle's
  bin removed). They are interpolated linearly to the time at which the particle
  crosses the plane.

### Wall work (`computedissipation`)

The power delivered by the moving lattice sites through the springs is

```
oml = Σ_lower  wallvel · K · (x_i − x_site,i)
omh = − Σ_upper wallvel · K · (x_i − x_site,i)
```

The power of the barostat piston is also computed, but it is not written to file.

### Profiles across the channel

The bins have width ≈ 0.1 σ and span y ∈ [−Ly, Ly] (`nbin` = 268). For each bin the
child records:

- the streaming velocity ⟨p<sub>x</sub>⟩;
- the temperature, ⟨p'²⟩/3;
- the number density.

All particles are included, walls and fluid.

---

## Parameters

### Interactions, thermostat, barostat

| Constant              | Default | Meaning |
|-----------------------|---------|---------|
| `sigmaW`, `sigmaF`    | 1       | Wall and fluid diameters |
| `epsW`, `epsF`        | 1       | Wall and fluid energies |
| `epsWF`               | 0.4     | Wall–fluid energy |
| `rWW`, `rWF`, `rFF`   | 2.2 σ   | Cutoffs |
| `T`                   | 1       | Wall temperature |
| `K`                   | 150     | Spring constant of the wall tethers |
| `QVV`, `QRK`          | 10, 0.4 | Nosé–Hoover masses for velocity Verlet and RK4 |
| `M`                   | 1       | Barostat piston mass |
| `dfree`               | 3       | Degrees of freedom removed when scaling the initial momenta |
| `gammadot` (constant) | 0       | Shear rate of the mother. The child's shear rate is the command-line argument |

### Geometry

| Constant               | Default | Meaning |
|------------------------|---------|---------|
| `tLx`, `tLy`, `tLz`    | 30, 10, 30 | Target box size. `tLy` is the width of the fluid slab |
| `rho_f`, `rho_w`       | 0.7, 0.8 | Fluid and wall densities |
| `wlayer`               | 2       | Wall layers. Must be even |
| `wcell`                | 1.710   | fcc cell of the walls |
| `nwx`, `nwy`, `nwz`    | 18, 1, 18 | Wall fcc cells in x, y, z |
| `Lx`, `Lz`             | 30.78   | Box lengths: whole numbers of wall cells |
| `Ly`                   | 13.42   | Fluid slab + both walls |
| `nwall` (= `nlatt`)    | 1296    | Particles (and lattice sites) per wall |
| `nfluid`               | 6632    | Fluid particles (multiple of 4) |
| `nxside`, `nyside`, `nzside` | 18, 6, 18 | Fluid fcc cells |
| `cellf`                | 1.710   | Fluid fcc cell |

### Cell list, profiles, planes

| Constant          | Default | Meaning |
|-------------------|---------|---------|
| `dcell`           | 2.2     | Minimum cell size (= largest cutoff) |
| `clx`, `cly`, `clz` | 13, 12, 13 | Cells. The grid spans y ∈ [−Ly, Ly] |
| `dcx`, `dcy`, `dcz` | 2.368, 2.237, 2.368 | Actual cell sizes |
| `dbin`, `bin`     | 0.1, 0.1001 | Target and actual profile bin width |
| `nbin`            | 268     | Profile bins |
| `vbin`            | 94.88   | Bin volume |
| `nplane`          | 7       | Method-of-planes planes |

### Time stepping

| Constant  | Default | Meaning |
|-----------|---------|---------|
| `dt`      | 0.002   | Time step |
| `warmup`  | 400000  | Mother: equilibration steps |
| `delay`   | 1000    | Mother: steps between snapshots |
| `simtime` | 5000    | Child: steps per mapped trajectory |

---

## Data structures

### `Phase`: state of the system

| Fields | Meaning |
|--------|---------|
| `x, y, z, px, py, pz` | Particle positions and momenta |
| `xl, yl, zl` | Lattice sites of the lower wall |
| `xh, yh, zh, pxh, pyh, pzh` | Lattice sites of the upper wall and their momenta. Only `yh`, `pyh` evolve (barostat) |
| `alpha[2]` | Nosé–Hoover variables: [0] lower, [1] upper wall |
| `gammadot`, `wallvel` | Shear rate and velocity of the lattice sites |
| `PE`, `FP` | Imposed pressure and force per lattice site |
| `x0…pz4`, `alpha0…4`, `yh0…4`, `pyh0…4` | RK4: state at the start of the step and stage increments |
| `plist` | Cell list: next particle in the same cell |
| `pxth, pyth, yp, pxthp, pythp` | Method of planes: peculiar momenta now and at the previous step, y at the previous step |
| `integrator, c_pressure, mapping, index, nchild` | Run settings |
| `ypl[nplane]` | Plane positions |
| `ptype`, `seff`, `eeff`, `reff` | Particle types and the interaction table of each pair of types |

### Other structures

| Structure  | Fields | Used by |
|------------|--------|---------|
| `Nlist`    | `atom[totlist]`: first particle of each cell | both |
| `Interact` | `fx, fy, fz`; `dalpha[2]` (dα/dt); `fyh` (barostat force) | both |
| `Thermo`   | `oml, omh` (wall work), `om2l, om2h` (always 0), `barostat`, `omtime`. One value per time step | child |
| `Pressure` | `mop_pxy, mop_pyy` (configurational), `mop_kpxy, mop_kpyy` (kinetic), per step and plane; `vstream, temp, n`, per step and bin; work arrays of the current step; `mop_time` | child |
| `RunParameters` | Command-line arguments | `mother`, `child` |

All structures are allocated with `new X()`, so every field starts at zero.

---

## Function reference

### Workflow

#### `int main(int argc, char **argv)`
Reads the six command-line arguments, sets `c_pressure` to 0 unless it is 1, and calls
`mother`.

#### `void mother(const RunParameters &par)`
Runs the mother:

1. creates its own `Phase`, `Interact` and `Nlist`;
2. builds the initial state;
3. runs `warmup` steps;
4. then integrates, and whenever the step counter is a multiple of `delay`:
   - saves a snapshot with `printphase`;
   - runs `child` on that snapshot;
5. stops after `nchild` snapshots.

#### `void child(const RunParameters &par, int k)`
Runs the four mapped trajectories from snapshot k:

1. creates its own structures (separate from the mother's);
2. opens `outputXX<c>_child<k>.bin` and `velprofile_child<k>.bin`;
3. for each mapping: reloads the snapshot, applies the mapping, resets the
   accumulators, and runs `simtime` steps. At every step it places the planes and calls
   `MOP_pressure` and `computedissipation`, then integrates;
4. writes the results of each mapping with `printonfile`.

### Mother routines

#### `void rand0(Phase *phase)`
Seeds `rand()` with the seed on line `index` of `seed.txt` and prints the seed. If the
file is missing, the seed comes from the clock.

#### `int randomInteger(int lower, int upper)`
Returns a uniform random integer in [lower, upper].

#### `void setLattice(Phase *phase)`
Builds the initial positions (fcc fluid between two fcc walls, see [System](#system)):

- centres the box in x and z, and the fluid in y;
- places the walls below and above the fluid;
- puts the lattice sites at the initial wall positions;
- sets `wallvel = 0`;
- calls `setparticletypes`.

#### `void randomMomenta(Phase *phase, int min, int max)`
Sets the initial momenta:

1. draws random integer momenta in [min, max];
2. removes the total momentum;
3. rescales to temperature `T` with 3·npart − `dfree` degrees of freedom;
4. sets the barostat momenta to zero.

#### `void printphase(Phase *phase, int i)`
Writes a snapshot to `child<i>.bin`. The format is given in
[Input and output files](#input-and-output-files).

### Child routines

#### `void importinstate(Phase *phase, int k)`
Reads `child<k>.bin`, sets `wallvel = gammadot·DL`, and calls `setparticletypes`.

#### `void mappings(Phase *phase)`
Applies mapping `phase->mapping` (see the table in
[Mother/child workflow and mappings](#motherchild-workflow-and-mappings)).

#### `void MOP_pressure(Phase *phase, Pressure *pressure, Nlist *nlist)`
Computes all child measurements for the current step:

- **Configurational term.** Rebuilds the cell list. For each plane, it visits only the
  cell rows next to the plane, and adds 4·F·sgn(r<sub>y</sub>) for every pair that
  interacts across the plane. The factor 4 is divided out in `printonfile`.
- **Profiles.** Computes the streaming velocity, the temperature and the number of
  particles in each bin, and the peculiar momenta of every particle.
- **Kinetic term.** For every particle that crossed a plane since the previous step,
  adds p'<sub>x</sub>p'<sub>y</sub> and p'<sub>y</sub>p'<sub>y</sub>, interpolated to
  the crossing time.
- **Storage.** Stores the profiles of this step and saves y and p' for the next
  crossing check.

#### `void computedissipation(Phase *phase, Thermo *thermo)`
Stores the wall work `oml`, `omh` and the barostat power of the current step.

#### `void printonfile(Phase *phase, Pressure *pressure, Thermo *thermo, FILE *writemop, FILE *velprof)`
Appends the results of one mapped trajectory to the two output files. The format is
given in [Input and output files](#input-and-output-files).

### Shared routines

#### `void pbcx(double *x)`, `void pbcz(double *z)`
Wrap an x (or z) coordinate or separation into [−L/2, L/2).

#### `void setparticletypes(Phase *phase)`
Sets `ptype` (0 = wall, 1 = fluid) and the tables `seff`, `eeff`, `reff`.

#### `void neighbourlist(Nlist *nlist, Phase *phase)`
Builds the cell list. Each cell holds a linked list (`atom` = head, `plist` = next).
Particles are appended at the tail, so each list is in increasing index.

#### `void wcaforce(Phase *phase, Interact *interact, Nlist *nlist)`
Pair forces through the cell list. It visits the 27 neighbouring cells, periodic in x
and z and clipped in y, and takes each pair once (`pn < pc`).

#### `void springforce(Phase *phase, Interact *interact)`
Adds the spring forces −K (q − site) to the wall particles. It also accumulates the
mean y reaction on the upper-wall sites in `fyh`.

#### `void thermostatforce(Phase *phase, Interact *interact)`
Computes dα/dt for each wall, as in [Equations of motion](#equations-of-motion).

#### `void force(Phase *phase, Interact *interact, Nlist *nlist)`
Computes all forces at the current state:

1. zeroes the forces, dα/dt and `fyh`;
2. calls `wcaforce`, `springforce` and `thermostatforce`;
3. sets fyh = c_pressure·(fyh − FP).

#### `void VelVerlet(Phase *phase, Interact *interact, Nlist *nlist)`
#### `void RK4(Phase *phase, Interact *interact, Nlist *nlist)`
One time step with each integrator (see [Integrators](#integrators)).

---

## Input and output files

All binary files contain native-endian `double`s.

### Input: `seed.txt`
One integer seed per line. The mother reads line `index`.

### Snapshot: `child<k>.bin` (written by the mother, read by the child)

| Block | Content | Doubles |
|-------|---------|---------|
| particles    | x y z px py pz for each of the `npart` particles | 6·npart |
| lower sites  | xl yl zl for each of the `nlatt` sites | 3·nlatt |
| upper sites  | xh yh zh pxh pyh pzh for each of the `nlatt` sites | 6·nlatt |
| thermostats  | α<sub>0</sub> α<sub>1</sub> | 2 |

With the default parameters a snapshot is 67 010 doubles (≈ 0.5 MB).

### Pressure and work: `outputVV<c>_child<k>.bin` or `outputRK<c>_child<k>.bin`

`<c>` is `c_pressure`. The file holds the 4 mappings one after another, each with
`simtime` records of `2·nplane + 2` doubles:

```
P_xy(plane 0) P_yy(plane 0) … P_xy(plane 6) P_yy(plane 6)  oml  omh
```

Shape: `[4][simtime][2·nplane + 2]` = [4][5000][16] with the default parameters.

### Profiles: `velprofile_child<k>.bin`

The file holds the 4 mappings one after another, each with `simtime × nbin` records:

```
streaming velocity   temperature   number density
```

Shape: `[4][simtime][nbin][3]` = [4][5000][268][3] with the default parameters
(≈ 129 MB per child).

### Reading the outputs (Python)

```python
import numpy as np
nplane, nbin, simtime = 7, 268, 5000

mop = np.fromfile("outputVV0_child0.bin").reshape(4, simtime, 2*nplane + 2)
Pxy = mop[:, :, 0:2*nplane:2]      # [mapping, step, plane]
Pyy = mop[:, :, 1:2*nplane:2]
oml, omh = mop[:, :, -2], mop[:, :, -1]

prof = np.fromfile("velprofile_child0.bin").reshape(4, simtime, nbin, 3)
vx, temp, dens = prof[..., 0], prof[..., 1], prof[..., 2]
y = -13.42 + (np.arange(nbin) + 0.5) * 2*13.42/nbin   # bin centres (default Ly)
```

---

## Notes and limitations

- **`wlayer` must be even.** `nwy = ceil(wlayer / 2)` uses integer division, and
  `nwall = 2·nwx·nwz·wlayer` agrees with the number of wall particles that
  `setLattice` places (4·nwx·nwz·nwy) only when `wlayer` is even. To allow odd values,
  use `nwy = (wlayer + 1) / 2` and `nwall = 4·nwx·nwz·nwy`.
- **Fluid density.** The number of fluid particles is rounded up to a multiple of 4,
  so the actual density (`rho_t`) is slightly above `rho_f`.
- **No bounds checks on cells and bins.** The cell list and the profiles cover
  y ∈ [−Ly, Ly], about twice the size of the system. A particle outside that range
  would write outside the arrays.
- **Memory.** The child stores the profiles of every step in memory: 3 × nbin ×
  simtime doubles, about 32 MB with the default parameters.
- **Unused fields.** `om2l`, `om2h` are always 0, and the barostat power is computed
  but not written.
- **Reproducibility.** The mother's random seed comes from `seed.txt`. If the file is
  missing, the clock is used and runs are not reproducible.
