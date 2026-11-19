/*
 * ============================================================================
 *  Mother / child NEMD workflow for a fluid confined between two walls
 * ============================================================================
 *
 *  Purpose
 *    Generate equilibrium starting states with a long "mother" trajectory and,
 *    from each of them, run short sheared "child" trajectories. For every
 *    starting state, the child runs the 4 phase-space mappings used in
 *    transient-time correlation function (TTCF) calculations:
 *      1. identity
 *      2. time reversal                       (p -> -p)
 *      3. reflection in x                     (x -> -x, px -> -px)
 *      4. time reversal + reflection in x     (x -> -x, py -> -py, pz -> -pz)
 *    During each child trajectory it records the shear and normal pressure
 *    on nplane planes (method of planes), the work done by the walls, and the
 *    velocity / temperature / density profiles across the channel.
 *
 *  System (3D)
 *    - Fluid of WCA-like (LJ cut at 2.2 sigma) particles, periodic in x and z,
 *      confined in y between two walls.
 *    - Each wall: an fcc slab of particles tethered by harmonic springs
 *      (constant K) to lattice sites. To shear the fluid, the sites of the
 *      upper ("high", h) wall move at +wallvel and those of the lower
 *      ("low", l) wall at -wallvel along x.
 *    - Each wall has its own Nose-Hoover thermostat (variables alpha[0] for
 *      the lower wall, alpha[1] for the upper wall). The fluid is not
 *      thermostatted directly.
 *    - Optional constant pressure (c_pressure = 1): the lattice sites of the
 *      upper wall move rigidly in y under the imposed pressure PE (a
 *      piston of mass M).
 *
 *  Particle indexing
 *    [0, nwall)                   lower wall
 *    [nwall, nwall+nfluid)        fluid
 *    [nwall+nfluid, npart)        upper wall
 *
 *  Workflow
 *    mother(): build the initial state, run `warmup` equilibrium steps, then
 *      every `delay` steps save a snapshot (child<k>.bin) and run child() on
 *      it. The child uses its own copies of all data structures, so after it
 *      returns the mother resumes exactly where it stopped. Repeat until
 *      nchild snapshots have been processed.
 *    child(): load snapshot k and, for each of the 4 mappings, apply the
 *      mapping and run `simtime` sheared steps. The results of all 4
 *      mappings go to
 *        outputVV<c>_child<k>.bin / outputRK<c>_child<k>.bin   (pressures, work)
 *        velprofile_child<k>.bin                               (profiles)
 *      where <c> is c_pressure.
 *
 *  Usage
 *    ./motherchild integrator c_pressure nchild PE index gammadot
 *
 *      integrator : 1 = velocity Verlet, any other value = RK4
 *      c_pressure : 1 = constant pressure (moving upper wall), otherwise 0
 *      nchild     : number of snapshots / child runs
 *      PE         : imposed pressure (used when c_pressure = 1)
 *      index      : line of seed.txt holding the random seed of the mother
 *      gammadot   : shear rate of the child runs
 *
 *  Input file
 *    seed.txt : one seed per line; the mother uses line `index`. If the file
 *               is missing, the seed is taken from the clock.
 *
 *  Compile:  g++ -O3 -march=native motherchild.cpp -o motherchild
 *
 *  See DOCUMENTATION.md for the output file formats.
 * ============================================================================
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;


/* ============================================================================
 *  Interaction parameters
 * ========================================================================== */

const double sigmaW = 1;  // wall particle diameter
const double sigmaF = 1;  // fluid particle diameter
const double epsW   = 1;  // wall energy
const double epsF   = 1;  // fluid energy

const double sigmaWW = 0.5 * (sigmaW + sigmaW);  // pair diameters (Lorentz rule)
const double sigmaWF = 0.5 * (sigmaW + sigmaF);
const double sigmaFF = 0.5 * (sigmaF + sigmaF);

const double epsWW = pow(epsW * epsW, 0.5);      // pair energies (Berthelot rule)
const double epsWF = 0.4;                        //   wall-fluid set explicitly
const double epsFF = pow(epsF * epsF, 0.5);

const double rWW = 2.2 * sigmaWW;                // cutoffs
const double rWF = 2.2 * sigmaWF;
const double rFF = 2.2 * sigmaFF;


/* ============================================================================
 *  Thermostat, springs, barostat
 * ========================================================================== */

const double T     = 1;     // wall temperature
const double K     = 150;   // spring constant tethering the wall particles
const int    dfree = 3;     // dof removed when the initial momenta are scaled (total momentum)
const double QRK   = 0.4;   // Nose-Hoover thermostat mass, RK4 integrator
const double QVV   = 10;    // Nose-Hoover thermostat mass, velocity Verlet integrator
const double M     = 1;     // mass of the barostat piston (upper wall lattice)

const double gammadot = 0;  // shear rate of the mother trajectory (equilibrium)


/* ============================================================================
 *  Geometry
 * ========================================================================== */

const double tLx = 30;   // target box length in x
const double tLy = 10;   // target width of the fluid slab in y
const double tLz = 30;   // target box length in z

const double rho_f = 0.7;   // fluid density (reduced)
const double rho_w = 0.8;   // wall density (reduced)

const double rho_wr = rho_w / (sigmaW * sigmaW * sigmaW);  // densities in units of the box
const double rho_fr = rho_f / (sigmaF * sigmaF * sigmaF);

const int    wlayer = 2;                              // wall layers (must be even: 2 layers per fcc cell)
const double wcell  = pow(4 / rho_wr, 1. / 3);        // fcc cell of the walls (4 particles)

const int nwx = ceil(tLx / wcell);                    // wall fcc cells along x
const int nwy = ceil(wlayer / 2);                     // wall fcc cells along y
const int nwz = ceil(tLz / wcell);                    // wall fcc cells along z

const double Lx = nwx * wcell;                        // box length in x (whole number of wall cells)
const double Lz = nwz * wcell;                        // box length in z
const double Ly = tLy + 2 * nwy * wcell;              // fluid slab + both walls

const int nwall = 2 * nwx * nwz * wlayer;             // particles per wall
const int nlatt = nwall;                              // lattice sites per wall

const double tNf   = Lx * tLy * Lz * rho_fr;          // target number of fluid particles
const int    lNf   = 4 * ceil(tNf / 4);               // rounded up to a multiple of 4 (fcc)
const double rho_t = lNf / (Lx * tLy * Lz);           // resulting fluid density

const double tcellf = pow(4 / rho_t, 1. / 3);         // target fcc cell of the fluid

const int nxside = ceil(Lx / tcellf);                 // fluid fcc cells along x, y, z
const int nyside = ceil(tLy / tcellf);
const int nzside = ceil(Lz / tcellf);

const double cellf = Lx / nxside;                     // fluid fcc cell size

const int ttNf = nxside * nyside * nzside * 4;        // sites of the fluid lattice (>= nfluid)

const int nfluid = round(lNf);                        // fluid particles
const int npart  = 2 * nwall + nfluid;                // total particles

const double Lx2 = Lx / 2;
const double Ly2 = Ly / 2;
const double Lz2 = Lz / 2;


/* ============================================================================
 *  Cell list (forces and method of planes)
 *  The grid spans y in [-Ly, Ly] so that the walls fit even if they move.
 * ========================================================================== */

const double dcell = rFF;                  // minimum cell size = largest cutoff
const int    clx   = floor(Lx / dcell);    // cells along x, y, z
const int    cly   = floor(2 * Ly / dcell);
const int    clz   = floor(Lz / dcell);

const double dcx = Lx / clx;               // actual cell sizes
const double dcy = 2 * Ly / cly;
const double dcz = Lz / clz;

const int totlist = clx * cly * clz;       // number of cells


/* ============================================================================
 *  Profiles and method of planes
 * ========================================================================== */

const double dbin = sigmaFF / 10;          // target bin width of the profiles
const int    nbin = 2 * Ly / dbin;         // bins over y in [-Ly, Ly]
const double bin  = 2 * Ly / nbin;         // actual bin width
const double vbin = bin * Lx * Lz;         // bin volume

const int nplane = 7;                      // planes for the method-of-planes pressure


/* ============================================================================
 *  Time stepping
 * ========================================================================== */

const double dt = 0.002;     // time step

const int warmup  = 400000;  // mother: equilibration steps before the first snapshot
const int delay   = 1000;    // mother: steps between two snapshots
const int simtime = 5000;    // child: steps of each mapped trajectory


/* ============================================================================
 *  Data structures
 * ========================================================================== */

// State of the system and work arrays of the integrators
typedef struct
{
    // particle positions and momenta
    double x[npart];
    double y[npart];
    double z[npart];

    double px[npart];
    double py[npart];
    double pz[npart];

    // lattice sites of the upper (h) wall and their momenta (only pyh is used: barostat)
    double xh[nlatt];
    double yh[nlatt];
    double zh[nlatt];

    double pxh[nlatt];
    double pyh[nlatt];
    double pzh[nlatt];

    // lattice sites of the lower (l) wall
    double xl[nlatt];
    double yl[nlatt];
    double zl[nlatt];

    double alpha[2];      // Nose-Hoover variables: [0] lower wall, [1] upper wall

    double gammadot;      // shear rate
    double wallvel;       // velocity of the wall lattice sites

    double PE;            // imposed pressure
    double FP;            // corresponding force per lattice site, PE*Lx*Lz/nwall

    // RK4: state at the start of the step (0) and the four stage increments (1-4)
    double x0[npart], y0[npart], z0[npart], px0[npart], py0[npart], pz0[npart];
    double x1[npart], y1[npart], z1[npart], px1[npart], py1[npart], pz1[npart];
    double x2[npart], y2[npart], z2[npart], px2[npart], py2[npart], pz2[npart];
    double x3[npart], y3[npart], z3[npart], px3[npart], py3[npart], pz3[npart];
    double x4[npart], y4[npart], z4[npart], px4[npart], py4[npart], pz4[npart];

    double alpha0[2], alpha1[2], alpha2[2], alpha3[2], alpha4[2];

    double yh0[nlatt], yh1[nlatt], yh2[nlatt], yh3[nlatt], yh4[nlatt];
    double pyh0[nlatt], pyh1[nlatt], pyh2[nlatt], pyh3[nlatt], pyh4[nlatt];

    int plist[npart];     // cell list: next particle in the same cell (-1 = end)

    // method of planes, kinetic term
    double pxth[npart];   // peculiar momenta at the current step
    double pyth[npart];
    double yp[npart];     // y at the previous step
    double pxthp[npart];  // peculiar momenta at the previous step
    double pythp[npart];

    int integrator;       // 1 = velocity Verlet, otherwise RK4
    int c_pressure;       // 1 = constant pressure
    int mapping;          // child: current mapping (1-4)
    int index;            // mother: line of seed.txt
    int nchild;           // mother: number of snapshots

    double ypl[nplane];   // y of the method-of-planes planes

    int ptype[npart];     // 0 = wall, 1 = fluid

    double seff[2][2];    // sigma, epsilon and cutoff of each pair of types
    double eeff[2][2];
    double reff[2][2];
} Phase;

// Cell list heads
typedef struct
{
    int atom[totlist];    // first particle of each cell (-1 = empty)
} Nlist;

// Forces
typedef struct
{
    double fx[npart];
    double fy[npart];
    double fz[npart];

    double dalpha[2];     // time derivatives of the thermostat variables
    double fyh;           // net y force on each upper-wall lattice site (barostat)
} Interact;

// Child: work done by the walls (one value per time step)
typedef struct
{
    double oml[simtime];      // lower wall: power from the lattice-site motion
    double omh[simtime];      // upper wall
    double om2l[simtime];     // reserved (always 0)
    double om2h[simtime];
    double barostat[simtime]; // power of the barostat piston

    int omtime;               // current time index
} Thermo;

// Child: method-of-planes pressure and profiles (one record per time step)
typedef struct
{
    double vstream[nbin * simtime];   // streaming velocity profile
    double temp[nbin * simtime];      // temperature profile
    double n[nbin * simtime];         // number of particles per bin

    double mop_vxstream[nbin];        // profiles at the current step
    double mop_vystream[nbin];
    double mop_T[nbin];
    double mop_n[nbin];

    double mop_pxy[simtime * nplane]; // configurational pressure P_xy, P_yy on each plane
    double mop_pyy[simtime * nplane];
    double mop_kpxy[simtime * nplane];// kinetic pressure P_xy, P_yy on each plane
    double mop_kpyy[simtime * nplane];

    int mop_time;                     // current time index
} Pressure;

// Command-line parameters shared by mother() and child()
struct RunParameters
{
    int    integrator;
    int    c_pressure;
    int    nchild;
    double PE;
    int    index;
    double gammadot;
};


/* ============================================================================
 *  Function declarations
 * ========================================================================== */

// shared
void pbcx(double *x);
void pbcz(double *z);
void setparticletypes(Phase *phase);
void neighbourlist(Nlist *nlist, Phase *phase);
void wcaforce(Phase *phase, Interact *interact, Nlist *nlist);
void springforce(Phase *phase, Interact *interact);
void thermostatforce(Phase *phase, Interact *interact);
void force(Phase *phase, Interact *interact, Nlist *nlist);
void RK4(Phase *phase, Interact *interact, Nlist *nlist);
void VelVerlet(Phase *phase, Interact *interact, Nlist *nlist);

// mother
void mother(const RunParameters &par);
void rand0(Phase *phase);
int  randomInteger(int lower, int upper);
void setLattice(Phase *phase);
void randomMomenta(Phase *phase, int min, int max);
void printphase(Phase *phase, int i);

// child
void child(const RunParameters &par, int k);
void importinstate(Phase *phase, int k);
void mappings(Phase *phase);
void MOP_pressure(Phase *phase, Pressure *pressure, Nlist *nlist);
void computedissipation(Phase *phase, Thermo *thermo);
void printonfile(Phase *phase, Pressure *pressure, Thermo *thermo, FILE *writemop, FILE *velprof);


/* ============================================================================
 *  Main
 * ========================================================================== */

int main(int argc, char **argv)
{
    if (argc < 7)
    {
        cout << "usage: " << argv[0] << " integrator c_pressure nchild PE index gammadot" << endl;
        return 1;
    }

    RunParameters par;

    stringstream ss1(argv[1]);
    ss1 >> par.integrator;

    stringstream ss2(argv[2]);
    ss2 >> par.c_pressure;

    stringstream ss3(argv[3]);
    ss3 >> par.nchild;

    stringstream ss4(argv[4]);
    ss4 >> par.PE;

    stringstream ss5(argv[5]);
    ss5 >> par.index;

    stringstream ss6(argv[6]);
    ss6 >> par.gammadot;

    if (par.c_pressure != 1)
    {
        par.c_pressure = 0;
    }

    mother(par);   // runs the mother and, after every snapshot, the child

    return 0;
}


/* ============================================================================
 *  MOTHER
 * ========================================================================== */

// Equilibrium trajectory. After `warmup` steps, every `delay` steps it saves a
// snapshot (child<counter>.bin) and runs child() on it, then resumes.
void mother(const RunParameters &par)
{
    Phase    *phase    = new Phase();     // () = zero-initialised
    Interact *interact = new Interact();
    Nlist    *nlist    = new Nlist();

    phase->integrator = par.integrator;
    phase->c_pressure = par.c_pressure;
    phase->nchild     = par.nchild;
    phase->PE         = par.PE;
    phase->index      = par.index;

    phase->FP = phase->PE * Lx * Lz / nwall;   // pressure -> force per lattice site

    setLattice(phase);              // initial positions (fcc fluid between fcc walls)
    randomMomenta(phase, 0, 100);   // initial momenta at temperature T

    int counter = 0;                // snapshots taken so far
    long long int i;

    // equilibration
    for (i = 0; i < warmup; i++)
    {
        if (phase->integrator == 1)
        {
            VelVerlet(phase, interact, nlist);
        }
        else
        {
            RK4(phase, interact, nlist);
        }
    }

    // production: snapshot + child every `delay` steps
    while (counter < phase->nchild)
    {
        if (phase->integrator == 1)
        {
            VelVerlet(phase, interact, nlist);
        }
        else
        {
            RK4(phase, interact, nlist);
        }

        if (i % delay == 0)
        {
            printphase(phase, counter);   // snapshot -> child<counter>.bin
            child(par, counter);          // run the 4 mapped trajectories from it
            counter++;
        }

        i++;
    }

    delete phase;
    delete interact;
    delete nlist;
}

// Seed rand() with the seed on line `index` of seed.txt (clock if the file is missing)
void rand0(Phase *phase)
{
    int seed;
    int i = 1;
    ifstream fseed("seed.txt");

    if (!fseed.is_open())
    {
        cout << "There was a problem opening the input file\n";
        srand((int)time(NULL));
    }
    else
    {
        // skip the first index-1 seeds
        while (i < phase->index)
        {
            fseed >> seed;
            i++;
        }
        fseed >> seed;
        cout << seed << endl;
        srand(seed);
    }

    fseed.close();
}

// Uniform random integer in [lower, upper]
int randomInteger(int lower, int upper)
{
    double f;
    f = (double)rand() / ((double)RAND_MAX + 1);
    return (lower + (int)(f * (upper - lower + 1)));
}

// Initial positions.
// Fluid: fcc lattice of nxside x nyside x nzside cells (ttNf sites); the
//        ttNf - nfluid sites in the middle of the list are left empty.
// Walls: fcc slabs of nwx x nwy x nwz cells, placed just below (lower wall)
//        and just above (upper wall) the fluid, at a distance
//        (sigmaFF + sigmaWW)/2 from the last fluid layer.
// The lattice sites of the wall springs start at the initial wall positions.
void setLattice(Phase *phase)
{
    int i, j, k, c;

    double tx[ttNf];   // fluid lattice sites
    double ty[ttNf];
    double tz[ttNf];

    // --- fluid fcc lattice: two half-layers per cell along y ---
    c = 0;
    for (i = 1; i <= nyside; i++)
    {
        for (j = 1; j <= nzside; j++)
        {
            for (k = 1; k <= nxside; k++)
            {
                tx[c] = (k - 1) * cellf;
                ty[c] = (i - 1) * cellf;
                tz[c] = (j - 1) * cellf;
                c++;

                tx[c] = (k - 0.5) * cellf;
                ty[c] = (i - 1) * cellf;
                tz[c] = (j - 0.5) * cellf;
                c++;
            }
        }

        for (j = 1; j <= nzside; j++)
        {
            for (k = 1; k <= nxside; k++)
            {
                tx[c] = (k - 0.5) * cellf;
                ty[c] = (i - 0.5) * cellf;
                tz[c] = (j - 1) * cellf;
                c++;

                tx[c] = (k - 1) * cellf;
                ty[c] = (i - 0.5) * cellf;
                tz[c] = (j - 0.5) * cellf;
                c++;
            }
        }
    }

    // --- keep nfluid sites: drop the dN sites in the middle of the list ---
    int dN = ttNf - nfluid;

    int xmin = floor((ttNf - dN) / 2);
    int xmax = floor((ttNf + dN) / 2);

    int cc = nwall;

    for (c = 0; c < ttNf; c++)
    {
        if (c < xmin || c >= xmax)
        {
            phase->x[cc] = tx[c];
            phase->y[cc] = ty[c];
            phase->z[cc] = tz[c];
            cc++;
        }
    }

    // --- lower wall: fcc slab, particles [0, nwall) ---
    c = 0;
    for (i = 1; i <= nwy; i++)
    {
        for (j = 1; j <= nwz; j++)
        {
            for (k = 1; k <= nwx; k++)
            {
                phase->x[c] = (k - 1) * wcell;
                phase->y[c] = (i - 1) * wcell;
                phase->z[c] = (j - 1) * wcell;
                c++;

                phase->x[c] = (k - 0.5) * wcell;
                phase->y[c] = (i - 1) * wcell;
                phase->z[c] = (j - 0.5) * wcell;
                c++;
            }
        }

        for (j = 1; j <= nwz; j++)
        {
            for (k = 1; k <= nwx; k++)
            {
                phase->x[c] = (k - 0.5) * wcell;
                phase->y[c] = (i - 0.5) * wcell;
                phase->z[c] = (j - 1) * wcell;
                c++;

                phase->x[c] = (k - 1) * wcell;
                phase->y[c] = (i - 0.5) * wcell;
                phase->z[c] = (j - 0.5) * wcell;
                c++;
            }
        }
    }

    // --- upper wall: same slab, particles [nwall+nfluid, npart) ---
    c = nwall + nfluid;
    for (i = 1; i <= nwy; i++)
    {
        for (j = 1; j <= nwz; j++)
        {
            for (k = 1; k <= nwx; k++)
            {
                phase->x[c] = (k - 1) * wcell;
                phase->y[c] = (i - 1) * wcell;
                phase->z[c] = (j - 1) * wcell;
                c++;

                phase->x[c] = (k - 0.5) * wcell;
                phase->y[c] = (i - 1) * wcell;
                phase->z[c] = (j - 0.5) * wcell;
                c++;
            }
        }

        for (j = 1; j <= nwz; j++)
        {
            for (k = 1; k <= nwx; k++)
            {
                phase->x[c] = (k - 0.5) * wcell;
                phase->y[c] = (i - 0.5) * wcell;
                phase->z[c] = (j - 1) * wcell;
                c++;

                phase->x[c] = (k - 1) * wcell;
                phase->y[c] = (i - 0.5) * wcell;
                phase->z[c] = (j - 0.5) * wcell;
                c++;
            }
        }
    }

    // --- centre the box in x and z ---
    for (i = 0; i < npart; i++)
    {
        phase->x[i] += -Lx2;
        phase->z[i] += -Lz2;
    }

    // --- centre the fluid in y ---
    double L2 = (ty[ttNf - 1] - ty[0]) / 2;

    for (i = nwall; i < nwall + nfluid; i++)
    {
        phase->y[i] += -L2;
    }

    // --- place the walls below and above the fluid ---
    for (i = 0; i < nwall; i++)
    {
        phase->y[i] += phase->y[nwall] - phase->y[nwall - 1] - (sigmaFF + sigmaWW) / 2;
    }

    for (i = nwall + nfluid; i < npart; i++)
    {
        phase->y[i] += phase->y[nwall + nfluid - 1] + (sigmaFF + sigmaWW) / 2;
    }

    // --- lattice sites of the wall springs ---
    for (i = 0; i < nlatt; i++)
    {
        phase->xl[i] = phase->x[i];
        phase->yl[i] = phase->y[i];
        phase->zl[i] = phase->z[i];

        phase->xh[i] = phase->x[i + nwall + nfluid];
        phase->yh[i] = phase->y[i + nwall + nfluid];
        phase->zh[i] = phase->z[i + nwall + nfluid];
    }

    phase->wallvel = phase->y[nwall + nfluid] * gammadot;   // 0 for the equilibrium mother

    setparticletypes(phase);
}

// Initial momenta: random integers in [min, max], zero total momentum,
// rescaled to temperature T with 3*npart - dfree degrees of freedom.
// The barostat momenta start at 0.
void randomMomenta(Phase *phase, int min, int max)
{
    double pxsum = 0, pysum = 0, pzsum = 0, kineticTemp = 0, scale;
    int i;

    rand0(phase);

    for (i = 0; i < npart; i++)
    {
        phase->px[i] = (double)randomInteger(min, max);
        phase->py[i] = (double)randomInteger(min, max);
        phase->pz[i] = (double)randomInteger(min, max);
    }

    // remove the total momentum
    for (i = 0; i < npart; i++)
    {
        pxsum += phase->px[i];
        pysum += phase->py[i];
        pzsum += phase->pz[i];
    }

    pxsum /= npart;
    pysum /= npart;
    pzsum /= npart;

    for (i = 0; i < npart; i++)
    {
        phase->px[i] -= pxsum;
        phase->py[i] -= pysum;
        phase->pz[i] -= pzsum;
    }

    // rescale to the target temperature
    for (i = 0; i < npart; i++)
    {
        kineticTemp += phase->px[i] * phase->px[i] + phase->py[i] * phase->py[i] + phase->pz[i] * phase->pz[i];
    }

    scale = sqrt((3 * npart - dfree) * T / kineticTemp);

    for (i = 0; i < npart; i++)
    {
        phase->px[i] *= scale;
        phase->py[i] *= scale;
        phase->pz[i] *= scale;
    }

    for (i = 0; i < nlatt; i++)
    {
        phase->pyh[i] = 0;
    }
}

// Save a snapshot to child<i>.bin (binary, doubles):
//   x y z px py pz of every particle,
//   xl yl zl of every lower-wall site,
//   xh yh zh pxh pyh pzh of every upper-wall site,
//   alpha[0] alpha[1]
void printphase(Phase *phase, int i)
{
    string filename = "child" + to_string(i) + ".bin";

    FILE *writephase = fopen(filename.c_str(), "wb");

    for (i = 0; i < npart; i++)
    {
        fwrite(&phase->x[i], sizeof(double), 1, writephase);
        fwrite(&phase->y[i], sizeof(double), 1, writephase);
        fwrite(&phase->z[i], sizeof(double), 1, writephase);

        fwrite(&phase->px[i], sizeof(double), 1, writephase);
        fwrite(&phase->py[i], sizeof(double), 1, writephase);
        fwrite(&phase->pz[i], sizeof(double), 1, writephase);
    }

    for (i = 0; i < nlatt; i++)
    {
        fwrite(&phase->xl[i], sizeof(double), 1, writephase);
        fwrite(&phase->yl[i], sizeof(double), 1, writephase);
        fwrite(&phase->zl[i], sizeof(double), 1, writephase);
    }

    for (i = 0; i < nlatt; i++)
    {
        fwrite(&phase->xh[i], sizeof(double), 1, writephase);
        fwrite(&phase->yh[i], sizeof(double), 1, writephase);
        fwrite(&phase->zh[i], sizeof(double), 1, writephase);
        fwrite(&phase->pxh[i], sizeof(double), 1, writephase);
        fwrite(&phase->pyh[i], sizeof(double), 1, writephase);
        fwrite(&phase->pzh[i], sizeof(double), 1, writephase);
    }

    fwrite(&phase->alpha[0], sizeof(double), 1, writephase);
    fwrite(&phase->alpha[1], sizeof(double), 1, writephase);

    fclose(writephase);
}


/* ============================================================================
 *  CHILD
 * ========================================================================== */

// Run the 4 mapped sheared trajectories from snapshot k, each for `simtime`
// steps, recording the method-of-planes pressure, the wall work and the profiles.
void child(const RunParameters &par, int k)
{
    Phase    *phase    = new Phase();      // separate from the mother's structures
    Interact *interact = new Interact();
    Nlist    *nlist    = new Nlist();
    Pressure *pressure = new Pressure();
    Thermo   *thermo   = new Thermo();

    phase->integrator = par.integrator;
    phase->c_pressure = par.c_pressure;
    phase->gammadot   = par.gammadot;
    phase->PE         = par.PE;

    phase->FP = phase->PE * Lx * Lz / nwall;

    long long int i, j, im;

    // output files of this child
    string name = (phase->integrator == 1) ? "outputVV" : "outputRK";
    string filename = name + to_string(phase->c_pressure) + "_child" + to_string(k) + ".bin";

    FILE *writemop = fopen(filename.c_str(), "wb");
    FILE *velprof  = fopen(("velprofile_child" + to_string(k) + ".bin").c_str(), "wb");

    for (im = 1; im <= 4; im++)
    {
        pressure->mop_time = 0;
        thermo->omtime = 0;
        phase->mapping = im;

        importinstate(phase, k);   // reload the snapshot
        mappings(phase);           // apply mapping im

        // reset the accumulators of this trajectory
        for (i = 0; i < simtime * nplane; i++)
        {
            pressure->mop_pxy[i] = 0;
            pressure->mop_pyy[i] = 0;

            pressure->mop_kpxy[i] = 0;
            pressure->mop_kpyy[i] = 0;
        }

        for (i = 0; i < simtime; i++)
        {
            thermo->oml[i] = 0;
            thermo->omh[i] = 0;

            thermo->om2l[i] = 0;
            thermo->om2h[i] = 0;

            thermo->barostat[i] = 0;
        }

        for (i = 0; i < simtime; i++)
        {
            // method-of-planes planes (they follow the walls):
            //   0          : between the two layers of the lower wall
            //   1          : fluid side of the lower wall
            //   nplane-2   : fluid side of the upper wall
            //   nplane-1   : between the two layers of the upper wall
            //   2..nplane-3: equally spaced across the fluid
            phase->ypl[0] = (phase->yl[nlatt / 2 - 1] + phase->yl[nlatt / 2]) / 2;
            phase->ypl[1] = phase->yl[nlatt - 1] + sigmaWW / 2;

            phase->ypl[nplane - 2] = phase->yh[0] - sigmaWW / 2;
            phase->ypl[nplane - 1] = (phase->yh[nlatt / 2 - 1] + phase->yh[nlatt / 2]) / 2;

            for (j = 2; j < nplane - 2; j++)
            {
                phase->ypl[j] = (phase->ypl[nplane - 2] - phase->ypl[1]) / (nplane - 3) * (j - 1) + phase->ypl[1];
            }

            MOP_pressure(phase, pressure, nlist);   // pressures and profiles at this step
            computedissipation(phase, thermo);      // wall work at this step

            if (phase->integrator == 1)
            {
                VelVerlet(phase, interact, nlist);
            }
            else
            {
                RK4(phase, interact, nlist);
            }
        }

        printonfile(phase, pressure, thermo, writemop, velprof);   // append this mapping
    }

    fclose(writemop);
    fclose(velprof);

    delete phase;
    delete interact;
    delete nlist;
    delete pressure;
    delete thermo;
}

// Load snapshot child<k>.bin (format of printphase), set the wall velocity
// from the shear rate and the gap between the walls, and set the particle types.
void importinstate(Phase *phase, int k)
{
    int i;
    string filename = "child" + to_string(k) + ".bin";
    FILE *readphase = fopen(filename.c_str(), "rb");

    for (i = 0; i < npart; i++)
    {
        fread(&phase->x[i], sizeof(double), 1, readphase);
        fread(&phase->y[i], sizeof(double), 1, readphase);
        fread(&phase->z[i], sizeof(double), 1, readphase);

        fread(&phase->px[i], sizeof(double), 1, readphase);
        fread(&phase->py[i], sizeof(double), 1, readphase);
        fread(&phase->pz[i], sizeof(double), 1, readphase);
    }

    for (i = 0; i < nlatt; i++)
    {
        fread(&phase->xl[i], sizeof(double), 1, readphase);
        fread(&phase->yl[i], sizeof(double), 1, readphase);
        fread(&phase->zl[i], sizeof(double), 1, readphase);
    }

    for (i = 0; i < nlatt; i++)
    {
        fread(&phase->xh[i], sizeof(double), 1, readphase);
        fread(&phase->yh[i], sizeof(double), 1, readphase);
        fread(&phase->zh[i], sizeof(double), 1, readphase);
        fread(&phase->pxh[i], sizeof(double), 1, readphase);
        fread(&phase->pyh[i], sizeof(double), 1, readphase);
        fread(&phase->pzh[i], sizeof(double), 1, readphase);
    }

    fread(&phase->alpha[0], sizeof(double), 1, readphase);
    fread(&phase->alpha[1], sizeof(double), 1, readphase);

    fclose(readphase);

    // half the distance between the innermost site layers of the two walls
    double DL = (phase->yh[0] - phase->yl[nlatt - 1]) / 2;

    phase->wallvel = DL * phase->gammadot;

    setparticletypes(phase);
}

// Phase-space mappings applied to the starting state (TTCF)
//   1: identity
//   2: time reversal               p -> -p (particles and barostat)
//   3: reflection in x             x -> -x, px -> -px (particles and wall sites)
//   4: reflection in x + time rev. x -> -x, py -> -py, pz -> -pz
void mappings(Phase *phase)
{
    int i;

    if (phase->mapping == 2)
    {
        for (i = 0; i < npart; i++)
        {
            phase->px[i] = -phase->px[i];
            phase->py[i] = -phase->py[i];
            phase->pz[i] = -phase->pz[i];
        }
        for (i = 0; i < nlatt; i++)
        {
            phase->pxh[i] = -phase->pxh[i];
            phase->pyh[i] = -phase->pyh[i];
            phase->pzh[i] = -phase->pzh[i];
        }
    }

    if (phase->mapping == 3)
    {
        for (i = 0; i < npart; i++)
        {
            phase->x[i] = -phase->x[i];
            phase->px[i] = -phase->px[i];
        }

        for (i = 0; i < nlatt; i++)
        {
            phase->xl[i] = -phase->xl[i];
            phase->xh[i] = -phase->xh[i];

            phase->pxh[i] = -phase->pxh[i];
        }
    }

    if (phase->mapping == 4)
    {
        for (i = 0; i < npart; i++)
        {
            phase->x[i] = -phase->x[i];
        }

        for (i = 0; i < nlatt; i++)
        {
            phase->xl[i] = -phase->xl[i];
            phase->xh[i] = -phase->xh[i];
        }

        for (i = 0; i < npart; i++)
        {
            phase->py[i] = -phase->py[i];
            phase->pz[i] = -phase->pz[i];
        }

        for (i = 0; i < nlatt; i++)
        {
            phase->pyh[i] = -phase->pyh[i];
        }
    }
}

// Method-of-planes pressure on the nplane planes, and the profiles, at the
// current step t = pressure->mop_time.
//   Configurational term: sum over the pairs that straddle the plane of
//     4 * F * sign(ry)   (the factor 4 is divided out in printonfile).
//   Kinetic term: sum over the particles that crossed the plane during the
//     last step of the peculiar momenta v_x*v_y and v_y*v_y, linearly
//     interpolated to the crossing time.
//   Profiles: streaming velocity, temperature and number of particles per bin.
void MOP_pressure(Phase *phase, Pressure *pressure, Nlist *nlist)
{
    int i, j, ip;
    double rx, ry, rz, rsi, r6, rsq;
    double rm, sm, em;
    double F;
    double sg;

    int cx, cy, cz, ic, cp;
    int nx, ny, nz, nx1, ny1, nz1, in;
    int pc, pn;

    int t = pressure->mop_time;

    neighbourlist(nlist, phase);

    // ------------------------------------------------ configurational term
    for (ip = 0; ip < nplane; ip++)
    {
        cp = (phase->ypl[ip] + Ly) / dcy;     // cell row containing the plane

        for (cx = 0; cx < clx; cx++)
        {
            for (cy = cp - 1; cy <= cp + 1; cy++)   // only rows next to the plane
            {
                for (cz = 0; cz < clz; cz++)
                {
                    ic = cx * cly * clz + cy * clz + cz;

                    if (nlist->atom[ic] != -1)
                    {
                        // neighbouring cells, periodic in x and z
                        for (nx = cx - 1; nx <= cx + 1; nx++)
                        {
                            for (ny = cy - 1; ny <= cy + 1; ny++)
                            {
                                for (nz = cz - 1; nz <= cz + 1; nz++)
                                {
                                    nx1 = nx;
                                    ny1 = ny;
                                    nz1 = nz;
                                    if (nx1 == -1)  nx1 += clx;
                                    if (nx1 == clx) nx1 -= clx;
                                    if (nz1 == -1)  nz1 += clz;
                                    if (nz1 == clz) nz1 -= clz;

                                    in = nx1 * cly * clz + ny1 * clz + nz1;

                                    // pairs (pc, pn) with pn < pc: each pair once
                                    pc = nlist->atom[ic];
                                    while (pc != -1)
                                    {
                                        pn = nlist->atom[in];
                                        while (pn < pc && pn != -1)
                                        {
                                            rx = phase->x[pc] - phase->x[pn];
                                            ry = phase->y[pc] - phase->y[pn];
                                            rz = phase->z[pc] - phase->z[pn];

                                            // only pairs on opposite sides of the plane
                                            if ((phase->y[pc] - phase->ypl[ip]) * (phase->y[pn] - phase->ypl[ip]) < 0)
                                            {
                                                pbcx(&rx);
                                                pbcz(&rz);

                                                rsq = rx * rx + ry * ry + rz * rz;
                                                rm = phase->reff[phase->ptype[pc]][phase->ptype[pn]];

                                                if (rsq < rm * rm)
                                                {
                                                    sm = phase->seff[phase->ptype[pc]][phase->ptype[pn]];
                                                    em = phase->eeff[phase->ptype[pc]][phase->ptype[pn]];

                                                    rsi = sm * sm / rsq;
                                                    r6  = rsi * rsi * rsi;
                                                    rsi = 1. / rsq;

                                                    F = 24. * em * rsi * r6 * (2. * r6 - 1.);   // force / r

                                                    sg = copysignf(1.0, ry);

                                                    pressure->mop_pxy[t * nplane + ip] += 4 * rx * F * sg;
                                                    pressure->mop_pyy[t * nplane + ip] += 4 * ry * F * sg;
                                                }
                                            }
                                            pn = phase->plist[pn];
                                        }
                                        pc = phase->plist[pc];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // ------------------------------------------------ profiles
    double t0, vx0, vy0;

    for (i = 0; i < nbin; i++)
    {
        pressure->mop_vxstream[i] = 0;
        pressure->mop_vystream[i] = 0;
        pressure->mop_T[i] = 0;
        pressure->mop_n[i] = 0;
    }

    // streaming velocity per bin
    for (i = 0; i < npart; i++)
    {
        ip = (phase->y[i] + Ly) / bin;
        pressure->mop_vxstream[ip] += phase->px[i];
        pressure->mop_vystream[ip] += phase->py[i];
        pressure->mop_n[ip]++;
    }

    for (i = 0; i < nbin; i++)
    {
        if (pressure->mop_n[i] != 0)
        {
            pressure->mop_vxstream[i] /= pressure->mop_n[i];
            pressure->mop_vystream[i] /= pressure->mop_n[i];
        }
    }

    // peculiar momenta and temperature per bin
    for (i = 0; i < npart; i++)
    {
        ip = (phase->y[i] + Ly) / bin;
        phase->pxth[i] = phase->px[i] - pressure->mop_vxstream[ip];
        phase->pyth[i] = phase->py[i] - pressure->mop_vystream[ip];
        pressure->mop_T[ip] += phase->pxth[i] * phase->pxth[i] + phase->pyth[i] * phase->pyth[i] + phase->pz[i] * phase->pz[i];
    }

    for (i = 0; i < nbin; i++)
    {
        if (pressure->mop_n[i] != 0)
        {
            pressure->mop_T[i] /= (3 * pressure->mop_n[i]);
        }
    }

    // ------------------------------------------------ kinetic term
    // particles that crossed a plane between the previous step and this one
    if (t > 0)
    {
        double tin = (t - 1) * dt;   // time of the previous step

        for (ip = 0; ip < nplane; ip++)
        {
            for (i = 0; i < npart; i++)
            {
                if ((phase->yp[i] - phase->ypl[ip]) * (phase->y[i] - phase->ypl[ip]) <= 0)
                {
                    // crossing time and momenta, by linear interpolation
                    t0  = (phase->yp[i] - phase->ypl[ip]) / ((phase->yp[i] - phase->ypl[ip]) - (phase->y[i] - phase->ypl[ip])) * dt + tin;
                    vx0 = phase->pxthp[i] + (phase->pxth[i] - phase->pxthp[i]) / dt * (t0 - tin);
                    vy0 = phase->pythp[i] + (phase->pyth[i] - phase->pythp[i]) / dt * (t0 - tin);

                    pressure->mop_kpxy[t * nplane + ip] += vx0 * vy0;
                    pressure->mop_kpyy[t * nplane + ip] += vy0 * vy0;
                }
            }
        }
    }

    // keep this step for the next crossing check
    for (j = 0; j < npart; j++)
    {
        phase->yp[j] = phase->y[j];
        phase->pxthp[j] = phase->pxth[j];
        phase->pythp[j] = phase->pyth[j];
    }

    // store the profiles of this step
    for (i = 0; i < nbin; i++)
    {
        pressure->vstream[t * nbin + i] = pressure->mop_vxstream[i];
        pressure->temp[t * nbin + i] = pressure->mop_T[i];
        pressure->n[t * nbin + i] = pressure->mop_n[i];
    }

    pressure->mop_time++;
}

// Power delivered by the moving lattice sites through the springs, for each
// wall: sum of wallvel * K * (x - x_site), with the sign of the site velocity.
// Also stores the power of the barostat piston.
void computedissipation(Phase *phase, Thermo *thermo)
{
    int t = thermo->omtime;
    int i;

    double ol = 0;
    double oh = 0;
    double o2l = 0;
    double o2h = 0;

    double rx;

    // lower wall (sites move at -wallvel)
    for (i = 0; i < nwall; i++)
    {
        rx = phase->x[i] - phase->xl[i];
        pbcx(&rx);
        ol += phase->wallvel * K * rx;
    }

    // upper wall (sites move at +wallvel)
    for (i = 0; i < nwall; i++)
    {
        rx = phase->x[i + nwall + nfluid] - phase->xh[i];
        pbcx(&rx);
        oh += -phase->wallvel * K * rx;
    }

    thermo->oml[t] = ol;
    thermo->omh[t] = oh;

    thermo->om2l[t] = o2l;
    thermo->om2h[t] = o2h;

    thermo->barostat[t] = nlatt * phase->pyh[0] / M * phase->FP;

    thermo->omtime++;
}

// Append the results of one mapped trajectory.
//   writemop: for each step, for each plane: P_xy, P_yy (configurational +
//             kinetic, per unit area); then oml, omh.
//             Record size per step: 2*nplane + 2 doubles.
//   velprof:  for each step, for each bin: streaming velocity, temperature,
//             number density. Record size per step: 3*nbin doubles.
void printonfile(Phase *phase, Pressure *pressure, Thermo *thermo, FILE *writemop, FILE *velprof)
{
    int i, j;
    double PP;

    for (i = 0; i < simtime; i++)
    {
        for (j = 0; j < nplane; j++)
        {
            PP = pressure->mop_pxy[i * nplane + j] / (4 * Lx * Lz) + pressure->mop_kpxy[i * nplane + j] / (dt * Lx * Lz);
            fwrite(&PP, sizeof(double), 1, writemop);
            PP = pressure->mop_pyy[i * nplane + j] / (4 * Lx * Lz) + pressure->mop_kpyy[i * nplane + j] / (dt * Lx * Lz);
            fwrite(&PP, sizeof(double), 1, writemop);
        }

        fwrite(&thermo->oml[i], sizeof(double), 1, writemop);
        fwrite(&thermo->omh[i], sizeof(double), 1, writemop);
    }

    for (i = 0; i < nbin * simtime; i++)
    {
        PP = pressure->vstream[i];
        fwrite(&PP, sizeof(double), 1, velprof);
        PP = pressure->temp[i];
        fwrite(&PP, sizeof(double), 1, velprof);
        PP = pressure->n[i] / vbin;
        fwrite(&PP, sizeof(double), 1, velprof);
    }
}


/* ============================================================================
 *  SHARED: boundary conditions, particle types, cell list
 * ========================================================================== */

// Wrap an x coordinate (or separation) into [-Lx/2, Lx/2)
void pbcx(double *x)
{
    if (*x <= -Lx2)
    {
        *x += Lx;
    }
    else if (*x >= Lx2)
    {
        *x -= Lx;
    }
}

// Wrap a z coordinate (or separation) into [-Lz/2, Lz/2)
void pbcz(double *z)
{
    if (*z <= -Lz2)
    {
        *z += Lz;
    }
    else if (*z >= Lz2)
    {
        *z -= Lz;
    }
}

// Particle types (0 = wall, 1 = fluid) and the interaction table of each pair of types
void setparticletypes(Phase *phase)
{
    int i;

    for (i = 0; i < nwall; i++)
    {
        phase->ptype[i] = 0;
    }
    for (i = nwall; i < nwall + nfluid; i++)
    {
        phase->ptype[i] = 1;
    }
    for (i = nwall + nfluid; i < npart; i++)
    {
        phase->ptype[i] = 0;
    }

    phase->seff[0][0] = sigmaWW;
    phase->seff[1][0] = sigmaWF;
    phase->seff[0][1] = sigmaWF;
    phase->seff[1][1] = sigmaFF;

    phase->eeff[0][0] = epsWW;
    phase->eeff[1][0] = epsWF;
    phase->eeff[0][1] = epsWF;
    phase->eeff[1][1] = epsFF;

    phase->reff[0][0] = rWW;
    phase->reff[1][0] = rWF;
    phase->reff[0][1] = rWF;
    phase->reff[1][1] = rFF;
}

// Build the cell list. Each cell holds a linked list (nlist->atom = head,
// phase->plist = next); particles are appended at the tail, so each list is
// in increasing particle index.
void neighbourlist(Nlist *nlist, Phase *phase)
{
    int i;
    int ix, iy, iz;
    int ind, indc;

    for (i = 0; i < totlist; i++)
    {
        nlist->atom[i] = -1;
    }
    for (i = 0; i < npart; i++)
    {
        phase->plist[i] = -1;
    }

    for (i = 0; i < npart; i++)
    {
        // cell of particle i (y measured from -Ly)
        ix = (phase->x[i] + Lx2) / dcx;
        iy = (phase->y[i] + Ly) / dcy;
        iz = (phase->z[i] + Lz2) / dcz;

        indc = ix * cly * clz + iy * clz + iz;

        ind = nlist->atom[indc];

        if (ind == -1)
        {
            nlist->atom[indc] = i;          // first particle of the cell
        }
        else
        {
            while (phase->plist[ind] != -1) // walk to the tail and append
            {
                ind = phase->plist[ind];
            }
            phase->plist[ind] = i;
        }
    }
}


/* ============================================================================
 *  SHARED: forces
 * ========================================================================== */

// Pair forces (LJ cut at reff, with sigma and epsilon from the type table)
// using the cell list. Each pair is computed once (pn < pc) over the 27
// neighbouring cells; periodic in x and z, not in y.
void wcaforce(Phase *phase, Interact *interact, Nlist *nlist)
{
    double rx, ry, rz, rsi, r6, rsq;
    double rm, sm, em;
    double F;

    int cx, cy, cz, ic;
    int nx, ny, nz, nx1, ny1, nz1, in;
    int pc, pn;

    neighbourlist(nlist, phase);

    for (cx = 0; cx < clx; cx++)
    {
        for (cy = 0; cy < cly; cy++)
        {
            for (cz = 0; cz < clz; cz++)
            {
                ic = cx * cly * clz + cy * clz + cz;

                if (nlist->atom[ic] != -1)
                {
                    // neighbouring cells: periodic in x and z, clipped in y
                    for (nx = cx - 1; nx <= cx + 1; nx++)
                    {
                        for (ny = max(cy - 1, 0); ny <= min(cy + 1, cly - 1); ny++)
                        {
                            for (nz = cz - 1; nz <= cz + 1; nz++)
                            {
                                nx1 = nx;
                                ny1 = ny;
                                nz1 = nz;
                                if (nx1 == -1)  nx1 += clx;
                                if (nx1 == clx) nx1 -= clx;
                                if (nz1 == -1)  nz1 += clz;
                                if (nz1 == clz) nz1 -= clz;

                                in = nx1 * cly * clz + ny1 * clz + nz1;

                                pc = nlist->atom[ic];
                                while (pc != -1)
                                {
                                    pn = nlist->atom[in];
                                    while (pn < pc && pn != -1)
                                    {
                                        rx = phase->x[pc] - phase->x[pn];
                                        ry = phase->y[pc] - phase->y[pn];
                                        rz = phase->z[pc] - phase->z[pn];

                                        pbcx(&rx);
                                        pbcz(&rz);

                                        rsq = rx * rx + ry * ry + rz * rz;

                                        rm = phase->reff[phase->ptype[pc]][phase->ptype[pn]];

                                        if (rsq < rm * rm)
                                        {
                                            sm = phase->seff[phase->ptype[pc]][phase->ptype[pn]];
                                            em = phase->eeff[phase->ptype[pc]][phase->ptype[pn]];

                                            rsi = sm * sm / rsq;
                                            r6  = rsi * rsi * rsi;
                                            rsi = 1. / rsq;

                                            F = 24. * em * rsi * r6 * (2. * r6 - 1.);   // force / r

                                            interact->fx[pc] += rx * F;
                                            interact->fy[pc] += ry * F;
                                            interact->fz[pc] += rz * F;
                                            interact->fx[pn] -= rx * F;
                                            interact->fy[pn] -= ry * F;
                                            interact->fz[pn] -= rz * F;
                                        }
                                        pn = phase->plist[pn];
                                    }
                                    pc = phase->plist[pc];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// Harmonic springs tethering every wall particle to its lattice site.
// Also accumulates in fyh the mean reaction force of the upper-wall springs on
// their sites (used by the barostat).
void springforce(Phase *phase, Interact *interact)
{
    int i;
    double rx, ry, rz;

    // lower wall
    for (i = 0; i < nwall; i++)
    {
        rx = phase->x[i] - phase->xl[i];
        ry = phase->y[i] - phase->yl[i];
        rz = phase->z[i] - phase->zl[i];

        pbcx(&rx);
        pbcz(&rz);

        interact->fx[i] -= rx * K;
        interact->fy[i] -= ry * K;
        interact->fz[i] -= rz * K;
    }

    // upper wall
    for (i = 0; i < nwall; i++)
    {
        rx = phase->x[i + nfluid + nwall] - phase->xh[i];
        ry = phase->y[i + nfluid + nwall] - phase->yh[i];
        rz = phase->z[i + nfluid + nwall] - phase->zh[i];

        pbcx(&rx);
        pbcz(&rz);

        interact->fx[i + nfluid + nwall] -= rx * K;
        interact->fy[i + nfluid + nwall] -= ry * K;
        interact->fz[i + nfluid + nwall] -= rz * K;

        interact->fyh += ry * K;   // reaction on the site
    }

    interact->fyh /= nlatt;
}

// Nose-Hoover thermostats of the two walls:
//   d(alpha)/dt = (sum p'^2 - 3*nwall*T) / Q
// with p' the wall momenta relative to the mean x momentum of the wall.
// Q = QVV (velocity Verlet) or QRK (RK4).
void thermostatforce(Phase *phase, Interact *interact)
{
    int i;
    double pmh, pml, sumh, suml;

    pmh = 0;
    pml = 0;
    sumh = 0;
    suml = 0;

    // mean x momentum of each wall
    for (i = 0; i < nwall; i++)
    {
        pml += phase->px[i];
        pmh += phase->px[i + nfluid + nwall];
    }
    pml /= nwall;
    pmh /= nwall;

    // peculiar kinetic energy (x2) of each wall
    for (i = 0; i < nwall; i++)
    {
        suml += (phase->px[i] - pml) * (phase->px[i] - pml) + phase->py[i] * phase->py[i] + phase->pz[i] * phase->pz[i];
        sumh += (phase->px[i + nfluid + nwall] - pmh) * (phase->px[i + nfluid + nwall] - pmh) + phase->py[i + nfluid + nwall] * phase->py[i + nfluid + nwall] + phase->pz[i + nfluid + nwall] * phase->pz[i + nfluid + nwall];
    }

    if (phase->integrator == 1)
    {
        interact->dalpha[0] = (suml - 3 * nwall * T) / QVV;
        interact->dalpha[1] = (sumh - 3 * nwall * T) / QVV;
    }
    else
    {
        interact->dalpha[0] = (suml - 3 * nwall * T) / QRK;
        interact->dalpha[1] = (sumh - 3 * nwall * T) / QRK;
    }
}

// All forces at the current state:
//   pair forces, wall springs, thermostat derivatives, and the net force on the
//   upper-wall sites (mean spring reaction minus the imposed FP, switched off
//   when c_pressure = 0).
void force(Phase *phase, Interact *interact, Nlist *nlist)
{
    int i;

    for (i = 0; i < npart; i++)
    {
        interact->fx[i] = 0.;
        interact->fy[i] = 0.;
        interact->fz[i] = 0.;
    }

    interact->dalpha[0] = 0.;
    interact->dalpha[1] = 0.;

    interact->fyh = 0.;

    wcaforce(phase, interact, nlist);
    springforce(phase, interact);
    thermostatforce(phase, interact);

    interact->fyh -= phase->FP;
    interact->fyh *= phase->c_pressure;
}


/* ============================================================================
 *  SHARED: integrators
 * ========================================================================== */

// One step of classical 4th-order Runge-Kutta for positions, momenta,
// thermostat variables and the barostat (y and py of the upper-wall sites).
// The wall sites move by wallvel*dt/2 before stage 2 and again before stage 4.
void RK4(Phase *phase, Interact *interact, Nlist *nlist)
{
    int i;

    // ---------------------------------------------------- state at the start
    for (i = 0; i < npart; i++)
    {
        phase->x0[i] = phase->x[i];
        phase->y0[i] = phase->y[i];
        phase->z0[i] = phase->z[i];

        phase->px0[i] = phase->px[i];
        phase->py0[i] = phase->py[i];
        phase->pz0[i] = phase->pz[i];
    }

    phase->alpha0[0] = phase->alpha[0];
    phase->alpha0[1] = phase->alpha[1];

    for (i = 0; i < nlatt; i++)
    {
        phase->yh0[i] = phase->yh[i];
        phase->pyh0[i] = phase->pyh[i];
    }

    // ---------------------------------------------------- stage 1
    force(phase, interact, nlist);

    for (i = 0; i < npart; i++)
    {
        phase->x1[i] = dt * phase->px[i];
        phase->y1[i] = dt * phase->py[i];
        phase->z1[i] = dt * phase->pz[i];

        phase->px1[i] = dt * interact->fx[i];
        phase->py1[i] = dt * interact->fy[i];
        phase->pz1[i] = dt * interact->fz[i];
    }

    phase->alpha1[0] = dt * interact->dalpha[0];
    phase->alpha1[1] = dt * interact->dalpha[1];

    for (i = 0; i < nlatt; i++)
    {
        phase->yh1[i] = dt * phase->pyh[i];
        phase->pyh1[i] = dt * interact->fyh;
    }

    // midpoint state for stage 2
    for (i = 0; i < npart; i++)
    {
        phase->x[i] = phase->x0[i] + 0.5 * phase->x1[i];
        phase->y[i] = phase->y0[i] + 0.5 * phase->y1[i];
        phase->z[i] = phase->z0[i] + 0.5 * phase->z1[i];

        phase->px[i] = phase->px0[i] + 0.5 * phase->px1[i];
        phase->py[i] = phase->py0[i] + 0.5 * phase->py1[i];
        phase->pz[i] = phase->pz0[i] + 0.5 * phase->pz1[i];
    }
    phase->alpha[0] = phase->alpha0[0] + 0.5 * phase->alpha1[0];
    phase->alpha[1] = phase->alpha0[1] + 0.5 * phase->alpha1[1];

    for (i = 0; i < npart; i++)
    {
        pbcx(&phase->x[i]);
        pbcz(&phase->z[i]);
    }

    // wall sites to the midpoint
    for (i = 0; i < nlatt; i++)
    {
        phase->xl[i] -= 0.5 * phase->wallvel * dt;
        phase->xh[i] += 0.5 * phase->wallvel * dt;
    }

    for (i = 0; i < nlatt; i++)
    {
        pbcx(&phase->xl[i]);
        pbcx(&phase->xh[i]);
    }

    for (i = 0; i < nlatt; i++)
    {
        phase->yh[i] = phase->yh0[i] + 0.5 * phase->yh1[i];
        phase->pyh[i] = phase->pyh0[i] + 0.5 * phase->pyh1[i];
    }

    // ---------------------------------------------------- stage 2
    force(phase, interact, nlist);

    for (i = 0; i < npart; i++)
    {
        phase->x2[i] = dt * phase->px[i];
        phase->y2[i] = dt * phase->py[i];
        phase->z2[i] = dt * phase->pz[i];

        phase->px2[i] = dt * interact->fx[i];
        phase->py2[i] = dt * interact->fy[i];
        phase->pz2[i] = dt * interact->fz[i];
    }

    phase->alpha2[0] = dt * interact->dalpha[0];
    phase->alpha2[1] = dt * interact->dalpha[1];

    for (i = 0; i < nlatt; i++)
    {
        phase->yh2[i] = dt * phase->pyh[i];
        phase->pyh2[i] = dt * interact->fyh;
    }

    // midpoint state for stage 3 (wall sites stay at the midpoint)
    for (i = 0; i < npart; i++)
    {
        phase->x[i] = phase->x0[i] + 0.5 * phase->x2[i];
        phase->y[i] = phase->y0[i] + 0.5 * phase->y2[i];
        phase->z[i] = phase->z0[i] + 0.5 * phase->z2[i];

        phase->px[i] = phase->px0[i] + 0.5 * phase->px2[i];
        phase->py[i] = phase->py0[i] + 0.5 * phase->py2[i];
        phase->pz[i] = phase->pz0[i] + 0.5 * phase->pz2[i];
    }

    phase->alpha[0] = phase->alpha0[0] + 0.5 * phase->alpha2[0];
    phase->alpha[1] = phase->alpha0[1] + 0.5 * phase->alpha2[1];

    for (i = 0; i < npart; i++)
    {
        pbcx(&phase->x[i]);
        pbcz(&phase->z[i]);
    }

    for (i = 0; i < nlatt; i++)
    {
        phase->yh[i] = phase->yh0[i] + 0.5 * phase->yh2[i];
        phase->pyh[i] = phase->pyh0[i] + 0.5 * phase->pyh2[i];
    }

    // ---------------------------------------------------- stage 3
    force(phase, interact, nlist);

    for (i = 0; i < npart; i++)
    {
        phase->x3[i] = dt * phase->px[i];
        phase->y3[i] = dt * phase->py[i];
        phase->z3[i] = dt * phase->pz[i];

        phase->px3[i] = dt * interact->fx[i];
        phase->py3[i] = dt * interact->fy[i];
        phase->pz3[i] = dt * interact->fz[i];
    }

    phase->alpha3[0] = dt * interact->dalpha[0];
    phase->alpha3[1] = dt * interact->dalpha[1];

    for (i = 0; i < nlatt; i++)
    {
        phase->yh3[i] = dt * phase->pyh[i];
        phase->pyh3[i] = dt * interact->fyh;
    }

    // end-point state for stage 4
    for (i = 0; i < npart; i++)
    {
        phase->x[i] = phase->x0[i] + phase->x3[i];
        phase->y[i] = phase->y0[i] + phase->y3[i];
        phase->z[i] = phase->z0[i] + phase->z3[i];

        phase->px[i] = phase->px0[i] + phase->px3[i];
        phase->py[i] = phase->py0[i] + phase->py3[i];
        phase->pz[i] = phase->pz0[i] + phase->pz3[i];
    }

    phase->alpha[0] = phase->alpha0[0] + phase->alpha3[0];
    phase->alpha[1] = phase->alpha0[1] + phase->alpha3[1];

    for (i = 0; i < npart; i++)
    {
        pbcx(&phase->x[i]);
        pbcz(&phase->z[i]);
    }

    // wall sites to the end of the step
    for (i = 0; i < nlatt; i++)
    {
        phase->xl[i] -= 0.5 * phase->wallvel * dt;
        phase->xh[i] += 0.5 * phase->wallvel * dt;
    }

    for (i = 0; i < nlatt; i++)
    {
        pbcx(&phase->xl[i]);
        pbcx(&phase->xh[i]);
    }

    for (i = 0; i < nlatt; i++)
    {
        phase->yh[i] = phase->yh0[i] + phase->yh3[i];
        phase->pyh[i] = phase->pyh0[i] + phase->pyh3[i];
    }

    // ---------------------------------------------------- stage 4
    force(phase, interact, nlist);

    for (i = 0; i < npart; i++)
    {
        phase->x4[i] = dt * phase->px[i];
        phase->y4[i] = dt * phase->py[i];
        phase->z4[i] = dt * phase->pz[i];

        phase->px4[i] = dt * interact->fx[i];
        phase->py4[i] = dt * interact->fy[i];
        phase->pz4[i] = dt * interact->fz[i];
    }

    phase->alpha4[0] = dt * interact->dalpha[0];
    phase->alpha4[1] = dt * interact->dalpha[1];

    for (i = 0; i < nlatt; i++)
    {
        phase->yh4[i] = dt * phase->pyh[i];
        phase->pyh4[i] = dt * interact->fyh;
    }

    // ---------------------------------------------------- RK4 update
    for (i = 0; i < npart; i++)
    {
        phase->x[i] = phase->x0[i] + (phase->x1[i] + 2 * phase->x2[i] + 2 * phase->x3[i] + phase->x4[i]) / 6;
        phase->y[i] = phase->y0[i] + (phase->y1[i] + 2 * phase->y2[i] + 2 * phase->y3[i] + phase->y4[i]) / 6;
        phase->z[i] = phase->z0[i] + (phase->z1[i] + 2 * phase->z2[i] + 2 * phase->z3[i] + phase->z4[i]) / 6;

        phase->px[i] = phase->px0[i] + (phase->px1[i] + 2 * phase->px2[i] + 2 * phase->px3[i] + phase->px4[i]) / 6;
        phase->py[i] = phase->py0[i] + (phase->py1[i] + 2 * phase->py2[i] + 2 * phase->py3[i] + phase->py4[i]) / 6;
        phase->pz[i] = phase->pz0[i] + (phase->pz1[i] + 2 * phase->pz2[i] + 2 * phase->pz3[i] + phase->pz4[i]) / 6;
    }

    phase->alpha[0] = phase->alpha0[0] + (phase->alpha1[0] + 2 * phase->alpha2[0] + 2 * phase->alpha3[0] + phase->alpha4[0]) / 6;
    phase->alpha[1] = phase->alpha0[1] + (phase->alpha1[1] + 2 * phase->alpha2[1] + 2 * phase->alpha3[1] + phase->alpha4[1]) / 6;

    for (i = 0; i < npart; i++)
    {
        pbcx(&phase->x[i]);
        pbcz(&phase->z[i]);
    }

    for (i = 0; i < nlatt; i++)
    {
        phase->yh[i] = phase->yh0[i] + (phase->yh1[i] + 2 * phase->yh2[i] + 2 * phase->yh3[i] + phase->yh4[i]) / 6;
        phase->pyh[i] = phase->pyh0[i] + (phase->pyh1[i] + 2 * phase->pyh2[i] + 2 * phase->pyh3[i] + phase->pyh4[i]) / 6;
    }
}

// One step of velocity Verlet with Nose-Hoover thermostats on the walls:
//   1. positions to t+dt (walls include the friction -alpha*p)
//   2. wall sites: x by +/-wallvel*dt, upper sites in y by the barostat
//   3. momenta half step with the old forces
//   4. alpha: two half steps around the force evaluation at t+dt
//   5. momenta second half step; for the walls the friction term is
//      treated implicitly:  p = (p + dt/2 F) / (1 + dt/2 alpha)
void VelVerlet(Phase *phase, Interact *interact, Nlist *nlist)
{
    int i;

    force(phase, interact, nlist);   // forces at t

    // ---------------------------------------------------- positions
    // lower wall
    for (i = 0; i < nwall; i++)
    {
        phase->x[i] += dt * phase->px[i] + 0.5 * dt * dt * (interact->fx[i] - phase->alpha[0] * phase->px[i]);
        phase->y[i] += dt * phase->py[i] + 0.5 * dt * dt * (interact->fy[i] - phase->alpha[0] * phase->py[i]);
        phase->z[i] += dt * phase->pz[i] + 0.5 * dt * dt * (interact->fz[i] - phase->alpha[0] * phase->pz[i]);
    }

    // fluid
    for (i = nwall; i < nwall + nfluid; i++)
    {
        phase->x[i] += dt * phase->px[i] + 0.5 * dt * dt * (interact->fx[i]);
        phase->y[i] += dt * phase->py[i] + 0.5 * dt * dt * (interact->fy[i]);
        phase->z[i] += dt * phase->pz[i] + 0.5 * dt * dt * (interact->fz[i]);
    }

    // upper wall
    for (i = nwall + nfluid; i < npart; i++)
    {
        phase->x[i] += dt * phase->px[i] + 0.5 * dt * dt * (interact->fx[i] - phase->alpha[1] * phase->px[i]);
        phase->y[i] += dt * phase->py[i] + 0.5 * dt * dt * (interact->fy[i] - phase->alpha[1] * phase->py[i]);
        phase->z[i] += dt * phase->pz[i] + 0.5 * dt * dt * (interact->fz[i] - phase->alpha[1] * phase->pz[i]);
    }

    // ---------------------------------------------------- wall sites
    for (i = 0; i < nlatt; i++)
    {
        phase->yh[i] += dt * phase->pyh[i] + 0.5 * dt * dt * interact->fyh;   // barostat
    }

    for (i = 0; i < nlatt; i++)
    {
        phase->xl[i] -= phase->wallvel * dt;
        phase->xh[i] += phase->wallvel * dt;
    }

    for (i = 0; i < npart; i++)
    {
        pbcx(&phase->x[i]);
        pbcz(&phase->z[i]);
    }

    for (i = 0; i < nlatt; i++)
    {
        pbcx(&phase->xl[i]);
        pbcx(&phase->xh[i]);
    }

    // ---------------------------------------------------- momenta, first half
    for (i = 0; i < nwall; i++)
    {
        phase->px[i] += 0.5 * dt * (interact->fx[i] - phase->alpha[0] * phase->px[i]);
        phase->py[i] += 0.5 * dt * (interact->fy[i] - phase->alpha[0] * phase->py[i]);
        phase->pz[i] += 0.5 * dt * (interact->fz[i] - phase->alpha[0] * phase->pz[i]);
    }

    for (i = nwall; i < nwall + nfluid; i++)
    {
        phase->px[i] += 0.5 * dt * (interact->fx[i]);
        phase->py[i] += 0.5 * dt * (interact->fy[i]);
        phase->pz[i] += 0.5 * dt * (interact->fz[i]);
    }

    for (i = nwall + nfluid; i < npart; i++)
    {
        phase->px[i] += 0.5 * dt * (interact->fx[i] - phase->alpha[1] * phase->px[i]);
        phase->py[i] += 0.5 * dt * (interact->fy[i] - phase->alpha[1] * phase->py[i]);
        phase->pz[i] += 0.5 * dt * (interact->fz[i] - phase->alpha[1] * phase->pz[i]);
    }

    for (i = 0; i < nlatt; i++)
    {
        phase->pyh[i] += 0.5 * dt * interact->fyh;
    }

    // ---------------------------------------------------- thermostat, first half
    phase->alpha[0] += 0.5 * dt * interact->dalpha[0];
    phase->alpha[1] += 0.5 * dt * interact->dalpha[1];

    force(phase, interact, nlist);   // forces at t+dt

    // ---------------------------------------------------- thermostat, second half
    phase->alpha[0] += 0.5 * dt * interact->dalpha[0];
    phase->alpha[1] += 0.5 * dt * interact->dalpha[1];

    // ---------------------------------------------------- momenta, second half
    for (i = 0; i < nwall; i++)
    {
        phase->px[i] = (phase->px[i] + 0.5 * dt * interact->fx[i]) / (1 + 0.5 * dt * phase->alpha[0]);
        phase->py[i] = (phase->py[i] + 0.5 * dt * interact->fy[i]) / (1 + 0.5 * dt * phase->alpha[0]);
        phase->pz[i] = (phase->pz[i] + 0.5 * dt * interact->fz[i]) / (1 + 0.5 * dt * phase->alpha[0]);
    }

    for (i = nwall; i < nwall + nfluid; i++)
    {
        phase->px[i] += 0.5 * dt * (interact->fx[i]);
        phase->py[i] += 0.5 * dt * (interact->fy[i]);
        phase->pz[i] += 0.5 * dt * (interact->fz[i]);
    }

    for (i = nwall + nfluid; i < npart; i++)
    {
        phase->px[i] = (phase->px[i] + 0.5 * dt * interact->fx[i]) / (1 + 0.5 * dt * phase->alpha[1]);
        phase->py[i] = (phase->py[i] + 0.5 * dt * interact->fy[i]) / (1 + 0.5 * dt * phase->alpha[1]);
        phase->pz[i] = (phase->pz[i] + 0.5 * dt * interact->fz[i]) / (1 + 0.5 * dt * phase->alpha[1]);
    }

    for (i = 0; i < nlatt; i++)
    {
        phase->pyh[i] += 0.5 * dt * interact->fyh;
    }
}
