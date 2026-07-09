//! Stage 3 + 5: the free-surface POUR validator and the deliverable
//! render. A tank with a weir lip pours under a rotating-gravity TILT
//! SCHEDULE (the moving frame fs-lbm's G3 gate pins); the outcome
//! record carries the STRICT mass-ledger drift, the Plateau–Rayleigh
//! fragment count of the free jet, and the dribble indicator (wetted
//! outer-lip cells after the pour) — computed under BOTH contact
//! models so the bracketing band is a first-class output, per the
//! plan's honest handling of contact-line fidelity. The render binds
//! fs-render's Woodcock tracker to the simulation's OWN mass buffer:
//! the marketing shot and the physics are the same bytes.

use fs_lbm::freesurface::{ContactModel, FreeSurface};
use fs_lbm::rheology::{Rheology, update_tau};
use fs_lbm::{Cell, Grid};
use fs_render::volumes::{MajorantGrid, VolumeGrid, render_transmittance};

/// The pour fixture: a tank with an interior weir (the lip) the fluid
/// crosses as gravity tilts.
pub struct PourRig {
    /// Lattice width.
    pub nx: usize,
    /// Lattice height.
    pub ny: usize,
    /// Weir column (the lip sits at this x).
    pub lip_x: usize,
    /// Weir height (cells).
    pub lip_h: usize,
    /// Gravity magnitude.
    pub g0: f64,
    /// Tilt schedule: final angle (radians) reached linearly over the
    /// run (rotating gravity — the fs-scenario moving-frame stand-in).
    pub tilt_final: f64,
    /// Steps.
    pub steps: u32,
}

impl Default for PourRig {
    fn default() -> Self {
        PourRig {
            nx: 48,
            ny: 28,
            lip_x: 26,
            lip_h: 10,
            // Gravity-wave timescale must beat the run length: at
            // 1.2e-4 the slosh barely crossed the chamber in 700 steps
            // and NOTHING poured (measured); 6e-4 pours decisively.
            g0: 6e-4,
            tilt_final: 0.7,
            steps: 900,
        }
    }
}

/// The validator's outcome record — the flagship's evidence row.
#[derive(Debug, Clone)]
pub struct PourOutcome {
    /// Worst relative mass-ledger drift over the run.
    pub mass_drift: f64,
    /// Free-surface fragments at the end (Plateau–Rayleigh score).
    pub fragments: usize,
    /// Fluid mass that crossed the lip (the pour actually poured).
    pub poured_mass: f64,
    /// Dribble indicator: wet cells clinging to the outer lip wall
    /// after the pour (the contact-line proxy).
    pub dribble_cells: usize,
    /// The final mass field (nx × ny) — the SAME bytes the render
    /// binds (kept for the deliverable path).
    pub mass_field: Vec<f64>,
}

/// Run one pour under a contact model and a Carreau fluid.
///
/// # Panics
/// Only on fs-lbm programmer contracts (fixture-scale).
#[must_use]
pub fn run_pour(rig: &PourRig, contact: ContactModel, law: Rheology) -> PourOutcome {
    let (nx, ny) = (rig.nx, rig.ny);
    let mut grid = Grid::uniform(nx, ny, 0.55);
    grid.periodic_x = false;
    grid.periodic_y = false;
    grid.g = [0.0, -rig.g0];
    for i in 0..nx * ny {
        grid.flags[i] = Cell::Gas;
    }
    // Walls: floor, ceiling, both sides, and the weir column up to
    // lip_h.
    for x in 0..nx {
        let b = grid.idx(x, 0);
        grid.flags[b] = Cell::Wall;
        let t = grid.idx(x, ny - 1);
        grid.flags[t] = Cell::Wall;
    }
    for y in 0..ny {
        let l = grid.idx(0, y);
        grid.flags[l] = Cell::Wall;
        let r = grid.idx(nx - 1, y);
        grid.flags[r] = Cell::Wall;
    }
    for y in 1..=rig.lip_h.min(ny - 2) {
        let w = grid.idx(rig.lip_x, y);
        grid.flags[w] = Cell::Wall;
    }
    // Fluid: fill the left chamber to just under the lip.
    for y in 1..rig.lip_h.min(ny - 2) {
        for x in 1..rig.lip_x {
            let i = grid.idx(x, y);
            grid.flags[i] = Cell::Fluid;
        }
    }
    let mut sim = FreeSurface::new(grid, 0.0, contact);
    let m0 = sim.ledger_mass();
    let mut mass_drift = 0.0f64;
    for step in 0..rig.steps {
        // Tilt schedule: rotate gravity linearly to tilt_final.
        let theta = rig.tilt_final * f64::from(step) / f64::from(rig.steps);
        sim.grid.g = [rig.g0 * theta.sin(), -rig.g0 * theta.cos()];
        // Carreau/power-law local viscosity adaptation (stage-3 band).
        let _ = update_tau(&mut sim.grid, law);
        sim.step();
        mass_drift = mass_drift.max(((sim.ledger_mass() - m0) / m0).abs());
    }
    // Poured mass: fluid/interface mass right of the weir.
    let mut poured = 0.0f64;
    let mut mass_field = vec![0.0f64; nx * ny];
    for y in 0..ny {
        for x in 0..nx {
            let i = sim.grid.idx(x, y);
            let m = match sim.grid.flags[i] {
                Cell::Fluid => sim.grid.f[i].iter().sum::<f64>(),
                Cell::Interface => sim.mass[i],
                _ => 0.0,
            };
            mass_field[i] = m;
            if x > rig.lip_x {
                poured += m;
            }
        }
    }
    // Dribble: wet cells hugging the outer face of the weir below the
    // lip crest (fluid that clung and crept down instead of jetting).
    let mut dribble = 0usize;
    for y in 1..rig.lip_h.min(ny - 2) {
        let i = sim.grid.idx(rig.lip_x + 1, y);
        if matches!(sim.grid.flags[i], Cell::Fluid | Cell::Interface) {
            dribble += 1;
        }
    }
    PourOutcome {
        mass_drift,
        fragments: sim.fragment_count(),
        poured_mass: poured,
        dribble_cells: dribble,
        mass_field,
    }
}

/// The deliverable: render the pour's mass field with the Woodcock
/// tracker, bound ZERO-COPY to the outcome's own buffer. Returns the
/// image (res × res transmittance).
#[must_use]
pub fn render_pour(outcome: &PourOutcome, nx: usize, ny: usize, res: usize) -> Vec<f64> {
    let grid = VolumeGrid::new(
        [nx, ny, 1],
        &outcome.mass_field,
        [0.0, 0.0, 0.0],
        [1.0, 1.0, f64::from(u32::try_from(ny).expect("small"))],
    );
    let majorant = MajorantGrid::build(&grid, 8);
    render_transmittance(&grid, &majorant, res, 24, 0x7E55E1)
}
