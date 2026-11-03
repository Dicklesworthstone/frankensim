//! Stage 3 + 5: the free-surface POUR validator and the deliverable
//! render. A tank with a weir lip pours in a prescribed rotating frame,
//! including Euler, centrifugal and Coriolis acceleration; the outcome
//! record carries the STRICT mass-ledger drift, the connected-component
//! count of the 2-D liquid, and the dribble indicator (wetted
//! outer-lip cells after the pour) — computed under BOTH contact
//! models so the bracketing band is a first-class output, per the
//! plan's honest handling of contact-line fidelity. The render binds
//! fs-render's Woodcock tracker to the simulation's OWN mass buffer:
//! the marketing shot and the physics are the same bytes.

use fs_lbm::core2::frame::RigidFrame2;
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
    /// Liquid/gas surface tension in lattice units; zero explicitly omits it.
    /// Convert SI values with `sigma_lu = sigma_si * dt^2 / (rho_si * dx^3)`.
    /// This fixture coefficient is not a sourced fluid-property calibration.
    pub surface_tension: f64,
    /// Final clockwise tilt (radians). The physical run uses a quintic
    /// rest-to-rest schedule with zero endpoint angular acceleration.
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
            surface_tension: 0.002,
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
    /// Connected liquid components in this 2-D fixture. This is not a
    /// quantitative Plateau–Rayleigh breakup model for a cylindrical jet.
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
    run_pour_in_frame(rig, contact, law, [0.5 * (rig.nx as f64 - 1.0), 0.5])
}

/// Pour in a prescribed rigid frame rotating about an explicit lattice pivot.
///
/// Cell centres have integer coordinates. The default [`run_pour`] pivot is
/// the middle of the halfway-bounce-back floor. Gravity is transformed into
/// this frame; native D2Q9 forcing includes centrifugal, Euler and Coriolis
/// acceleration. Tilt and its derivatives come from the same quintic schedule,
/// sampled at step midpoints, starting and finishing with zero spin and angular
/// acceleration. The stationary lattice walls are the co-moving vessel walls.
///
/// The prescribed support supplies the vessel motion; no coupled finite-inertia
/// vessel, gas motion, 3-D pour or quantitative jet-breakup model is claimed.
///
/// # Panics
/// Nonfinite pivot/motion or the original native flow admission contracts.
#[must_use]
pub fn run_pour_in_frame(
    rig: &PourRig,
    contact: ContactModel,
    law: Rheology,
    pivot: [f64; 2],
) -> PourOutcome {
    assert!(
        pivot.into_iter().all(f64::is_finite),
        "finite pour pivot required"
    );
    run_pour_motion(rig, contact, law, Some(pivot))
}

/// Historical rotating-gravity smoke model, retained for explicit comparisons.
///
/// The linear gravity-angle ramp omits frame inertia and is not the momentum
/// equation of a physically rotating vessel. Existing cross-ISA smoke goldens
/// can select this model explicitly; physical pours use [`run_pour`].
///
/// # Panics
/// The original native flow admission contracts.
#[must_use]
pub fn run_pour_rotating_gravity(
    rig: &PourRig,
    contact: ContactModel,
    law: Rheology,
) -> PourOutcome {
    run_pour_motion(rig, contact, law, None)
}

fn tilt_kinematics(rig: &PourRig, time: f64) -> [f64; 3] {
    let duration = f64::from(rig.steps);
    let s = time / duration;
    let shape = s * s * s * (10.0 + s * (-15.0 + 6.0 * s));
    let speed = 30.0 * s * s * (1.0 - s) * (1.0 - s);
    let acceleration = 60.0 * s * (1.0 - s) * (1.0 - 2.0 * s);
    [
        rig.tilt_final * shape,
        rig.tilt_final * speed / duration,
        rig.tilt_final * acceleration / duration / duration,
    ]
}

fn run_pour_motion(
    rig: &PourRig,
    contact: ContactModel,
    law: Rheology,
    pivot: Option<[f64; 2]>,
) -> PourOutcome {
    let (nx, ny) = (rig.nx, rig.ny);
    let mut grid = Grid::uniform(nx, ny, 0.55);
    grid.periodic_x = false;
    grid.periodic_y = false;
    // Start the physical run at rest, before applying any half-source frame.
    // The historical comparison retains its original initialization exactly.
    if pivot.is_none() {
        grid.g = [0.0, -rig.g0];
    }
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
    let mut sim = FreeSurface::new(grid, rig.surface_tension, contact);
    let m0 = sim.ledger_mass();
    let mut mass_drift = 0.0f64;
    for step in 0..rig.steps {
        if let Some(origin) = pivot {
            let [theta, speed, acceleration] = tilt_kinematics(rig, f64::from(step) + 0.5);
            // Clockwise tilt is negative counterclockwise frame rotation.
            RigidFrame2 {
                origin,
                angular_velocity: -speed,
                angular_acceleration: -acceleration,
                origin_acceleration: [0.0; 2],
            }
            .prepare_step(
                &mut sim.grid,
                [
                    rig.g0 * fs_math::det::sin(theta),
                    -rig.g0 * fs_math::det::cos(theta),
                ],
            )
            .expect("pour rigid-frame forcing");
        } else {
            // Preserve the explicit historical model's deterministic arithmetic.
            let theta = rig.tilt_final * f64::from(step) / f64::from(rig.steps);
            sim.grid.g = [
                rig.g0 * fs_math::det::sin(theta),
                -rig.g0 * fs_math::det::cos(theta),
            ];
        }
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

#[cfg(test)]
mod capillary_tests {
    use super::*;

    #[test]
    fn rigid_tilt_schedule_starts_and_finishes_at_rest() {
        let rig = PourRig::default();
        assert_eq!(tilt_kinematics(&rig, 0.0), [0.0; 3]);
        assert_eq!(
            tilt_kinematics(&rig, f64::from(rig.steps)),
            [rig.tilt_final, 0.0, -0.0]
        );
        let middle = tilt_kinematics(&rig, 0.5 * f64::from(rig.steps));
        assert!((middle[0] - 0.5 * rig.tilt_final).abs() < 1e-15);
        assert!(middle[1] > 0.0);
        assert_eq!(middle[2], 0.0);
    }

    #[test]
    fn original_flow_responds_to_rigid_tilt_pivot_even_without_gravity() {
        let rig = PourRig {
            g0: 0.0,
            steps: 300,
            tilt_final: 0.5,
            surface_tension: 0.0,
            ..PourRig::default()
        };
        let law = Rheology::Newtonian { nu: 0.0167 };
        let left = run_pour_in_frame(&rig, ContactModel::Neutral, law, [0.5, 0.5]);
        let right = run_pour_in_frame(&rig, ContactModel::Neutral, law, [46.5, 0.5]);
        let gravity_only = run_pour_rotating_gravity(&rig, ContactModel::Neutral, law);
        for result in [&left, &right, &gravity_only] {
            assert!(result.mass_drift < 1e-10, "{}", result.mass_drift);
            assert!(result.mass_field.iter().all(|v| v.is_finite()));
        }
        let difference = |a: &PourOutcome, b: &PourOutcome| {
            a.mass_field
                .iter()
                .zip(&b.mass_field)
                .map(|(x, y)| (x - y).abs())
                .sum::<f64>()
        };
        assert!(difference(&left, &right) > 1e-3);
        assert!(difference(&left, &gravity_only) > 1e-3);
        // This response cannot come from rotating gravity: it is identically
        // zero in all three runs. The only driving mechanism is frame inertia.
    }

    #[test]
    fn capillary_pour_responds_to_tension_and_wetting_without_losing_mass() {
        let rig = PourRig::default();
        assert!(rig.surface_tension > 0.0);
        let liquid = Rheology::Newtonian { nu: 0.0167 };
        let neutral = run_pour(&rig, ContactModel::Neutral, liquid);
        let wetting = run_pour(&rig, ContactModel::Wetting, liquid);
        let no_tension = run_pour(
            &PourRig {
                surface_tension: 0.0,
                ..rig
            },
            ContactModel::Neutral,
            liquid,
        );
        for result in [&neutral, &wetting, &no_tension] {
            assert!(
                result.mass_drift < 1e-10,
                "mass drift {}",
                result.mass_drift
            );
            assert!(
                result.poured_mass > 1.0,
                "poured mass {}",
                result.poured_mass
            );
            assert!(result.mass_field.iter().all(|mass| mass.is_finite()));
        }
        let difference = |a: &PourOutcome, b: &PourOutcome| {
            a.mass_field
                .iter()
                .zip(&b.mass_field)
                .map(|(x, y)| (x - y).abs())
                .sum::<f64>()
        };
        assert!(difference(&neutral, &wetting) > 1e-8);
        assert!(difference(&neutral, &no_tension) > 1e-8);
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
