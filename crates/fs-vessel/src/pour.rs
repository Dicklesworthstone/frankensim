//! Stage 3 + 5: the free-surface POUR validator and the deliverable
//! render. A tank with a weir lip pours in a prescribed rotating frame,
//! including Euler, centrifugal and Coriolis acceleration; the outcome
//! record carries the STRICT mass-ledger drift, the connected-component
//! count of the 2-D liquid, and the dribble indicator (wetted
//! outer-lip cells after the pour) — computed under BOTH contact
//! models so the bracketing band is a first-class output, per the
//! plan's honest handling of contact-line fidelity. The render binds
//! fs-render's exact slab absorption to the simulation's OWN mass buffer:
//! the marketing shot and the physics are the same bytes.

use fs_lbm::core2::frame::{FrameForcingError2, RigidFrame2};
use fs_lbm::freesurface::{ContactModel, FreeSurface};
use fs_lbm::rheology::{Rheology, TAU_CAP, TAU_FLOOR, update_tau};
use fs_lbm::{Cell, Grid};
use fs_render::volumes::{AbsorptionImageSettings, VolumeGrid, render_absorption_transmittance};

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

impl PourRig {
    /// Prepare one physical pour step on an existing native liquid grid.
    ///
    /// Native and browser consumers share this rest-to-rest tilt schedule and
    /// its midpoint gravity, centrifugal, Euler and Coriolis forcing. `step`
    /// is zero-based and `pivot` uses native lattice cell-centre coordinates.
    /// Preparation preserves the accepted physical momentum and stress; the
    /// caller then advances the original flow owner exactly once.
    ///
    /// Refuses an exhausted/empty schedule or invalid motion/grid before any
    /// accepted population or force changes. Motion is prescribed; this does
    /// not evolve finite-inertia vessel/support dynamics.
    pub fn prepare_step(
        &self,
        grid: &mut Grid,
        step: u32,
        pivot: [f64; 2],
    ) -> Result<(), FrameForcingError2> {
        if step >= self.steps {
            return Err(FrameForcingError2::InvalidMotion);
        }
        let [theta, speed, acceleration] = tilt_kinematics(self, f64::from(step) + 0.5);
        RigidFrame2 {
            origin: pivot,
            angular_velocity: -speed,
            angular_acceleration: -acceleration,
            origin_acceleration: [0.0; 2],
        }
        .prepare_step(
            grid,
            [
                self.g0 * fs_math::det::sin(theta),
                -self.g0 * fs_math::det::cos(theta),
            ],
        )
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

/// A pour that cannot report a trajectory for the supplied fluid and motion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PourError {
    /// A zero-step run supplies no evolved pour to validate.
    EmptySchedule,
    /// The original frame owner refused before advancing this step.
    FrameForcing {
        /// Zero-based step whose forcing could not be prepared.
        step: u32,
        /// Original native refusal, without a substitute motion.
        source: FrameForcingError2,
    },
    /// The native rheology update would substitute a clipped viscosity.
    RheologyClipped {
        /// Zero-based step refused before collision, streaming or VOF transport.
        step: u32,
        /// Fluid/interface cells below the representable relaxation floor.
        floored_cells: u32,
        /// Fluid/interface cells above the representable relaxation ceiling.
        capped_cells: u32,
    },
}

impl core::fmt::Display for PourError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptySchedule => f.write_str("a pour requires at least one time step"),
            Self::FrameForcing { step, source } => {
                write!(f, "pour step {step}: {source}")
            }
            Self::RheologyClipped {
                step,
                floored_cells,
                capped_cells,
            } => write!(
                f,
                "pour step {step}: {floored_cells} cells below tau={TAU_FLOOR} and \
                 {capped_cells} cells above tau={TAU_CAP}; rescale the lattice to \
                 represent the supplied fluid instead of clipping its viscosity"
            ),
        }
    }
}

impl std::error::Error for PourError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::FrameForcing { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// Advance the original liquid owner only if its local law is representable.
///
/// Native and browser pours share this check at every step, because a
/// shear-dependent law can leave the lattice range after starting inside it.
/// Refusal leaves physical populations and VOF transport unadvanced. The
/// relaxation buffer holds the rejected proposal; callers must discard that
/// trial rather than continuing with its clipped viscosity. There is no added
/// population pass, scratch allocation or change to admitted-step arithmetic.
///
/// # Errors
/// Returns the step and native floor/cap counts whenever any live cell would
/// use a viscosity different from the supplied constitutive law.
///
/// # Panics
/// Invalid constitutive parameters or the original native flow contracts.
pub fn advance_pour_liquid(
    sim: &mut FreeSurface,
    law: Rheology,
    step: u32,
) -> Result<(), PourError> {
    let stats = update_tau(&mut sim.grid, law);
    if stats.floored != 0 || stats.capped != 0 {
        return Err(PourError::RheologyClipped {
            step,
            floored_cells: stats.floored,
            capped_cells: stats.capped,
        });
    }
    sim.step();
    Ok(())
}

/// Run one pour under a contact model and a Carreau fluid.
///
/// # Errors
/// Refuses empty motion, invalid frame forcing or any clipped local viscosity.
/// A refused run returns no mass field or successful validation outcome.
///
/// # Panics
/// Only on fs-lbm programmer contracts (fixture-scale).
pub fn run_pour(
    rig: &PourRig,
    contact: ContactModel,
    law: Rheology,
) -> Result<PourOutcome, PourError> {
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
/// # Errors
/// Empty motion, nonfinite pivot/motion, or clipped local viscosity. No
/// successful pour outcome is returned for a numerically substituted fluid.
///
/// # Panics
/// The original native flow admission contracts.
pub fn run_pour_in_frame(
    rig: &PourRig,
    contact: ContactModel,
    law: Rheology,
    pivot: [f64; 2],
) -> Result<PourOutcome, PourError> {
    run_pour_motion(rig, contact, law, Some(pivot))
}

/// Historical rotating-gravity smoke model, retained for explicit comparisons.
///
/// The linear gravity-angle ramp omits frame inertia and is not the momentum
/// equation of a physically rotating vessel. Existing cross-ISA smoke goldens
/// can select this model explicitly; physical pours use [`run_pour`].
///
/// This comparison still refuses clipped rheology: retaining its historical
/// motion does not permit silently substituting a different fluid.
///
/// # Errors
/// An empty schedule or any locally clipped viscosity.
///
/// # Panics
/// The original native flow admission contracts.
pub fn run_pour_rotating_gravity(
    rig: &PourRig,
    contact: ContactModel,
    law: Rheology,
) -> Result<PourOutcome, PourError> {
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
) -> Result<PourOutcome, PourError> {
    if rig.steps == 0 {
        return Err(PourError::EmptySchedule);
    }
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
            rig.prepare_step(&mut sim.grid, step, origin)
                .map_err(|source| PourError::FrameForcing { step, source })?;
        } else {
            // Preserve the explicit historical model's deterministic arithmetic.
            let theta = rig.tilt_final * f64::from(step) / f64::from(rig.steps);
            sim.grid.g = [
                rig.g0 * fs_math::det::sin(theta),
                -rig.g0 * fs_math::det::cos(theta),
            ];
        }
        advance_pour_liquid(&mut sim, law, step)?;
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
    Ok(PourOutcome {
        mass_drift,
        fragments: sim.fragment_count(),
        poured_mass: poured,
        dribble_cells: dribble,
        mass_field,
    })
}

#[cfg(test)]
mod capillary_tests {
    use super::*;

    /// G0: the complete validator rejects both directions of constitutive
    /// substitution, including the explicit historical motion comparison.
    #[test]
    fn pours_refuse_viscosities_outside_the_native_lattice_range() {
        let rig = PourRig {
            nx: 12,
            ny: 8,
            lip_x: 6,
            lip_h: 4,
            g0: 0.0,
            surface_tension: 0.0,
            tilt_final: 0.0,
            steps: 4,
        };
        // Exactly 5 columns by 3 rows are wet initially. Interface insertion
        // changes their flags but not the count requiring a constitutive law.
        for (nu, floored_cells, capped_cells) in [(1e-5, 15, 0), (2.0, 0, 15)] {
            let law = Rheology::Newtonian { nu };
            let expected = PourError::RheologyClipped {
                step: 0,
                floored_cells,
                capped_cells,
            };
            assert_eq!(
                run_pour(&rig, ContactModel::Neutral, law).unwrap_err(),
                expected
            );
            assert_eq!(
                run_pour_in_frame(&rig, ContactModel::Wetting, law, [0.5, 0.5]).unwrap_err(),
                expected,
            );
            assert_eq!(
                run_pour_rotating_gravity(&rig, ContactModel::Neutral, law).unwrap_err(),
                expected,
            );
            assert!(expected.to_string().contains("rescale the lattice"));
        }
    }

    /// G0/G3: zero-shear admission cannot authorize a later clipped response.
    #[test]
    fn transient_shear_thinning_refuses_when_the_actual_flow_reaches_the_floor() {
        let rig = PourRig {
            steps: 300,
            surface_tension: 0.0,
            ..PourRig::default()
        };
        let law = Rheology::Carreau {
            nu0: 0.0167,
            nu_inf: 1e-6,
            lambda: 10_000.0,
            n: 0.1,
        };
        assert!(0.5 + 3.0 * law.viscosity(0.0) > TAU_FLOOR);
        let error = run_pour(&rig, ContactModel::Neutral, law).unwrap_err();
        let PourError::RheologyClipped {
            step,
            floored_cells,
            capped_cells,
        } = error
        else {
            panic!("expected constitutive refusal, got {error:?}");
        };
        assert!(step > 0 && step < rig.steps, "step {step}");
        assert!(floored_cells > 0);
        assert_eq!(capped_cells, 0);
    }

    /// G0: a rejected viscosity proposal never drives collision or mass flux.
    #[test]
    fn rheology_refusal_precedes_population_and_free_surface_transport() {
        let mut grid = Grid::uniform(6, 5, 0.8);
        for y in 3..5 {
            for x in 0..6 {
                let cell = grid.idx(x, y);
                grid.flags[cell] = Cell::Gas;
            }
        }
        let mut sim = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
        let before = sim.clone();
        assert_eq!(
            advance_pour_liquid(&mut sim, Rheology::Newtonian { nu: 2.0 }, 7),
            Err(PourError::RheologyClipped {
                step: 7,
                floored_cells: 0,
                capped_cells: 18,
            }),
        );
        assert_eq!(sim.grid.f, before.grid.f);
        assert_eq!(sim.grid.flags, before.grid.flags);
        assert_eq!(sim.mass, before.mass);
        assert_eq!(sim.ledger_mass(), before.ledger_mass());
        assert_eq!(sim.conversions.to_fluid, before.conversions.to_fluid);
        assert_eq!(sim.conversions.to_gas, before.conversions.to_gas);
        assert_eq!(
            sim.conversions.gas_to_interface,
            before.conversions.gas_to_interface
        );
        assert_eq!(
            sim.conversions.fluid_to_interface,
            before.conversions.fluid_to_interface
        );
    }

    #[test]
    fn empty_or_invalid_motion_does_not_return_a_successful_pour() {
        let empty = PourRig {
            steps: 0,
            ..PourRig::default()
        };
        let law = Rheology::Newtonian { nu: 0.0167 };
        assert_eq!(
            run_pour(&empty, ContactModel::Neutral, law).unwrap_err(),
            PourError::EmptySchedule,
        );
        assert_eq!(
            run_pour_in_frame(
                &PourRig::default(),
                ContactModel::Neutral,
                law,
                [f64::NAN, 0.0]
            )
            .unwrap_err(),
            PourError::FrameForcing {
                step: 0,
                source: FrameForcingError2::InvalidMotion
            },
        );
    }

    #[test]
    fn pour_image_retains_thin_absorption_without_collision_noise() {
        let outcome = PourOutcome {
            mass_drift: 0.0,
            fragments: 1,
            poured_mass: 0.0,
            dribble_cells: 0,
            mass_field: vec![0.0, 1e-9, 0.5, 1.0],
        };
        let image = render_pour(&outcome, 2, 2, 2);
        for (actual, mass) in image.iter().zip(&outcome.mass_field) {
            assert!((actual - (-mass * 2.0).exp()).abs() < 2e-15);
        }
        assert!(image[1] < 1.0);
        assert_eq!(image, render_pour(&outcome, 2, 2, 2));
    }

    #[test]
    fn shared_pour_preparation_preserves_momentum_then_delivers_native_impulse() {
        let rig = PourRig {
            g0: 0.0,
            steps: 300,
            ..PourRig::default()
        };
        let mut grid = Grid::uniform(8, 6, 0.8);
        grid.f.fill(fs_lbm::equilibrium(0.73, 0.01, -0.02));
        let initial = grid.moments(0);
        rig.prepare_step(&mut grid, 100, [0.5, 0.5]).unwrap();
        assert!(
            grid.fext
                .iter()
                .any(|force| force.iter().any(|v| v.abs() > 1e-6))
        );
        let mut collided = grid.clone();
        grid.collide_into(&mut collided.f);
        for i in 0..grid.f.len() {
            let prepared = grid.moments(i);
            let after = collided.moments(i);
            for axis in 0..2 {
                let original = initial.rho * initial.u[axis];
                let held = prepared.rho * prepared.u[axis];
                let delivered = after.rho * after.u[axis] - held;
                assert!((held - original).abs() < 3e-15);
                assert!((delivered - grid.fext[i][axis]).abs() < 3e-15);
            }
        }
    }

    #[test]
    fn refused_pour_schedule_leaves_native_state_untouched() {
        let mut rig = PourRig::default();
        let mut grid = Grid::uniform(4, 4, 0.8);
        let before = grid.clone();
        assert_eq!(
            rig.prepare_step(&mut grid, rig.steps, [0.5, 0.5]),
            Err(FrameForcingError2::InvalidMotion)
        );
        rig.steps = 0;
        assert_eq!(
            rig.prepare_step(&mut grid, 0, [0.5, 0.5]),
            Err(FrameForcingError2::InvalidMotion)
        );
        assert_eq!(grid.f, before.f);
        assert_eq!(grid.fext, before.fext);
        assert_eq!(grid.g, before.g);
    }

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
        let left = run_pour_in_frame(&rig, ContactModel::Neutral, law, [0.5, 0.5]).unwrap();
        let right = run_pour_in_frame(&rig, ContactModel::Neutral, law, [46.5, 0.5]).unwrap();
        let gravity_only = run_pour_rotating_gravity(&rig, ContactModel::Neutral, law).unwrap();
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
        let neutral = run_pour(&rig, ContactModel::Neutral, liquid).unwrap();
        let wetting = run_pour(&rig, ContactModel::Wetting, liquid).unwrap();
        let no_tension = run_pour(
            &PourRig {
                surface_tension: 0.0,
                ..rig
            },
            ContactModel::Neutral,
            liquid,
        )
        .unwrap();
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

/// Render exact slab transmittance, borrowing the outcome's mass field.
///
/// The fixture maps lattice mass to absorption per lattice length and uses
/// `ny` as the displayed slab thickness. This is an absorption visualization,
/// not a calibrated liquid optical material or refracting 3-D free surface.
///
/// # Panics
/// Invalid dimensions, field values or an unrepresentable image allocation.
#[must_use]
pub fn render_pour(outcome: &PourOutcome, nx: usize, ny: usize, res: usize) -> Vec<f64> {
    let grid = VolumeGrid::new(
        [nx, ny, 1],
        &outcome.mass_field,
        [0.0, 0.0, 0.0],
        [1.0, 1.0, f64::from(u32::try_from(ny).expect("small"))],
    );
    let pixels = res.checked_mul(res).expect("render resolution fits usize");
    render_absorption_transmittance(
        &grid,
        AbsorptionImageSettings {
            resolution: [res; 2],
            max_pixels: pixels,
            max_grid_cells: outcome.mass_field.len(),
            max_cell_visits: u64::try_from(pixels).expect("render work fits u64"),
        },
    )
    .expect("valid finite pour absorption image")
}
