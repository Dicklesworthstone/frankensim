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

/// Original sourced liquid and SI scaling on this same native pour owner.
pub mod material;

use fs_lbm::core2::frame::{FrameForcingError2, RigidFrame2};
use fs_lbm::freesurface::{ContactModel, FreeSurface};
use fs_lbm::rheology::{Rheology, TAU_CAP, TAU_FLOOR, update_tau};
use fs_lbm::{Cell, Grid};
use fs_render::volumes::{
    AbsorptionImageSettings, DensityAbsorption, DensityAbsorptionSpectrum, DvrError,
    SpectralAbsorptionImageSettings, SpectralOpacityImage, VolumeGrid,
    render_absorption_transmittance, render_density_spectral_box,
    render_density_spectral_transmittance_box, render_density_transmittance_box,
};

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

/// A monochromatic SI absorption view of an accepted two-dimensional pour.
///
/// The original field contains cell-averaged liquid mass per lattice cell:
/// full-cell lattice density in liquid cells, volume-fraction-weighted density
/// in interface cells, and zero in gas/wall cells. Multiplying it by
/// `density_unit_kg_m3` supplies the liquid mass density for the optical slab.
///
/// The view explicitly extrudes this two-dimensional field through
/// `slab_depth_m`. Interface cells use a homogeneous density along that ray;
/// the mass map does not determine a sharp interface or subcell coverage.
/// This computes liquid absorption only. It does not infer wall occlusion,
/// surface reflection/refraction, scattering, thermal emission or a spectrum.
/// The caller supplies the coefficient at the stated wavelength; changing only
/// the wavelength label does not invent a new material coefficient.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PourOpticalView {
    /// Positive vacuum wavelength at which the coefficient was supplied [nm].
    pub wavelength_nm: f64,
    /// Intensity mass absorption at `wavelength_nm` [m²/kg], nonnegative.
    pub mass_absorption_m2_kg: f64,
    /// Physical density represented by one unit of cell-average lattice mass
    /// [kg/m³]. This is a density conversion, not kilograms per grid cell.
    pub density_unit_kg_m3: f64,
    /// Positive physical spacing of the square in-plane lattice cells [m].
    pub lattice_spacing_m: f64,
    /// Positive out-of-plane path length [m], independent of grid resolution.
    pub slab_depth_m: f64,
}

/// Wavelength-resolved SI absorption of the same accepted pour mass field.
///
/// The caller supplies one to sixteen positive, strictly increasing vacuum
/// wavelengths and their nonnegative mass absorption coefficients at the
/// material state being observed. These are monochromatic samples: no spectral
/// interpolation, wavelength-bin widths, illumination or detector response is
/// inferred. The spectrum and the original mass field are borrowed.
///
/// As for [`PourOpticalView`], the two-dimensional mass field is extruded
/// through an explicit slab depth. Each wavelength's transmitted intensity is
/// averaged over the exact pixel footprint, preserving partial liquid coverage
/// independently of its absorption strength. This includes bulk absorption,
/// not surface reflection/refraction, wall occlusion, scattering or emission.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PourSpectralOpticalView<'a> {
    /// Positive strictly increasing vacuum wavelength samples [nm].
    pub wavelengths_nm: &'a [f64],
    /// Intensity mass absorption at each corresponding wavelength [m²/kg].
    pub mass_absorption_m2_kg: &'a [f64],
    /// Physical density represented by one unit of cell-average mass [kg/m³].
    pub density_unit_kg_m3: f64,
    /// Positive physical spacing of the square in-plane lattice cells [m].
    pub lattice_spacing_m: f64,
    /// Positive out-of-plane path length [m], independent of grid resolution.
    pub slab_depth_m: f64,
}

/// A physical pour view that cannot produce an admitted transmittance image.
#[derive(Debug)]
pub enum PourOpticalError {
    /// The declared vacuum wavelength is nonfinite or nonpositive.
    InvalidWavelength,
    /// The physical cell spacing or extruded path length is invalid.
    InvalidGeometry,
    /// The original bounded renderer refused the field, conversion or budget.
    Transport(DvrError),
}

impl core::fmt::Display for PourOpticalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidWavelength => {
                f.write_str("pour absorption requires a positive finite vacuum wavelength")
            }
            Self::InvalidGeometry => {
                f.write_str("pour absorption requires positive finite SI spacing and slab depth")
            }
            Self::Transport(error) => write!(f, "pour absorption transport refused: {error:?}"),
        }
    }
}

impl std::error::Error for PourOpticalError {}

impl PourOpticalView {
    /// Render a borrowed accepted mass field with explicit work/memory limits.
    ///
    /// Each original column transmits
    /// `exp(-mass * density_unit_kg_m3 * mass_absorption_m2_kg * slab_depth_m)`.
    /// Each output pixel averages that transmission over its complete rectangular
    /// footprint. Narrow liquid features retain their actual projected coverage;
    /// density or optical depth is never averaged before exponentiation.
    /// The original mass buffer is borrowed without allocating a converted
    /// density/extinction field or a supersampling image.
    ///
    /// Invalid geometry, samples, conversion or budget returns no partial image.
    pub fn render_mass_field(
        self,
        mass_field: &[f64],
        dimensions: [usize; 2],
        settings: AbsorptionImageSettings,
    ) -> Result<Vec<f64>, PourOpticalError> {
        if !self.wavelength_nm.is_finite() || self.wavelength_nm <= 0.0 {
            return Err(PourOpticalError::InvalidWavelength);
        }
        if !self.lattice_spacing_m.is_finite()
            || self.lattice_spacing_m <= 0.0
            || !self.slab_depth_m.is_finite()
            || self.slab_depth_m <= 0.0
        {
            return Err(PourOpticalError::InvalidGeometry);
        }
        // Struct construction lets the shared validator reject dimension
        // overflow/mismatched storage without VolumeGrid::new's panic.
        let grid = VolumeGrid {
            dims: [dimensions[0], dimensions[1], 1],
            data: mass_field,
            origin: [0.0; 3],
            cell: [
                self.lattice_spacing_m,
                self.lattice_spacing_m,
                self.slab_depth_m,
            ],
        };
        render_density_transmittance_box(
            &grid,
            settings,
            DensityAbsorption {
                density_unit_kg_m3: self.density_unit_kg_m3,
                length_unit_m: 1.0,
                mass_absorption_m2_kg: self.mass_absorption_m2_kg,
            },
        )
        .map_err(PourOpticalError::Transport)
    }
}

impl<'a> PourSpectralOpticalView<'a> {
    /// Render per-wavelength transmission, opacity and effective optical depth.
    ///
    /// The image retains the supplied wavelengths. Samples are pixel-major,
    /// with wavelengths varying fastest: `(y * width + x) * bands + band`.
    /// Geometry and density traversal are shared across wavelengths; no
    /// converted density field or per-wavelength simulation copy is allocated.
    /// Invalid input or exhausted image/spectral work budgets returns no image.
    pub fn render_mass_field(
        self,
        mass_field: &[f64],
        dimensions: [usize; 2],
        settings: SpectralAbsorptionImageSettings,
    ) -> Result<SpectralOpacityImage, PourOpticalError> {
        let grid = self.grid(mass_field, dimensions)?;
        render_density_spectral_box(&grid, settings, self.absorption())
            .map_err(PourOpticalError::Transport)
    }

    /// Render only direct intensity transmittance in pixel-major spectral order.
    ///
    /// This writes the final `pixels * bands` values without allocating the
    /// full opacity/depth image first. Each wavelength is filtered after
    /// exponentiation; neither density nor optical depth is spatially averaged.
    pub fn render_transmittance(
        self,
        mass_field: &[f64],
        dimensions: [usize; 2],
        settings: SpectralAbsorptionImageSettings,
    ) -> Result<Vec<f64>, PourOpticalError> {
        let grid = self.grid(mass_field, dimensions)?;
        render_density_spectral_transmittance_box(&grid, settings, self.absorption())
            .map_err(PourOpticalError::Transport)
    }

    fn absorption(self) -> DensityAbsorptionSpectrum<'a> {
        DensityAbsorptionSpectrum {
            wavelengths_nm: self.wavelengths_nm,
            mass_absorption_m2_kg: self.mass_absorption_m2_kg,
            density_unit_kg_m3: self.density_unit_kg_m3,
            length_unit_m: 1.0,
        }
    }

    fn grid<'m>(
        self,
        mass_field: &'m [f64],
        dimensions: [usize; 2],
    ) -> Result<VolumeGrid<'m>, PourOpticalError> {
        if !self.lattice_spacing_m.is_finite()
            || self.lattice_spacing_m <= 0.0
            || !self.slab_depth_m.is_finite()
            || self.slab_depth_m <= 0.0
        {
            return Err(PourOpticalError::InvalidGeometry);
        }
        // Shared transport admission checks the original shape and samples
        // without VolumeGrid::new's panic on invalid dimensions.
        Ok(VolumeGrid {
            dims: [dimensions[0], dimensions[1], 1],
            data: mass_field,
            origin: [0.0; 3],
            cell: [
                self.lattice_spacing_m,
                self.lattice_spacing_m,
                self.slab_depth_m,
            ],
        })
    }
}

impl PourOutcome {
    /// Render this accepted pour using explicit SI material/view parameters.
    ///
    /// `dimensions` must be the original simulation shape. This uses the same
    /// borrowed mass buffer as the authored preview, with an independently
    /// supplied physical slab thickness and absorption calibration.
    pub fn render_transmittance(
        &self,
        dimensions: [usize; 2],
        view: PourOpticalView,
        settings: AbsorptionImageSettings,
    ) -> Result<Vec<f64>, PourOpticalError> {
        view.render_mass_field(&self.mass_field, dimensions, settings)
    }

    /// Observe this accepted mass field at explicitly supplied wavelengths.
    ///
    /// The original simulation dimensions, SI calibration and spectral work
    /// limit remain explicit. Output is pixel-major with wavelengths fastest.
    pub fn render_spectral_transmittance(
        &self,
        dimensions: [usize; 2],
        view: PourSpectralOpticalView<'_>,
        settings: SpectralAbsorptionImageSettings,
    ) -> Result<Vec<f64>, PourOpticalError> {
        view.render_transmittance(&self.mass_field, dimensions, settings)
    }
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
    run_pour_motion_observed(rig, contact, law, pivot, |_, _| Ok::<(), PourError>(()))
}

// One integration path for legacy fixtures and physically bound material pours.
// An observer sees only completed native states; it cannot alter the trajectory.
fn run_pour_motion_observed<E: From<PourError>>(
    rig: &PourRig,
    contact: ContactModel,
    law: Rheology,
    pivot: Option<[f64; 2]>,
    observe: impl FnMut(&FreeSurface, u32) -> Result<(), E>,
) -> Result<PourOutcome, E> {
    run_pour_motion_integrated(
        rig,
        contact,
        pivot,
        |sim, step| advance_pour_liquid(sim, law, step).map_err(E::from),
        observe,
    )
}

// Constitutive bindings provide one native advance; geometry, frame motion,
// endpoint observation and outcome construction remain shared by every pour.
fn run_pour_motion_integrated<E: From<PourError>>(
    rig: &PourRig,
    contact: ContactModel,
    pivot: Option<[f64; 2]>,
    mut advance: impl FnMut(&mut FreeSurface, u32) -> Result<(), E>,
    mut observe: impl FnMut(&FreeSurface, u32) -> Result<(), E>,
) -> Result<PourOutcome, E> {
    if rig.steps == 0 {
        return Err(PourError::EmptySchedule.into());
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
    observe(&sim, 0)?;
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
        advance(&mut sim, step)?;
        observe(&sim, step + 1)?;
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
/// Use [`PourOutcome::render_transmittance`] with [`PourOpticalView`] for a
/// supplied SI absorption coefficient, density conversion and physical depth.
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

#[cfg(test)]
mod physical_opacity_tests {
    use super::*;

    fn view() -> PourOpticalView {
        // Authored analytic fixture, not measured liquid spectroscopy.
        PourOpticalView {
            wavelength_nm: 550.0,
            mass_absorption_m2_kg: 0.4,
            density_unit_kg_m3: 1000.0,
            lattice_spacing_m: 1e-3,
            slab_depth_m: 2e-3,
        }
    }

    fn settings(cells: usize) -> AbsorptionImageSettings {
        AbsorptionImageSettings {
            resolution: [2; 2],
            max_pixels: 4,
            max_grid_cells: cells,
            max_cell_visits: cells as u64,
        }
    }

    fn close(actual: f64, expected: f64) {
        assert!((actual - expected).abs() < 8e-15, "{actual} != {expected}");
    }

    /// G1: the accepted full/interface-cell masses use the supplied SI column.
    #[test]
    fn physical_pour_absorption_uses_density_coefficient_and_viewing_depth() {
        let outcome = PourOutcome {
            mass_drift: 0.0,
            fragments: 1,
            poured_mass: 0.0,
            dribble_cells: 0,
            mass_field: vec![0.0, 0.25, 0.5, 1.0],
        };
        let optics = view();
        let image = outcome
            .render_transmittance([2, 2], optics, settings(4))
            .unwrap();
        for (&transmittance, &mass) in image.iter().zip(&outcome.mass_field) {
            close(transmittance, (-mass * 0.8).exp());
        }
        let deeper = PourOpticalView {
            slab_depth_m: 2.0 * optics.slab_depth_m,
            ..optics
        }
        .render_mass_field(&outcome.mass_field, [2, 2], settings(4))
        .unwrap();
        for (&twice, &once) in deeper.iter().zip(&image) {
            close(twice, once * once);
        }
        let red = PourOpticalView {
            wavelength_nm: 700.0,
            mass_absorption_m2_kg: 0.1,
            ..optics
        }
        .render_mass_field(&outcome.mass_field, [2, 2], settings(4))
        .unwrap();
        close(red[3], (-0.2_f64).exp());
        assert!(red[3] > image[3]);

        let vacuum = PourOpticalView {
            mass_absorption_m2_kg: 0.0,
            ..optics
        }
        .render_mass_field(&outcome.mass_field, [2, 2], settings(4))
        .unwrap();
        assert_eq!(vacuum, vec![1.0; 4]);
    }

    /// G3: lateral refinement and density-unit choices cannot change opacity.
    #[test]
    fn physical_pour_opacity_is_independent_of_lattice_resolution_and_mass_units() {
        let coarse = [0.0, 0.25, 0.5, 1.0];
        let optics = view();
        let reference = optics.render_mass_field(&coarse, [2, 2], settings(4)).unwrap();
        let fine: Vec<_> = (0..4)
            .flat_map(|y| (0..4).map(move |x| coarse[x / 2 + 2 * (y / 2)]))
            .collect();
        let refined = PourOpticalView {
            lattice_spacing_m: 0.5 * optics.lattice_spacing_m,
            ..optics
        }
        .render_mass_field(&fine, [4, 4], settings(16))
        .unwrap();
        assert_eq!(refined, reference);
        let rescaled = PourOpticalView {
            density_unit_kg_m3: 0.5 * optics.density_unit_kg_m3,
            ..optics
        }
        .render_mass_field(&coarse.map(|mass| 2.0 * mass), [2, 2], settings(4))
        .unwrap();
        assert_eq!(rescaled, reference);
    }

    /// G1: a subpixel stream contributes its area even when the centre misses it.
    #[test]
    fn physical_pour_filters_transmission_before_downsampling() {
        let optics = PourOpticalView {
            mass_absorption_m2_kg: 1000.0,
            density_unit_kg_m3: 1.0,
            slab_depth_m: 1.0,
            ..view()
        };
        let mass = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
        let limits = AbsorptionImageSettings {
            resolution: [1, 1], max_pixels: 1,
            max_grid_cells: 8, max_cell_visits: 8,
        };
        assert_eq!(optics.render_mass_field(&mass, [4, 2], limits).unwrap(), vec![0.75]);
        let refined: Vec<_> = (0..4)
            .flat_map(|y| (0..8).map(move |x| mass[x / 2 + 4 * (y / 2)]))
            .collect();
        let refined_limits = AbsorptionImageSettings {
            max_grid_cells: 32, max_cell_visits: 32, ..limits
        };
        assert_eq!(
            PourOpticalView { lattice_spacing_m: 0.5 * optics.lattice_spacing_m, ..optics }
                .render_mass_field(&refined, [8, 4], refined_limits).unwrap(),
            vec![0.75]
        );
        assert!(optics.render_mass_field(
            &mass, [4, 2], AbsorptionImageSettings { max_cell_visits: 7, ..limits },
        ).is_err());
    }

    /// G1: a coarse pixel preserves separate transparent and absorbing columns.
    #[test]
    fn physical_pour_coarse_pixels_do_not_use_the_mean_density() {
        let mass = [0.0, 0.25, 0.5, 1.0];
        let limits = AbsorptionImageSettings {
            resolution: [1, 1], max_pixels: 1,
            max_grid_cells: 4, max_cell_visits: 4,
        };
        let actual = view().render_mass_field(&mass, [2, 2], limits).unwrap()[0];
        let expected = mass.iter().map(|rho| (-0.8_f64 * rho).exp()).sum::<f64>() / 4.0;
        close(actual, expected);
        assert!((actual - (-0.8_f64 * 0.4375).exp()).abs() > 0.02);
    }

    /// G0: invalid source/view input and bounded work never return an image.
    #[test]
    fn physical_pour_refuses_invalid_calibration_shape_and_budget() {
        let mass = [0.0, 0.25, 0.5, 1.0];
        let optics = view();
        for invalid in [
            PourOpticalView { wavelength_nm: f64::NAN, ..optics },
            PourOpticalView { wavelength_nm: 0.0, ..optics },
            PourOpticalView { lattice_spacing_m: 0.0, ..optics },
            PourOpticalView { slab_depth_m: f64::INFINITY, ..optics },
            PourOpticalView { density_unit_kg_m3: -1.0, ..optics },
            PourOpticalView { mass_absorption_m2_kg: -1.0, ..optics },
        ] {
            assert!(invalid.render_mass_field(&mass, [2, 2], settings(4)).is_err());
        }
        assert!(optics.render_mass_field(&mass, [3, 2], settings(4)).is_err());
        assert!(optics.render_mass_field(&mass, [usize::MAX, 2], settings(4)).is_err());
        assert!(optics.render_mass_field(&[-0.1, 0.0, 0.0, 0.0], [2, 2], settings(4)).is_err());
        let limited = AbsorptionImageSettings { max_pixels: 3, ..settings(4) };
        assert!(optics.render_mass_field(&mass, [2, 2], limited).is_err());
    }

    /// G1/G3: each wavelength keeps the filament's area and the SI slab law.
    #[test]
    fn spectral_pour_filters_each_wavelength_and_matches_the_monochromatic_view() {
        let mass = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
        let wavelengths = [450.0, 550.0, 650.0];
        let coefficients = [0.0, 0.4, 1000.0];
        let optics = PourSpectralOpticalView {
            wavelengths_nm: &wavelengths,
            mass_absorption_m2_kg: &coefficients,
            density_unit_kg_m3: 1000.0,
            lattice_spacing_m: 1e-3,
            slab_depth_m: 2e-3,
        };
        let limits = SpectralAbsorptionImageSettings {
            image: AbsorptionImageSettings {
                resolution: [1, 1],
                max_pixels: 1,
                max_grid_cells: 8,
                max_cell_visits: 8,
            },
            max_spectral_evaluations: 24,
        };
        let image = optics.render_mass_field(&mass, [4, 2], limits).unwrap();
        assert_eq!(image.resolution(), [1, 1]);
        assert_eq!(image.wavelengths_nm(), &wavelengths);
        assert_eq!(image.samples().len(), 3);
        let expected = [1.0, 0.75 + 0.25 * (-0.8_f64).exp(), 0.75];
        for (sample, expected) in image.samples().iter().zip(expected) {
            close(sample.transmittance(), expected);
            close(sample.opacity(), 1.0 - expected);
            close(sample.optical_depth(), -expected.ln());
        }
        let direct = optics.render_transmittance(&mass, [4, 2], limits).unwrap();
        assert_eq!(
            direct,
            image
                .samples()
                .iter()
                .map(|sample| sample.transmittance())
                .collect::<Vec<_>>()
        );
        for (band, &wavelength_nm) in wavelengths.iter().enumerate() {
            let single = PourSpectralOpticalView {
                wavelengths_nm: &wavelengths[band..=band],
                mass_absorption_m2_kg: &coefficients[band..=band],
                ..optics
            }
            .render_transmittance(&mass, [4, 2], limits)
            .unwrap();
            let monochromatic = PourOpticalView {
                wavelength_nm,
                mass_absorption_m2_kg: coefficients[band],
                density_unit_kg_m3: optics.density_unit_kg_m3,
                lattice_spacing_m: optics.lattice_spacing_m,
                slab_depth_m: optics.slab_depth_m,
            }
            .render_mass_field(&mass, [4, 2], limits.image)
            .unwrap();
            close(single[0], monochromatic[0]);
            close(single[0], direct[band]);
        }
        let outcome = PourOutcome {
            mass_drift: 0.0,
            fragments: 1,
            poured_mass: 0.0,
            dribble_cells: 0,
            mass_field: mass.to_vec(),
        };
        assert_eq!(
            outcome
                .render_spectral_transmittance([4, 2], optics, limits)
                .unwrap(),
            direct
        );
    }

    /// G0: the complete spectrum and all wavelengths' work must be admitted.
    #[test]
    fn spectral_pour_refuses_invalid_spectra_geometry_and_total_work() {
        let mass = [0.0, 0.25, 0.5, 1.0];
        let optics = PourSpectralOpticalView {
            wavelengths_nm: &[450.0, 650.0],
            mass_absorption_m2_kg: &[0.4, 0.1],
            density_unit_kg_m3: 1000.0,
            lattice_spacing_m: 1e-3,
            slab_depth_m: 2e-3,
        };
        let limits = SpectralAbsorptionImageSettings {
            image: settings(4),
            max_spectral_evaluations: 8,
        };
        let insufficient = SpectralAbsorptionImageSettings {
            max_spectral_evaluations: 7,
            ..limits
        };
        assert!(matches!(
            optics.render_transmittance(&mass, [2, 2], insufficient),
            Err(PourOpticalError::Transport(
                DvrError::SpectralEvaluationBudgetExceeded {
                    requested: 8,
                    limit: 7
                }
            ))
        ));
        assert!(
            optics
                .render_mass_field(&mass, [2, 2], insufficient)
                .is_err()
        );
        let invalid: [(&[f64], &[f64]); 6] = [
            (&[], &[]),
            (&[450.0], &[0.4, 0.1]),
            (&[650.0, 450.0], &[0.4, 0.1]),
            (&[450.0, 450.0], &[0.4, 0.1]),
            (&[450.0, f64::NAN], &[0.4, 0.1]),
            (&[450.0, 650.0], &[0.4, -0.1]),
        ];
        for (wavelengths_nm, mass_absorption_m2_kg) in invalid {
            let invalid = PourSpectralOpticalView {
                wavelengths_nm,
                mass_absorption_m2_kg,
                ..optics
            };
            assert!(invalid.render_transmittance(&mass, [2, 2], limits).is_err());
        }
        assert!(matches!(
            PourSpectralOpticalView {
                slab_depth_m: 0.0,
                ..optics
            }
            .render_mass_field(&mass, [2, 2], limits),
            Err(PourOpticalError::InvalidGeometry)
        ));
        assert!(optics.render_transmittance(&mass, [3, 2], limits).is_err());
    }
}
