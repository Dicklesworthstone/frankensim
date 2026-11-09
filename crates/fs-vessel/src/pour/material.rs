//! Source-bound SI pours through the original rotating-frame weir owner.
//!
//! One frozen Newtonian liquid supplies density, viscosity and capillarity to
//! the native material free-surface planner. That plan drives the existing
//! PourRig and its unchanged D2Q9 integration, contact model and mass ledger.
//! No apparent viscosity is frozen from a source requiring evolving rheology.
//!
//! Coordinates start at the lower-left corner of cell (0,0), so its centre is
//! at (dx/2,dx/2). The flow is two-dimensional: masses are kg per metre of
//! out-of-plane depth, never kilograms for an inferred vessel depth. Geometry,
//! support motion and contact-model choice remain explicit scenario inputs.
//! Speeds are relative to the rotating vessel, the frame used by the native
//! lattice Mach admission.

/// Evolving source viscosity on the same physical pour trajectory.
pub mod rheology;

use fs_exec::CancelGate;
use fs_lbm::film::{FilmError, material::MaterialFilmLiquid};
use fs_lbm::freesurface::{ContactModel, FreeSurface, FreeSurfacePressureError2};
use fs_lbm::material::MaterialLbmError;
use fs_lbm::material::free_surface::rheology::FreeSurfaceRheologyError;
use fs_lbm::material::free_surface::{MaterialFreeSurfaceScaling, plan_free_surface_from_material};
use fs_lbm::rheology::{Rheology, TAU_CAP, TAU_FLOOR};
use fs_lbm::{Cell, Grid};
use fs_material::liquid::surface::ResolvedLiquidSurfaceState;

use super::{PourError, PourOutcome, PourRig, advance_pour_liquid, run_pour_motion_integrated};

/// Grid-aligned weir geometry and physical motion, with numerical limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialPourSpec {
    /// Includes one outer wall row/column at each boundary.
    pub cells: [usize; 2],
    /// Interior weir column and its height in cells, respectively.
    pub lip: [usize; 2],
    pub cell_size_m: f64,
    pub timestep_s: f64,
    /// Nonnegative downward gravity before the vessel rotates [m/s²].
    pub gravity_m_s2: f64,
    /// Final clockwise angle of the native quintic rest-to-rest schedule [rad].
    pub tilt_final_radians: f64,
    /// Rotation pivot in the physical coordinate system described above [m].
    pub pivot_m: [f64; 2],
    /// Positive number of original native steps; duration is steps * timestep.
    pub steps: u32,
    /// Allocation limit, including liquid, gas, weir and exterior wall cells.
    pub maximum_cells: usize,
    /// Bound on actual liquid/interface speed in the rotating vessel frame [m/s].
    /// Checked at every endpoint; the native planner also admits its Mach number.
    pub maximum_speed_m_s: f64,
    /// Bound on abs(rho_lattice - 1), strictly between zero and one.
    /// This bounds numerical compressibility, not a physical liquid EOS.
    pub maximum_density_excursion: f64,
}

#[derive(Debug)]
pub enum MaterialPourError {
    Invalid(&'static str),
    Liquid(FilmError),
    Scaling(MaterialLbmError),
    Pour(PourError),
    Rheology(FreeSurfaceRheologyError),
    Cancelled,
    /// An actual wet-cell observation exceeded a declared physical envelope.
    State {
        step: u32,
        cell: usize,
        reason: &'static str,
    },
}

impl core::fmt::Display for MaterialPourError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "material pour: {reason}"),
            Self::Liquid(error) => write!(f, "material pour liquid: {error}"),
            Self::Scaling(error) => write!(f, "material pour scaling: {error}"),
            Self::Pour(error) => write!(f, "material pour: {error}"),
            Self::Rheology(error) => write!(f, "material pour rheology: {error}"),
            Self::Cancelled => f.write_str("material pour cancelled"),
            Self::State { step, cell, reason } => {
                write!(f, "material pour after step {step}, cell {cell}: {reason}")
            }
        }
    }
}
impl std::error::Error for MaterialPourError {}
impl From<PourError> for MaterialPourError {
    fn from(error: PourError) -> Self {
        Self::Pour(error)
    }
}
impl From<FreeSurfaceRheologyError> for MaterialPourError {
    fn from(error: FreeSurfaceRheologyError) -> Self {
        Self::Rheology(error)
    }
}
impl From<MaterialLbmError> for MaterialPourError {
    fn from(error: MaterialLbmError) -> Self {
        Self::Scaling(error)
    }
}

/// Admitted immutable source, physical scale and original native pour geometry.
/// Preparation performs no material reselection, allocation or physical step.
pub struct MaterialPour<'a> {
    scaling: MaterialFreeSurfaceScaling<'a>,
    spec: MaterialPourSpec,
    rig: PourRig,
    pivot: [f64; 2],
    mass_per_depth_unit_kg_m: f64,
}

/// Physical observations of the same native trajectory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialPourReport {
    pub duration_s: f64,
    pub initial_mass_kg_per_m: f64,
    pub final_mass_kg_per_m: f64,
    /// Located liquid to the right of the weir, from the native mass field.
    pub poured_mass_kg_per_m: f64,
    /// Maximum observed liquid speed relative to the rotating vessel [m/s].
    pub peak_speed_m_s: f64,
    pub maximum_density_excursion: f64,
}

/// Accepted source identity, physical observations and the original render field.
pub struct MaterialPourOutcome<'a> {
    scaling: MaterialFreeSurfaceScaling<'a>,
    native: PourOutcome,
    report: MaterialPourReport,
    mass_per_depth_unit_kg_m: f64,
}

impl<'a> MaterialPour<'a> {
    /// Admit one actual liquid source and a representable physical weir run.
    ///
    /// The source must be frozen Newtonian at its retained conditions.
    /// Temperature, species, pressure-dependent or shear-dependent evolution
    /// requires its own constitutive/thermal owner and is refused here.
    /// Native capillary timestep, reconstruction, Mach and relaxation checks
    /// remain authoritative; no viscosity or surface tension is clipped.
    pub fn try_new(
        material: &'a ResolvedLiquidSurfaceState,
        spec: MaterialPourSpec,
    ) -> Result<Self, MaterialPourError> {
        Self::prepare(spec, || {
            MaterialFilmLiquid::try_new(material).map_err(MaterialPourError::Liquid)?;
            Ok(plan_free_surface_from_material(
                material,
                spec.cell_size_m,
                spec.timestep_s,
                spec.maximum_speed_m_s,
            )?)
        })
    }

    // Both constitutive paths share geometry, SI units, motion and envelope
    // admission. Only the explicitly selected source owner supplies scaling.
    fn prepare(
        spec: MaterialPourSpec,
        plan: impl FnOnce() -> Result<MaterialFreeSurfaceScaling<'a>, MaterialPourError>,
    ) -> Result<Self, MaterialPourError> {
        let [nx, ny] = spec.cells;
        let [lip_x, lip_h] = spec.lip;
        let cells = nx.checked_mul(ny).ok_or(MaterialPourError::Invalid(
            "cell count cannot be represented",
        ))?;
        if nx < 5
            || ny < 4
            || lip_x < 2
            || lip_x > nx - 3
            || lip_h < 2
            || lip_h > ny - 2
            || cells > spec.maximum_cells
            || spec.steps == 0
        {
            return Err(MaterialPourError::Invalid(
                "a bounded nonempty grid must leave liquid below the interior weir and space on both sides",
            ));
        }
        if !spec.gravity_m_s2.is_finite()
            || spec.gravity_m_s2 < 0.0
            || !spec.tilt_final_radians.is_finite()
            || spec.pivot_m.iter().any(|value| !value.is_finite())
            || !spec.maximum_speed_m_s.is_finite()
            || spec.maximum_speed_m_s <= 0.0
            || !spec.maximum_density_excursion.is_finite()
            || spec.maximum_density_excursion <= 0.0
            || spec.maximum_density_excursion >= 1.0
        {
            return Err(MaterialPourError::Invalid(
                "invalid physical motion or liquid envelope",
            ));
        }
        let scaling = plan()?;
        let material = scaling.material();
        if !(TAU_FLOOR..=TAU_CAP).contains(&scaling.flow().tau) {
            return Err(MaterialPourError::Invalid(
                "source viscosity exceeds the native pour relaxation range; rescale dx/dt",
            ));
        }
        let gravity = scaling.gravity_lattice(spec.gravity_m_s2)?;
        if spec.gravity_m_s2 != 0.0 && gravity == 0.0 {
            return Err(MaterialPourError::Invalid(
                "nonzero physical gravity disappeared in conversion",
            ));
        }
        let pivot = spec.pivot_m.map(|value| value / spec.cell_size_m - 0.5);
        if pivot.iter().any(|value| !value.is_finite())
            || !positive(f64::from(spec.steps) * spec.timestep_s)
            || !positive(nx as f64 * spec.cell_size_m)
            || !positive(ny as f64 * spec.cell_size_m)
        {
            return Err(MaterialPourError::Invalid(
                "physical coordinates or duration cannot be represented",
            ));
        }
        let mass_per_depth_unit_kg_m =
            (material.density_kg_m3() * spec.cell_size_m) * spec.cell_size_m;
        if !positive(mass_per_depth_unit_kg_m) || !positive(mass_per_depth_unit_kg_m * cells as f64)
        {
            return Err(MaterialPourError::Invalid(
                "two-dimensional mass scale cannot be represented",
            ));
        }
        let rig = PourRig {
            nx,
            ny,
            lip_x,
            lip_h,
            g0: gravity,
            surface_tension: scaling.surface_tension_lattice(),
            tilt_final: spec.tilt_final_radians,
            steps: spec.steps,
        };
        Ok(Self {
            scaling,
            spec,
            rig,
            pivot,
            mass_per_depth_unit_kg_m,
        })
    }

    pub const fn scaling(&self) -> &MaterialFreeSurfaceScaling<'a> {
        &self.scaling
    }
    pub const fn spec(&self) -> &MaterialPourSpec {
        &self.spec
    }
    pub const fn native_rig(&self) -> &PourRig {
        &self.rig
    }
    /// Original integer-centred lattice coordinate convention.
    pub const fn pivot_lattice(&self) -> [f64; 2] {
        self.pivot
    }

    /// Run the original weir solver with actual source viscosity and tension.
    ///
    /// Vessel-frame velocities and numerical density excursions are checked at
    /// initialization and after every complete native step. Cancellation is
    /// checked at these same boundaries and before publishing the outcome.
    /// A late refusal returns no partial field and leaves this preparation and
    /// its original material unchanged. Bounds are sampled trajectory limits,
    /// not substep, interface-resolution or laminar-transition certificates.
    ///
    /// ContactModel remains the native explicit numerical bracket; this API
    /// does not turn it into a measured contact-angle or 3-D breakup model.
    pub fn run(
        &self,
        contact: ContactModel,
        gate: &CancelGate,
    ) -> Result<MaterialPourOutcome<'a>, MaterialPourError> {
        self.run_with(contact, gate, |sim, step| {
            advance_pour_liquid(
                sim,
                Rheology::Newtonian {
                    nu: self.scaling.flow().viscosity,
                },
                step,
            )
            .map_err(Into::into)
        })
    }

    /// Start the original weir pour with its liquid pressure balancing gravity.
    ///
    /// The existing below-crest liquid mask occupies `lip[1] - 1` complete
    /// layers. Its flat free surface lies halfway between the top wet and first
    /// gas cell centres. Interface cells start fully filled, unlike the
    /// half-filled boundary row in [`Self::run`]. The reported initial mass
    /// includes this full-layer geometry and the bounded numerical density
    /// increase with depth; this is not a mass-preserving pressure projection.
    ///
    /// The positive `free_surface_pressure_pa` is the gas pressure datum and is
    /// retained in the accepted outcome's scaling. Neutral contact is implicit
    /// and surface tension must be zero. The native rest profile balances
    /// collision, streaming and gas reconstruction under the actual `rho * g`
    /// force. After initialization, the same rotating-frame schedule and source
    /// viscosity advance the original solver; no relaxation steps are inserted.
    ///
    /// Density variation remains weakly compressible numerical pressure, not
    /// a liquid equation of state. The existing speed and density-excursion
    /// limits apply at initialization and every completed step. Physical
    /// absolute pressure must also remain positive and representable.
    ///
    /// # Errors
    /// Refuses nonzero capillarity, an invalid pressure datum, an unresolved
    /// hydrostatic profile, cancellation, or any original pour/trajectory limit.
    /// No existing source or prepared pour is changed by a failed run.
    pub fn run_hydrostatic(
        &self,
        free_surface_pressure_pa: f64,
        gate: &CancelGate,
    ) -> Result<MaterialPourOutcome<'a>, MaterialPourError> {
        self.run_hydrostatic_with(free_surface_pressure_pa, gate, |sim, step| {
            advance_pour_liquid(
                sim,
                Rheology::Newtonian {
                    nu: self.scaling.flow().viscosity,
                },
                step,
            )
            .map_err(Into::into)
        })
    }

    fn run_hydrostatic_with(
        &self,
        free_surface_pressure_pa: f64,
        gate: &CancelGate,
        advance: impl FnMut(&mut FreeSurface, u32) -> Result<(), MaterialPourError>,
    ) -> Result<MaterialPourOutcome<'a>, MaterialPourError> {
        checkpoint(gate)?;
        if self.rig.surface_tension != 0.0 {
            return Err(MaterialPourError::Invalid(
                "hydrostatic pour initialization requires zero surface tension",
            ));
        }
        let scaling = self
            .scaling
            .with_pressure_reference(free_surface_pressure_pa)?;
        self.run_initialized_with(
            scaling,
            gate,
            |mut grid| {
                // This is a fresh raw-equilibrium grid, not an accepted flow.
                // The native initializer supplies the physical-rest Guo frame.
                grid.g = [0.0, -self.rig.g0];
                FreeSurface::new_hydrostatic_column(
                    grid,
                    self.rig.lip_h - 1,
                    self.spec.maximum_density_excursion,
                    gate,
                )
                .map_err(|error| match error {
                    FreeSurfacePressureError2::Cancelled => MaterialPourError::Cancelled,
                    error => MaterialPourError::Scaling(error.into()),
                })
            },
            advance,
        )
    }

    fn run_with(
        &self,
        contact: ContactModel,
        gate: &CancelGate,
        advance: impl FnMut(&mut FreeSurface, u32) -> Result<(), MaterialPourError>,
    ) -> Result<MaterialPourOutcome<'a>, MaterialPourError> {
        self.run_initialized_with(
            self.scaling,
            gate,
            |grid| Ok(FreeSurface::new(grid, self.rig.surface_tension, contact)),
            advance,
        )
    }

    fn run_initialized_with(
        &self,
        scaling: MaterialFreeSurfaceScaling<'a>,
        gate: &CancelGate,
        initialize: impl FnOnce(Grid) -> Result<FreeSurface, MaterialPourError>,
        advance: impl FnMut(&mut FreeSurface, u32) -> Result<(), MaterialPourError>,
    ) -> Result<MaterialPourOutcome<'a>, MaterialPourError> {
        checkpoint(gate)?;
        let mut report = MaterialPourReport {
            duration_s: f64::from(self.spec.steps) * self.spec.timestep_s,
            initial_mass_kg_per_m: 0.0,
            final_mass_kg_per_m: 0.0,
            poured_mass_kg_per_m: 0.0,
            peak_speed_m_s: 0.0,
            maximum_density_excursion: 0.0,
        };
        let native = run_pour_motion_integrated(
            &self.rig,
            Some(self.pivot),
            initialize,
            advance,
            |sim, step| {
                checkpoint(gate)?;
                self.observe(sim, step, &scaling, &mut report)?;
                Ok::<(), MaterialPourError>(())
            },
        )?;
        report.poured_mass_kg_per_m = mass_si(native.poured_mass, self.mass_per_depth_unit_kg_m)?;
        for &mass in &native.mass_field {
            mass_si(mass, self.mass_per_depth_unit_kg_m)?;
        }
        checkpoint(gate)?;
        Ok(MaterialPourOutcome {
            scaling,
            native,
            report,
            mass_per_depth_unit_kg_m: self.mass_per_depth_unit_kg_m,
        })
    }

    fn observe(
        &self,
        sim: &FreeSurface,
        step: u32,
        scaling: &MaterialFreeSurfaceScaling<'_>,
        report: &mut MaterialPourReport,
    ) -> Result<(), MaterialPourError> {
        for (cell, flag) in sim.grid.flags.iter().enumerate() {
            if !matches!(flag, Cell::Fluid | Cell::Interface) {
                continue;
            }
            let rho: f64 = sim.grid.f[cell].iter().sum();
            if !positive(rho) {
                return Err(MaterialPourError::State {
                    step,
                    cell,
                    reason: "wet-cell density is not positive and finite",
                });
            }
            if scaling.reference_pressure_pa().is_some() {
                scaling.absolute_pressure_si(rho)?;
            }
            let moments = sim.grid.moments(cell);
            let speed =
                fs_math::det::sqrt(moments.u[0] * moments.u[0] + moments.u[1] * moments.u[1]);
            let speed = self.scaling.velocity_si(speed)?;
            if speed > self.spec.maximum_speed_m_s {
                return Err(MaterialPourError::State {
                    step,
                    cell,
                    reason: "actual vessel-frame liquid speed exceeds maximum_speed_m_s",
                });
            }
            let excursion = (rho - 1.0).abs();
            if excursion > self.spec.maximum_density_excursion {
                return Err(MaterialPourError::State {
                    step,
                    cell,
                    reason: "actual liquid density exceeds maximum_density_excursion",
                });
            }
            report.peak_speed_m_s = report.peak_speed_m_s.max(speed);
            report.maximum_density_excursion = report.maximum_density_excursion.max(excursion);
        }
        let mass = mass_si(sim.ledger_mass(), self.mass_per_depth_unit_kg_m)?;
        if step == 0 {
            report.initial_mass_kg_per_m = mass;
        }
        report.final_mass_kg_per_m = mass;
        Ok(())
    }
}

impl<'a> MaterialPourOutcome<'a> {
    pub const fn scaling(&self) -> &MaterialFreeSurfaceScaling<'a> {
        &self.scaling
    }
    pub const fn native(&self) -> &PourOutcome {
        &self.native
    }
    pub const fn report(&self) -> &MaterialPourReport {
        &self.report
    }
    /// Convert lazily without copying or replacing the native render field.
    pub fn mass_field_kg_per_m(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.native
            .mass_field
            .iter()
            .map(|mass| mass * self.mass_per_depth_unit_kg_m)
    }
}

fn checkpoint(gate: &CancelGate) -> Result<(), MaterialPourError> {
    if gate.is_requested() {
        Err(MaterialPourError::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn mass_si(mass: f64, unit: f64) -> Result<f64, MaterialPourError> {
    let result = mass * unit;
    if !mass.is_finite() || !result.is_finite() || (mass != 0.0 && result == 0.0) {
        Err(MaterialPourError::Invalid(
            "native mass cannot be represented in kg per metre of depth",
        ))
    } else {
        Ok(result)
    }
}
