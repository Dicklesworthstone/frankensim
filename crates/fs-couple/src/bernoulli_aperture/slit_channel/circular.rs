//! Axisymmetric unsteady no-slip flow in a fixed circular passage.
//!
//! Annular finite volumes resolve rho dv/dt = dp/L + mu/r d(r dv/dr)/dr.
//! Annulus areas carry pressure force, mass and flow; cylindrical face areas
//! carry shear. Zero area at r=0 imposes symmetry without a fictitious inner
//! wall. The outer stationary wall lies half a radial cell from its centre.
//! The SAME momentum/time/work owner also serves parallel-wall slits.
//!
//! This is the incompressible, fully developed, laminar Newtonian reduction at
//! fixed material state. It excludes entrance development, axial taper, moving
//! walls, separation, turbulence and thermodynamic compression; a whole oboe
//! staple is not qualified by a cylindrical fixture. Space/time refinement is
//! required separately from validation of that physical reduction.
//!
//! Steady reference: C. E. Brennen, *An Internet Book on Fluid Dynamics*,
//! <https://brennen.caltech.edu/fluidbook/basicfluiddynamics/Navierstokesexactsolutions/poiseuilleflow.pdf>.

use super::{AcousticRealizeError, ChannelGasUpdate, GasState, MomentumProfile, SlitChannelLoad,
    SlitChannelStep, SlitExit, MAX_TRANSVERSE_CELLS, invalid};

/// Fixed cylindrical geometry and radial/time discretization in SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CircularChannelSpec {
    /// Distance between supplied static-pressure planes [m].
    pub length_m: f64,
    /// Actual circular radius [m], not a hydraulic-radius substitution.
    pub radius_m: f64,
    /// Uniform radial annuli in 2..=4096, independently refinable.
    pub radial_cells: usize,
    /// Backward-Euler period [s], independently refinable.
    pub sample_period_s: f64,
    /// Explicit outlet datum, with the original optional kinetic-head closure.
    pub exit: SlitExit,
}

/// Material-bound circular field, with no editable density/viscosity overrides.
#[derive(Debug, Clone)]
pub struct CircularChannel {
    gas: GasState,
    spec: CircularChannelSpec,
    profile: MomentumProfile,
}

impl CircularChannel {
    /// Construct a quiescent circular passage. Geometry and complete gas state
    /// determine every inertial and viscous coefficient before any sample runs.
    pub fn from_gas(gas: &GasState, spec: CircularChannelSpec)
        -> Result<Self, AcousticRealizeError>
    {
        gas.validate_coherence().map_err(|_| invalid("circular channel needs one coherent material gas"))?;
        if ![spec.length_m, spec.radius_m, spec.sample_period_s].iter()
            .all(|x| x.is_finite() && *x > 0.0)
            || !(2..=MAX_TRANSVERSE_CELLS).contains(&spec.radial_cells) {
            return Err(invalid("circular channel needs positive SI geometry/time and 2..=4096 radial cells"));
        }
        let n = spec.radial_cells;
        let dr = spec.radius_m / n as f64;
        let ring = core::f64::consts::PI * dr * dr;
        let shear = 2.0 * core::f64::consts::PI * gas.dynamic_viscosity * spec.length_m;
        let areas = (0..n).map(|i| ring * (2 * i + 1) as f64).collect();
        let mut conductances: Vec<_> = (0..=n).map(|i| shear * i as f64).collect();
        conductances[n] *= 2.0;
        let profile = MomentumProfile::new(gas.density, spec.length_m,
            spec.sample_period_s, spec.exit, areas, conductances)?;
        Ok(Self { gas: *gas, spec, profile })
    }
    /// Complete source gas driving inertia, shear and coupled admissibility.
    #[must_use] pub const fn gas(&self) -> &GasState { &self.gas }
    /// Original physical cylinder and numerical refinement counts.
    #[must_use] pub const fn spec(&self) -> &CircularChannelSpec { &self.spec }
    /// Accepted axial velocities at radial centres, from the axis to the wall.
    #[must_use] pub fn velocities_m_s(&self) -> &[f64] { self.profile.velocities() }
    /// Last accepted storage, work and separately measured losses.
    #[must_use] pub const fn last_step(&self) -> SlitChannelStep { self.profile.last_step() }
    /// Current kinetic storage, including an accepted zero-duration gas update.
    #[must_use] pub const fn stored_energy_j(&self) -> f64 { self.profile.stored_energy_j() }
    /// Replace the complete source gas at fixed geometry while retaining every
    /// annulus's axial momentum. New viscosity drives the original no-slip
    /// operator. The signed kinetic storage change is parameter work, never
    /// pressure work, heat or numerical dissipation. No physical sample runs;
    /// `last_step()` remains the original historical exchange.
    ///
    /// The finite nonnegative allowance bounds absolute parameter work. All
    /// gas, coefficient, transfer and cancellation checks precede publication.
    /// This prescribed intervention does not solve mass, heat or species flow.
    pub fn rebind_gas_controlled(
        &mut self, gas: &GasState, maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ChannelGasUpdate, AcousticRealizeError> {
        super::gas::admit(gas, maximum_absolute_parameter_work_j, cancelled)?;
        if gas == &self.gas {
            super::gas::poll(cancelled)?;
            return Ok(ChannelGasUpdate::unchanged(self.stored_energy_j()));
        }
        let mut candidate = Self::from_gas(gas, self.spec)?;
        let update = candidate.profile.transfer_material_history(
            &self.profile, maximum_absolute_parameter_work_j, cancelled,
        )?;
        super::gas::poll(cancelled)?;
        *self = candidate;
        Ok(update)
    }
    /// Exact affine boundary of the next implicit field step, retaining inertia.
    pub fn pressure_load(&self) -> Result<SlitChannelLoad, AcousticRealizeError> {
        self.profile.pressure_load()
    }
    /// Candidate field at shared through-flow, leaving original state unchanged.
    pub fn preview_at_flow(&self, flow: f64) -> Result<Self, AcousticRealizeError> {
        let pressure = self.pressure_load()?.pressure_drop_pa(flow)?;
        let mut trial = self.clone();
        trial.step(pressure)?;
        Ok(trial)
    }
    /// Steady response of the actual radial mesh, not a replacement for dynamics.
    pub fn steady_bulk_speed_m_s(&self, pressure: f64) -> Result<f64, AcousticRealizeError> {
        self.profile.steady_bulk_speed(pressure)
    }
    /// Advance one signed physical pressure interval with atomic acceptance.
    pub fn step(&mut self, pressure: f64) -> Result<SlitChannelStep, AcousticRealizeError> {
        self.profile.step(pressure)
    }
}
