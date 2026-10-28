//! Unsteady, fully developed flow between fixed parallel no-slip walls.
//!
//! Finite volumes solve rho dv/dt = dp/L + mu d²v/dy² with cell-centred velocities
//! across the physical gap, not a prescribed parabolic profile or fitted
//! inertance. One retained GasState supplies density AND viscosity. Refinement
//! resolves the evolving viscous layers; backward-Euler numerical dissipation
//! remains separate from wall heat.
//!
//! An optional ideal uniform exit kinetic head is solved implicitly with EVERY
//! cell. This terminal approximation does not resolve axial development,
//! nonuniform exit energy flux, separation or the jet. The channel assumes an
//! incompressible Newtonian gas at fixed material state and geometry. Side walls,
//! compressibility, heat transport, moving walls and turbulence are excluded.
//!
//! Reference limits: NETL MFiX FLD01 plane Poiseuille verification:
//! <https://mfix.netl.doe.gov/doc/vvuq-manual/main/html/fluid/fld-01.html>.
//! Unsteady momentum equation and eigenfunction solution:
//! <https://courses.washington.edu/me507/Lecture_Notes/Lecture10_ExactSolutions_UnsteadyPoiseuille.pdf>.

use crate::acoustic_realize::AcousticRealizeError;
use fs_material::gas::GasState;
mod profile;
/// Axisymmetric passages using the same no-slip momentum and work owner.
pub mod circular;
use profile::MomentumProfile;

/// Maximum transverse cells. All hot-step storage is allocated at construction.
pub const MAX_TRANSVERSE_CELLS: usize = 4096;

/// Pressure datum at the end of the fully developed channel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SlitExit {
    /// The supplied pressure difference acts across the channel itself.
    StaticPressure,
    /// Subtract ideal uniform exit kinetic head rho v_bulk |v_bulk|/2.
    /// This explicit terminal approximation is not resolved outflow.
    UniformKineticHead,
}

/// Fixed physical geometry and independent numerical resolution, in SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlitChannelSpec {
    /// Streamwise distance between pressure planes [m].
    pub length_m: f64,
    /// Spanwise width [m]; lateral wall effects are excluded.
    pub width_m: f64,
    /// Distance between the stationary no-slip walls [m].
    pub height_m: f64,
    /// Finite-volume cells across the gap, in 2..=4096.
    pub transverse_cells: usize,
    /// Backward-Euler sample period [s], independent of spatial resolution.
    pub sample_period_s: f64,
    /// Explicit exit pressure datum.
    pub exit: SlitExit,
}

/// An accepted pressure/flow exchange and its discrete physical energy balance.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SlitChannelStep {
    /// Signed endpoint bulk velocity [m/s].
    pub bulk_speed_m_s: f64,
    /// Signed endpoint volume flow [m³/s].
    pub flow_m3_s: f64,
    /// Stored resolved fluid kinetic energy [J].
    pub kinetic_energy_j: f64,
    /// Change in resolved fluid kinetic energy [J].
    pub stored_delta_j: f64,
    /// Supplied pressure times accepted volume [J].
    pub pressure_work_j: f64,
    /// Irreversible no-slip wall-friction heat [J].
    pub viscous_dissipation_j: f64,
    /// Work passed to the exit terminal [J], nonnegative for either flow sign.
    /// It is separate from wall heat and the downstream acoustic mouth ledger.
    pub exit_work_j: f64,
    /// Backward-Euler numerical dissipation [J], not physical heat.
    pub numerical_dissipation_j: f64,
    /// Pressure work minus storage change, wall heat, exit and numerical work [J].
    pub energy_residual_j: f64,
}

/// Exact affine pressure/flow boundary of a fixed channel's next discrete step.
/// Positive volume flow and pressure drop point along the supplied channel.
/// The history term contains the retained velocity field; it is not a pressure
/// source supplied by the caller. It must remain inside the channel energy audit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlitChannelLoad {
    /// Positive coefficient Z in dp = Z Q - history [Pa s/m³].
    pub impedance_pa_s_m3: f64,
    /// Retained inertial history pressure [Pa], signed with the previous flow.
    pub history_pressure_pa: f64,
}
impl SlitChannelLoad {
    /// Pressure drop producing the requested flow in the same implicit step.
    pub fn pressure_drop_pa(self, flow_m3_s: f64) -> Result<f64, AcousticRealizeError> {
        let pressure = self.impedance_pa_s_m3 * flow_m3_s - self.history_pressure_pa;
        if !self.impedance_pa_s_m3.is_finite() || self.impedance_pa_s_m3 <= 0.0
            || !self.history_pressure_pa.is_finite() || !flow_m3_s.is_finite() || !pressure.is_finite() {
            return Err(invalid("unsteady slit requested flow gives nonfinite pressure"));
        }
        Ok(pressure)
    }
}

/// Stateful transverse momentum field with an allocation-free implicit step.
#[derive(Debug, Clone)]
pub struct SlitChannel {
    gas: GasState,
    spec: SlitChannelSpec,
    profile: MomentumProfile,
}

/// Explicit cross-section of one fixed material-gas passage. Both variants
/// share the same finite-volume momentum, pressure condensation and energy
/// owner; their actual geometry is retained instead of a hydraulic-diameter fit.
#[derive(Debug, Clone)]
pub enum FixedChannel {
    /// Parallel, planar, no-slip walls.
    Slit(SlitChannel),
    /// Axisymmetric circular no-slip wall, with symmetry at the axis.
    Circular(circular::CircularChannel),
}
impl From<SlitChannel> for FixedChannel {
    fn from(value: SlitChannel) -> Self { Self::Slit(value) }
}
impl From<circular::CircularChannel> for FixedChannel {
    fn from(value: circular::CircularChannel) -> Self { Self::Circular(value) }
}
impl From<&SlitChannel> for FixedChannel {
    fn from(value: &SlitChannel) -> Self { Self::Slit(value.clone()) }
}
impl From<&circular::CircularChannel> for FixedChannel {
    fn from(value: &circular::CircularChannel) -> Self { Self::Circular(value.clone()) }
}
impl From<&FixedChannel> for FixedChannel {
    fn from(value: &FixedChannel) -> Self { value.clone() }
}
impl FixedChannel {
    /// Original planar field when the supplied geometry is a slit.
    pub const fn slit(&self) -> Option<&SlitChannel> {
        match self { Self::Slit(v) => Some(v), Self::Circular(_) => None }
    }
    /// Original annular field when the supplied geometry is circular.
    pub const fn circular(&self) -> Option<&circular::CircularChannel> {
        match self { Self::Circular(v) => Some(v), Self::Slit(_) => None }
    }
    /// Original complete material state; no caller-owned density override.
    pub const fn gas(&self) -> &GasState {
        match self { Self::Slit(v) => v.gas(), Self::Circular(v) => v.gas() }
    }
    /// One unchanged physical integration period [s].
    pub const fn sample_period_s(&self) -> f64 {
        match self { Self::Slit(v) => v.spec.sample_period_s, Self::Circular(v) => v.spec().sample_period_s }
    }
    /// Original accepted transverse velocity field.
    pub fn velocities_m_s(&self) -> &[f64] {
        match self { Self::Slit(v) => v.velocities_m_s(), Self::Circular(v) => v.velocities_m_s() }
    }
    /// Original accepted fluid work and storage account.
    pub const fn last_step(&self) -> SlitChannelStep {
        match self { Self::Slit(v) => v.last_step(), Self::Circular(v) => v.last_step() }
    }
    /// Exact next-step affine pressure boundary of the actual geometry.
    pub fn pressure_load(&self) -> Result<SlitChannelLoad, AcousticRealizeError> {
        match self { Self::Slit(v) => v.pressure_load(), Self::Circular(v) => v.pressure_load() }
    }
    /// Candidate original field at the common through-flow; no state publishes.
    pub fn preview_at_flow(&self, flow: f64) -> Result<Self, AcousticRealizeError> {
        match self {
            Self::Slit(v) => v.preview_at_flow(flow).map(Self::Slit),
            Self::Circular(v) => v.preview_at_flow(flow).map(Self::Circular),
        }
    }
    /// Advance with the shared implicit field and work acceptance.
    pub fn step(&mut self, pressure: f64) -> Result<SlitChannelStep, AcousticRealizeError> {
        match self { Self::Slit(v) => v.step(pressure), Self::Circular(v) => v.step(pressure) }
    }
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

impl SlitChannel {
    /// Build a channel at rest from a complete, coherent gas state.
    ///
    /// # Errors
    /// Invalid material state, geometry, resolution, or unrepresentable matrices.
    pub fn from_gas(gas: &GasState, spec: SlitChannelSpec) -> Result<Self, AcousticRealizeError> {
        gas.validate_coherence().map_err(|_| invalid("unsteady slit requires one coherent material gas"))?;
        if ![spec.length_m, spec.width_m, spec.height_m, spec.sample_period_s]
            .iter().all(|x| x.is_finite() && *x > 0.0)
            || !(2..=MAX_TRANSVERSE_CELLS).contains(&spec.transverse_cells)
        {
            return Err(invalid("unsteady slit requires positive SI dimensions/time and 2..=4096 transverse cells"));
        }
        let n = spec.transverse_cells;
        let dy = spec.height_m / n as f64;
        let area = spec.width_m * dy;
        let conductance = gas.dynamic_viscosity * spec.length_m * spec.width_m / dy;
        if ![dy, area, conductance]
            .iter().all(|x| x.is_finite() && *x > 0.0)
        {
            return Err(invalid("unsteady slit coefficients cannot be represented"));
        }
        let mut faces = vec![conductance; n + 1];
        // Both stationary walls lie half a cell from the adjacent centre.
        faces[0] *= 2.0;
        faces[n] *= 2.0;
        let profile = MomentumProfile::new(gas.density, spec.length_m,
            spec.sample_period_s, spec.exit, vec![area; n], faces)?;
        Ok(Self { gas: *gas, spec, profile })
    }

    /// Unchanged source material driving storage and wall friction.
    #[must_use]
    pub const fn gas(&self) -> &GasState { &self.gas }

    /// Original geometry, mesh and time step.
    #[must_use]
    pub const fn spec(&self) -> &SlitChannelSpec { &self.spec }

    /// Actual cell-centre velocities from the lower to the upper wall [m/s].
    #[must_use]
    pub fn velocities_m_s(&self) -> &[f64] { self.profile.velocities() }

    /// Last accepted exchange, all zero before the first step.
    #[must_use]
    pub const fn last_step(&self) -> SlitChannelStep { self.profile.last_step() }

    /// Exact next-step affine load of this fixed, no-slip channel.
    ///
    /// Only StaticPressure has a linear pressure port. A kinetic-head terminal
    /// is nonlinear and refuses this condensation instead of being dropped.
    /// The read-only evaluation allocates one transverse work vector. It changes
    /// no accepted fluid state, and uses the SAME factored diffusion operator as
    /// step(). The positive impedance is constant at fixed geometry/gas/clock;
    /// the history pressure must be recomputed from every accepted field.
    pub fn pressure_load(&self) -> Result<SlitChannelLoad, AcousticRealizeError> {
        self.profile.pressure_load()
    }

    /// Prepare the next channel state at a requested signed volume flow.
    ///
    /// The original state stays unchanged. A coupled owner can publish the
    /// returned channel only after its other participants and complete energy
    /// account pass. The returned last_step() reports the actual floating-point
    /// flow and physical/numerical work from the original stepper; a coupled
    /// owner must include its flow residual in the common admission tolerance.
    /// This preview allocates a complete trial field; no real-time claim follows.
    pub fn preview_at_flow(&self, flow_m3_s: f64) -> Result<Self, AcousticRealizeError> {
        let pressure = self.pressure_load()?.pressure_drop_pa(flow_m3_s)?;
        let mut trial = self.clone();
        trial.step(pressure)?;
        Ok(trial)
    }

    /// Exact steady bulk response of THIS spatial discretization [m/s].
    /// It is useful for cold admission, but never replaces the live field.
    pub fn steady_bulk_speed_m_s(&self, pressure_pa: f64) -> Result<f64, AcousticRealizeError> {
        self.profile.steady_bulk_speed(pressure_pa)
    }

    /// Advance one physical pressure interval, allowing signed pressure and flow.
    /// Failure leaves the accepted velocity field and last exchange unchanged.
    ///
    /// The terminal law makes the mean response a scalar monotone quadratic.
    /// Its analytic solution feeds the complete factored transverse field; no
    /// fitted flow coefficient, nonlinear iteration or one-sample lag appears.
    pub fn step(&mut self, pressure_pa: f64) -> Result<SlitChannelStep, AcousticRealizeError> {
        self.profile.step(pressure_pa)
    }
}
