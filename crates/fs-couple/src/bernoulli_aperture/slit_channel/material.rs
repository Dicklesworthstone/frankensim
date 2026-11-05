//! Source-bearing incompressible Newtonian channels for liquids and gases.
//!
//! One original FluidTransport owner supplies the coherent density and dynamic
//! viscosity. Geometry is assembled by the SAME planar, circular or annular
//! constructor as the gas channel, and MomentumProfile retains every transverse
//! velocity, pressure port, wall impulse and physical/numerical energy account.
//! No synthetic gas, hydraulic-diameter substitution, fitted resistance or
//! resampling of a material source occurs during stepping.
//!
//! A source change is a cold, immutable proposal at fixed geometry and EXACTLY
//! fixed density. It preserves every cell velocity, kinetic energy and historical
//! exchange while changing future viscous diffusion. Density changes require a
//! separate mass/parameter-work owner and refuse here. Temperature, heat
//! conduction, species transport, phase changes and transition prediction remain
//! with their physical owners; a resolved transport snapshot does not evolve them.
//! Viscosity is uniform and frozen within each step. An apparent-viscosity
//! snapshot does not turn this linear owner into a shear-dependent rheology solve.

use core::fmt;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

use super::{
    AcousticRealizeError, MomentumProfile, SlitChannel, SlitChannelLoad, SlitChannelSpec,
    SlitChannelStep,
    annular::{AnnularChannel, AnnularChannelSpec},
    circular::{CircularChannel, CircularChannelSpec},
    harmonic::ChannelContinuousOperator,
};

/// Original geometry, resolution, time interval and pressure datum.
/// Each variant uses its own physical cell areas and face shear coefficients.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MaterialChannelGeometry {
    Slit(SlitChannelSpec),
    Circular(CircularChannelSpec),
    Annular(AnnularChannelSpec),
}

impl MaterialChannelGeometry {
    pub const fn sample_period_s(self) -> f64 {
        match self {
            Self::Slit(spec) => spec.sample_period_s,
            Self::Circular(spec) => spec.sample_period_s,
            Self::Annular(spec) => spec.sample_period_s,
        }
    }

    fn profile(
        self,
        transport: FluidTransportState,
    ) -> Result<MomentumProfile, AcousticRealizeError> {
        let density = transport.density_kg_m3();
        let viscosity = transport.dynamic_viscosity_pa_s();
        match self {
            Self::Slit(spec) => SlitChannel::profile_from_transport(density, viscosity, spec),
            Self::Circular(spec) => {
                CircularChannel::profile_from_transport(density, viscosity, spec)
            }
            Self::Annular(spec) => AnnularChannel::profile_from_transport(density, viscosity, spec),
        }
    }
}

/// Original material failure or original channel refusal, without erased source diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialChannelError {
    Material(MaterialError),
    Channel(AcousticRealizeError),
    /// This fixed-volume owner cannot change fluid mass or pay density parameter work.
    DensityChange {
        previous_kg_m3: f64,
        proposed_kg_m3: f64,
    },
}

impl From<MaterialError> for MaterialChannelError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<AcousticRealizeError> for MaterialChannelError {
    fn from(error: AcousticRealizeError) -> Self {
        Self::Channel(error)
    }
}
impl fmt::Display for MaterialChannelError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Material(error) => write!(formatter, "channel material: {error}"),
            Self::Channel(error) => write!(formatter, "material channel: {error}"),
            Self::DensityChange {
                previous_kg_m3,
                proposed_kg_m3,
            } => write!(
                formatter,
                "fixed-density channel cannot change density from {previous_kg_m3} to {proposed_kg_m3} kg/m3"
            ),
        }
    }
}
impl std::error::Error for MaterialChannelError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Material(error) => Some(error),
            Self::Channel(error) => Some(error),
            Self::DensityChange { .. } => None,
        }
    }
}

/// Original source, evaluated transport and resolved no-slip momentum history.
/// Borrowed material-card owners work through FluidTransport for &T; stepping
/// does not require Clone, does not query the source, and adds no allocation.
#[derive(Debug, Clone)]
pub struct MaterialChannel<M> {
    source: M,
    transport: FluidTransportState,
    geometry: MaterialChannelGeometry,
    profile: MomentumProfile,
}

impl<M: FluidTransport> MaterialChannel<M> {
    /// Resolve the original source ONCE and construct a quiescent channel.
    /// Material-card receipts, state-point identity and constitutive choices
    /// remain on source(). Raw FluidTransportState input remains explicitly raw.
    pub fn try_new(
        source: M,
        geometry: MaterialChannelGeometry,
    ) -> Result<Self, MaterialChannelError> {
        let transport = source.transport_state()?;
        Self::from_transport(source, transport, geometry)
    }

    fn from_transport(
        source: M,
        transport: FluidTransportState,
        geometry: MaterialChannelGeometry,
    ) -> Result<Self, MaterialChannelError> {
        let profile = geometry.profile(transport)?;
        Ok(Self {
            source,
            transport,
            geometry,
            profile,
        })
    }

    pub const fn source(&self) -> &M {
        &self.source
    }
    pub const fn transport(&self) -> FluidTransportState {
        self.transport
    }
    pub const fn geometry(&self) -> MaterialChannelGeometry {
        self.geometry
    }
    pub const fn sample_period_s(&self) -> f64 {
        self.geometry.sample_period_s()
    }
    pub fn velocities_m_s(&self) -> &[f64] {
        self.profile.velocities()
    }
    pub const fn last_step(&self) -> SlitChannelStep {
        self.profile.last_step()
    }
    pub const fn stored_energy_j(&self) -> f64 {
        self.profile.stored_energy_j()
    }

    /// Propose another original source with unchanged physical density.
    /// This is a zero-duration intervention: no pressure/wall work, heat, clock
    /// advance or flow-profile reset is fabricated. Coefficients may change,
    /// but temperature or composition evolution is not solved by this method.
    pub fn with_source<N: FluidTransport>(
        &self,
        source: N,
    ) -> Result<MaterialChannel<N>, MaterialChannelError> {
        self.with_source_controlled(source, &mut || false)
    }

    /// Cancellable material continuation. The old source and accepted field
    /// remain untouched on every refusal, including final cancellation.
    /// Construction is bounded by the original 4096-cell limit, and the shared
    /// history transfer polls each 64 cells. Cold construction allocates.
    pub fn with_source_controlled<N: FluidTransport>(
        &self,
        source: N,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MaterialChannel<N>, MaterialChannelError> {
        super::gas::poll(cancelled)?;
        let transport = source.transport_state()?;
        if transport.density_kg_m3() != self.transport.density_kg_m3() {
            return Err(MaterialChannelError::DensityChange {
                previous_kg_m3: self.transport.density_kg_m3(),
                proposed_kg_m3: transport.density_kg_m3(),
            });
        }
        super::gas::poll(cancelled)?;
        let mut candidate = MaterialChannel::from_transport(source, transport, self.geometry)?;
        candidate
            .profile
            .transfer_material_history(&self.profile, 0.0, cancelled)?;
        super::gas::poll(cancelled)?;
        Ok(candidate)
    }

    /// Borrow the original continuous cell equations for stationary-boundary
    /// harmonic analysis. The nonlinear kinetic-head terminal refuses this view.
    pub fn continuous_operator(
        &self,
    ) -> Result<ChannelContinuousOperator<'_>, MaterialChannelError> {
        Ok(self.profile.continuous_operator()?)
    }

    /// Original next-step affine pressure/flow relation with stationary walls.
    /// The retained history is condensed, never replaced by a steady profile.
    pub fn pressure_load(&self) -> Result<SlitChannelLoad, MaterialChannelError> {
        Ok(self.profile.pressure_load()?)
    }

    /// Original affine pressure port with prescribed axial boundary drive.
    /// Wall work remains separate in the accepted step's mechanical account.
    pub fn pressure_load_with_wall_velocities(
        &self,
        walls_m_s: [f64; 2],
    ) -> Result<SlitChannelLoad, MaterialChannelError> {
        Ok(self.profile.pressure_load_with_walls(walls_m_s)?)
    }

    pub fn inertial_impedance_pa_s_m3(&self) -> f64 {
        self.profile.inertial_impedance()
    }

    /// Stationary-wall response of the original mesh, used only for cold
    /// admission. Live startup and reversal always retain the resolved field.
    pub fn steady_bulk_speed_m_s(&self, pressure_pa: f64) -> Result<f64, MaterialChannelError> {
        Ok(self.profile.steady_bulk_speed(pressure_pa)?)
    }

    /// Advance the SAME finite-volume momentum owner with stationary walls.
    /// Failure preserves all accepted velocities, energy and historical work.
    pub fn step(&mut self, pressure_pa: f64) -> Result<SlitChannelStep, MaterialChannelError> {
        Ok(self.profile.step(pressure_pa)?)
    }

    /// Resolve pressure plus prescribed axial wall motion. Walls are lower/upper
    /// for a slit, inner/outer for an annulus, and [0, outer] for a circle.
    /// The circular symmetry axis cannot supply traction or mechanical work.
    pub fn step_with_wall_velocities(
        &mut self,
        pressure_pa: f64,
        walls_m_s: [f64; 2],
    ) -> Result<SlitChannelStep, MaterialChannelError> {
        Ok(self.profile.step_with_walls(pressure_pa, walls_m_s)?)
    }
}

impl<M: FluidTransport + Clone> MaterialChannel<M> {
    /// Complete cold trial at a specified volume flow, with no publication.
    /// Clones the retained source and field without re-evaluating transport.
    /// Borrow the original source at construction if copying it is undesirable.
    pub fn preview_at_flow(&self, flow_m3_s: f64) -> Result<Self, MaterialChannelError> {
        self.preview_at_flow_with_wall_velocities(flow_m3_s, [0.0; 2])
    }

    /// Condense and advance the same wall-driven trial. A coupled owner must
    /// admit its actual flow residual and separate wall work before publication.
    pub fn preview_at_flow_with_wall_velocities(
        &self,
        flow_m3_s: f64,
        walls_m_s: [f64; 2],
    ) -> Result<Self, MaterialChannelError> {
        let pressure = self
            .pressure_load_with_wall_velocities(walls_m_s)?
            .pressure_drop_pa(flow_m3_s)?;
        let mut candidate = self.clone();
        candidate.step_with_wall_velocities(pressure, walls_m_s)?;
        Ok(candidate)
    }
}
