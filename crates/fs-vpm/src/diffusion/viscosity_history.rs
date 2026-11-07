//! Explicitly interpolated, spatially homogeneous time-dependent viscosity.
//!
//! The existing PSE pair law is linear in viscosity. On fixed sites its
//! coefficient exposure is therefore the physical integral of `nu(t)`. Moving
//! Strang evolution must use the separate exposures before and after advection;
//! replacing both halves by one full-step mean changes that ordering.

use super::{PseDiffusion, PseKernel};
use crate::{VpmError, reserved_vec};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

mod profile;
pub use profile::{ViscosityHistory, ViscosityHistoryError, ViscosityInterval, ViscosityKnot};

/// Original history, material or native numerical refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum PseViscosityError {
    /// Invalid support or an unrepresentable interpolated quantity.
    History(ViscosityHistoryError),
    /// Original fluid source admission failure.
    Material(MaterialError),
    /// Original native kernel or allocation refusal.
    Diffusion(VpmError),
    /// This incompressible homogeneous owner does not evolve density or areas.
    DensityChanged { sample: usize },
}

impl core::fmt::Display for PseViscosityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::History(error) => error.fmt(f),
            Self::Material(error) => write!(f, "viscosity history material: {error}"),
            Self::Diffusion(error) => error.fmt(f),
            Self::DensityChanged { sample } => {
                write!(f, "viscosity history density changes at sample {sample}")
            }
        }
    }
}
impl std::error::Error for PseViscosityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::History(error) => Some(error),
            Self::Material(error) => Some(error),
            Self::Diffusion(error) => Some(error),
            Self::DensityChanged { .. } => None,
        }
    }
}
impl From<ViscosityHistoryError> for PseViscosityError {
    fn from(error: ViscosityHistoryError) -> Self {
        Self::History(error)
    }
}
impl From<VpmError> for PseViscosityError {
    fn from(error: VpmError) -> Self {
        Self::Diffusion(error)
    }
}

/// One authored homogeneous `nu(t)` interpolant with a fixed spatial kernel.
///
/// The interpolation is piecewise linear in kinematic viscosity, not inferred
/// interpolation of temperature, pressure, a material's constitutive law or
/// source receipts. Knot cuts and exact half-interval means remove coefficient
/// quadrature error for this interpolant; mechanical time and cloud refinement
/// are still required. No spatial viscosity gradient or density evolution is
/// supplied. Zero viscosity is allowed only as an explicit authored value.
#[derive(Debug, Clone, Copy)]
pub struct PseViscosityHistory<'a> {
    profile: ViscosityHistory<'a>,
    template: PseDiffusion,
}

impl<'a> PseViscosityHistory<'a> {
    /// Admit the declared linear coefficient history and original PSE kernel.
    ///
    /// # Errors
    /// Refuses unordered/nonfinite samples and invalid native kernel width.
    pub fn piecewise_linear(
        knots: &'a [ViscosityKnot],
        kernel_width_m: f64,
        kernel: PseKernel,
    ) -> Result<Self, PseViscosityError> {
        Ok(Self {
            profile: ViscosityHistory::piecewise_linear(knots)?,
            template: PseDiffusion::with_kernel(0.0, kernel_width_m, kernel)?,
        })
    }

    /// The admitted authored coefficient history.
    #[must_use]
    pub const fn profile(self) -> ViscosityHistory<'a> {
        self.profile
    }

    /// Original native spatial kernel family.
    #[must_use]
    pub const fn kernel(self) -> PseKernel {
        self.template.kernel()
    }

    /// Gaussian width or compact support radius [m].
    #[must_use]
    pub const fn kernel_width_m(self) -> f64 {
        self.template.kernel_width_m()
    }

    /// Native law using the exactly integrated mean on any supported window.
    /// This is suitable for one fixed-site sweep. Moving Strang consumers must
    /// request each physical half separately, as the native time-dependent
    /// conservative owner does automatically.
    pub fn mean_law(self, start_s: f64, end_s: f64) -> Result<PseDiffusion, PseViscosityError> {
        Ok(PseDiffusion::with_kernel(
            self.profile.mean_viscosity(start_s, end_s)?,
            self.kernel_width_m(),
            self.kernel(),
        )?)
    }

    pub(crate) fn half_laws(
        self,
        interval: ViscosityInterval,
    ) -> Result<[PseDiffusion; 2], PseViscosityError> {
        let [first, second] = interval.half_mean_viscosities_m2_s;
        Ok([
            PseDiffusion::with_kernel(first, self.kernel_width_m(), self.kernel())?,
            PseDiffusion::with_kernel(second, self.kernel_width_m(), self.kernel())?,
        ])
    }
}

/// Original fluid sources retained beside their complete evaluated snapshots.
///
/// Every source is evaluated once at construction. The caller explicitly
/// chooses piecewise-linear interpolation of the resulting `mu/rho` values;
/// this is not a claim that the source constitutive law is exactly linear in
/// time. A common density is required because this incompressible vortex owner
/// transports fixed quadrature areas. Sources remain borrowed so their own
/// receipts, state points and applicability limits are not replaced by raw
/// coefficients. Raw transport snapshots retain only their raw authority.
#[derive(Debug)]
pub struct MaterialPseViscosityHistory<'m, F: FluidTransport> {
    sources: &'m [F],
    transport: Vec<FluidTransportState>,
    knots: Vec<ViscosityKnot>,
    template: PseDiffusion,
}

impl<'m, F: FluidTransport> MaterialPseViscosityHistory<'m, F> {
    /// Explicitly interpolate the original source snapshots at authored times.
    /// A slice of `&dyn FluidTransport` can hold different source-owner types.
    ///
    /// # Errors
    /// Refuses mismatched/invalid times, original material errors, changing
    /// density, invalid kernel width or failed fallible reservation. No binding
    /// escapes until all snapshots have been admitted.
    pub fn piecewise_linear(
        sources: &'m [F],
        times_s: &[f64],
        kernel_width_m: f64,
        kernel: PseKernel,
    ) -> Result<Self, PseViscosityError> {
        if sources.len() != times_s.len() {
            return Err(VpmError::LengthMismatch {
                field: "viscosity source times",
                expected: sources.len(),
                actual: times_s.len(),
            }
            .into());
        }
        let template = PseDiffusion::with_kernel(0.0, kernel_width_m, kernel)?;
        let mut knots = reserved_vec("viscosity history times", times_s.len())?;
        for &time_s in times_s {
            knots.push(ViscosityKnot {
                time_s,
                kinematic_viscosity_m2_s: 0.0,
            });
        }
        ViscosityHistory::piecewise_linear(&knots)?;
        let mut transport: Vec<FluidTransportState> =
            reserved_vec("viscosity history snapshots", sources.len())?;
        for (index, source) in sources.iter().enumerate() {
            let state = source
                .transport_state()
                .map_err(PseViscosityError::Material)?;
            if transport
                .first()
                .is_some_and(|first| first.density_kg_m3() != state.density_kg_m3())
            {
                return Err(PseViscosityError::DensityChanged { sample: index });
            }
            knots[index].kinematic_viscosity_m2_s = state.kinematic_viscosity_m2_s();
            transport.push(state);
        }
        Ok(Self {
            sources,
            transport,
            knots,
            template,
        })
    }

    /// Original borrowed source owners in authored time order.
    #[must_use]
    pub const fn sources(&self) -> &'m [F] {
        self.sources
    }

    /// Complete original snapshots used for the coefficient knots.
    #[must_use]
    pub fn transport_states(&self) -> &[FluidTransportState] {
        &self.transport
    }

    /// Borrow the admitted interpolant, retaining this source binding beside
    /// the returned numerical run whenever provenance is required.
    #[must_use]
    pub fn history(&self) -> PseViscosityHistory<'_> {
        PseViscosityHistory {
            // Construction has admitted these private immutable vectors and
            // FluidTransportState can only supply positive finite viscosity.
            profile: ViscosityHistory::from_admitted(&self.knots),
            template: self.template,
        }
    }
}
