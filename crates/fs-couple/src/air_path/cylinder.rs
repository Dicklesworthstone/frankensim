//! Causal transverse loading from the unsteady Stokes cylinder continuum.
//!
//! rho u_t = -grad p + mu laplacian u, div u = 0, no slip on a circular
//! cylinder and rest at infinity give, per unit length,
//! Z(s) = pi rho a² s + 4 pi mu z K1(z)/K0(z), z=a sqrt(rho s/mu).
//! See Berg, Yao & Panzarella, NASA/CR-2007-215050, equation (7):
//! <https://ntrs.nasa.gov/api/citations/20080006726/downloads/20080006726.pdf>.
//!
//! Across the negative real s axis the Bessel Wronskian gives a POSITIVE
//! spectral density 8 mu/[pi (J0(x)²+Y0(x)²)] per d log(lambda), where
//! x=a sqrt(rho lambda/mu). A midpoint quadrature therefore produces physical
//! mass/dashpot branches, not fitted poles or per-mode damping coefficients:
//! Z_j(s)=c_j s/(s+lambda_j), m_j=c_j/lambda_j.
//!
//! Resolution has three independent axes: two spectral cutoffs and log-grid
//! intervals. The low tail's exact integrated weight is collapsed to resistance;
//! its error is at most R_low*lambda_min/omega at a positive frequency omega.
//! The high tail is omitted. Neither cutoff nor a grid count certifies a whole
//! performance: compare complex impedance and actual coupled histories on a
//! declared frequency/time band while refining all three axes.
//!
//! This is linear, incompressible, unconfined two-dimensional flow around a long
//! circular cylinder. It omits convection, end/axial interactions, nearby wires,
//! rarefaction, rotary drag and compressible sound radiation. A string's material
//! loss is separate; the older sampled cylinder resistance must not be added.
//! Gas and single-phase liquid sources use the SAME Newtonian continuum through
//! fs-material's checked transport view. The original source is retained, not
//! replaced by a synthetic gas or a second material card. Coefficients are frozen
//! at preparation; evolving temperature/composition and non-Newtonian histories
//! require their own coupled update, not mutation of an already loaded network.

use super::invalid;
use crate::acoustic_realize::AcousticRealizeError;
use fs_aeroac::bessel::{j0, y0};
use fs_exec::CancelGate;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};

/// Explicit numerical resolution of the SAME continuum, not material data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CylinderFlowResolution {
    /// Smallest dimensionless diffusion wavenumber times physical radius.
    /// Finite in (0,0.5]; the excluded slow tail is integrated analytically.
    pub minimum_wavenumber_radius: f64,
    /// Largest dimensionless diffusion wavenumber times physical radius.
    /// Greater than the minimum and at most 3000 (the Bessel owner's domain).
    pub maximum_wavenumber_radius: f64,
    /// Positive midpoint intervals in log(wavenumber), in 1..=4096.
    /// The consumer must additionally admit the COMPLETE resulting state count.
    pub intervals: usize,
}
impl CylinderFlowResolution {
    /// Validate declared resolution without silently clipping any axis.
    pub fn validate(self) -> Result<(), AcousticRealizeError> {
        if !self.minimum_wavenumber_radius.is_finite()
            || self.minimum_wavenumber_radius <= 0.0
            || self.minimum_wavenumber_radius > 0.5
            || !self.maximum_wavenumber_radius.is_finite()
            || self.maximum_wavenumber_radius <= self.minimum_wavenumber_radius
            || self.maximum_wavenumber_radius > 3000.0
            || !(1..=4096).contains(&self.intervals)
        {
            return Err(invalid("cylinder flow needs 0 < X_MIN <= 0.5, X_MIN < X_MAX <= 3000 and 1..=4096 intervals"));
        }
        Ok(())
    }
}

/// One positive quadrature branch of the original fluid continuum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CylinderFlowBranch {
    /// Exponential fluid relaxation rate [1/s].
    pub rate_per_s: f64,
    /// Resistance per cylinder length [N s/m²].
    pub resistance_per_length_n_s_m2: f64,
    /// Actual branch inertia per cylinder length [kg/m], c/rate.
    pub mass_per_length_kg_m: f64,
}

/// Source-bound fluid loading before any mechanical normalization.
/// The default source keeps existing gas consumers unchanged; resolved liquid,
/// gas-mixture and moist-air owners retain their own complete query receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct CylinderFlow<F = GasState> {
    fluid: F,
    transport: FluidTransportState,
    radius_m: f64,
    resolution: CylinderFlowResolution,
    added_mass_per_length_kg_m: f64,
    low_tail_resistance_per_length_n_s_m2: f64,
    minimum_rate_per_s: f64,
    branches: Vec<CylinderFlowBranch>,
}
impl CylinderFlow<GasState> {
    /// Derive potential-flow inertia and every viscous memory coefficient from
    /// the SAME gas and actual external circular radius. Cold construction polls
    /// cancellation between quadrature cells. No pressure or viscosity default.
    pub fn from_gas(
        radius_m: f64,
        gas: GasState,
        resolution: CylinderFlowResolution,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        Self::from_fluid(radius_m, gas, resolution, gate)
    }

    /// Original coherent material gas state.
    #[must_use] pub const fn gas(&self) -> &GasState { &self.fluid }
}
impl<F: FluidTransport> CylinderFlow<F> {
    /// Prepare the original Stokes operator from one existing fluid source.
    ///
    /// `ResolvedLiquidState`, resolved gases/mixtures and explicit transport-only
    /// tuples all use this path. The checked snapshot is obtained exactly once;
    /// the original source, including any material/query receipts, remains on
    /// [`Self::fluid`]. A raw tuple does not gain material-card authority.
    /// Only density and dynamic viscosity enter this incompressible operator.
    /// Heat capacity/conductivity stay with the shared transport state; they do
    /// not invent thermal damping, sound speed, phase stability or cavitation.
    /// The caller must admit constant Newtonian transport over the motion.
    pub fn from_fluid(
        radius_m: f64,
        fluid: F,
        resolution: CylinderFlowResolution,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
        resolution.validate()?;
        if !radius_m.is_finite() || radius_m <= 0.0 {
            return Err(invalid("cylinder flow requires positive finite radius"));
        }
        let transport = fluid.transport_state().map_err(|error|
            AcousticRealizeError::Ambient(format!("cylinder fluid transport: {error}")))?;
        if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
        let density = transport.density_kg_m3();
        let viscosity = transport.dynamic_viscosity_pa_s();
        // Retain the original evaluation order for existing gas consumers.
        let diffusion_rate = viscosity / density / radius_m / radius_m;
        let added_mass = core::f64::consts::PI * radius_m * radius_m * density;
        let x_min = resolution.minimum_wavenumber_radius;
        let minimum_rate = diffusion_rate * x_min * x_min;
        // d arg(J0+iY0)/dx = 2/[pi*x*(J0²+Y0²)]. The x->0 phase is -pi/2.
        let low_tail = 8.0 * viscosity * det::atan2(j0(x_min), -y0(x_min));
        let log_min = det::ln(x_min);
        let step = (det::ln(resolution.maximum_wavenumber_radius) - log_min)
            / resolution.intervals as f64;
        if [diffusion_rate, added_mass, minimum_rate, low_tail, step].iter()
            .any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(invalid("cylinder fluid inertia, tail or diffusion grid is unrepresentable"));
        }
        let mut branches = Vec::with_capacity(resolution.intervals);
        for cell in 0..resolution.intervals {
            if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
            let x = det::exp(log_min + (cell as f64 + 0.5) * step);
            let denominator = j0(x).powi(2) + y0(x).powi(2);
            // d log(lambda) = 2 d log(x).
            let resistance = (16.0 / core::f64::consts::PI) * viscosity
                * step / denominator;
            let rate = diffusion_rate * x * x;
            let mass = resistance / rate;
            if [x, denominator, resistance, rate, mass].iter()
                .any(|v| !v.is_finite() || *v <= 0.0) {
                return Err(invalid("cylinder fluid memory branch is not finite and passive"));
            }
            branches.push(CylinderFlowBranch {
                rate_per_s: rate,
                resistance_per_length_n_s_m2: resistance,
                mass_per_length_kg_m: mass,
            });
        }
        if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
        Ok(Self { fluid, transport, radius_m, resolution, added_mass_per_length_kg_m: added_mass,
            low_tail_resistance_per_length_n_s_m2: low_tail,
            minimum_rate_per_s: minimum_rate, branches })
    }
}
impl<F> CylinderFlow<F> {
    /// Original input owner, with its source identity and complete conditions.
    #[must_use] pub const fn fluid(&self) -> &F { &self.fluid }
    /// Immutable checked coefficients used by every branch of this operator.
    #[must_use] pub const fn transport(&self) -> &FluidTransportState { &self.transport }
    /// Original physical external radius [m].
    #[must_use] pub const fn radius_m(&self) -> f64 { self.radius_m }
    /// Original three resolution axes.
    #[must_use] pub const fn resolution(&self) -> CylinderFlowResolution { self.resolution }
    /// Instantaneous potential-flow inertia [kg/m], excluding memory masses.
    #[must_use] pub const fn added_mass_per_length_kg_m(&self) -> f64 { self.added_mass_per_length_kg_m }
    /// Collapsed slow spectral tail's exact resistance [N s/m²].
    #[must_use] pub const fn low_tail_resistance_per_length_n_s_m2(&self) -> f64 {
        self.low_tail_resistance_per_length_n_s_m2
    }
    /// Slow cutoff [1/s], for the explicit low-tail error estimate.
    #[must_use] pub const fn minimum_rate_per_s(&self) -> f64 { self.minimum_rate_per_s }
    /// Positive mechanical branches in original quadrature order.
    #[must_use] pub fn branches(&self) -> &[CylinderFlowBranch] { &self.branches }

    /// Discretized force/velocity impedance [N s/m²], convention exp(+i omega t).
    /// Evaluates the actual passive operator; this is not an error certificate.
    pub fn impedance_per_length(&self, omega_rad_s: f64) -> Result<C64, AcousticRealizeError> {
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 {
            return Err(invalid("cylinder fluid impedance requires positive finite angular frequency"));
        }
        let mut z = C64 { re: self.low_tail_resistance_per_length_n_s_m2,
            im: self.added_mass_per_length_kg_m * omega_rad_s };
        for branch in &self.branches {
            // Scale before squaring to avoid overflow in rate²+omega².
            let scale = omega_rad_s.max(branch.rate_per_s);
            let w = omega_rad_s / scale;
            let rate = branch.rate_per_s / scale;
            let denominator = w*w + rate*rate;
            z.re += branch.resistance_per_length_n_s_m2 * w*w / denominator;
            z.im += branch.resistance_per_length_n_s_m2 * w*rate / denominator;
        }
        if !z.re.is_finite() || !z.im.is_finite() {
            return Err(invalid("cylinder fluid impedance is unrepresentable"));
        }
        Ok(z)
    }
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod transport_tests;
