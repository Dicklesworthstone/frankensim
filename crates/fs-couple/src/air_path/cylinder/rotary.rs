//! Causal torque from the rotational Stokes cylinder continuum.
//!
//! Angular motion has a different boundary condition from transverse motion:
//! rho*u_t = mu*(u_rr + u_r/r - u/r^2), u(a)=a*Omega, u(infinity)=0.
//! Its Laplace solution is u(r)=a*Omega*K1(q*r)/K1(q*a), q=sqrt(rho*s/mu).
//! Integrating mu*(u_r-u/r) over the circumference gives resisting torque per
//! length Z(s)*Omega, with Z(s)=2*pi*mu*a^2*(2+z*K0(z)/K1(z)), z=q*a.
//! The exact steady resistance is 4*pi*mu*a^2; there is no instantaneous
//! potential-flow rotary inertia. All remaining inertia belongs to the fluid
//! diffusion history, never to the solid material's polar moment.
//!
//! The Bessel Wronskian (NIST DLMF 10.5.2) gives positive spectral weight
//! 4*mu*a^2/[pi*(J1(x)^2+Y1(x)^2)] per d(log lambda), x=a*sqrt(rho*lambda/mu).
//! Each quadrature cell is therefore an actual inertia/dashpot branch,
//! Z_j(s)=c_j*s/(s+lambda_j), I_j=c_j/lambda_j. The slow tail is integrated
//! by the phase of J1+iY1 and collapsed to resistance, with impedance error
//! at most R_low*lambda_min/omega. Refine cutoffs and grid independently.
//! <https://dlmf.nist.gov/10.5.E2>; the independent modified-Bessel integral
//! oracle is <https://dlmf.nist.gov/10.32.E9>.
//!
//! This is infinitesimal, unconfined, incompressible circular-cylinder rotation
//! in a frozen Newtonian fluid. It does not add end effects, nonlinear swirl,
//! neighboring-string flow, compressible radiation or thermal feedback.

use super::{
    AcousticRealizeError, C64, CancelGate, CylinderFlowResolution, FluidTransport,
    FluidTransportState, GasState, det, invalid,
};
use fs_aeroac::bessel::{j1, y1};

/// One positive rotational diffusion branch, before mechanical normalization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotaryCylinderBranch {
    /// Fluid relaxation rate [1/s].
    pub rate_per_s: f64,
    /// Torque/angular-velocity resistance per length [N m s / m = N s].
    pub resistance_per_length_n_s: f64,
    /// Polar fluid inertia per length [kg m^2 / m = kg m], c/rate.
    pub inertia_per_length_kg_m: f64,
}

/// The original fluid and its passive rotational continuum discretization.
#[derive(Clone, Debug, PartialEq)]
pub struct RotaryCylinderFlow<F = GasState> {
    fluid: F,
    transport: FluidTransportState,
    radius_m: f64,
    resolution: CylinderFlowResolution,
    steady_resistance_per_length_n_s: f64,
    low_tail_resistance_per_length_n_s: f64,
    minimum_rate_per_s: f64,
    branches: Vec<RotaryCylinderBranch>,
}
impl RotaryCylinderFlow<GasState> {
    /// Use the exact gas state already shared by the structure and acoustics.
    pub fn from_gas(radius_m: f64, gas: GasState, resolution: CylinderFlowResolution,
        gate: &CancelGate) -> Result<Self, AcousticRealizeError>
    {
        Self::from_fluid(radius_m, gas, resolution, gate)
    }
    #[must_use]
    pub const fn gas(&self) -> &GasState { &self.fluid }
}
impl<F: FluidTransport> RotaryCylinderFlow<F> {
    /// Compile the same circular diffusion equation from checked fluid data.
    /// The complete original owner/receipts survive. No density, viscosity,
    /// rotary material loss or resolution default is supplied.
    pub fn from_fluid(radius_m: f64, fluid: F, resolution: CylinderFlowResolution,
        gate: &CancelGate) -> Result<Self, AcousticRealizeError>
    {
        if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
        resolution.validate()?;
        if !radius_m.is_finite() || radius_m <= 0.0 {
            return Err(invalid("rotary cylinder flow requires positive finite radius"));
        }
        let transport = fluid.transport_state().map_err(|error|
            AcousticRealizeError::Ambient(format!("rotary cylinder fluid transport: {error}")))?;
        let viscosity = transport.dynamic_viscosity_pa_s();
        let diffusion = viscosity / transport.density_kg_m3() / radius_m / radius_m;
        let scale = viscosity * radius_m * radius_m;
        let steady = 4.0 * core::f64::consts::PI * scale;
        let x_min = resolution.minimum_wavenumber_radius;
        let minimum_rate = diffusion * x_min * x_min;
        let low_tail = 4.0 * scale * det::atan2(j1(x_min), -y1(x_min));
        let log_min = det::ln(x_min);
        let step = (det::ln(resolution.maximum_wavenumber_radius) - log_min)
            / resolution.intervals as f64;
        if [diffusion, scale, steady, minimum_rate, low_tail, step].iter()
            .any(|value| !value.is_finite() || *value <= 0.0) {
            return Err(invalid("rotary cylinder transport, tail or diffusion grid is unrepresentable"));
        }
        let mut branches = Vec::with_capacity(resolution.intervals);
        for cell in 0..resolution.intervals {
            if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
            let x = det::exp(log_min + (cell as f64 + 0.5) * step);
            let denominator = j1(x).powi(2) + y1(x).powi(2);
            // d(log lambda) = 2*d(log x); this is the order-one rotary spectrum.
            let resistance = (8.0 / core::f64::consts::PI) * scale * step / denominator;
            let rate = diffusion * x * x;
            let inertia = resistance / rate;
            if [x, denominator, resistance, rate, inertia].iter()
                .any(|value| !value.is_finite() || *value <= 0.0) {
                return Err(invalid("rotary cylinder memory is not finite and passive"));
            }
            branches.push(RotaryCylinderBranch {
                rate_per_s: rate, resistance_per_length_n_s: resistance,
                inertia_per_length_kg_m: inertia,
            });
        }
        if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
        Ok(Self { fluid, transport, radius_m, resolution,
            steady_resistance_per_length_n_s: steady,
            low_tail_resistance_per_length_n_s: low_tail,
            minimum_rate_per_s: minimum_rate, branches })
    }
}
impl<F> RotaryCylinderFlow<F> {
    #[must_use] pub const fn fluid(&self) -> &F { &self.fluid }
    #[must_use] pub const fn transport(&self) -> &FluidTransportState { &self.transport }
    #[must_use] pub const fn radius_m(&self) -> f64 { self.radius_m }
    #[must_use] pub const fn resolution(&self) -> CylinderFlowResolution { self.resolution }
    /// Exact stationary-fluid resistance, independent of fluid density [N s].
    #[must_use] pub const fn steady_resistance_per_length_n_s(&self) -> f64 {
        self.steady_resistance_per_length_n_s
    }
    /// Collapsed slow-tail weight, excluding the separate exact steady term [N s].
    #[must_use] pub const fn low_tail_resistance_per_length_n_s(&self) -> f64 {
        self.low_tail_resistance_per_length_n_s
    }
    #[must_use] pub const fn minimum_rate_per_s(&self) -> f64 { self.minimum_rate_per_s }
    #[must_use] pub fn branches(&self) -> &[RotaryCylinderBranch] { &self.branches }

    /// Actual discrete torque/angular-velocity impedance per length [N s],
    /// exp(+i*omega*t). This evaluates the prepared approximation, not an error
    /// certificate; its nonzero low cutoff is not the exact DC limit.
    pub fn impedance_per_length(&self, omega_rad_s: f64) -> Result<C64, AcousticRealizeError> {
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 {
            return Err(invalid("rotary cylinder impedance needs positive finite angular frequency"));
        }
        let mut z = C64 { re: self.steady_resistance_per_length_n_s
            + self.low_tail_resistance_per_length_n_s, im: 0.0 };
        for branch in &self.branches {
            let scale = omega_rad_s.max(branch.rate_per_s);
            let w = omega_rad_s / scale;
            let rate = branch.rate_per_s / scale;
            let denominator = w*w + rate*rate;
            z.re += branch.resistance_per_length_n_s * w*w / denominator;
            z.im += branch.resistance_per_length_n_s * w*rate / denominator;
        }
        if !z.re.is_finite() || !z.im.is_finite() {
            return Err(invalid("rotary cylinder impedance is unrepresentable"));
        }
        Ok(z)
    }
}

#[cfg(test)]
mod tests;
