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

/// Source-defined elastic frames with complete exterior-fluid inertia and memory.
pub mod frame;

/// The same Newtonian exterior with rotational no-slip boundary conditions.
pub mod rotary;

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

/// A translating cylindrical body and its original fluid-memory coordinates.
///
/// This borrows the exact prepared fluid/source owner and supplies a native pHS
/// component, not another integrator. State order is `[x, p, p_0, ...]`: physical
/// displacement [m], momentum of the solid PLUS instantaneous potential-flow
/// inertia [kg m/s], then each spectral fluid branch momentum [kg m/s].
/// The single port is applied force [N] / body velocity [m/s].
///
/// Uniform transverse translation along `length_m` uses the same long-cylinder
/// local-flow approximation as the impedance. End effects, confinement,
/// nonlinear drag, radiation, thermal feedback and intrinsic solid loss are not
/// added. The spring is an explicit external linear support, not an inferred
/// material law. Retain every fluid state when composing, stepping or resuming.
pub struct CylinderBody<'a, F = GasState> {
    flow: &'a CylinderFlow<F>,
    length_m: f64,
    dry_mass_kg: f64,
    moving_mass_kg: f64,
    memory_masses_kg: Vec<f64>,
    system: fs_phs::PortHamiltonian,
}

impl<F> CylinderFlow<F> {
    /// Couple the SAME fluid operator to physical transverse body motion.
    ///
    /// `dry_mass_kg` excludes ALL fluid inertia. Added potential-flow mass is
    /// included once; each remaining spectral mass is an independent dynamical
    /// coordinate linked by its original dashpot. The low tail acts against the
    /// stationary exterior. No old sampled fluid damping may be added again.
    ///
    /// `max_states` bounds the COMPLETE body and memory dimension before dense
    /// matrix allocation. Native pHS admission/solvers own their workspace and
    /// time-discretization controls. Construction polls around native admission
    /// and between branches, not within the native dense eigensolver.
    pub fn transverse_body(
        &self,
        length_m: f64,
        dry_mass_kg: f64,
        support_stiffness_n_m: f64,
        max_states: usize,
        gate: &CancelGate,
    ) -> Result<CylinderBody<'_, F>, AcousticRealizeError> {
        if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
        let n = self.branches.len().checked_add(2)
            .ok_or_else(|| invalid("cylinder body state count overflow"))?;
        if n > max_states {
            return Err(invalid("cylinder body and ALL fluid memory exceed the state budget"));
        }
        let square = n.checked_mul(n)
            .ok_or_else(|| invalid("cylinder body matrix dimension overflow"))?;
        let added_mass = length_m * self.added_mass_per_length_kg_m;
        let moving_mass = dry_mass_kg + added_mass;
        let low_resistance = length_m * self.low_tail_resistance_per_length_n_s_m2;
        if [length_m, dry_mass_kg, added_mass, moving_mass, 1.0 / moving_mass, low_resistance]
            .iter().any(|v| !v.is_finite() || *v <= 0.0)
            || moving_mass <= dry_mass_kg
            || !support_stiffness_n_m.is_finite() || support_stiffness_n_m < 0.0 {
            return Err(invalid("cylinder body requires representable positive geometry/inertia and nonnegative support stiffness"));
        }
        let mut q = vec![0.0; square];
        let mut j = vec![0.0; square];
        let mut r = vec![0.0; square];
        let mut g = vec![0.0; n];
        let mut memory_masses = Vec::with_capacity(self.branches.len());
        q[0] = support_stiffness_n_m;
        q[n + 1] = 1.0 / moving_mass;
        j[1] = 1.0;
        j[n] = -1.0;
        g[1] = 1.0;
        r[n + 1] = low_resistance;
        for (branch_index, branch) in self.branches.iter().enumerate() {
            if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
            let index = branch_index + 2;
            let mass = length_m * branch.mass_per_length_kg_m;
            let resistance = length_m * branch.resistance_per_length_n_s_m2;
            if [mass, 1.0 / mass, resistance].iter().any(|v| !v.is_finite() || *v <= 0.0) {
                return Err(invalid("physical cylinder fluid mass or resistance is unrepresentable"));
            }
            memory_masses.push(mass);
            q[index * n + index] = 1.0 / mass;
            // Each c*(v-v_j)^2 is one PSD resistance block. Its equal/opposite
            // forces advance both momenta; memory storage is never counted as heat.
            r[n + 1] += resistance;
            r[n + index] = -resistance;
            r[index * n + 1] = -resistance;
            r[index * n + index] = resistance;
        }
        if !r[n + 1].is_finite() {
            return Err(invalid("summed cylinder fluid resistance is unrepresentable"));
        }
        let storage = fs_phs::QuadraticStorage::new(q, n)
            .map_err(AcousticRealizeError::AlgebraicSolve)?;
        if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
        let system = fs_phs::PortHamiltonian::new(n, 1, j, r, g, Box::new(storage))
            .map_err(AcousticRealizeError::AlgebraicSolve)?;
        if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
        Ok(CylinderBody { flow: self, length_m, dry_mass_kg,
            moving_mass_kg: moving_mass, memory_masses_kg: memory_masses, system })
    }
}

impl<F> CylinderBody<'_, F> {
    /// Original typed source and unchanged quadrature/resolution.
    #[must_use] pub const fn flow(&self) -> &CylinderFlow<F> { self.flow }
    /// Uniformly immersed cylinder length [m].
    #[must_use] pub const fn length_m(&self) -> f64 { self.length_m }
    /// Original solid/body mass, excluding fluid [kg].
    #[must_use] pub const fn dry_mass_kg(&self) -> f64 { self.dry_mass_kg }
    /// Solid plus instantaneous potential-flow inertia, excluding memory [kg].
    #[must_use] pub const fn moving_mass_kg(&self) -> f64 { self.moving_mass_kg }
    /// Independent spectral masses, in the original quadrature order [kg].
    #[must_use] pub fn memory_masses_kg(&self) -> &[f64] { &self.memory_masses_kg }
    /// Native force/velocity component. Existing pHS steppers and composers
    /// advance this COMPLETE state; there is no private fluid time stepper.
    #[must_use] pub const fn system(&self) -> &fs_phs::PortHamiltonian { &self.system }

    /// Convert physical initial motion to the native conserved state order.
    /// Fluid velocities must be explicit, including a declared quiescent start;
    /// this method does not erase an evolving or resumed fluid history.
    pub fn state_from_motion(
        &self,
        displacement_m: f64,
        velocity_m_s: f64,
        fluid_velocities_m_s: &[f64],
    ) -> Result<Vec<f64>, AcousticRealizeError> {
        if !displacement_m.is_finite() || !velocity_m_s.is_finite()
            || fluid_velocities_m_s.len() != self.memory_masses_kg.len()
            || fluid_velocities_m_s.iter().any(|v| !v.is_finite()) {
            return Err(invalid("cylinder initial motion requires every finite physical fluid velocity"));
        }
        let momentum = |mass: f64, velocity: f64| {
            let value = mass * velocity;
            if !value.is_finite() || velocity != 0.0 && value == 0.0 {
                Err(invalid("cylinder momentum is unrepresentable"))
            } else { Ok(value) }
        };
        let mut state = Vec::with_capacity(self.system.state_dim());
        state.push(displacement_m);
        state.push(momentum(self.moving_mass_kg, velocity_m_s)?);
        for (&mass, &velocity) in self.memory_masses_kg.iter().zip(fluid_velocities_m_s) {
            state.push(momentum(mass, velocity)?);
        }
        if !self.system.hamiltonian(&state).is_finite() {
            return Err(invalid("cylinder initial stored energy is unrepresentable"));
        }
        Ok(state)
    }
}

#[cfg(test)]
mod mechanical_tests {
    use super::*;

    fn flow(density: f64, viscosity: f64) -> CylinderFlow<FluidTransportState> {
        CylinderFlow::from_fluid(0.2,
            FluidTransportState::try_new(density, 4000.0, 0.5, viscosity).unwrap(),
            CylinderFlowResolution { minimum_wavenumber_radius: 0.01,
                maximum_wavenumber_radius: 2.0, intervals: 4 },
            &CancelGate::new_clock_free()).unwrap()
    }
    fn near(a: f64, b: f64) {
        assert!((a-b).abs() < 1e-10 * (1.0 + a.abs() + b.abs()), "{a:e} != {b:e}");
    }
    fn rhs(body: &CylinderBody<'_, FluidTransportState>, x: &[f64], force: f64) -> Vec<f64> {
        let e = body.system().effort(x);
        let (j, r, g) = body.system().structure();
        let n = x.len();
        (0..n).map(|i| (0..n).map(|k| (j[i*n+k]-r[i*n+k])*e[k]).sum::<f64>()
            + g[i]*force).collect()
    }

    #[test]
    fn native_component_has_original_reciprocal_forces_energy_and_heat() {
        let flow = flow(2.0, 0.1);
        let body = flow.transverse_body(0.7, 0.6, 4.0, 6, &CancelGate::new_clock_free()).unwrap();
        let fluid_v = [0.1, -0.2, 0.25, -0.1];
        let (x, v, force) = (0.02, 0.3, 0.8);
        let state = body.state_from_motion(x, v, &fluid_v).unwrap();
        let dx = rhs(&body, &state, force);
        near(body.moving_mass_kg(), 0.6 + 0.7 * core::f64::consts::PI * 2.0 * 0.2 * 0.2);
        near(dx[0], v);
        near(body.system().output(&state)[0], v);
        let mut body_force = force - 4.0*x;
        let low = body.length_m() * flow.low_tail_resistance_per_length_n_s_m2();
        body_force -= low*v;
        let mut heat = low*v*v;
        let mut energy = 0.5 * 4.0*x*x + 0.5 * body.moving_mass_kg()*v*v;
        for (k, branch) in flow.branches().iter().enumerate() {
            let c = 0.7 * branch.resistance_per_length_n_s_m2;
            let transfer = c * (v - fluid_v[k]);
            near(dx[k+2], transfer);
            body_force -= transfer;
            heat += c * (v-fluid_v[k]).powi(2);
            energy += 0.5 * body.memory_masses_kg()[k] * fluid_v[k].powi(2);
        }
        near(dx[1], body_force);
        near(dx[1..].iter().sum(), force - 4.0*x - low*v);
        near(body.system().hamiltonian(&state), energy);
        let power: f64 = body.system().effort(&state).iter().zip(&dx).map(|(e,d)| e*d).sum();
        near(power, force*v - heat);
        assert!(heat > 0.0);
        assert!(core::ptr::eq(body.flow(), &flow));
    }

    #[test]
    fn actual_native_step_preserves_fluid_history_and_closes_physical_work() {
        let flow = flow(2.0, 0.1);
        let body = flow.transverse_body(0.7, 0.6, 4.0, 6, &CancelGate::new_clock_free()).unwrap();
        let mut x = body.state_from_motion(0.02, 0.3, &[0.1, -0.2, 0.25, -0.1]).unwrap();
        let mut cumulative_heat = 0.0;
        let initial = body.system().hamiltonian(&x);
        let mut cumulative_work = 0.0;
        for _ in 0..50 {
            let (force, dt) = (0.03, 0.002);
            let record = fs_phs::step(body.system(), &x, &[force], dt).unwrap();
            let mid: Vec<_> = x.iter().zip(&record.x).map(|(a,b)| 0.5*(a+b)).collect();
            let v = mid[1] / body.moving_mass_kg();
            let mut heat = flow.low_tail_resistance_per_length_n_s_m2() * 0.7 * v*v;
            for (k, branch) in flow.branches().iter().enumerate() {
                let fluid_v = mid[k+2] / body.memory_masses_kg()[k];
                heat += 0.7 * branch.resistance_per_length_n_s_m2 * (v-fluid_v).powi(2);
            }
            let expected_rhs = rhs(&body, &mid, force);
            for k in 0..x.len() { near((record.x[k]-x[k])/dt, expected_rhs[k]); }
            cumulative_heat += dt*heat;
            cumulative_work += force*(record.x[0]-x[0]);
            x = record.x;
        }
        assert!(cumulative_heat > 0.0);
        near(body.system().hamiltonian(&x)-initial+cumulative_heat, cumulative_work);
        assert!(x[2..].iter().any(|p| *p != 0.0));
    }

    #[test]
    fn harmonic_elimination_recovers_same_complex_cylinder_impedance() {
        let flow = flow(2.0, 0.1);
        let body = flow.transverse_body(0.7, 0.6, 4.0, 6, &CancelGate::new_clock_free()).unwrap();
        let (_, r, _) = body.system().structure();
        let n = body.system().state_dim();
        for omega in [0.1, 1.0, 10.0] {
            // Independently eliminate every real fluid momentum equation at
            // unit harmonic body velocity. This includes both coupling signs.
            let mut real = r[n+1];
            let mut imag = omega*body.moving_mass_kg() - 4.0/omega;
            for k in 2..n {
                let c = r[k*n+k];
                let wm = omega*body.memory_masses_kg()[k-2];
                let denominator = c*c + wm*wm;
                let product = r[n+k]*r[k*n+1];
                real -= product*c/denominator;
                imag += product*wm/denominator;
            }
            let z = flow.impedance_per_length(omega).unwrap();
            near(real, 0.7*z.re);
            near(imag, omega*0.6 - 4.0/omega + 0.7*z.im);
        }
    }

    #[test]
    fn source_data_changes_native_motion_not_only_the_impedance_diagnostic() {
        let a = flow(2.0, 0.1);
        let b = flow(8.0, 0.1);
        let c = flow(2.0, 0.4);
        let gate = CancelGate::new_clock_free();
        let a = a.transverse_body(0.7, 0.6, 0.0, 6, &gate).unwrap();
        let b = b.transverse_body(0.7, 0.6, 0.0, 6, &gate).unwrap();
        let c = c.transverse_body(0.7, 0.6, 0.0, 6, &gate).unwrap();
        let zero = [0.0; 6];
        near(rhs(&a, &zero, 1.0)[1] / a.moving_mass_kg(), 1.0 / a.moving_mass_kg());
        assert!(b.moving_mass_kg() > a.moving_mass_kg());
        near(c.moving_mass_kg(), a.moving_mass_kg());
        let ra = fs_phs::step(a.system(), &zero, &[1.0], 0.01).unwrap();
        let rb = fs_phs::step(b.system(), &zero, &[1.0], 0.01).unwrap();
        let rc = fs_phs::step(c.system(), &zero, &[1.0], 0.01).unwrap();
        assert!(rb.x[0] < ra.x[0]);
        assert!(rc.x[0] < ra.x[0]);
    }

    #[test]
    fn complete_capacity_and_nonphysical_initial_motion_refuse() {
        let flow = flow(2.0, 0.1);
        let gate = CancelGate::new_clock_free();
        assert!(flow.transverse_body(0.7, 0.6, 4.0, 5, &gate).is_err());
        for value in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
            assert!(flow.transverse_body(value, 0.6, 4.0, 6, &gate).is_err());
            assert!(flow.transverse_body(0.7, value, 4.0, 6, &gate).is_err());
        }
        assert!(flow.transverse_body(0.7, 0.6, -1.0, 6, &gate).is_err());
        let body = flow.transverse_body(0.7, 0.6, 0.0, 6, &gate).unwrap();
        assert!(body.state_from_motion(0.0, 0.0, &[0.0; 3]).is_err());
        assert!(body.state_from_motion(0.0, f64::NAN, &[0.0; 4]).is_err());
        assert!(body.state_from_motion(0.0, 0.0, &[f64::INFINITY; 4]).is_err());
        assert_eq!(body.state_from_motion(0.0, 0.0, &[0.0; 4]).unwrap(), vec![0.0; 6]);
    }
}
