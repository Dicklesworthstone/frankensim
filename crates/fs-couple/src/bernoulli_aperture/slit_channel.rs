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
use fs_math::det;

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

/// Stateful transverse momentum field with an allocation-free implicit step.
#[derive(Debug, Clone)]
pub struct SlitChannel {
    gas: GasState,
    spec: SlitChannelSpec,
    cell_mass_kg: f64,
    viscous_factor: f64,
    area_m2: f64,
    inverse_diagonal: Vec<f64>,
    upper_ratio: Vec<f64>,
    pressure_response: Vec<f64>,
    mean_pressure_response: f64,
    velocity_m_s: Vec<f64>,
    trial: Vec<f64>,
    last: SlitChannelStep,
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

// Thomas elimination for the factored symmetric diffusion matrix.
// Both arrays are cold allocations; no state is shared between channels.
fn solve(inverse: &[f64], upper: &[f64], rhs: &mut [f64]) {
    rhs[0] *= inverse[0];
    for i in 1..rhs.len() {
        rhs[i] = rhs[i] * inverse[i] + upper[i] * rhs[i - 1];
    }
    for i in (0..rhs.len() - 1).rev() {
        rhs[i] += upper[i] * rhs[i + 1];
    }
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
        let r = gas.dynamic_viscosity / gas.density * spec.sample_period_s / (dy * dy);
        let pressure_scale = spec.sample_period_s / (gas.density * spec.length_m);
        let area_m2 = spec.width_m * spec.height_m;
        let cell_mass_kg = gas.density * spec.length_m * spec.width_m * dy;
        let viscous_factor = gas.dynamic_viscosity * spec.length_m * spec.width_m / dy;
        if ![dy, r, pressure_scale, area_m2, cell_mass_kg, viscous_factor]
            .iter().all(|x| x.is_finite() && *x > 0.0)
        {
            return Err(invalid("unsteady slit coefficients cannot be represented"));
        }
        let mut inverse_diagonal = vec![0.0; n];
        let mut upper_ratio = vec![0.0; n];
        for i in 0..n {
            // A wall is half a cell from its centre: 2r at the wall plus r
            // at the interior face. Replacing 3r by 2r would move the wall.
            let diagonal = 1.0 + if i == 0 || i + 1 == n { 3.0 * r } else { 2.0 * r };
            let pivot = diagonal - if i == 0 { 0.0 } else { r * upper_ratio[i - 1] };
            inverse_diagonal[i] = 1.0 / pivot;
            upper_ratio[i] = r * inverse_diagonal[i];
            if !pivot.is_finite() || pivot <= 0.0 || !upper_ratio[i].is_finite()
                || !inverse_diagonal[i].is_finite() || inverse_diagonal[i] <= 0.0
            {
                return Err(invalid("unsteady slit diffusion factorization is not finite positive definite"));
            }
        }
        let mut pressure_response = vec![pressure_scale; n];
        solve(&inverse_diagonal, &upper_ratio, &mut pressure_response);
        let mean_pressure_response = pressure_response.iter().sum::<f64>() / n as f64;
        if !pressure_response.iter().all(|v| v.is_finite() && *v > 0.0)
            || !mean_pressure_response.is_finite() || mean_pressure_response <= 0.0
        {
            return Err(invalid("unsteady slit pressure response cannot be represented"));
        }
        Ok(Self { gas: *gas, spec, cell_mass_kg, viscous_factor, area_m2,
            inverse_diagonal, upper_ratio, pressure_response, mean_pressure_response,
            velocity_m_s: vec![0.0; n], trial: vec![0.0; n], last: SlitChannelStep::default() })
    }

    /// Unchanged source material driving storage and wall friction.
    #[must_use]
    pub const fn gas(&self) -> &GasState { &self.gas }

    /// Original geometry, mesh and time step.
    #[must_use]
    pub const fn spec(&self) -> &SlitChannelSpec { &self.spec }

    /// Actual cell-centre velocities from the lower to the upper wall [m/s].
    #[must_use]
    pub fn velocities_m_s(&self) -> &[f64] { &self.velocity_m_s }

    /// Last accepted exchange, all zero before the first step.
    #[must_use]
    pub const fn last_step(&self) -> SlitChannelStep { self.last }

    /// Exact steady bulk response of THIS spatial discretization [m/s].
    /// It is useful for cold admission, but never replaces the live field.
    pub fn steady_bulk_speed_m_s(&self, pressure_pa: f64) -> Result<f64, AcousticRealizeError> {
        if !pressure_pa.is_finite() { return Err(invalid("unsteady slit pressure must be finite")); }
        // Cell-centred no-slip FV conductance differs from continuum Poiseuille
        // by 1+2/N². Retain it so the admitted pressure agrees with the actual mesh.
        let n = self.spec.transverse_cells as f64;
        let resistance = 12.0 * self.gas.dynamic_viscosity * self.spec.length_m
            / (self.spec.height_m * self.spec.height_m * (1.0 + 2.0 / (n * n)));
        let speed = match self.spec.exit {
            SlitExit::StaticPressure => pressure_pa / resistance,
            SlitExit::UniformKineticHead => 2.0 * pressure_pa / (resistance
                + det::hypot(resistance, det::sqrt(2.0 * self.gas.density * pressure_pa.abs()))),
        };
        if !resistance.is_finite() || resistance <= 0.0 || !speed.is_finite() {
            return Err(invalid("unsteady slit steady pressure response cannot be represented"));
        }
        Ok(speed)
    }

    /// Advance one physical pressure interval, allowing signed pressure and flow.
    /// Failure leaves the accepted velocity field and last exchange unchanged.
    ///
    /// The terminal law makes the mean response a scalar monotone quadratic.
    /// Its analytic solution feeds the complete factored transverse field; no
    /// fitted flow coefficient, nonlinear iteration or one-sample lag appears.
    pub fn step(&mut self, pressure_pa: f64) -> Result<SlitChannelStep, AcousticRealizeError> {
        if !pressure_pa.is_finite() { return Err(invalid("unsteady slit pressure must be finite")); }
        self.trial.copy_from_slice(&self.velocity_m_s);
        solve(&self.inverse_diagonal, &self.upper_ratio, &mut self.trial);
        let base = self.trial.iter().sum::<f64>() / self.trial.len() as f64;
        let target = base + self.mean_pressure_response * pressure_pa;
        let speed = match self.spec.exit {
            SlitExit::StaticPressure => target,
            SlitExit::UniformKineticHead => 2.0 * target / (1.0 + det::sqrt(
                1.0 + 2.0 * self.mean_pressure_response * self.gas.density * target.abs())),
        };
        let exit_pressure = match self.spec.exit {
            SlitExit::StaticPressure => 0.0,
            SlitExit::UniformKineticHead => 0.5 * self.gas.density * speed * speed.abs(),
        };
        let channel_pressure = pressure_pa - exit_pressure;
        if ![target, speed, exit_pressure, channel_pressure].iter().all(|x| x.is_finite()) {
            return Err(invalid("unsteady slit implicit pressure/flow solve is not finite"));
        }
        let mut sum = 0.0;
        let mut energy_sum = 0.0;
        let mut delta_sum = 0.0;
        let mut numerical_sum = 0.0;
        let mut shear_sum = 0.0;
        for i in 0..self.trial.len() {
            let old = self.velocity_m_s[i];
            self.trial[i] += self.pressure_response[i] * channel_pressure;
            let v = self.trial[i];
            sum += v;
            energy_sum += v * v;
            delta_sum += (v - old) * (v + old);
            numerical_sum += (v - old) * (v - old);
            if i != 0 { shear_sum += (v - self.trial[i - 1]).powi(2); }
        }
        shear_sum += 2.0 * (self.trial[0].powi(2) + self.trial[self.trial.len() - 1].powi(2));
        let bulk = sum / self.trial.len() as f64;
        let flow = self.area_m2 * bulk;
        let dt = self.spec.sample_period_s;
        let mut step = SlitChannelStep {
            bulk_speed_m_s: bulk, flow_m3_s: flow,
            kinetic_energy_j: 0.5 * self.cell_mass_kg * energy_sum,
            stored_delta_j: 0.5 * self.cell_mass_kg * delta_sum,
            pressure_work_j: pressure_pa * flow * dt,
            viscous_dissipation_j: self.viscous_factor * shear_sum * dt,
            exit_work_j: exit_pressure * flow * dt,
            numerical_dissipation_j: 0.5 * self.cell_mass_kg * numerical_sum,
            energy_residual_j: 0.0,
        };
        step.energy_residual_j = step.pressure_work_j - step.stored_delta_j
            - step.viscous_dissipation_j - step.exit_work_j - step.numerical_dissipation_j;
        let scale = step.pressure_work_j.abs() + step.stored_delta_j.abs()
            + step.viscous_dissipation_j + step.exit_work_j + step.numerical_dissipation_j;
        if ![bulk, flow, step.kinetic_energy_j, step.stored_delta_j, step.pressure_work_j,
            step.viscous_dissipation_j, step.exit_work_j, step.numerical_dissipation_j,
            step.energy_residual_j, scale].iter().all(|x| x.is_finite())
            || step.exit_work_j < 0.0
            || step.energy_residual_j.abs() > 4096.0 * f64::EPSILON * scale
        {
            return Err(invalid("unsteady slit step does not close finite dissipation and kinetic storage"));
        }
        core::mem::swap(&mut self.velocity_m_s, &mut self.trial);
        self.last = step;
        Ok(step)
    }
}
