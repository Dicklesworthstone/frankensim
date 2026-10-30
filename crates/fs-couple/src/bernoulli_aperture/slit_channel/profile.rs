//! Shared finite-volume momentum, pressure condensation and work accounting.
//! Geometry contributes cell areas and no-slip face conductances, never a
//! separately fitted resistance. In cell i, M_i dv_i/dt = A_i dp - (Gv)_i.
//! Symmetric G has one nonnegative conductance per face; its quadratic form
//! is exactly the wall heat used below. The same owner serves planar and
//! axisymmetric passages, including pressure reversal and retained inertia.

use super::{AcousticRealizeError, ChannelGasUpdate, SlitChannelLoad, SlitChannelStep, SlitExit, invalid};
use fs_math::det;

#[derive(Debug, Clone)]
pub(super) struct MomentumProfile {
    density: f64,
    dt: f64,
    exit: SlitExit,
    area: Vec<f64>,
    total_area: f64,
    mass: Vec<f64>,
    conductance: Vec<f64>,
    inverse_diagonal: Vec<f64>,
    lower_ratio: Vec<f64>,
    upper_ratio: Vec<f64>,
    pressure_response: Vec<f64>,
    mean_pressure_response: f64,
    steady_mean_response: f64,
    steady_midpoint_response: f64,
    velocity: Vec<f64>,
    trial: Vec<f64>,
    energy_j: f64,
    last: SlitChannelStep,
}

// Cell centres straddle the geometric midpoint for an even mesh; an odd mesh
// has a centre exactly there. This is linear reconstruction of the original
// field, not a parabolic profile or a rescaling of its bulk velocity.
fn midpoint_sample(values: &[f64]) -> f64 {
    let upper = values.len() / 2;
    if values.len() % 2 == 1 { values[upper] }
    else { 0.5 * values[upper - 1] + 0.5 * values[upper] }
}

fn factor(mass: &[f64], conductance: &[f64], dt: f64)
    -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), AcousticRealizeError>
{
    let mut inverse = vec![0.0; mass.len()];
    let mut lower = vec![0.0; mass.len()];
    let mut upper = vec![0.0; mass.len()];
    for i in 0..mass.len() {
        let down = dt * conductance[i];
        let up = dt * conductance[i + 1];
        let diagonal = mass[i] + down + up;
        let pivot = diagonal - if i == 0 { 0.0 } else { down * upper[i - 1] };
        inverse[i] = 1.0 / pivot;
        lower[i] = down * inverse[i];
        upper[i] = up * inverse[i];
        if ![pivot, inverse[i], lower[i], upper[i]].iter().all(|v| v.is_finite())
            || pivot <= 0.0 || inverse[i] <= 0.0 {
            return Err(invalid("channel momentum factorization is not finite positive definite"));
        }
    }
    Ok((inverse, lower, upper))
}

fn solve(inverse: &[f64], lower: &[f64], upper: &[f64], rhs: &mut [f64]) {
    rhs[0] *= inverse[0];
    for i in 1..rhs.len() { rhs[i] = rhs[i] * inverse[i] + lower[i] * rhs[i - 1]; }
    for i in (0..rhs.len() - 1).rev() { rhs[i] += upper[i] * rhs[i + 1]; }
}

impl MomentumProfile {
    pub(super) fn new(density: f64, length: f64, dt: f64, exit: SlitExit,
        area: Vec<f64>, conductance: Vec<f64>) -> Result<Self, AcousticRealizeError>
    {
        let n = area.len();
        if n < 2 || n > super::MAX_TRANSVERSE_CELLS || conductance.len() != n + 1
            || ![density, length, dt].iter().all(|x| x.is_finite() && *x > 0.0)
            || area.iter().any(|x| !x.is_finite() || *x <= 0.0)
            || conductance.iter().any(|x| !x.is_finite() || *x < 0.0)
            || conductance[1..].iter().any(|x| *x <= 0.0)
        { return Err(invalid("channel momentum needs positive physical cells and connected no-slip faces")); }
        let total_area = area.iter().sum::<f64>();
        let mass: Vec<_> = area.iter().map(|a| density * length * a).collect();
        if !total_area.is_finite() || mass.iter().any(|m| !m.is_finite() || *m <= 0.0) {
            return Err(invalid("channel momentum geometry or inertial mass is not representable"));
        }
        let (inverse_diagonal, lower_ratio, upper_ratio) = factor(&mass, &conductance, dt)?;
        let mut pressure_response: Vec<_> = area.iter().map(|a| dt * a).collect();
        solve(&inverse_diagonal, &lower_ratio, &upper_ratio, &mut pressure_response);
        let mean_pressure_response = area.iter().zip(&pressure_response)
            .map(|(a, v)| a * v).sum::<f64>() / total_area;
        // Steady response of the actual mesh, using the SAME no-slip operator.
        let (steady_inverse, steady_lower, steady_upper) = factor(&vec![0.0; n], &conductance, 1.0)?;
        let mut steady = area.clone();
        solve(&steady_inverse, &steady_lower, &steady_upper, &mut steady);
        let steady_mean_response = area.iter().zip(&steady).map(|(a, v)| a * v).sum::<f64>() / total_area;
        let steady_midpoint_response = midpoint_sample(&steady);
        if ![mean_pressure_response, steady_mean_response, steady_midpoint_response]
            .iter().all(|x| x.is_finite() && *x > 0.0)
            || pressure_response.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(invalid("channel pressure response cannot be represented"));
        }
        Ok(Self { density, dt, exit, area, total_area, mass, conductance,
            inverse_diagonal, lower_ratio, upper_ratio, pressure_response,
            mean_pressure_response, steady_mean_response, steady_midpoint_response,
            velocity: vec![0.0; n], trial: vec![0.0; n], energy_j: 0.0,
            last: SlitChannelStep::default() })
    }

    pub(super) fn velocities(&self) -> &[f64] { &self.velocity }
    pub(super) fn midpoint_speed(&self) -> f64 { midpoint_sample(&self.velocity) }

    // The same no-slip operator supplies the steady shape and bulk response.
    // The nonlinear exit changes its amplitude, never this normalized shape.
    pub(super) fn steady_midpoint_speed(&self, pressure: f64)
        -> Result<f64, AcousticRealizeError>
    {
        let speed = self.steady_bulk_speed(pressure)?
            * (self.steady_midpoint_response / self.steady_mean_response);
        if !speed.is_finite() {
            return Err(invalid("channel steady midpoint velocity cannot be represented"));
        }
        Ok(speed)
    }
    pub(super) const fn last_step(&self) -> SlitChannelStep { self.last }
    pub(super) const fn stored_energy_j(&self) -> f64 { self.energy_j }

    // This profile has already been assembled with the new complete gas at
    // identical physical geometry. Preserve each cell's axial momentum rather
    // than its speed: M_new v_new = M_old v_old. The original common density
    // ratio avoids introducing cell-dependent rounding into a uniform change.
    pub(super) fn transfer_material_history(
        &mut self, previous: &Self, maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ChannelGasUpdate, AcousticRealizeError> {
        super::gas::poll(cancelled)?;
        if self.area != previous.area || self.dt.to_bits() != previous.dt.to_bits()
            || self.exit != previous.exit || self.mass.len() != previous.mass.len()
        {
            return Err(invalid("channel gas transfer requires the same geometry, cells, clock and exit"));
        }
        let ratio = previous.density / self.density;
        if !ratio.is_finite() || ratio <= 0.0 {
            return Err(invalid("channel gas transfer has unrepresentable density ratio"));
        }
        let mut twice_energy = 0.0;
        for i in 0..self.velocity.len() {
            if i % 64 == 0 { super::gas::poll(cancelled)?; }
            let old = previous.velocity[i];
            let velocity = old * ratio;
            if !velocity.is_finite() || (old != 0.0 && velocity == 0.0) {
                return Err(invalid("channel gas transfer lost finite cell momentum"));
            }
            self.velocity[i] = velocity;
            twice_energy += self.mass[i] * velocity * velocity;
        }
        let after = 0.5 * twice_energy;
        let work = after - previous.energy_j;
        if !after.is_finite() || !work.is_finite()
            || work.abs() > maximum_absolute_parameter_work_j
        {
            return Err(invalid("channel gas transfer exceeds finite kinetic parameter-work allowance"));
        }
        self.energy_j = after;
        self.last = previous.last;
        Ok(ChannelGasUpdate {
            kinetic_energy_before_j: previous.energy_j,
            kinetic_energy_after_j: after,
            parameter_work_j: work,
        })
    }

    fn free_profile(&self) -> Vec<f64> {
        let mut free: Vec<_> = self.mass.iter().zip(&self.velocity).map(|(m, v)| m * v).collect();
        solve(&self.inverse_diagonal, &self.lower_ratio, &self.upper_ratio, &mut free);
        free
    }
    fn mean(&self, velocity: &[f64]) -> f64 {
        self.area.iter().zip(velocity).map(|(a, v)| a * v).sum::<f64>() / self.total_area
    }

    pub(super) fn pressure_load(&self) -> Result<SlitChannelLoad, AcousticRealizeError> {
        if self.exit != SlitExit::StaticPressure {
            return Err(invalid("affine channel condensation requires the actual static-pressure terminal"));
        }
        let load = SlitChannelLoad {
            impedance_pa_s_m3: 1.0 / (self.total_area * self.mean_pressure_response),
            history_pressure_pa: self.mean(&self.free_profile()) / self.mean_pressure_response,
        };
        if !load.impedance_pa_s_m3.is_finite() || load.impedance_pa_s_m3 <= 0.0
            || !load.history_pressure_pa.is_finite() {
            return Err(invalid("unsteady channel affine pressure boundary cannot be represented"));
        }
        Ok(load)
    }

    pub(super) fn steady_bulk_speed(&self, pressure: f64) -> Result<f64, AcousticRealizeError> {
        if !pressure.is_finite() { return Err(invalid("unsteady channel pressure must be finite")); }
        let resistance = 1.0 / self.steady_mean_response;
        if !resistance.is_finite() || resistance <= 0.0 {
            return Err(invalid("channel steady resistance cannot be represented"));
        }
        let speed = match self.exit {
            SlitExit::StaticPressure => pressure / resistance,
            SlitExit::UniformKineticHead => 2.0 * pressure / (resistance
                + det::hypot(resistance, det::sqrt(2.0 * self.density * pressure.abs()))),
        };
        if !speed.is_finite() { return Err(invalid("channel steady response cannot be represented")); }
        Ok(speed)
    }

    pub(super) fn step(&mut self, pressure: f64) -> Result<SlitChannelStep, AcousticRealizeError> {
        if !pressure.is_finite() { return Err(invalid("unsteady channel pressure must be finite")); }
        for ((rhs, mass), velocity) in self.trial.iter_mut().zip(&self.mass).zip(&self.velocity) {
            *rhs = mass * velocity;
        }
        solve(&self.inverse_diagonal, &self.lower_ratio, &self.upper_ratio, &mut self.trial);
        let target = self.mean(&self.trial) + self.mean_pressure_response * pressure;
        let speed = match self.exit {
            SlitExit::StaticPressure => target,
            SlitExit::UniformKineticHead => 2.0 * target / (1.0 + det::sqrt(
                1.0 + 2.0 * self.mean_pressure_response * self.density * target.abs())),
        };
        let exit_pressure = match self.exit {
            SlitExit::StaticPressure => 0.0,
            SlitExit::UniformKineticHead => 0.5 * self.density * speed * speed.abs(),
        };
        let channel_pressure = pressure - exit_pressure;
        if ![target, speed, exit_pressure, channel_pressure].iter().all(|x| x.is_finite()) {
            return Err(invalid("channel implicit pressure/flow solve is not finite"));
        }
        let mut energy = 0.0;
        let mut delta = 0.0;
        let mut numerical = 0.0;
        let mut shear = 0.0;
        for i in 0..self.trial.len() {
            let old = self.velocity[i];
            self.trial[i] += self.pressure_response[i] * channel_pressure;
            let v = self.trial[i];
            energy += self.mass[i] * v * v;
            delta += self.mass[i] * (v - old) * (v + old);
            numerical += self.mass[i] * (v - old).powi(2);
            if i != 0 { shear += self.conductance[i] * (v - self.trial[i - 1]).powi(2); }
        }
        let n = self.trial.len();
        shear += self.conductance[0] * self.trial[0].powi(2)
            + self.conductance[n] * self.trial[n - 1].powi(2);
        let bulk = self.mean(&self.trial);
        let flow = self.total_area * bulk;
        let mut step = SlitChannelStep { bulk_speed_m_s: bulk, flow_m3_s: flow,
            kinetic_energy_j: 0.5 * energy, stored_delta_j: 0.5 * delta,
            pressure_work_j: pressure * flow * self.dt,
            viscous_dissipation_j: shear * self.dt,
            exit_work_j: exit_pressure * flow * self.dt,
            numerical_dissipation_j: 0.5 * numerical, energy_residual_j: 0.0 };
        step.energy_residual_j = step.pressure_work_j - step.stored_delta_j
            - step.viscous_dissipation_j - step.exit_work_j - step.numerical_dissipation_j;
        let scale = step.pressure_work_j.abs() + step.stored_delta_j.abs()
            + step.viscous_dissipation_j + step.exit_work_j + step.numerical_dissipation_j;
        if ![bulk, flow, step.kinetic_energy_j, step.stored_delta_j, step.pressure_work_j,
            step.viscous_dissipation_j, step.exit_work_j, step.numerical_dissipation_j,
            step.energy_residual_j, scale].iter().all(|x| x.is_finite())
            || step.exit_work_j < 0.0
            || step.energy_residual_j.abs() > 4096.0 * f64::EPSILON * scale {
            return Err(invalid("channel step does not close finite dissipation and kinetic storage"));
        }
        core::mem::swap(&mut self.velocity, &mut self.trial);
        self.energy_j = step.kinetic_energy_j;
        self.last = step;
        Ok(step)
    }
}
