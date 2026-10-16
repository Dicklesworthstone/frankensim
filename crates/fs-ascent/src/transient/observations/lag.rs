//! First-order sensor dynamics on the existing transient calibration path.
//!
//! For each channel, tau * dz/dt = q(t,x,p) - z, and reading = z + bias.
//! `q` is the INNER model's sensor prediction; lag channels do not feed back
//! into the physical dynamics. Unmapped outputs remain instantaneous. This is
//! response lag, NOT a timestamp shift, transport delay, or deconvolution.
//! The model is a declared approximation, not proof of instrument calibration.
//! First-order lag conventions: MathWorks PS Transfer Function documentation,
//! https://www.mathworks.com/help/simscape/ref/pstransferfunction.html .
//!
//! The existing RK45 recording advances [physical state, sensor states]; its
//! error controller and checkpoint/adjoint include every component. Small tau
//! may make this explicit system stiff: use suitable budgets/tolerances; no
//! stiffness-independent cost or equivalence to the IMEX lane is claimed.

use super::{OdeVjp, SensorFamily, SensorModel};
use std::sync::Arc;

/// A constant or an existing decision coordinate, in the model's units. This
/// adapter does not append decision variables or infer parameter scaling.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LagValue { Fixed(f64), Parameter(usize) }
impl LagValue {
    fn resolve(self, point: &[f64]) -> Result<f64, String> {
        let value = match self {
            Self::Fixed(value) => value,
            Self::Parameter(i) => *point.get(i).ok_or("unknown lag parameter")?,
        };
        if value.is_finite() { Ok(value) } else { Err("non-finite lag value".into()) }
    }
    fn bounds(self, bounds: &[[f64; 2]]) -> Result<[f64; 2], String> {
        match self {
            Self::Fixed(value) if value.is_finite() => Ok([value, value]),
            Self::Parameter(i) => bounds.get(i).copied().ok_or_else(|| "unknown lag parameter".into()),
            _ => Err("non-finite lag value".into()),
        }
    }
    fn accumulate(self, seed: f64, out: &mut [f64]) {
        if let Self::Parameter(i) = self { out[i] += seed; }
    }
}

/// Initial sensor INTERNAL state (before adding readout bias).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LagInitial {
    Value(LagValue),
    /// q(initial_time, x_initial(p), p); both dependencies are differentiated.
    Equilibrium,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LagSensor {
    /// External channel to replace/add. Channels not listed pass to the inner
    /// model unchanged. IDs are unique here; an inner ID may be shadowed.
    pub channel: u64,
    /// Always resolves against the inner model, even if the ID is shadowed.
    pub source: u64,
    /// Positive time constant in the same time unit as the integration clock.
    pub time_constant: LagValue,
    pub initial: LagInitial,
    /// Added AFTER the lag, in the signal unit; use zero when absent. A bias
    /// already in the source signal is filtered and must not be counted twice.
    pub bias: LagValue,
}

/// Reusable `SensorFamily` adapter, composable with `ObservedFamily` and
/// `CorrelatedFamily`, and hence ordinary and multi-experiment fitting.
/// `initial_time` must equal the caller's recording start; like initial_values,
/// it is a declaration, not a change to that clock. Clone checkpoints retain
/// the sensor state: do not reinitialize sensors at observation boundaries.
pub struct LaggedFamily<F> {
    inner: F,
    sensors: Arc<[LagSensor]>,
    initial_time: f64,
    max_components: usize,
}
impl<F: SensorFamily> LaggedFamily<F> {
    /// Caps are checked before owned descriptor allocation. max_components
    /// caps n + sensor_count, not model-owned memory or the full RK workspace.
    /// Derivative actions allocate at most 2n+p extra scalar scratch; all
    /// expensive inner callbacks remain responsible for their own resources.
    pub fn new(inner: F, sensors: &[LagSensor], initial_time: f64,
               max_sensors: usize, max_components: usize) -> Result<Self, String> {
        if sensors.is_empty() || sensors.len() > max_sensors || sensors.len() >= max_components
            || !initial_time.is_finite()
        { return Err("invalid sensor count, component budget or initial time".into()); }
        let bounds = inner.bounds();
        if bounds.iter().any(|b| !b[0].is_finite() || !b[1].is_finite() || b[0] > b[1]) {
            return Err("invalid inner parameter bounds".into());
        }
        let mut channels = Vec::new();
        channels.try_reserve_exact(sensors.len()).map_err(|_| "lag descriptor allocation refused")?;
        for sensor in sensors {
            let tau = sensor.time_constant.bounds(bounds)?;
            if tau[0] <= 0.0 || !tau[0].recip().is_finite() {
                return Err("lag time constant must stay positive with a finite reciprocal over its box".into());
            }
            sensor.bias.bounds(bounds)?;
            if let LagInitial::Value(value) = sensor.initial { value.bounds(bounds)?; }
            channels.push(sensor.channel);
        }
        channels.sort_unstable();
        if channels.windows(2).any(|ids| ids[0] == ids[1]) { return Err("duplicate lag output channel".into()); }
        Ok(Self { inner, sensors: Arc::from(sensors), initial_time, max_components })
    }
}

/// An immutable parameter point and its augmented physical/sensor dynamics.
pub struct LaggedModel<M> {
    inner: M,
    sensors: Arc<[LagSensor]>,
    initial_time: f64,
    n: usize,
    p: usize,
    initial: Vec<f64>,
    taus: Vec<f64>,
    biases: Vec<f64>,
}
fn finite(values: &[f64]) -> bool { values.iter().all(|v| v.is_finite()) }
fn add(total: &mut [f64], partial: &[f64]) -> Result<(), String> {
    for (a, b) in total.iter_mut().zip(partial) { *a += b; }
    if finite(total) { Ok(()) } else { Err("non-finite lag derivative accumulation".into()) }
}
impl<F: SensorFamily> SensorFamily for LaggedFamily<F> {
    type Model = LaggedModel<F::Model>;
    fn bounds(&self) -> &[[f64; 2]] { self.inner.bounds() }
    fn instantiate(&self, point: &[f64]) -> Result<Self::Model, String> {
        if point.len() != self.bounds().len() || point.iter().zip(self.bounds())
            .any(|(p,b)| !p.is_finite() || *p < b[0] || *p > b[1])
        { return Err("lag decision point is outside its declared parameter box".into()); }
        let inner = self.inner.instantiate(point)?;
        let (n, p) = (inner.dimension(), inner.parameter_count());
        let total = n.checked_add(self.sensors.len()).ok_or("lag state dimension overflow")?;
        if n == 0 || total > self.max_components || p != point.len()
            || inner.initial_values().len() != n || !finite(inner.initial_values())
        { return Err("invalid inner state/parameter dimensions or augmented component limit".into()); }
        let mut initial = Vec::new();
        initial.try_reserve_exact(total).map_err(|_| "lag state allocation refused")?;
        initial.extend_from_slice(inner.initial_values());
        let mut taus = Vec::with_capacity(self.sensors.len());
        let mut biases = Vec::with_capacity(self.sensors.len());
        for sensor in self.sensors.iter() {
            let tau = sensor.time_constant.resolve(point)?;
            if tau <= 0.0 || !tau.recip().is_finite() { return Err("invalid lag time constant".into()); }
            let z = match sensor.initial {
                LagInitial::Value(value) => value.resolve(point)?,
                LagInitial::Equilibrium => inner.predict(sensor.source, self.initial_time, &initial[..n])?,
            };
            let bias = sensor.bias.resolve(point)?;
            if !z.is_finite() || !(z + bias).is_finite() { return Err("non-finite initial sensor output".into()); }
            initial.push(z); taus.push(tau); biases.push(bias);
        }
        Ok(LaggedModel { inner, sensors: self.sensors.clone(), initial_time: self.initial_time,
            n, p, initial, taus, biases })
    }
}
impl<M: SensorModel> LaggedModel<M> {
    pub fn physical_dimension(&self) -> usize { self.n }
    fn state_ok(&self, time: f64, state: &[f64]) -> bool {
        time.is_finite() && time >= self.initial_time && state.len() == self.initial.len() && finite(state)
            && self.inner.dimension() == self.n && self.inner.parameter_count() == self.p
    }
    fn output(&self, channel: u64) -> Option<usize> { self.sensors.iter().position(|s| s.channel == channel) }
    fn source_vjp(&self, i: usize, t: f64, x: &[f64], seed: f64,
                  xb: &mut [f64], pb: &mut [f64]) -> Result<(), String> {
        if !seed.is_finite() { return Err("non-finite lag source cotangent".into()); }
        xb.fill(f64::NAN); pb.fill(f64::NAN);
        self.inner.prediction_vjp(self.sensors[i].source, t, x, seed, xb, pb)?;
        if finite(xb) && finite(pb) { Ok(()) } else { Err("non-finite or unwritten lag source derivative".into()) }
    }
}
impl<M: SensorModel> OdeVjp for LaggedModel<M> {
    fn dimension(&self) -> usize { self.initial.len() }
    fn parameter_count(&self) -> usize { self.p }
    fn rhs(&self, time: f64, state: &[f64], out: &mut [f64]) {
        out.fill(f64::NAN);
        if !self.state_ok(time, state) || out.len() != self.dimension() { return; }
        self.inner.rhs(time, &state[..self.n], &mut out[..self.n]);
        for (i, sensor) in self.sensors.iter().enumerate() {
            // OdeVjp's RHS cannot return a typed model error. Leave poison on
            // failure so production integration returns NonFiniteRhs instead.
            if let Ok(q) = self.inner.predict(sensor.source, time, &state[..self.n]) {
                out[self.n+i] = (q - state[self.n+i]) / self.taus[i];
            }
        }
    }
    fn rhs_vjp(&self, time: f64, state: &[f64], seed: &[f64],
               xb: &mut [f64], pb: &mut [f64]) -> Result<(), String> {
        xb.fill(f64::NAN); pb.fill(f64::NAN);
        if !self.state_ok(time, state) || seed.len() != self.dimension() || !finite(seed)
            || xb.len() != self.dimension() || pb.len() != self.p
        { return Err("invalid lag RHS cotangent shape/value".into()); }
        self.inner.rhs_vjp(time, &state[..self.n], &seed[..self.n], &mut xb[..self.n], pb)?;
        if !finite(&xb[..self.n]) || !finite(pb) { return Err("non-finite or unwritten physical derivative".into()); }
        let (mut dx, mut dp) = (vec![0.0; self.n], vec![0.0; self.p]);
        for (i, sensor) in self.sensors.iter().enumerate() {
            let q = self.inner.predict(sensor.source, time, &state[..self.n])?;
            let speed = (q - state[self.n+i]) / self.taus[i];
            let b = seed[self.n+i] / self.taus[i];
            if !speed.is_finite() || !b.is_finite() { return Err("non-finite lag dynamics/derivative".into()); }
            self.source_vjp(i, time, &state[..self.n], b, &mut dx, &mut dp)?;
            add(&mut xb[..self.n], &dx)?; add(pb, &dp)?;
            xb[self.n+i] = -b;
            sensor.time_constant.accumulate(-b * speed, pb);
        }
        if finite(xb) && finite(pb) { Ok(()) } else { Err("non-finite lag RHS pullback".into()) }
    }
}
impl<M: SensorModel> SensorModel for LaggedModel<M> {
    fn initial_values(&self) -> &[f64] { &self.initial }
    fn initial_vjp(&self, seed: &[f64], pb: &mut [f64]) -> Result<(), String> {
        pb.fill(f64::NAN);
        if seed.len() != self.dimension() || pb.len() != self.p || !finite(seed) {
            return Err("invalid lag initial cotangent".into());
        }
        let mut physical = seed[..self.n].to_vec();
        let (mut dx, mut dp) = (vec![0.0; self.n], vec![0.0; self.p]);
        pb.fill(0.0);
        for (i, sensor) in self.sensors.iter().enumerate() {
            let b = seed[self.n+i];
            match sensor.initial {
                LagInitial::Value(value) => value.accumulate(b, pb),
                LagInitial::Equilibrium => {
                    self.source_vjp(i, self.initial_time, &self.initial[..self.n], b, &mut dx, &mut dp)?;
                    add(&mut physical, &dx)?; add(pb, &dp)?;
                }
            }
        }
        dp.fill(f64::NAN);
        self.inner.initial_vjp(&physical, &mut dp)?;
        if !finite(&dp) { return Err("non-finite or unwritten physical initial derivative".into()); }
        add(pb, &dp)
    }
    fn predict(&self, channel: u64, time: f64, state: &[f64]) -> Result<f64, String> {
        if !self.state_ok(time, state) { return Err("invalid lag prediction state/time".into()); }
        let value = match self.output(channel) {
            Some(i) => state[self.n+i] + self.biases[i],
            None => self.inner.predict(channel, time, &state[..self.n])?,
        };
        if value.is_finite() { Ok(value) } else { Err("non-finite lag sensor prediction".into()) }
    }
    fn prediction_vjp(&self, channel: u64, time: f64, state: &[f64], seed: f64,
                      xb: &mut [f64], pb: &mut [f64]) -> Result<(), String> {
        xb.fill(f64::NAN); pb.fill(f64::NAN);
        if !self.state_ok(time, state) || !seed.is_finite()
            || xb.len() != self.dimension() || pb.len() != self.p
        { return Err("invalid lag sensor cotangent".into()); }
        match self.output(channel) {
            Some(i) => {
                xb.fill(0.0); pb.fill(0.0); xb[self.n+i] = seed;
                self.sensors[i].bias.accumulate(seed, pb);
            }
            None => {
                self.inner.prediction_vjp(channel, time, &state[..self.n], seed, &mut xb[..self.n], pb)?;
                xb[self.n..].fill(0.0);
            }
        }
        if finite(xb) && finite(pb) { Ok(()) } else { Err("non-finite or unwritten prediction derivative".into()) }
    }
}

#[cfg(test)]
mod tests;
