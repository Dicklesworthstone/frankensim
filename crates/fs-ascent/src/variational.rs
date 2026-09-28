//! State-form weak-constraint variational assimilation with matched adjoints.
//!
//! Every time-frame state is a decision variable. Minimize the sum of a
//! background penalty, per-interval model-error penalties, and observation
//! penalties. The model error is x[k+1] - M[k](x[k]), NOT a sensor residual.
//! A model supplies its discrete forecast and a pullback from the SAME tape.
//! No state covariance, numerical Jacobian, or new optimizer is constructed.
//!
//! Background, model and observation errors here are independent, fixed,
//! diagonal Gaussian declarations. Model standard deviations describe an
//! INTERVAL defect, not a continuous-time noise intensity. They must be
//! positive; zero model error requires a different, strong-constraint problem.
//! Controls are dimensionless: x[i] = reference[i] + scale[i] * control[i].
//! See Tremolet, ECMWF TM 520 (2007), and the state-form weak-constraint
//! objective in Freitag, GAMM-Mitteilungen 43, e202000014 (2020).
//!
//! This is a numerical MAP objective, not a posterior distribution, calibrated
//! uncertainty, a certified adjoint, or validation of the physical model. The
//! caller owns units, model/data identity and the correctness of derivatives.

/// One independent scalar observation. IDs strictly increase in the supplied
/// array; frame selects an exact supplied time, not a nearest-time lookup.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindowObservation {
    pub id: u64,
    pub frame: usize,
    pub channel: u64,
    pub value: f64,
    pub sigma: f64,
}

/// Borrowed declarations copied once by `WeakConstraintWindow::new`.
/// Frame-major arrays reference/scale have frames*dimension entries; model_std
/// has (frames-1)*dimension. Background arrays have dimension entries.
pub struct WindowInputs<'a> {
    pub times: &'a [f64],
    pub background: &'a [f64],
    pub background_std: &'a [f64],
    pub model_std: &'a [f64],
    pub reference: &'a [f64],
    pub scale: &'a [f64],
    pub observations: &'a [WindowObservation],
}

#[derive(Debug, Clone, Copy)]
pub struct WindowLimits {
    pub max_frames: usize,
    pub max_dimension: usize,
    pub max_observations: usize,
    /// Owned scalar/ID slots: times + 2*dimension + model_std + 2*controls
    /// + 5*observations. Vec descriptors and allocator metadata are excluded.
    pub max_owned_components: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindowStage { Forecast, Pullback, Observation }

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WindowError {
    Invalid(&'static str),
    Limit { required: usize, limit: usize },
    EvaluationLimit,
    ModelCallLimit,
    Allocation,
    NonFinite(&'static str),
    Model { stage: WindowStage, index: usize, message: String },
    Cancelled,
}
impl std::fmt::Display for WindowError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "weak-constraint assimilation: {self:?}")
    }
}
impl std::error::Error for WindowError {}

/// An immutable forward model and observation map. A tape may borrow the
/// model, but must own any initial state it needs. `forecast_vjp` consumes that
/// exact tape; do not recompute a different mesh or differentiate solver stops.
/// All output slices must be overwritten, including zero derivatives.
/// Long callbacks must poll the provided latched cancellation check and bound
/// their own work/tape storage. Callback internals are NOT covered by our cap.
pub trait WindowModel {
    type Tape<'a> where Self: 'a;
    fn dimension(&self) -> usize;
    #[allow(clippy::too_many_arguments)]
    fn forecast<'a>(
        &'a self, interval: usize, start: f64, end: f64, state: &[f64],
        output: &mut [f64], cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Self::Tape<'a>, String>;
    fn forecast_vjp(
        &self, tape: Self::Tape<'_>, seed: &[f64], state_bar: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(), String>;
    /// Return prediction and dh/dstate, WITHOUT noise scaling or loss terms.
    fn observe(
        &self, channel: u64, time: f64, state: &[f64], state_bar: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<f64, String>;
}

/// Non-cloneable cumulative work allowance. Failed evaluations and callbacks
/// stay charged; copying a numerical checkpoint does not copy this authority.
#[derive(Debug)]
pub struct WindowControl {
    max_evaluations: usize,
    max_calls: usize,
    workspace: usize,
    evaluations: usize,
    calls: usize,
}
impl WindowControl {
    pub fn new(max_evaluations: usize, max_model_calls: usize, max_workspace_components: usize) -> Self {
        Self { max_evaluations, max_calls: max_model_calls, workspace: max_workspace_components,
            evaluations: 0, calls: 0 }
    }
    pub fn evaluations(&self) -> usize { self.evaluations }
    pub fn model_calls(&self) -> usize { self.calls }
    /// Only increase ceilings, never reset already spent work.
    pub fn extend(&mut self, evaluations: usize, calls: usize, workspace: usize) -> Result<(), WindowError> {
        if evaluations < self.max_evaluations || calls < self.max_calls || workspace < self.workspace {
            return Err(WindowError::Invalid("resource ceilings may only increase"));
        }
        self.max_evaluations = evaluations; self.max_calls = calls; self.workspace = workspace;
        Ok(())
    }
    fn begin(&mut self, calls: usize, required: usize) -> Result<(), WindowError> {
        self.admit_workspace(required)?;
        if self.evaluations == self.max_evaluations { return Err(WindowError::EvaluationLimit); }
        if calls > self.max_calls - self.calls { return Err(WindowError::ModelCallLimit); }
        self.evaluations += 1;
        Ok(())
    }
    fn admit_workspace(&self, required: usize) -> Result<(), WindowError> {
        if required > self.workspace { Err(WindowError::Limit { required, limit: self.workspace }) }
        else { Ok(()) }
    }
    fn charge(&mut self) -> Result<(), WindowError> {
        if self.calls == self.max_calls { return Err(WindowError::ModelCallLimit); }
        self.calls += 1; Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct WeakConstraintWindow {
    dimension: usize,
    times: Vec<f64>,
    background: Vec<f64>,
    background_std: Vec<f64>,
    model_std: Vec<f64>,
    reference: Vec<f64>,
    scale: Vec<f64>,
    observations: Vec<WindowObservation>,
}

/// Complete objective and gradient at ONE control point. Frame-major states
/// and physical model defects are retained together. Gradients use the declared
/// dimensionless control coordinates, not the physical state coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowEvaluation {
    states: Vec<f64>,
    gradient: Vec<f64>,
    model_errors: Vec<f64>,
    background_cost: f64,
    model_cost: f64,
    observation_cost: f64,
    value: f64,
}
impl WindowEvaluation {
    pub fn states(&self) -> &[f64] { &self.states }
    pub fn gradient(&self) -> &[f64] { &self.gradient }
    pub fn model_errors(&self) -> &[f64] { &self.model_errors }
    pub fn background_cost(&self) -> f64 { self.background_cost }
    pub fn model_cost(&self) -> f64 { self.model_cost }
    pub fn observation_cost(&self) -> f64 { self.observation_cost }
    pub fn value(&self) -> f64 { self.value }
}

fn finite(value: f64, label: &'static str) -> Result<f64, WindowError> {
    if value.is_finite() { Ok(value) } else { Err(WindowError::NonFinite(label)) }
}
fn poll<C: FnMut() -> bool>(cancelled: &mut C) -> Result<(), WindowError> {
    if cancelled() { Err(WindowError::Cancelled) } else { Ok(()) }
}
fn copy<T: Clone>(values: &[T]) -> Result<Vec<T>, WindowError> {
    let mut out = Vec::new(); out.try_reserve_exact(values.len()).map_err(|_| WindowError::Allocation)?;
    out.extend_from_slice(values); Ok(out)
}
fn zeros(n: usize) -> Result<Vec<f64>, WindowError> {
    let mut out = Vec::new(); out.try_reserve_exact(n).map_err(|_| WindowError::Allocation)?;
    out.resize(n, 0.0); Ok(out)
}
fn extent(parts: &[Option<usize>]) -> Result<usize, WindowError> {
    parts.iter().try_fold(0_usize, |s, p| p.and_then(|v| s.checked_add(v)))
        .ok_or(WindowError::Invalid("window extent overflow"))
}
fn checked_slice<C: FnMut() -> bool>(values: &[f64], positive: bool, cancelled: &mut C) -> Result<(), WindowError> {
    for chunk in values.chunks(256) {
        poll(cancelled)?;
        if chunk.iter().any(|v| !v.is_finite() || (positive && *v <= 0.0)) {
            return Err(WindowError::Invalid("expected finite values and strictly positive standard deviations/scales"));
        }
    }
    Ok(())
}
#[derive(Default)]
struct Cost { sum: f64, correction: f64 }
impl Cost {
    fn add(&mut self, value: f64) -> Result<(), WindowError> {
        let next = finite(self.sum + value, "cost sum")?;
        self.correction += if self.sum.abs() >= value.abs() { (self.sum-next)+value } else { (value-next)+self.sum };
        self.sum = next; finite(self.correction, "cost compensation")?; Ok(())
    }
    fn value(&self) -> Result<f64, WindowError> { finite(self.sum + self.correction, "objective") }
}

impl WeakConstraintWindow {
    pub fn new<C: FnMut() -> bool>(input: WindowInputs<'_>, limits: WindowLimits, cancelled: &mut C)
        -> Result<Self, WindowError>
    {
        poll(cancelled)?;
        let n = input.background.len(); let frames = input.times.len();
        if n == 0 || n > limits.max_dimension || frames == 0 || frames > limits.max_frames
            || input.observations.len() > limits.max_observations {
            return Err(WindowError::Invalid("nonempty bounded state/window required"));
        }
        let p = n.checked_mul(frames).ok_or(WindowError::Invalid("control extent overflow"))?;
        let required = extent(&[Some(frames), n.checked_mul(2), Some(p-n), p.checked_mul(2),
            input.observations.len().checked_mul(5)])?;
        if required > limits.max_owned_components { return Err(WindowError::Limit { required, limit: limits.max_owned_components }); }
        if input.background_std.len()!=n || input.model_std.len()!=p-n || input.reference.len()!=p || input.scale.len()!=p {
            return Err(WindowError::Invalid("window array dimensions differ"));
        }
        checked_slice(input.times, false, cancelled)?;
        if input.times.windows(2).any(|t| t[1] <= t[0] || !(t[1]-t[0]).is_finite()) {
            return Err(WindowError::Invalid("frame times must strictly increase with finite intervals"));
        }
        for values in [input.background, input.reference] { checked_slice(values, false, cancelled)?; }
        for values in [input.background_std, input.model_std, input.scale] { checked_slice(values, true, cancelled)?; }
        let mut last = None;
        for reading in input.observations {
            poll(cancelled)?;
            if reading.frame >= frames || !reading.value.is_finite() || !reading.sigma.is_finite() || reading.sigma <= 0.0
                || last.is_some_and(|id| reading.id <= id) {
                return Err(WindowError::Invalid("invalid observation frame, value, noise or increasing ID order"));
            }
            last = Some(reading.id);
        }
        let window = Self { dimension:n, times:copy(input.times)?, background:copy(input.background)?,
            background_std:copy(input.background_std)?, model_std:copy(input.model_std)?,
            reference:copy(input.reference)?, scale:copy(input.scale)?, observations:copy(input.observations)? };
        poll(cancelled)?; Ok(window)
    }
    pub fn dimension(&self) -> usize { self.dimension }
    pub fn times(&self) -> &[f64] { &self.times }
    pub fn control_dimension(&self) -> usize { self.reference.len() }
    pub fn observations(&self) -> &[WindowObservation] { &self.observations }
    /// 3*frames*dimension + 2*dimension live scalar slots, including outputs.
    /// Excludes the input controls/window, callback tapes and optimizer storage.
    pub fn workspace_components(&self) -> Result<usize, WindowError> {
        extent(&[self.control_dimension().checked_mul(3), self.dimension.checked_mul(2)])
    }
    pub fn model_calls_per_evaluation(&self) -> Result<usize, WindowError> {
        extent(&[(self.times.len()-1).checked_mul(2), Some(self.observations.len())])
    }

    /// Evaluate the state-form weak-constraint objective and its exact supplied
    /// discrete pullbacks. Each interval is forecast once and pulled back once;
    /// observations each use one prediction/gradient callback. Only one model
    /// tape is live. Failures return no partial state history or gradient.
    pub fn evaluate<M: WindowModel, C: FnMut() -> bool>(
        &self, model: &M, controls: &[f64], control: &mut WindowControl, cancelled: &mut C,
    ) -> Result<WindowEvaluation, WindowError> {
        poll(cancelled)?;
        let n = self.dimension; let p = self.control_dimension();
        if model.dimension()!=n || controls.len()!=p { return Err(WindowError::Invalid("model or control dimensions differ")); }
        checked_slice(controls, false, cancelled)?;
        control.begin(self.model_calls_per_evaluation()?, self.workspace_components()?)?;
        let mut states = zeros(p)?; let mut gradient = zeros(p)?; let mut errors = zeros(p-n)?;
        let mut forecast = zeros(n)?; let mut seed = zeros(n)?; let mut bar = zeros(n)?;
        for i in 0..p {
            if i%256==0 { poll(cancelled)?; }
            states[i] = finite(self.scale[i].mul_add(controls[i], self.reference[i]), "physical state")?;
        }
        let (mut background, mut dynamics, mut data) = (Cost::default(), Cost::default(), Cost::default());
        for j in 0..n {
            if j%256==0 { poll(cancelled)?; }
            let residual = finite((states[j]-self.background[j])/self.background_std[j], "background residual")?;
            background.add(finite((0.5*residual)*residual, "background cost")?)?;
            gradient[j] = finite(residual/self.background_std[j], "background derivative")?;
        }
        for k in 0..self.times.len()-1 {
            poll(cancelled)?; forecast.fill(f64::NAN); control.charge()?;
            let mut stopped = false; let mut check = || { stopped |= cancelled(); stopped };
            let result = model.forecast(k,self.times[k],self.times[k+1],&states[k*n..(k+1)*n],&mut forecast,&mut check);
            if check() { return Err(WindowError::Cancelled); }
            let tape = result.map_err(|message| WindowError::Model { stage:WindowStage::Forecast,index:k,message })?;
            for j in 0..n {
                if j%256==0 { poll(cancelled)?; }
                finite(forecast[j], "unwritten or non-finite forecast")?;
                let i=k*n+j; let sigma=self.model_std[i];
                errors[i] = finite(states[i+n]-forecast[j], "model defect")?;
                let residual = finite(errors[i]/sigma, "standardized model defect")?;
                dynamics.add(finite((0.5*residual)*residual,"model cost")?)?;
                seed[j] = finite(residual/sigma, "model derivative")?;
                gradient[i+n] = finite(gradient[i+n]+seed[j], "arrival derivative")?;
            }
            poll(cancelled)?; bar.fill(f64::NAN); control.charge()?;
            let mut stopped=false; let mut check=|| { stopped |= cancelled(); stopped };
            let result=model.forecast_vjp(tape,&seed,&mut bar,&mut check);
            if check() { return Err(WindowError::Cancelled); }
            result.map_err(|message| WindowError::Model { stage:WindowStage::Pullback,index:k,message })?;
            for j in 0..n {
                if j%256==0 { poll(cancelled)?; }
                finite(bar[j],"unwritten or non-finite model derivative")?;
                gradient[k*n+j]=finite(gradient[k*n+j]-bar[j],"departure derivative")?;
            }
        }
        for (index, reading) in self.observations.iter().enumerate() {
            poll(cancelled)?; bar.fill(f64::NAN); control.charge()?;
            let start=reading.frame*n;
            let mut stopped=false; let mut check=|| { stopped |= cancelled(); stopped };
            let result=model.observe(reading.channel,self.times[reading.frame],&states[start..start+n],&mut bar,&mut check);
            if check() { return Err(WindowError::Cancelled); }
            let prediction=result.map_err(|message| WindowError::Model { stage:WindowStage::Observation,index,message })?;
            finite(prediction,"observation prediction")?;
            let residual=finite((prediction-reading.value)/reading.sigma,"observation residual")?;
            data.add(finite((0.5*residual)*residual,"observation cost")?)?;
            let weight=finite(residual/reading.sigma,"observation derivative scale")?;
            for j in 0..n {
                if j%256==0 { poll(cancelled)?; }
                finite(bar[j],"unwritten or non-finite observation derivative")?;
                gradient[start+j]=finite(weight.mul_add(bar[j],gradient[start+j]),"observation derivative")?;
            }
        }
        for (i,g) in gradient.iter_mut().enumerate() {
            if i%256==0 { poll(cancelled)?; }
            *g=finite(*g*self.scale[i],"scaled control derivative")?;
        }
        let background_cost=background.value()?; let model_cost=dynamics.value()?; let observation_cost=data.value()?;
        let mut total=Cost::default(); for cost in [background_cost,model_cost,observation_cost] { total.add(cost)?; }
        poll(cancelled)?;
        Ok(WindowEvaluation { states,gradient,model_errors:errors,background_cost,model_cost,observation_cost,value:total.value()? })
    }
}

#[cfg(test)]
#[path = "variational/tests.rs"]
mod tests;
