//! Matrix-free native structural adapter with cached original heat/volume trials.
use super::*;
use fs_time::galpha::{SecondOrderOperatorWeights, SecondOrderProblem};
use std::cell::{Cell, RefCell};

pub(super) struct ThermalOperator<'a, 'p, C> {
    model: &'a GasMembrane<'p>,
    masses: &'a [f64],
    initial: &'a GasMembraneThermalEvaluation,
    old_q: &'a [f64],
    heat: f64,
    exterior_heat: Option<f64>,
    controls: GasMembraneThermalControls,
    cancel: RefCell<C>,
    stop: RefCell<Option<GasMembraneMotionError>>,
    failure: RefCell<Option<GasMembraneMotionError>>,
    cache: RefCell<Option<gas::Trial>>,
    evaluations: Cell<usize>,
    gas_trials: Cell<usize>,
}
impl<'a, 'p, C: FnMut() -> bool> ThermalOperator<'a, 'p, C> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        model: &'a GasMembrane<'p>,
        masses: &'a [f64],
        initial: &'a GasMembraneThermalEvaluation,
        old_q: &'a [f64],
        heat: f64,
        exterior_heat: Option<f64>,
        controls: GasMembraneThermalControls,
        cancelled: C,
    ) -> Self {
        Self {
            model,
            masses,
            initial,
            old_q,
            heat,
            exterior_heat,
            controls,
            cancel: RefCell::new(cancelled),
            stop: RefCell::new(None),
            failure: RefCell::new(None),
            cache: RefCell::new(None),
            evaluations: Cell::new(0),
            gas_trials: Cell::new(0),
        }
    }
    pub(super) fn evaluations(&self) -> usize {
        self.evaluations.get()
    }
    pub(super) fn gas_trials(&self) -> usize {
        self.gas_trials.get()
    }
    pub(super) fn take_stop(&self) -> Option<GasMembraneMotionError> {
        self.stop.borrow_mut().take()
    }
    pub(super) fn take_failure(&self) -> Option<GasMembraneMotionError> {
        self.failure.borrow_mut().take()
    }
    pub(super) fn cancelled(&self) -> bool {
        if self.stop.borrow().is_some() {
            return true;
        }
        if (self.cancel.borrow_mut())() {
            *self.stop.borrow_mut() = Some(GasMembraneMotionError::Cancelled);
            return true;
        }
        false
    }
    fn begin(&self) -> bool {
        if self.cancelled() {
            return false;
        }
        if self.evaluations.get() >= self.controls.maximum_evaluations {
            *self.stop.borrow_mut() = Some(GasMembraneMotionError::EvaluationBudget);
            return false;
        }
        self.evaluations.set(self.evaluations.get() + 1);
        true
    }
    fn remember(&self, error: GasMembraneMotionError) {
        if matches!(error, GasMembraneMotionError::Cancelled) {
            *self.stop.borrow_mut() = Some(GasMembraneMotionError::Cancelled);
        } else {
            // Source/geometry trial refusals are recoverable by native backtracking.
            *self.failure.borrow_mut() = Some(error);
        }
    }
    fn matches(&self, geometry: &gas::Geometry, q: &[f64]) -> bool {
        q.len() == self.masses.len()
            && self
                .model
                .free_dofs
                .iter()
                .zip(q)
                .all(|(&i, &value)| geometry.displacement[i / 3][i % 3] == value)
    }
    fn fresh(&self, stage: &[f64], endpoint: &[f64]) -> Result<gas::Trial, GasMembraneMotionError> {
        self.gas_trials.set(self.gas_trials.get() + 1);
        gas::Trial::new(
            self.model,
            self.initial,
            stage,
            endpoint,
            self.heat,
            self.exterior_heat,
            self.controls,
            &mut || self.cancelled(),
        )
    }
    fn prepare(&self, q: &[f64]) -> bool {
        if !self.begin() {
            return false;
        }
        if q.len() != self.masses.len() || q.iter().any(|x| !x.is_finite()) {
            return false;
        }
        if self
            .cache
            .borrow()
            .as_ref()
            .is_some_and(|trial| self.matches(&trial.stage, q))
        {
            return true;
        }
        let endpoint: Vec<_> = self
            .old_q
            .iter()
            .zip(q)
            .map(|(&old, &stage)| old + 2.0 * (stage - old))
            .collect();
        match self.fresh(q, &endpoint) {
            Ok(trial) => {
                *self.cache.borrow_mut() = Some(trial);
                true
            }
            Err(error) => {
                self.remember(error);
                false
            }
        }
    }
    pub(super) fn final_trial(
        &self,
        endpoint: &[f64],
    ) -> Result<gas::Trial, GasMembraneMotionError> {
        if !self.begin() {
            return Err(self.take_stop().expect("stopped operator"));
        }
        let stage: Vec<_> = self
            .old_q
            .iter()
            .zip(endpoint)
            .map(|(&old, &end)| old + 0.5 * (end - old))
            .collect();
        self.fresh(&stage, endpoint)
    }
}

impl<C: FnMut() -> bool> SecondOrderProblem for ThermalOperator<'_, '_, C> {
    fn dimension(&self) -> usize {
        self.masses.len()
    }
    fn validate_state(&self, q: &[f64], v: &[f64]) -> Result<(), String> {
        if !self.begin() {
            return Err("thermal motion cancelled or work exhausted".into());
        }
        if v.len() != self.dimension() || v.iter().any(|x| !x.is_finite()) {
            return Err("finite dimension-matched thermal motion velocity required".into());
        }
        if self
            .cache
            .borrow()
            .as_ref()
            .is_some_and(|trial| self.matches(&trial.stage, q) || self.matches(&trial.endpoint, q))
        {
            return Ok(());
        }
        // This native callback has no initial/stage/endpoint tag. Admit actual
        // geometry here; force preparation admits BOTH thermodynamic states,
        // and final_trial freshly rechecks the accepted endpoint before publish.
        match gas::geometry(self.model, q, &mut || self.cancelled()) {
            Ok(_) => Ok(()),
            Err(error) => {
                self.remember(error);
                Err("original membrane geometry/volume refused thermal motion".into())
            }
        }
    }
    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if self.cancelled() || input.len() != self.dimension() || output.len() != self.dimension() {
            return;
        }
        for ((out, value), mass) in output.iter_mut().zip(input).zip(self.masses) {
            *out = value * mass;
        }
    }
    fn damping_apply(&self, _input: &[f64], output: &mut [f64]) {
        output.fill(if self.cancelled() { f64::NAN } else { 0.0 });
    }
    fn internal_force(&self, q: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if output.len() != self.dimension() || !self.prepare(q) {
            return;
        }
        let cache = self.cache.borrow();
        let trial = cache.as_ref().expect("prepared heat/volume trial");
        for (out, &i) in output.iter_mut().zip(&self.model.free_dofs) {
            *out = trial.stage.internal[i / 3][i % 3]
                - trial.net_pressure_pa * trial.stage.gradient[i / 3][i % 3];
        }
    }
    fn tangent_apply(&self, q: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if output.len() != self.dimension() || !self.prepare(q) {
            return;
        }
        let cache = self.cache.borrow();
        match cache.as_ref().expect("prepared thermal trial").tangent(
            self.model,
            direction,
            &mut || self.cancelled(),
        ) {
            Ok(action) => output.copy_from_slice(&action),
            Err(error) => self.remember(error),
        }
    }
    fn preconditioner_apply(
        &self,
        _q: &[f64],
        weights: SecondOrderOperatorWeights,
        _outer: usize,
        _inner: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        output.fill(f64::NAN);
        if self.cancelled()
            || residual.len() != self.dimension()
            || output.len() != self.dimension()
        {
            return;
        }
        for ((out, value), mass) in output.iter_mut().zip(residual).zip(self.masses) {
            let scale = weights.mass * mass;
            if scale.is_finite() && scale > 0.0 {
                *out = value / scale;
            }
        }
    }
}
