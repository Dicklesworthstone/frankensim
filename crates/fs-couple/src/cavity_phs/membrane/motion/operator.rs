//! Fallible physical adapter; the original structural owner performs all steps.
use super::*;
use fs_time::galpha::{SecondOrderOperatorWeights, SecondOrderProblem};
use std::cell::{Cell, RefCell};

pub(super) struct Operator<'a, 'p, C> {
    model: &'a GasMembrane<'p>,
    masses: &'a [f64],
    gas: GasEntropySolve,
    maximum: usize,
    zero: Vec<[f64; 3]>,
    cancel: RefCell<C>,
    stop: RefCell<Option<GasMembraneMotionError>>,
    failure: RefCell<Option<GasMembraneError>>,
    cache: RefCell<Option<Trial>>,
    evaluations: Cell<usize>,
    gas_trials: Cell<usize>,
}
impl<'a, 'p, C: FnMut() -> bool> Operator<'a, 'p, C> {
    pub(super) fn new(
        model: &'a GasMembrane<'p>,
        masses: &'a [f64],
        controls: GasMembraneMotionControls,
        cancelled: C,
    ) -> Self {
        Self {
            model,
            masses,
            gas: controls.gas,
            maximum: controls.maximum_evaluations,
            zero: vec![[0.0; 3]; model.reference_m.len()],
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
    pub(super) fn take_failure(&self) -> Option<GasMembraneError> {
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
        if self.evaluations.get() >= self.maximum {
            *self.stop.borrow_mut() = Some(GasMembraneMotionError::EvaluationBudget);
            return false;
        }
        self.evaluations.set(self.evaluations.get() + 1);
        true
    }
    fn remember(&self, error: GasMembraneError) {
        if matches!(error, GasMembraneError::Cancelled) {
            if self.stop.borrow().is_none() {
                *self.stop.borrow_mut() = Some(GasMembraneMotionError::Cancelled);
            }
        } else {
            // An inadmissible stage can recover during native backtracking.
            // Preserve its cause only if the complete solve ultimately fails.
            *self.failure.borrow_mut() = Some(error);
        }
    }
    fn fresh(&self, q: &[f64], loads: &[[f64; 3]]) -> Result<Trial, GasMembraneError> {
        self.gas_trials.set(self.gas_trials.get() + 1);
        self.model
            .trial(&self.model.expand(q, false), loads, self.gas, &mut || {
                self.cancelled()
            })
    }
    fn prepare(&self, q: &[f64]) -> bool {
        if !self.begin() {
            return false;
        }
        if q.len() != self.masses.len() || q.iter().any(|q| !q.is_finite()) {
            self.remember(GasMembraneError::Invalid(
                "finite dimension-matched motion coordinates required",
            ));
            return false;
        }
        let matches = self.cache.borrow().as_ref().is_some_and(|trial| {
            self.model
                .free_dofs
                .iter()
                .zip(q)
                .all(|(&i, &q)| trial.value.displacement_m[i / 3][i % 3] == q)
        });
        if !matches {
            match self.fresh(q, &self.zero) {
                Ok(trial) => *self.cache.borrow_mut() = Some(trial),
                Err(error) => {
                    self.remember(error);
                    return false;
                }
            }
        }
        true
    }
    pub(super) fn final_evaluation(
        &self,
        q: &[f64],
        loads: &[[f64; 3]],
    ) -> Result<GasMembraneEvaluation, GasMembraneMotionError> {
        if !self.begin() {
            return Err(self.take_stop().expect("stopped operator"));
        }
        Ok(self.fresh(q, loads)?.value)
    }
}

impl<C: FnMut() -> bool> SecondOrderProblem for Operator<'_, '_, C> {
    fn dimension(&self) -> usize {
        self.masses.len()
    }
    fn validate_state(&self, q: &[f64], v: &[f64]) -> Result<(), String> {
        if v.len() != self.dimension() || v.iter().any(|v| !v.is_finite()) {
            return Err("finite dimension-matched motion velocity required".into());
        }
        if self.prepare(q) {
            Ok(())
        } else {
            Err(
                "original membrane/gas stage or endpoint refused; physical diagnostic retained"
                    .into(),
            )
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
        let trial = cache.as_ref().expect("prepared trial");
        // Trial has ZERO dead force. M*a + this residual = native forcing;
        // pressure is already present and is not applied a second time.
        for (out, &i) in output.iter_mut().zip(&self.model.free_dofs) {
            *out = trial.value.residual_n[i / 3][i % 3];
        }
    }
    fn tangent_apply(&self, q: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if output.len() != self.dimension() || !self.prepare(q) {
            return;
        }
        let cache = self.cache.borrow();
        match self.model.tangent_at(
            cache.as_ref().expect("prepared trial"),
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
        // Original mass diagonal supplies a cheap right preconditioner; the
        // complete nonlinear tangent remains in every native Krylov action.
        for ((out, value), mass) in output.iter_mut().zip(residual).zip(self.masses) {
            let scale = weights.mass * mass;
            if scale.is_finite() && scale > 0.0 {
                *out = value / scale;
            }
        }
    }
}
