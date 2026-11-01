//! Numerical adapters for the single geometry event walk.

use super::*;

// The driver only sees the optical prefix. An augmented stepper retains its
// derivative columns and includes their local errors in the same step decision.
pub(crate) trait EventStepper {
    type Trial;

    fn state(&self) -> State;
    fn initial_sample(&self) -> GradientIndexSample;
    fn final_sample(&self) -> GradientIndexSample;
    fn is_uniform(&self) -> bool;
    fn trial_state(trial: &Self::Trial) -> State;
    fn trial(
        &self,
        step: f64,
        new_distance: f64,
        guard: &mut impl FnMut(Point3, Vec3) -> bool,
        evaluations: &mut u32,
    ) -> Result<Self::Trial, GradientIndexError>;
    fn normalized_error(
        &self,
        trial: &Self::Trial,
        step: f64,
        options: GradientIndexOptions,
    ) -> Result<f64, GradientIndexError>;
    fn accept(&mut self, trial: Self::Trial) -> Result<(), GradientIndexError>;
}

pub(super) struct OpticalStepper<'a, F: GradientIndexField + ?Sized> {
    field: &'a F,
    ray: Ray,
    wavelength_nm: f64,
    initial: GradientIndexSample,
    sample: GradientIndexSample,
    uniform: bool,
    track_group: bool,
    state: State,
    derivative: State,
}

impl<'a, F: GradientIndexField + ?Sized> OpticalStepper<'a, F> {
    pub(super) fn new(
        field: &'a F,
        ray: Ray,
        axis: Vec3,
        wavelength_nm: f64,
    ) -> Result<Self, GradientIndexError> {
        let initial = field.sample(ray.origin, wavelength_nm)?;
        let uniform = field.is_spatially_uniform();
        if uniform && initial.phase_gradient_per_m != Vec3::new(0.0, 0.0, 0.0) {
            return Err(GradientIndexError::InvalidSample);
        }
        let track_group = initial.group_index.is_some();
        let mut state = [0.0; 12];
        put_vec(&mut state, 3, ray.dir);
        put_vec(&mut state, 6, axis);
        Ok(Self {
            field,
            ray,
            wavelength_nm,
            initial,
            sample: initial,
            uniform,
            track_group,
            derivative: rhs(&state, initial, track_group)?,
            state,
        })
    }
}

impl<F: GradientIndexField + ?Sized> EventStepper for OpticalStepper<'_, F> {
    type Trial = Trial;

    fn state(&self) -> State {
        self.state
    }
    fn initial_sample(&self) -> GradientIndexSample {
        self.initial
    }
    fn final_sample(&self) -> GradientIndexSample {
        self.sample
    }
    fn is_uniform(&self) -> bool {
        self.uniform
    }
    fn trial_state(trial: &Self::Trial) -> State {
        trial.state
    }
    fn trial(
        &self,
        step: f64,
        new_distance: f64,
        guard: &mut impl FnMut(Point3, Vec3) -> bool,
        evaluations: &mut u32,
    ) -> Result<Trial, GradientIndexError> {
        if !self.uniform {
            return trial_step_guarded(
                self.field,
                self.ray.origin,
                self.wavelength_nm,
                &self.state,
                self.derivative,
                step,
                self.track_group,
                guard,
                evaluations,
            );
        }
        let mut next = self.state;
        put_vec(&mut next, 0, self.ray.dir.scale(new_distance));
        next[9] = self.initial.phase_index * new_distance;
        next[10] = self.initial.group_index.unwrap_or(0.0) * new_distance;
        next[11] = self.initial.absorption_per_m * new_distance;
        let point = self.ray.origin.offset(get_vec(&next, 0));
        if next.iter().any(|value| !value.is_finite()) || !finite_point(point) {
            return Err(GradientIndexError::NumericalRange);
        }
        if !guard(point, self.ray.dir) {
            return Err(GradientIndexError::BoundaryClearanceExceeded);
        }
        Ok(Trial {
            state: next,
            error: [0.0; 12],
            projection_error: 0.0,
            sample: self.initial,
        })
    }
    fn normalized_error(
        &self,
        trial: &Trial,
        step: f64,
        options: GradientIndexOptions,
    ) -> Result<f64, GradientIndexError> {
        normalized_error(&self.state, trial, step, options)
    }
    fn accept(&mut self, trial: Trial) -> Result<(), GradientIndexError> {
        self.state = trial.state;
        self.sample = trial.sample;
        self.derivative = rhs(&self.state, self.sample, self.track_group)?;
        Ok(())
    }
}
