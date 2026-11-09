//! Signed field flux and reversible storage on one common temporal factor.
//!
//! For native EM amplitudes A exp(+i omega t), dA/dt = +i omega A.
//! A proper covariance permits removing one common carrier: Re(i omega_ref
//! |A|²)=0 exactly. Actual real fields require the full derivative frequency.
//! Storage is value times derivative, with its material coefficient; it is
//! neither a positive power nor irreversible heat. Reservoir fields keep the
//! original P/sqrt(Delta) and P_dot/(omega0 sqrt(Delta)) normalization.

use super::*;

type Modes = [[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];

pub(super) struct ProjectedLayerField {
    value: [Modes; 6],
    derivative: [Modes; 6],
}

impl PositiveLayerMoments<'_> {
    fn project_value_and_derivative(
        &self,
        field: impl Fn(usize) -> [Dd; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Modes, Modes), MomentError> {
        if cancelled() {
            return Err(MomentError::Cancelled);
        }
        let mut values = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
        let mut derivatives = values;
        for i in 0..self.frequencies {
            let [real, imaginary] = field(i);
            let omega = self.derivative_frequencies[i];
            // Conjugate the native positive-time value AND its +i*omega
            // derivative into the shared negative-time phase convention.
            values[i] = [real, -imaginary];
            derivatives[i] = [-imaginary * omega, -real * omega];
        }
        let mut value = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
        let mut derivative = value;
        let count = self.factor.mode_count();
        self.factor
            .project_parts_dd(&values[..self.frequencies], &mut value[..count])
            .map_err(MomentError::Product)?;
        self.factor
            .project_parts_dd(&derivatives[..self.frequencies], &mut derivative[..count])
            .map_err(MomentError::Product)?;
        if cancelled() {
            return Err(MomentError::Cancelled);
        }
        Ok((value, derivative))
    }

    fn cross_modes(
        &self,
        left: &Modes,
        right: &Modes,
        material_scale: Dd,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Dd, MomentError> {
        if cancelled() {
            return Err(MomentError::Cancelled);
        }
        let scale = material_scale
            * if self.proper {
                Dd::from_f64(0.5)
            } else {
                Dd::ONE
            };
        let mut moment = Dd::ZERO;
        for mode in 0..self.factor.mode_count() {
            moment = moment
                + self
                    .factor
                    .mode_cross_parts_dd_scaled(mode, left[mode], right[mode], scale)
                    .map_err(MomentError::Product)?[0];
        }
        if !moment.to_f64().is_finite() {
            return Err(MomentError::Unrepresentable);
        }
        Ok(moment)
    }

    pub(super) fn project_field(
        &self,
        field: impl Fn(usize) -> [C64; 6],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ProjectedLayerField, MomentError> {
        let mut projected = ProjectedLayerField {
            value: [[[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES]; 6],
            derivative: [[[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES]; 6],
        };
        for c in 0..6 {
            (projected.value[c], projected.derivative[c]) = self.project_value_and_derivative(
                |i| {
                    let value = field(i)[c];
                    [Dd::from_f64(value.re), Dd::from_f64(value.im)]
                },
                cancelled,
            )?;
        }
        Ok(projected)
    }

    pub(super) fn reservoir_moments(
        &self,
        field: impl Fn(usize) -> [C64; 3],
        transfer: impl Fn(usize) -> PoleTransfer,
        inertial: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Dd, Dd), MomentError> {
        let mut heat = Dd::ZERO;
        let mut storage = Dd::ZERO;
        for c in 0..3 {
            heat =
                heat + self.square(|i| filtered_parts(field(i)[c], transfer(i).heat), cancelled)?;
            let (value, derivative) = self.project_value_and_derivative(
                |i| filtered_parts(field(i)[c], transfer(i).storage),
                cancelled,
            )?;
            storage = storage + self.cross_modes(&derivative, &value, Dd::ONE, cancelled)?;
            if inertial {
                let (value, derivative) = self.project_value_and_derivative(
                    |i| filtered_parts(field(i)[c], transfer(i).kinetic_storage),
                    cancelled,
                )?;
                storage = storage + self.cross_modes(&derivative, &value, Dd::ONE, cancelled)?;
            }
        }
        if !heat.to_f64().is_finite() || heat.to_f64() < 0.0 || !storage.to_f64().is_finite() {
            return Err(MomentError::Unrepresentable);
        }
        Ok((heat, storage))
    }
}

impl ProjectedLayerField {
    pub(super) fn moments(
        &self,
        factor: &PositiveLayerMoments<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Dd; 9], MomentError> {
        let mut sum = [Dd::ZERO; 9];
        for c in 0..6 {
            sum[c] = factor.cross_modes(&self.value[c], &self.value[c], Dd::ONE, cancelled)?;
        }
        for c in 0..3 {
            let (u, v) = ((c + 1) % 3, (c + 2) % 3);
            sum[c + 6] =
                factor.cross_modes(&self.value[u], &self.value[v + 3], Dd::ONE, cancelled)?
                    - factor.cross_modes(&self.value[v], &self.value[u + 3], Dd::ONE, cancelled)?;
        }
        if sum.iter().any(|value| !value.to_f64().is_finite()) {
            return Err(MomentError::Unrepresentable);
        }
        Ok(sum)
    }

    pub(super) fn storage_rate(
        &self,
        factor: &PositiveLayerMoments<'_>,
        epsilon: f64,
        mu: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Dd, MomentError> {
        let mut rate = Dd::ZERO;
        for c in 0..6 {
            rate = rate
                + factor.cross_modes(
                    &self.derivative[c],
                    &self.value[c],
                    Dd::from_f64(if c < 3 { epsilon } else { mu }),
                    cancelled,
                )?;
        }
        if !rate.to_f64().is_finite() {
            return Err(MomentError::Unrepresentable);
        }
        Ok(rate)
    }
}
