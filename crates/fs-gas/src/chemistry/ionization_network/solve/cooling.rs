//! Optically thin electron–ion free-free cooling on the original atomic bank.
//!
//! The supplied law is frequency integrated. Its emitted energy is explicitly
//! assumed to escape; no photon-number moment or absorption spectrum is guessed.
//! Atomic bound-free photon accounts remain separate. At fixed populations,
//! Ce dTe/dt = -C sqrt(Te), with C = P/sqrt(Te). Freezing C at a predictor
//! midpoint gives the exact constant-Gaunt cooling solution and a second-order
//! approximation when the original Gaunt factor varies with temperature.
//! Every predictor and endpoint must still pass the original source's domains.

use super::*;
use fs_math::det::sqrt;
use fs_thermochem::plasma_radiation::{BremsstrahlungErrorV1, BremsstrahlungEvaluationV1};

impl Problem<'_> {
    /// Bind the exact bank and evaluate initial source admission even when the
    /// requested interval is zero or the initial electron population is absent.
    pub(super) fn prepare_cooling(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), IonizationNetworkError> {
        let Some(source) = self.cooling else {
            return Ok(());
        };
        if self.edges.is_empty() {
            return Err(IonizationNetworkError::Invalid(
                "free-free cooling requires an admitted atomic electron source",
            ));
        }
        if source.thermodynamics() != self.initial.thermodynamics()
            || self.electron != Some(source.electron_index())
        {
            return Err(IonizationNetworkError::Invalid(
                "free-free cooling requires the original bank and electron",
            ));
        }
        source
            .evaluate(self.initial, cancelled)
            .map_err(source_error)?;
        Ok(())
    }

    /// One full source evaluation is one globally bounded rate query, including
    /// work whose trial will subsequently be rejected or cancelled.
    fn cooling_stage(
        &self,
        source: &ElectronIonBremsstrahlungV1,
        amounts: &[Amount],
        th: f64,
        te: f64,
        balance: &mut IonizationNetworkBalance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BremsstrahlungEvaluationV1, IonizationNetworkError> {
        if balance.rate_evaluations >= self.controls.maximum_rate_evaluations {
            return Err(IonizationNetworkError::StepBudget);
        }
        balance.rate_evaluations += 1;
        balance.free_free_evaluations += 1;
        source
            .evaluate_point(
                amounts,
                self.initial.volume(),
                Temperature::new(th),
                Temperature::new(te),
                cancelled,
            )
            .map_err(source_error)
    }

    /// A symmetric outer composition applies this to the complete shared
    /// electron pool once per interval, independently of the atomic edge count.
    pub(super) fn cooling_step(
        &self,
        point: &mut Point,
        dt: f64,
        balance: &mut IonizationNetworkBalance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), IonizationNetworkError> {
        let Some(source) = self.cooling else {
            return Ok(());
        };
        if dt == 0.0 {
            return Ok(());
        }
        poll(cancelled)?;
        let mut amounts = [Amount::new(0.0); MAX_EQUILIBRIUM_SPECIES_V1];
        for (target, &value) in amounts.iter_mut().zip(&point.amounts) {
            *target = Amount::new(value);
        }
        let amounts = &amounts[..point.amounts.len()];
        let start = self.cooling_stage(source, amounts, point.th, point.te, balance, cancelled)?;
        let coefficient = start.total_power_over_sqrt_temperature_w_k_sqrt();
        if coefficient == 0.0 {
            return Ok(());
        }
        let electron = point.amounts[self.electron.expect("admitted cooling electron")];
        let capacity = CV * electron;
        if !positive(capacity) {
            return Err(IonizationNetworkError::Arithmetic(
                "free-free electron heat capacity",
            ));
        }
        let (predicted, _) = frozen_step(point.te, capacity, 0.5 * dt, coefficient)?;
        self.admit_temperature(predicted)?;
        let middle =
            self.cooling_stage(source, amounts, point.th, predicted, balance, cancelled)?;
        let (next, loss) = frozen_step(
            point.te,
            capacity,
            dt,
            middle.total_power_over_sqrt_temperature_w_k_sqrt(),
        )?;
        self.admit_temperature(next)?;
        // The actual endpoint may be cooler than the predictor. In particular,
        // degeneracy, stationary-ion and each populated ion's Gaunt bounds
        // must be checked here, without clipping or source extrapolation.
        self.cooling_stage(source, amounts, point.th, next, balance, cancelled)?;
        let escaped = point.free_free + Dd::from_f64(loss);
        if !nonnegative(escaped.to_f64()) {
            return Err(IonizationNetworkError::Arithmetic(
                "cumulative free-free escape",
            ));
        }
        point.te = next;
        point.free_free = escaped;
        Ok(())
    }
}

/// For constant C, sqrt(Te_new) = sqrt(Te_old) - C dt/(2 Ce). Compute
/// the emitted energy as C dt [sqrt(Te_old) - C dt/(4 Ce)], avoiding
/// subtraction of nearly equal thermal energies in a weak-cooling interval.
fn frozen_step(
    temperature: f64,
    capacity: f64,
    dt: f64,
    coefficient: f64,
) -> Result<(f64, f64), IonizationNetworkError> {
    if !positive(temperature)
        || !positive(capacity)
        || !nonnegative(dt)
        || !nonnegative(coefficient)
    {
        return Err(IonizationNetworkError::Arithmetic(
            "free-free frozen cooling coefficients",
        ));
    }
    if coefficient == 0.0 || dt == 0.0 {
        return Ok((temperature, 0.0));
    }
    let root = sqrt(temperature);
    let decrement = 0.5 * coefficient / capacity * dt;
    if !nonnegative(decrement) || decrement >= root {
        return Err(IonizationNetworkError::OutsideDomain);
    }
    let loss = coefficient * dt * (root - 0.5 * decrement);
    // Preserve small decrements relative to the original Te, and use the
    // positive square near strong cooling where Te - loss/Ce would cancel.
    let next = if decrement <= 0.5 * root {
        temperature - decrement * (2.0 * root - decrement)
    } else {
        let remaining = root - decrement;
        remaining * remaining
    };
    if !nonnegative(loss) || !positive(next) {
        return Err(IonizationNetworkError::OutsideDomain);
    }
    Ok((next, loss))
}

fn source_error(error: BremsstrahlungErrorV1) -> IonizationNetworkError {
    match error {
        BremsstrahlungErrorV1::Cancelled => IonizationNetworkError::Cancelled,
        BremsstrahlungErrorV1::Domain(_) => IonizationNetworkError::OutsideDomain,
        error => IonizationNetworkError::Bremsstrahlung(error),
    }
}
