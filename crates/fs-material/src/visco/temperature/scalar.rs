//! The same calibrated time shift in existing scalar constitutive consumers.

use super::{TemperatureShift, TemperatureShiftState, invalid, positive, scaled};
use crate::visco::{
    FractionalZener, GeneralizedMaxwell, LoweredModel, ModifiedFractionalZener, ViscoError,
};

/// Existing constitutive law at one temperature, retaining its exact calibration.
/// The original numerical owner supplies harmonic response, history and heat;
/// this object does not introduce a second evolution or fitting algorithm.
#[derive(Clone, Debug, PartialEq)]
pub struct ShiftedViscoelastic<M> {
    pub(super) law: M,
    pub(super) state: TemperatureShiftState,
}

impl<M> ShiftedViscoelastic<M> {
    /// Original constitutive implementation with shifted viscous kinetics.
    pub const fn law(&self) -> &M {
        &self.law
    }

    /// Actual temperature, shift factor, derivative and calibration declaration.
    pub const fn state(&self) -> &TemperatureShiftState {
        &self.state
    }
}

impl TemperatureShift {
    /// Shift the existing scalar Maxwell spectrum, preserving branch order,
    /// every elastic modulus and its exact strain-ramp/heat implementation.
    /// Apply this to the original reference-temperature law on every call.
    ///
    /// The scalar law carries no frequency support or source identity. The
    /// caller must retain those qualifications; a lowered source approximation
    /// should use [Self::shift_lowered] so its support travels with the law.
    /// No branch fitting, reduction or per-step allocation is introduced.
    ///
    /// # Errors
    /// Invalid original public fields, temperature outside calibration, or
    /// overflow/underflow of a positive shifted relaxation time.
    pub fn shift_maxwell(
        &self,
        reference: &GeneralizedMaxwell,
        temperature_k: f64,
    ) -> Result<ShiftedViscoelastic<GeneralizedMaxwell>, ViscoError> {
        let state = self.evaluate(temperature_k)?;
        // These fields are public: admission cannot assume an earlier new().
        let mut law = GeneralizedMaxwell::new(reference.e_inf, reference.terms.clone())?;
        for (_, tau) in &mut law.terms {
            *tau = scaled(*tau, state.factor())?;
        }
        Ok(ShiftedViscoelastic { law, state })
    }

    /// Shift the original fractional Zener relaxation time. The same harmonic
    /// and automatic-differentiation owner remains in use; elastic moduli and
    /// fractional order are unchanged. Source frequency support stays with
    /// the caller because [FractionalZener] does not contain that metadata.
    pub fn shift_fractional_zener(
        &self,
        reference: &FractionalZener,
        temperature_k: f64,
    ) -> Result<ShiftedViscoelastic<FractionalZener>, ViscoError> {
        let state = self.evaluate(temperature_k)?;
        let reference = FractionalZener::new(
            reference.e0,
            reference.e_inf,
            reference.alpha,
            reference.tau,
        )?;
        let law = FractionalZener::new(
            reference.e0,
            reference.e_inf,
            reference.alpha,
            scaled(reference.tau, state.factor())?,
        )?;
        Ok(ShiftedViscoelastic { law, state })
    }

    /// Shift the bounded two-exponent Zener law with its complete original
    /// frequency band. Both exponents and both elastic fit scales are retained.
    /// Existing checked evaluation and positive-Prony lowering use the shifted
    /// law directly; nothing is refitted to temperature-dependent output.
    pub fn shift_modified_zener(
        &self,
        reference: &ModifiedFractionalZener,
        temperature_k: f64,
    ) -> Result<ShiftedViscoelastic<ModifiedFractionalZener>, ViscoError> {
        let state = self.evaluate(temperature_k)?;
        let law = ModifiedFractionalZener::new(
            reference.equilibrium_modulus_pa(),
            reference.fit_modulus_pa(),
            reference.alpha(),
            reference.beta(),
            scaled(reference.relaxation_time_s(), state.factor())?,
            shifted_band(reference.band_hz(), state.factor())?,
        )?;
        Ok(ShiftedViscoelastic { law, state })
    }

    /// Shift an already lowered Prony approximation without repeating NNLS.
    ///
    /// The source and approximation use the SAME frequency map:
    /// omega_ref = a_T * omega. The original measured fit errors and sample
    /// count are retained, now referring to that mapped sample grid. They are
    /// not a new material validation or an inter-sample error bound.
    /// Checked evaluation refuses outside [omega_lo/a_T, omega_hi/a_T].
    ///
    /// # Errors
    /// Invalid original fit metadata or unrepresentable transformed support,
    /// in addition to the original scalar-law and temperature refusals.
    pub fn shift_lowered(
        &self,
        reference: &LoweredModel,
        temperature_k: f64,
    ) -> Result<ShiftedViscoelastic<LoweredModel>, ViscoError> {
        if reference.verification_points < 2
            || [
                reference.sup_rel_err,
                reference.sup_rel_storage_err,
                reference.sup_rel_loss_err,
            ]
            .iter()
            .any(|error| !error.is_finite() || *error < 0.0)
        {
            return Err(invalid("original lowered-model fit metadata is invalid"));
        }
        let shifted = self.shift_maxwell(&reference.model, temperature_k)?;
        let band = shifted_band(reference.band, shifted.state.factor())?;
        Ok(ShiftedViscoelastic {
            law: LoweredModel {
                model: shifted.law,
                band,
                sup_rel_err: reference.sup_rel_err,
                sup_rel_storage_err: reference.sup_rel_storage_err,
                sup_rel_loss_err: reference.sup_rel_loss_err,
                verification_points: reference.verification_points,
            },
            state: shifted.state,
        })
    }
}

fn shifted_band(band: (f64, f64), factor: f64) -> Result<(f64, f64), ViscoError> {
    let next = (band.0 / factor, band.1 / factor);
    if !positive(band.0)
        || !positive(band.1)
        || band.0 >= band.1
        || !positive(next.0)
        || !positive(next.1)
        || next.0 >= next.1
    {
        return Err(invalid(
            "shifted scalar frequency support is unrepresentable",
        ));
    }
    Ok(next)
}
