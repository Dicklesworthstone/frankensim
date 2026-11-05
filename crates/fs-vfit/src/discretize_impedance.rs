//! Passive preparation of computed scalar driving-point impedances.
//!
//! The complete continuous model and the actual rounded digital coefficients
//! must both pass global admission before a filter is returned. Identification
//! accuracy and passivity-repair distortion remain separate observations.

use super::{
    DigitalFilter, DiscretizeError, RealizeError, bilinear, check_digital_passivity,
    validated_bilinear_parameters,
};
use crate::model::RationalModel;
use crate::passivity::{
    GlobalPassivityBudget, GlobalPassivityReport, GlobalPassivityVerdict, GlobalRepairReport,
    repair_passivity_global,
};
use crate::vf::{FitOptions, FitReport, vector_fit};
use fs_math::c64::C64;

/// Errors against the original samples, in the impedance's physical units.
///
/// These are sampled observations, not bounds between samples. The digital
/// error includes bilinear frequency warping as well as coefficient rounding;
/// the continuous error includes identification and passivity repair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImpedanceSampleErrors {
    /// Maximum absolute error of the final repaired continuous impedance.
    pub continuous_max_abs: f64,
    /// Maximum absolute error of the emitted filter at the supplied frequencies.
    pub digital_max_abs: f64,
}

/// A passive runtime filter with the physical approximation that produced it.
///
/// Reports describe these coefficients at return. Mutating the public model
/// or filter requires a new passivity check; these records are not authority
/// for arbitrary later edits. Preparation allocates; ordinary filter stepping
/// has no additional work.
#[derive(Debug, Clone)]
pub struct ImpedanceRealization {
    /// The admitted digital impedance, ready for its original runtime owner.
    pub filter: DigitalFilter,
    /// The final continuous impedance, after any residue-only repair.
    pub model: RationalModel,
    /// Identification diagnostics BEFORE passivity repair.
    pub fit: FitReport,
    /// Actual residue perturbation and full-frequency continuous certificate.
    pub repair: GlobalRepairReport,
    /// Full-circle admission of the rounded coefficients in filter.
    pub digital_certificate: GlobalPassivityReport,
    /// Errors of the final models against the original input samples.
    pub sample_errors: ImpedanceSampleErrors,
}

/// Fit computed impedance samples, enforce global passivity and discretize.
///
/// Stable poles or a clean frequency grid alone do not admit a feedback load.
/// Failed repair, an unresolved global proof, or a nonpassive rounded filter
/// refuses. No raw-fit fallback, global gain reduction or artificial resistance
/// is substituted. Uses the default finite passivity work allowance.
///
/// Conjugate data using the opposite time-harmonic convention before calling.
///
/// # Errors
/// Identification, passivity, sample-clock, Nyquist or representation failure.
pub fn realize_tabulated_impedance(
    omega: &[f64],
    z: &[C64],
    t_s: f64,
    opts: &FitOptions,
    prewarp: f64,
) -> Result<DigitalFilter, RealizeError> {
    realize_tabulated_impedance_with_report(
        omega,
        z,
        t_s,
        opts,
        prewarp,
        GlobalPassivityBudget::default(),
    )
    .map(|realized| realized.filter)
}

/// Prepare an impedance while retaining its final error and repair reports.
///
/// The supplied finite budget applies to each global check during the bounded
/// residue-repair rounds and to the final digital check. Poles, direct
/// resistance and inertance remain fixed during repair. Input sample order
/// does not change the declared band, and every requested digital frequency
/// must be strictly below Nyquist.
///
/// These certificates establish passivity of the represented operators. They
/// do not prove that the supplied physical data or the fit accuracy are valid;
/// the caller owns its physical approximation tolerance.
///
/// # Errors
/// As realize_tabulated_impedance; every refusal happens during preparation,
/// before a runtime filter or accepted filter history is published.
pub fn realize_tabulated_impedance_with_report(
    omega: &[f64],
    z: &[C64],
    t_s: f64,
    opts: &FitOptions,
    prewarp: f64,
    budget: GlobalPassivityBudget,
) -> Result<ImpedanceRealization, RealizeError> {
    let (_, nyquist, _) =
        validated_bilinear_parameters(t_s, prewarp).map_err(RealizeError::Discretize)?;
    for &frequency in omega {
        if frequency.is_finite() && frequency >= nyquist {
            return Err(RealizeError::Discretize(DiscretizeError::BeyondNyquist {
                omega: frequency,
                nyquist,
            }));
        }
    }
    let fit = vector_fit(omega, z, opts).map_err(RealizeError::Fit)?;
    let lo = omega.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = omega.iter().copied().fold(0.0_f64, f64::max);
    let (model, repair) =
        repair_passivity_global(&fit.model, (lo, hi), budget).map_err(RealizeError::Passivity)?;
    let filter = bilinear(&model, t_s, prewarp).map_err(RealizeError::Discretize)?;
    let digital_certificate =
        check_digital_passivity(&filter, budget).map_err(RealizeError::Discretize)?;
    if digital_certificate.verdict != GlobalPassivityVerdict::Certified {
        return Err(RealizeError::UncertifiedDigital {
            verdict: digital_certificate.verdict,
            frequency_rad_s: digital_certificate.violation_frequency_rad_s,
            reason: digital_certificate.unresolved,
        });
    }
    let sample_errors = final_sample_errors(omega, z, &model, &filter)?;
    Ok(ImpedanceRealization {
        filter,
        model,
        fit: fit.report,
        repair,
        digital_certificate,
        sample_errors,
    })
}

fn final_sample_errors(
    omega: &[f64],
    z: &[C64],
    model: &RationalModel,
    filter: &DigitalFilter,
) -> Result<ImpedanceSampleErrors, RealizeError> {
    let mut errors = ImpedanceSampleErrors {
        continuous_max_abs: 0.0,
        digital_max_abs: 0.0,
    };
    for (index, (&frequency, &expected)) in omega.iter().zip(z).enumerate() {
        let continuous = (model.eval_iw(frequency) - expected).abs();
        let digital = (filter.eval(frequency) - expected).abs();
        if !continuous.is_finite() || !digital.is_finite() {
            return Err(RealizeError::Discretize(
                DiscretizeError::NonFiniteRuntimeValue {
                    field: "impedance_sample_error",
                    index: Some(index),
                },
            ));
        }
        errors.continuous_max_abs = errors.continuous_max_abs.max(continuous);
        errors.digital_max_abs = errors.digital_max_abs.max(digital);
    }
    Ok(errors)
}
