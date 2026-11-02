//! One physical fundamental from all recorded Cartesian force components.

use fs_fft::harmonic::{
    HarmonicFitError, HarmonicFitSummary, HarmonicFrequencyBudget, HarmonicFrequencyError,
    refine_harmonic_frequency_channels,
};
use fs_math::c64::C64;

/// Numerical diagnostics accompanying a refined physical force spectrum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ForceFrequencyRefinement {
    /// One locally refined fundamental shared by all recorded axes [Hz].
    pub fundamental_hz: f64,
    /// Common-frequency QR candidates evaluated, including rejected trials.
    pub evaluations: usize,
    /// Remaining projected frequency correction [Hz], not an uncertainty bound.
    pub last_correction_hz: f64,
}

pub(super) fn refine_force_harmonics<const D: usize>(
    force: &[[f64; D]],
    timestep_s: f64,
    initial_fundamental_hz: f64,
    harmonics: &[usize],
    budget: HarmonicFrequencyBudget,
) -> Result<
    (
        Vec<[C64; D]>,
        [HarmonicFitSummary; D],
        ForceFrequencyRefinement,
    ),
    HarmonicFrequencyError,
> {
    let invalid = |what| HarmonicFrequencyError::InvalidParameter { what };
    if D == 0 || force.is_empty() {
        return Err(invalid(
            "force refinement needs nonempty recorded components",
        ));
    }
    // Reserve the structure-of-arrays input and final force/order arrays in
    // addition to the shared estimator's complete internal peak. The original
    // accepted force record is caller-owned and excluded from this budget.
    let extra = force
        .len()
        .checked_mul(D)
        .and_then(|n| n.checked_mul(8))
        .and_then(|n| {
            harmonics
                .len()
                .checked_mul(
                    D.checked_mul(16)?
                        .checked_add(core::mem::size_of::<usize>())?,
                )
                .and_then(|k| n.checked_add(k))
        })
        .ok_or_else(|| invalid("force refinement workspace size overflows"))?;
    if extra > budget.fit.max_workspace_bytes {
        return Err(HarmonicFrequencyError::Fit(
            HarmonicFitError::BudgetExceeded {
                required_bytes: extra,
                max_bytes: budget.fit.max_workspace_bytes,
            },
        ));
    }
    let channels: [Vec<f64>; D] =
        core::array::from_fn(|axis| force.iter().map(|row| row[axis]).collect());
    let views: [&[f64]; D] = core::array::from_fn(|axis| channels[axis].as_slice());
    let mut inner_budget = budget;
    inner_budget.fit.max_workspace_bytes -= extra;
    let fit = refine_harmonic_frequency_channels(
        &views,
        timestep_s,
        initial_fundamental_hz,
        harmonics,
        inner_budget,
    )
    .map_err(|error| match error {
        HarmonicFrequencyError::Fit(HarmonicFitError::BudgetExceeded {
            required_bytes, ..
        }) => match required_bytes.checked_add(extra) {
            Some(required_bytes) => HarmonicFrequencyError::Fit(HarmonicFitError::BudgetExceeded {
                required_bytes,
                max_bytes: budget.fit.max_workspace_bytes,
            }),
            None => invalid("force refinement workspace size overflows"),
        },
        other => other,
    })?;
    let phasors = (0..harmonics.len())
        .map(|tone| {
            core::array::from_fn(|axis| {
                let value = fit.fits[axis].phasors[tone];
                C64::new(value.re, value.im)
            })
        })
        .collect();
    let summary = core::array::from_fn(|axis| HarmonicFitSummary {
        dc: fit.fits[axis].dc,
        residual_rms: fit.fits[axis].residual_rms,
    });
    Ok((
        phasors,
        summary,
        ForceFrequencyRefinement {
            fundamental_hz: fit.fundamental_hz,
            evaluations: fit.evaluations,
            last_correction_hz: fit.last_correction_hz,
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::JetForceHarmonics;
    use fs_fft::harmonic::HarmonicFitBudget;
    use fs_material::fluid_acoustics::FluidAcousticState;

    fn budget() -> HarmonicFrequencyBudget {
        HarmonicFrequencyBudget {
            frequency_bounds_hz: [211.0, 215.0],
            absolute_tolerance_hz: 1e-8,
            max_evaluations: 24,
            fit: HarmonicFitBudget {
                max_workspace_bytes: 4 * 1024 * 1024,
                min_reciprocal_condition: 1e-10,
            },
        }
    }

    #[test]
    fn g1_refined_vector_force_retains_missing_fundamental_phase_and_radiation() {
        let f0 = 213.17;
        let dt = 1.0 / 8000.0;
        let expected = [
            [C64::ZERO; 2],
            [C64::new(1e-5, -2e-5), C64::new(0.7, -0.2)],
            [C64::ZERO, C64::new(0.17, 0.31)],
        ];
        let dc = [0.4, -0.8];
        let force: Vec<_> = (0..1021)
            .map(|i| {
                let mut row = dc;
                for (h, amplitude) in expected.iter().enumerate() {
                    let (s, c) =
                        (core::f64::consts::TAU * f0 * (h + 1) as f64 * i as f64 * dt).sin_cos();
                    for axis in 0..2 {
                        row[axis] += amplitude[axis].re * c + amplitude[axis].im * s;
                    }
                }
                row
            })
            .collect();
        let (phasors, summary, refined) =
            refine_force_harmonics(&force, dt, 213.7, &[1, 2, 3], budget()).unwrap();
        assert!((refined.fundamental_hz - f0).abs() < 1e-8);
        assert!(refined.evaluations > 1 && refined.evaluations <= 24);
        assert!(refined.last_correction_hz.abs() <= budget().absolute_tolerance_hz);
        let fitted = JetForceHarmonics {
            acoustic: FluidAcousticState::try_new(1.7, 277.0).unwrap(),
            harmonics: vec![1, 2, 3],
            fundamental_hz: refined.fundamental_hz,
            force_n_m: phasors,
            dc_force_n_m: summary.map(|s| s.dc),
            residual_rms_n_m: summary.map(|s| s.residual_rms),
        };
        for axis in 0..2 {
            assert!((fitted.dc_force_n_m[axis] - dc[axis]).abs() < 1e-9);
            assert!(fitted.residual_rms_n_m[axis] < 1e-9);
        }
        for (h, amplitude) in expected.iter().enumerate() {
            let observed = fitted
                .observe_harmonic(h + 1, [0.4, 0.3], [0.0, 0.0])
                .unwrap();
            assert!((observed.frequency_hz - (h + 1) as f64 * f0).abs() < 3e-8);
            for axis in 0..2 {
                assert!((observed.force_n_m[axis] - amplitude[axis]).abs() < 1e-9);
            }
            let reference = crate::curle2d::dipole_pressure(
                *amplitude,
                core::f64::consts::TAU * (h + 1) as f64 * f0 / 277.0,
                [0.4, 0.3],
                [0.0, 0.0],
            )
            .unwrap();
            assert!((observed.field.pressure - reference.pressure).abs() < 1e-8);
            assert_eq!(observed.field.scope, crate::SCOPE_STATEMENT);
        }
    }

    #[test]
    fn g0_refinement_counts_transposed_record_and_refuses_nonfinite_axes() {
        let force: Vec<_> = (0..128)
            .map(|i| {
                [
                    1.0,
                    (core::f64::consts::TAU * 213.17 * i as f64 / 8000.0).cos(),
                ]
            })
            .collect();
        let extra = 128 * 2 * 8 + 2 * 16 + core::mem::size_of::<usize>();
        let mut limits = budget();
        limits.fit.max_workspace_bytes = extra - 1;
        assert!(matches!(
            refine_force_harmonics(&force, 1.0 / 8000.0, 213.7, &[1], limits),
            Err(HarmonicFrequencyError::Fit(HarmonicFitError::BudgetExceeded {
                required_bytes,
                max_bytes,
            })) if required_bytes == extra && max_bytes == extra - 1
        ));
        limits.fit.max_workspace_bytes = extra + 1;
        assert!(matches!(
            refine_force_harmonics(&force, 1.0 / 8000.0, 213.7, &[1], limits),
            Err(HarmonicFrequencyError::Fit(HarmonicFitError::BudgetExceeded {
                required_bytes,
                max_bytes,
            })) if required_bytes > extra + 1 && max_bytes == extra + 1
        ));
        let mut bad = force;
        bad[25][1] = f64::NAN;
        assert!(matches!(
            refine_force_harmonics(&bad, 1.0 / 8000.0, 213.7, &[1], budget()),
            Err(HarmonicFrequencyError::Fit(
                HarmonicFitError::NonFinite { .. }
            ))
        ));
        assert!(
            refine_force_harmonics(&[[1.0, 2.0]; 128], 1.0 / 8000.0, 213.7, &[1], budget())
                .is_err()
        );
    }
}
