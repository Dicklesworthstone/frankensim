use super::super::{ComplexSum, PreparedFrequency, PreparedLine};
use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

#[test]
fn triangular_native_hz_and_radian_nulls_keep_all_signed_field_moments() {
    let duration = 1e-12_f64;
    let window =
        HarmonicWindow::with_profile(-0.5 * duration, duration, HarmonicWindowProfile::Triangular)
            .unwrap();
    for angular in [false, true] {
        for proper in [false, true] {
            let carrier = if proper { 2.0_f64.powi(49) } else { 1.0 };
            let field = PreparedHarmonicField::<2> {
                lines: [1.0, -2.0, 1.0]
                    .into_iter()
                    .enumerate()
                    .map(|(i, amplitude)| {
                        let component = |gain| ComplexSum {
                            re: Dd::from_f64(amplitude * gain),
                            im: Dd::ZERO,
                        };
                        PreparedLine {
                            frequency: if angular {
                                PreparedFrequency::Angular(carrier + i as f64)
                            } else {
                                PreparedFrequency::Hertz(carrier + i as f64)
                            },
                            pressure: component(1.0),
                            velocity: [component(-2.0), component(0.5)],
                        }
                    })
                    .collect(),
            };
            let mut sums = FieldMoments::new();
            assert!(
                accumulate(
                    &field,
                    window,
                    if proper {
                        FactorKind::Proper
                    } else {
                        FactorKind::Deterministic
                    },
                    !proper,
                    true,
                    &mut sums,
                    &mut || false
                )
                .unwrap()
            );
            let span = duration * if angular { 1.0 } else { core::f64::consts::TAU };
            let square = span.powi(4) / if proper { 480.0 } else { 240.0 };
            let mean = if proper { 0.0 } else { -span.powi(2) / 24.0 };
            for (actual, expected) in [
                sums.pressure_mean,
                sums.pressure_square,
                sums.intensity[0],
                sums.velocity[0][1],
            ]
            .into_iter()
            .zip([mean, square, -2.0 * square, -square])
            {
                assert!(
                    (actual.to_f64() - expected).abs() < 2e-9 * expected.abs().max(1e-300),
                    "{actual:?} != {expected:e}"
                );
            }
        }
    }
}
