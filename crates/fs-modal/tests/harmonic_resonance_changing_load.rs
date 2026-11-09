//! G1/G3: material and attachment changes move actual peaks, linewidths and Q.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{
    HarmonicTransferParameterJet, ResonanceControls,
    refine_half_power_bandwidth_with_parameter_derivative,
    refine_transfer_extremum_with_parameter_derivative,
};
use fs_modal::harmonic::{
    HarmonicFrequencyControls, HarmonicFrequencyParameterControls, HarmonicPortControls,
    HarmonicPortFrequencyParameterDirection, HarmonicPortParameterDirection,
    HarmonicResponseControls, HarmonicResponseError, PreparedHarmonicResponse,
};
use fs_modal::quadratic_modes::QuadraticParameterDirection;

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 73,
                kernel_id: 27,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 4,
        maximum_payload_bytes: 1 << 20,
        maximum_relative_residual: 2e-11,
    }
}
fn port_controls() -> HarmonicPortControls {
    HarmonicPortControls {
        maximum_ports: 4,
        maximum_projection_terms: 1_000_000,
        maximum_payload_bytes: 1 << 20,
    }
}
fn query_controls() -> ResonanceControls {
    ResonanceControls {
        absolute_angular_tolerance_rad_s: 1e-12,
        relative_angular_tolerance: 1e-12,
        maximum_evaluations: 300,
    }
}

/// Independent exact roots of the mechanical displacement half-height equation.
fn extrema(m: f64, c: f64, k: f64, dm: f64, dc: f64, dk: f64) -> [f64; 6] {
    let a = k / m - c * c / (2.0 * m * m);
    let da = dk / m - k * dm / (m * m) - c * dc / (m * m) + c * c * dm / (m * m * m);
    let r = k / m - c * c / (4.0 * m * m);
    let dr = dk / m - k * dm / (m * m) - c * dc / (2.0 * m * m) + c * c * dm / (2.0 * m * m * m);
    let u = c / m * r.sqrt();
    let du = (dc / m - c * dm / (m * m)) * r.sqrt() + c / m * dr / (2.0 * r.sqrt());
    let peak = a.sqrt();
    let lower = (a - u).sqrt();
    let upper = (a + u).sqrt();
    [
        peak,
        da / (2.0 * peak),
        lower,
        (da - du) / (2.0 * lower),
        upper,
        (da + du) / (2.0 * upper),
    ]
}

#[test]
fn changing_load_mass_loss_spring_and_attachment_follow_the_complete_physical_line_shape() {
    with_cx(|cx| {
        let (m, c, k, b, lm, lc, lk) = (
            1.3_f64, 0.12_f64, 7.0_f64, 1.4_f64, 0.25_f64, 0.08_f64, 0.7_f64,
        );
        let (effective_m, effective_c, effective_k) =
            (m + b * b * lm, c + b * b * lc, k + b * b * lk);
        let force = [C64::new(1.3, -0.7)];
        let zero = [0.0];
        for (db, dlm, dlc, dlk) in [
            (0.0, 1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0, 0.0),
            (0.0, 0.0, 0.0, 1.0),
            (1.0, 0.0, 0.0, 0.0),
        ] {
            let direction = QuadraticParameterDirection {
                mass: &zero,
                damping: &zero,
                stiffness: &zero,
            };
            let prepare = |w: f64| {
                let z = C64::new(lc, lk / w - w * lm);
                PreparedHarmonicResponse::try_with_port_impedance(
                    &[m],
                    &[c],
                    &[k],
                    1,
                    &[b],
                    1,
                    &[z],
                    w,
                    controls(),
                    port_controls(),
                    cx,
                )
            };
            let evaluate = |w: f64| {
                let zw = C64::new(0.0, -lk / (w * w) - lm);
                prepare(w)?
                    .solve_with_frequency_derivative(
                        &force,
                        None,
                        Some(&[zw]),
                        HarmonicFrequencyControls {
                            maximum_additional_payload_bytes: 1 << 20,
                            maximum_work: 1_000_000,
                        },
                        cx,
                    )?
                    .displacement_transfer(&[1.0])
            };
            let evaluate_jet =
                |w: f64| -> Result<HarmonicTransferParameterJet, HarmonicResponseError> {
                    let zw = C64::new(0.0, -lk / (w * w) - lm);
                    let zww = C64::new(0.0, 2.0 * lk / (w * w * w));
                    let zt = C64::new(dlc, dlk / w - w * dlm);
                    let zwt = C64::new(0.0, -dlk / (w * w) - dlm);
                    let response = prepare(w)?.solve_with_frequency_parameter_and_port_derivative(
                        &force,
                        direction,
                        None,
                        None,
                        None,
                        None,
                        &[zw],
                        &[zww],
                        HarmonicPortFrequencyParameterDirection {
                            parameter: HarmonicPortParameterDirection {
                                force_map: &[db],
                                impedance: &[zt],
                            },
                            impedance_angular_parameter_derivative: &[zwt],
                        },
                        HarmonicFrequencyParameterControls {
                            maximum_additional_payload_bytes: 1 << 20,
                            maximum_work: 1_000_000,
                        },
                        cx,
                    )?;
                    let power_scale = 1.0
                        + response.average_input_power_mixed_derivative.abs()
                        + response.average_viscous_power_mixed_derivative.abs()
                        + response
                            .port_load
                            .as_ref()
                            .unwrap()
                            .average_absorbed_power_mixed_derivative
                            .abs();
                    assert!(
                        response.power_balance_residual_mixed_derivative.abs()
                            < 3e-10 * power_scale
                    );
                    response.displacement_transfer_parameter(&[1.0])
                };
            let peak = refine_transfer_extremum_with_parameter_derivative(
                (1.8, 2.5),
                query_controls(),
                cx,
                evaluate,
                evaluate_jet,
            )
            .unwrap();
            let width = refine_half_power_bandwidth_with_parameter_derivative(
                &peak,
                (1.0, peak.extremum.angular_frequency_rad_s),
                (peak.extremum.angular_frequency_rad_s, 3.3),
                query_controls(),
                cx,
                evaluate,
                evaluate_jet,
            )
            .unwrap();
            let dm = 2.0 * b * db * lm + b * b * dlm;
            let dc = 2.0 * b * db * lc + b * b * dlc;
            let dk = 2.0 * b * db * lk + b * b * dlk;
            let [w, dw, lower, dlower, upper, dupper] =
                extrema(effective_m, effective_c, effective_k, dm, dc, dk);
            let bandwidth = upper - lower;
            let bandwidth_derivative = dupper - dlower;
            let quality_derivative =
                dw / bandwidth - w * bandwidth_derivative / (bandwidth * bandwidth);
            assert!((peak.extremum.angular_frequency_rad_s - w).abs() < 4e-12);
            assert!((peak.angular_frequency_derivative_rad_s - dw).abs() < 3e-9);
            assert!((width.bandwidth.angular_bandwidth_rad_s - bandwidth).abs() < 6e-12);
            assert!((width.lower_angular_frequency_derivative_rad_s - dlower).abs() < 3e-9);
            assert!((width.upper_angular_frequency_derivative_rad_s - dupper).abs() < 3e-9);
            assert!((width.angular_bandwidth_derivative_rad_s - bandwidth_derivative).abs() < 6e-9);
            assert!((width.response_quality_factor_derivative - quality_derivative).abs() < 1e-7);
        }
    });
}
