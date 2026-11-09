//! Physical port receivers retain the actual transmitted-force/flow line shape.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{
    HarmonicReceiverParameterJet, HarmonicTransferParameterJet, ResonanceControls,
    refine_half_power_bandwidth_with_parameter_derivative, refine_transfer_extremum,
    refine_transfer_extremum_with_parameter_derivative,
};
use fs_modal::harmonic::{
    HarmonicFrequencyControls, HarmonicFrequencyParameterControls,
    HarmonicFrequencyParameterResponse, HarmonicPortControls,
    HarmonicPortFrequencyParameterDirection, HarmonicPortParameterDirection,
    HarmonicResponseControls, HarmonicResponseError, PreparedHarmonicResponse,
};
use fs_modal::quadratic_modes::QuadraticParameterDirection;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 101,
                kernel_id: 37,
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
        maximum_projection_terms: 100_000,
        maximum_payload_bytes: 1 << 20,
    }
}

fn frequency_controls() -> HarmonicFrequencyControls {
    HarmonicFrequencyControls {
        maximum_additional_payload_bytes: 1 << 20,
        maximum_work: 1_000_000,
    }
}

fn jet_controls() -> HarmonicFrequencyParameterControls {
    HarmonicFrequencyParameterControls {
        maximum_additional_payload_bytes: 1 << 20,
        maximum_work: 1_000_000,
    }
}

fn query_controls() -> ResonanceControls {
    ResonanceControls {
        absolute_angular_tolerance_rad_s: 1e-12,
        relative_angular_tolerance: 1e-12,
        maximum_evaluations: 300,
    }
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * expected.abs().max(1.0),
        "{actual:?} != {expected:?}"
    );
}

fn close_real(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * expected.abs().max(1.0),
        "{actual:?} != {expected:?}"
    );
}

fn damper_owner(w: f64, cx: &Cx<'_>) -> Result<PreparedHarmonicResponse, HarmonicResponseError> {
    PreparedHarmonicResponse::try_with_port_impedance(
        &[1.0],
        &[0.15],
        &[4.0],
        1,
        &[1.3],
        1,
        &[C64::from_re(0.4)],
        w,
        controls(),
        port_controls(),
        cx,
    )
}

fn damper_jet(
    w: f64,
    cx: &Cx<'_>,
) -> Result<HarmonicFrequencyParameterResponse, HarmonicResponseError> {
    damper_owner(w, cx)?.solve_with_frequency_parameter_and_port_derivative(
        &[C64::new(1.2, -0.7)],
        QuadraticParameterDirection {
            mass: &[0.0],
            damping: &[0.0],
            stiffness: &[0.0],
        },
        None,
        None,
        None,
        None,
        &[C64::ZERO],
        &[C64::ZERO],
        HarmonicPortFrequencyParameterDirection {
            parameter: HarmonicPortParameterDirection {
                force_map: &[0.0],
                impedance: &[C64::ONE],
            },
            impedance_angular_parameter_derivative: &[C64::ZERO],
        },
        jet_controls(),
        cx,
    )
}

#[test]
fn transmitted_force_resonance_and_width_follow_physical_load_loss() {
    with_cx(|cx| {
        let (b, r, c0) = (1.3_f64, 0.4_f64, 0.15_f64);
        let c = c0 + b * b * r;
        let force = C64::new(1.2, -0.7);
        let first = |w: f64| {
            damper_owner(w, cx)?
                .solve_with_frequency_derivative(
                    &[force],
                    None,
                    Some(&[C64::ZERO]),
                    frequency_controls(),
                    cx,
                )?
                .port_effort_transfer(&[1.0])
        };
        let full = |w: f64| damper_jet(w, cx)?.port_effort_transfer_parameter(&[1.0]);
        let peak = refine_transfer_extremum_with_parameter_derivative(
            (1.5, 2.5),
            query_controls(),
            cx,
            first,
            full,
        )
        .unwrap();
        let width = refine_half_power_bandwidth_with_parameter_derivative(
            &peak,
            (0.5, peak.extremum.angular_frequency_rad_s),
            (peak.extremum.angular_frequency_rad_s, 4.0),
            query_controls(),
            cx,
            first,
            full,
        )
        .unwrap();
        close_real(peak.extremum.angular_frequency_rad_s, 2.0);
        close_real(peak.angular_frequency_derivative_rad_s, 0.0);
        close(peak.extremum.sample.amplitude, force.scale(r * b / c));
        close(peak.amplitude_derivative, force.scale(b * c0 / (c * c)));
        close_real(peak.log_amplitude_derivative.unwrap(), 1.0 / r - b * b / c);
        close_real(peak.phase_derivative.unwrap(), 0.0);
        let root = (c * c + 16.0).sqrt();
        close_real(
            width.bandwidth.lower_angular_frequency_rad_s,
            0.5 * (root - c),
        );
        close_real(
            width.bandwidth.upper_angular_frequency_rad_s,
            0.5 * (root + c),
        );
        close_real(width.bandwidth.angular_bandwidth_rad_s, c);
        close_real(width.angular_bandwidth_derivative_rad_s, b * b);
        close_real(
            width.response_quality_factor_derivative,
            -2.0 * b * b / (c * c),
        );

        let displacement_peak = refine_transfer_extremum((1.5, 2.5), query_controls(), cx, |w| {
            damper_owner(w, cx)?
                .solve_with_frequency_derivative(
                    &[force],
                    None,
                    Some(&[C64::ZERO]),
                    frequency_controls(),
                    cx,
                )?
                .displacement_transfer(&[1.0])
        })
        .unwrap();
        close_real(
            displacement_peak.angular_frequency_rad_s,
            (4.0 - c * c / 2.0).sqrt(),
        );
        assert!(
            peak.extremum.angular_frequency_rad_s - displacement_peak.angular_frequency_rad_s
                > 0.08
        );

        // Independent original-pencil quotient verifies coherent rows and
        // frequency-dependent receiver phase, without conjugating weights.
        let w = 1.7;
        let d = C64::new(4.0 - w * w, -w * c);
        let q = force / d;
        let qw = force * C64::new(2.0 * w, c) / (d * d);
        let flow = C64::new(0.0, -w) * q.scale(b);
        let flow_w = C64::new(0.0, -1.0) * (q + qw.scale(w)).scale(b);
        let jet = damper_jet(w, cx).unwrap();
        let response = &jet.frequency_response.frequency_response;
        let weight = C64::new(0.7, -0.2);
        let weight_w = C64::new(-0.15, 0.4);
        close(response.port_flow_transfer(&[1.0]).unwrap().amplitude, flow);
        close(
            response
                .port_flow_transfer_complex(&[weight])
                .unwrap()
                .amplitude,
            weight * flow,
        );
        close(
            response
                .port_effort_transfer_complex(&[weight])
                .unwrap()
                .amplitude,
            (weight * flow).scale(r),
        );
        let observed_flow = response
            .port_flow_transfer_complex_with_derivative(&[weight], &[weight_w])
            .unwrap();
        let observed_effort = response
            .port_effort_transfer_complex_with_derivative(&[weight], &[weight_w])
            .unwrap();
        close(
            observed_flow.angular_derivative,
            weight_w * flow + weight * flow_w,
        );
        close(
            observed_effort.angular_derivative,
            (weight_w * flow + weight * flow_w).scale(r),
        );
        let fixed_flow = jet.port_flow_transfer_parameter(&[1.0]).unwrap();
        let complex_flow = jet.port_flow_transfer_parameter_complex(&[weight]).unwrap();
        let complex_effort = jet
            .port_effort_transfer_parameter_complex(&[weight])
            .unwrap();
        close(
            complex_flow.parameter_derivative,
            weight * fixed_flow.parameter_derivative,
        );
        close(
            complex_flow.angular_parameter_derivative,
            weight * fixed_flow.angular_parameter_derivative,
        );
        // e_theta = u + R*u_theta because the parameter is the physical damper R.
        close(
            complex_effort.parameter_derivative,
            weight * (flow + fixed_flow.parameter_derivative.scale(r)),
        );
        close(
            complex_effort.angular_parameter_derivative,
            weight * (flow_w + fixed_flow.angular_parameter_derivative.scale(r)),
        );
    });
}

fn maxwell_law(w: f64) -> (C64, C64, C64) {
    let (kb, tau) = (1.2_f64, 0.4_f64);
    let s = C64::new(1.0, -w * tau);
    (
        C64::from_re(0.1) + C64::from_re(kb * tau) / s,
        C64::new(0.0, kb * tau * tau) / (s * s),
        C64::from_re(-2.0 * kb * tau * tau * tau) / (s * s * s),
    )
}

fn maxwell_owner(w: f64, cx: &Cx<'_>) -> Result<PreparedHarmonicResponse, HarmonicResponseError> {
    PreparedHarmonicResponse::try_with_port_impedance(
        &[1.1],
        &[0.1],
        &[4.3],
        1,
        &[1.5],
        1,
        &[maxwell_law(w).0],
        w,
        controls(),
        port_controls(),
        cx,
    )
}

fn relabeled_port_jet(
    w: f64,
    cx: &Cx<'_>,
) -> Result<HarmonicFrequencyParameterResponse, HarmonicResponseError> {
    let (z, zw, zww) = maxwell_law(w);
    maxwell_owner(w, cx)?.solve_with_frequency_parameter_and_port_derivative(
        &[C64::new(1.2, -0.7)],
        QuadraticParameterDirection {
            mass: &[0.0],
            damping: &[0.0],
            stiffness: &[0.0],
        },
        None,
        None,
        None,
        None,
        &[zw],
        &[zww],
        HarmonicPortFrequencyParameterDirection {
            parameter: HarmonicPortParameterDirection {
                force_map: &[1.5],
                impedance: &[z.scale(-2.0)],
            },
            impedance_angular_parameter_derivative: &[zw.scale(-2.0)],
        },
        jet_controls(),
        cx,
    )
}

fn physical_force_jet(
    w: f64,
    cx: &Cx<'_>,
) -> Result<HarmonicTransferParameterJet, HarmonicResponseError> {
    let response = relabeled_port_jet(w, cx)?;
    HarmonicTransferParameterJet::project_complex_with_derivatives(
        HarmonicReceiverParameterJet {
            weights: &[C64::from_re(1.5)],
            angular_derivative: None,
            angular_second_derivative: None,
            parameter_derivative: Some(&[C64::from_re(1.5)]),
            angular_parameter_derivative: None,
        },
        response.port_effort_field_parameter_jet()?,
    )
}

#[test]
fn port_coordinate_change_preserves_mechanical_force_velocity_and_resonance() {
    with_cx(|cx| {
        // B(theta)=exp(theta)*B, Z(w,theta)=exp(-2*theta)*Z(w).
        // The full mechanical load is invariant. Flow and effort carry
        // reciprocal units; physical receiver weights must change with them.
        for w in [0.0, 0.7, 2.0, 6.0] {
            let response = relabeled_port_jet(w, cx).unwrap();
            let first = &response.frequency_response.frequency_response;
            close(
                response.parameter_response.displacement_derivative[0],
                C64::ZERO,
            );
            close(response.displacement_mixed_derivative[0], C64::ZERO);
            let flow = response.port_flow_field_parameter_jet().unwrap();
            let effort = response.port_effort_field_parameter_jet().unwrap();
            close(flow.parameter_derivative[0], flow.values[0]);
            close(
                flow.angular_parameter_derivative[0],
                flow.angular_derivative[0],
            );
            close(effort.parameter_derivative[0], -effort.values[0]);
            close(
                effort.angular_parameter_derivative[0],
                -effort.angular_derivative[0],
            );
            let force = physical_force_jet(w, cx).unwrap();
            close(force.sample.amplitude, effort.values[0].scale(1.5));
            close(force.parameter_derivative, C64::ZERO);
            close(force.angular_parameter_derivative, C64::ZERO);
            let velocity = HarmonicTransferParameterJet::project_complex_with_derivatives(
                HarmonicReceiverParameterJet {
                    weights: &[C64::from_re(1.0 / 1.5)],
                    angular_derivative: None,
                    angular_second_derivative: None,
                    parameter_derivative: Some(&[C64::from_re(-1.0 / 1.5)]),
                    angular_parameter_derivative: None,
                },
                flow,
            )
            .unwrap();
            close(velocity.sample.amplitude, first.response.velocity[0]);
            close(velocity.parameter_derivative, C64::ZERO);
            close(velocity.angular_parameter_derivative, C64::ZERO);
            close_real(
                response
                    .port_load
                    .as_ref()
                    .unwrap()
                    .average_absorbed_power_mixed_derivative,
                0.0,
            );
            close_real(response.power_balance_residual_mixed_derivative, 0.0);
        }
        let peak = refine_transfer_extremum_with_parameter_derivative(
            (1.5, 3.0),
            query_controls(),
            cx,
            |w| {
                let (_, zw, _) = maxwell_law(w);
                maxwell_owner(w, cx)?
                    .solve_with_frequency_derivative(
                        &[C64::new(1.2, -0.7)],
                        None,
                        Some(&[zw]),
                        frequency_controls(),
                        cx,
                    )?
                    .port_effort_transfer(&[1.5])
            },
            |w| physical_force_jet(w, cx),
        )
        .unwrap();
        close_real(peak.angular_frequency_derivative_rad_s, 0.0);
        close(peak.amplitude_derivative, C64::ZERO);
        close_real(peak.log_amplitude_derivative.unwrap(), 0.0);
        close_real(peak.phase_derivative.unwrap(), 0.0);
    });
}

#[test]
fn port_receivers_refuse_missing_physical_jets_and_invalid_weights() {
    with_cx(|cx| {
        let owner =
            PreparedHarmonicResponse::try_new(&[1.0], &[0.15], &[4.0], 1, 1.7, controls(), cx)
                .unwrap();
        let response = owner
            .solve_with_frequency_parameter_derivative(
                &[C64::ONE],
                QuadraticParameterDirection {
                    mass: &[0.0],
                    damping: &[0.0],
                    stiffness: &[0.0],
                },
                None,
                None,
                None,
                None,
                jet_controls(),
                cx,
            )
            .unwrap();
        let first = &response.frequency_response.frequency_response;
        assert!(first.port_flow_transfer(&[1.0]).is_err());
        assert!(first.port_effort_transfer_complex(&[C64::ONE]).is_err());
        assert!(response.port_flow_field_parameter_jet().is_err());
        assert!(response.port_effort_transfer_parameter(&[1.0]).is_err());

        let mut loaded = damper_jet(1.7, cx).unwrap();
        assert!(loaded.port_flow_transfer_parameter(&[]).is_err());
        assert!(
            loaded
                .port_effort_transfer_parameter_complex(&[C64::new(f64::NAN, 0.0)])
                .is_err()
        );
        assert!(
            loaded
                .frequency_response
                .frequency_response
                .port_effort_transfer_complex_with_derivative(&[C64::ONE], &[])
                .is_err()
        );
        loaded.port_load = None;
        assert!(loaded.port_effort_field_parameter_jet().is_err());
    });
}
