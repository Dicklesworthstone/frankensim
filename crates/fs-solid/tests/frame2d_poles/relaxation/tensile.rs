//! G1/G3/G4: actual per-layer Maxwell stress, constrained memory and refusals.
//! Manufactured coefficients test the native history observation, not rupture.

use super::*;
use fs_material::section::strength::AxialTensileCapacity;
use fs_phs::linear_exact::DenseExponentialBudget;
use fs_solid::frame2d::coupled::transient::FrameTransientControls;
use fs_solid::frame2d::relaxation::transient::FrameRelaxationInitialHistory;
use fs_solid::frame2d::strength::FrameTensileError;

const BROADBAND: (f64, f64) = (0.0, f64::MAX);

fn transient_controls() -> FrameTransientControls {
    FrameTransientControls {
        exact: DenseExponentialBudget {
            max_dimension: 128,
            max_dense_bytes: 64 << 20,
            max_multiply_terms: 1 << 42,
        },
        maximum_relative_balance_error: 1e-10,
    }
}

fn capacity(frame: &Frame2d) -> AxialTensileCapacity {
    frame
        .material_section(0)
        .unwrap()
        .axial_tensile_capacity("tensile_limit")
        .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn g1_equal_motion_keeps_distinct_layer_histories_and_signed_stress_relief() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let (frame, selections) = frame(
            &[&[(2.0, 0.5)], &[(20.0, 2.0)]],
            BROADBAND,
            true,
            (1.0, 0.0),
        );
        let capacities = [capacity(&frame)];
        let prepared = frame
            .prepare_relaxing_resonance(NO_RAYLEIGH, &selections, None, controls(), cx)
            .unwrap();
        let mut runtime = prepared
            .prepare_transient(0.1, transient_controls(), cx)
            .unwrap();
        let mut displacement = [0.0; 6];
        displacement[3] = 0.08;
        runtime
            .initialize(
                0.3,
                &displacement,
                &[0.0; 6],
                FrameRelaxationInitialHistory::Relaxed,
                &[],
                cx,
            )
            .unwrap();
        let relaxed = runtime.assess_tensile_strength(&capacities, 4, cx).unwrap();
        assert_eq!(relaxed.time_s, 0.3);
        assert!(!relaxed.reaches_limit());
        for sample in relaxed.samples() {
            near(sample.mean_normal_stress_pa, 8.0 * 0.08);
        }
        runtime
            .initialize(
                0.3,
                &displacement,
                &[0.0; 6],
                FrameRelaxationInitialHistory::Unrelaxed,
                &[],
                cx,
            )
            .unwrap();
        let before = runtime.canonical_state().to_vec();
        let energy = runtime.report().motion.total_energy_j;
        let unrelaxed = runtime.assess_tensile_strength(&capacities, 4, cx).unwrap();
        assert_eq!(unrelaxed.samples().len(), 4);
        assert!(unrelaxed.reaches_limit());
        assert_eq!(unrelaxed.governing_sample().unwrap().layer, 1);
        assert_eq!(unrelaxed.governing_sample().unwrap().end, 0);
        for sample in unrelaxed.samples() {
            let expected = (8.0 + [2.0, 20.0][sample.layer]) * 0.08;
            near(sample.mean_normal_stress_pa, expected);
            near(sample.assessment.maximum_axial_stress_pa, expected);
            near(sample.assessment.tensile_utilization, expected);
            assert!(core::ptr::eq(sample.capacity, &capacities[0]));
        }
        assert_eq!(runtime.canonical_state(), before);
        assert_eq!(runtime.report().motion.total_energy_j, energy);

        // Oppositely signed old arm stress relieves the SAME positive current
        // displacement. An absolute arm envelope would falsely report failure.
        let extensions = [-0.32, 0.0, 0.0, -0.032, 0.0, 0.0];
        runtime
            .initialize(
                0.3,
                &displacement,
                &[0.0; 6],
                FrameRelaxationInitialHistory::ElasticExtensions(&extensions),
                &[],
                cx,
            )
            .unwrap();
        let relieved = runtime.assess_tensile_strength(&capacities, 4, cx).unwrap();
        assert!(!relieved.reaches_limit());
        for sample in relieved.samples() {
            near(sample.mean_normal_stress_pa, 0.0);
        }
        assert!(runtime.report().material_memory_energy_j > 0.0);
    });
}

#[test]
fn g1_constrained_bending_memory_has_tensile_stress_and_decays_at_its_original_rate() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let (frame, selections) = frame(&[&[(10.0, 0.5)]], BROADBAND, true, (1.0, 0.0));
        let capacities = [capacity(&frame)];
        let prepared = frame
            .prepare_relaxing_resonance(NO_RAYLEIGH, &selections, None, controls(), cx)
            .unwrap();
        let duration = 1.0;
        let mut runtime = prepared
            .prepare_transient(duration, transient_controls(), cx)
            .unwrap();
        runtime
            .initialize(
                0.0,
                &[0.0; 6],
                &[0.0; 6],
                FrameRelaxationInitialHistory::ElasticExtensions(&[0.0, 0.2, -0.3]),
                &[],
                cx,
            )
            .unwrap();
        let first = runtime.assess_tensile_strength(&capacities, 2, cx).unwrap();
        // sigma(y) = delta_E * (3*(1-2*xi)*0.2 - 0.3) * y / L^2.
        near(first.samples()[0].normal_stress_gradient_pa_per_m, 3.0);
        near(first.samples()[1].normal_stress_gradient_pa_per_m, -9.0);
        near(first.samples()[0].assessment.tensile_utilization, 0.6);
        near(first.samples()[1].assessment.tensile_utilization, 1.8);
        assert!(first.reaches_limit());
        assert_eq!(first.governing_sample().unwrap().end, 1);
        assert_eq!(runtime.report().motion.displacement, vec![0.0; 6]);
        runtime.step(&[0.0; 6], cx).unwrap();
        let last = runtime.assess_tensile_strength(&capacities, 2, cx).unwrap();
        let decay = (-duration / 0.5).exp();
        for (before, after) in first.samples().iter().zip(last.samples()) {
            near(
                after.assessment.tensile_utilization,
                before.assessment.tensile_utilization * decay,
            );
        }
        assert!(!last.reaches_limit());
        assert_eq!(runtime.report().motion.displacement, vec![0.0; 6]);
        assert_eq!(last.time_s, duration);
    });
}

#[test]
fn g3_timoshenko_memory_keeps_original_curvature_and_distinct_shear_channels() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let (frame, longitudinal, transverse) = super::transverse::shear_frame(true);
        let capacities = [capacity(&frame)];
        let prepared = frame
            .prepare_relaxing_resonance_with_transverse_shear(
                NO_RAYLEIGH,
                &longitudinal,
                &transverse,
                None,
                controls(),
                cx,
            )
            .unwrap();
        let mut runtime = prepared
            .prepare_transient(0.1, transient_controls(), cx)
            .unwrap();
        let length = super::transverse::LENGTH;
        let correction = super::transverse::CORRECTION;
        let phi =
            12.0 * YOUNG * (RADIUS * RADIUS / 4.0) / (correction * (YOUNG / 3.0) * length * length);
        let beta = 1.0 / (1.0 + phi);
        let mut reference = None;
        for shear_history in [0.0, 0.7] {
            runtime
                .initialize(
                    0.0,
                    &[0.0; 6],
                    &[0.0; 6],
                    FrameRelaxationInitialHistory::ElasticExtensions(&[
                        0.0,
                        0.2,
                        -0.3,
                        shear_history,
                        -shear_history,
                    ]),
                    &[],
                    cx,
                )
                .unwrap();
            let result = runtime.assess_tensile_strength(&capacities, 2, cx).unwrap();
            for sample in result.samples() {
                let sign = if sample.end == 0 { 1.0 } else { -1.0 };
                let gradient = 0.9 * (sign * 3.0 * beta * 0.2 - 0.3) / length.powi(2);
                near(sample.normal_stress_gradient_pa_per_m, gradient);
                near(
                    sample.assessment.maximum_axial_stress_pa,
                    RADIUS * gradient.abs(),
                );
            }
            if let Some(previous) = &reference {
                // Normal-only observation must not reinterpret G memory as E.
                assert_eq!(&result, previous);
            }
            reference = Some(result);
        }
    });
}

#[test]
fn g4_original_sources_budget_and_cancellation_refuse_without_changing_history() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let (frame, selections) = frame(&[&[(10.0, 0.5)]], BROADBAND, true, (1.0, 0.0));
        let capacities = [capacity(&frame)];
        let prepared = frame
            .prepare_relaxing_resonance(NO_RAYLEIGH, &selections, None, controls(), cx)
            .unwrap();
        let mut runtime = prepared
            .prepare_transient(0.1, transient_controls(), cx)
            .unwrap();
        let before = runtime.canonical_state().to_vec();
        assert!(matches!(
            runtime.assess_tensile_strength(&capacities, 1, cx),
            Err(FrameTensileError::Budget)
        ));
        assert!(runtime.assess_tensile_strength(&[], 2, cx).is_err());
        let (foreign, _) = super::frame(&[&[(11.0, 0.5)]], BROADBAND, true, (1.0, 0.0));
        assert!(matches!(
            runtime.assess_tensile_strength(&[capacity(&foreign)], 2, cx),
            Err(FrameTensileError::SourceMismatch { element: 0 })
        ));
        // A tiny but nonzero canonical bending history can underflow its native
        // force diagnostic. Observation must refuse instead of claiming safety.
        let memory_state = prepared.coordinates()[1].native_state;
        let mut tiny = before.clone();
        tiny[memory_state] = f64::from_bits(1);
        runtime.restore_state(0.0, &tiny, cx).unwrap();
        assert_ne!(runtime.canonical_state()[memory_state], 0.0);
        assert_eq!(runtime.report().memory[1].opposing_force_n, 0.0);
        assert!(runtime.assess_tensile_strength(&capacities, 2, cx).is_err());
        assert_eq!(runtime.canonical_state(), tiny);
        runtime.restore_state(0.0, &before, cx).unwrap();
        gate.request();
        assert!(matches!(
            runtime.assess_tensile_strength(&capacities, 2, cx),
            Err(FrameTensileError::Cancelled)
        ));
        assert_eq!(runtime.canonical_state(), before);
    });
}
