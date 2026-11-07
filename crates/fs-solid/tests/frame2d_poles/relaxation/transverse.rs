//! Independent complex rigidity and constrained-memory decay on the original
//! Timoshenko interpolation. No complex-frequency shape refitting is used.

use super::*;
use fs_phs::linear_exact::DenseExponentialBudget;
use fs_solid::frame2d::coupled::transient::FrameTransientControls;
use fs_solid::frame2d::relaxation::FrameTransverseRelaxation;

pub(super) const LENGTH: f64 = 0.7;
pub(super) const CORRECTION: f64 = 0.8;

pub(super) fn shear_frame(
    guided_axial: bool,
) -> (
    Frame2d,
    Vec<FrameLongitudinalRelaxation>,
    Vec<FrameTransverseRelaxation>,
) {
    let (layer, longitudinal) = material_layer_with_transverse(
        RADIUS,
        &[(0.9, 0.4)],
        Some(&[(0.0, 0.07), (1.3, 0.23), (0.6, 0.91)]),
        (0.0, f64::MAX),
    );
    let section = ConcentricElasticSection::try_new(0.0, vec![layer])
        .unwrap()
        .with_shear_modulus("shear_modulus", CORRECTION)
        .unwrap();
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0);
    frame.add_node(LENGTH, 0.0);
    frame.fix(0, [true; 3]);
    frame.fix(
        1,
        if guided_axial {
            [false, true, true]
        } else {
            [true, false, true]
        },
    );
    let element = frame.add_material_shear_element(0, 1, 5, section).unwrap();
    (
        frame,
        vec![FrameLongitudinalRelaxation {
            element,
            layers: vec![longitudinal],
        }],
        vec![FrameTransverseRelaxation {
            element,
            layers: vec![LayerRelaxationProperties {
                equilibrium_modulus: "shear_equilibrium".into(),
                branches: (0..3)
                    .map(|i| (format!("shear_increment_{i}"), format!("shear_time_{i}")))
                    .collect(),
            }],
        }],
    )
}

#[test]
fn independent_shear_spectra_change_native_response_and_close_physical_power() {
    with_cx(&CancelGate::new(), |cx| {
        let (frame, longitudinal, transverse) = shear_frame(false);
        assert!(
            frame
                .prepare_relaxing_resonance(NO_RAYLEIGH, &longitudinal, None, controls(), cx)
                .is_err()
        );
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
        assert_eq!(prepared.coordinates().len(), 5);
        let shear = prepared.members()[0].transverse_source.as_ref().unwrap();
        assert_eq!(shear.branches().len(), 3);
        assert_eq!(shear.branches()[0].relaxing_rigidity_n, 0.0);
        let ei = YOUNG * area() * RADIUS.powi(2) / 4.0;
        let ga = CORRECTION * (YOUNG / 3.0) * area();
        let phi = 12.0 * ei / (ga * LENGTH.powi(2));
        let beta = 1.0 / (1.0 + phi);
        let shear_fraction = phi / (1.0 + phi);
        let equilibrium = 12.0 * ei * beta / LENGTH.powi(3);
        let native = prepared.coupled();
        // Exact integrals of original tip translation h=beta*x^2*(3-2*x)+s*x
        // and section spin theta=6*beta*x*(1-x)/L, including rotary inertia.
        let mass = RHO
            * area()
            * LENGTH
            * (shear_fraction.powi(2) / 3.0
                + 0.7 * shear_fraction * beta
                + 13.0 * beta.powi(2) / 35.0)
            + RHO * area() * RADIUS.powi(2) / 4.0 / LENGTH * 1.2 * beta.powi(2);
        for omega in [0.1, 1.0 / 0.91, 1.0 / 0.23, 20.0] {
            let s = C64::new(0.0, -omega);
            let e_arm = 12.0 * (0.9 / YOUNG) * ei * beta.powi(2) / LENGTH.powi(3);
            let mut dynamic = C64::from_re(equilibrium)
                + (s * s).scale(mass)
                + s.scale(e_arm) / (s + C64::from_re(1.0 / 0.4));
            for (delta, tau) in [(1.3, 0.23), (0.6, 0.91)] {
                let stiffness = CORRECTION * delta * area() * shear_fraction.powi(2) / LENGTH;
                dynamic = dynamic + s.scale(stiffness) / (s + C64::from_re(1.0 / tau));
            }
            let response = native
                .factor(omega, cx)
                .unwrap()
                .response(&[C64::ONE], cx)
                .unwrap();
            close(response.state[0], C64::ONE / dynamic);
            let memory = prepared.harmonic_memory(&response, cx).unwrap();
            assert_eq!(
                memory
                    .coordinates
                    .iter()
                    .filter(|c| c.coordinate.channel == FrameStrainChannel::TransverseShear)
                    .count(),
                2
            );
            let supplied = 0.5 * (s * response.state[0]).re;
            assert!(
                (supplied - memory.mean_dissipated_power_w).abs() < 2e-10 * supplied.abs().max(1.0)
            );
        }
        let mut too_small = controls();
        too_small.maximum_states = 7; // Actual 2+3+2 states fit; declared zero G arm still costs one.
        assert!(matches!(
            frame.prepare_relaxing_resonance_with_transverse_shear(
                NO_RAYLEIGH,
                &longitudinal,
                &transverse,
                None,
                too_small,
                cx,
            ),
            Err(FrameRelaxationError::Native(
                FrameCoupledError::Budget { .. }
            ))
        ));
        let mut duplicate = transverse.clone();
        duplicate.push(transverse[0].clone());
        assert!(
            frame
                .prepare_relaxing_resonance_with_transverse_shear(
                    NO_RAYLEIGH,
                    &longitudinal,
                    &duplicate,
                    None,
                    controls(),
                    cx,
                )
                .is_err()
        );
        let mut elastic_shear = transverse.clone();
        elastic_shear[0].layers[0].branches.clear();
        let e_only = frame
            .prepare_relaxing_resonance_with_transverse_shear(
                NO_RAYLEIGH,
                &longitudinal,
                &elastic_shear,
                None,
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(e_only.coordinates().len(), 3);
        assert!(
            e_only.members()[0]
                .transverse_source
                .as_ref()
                .unwrap()
                .branches()
                .is_empty()
        );
        let mut elastic_longitudinal = longitudinal.clone();
        elastic_longitudinal[0].layers[0].branches.clear();
        let g_only = frame
            .prepare_relaxing_resonance_with_transverse_shear(
                NO_RAYLEIGH,
                &elastic_longitudinal,
                &transverse,
                None,
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(g_only.coordinates().len(), 2);
        assert!(
            g_only
                .coordinates()
                .iter()
                .all(|c| c.channel == FrameStrainChannel::TransverseShear)
        );
    });
}

#[test]
fn constrained_shear_memories_decay_at_original_times_and_release_only_their_storage() {
    with_cx(&CancelGate::new(), |cx| {
        let (frame, longitudinal, transverse) = shear_frame(true);
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
        let duration = 0.31;
        let mut runtime = prepared
            .prepare_transient(
                duration,
                FrameTransientControls {
                    exact: DenseExponentialBudget {
                        max_dimension: 128,
                        max_dense_bytes: 64 << 20,
                        max_multiply_terms: 1 << 42,
                    },
                    maximum_relative_balance_error: 1e-10,
                },
                cx,
            )
            .unwrap();
        let mut initial = vec![0.0; prepared.coupled().dimension()];
        let mut expected_heat = 0.0;
        for coordinate in prepared
            .coordinates()
            .iter()
            .filter(|c| c.channel == FrameStrainChannel::TransverseShear)
        {
            let z = 0.02 * (coordinate.source_arm + 1) as f64;
            initial[coordinate.native_state] = z;
            expected_heat +=
                0.5 * z * z * (1.0 - (-2.0 * duration / coordinate.relaxation_time_s).exp());
        }
        runtime.restore_state(0.0, &initial, cx).unwrap();
        runtime.step(&[0.0; 6], cx).unwrap();
        for coordinate in prepared.coordinates() {
            let expected =
                initial[coordinate.native_state] * (-duration / coordinate.relaxation_time_s).exp();
            assert!((runtime.canonical_state()[coordinate.native_state] - expected).abs() < 2e-12);
        }
        assert_eq!(runtime.canonical_state()[0], 0.0);
        assert!((runtime.report().material_dissipation_j - expected_heat).abs() < 2e-12);
    });
}
