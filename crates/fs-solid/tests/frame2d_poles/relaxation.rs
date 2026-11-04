//! Source Maxwell storage/loss, original native poles, full strain coupling and
//! physical descriptor loading. Independent closed-form member equations.

use super::{NO_RAYLEIGH, RADIUS, RHO, YOUNG, guided, with_cx};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::relaxation::LayerRelaxationProperties;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_math::c64::C64;
use fs_phs::harmonic::PoleOptions;
use fs_phs::{DescriptorPortHamiltonian, QuadraticStorage};
use fs_qty::{Density, Pressure, Time};
use fs_solid::frame2d::{
    FrameMassFormulation,
    coupled::{FrameCoupledControls, FrameCoupledError, FrameDynamicPortLoad, FramePortDynamics},
    relaxation::{FrameLongitudinalRelaxation, FrameRelaxationError, FrameStrainChannel},
};
use fs_solid::{Frame2d, Geometry, Rayleigh};

fn controls() -> FrameCoupledControls {
    FrameCoupledControls {
        maximum_coordinates: 12,
        maximum_ports: 8,
        maximum_states: 32,
        maximum_payload_bytes: 8 << 20,
        maximum_work: 1 << 36,
        maximum_relative_residual: 1e-10,
    }
}
fn options() -> PoleOptions {
    PoleOptions {
        spectral_scale: 16.0,
        step_tolerance: 1e-13,
        residual_tolerance: 1e-13,
        max_iterations: 40,
        max_dimension: 32,
        max_dense_bytes: 8 << 20,
        max_work: 100_000_000,
    }
}
fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() < 2e-10 * (1.0 + expected.abs()),
        "actual {actual:?}, expected {expected:?}"
    );
}
fn area() -> f64 {
    core::f64::consts::PI * RADIUS * RADIUS
}

fn material_layer(
    radius: f64,
    arms: &[(f64, f64)],
    band: (f64, f64),
) -> (ConcentricMaterialLayer, LayerRelaxationProperties) {
    let mut definitions = vec![
        ("density".to_string(), RHO, Density::DIMS, false),
        ("young_modulus".to_string(), YOUNG, Pressure::DIMS, false),
        ("equilibrium".to_string(), YOUNG, Pressure::DIMS, true),
    ];
    let mut branches = Vec::new();
    for (index, &(delta, tau)) in arms.iter().enumerate() {
        let delta_name = format!("increment_{index}");
        let tau_name = format!("time_{index}");
        definitions.push((delta_name.clone(), delta, Pressure::DIMS, true));
        definitions.push((tau_name.clone(), tau, Time::DIMS, true));
        branches.push((delta_name, tau_name));
    }
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims, frequency_bound) in definitions {
        let mut validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
        if frequency_bound {
            validity = validity.with("omega", band.0, band.1);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name.clone(), dims),
                value: PropertyValue::Scalar { value, dims },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured longitudinal Maxwell frame fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured solid".into(),
            phase: "solid".into(),
            process: "source Maxwell fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with("omega", 0.5 * (band.0 + band.1))
        .unwrap();
    let material = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (
        ConcentricMaterialLayer {
            outer_radius_m: radius,
            material,
        },
        LayerRelaxationProperties {
            equilibrium_modulus: "equilibrium".into(),
            branches,
        },
    )
}

fn frame(
    spectra: &[&[(f64, f64)]],
    band: (f64, f64),
    guided_tip: bool,
    direction: (f64, f64),
) -> (Frame2d, Vec<FrameLongitudinalRelaxation>) {
    let mut layers = Vec::new();
    let mut properties = Vec::new();
    for (index, arms) in spectra.iter().enumerate() {
        // Equal-area original layers; distinct material spectra remain distinct.
        let radius = RADIUS * ((index + 1) as f64 / spectra.len() as f64).sqrt();
        let (layer, binding) = material_layer(radius, arms, band);
        layers.push(layer);
        properties.push(binding);
    }
    let source = ConcentricElasticSection::try_new(0.0, layers).unwrap();
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0);
    frame.add_node(direction.0, direction.1);
    frame.fix(0, [true; 3]);
    if guided_tip {
        frame.fix(1, [false, true, true]);
    }
    let element = frame.add_material_element(0, 1, 5, source).unwrap();
    (
        frame,
        vec![FrameLongitudinalRelaxation {
            element,
            layers: properties,
        }],
    )
}

#[test]
fn authored_maxwell_cubic_moves_all_native_rates_and_closes_material_heat() {
    with_cx(&CancelGate::new(), |cx| {
        let (frame, selections) = frame(&[&[(0.5, 2.0 / 3.0)]], (0.0, 1e8), true, (1.0, 0.0));
        let accepted = format!("{frame:?}");
        let prepared = frame
            .prepare_relaxing_resonance(
                Rayleigh { a0: 4.5, a1: 0.0 },
                &selections,
                None,
                controls(),
                cx,
            )
            .unwrap();
        let model = prepared.coupled();
        assert_eq!(model.free_dofs(), [3]);
        assert_eq!(prepared.memory_state_range(), 2..5);
        assert_eq!(prepared.boundary_state_range(), 5..5);
        assert_eq!(model.material_omega_band_rad_s(), Some((0.0, 1e8)));
        // M=2A, K_infinity=8A, C=9A, delta_K=.5A, tau=2/3.
        // Dividing the original characteristic polynomial by M gives
        // (s+1.5)*(s^2+4.5s+4)+.25s=(s+1)(s+2)(s+3).
        let mass = 2.0 * area();
        let stiffness = 8.0 * area();
        let loss = 9.0 * area();
        let delta = 0.5 * area();
        for rate in [-1.0, -2.0, -3.0] {
            let seed = [
                C64::ONE,
                C64::from_re(mass * rate),
                C64::from_re(delta.sqrt() * rate / (rate + 1.5)),
                C64::ZERO,
                C64::ZERO,
            ];
            let pole = model
                .resonance_pole(C64::from_re(rate + 0.02), &seed, &options(), cx)
                .unwrap();
            close(pole.laplace_pole_per_s, C64::from_re(rate));
            assert!(pole.relative_residual < 1e-10);
            assert!(pole.decay_power_residual_w.abs() < 1e-9);
        }
        for omega in [0.25, 1.5, 3.0, 20.0] {
            let s = C64::new(0.0, -omega);
            let maxwell = s.scale(delta) / (s + C64::from_re(1.5));
            let expected = C64::ONE
                / ((s * s).scale(mass) + s.scale(loss) + C64::from_re(stiffness) + maxwell);
            let response = model
                .factor(omega, cx)
                .unwrap()
                .response(&[C64::ONE], cx)
                .unwrap();
            close(response.state[0], expected);
            let memory = prepared.harmonic_memory(&response, cx).unwrap();
            assert_eq!(memory.coordinates.len(), 3);
            assert_eq!(
                memory.coordinates[0].coordinate.channel,
                FrameStrainChannel::Axial
            );
            close(memory.coordinates[0].opposing_force_n, maxwell * expected);
            close(memory.coordinates[0].extension_m, expected);
            for dark in &memory.coordinates[1..] {
                assert_eq!(dark.extension_m, C64::ZERO);
                assert_eq!(dark.mean_stored_energy_j, 0.0);
            }
            let velocity = s * expected;
            let energy = 0.25 * (stiffness * expected.norm_sq() + mass * velocity.norm_sq())
                + memory.mean_stored_energy_j;
            let heat = 0.5 * loss * velocity.norm_sq() + memory.mean_dissipated_power_w;
            assert!((response.mean_stored_energy_j - energy).abs() < 1e-11);
            assert!((response.dissipated_power_w - heat).abs() < 1e-11);
            assert!(response.power_residual_w.abs() < 1e-11);
        }
        // Dark constrained strain memories are real energy states. They remain
        // independently observable rather than being removed by a zero projection.
        let mut observation = model
            .factor(1.0, cx)
            .unwrap()
            .response(&[C64::ZERO], cx)
            .unwrap();
        observation.state[3] = C64::new(2.0, -1.0);
        let memory = prepared.harmonic_memory(&observation, cx).unwrap();
        assert_eq!(memory.coordinates[1].extension_m, C64::ZERO);
        assert_eq!(memory.coordinates[1].mean_stored_energy_j, 1.25);
        assert!((memory.coordinates[1].mean_dissipated_power_w - 3.75).abs() < 1e-14);
        assert_eq!(format!("{frame:?}"), accepted);
    });
}

#[test]
fn original_layers_keep_different_relaxation_times_and_zero_source_identity() {
    with_cx(&CancelGate::new(), |cx| {
        let (frame, selections) = frame(
            &[&[(0.0, 0.1), (0.5, 2.0 / 3.0)], &[(1.5, 0.2)]],
            (0.0, 100.0),
            true,
            (1.0, 0.0),
        );
        let prepared = frame
            .prepare_relaxing_resonance(NO_RAYLEIGH, &selections, None, controls(), cx)
            .unwrap();
        assert_eq!(prepared.coordinates().len(), 6);
        assert_eq!(prepared.members()[0].source.spectra()[0].terms.len(), 2);
        assert_eq!(prepared.coordinates()[0].source_arm, 1);
        assert_eq!(prepared.coordinates()[3].layer, 1);
        assert_eq!(prepared.coordinates()[0].relaxation_time_s, 2.0 / 3.0);
        assert_eq!(prepared.coordinates()[3].relaxation_time_s, 0.2);
        for omega in [0.2, 1.5, 5.0, 40.0] {
            let s = C64::new(0.0, -omega);
            let dynamic = s.scale(0.25 * area()) / (s + C64::from_re(1.5))
                + s.scale(0.75 * area()) / (s + C64::from_re(5.0));
            let expected =
                C64::ONE / ((s * s).scale(2.0 * area()) + C64::from_re(8.0 * area()) + dynamic);
            let response = prepared
                .coupled()
                .factor(omega, cx)
                .unwrap()
                .response(&[C64::ONE], cx)
                .unwrap();
            close(response.state[0], expected);
            let material = prepared.harmonic_memory(&response, cx).unwrap();
            close(
                material.coordinates[0].opposing_force_n + material.coordinates[3].opposing_force_n,
                dynamic * expected,
            );
            assert!((response.dissipated_power_w - material.mean_dissipated_power_w).abs() < 1e-11);
        }
    });
}

#[test]
fn both_bending_strains_preserve_full_mass_and_orientation_coupling() {
    with_cx(&CancelGate::new(), |cx| {
        for direction in [(1.0, 0.0), (0.6, 0.8)] {
            let (frame, selections) = frame(&[&[(2.0, 0.3)]], (0.0, 100.0), false, direction);
            let prepared = frame
                .prepare_relaxing_resonance(NO_RAYLEIGH, &selections, None, controls(), cx)
                .unwrap();
            assert_eq!(prepared.coupled().free_dofs(), [3, 4, 5]);
            let second_moment = core::f64::consts::PI * RADIUS.powi(4) / 4.0;
            // Exact consistent transverse AND rotary-inertia submatrices.
            let mass = RHO * area();
            let spin = RHO * second_moment;
            let m11 = mass * 156.0 / 420.0 + spin * 36.0 / 30.0;
            let m12 = -mass * 22.0 / 420.0 - spin * 3.0 / 30.0;
            let m22 = mass * 4.0 / 420.0 + spin * 4.0 / 30.0;
            for omega in [0.7, 1.2, 4.0, 15.0] {
                let s = C64::new(0.0, -omega);
                let ei = (C64::from_re(YOUNG) + s.scale(2.0) / (s + C64::from_re(1.0 / 0.3)))
                    .scale(second_moment);
                let d11 = (s * s).scale(m11) + ei.scale(12.0);
                let d12 = (s * s).scale(m12) - ei.scale(6.0);
                let d22 = (s * s).scale(m22) + ei.scale(4.0);
                let determinant = d11 * d22 - d12 * d12;
                let displacement = d22 / determinant;
                let rotation = -d12 / determinant;
                let force = [
                    C64::from_re(-direction.1),
                    C64::from_re(direction.0),
                    C64::ZERO,
                ];
                let response = prepared
                    .coupled()
                    .factor(omega, cx)
                    .unwrap()
                    .response(&force, cx)
                    .unwrap();
                close(response.state[0], displacement.scale(-direction.1));
                close(response.state[1], displacement.scale(direction.0));
                close(response.state[2], rotation);
                let material = prepared.harmonic_memory(&response, cx).unwrap();
                close(
                    material.coordinates[1].extension_m,
                    displacement.scale(-2.0) + rotation,
                );
                close(material.coordinates[2].extension_m, -rotation);
                assert!(
                    (material.mean_dissipated_power_w - response.dissipated_power_w).abs() < 1e-9
                );
            }
        }
    });
}

#[test]
fn maxwell_and_descriptor_boundary_share_original_constraints_and_work() {
    with_cx(&CancelGate::new(), |cx| {
        let (frame, selections) = frame(&[&[(0.5, 2.0 / 3.0)]], (0.0, 100.0), true, (1.0, 0.0));
        let added_mass = 3.0 * area();
        let drag = area();
        // z_dot=lambda, 0=-z/m+v, 0=-mu/drag+v. Opposing force=lambda+mu.
        let boundary = DescriptorPortHamiltonian::new(
            3,
            1,
            1,
            vec![0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / drag],
            vec![0.0, 1.0, 1.0],
            Box::new(QuadraticStorage::new(vec![1.0 / added_mass], 1).unwrap()),
        )
        .unwrap();
        let map = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let prepared = frame
            .prepare_relaxing_resonance(
                NO_RAYLEIGH,
                &selections,
                Some(FrameDynamicPortLoad {
                    force_map: &map,
                    dynamics: FramePortDynamics::Descriptor(&boundary),
                }),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(prepared.memory_state_range(), 2..5);
        assert_eq!(prepared.boundary_state_range(), 5..8);
        assert_eq!(prepared.coupled().load_algebraic_range(), 6..8);
        for omega in [0.3, 1.5, 3.0, 9.0] {
            let s = C64::new(0.0, -omega);
            let maxwell = s.scale(0.5 * area()) / (s + C64::from_re(1.5));
            let expected = C64::ONE
                / ((s * s).scale(5.0 * area())
                    + s.scale(drag)
                    + C64::from_re(8.0 * area())
                    + maxwell);
            let response = prepared
                .coupled()
                .factor(omega, cx)
                .unwrap()
                .response(&[C64::ONE], cx)
                .unwrap();
            close(response.state[0], expected);
            let velocity = s * expected;
            close(response.state[5], velocity.scale(added_mass));
            close(response.state[6], (s * velocity).scale(added_mass));
            close(response.state[7], velocity.scale(drag));
            let material = prepared.harmonic_memory(&response, cx).unwrap();
            let energy = 0.25
                * (8.0 * area() * expected.norm_sq() + 5.0 * area() * velocity.norm_sq())
                + material.mean_stored_energy_j;
            let heat = 0.5 * drag * velocity.norm_sq() + material.mean_dissipated_power_w;
            assert!((response.mean_stored_energy_j - energy).abs() < 1e-10);
            assert!((response.dissipated_power_w - heat).abs() < 1e-10);
            assert!(response.power_residual_w.abs() < 1e-10);
        }
    });
}

#[test]
fn zero_and_explicitly_empty_spectra_keep_the_original_elastic_frame_only() {
    with_cx(&CancelGate::new(), |cx| {
        for arms in [&[(0.0, 0.5)][..], &[][..]] {
            let (frame, selections) = frame(&[arms], (0.0, 20.0), true, (1.0, 0.0));
            let prepared = frame
                .prepare_relaxing_resonance(NO_RAYLEIGH, &selections, None, controls(), cx)
                .unwrap();
            assert!(prepared.coordinates().is_empty());
            assert_eq!(prepared.memory_state_range(), 2..2);
            assert_eq!(prepared.coupled().dimension(), 2);
            assert_eq!(
                prepared.members()[0].source.spectra()[0].terms.len(),
                arms.len()
            );
            let response = prepared
                .coupled()
                .factor(1.0, cx)
                .unwrap()
                .response(&[C64::ONE], cx)
                .unwrap();
            close(response.state[0], C64::from_re(1.0 / (6.0 * area())));
            assert_eq!(response.dissipated_power_w, 0.0);
            let memory = prepared.harmonic_memory(&response, cx).unwrap();
            assert!(memory.coordinates.is_empty());
            assert_eq!(memory.mean_stored_energy_j, 0.0);
            let rate = C64::new(0.0, 2.0);
            let pole = prepared
                .coupled()
                .resonance_pole(
                    C64::new(0.01, 2.01),
                    &[C64::ONE, rate.scale(2.0 * area())],
                    &options(),
                    cx,
                )
                .unwrap();
            close(pole.laplace_pole_per_s, rate);
        }
    });
}

#[test]
fn source_support_duplicate_loss_resources_and_cancellation_remain_binding() {
    let (frame, selections) = frame(&[&[(0.5, 2.0 / 3.0)]], (0.5, 3.0), true, (1.0, 0.0));
    let accepted = format!("{frame:?}");
    with_cx(&CancelGate::new(), |cx| {
        let prepared = frame
            .prepare_relaxing_resonance(
                Rayleigh { a0: 4.5, a1: 0.0 },
                &selections,
                None,
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            prepared.coupled().material_omega_band_rad_s(),
            Some((0.5, 3.0))
        );
        assert!(prepared.coupled().factor(1.0, cx).is_ok());
        for omega in [0.0, 0.25, 4.0, f64::NAN] {
            assert!(prepared.coupled().factor(omega, cx).is_err());
        }
        // The finite real-frequency band is provenance for complex poles, not
        // an Im(s) filter that would erase this original real decay mode.
        let seed = [
            C64::ONE,
            C64::from_re(-2.0 * area()),
            C64::from_re(-2.0 * (0.5 * area()).sqrt()),
            C64::ZERO,
            C64::ZERO,
        ];
        let pole = prepared
            .coupled()
            .resonance_pole(C64::from_re(-0.98), &seed, &options(), cx)
            .unwrap();
        close(pole.laplace_pole_per_s, C64::from_re(-1.0));
        let mut missing = selections.clone();
        missing[0].layers[0].branches[0].1 = "absent_source_time".into();
        assert!(matches!(
            frame.prepare_relaxing_resonance(NO_RAYLEIGH, &missing, None, controls(), cx),
            Err(FrameRelaxationError::Source(_))
        ));
        let duplicate = [selections[0].clone(), selections[0].clone()];
        assert!(
            frame
                .prepare_relaxing_resonance(NO_RAYLEIGH, &duplicate, None, controls(), cx)
                .is_err()
        );
        // Reject an unspecified second material loss before asking for selectors.
        assert!(matches!(
            guided(0.2).prepare_relaxing_resonance(NO_RAYLEIGH, &selections, None, controls(), cx),
            Err(FrameRelaxationError::Native(FrameCoupledError::Invalid(_)))
        ));
        for cap in [
            FrameCoupledControls {
                maximum_ports: 2,
                ..controls()
            },
            FrameCoupledControls {
                maximum_states: 4,
                ..controls()
            },
            FrameCoupledControls {
                maximum_payload_bytes: 1,
                ..controls()
            },
            FrameCoupledControls {
                maximum_work: 1,
                ..controls()
            },
        ] {
            assert!(matches!(
                frame.prepare_relaxing_resonance(NO_RAYLEIGH, &selections, None, cap, cx),
                Err(FrameRelaxationError::Native(
                    FrameCoupledError::Budget { .. }
                ))
            ));
        }
    });
    let gate = CancelGate::new();
    gate.request();
    with_cx(&gate, |cx| {
        assert!(matches!(
            frame.prepare_relaxing_resonance(NO_RAYLEIGH, &selections, None, controls(), cx),
            Err(FrameRelaxationError::Native(FrameCoupledError::Cancelled))
        ));
    });
    assert_eq!(format!("{frame:?}"), accepted);
}
