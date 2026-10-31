//! Independent continuum response and original contact consumer of fiber memory.
use super::*;
use fs_couple::fiber_ribbon::FiberRibbonError;
use fs_material::section::relaxation::{ConcentricPronySection, LayerRelaxationProperties};
use fs_math::c64::C64;
use fs_qty::Time;

fn source(
    rho: f64,
    young: f64,
    delta: f64,
    tau: f64,
    band: (f64, f64),
) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (key, dims, value) in [
        ("density", Density::DIMS, rho),
        ("young_modulus", Pressure::DIMS, young),
        (ETA, DynViscosity::DIMS, 0.0),
        (YOUNG_VISCOSITY_PROPERTY, DynViscosity::DIMS, 0.0),
        ("delta", Pressure::DIMS, delta),
        ("tau", Time::DIMS, tau),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 310.0)
                    .with("omega", band.0, band.1),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured layered fiber spectrum, not measured horsehair".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-fiber-spectrum".into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with("T", 293.15)
            .unwrap()
            .with("omega", 1000.0)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn specimen(
    cells: usize,
    delta: f64,
    tau: f64,
    band: (f64, f64),
    gate: &CancelGate,
) -> (TensionedFiberRibbon, ConcentricPronySection) {
    let section = ConcentricElasticSection::try_new(
        0.0,
        vec![
            ConcentricMaterialLayer {
                outer_radius_m: 20e-6,
                material: source(1300.0, 3e9, delta, tau, band),
            },
            ConcentricMaterialLayer {
                outer_radius_m: 40e-6,
                material: source(2000.0, 5e9, 0.5 * delta, 2.0 * tau, band),
            },
        ],
    )
    .unwrap();
    let ribbon = TensionedFiberRibbon::try_new(
        FiberRibbonGeometry {
            length_m: 0.6,
            width_m: 0.01,
            fiber_count: 180,
            tension_n: 40.0,
        },
        FiberRibbonResolution {
            cells,
            axial_modes: 2,
            transverse_modes: 2,
        },
        section
            .clone()
            .with_axial_viscosity(YOUNG_VISCOSITY_PROPERTY)
            .unwrap(),
        section.clone().with_kelvin_voigt_bending(ETA).unwrap(),
        48000,
        ModalAcousticTimeBudget::audible_reference(),
        gate,
    )
    .unwrap();
    let spectrum = section
        .with_prony_bending(&vec![
            LayerRelaxationProperties {
                equilibrium_modulus: "young_modulus".into(),
                branches: vec![("delta".into(), "tau".into())],
            };
            2
        ])
        .unwrap();
    (ribbon, spectrum)
}

fn relaxing(cells: usize, delta: f64, gate: &CancelGate) -> TensionedFiberRibbon {
    let (ribbon, spectrum) = specimen(cells, delta, 0.0003, (0.0, 1e8), gate);
    ribbon
        .with_material_relaxation(spectrum, 256, gate)
        .unwrap()
}

#[test]
fn hereditary_fiber_response_matches_complex_continuum_and_width_partition() {
    with_cx(|cx, gate| {
        for cells in [1, 3, 5] {
            let fiber = relaxing(cells, 1e9, gate);
            let source = fiber.material_relaxation().unwrap().source().clone();
            let expected_states = cells * 4 * 2;
            assert_eq!(
                fiber.material_relaxation().unwrap().memory_states(),
                expected_states
            );
            let bank = base(cx, gate, 7800.0)
                .with_tensioned_fiber_ribbon(fiber, 0.4, ribbon_limits(), gate)
                .unwrap();
            let system = bank.system();
            let first = bank.source().fiber_ribbon().unwrap().component_start();
            let start: usize = system.components()[..first]
                .iter()
                .map(|m| m.modes().len())
                .sum();
            assert_eq!(system.material_memory().len(), expected_states);
            let mu = source.elastic().line_density_kg_m();
            // Complex continuum modulus is evaluated before any Galerkin projection.
            // The two original layers keep their distinct relaxation times.
            for omega in [10.0, 300.0, 2500.0, 9000.0] {
                let mut ea = C64::ZERO;
                let mut ei = C64::ZERO;
                for (law, shape) in source
                    .spectra()
                    .iter()
                    .zip(source.elastic().layer_integrals())
                {
                    let (storage, loss) = law.modulus(omega);
                    let young = C64::new(storage, -loss);
                    ea = ea + young.scale(shape.area_m2());
                    ei = ei + young.scale(shape.second_moment_m4());
                }
                for mode in 0..4 {
                    let k = ((mode % 2) + 1) as f64 * core::f64::consts::PI / 0.6;
                    let dynamic = if mode < 2 {
                        ea.scale(k * k / mu)
                    } else {
                        ei.scale(k.powi(4) / mu) + C64::new((40.0 / 180.0) * k * k / mu, 0.0)
                    } - C64::new(omega * omega, 0.0);
                    let expected = C64::new(1e-4, 0.0) / dynamic;
                    let mut force = vec![C64::ZERO; system.mode_count()];
                    force[start + mode] = C64::new(1e-4, 0.0);
                    let response = system
                        .released_harmonic_response(omega, &force, 1_048_576, gate)
                        .unwrap();
                    let actual = response.dynamics.displacement_m_sqrt_kg[start + mode];
                    assert!(
                        (actual - expected).abs() < 2e-10 * expected.abs(),
                        "cells={cells}, mode={mode}, omega={omega}: {actual:?} versus {expected:?}"
                    );
                    let power = -0.5 * dynamic.im * omega * expected.abs().powi(2);
                    assert!(
                        (response.dynamics.material_dissipation_w - power).abs() < 2e-10 * power
                    );
                    assert!(response.dynamics.power_residual_w.abs() < 2e-10 * power);
                }
            }
            assert_eq!(system.samples_rendered(), 0);
            assert!(system.material_memory().iter().all(|z| *z == 0.0));
        }
    });
}

#[test]
fn material_fiber_history_changes_contact_and_survives_support_motion_and_callbacks() {
    with_cx(|cx, gate| {
        let make = |delta| {
            base(cx, gate, 7800.0)
                .with_tensioned_fiber_ribbon(relaxing(1, delta, gate), 0.4, ribbon_limits(), gate)
                .unwrap()
                .into_scheduled_ribbon(
                    fiber_load_schedule(0.03),
                    "bow".into(),
                    performance_config(),
                )
                .unwrap()
                .with_fiber_support_program(
                    support_program(&[(128, 0.02), (256, -0.02), (384, 0.0)]),
                    gate,
                )
                .unwrap()
        };
        let (mut a, mut b, mut elastic) = (make(1e9), make(1e9), make(0.0));
        let (mut pa, mut pb, mut pe) = (vec![0.0; 513], vec![0.0; 513], vec![0.0; 513]);
        let mut memory_loss = 0.0;
        let mut traction_difference = 0.0;
        for (p, e) in pa.iter_mut().zip(&mut pe) {
            a.block(core::slice::from_mut(p)).unwrap();
            elastic.block(core::slice::from_mut(e)).unwrap();
            let frame = a.bank().system().last_frame().unwrap();
            let reference = elastic.bank().system().last_frame().unwrap();
            traction_difference +=
                (frame.ports[0].traction_n - reference.ports[0].traction_n).abs();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            memory_loss += a
                .bank()
                .system()
                .coupled_frame()
                .unwrap()
                .material_dissipation_j;
        }
        for chunk in pb.chunks_mut(37) {
            b.block(chunk).unwrap();
        }
        assert_eq!(pa, pb);
        assert_eq!(
            a.bank().system().material_memory(),
            b.bank().system().material_memory()
        );
        assert_eq!(states(a.bank()), states(b.bank()));
        assert!(
            a.bank()
                .system()
                .material_memory()
                .iter()
                .any(|z| z.abs() > 1e-16)
        );
        assert!(memory_loss > 0.0 && traction_difference > 1e-8);
        assert!(pa.iter().zip(&pe).any(|(a, b)| (a - b).abs() > 1e-12));
        assert_eq!(
            a.bank()
                .source()
                .fiber_ribbon()
                .unwrap()
                .model()
                .material_relaxation()
                .unwrap()
                .source()
                .spectra()[0]
                .terms,
            vec![(1e9, 0.0003)]
        );
        let mut direct = base(cx, gate, 7800.0)
            .with_tensioned_fiber_ribbon(relaxing(1, 1e9, gate), 0.4, ribbon_limits(), gate)
            .unwrap();
        let zero = vec![0.0; direct.system().mode_count()];
        for _ in 0..16 {
            direct.step(&zero, &[drive(0.01, 0.03)]).unwrap();
        }
        let before_memory = direct.system().material_memory().to_vec();
        let before_states = states(&direct);
        gate.request();
        assert!(
            direct
                .step_under_gate(&zero, &[drive(0.01, 0.03)], gate)
                .is_err()
        );
        assert_eq!(direct.system().material_memory(), before_memory);
        assert_eq!(states(&direct), before_states);
    });
}

#[test]
fn fast_fiber_axial_and_bending_history_reaches_the_original_contact_owner() {
    with_cx(|cx, gate| {
        let tau = 1e-7;
        let (ribbon, spectrum) = specimen(1, 1e9, tau, (0.0, 1e8), gate);
        let fiber = ribbon.with_material_relaxation(spectrum, 256, gate).unwrap();
        let mut bank = base(cx, gate, 7800.0)
            .with_tensioned_fiber_ribbon(fiber, 0.4, ribbon_limits(), gate).unwrap();
        assert!(bank.system().sample_period_s() / tau > 200.0);
        assert!(bank.system().material_relaxation_branches().iter()
            .any(|branch| branch.relaxation_time_s == tau));
        let zero = vec![0.0; bank.system().mode_count()];
        let mut heat = 0.0;
        for _ in 0..32 {
            bank.step(&zero, &[drive(0.01, 0.03)]).unwrap();
            let frame = bank.system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            heat += bank.system().coupled_frame().unwrap().material_dissipation_j;
        }
        assert!(heat > 0.0);
        assert!(bank.system().material_memory().iter().any(|z| *z != 0.0));
    });
}

#[test]
fn full_strain_fiber_memory_composes_in_both_orders_and_retains_original_coordinates() {
    with_cx(|cx, gate| {
        // Both exact bending and exact nonlinear strain histories must remain
        // composable when each physical sample spans hundreds of decay times.
        let (ribbon, spectrum) = specimen(3, 1e9, 1e-7, (0.0, 1e8), gate);
        let strain_first = ribbon
            .clone()
            .with_axial_transverse_strain(strain_limits(), strain_config())
            .unwrap()
            .with_material_relaxation(spectrum.clone(), 30, gate)
            .unwrap();
        let memory_first = ribbon
            .with_material_relaxation(spectrum, 30, gate)
            .unwrap()
            .with_axial_transverse_strain(strain_limits(), strain_config())
            .unwrap();
        for source in [&strain_first, &memory_first] {
            let memory = source.material_relaxation().unwrap();
            assert!(memory.has_finite_strain());
            assert_eq!(memory.memory_states(), 30);
            assert_eq!(memory.channels().len(), 4);
            assert!(memory.channels().iter().all(|channel| channel.mode >= 2));
        }
        let mut a = base(cx, gate, 7800.0)
            .with_tensioned_fiber_ribbon(strain_first, 0.4, ribbon_limits(), gate)
            .unwrap();
        let mut b = base(cx, gate, 7800.0)
            .with_tensioned_fiber_ribbon(memory_first, 0.4, ribbon_limits(), gate)
            .unwrap();
        assert_eq!(a.system().material_memory().len(), 12);
        assert_eq!(
            a.system()
                .nonlinear_material_memory()
                .iter()
                .map(|m| m.len())
                .sum::<usize>(),
            18
        );
        let component = a.source().fiber_ribbon().unwrap().component_start();
        let first: usize = a.system().components()[..component]
            .iter()
            .map(|c| c.modes().len())
            .sum();
        let mut force = vec![0.0; a.system().mode_count()];
        // A physical force on the retained NORMAL fiber mode must also drive
        // its original axial material memory through geometric stretching.
        force[first + 2] = 5e-6;
        let drives = vec![drive(0.01, 0.03); 3];
        for _ in 0..64 {
            a.step(&force, &drives).unwrap();
            b.step(&force, &drives).unwrap();
            let frame = a.system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        assert_eq!(states(&a), states(&b));
        assert_eq!(a.system().material_memory(), b.system().material_memory());
        assert_eq!(
            a.system().nonlinear_material_memory(),
            b.system().nonlinear_material_memory()
        );
        assert!(
            a.system()
                .nonlinear_material_memory()
                .iter()
                .any(|memory| memory[0] > 0.0)
        );
        let before: Vec<_> = a
            .system()
            .nonlinear_material_memory()
            .iter()
            .map(|m| m.to_vec())
            .collect();
        let state = states(&a);
        let bending = a.system().material_memory().to_vec();
        gate.request();
        assert!(a.step_under_gate(&force, &drives, gate).is_err());
        assert_eq!(states(&a), state);
        assert_eq!(a.system().material_memory(), bending);
        assert_eq!(
            a.system()
                .nonlinear_material_memory()
                .iter()
                .map(|m| m.to_vec())
                .collect::<Vec<_>>(),
            before
        );
    });
}

#[test]
fn material_fiber_memory_refuses_incompatible_sources_clocks_and_incomplete_laws() {
    let gate = CancelGate::new_clock_free();
    let (ribbon, spectrum) = specimen(3, 1e9, 0.0003, (0.0, 1e8), &gate);
    assert!(
        ribbon
            .clone()
            .with_material_relaxation(spectrum.clone(), 23, &gate)
            .is_err()
    );
    let admitted = ribbon
        .clone()
        .with_material_relaxation(spectrum.clone(), 24, &gate)
        .unwrap();
    assert!(
        admitted
            .clone()
            .with_material_relaxation(spectrum.clone(), 24, &gate)
            .is_err()
    );
    // The linear 24-state budget omits the six mean-strain memories; neither
    // builder order may silently drop them when nonlinear strain is selected.
    assert!(
        admitted
            .with_axial_transverse_strain(strain_limits(), strain_config())
            .is_err()
    );
    assert!(
        ribbon
            .with_axial_transverse_strain(strain_limits(), strain_config())
            .unwrap()
            .with_material_relaxation(spectrum, 24, &gate)
            .is_err()
    );
    let (ribbon, _) = specimen(1, 1e9, 0.0003, (0.0, 1e8), &gate);
    let (_, different) = specimen(1, 2e9, 0.0003, (0.0, 1e8), &gate);
    assert!(
        ribbon
            .with_material_relaxation(different, 256, &gate)
            .is_err()
    );
    let (ribbon, spectrum) = specimen(1, 1e9, 0.0003, (1.0, 1e8), &gate);
    assert!(
        ribbon
            .with_material_relaxation(spectrum, 256, &gate)
            .is_err()
    );
    let (ribbon, spectrum) = specimen(1, 1e9, 0.0003, (0.0, 1e8), &gate);
    gate.request();
    assert!(matches!(
        ribbon.with_material_relaxation(spectrum, 256, &gate),
        Err(FiberRibbonError::Cancelled)
    ));
}
