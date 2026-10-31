//! Regional constitutive histories stay in the same string/body/gas owner.
//! Synthetic material fixtures test coupling, not measured instrument accuracy.
use super::*;
use fs_couple::air_path::cylinder::CylinderFlowResolution;
use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::LongitudinalStringResolution;
use fs_couple::plate_viscosity::relaxation::{InitialPlateMemory, PlateRelaxationLimits};
use fs_couple::render::schedule::force::coupled::contact::friction::nonlinear::NonlinearModalConfig;

fn properties() -> Vec<(String, String)> {
    vec![("maxwell_modulus".into(), "maxwell_time".into())]
}

fn limits() -> PlateRelaxationLimits {
    PlateRelaxationLimits {
        max_modes: 3,
        max_projection_terms: 1_000_000,
        max_memory_states: 6,
    }
}

fn memory_material(
    rho: f64,
    excess: f64,
    tau: f64,
    band: (f64, f64),
) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (key, dims, value, omega) in [
        ("density", Density::DIMS, rho, (0.0, 100_000.0)),
        ("young_modulus", Pressure::DIMS, 200e9, (0.0, 100_000.0)),
        ("equilibrium_young_modulus", Pressure::DIMS, 200e9, (0.0, 100_000.0)),
        ("poisson_ratio", Dims::NONE, 0.3, (0.0, 100_000.0)),
        ("maxwell_modulus", Pressure::DIMS, excess, band),
        ("maxwell_time", Dims([0, 0, 1, 0, 0, 0]), tau, band),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 310.0)
                    .with("omega", omega.0, omega.1),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                provenance: Provenance {
                    source: "synthetic regional Maxwell body".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
                observations: vec![],
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "regional Maxwell fixture".into(),
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

fn memory_body(right_tau: f64, band: (f64, f64)) -> SharedMaterialPlateChartBody {
    let mut body = body();
    body.reduction.damping_ratio = 0.0;
    body.specimen = body
        .specimen
        .with_region_material("left", &memory_material(7800.0, 20e9, 0.001, band))
        .unwrap()
        .with_region_material("right", &memory_material(3900.0, 60e9, right_tau, band))
        .unwrap();
    body
}

fn base(
    cx: &Cx<'_>,
    gate: &CancelGate,
    body: SharedMaterialPlateChartBody,
) -> ReciprocalStringBank {
    base_with_sources(cx, gate, body, sources())
}

fn base_with_sources(
    cx: &Cx<'_>, gate: &CancelGate, body: SharedMaterialPlateChartBody,
    strings: Vec<MaterialStringSource>,
) -> ReciprocalStringBank {
    let mut config = options();
    config.network.max_modes = 64;
    config.network.max_connections = 64;
    config.network.max_setup_terms = 1_000_000;
    let medium = gas(101_325.0);
    ReciprocalStringBank::from_material_plate_chart(
        cx,
        strings,
        medium,
        body,
        Some(cavity(medium)),
        law(),
        config,
        gate,
    )
    .unwrap()
}

fn memory(bank: ReciprocalStringBank, gate: &CancelGate) -> ReciprocalStringBank {
    bank.with_material_plate_relaxation(&properties(), InitialPlateMemory::Relaxed, limits(), gate)
        .unwrap()
}

#[test]
fn fast_regional_memory_reaches_the_original_string_body_cavity_trajectory() {
    with_cx(|cx, gate| {
        let tau = 1e-7;
        let mut bank = memory(base(cx, gate, memory_body(tau, (0.0, 100_000.0))), gate);
        assert!(bank.system().sample_period_s() / tau > 200.0);
        assert!(bank.system().material_relaxation_branches().iter()
            .any(|branch| branch.relaxation_time_s == tau));
        let zero = vec![0.0; bank.system().mode_count()];
        let mut material_heat = 0.0;
        let mut peak_pressure = 0.0_f64;
        for _ in 0..32 {
            let frame = *bank.step(&zero, drive()).unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            peak_pressure = peak_pressure.max(frame.observer_pressure_pa.abs());
            material_heat += bank.system().coupled_frame().unwrap().material_dissipation_j;
        }
        assert!(peak_pressure > 0.0);
        assert!(material_heat > 0.0);
        assert!(bank.system().material_memory().iter().any(|z| *z != 0.0));
    });
}

#[test]
fn original_string_and_regional_body_memories_survive_gas_rebuild_together() {
    use fs_scenario::acoustic::BendingRelaxationProperties;
    with_cx(|cx, gate| {
        let make = || {
            let strings = sources().into_iter().map(|source| {
                let MaterialStringSource::Homogeneous(specimen) = source else { panic!("fixture") };
                let mut geometry = specimen.string();
                geometry.kelvin_voigt_bending = None;
                let source = with_uniform_circular_material_state(geometry, specimen.radius_m(),
                    &memory_material(7800.0, 40e9, 0.002, (0.0, 100_000.0))).unwrap()
                    .with_prony_bending_loss(&[BendingRelaxationProperties {
                        modulus: "maxwell_modulus".into(), relaxation_time: "maxwell_time".into(),
                    }]).unwrap();
                MaterialStringSource::Homogeneous(source)
            }).collect();
            base_with_sources(cx, gate, memory_body(0.004, (0.0, 100_000.0)), strings)
        };
        let before = make();
        assert_eq!(before.system().material_memory().len(), 7);
        let original_history = before.system().material_memory().to_vec();
        let body_first = memory(before, gate);
        assert_eq!(&body_first.system().material_memory()[..7], original_history);
        let mut a = extend(body_first, false, gate);
        let mut b = memory(extend(make(), false, gate), gate);
        assert_eq!(a.system().mode_count(), 25);
        assert_eq!(a.system().material_memory().len(), 13);
        let image_a = a.source().material_chart_relaxation().unwrap();
        let image_b = b.source().material_chart_relaxation().unwrap();
        assert_eq!(image_a.source().chart_identity(), image_b.source().chart_identity());
        assert_eq!(image_a.coupling_rows_per_s(), image_b.coupling_rows_per_s());
        for (index, (left, right)) in a.system().material_relaxation_branches().iter()
            .zip(b.system().material_relaxation_branches()).enumerate() {
            assert_eq!(left.projection, right.projection);
            assert_eq!(left.stiffness, right.stiffness);
            assert_eq!(left.relaxation_time_s, right.relaxation_time_s);
            assert!(left.projection.iter().enumerate().all(|(i, &v)| v == 0.0
                || (i % 2 == 0 && if index < 7 { i/2 < 7 } else { (7..10).contains(&(i/2)) })));
        }
        let force = vec![0.0; a.system().mode_count()];
        let mut pressure = 0.0_f64;
        for _ in 0..48 {
            let fa = *a.step(&force, drive()).unwrap();
            let fb = *b.step(&force, drive()).unwrap();
            assert_eq!(fa, fb);
            assert!(fa.energy_residual_j.abs() <= fa.energy_tolerance_j);
            assert_eq!(a.system().material_memory(), b.system().material_memory());
            pressure = pressure.max(fa.observer_pressure_pa.abs());
        }
        assert!(pressure > 0.0);
        assert!(a.system().material_memory()[..7].iter().any(|z| *z != 0.0));
        assert!(a.system().material_memory()[7..].iter().any(|z| *z != 0.0));
        assert!(a.system().coupled_frame().unwrap().material_dissipation_j > 0.0);
        assert_eq!(a.source().gas(), a.source().material_cavity().unwrap().cavity.gas());
    });
}

fn extend(bank: ReciprocalStringBank, axial: bool, gate: &CancelGate) -> ReciprocalStringBank {
    if axial {
        bank.with_longitudinal_waves(
            vec![
                LongitudinalStringResolution {
                    axial_modes: vec![1],
                    maximum_slope: 0.2,
                    maximum_strain: 0.01,
                };
                2
            ],
            NonlinearModalConfig {
                max_iterations: 64,
                maximum_generalized_force: 1e5,
                force_absolute_tolerance: 1e-11,
                force_relative_tolerance: 1e-10,
            },
            gate,
        )
        .unwrap()
    } else {
        bank.with_cylinder_gas_loading(
            CylinderFlowResolution {
                minimum_wavenumber_radius: 1e-3,
                maximum_wavenumber_radius: 30.0,
                intervals: 2,
            },
            gate,
        )
        .unwrap()
    }
}

fn assert_same_body_and_memory(a: &ReciprocalStringBank, b: &ReciprocalStringBank) {
    let a_image = a.source().material_chart_relaxation().unwrap();
    let b_image = b.source().material_chart_relaxation().unwrap();
    assert_eq!(
        a_image.source().chart_identity(),
        b_image.source().chart_identity()
    );
    assert_eq!(
        a_image.source().chart_identity(),
        a.source()
            .material_chart()
            .unwrap()
            .specimen
            .specimen_identity()
    );
    assert_eq!(a_image.coupling_rows_per_s(), b_image.coupling_rows_per_s());
    assert_eq!(a_image.relaxation_times_s(), b_image.relaxation_times_s());
    assert_eq!(a_image.source().properties(), properties());
    assert_eq!(
        a_image.source().omega_bands_rad_s(),
        b_image.source().omega_bands_rad_s()
    );
    assert_eq!(
        a_image.coupling_rows_per_s().len(),
        6,
        "a regional strain history was discarded"
    );
    assert_eq!(a.system().material_memory(), b.system().material_memory());
    assert_eq!(a.system().mode_count(), b.system().mode_count());
    assert_eq!(
        a.system().acceleration_pressure_weights(),
        b.system().acceleration_pressure_weights()
    );
    assert_eq!(a.system().material_relaxation_branches().len(), 6);
    assert_eq!(b.system().material_relaxation_branches().len(), 6);
    for (a_branch, b_branch) in a
        .system()
        .material_relaxation_branches()
        .iter()
        .zip(b.system().material_relaxation_branches())
    {
        assert_eq!(a_branch.projection, b_branch.projection);
        assert_eq!(a_branch.stiffness, b_branch.stiffness);
        assert_eq!(a_branch.relaxation_time_s, b_branch.relaxation_time_s);
        assert_eq!(a_branch.projection.len(), 2 * a.system().mode_count());
        // The two unequal string bases occupy seven modes, followed by the
        // three original body modes. New fluid/axial states cannot steal a row.
        assert!(
            a_branch
                .projection
                .iter()
                .enumerate()
                .all(|(i, &value)| value == 0.0 || (i % 2 == 0 && (7..10).contains(&(i / 2))))
        );
    }
    assert_eq!(a.system().components().len(), b.system().components().len());
    for (a_component, b_component) in a.system().components().iter().zip(b.system().components()) {
        assert_eq!(a_component.modes(), b_component.modes());
        assert_eq!(a_component.states(), b_component.states());
    }
}

#[test]
fn gas_rebuild_and_axial_append_preserve_body_memory_in_either_attachment_order() {
    with_cx(|cx, gate| {
        for axial in [false, true] {
            let body = memory_body(0.004, (0.0, 100_000.0));
            let mut a = extend(memory(base(cx, gate, body.clone()), gate), axial, gate);
            let mut b = memory(extend(base(cx, gate, body), axial, gate), gate);
            assert_same_body_and_memory(&a, &b);
            assert!(a.system().mode_count() > 11);
            let zero = vec![0.0; a.system().mode_count()];
            let mut work = 0.0;
            let mut loss = 0.0;
            let mut pressure = 0.0_f64;
            for sample in 0..48 {
                let drive = if sample < 24 {
                    drive()
                } else {
                    FrictionDrive {
                        speed_m_s: 0.0,
                        normal_force_n: 0.0,
                    }
                };
                let fa = *a.step(&zero, drive).unwrap();
                let fb = *b.step(&zero, drive).unwrap();
                assert_eq!(
                    fa, fb,
                    "attachment order changed accepted sample {sample}, axial={axial}"
                );
                assert!(fa.energy_residual_j.abs() <= fa.energy_tolerance_j);
                pressure = pressure.max(fa.observer_pressure_pa.abs());
                work += fa.external_work_j + fa.drive_work_j;
                loss += fa.network_dissipation_j + fa.friction_dissipation_j;
                assert_same_body_and_memory(&a, &b);
            }
            assert!(pressure > 0.0);
            assert!(a.system().material_memory().iter().any(|z| *z != 0.0));
            assert!(a.system().material_energy_j().unwrap() > 0.0);
            assert!(
                a.system().components()[1]
                    .states()
                    .iter()
                    .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0),
                "the shared body lost its sympathetic string reaction"
            );
            assert!((a.system().total_energy_j().unwrap() + loss - work).abs() < 1e-10);
            if axial {
                assert_eq!(
                    a.longitudinal_observation(0).unwrap(),
                    b.longitudinal_observation(0).unwrap()
                );
                assert!(a.source().longitudinal_waves().is_some());
            } else {
                assert!(a.source().cylinder_gas_loading().is_some());
            }
        }
    });
}

#[test]
fn changing_only_one_regional_relaxation_time_changes_pressure_and_neighbor_motion() {
    with_cx(|cx, gate| {
        let mut a = memory(base(cx, gate, memory_body(0.001, (0.0, 100_000.0))), gate);
        let mut b = memory(base(cx, gate, memory_body(0.008, (0.0, 100_000.0))), gate);
        let a_image = a.source().material_chart_relaxation().unwrap();
        let b_image = b.source().material_chart_relaxation().unwrap();
        assert_ne!(
            a_image.source().chart_identity(),
            b_image.source().chart_identity()
        );
        assert_eq!(a_image.coupling_rows_per_s(), b_image.coupling_rows_per_s());
        assert_ne!(a_image.relaxation_times_s(), b_image.relaxation_times_s());
        assert_eq!(
            a.source().material_chart().unwrap().specimen.regions()[0]
                .input
                .material
                .identity(),
            b.source().material_chart().unwrap().specimen.regions()[0]
                .input
                .material
                .identity()
        );
        for (a, b) in a.system().components().iter().zip(b.system().components()) {
            assert_eq!(
                a.modes(),
                b.modes(),
                "changing memory changed equilibrium mass or stiffness"
            );
        }
        assert_eq!(
            a.system().acceleration_pressure_weights(),
            b.system().acceleration_pressure_weights()
        );
        let zero = vec![0.0; a.system().mode_count()];
        let mut pressure_difference = 0.0_f64;
        let mut neighbor_difference = 0.0_f64;
        for _ in 0..128 {
            let fa = *a.step(&zero, drive()).unwrap();
            let fb = *b.step(&zero, drive()).unwrap();
            assert!(fa.energy_residual_j.abs() <= fa.energy_tolerance_j);
            assert!(fb.energy_residual_j.abs() <= fb.energy_tolerance_j);
            pressure_difference =
                pressure_difference.max((fa.observer_pressure_pa - fb.observer_pressure_pa).abs());
            for (a, b) in a.system().components()[1]
                .states()
                .iter()
                .zip(b.system().components()[1].states())
            {
                neighbor_difference = neighbor_difference
                    .max((a.displacement_m_sqrt_kg - b.displacement_m_sqrt_kg).abs());
            }
        }
        assert!(
            pressure_difference > 1e-15,
            "regional material history never reached actual pressure"
        );
        assert!(
            neighbor_difference > 1e-20,
            "regional history did not react back through the physical bridge"
        );
        assert_ne!(a.system().material_memory(), b.system().material_memory());
    });
}

#[test]
fn duplicate_or_unsupported_sources_refuse_and_cancellation_keeps_accepted_memory() {
    with_cx(|cx, gate| {
        let body = memory_body(0.004, (0.0, 100_000.0));
        let duplicate = vec![properties()[0].clone(); 2];
        assert!(
            base(cx, gate, body.clone())
                .with_material_plate_relaxation(
                    &duplicate,
                    InitialPlateMemory::Relaxed,
                    limits(),
                    gate
                )
                .is_err()
        );
        assert!(
            memory(base(cx, gate, body.clone()), gate)
                .with_material_plate_relaxation(
                    &properties(),
                    InitialPlateMemory::Relaxed,
                    limits(),
                    gate
                )
                .is_err()
        );
        for band in [(10.0, 100_000.0), (0.0, 1000.0)] {
            // Every query coefficient is valid at omega=1000; neither band
            // supports the actual coupled owner's DC-to-instantaneous range.
            assert!(
                base(cx, gate, memory_body(0.004, band))
                    .with_material_plate_relaxation(
                        &properties(),
                        InitialPlateMemory::Relaxed,
                        limits(),
                        gate
                    )
                    .is_err()
            );
        }
        let stop = CancelGate::new_clock_free();
        stop.request();
        assert!(
            base(cx, gate, body.clone())
                .with_material_plate_relaxation(
                    &properties(),
                    InitialPlateMemory::Relaxed,
                    limits(),
                    &stop
                )
                .is_err()
        );
        let mut bank = memory(base(cx, gate, body.clone()), gate);
        let mut reference = memory(base(cx, gate, body), gate);
        let zero = vec![0.0; bank.system().mode_count()];
        for _ in 0..8 {
            assert_eq!(
                bank.step(&zero, drive()).unwrap(),
                reference.step(&zero, drive()).unwrap()
            );
        }
        let history = bank.system().material_memory().to_vec();
        assert!(history.iter().any(|z| *z != 0.0));
        let energy = bank.system().total_energy_j().unwrap();
        let count = bank.system().samples_rendered();
        assert!(bank.step_under_gate(&zero, drive(), &stop).is_err());
        assert_eq!(bank.system().samples_rendered(), count);
        assert_eq!(bank.system().material_memory(), history);
        assert_eq!(bank.system().total_energy_j().unwrap(), energy);
        assert_same_body_and_memory(&bank, &reference);
        assert_eq!(
            bank.step(&zero, drive()).unwrap(),
            reference.step(&zero, drive()).unwrap()
        );
        assert_same_body_and_memory(&bank, &reference);
    });
}
