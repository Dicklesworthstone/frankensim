//! Local source laws versus the independent full complex weak equation.
use super::*;
use fs_duct::modal::time_domain::{ModalTimeError, assemble_lossless_field};
use fs_matdb::QueryPoint;
use fs_material::gas_field::{AxialGasField, PrimitiveGasState};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::thermochemical_gas::packs::from_packs;
use fs_math::det;
use fs_thermochem::{Composition, CompositionBasis};

#[path = "../../../fs-couple/src/render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

fn field(length: f64, cellwise: bool) -> AxialGasField {
    let a = gas_support::packs("A", 0.028, 3.5, 200.0);
    let b = gas_support::packs("B", 0.044, 4.5, 200.0);
    let primitive = |t, f| PrimitiveGasState {
        temperature_k: t,
        composition: Composition::new(CompositionBasis::MoleFraction, vec![f, 1.0 - f]).unwrap(),
    };
    let model = from_packs(
        &[a.binding(), b.binding()],
        primitive(300.0, 0.8).composition,
        None,
    )
    .unwrap();
    let (stations, states) = if cellwise {
        (
            vec![0.0, length / 2.0, length],
            vec![
                primitive(300.0, 1.0),
                primitive(300.0, 1.0),
                primitive(600.0, 0.0),
            ],
        )
    } else {
        (
            vec![0.0, length],
            vec![primitive(300.0, 0.8), primitive(600.0, 0.2)],
        )
    };
    let constructor = if cellwise {
        AxialGasField::new_piecewise_constant
    } else {
        AxialGasField::new
    };
    constructor(
        stations,
        101325.0,
        states,
        model,
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap()
}

fn spatial_law(a: bool, b: bool) -> RelaxationAbsorption {
    let process = |theta| RelaxationProcess {
        reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
            theta_k: theta,
            degeneracy: 1,
        }]),
        rate: RelaxationRateLaw::PartnerIndependent {
            p_tau_ref_pa_s: 101325.0 * 2e-5,
            t_ref_k: 300.0,
            exponent: 1.2,
            t_min_k: 200.0,
            t_max_k: 1000.0,
        },
    };
    RelaxationAbsorption::new(vec![
        SpeciesRelaxation {
            species: SpeciesId::new("A").unwrap(),
            // Exactly equal times remain distinct; the last capacity underflows
            // physically to zero and remains in the original evaluation.
            processes: if a {
                vec![process(1000.0), process(1500.0), process(1e7)]
            } else {
                vec![]
            },
        },
        SpeciesRelaxation {
            species: SpeciesId::new("B").unwrap(),
            processes: if b {
                vec![RelaxationProcess {
                    reservoir: InternalReservoir::Rotational {
                        degrees_of_freedom: 2,
                    },
                    rate: RelaxationRateLaw::PartnerIndependent {
                        p_tau_ref_pa_s: 101325.0 * 1e-4,
                        t_ref_k: 300.0,
                        exponent: -0.7,
                        t_min_k: 200.0,
                        t_max_k: 1000.0,
                    },
                }]
            } else {
                vec![]
            },
        },
    ])
    .unwrap()
}

#[test]
fn local_caloric_fields_match_original_complex_compliance_and_physical_heat() {
    let original = field(0.035, false);
    let law = spatial_law(true, true);
    for flare in [0.0, 7.0] {
        let profile = RadiusProfile::exponential(0.05, flare, original.length_m()).unwrap();
        let image = assemble_lossless_field(
            &profile,
            &original,
            ModeCount(2),
            3,
            ModalEnd::Ports,
            124,
            || false,
        )
        .unwrap()
        .with_field_molecular_relaxation(law.clone(), budget(), || false)
        .unwrap();
        let source = image.source();
        let n = source.coordinates();
        assert_eq!(source.state_dimension(), 64);
        assert!(image.system().resistance_witness().is_some());
        assert_eq!(source.molecular_field_law(), Some(&law));
        assert_eq!(source.gas_field(), Some(&original));
        assert!(source.molecular_relaxation().is_none());
        for element in 0..3 {
            for point in 0..2 {
                let local = source
                    .molecular_relaxation_at_quadrature(element, point)
                    .unwrap();
                let expected = law
                    .evaluate(
                        &RelaxingMedium::from_thermochemical(
                            source.gas_at_quadrature(element, point).unwrap(),
                        )
                        .unwrap(),
                    )
                    .unwrap();
                assert_eq!(local, &expected);
                assert_eq!(local.processes().len(), 4);
                assert_eq!(
                    local.processes()[0].relaxation_time_s,
                    local.processes()[1].relaxation_time_s
                );
                assert_eq!(
                    source.molecular_process_indices_at_quadrature(element, point),
                    Some([0, 1, 3].as_slice())
                );
                assert_eq!(
                    source
                        .molecular_state_range_at_quadrature(element, point)
                        .unwrap()
                        .len(),
                    8
                );
                assert_eq!(
                    source
                        .molecular_heat_range_at_quadrature(element, point)
                        .unwrap()
                        .len(),
                    3
                );
            }
        }
        assert_ne!(
            source
                .molecular_relaxation_at_quadrature(0, 0)
                .unwrap()
                .processes()[0]
                .relaxation_time_s,
            source
                .molecular_relaxation_at_quadrature(2, 1)
                .unwrap()
                .processes()[0]
                .relaxation_time_s
        );
        let model = HarmonicModel::from_ode(image.system(), harmonic(), || false).unwrap();
        for omega in [100.0, 3000.0, 19000.0, 40000.0] {
            let drive = [
                C64::new(1e-5, 2e-6),
                C64::new(-3e-6, 4e-6),
                C64::new(2e-6, -1e-6),
                C64::new(1e-6, 3e-6),
            ];
            let response = model
                .factor(omega, || false)
                .unwrap()
                .response(&drive, || false)
                .unwrap();
            let mut dynamic: Vec<_> = source
                .stiffness()
                .iter()
                .copied()
                .map(C64::from_re)
                .collect();
            let mut local_heat = Vec::new();
            let h = profile.length() / 3.0;
            let offset = 0.5 / det::sqrt(3.0);
            // Assemble the original COMPLEX material law directly, independently
            // of all augmented states, frozen mass and pHS loss outer products.
            for element in 0..3 {
                for (point, fraction) in [0.5 - offset, 0.5 + offset].into_iter().enumerate() {
                    let x = (element as f64 + fraction) * h;
                    let radius = profile.radius(x).unwrap();
                    let volume = 0.5 * h * core::f64::consts::PI * radius * radius;
                    let local = law
                        .evaluate(
                            &RelaxingMedium::from_thermochemical(&original.state_at(x).unwrap())
                                .unwrap(),
                        )
                        .unwrap();
                    let compliance = local
                        .caloric_response(omega)
                        .unwrap()
                        .bulk_modulus_pa
                        .recip()
                        .scale(volume);
                    let shape = [1.0 - fraction, fraction];
                    for mode in 0..2 {
                        for a in 0..2 {
                            for b in 0..2 {
                                dynamic
                                    [((element + a) * 2 + mode) * n + (element + b) * 2 + mode] =
                                    dynamic
                                        [((element + a) * 2 + mode) * n + (element + b) * 2 + mode]
                                        - compliance.scale(omega * omega * shape[a] * shape[b]);
                            }
                        }
                    }
                    local_heat.push((element, point, shape, compliance.im));
                }
            }
            let mut expected = vec![C64::ZERO; n];
            expected[..2].copy_from_slice(&drive[..2]);
            expected[n - 2..].copy_from_slice(&drive[2..]);
            lu_complex(&dynamic, n).unwrap().solve(&mut expected);
            let pressure: Vec<_> = expected
                .iter()
                .map(|q| C64::new(0.0, -omega) * *q)
                .collect();
            for (port, index) in [0, 1, n - 2, n - 1].into_iter().enumerate() {
                close(response.output[port], pressure[index], 4e-8);
            }
            let heats = source
                .molecular_mean_heat_per_process_w(&response.state)
                .unwrap();
            assert_eq!(heats.len(), 18);
            let mut expected_total = 0.0;
            for (element, point, shape, imaginary_compliance) in local_heat {
                let expected = 0.5
                    * omega
                    * imaginary_compliance
                    * (0..2)
                        .map(|mode| {
                            (pressure[element * 2 + mode].scale(shape[0])
                                + pressure[(element + 1) * 2 + mode].scale(shape[1]))
                            .norm_sq()
                        })
                        .sum::<f64>();
                let actual = heats[source
                    .molecular_heat_range_at_quadrature(element, point)
                    .unwrap()]
                .iter()
                .sum();
                close(C64::from_re(actual), C64::from_re(expected), 5e-8);
                expected_total += expected;
            }
            close(
                C64::from_re(response.dissipated_power_w),
                C64::from_re(expected_total),
                5e-8,
            );
            assert!(response.power_residual_w.abs() < 2e-8 * expected_total);
            let real: Vec<_> = response.state.iter().map(|v| v.re).collect();
            let imaginary: Vec<_> = response.state.iter().map(|v| v.im).collect();
            let a = source.molecular_heat_per_process_w(&real).unwrap();
            let b = source.molecular_heat_per_process_w(&imaginary).unwrap();
            for i in 0..heats.len() {
                close(
                    C64::from_re(heats[i]),
                    C64::from_re(0.5 * (a[i] + b[i])),
                    2e-12,
                );
            }
        }
    }
}

#[test]
fn empty_and_locally_absent_reservoirs_keep_original_states_and_matrices() {
    let original = field(0.04, true);
    let profile = RadiusProfile::exponential(0.03, 2.0, 0.04).unwrap();
    let make = || {
        assemble_lossless_field(
            &profile,
            &original,
            ModeCount(1),
            2,
            ModalEnd::Ports,
            124,
            || false,
        )
        .unwrap()
    };
    let bare = make();
    let empty = make()
        .with_field_molecular_relaxation(spatial_law(false, false), budget(), || false)
        .unwrap();
    assert_eq!(bare.system().structure(), empty.system().structure());
    assert_eq!(
        bare.system().storage().quadratic_part().unwrap().matrix,
        empty.system().storage().quadratic_part().unwrap().matrix
    );
    assert_eq!(bare.source().mass(), empty.source().mass());
    assert_eq!(bare.source().inverse_mass(), empty.source().inverse_mass());
    assert_eq!(empty.source().state_dimension(), 6);
    let active = make()
        .with_field_molecular_relaxation(spatial_law(true, false), budget(), || false)
        .unwrap();
    assert_eq!(active.source().state_dimension(), 12);
    for point in 0..2 {
        assert_eq!(
            active
                .source()
                .molecular_state_range_at_quadrature(0, point)
                .unwrap()
                .len(),
            3
        );
        assert!(
            active
                .source()
                .molecular_state_range_at_quadrature(1, point)
                .unwrap()
                .is_empty()
        );
        assert!(
            active
                .source()
                .molecular_heat_range_at_quadrature(1, point)
                .unwrap()
                .is_empty()
        );
    }
    assert!(
        active
            .source()
            .molecular_relaxation_at_quadrature(2, 0)
            .is_none()
    );
    assert!(
        active
            .source()
            .molecular_relaxation_at_quadrature(0, 2)
            .is_none()
    );
    assert!(
        active
            .with_field_molecular_relaxation(spatial_law(true, false), budget(), || false)
            .is_err()
    );
}

#[test]
fn local_material_capacity_and_cancellation_fail_without_reduced_models() {
    let original = field(0.035, false);
    let profile = RadiusProfile::exponential(0.05, 0.0, 0.035).unwrap();
    let make = || {
        assemble_lossless_field(
            &profile,
            &original,
            ModeCount(2),
            3,
            ModalEnd::Ports,
            124,
            || false,
        )
        .unwrap()
    };
    for limits in [
        ModalMolecularBudget {
            maximum_states: 63,
            ..budget()
        },
        ModalMolecularBudget {
            maximum_payload_bytes: 1,
            ..budget()
        },
        ModalMolecularBudget {
            maximum_preparation_work: 1,
            ..budget()
        },
    ] {
        assert!(
            make()
                .with_field_molecular_relaxation(spatial_law(true, true), limits, || false)
                .is_err()
        );
    }
    let mut polls = 0;
    make()
        .with_field_molecular_relaxation(spatial_law(true, true), budget(), || {
            polls += 1;
            false
        })
        .unwrap();
    for limit in [1, 7, polls] {
        let mut seen = 0;
        let result =
            make().with_field_molecular_relaxation(spatial_law(true, true), budget(), || {
                seen += 1;
                seen == limit
            });
        assert!(matches!(
            result,
            Err(ModalTimeError::Duct(fs_duct::DuctError::Cancelled { .. }))
        ));
        assert_eq!(seen, limit);
    }
    let missing = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: SpeciesId::new("A").unwrap(),
        processes: vec![],
    }])
    .unwrap();
    assert!(matches!(
        make().with_field_molecular_relaxation(missing, budget(), || false),
        Err(ModalTimeError::Material(_))
    ));
    let mut entries = spatial_law(true, true).species().to_vec();
    for entry in &mut entries {
        for process in &mut entry.processes {
            if let RelaxationRateLaw::PartnerIndependent { t_max_k, .. } = &mut process.rate {
                *t_max_k = 400.0;
            }
        }
    }
    let limited = RelaxationAbsorption::new(entries).unwrap();
    assert!(matches!(
        make().with_field_molecular_relaxation(limited, budget(), || false),
        Err(ModalTimeError::Material(_))
    ));
    assert!(
        make()
            .with_field_molecular_relaxation(spatial_law(true, true), budget(), || false)
            .is_ok()
    );
    let unbound = assemble_lossless(
        &profile,
        &gas(),
        ModeCount(1),
        2,
        ModalEnd::Rigid,
        124,
        || false,
    )
    .unwrap();
    assert!(
        unbound
            .with_field_molecular_relaxation(spatial_law(true, true), budget(), || false)
            .is_err()
    );
}
