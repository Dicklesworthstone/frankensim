//! G1/G3/G4: original interaction kernels, interpolation and source support.

use super::*;
use crate::radiation::transport::{
    PhotonOrdinate, PhotonSlabBoundary, PhotonSlabSource, PhotonSourceGeometry,
    prepare_slab_response_with_sources,
};

#[test]
fn g1_affine_coefficients_match_an_independent_polynomial() {
    for (a, b) in [(0.0, 4.0), (4.0, 0.0), (1.0, 9.0), (9.0, 1.0)] {
        for index in 0..=64 {
            let t = 300.0 + index as f64;
            let expected = ((364.0 - t) * a + (t - 300.0) * b) / 64.0;
            assert_eq!(interpolate(a, b, 300.0, 364.0, t), Some(expected));
        }
    }
}

#[test]
fn g3_interpolation_retains_weak_branches_and_constant_extremes() {
    let temperature = 500.0_f64.next_down();
    let weak = interpolate(1.0, 0.0, 300.0, 500.0, temperature).unwrap();
    assert!(weak > 0.0);
    assert_eq!(weak, (500.0 - temperature) / 200.0);
    for coefficient in [0.0, f64::from_bits(1), 1.0, f64::MAX] {
        assert_eq!(
            interpolate(coefficient, coefficient, 300.0, 500.0, 400.0),
            Some(coefficient)
        );
    }
    assert_eq!(
        interpolate(0.0, f64::MAX, 300.0, 500.0, 400.0),
        Some(0.5 * f64::MAX)
    );
    assert_eq!(
        interpolate(f64::MAX, 0.0, 300.0, 500.0, 400.0),
        Some(0.5 * f64::MAX)
    );
}

#[test]
fn g4_interpolation_refuses_lost_positive_coefficients() {
    assert_eq!(
        interpolate(0.0, f64::from_bits(1), 300.0, 500.0, 301.0),
        None
    );
    assert_eq!(
        interpolate(0.0, f64::from_bits(1), 300.0, 500.0, 500.0),
        Some(f64::from_bits(1))
    );
    assert_eq!(interpolate(0.0, 0.0, 300.0, 500.0, 301.0), Some(0.0));
}

fn material(id: &str, hot: bool) -> PhotonMassMaterial {
    PhotonMassMaterial::new(
        id.into(),
        37,
        DecayRadiation::Gamma,
        vec![0.0, 1.0, 2.0],
        if hot {
            vec![
                PhotonMassCollision {
                    total_m2_kg: 4.0,
                    heat_m2_kg: 1.0,
                    redistribution_m2_kg: vec![3.0, 0.0],
                },
                PhotonMassCollision {
                    total_m2_kg: 8.0,
                    heat_m2_kg: 2.0,
                    redistribution_m2_kg: vec![1.0, 5.0],
                },
            ]
        } else {
            vec![
                PhotonMassCollision {
                    total_m2_kg: 2.0,
                    heat_m2_kg: 2.0,
                    redistribution_m2_kg: vec![0.0, 0.0],
                },
                PhotonMassCollision {
                    total_m2_kg: 4.0,
                    heat_m2_kg: 1.0,
                    redistribution_m2_kg: vec![2.0, 1.0],
                },
            ]
        },
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap()
}
fn knots() -> Vec<PhotonMassTemperatureKnot> {
    vec![
        PhotonMassTemperatureKnot {
            temperature_k: 300.0,
            material: material("cold original source", false),
        },
        PhotonMassTemperatureKnot {
            temperature_k: 500.0,
            material: material("hot original source", true),
        },
    ]
}
fn table() -> TemperaturePhotonMassMaterial {
    TemperaturePhotonMassMaterial::new(
        "authored temperature interpolation".into(),
        knots(),
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn g1_source_temperature_changes_heat_and_outgoing_energy_together() {
    let table = table();
    assert_eq!(table.material_id(), 37);
    assert_eq!(table.temperature_bounds_k(), (300.0, 500.0));
    assert_eq!(table.knots()[1].material.data_id(), "hot original source");
    assert_eq!(
        table.sample(300.0, || Ok(())).unwrap().data_id(),
        "cold original source"
    );
    let midpoint = table.sample(400.0, || Ok(())).unwrap();
    assert_eq!(midpoint.data_id(), table.data_id());
    assert_eq!(midpoint.groups()[0].total_m2_kg, 3.0);
    assert_eq!(midpoint.groups()[0].heat_m2_kg, 1.5);
    assert_eq!(midpoint.groups()[0].redistribution_m2_kg, [1.5, 0.0]);
    assert_eq!(midpoint.groups()[1].total_m2_kg, 6.0);
    assert_eq!(midpoint.groups()[1].heat_m2_kg, 1.5);
    assert_eq!(midpoint.groups()[1].redistribution_m2_kg, [1.5, 3.0]);
    let cell = midpoint.macroscopic_cell(91, 0.25, 3.0, || Ok(())).unwrap();
    assert_eq!(cell.groups[1].heat_per_m, 4.5);
    assert_eq!(cell.groups[1].redistribution_per_m, [4.5, 9.0]);
    assert_eq!(cell.groups[1].total_per_m, 18.0);
    assert_eq!(
        table.sample(500.0, || Ok(())).unwrap().groups(),
        material("hot", true).groups()
    );
}

#[test]
fn g4_temperature_support_identity_budgets_and_callbacks_are_retained() {
    let table = table();
    for temperature in [f64::NAN, 0.0, 299.9, 500.1] {
        assert!(table.sample(temperature, || Ok(())).is_err());
    }
    let mut wrong = knots();
    wrong[1].temperature_k = 300.0;
    assert!(
        TemperaturePhotonMassMaterial::new(
            "unordered".into(),
            wrong,
            PhotonSlabLimits::default(),
            || Ok(())
        )
        .is_err()
    );
    let mut wrong = knots();
    wrong[1].material = PhotonMassMaterial::new(
        "different material".into(),
        38,
        DecayRadiation::Gamma,
        vec![0.0, 1.0, 2.0],
        material("hot", true).groups().to_vec(),
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap();
    assert!(
        TemperaturePhotonMassMaterial::new(
            "mixed identity".into(),
            wrong,
            PhotonSlabLimits::default(),
            || Ok(())
        )
        .is_err()
    );
    assert!(matches!(
        TemperaturePhotonMassMaterial::new(
            "limited".into(),
            knots(),
            PhotonSlabLimits {
                max_values: 23,
                ..PhotonSlabLimits::default()
            },
            || Ok(())
        ),
        Err(RadiationError::Budget(_))
    ));
    let expected = RadiationError::InvalidInput("original interpolation callback");
    let mut calls = 0;
    let result = table.sample(400.0, || {
        calls += 1;
        if calls == 3 {
            Err(expected.clone())
        } else {
            Ok(())
        }
    });
    assert_eq!(result.unwrap_err(), expected);
    assert_eq!(
        table.sample(400.0, || Ok(())).unwrap().groups()[0].heat_m2_kg,
        1.5
    );
}

fn mixture_source(
    material_id: u64,
    points: &[(f64, f64)],
    scattering: bool,
) -> TemperaturePhotonMassMaterial {
    TemperaturePhotonMassMaterial::new(
        format!("manufactured temperature source {material_id}"),
        points
            .iter()
            .map(|&(temperature_k, x)| PhotonMassTemperatureKnot {
                temperature_k,
                material: PhotonMassMaterial::new(
                    format!("manufactured source {material_id} at {temperature_k} K"),
                    material_id,
                    DecayRadiation::Gamma,
                    vec![0.0, 1.0, 2.0],
                    if scattering {
                        vec![
                            PhotonMassCollision {
                                total_m2_kg: 4.0,
                                heat_m2_kg: x,
                                redistribution_m2_kg: vec![4.0 - x, 0.0],
                            },
                            PhotonMassCollision {
                                total_m2_kg: 8.0,
                                heat_m2_kg: 2.0 * x,
                                redistribution_m2_kg: vec![4.0 - x, 4.0 - x],
                            },
                        ]
                    } else {
                        [x, 2.0 * x]
                            .map(|coefficient| PhotonMassCollision {
                                total_m2_kg: coefficient,
                                heat_m2_kg: coefficient,
                                redistribution_m2_kg: vec![0.0, 0.0],
                            })
                            .to_vec()
                    },
                    PhotonSlabLimits::default(),
                    || Ok(()),
                )
                .unwrap(),
            })
            .collect(),
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

fn composition(scattering: bool) -> Vec<TemperaturePhotonMassConstituent> {
    vec![
        TemperaturePhotonMassConstituent {
            mass_fraction: 0.25,
            material: mixture_source(
                11,
                &[(200.0, 0.0), (300.0, 2.0), (400.0, 4.0), (600.0, 0.0)],
                scattering,
            ),
        },
        TemperaturePhotonMassConstituent {
            mass_fraction: 0.75,
            material: mixture_source(22, &[(300.0, 4.0), (500.0, 0.0), (700.0, 4.0)], scattering),
        },
    ]
}

const MIXTURE_BASIS: &str = "Manufactured independent kernels with common group weighting; \
    the supplied local deposition and outgoing energy apply in this mixture at every admitted kelvin.";

fn mixture(scattering: bool) -> TemperaturePhotonMassMaterial {
    TemperaturePhotonMassMaterial::from_mass_fractions(
        "manufactured temperature mixture".into(),
        99,
        composition(scattering),
        MIXTURE_BASIS.into(),
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

// Independently evaluated 1/4 A(T) + 3/4 B(T), including both changes of slope.
const MIXED_VALUES: [(f64, f64); 7] = [
    (300.0, 3.5),
    (350.0, 3.0),
    (400.0, 2.5),
    (450.0, 1.5),
    (500.0, 0.5),
    (550.0, 1.0),
    (600.0, 1.5),
];

#[test]
fn g1_temperature_mixture_preserves_original_knots_and_complete_energy_kernels() {
    let sources = composition(true);
    let table = TemperaturePhotonMassMaterial::from_mass_fractions(
        "manufactured temperature mixture".into(),
        99,
        sources.clone(),
        MIXTURE_BASIS.into(),
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(table.material_id(), 99);
    assert_eq!(table.temperature_bounds_k(), (300.0, 600.0));
    assert_eq!(
        table
            .knots()
            .iter()
            .map(|k| k.temperature_k)
            .collect::<Vec<_>>(),
        [300.0, 400.0, 500.0, 600.0]
    );
    let retained = table.mixture().unwrap();
    assert_eq!(retained.kernel_applicability, MIXTURE_BASIS);
    for (original, retained) in sources.iter().zip(&retained.constituents) {
        assert!(Arc::ptr_eq(&original.material.0, &retained.material.0));
        assert_eq!(original.mass_fraction, retained.mass_fraction);
    }
    for (temperature, x) in MIXED_VALUES {
        let sampled = table.sample(temperature, || Ok(())).unwrap();
        let rows = sampled.groups();
        assert_eq!(rows[0].total_m2_kg, 4.0);
        assert_eq!(rows[0].heat_m2_kg, x);
        assert_eq!(rows[0].redistribution_m2_kg, [4.0 - x, 0.0]);
        assert_eq!(rows[1].total_m2_kg, 8.0);
        assert_eq!(rows[1].heat_m2_kg, 2.0 * x);
        assert_eq!(rows[1].redistribution_m2_kg, [4.0 - x, 4.0 - x]);
    }
    let at_knot = table.sample(400.0, || Ok(())).unwrap();
    assert!(Arc::ptr_eq(&at_knot.0, &table.knots()[1].material.0));
    assert_eq!(at_knot.mixture().unwrap().constituents.len(), 2);
    // Interpolating just the mixed endpoints would erase the 500 K minimum.
    let endpoint_chord = 3.5 + (1.5 - 3.5) * (500.0 - 300.0) / (600.0 - 300.0);
    assert!(endpoint_chord - table.sample(500.0, || Ok(())).unwrap().groups()[0].heat_m2_kg > 1.0);
    for outside in [299.0, 601.0, f64::NAN] {
        assert!(table.sample(outside, || Ok(())).is_err());
    }
}

#[test]
fn g1_temperature_mixture_heating_follows_mass_thickness_and_constituent_breakpoints() {
    let table = mixture(false);
    for (temperature, coefficient) in MIXED_VALUES {
        for mu in [0.4, 1.0] {
            let density = 2.0;
            let depth = 0.3;
            let kernel = table.sample(temperature, || Ok(())).unwrap();
            let prepared = prepare_slab_response_with_sources(
                "manufactured temperature-mixture irradiation".into(),
                DecayRadiation::Gamma,
                table.group_edges_j().to_vec(),
                vec![
                    kernel
                        .macroscopic_cell(91, depth, density, || Ok(()))
                        .unwrap(),
                ],
                vec![PhotonSlabSource {
                    source_id: 17,
                    geometry: PhotonSourceGeometry::Collimated {
                        boundary: PhotonSlabBoundary::Left,
                        mu,
                    },
                }],
                vec![PhotonOrdinate {
                    mu: 1.0,
                    weight: 1.0,
                }],
                PhotonSlabLimits::default(),
                || Ok(()),
            )
            .unwrap();
            let result = prepared
                .response
                .apply(&[(17, &[5.0, 2.0])], 1e-12, || Ok(()))
                .unwrap();
            let tau = coefficient * density * depth / mu;
            let heat = -5.0 * (-tau).exp_m1() - 2.0 * (-2.0 * tau).exp_m1();
            let escaped = 5.0 * (-tau).exp() + 2.0 * (-2.0 * tau).exp();
            assert!((result.deposited_j - heat).abs() < 2e-13);
            assert!((result.escaped_j - escaped).abs() < 2e-13);
            assert!((result.receiver_energy_j[&91] - heat).abs() < 2e-13);
        }
    }
}

#[test]
fn g3_temperature_mixture_repartition_preserves_the_same_piecewise_law() {
    let original = mixture(true);
    let mut sources = composition(true);
    sources[0].mass_fraction = 0.125;
    sources.push(sources[0].clone());
    sources.reverse();
    let repartitioned = TemperaturePhotonMassMaterial::from_mass_fractions(
        "same temperature-dependent physical mixture".into(),
        99,
        sources,
        MIXTURE_BASIS.into(),
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap();
    for temperature in [300.0, 325.0, 400.0, 475.0, 500.0, 575.0, 600.0] {
        assert_eq!(
            original.sample(temperature, || Ok(())).unwrap().groups(),
            repartitioned
                .sample(temperature, || Ok(()))
                .unwrap()
                .groups()
        );
    }
    assert_eq!(repartitioned.knots().len(), 4);
}

#[test]
fn g4_temperature_mixture_refuses_incomplete_composition_and_unsupported_data() {
    let admit = |sources, basis: &str, limits| {
        TemperaturePhotonMassMaterial::from_mass_fractions(
            "manufactured temperature mixture".into(),
            99,
            sources,
            basis.into(),
            limits,
            || Ok(()),
        )
    };
    for fraction in [-0.1, 0.1, 1.1, f64::NAN, f64::INFINITY] {
        let mut sources = composition(false);
        sources[0].mass_fraction = fraction;
        assert!(matches!(
            admit(sources, MIXTURE_BASIS, PhotonSlabLimits::default()),
            Err(RadiationError::InvalidInput(_))
        ));
    }
    assert!(admit(Vec::new(), MIXTURE_BASIS, PhotonSlabLimits::default()).is_err());
    assert!(admit(composition(false), " ", PhotonSlabLimits::default()).is_err());
    for minimum in [600.0, 650.0] {
        let mut sources = composition(false);
        sources[1].material = mixture_source(22, &[(minimum, 1.0), (700.0, 2.0)], false);
        assert!(matches!(
            admit(sources, MIXTURE_BASIS, PhotonSlabLimits::default()),
            Err(RadiationError::InvalidInput(_))
        ));
    }
    let mut sources = composition(false);
    sources[1].material = TemperaturePhotonMassMaterial::new(
        "manufactured mismatched energy grid".into(),
        [300.0, 600.0]
            .map(|temperature_k| PhotonMassTemperatureKnot {
                temperature_k,
                material: PhotonMassMaterial::new(
                    "different grid".into(),
                    22,
                    DecayRadiation::Gamma,
                    vec![0.0, 1.0, 3.0],
                    sources[1].material.knots()[0].material.groups().to_vec(),
                    PhotonSlabLimits::default(),
                    || Ok(()),
                )
                .unwrap(),
            })
            .to_vec(),
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap();
    assert!(matches!(
        admit(sources, MIXTURE_BASIS, PhotonSlabLimits::default()),
        Err(RadiationError::InvalidInput(_))
    ));
    let table = mixture(false);
    let limits = PhotonSlabLimits {
        max_values: table.value_count(),
        ..PhotonSlabLimits::default()
    };
    assert!(admit(composition(false), MIXTURE_BASIS, limits).is_ok());
    for limits in [
        PhotonSlabLimits {
            max_values: table.value_count() - 1,
            ..limits
        },
        PhotonSlabLimits {
            max_work: 1,
            ..limits
        },
    ] {
        assert!(matches!(
            admit(composition(false), MIXTURE_BASIS, limits),
            Err(RadiationError::Budget(_))
        ));
    }
    let expected = RadiationError::InvalidInput("original temperature-mixture cancellation");
    let mut calls = 0;
    let result = TemperaturePhotonMassMaterial::from_mass_fractions(
        "cancelled mixture".into(),
        99,
        composition(false),
        MIXTURE_BASIS.into(),
        PhotonSlabLimits::default(),
        || {
            calls += 1;
            if calls == 30 {
                Err(expected.clone())
            } else {
                Ok(())
            }
        },
    );
    assert_eq!(result.unwrap_err(), expected);
    assert_eq!(
        table.sample(500.0, || Ok(())).unwrap().groups()[0].heat_m2_kg,
        0.5
    );
}
