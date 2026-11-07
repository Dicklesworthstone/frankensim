//! G1/G3/G4: physical density conversion, optical-depth scaling and refusal.

use super::*;
use crate::radiation::transport::{
    PhotonOrdinate, PhotonSlabBoundary, PhotonSlabSource, PhotonSourceGeometry, prepare_slab_response,
    prepare_slab_response_with_sources,
};

fn rows() -> Vec<PhotonMassCollision> {
    vec![
        PhotonMassCollision {
            total_m2_kg: 2.0,
            heat_m2_kg: 0.5,
            redistribution_m2_kg: vec![1.5, 0.0],
        },
        PhotonMassCollision {
            total_m2_kg: 4.0,
            heat_m2_kg: 1.0,
            redistribution_m2_kg: vec![2.0, 1.0],
        },
    ]
}

fn material() -> PhotonMassMaterial {
    PhotonMassMaterial::new(
        "manufactured mass collision data".into(),
        37,
        DecayRadiation::Gamma,
        vec![0.0, 1.0, 2.0],
        rows(),
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn g1_original_density_scales_distinct_heat_and_photon_energy_coefficients() {
    let material = material();
    let cell = material.macroscopic_cell(91, 0.25, 3.0, || Ok(())).unwrap();
    assert_eq!(material.data_id(), "manufactured mass collision data");
    assert_eq!(material.material_id(), 37);
    assert_eq!(material.radiation(), DecayRadiation::Gamma);
    assert_eq!(material.group_edges_j(), &[0.0, 1.0, 2.0]);
    assert_eq!(cell.cell_id, 91);
    assert_eq!(cell.thickness_m, 0.25);
    assert_eq!(cell.groups[0].total_per_m, 6.0);
    assert_eq!(cell.groups[0].heat_per_m, 1.5);
    assert_eq!(cell.groups[0].redistribution_per_m, [4.5, 0.0]);
    assert_eq!(cell.groups[1].total_per_m, 12.0);
    assert_eq!(cell.groups[1].heat_per_m, 3.0);
    assert_eq!(cell.groups[1].redistribution_per_m, [6.0, 3.0]);
    assert_eq!(material.groups(), rows());
}

#[test]
fn g3_density_depth_reciprocal_scaling_preserves_transport_heat_and_escape() {
    let material = material();
    let solve = |density, depth| {
        let cell = material
            .macroscopic_cell(91, depth, density, || Ok(()))
            .unwrap();
        prepare_slab_response(
            "same physical optical thickness".into(),
            DecayRadiation::Gamma,
            material.group_edges_j().to_vec(),
            vec![cell],
            vec![91],
            vec![PhotonOrdinate {
                mu: 1.0,
                weight: 1.0,
            }],
            PhotonSlabLimits::default(),
            || Ok(()),
        )
        .unwrap()
        .response
        .apply(&[(91, &[2.0, 3.0])], 1e-12, || Ok(()))
        .unwrap()
    };
    let a = solve(3.0, 0.25);
    let b = solve(1.5, 0.5);
    assert_eq!(a, b);
    assert!(a.deposited_j > 0.0 && a.deposited_j < 5.0);
    assert!(a.escaped_j > 0.0 && a.escaped_j < 5.0);
    assert!((a.deposited_j + a.escaped_j - 5.0).abs() < 2e-14);
}

#[test]
fn g4_incomplete_kernels_scaling_loss_limits_and_cancellation_refuse() {
    let admit = |rows, limits| {
        PhotonMassMaterial::new(
            "manufactured refusal data".into(),
            37,
            DecayRadiation::Gamma,
            vec![0.0, 1.0, 2.0],
            rows,
            limits,
            || Ok(()),
        )
    };
    let mut incomplete = rows();
    incomplete[0].heat_m2_kg = 0.25;
    assert!(admit(incomplete, PhotonSlabLimits::default()).is_err());
    let mut upscatter = rows();
    upscatter[0].redistribution_m2_kg = vec![1.0, 0.5];
    assert!(admit(upscatter, PhotonSlabLimits::default()).is_err());
    assert!(matches!(
        admit(
            rows(),
            PhotonSlabLimits {
                max_values: 10,
                ..PhotonSlabLimits::default()
            }
        ),
        Err(RadiationError::Budget(_))
    ));
    let material = material();
    assert_eq!(
        material
            .macroscopic_cell(91, 1.0, f64::MAX, || Ok(()))
            .unwrap_err(),
        RadiationError::Arithmetic
    );
    let tiny = PhotonMassMaterial::new(
        "manufactured subnormal material".into(),
        1,
        DecayRadiation::Gamma,
        vec![0.0, 1.0],
        vec![PhotonMassCollision {
            total_m2_kg: f64::from_bits(1),
            heat_m2_kg: f64::from_bits(1),
            redistribution_m2_kg: vec![0.0],
        }],
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(
        tiny.macroscopic_cell(91, 1.0, 0.5, || Ok(())).unwrap_err(),
        RadiationError::Arithmetic
    );
    let expected = RadiationError::InvalidInput("original material callback failure");
    let mut calls = 0;
    let error = material
        .macroscopic_cell(91, 0.25, 3.0, || {
            calls += 1;
            if calls == 4 {
                Err(expected.clone())
            } else {
                Ok(())
            }
        })
        .unwrap_err();
    assert_eq!(error, expected);
    assert_eq!(material.groups(), rows());
    assert_eq!(
        material
            .macroscopic_cell(91, 0.25, 3.0, || Ok(()))
            .unwrap()
            .groups[1]
            .heat_per_m,
        3.0
    );
}

fn constituent(mass_fraction: f64, material: PhotonMassMaterial) -> PhotonMassConstituent {
    PhotonMassConstituent {
        mass_fraction,
        material,
    }
}

fn mix(constituents: Vec<PhotonMassConstituent>) -> PhotonMassMaterial {
    PhotonMassMaterial::from_mass_fractions(
        "manufactured homogeneous mixture".into(),
        900,
        constituents,
        "Independent manufactured interactions with common spectral weighting; the supplied local heat and outgoing photon kernels apply in this homogeneous target.".into(),
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

fn mono_material(id: u64, heat: f64, scatter: f64) -> PhotonMassMaterial {
    PhotonMassMaterial::new(
        format!("manufactured constituent {id}"),
        id,
        DecayRadiation::Gamma,
        vec![0.0, 1.0],
        vec![PhotonMassCollision {
            total_m2_kg: heat + scatter,
            heat_m2_kg: heat,
            redistribution_m2_kg: vec![scatter],
        }],
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

fn beam_deposition(
    material: &PhotonMassMaterial,
    density: f64,
    thickness: f64,
    mu: f64,
    energy: &[f64],
) -> crate::radiation::RadiationDeposition {
    let cell = material
        .macroscopic_cell(91, thickness, density, || Ok(()))
        .unwrap();
    prepare_slab_response_with_sources(
        "manufactured mixture beam calorimetry".into(),
        material.radiation(),
        material.group_edges_j().to_vec(),
        vec![cell],
        vec![PhotonSlabSource {
            source_id: 99,
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
    .unwrap()
    .response
    .apply(&[(99, energy)], 1e-12, || Ok(()))
    .unwrap()
}

#[test]
fn g1_mixture_preserves_original_mass_sources_and_complete_energy_kernels() {
    let first = material();
    let second = PhotonMassMaterial::new(
        "manufactured second constituent".into(),
        38,
        DecayRadiation::Gamma,
        vec![0.0, 1.0, 2.0],
        vec![
            PhotonMassCollision {
                total_m2_kg: 4.0,
                heat_m2_kg: 1.0,
                redistribution_m2_kg: vec![3.0, 0.0],
            },
            PhotonMassCollision {
                total_m2_kg: 8.0,
                heat_m2_kg: 3.0,
                redistribution_m2_kg: vec![4.0, 1.0],
            },
        ],
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap();
    let mixed = mix(vec![
        constituent(0.25, first.clone()),
        constituent(0.75, second.clone()),
    ]);
    assert_eq!(mixed.material_id(), 900);
    assert_eq!(mixed.groups()[0].total_m2_kg, 3.5);
    assert_eq!(mixed.groups()[0].heat_m2_kg, 0.875);
    assert_eq!(mixed.groups()[0].redistribution_m2_kg, [2.625, 0.0]);
    assert_eq!(mixed.groups()[1].total_m2_kg, 7.0);
    assert_eq!(mixed.groups()[1].heat_m2_kg, 2.5);
    assert_eq!(mixed.groups()[1].redistribution_m2_kg, [3.5, 1.0]);
    let retained = mixed.mixture().unwrap();
    assert!(!retained.kernel_applicability.is_empty());
    assert_eq!(retained.constituents[0].mass_fraction, 0.25);
    assert_eq!(retained.constituents[1].mass_fraction, 0.75);
    assert!(Arc::ptr_eq(&retained.constituents[0].material.0, &first.0));
    assert!(Arc::ptr_eq(&retained.constituents[1].material.0, &second.0));
    assert!(first.mixture().is_none());

    // Both same-group scattering and downscatter reach the existing transport
    // solve; local heat remains distinct from total photon removal.
    let result = beam_deposition(&mixed, 2.0, 0.4, 0.6, &[0.0, 7.0]);
    assert!(result.deposited_j > 0.0 && result.deposited_j < 7.0);
    assert!((result.deposited_j + result.escaped_j - 7.0).abs() < 5e-14);
    assert_eq!(result.receiver_energy_j[&91], result.deposited_j);
}

#[test]
fn g1_mixed_absorbers_follow_mass_thickness_and_oblique_beer_lambert_heat() {
    let mixed = mix(vec![
        constituent(0.75, mono_material(1, 0.4, 0.0)),
        constituent(0.25, mono_material(2, 3.6, 0.0)),
    ]);
    for thickness in [1e-12, 0.25, 3.0] {
        for mu in [1.0, 0.6] {
            let result = beam_deposition(&mixed, 2.0, thickness, mu, &[7.0]);
            // Independent constituent column masses determine the exponent.
            let tau = (0.75 * 0.4 + 0.25 * 3.6) * 2.0 * thickness / mu;
            let heat = 7.0 * (-(-tau).exp_m1());
            let escaped = 7.0 * (-tau).exp();
            assert!((result.deposited_j - heat).abs() < 6e-14 * heat);
            assert!((result.escaped_j - escaped).abs() < 6e-14 * escaped);
            assert!((result.deposited_j + result.escaped_j - 7.0).abs() < 6e-14);
        }
    }
}

#[test]
fn g1_mixture_scattering_reabsorbs_photons_without_treating_extinction_as_heat() {
    let mixed = mix(vec![
        constituent(0.25, mono_material(1, 2.0, 0.0)),
        constituent(0.75, mono_material(2, 0.0, 2.0)),
    ]);
    let result = beam_deposition(&mixed, 1.0, 0.8, 0.6, &[7.0]);
    let tau = 1.6_f64;
    let first_collision = -(-tau / 0.6).exp_m1();
    let uniform_collision = 1.0 + (-tau).exp_m1() / tau;
    // The mixed medium has 1/4 heat and 3/4 elastic-scattering probability
    // per collision. Sum its infinite repeated-collision geometric series.
    let expected = 7.0 * first_collision * 0.25 / (1.0 - 0.75 * uniform_collision);
    assert!((result.deposited_j - expected).abs() < 5e-14);
    assert!((result.escaped_j - (7.0 - expected)).abs() < 5e-14);
    assert!(result.deposited_j < 7.0 * first_collision);
}

#[test]
fn g3_mixture_repartition_and_partial_density_preserve_the_physical_kernel() {
    let a = mono_material(1, 2.0, 0.0);
    let b = mono_material(2, 0.0, 2.0);
    let direct = mix(vec![
        constituent(0.25, a.clone()),
        constituent(0.75, b.clone()),
    ]);
    let half = mix(vec![
        constituent(0.5, a.clone()),
        constituent(0.5, b.clone()),
    ]);
    let regrouped = mix(vec![constituent(0.5, half), constituent(0.5, b.clone())]);
    assert_eq!(direct.groups(), regrouped.groups());
    let cell = direct.macroscopic_cell(91, 0.8, 4.0, || Ok(())).unwrap();
    let partial_a = a.macroscopic_cell(91, 0.8, 1.0, || Ok(())).unwrap();
    let partial_b = b.macroscopic_cell(91, 0.8, 3.0, || Ok(())).unwrap();
    assert_eq!(
        cell.groups[0].total_per_m,
        partial_a.groups[0].total_per_m + partial_b.groups[0].total_per_m
    );
    assert_eq!(
        cell.groups[0].heat_per_m,
        partial_a.groups[0].heat_per_m + partial_b.groups[0].heat_per_m
    );
    assert_eq!(
        cell.groups[0].redistribution_per_m[0],
        partial_a.groups[0].redistribution_per_m[0] + partial_b.groups[0].redistribution_per_m[0]
    );
}

#[test]
fn g4_mixture_rejects_unowned_physics_composition_grid_and_resource_mismatches() {
    let source = mono_material(1, 2.0, 0.0);
    let admit = |parts, applicability: &str, limits| {
        PhotonMassMaterial::from_mass_fractions(
            "manufactured refusal".into(),
            900,
            parts,
            applicability.into(),
            limits,
            || Ok(()),
        )
    };
    for fraction in [-0.1, 0.0, 0.5, 1.1, f64::NAN, f64::INFINITY] {
        assert!(
            admit(
                vec![constituent(fraction, source.clone())],
                "declared basis",
                PhotonSlabLimits::default()
            )
            .is_err()
        );
    }
    assert!(admit(Vec::new(), "declared basis", PhotonSlabLimits::default()).is_err());
    assert!(
        admit(
            vec![constituent(1.0, source.clone())],
            " ",
            PhotonSlabLimits::default()
        )
        .is_err()
    );
    for (radiation, edges) in [
        (DecayRadiation::XRay, vec![0.0, 1.0]),
        (DecayRadiation::Gamma, vec![0.0, 2.0]),
    ] {
        let incompatible = PhotonMassMaterial::new(
            "other basis".into(),
            2,
            radiation,
            edges,
            source.groups().to_vec(),
            PhotonSlabLimits::default(),
            || Ok(()),
        )
        .unwrap();
        assert!(
            admit(
                vec![
                    constituent(0.5, source.clone()),
                    constituent(0.5, incompatible)
                ],
                "declared basis",
                PhotonSlabLimits::default()
            )
            .is_err()
        );
    }
    for limits in [
        PhotonSlabLimits {
            max_values: 1,
            ..PhotonSlabLimits::default()
        },
        PhotonSlabLimits {
            max_work: 1,
            ..PhotonSlabLimits::default()
        },
        PhotonSlabLimits {
            max_groups: 0,
            ..PhotonSlabLimits::default()
        },
    ] {
        assert!(matches!(
            admit(
                vec![constituent(1.0, source.clone())],
                "declared basis",
                limits
            ),
            Err(RadiationError::Budget(_))
        ));
    }
    let expected = RadiationError::InvalidInput("cancelled mixture preparation");
    let mut calls = 0;
    let result = PhotonMassMaterial::from_mass_fractions(
        "cancelled mixture".into(),
        900,
        vec![constituent(1.0, source.clone())],
        "declared basis".into(),
        PhotonSlabLimits::default(),
        || {
            calls += 1;
            if calls == 5 {
                Err(expected.clone())
            } else {
                Ok(())
            }
        },
    );
    assert_eq!(result.unwrap_err(), expected);
    assert_eq!(source.groups()[0].heat_m2_kg, 2.0);
}

#[test]
fn g4_mixture_retains_representable_extremes_and_refuses_lost_positive_terms() {
    for value in [f64::from_bits(1), 1.0, f64::MAX] {
        let source = mono_material(1, value, 0.0);
        assert_eq!(
            mix(vec![constituent(1.0, source)]).groups()[0].total_m2_kg,
            value
        );
    }
    let source = mono_material(1, f64::MAX, 0.0);
    assert_eq!(
        mix(vec![
            constituent(0.5, source.clone()),
            constituent(0.5, source)
        ])
        .groups()[0]
            .total_m2_kg,
        f64::MAX
    );
    let tiny = mono_material(1, f64::from_bits(1), 0.0);
    let normal = mono_material(2, 1.0, 0.0);
    let refused = PhotonMassMaterial::from_mass_fractions(
        "positive term underflow".into(),
        900,
        vec![
            constituent(0.5, tiny.clone()),
            constituent(0.5, normal.clone()),
        ],
        "declared basis".into(),
        PhotonSlabLimits::default(),
        || Ok(()),
    );
    assert_eq!(refused.unwrap_err(), RadiationError::Arithmetic);
    assert_eq!(
        mix(vec![constituent(0.0, tiny), constituent(1.0, normal)]).groups()[0].total_m2_kg,
        1.0
    );
}
