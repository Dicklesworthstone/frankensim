//! G1/G3/G4: physical density conversion, optical-depth scaling and refusal.

use super::*;
use crate::radiation::transport::{PhotonOrdinate, prepare_slab_response};

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
