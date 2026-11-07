//! G1/G3/G4: original interaction kernels, interpolation and source support.

use super::*;

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
