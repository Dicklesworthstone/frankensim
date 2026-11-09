use super::*;

fn admitted(reynolds: f64, schmidt: f64) -> Result<ForcedSphereFilm, ForcedSphereFilmError> {
    evaluate(
        ForcedSphereMassTransfer::Frossling1938,
        ForcedSphereGeometry {
            diameter_m: 1.0,
            relative_speed_m_s: reynolds * schmidt,
        },
        schmidt,
        1.0,
    )
}

#[test]
fn original_wind_factor_sets_diameter_sherwood_and_dimensional_conductance() {
    let result = evaluate(
        ForcedSphereMassTransfer::Frossling1938,
        ForcedSphereGeometry {
            diameter_m: 0.002,
            relative_speed_m_s: 0.5,
        },
        1e-5,
        1e-5,
    )
    .unwrap();
    // Re=100, Sc=1: original wind factor 1+.276*10=3.76.
    assert!((result.reynolds() - 100.0).abs() < 2e-14);
    assert_eq!(result.schmidt(), 1.0);
    assert!((result.sherwood_number() - 7.52).abs() < 2e-14);
    assert!((result.mass_transfer_coefficient_m_s() - 0.0376).abs() < 2e-16);
    assert_eq!(result.diffusivity_m2_s(), 1e-5);
    assert_eq!(result.kinematic_viscosity_m2_s(), 1e-5);
    assert_eq!(result.model(), ForcedSphereMassTransfer::Frossling1938);
    assert_eq!(result.geometry().diameter_m, 0.002);

    // At fixed groups the transport coefficient scales as D/d, not as area.
    let scaled = evaluate(
        ForcedSphereMassTransfer::Frossling1938,
        ForcedSphereGeometry {
            diameter_m: 0.004,
            relative_speed_m_s: 0.25,
        },
        1e-5,
        1e-5,
    )
    .unwrap();
    assert_eq!(scaled.sherwood_number(), result.sherwood_number());
    assert_eq!(
        scaled.mass_transfer_coefficient_m_s(),
        0.5 * result.mass_transfer_coefficient_m_s(),
    );
}

#[test]
fn admitted_original_range_refuses_extrapolation_and_never_switches_to_stagnant() {
    for reynolds in [2.0, 100.0, 800.0] {
        assert!(admitted(reynolds, 1.0).is_ok());
    }
    for schmidt in [0.6, 1.0, 2.6] {
        assert!(admitted(100.0, schmidt).is_ok());
    }
    for (reynolds, schmidt) in [(1.99, 1.0), (801.0, 1.0), (100.0, 0.59), (100.0, 2.61)] {
        assert!(matches!(
            admitted(reynolds, schmidt),
            Err(ForcedSphereFilmError::OutsideCorrelation { .. }),
        ));
    }
    assert!(admitted(0.0, 1.0).is_err());
    assert!(admitted(-1.0, 1.0).is_err());
}

#[test]
fn nonfinite_geometry_and_transport_do_not_generate_a_film() {
    let geometry = ForcedSphereGeometry { diameter_m: 0.001, relative_speed_m_s: 1.0 };
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(evaluate(
            ForcedSphereMassTransfer::Frossling1938,
            ForcedSphereGeometry { diameter_m: invalid, ..geometry },
            1e-5,
            1e-5,
        ).is_err());
        assert!(evaluate(
            ForcedSphereMassTransfer::Frossling1938,
            ForcedSphereGeometry { relative_speed_m_s: invalid, ..geometry },
            1e-5,
            1e-5,
        ).is_err());
        assert!(evaluate(
            ForcedSphereMassTransfer::Frossling1938, geometry, invalid, 1e-5,
        ).is_err());
        assert!(evaluate(
            ForcedSphereMassTransfer::Frossling1938, geometry, 1e-5, invalid,
        ).is_err());
    }
}
