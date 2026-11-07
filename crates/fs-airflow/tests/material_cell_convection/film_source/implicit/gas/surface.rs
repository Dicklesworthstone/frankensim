//! Manufactured original T/p sources consumed by the native exterior owner.
use super::*;
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};

pub(super) fn law(pressure: f64, corner: bool) -> IsobaricGasLaw {
    let calorics =
        thermochemical_source::model_transport_with_heat_capacity_slope(1.0, 1.0 / 512.0);
    let mut claims = ClaimSet::new();
    for (name, dims, base, thermal, pressure_slope) in [
        (
            DYNAMIC_VISCOSITY_PROPERTY,
            DYNAMIC_VISCOSITY_DIMS,
            2.0_f64.powi(-16),
            2.0_f64.powi(-25),
            2.0_f64.powi(-36),
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            2.0_f64.powi(-6),
            2.0_f64.powi(-15),
            2.0_f64.powi(-24),
        ),
    ] {
        let mut samples: Vec<_> = [
            [280.0, 65536.0],
            [320.0, 65536.0],
            [280.0, 131072.0],
            [320.0, 131072.0],
        ]
        .into_iter()
        .map(|[t, p]| {
            (
                [t, p],
                base + thermal * (t - 280.0) + pressure_slope * (p - 65536.0),
            )
        })
        .collect();
        if corner && name == DYNAMIC_VISCOSITY_PROPERTY {
            samples[3].1 += 2.0_f64.powi(-18);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Surface2D {
                    axes: ["temperature".into(), "pressure".into()],
                    axis_dims: [Temperature::DIMS, Pressure::DIMS],
                    samples,
                    triangles: vec![[0, 1, 2], [1, 3, 2]],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        thermochemical_source::typed(QuantityKind::AbsoluteTemperature),
                        280.0,
                        320.0,
                    )
                    .with_quantity(
                        "pressure",
                        thermochemical_source::typed(QuantityKind::Pressure),
                        65536.0,
                        131072.0,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured T/p transport; not measured".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let source = SpeciesTransportMaterial {
        species: calorics.transport_materials()[0].species.clone(),
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: "A".into(),
                phase: "gas".into(),
                process: "manufactured surface".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    let model =
        ThermochemicalGasModel::new(calorics.thermodynamics().clone(), vec![source], None).unwrap();
    IsobaricGasLaw::try_from_model(model, thermochemical_source::point(300.0, pressure)).unwrap()
}

#[test]
fn surface_transport_at_pressure_edges_drives_original_film_coefficient_and_tangent() {
    let body = curve();
    for pressure in [65536.0, 98304.0, 131072.0] {
        let gas = law(pressure, false);
        let saved = gas.reference().clone();
        // Construction requests derivatives at exactly the original 280/320 K
        // source endpoints, including both corners of each fixed-pressure edge.
        let source =
            ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &gas, settings(plate()))
                .unwrap();
        for wall in [290.0, 300.0, 310.0] {
            let film = 0.5 * (wall + 300.0);
            let original = gas.at_temperature(Temperature::new(film)).unwrap();
            let expected = evaluate_flat_plate_section_from_material(
                CorrelationId::FlatPlateLaminarAverage,
                &original,
                1.0,
                Length::new(0.025),
                Length::new(0.1),
            )
            .unwrap();
            let actual = source.evaluate(wall, 300.0).unwrap();
            near(
                actual.coefficient.value.value(),
                expected.coefficient().value.value(),
                1e-11,
            );
            assert_eq!(actual.transport, original.transport_state().unwrap());
            let state = original.state();
            // h ~ rho^(1/2) mu^(-1/6) Cp^(1/3) k^(2/3); all pure
            // source slopes are authored above. The extra 1/2 is dTfilm/dTs.
            let logarithmic = -0.5 / film - 2.0_f64.powi(-25) / (6.0 * state.dynamic_viscosity)
                + (1.0 / 512.0) / (3.0 * (3.5 + film / 512.0))
                + 2.0 * 2.0_f64.powi(-15) / (3.0 * state.thermal_conductivity);
            near(
                actual.coefficient_surface_derivative_w_m2_k2,
                0.5 * actual.coefficient.value.value() * logarithmic,
                1e-12,
            );
            derivative(&source, wall, 300.0);
        }
        assert_eq!(gas.reference(), &saved);
    }
}

#[test]
fn active_film_queries_refuse_an_original_surface_corner_without_changing_the_gas() {
    let body = curve();
    let gas = law(98304.0, true);
    let saved = gas.reference().clone();
    let source =
        ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &gas, settings(plate()))
            .unwrap();
    assert!(source.evaluate(299.5, 300.0).is_ok());
    assert!(source.evaluate(300.5, 300.0).is_ok());
    assert!(source.evaluate(300.0, 300.0).is_err());
    assert_eq!(gas.reference(), &saved);
}
