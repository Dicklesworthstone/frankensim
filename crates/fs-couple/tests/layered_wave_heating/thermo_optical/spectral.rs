//! Finite-band source feedback through the existing coherent/thermal owners.

use super::*;
use fs_couple::layered_wave_heating::{SpectralThermoOpticalSource, ThermoOpticalSource};
use fs_em::plane_wave::SpectralOpticalResponse;
use fs_material::state_point::{
    SPECTRAL_CONDUCTOR_ETA_PROPERTIES, SPECTRAL_CONDUCTOR_K_PROPERTIES,
    SPECTRAL_CONDUCTOR_K_QUANTITY, SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
    SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY, SPECTRAL_DIELECTRIC_IOR_PROPERTIES,
    SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES, SpectralConductorStatePoint,
    SpectralOpticalStatePoint, SpectralTabulatedDielectricStatePoint,
    resolve_spectral_conductor_state_point, resolve_spectral_tabulated_dielectric_state_point,
};

const IR_OMEGA: f64 = core::f64::consts::TAU * LIGHT_SPEED / 11e-6;
const MAX_GAP_M: f64 = 0.751e-6;

fn infrared_card(dielectric: bool, upper: f64, widen_gap: bool) -> MaterialCard {
    let base = card(
        if dielectric {
            "infrared glass"
        } else {
            "infrared metal"
        },
        1.0,
        0.0,
    );
    let mut claims = base.claims().clone();
    for index in 0..9 {
        let wavelength = (8.0 + 0.75 * index as f64) * 1e-6;
        let wavelengths = if widen_gap && index == 8 {
            [wavelength - 0.5e-6, wavelength + 1.5e-6]
        } else {
            [wavelength; 2]
        };
        let (n_name, n, loss_name, loss_quantity, loss) = if dielectric {
            (
                SPECTRAL_DIELECTRIC_IOR_PROPERTIES[index],
                [1.4, 1.56],
                SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[index],
                SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
                [1e4, 4e4],
            )
        } else {
            (
                SPECTRAL_CONDUCTOR_ETA_PROPERTIES[index],
                [0.4, 1.6],
                SPECTRAL_CONDUCTOR_K_PROPERTIES[index],
                SPECTRAL_CONDUCTOR_K_QUANTITY,
                [2.8, 1.2],
            )
        };
        for (name, quantity, values) in [
            (
                SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[index],
                QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                wavelengths,
            ),
            (n_name, QuantitySpec::dimensional(Dims::NONE), n),
            (loss_name, loss_quantity, loss),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::with_quantity(name, quantity),
                    value: PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: fs_qty::Temperature::DIMS,
                        knots: vec![(280.0, values[0]), (360.0, values[1])],
                        dims: quantity.dims(),
                    },
                    validity: ValidityDomain::unconstrained()
                        .with_quantity(
                            "T",
                            super::quantity(QuantityKind::AbsoluteTemperature),
                            280.0,
                            upper,
                        )
                        .with_quantity(
                            "pressure",
                            super::quantity(QuantityKind::Pressure),
                            90_000.0,
                            200_000.0,
                        ),
                    interpolation: InterpolationPolicy::LinearInside,
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: "authored finite-band thermo-optical regression".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

fn infrared_source(
    card: &MaterialCard,
    dielectric: bool,
    include_mu: bool,
) -> SpectralOpticalStatePoint {
    let source = if dielectric {
        SpectralOpticalStatePoint::TabulatedDielectric(
            resolve_spectral_tabulated_dielectric_state_point(
                card,
                &at(300.0, 100_000.0),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap(),
        )
    } else {
        SpectralOpticalStatePoint::Conductor(
            resolve_spectral_conductor_state_point(
                card,
                &at(300.0, 100_000.0),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap(),
        )
    };
    if !include_mu {
        return source;
    }
    let complete = source
        .resolved()
        .with_additional_requirements(
            card,
            &[ScalarPropertyRequirement::try_new(
                MU_NAME,
                MU,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
    match source {
        SpectralOpticalStatePoint::Conductor(_) => SpectralOpticalStatePoint::Conductor(
            SpectralConductorStatePoint::try_from_material(complete).unwrap(),
        ),
        SpectralOpticalStatePoint::TabulatedDielectric(_) => {
            SpectralOpticalStatePoint::TabulatedDielectric(
                SpectralTabulatedDielectricStatePoint::try_from_material(complete).unwrap(),
            )
        }
    }
}

#[test]
fn infrared_mixed_sources_refresh_each_native_cell_and_feed_conservative_heat_and_next_step() {
    let glass_card = infrared_card(true, 360.0, false);
    let metal_card = infrared_card(false, 360.0, false);
    let glass =
        SpectralThermoOpticalSource::try_new(infrared_source(&glass_card, true, true), MAX_GAP_M)
            .unwrap();
    let metal =
        SpectralThermoOpticalSource::try_new(infrared_source(&metal_card, false, true), MAX_GAP_M)
            .unwrap();
    let outside = exterior();
    let native = body(&[&glass_card, &metal_card], &[4e-6, 0.3e-6], &[2, 1], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let initial = initial_at(&bound, &[290.0, 325.0, 305.0]);
    let before = initial.cells().clone();
    let original_cells = [
        (&glass, &glass_card),
        (&glass, &glass_card),
        (&metal, &metal_card),
    ];
    let initial_temperatures: Vec<_> = bound
        .cell_ids()
        .into_iter()
        .map(|id| initial.temperature_k(id).unwrap())
        .collect();
    let refreshed: Vec<_> = initial_temperatures
        .iter()
        .zip(original_cells)
        .map(|(&t, (source, card))| source.at_query_point(card, &at(t, 100_000.0)).unwrap())
        .collect();
    let responses: Vec<_> = refreshed
        .iter()
        .map(|source| {
            DielectricResponse::SpectralOptical(
                SpectralOpticalResponse::try_new(source.source().into(), MAX_GAP_M).unwrap(),
            )
        })
        .collect();
    let reference = LayeredWaveHeating::try_new_cellwise(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        IR_OMEGA,
        &responses,
        || false,
    )
    .unwrap();
    let drive = LayeredWaveDrive::Coherent([C64::new(1e6, 0.0), C64::new(2e5, -1e5)]);
    let mut powers = [0.0; 3];
    reference
        .prepare(0.9, Polarization::Tm, || false)
        .unwrap()
        .cell_power_into(drive, &mut powers, || false)
        .unwrap();
    assert!(powers.iter().all(|p| *p > 0.0));
    let owner = ThermoOpticalWaveHeating::try_new(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        IR_OMEGA,
        0.9,
        Polarization::Tm,
        &[
            ThermoOpticalLayer {
                source: &glass,
                card: &glass_card,
            },
            ThermoOpticalLayer {
                source: &metal,
                card: &metal_card,
            },
        ],
        || false,
    )
    .unwrap();
    let mut workspace = LayeredThermalWorkspace::new(3).unwrap();
    let controls = LayeredThermalSolveControls {
        cell_tolerance_j: 1e-11,
        total_tolerance_j: 1e-11,
        ..LayeredThermalSolveControls::default()
    };
    let dt = 2e-10;
    let result = owner
        .advance_to(
            &initial,
            dt,
            drive,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(result.applied_optical_states, refreshed);
    let temperatures: Vec<_> = bound
        .cell_ids()
        .into_iter()
        .map(|id| result.step.thermal.temperature_k(id).unwrap())
        .collect();
    let widths = [2e-6, 2e-6, 0.3e-6];
    for (i, id) in bound.cell_ids().into_iter().enumerate() {
        // Independent local constant-capacity backward-Euler heat balance,
        // including both conductive neighbors at the solved endpoint.
        let storage = RHO * AREA * widths[i] * CP * (temperatures[i] - initial_temperatures[i]);
        let mut expected_heat = powers[i] * dt;
        for neighbor in [i.checked_sub(1), (i + 1 < 3).then_some(i + 1)]
            .into_iter()
            .flatten()
        {
            expected_heat += AREA * K * dt / (0.5 * (widths[i] + widths[neighbor]))
                * (temperatures[neighbor] - temperatures[i]);
        }
        near(storage, expected_heat, 3e-10);
        near(result.step.cell_energy_j[i].1, powers[i] * dt, 1e-11);
        assert_eq!(result.step.cell_energy_j[i].0, id);
        let (original, card) = original_cells[i];
        assert_eq!(
            result.endpoint_optical_states[i],
            original
                .at_query_point(card, &at(temperatures[i], 100_000.0))
                .unwrap()
        );
        assert_eq!(
            result.endpoint_optical_states[i].maximum_knot_span_m(),
            MAX_GAP_M
        );
    }
    near(
        result.step.heat_report.energy.stored_energy_change_j,
        result.step.optical.absorbed_j,
        3e-10,
    );
    assert!(matches!(
        result.endpoint_optical_states[0].source(),
        SpectralOpticalStatePoint::TabulatedDielectric(_)
    ));
    assert!(matches!(
        result.endpoint_optical_states[2].source(),
        SpectralOpticalStatePoint::Conductor(_)
    ));
    let next = owner
        .advance_to(
            &result.step.thermal,
            2.0 * dt,
            drive,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(next.applied_optical_states, result.endpoint_optical_states);
    assert_ne!(next.applied_optical_states, result.applied_optical_states);
    assert_eq!(initial.cells(), &before);
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn infrared_endpoint_source_exit_refuses_after_native_heat_without_publishing_state() {
    let card = infrared_card(true, 300.01, false);
    let source =
        SpectralThermoOpticalSource::try_new(infrared_source(&card, true, true), MAX_GAP_M)
            .unwrap();
    let outside = exterior();
    let native = body(&[&card], &[2e-6], &[1], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let initial = initial_at(&bound, &[300.0]);
    let before = initial.cells().clone();
    let owner = ThermoOpticalWaveHeating::try_new(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        IR_OMEGA,
        1.0,
        Polarization::Te,
        &[ThermoOpticalLayer {
            source: &source,
            card: &card,
        }],
        || false,
    )
    .unwrap();
    let responses = [source.wave_response(IR_OMEGA).unwrap()];
    let frozen = LayeredWaveHeating::try_new(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        IR_OMEGA,
        &responses,
        || false,
    )
    .unwrap();
    let drive = LayeredWaveDrive::Coherent([C64::new(1.5e6, 0.0), C64::ZERO]);
    let dt = 1e-9;
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    let controls = LayeredThermalSolveControls::default();
    let candidate = frozen
        .prepare(1.0, Polarization::Te, || false)
        .unwrap()
        .advance_to(
            &initial,
            dt,
            drive,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    let t = candidate
        .thermal
        .temperature_k(bound.cell_ids()[0])
        .unwrap();
    assert!(t > 300.01 && t < 360.0);
    assert!(matches!(
        owner.advance_to(
            &initial,
            dt,
            drive,
            conduction(),
            controls,
            &mut workspace,
            || false
        ),
        Err(LayeredWaveHeatingError::Material(_))
    ));
    assert_eq!(initial.cells(), &before);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(source.resolved().state_coordinate("T"), Some(300.0));
}

#[test]
fn infrared_constructor_and_refreshed_response_retain_gap_mu_and_frequency_admission() {
    let card = infrared_card(false, 360.0, true);
    let original = infrared_source(&card, false, true);
    for maximum in [0.0, f64::NAN, f64::INFINITY, 0.5e-6] {
        assert!(SpectralThermoOpticalSource::try_new(original.clone(), maximum).is_err());
    }
    assert!(
        SpectralThermoOpticalSource::try_new(infrared_source(&card, false, false), MAX_GAP_M)
            .is_err()
    );
    let source = SpectralThermoOpticalSource::try_new(original, MAX_GAP_M).unwrap();
    assert!(source.wave_response(IR_OMEGA).is_ok());
    assert!(source.wave_response(OPTICAL_OMEGA).is_err());
    let hot = source.at_query_point(&card, &at(330.0, 100_000.0)).unwrap();
    assert_eq!(hot.maximum_knot_span_m(), MAX_GAP_M);
    assert!(
        hot.wave_response(IR_OMEGA).is_err(),
        "wavelength refresh must not loosen the original gap policy"
    );
    assert_eq!(source.resolved().state_coordinate("T"), Some(300.0));
}
