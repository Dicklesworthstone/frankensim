use super::*;
use crate::thermal_surface::heating::{
    CoatingSurfaceSource, StatefulMaterialSlabIlluminationSource, SurfaceHeatingBudget,
    SurfaceHeatingTolerance, SurfaceIncidence, SurfaceOpticalHeating,
};
use crate::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};

const WAVELENGTHS: [f64; 2] = [9500.0, 11_000.0];
const INCIDENT: f64 = 1.25;
const EXIT: f64 = 1.05;
const COSINE: f64 = 0.6;
const PROPERTIES: usize = 29;

fn segment(start: f64, end: f64, first: [f64; 2], last: [f64; 2]) -> SpectralPowerSegment<2> {
    SpectralPowerSegment {
        start_s: start,
        end_s: end,
        start_power_w: first,
        end_power_w: last,
    }
}

fn illumination(
    id: u64,
    optics: Arc<MaterialSlabLaw>,
    segments: &[SpectralPowerSegment<2>],
) -> StatefulMaterialSlabIlluminationSource<2> {
    StatefulMaterialSlabIlluminationSource {
        source_id: id,
        receiver: receiver(),
        optics,
        thickness_m: DEPTH,
        incidence: SurfaceIncidence {
            wavelengths_nm: WAVELENGTHS,
            incident_eta: [INCIDENT; 2],
            incidence_cosine: COSINE,
            // Separate pure s and p lanes expose mixing and interface errors.
            s_power_fraction: [1.0, 0.0],
        },
        transmitted_eta: [EXIT; 2],
        source_identity: ContentHash([215; 32]),
        power: SpectralPowerSchedule::try_new(segments, 2).unwrap(),
    }
}

fn observer(
    thermal: &DepositionSimulation<'_>,
    sources: &[StatefulMaterialSlabIlluminationSource<2>],
    properties: usize,
) -> Result<SurfaceOpticalHeating<2>, SurfaceHeatingError> {
    let inputs: Vec<_> = sources
        .iter()
        .map(CoatingSurfaceSource::MaterialIncoherentSlab)
        .collect();
    SurfaceOpticalHeating::try_new_coating_stacks(
        &inputs,
        thermal,
        SurfaceHeatingBudget {
            max_sources: 2,
            max_total_power_segments: 4,
            max_source_evaluations: 128,
        },
        properties,
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e-10,
        },
        || Ok(()),
    )
}

fn native_with_background(temperature: f64, watts: f64) -> DepositionSimulation<'static> {
    let initial = native(temperature);
    let mut cell = initial.cells()[&CELL].clone();
    cell.boundary.source_w = watts;
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        }],
        &[(MATERIAL, &SENSIBLE)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0
}

fn slab_power(wavelength: f64, temperature: f64, perpendicular: bool) -> [f64; 3] {
    // Independent asymmetric two-interface geometric series. Normal flux
    // transmission is 1-R at each lossless interface, with no radiance eta^2.
    let n = schott::phase(wavelength, temperature);
    let transverse2 = INCIDENT * INCIDENT * (1.0 - COSINE * COSINE);
    let ct = (1.0 - transverse2 / (n * n)).sqrt();
    let ce = (1.0 - transverse2 / (EXIT * EXIT)).sqrt();
    let (r0, r1) = if perpendicular {
        (
            ((INCIDENT * COSINE - n * ct) / (INCIDENT * COSINE + n * ct)).powi(2),
            ((n * ct - EXIT * ce) / (n * ct + EXIT * ce)).powi(2),
        )
    } else {
        (
            ((n * COSINE - INCIDENT * ct) / (n * COSINE + INCIDENT * ct)).powi(2),
            ((EXIT * ct - n * ce) / (EXIT * ct + n * ce)).powi(2),
        )
    };
    let a = (-absorption(wavelength, 300.0) * DEPTH / ct).exp();
    let denominator = 1.0 - r0 * r1 * a * a;
    let reflected = r0 + (1.0 - r0).powi(2) * r1 * a * a / denominator;
    let transmitted = (1.0 - r0) * (1.0 - r1) * a / denominator;
    let absorbed = (1.0 - r0) * (1.0 - a) * (1.0 + r1 * a) / denominator;
    near(reflected + transmitted + absorbed, 1.0);
    [reflected, transmitted, absorbed]
}

#[test]
fn oblique_asymmetric_illumination_feeds_back_refraction_without_heating_transmitted_power() {
    let thermal = native(400.0);
    let power = [40.0, 60.0];
    let source = illumination(
        1,
        schott::schott_law(600.0, 600.0),
        &[segment(0.0, 0.2, power, power)],
    );
    let input = observer(&thermal, &[source], PROPERTIES).unwrap();
    let mut expected_temperature = 400.0;
    let mut expected = [[0.0; 2]; 3];
    for _ in 0..2 {
        let mut heat = 0.0;
        for lane in 0..2 {
            let fractions = slab_power(WAVELENGTHS[lane], expected_temperature, lane == 0);
            for channel in 0..3 {
                expected[channel][lane] += 0.1 * power[lane] * fractions[channel];
            }
            heat += 0.1 * power[lane] * fractions[2];
        }
        expected_temperature += heat;
    }
    let (advanced, observed, audit) = thermal
        .advance_observed(&input, 0.2, config(), || Ok(()))
        .unwrap();
    near(advanced.temperature_k(CELL).unwrap(), expected_temperature);
    for lane in 0..2 {
        near(observed.energy().incident_j[lane], 0.2 * power[lane]);
        near(observed.energy().reflected_j[lane], expected[0][lane]);
        near(observed.energy().transmitted_j[lane], expected[1][lane]);
        near(observed.energy().absorbed_j[lane], expected[2][lane]);
        assert!(observed.energy().transmitted_j[lane] > 0.0);
    }
    let heat: f64 = expected[2].iter().sum();
    near(observed.energy().thermal_input_j, heat);
    near(audit.external_heat_j, heat);
    near(advanced.energy_j() - thermal.energy_j(), heat);
    assert_eq!(audit.substeps, 2);
    assert_eq!(observed.response_evaluations(), 1);
    let frozen_temperature = 400.0
        + (0..2)
            .map(|lane| 0.2 * power[lane] * slab_power(WAVELENGTHS[lane], 400.0, lane == 0)[2])
            .sum::<f64>();
    assert!((expected_temperature - frozen_temperature).abs() > 1e-4);
    assert_eq!(thermal.temperature_k(CELL).unwrap(), 400.0);
    assert_eq!(input.energy().thermal_input_j, 0.0);
}

#[test]
fn shared_law_pulses_keep_exact_fluence_source_order_and_dark_packet_reuse() {
    let law = schott::schott_law(600.0, 600.0);
    let pulse = [segment(0.05, 0.15, [0.0, 40.0], [40.0, 0.0])];
    let first = illumination(10, Arc::clone(&law), &pulse);
    let second = illumination(20, law, &pulse);
    let thermal = native(400.0);
    let forward = observer(&thermal, &[first.clone(), second.clone()], PROPERTIES).unwrap();
    let backward = observer(&thermal, &[second, first], PROPERTIES).unwrap();
    assert_eq!(forward.model_identity(), backward.model_identity());
    let step = StepConfig {
        max_dt_s: 0.04,
        ..config()
    };
    let (dark, dark_history, _) = thermal
        .advance_observed(&forward, 0.05, step, || Ok(()))
        .unwrap();
    assert_eq!(dark.temperature_k(CELL).unwrap(), 400.0);
    assert_eq!(dark_history.energy().incident_j, [0.0; 2]);
    assert_eq!(dark_history.response_evaluations(), 0);
    let (heated, observed, _) = thermal
        .advance_observed(&forward, 0.2, step, || Ok(()))
        .unwrap();
    let (reordered, reordered_history, _) = thermal
        .advance_observed(&backward, 0.2, step, || Ok(()))
        .unwrap();
    assert_eq!(heated.temperature_k(CELL), reordered.temperature_k(CELL));
    assert_eq!(observed.energy(), reordered_history.energy());
    for incident in observed.energy().incident_j {
        near(incident, 4.0);
    }
    // Three illuminated native intervals need two updated optical packets;
    // both sources share each packet, and the trailing darkness needs none.
    assert_eq!(observed.response_evaluations(), 2);
    assert!(observed.energy().thermal_input_j > 0.0);
    assert!(observed.energy().thermal_input_j < 8.0);
    near(
        heated.energy_j() - thermal.energy_j(),
        observed.energy().thermal_input_j,
    );
}

#[test]
fn unsupported_hot_and_dark_states_cancel_without_publishing_heat_or_packets() {
    let thermal = native(400.0);
    let limited = illumination(
        1,
        schott::schott_law(600.0, 400.0),
        &[segment(0.0, 0.2, [40.0, 60.0], [40.0, 60.0])],
    );
    let input = observer(&thermal, &[limited.clone()], PROPERTIES).unwrap();
    assert!(
        thermal
            .advance_observed(&input, 0.1, config(), || Ok(()))
            .is_err()
    );
    assert!(matches!(
        thermal.advance_observed(&input, 0.1, config(), || Err(DepositionError::Cancelled)),
        Err(SurfaceHeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert!(matches!(
        observer(&thermal, &[limited], PROPERTIES - 1),
        Err(SurfaceHeatingError::Budget)
    ));
    assert_eq!(thermal.temperature_k(CELL).unwrap(), 400.0);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(input.time_s(), 0.0);
    assert_eq!(input.energy().thermal_input_j, 0.0);
    assert_eq!(input.response_evaluations(), 0);

    // Current-claim mode may move optical support as mechanics heats, even
    // with zero illumination. Every accepted endpoint must still be admitted.
    let (card, original) = fixture(|claims| {
        for claim in claims {
            if SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES.contains(&claim.key.name()) {
                if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                    knots[1].1 += 6e-6;
                    knots[2].1 += 12e-6;
                }
            }
        }
    });
    let mut dark_source = illumination(2, law(card, original), &[]);
    dark_source.incidence.wavelengths_nm[0] = 8001.0;
    let externally_heated = native_with_background(300.0, 100.0);
    let dark = observer(&externally_heated, &[dark_source], PROPERTIES).unwrap();
    assert!(
        externally_heated
            .advance_observed(&dark, 0.1, config(), || Ok(()))
            .is_err()
    );
    assert_eq!(externally_heated.temperature_k(CELL).unwrap(), 300.0);
    assert_eq!(externally_heated.time_s(), 0.0);
    assert_eq!(dark.time_s(), 0.0);
    assert_eq!(dark.energy().incident_j, [0.0; 2]);
    assert_eq!(dark.energy().thermal_input_j, 0.0);
    assert_eq!(dark.response_evaluations(), 0);
}
