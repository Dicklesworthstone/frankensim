//! G1/G3: independently known film heat enters the native thermal cells once.

use super::*;
use fs_blake3::ContentHash;
use fs_render::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, CoherentCoatingStack,
    ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource,
};
use fs_render::thermal_surface::heating::{
    CoatedSurfaceIlluminationSource, CoatedSurfacePowerResponse, SurfaceHeatingBudget,
    SurfaceHeatingError, SurfaceHeatingReceiver, SurfaceHeatingTolerance, SurfaceOpticalHeating,
};
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};

const WAVELENGTHS: [f64; 2] = [480.0, 680.0];
const ETA: f64 = 1.5;
const K: f64 = 0.3;

fn optical_source() -> ConductorSource {
    ConductorSource::try_new(ContentHash([207; 32]), ConductorDataStatus::Representative).unwrap()
}
fn receiver(cell: u64) -> SurfaceHeatingReceiver {
    SurfaceHeatingReceiver {
        thermal_cell: CellId(cell),
        material_id: MaterialId(7),
        temperature_bounds_k: (300.0, 600.0),
    }
}
fn source(id: u64, cells: [u64; 3]) -> CoatedSurfaceIlluminationSource<2> {
    let substrate = ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + i as f64 * 50.0, ETA, K).unwrap()
        }),
        optical_source(),
    )
    .unwrap();
    let film = |depth_nm: f64| {
        CoherentAbsorbingCoating::try_new(
            core::array::from_fn(|i| {
                CoatingIorSample::try_new(380.0 + i as f64 * 50.0, ETA, K).unwrap()
            }),
            depth_nm * 1e-9,
            optical_source(),
        )
        .unwrap()
    };
    let first = film(80.0);
    let second = film(120.0);
    let layers = [
        CoatingLayer::Absorbing(&first),
        CoatingLayer::Absorbing(&second),
    ];
    let response = CoatedSurfacePowerResponse::from_coated_conductor(
        WAVELENGTHS,
        &substrate,
        CoherentCoatingStack::try_new(&layers).unwrap(),
        [1.0; 2],
        1.0,
        [0.0, 1.0],
    )
    .unwrap();
    CoatedSurfaceIlluminationSource {
        source_id: id,
        film_receivers: vec![receiver(cells[0]), receiver(cells[1])],
        substrate_receiver: receiver(cells[2]),
        response,
        power: SpectralPowerSchedule::try_new(
            &[SpectralPowerSegment {
                start_s: 0.25,
                end_s: 1.25,
                start_power_w: [8.0, 12.0],
                end_power_w: [16.0, 4.0],
            }],
            1,
        )
        .unwrap(),
    }
}
fn budget() -> SurfaceHeatingBudget {
    SurfaceHeatingBudget {
        max_sources: 4,
        max_total_power_segments: 4,
        max_source_evaluations: 1000,
    }
}
fn heating(
    thermal: &DepositionSimulation<'_>,
    sources: &[CoatedSurfaceIlluminationSource<2>],
) -> SurfaceOpticalHeating<2> {
    SurfaceOpticalHeating::try_new_coated(
        sources,
        thermal,
        budget(),
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn matched_films_heat_distinct_native_cells_by_beer_lambert_not_uniform_partition() {
    let original = thermal(
        vec![
            cell(11, 1.0, 1.0, 0.0),
            cell(12, 2.0, 1.0, 0.0),
            cell(13, 3.0, 1.0, 3.0),
        ],
        &SENSIBLE,
        vec![],
    );
    let input = heating(&original, &[source(1, [11, 12, 13])]);
    let (next, accepted, audit) = original
        .advance_observed(&input, 2.0, config(0.1), || Ok(()))
        .unwrap();
    let mut expected = [0.0; 3];
    let absorbed = 4.0 * ETA / ((1.0 + ETA).powi(2) + K * K);
    for (lane, (wavelength, incident)) in WAVELENGTHS.into_iter().zip([12.0, 8.0]).enumerate() {
        // All interior indices coincide: there are no internal reflections.
        // Entrance Fresnel transmission followed by Beer-Lambert is independent
        // of the production coherent layer integrals and thermal stepper.
        let tau1 = 4.0 * core::f64::consts::PI * K * 80.0 / wavelength;
        let tau2 = 4.0 * core::f64::consts::PI * K * 120.0 / wavelength;
        expected[0] += incident * absorbed * -(-tau1).exp_m1();
        expected[1] += incident * absorbed * (-tau1).exp() * -(-tau2).exp_m1();
        expected[2] += incident * absorbed * (-(tau1 + tau2)).exp();
        close(accepted.energy().incident_j[lane], incident);
        close(
            accepted.energy().reflected_j[lane],
            incident * (1.0 - absorbed),
        );
        close(accepted.energy().absorbed_j[lane], incident * absorbed);
        assert_eq!(accepted.energy().transmitted_j[lane], 0.0);
    }
    for (index, id) in [11, 12, 13].into_iter().enumerate() {
        let preexisting_heat = if id == 13 { 6.0 } else { 0.0 };
        close(
            next.cells()[&CellId(id)].specific_h_j_kg * (index + 1) as f64,
            expected[index] + preexisting_heat,
        );
        close(
            next.temperature_k(CellId(id)).unwrap(),
            300.0 + (expected[index] + preexisting_heat) / (4.0 * (index + 1) as f64),
        );
    }
    close(accepted.energy().thermal_input_j, expected.iter().sum());
    close(audit.external_heat_j, expected.iter().sum::<f64>() + 6.0);
    assert_ne!(expected[0], expected[1]);
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(input.time_s(), 0.0);
}

#[test]
fn shared_receivers_and_source_permutation_preserve_single_incident_account() {
    let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
    let a = source(7, [11; 3]);
    let b = source(3, [11; 3]);
    let first = heating(&original, &[a.clone(), b.clone()]);
    let second = heating(&original, &[b, a]);
    assert_eq!(first.model_identity(), second.model_identity());
    let (a, a_heat, _) = original
        .advance_observed(&first, 2.0, config(0.125), || Ok(()))
        .unwrap();
    let (b, b_heat, _) = original
        .advance_observed(&second, 2.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(a.cells(), b.cells());
    assert_eq!(a_heat.energy(), b_heat.energy());
    assert_eq!(a_heat.source_count(), 2);
    close(a_heat.energy().incident_j[0], 24.0);
    close(a_heat.energy().incident_j[1], 16.0);
    close(
        a.cells()[&CellId(11)].specific_h_j_kg,
        a_heat.energy().thermal_input_j,
    );
}

#[test]
fn every_layer_binding_and_temperature_domain_is_enforced_transactionally() {
    let original = thermal(
        vec![
            cell(11, 1.0, 1.0, 0.0),
            cell(12, 1.0, 1.0, 0.0),
            cell(13, 1.0, 1.0, 0.0),
        ],
        &SENSIBLE,
        vec![],
    );
    let valid = source(1, [11, 12, 13]);
    let admit = |s| {
        SurfaceOpticalHeating::try_new_coated(
            &[s],
            &original,
            budget(),
            SurfaceHeatingTolerance {
                power_roundoff_w: 1e-10,
            },
            || Ok(()),
        )
    };
    let mut missing = valid.clone();
    missing.film_receivers.pop();
    assert!(matches!(
        admit(missing),
        Err(SurfaceHeatingError::Invalid(_))
    ));
    let mut foreign = valid.clone();
    foreign.film_receivers[0].material_id = MaterialId(8);
    assert!(matches!(
        admit(foreign),
        Err(SurfaceHeatingError::Invalid(_))
    ));
    let mut inactive = valid.clone();
    inactive.film_receivers[0].thermal_cell = CellId(99);
    assert!(matches!(
        admit(inactive),
        Err(SurfaceHeatingError::Invalid(_))
    ));
    let mut limited = valid.clone();
    limited.film_receivers[0].temperature_bounds_k = (300.0, 300.01);
    let input = admit(limited).unwrap();
    assert!(matches!(
        original.advance_observed(&input, 2.0, config(0.1), || Ok(())),
        Err(SurfaceHeatingError::TemperatureDomain {
            thermal_cell: CellId(11),
            ..
        })
    ));
    assert_eq!(input.time_s(), 0.0);
    assert_eq!(input.energy().thermal_input_j, 0.0);
    assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 0.0);
    let input = admit(valid).unwrap();
    let mut checkpoints = 0;
    original
        .advance_observed(&input, 2.0, config(0.1), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        original.advance_observed(&input, 2.0, config(0.1), || {
            calls += 1;
            if calls == checkpoints - 2 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(SurfaceHeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(input.time_s(), 0.0);
    assert_eq!(original.time_s(), 0.0);
}
