//! G1/G3: finite-film heat is separate from outgoing transmitted radiation.

#[path = "transmitting_coating/refinement.rs"]
mod refinement;
#[path = "transmitting_coating/state.rs"]
mod state;

use super::*;
use fs_blake3::ContentHash;
use fs_render::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, CoherentDielectricCoating,
    ConductorDataStatus, ConductorSource,
};
use fs_render::dielectric::{CauchyIor, DielectricCoatingStack};
use fs_render::thermal_surface::heating::{
    SurfaceHeatingBudget, SurfaceHeatingError, SurfaceHeatingReceiver, SurfaceHeatingTolerance,
    SurfaceOpticalHeating, SurfacePowerResponse, TransmittingSurfaceIlluminationSource,
    TransmittingSurfacePowerResponse,
};
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};

const WAVELENGTHS: [f64; 2] = [480.0, 680.0];
const POLARIZATION: [f64; 2] = [0.2, 0.8];
const POWER: [f64; 2] = [8.0, 12.0];
const SOURCE: ContentHash = ContentHash([242; 32]);

fn optical_source() -> ConductorSource {
    ConductorSource::try_new(ContentHash([241; 32]), ConductorDataStatus::Representative).unwrap()
}
fn film(eta: f64, k: f64, depth_nm: f64) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new(
        core::array::from_fn(|i| {
            CoatingIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()
        }),
        depth_nm * 1e-9,
        optical_source(),
    )
    .unwrap()
}
fn clear_film() -> CoherentDielectricCoating {
    CoherentDielectricCoating::try_new(
        CauchyIor::try_new(1.4, 0.004, 0.0).unwrap(),
        80e-9,
        optical_source(),
    )
    .unwrap()
}
fn receiver(id: u64) -> SurfaceHeatingReceiver {
    SurfaceHeatingReceiver {
        thermal_cell: CellId(id),
        material_id: MaterialId(7),
        temperature_bounds_k: (300.0, 600.0),
    }
}
fn power() -> SpectralPowerSchedule<2> {
    SpectralPowerSchedule::try_new(
        &[SpectralPowerSegment {
            start_s: 0.0,
            end_s: 1.0,
            start_power_w: POWER,
            end_power_w: POWER,
        }],
        1,
    )
    .unwrap()
}
fn response(
    layers: &[CoatingLayer<'_>],
    eta_i: f64,
    eta_t: f64,
    cosine: f64,
    reverse: bool,
) -> TransmittingSurfacePowerResponse<2> {
    TransmittingSurfacePowerResponse::from_dielectric_stack(
        WAVELENGTHS,
        layers,
        [eta_i; 2],
        [eta_t; 2],
        cosine,
        reverse,
        POLARIZATION,
        SOURCE,
    )
    .unwrap()
}
fn source(
    response: TransmittingSurfacePowerResponse<2>,
    ids: &[u64],
) -> TransmittingSurfaceIlluminationSource<2> {
    TransmittingSurfaceIlluminationSource {
        source_id: 1,
        film_receivers: ids.iter().copied().map(receiver).collect(),
        response,
        power: power(),
    }
}
fn budget() -> SurfaceHeatingBudget {
    SurfaceHeatingBudget {
        max_sources: 4,
        max_total_power_segments: 4,
        max_source_evaluations: 128,
    }
}
fn tolerance() -> SurfaceHeatingTolerance {
    SurfaceHeatingTolerance {
        power_roundoff_w: 1e-10,
    }
}
fn heating(
    native: &DepositionSimulation<'_>,
    sources: &[TransmittingSurfaceIlluminationSource<2>],
) -> SurfaceOpticalHeating<2> {
    SurfaceOpticalHeating::try_new_transmitting(sources, native, budget(), tolerance(), || Ok(()))
        .unwrap()
}

#[test]
fn asymmetric_forward_reverse_and_tir_heat_authored_films_only() {
    let front = film(1.3, 0.05, 45.0);
    let clear = clear_film();
    let back = film(2.2, 0.4, 120.0);
    let layers = [
        CoatingLayer::Absorbing(&front),
        CoatingLayer::Dielectric(&clear),
        CoatingLayer::Absorbing(&back),
    ];
    let mirrored = [
        CoatingLayer::Absorbing(&back),
        CoatingLayer::Dielectric(&clear),
        CoatingLayer::Absorbing(&front),
    ];
    let native = thermal(
        vec![
            cell(11, 1.0, 1.0, 0.0),
            cell(12, 2.0, 1.0, 0.0),
            cell(13, 3.0, 1.0, 0.0),
            cell(14, 2.0, 1.0, 3.0),
        ],
        &SENSIBLE,
        vec![],
    );
    let reverse_cosine = (1.0_f64 - (1.0 / 1.5_f64).powi(2) * (1.0 - 0.8_f64.powi(2))).sqrt();
    for (reverse, eta_i, eta_t, cosine) in [
        (false, 1.0, 1.5, 0.8),
        (true, 1.5, 1.0, reverse_cosine),
        (false, 1.6, 1.0, 0.5),
    ] {
        let compiled = response(&layers, eta_i, eta_t, cosine, reverse);
        let aggregate = SurfacePowerResponse::from_dielectric_stack(
            WAVELENGTHS,
            &layers,
            [eta_i; 2],
            [eta_t; 2],
            cosine,
            reverse,
            POLARIZATION,
            SOURCE,
        )
        .unwrap();
        assert_eq!(*compiled.total(), aggregate);
        let input = heating(&native, &[source(compiled, &[11, 12, 13])]);
        let (hot, accepted, audit) = native
            .advance_observed(&input, 1.0, config(0.125), || Ok(()))
            .unwrap();
        let mut expected = [0.0; 3];
        for lane in 0..2 {
            let original = DielectricCoatingStack::try_mixed(&layers)
                .unwrap()
                .layer_power_balance(WAVELENGTHS[lane], reverse, cosine, eta_i, eta_t)
                .unwrap();
            for i in 0..3 {
                let a = original.layer_absorbed()[i];
                expected[i] +=
                    POWER[lane] * (POLARIZATION[lane] * a[0] + (1.0 - POLARIZATION[lane]) * a[1]);
            }
            if reverse {
                let reversed = DielectricCoatingStack::try_mixed(&mirrored)
                    .unwrap()
                    .layer_power_balance(WAVELENGTHS[lane], false, cosine, eta_i, eta_t)
                    .unwrap();
                for i in 0..3 {
                    assert_eq!(
                        original.layer_absorbed()[i],
                        reversed.layer_absorbed()[2 - i]
                    );
                }
            }
            let reflected = original.fresnel.reflection.0[0][0]
                + (2.0 * POLARIZATION[lane] - 1.0) * original.fresnel.reflection.0[0][1];
            let transmitted = original.fresnel.transmission.0[0][0]
                + (2.0 * POLARIZATION[lane] - 1.0) * original.fresnel.transmission.0[0][1];
            close(accepted.energy().reflected_j[lane], POWER[lane] * reflected);
            close(
                accepted.energy().transmitted_j[lane],
                POWER[lane] * transmitted,
            );
            close(accepted.energy().incident_j[lane], POWER[lane]);
            close(
                accepted.energy().reflected_j[lane]
                    + accepted.energy().transmitted_j[lane]
                    + accepted.energy().absorbed_j[lane],
                POWER[lane],
            );
            if eta_i == 1.6 {
                assert_eq!(accepted.energy().transmitted_j[lane], 0.0);
                assert!(accepted.energy().absorbed_j[lane] > 0.0);
            }
        }
        for i in 0..3 {
            close(
                hot.cells()[&CellId(11 + i as u64)].specific_h_j_kg * (i + 1) as f64,
                expected[i],
            );
        }
        assert_eq!(expected[1], 0.0);
        assert!((expected[0] - expected[2]).abs() > 1e-3);
        close(hot.cells()[&CellId(14)].specific_h_j_kg, 1.5);
        close(accepted.energy().thermal_input_j, expected.iter().sum());
        close(audit.external_heat_j, expected.iter().sum::<f64>() + 3.0);
        assert_eq!(hot.cells().len(), 4);
    }
}

#[test]
fn weak_film_work_below_rt_roundoff_enters_each_native_enthalpy() {
    let front = film(1.0, 1e-20, 50.0);
    let back = film(1.0, 1e-20, 150.0);
    let layers = [
        CoatingLayer::Absorbing(&front),
        CoatingLayer::Absorbing(&back),
    ];
    let native = thermal(
        vec![cell(11, 1.0, 1.0, 0.0), cell(12, 2.0, 1.0, 0.0)],
        &SENSIBLE,
        vec![],
    );
    let input = heating(
        &native,
        &[source(response(&layers, 1.0, 1.0, 1.0, false), &[11, 12])],
    );
    let (hot, accepted, _) = native
        .advance_observed(&input, 1.0, config(0.25), || Ok(()))
        .unwrap();
    let mut total = 0.0;
    for (id, mass, depth_nm) in [(11, 1.0, 50.0), (12, 2.0, 150.0)] {
        // Matched real indices, infinitesimal loss: entrance and interior
        // reflection corrections are O(k^2), while work is 4*pi*k*d/lambda.
        let expected: f64 = (0..2)
            .map(|lane| {
                POWER[lane] * 4.0 * core::f64::consts::PI * 1e-20 * depth_nm / WAVELENGTHS[lane]
            })
            .sum();
        let stored = hot.cells()[&CellId(id)].specific_h_j_kg * mass;
        assert!(stored > 0.0);
        assert!(
            (stored / expected - 1.0).abs() < 2e-12,
            "{stored:e} != {expected:e}"
        );
        total += expected;
    }
    assert!((accepted.energy().thermal_input_j / total - 1.0).abs() < 2e-12);
    for lane in 0..2 {
        assert!(accepted.energy().transmitted_j[lane] > POWER[lane] * (1.0 - 1e-13));
        assert!(accepted.energy().absorbed_j[lane] > 0.0);
    }
}

#[test]
fn bare_and_dispersive_lossless_interfaces_do_not_invent_heat_receivers() {
    let empty = thermal(vec![], &SENSIBLE, vec![]);
    assert!(empty.cells().is_empty());
    let mut bare = source(response(&[], 1.0, 1.5, 1.0, false), &[]);
    bare.power = SpectralPowerSchedule::try_new(
        &[
            SpectralPowerSegment {
                start_s: 0.125,
                end_s: 0.375,
                start_power_w: POWER,
                end_power_w: POWER,
            },
            SpectralPowerSegment {
                start_s: 0.625,
                end_s: 0.875,
                start_power_w: POWER,
                end_power_w: [0.0; 2],
            },
        ],
        2,
    )
    .unwrap();
    let input = heating(&empty, &[bare]);
    let mut checkpoints = 0;
    let (advanced, accepted, audit) = empty
        .advance_observed(&input, 1.0, config(0.4), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    assert!(advanced.cells().is_empty());
    assert!(accepted.sources_w().is_empty());
    assert_eq!(accepted.energy().thermal_input_j, 0.0);
    assert_eq!(accepted.energy().absorbed_j, [0.0; 2]);
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(audit.substeps, 5);
    assert_eq!(accepted.time_s(), 1.0);
    close(accepted.energy().duration_s, 1.0);
    for lane in 0..2 {
        // One rectangular pulse, one ramp, and dark gaps: integral = 3P/8.
        close(accepted.energy().incident_j[lane], 0.375 * POWER[lane]);
        close(
            accepted.energy().reflected_j[lane],
            0.04 * 0.375 * POWER[lane],
        );
        close(
            accepted.energy().transmitted_j[lane],
            0.96 * 0.375 * POWER[lane],
        );
    }
    assert!(matches!(
        empty.advance_observed_limited(&input, 1.0, config(0.4), 2, || Ok(())),
        Err(SurfaceHeatingError::Thermal(DepositionError::Budget(_)))
    ));
    let mut calls = 0;
    assert!(matches!(
        empty.advance_observed(&input, 1.0, config(0.4), || {
            calls += 1;
            if calls == checkpoints - 2 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(SurfaceHeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(empty.time_s(), 0.0);
    assert_eq!(input.time_s(), 0.0);
    assert_eq!(input.energy().incident_j, [0.0; 2]);
    let (idle, (), idle_audit) = empty
        .advance_observed_limited(&(), 1.0, config(0.4), 0, || Ok(()))
        .unwrap();
    assert_eq!(idle.time_s(), 1.0);
    assert_eq!(idle_audit.substeps, 0);

    let clear = clear_film();
    let layers = [CoatingLayer::Dielectric(&clear)];
    let compiled = response(&layers, 1.0, 1.5, 0.8, false);
    assert_eq!(compiled.film_absorbed_fractions(), &[[0.0; 2]]);
    assert_ne!(
        compiled.total().reflected_fraction()[0],
        compiled.total().reflected_fraction()[1]
    );
    let native = thermal(vec![cell(11, 1.0, 1.0, 2.0)], &SENSIBLE, vec![]);
    let input = heating(&native, &[source(compiled, &[11])]);
    let (hot, accepted, _) = native
        .advance_observed(&input, 1.0, config(0.25), || Ok(()))
        .unwrap();
    assert_eq!(hot.cells().len(), 1);
    close(hot.cells()[&CellId(11)].specific_h_j_kg, 2.0);
    assert_eq!(accepted.energy().thermal_input_j, 0.0);
    assert!(accepted.sources_w().is_empty());
}

#[test]
fn every_finite_film_binding_budget_and_domain_is_transactional() {
    let film = film(1.8, 0.3, 90.0);
    let layers = [CoatingLayer::Absorbing(&film)];
    let native = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
    let valid = source(response(&layers, 1.0, 1.5, 0.8, false), &[11]);
    let admit = |sources: &[TransmittingSurfaceIlluminationSource<2>], limit| {
        SurfaceOpticalHeating::try_new_transmitting(sources, &native, limit, tolerance(), || Ok(()))
    };
    let mut missing = valid.clone();
    missing.film_receivers.clear();
    assert!(admit(&[missing], budget()).is_err());
    let mut foreign = valid.clone();
    foreign.film_receivers[0].material_id = MaterialId(8);
    assert!(admit(&[foreign], budget()).is_err());
    let mut inactive = valid.clone();
    inactive.film_receivers[0].thermal_cell = CellId(99);
    assert!(admit(&[inactive], budget()).is_err());
    assert!(admit(&[valid.clone(), valid.clone()], budget()).is_err());
    let mut limited_budget = budget();
    limited_budget.max_total_power_segments = 0;
    assert!(matches!(
        admit(&[valid.clone()], limited_budget),
        Err(SurfaceHeatingError::Budget)
    ));
    let mut limited = valid.clone();
    limited.film_receivers[0].temperature_bounds_k = (300.0, 300.01);
    let narrow = admit(&[limited], budget()).unwrap();
    assert!(matches!(
        native.advance_observed(&narrow, 1.0, config(0.125), || Ok(())),
        Err(SurfaceHeatingError::TemperatureDomain {
            thermal_cell: CellId(11),
            ..
        })
    ));
    let input = admit(&[valid], budget()).unwrap();
    let mut checkpoints = 0;
    let (expected, energy, _) = native
        .advance_observed(&input, 1.0, config(0.125), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        native.advance_observed(&input, 1.0, config(0.125), || {
            calls += 1;
            if calls == checkpoints - 2 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(SurfaceHeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(native.time_s(), 0.0);
    assert_eq!(native.cells()[&CellId(11)].specific_h_j_kg, 0.0);
    assert_eq!(input.time_s(), 0.0);
    assert_eq!(input.energy().thermal_input_j, 0.0);
    let (retry, retried, _) = native
        .advance_observed(&input, 1.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(retry.cells(), expected.cells());
    assert_eq!(retried.energy(), energy.energy());
}
