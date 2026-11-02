//! Native transmitting-film state, shared packets and mixed source accounting.

use super::*;
use fs_render::conductor::{CoherentCoatingStack, ConductorIorSample, ConductorOptics};
use fs_render::thermal_surface::heating::feedback::SurfaceConductorLaw;
use fs_render::thermal_surface::heating::{
    CoatedSurfaceIlluminationSource, CoatedSurfacePowerResponse, CoatingOpticsModel,
    CoatingSurfaceSource, StatefulCoatedSurfaceIlluminationSource, StatefulCoatingLayer,
    StatefulTransmittingSurfaceIlluminationSource, SubstrateOpticsModel, SurfaceCoatingLaw,
    SurfaceIncidence,
};
use std::sync::Arc;

fn incidence() -> SurfaceIncidence<2> {
    SurfaceIncidence {
        wavelengths_nm: WAVELENGTHS,
        incident_eta: [1.0; 2],
        incidence_cosine: 0.8,
        s_power_fraction: POLARIZATION,
    }
}
fn constant_power(watts: [f64; 2]) -> SpectralPowerSchedule<2> {
    SpectralPowerSchedule::try_new(
        &[SpectralPowerSegment {
            start_s: 0.0,
            end_s: 1.0,
            start_power_w: watts,
            end_power_w: watts,
        }],
        1,
    )
    .unwrap()
}
fn transmitting(
    id: u64,
    films: Vec<StatefulCoatingLayer>,
    watts: [f64; 2],
) -> StatefulTransmittingSurfaceIlluminationSource<2> {
    StatefulTransmittingSurfaceIlluminationSource {
        source_id: id,
        films,
        incidence: incidence(),
        transmitted_eta: [1.5; 2],
        reverse: false,
        source_identity: SOURCE,
        power: constant_power(watts),
    }
}
fn admit(
    native: &DepositionSimulation<'_>,
    sources: &[CoatingSurfaceSource<'_, 2>],
    rows: usize,
) -> Result<SurfaceOpticalHeating<2>, SurfaceHeatingError> {
    SurfaceOpticalHeating::try_new_coating_stacks(
        sources,
        native,
        budget(),
        rows,
        tolerance(),
        || Ok(()),
    )
}
fn metal() -> ConductorOptics {
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.5, 3.0).unwrap()
        }),
        optical_source(),
    )
    .unwrap()
}
fn opaque_response(
    layers: &[CoatingLayer<'_>],
    substrate: &ConductorOptics,
) -> CoatedSurfacePowerResponse<2> {
    CoatedSurfacePowerResponse::from_coated_conductor(
        WAVELENGTHS,
        substrate,
        CoherentCoatingStack::try_new(layers).unwrap(),
        [1.0; 2],
        0.8,
        POLARIZATION,
    )
    .unwrap()
}

#[test]
fn all_four_source_kinds_compose_with_static_bits_and_canonical_receiver_maps() {
    let coating = film(1.8, 0.2, 90.0);
    let clear = clear_film();
    let auxiliary = film(2.0, 0.3, 30.0);
    let opaque_film = film(1.5, 0.1, 60.0);
    let substrate = metal();
    let film_law = Arc::new(
        SurfaceCoatingLaw::try_temperature(
            optical_source(),
            &[(300.0, coating), (600.0, coating)],
            2,
        )
        .unwrap(),
    );
    let substrate_law = Arc::new(
        SurfaceConductorLaw::try_temperature(
            optical_source(),
            &[(300.0, substrate), (600.0, substrate)],
            2,
        )
        .unwrap(),
    );
    let varying_transmission = transmitting(
        1,
        vec![
            StatefulCoatingLayer {
                receiver: receiver(11),
                optics: CoatingOpticsModel::Law(film_law),
            },
            StatefulCoatingLayer {
                receiver: receiver(12),
                optics: CoatingOpticsModel::Lossless(clear),
            },
        ],
        POWER,
    );
    let frozen_transmission = TransmittingSurfaceIlluminationSource {
        source_id: 1,
        film_receivers: vec![receiver(11), receiver(12)],
        response: response(
            &[
                CoatingLayer::Absorbing(&coating),
                CoatingLayer::Dielectric(&clear),
            ],
            1.0,
            1.5,
            0.8,
            false,
        ),
        power: power(),
    };
    let fixed_transmission = TransmittingSurfaceIlluminationSource {
        source_id: 2,
        film_receivers: vec![receiver(13)],
        response: response(&[CoatingLayer::Absorbing(&auxiliary)], 1.0, 1.5, 0.8, false),
        power: power(),
    };
    let varying_opaque = StatefulCoatedSurfaceIlluminationSource {
        source_id: 3,
        films: vec![StatefulCoatingLayer {
            receiver: receiver(14),
            optics: CoatingOpticsModel::Fixed(opaque_film),
        }],
        substrate_receiver: receiver(15),
        substrate: SubstrateOpticsModel::Law(substrate_law),
        incidence: incidence(),
        power: power(),
    };
    let frozen_opaque = CoatedSurfaceIlluminationSource {
        source_id: 3,
        film_receivers: vec![receiver(14)],
        substrate_receiver: receiver(15),
        response: opaque_response(&[CoatingLayer::Absorbing(&opaque_film)], &substrate),
        power: power(),
    };
    let fixed_opaque = CoatedSurfaceIlluminationSource {
        source_id: 4,
        film_receivers: vec![],
        substrate_receiver: receiver(15),
        response: opaque_response(&[], &substrate),
        power: power(),
    };
    let native = thermal(
        (11..=15)
            .map(|id| cell(id, (id - 10) as f64, 1.0, 0.0))
            .collect(),
        &SENSIBLE,
        vec![],
    );
    let frozen = admit(
        &native,
        &[
            CoatingSurfaceSource::FixedOpaque(&fixed_opaque),
            CoatingSurfaceSource::FixedTransmitting(&frozen_transmission),
            CoatingSurfaceSource::FixedOpaque(&frozen_opaque),
            CoatingSurfaceSource::FixedTransmitting(&fixed_transmission),
        ],
        0,
    )
    .unwrap();
    let mut inputs = [
        CoatingSurfaceSource::FixedOpaque(&fixed_opaque),
        CoatingSurfaceSource::StatefulTransmitting(&varying_transmission),
        CoatingSurfaceSource::StatefulOpaque(&varying_opaque),
        CoatingSurfaceSource::FixedTransmitting(&fixed_transmission),
    ];
    let dynamic = admit(&native, &inputs, 4).unwrap();
    inputs.reverse();
    let permuted = admit(&native, &inputs, 4).unwrap();
    assert_eq!(dynamic.model_identity(), permuted.model_identity());
    assert!(matches!(
        admit(&native, &inputs, 3),
        Err(SurfaceHeatingError::Budget)
    ));
    let (frozen_hot, frozen_heat, _) = native
        .advance_observed(&frozen, 1.0, config(0.125), || Ok(()))
        .unwrap();
    let (dynamic_hot, dynamic_heat, _) = native
        .advance_observed(&dynamic, 1.0, config(0.125), || Ok(()))
        .unwrap();
    let (permuted_hot, permuted_heat, _) = native
        .advance_observed(&permuted, 1.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(dynamic_hot.cells(), frozen_hot.cells());
    assert_eq!(dynamic_heat.energy(), frozen_heat.energy());
    assert_eq!(dynamic_hot.cells(), permuted_hot.cells());
    assert_eq!(dynamic_heat.energy(), permuted_heat.energy());
    assert_eq!(dynamic_heat.response_evaluations(), 0);
    assert_eq!(dynamic_heat.source_count(), 4);
    assert_eq!(dynamic_hot.cells()[&CellId(12)].specific_h_j_kg, 0.0);
    for lane in 0..2 {
        close(dynamic_heat.energy().incident_j[lane], 4.0 * POWER[lane]);
        assert!(dynamic_heat.energy().transmitted_j[lane] > 0.0);
    }
    close(
        dynamic_hot.energy_j(),
        dynamic_heat.energy().thermal_input_j,
    );
}

struct IdentifiedLatent;
impl ThermalLaw for IdentifiedLatent {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        LATENT.temperature_k(h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        LATENT.maximum_temperature_slope()
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([244; 32])
    }
}
static IDENTIFIED_LATENT: IdentifiedLatent = IdentifiedLatent;
fn latent_film(h: f64) -> CoherentAbsorbingCoating {
    film(1.4 + 0.01 * (h - 10.0), 0.02 + 0.0078 * (h - 10.0), 90.0)
}
fn latent_watts(h: f64, reverse: bool, watts: [f64; 2]) -> f64 {
    let coating = latent_film(h);
    let layers = [CoatingLayer::Absorbing(&coating)];
    let stack = DielectricCoatingStack::try_mixed(&layers).unwrap();
    let (ni, nt, cosine) = if reverse {
        (1.5, 1.0, 0.5)
    } else {
        (1.0, 1.5, 0.8)
    };
    (0..2)
        .map(|lane| {
            let work = stack
                .layer_power_balance(WAVELENGTHS[lane], reverse, cosine, ni, nt)
                .unwrap();
            let a = work.layer_absorbed()[0];
            watts[lane] * (POLARIZATION[lane] * a[0] + (1.0 - POLARIZATION[lane]) * a[1])
        })
        .sum()
}

#[test]
fn latent_reverse_and_forward_packets_share_only_identical_bindings_and_roll_back() {
    let mut cells = vec![cell(11, 1.0, 1.0, 0.0), cell(12, 1.0, 1.0, 0.0)];
    cells[0].specific_h_j_kg = 30.0;
    cells[1].specific_h_j_kg = 90.0;
    let native = thermal(cells, &IDENTIFIED_LATENT, vec![]);
    let law = Arc::new(
        SurfaceCoatingLaw::try_enthalpy(
            optical_source(),
            ContentHash([244; 32]),
            (310.0, 310.0),
            &[(10.0, latent_film(10.0)), (110.0, latent_film(110.0))],
            2,
        )
        .unwrap(),
    );
    let make = |id, cell, watts| {
        transmitting(
            id,
            vec![StatefulCoatingLayer {
                receiver: receiver(cell),
                optics: CoatingOpticsModel::Law(Arc::clone(&law)),
            }],
            watts,
        )
    };
    let mut reverse = make(4, 11, [1.0, 1.0]);
    reverse.reverse = true;
    reverse.incidence.incident_eta = [1.5; 2];
    reverse.incidence.incidence_cosine = 0.5;
    reverse.transmitted_eta = [1.0; 2];
    let sources = [
        make(3, 11, [3.0, 2.0]),
        make(1, 11, [2.0, 1.0]),
        make(2, 12, [5.0, 3.0]),
        reverse,
    ];
    let mut inputs: Vec<_> = sources
        .iter()
        .map(CoatingSurfaceSource::StatefulTransmitting)
        .collect();
    let input = admit(&native, &inputs, 2).unwrap();
    inputs.reverse();
    let permuted = admit(&native, &inputs, 2).unwrap();
    assert_eq!(input.model_identity(), permuted.model_identity());
    assert!(matches!(
        admit(&native, &inputs, 1),
        Err(SurfaceHeatingError::Budget)
    ));
    let watts = |id, h| {
        latent_watts(h, false, [5.0, 3.0])
            + if id == 11 {
                latent_watts(h, true, [1.0, 1.0])
            } else {
                0.0
            }
    };
    let (first, accepted, _) = native
        .advance_observed(&input, 0.125, config(0.125), || Ok(()))
        .unwrap();
    for (id, initial) in [(11, 30.0), (12, 90.0)] {
        close(
            first.cells()[&CellId(id)].specific_h_j_kg,
            initial + 0.125 * watts(id, initial),
        );
        close(first.temperature_k(CellId(id)).unwrap(), 310.0);
    }
    assert!(
        (latent_watts(30.0, false, [5.0, 3.0]) - latent_watts(90.0, false, [5.0, 3.0])).abs() > 0.1
    );
    let (second, refreshed, _) = first
        .advance_observed(&accepted, 0.25, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(refreshed.response_evaluations(), 3);
    for id in [11, 12] {
        let h = first.cells()[&CellId(id)].specific_h_j_kg;
        close(
            second.cells()[&CellId(id)].specific_h_j_kg,
            h + 0.125 * watts(id, h),
        );
        close(second.temperature_k(CellId(id)).unwrap(), 310.0);
    }
    let mut checkpoints = 0;
    let (whole, complete, _) = native
        .advance_observed(&input, 1.0, config(0.125), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    assert_eq!(complete.response_evaluations(), 21);
    let (middle, partial, _) = native
        .advance_observed(&permuted, 0.5, config(0.125), || Ok(()))
        .unwrap();
    let (split, continued, _) = middle
        .advance_observed(&partial, 1.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(whole.cells(), split.cells());
    assert_eq!(complete.energy(), continued.energy());
    close(complete.energy().incident_j.iter().sum(), 18.0);
    close(
        whole.energy_j() - native.energy_j(),
        complete.energy().thermal_input_j,
    );
    assert!(complete.energy().transmitted_j.iter().all(|&x| x > 0.0));
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
    assert_eq!(input.time_s(), 0.0);
    assert_eq!(input.energy().thermal_input_j, 0.0);
    assert_eq!(input.response_evaluations(), 0);
    let (retry, retried, _) = native
        .advance_observed(&input, 1.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(retry.cells(), whole.cells());
    assert_eq!(retried.energy(), complete.energy());
    assert_eq!(retried.response_evaluations(), 21);
}

#[test]
fn dark_transmitting_film_domains_and_caloric_identity_remain_binding() {
    let first = film(1.5, 0.1, 90.0);
    let last = film(2.0, 0.5, 90.0);
    let law = Arc::new(
        SurfaceCoatingLaw::try_temperature(optical_source(), &[(300.0, first), (300.5, last)], 2)
            .unwrap(),
    );
    let native = thermal(
        vec![cell(11, 1.0, 1.0, 8.0), cell(12, 1.0, 1.0, 0.0)],
        &SENSIBLE,
        vec![],
    );
    let mut source = transmitting(
        1,
        vec![
            StatefulCoatingLayer {
                receiver: receiver(11),
                optics: CoatingOpticsModel::Law(law),
            },
            StatefulCoatingLayer {
                receiver: receiver(12),
                optics: CoatingOpticsModel::Lossless(clear_film()),
            },
        ],
        [0.0; 2],
    );
    source.power = SpectralPowerSchedule::try_new(&[], 0).unwrap();
    let input = admit(
        &native,
        &[CoatingSurfaceSource::StatefulTransmitting(&source)],
        2,
    )
    .unwrap();
    let (partial, accepted, audit) = native
        .advance_observed(&input, 0.25, config(0.125), || Ok(()))
        .unwrap();
    close(partial.temperature_k(CellId(11)).unwrap(), 300.5);
    close(audit.external_heat_j, 2.0);
    assert_eq!(accepted.response_evaluations(), 0);
    assert_eq!(accepted.energy().incident_j, [0.0; 2]);
    assert_eq!(accepted.energy().thermal_input_j, 0.0);
    assert!(matches!(
        partial.advance_observed(&accepted, 0.5, config(0.125), || Ok(())),
        Err(SurfaceHeatingError::TemperatureDomain {
            thermal_cell: CellId(11),
            ..
        })
    ));
    assert!(matches!(
        native.advance_observed(&input, 1.0, config(0.125), || Ok(())),
        Err(SurfaceHeatingError::TemperatureDomain {
            thermal_cell: CellId(11),
            ..
        })
    ));
    assert_eq!(native.time_s(), 0.0);
    assert_eq!(native.cells()[&CellId(11)].specific_h_j_kg, 0.0);
    assert_eq!(input.time_s(), 0.0);
    assert_eq!(partial.time_s(), 0.25);
    assert_eq!(accepted.time_s(), 0.25);
    assert_eq!(accepted.response_evaluations(), 0);
    let mut unsupported_clear = source.clone();
    unsupported_clear.films[1].receiver.temperature_bounds_k = (301.0, 600.0);
    assert!(matches!(
        admit(
            &native,
            &[CoatingSurfaceSource::StatefulTransmitting(
                &unsupported_clear
            )],
            2
        ),
        Err(SurfaceHeatingError::TemperatureDomain {
            thermal_cell: CellId(12),
            ..
        })
    ));
    let mut foreign = source;
    foreign.films[0].optics = CoatingOpticsModel::Law(Arc::new(
        SurfaceCoatingLaw::try_enthalpy(
            optical_source(),
            ContentHash([244; 32]),
            (300.0, 600.0),
            &[(0.0, first), (1000.0, last)],
            2,
        )
        .unwrap(),
    ));
    assert!(
        admit(
            &native,
            &[CoatingSurfaceSource::StatefulTransmitting(&foreign)],
            2
        )
        .is_err()
    );
}
