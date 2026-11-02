//! Native state, cache and physical port ownership for transmitting pairs.
use super::transmitting_pair::{ir_film, pair, vacuum};
use super::*;
use fs_render::thermal_surface::heating::coating_exchange::transmitting_pair::*;

fn limits() -> NativeCoatingTransmittingBudget {
    NativeCoatingTransmittingBudget {
        max_pairs: 4,
        max_total_optical_knots: 4,
        max_source_evaluations: 10_000,
        max_optical_evaluations: 1_000_000,
        max_channel_spectral_nodes: 4 * 6 * 8,
    }
}
fn bodies(temperatures: [f64; 2]) -> Vec<DepositedCell> {
    temperatures
        .into_iter()
        .enumerate()
        .map(|(side, t)| {
            let mut body = specimen(t, 1.0);
            body.id = CellId(11 + side as u64);
            body.boundary.ambient_k = if side == 0 { 260.0 } else { 290.0 };
            body
        })
        .collect()
}
fn fixed_pair() -> CoatingTransmittingPair {
    let mut p = pair();
    for (side, surface) in p.surfaces.iter_mut().enumerate() {
        surface.stack = IsothermalCoatingStack::Transmitting {
            films: vec![CoatingOpticsModel::Fixed(ir_film(
                1.4 + side as f64,
                0.1,
                2e-6,
            ))],
            opposite_ior: vacuum(),
            reverse: side == 1,
        };
    }
    p
}
fn prepare(
    native: &DepositionSimulation<'_>,
    pair: &CoatingTransmittingPair,
) -> NativeCoatingTransmittingExchange {
    NativeCoatingTransmittingExchange::try_new(native, &[pair.clone()], limits(), || Ok(()))
        .unwrap()
}

#[test]
fn transmitted_bath_radiation_never_heats_transparent_bodies_and_both_faces_occupy_area() {
    let native = thermal(bodies([340.0, 280.0]), &HOT, vec![]);
    let mut p = fixed_pair();
    for side in &mut p.surfaces {
        if let IsothermalCoatingStack::Transmitting { films, .. } = &mut side.stack {
            films.clear();
        }
    }
    p.area_m2 = 3.0; // Exactly two faces fill each original 6 m² cell surface.
    let initial = prepare(&native, &p);
    assert_eq!(
        initial.sources_w().values().copied().collect::<Vec<_>>(),
        vec![0.0, 0.0]
    );
    assert_eq!(initial.heat_transfers()[0].power_w, 0.0);
    let throughput = initial.rates(17).unwrap().bath_transfer.outward_power_w;
    let expected: f64 = (0..8)
        .map(|w| {
            let lambda = 4000.0 + (w as f64 + 0.5) * 2000.0;
            std::f64::consts::PI
                * p.area_m2
                * 2000.0
                * (planck(430.0, lambda) - planck(240.0, lambda))
        })
        .sum();
    near(throughput, expected, 3e-12);
    let (after, observed, audit) = native
        .advance_observed(&initial, 0.1, config(0.05), || Ok(()))
        .unwrap();
    assert_eq!(after.cells(), native.cells());
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(observed.energy().thermal_external_j, 0.0);
    assert_eq!(observed.energy().by_pair_first_to_second_j[&17], 0.0);
    near(
        observed.energy().by_pair_bath_transfer_j[&17],
        throughput * 0.1,
        2e-12,
    );
    assert_eq!(observed.work().optical_evaluations, 0);

    let mut excess = p.clone();
    excess.area_m2 += 1e-6;
    assert!(matches!(
        NativeCoatingTransmittingExchange::try_new(&native, &[excess], limits(), || Ok(())),
        Err(NativeCoatingExchangeError::Invalid(_))
    ));
    let contacted = thermal(
        bodies([340.0, 280.0]),
        &HOT,
        vec![TimedEvent {
            time_s: 0.001,
            event: Event::Contact(ThermalLink {
                id: 1,
                a: CellId(11),
                b: CellId(12),
                area_m2: 0.1,
                conductance_w_k: 0.0,
            }),
        }],
    );
    let before = prepare(&contacted, &p);
    assert!(matches!(
        contacted.advance_observed(&before, 0.001, config(1.0), || Ok(())),
        Err(NativeCoatingExchangeError::Invalid(_))
    ));
    assert!(contacted.links().is_empty());
    assert_eq!(before.time_s(), 0.0);
}

#[test]
fn fixed_optics_reuse_all_six_channels_across_bath_events_and_split_continuation() {
    let cells = bodies([340.0, 280.0]);
    let mut warm = cells[0].boundary;
    warm.ambient_k = 420.0;
    let native = thermal(
        cells,
        &HOT,
        vec![TimedEvent {
            time_s: 0.05,
            event: Event::Boundary(CellId(11), warm),
        }],
    );
    let mut p = fixed_pair();
    p.baths = [
        CoatingExchangeBath::NativeAmbient,
        CoatingExchangeBath::NativeAmbient,
    ];
    let initial = prepare(&native, &p);
    let optics = [
        initial.emitter(17, 0).unwrap().optical_identity(),
        initial.emitter(17, 1).unwrap().optical_identity(),
    ];
    let (whole, observed, audit) = native
        .advance_observed(&initial, 0.1, config(0.025), || Ok(()))
        .unwrap();
    let (middle, partial, _) = native
        .advance_observed(&initial, 0.05, config(0.025), || Ok(()))
        .unwrap();
    let (split, resumed, _) = middle
        .advance_observed(&partial, 0.1, config(0.025), || Ok(()))
        .unwrap();
    assert!(initial.rates(17).unwrap().first_baths[0].outward_power_w > 0.0);
    assert!(partial.rates(17).unwrap().first_baths[0].outward_power_w < 0.0);
    // The left bath reaches the second body through the first transmitting body.
    assert!(initial.rates(17).unwrap().second_baths[0].outward_power_w > 0.0);
    assert!(partial.rates(17).unwrap().second_baths[0].outward_power_w < 0.0);
    assert_eq!(observed.work().optical_evaluations, 0);
    near(whole.energy_j(), split.energy_j(), 2e-13);
    near(
        observed.energy().thermal_external_j,
        resumed.energy().thermal_external_j,
        2e-13,
    );
    near(
        observed.energy().thermal_external_j,
        audit.external_heat_j,
        2e-13,
    );
    near(observed.energy().balance_residual_j, 0.0, 1e-13);
    for side in 0..2 {
        let id = CellId(11 + side as u64);
        assert_eq!(
            observed.emitter(17, side).unwrap().optical_identity(),
            optics[side]
        );
        assert_eq!(
            observed.emitter(17, side).unwrap().film_temperature_k(0),
            Some(whole.temperature_k(id).unwrap())
        );
        near(
            10.0 * (whole.temperature_k(id).unwrap() - native.temperature_k(id).unwrap()),
            observed.energy().by_cell_internal_j[&id] + observed.energy().by_cell_external_j[&id],
            2e-11,
        );
    }
    // Canonical order and summed area give the same physical source.
    let mut a = p.clone();
    a.area_m2 *= 0.5;
    let mut b = a.clone();
    b.pair_id = 18;
    let ordered = NativeCoatingTransmittingExchange::try_new(
        &native,
        &[a.clone(), b.clone()],
        limits(),
        || Ok(()),
    )
    .unwrap();
    let reverse =
        NativeCoatingTransmittingExchange::try_new(&native, &[b, a], limits(), || Ok(())).unwrap();
    assert_eq!(ordered.identity(), reverse.identity());
    for (&id, &watts) in initial.sources_w() {
        near(ordered.sources_w()[&id], watts, 1e-14);
    }
}

struct Latent;
impl ThermalLaw for Latent {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && (10.0..=110.0).contains(&h) {
            Ok(310.0)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.0
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([239; 32])
    }
}

#[test]
fn independent_enthalpy_optics_refresh_on_a_shared_latent_temperature_plateau() {
    let mut cells = bodies([310.0, 310.0]);
    for (side, c) in cells.iter_mut().enumerate() {
        c.mass_kg = 0.001;
        c.specific_h_j_kg = 70.0 - side as f64 * 25.0;
    }
    let native = thermal(cells, &Latent, vec![]);
    let law = Arc::new(
        SurfaceCoatingLaw::try_enthalpy(
            authority(),
            ContentHash([239; 32]),
            (310.0, 310.0),
            &[
                (10.0, ir_film(1.5, 0.01, 2e-6)),
                (110.0, ir_film(1.5, 0.6, 2e-6)),
            ],
            2,
        )
        .unwrap(),
    );
    let mut p = fixed_pair();
    p.area_m2 = 0.0001;
    p.baths = [CoatingExchangeBath::Blackbody {
        temperature_k: 250.0,
    }; 2];
    for side in &mut p.surfaces {
        side.receiver.temperature_bounds_k = (310.0, 310.0);
        if let IsothermalCoatingStack::Transmitting { films, .. } = &mut side.stack {
            *films = vec![CoatingOpticsModel::Law(Arc::clone(&law))];
        }
    }
    let mut budget = limits();
    budget.max_total_optical_knots = 2;
    let initial =
        NativeCoatingTransmittingExchange::try_new(&native, &[p.clone()], budget, || Ok(()))
            .unwrap();
    let (after, observed, audit) = native
        .advance_observed(&initial, 0.04, config(0.005), || Ok(()))
        .unwrap();
    for side in 0..2 {
        let id = CellId(11 + side as u64);
        assert_eq!(after.temperature_k(id).unwrap(), 310.0);
        assert!(after.cells()[&id].specific_h_j_kg < native.cells()[&id].specific_h_j_kg);
        assert_ne!(
            observed.emitter(17, side).unwrap().optical_identity(),
            initial.emitter(17, side).unwrap().optical_identity()
        );
    }
    assert_eq!(observed.energy().by_pair_first_to_second_j[&17], 0.0);
    assert!(observed.work().optical_evaluations >= 8 * 4 * 8 * 4);
    near(
        after.energy_j() - native.energy_j(),
        audit.external_heat_j,
        1e-12,
    );
    budget.max_total_optical_knots = 1;
    assert!(matches!(
        NativeCoatingTransmittingExchange::try_new(&native, &[p], budget, || Ok(())),
        Err(NativeCoatingExchangeError::Budget)
    ));
}

#[test]
fn retained_work_limits_and_late_cancellation_preserve_the_accepted_pair() {
    let native = thermal(bodies([340.0, 280.0]), &HOT, vec![]);
    let p = fixed_pair();
    for budget in [
        NativeCoatingTransmittingBudget {
            max_channel_spectral_nodes: 47,
            ..limits()
        },
        NativeCoatingTransmittingBudget {
            max_optical_evaluations: 127,
            ..limits()
        },
    ] {
        assert!(matches!(
            NativeCoatingTransmittingExchange::try_new(&native, &[p.clone()], budget, || Ok(())),
            Err(NativeCoatingExchangeError::Budget)
        ));
    }
    assert!(
        NativeCoatingTransmittingExchange::try_new(
            &native,
            &[p.clone(), p.clone()],
            limits(),
            || Ok(())
        )
        .is_err()
    );
    let initial = prepare(&native, &p);
    let limited = NativeCoatingTransmittingExchange::try_new(
        &native,
        &[p.clone()],
        NativeCoatingTransmittingBudget {
            max_source_evaluations: 1,
            ..limits()
        },
        || Ok(()),
    )
    .unwrap();
    assert!(matches!(
        native.advance_observed(&limited, 0.01, config(0.01), || Ok(())),
        Err(NativeCoatingExchangeError::Budget)
    ));
    let mut calls = 0;
    let accepted = native
        .advance_observed(&initial, 0.01, config(0.01), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    let mut cancelled = 0;
    assert!(matches!(
        native.advance_observed(&initial, 0.01, config(0.01), || {
            cancelled += 1;
            if cancelled == calls - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(NativeCoatingExchangeError::Thermal(
            DepositionError::Cancelled
        ))
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(native.time_s(), 0.0);
    assert_eq!(
        initial.energy(),
        &NativeCoatingTransmittingEnergy::default()
    );
    let retry = native
        .advance_observed(&initial, 0.01, config(0.01), || Ok(()))
        .unwrap();
    assert_eq!(
        retry.0.energy_j().to_bits(),
        accepted.0.energy_j().to_bits()
    );
    assert_eq!(retry.1.energy(), accepted.1.energy());
}
