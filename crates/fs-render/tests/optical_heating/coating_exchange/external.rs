use super::*;
use fs_render::dielectric::{DielectricCoatingStack, DielectricIorSample, TabulatedIor};
use fs_render::thermal_surface::heating::StatefulCoatingLayer;
use fs_render::thermal_surface::heating::coating_exchange::external::*;

fn ir_film(n: f64, k: f64) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new_in_band(
        std::array::from_fn(|i| {
            CoatingIorSample::try_new_in_band(4000.0 + 2000.0 * i as f64, n, k).unwrap()
        }),
        2e-6,
        2000.0,
        authority(),
    )
    .unwrap()
}
fn ir_metal() -> ConductorOptics {
    ConductorOptics::try_new_in_band(
        std::array::from_fn(|i| {
            ConductorIorSample::try_new_in_band(4000.0 + 2000.0 * i as f64, 0.7, 2.5).unwrap()
        }),
        2000.0,
        authority(),
    )
    .unwrap()
}
fn vacuum() -> DielectricIor {
    TabulatedIor::try_new_in_band(
        std::array::from_fn(|i| {
            DielectricIorSample::try_new_in_band(4000.0 + 2000.0 * i as f64, 1.0).unwrap()
        }),
        2000.0,
    )
    .unwrap()
    .into()
}
fn body(id: u64, t: f64, mass: f64) -> DepositedCell {
    let mut c = specimen(t, mass);
    c.id = CellId(id);
    c
}
fn receiver(id: u64) -> SurfaceHeatingReceiver {
    SurfaceHeatingReceiver {
        thermal_cell: CellId(id),
        material_id: MaterialId(7),
        temperature_bounds_k: (100.0, 5000.0),
    }
}
fn layer(id: u64, optics: CoatingOpticsModel) -> StatefulCoatingLayer {
    StatefulCoatingLayer {
        receiver: receiver(id),
        optics,
    }
}
fn limits() -> NativeCoatingExternalBudget {
    NativeCoatingExternalBudget {
        exchange: budget(),
        max_channel_spectral_nodes: 4096,
    }
}
fn prepare(
    native: &DepositionSimulation<'_>,
    sources: &[CoatingExternalBathSource],
) -> NativeCoatingExternalBathExchange {
    NativeCoatingExternalBathExchange::try_new(native, sources, limits(), || Ok(())).unwrap()
}
fn opaque() -> CoatingExternalBathSource {
    CoatingExternalBathSource {
        source_id: 10,
        source_identity: ContentHash([251; 32]),
        surface_cell: CellId(11),
        area_m2: 0.001,
        stack: ExternalCoatingStack::Opaque {
            films: vec![
                layer(11, CoatingOpticsModel::Fixed(ir_film(1.5, 0.15))),
                layer(12, CoatingOpticsModel::Fixed(ir_film(2.1, 0.3))),
            ],
            substrate_receiver: receiver(13),
            substrate: SubstrateOpticsModel::Fixed(ir_metal()),
        },
        bath: CoatingExchangeBath::Blackbody {
            temperature_k: 300.0,
        },
        quadrature: CoatingExchangeQuadrature {
            wavelength_bins: 16,
            cosine_bins: 8,
        },
    }
}

// Independent integration of the original constitutive/field work and SI
// Planck law. No CoatingExternalBathResponse, emission radiance or native heat
// helper participates in this reference. The low-level optical tests supply
// the independent coherent-field witness for those original absorptances.
fn reference(
    films: &[CoherentAbsorbingCoating],
    substrate: Option<ConductorOptics>,
    reverse: bool,
    temperatures: &[f64],
    bath: Option<f64>,
    area: f64,
    q: CoatingExchangeQuadrature,
) -> Vec<f64> {
    let layers: Vec<_> = films.iter().map(CoatingLayer::Absorbing).collect();
    let mut rates = vec![0.0; films.len() + usize::from(substrate.is_some())];
    let width = 16000.0 / q.wavelength_bins as f64;
    for i in 0..q.wavelength_bins {
        let lambda = 4000.0 + (i as f64 + 0.5) * width;
        for j in 0..q.cosine_bins {
            let mu = (j as f64 + 0.5) / q.cosine_bins as f64;
            let absorbed: Vec<_> = if let Some(substrate) = substrate {
                let work = CoherentCoatingStack::try_new(&layers)
                    .unwrap()
                    .layer_power_balance(&substrate, lambda, 1.0, mu)
                    .unwrap();
                work.layer_absorbed()
                    .iter()
                    .copied()
                    .chain([work.substrate_absorbed])
                    .collect()
            } else {
                DielectricCoatingStack::try_mixed_in_band(&layers, (4000.0, 20000.0))
                    .unwrap()
                    .layer_power_balance(lambda, reverse, mu, 1.0, 1.0)
                    .unwrap()
                    .layer_absorbed()
                    .to_vec()
            };
            for (k, a) in absorbed.into_iter().enumerate() {
                rates[k] += area * std::f64::consts::PI * mu * width / q.cosine_bins as f64
                    * (a[0] + a[1])
                    * (planck(temperatures[k], lambda) - bath.map_or(0.0, |t| planck(t, lambda)));
            }
        }
    }
    rates
}

#[test]
fn independent_film_baths_reach_original_stores_and_charge_only_the_physical_surface() {
    let first = body(11, 450.0, 0.001);
    let mut buried = body(12, 250.0, 0.002);
    let mut substrate = body(13, 320.0, 0.003);
    // These complete native faces are occupied by zero-conductance reservoirs.
    // Exterior radiation still reaches the buried material through the films.
    for c in [&mut buried, &mut substrate] {
        c.boundary.reservoir = Some(Reservoir {
            temperature_k: 300.0,
            conductance_w_k: 0.0,
            area_m2: c.geometry.surface_area_m2,
        });
    }
    let native = thermal(vec![first, buried, substrate], &HOT, vec![]);
    let physical = opaque();
    let owner = prepare(&native, &[physical.clone()]);
    let expected = reference(
        &[ir_film(1.5, 0.15), ir_film(2.1, 0.3)],
        Some(ir_metal()),
        false,
        &[450.0, 250.0, 320.0],
        Some(300.0),
        physical.area_m2,
        physical.quadrature,
    );
    assert!(expected[0] > 0.0 && expected[1] < 0.0 && expected[2] > 0.0);
    for (i, watts) in expected.iter().enumerate() {
        near(owner.sources_w()[&CellId(11 + i as u64)], -*watts, 3e-12);
        assert!(owner.conductance_bounds_w_k()[&CellId(11 + i as u64)] > 0.0);
    }
    assert_eq!(
        owner.work().optical_evaluations,
        16 * 8,
        "one solve serves all channels"
    );
    let (after, observed, audit) = native
        .advance_observed(&owner, 0.02, config(0.005), || Ok(()))
        .unwrap();
    near(
        after.energy_j() - native.energy_j(),
        audit.external_heat_j,
        1e-11,
    );
    near(audit.external_heat_j, -observed.energy().outward_j, 1e-11);
    near(observed.energy().balance_residual_j, 0.0, 1e-13);
    for id in [11, 12, 13] {
        let a = &after.cells()[&CellId(id)];
        let b = &native.cells()[&CellId(id)];
        near(
            a.mass_kg * (a.specific_h_j_kg - b.specific_h_j_kg),
            observed.energy_by_cell_j()[&CellId(id)],
            1e-12,
        );
    }
    assert_eq!(
        observed.work().optical_evaluations,
        0,
        "changing temperatures reuse fixed optical coefficients"
    );
    assert_eq!(owner.time_s(), 0.0);

    // Repeated receivers are one real thermal store and recover the original
    // isothermal owner instead of creating three independent source accounts.
    let mut same = physical;
    if let ExternalCoatingStack::Opaque {
        films,
        substrate_receiver,
        ..
    } = &mut same.stack
    {
        for film in films {
            film.receiver = receiver(11);
        }
        *substrate_receiver = receiver(11);
    }
    let resolved = prepare(&native, &[same.clone()]);
    let old = exchange(
        &native,
        &[CoatingExchangeSource {
            source_id: same.source_id,
            source_identity: same.source_identity,
            receiver: receiver(11),
            area_m2: same.area_m2,
            bath: same.bath,
            quadrature: same.quadrature,
            stack: IsothermalCoatingStack::Opaque {
                films: vec![
                    CoatingOpticsModel::Fixed(ir_film(1.5, 0.15)),
                    CoatingOpticsModel::Fixed(ir_film(2.1, 0.3)),
                ],
                substrate: SubstrateOpticsModel::Fixed(ir_metal()),
            },
        }],
    );
    assert_eq!(resolved.sources_w().len(), 1);
    near(
        resolved.sources_w()[&CellId(11)],
        old.sources_w()[&CellId(11)],
        3e-13,
    );
    near(
        resolved.conductance_bounds_w_k()[&CellId(11)],
        old.conductance_bounds_w_k()[&CellId(11)],
        3e-13,
    );
}

fn changing_film(t: f64) -> CoherentAbsorbingCoating {
    ir_film(1.5, 0.08 + 0.002 * (t - 250.0))
}
fn transmitting_pair() -> [CoatingExternalBathSource; 2] {
    let law = Arc::new(
        SurfaceCoatingLaw::try_temperature(
            authority(),
            &[(250.0, changing_film(250.0)), (450.0, changing_film(450.0))],
            2,
        )
        .unwrap(),
    );
    std::array::from_fn(|side| CoatingExternalBathSource {
        source_id: 10 + side as u64,
        source_identity: ContentHash([252; 32]),
        surface_cell: CellId(11 + side as u64),
        area_m2: 0.001,
        stack: ExternalCoatingStack::Transmitting {
            films: vec![
                layer(11, CoatingOpticsModel::Law(Arc::clone(&law))),
                layer(12, CoatingOpticsModel::Fixed(ir_film(1.8, 0.0))),
            ],
            opposite_ior: vacuum(),
            reverse: side == 1,
        },
        bath: CoatingExchangeBath::Blackbody {
            temperature_k: if side == 0 { 420.0 } else { 250.0 },
        },
        quadrature: CoatingExchangeQuadrature {
            wavelength_bins: 8,
            cosine_bins: 4,
        },
    })
}

#[test]
fn complete_two_bath_single_absorber_preserves_reverse_mapping_and_refines_native_feedback() {
    let native = thermal(
        vec![body(11, 330.0, 0.001), body(12, 370.0, 0.002)],
        &HOT,
        vec![],
    );
    let ports = transmitting_pair();
    let mut budget = limits();
    budget.exchange.max_total_optical_knots = 2;
    let initial =
        NativeCoatingExternalBathExchange::try_new(&native, &ports, budget, || Ok(())).unwrap();
    let shuffled = NativeCoatingExternalBathExchange::try_new(
        &native,
        &[ports[1].clone(), ports[0].clone()],
        budget,
        || Ok(()),
    )
    .unwrap();
    assert_eq!(initial.identity(), shuffled.identity());
    let outward = |t: f64| -> f64 {
        ports
            .iter()
            .enumerate()
            .map(|(side, port)| {
                let bath = if side == 0 { 420.0 } else { 250.0 };
                reference(
                    &[changing_film(t), ir_film(1.8, 0.0)],
                    None,
                    side == 1,
                    &[t, 370.0],
                    Some(bath),
                    port.area_m2,
                    port.quadrature,
                )[0]
            })
            .sum()
    };
    near(initial.sources_w()[&CellId(11)], -outward(330.0), 3e-12);
    assert_eq!(initial.sources_w()[&CellId(12)], 0.0);
    assert!(initial.rates(10).unwrap().film_rates()[0].outward_power_w < 0.0);
    assert!(initial.rates(11).unwrap().film_rates()[0].outward_power_w > 0.0);
    for id in [10, 11] {
        assert_eq!(
            initial.rates(id).unwrap().film_rates()[1].outward_power_w,
            0.0
        );
    }
    let duration = 0.1;
    let step = duration / 128.0;
    let mut reference_t = 330.0;
    let derivative = |t| -outward(t) / 0.01;
    for _ in 0..128 {
        let a = derivative(reference_t);
        let b = derivative(reference_t + 0.5 * step * a);
        let c = derivative(reference_t + 0.5 * step * b);
        let d = derivative(reference_t + step * c);
        reference_t += step * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    let mut errors = Vec::new();
    for dt in [0.025, 0.0125, 0.00625] {
        let (after, observed, audit) = native
            .advance_observed(&initial, duration, config(dt), || Ok(()))
            .unwrap();
        let t = after.temperature_k(CellId(11)).unwrap();
        errors.push((t - reference_t).abs());
        assert_eq!(after.temperature_k(CellId(12)).unwrap(), 370.0);
        assert_eq!(observed.energy_by_cell_j()[&CellId(12)], 0.0);
        near(
            0.01 * (t - 330.0),
            observed.energy_by_cell_j()[&CellId(11)],
            1e-12,
        );
        near(
            audit.external_heat_j,
            observed.energy().thermal_external_j,
            1e-12,
        );
        near(observed.energy().balance_residual_j, 0.0, 1e-13);
        assert_ne!(
            observed.emitter(10).unwrap().optical_identity(),
            initial.emitter(10).unwrap().optical_identity()
        );
        assert!(observed.work().optical_evaluations >= 2 * 32 * 4);
        assert_eq!(observed.emitter(11).unwrap().film_temperature_k(0), Some(t));
    }
    assert!(
        errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1],
        "{errors:?}"
    );
    budget.exchange.max_total_optical_knots = 1;
    assert!(matches!(
        NativeCoatingExternalBathExchange::try_new(&native, &ports, budget, || Ok(())),
        Err(NativeCoatingExchangeError::Budget)
    ));
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
        Some([253; 32])
    }
}

#[test]
fn layer_enthalpy_and_surface_bath_events_refresh_on_the_same_native_clock() {
    let mut a = cell(11, 0.001, 1.0, 0.0);
    a.specific_h_j_kg = 70.0;
    let mut b = a.clone();
    b.id = CellId(12);
    b.specific_h_j_kg = 45.0;
    a.boundary.ambient_k = 250.0;
    let mut warm = a.boundary;
    warm.ambient_k = 350.0;
    let native = thermal(
        vec![a, b],
        &Latent,
        vec![TimedEvent {
            time_s: 0.02,
            event: Event::Boundary(CellId(11), warm),
        }],
    );
    let law = Arc::new(
        SurfaceCoatingLaw::try_enthalpy(
            authority(),
            ContentHash([253; 32]),
            (310.0, 310.0),
            &[(10.0, ir_film(1.5, 0.01)), (110.0, ir_film(1.5, 0.6))],
            2,
        )
        .unwrap(),
    );
    let mut films = vec![
        layer(11, CoatingOpticsModel::Law(Arc::clone(&law))),
        layer(12, CoatingOpticsModel::Law(law)),
    ];
    for film in &mut films {
        film.receiver.temperature_bounds_k = (310.0, 310.0);
    }
    let port = CoatingExternalBathSource {
        source_id: 1,
        source_identity: ContentHash([254; 32]),
        surface_cell: CellId(11),
        area_m2: 0.0001,
        stack: ExternalCoatingStack::Transmitting {
            films,
            opposite_ior: vacuum(),
            reverse: false,
        },
        bath: CoatingExchangeBath::NativeAmbient,
        quadrature: CoatingExchangeQuadrature {
            wavelength_bins: 8,
            cosine_bins: 4,
        },
    };
    let mut budget = limits();
    budget.exchange.max_total_optical_knots = 2;
    let initial =
        NativeCoatingExternalBathExchange::try_new(&native, &[port], budget, || Ok(())).unwrap();
    let (whole, observed, audit) = native
        .advance_observed(&initial, 0.04, config(0.005), || Ok(()))
        .unwrap();
    let (middle, partial, _) = native
        .advance_observed(&initial, 0.02, config(0.005), || Ok(()))
        .unwrap();
    let (split, resumed, _) = middle
        .advance_observed(&partial, 0.04, config(0.005), || Ok(()))
        .unwrap();
    assert!(initial.sources_w()[&CellId(11)] < 0.0 && partial.sources_w()[&CellId(11)] > 0.0);
    for id in [11, 12] {
        assert_eq!(whole.temperature_k(CellId(id)).unwrap(), 310.0);
        near(
            whole.cells()[&CellId(id)].specific_h_j_kg,
            split.cells()[&CellId(id)].specific_h_j_kg,
            1e-13,
        );
        near(
            observed.energy_by_cell_j()[&CellId(id)],
            resumed.energy_by_cell_j()[&CellId(id)],
            1e-13,
        );
        assert_ne!(
            whole.cells()[&CellId(id)].specific_h_j_kg,
            native.cells()[&CellId(id)].specific_h_j_kg
        );
    }
    assert_ne!(
        observed.emitter(1).unwrap().optical_identity(),
        initial.emitter(1).unwrap().optical_identity()
    );
    near(
        whole.energy_j() - native.energy_j(),
        audit.external_heat_j,
        1e-12,
    );
    near(audit.external_heat_j, -observed.energy().outward_j, 1e-12);
}

#[test]
fn exterior_channel_budgets_cancellation_and_surface_coverage_refuse_atomically() {
    let bodies = || {
        vec![
            body(11, 450.0, 0.001),
            body(12, 250.0, 0.002),
            body(13, 320.0, 0.003),
        ]
    };
    let native = thermal(bodies(), &HOT, vec![]);
    let port = opaque();
    let owner = prepare(&native, &[port.clone()]);
    let mut budget = limits();
    budget.max_channel_spectral_nodes = 3 * 16 - 1;
    assert!(matches!(
        NativeCoatingExternalBathExchange::try_new(&native, &[port.clone()], budget, || Ok(())),
        Err(NativeCoatingExchangeError::Budget)
    ));
    let mut duplicate = port.clone();
    duplicate.source_id = 11;
    budget.max_channel_spectral_nodes = 2 * 3 * 16 - 1;
    assert!(matches!(
        NativeCoatingExternalBathExchange::try_new(
            &native,
            &[port.clone(), duplicate],
            budget,
            || Ok(())
        ),
        Err(NativeCoatingExchangeError::Budget)
    ));
    assert!(
        NativeCoatingExternalBathExchange::try_new(
            &native,
            &[port.clone(), port.clone()],
            limits(),
            || Ok(())
        )
        .is_err()
    );
    budget = limits();
    budget.exchange.max_optical_evaluations = 16 * 8 - 1;
    assert!(matches!(
        NativeCoatingExternalBathExchange::try_new(&native, &[port.clone()], budget, || Ok(())),
        Err(NativeCoatingExchangeError::Budget)
    ));
    budget = limits();
    budget.exchange.max_source_evaluations = 1;
    let limited =
        NativeCoatingExternalBathExchange::try_new(&native, &[port.clone()], budget, || Ok(()))
            .unwrap();
    assert!(matches!(
        native.advance_observed(&limited, 0.01, config(0.01), || Ok(())),
        Err(NativeCoatingExchangeError::Budget)
    ));
    let mut calls = 0;
    assert!(matches!(
        native.advance_observed(&owner, 0.1, config(0.01), || {
            calls += 1;
            if calls == 150 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(NativeCoatingExchangeError::Thermal(
            DepositionError::Cancelled
        ))
    ));
    assert_eq!(owner.time_s(), 0.0);
    assert_eq!(owner.energy().duration_s, 0.0);
    assert_eq!(native.time_s(), 0.0);
    let accepted = native
        .advance_observed(&owner, 0.02, config(0.01), || Ok(()))
        .unwrap();
    let repeated = native
        .advance_observed(&owner, 0.02, config(0.01), || Ok(()))
        .unwrap();
    assert_eq!(accepted.1.energy(), repeated.1.energy());
    let covered = thermal(
        bodies(),
        &HOT,
        vec![TimedEvent {
            time_s: 0.01,
            event: Event::Contact(ThermalLink {
                id: 9,
                a: CellId(11),
                b: CellId(12),
                area_m2: 6.0,
                conductance_w_k: 0.0,
            }),
        }],
    );
    let before = prepare(&covered, &[port]);
    assert!(
        covered
            .advance_observed(&before, 0.02, config(0.005), || Ok(()))
            .is_err()
    );
    assert_eq!(before.time_s(), 0.0);
    assert_eq!(covered.time_s(), 0.0);
}
