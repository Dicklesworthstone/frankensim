use super::*;
use fs_render::thermal_surface::heating::coating_exchange::gap::*;
use fs_render::thermal_surface::heating::feedback::SurfaceConductorLaw;

fn ir_metal(eta: f64, k: f64) -> ConductorOptics {
    ConductorOptics::try_new_in_band(
        std::array::from_fn(|i| {
            ConductorIorSample::try_new_in_band(4000.0 + 2000.0 * i as f64, eta, k).unwrap()
        }),
        2000.0,
        authority(),
    )
    .unwrap()
}
fn ir_film(k: f64) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new_in_band(
        std::array::from_fn(|i| {
            CoatingIorSample::try_new_in_band(4000.0 + 2000.0 * i as f64, 1.5, k).unwrap()
        }),
        2e-6,
        2000.0,
        authority(),
    )
    .unwrap()
}
fn body(id: u64, t: f64, mass: f64) -> DepositedCell {
    let mut c = specimen(t, mass);
    c.id = CellId(id);
    c
}
fn side(id: u64) -> CoatingGapSurface {
    CoatingGapSurface {
        source_identity: ContentHash([id as u8; 32]),
        receiver: SurfaceHeatingReceiver {
            thermal_cell: CellId(id),
            material_id: MaterialId(7),
            temperature_bounds_k: (100.0, 5000.0),
        },
        stack: IsothermalCoatingStack::Opaque {
            films: vec![],
            substrate: SubstrateOpticsModel::Fixed(ir_metal(1.0, 1e-20)),
        },
    }
}
fn link(id: u64, a: u64, b: u64) -> CoatingGapLink {
    CoatingGapLink {
        link_id: id,
        surfaces: [side(a), side(b)],
        area_m2: 0.001,
        quadrature: CoatingExchangeQuadrature {
            wavelength_bins: 32,
            cosine_bins: 4,
        },
    }
}
fn limits() -> NativeCoatingGapBudget {
    NativeCoatingGapBudget {
        max_links: 8,
        max_total_optical_knots: 16,
        max_source_evaluations: 10_000,
        max_optical_evaluations: 1_000_000,
    }
}
fn prepare(
    native: &DepositionSimulation<'_>,
    links: &[CoatingGapLink],
) -> NativeCoatingGapExchange {
    NativeCoatingGapExchange::try_new(native, links, limits(), || Ok(())).unwrap()
}
fn black_power(a: f64, b: f64) -> f64 {
    (0..32)
        .map(|i| {
            let wavelength = 4000.0 + (i as f64 + 0.5) * 500.0;
            std::f64::consts::PI * 0.001 * 500.0 * (planck(a, wavelength) - planck(b, wavelength))
        })
        .sum()
}

#[test]
fn native_gap_is_internal_conservative_and_uses_both_original_caloric_stores() {
    for (first, second) in [(350.0, 250.0), (250.0, 350.0), (300.0, 300.0)] {
        let native = thermal(
            vec![body(11, first, 0.001), body(12, second, 0.002)],
            &HOT,
            vec![],
        );
        let owner = prepare(&native, &[link(1, 11, 12)]);
        near(
            owner.heat_transfers()[0].power_w,
            black_power(first, second),
            3e-12,
        );
        assert_eq!(owner.work().optical_evaluations, 2 * 32 * 4);
        let power = native.power_balance_observed(&owner, || Ok(())).unwrap();
        assert_eq!(power.external_power_w, 0.0);
        near(
            power.net_power_w[&CellId(11)],
            -owner.heat_transfers()[0].power_w,
            1e-14,
        );
        near(
            power.net_power_w[&CellId(12)],
            owner.heat_transfers()[0].power_w,
            1e-14,
        );
        let (after, observed, audit) = native
            .advance_observed(&owner, 5.0, config(1e6), || Ok(()))
            .unwrap();
        assert_eq!(audit.external_heat_j, 0.0);
        near(after.energy_j(), native.energy_j(), 1e-12);
        near(observed.energy().balance_residual_j, 0.0, 1e-13);
        let exchanged = observed.energy().by_link_first_to_second_j[&1];
        near(
            10.0 * 0.001 * (after.temperature_k(CellId(11)).unwrap() - first),
            -exchanged,
            1e-11,
        );
        near(
            10.0 * 0.002 * (after.temperature_k(CellId(12)).unwrap() - second),
            exchanged,
            1e-11,
        );
        for id in [11, 12] {
            let t = after.temperature_k(CellId(id)).unwrap();
            assert!(t >= first.min(second) - 1e-10 && t <= first.max(second) + 1e-10);
            assert_eq!(
                observed
                    .emitter(1, usize::from(id == 12))
                    .unwrap()
                    .substrate_temperature_k(),
                Some(t)
            );
        }
        assert_eq!(
            observed.work().optical_evaluations,
            0,
            "fixed optics reuse the joint response"
        );
        if first == second {
            assert_eq!(exchanged, 0.0);
        }
        assert_eq!(owner.time_s(), 0.0);
        assert_eq!(native.time_s(), 0.0);
    }
}

fn temperature_film(t: f64) -> CoherentAbsorbingCoating {
    ir_film(0.05 + 0.0006 * (t - 200.0))
}
fn temperature_metal(t: f64) -> ConductorOptics {
    ir_metal(0.5 + 0.001 * (t - 200.0), 2.5 + 0.002 * (t - 200.0))
}
fn dynamic_link() -> CoatingGapLink {
    let mut pair = link(1, 11, 12);
    pair.area_m2 = 0.01;
    pair.quadrature.wavelength_bins = 8;
    pair.surfaces[0].stack = IsothermalCoatingStack::Opaque {
        films: vec![CoatingOpticsModel::Law(Arc::new(
            SurfaceCoatingLaw::try_temperature(
                authority(),
                &[
                    (200.0, temperature_film(200.0)),
                    (500.0, temperature_film(500.0)),
                ],
                2,
            )
            .unwrap(),
        ))],
        substrate: SubstrateOpticsModel::Fixed(ir_metal(0.5, 3.0)),
    };
    pair.surfaces[1].stack = IsothermalCoatingStack::Opaque {
        films: vec![],
        substrate: SubstrateOpticsModel::Law(Arc::new(
            SurfaceConductorLaw::try_temperature(
                authority(),
                &[
                    (200.0, temperature_metal(200.0)),
                    (500.0, temperature_metal(500.0)),
                ],
                2,
            )
            .unwrap(),
        )),
    };
    pair
}
// Independent original coherent fields + direct SI Planck and a geometric
// reflection-series resistance, evaluated at the declared quadrature nodes.
fn pair_power(first: f64, second: f64) -> f64 {
    let film = temperature_film(first);
    let layers = [CoatingLayer::Absorbing(&film)];
    let first_stack = CoherentCoatingStack::try_new(&layers).unwrap();
    let second_stack = CoherentCoatingStack::try_new(&[]).unwrap();
    let mut power = 0.0;
    for w in 0..8 {
        let wavelength = 4000.0 + (w as f64 + 0.5) * 2000.0;
        for a in 0..4 {
            let mu = (a as f64 + 0.5) / 4.0;
            let left = first_stack
                .layer_power_balance(&ir_metal(0.5, 3.0), wavelength, 1.0, mu)
                .unwrap()
                .total
                .absorbed;
            let right = second_stack
                .layer_power_balance(&temperature_metal(second), wavelength, 1.0, mu)
                .unwrap()
                .total
                .absorbed;
            for pol in 0..2 {
                let exchange = 1.0 / (1.0 / left[pol] + 1.0 / right[pol] - 1.0);
                power += std::f64::consts::PI * 0.01 * 2000.0 / 4.0
                    * mu
                    * exchange
                    * (planck(first, wavelength) - planck(second, wavelength));
            }
        }
    }
    power
}

#[test]
fn independent_optical_states_and_native_conduction_refine_to_a_two_body_rk4_reference() {
    let mut first = body(11, 340.0, 0.001);
    first.boundary.source_w = 0.02;
    let native = thermal(
        vec![first, body(12, 280.0, 0.002)],
        &HOT,
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 99,
                a: CellId(11),
                b: CellId(12),
                area_m2: 0.01,
                conductance_w_k: 0.001,
            }),
        }],
    );
    let owner = prepare(&native, &[dynamic_link()]);
    near(
        owner.heat_transfers()[0].power_w,
        pair_power(340.0, 280.0),
        2e-11,
    );
    let duration = 0.04;
    let step = duration / 256.0;
    let derivative = |x: [f64; 3]| {
        let radiation = pair_power(x[0], x[1]);
        let conduction = 0.001 * (x[0] - x[1]);
        [
            (0.02 - radiation - conduction) / 0.01,
            (radiation + conduction) / 0.02,
            radiation,
        ]
    };
    let mut reference = [340.0, 280.0, 0.0];
    for _ in 0..256 {
        let a = derivative(reference);
        let b = derivative(std::array::from_fn(|i| reference[i] + step * 0.5 * a[i]));
        let c = derivative(std::array::from_fn(|i| reference[i] + step * 0.5 * b[i]));
        let d = derivative(std::array::from_fn(|i| reference[i] + step * c[i]));
        for i in 0..3 {
            reference[i] += step * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    let mut errors = Vec::new();
    for dt in [0.01, 0.005, 0.0025] {
        let (after, observed, audit) = native
            .advance_observed(&owner, duration, config(dt), || Ok(()))
            .unwrap();
        let error = (after.temperature_k(CellId(11)).unwrap() - reference[0]).abs()
            + (after.temperature_k(CellId(12)).unwrap() - reference[1]).abs()
            + (observed.energy().by_link_first_to_second_j[&1] - reference[2]).abs();
        errors.push(error);
        near(audit.external_heat_j, 0.02 * duration, 1e-13);
        near(
            after.energy_j() - native.energy_j(),
            audit.external_heat_j,
            1e-11,
        );
        assert!(observed.work().optical_evaluations > 2 * 8 * 4);
    }
    assert!(
        errors[1] < 0.65 * errors[0] && errors[2] < 0.65 * errors[1],
        "{errors:?}"
    );
}

struct GapLatent;
impl ThermalLaw for GapLatent {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(0.0..=1000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if h < 10.0 {
            300.0 + h
        } else if h <= 110.0 {
            310.0
        } else {
            310.0 + h - 110.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([244; 32])
    }
}

#[test]
fn gap_network_retains_latent_optics_canonical_order_events_and_continuation() {
    let mut hot = cell(11, 0.001, 1.0, 0.0);
    hot.specific_h_j_kg = 200.0;
    let mut latent = cell(12, 0.001, 1.0, 0.0);
    latent.specific_h_j_kg = 70.0;
    let mut later = latent.boundary;
    later.source_w = 0.01;
    let native = thermal(
        vec![hot, latent, cell(13, 0.001, 1.0, 0.0)],
        &GapLatent,
        vec![TimedEvent {
            time_s: 0.02,
            event: Event::Boundary(CellId(12), later),
        }],
    );
    let law = Arc::new(
        SurfaceCoatingLaw::try_enthalpy(
            authority(),
            ContentHash([244; 32]),
            (310.0, 310.0),
            &[(10.0, ir_film(0.01)), (110.0, ir_film(0.6))],
            2,
        )
        .unwrap(),
    );
    let mut a = link(1, 11, 12);
    let mut b = link(2, 12, 13);
    for s in [&mut a.surfaces[1], &mut b.surfaces[0]] {
        s.receiver.temperature_bounds_k = (310.0, 310.0);
        s.stack = IsothermalCoatingStack::Opaque {
            films: vec![CoatingOpticsModel::Law(Arc::clone(&law))],
            substrate: SubstrateOpticsModel::Fixed(ir_metal(0.5, 3.0)),
        };
    }
    let mut budget = limits();
    budget.max_total_optical_knots = 2;
    let ordered =
        NativeCoatingGapExchange::try_new(&native, &[a.clone(), b.clone()], budget, || Ok(()))
            .unwrap();
    let shuffled =
        NativeCoatingGapExchange::try_new(&native, &[b.clone(), a.clone()], budget, || Ok(()))
            .unwrap();
    assert_eq!(ordered.identity(), shuffled.identity());
    let old_optics = ordered.emitter(1, 1).unwrap().optical_identity();
    let (whole, observed, audit) = native
        .advance_observed(&ordered, 0.04, config(0.005), || Ok(()))
        .unwrap();
    let (middle, partial, _) = native
        .advance_observed(&shuffled, 0.02, config(0.005), || Ok(()))
        .unwrap();
    let (split, resumed, _) = middle
        .advance_observed(&partial, 0.04, config(0.005), || Ok(()))
        .unwrap();
    for id in [11, 12, 13] {
        near(
            whole.cells()[&CellId(id)].specific_h_j_kg,
            split.cells()[&CellId(id)].specific_h_j_kg,
            1e-12,
        );
    }
    assert_eq!(whole.temperature_k(CellId(12)).unwrap(), 310.0);
    assert_ne!(
        observed.emitter(1, 1).unwrap().optical_identity(),
        old_optics
    );
    for id in [1, 2] {
        near(
            observed.energy().by_link_first_to_second_j[&id],
            resumed.energy().by_link_first_to_second_j[&id],
            1e-12,
        );
    }
    near(audit.external_heat_j, 0.01 * 0.02, 1e-13);
    near(
        whole.energy_j() - native.energy_j(),
        audit.external_heat_j,
        1e-12,
    );
    near(observed.energy().balance_residual_j, 0.0, 1e-13);
    budget.max_total_optical_knots = 1;
    assert!(matches!(
        NativeCoatingGapExchange::try_new(&native, &[a, b], budget, || Ok(())),
        Err(NativeCoatingExchangeError::Budget)
    ));
}

#[test]
fn gap_domains_coverage_work_and_late_cancellation_are_transactional() {
    let native = thermal(
        vec![body(11, 350.0, 0.001), body(12, 250.0, 0.002)],
        &HOT,
        vec![],
    );
    let valid = link(1, 11, 12);
    let owner = prepare(&native, &[valid.clone()]);
    for bad in [link(1, 11, 11), {
        let mut x = valid.clone();
        x.area_m2 = 6.0001;
        x
    }] {
        assert!(NativeCoatingGapExchange::try_new(&native, &[bad], limits(), || Ok(())).is_err());
    }
    assert!(
        NativeCoatingGapExchange::try_new(
            &native,
            &[valid.clone(), valid.clone()],
            limits(),
            || Ok(())
        )
        .is_err()
    );
    let mut budget = limits();
    budget.max_optical_evaluations = 2 * 32 * 4 - 1;
    assert!(matches!(
        NativeCoatingGapExchange::try_new(&native, &[valid.clone()], budget, || Ok(())),
        Err(NativeCoatingExchangeError::Budget)
    ));
    budget = limits();
    budget.max_source_evaluations = 1;
    let limited =
        NativeCoatingGapExchange::try_new(&native, &[valid.clone()], budget, || Ok(())).unwrap();
    assert!(matches!(
        native.advance_observed(&limited, 0.01, config(0.01), || Ok(())),
        Err(NativeCoatingExchangeError::Budget)
    ));
    let births = [body(11, 350.0, 0.001), body(12, 250.0, 0.002)]
        .into_iter()
        .map(|c| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(c),
        })
        .collect();
    let no_transfers = DepositionSimulation::new(
        births,
        &[(MaterialId(7), &HOT as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 2,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(0.0, config(0.01), || Ok(()))
    .unwrap()
    .0;
    let source = prepare(&no_transfers, &[valid.clone()]);
    assert!(matches!(
        no_transfers.advance_observed(&source, 0.01, config(0.01), || Ok(())),
        Err(NativeCoatingExchangeError::Thermal(
            DepositionError::Budget(_)
        ))
    ));
    let mut calls = 0;
    assert!(matches!(
        native.advance_observed(&owner, 0.1, config(0.01), || {
            calls += 1;
            if calls == 100 {
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
    assert_eq!(native.time_s(), 0.0);
    assert_eq!(owner.energy().duration_s, 0.0);
    let a = native
        .advance_observed(&owner, 0.02, config(0.01), || Ok(()))
        .unwrap();
    let b = native
        .advance_observed(&owner, 0.02, config(0.01), || Ok(()))
        .unwrap();
    assert_eq!(a.1.energy(), b.1.energy());
    assert_eq!(a.0.energy_j(), b.0.energy_j());
    let contacted = thermal(
        vec![body(11, 350.0, 0.001), body(12, 250.0, 0.002)],
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
    let pending = prepare(&contacted, &[valid]);
    assert!(matches!(
        contacted.advance_observed(&pending, 0.02, config(0.005), || Ok(())),
        Err(NativeCoatingExchangeError::Invalid(_))
    ));
    assert_eq!(contacted.time_s(), 0.0);
    assert_eq!(pending.time_s(), 0.0);
}
