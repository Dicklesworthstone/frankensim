//! G1/G3/G4: independent multi-gap light balances drive original thermal stores.
use super::transmitting_pair::{ir_film, vacuum};
use super::*;
use fs_render::dielectric::DielectricCoatingStack;
use fs_render::thermal_surface::heating::coating_exchange::radiative_stack::*;

const AREA: f64 = 0.003;
const BATHS: [f64; 2] = [430.0, 240.0];
const REVERSE: [bool; 3] = [false, true, false];
const INITIAL: [f64; 3] = [340.0, 280.0, 375.0];
const MASSES: [f64; 3] = [0.001, 0.0016, 0.0022];

fn evolving_film(body: usize, t: f64) -> CoherentAbsorbingCoating {
    let i = body as f64;
    ir_film(
        1.3 + 0.35 * i + 0.0003 * (t - 250.0),
        0.04 + 0.06 * i + (0.0004 + i * 0.0001) * (t - 250.0),
        (1.2 + 0.3 * i) * 1e-6,
    )
}

fn outer_film(body: usize) -> CoherentAbsorbingCoating {
    ir_film(2.3 - 0.4 * body as f64, 0.025 + 0.05 * body as f64, 2.1e-6)
}

fn stack(dynamic: bool) -> CoatingRadiativeStack {
    CoatingRadiativeStack {
        stack_id: 31,
        bodies: (0..3)
            .map(|body| CoatingRadiativeBody {
                source_identity: ContentHash([210 + body as u8; 32]),
                receiver: SurfaceHeatingReceiver {
                    thermal_cell: CellId(11 + body as u64),
                    material_id: MaterialId(7),
                    temperature_bounds_k: (250.0, 450.0),
                },
                stack: IsothermalCoatingStack::Transmitting {
                    films: vec![
                        if dynamic {
                            CoatingOpticsModel::Law(Arc::new(
                                SurfaceCoatingLaw::try_temperature(
                                    authority(),
                                    &[
                                        (250.0, evolving_film(body, 250.0)),
                                        (450.0, evolving_film(body, 450.0)),
                                    ],
                                    2,
                                )
                                .unwrap(),
                            ))
                        } else {
                            CoatingOpticsModel::Fixed(evolving_film(body, INITIAL[body]))
                        },
                        CoatingOpticsModel::Fixed(outer_film(body)),
                    ],
                    opposite_ior: vacuum(),
                    reverse: REVERSE[body],
                },
            })
            .collect(),
        area_m2: AREA,
        baths: BATHS.map(|temperature_k| CoatingExchangeBath::Blackbody { temperature_k }),
        quadrature: CoatingExchangeQuadrature {
            wavelength_bins: 8,
            cosine_bins: 4,
        },
    }
}

fn cells(masses: [f64; 3]) -> Vec<DepositedCell> {
    INITIAL
        .into_iter()
        .enumerate()
        .map(|(body, t)| {
            let mut c = specimen(t, masses[body]);
            c.id = CellId(11 + body as u64);
            c.boundary.ambient_k = [260.0, 400.0, 290.0][body];
            c
        })
        .collect()
}

fn limits() -> NativeCoatingRadiativeStackBudget {
    NativeCoatingRadiativeStackBudget {
        max_stacks: 3,
        max_total_bodies: 9,
        max_total_optical_knots: 6,
        max_source_evaluations: 10_000,
        max_optical_evaluations: 1_000_000,
        max_channel_spectral_nodes: 3 * 10 * 8,
    }
}

fn prepare(
    native: &DepositionSimulation<'_>,
    stack: &CoatingRadiativeStack,
) -> NativeCoatingRadiativeStack {
    NativeCoatingRadiativeStack::try_new(native, &[stack.clone()], limits(), || Ok(())).unwrap()
}

// Four independent unknowns are the two directional intensities in each gap.
// This deliberately solves the unreduced transport equations by pivoted dense
// elimination; no production exchange coefficient or elimination is reused.
fn solve(mut a: [[f64; 5]; 4]) -> [f64; 4] {
    for column in 0..4 {
        let pivot = (column..4)
            .max_by(|&i, &j| a[i][column].abs().total_cmp(&a[j][column].abs()))
            .unwrap();
        a.swap(column, pivot);
        assert!(a[column][column].abs() > 1e-10);
        let diagonal = a[column][column];
        for j in column..5 {
            a[column][j] /= diagonal;
        }
        for row in 0..4 {
            if row != column {
                let factor = a[row][column];
                for j in column..5 {
                    a[row][j] -= factor * a[column][j];
                }
            }
        }
    }
    std::array::from_fn(|i| a[i][4])
}

fn absorbed_power(temperatures: [f64; 3]) -> [f64; 5] {
    let films =
        std::array::from_fn::<_, 3, _>(|i| [evolving_film(i, temperatures[i]), outer_film(i)]);
    let layers = films
        .each_ref()
        .map(|f| f.each_ref().map(CoatingLayer::Absorbing));
    let optical = layers
        .each_ref()
        .map(|l| DielectricCoatingStack::try_mixed_in_band(l, (4000.0, 20000.0)).unwrap());
    let mut watts = [0.0; 5];
    for w in 0..8 {
        let lambda = 4000.0 + (w as f64 + 0.5) * 2000.0;
        let body = temperatures.map(|t| planck(t, lambda));
        let bath = BATHS.map(|t| planck(t, lambda));
        for angular in 0..4 {
            let mu = (angular as f64 + 0.5) / 4.0;
            let ports = std::array::from_fn::<_, 3, _>(|i| {
                std::array::from_fn::<_, 2, _>(|side| {
                    optical[i]
                        .layer_power_balance(lambda, REVERSE[i] ^ (side == 1), mu, 1.0, 1.0)
                        .unwrap()
                })
            });
            for polarization in 0..2 {
                let sign = if polarization == 0 { 1.0 } else { -1.0 };
                let r = ports.each_ref().map(|p| {
                    p.each_ref()
                        .map(|w| w.fresnel.reflection.0[0][0] + sign * w.fresnel.reflection.0[0][1])
                });
                let t = ports.each_ref().map(|p| {
                    p.each_ref().map(|w| {
                        w.fresnel.transmission.0[0][0] + sign * w.fresnel.transmission.0[0][1]
                    })
                });
                let a = ports
                    .each_ref()
                    .map(|p| p.each_ref().map(|w| w.absorbed[polarization]));
                // x0 = TL0*BL + RR0*y0 + AR0*B0
                // y0 = RL1*x0 + TR1*y1 + AL1*B1
                // x1 = TL1*x0 + RR1*y1 + AR1*B1
                // y1 = RL2*x1 + TR2*BR + AL2*B2
                let gap = solve([
                    [
                        1.0,
                        -r[0][1],
                        0.0,
                        0.0,
                        t[0][0] * bath[0] + a[0][1] * body[0],
                    ],
                    [-r[1][0], 1.0, 0.0, -t[1][1], a[1][0] * body[1]],
                    [-t[1][0], 0.0, 1.0, -r[1][1], a[1][1] * body[1]],
                    [
                        0.0,
                        0.0,
                        -r[2][0],
                        1.0,
                        t[2][1] * bath[1] + a[2][0] * body[2],
                    ],
                ]);
                let incoming = [[bath[0], gap[1]], [gap[0], gap[3]], [gap[2], bath[1]]];
                let body_heat = std::array::from_fn::<_, 3, _>(|i| {
                    a[i][0] * incoming[i][0] + a[i][1] * incoming[i][1]
                        - (a[i][0] + a[i][1]) * body[i]
                });
                let bath_heat = [
                    r[0][0] * bath[0] + t[0][1] * gap[1] + a[0][0] * body[0] - bath[0],
                    t[2][0] * gap[2] + r[2][1] * bath[1] + a[2][1] * body[2] - bath[1],
                ];
                // SI Planck includes both polarizations; /2 times azimuthal
                // 2*pi leaves this pi-weight for each independently solved mode.
                let measure = std::f64::consts::PI * AREA * mu * 2000.0 / 4.0;
                for (i, heat) in body_heat.into_iter().chain(bath_heat).enumerate() {
                    watts[i] += measure * heat;
                }
            }
        }
    }
    watts
}

#[test]
fn three_body_native_radiation_matches_independent_intensities_and_thermal_refinement() {
    let mut cells = cells(MASSES);
    cells[0].boundary.source_w = 0.025;
    let native = thermal(
        cells,
        &HOT,
        (0..2)
            .map(|i| TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 70 + i,
                    a: CellId(11 + i),
                    b: CellId(12 + i),
                    area_m2: 0.01,
                    conductance_w_k: [0.001, 0.0007][i as usize],
                }),
            })
            .collect(),
    );
    let initial = prepare(&native, &stack(true));
    let heat = absorbed_power(INITIAL);
    near(heat.iter().sum(), 0.0, 3e-12);
    let conduction = [
        0.001 * (INITIAL[0] - INITIAL[1]),
        0.0007 * (INITIAL[1] - INITIAL[2]),
    ];
    let expected = [
        heat[0] + 0.025 - conduction[0],
        heat[1] + conduction[0] - conduction[1],
        heat[2] + conduction[1],
    ];
    let power = native.power_balance_observed(&initial, || Ok(())).unwrap();
    for (i, expected) in expected.into_iter().enumerate() {
        near(power.net_power_w[&CellId(11 + i as u64)], expected, 5e-12);
    }
    near(power.external_power_w, -heat[3] - heat[4] + 0.025, 5e-12);
    assert_eq!(initial.work().optical_evaluations, 6 * 8 * 4);
    assert_eq!(initial.heat_transfers().len(), 3);
    let distant = initial
        .rates(31)
        .unwrap()
        .body_pairs()
        .iter()
        .find(|edge| (edge.first, edge.second) == (0, 2))
        .unwrap();
    assert!(distant.exchange.secant_w_k > 0.0);

    let duration = 0.04;
    let step = duration / 64.0;
    let derivative = |x: [f64; 4]| {
        let heat = absorbed_power([x[0], x[1], x[2]]);
        let q01 = 0.001 * (x[0] - x[1]);
        let q12 = 0.0007 * (x[1] - x[2]);
        [
            (heat[0] + 0.025 - q01) / (10.0 * MASSES[0]),
            (heat[1] + q01 - q12) / (10.0 * MASSES[1]),
            (heat[2] + q12) / (10.0 * MASSES[2]),
            -heat[3] - heat[4] + 0.025,
        ]
    };
    let mut reference = [INITIAL[0], INITIAL[1], INITIAL[2], 0.0];
    for _ in 0..64 {
        let a = derivative(reference);
        let b = derivative(std::array::from_fn(|i| reference[i] + step * 0.5 * a[i]));
        let c = derivative(std::array::from_fn(|i| reference[i] + step * 0.5 * b[i]));
        let d = derivative(std::array::from_fn(|i| reference[i] + step * c[i]));
        for i in 0..4 {
            reference[i] += step * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    let mut errors = Vec::new();
    for dt in [0.01, 0.005, 0.0025] {
        let (after, observed, audit) = native
            .advance_observed(&initial, duration, config(dt), || Ok(()))
            .unwrap();
        errors.push(
            (0..3)
                .map(|i| (after.temperature_k(CellId(11 + i as u64)).unwrap() - reference[i]).abs())
                .sum::<f64>()
                + (audit.external_heat_j - reference[3]).abs(),
        );
        near(
            after.energy_j() - native.energy_j(),
            audit.external_heat_j,
            3e-11,
        );
        near(
            observed.energy().thermal_external_j + 0.025 * duration,
            audit.external_heat_j,
            3e-12,
        );
        near(observed.energy().balance_residual_j, 0.0, 2e-13);
        assert_eq!(observed.energy().by_stack_pair_j.len(), 3);
        assert!(observed.energy().by_stack_bath_transfer_j[&31] > 0.0);
        assert!(observed.work().optical_evaluations > 6 * 8 * 4);
        for i in 0..3 {
            assert_eq!(
                observed.emitter(31, i).unwrap().film_temperature_k(0),
                Some(after.temperature_k(CellId(11 + i as u64)).unwrap())
            );
        }
    }
    assert!(
        errors[1] < 0.65 * errors[0] && errors[2] < 0.65 * errors[1],
        "{errors:?}"
    );
}

#[test]
fn transparent_intermediate_body_passes_heat_and_opaque_shield_blocks_distant_channels() {
    let native = thermal(cells([1.0; 3]), &HOT, vec![]);
    let mut s = stack(false);
    if let IsothermalCoatingStack::Transmitting { films, .. } = &mut s.bodies[1].stack {
        films.clear();
    }
    let transparent = prepare(&native, &s);
    assert_eq!(transparent.sources_w()[&CellId(12)], 0.0);
    assert_eq!(
        transparent.rates(31).unwrap().body_outward_power_w()[1],
        0.0
    );
    let mut two = s.clone();
    two.bodies.remove(1);
    let direct = prepare(&native, &two);
    let pass = native
        .power_balance_observed(&transparent, || Ok(()))
        .unwrap();
    let absent = native.power_balance_observed(&direct, || Ok(())).unwrap();
    for id in [CellId(11), CellId(12), CellId(13)] {
        near(pass.net_power_w[&id], absent.net_power_w[&id], 3e-13);
    }
    let (after, observed, _) = native
        .advance_observed(&transparent, 0.01, config(0.01), || Ok(()))
        .unwrap();
    assert_eq!(after.cells()[&CellId(12)], native.cells()[&CellId(12)]);
    assert_eq!(observed.energy().by_cell_external_j[&CellId(12)], 0.0);
    if let IsothermalCoatingStack::Transmitting { films, .. } = &mut s.bodies[1].stack {
        *films = vec![CoatingOpticsModel::Fixed(ir_film(1.8, 30.0, 1e-4))];
    }
    let shielded = prepare(&native, &s);
    let rates = shielded.rates(31).unwrap();
    let distant = rates
        .body_pairs()
        .iter()
        .find(|e| (e.first, e.second) == (0, 2))
        .unwrap();
    assert_eq!(distant.exchange.outward_power_w, 0.0);
    assert_eq!(distant.exchange.secant_w_k, 0.0);
    assert_eq!(rates.body_baths()[0][1].outward_power_w, 0.0);
    assert_eq!(rates.body_baths()[2][0].outward_power_w, 0.0);
    assert_eq!(rates.bath_transfer.outward_power_w, 0.0);
    assert!(
        rates
            .body_pairs()
            .iter()
            .filter(|e| e.exchange.secant_w_k > 0.0)
            .count()
            == 2
    );
}

#[test]
fn endpoint_baths_events_and_continuation_reuse_fixed_optics_without_duplicate_heat() {
    let cells = cells([1.0; 3]);
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
    let mut s = stack(false);
    s.baths = [CoatingExchangeBath::NativeAmbient; 2];
    let initial = prepare(&native, &s);
    let optics: Vec<_> = (0..3)
        .map(|i| initial.emitter(31, i).unwrap().optical_identity())
        .collect();
    let (whole, observed, audit) = native
        .advance_observed(&initial, 0.1, config(0.025), || Ok(()))
        .unwrap();
    let (middle, partial, _) = native
        .advance_observed(&initial, 0.05, config(0.025), || Ok(()))
        .unwrap();
    let (split, resumed, _) = middle
        .advance_observed(&partial, 0.1, config(0.025), || Ok(()))
        .unwrap();
    for i in 0..3 {
        assert!(initial.rates(31).unwrap().body_baths()[i][0].outward_power_w > 0.0);
        assert!(partial.rates(31).unwrap().body_baths()[i][0].outward_power_w < 0.0);
        assert_eq!(
            observed.emitter(31, i).unwrap().optical_identity(),
            optics[i]
        );
        let id = CellId(11 + i as u64);
        near(
            10.0 * (whole.temperature_k(id).unwrap() - native.temperature_k(id).unwrap()),
            observed.energy().by_cell_internal_j[&id] + observed.energy().by_cell_external_j[&id],
            3e-11,
        );
    }
    assert_eq!(observed.work().optical_evaluations, 0);
    near(whole.energy_j(), split.energy_j(), 3e-13);
    near(
        observed.energy().thermal_external_j,
        audit.external_heat_j,
        3e-13,
    );
    near(
        observed.energy().thermal_external_j,
        resumed.energy().thermal_external_j,
        3e-13,
    );
    near(observed.energy().balance_residual_j, 0.0, 2e-13);
    let mut a = s.clone();
    a.area_m2 *= 0.5;
    let mut b = a.clone();
    b.stack_id = 32;
    let ordered =
        NativeCoatingRadiativeStack::try_new(&native, &[a.clone(), b.clone()], limits(), || Ok(()))
            .unwrap();
    let reversed =
        NativeCoatingRadiativeStack::try_new(&native, &[b, a], limits(), || Ok(())).unwrap();
    assert_eq!(ordered.identity(), reversed.identity());
    for (&id, &watts) in initial.sources_w() {
        near(ordered.sources_w()[&id], watts, 2e-14);
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
        Some([209; 32])
    }
}

#[test]
fn native_enthalpy_changes_refresh_all_shield_optics_on_a_latent_plateau() {
    let mut cells = cells([0.001; 3]);
    for (i, cell) in cells.iter_mut().enumerate() {
        cell.specific_h_j_kg = 80.0 - 20.0 * i as f64;
    }
    let native = thermal(cells, &Latent, vec![]);
    let law = Arc::new(
        SurfaceCoatingLaw::try_enthalpy(
            authority(),
            ContentHash([209; 32]),
            (310.0, 310.0),
            &[
                (10.0, ir_film(1.5, 0.01, 2e-6)),
                (110.0, ir_film(1.5, 0.6, 2e-6)),
            ],
            2,
        )
        .unwrap(),
    );
    let mut s = stack(false);
    s.area_m2 = 0.0001;
    s.baths = [CoatingExchangeBath::Blackbody {
        temperature_k: 250.0,
    }; 2];
    for body in &mut s.bodies {
        body.receiver.temperature_bounds_k = (310.0, 310.0);
        if let IsothermalCoatingStack::Transmitting { films, .. } = &mut body.stack {
            *films = vec![CoatingOpticsModel::Law(Arc::clone(&law))];
        }
    }
    let initial = NativeCoatingRadiativeStack::try_new(
        &native,
        &[s],
        NativeCoatingRadiativeStackBudget {
            max_total_optical_knots: 2,
            ..limits()
        },
        || Ok(()),
    )
    .unwrap();
    let (after, observed, audit) = native
        .advance_observed(&initial, 0.04, config(0.005), || Ok(()))
        .unwrap();
    for i in 0..3 {
        let id = CellId(11 + i as u64);
        assert_eq!(after.temperature_k(id).unwrap(), 310.0);
        assert!(after.cells()[&id].specific_h_j_kg < native.cells()[&id].specific_h_j_kg);
        assert_ne!(
            initial.emitter(31, i).unwrap().optical_identity(),
            observed.emitter(31, i).unwrap().optical_identity()
        );
    }
    assert!(
        observed
            .energy()
            .by_stack_pair_j
            .values()
            .all(|&j| j == 0.0)
    );
    assert!(observed.work().optical_evaluations >= 8 * 6 * 8 * 4);
    near(
        after.energy_j() - native.energy_j(),
        audit.external_heat_j,
        2e-12,
    );
}

#[test]
fn radiative_stack_coverage_budgets_and_late_cancellation_preserve_native_state() {
    let native = thermal(cells([1.0; 3]), &HOT, vec![]);
    let s = stack(false);
    for budget in [
        NativeCoatingRadiativeStackBudget {
            max_total_bodies: 2,
            ..limits()
        },
        NativeCoatingRadiativeStackBudget {
            max_channel_spectral_nodes: 79,
            ..limits()
        },
        NativeCoatingRadiativeStackBudget {
            max_optical_evaluations: 191,
            ..limits()
        },
    ] {
        assert!(matches!(
            NativeCoatingRadiativeStack::try_new(&native, &[s.clone()], budget, || Ok(())),
            Err(NativeCoatingExchangeError::Budget)
        ));
    }
    let mut too_large = s.clone();
    too_large.area_m2 = 3.0001;
    assert!(
        NativeCoatingRadiativeStack::try_new(&native, &[too_large], limits(), || Ok(())).is_err()
    );
    let contacted = thermal(
        cells([1.0; 3]),
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
    let mut whole_faces = s.clone();
    whole_faces.area_m2 = 3.0;
    let contact_owner = prepare(&contacted, &whole_faces);
    assert!(
        contacted
            .advance_observed(&contact_owner, 0.001, config(1.0), || Ok(()))
            .is_err()
    );
    assert!(contacted.links().is_empty());

    let initial = prepare(&native, &s);
    let limited = NativeCoatingRadiativeStack::try_new(
        &native,
        &[s],
        NativeCoatingRadiativeStackBudget {
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
    assert_eq!(native.time_s(), 0.0);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.energy().duration_s, 0.0);
    let retry = native
        .advance_observed(&initial, 0.01, config(0.01), || Ok(()))
        .unwrap();
    assert_eq!(
        retry.0.energy_j().to_bits(),
        accepted.0.energy_j().to_bits()
    );
    assert_eq!(retry.1.energy(), accepted.1.energy());
}
