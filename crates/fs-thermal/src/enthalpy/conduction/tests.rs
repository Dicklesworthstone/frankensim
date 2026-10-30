use super::*;
use crate::observers::{CoupledObservers, ObserverCouplingError};
use std::cell::Cell;

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if positive(h) {
            Ok(h)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.temperature_k(a)?;
        self.temperature_k(b)?;
        Ok(Some(((b - a) / a).ln_1p()))
    }
}
static SENSIBLE: Sensible = Sensible;

struct Constant(f64);
impl ConductivityLaw for Constant {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(self.0)
    }
}
static UNIT_K: Constant = Constant(1.0);

struct LinearK<'a>(&'a Cell<usize>);
impl ConductivityLaw for LinearK<'_> {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.0.set(self.0.get() + 1);
        if positive(h) {
            Ok(h / 100.0)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
}

fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-9,
    }
}
fn boundary(t: Option<f64>) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        source_w: 0.0,
        reservoir: t.map(|temperature_k| Reservoir {
            area_m2: 1.0,
            conductance_w_k: 0.0,
            temperature_k,
        }),
    }
}
fn cells(h: &[f64], width: f64, cold: Option<f64>) -> Vec<DepositedCell> {
    h.iter()
        .enumerate()
        .map(|(i, &specific_h_j_kg)| DepositedCell {
            id: CellId(i as u64),
            body_id: 1,
            segment_id: 1,
            material_id: MaterialId(1),
            mass_kg: width,
            specific_h_j_kg,
            geometry: CellGeometry {
                start_m: [i as f64 * width, 0.0, 0.0],
                end_m: [(i + 1) as f64 * width, 0.0, 0.0],
                volume_m3: width,
                surface_area_m2: 2.0,
            },
            boundary: boundary(if i == 0 { cold } else { None }),
        })
        .collect()
}
fn simulation<'a>(
    cells: Vec<DepositedCell>,
    law: &'a dyn ThermalLaw,
    mut later: Vec<TimedEvent>,
    max_substeps: usize,
) -> DepositionSimulation<'a> {
    let n = cells.len();
    let mut schedule: Vec<_> = cells
        .into_iter()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    for i in 1..n {
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: i as u64,
                a: CellId((i - 1) as u64),
                b: CellId(i as u64),
                area_m2: 1.0,
                conductance_w_k: 0.0,
            }),
        });
    }
    schedule.append(&mut later);
    let event_count = schedule.len();
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: n,
            max_links: n,
            max_events: event_count,
            max_substeps,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0
}
fn link_geometry(n: usize, width: f64) -> Vec<ConductiveLink> {
    (1..n)
        .map(|i| ConductiveLink {
            link_id: i as u64,
            half_distances_m: [0.5 * width; 2],
            contact_resistance_m2_k_w: 0.0,
        })
        .collect()
}
fn face(width: f64, resistance: f64) -> ConductiveReservoir {
    ConductiveReservoir {
        cell_id: CellId(0),
        half_distance_m: 0.5 * width,
        surface_resistance_m2_k_w: resistance,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b} ({tolerance})");
}

#[test]
fn geometric_resistances_query_each_cell_once_and_keep_original_heat_owner() {
    let count = Cell::new(0);
    let conductivity = LinearK(&count);
    let initial = simulation(
        cells(&[300.0, 350.0, 400.0], 1.0, Some(600.0)),
        &SENSIBLE,
        vec![],
        100,
    );
    let mut links = link_geometry(3, 1.0);
    links[1].contact_resistance_m2_k_w = 0.1;
    let observer = ConductionObserver::new(
        &initial,
        &[(MaterialId(1), &conductivity)],
        &links,
        &[face(1.0, 0.2)],
        || Ok(()),
    )
    .unwrap();
    assert_eq!(count.get(), 3);
    count.set(0);
    let rates = initial
        .power_balance_observed(&observer, || Ok(()))
        .unwrap();
    assert_eq!(
        count.get(),
        3,
        "face count must not multiply material queries"
    );
    let g01 = 1.0 / (0.5 / 3.0 + 0.5 / 3.5);
    let g12 = 1.0 / (0.5 / 3.5 + 0.5 / 4.0 + 0.1);
    let bath = 300.0 / (0.5 / 3.0 + 0.2);
    near(rates.net_power_w[&CellId(0)], 50.0 * g01 + bath, 1e-12);
    near(rates.net_power_w[&CellId(1)], 50.0 * (g12 - g01), 1e-12);
    near(rates.net_power_w[&CellId(2)], -50.0 * g12, 1e-12);
    near(rates.external_power_w, bath, 1e-12);
    assert!(
        initial
            .links()
            .values()
            .all(|link| link.conductance_w_k == 0.0)
    );
    assert_eq!(
        initial.cells()[&CellId(0)]
            .boundary
            .reservoir
            .unwrap()
            .conductance_w_k,
        0.0
    );
}

#[test]
fn surface_reservoir_uses_half_cell_cfl_actual_bath_entropy_and_exact_event_time() {
    let conductivity = Constant(4.0);
    let initial = simulation(cells(&[300.0], 1.0, Some(600.0)), &SENSIBLE, vec![], 100);
    let observer = ConductionObserver::new(
        &initial,
        &[(MaterialId(1), &conductivity)],
        &[],
        &[face(1.0, 0.0)],
        || Ok(()),
    )
    .unwrap();
    let coupled = CoupledObservers::new(observer.clone(), ());
    let (next, _, energy, entropy) = initial
        .advance_observed_with_entropy(&coupled, 0.1, config(100.0), || Ok(()))
        .unwrap();
    assert_eq!(energy.substeps, 1);
    near(next.temperature_k(CellId(0)).unwrap(), 540.0, 1e-12);
    near(energy.external_heat_j, 240.0, 1e-12);
    near(entropy.reservoir_entropy_input_j_k, 240.0 / 600.0, 1e-14);
    near(entropy.conduction_generation_j_k, 0.4, 1e-14);
    near(entropy.thermal_storage_change_j_k, 1.8_f64.ln(), 1e-14);
    assert!(!entropy.has_unclassified_sources);
    assert_eq!(
        next.cells()[&CellId(0)]
            .boundary
            .reservoir
            .unwrap()
            .conductance_w_k,
        0.0
    );
    let advected = crate::enthalpy::advection::Advected::new(observer.clone(), vec![], 0).unwrap();
    let (same, _, same_energy, same_entropy) = initial
        .advance_observed_with_entropy(&advected, 0.1, config(100.0), || Ok(()))
        .unwrap();
    assert_eq!(same.cells(), next.cells());
    assert_eq!(same_energy, energy);
    assert_eq!(same_entropy, entropy);
    let conflict = CoupledObservers::new(observer.clone(), observer);
    assert!(matches!(
        initial.power_balance_observed(&conflict, || Ok(())),
        Err(ObserverCouplingError::ConflictingReservoirConductance(
            CellId(0)
        ))
    ));

    let event = TimedEvent {
        time_s: 0.05,
        event: Event::Boundary(CellId(0), boundary(Some(250.0))),
    };
    let initial = simulation(
        cells(&[300.0], 1.0, Some(600.0)),
        &SENSIBLE,
        vec![event],
        100,
    );
    let observer = ConductionObserver::new(
        &initial,
        &[(MaterialId(1), &conductivity)],
        &[],
        &[face(1.0, 0.0)],
        || Ok(()),
    )
    .unwrap();
    let (next, _, energy, entropy) = initial
        .advance_observed_with_entropy(&observer, 0.1, config(100.0), || Ok(()))
        .unwrap();
    assert_eq!((energy.substeps, energy.events), (2, 1));
    near(next.temperature_k(CellId(0)).unwrap(), 352.0, 1e-12);
    near(energy.external_heat_j, 120.0 - 68.0, 1e-12);
    near(
        entropy.reservoir_entropy_input_j_k,
        120.0 / 600.0 - 68.0 / 250.0,
        1e-14,
    );
    assert!(!entropy.has_unclassified_sources);
}

#[test]
fn nonlinear_conductivity_evolves_on_actual_substeps_and_refines_to_logistic_solution() {
    // cp=m=A=1, d=1/2, k(T)=T/100 and Tb=600 give
    // dT/dt=T(600-T)/50; T(0)=300 => T=600/(1+exp(-12t)).
    let exact = 600.0 / (1.0 + (-1.2_f64).exp());
    let count = Cell::new(0);
    let law = LinearK(&count);
    let initial = simulation(cells(&[300.0], 1.0, Some(600.0)), &SENSIBLE, vec![], 1000);
    let observer = ConductionObserver::new(
        &initial,
        &[(MaterialId(1), &law)],
        &[],
        &[face(1.0, 0.0)],
        || Ok(()),
    )
    .unwrap();
    let mut errors = Vec::new();
    for dt in [0.005, 0.0025] {
        count.set(0);
        let (next, _, energy, entropy) = initial
            .advance_observed_with_entropy(&observer, 0.1, config(dt), || Ok(()))
            .unwrap();
        assert_eq!(count.get(), energy.substeps);
        let t = next.temperature_k(CellId(0)).unwrap();
        errors.push((t - exact).abs());
        near(energy.external_heat_j, t - 300.0, 1e-10);
        near(
            entropy.reservoir_entropy_input_j_k,
            (t - 300.0) / 600.0,
            1e-13,
        );
        assert!(!entropy.has_unclassified_sources);
    }
    assert!((1.8..2.2).contains(&(errors[0] / errors[1])), "{errors:?}");
}

#[test]
fn entropy_wrapper_forwards_each_rate_hook_once_and_keeps_internal_transfers_unclassified() {
    #[derive(Clone, Default)]
    struct Rates([usize; 6]);
    impl DepositionObserver for Rates {
        type Error = DepositionError;
        fn sources(
            &mut self,
            _: &DepositionSimulation<'_>,
            _: &mut dyn FnMut() -> Result<(), DepositionError>,
        ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
            self.0[0] += 1;
            Ok([(CellId(0), 5.0)].into())
        }
        fn conductances(
            &mut self,
            _: &DepositionSimulation<'_>,
            _: &mut dyn FnMut() -> Result<(), DepositionError>,
        ) -> Result<BTreeMap<u64, f64>, Self::Error> {
            self.0[1] += 1;
            Ok([(1, 2.0)].into())
        }
        fn reservoir_conductances(
            &mut self,
            _: &DepositionSimulation<'_>,
            _: &mut dyn FnMut() -> Result<(), DepositionError>,
        ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
            self.0[2] += 1;
            Ok([(CellId(0), 3.0)].into())
        }
        fn heat_transfers(
            &mut self,
            _: &DepositionSimulation<'_>,
            _: &mut dyn FnMut() -> Result<(), DepositionError>,
        ) -> Result<Vec<HeatTransfer>, Self::Error> {
            self.0[3] += 1;
            Ok(vec![HeatTransfer {
                a: CellId(0),
                b: CellId(1),
                power_w: 10.0,
            }])
        }
        fn requires_step_interval(&self) -> bool {
            true
        }
        fn step_sources(
            &mut self,
            _: &DepositionSimulation<'_>,
            _: f64,
            _: f64,
            _: f64,
            _: &mut dyn FnMut() -> Result<(), DepositionError>,
        ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
            self.0[4] += 1;
            Ok([(CellId(1), 7.0)].into())
        }
        fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
            Ok(())
        }
        fn step(
            &mut self,
            _: f64,
            _: f64,
            _: &BTreeMap<CellId, f64>,
            _: &DepositionSimulation<'_>,
            _: &mut dyn FnMut() -> Result<(), DepositionError>,
        ) -> Result<(), Self::Error> {
            self.0[5] += 1;
            Ok(())
        }
    }
    let mut body = cells(&[300.0, 400.0], 1.0, Some(600.0));
    body[0].boundary.source_w = 3.0;
    let initial = simulation(body, &SENSIBLE, vec![], 100);
    let observer = Rates::default();
    let dt = 0.001;
    let (next, rates, energy, entropy) = initial
        .advance_observed_with_entropy(&observer, dt, config(dt), || Ok(()))
        .unwrap();
    assert_eq!(rates.0, [1; 6]);
    assert_eq!(observer.0, [0; 6]);
    near(next.temperature_k(CellId(0)).unwrap(), 301.098, 1e-12);
    near(next.temperature_k(CellId(1)).unwrap(), 399.817, 1e-12);
    near(
        energy.external_heat_j,
        dt * (900.0 + 3.0 + 5.0 + 7.0),
        1e-14,
    );
    near(
        entropy.reservoir_entropy_input_j_k,
        dt * 900.0 / 600.0,
        1e-14,
    );
    near(
        entropy.conduction_generation_j_k,
        dt * (200.0 / 1200.0 + 900.0 / 600.0),
        1e-14,
    );
    near(
        entropy.unclassified_source_at_cell_j_k,
        dt * (-2.0 / 300.0 + 17.0 / 400.0),
        1e-14,
    );
    assert!(entropy.has_unclassified_sources);
}

struct Stefan;
impl ThermalLaw for Stefan {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = if h < 0.0 {
            300.0 + h
        } else if h <= 1.0 {
            300.0
        } else {
            299.0 + h
        };
        if positive(t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let entropy = |h: f64| -> Result<f64, DepositionError> {
            let t = self.temperature_k(h)?;
            Ok(if h < 0.0 {
                (t / 300.0).ln()
            } else if h <= 1.0 {
                h / 300.0
            } else {
                1.0 / 300.0 + (t / 300.0).ln()
            })
        };
        Ok(Some(entropy(b)? - entropy(a)?))
    }
}
static STEFAN: Stefan = Stefan;

#[test]
fn computed_freezing_front_and_heat_refine_to_one_phase_stefan_similarity() {
    // Same one-phase similarity family as the existing fs-conduction Stefan
    // example: rho=cp=k=L=1, lambda=1/2. A cold prescribed FACE and uniform
    // initial liquid at Tm are the only inputs; no front is imposed in time.
    // DeltaT=lambda*sqrt(pi)*exp(lambda²)*erf(lambda) gives s(t)=sqrt(t),
    // and independently integrated boundary heat Q(t)=-exp(1/4)*sqrt(t).
    let lambda = 0.5_f64;
    let delta =
        lambda * std::f64::consts::PI.sqrt() * (lambda * lambda).exp() * fs_math::det::erf(lambda);
    let time = 0.16_f64;
    let front = time.sqrt();
    let heat = -(0.25_f64).exp() * time.sqrt();
    let mut errors = Vec::new();
    for n in [24, 48, 96] {
        let dx = 1.0 / n as f64;
        let initial = simulation(
            cells(&vec![1.0; n], dx, Some(300.0 - delta)),
            &STEFAN,
            vec![],
            100_000,
        );
        let observer = ConductionObserver::new(
            &initial,
            &[(MaterialId(1), &UNIT_K)],
            &link_geometry(n, dx),
            &[face(dx, 0.0)],
            || Ok(()),
        )
        .unwrap();
        let (next, _, audit) = initial
            .advance_observed(&observer, time, config(0.2 * dx * dx), || Ok(()))
            .unwrap();
        let frozen_length = next
            .cells()
            .values()
            .map(|c| dx * (1.0 - c.specific_h_j_kg.clamp(0.0, 1.0)))
            .sum::<f64>();
        let mut temperature_squared = 0.0;
        for (i, cell) in next.cells().values().enumerate() {
            let x = (i as f64 + 0.5) * dx;
            let exact = if x < front {
                300.0 - delta
                    + delta * fs_math::det::erf(x / (2.0 * time.sqrt())) / fs_math::det::erf(lambda)
            } else {
                300.0
            };
            let actual = next.temperature_k(cell.id).unwrap();
            temperature_squared += dx * (actual - exact).powi(2);
            assert!((300.0 - delta - 1e-12..=300.0).contains(&actual));
        }
        let row = [
            (frozen_length - front).abs(),
            (audit.external_heat_j - heat).abs(),
            temperature_squared.sqrt(),
        ];
        println!(
            "Stefan cells={n}: frozen_length={frozen_length:.12}, heat={:.12}, front/heat/RMS errors={row:?}, closure={:.3e}",
            audit.external_heat_j, audit.residual_j
        );
        near(audit.residual_j, 0.0, 1e-9);
        errors.push(row);
    }
    for metric in 0..3 {
        assert!(
            errors[2][metric] < errors[1][metric] && errors[1][metric] < errors[0][metric],
            "{errors:?}"
        );
    }
    assert!(
        errors[2][0] < 0.002 && errors[2][1] < 0.003 && errors[2][2] < 0.002,
        "{errors:?}"
    );
}

#[test]
fn unsupported_geometry_source_cancellation_and_budget_leave_original_state_unchanged() {
    let initial = simulation(
        cells(&[300.0, 350.0], 1.0, Some(600.0)),
        &SENSIBLE,
        vec![],
        2,
    );
    let bindings = link_geometry(2, 1.0);
    let observer = ConductionObserver::new(
        &initial,
        &[(MaterialId(1), &UNIT_K)],
        &bindings,
        &[face(1.0, 0.0)],
        || Ok(()),
    )
    .unwrap();
    let before = initial.cells().clone();
    assert!(matches!(
        initial.advance_observed_with_entropy(&observer, 1.0, config(0.01), || Ok(())),
        Err(DepositionError::Budget("thermal substeps"))
    ));
    assert_eq!(initial.cells(), &before);
    assert!(matches!(
        initial.advance_observed_with_entropy(&observer, 1.0, config(0.01), || Err(
            DepositionError::Cancelled
        )),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(initial.cells(), &before);
    let invalid = Constant(f64::NAN);
    assert!(matches!(
        ConductionObserver::new(
            &initial,
            &[(MaterialId(1), &invalid)],
            &bindings,
            &[],
            || Ok(())
        ),
        Err(DepositionError::MaterialDomain)
    ));
    let mut changed = initial.clone();
    changed.links.get_mut(&1).unwrap().area_m2 = 0.5;
    assert!(
        changed
            .power_balance_observed(&observer, || Ok(()))
            .is_err()
    );
    changed = initial.clone();
    changed.cells.get_mut(&CellId(0)).unwrap().mass_kg = 2.0;
    assert!(
        changed
            .power_balance_observed(&observer, || Ok(()))
            .is_err()
    );
    let insulation = Constant(0.0);
    let observer = ConductionObserver::new(
        &initial,
        &[(MaterialId(1), &insulation)],
        &bindings,
        &[face(1.0, 0.0)],
        || Ok(()),
    )
    .unwrap();
    let (next, _, audit) = initial
        .advance_observed(&observer, 1.0, config(1.0), || Ok(()))
        .unwrap();
    assert_eq!(next.cells(), &before);
    assert_eq!(audit.external_heat_j, 0.0);
}
