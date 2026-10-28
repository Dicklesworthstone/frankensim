use super::*;
use crate::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};

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
static LAW: Sensible = Sensible;
fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-9,
    }
}
fn cell(id: u64, h: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: 1,
        material_id: MaterialId(1),
        mass_kg: 1.0,
        specific_h_j_kg: h,
        geometry: CellGeometry {
            start_m: [0.0, 0.0, id as f64],
            end_m: [1.0, 0.0, id as f64],
            volume_m3: 1.0,
            surface_area_m2: 1.0,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}
fn simulation<'a>(
    cells: Vec<DepositedCell>,
    linked: bool,
    law: &'a dyn ThermalLaw,
) -> DepositionSimulation<'a> {
    let mut events: Vec<_> = cells
        .into_iter()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    if linked {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 1,
                a: CellId(1),
                b: CellId(2),
                area_m2: 0.5,
                conductance_w_k: 1.0,
            }),
        });
    }
    let sim = DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: 10_000,
        },
    )
    .unwrap();
    sim.advance_to(0.0, config(0.1), || Ok(())).unwrap().0
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-12, "{a} != {b}");
}

// G1: exact two-body state entropy and original explicit trajectory.
#[test]
fn isolated_exchange_retains_energy_and_separates_explicit_entropy_defect() {
    let sim = simulation(vec![cell(1, 300.0), cell(2, 600.0)], true, &LAW);
    let (plain, plain_energy) = sim.advance_to(1.0, config(0.125), || Ok(())).unwrap();
    let (next, energy, entropy) = sim
        .advance_to_with_entropy(1.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(next.cells(), plain.cells());
    assert_eq!(energy, plain_energy);
    close(next.energy_j(), 900.0);
    let expected = (next.temperature_k(CellId(1)).unwrap() / 300.0).ln()
        + (next.temperature_k(CellId(2)).unwrap() / 600.0).ln();
    close(entropy.thermal_storage_change_j_k, expected);
    assert!(entropy.thermal_storage_change_j_k > 0.0);
    assert!(entropy.conduction_generation_j_k > 0.0);
    assert!(entropy.temporal_defect_j_k < 0.0);
    assert_eq!(entropy.reservoir_entropy_input_j_k, 0.0);
    assert!(!entropy.has_unclassified_sources);
    let (_, _, fine) = sim
        .advance_to_with_entropy(1.0, config(0.0625), || Ok(()))
        .unwrap();
    assert!((1.6..2.6).contains(&(entropy.temporal_defect_j_k / fine.temporal_defect_j_k)));
}

// G1: same endpoint watts, exact state entropy, positive BE temporal defect.
#[test]
fn implicit_exchange_uses_endpoint_temperature_and_original_residual() {
    let sim = simulation(vec![cell(1, 300.0), cell(2, 600.0)], true, &LAW);
    let step = BackwardEulerStep::new(&sim, 1.0, || Ok(())).unwrap();
    let resolve = |_: &DepositionSimulation<'_>| {
        Ok::<_, DepositionError>(StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        })
    };
    let trial = step
        .evaluate_with_entropy(&[100.0, -100.0], resolve, || Ok(()))
        .unwrap();
    let entropy = *trial.entropy_audit().unwrap();
    close(entropy.thermal_storage_change_j_k, (10.0_f64 / 9.0).ln());
    close(entropy.conduction_generation_j_k, 0.05);
    assert!(entropy.temporal_defect_j_k > 0.0);
    assert!(!entropy.has_unclassified_sources);
    let (next, (), _) = trial.accept(1e-10, 1e-10, || Ok(())).unwrap();
    close(next.temperature_k(CellId(1)).unwrap(), 400.0);
    close(next.temperature_k(CellId(2)).unwrap(), 500.0);
    assert!(
        step.evaluate(&[100.0, -100.0], resolve, || Ok(()))
            .unwrap()
            .entropy_audit()
            .is_none()
    );
}

#[test]
fn ambient_and_support_use_actual_reservoir_entropy() {
    for support in [false, true] {
        let mut c = cell(1, 300.0);
        if support {
            c.boundary.reservoir = Some(Reservoir {
                area_m2: 0.5,
                conductance_w_k: 1.0,
                temperature_k: 600.0,
            });
        } else {
            c.boundary.ambient_k = 600.0;
            c.boundary.htc_w_m2_k = 1.0;
        }
        let sim = simulation(vec![c], false, &LAW);
        let (next, energy, entropy) = sim
            .advance_to_with_entropy(0.25, config(0.25), || Ok(()))
            .unwrap();
        close(energy.external_heat_j, 75.0);
        close(next.temperature_k(CellId(1)).unwrap(), 375.0);
        close(entropy.reservoir_entropy_input_j_k, 75.0 / 600.0);
        close(
            entropy.conduction_generation_j_k,
            75.0 * (1.0 / 300.0 - 1.0 / 600.0),
        );
        close(entropy.thermal_storage_change_j_k, 1.25_f64.ln());
        assert!(entropy.temporal_defect_j_k < 0.0);
    }
}

#[test]
fn equilibrium_has_no_entropy_production_or_numerical_defect() {
    let sim = simulation(vec![cell(1, 300.0), cell(2, 300.0)], true, &LAW);
    let (_, _, entropy) = sim
        .advance_to_with_entropy(1.0, config(0.1), || Ok(()))
        .unwrap();
    assert_eq!(entropy, EntropyAudit::default());
}

#[test]
fn latent_storage_and_unclassified_source_are_each_counted_once() {
    struct Plateau;
    impl ThermalLaw for Plateau {
        fn temperature_k(&self, _: f64) -> Result<f64, DepositionError> {
            Ok(300.0)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            0.0
        }
        fn specific_entropy_change_j_kg_k(
            &self,
            a: f64,
            b: f64,
        ) -> Result<Option<f64>, DepositionError> {
            Ok(Some((b - a) / 300.0))
        }
    }
    let mut c = cell(1, 100.0);
    c.boundary.source_w = 6.0;
    let sim = simulation(vec![c], false, &Plateau);
    let (next, energy, entropy) = sim
        .advance_to_with_entropy(5.0, config(5.0), || Ok(()))
        .unwrap();
    close(next.temperature_k(CellId(1)).unwrap(), 300.0);
    close(energy.external_heat_j, 30.0);
    close(entropy.thermal_storage_change_j_k, 0.1);
    close(entropy.unclassified_source_at_cell_j_k, 0.1);
    close(entropy.temporal_defect_j_k, 0.0);
    assert_eq!(entropy.conduction_generation_j_k, 0.0);
    assert!(entropy.has_unclassified_sources);
}

#[test]
fn unavailable_entropy_and_cancellation_never_publish_partial_states() {
    struct Missing;
    impl ThermalLaw for Missing {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            LAW.temperature_k(h)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            1.0
        }
    }
    let sim = simulation(vec![cell(1, 300.0), cell(2, 600.0)], true, &Missing);
    assert!(sim.advance_to(0.25, config(0.125), || Ok(())).is_ok());
    assert!(matches!(
        sim.advance_to_with_entropy(0.25, config(0.125), || Ok(())),
        Err(DepositionError::InvalidInput(
            "thermal law does not provide entropy change"
        ))
    ));
    close(sim.temperature_k(CellId(1)).unwrap(), 300.0);
    let sim = simulation(vec![cell(1, 300.0), cell(2, 600.0)], true, &LAW);
    let mut polls = 0;
    assert!(matches!(
        sim.advance_to_with_entropy(1.0, config(0.125), || {
            polls += 1;
            if polls > 8 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    ));
    close(sim.temperature_k(CellId(1)).unwrap(), 300.0);
    assert_eq!(sim.time_s(), 0.0);
}

#[test]
fn deferred_birth_counts_only_entropy_after_each_cells_arrival() {
    let mut first = cell(1, 300.0);
    let mut second = cell(2, 300.0);
    first.boundary.source_w = 6.0;
    second.boundary.source_w = 6.0;
    let sim = DepositionSimulation::new(
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(first),
            },
            TimedEvent {
                time_s: 0.5,
                event: Event::Birth(second),
            },
        ],
        &[(MaterialId(1), &LAW as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: 100,
        },
    )
    .unwrap();
    let (next, energy, entropy) = sim
        .advance_to_with_entropy(1.0, config(0.25), || Ok(()))
        .unwrap();
    close(next.temperature_k(CellId(1)).unwrap(), 306.0);
    close(next.temperature_k(CellId(2)).unwrap(), 303.0);
    close(energy.external_heat_j, 9.0);
    close(
        entropy.thermal_storage_change_j_k,
        (306.0_f64 / 300.0).ln() + (303.0_f64 / 300.0).ln(),
    );
    assert_eq!(entropy.conduction_generation_j_k, 0.0);
    assert_eq!(entropy.reservoir_entropy_input_j_k, 0.0);
    assert!(entropy.has_unclassified_sources);
}

#[test]
fn implicit_trial_audits_resolved_sources_before_restoring_authored_inputs() {
    let sim = simulation(vec![cell(1, 300.0)], false, &LAW);
    let step = BackwardEulerStep::new(&sim, 1.0, || Ok(())).unwrap();
    let trial = step
        .evaluate_with_entropy(
            &[10.0],
            |_| {
                Ok::<_, DepositionError>(StepCoefficients {
                    sources_w: BTreeMap::from([(CellId(1), 10.0)]),
                    conductances_w_k: BTreeMap::new(),
                    payload: (),
                })
            },
            || Ok(()),
        )
        .unwrap();
    let entropy = *trial.entropy_audit().unwrap();
    close(entropy.thermal_storage_change_j_k, (31.0_f64 / 30.0).ln());
    close(entropy.unclassified_source_at_cell_j_k, 10.0 / 310.0);
    assert_eq!(entropy.conduction_generation_j_k, 0.0);
    assert_eq!(entropy.reservoir_entropy_input_j_k, 0.0);
    assert!(entropy.has_unclassified_sources);
    assert!(entropy.temporal_defect_j_k > 0.0);
    let (next, (), energy) = trial.accept(1e-10, 1e-10, || Ok(())).unwrap();
    close(energy.external_heat_j, 10.0);
    assert_eq!(next.cells[&CellId(1)].boundary.source_w, 0.0);
}
