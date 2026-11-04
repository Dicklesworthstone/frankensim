//! G1 network calorimetry, original deposition events and atomic solver refusal.
use fs_solver::NewtonKrylovConfig;
use fs_thermal::crystallization::{
    CrystallizationAdaptivity, CrystallizationCaloric, CrystallizationConfig, CrystallizationError,
    CrystallizationLaw,
};
use fs_thermal::enthalpy::backward_euler::{BoundaryStepCoefficients, StepCoefficients};
use fs_thermal::enthalpy::crystallization::CrystallizationMaterial;
use fs_thermal::enthalpy::crystallization::implicit::{
    CrystallizationImplicitControls, CrystallizationImplicitError,
};
use fs_thermal::enthalpy::crystallization::{
    CrystallizationControls, CrystallizationSimulation, CrystallizationTimelineError,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink,
    ThermalStorageBasis, TimedEvent,
};
use std::collections::BTreeMap;

struct Branch(f64);
impl CrystallizationCaloric for Branch {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        let mut identity = [0x91; 32];
        identity[..8].copy_from_slice(&self.0.to_bits().to_le_bytes());
        Some(identity)
    }
    fn specific_enthalpy_j_kg(&self, temperature: f64) -> Result<f64, CrystallizationError> {
        if !(200.0..=600.0).contains(&temperature) {
            return Err(CrystallizationError::Invalid("original test branch domain"));
        }
        Ok(1000.0 * (temperature - 300.0) + self.0)
    }
    fn minimum_heat_capacity_j_kg_k(&self) -> Option<f64> {
        Some(1000.0)
    }
}
impl ThermalLaw for Branch {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        CrystallizationCaloric::caloric_identity(self)
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let temperature = 300.0 + (h - self.0) / 1000.0;
        if (200.0..=600.0).contains(&temperature) {
            Ok(temperature)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.001
    }
}
static PARENT: Branch = Branch(0.0);
static CRYSTAL: Branch = Branch(-100_000.0);
struct Bath;
impl ThermalLaw for Bath {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn temperature_k(&self, u: f64) -> Result<f64, DepositionError> {
        let temperature = 300.0 + u / 2000.0;
        if (200.0..=600.0).contains(&temperature) {
            Ok(temperature)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.0005
    }
}
static BATH: Bath = Bath;
fn rate(t: f64) -> Result<f64, CrystallizationError> {
    Ok(if t <= 480.0 {
        0.1
    } else if t < 500.0 {
        0.005 * (500.0 - t)
    } else {
        0.0
    })
}
fn law() -> CrystallizationLaw<'static> {
    CrystallizationLaw::new(
        &PARENT,
        &CRYSTAL,
        &rate,
        CrystallizationConfig {
            identity: [0x92; 32],
            avrami_exponent: 1.0,
            maximum_crystalline_mass_fraction: 0.8,
            temperature_range_k: [200.0, 600.0],
            melting_temperature_k: 500.0,
        },
    )
    .unwrap()
}
fn cell(id: u64, material: u64, h: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(material),
        mass_kg: 1.0,
        specific_h_j_kg: h,
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64 + 0.1, 0.0, 0.0],
            volume_m3: 0.001,
            surface_area_m2: 1.0,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            source_w: 0.0,
            reservoir: None,
        },
    }
}
fn link(id: u64, a: u64, b: u64, g: f64) -> Event {
    Event::Contact(ThermalLink {
        id,
        a: CellId(a),
        b: CellId(b),
        area_m2: 0.1,
        conductance_w_k: g,
    })
}
fn simulation<'a>(
    law: &'a CrystallizationLaw<'a>,
    later: Vec<TimedEvent>,
) -> CrystallizationSimulation<'a> {
    let mut events = vec![
        Event::Birth(cell(1, 1, 100_000.0)),
        Event::Birth(cell(2, 2, 0.0)),
        Event::Birth(cell(3, 2, 0.0)),
        link(10, 1, 2, 4000.0),
        link(11, 1, 3, 2000.0),
    ]
    .into_iter()
    .map(|event| TimedEvent { time_s: 0.0, event })
    .collect::<Vec<_>>();
    events.extend(later);
    let thermal = DepositionSimulation::new(
        events,
        &[
            (MaterialId(1), &PARENT as &dyn ThermalLaw),
            (MaterialId(2), &BATH as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 20,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-6,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    CrystallizationSimulation::new(
        thermal,
        &[CrystallizationMaterial {
            material_id: MaterialId(1),
            amorphous: &PARENT,
            law,
            maximum_temperature_slope: 0.001,
        }],
        || Ok(()),
    )
    .unwrap()
    .with_adaptive_phase_integration(MaterialId(1), CrystallizationAdaptivity::default())
    .unwrap()
}
fn controls() -> CrystallizationImplicitControls {
    CrystallizationImplicitControls {
        phase: CrystallizationControls {
            maximum_dt_s: 1.0,
            cell_energy_tolerance_j: 1e-4,
            total_energy_tolerance_j: 1e-3,
            ..CrystallizationControls::default()
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-7,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        difference_step_j: 1.0,
        ..CrystallizationImplicitControls::default()
    }
}
fn empty(
    _: &CrystallizationSimulation<'_>,
) -> Result<BoundaryStepCoefficients<()>, CrystallizationTimelineError> {
    Ok(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    }
    .into())
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

#[test]
fn original_network_solves_latent_growth_with_two_finite_internal_energy_baths() {
    let law = law();
    let old = simulation(&law, vec![]);
    let (next, (), report) = old
        .solve_implicit_step(1.0, controls(), empty, || Ok(()))
        .unwrap();
    // Exact JMAK conversion on the original constant-rate interval, followed by
    // the independently solved three-body backward-Euler heat balance.
    let x = 0.8 * (1.0 - (-0.1_f64).exp());
    let difference = 0.3 * (100.0 + 100.0 * x);
    close(
        next.states()[&CellId(1)].crystalline_mass_fraction(),
        x,
        1e-10,
    );
    close(
        next.temperature_k(CellId(1)).unwrap(),
        300.0 + difference,
        1e-6,
    );
    close(
        next.temperature_k(CellId(2)).unwrap(),
        300.0 + 2.0 * difference / 3.0,
        1e-6,
    );
    close(
        next.temperature_k(CellId(3)).unwrap(),
        300.0 + difference / 2.0,
        1e-6,
    );
    assert_eq!(next.states().len(), 1);
    assert_eq!(
        next.material_storage_basis(MaterialId(2)).unwrap(),
        ThermalStorageBasis::InternalEnergy
    );
    close(next.energy_j(), old.energy_j(), 1e-4);
    assert_eq!(report.phase.thermal.substeps, 1);
    assert_eq!(next.time_s(), 1.0);
    assert!(report.residual_evaluations > 0 && report.newton_iterations > 0);
    assert!(
        report.phase.maximum_cell_energy_residual_j <= controls().phase.cell_energy_tolerance_j
    );
    assert!(report.phase.absolute_energy_residual_j <= controls().phase.total_energy_tolerance_j);
    let (replayed, (), replay) = old
        .solve_implicit_step(1.0, controls(), empty, || Ok(()))
        .unwrap();
    assert_eq!(next.states(), replayed.states());
    assert_eq!(report, replay);
    assert_eq!(old.states()[&CellId(1)].exposure(), 0.0);
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn implicit_march_uses_original_birth_contact_and_separation_events_once() {
    let law = law();
    let old = simulation(
        &law,
        vec![
            TimedEvent {
                time_s: 0.5,
                event: Event::Birth(cell(4, 1, 120_000.0)),
            },
            TimedEvent {
                time_s: 0.5,
                event: link(12, 4, 2, 2500.0),
            },
            TimedEvent {
                time_s: 0.5,
                event: Event::Separation(10),
            },
        ],
    );
    let (whole, report) = old.advance_implicit_to(1.0, controls(), || Ok(())).unwrap();
    let (middle, _) = old.advance_implicit_to(0.5, controls(), || Ok(())).unwrap();
    let (split, _) = middle
        .advance_implicit_to(1.0, controls(), || Ok(()))
        .unwrap();
    assert_eq!(whole.states(), split.states());
    for (id, cell) in whole.cells() {
        assert_eq!(cell.specific_h_j_kg, split.cells()[id].specific_h_j_kg);
    }
    assert_eq!(whole.mass_kg(), 4.0);
    close(whole.states()[&CellId(1)].exposure(), 0.1, 1e-12);
    close(whole.states()[&CellId(4)].exposure(), 0.05, 1e-12);
    assert_eq!(report.phase.thermal.events, 3);
    assert_eq!(report.phase.thermal.substeps, 2);
    assert_eq!(report.phase.thermal.incoming_mass_kg, 1.0);
    assert_eq!(report.phase.thermal.incoming_energy_j, 120_000.0);
    close(report.phase.thermal.stored_energy_change_j, 120_000.0, 1e-3);
    assert!(!whole.links().contains_key(&10));
    assert!(whole.links().contains_key(&12));
    assert_eq!(old.pending_events().len(), 3);
    assert_eq!(old.cells().len(), 3);
}

#[test]
fn original_source_evaluation_budget_and_final_cancellation_publish_no_history() {
    let law = law();
    let old = simulation(&law, vec![]);
    let source_failure = old.solve_implicit_step(
        1.0,
        controls(),
        |candidate| {
            assert!(candidate.states()[&CellId(1)].exposure() > 0.0);
            Err::<BoundaryStepCoefficients<()>, _>(
                DepositionError::InvalidInput("original interface unavailable").into(),
            )
        },
        || Ok(()),
    );
    assert!(matches!(
        source_failure,
        Err(CrystallizationImplicitError::Physical(
            CrystallizationTimelineError::Thermal(DepositionError::InvalidInput(
                "original interface unavailable"
            ))
        ))
    ));
    assert!(matches!(
        old.solve_implicit_step(
            1.0,
            CrystallizationImplicitControls {
                maximum_evaluations: 1,
                ..controls()
            },
            empty,
            || Ok(())
        ),
        Err(CrystallizationImplicitError::Physical(
            CrystallizationTimelineError::Thermal(DepositionError::Budget(_))
        ))
    ));
    let mut calls = 0;
    old.advance_implicit_to(1.0, controls(), || {
        calls += 1;
        Ok(())
    })
    .unwrap();
    let mut visited = 0;
    let refused = old.advance_implicit_to(1.0, controls(), || {
        visited += 1;
        if visited == calls {
            Err(DepositionError::InvalidInput("stop at final publication"))
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        refused,
        Err(CrystallizationImplicitError::Physical(
            CrystallizationTimelineError::Thermal(DepositionError::InvalidInput(
                "stop at final publication"
            ))
        ))
    ));
    assert_eq!(visited, calls);
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.states()[&CellId(1)].exposure(), 0.0);
    assert_eq!(old.energy_j(), 100_000.0);
}

#[test]
fn state_dependent_network_sources_continue_through_the_same_atomic_march() {
    let law = law();
    let old = simulation(&law, vec![]);
    let mut settings = controls();
    settings.phase.maximum_dt_s = 0.5;
    let source = |candidate: &CrystallizationSimulation<'_>| {
        let t = candidate.temperature_k(CellId(1))?;
        let bath = candidate.temperature_k(CellId(2))?;
        let power = 4.0 * (360.0 - candidate.temperature_k(CellId(3))?);
        Ok::<_, CrystallizationTimelineError>(
            StepCoefficients {
                sources_w: BTreeMap::from([(CellId(3), power)]),
                conductances_w_k: BTreeMap::from([(10, 100.0 + 0.1 * (t - bath).powi(2))]),
                payload: (
                    candidate.time_s(),
                    power,
                    candidate.states()[&CellId(1)].exposure(),
                ),
            }
            .into(),
        )
    };
    let (whole, last, report) = old
        .advance_implicit_with_coefficients_to(1.0, settings, source, || Ok(()))
        .unwrap();
    let (middle, first_payload, _) = old
        .solve_implicit_step(0.5, settings, source, || Ok(()))
        .unwrap();
    let (split, second_payload, _) = middle
        .solve_implicit_step(1.0, settings, source, || Ok(()))
        .unwrap();
    assert_eq!(whole.states(), split.states());
    for (id, cell) in whole.cells() {
        assert_eq!(cell.specific_h_j_kg, split.cells()[id].specific_h_j_kg);
    }
    assert_eq!(last, Some(second_payload));
    assert_eq!(second_payload.0, 1.0);
    close(second_payload.2, 0.1, 1e-12);
    close(
        report.phase.thermal.external_heat_j,
        0.5 * (first_payload.1 + second_payload.1),
        1e-6,
    );
    close(
        report.phase.thermal.stored_energy_change_j,
        report.phase.thermal.external_heat_j,
        1e-3,
    );
    let refused = old.advance_implicit_with_coefficients_to(
        1.0,
        settings,
        |candidate| {
            if candidate.time_s() > 0.5 {
                Err(DepositionError::InvalidInput("later original source refused").into())
            } else {
                source(candidate)
            }
        },
        || Ok(()),
    );
    assert!(matches!(
        refused,
        Err(CrystallizationImplicitError::Physical(
            CrystallizationTimelineError::Thermal(DepositionError::InvalidInput(
                "later original source refused"
            ))
        ))
    ));
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.states()[&CellId(1)].exposure(), 0.0);
}
