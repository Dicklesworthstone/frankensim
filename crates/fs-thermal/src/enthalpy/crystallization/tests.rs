//! G0/G1 manufactured calorimetry and G3/G4 timeline/cancellation checks.
#[path = "tests/adaptive.rs"]
mod adaptive;
#[path = "tests/remelting.rs"]
mod remelting;

use super::*;
use crate::crystallization::{CrystallizationCaloric, CrystallizationConfig, CrystallizationRate};
use crate::enthalpy::{CellGeometry, DepositionBudget, StepConfig, ThermalLink, TimedEvent};

struct Branch {
    heat_capacity_j_kg_k: f64,
    offset: f64,
}
impl Branch {
    fn identity(&self) -> [u8; 32] {
        let mut identity = [0xCB; 32];
        identity[..8].copy_from_slice(&self.heat_capacity_j_kg_k.to_bits().to_le_bytes());
        identity[8..16].copy_from_slice(&self.offset.to_bits().to_le_bytes());
        identity
    }
}
impl ThermalLaw for Branch {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some(self.identity())
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + (h - self.offset) / self.heat_capacity_j_kg_k;
        if (250.0..=500.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        self.heat_capacity_j_kg_k.recip()
    }
}
impl CrystallizationCaloric for Branch {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some(self.identity())
    }
    fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, CrystallizationError> {
        if !(250.0..=500.0).contains(&t) {
            return Err(CrystallizationError::Invalid("manufactured branch domain"));
        }
        Ok(self.heat_capacity_j_kg_k * (t - 300.0) + self.offset)
    }
}
static AMORPHOUS: Branch = Branch {
    heat_capacity_j_kg_k: 1000.0,
    offset: 0.0,
};
static CRYSTALLINE: Branch = Branch {
    heat_capacity_j_kg_k: 1000.0,
    offset: -100_000.0,
};

/// Numerically identical, but deliberately omits caloric identity on BOTH traits.
struct UnidentifiedBranch;
impl ThermalLaw for UnidentifiedBranch {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        ThermalLaw::temperature_k(&AMORPHOUS, h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        ThermalLaw::maximum_temperature_slope(&AMORPHOUS)
    }
}
impl CrystallizationCaloric for UnidentifiedBranch {
    fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, CrystallizationError> {
        CrystallizationCaloric::specific_enthalpy_j_kg(&AMORPHOUS, t)
    }
}
static UNIDENTIFIED: UnidentifiedBranch = UnidentifiedBranch;

/// Distinct, noncrystallizing finite mold: cp=500 J/(kg K).
struct Mold;
impl ThermalLaw for Mold {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let temperature = 300.0 + h / 500.0;
        if (250.0..=800.0).contains(&temperature) {
            Ok(temperature)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / 500.0
    }
}
static MOLD: Mold = Mold;

struct Undercooling;
impl CrystallizationRate for Undercooling {
    fn rate_s_inv(&self, t: f64) -> Result<f64, CrystallizationError> {
        Ok(0.002 * (500.0 - t).max(0.0))
    }
}
static RATE: Undercooling = Undercooling;

fn law() -> CrystallizationLaw<'static> {
    CrystallizationLaw::new(
        &AMORPHOUS,
        &CRYSTALLINE,
        &RATE,
        CrystallizationConfig {
            identity: [37; 32],
            avrami_exponent: 2.0,
            maximum_crystalline_mass_fraction: 0.8,
            temperature_range_k: [250.0, 500.0],
            melting_temperature_k: 500.0,
        },
    )
    .unwrap()
}
fn boundary(htc: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: htc,
        source_w: 0.0,
        reservoir: None,
    }
}
fn thermal(temperatures: &[f64], htc: f64, conductance: f64) -> DepositionSimulation<'static> {
    let mut events = Vec::new();
    for (i, &t) in temperatures.iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64 + 10),
                segment_id: 81,
                body_id: 7,
                material_id: MaterialId(3),
                mass_kg: 0.01,
                specific_h_j_kg: CrystallizationCaloric::specific_enthalpy_j_kg(&AMORPHOUS, t)
                    .unwrap(),
                geometry: CellGeometry {
                    start_m: [i as f64 * 0.01, 0.0, 0.0],
                    end_m: [(i + 1) as f64 * 0.01, 0.0, 0.0],
                    volume_m3: 1e-5,
                    surface_area_m2: 0.01,
                },
                boundary: boundary(htc),
            }),
        });
    }
    for i in 1..temperatures.len() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: i as u64,
                a: CellId(i as u64 + 9),
                b: CellId(i as u64 + 10),
                area_m2: 0.001,
                conductance_w_k: conductance,
            }),
        });
    }
    activate(events, &[(MaterialId(3), &AMORPHOUS)])
}
fn activate(
    events: Vec<TimedEvent>,
    materials: &[(MaterialId, &'static dyn ThermalLaw)],
) -> DepositionSimulation<'static> {
    DepositionSimulation::new(
        events,
        materials,
        DepositionBudget {
            max_cells: 16,
            max_links: 16,
            max_events: 32,
            max_substeps: 100_000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn finite_bath() -> DepositionSimulation<'static> {
    let specimen = DepositedCell {
        id: CellId(10),
        segment_id: 81,
        body_id: 7,
        material_id: MaterialId(3),
        mass_kg: 0.01,
        specific_h_j_kg: 80_000.0,
        geometry: CellGeometry {
            start_m: [0.0, 0.0, 0.0],
            end_m: [0.01, 0.0, 0.0],
            volume_m3: 1e-5,
            surface_area_m2: 0.01,
        },
        boundary: boundary(0.0),
    };
    let mold = DepositedCell {
        id: CellId(11),
        segment_id: 82,
        body_id: 8,
        material_id: MaterialId(4),
        mass_kg: 0.05,
        specific_h_j_kg: 0.0,
        geometry: CellGeometry {
            start_m: [0.01, 0.0, 0.0],
            end_m: [0.02, 0.0, 0.0],
            volume_m3: 1e-5,
            surface_area_m2: 0.01,
        },
        boundary: boundary(0.0),
    };
    let events = vec![
        TimedEvent {
            time_s: 0.0,
            event: Event::Birth(specimen),
        },
        TimedEvent {
            time_s: 0.0,
            event: Event::Birth(mold),
        },
        TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 1,
                a: CellId(10),
                b: CellId(11),
                area_m2: 0.001,
                conductance_w_k: 1.25,
            }),
        },
    ];
    activate(
        events,
        &[(MaterialId(3), &AMORPHOUS), (MaterialId(4), &MOLD)],
    )
}
fn material<'a>(law: &'a CrystallizationLaw<'a>) -> CrystallizationMaterial<'a> {
    CrystallizationMaterial {
        material_id: MaterialId(3),
        amorphous: &AMORPHOUS,
        law,
        maximum_temperature_slope: 0.001,
    }
}
fn controls(dt: f64) -> CrystallizationControls {
    CrystallizationControls {
        maximum_dt_s: dt,
        caloric: CrystallizationSolve {
            maximum_iterations: 128,
            specific_enthalpy_tolerance_j_kg: 1e-7,
            temperature_tolerance_k: 1e-10,
        },
        cell_energy_tolerance_j: 1e-7,
        total_energy_tolerance_j: 1e-4,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

#[test]
fn spatial_cooling_has_recalescence_on_one_original_energy_account() {
    let law = law();
    let original = thermal(&[380.0, 390.0], 20.0, 0.5);
    let initial_h = original.energy_j();
    let initial_mass = original.mass_kg();
    let original_cells = original.cells.clone();
    let start = CrystallizationSimulation::new(original, &[material(&law)], || Ok(())).unwrap();
    let (next, audit) = start.advance_to(5.0, controls(0.025), || Ok(())).unwrap();
    assert!(audit.thermal.external_heat_j < 0.0);
    assert!(next.temperature_k(CellId(10)).unwrap() > 380.0);
    assert!(next.temperature_k(CellId(11)).unwrap() > 390.0);
    assert!(audit.crystallized_mass_kg > 0.0);
    assert_eq!(next.mass_kg(), initial_mass);
    close(
        next.energy_j() - initial_h,
        audit.thermal.external_heat_j,
        1e-8,
    );
    close(audit.thermal.residual_j, 0.0, 1e-8);
    for (&id, cell) in next.cells() {
        let state = next.states()[&id];
        let original = &original_cells[&id];
        assert_eq!(
            (cell.segment_id, cell.body_id, cell.material_id),
            (original.segment_id, original.body_id, original.material_id)
        );
        assert_eq!(cell.mass_kg, original.mass_kg);
        assert_eq!(cell.geometry, original.geometry);
        assert_eq!(cell.specific_h_j_kg, state.specific_enthalpy_j_kg());
        assert!(state.crystalline_mass_fraction() > 0.0 && state.crystalline_mass_fraction() < 0.8);
        // Independent constant-Cp mixture: h=cp(T-300)-L*X.
        close(
            state.specific_enthalpy_j_kg(),
            1000.0 * (state.temperature_k() - 300.0)
                - 100_000.0 * state.crystalline_mass_fraction(),
            1e-6,
        );
    }
    assert_eq!(start.time_s(), 0.0);
    assert_eq!(start.states()[&CellId(10)].crystalline_mass_fraction(), 0.0);
}

#[test]
fn adiabatic_crystallization_preserves_enthalpy_and_mass_exactly() {
    let law = law();
    let start = CrystallizationSimulation::new(
        thermal(&[380.0, 410.0], 0.0, 0.0),
        &[material(&law)],
        || Ok(()),
    )
    .unwrap();
    let (next, audit) = start.advance_to(10.0, controls(0.05), || Ok(())).unwrap();
    assert_eq!(next.energy_j(), start.energy_j());
    assert_eq!(next.mass_kg(), start.mass_kg());
    assert_eq!(audit.thermal.external_heat_j, 0.0);
    assert_eq!(audit.thermal.stored_energy_change_j, 0.0);
    for (&id, state) in next.states() {
        close(
            next.temperature_k(id).unwrap() - start.temperature_k(id).unwrap(),
            100.0 * state.crystalline_mass_fraction(),
            1e-8,
        );
    }
}

#[test]
fn spatial_timestep_refinement_converges_without_changing_heat_capacity() {
    let law = law();
    let start = CrystallizationSimulation::new(
        thermal(&[370.0, 420.0], 35.0, 1.0),
        &[material(&law)],
        || Ok(()),
    )
    .unwrap();
    let run = |dt| start.advance_to(2.0, controls(dt), || Ok(())).unwrap().0;
    let reference = run(1.0 / 2048.0);
    let error = |state: &CrystallizationSimulation<'_>| {
        state
            .states()
            .iter()
            .map(|(id, s)| (s.temperature_k() - reference.states()[id].temperature_k()).abs())
            .sum::<f64>()
    };
    let errors = [0.125, 0.0625, 0.03125].map(|dt| error(&run(dt)));
    assert!(
        errors[0] > 1.7 * errors[1] && errors[1] > 1.7 * errors[2],
        "{errors:?}"
    );
}

#[test]
fn finite_mold_exchanges_conserved_heat_and_refines_without_a_fake_phase_state() {
    let law = law();
    let original = finite_bath();
    let initial_energy = original.energy_j();
    let initial_mass = original.mass_kg();
    let original_mold = original.cells()[&CellId(11)].clone();
    let start = CrystallizationSimulation::new(original, &[material(&law)], || Ok(())).unwrap();
    assert_eq!(start.states().len(), 1);
    assert!(!start.states().contains_key(&CellId(11)));
    assert_eq!(start.temperature_k(CellId(11)).unwrap(), 300.0);
    assert!(start.temperature_k(CellId(999)).is_err());

    let (next, audit) = start.advance_to(4.0, controls(0.0625), || Ok(())).unwrap();
    assert_eq!(next.mass_kg(), initial_mass);
    assert_eq!(audit.thermal.external_heat_j, 0.0);
    close(next.energy_j(), initial_energy, 1e-8);
    close(audit.thermal.stored_energy_change_j, 0.0, 1e-8);
    close(audit.thermal.residual_j, 0.0, 1e-8);
    assert!(audit.crystallized_mass_kg > 0.0);
    close(
        audit.crystallized_mass_kg,
        next.cells()[&CellId(10)].mass_kg * next.states()[&CellId(10)].crystalline_mass_fraction(),
        1e-14,
    );
    assert_eq!(next.states().len(), 1);
    assert!(!next.states().contains_key(&CellId(11)));
    let mold = &next.cells()[&CellId(11)];
    assert_eq!(mold.material_id, original_mold.material_id);
    assert_eq!(mold.body_id, original_mold.body_id);
    assert_eq!(mold.segment_id, original_mold.segment_id);
    assert_eq!(mold.geometry, original_mold.geometry);
    assert_eq!(mold.mass_kg, original_mold.mass_kg);
    assert!(mold.specific_h_j_kg > original_mold.specific_h_j_kg);
    assert!(next.temperature_k(CellId(11)).unwrap() > 300.0);
    // Independent finite-bath calorimetry: its heat gain is m*cp*delta T.
    close(
        mold.mass_kg * mold.specific_h_j_kg,
        0.05 * 500.0 * (next.temperature_k(CellId(11)).unwrap() - 300.0),
        1e-10,
    );
    close(
        next.cells()[&CellId(10)].mass_kg
            * (start.cells()[&CellId(10)].specific_h_j_kg
                - next.cells()[&CellId(10)].specific_h_j_kg),
        mold.mass_kg * mold.specific_h_j_kg,
        1e-8,
    );

    let reference = start
        .advance_to(4.0, controls(1.0 / 1024.0), || Ok(()))
        .unwrap()
        .0;
    let errors = [0.25, 0.125, 0.0625].map(|dt| {
        let candidate = start.advance_to(4.0, controls(dt), || Ok(())).unwrap().0;
        [CellId(10), CellId(11)]
            .into_iter()
            .map(|id| {
                (candidate.temperature_k(id).unwrap() - reference.temperature_k(id).unwrap()).abs()
            })
            .sum::<f64>()
    });
    assert!(
        errors[0] > 1.7 * errors[1] && errors[1] > 1.7 * errors[2],
        "{errors:?}"
    );
    assert_eq!(start.time_s(), 0.0);
    assert_eq!(start.cells()[&CellId(11)], original_mold);
}

#[test]
fn ordinary_mold_caloric_domain_failure_is_atomic() {
    let law = law();
    let mut start =
        CrystallizationSimulation::new(finite_bath(), &[material(&law)], || Ok(())).unwrap();
    let mut heated = boundary(0.0);
    heated.source_w = 1e12;
    start.set_boundary(CellId(11), heated, || Ok(())).unwrap();
    let original_cells = start.cells().clone();
    let original_states = start.states().clone();
    assert!(matches!(
        start.advance_to(1.0, controls(0.1), || Ok(())),
        Err(CrystallizationTimelineError::Thermal(
            DepositionError::MaterialDomain
        )),
    ));
    assert_eq!(start.cells(), &original_cells);
    assert_eq!(start.states(), &original_states);
    assert_eq!(start.time_s(), 0.0);
}

#[test]
fn cloned_checkpoint_preserves_exposure_and_restarts_bit_for_bit() {
    let law = law();
    let start = CrystallizationSimulation::new(
        thermal(&[380.0, 390.0], 20.0, 0.5),
        &[material(&law)],
        || Ok(()),
    )
    .unwrap();
    let all = start
        .advance_to(4.0, controls(0.0625), || Ok(()))
        .unwrap()
        .0;
    let half = start
        .advance_to(2.0, controls(0.0625), || Ok(()))
        .unwrap()
        .0;
    let resumed = half
        .clone()
        .advance_to(4.0, controls(0.0625), || Ok(()))
        .unwrap()
        .0;
    assert_eq!(all.states(), resumed.states());
    assert_eq!(all.cells(), resumed.cells());
    assert_eq!(all.time_s(), resumed.time_s());
}

#[test]
fn cancellation_and_domain_failures_leave_all_original_state_unchanged() {
    let law = law();
    let mut start = CrystallizationSimulation::new(
        thermal(&[380.0, 390.0], 20.0, 0.5),
        &[material(&law)],
        || Ok(()),
    )
    .unwrap();
    let original_states = start.states.clone();
    let original_cells = start.cells().clone();
    let mut calls = 0;
    let failed = start.advance_to(4.0, controls(0.0625), || {
        calls += 1;
        if calls == 80 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        failed,
        Err(CrystallizationTimelineError::Thermal(
            DepositionError::Cancelled
        ))
    ));
    assert_eq!(start.states, original_states);
    assert_eq!(start.cells(), &original_cells);
    assert_eq!(start.time_s(), 0.0);
    let mut heated = boundary(0.0);
    heated.source_w = 1e12;
    start.set_boundary(CellId(10), heated, || Ok(())).unwrap();
    assert!(start.advance_to(1.0, controls(0.1), || Ok(())).is_err());
    assert_eq!(start.states, original_states);
    assert_eq!(start.time_s(), 0.0);
}

#[test]
fn attachment_preserves_an_original_account_despite_caloric_roundoff() {
    let law = law();
    let mut original = thermal(&[380.0], 0.0, 0.0);
    let initial_h = 75_000.123_456_789;
    original.cells.get_mut(&CellId(10)).unwrap().specific_h_j_kg = initial_h;
    let original_energy = original.energy_j();
    let round_tripped = CrystallizationCaloric::specific_enthalpy_j_kg(
        &AMORPHOUS,
        original.temperature_k(CellId(10)).unwrap(),
    )
    .unwrap();
    assert_ne!(round_tripped, initial_h);
    let adopted = CrystallizationSimulation::new(original, &[material(&law)], || Ok(())).unwrap();
    assert_eq!(adopted.energy_j(), original_energy);
    assert_eq!(
        adopted.states()[&CellId(10)].specific_enthalpy_j_kg(),
        initial_h
    );
    assert_eq!(adopted.cells()[&CellId(10)].specific_h_j_kg, initial_h);
}

#[test]
fn wrong_source_binding_and_exhausted_work_budget_refuse() {
    let law = law();
    let other = Branch {
        heat_capacity_j_kg_k: 1000.0,
        offset: 0.0,
    };
    let mut wrong = material(&law);
    wrong.amorphous = &other;
    assert!(
        CrystallizationSimulation::new(thermal(&[380.0], 0.0, 0.0), &[wrong], || Ok(())).is_err()
    );
    let mut original = thermal(&[380.0], 0.0, 0.0);
    original.budget.max_substeps = 1;
    let start = CrystallizationSimulation::new(original, &[material(&law)], || Ok(())).unwrap();
    assert!(matches!(
        start.advance_to(2.0, controls(0.1), || Ok(())),
        Err(CrystallizationTimelineError::Thermal(
            DepositionError::Budget(_)
        ))
    ));
    assert_eq!(start.time_s(), 0.0);
}

#[test]
fn crossing_parent_calorics_refuse_even_with_zero_crystallization_rate() {
    let crossing = Branch {
        heat_capacity_j_kg_k: 2000.0,
        offset: 0.0,
    };
    let zero_rate = |_: f64| Ok(0.0);
    let wrong_parent =
        CrystallizationLaw::new(&crossing, &CRYSTALLINE, &zero_rate, law().config()).unwrap();
    let original = thermal(&[300.0], 0.0, 0.0);
    // A one-point attachment check cannot distinguish these different laws.
    assert_eq!(original.temperature_k(CellId(10)).unwrap(), 300.0);
    assert_eq!(original.cells()[&CellId(10)].specific_h_j_kg, 0.0);
    assert_eq!(
        wrong_parent
            .state_at_temperature(300.0, 0.0)
            .unwrap()
            .specific_enthalpy_j_kg(),
        0.0
    );
    assert_ne!(
        ThermalLaw::caloric_identity(&AMORPHOUS),
        wrong_parent.parent_caloric_identity(),
    );
    assert!(
        CrystallizationSimulation::new(original, &[material(&wrong_parent)], || Ok(())).is_err()
    );
}

#[test]
fn attachment_requires_identified_original_and_kinetic_parent_calorics() {
    let identified = law();
    let zero_rate = |_: f64| Ok(0.0);
    let unidentified =
        CrystallizationLaw::new(&UNIDENTIFIED, &CRYSTALLINE, &zero_rate, identified.config())
            .unwrap();
    let cell = thermal(&[300.0], 0.0, 0.0).cells()[&CellId(10)].clone();
    let cases: [(&'static dyn ThermalLaw, &CrystallizationLaw<'_>); 3] = [
        (&UNIDENTIFIED, &identified),
        (&AMORPHOUS, &unidentified),
        (&UNIDENTIFIED, &unidentified),
    ];
    for (original, law) in cases {
        let thermal = activate(
            vec![TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell.clone()),
            }],
            &[(MaterialId(3), original)],
        );
        let mut binding = material(law);
        binding.amorphous = original;
        assert!(CrystallizationSimulation::new(thermal, &[binding], || Ok(())).is_err());
    }
}

fn scheduled(
    events: Vec<TimedEvent>,
    materials: &[(MaterialId, &'static dyn ThermalLaw)],
) -> DepositionSimulation<'static> {
    DepositionSimulation::new(
        events,
        materials,
        DepositionBudget {
            max_cells: 16,
            max_links: 16,
            max_events: 32,
            max_substeps: 100_000,
        },
    )
    .unwrap()
}

#[test]
fn delayed_births_keep_their_own_kinetic_age_and_conserve_incoming_enthalpy() {
    // K is constant over every temperature reached by this calorimetry test,
    // then continuously vanishes at the original 500 K melting boundary.
    let rate = |t: f64| Ok(0.025 * (500.0 - t).clamp(0.0, 10.0));
    let law = CrystallizationLaw::new(&AMORPHOUS, &CRYSTALLINE, &rate, law().config()).unwrap();
    let template = thermal(&[380.0, 380.0], 0.0, 0.0);
    let events = [2.0, 4.0]
        .into_iter()
        .enumerate()
        .map(|(i, time_s)| TimedEvent {
            time_s,
            event: Event::Birth(template.cells()[&CellId(i as u64 + 10)].clone()),
        })
        .collect();
    let start = CrystallizationSimulation::new(
        scheduled(events, &[(MaterialId(3), &AMORPHOUS)]),
        &[material(&law)],
        || Ok(()),
    )
    .unwrap();
    let (empty, idle_audit) = start.advance_to(1.0, controls(0.125), || Ok(())).unwrap();
    assert!(empty.cells().is_empty());
    assert_eq!(idle_audit.thermal, empty_audit());
    let (born, born_audit) = start.advance_to(2.0, controls(0.125), || Ok(())).unwrap();
    assert_eq!(born_audit.thermal.substeps, 0);
    assert_eq!(born_audit.thermal.events, 1);
    assert_eq!(born_audit.thermal.incoming_mass_kg, 0.01);
    assert_eq!(born_audit.thermal.incoming_energy_j, 800.0);
    assert_eq!(born_audit.thermal.stored_energy_change_j, 800.0);
    assert_eq!(born.states()[&CellId(10)].exposure(), 0.0);
    assert!(!born.cells().contains_key(&CellId(11)));

    let (middle, _) = start.advance_to(4.0, controls(0.125), || Ok(())).unwrap();
    assert_eq!(middle.states()[&CellId(10)].exposure(), 0.5);
    assert_eq!(middle.states()[&CellId(11)].exposure(), 0.0);
    let (next, audit) = start.advance_to(6.0, controls(0.125), || Ok(())).unwrap();
    for (id, age) in [(CellId(10), 4.0_f64), (CellId(11), 2.0_f64)] {
        let exposure = 0.25 * age;
        let fraction = 0.8 * -(-(exposure * exposure)).exp_m1();
        assert_eq!(next.states()[&id].exposure(), exposure);
        close(
            next.states()[&id].crystalline_mass_fraction(),
            fraction,
            2e-15,
        );
        close(
            next.temperature_k(id).unwrap(),
            380.0 + 100.0 * fraction,
            1e-8,
        );
        assert_eq!(next.cells()[&id].specific_h_j_kg, 80_000.0);
    }
    assert_eq!(audit.thermal.events, 2);
    assert_eq!(audit.thermal.incoming_mass_kg, 0.02);
    assert_eq!(audit.thermal.incoming_energy_j, 1600.0);
    assert_eq!(audit.thermal.stored_energy_change_j, 1600.0);
    assert_eq!(audit.thermal.external_heat_j, 0.0);
    assert_eq!(audit.thermal.residual_j, 0.0);
    close(
        audit.crystallized_mass_kg,
        next.states()
            .values()
            .map(|state| 0.01 * state.crystalline_mass_fraction())
            .sum(),
        1e-16,
    );
    let resumed = middle
        .clone()
        .advance_to(6.0, controls(0.125), || Ok(()))
        .unwrap()
        .0;
    assert_eq!(next.cells(), resumed.cells());
    assert_eq!(next.states(), resumed.states());
    assert!(next.pending_events().is_empty());
    assert_eq!(start.pending_events().len(), 2);
    assert!(start.cells().is_empty());
}

#[test]
fn scheduled_contact_separation_and_baths_use_the_original_event_intervals() {
    let law = law();
    let template = finite_bath();
    let events = vec![
        TimedEvent {
            time_s: 0.0,
            event: Event::Birth(template.cells()[&CellId(10)].clone()),
        },
        TimedEvent {
            time_s: 0.0,
            event: Event::Birth(template.cells()[&CellId(11)].clone()),
        },
        TimedEvent {
            time_s: 1.0,
            event: Event::Contact(template.links()[&1].clone()),
        },
        TimedEvent {
            time_s: 2.0,
            event: Event::Separation(1),
        },
        TimedEvent {
            time_s: 2.0,
            event: Event::Boundary(CellId(10), boundary(20.0)),
        },
    ];
    let start = CrystallizationSimulation::new(
        scheduled(
            events,
            &[(MaterialId(3), &AMORPHOUS), (MaterialId(4), &MOLD)],
        ),
        &[material(&law)],
        || Ok(()),
    )
    .unwrap();
    let (contact, first) = start.advance_to(1.0, controls(0.125), || Ok(())).unwrap();
    assert_eq!(first.thermal.events, 3);
    assert_eq!(contact.links().len(), 1);
    assert_eq!(contact.temperature_k(CellId(11)).unwrap(), 300.0);
    assert!(contact.temperature_k(CellId(10)).unwrap() > 380.0);
    assert!(!contact.states().contains_key(&CellId(11)));
    let (separated, second) = contact.advance_to(2.0, controls(0.125), || Ok(())).unwrap();
    assert!(separated.links().is_empty());
    assert_eq!(second.thermal.events, 2);
    assert_eq!(second.thermal.external_heat_j, 0.0);
    close(separated.energy_j(), contact.energy_j(), 1e-8);
    assert!(separated.temperature_k(CellId(11)).unwrap() > 300.0);
    let mold_h = separated.cells()[&CellId(11)].specific_h_j_kg;
    let (cooled, third) = separated
        .advance_to(2.5, controls(0.125), || Ok(()))
        .unwrap();
    assert_eq!(cooled.cells()[&CellId(11)].specific_h_j_kg, mold_h);
    assert!(third.thermal.external_heat_j < 0.0);
    close(
        cooled.energy_j() - separated.energy_j(),
        third.thermal.external_heat_j,
        1e-8,
    );
    let (whole, audit) = start.advance_to(2.5, controls(0.125), || Ok(())).unwrap();
    assert_eq!(whole.cells(), cooled.cells());
    assert_eq!(whole.states(), cooled.states());
    assert_eq!(audit.thermal.events, 5);
    close(
        whole.energy_j(),
        audit.thermal.incoming_energy_j + audit.thermal.external_heat_j,
        1e-8,
    );
    close(audit.thermal.incoming_mass_kg, 0.06, 1e-17);
    assert!(audit.crystallized_mass_kg > 0.0);
}

#[test]
fn failed_scheduled_birth_and_cancellation_preserve_phase_history_and_cursor() {
    let full = law();
    let law = CrystallizationLaw::new(
        &AMORPHOUS,
        &CRYSTALLINE,
        &RATE,
        CrystallizationConfig {
            temperature_range_k: [250.0, 450.0],
            ..full.config()
        },
    )
    .unwrap();
    let template = thermal(&[380.0, 490.0], 0.0, 0.0);
    let events = [0.0, 1.0]
        .into_iter()
        .enumerate()
        .map(|(i, time_s)| TimedEvent {
            time_s,
            event: Event::Birth(template.cells()[&CellId(i as u64 + 10)].clone()),
        })
        .collect();
    let start = CrystallizationSimulation::new(
        scheduled(events, &[(MaterialId(3), &AMORPHOUS)]),
        &[material(&law)],
        || Ok(()),
    )
    .unwrap()
    .advance_to(0.0, controls(0.125), || Ok(()))
    .unwrap()
    .0;
    assert!(matches!(
        start.advance_to(2.0, controls(0.125), || Ok(())),
        Err(CrystallizationTimelineError::Material(
            CrystallizationError::TemperatureOutOfRange { .. }
        ))
    ));
    assert_eq!(start.time_s(), 0.0);
    assert_eq!(start.cells().len(), 1);
    assert_eq!(start.states()[&CellId(10)].exposure(), 0.0);
    assert_eq!(start.pending_events().len(), 1);
    let mut polls = 0;
    assert!(matches!(
        start.advance_to(0.5, controls(0.125), || {
            polls += 1;
            if polls == 80 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(CrystallizationTimelineError::Thermal(
            DepositionError::Cancelled
        ))
    ));
    assert_eq!(start.time_s(), 0.0);
    assert_eq!(start.states()[&CellId(10)].exposure(), 0.0);
    assert_eq!(start.pending_events()[0].time_s, 1.0);
}

#[test]
fn schedule_cannot_renew_substep_budget_and_due_events_precede_manual_boundaries() {
    let law = law();
    let template = thermal(&[380.0, 390.0], 0.0, 0.0);
    let events = [0.0, 0.125]
        .into_iter()
        .enumerate()
        .map(|(i, time_s)| TimedEvent {
            time_s,
            event: Event::Birth(template.cells()[&CellId(i as u64 + 10)].clone()),
        })
        .collect();
    let mut thermal = scheduled(events, &[(MaterialId(3), &AMORPHOUS)]);
    thermal.budget.max_substeps = 1;
    let mut start = CrystallizationSimulation::new(thermal, &[material(&law)], || Ok(())).unwrap();
    assert!(
        start
            .set_boundary(CellId(10), boundary(10.0), || Ok(()))
            .is_err()
    );
    assert!(matches!(
        start.advance_to(0.25, controls(0.125), || Ok(())),
        Err(CrystallizationTimelineError::Thermal(
            DepositionError::Budget(_)
        ))
    ));
    assert_eq!(start.time_s(), 0.0);
    assert!(start.cells().is_empty());
    assert_eq!(start.pending_events().len(), 2);
    let mut born = start.advance_to(0.0, controls(0.125), || Ok(())).unwrap().0;
    born.set_boundary(CellId(10), boundary(10.0), || Ok(()))
        .unwrap();
    assert_eq!(born.cells()[&CellId(10)].boundary.htc_w_m2_k, 10.0);
}

#[test]
fn superheated_scheduled_feed_cools_into_crystallization_against_a_finite_mold() {
    // Identified manufactured affine calorics extend to 600 K on BOTH
    // original branches. No lower-domain source is silently extrapolated.
    struct HotBranch(f64);
    impl HotBranch {
        fn identity(&self) -> [u8; 32] {
            let mut id = [0x61; 32];
            id[..8].copy_from_slice(&self.0.to_bits().to_le_bytes());
            id
        }
    }
    impl ThermalLaw for HotBranch {
        fn caloric_identity(&self) -> Option<[u8; 32]> {
            Some(self.identity())
        }
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            let t = 300.0 + (h - self.0) / 1000.0;
            if (250.0..=600.0).contains(&t) {
                Ok(t)
            } else {
                Err(DepositionError::MaterialDomain)
            }
        }
        fn maximum_temperature_slope(&self) -> f64 {
            0.001
        }
    }
    impl CrystallizationCaloric for HotBranch {
        fn caloric_identity(&self) -> Option<[u8; 32]> {
            Some(self.identity())
        }
        fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, CrystallizationError> {
            if (250.0..=600.0).contains(&t) {
                Ok(1000.0 * (t - 300.0) + self.0)
            } else {
                Err(CrystallizationError::Invalid("hot caloric source domain"))
            }
        }
        fn minimum_heat_capacity_j_kg_k(&self) -> Option<f64> {
            Some(1000.0)
        }
    }
    let amorphous = HotBranch(0.0);
    let crystalline = HotBranch(-100_000.0);
    let law = CrystallizationLaw::new(
        &amorphous,
        &crystalline,
        &RATE,
        CrystallizationConfig {
            identity: [91; 32],
            temperature_range_k: [250.0, 600.0],
            ..law().config()
        },
    )
    .unwrap();
    let template = finite_bath();
    let mut feed = template.cells()[&CellId(10)].clone();
    feed.specific_h_j_kg = 250_000.0; // Original amorphous chart at 550 K.
    let original = DepositionSimulation::new(
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(template.cells()[&CellId(11)].clone()),
            },
            TimedEvent {
                time_s: 1.0,
                event: Event::Birth(feed),
            },
            TimedEvent {
                time_s: 1.0,
                event: Event::Contact(template.links()[&1].clone()),
            },
        ],
        &[(MaterialId(3), &amorphous), (MaterialId(4), &MOLD)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 1000,
        },
    )
    .unwrap();
    let start = CrystallizationSimulation::new(
        original,
        &[CrystallizationMaterial {
            material_id: MaterialId(3),
            amorphous: &amorphous,
            law: &law,
            maximum_temperature_slope: 0.001,
        }],
        || Ok(()),
    )
    .unwrap();
    let (hot, _) = start.advance_to(2.0, controls(0.125), || Ok(())).unwrap();
    assert!(hot.temperature_k(CellId(10)).unwrap() > 500.0);
    assert_eq!(hot.states()[&CellId(10)].exposure(), 0.0);
    assert_eq!(hot.states()[&CellId(10)].crystalline_mass_fraction(), 0.0);
    let (cooled, audit) = start.advance_to(20.0, controls(0.125), || Ok(())).unwrap();
    assert!(cooled.temperature_k(CellId(10)).unwrap() < 500.0);
    assert!(cooled.temperature_k(CellId(11)).unwrap() > 300.0);
    assert!(cooled.states()[&CellId(10)].crystalline_mass_fraction() > 0.0);
    assert_eq!(audit.thermal.incoming_energy_j, 2500.0);
    assert_eq!(audit.thermal.external_heat_j, 0.0);
    close(cooled.energy_j(), 2500.0, 1e-7);
    close(cooled.mass_kg(), 0.06, 1e-17);
    assert!(!cooled.states().contains_key(&CellId(11)));
    assert!(start.cells().is_empty());
}
