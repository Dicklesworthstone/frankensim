use super::*;

#[test]
fn g1_heat_impulses_use_original_latent_storage_without_aging_or_recounting_sources() {
    let law = ReferencedPhase { datum: 0.0 };
    let initial = prepare(&law, vec![birth(1, 0.0, 80.0, 7.0), birth(2, 0.0, 140.0, 3.0)]);
    let heat = BTreeMap::from([(CellId(1), 24.0), (CellId(2), -24.0)]);
    let (next, audit) = initial.apply_heat_j(&heat, 0.0, 0.0, || Ok(())).unwrap();
    assert_eq!(next.cells()[&CellId(1)].specific_h_j_kg, 104.0);
    assert_eq!(next.cells()[&CellId(2)].specific_h_j_kg, 116.0);
    assert_eq!(next.temperature_k(CellId(1)).unwrap(), 364.0);
    assert_eq!(next.temperature_k(CellId(2)).unwrap(), 364.0);
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(audit.stored_energy_change_j, 0.0);
    assert_eq!(audit.residual_j, 0.0);
    assert_eq!((audit.events, audit.substeps), (0, 0));
    assert_eq!((next.time_s, next.cursor), (initial.time_s, initial.cursor));
    for id in [CellId(1), CellId(2)] {
        let mut expected = initial.cells()[&id].clone();
        expected.specific_h_j_kg = next.cells()[&id].specific_h_j_kg;
        assert_eq!(next.cells()[&id], expected);
    }
    // Authored watts remain active only when physical time subsequently advances.
    let (_, later) = next.advance_to(1.0, config(), || Ok(())).unwrap();
    assert_eq!(later.external_heat_j, 10.0);
    assert_eq!(initial.cells()[&CellId(1)].specific_h_j_kg, 80.0);
}

#[test]
fn g0_heat_impulse_balance_resolves_cell_changes_below_aggregate_energy_precision() {
    let datum = 2.0_f64.powi(26);
    let quantum_j = 2.0_f64.powi(-26);
    let law = ReferencedPhase { datum };
    let initial = prepare(&law, vec![
        birth(1, 0.0, datum + 80.0, 0.0), birth(2, 0.0, datum + 80.0, 0.0),
    ]);
    let (next, audit) = initial.apply_heat_j(
        &BTreeMap::from([(CellId(1), quantum_j)]), 0.0, 0.0, || Ok(()),
    ).unwrap();
    assert_eq!(next.energy_j(), initial.energy_j());
    assert_eq!(audit.external_heat_j, quantum_j);
    assert_eq!(audit.stored_energy_change_j, quantum_j);
    assert_eq!(audit.residual_j, 0.0);
    assert_eq!(next.temperature_k(CellId(1)).unwrap(), 364.0);
}

#[test]
fn g0_heat_impulses_gate_local_and_aggregate_absolute_storage_errors() {
    let datum = 2.0_f64.powi(26);
    let lost_j = 2.0_f64.powi(-28);
    let law = ReferencedPhase { datum };
    let initial = prepare(&law, vec![
        birth(1, 0.0, datum + 80.0, 0.0), birth(2, 0.0, datum + 80.0, 0.0),
    ]);
    let heat = BTreeMap::from([(CellId(1), lost_j), (CellId(2), -lost_j)]);
    // Both signs disappear from h. Zero global signed error cannot excuse them.
    assert!(matches!(initial.apply_heat_j(&heat, lost_j, lost_j, || Ok(())),
        Err(DepositionError::EnergyBalance(_))));
    assert!(matches!(initial.apply_heat_j(&heat, 0.5 * lost_j, 8.0 * lost_j, || Ok(())),
        Err(DepositionError::EnergyBalance(_))));
    let (next, audit) = initial.apply_heat_j(&heat, lost_j, 2.0 * lost_j, || Ok(())).unwrap();
    assert_eq!(next.cells(), initial.cells());
    assert_eq!(audit.residual_j, 0.0);
}

#[test]
fn g4_heat_impulses_refuse_domain_errors_unknown_cells_and_late_cancellation_atomically() {
    let law = ReferencedPhase { datum: 0.0 };
    let initial = prepare(&law, vec![birth(1, 0.0, 80.0, 0.0), birth(2, 0.0, 80.0, 0.0)]);
    for heat in [
        BTreeMap::from([(CellId(1), 2.0), (CellId(2), 200.0)]),
        BTreeMap::from([(CellId(9), 0.0)]),
        BTreeMap::from([(CellId(1), f64::INFINITY)]),
    ] {
        assert!(initial.apply_heat_j(&heat, 0.0, 0.0, || Ok(())).is_err());
    }
    let heat = BTreeMap::from([(CellId(1), 2.0)]);
    for allowed in [0, 1, 2] {
        let mut calls = 0;
        assert!(matches!(initial.apply_heat_j(&heat, 0.0, 0.0, || {
            calls += 1;
            if calls > allowed { Err(DepositionError::Cancelled) } else { Ok(()) }
        }), Err(DepositionError::Cancelled)));
    }
    for tolerance in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(initial.apply_heat_j(&heat, tolerance, 0.0, || Ok(())).is_err());
        assert!(initial.apply_heat_j(&heat, 0.0, tolerance, || Ok(())).is_err());
    }
    assert_eq!(initial.cells()[&CellId(1)].specific_h_j_kg, 80.0);
    assert_eq!(initial.cells()[&CellId(2)].specific_h_j_kg, 80.0);
    let (same, audit) = initial.apply_heat_j(&BTreeMap::new(), 0.0, 0.0, || Ok(())).unwrap();
    assert_eq!(same.cells(), initial.cells());
    assert_eq!(audit.external_heat_j, 0.0);
}

#[test]
fn g0_heat_impulses_require_same_time_births_and_boundary_events_to_be_processed_first() {
    let law = ReferencedPhase { datum: 0.0 };
    let prepared = prepare(&law, vec![birth(1, 0.0, 80.0, 0.0)]);
    let pending = DepositionSimulation::new(
        vec![birth(1, 0.0, 80.0, 0.0)], &[(MaterialId(1), &law)], prepared.budget,
    ).unwrap();
    assert!(matches!(pending.apply_heat_j(&BTreeMap::new(), 0.0, 0.0, || Ok(())),
        Err(DepositionError::InvalidInput(_))));
    assert_eq!(pending.cursor, 0);
}

struct ReferencedPhase {
    datum: f64,
}
impl ThermalLaw for ReferencedPhase {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let sensible_and_latent = h - self.datum;
        if !(0.0..=256.0).contains(&sensible_and_latent) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if sensible_and_latent < 64.0 {
            300.0 + sensible_and_latent
        } else if sensible_and_latent <= 128.0 {
            364.0
        } else {
            364.0 + sensible_and_latent - 128.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}

fn config() -> StepConfig {
    StepConfig {
        max_dt_s: 1.0,
        energy_tolerance_j: 0.0,
    }
}

fn birth(id: u64, time_s: f64, h: f64, source_w: f64) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Birth(DepositedCell {
            id: CellId(id),
            segment_id: id,
            body_id: 1,
            material_id: MaterialId(1),
            mass_kg: 1.0,
            specific_h_j_kg: h,
            geometry: CellGeometry {
                start_m: [0.0, 0.0, id as f64],
                end_m: [1.0, 0.0, id as f64],
                volume_m3: 0.001,
                surface_area_m2: 1.0,
            },
            boundary: Boundary {
                ambient_k: 300.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w,
            },
        }),
    }
}

fn prepare(law: &dyn ThermalLaw, events: Vec<TimedEvent>) -> DepositionSimulation<'_> {
    let budget = DepositionBudget {
        max_cells: 8,
        max_links: 8,
        max_events: 16,
        max_substeps: 16,
    };
    DepositionSimulation::new(events, &[(MaterialId(1), law)], budget)
        .unwrap()
        .advance_to(0.0, config(), || Ok(()))
        .unwrap()
        .0
}

#[test]
fn g0_latent_heat_balance_is_invariant_to_a_large_enthalpy_reference() {
    // One cell increment is exactly representable, but the aggregate of two
    // reference-loaded cells has twice the ULP and loses the entire transfer.
    let quantum_j = 2.0_f64.powi(-26);
    for datum in [0.0, 2.0_f64.powi(26)] {
        let law = ReferencedPhase { datum };
        let initial = prepare(
            &law,
            vec![
                birth(1, 0.0, datum + 80.0, quantum_j),
                birth(2, 0.0, datum + 80.0, 0.0),
            ],
        );
        let (next, audit) = initial.advance_to(1.0, config(), || Ok(())).unwrap();
        assert_eq!(
            next.cells()[&CellId(1)].specific_h_j_kg - datum,
            80.0 + quantum_j
        );
        assert_eq!(next.cells()[&CellId(2)].specific_h_j_kg - datum, 80.0);
        assert_eq!(next.temperature_k(CellId(1)).unwrap(), 364.0);
        assert_eq!(audit.stored_energy_change_j, quantum_j);
        assert_eq!(audit.external_heat_j, quantum_j);
        assert_eq!(audit.residual_j, 0.0);
        if datum != 0.0 {
            assert_eq!(next.energy_j(), initial.energy_j());
        }
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.cells()[&CellId(1)].specific_h_j_kg, datum + 80.0);
    }
}

#[test]
fn g0_new_material_keeps_its_incoming_enthalpy_separate_from_subsequent_heat() {
    let datum = 2.0_f64.powi(26);
    let quantum_j = 2.0_f64.powi(-26);
    let law = ReferencedPhase { datum };
    let initial = prepare(
        &law,
        vec![
            birth(1, 0.0, datum + 80.0, 2.0 * quantum_j),
            birth(2, 0.5, datum + 80.0, 2.0 * quantum_j),
        ],
    );
    let (next, audit) = initial.advance_to(1.0, config(), || Ok(())).unwrap();
    assert_eq!(
        next.cells()[&CellId(1)].specific_h_j_kg - datum,
        80.0 + 2.0 * quantum_j
    );
    assert_eq!(
        next.cells()[&CellId(2)].specific_h_j_kg - datum,
        80.0 + quantum_j
    );
    assert_eq!(audit.incoming_mass_kg, 1.0);
    assert_eq!(audit.incoming_energy_j, datum + 80.0);
    assert_eq!(audit.external_heat_j, 3.0 * quantum_j);
    assert_eq!(audit.residual_j, 0.0);
    assert_eq!(audit.events, 1);
    assert_eq!(audit.substeps, 2);
}

#[derive(Clone)]
struct Exchange(f64);
impl DepositionObserver for Exchange {
    type Error = DepositionError;
    fn heat_transfers(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        Ok(vec![HeatTransfer {
            a: CellId(1),
            b: CellId(2),
            power_w: self.0,
        }])
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
        Ok(())
    }
}

#[test]
fn g0_canceling_but_unrepresentable_transfers_cannot_pass_the_global_balance() {
    let datum = 2.0_f64.powi(26);
    let law = ReferencedPhase { datum };
    let initial = prepare(
        &law,
        vec![
            birth(1, 0.0, datum + 80.0, 0.0),
            birth(2, 0.0, datum + 80.0, 0.0),
        ],
    );
    let result = initial.advance_observed(
        &Exchange(2.0_f64.powi(-28)),
        1.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 0.0,
        },
        || Ok(()),
    );
    assert!(matches!(result, Err(DepositionError::Arithmetic)));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells()[&CellId(1)].specific_h_j_kg, datum + 80.0);
    assert_eq!(initial.cells()[&CellId(2)].specific_h_j_kg, datum + 80.0);
}

#[test]
fn g0_roundoff_sized_transfers_respect_the_declared_energy_tolerance() {
    let datum = 2.0_f64.powi(26);
    let power = 2.0_f64.powi(-28);
    let law = ReferencedPhase { datum };
    let initial = prepare(
        &law,
        vec![
            birth(1, 0.0, datum + 80.0, 0.0),
            birth(2, 0.0, datum + 80.0, 0.0),
        ],
    );
    let (next, _, audit) = initial
        .advance_observed(
            &Exchange(power),
            1.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 8.0 * power,
            },
            || Ok(()),
        )
        .unwrap();
    assert_eq!(next.time_s(), 1.0);
    assert_eq!(next.cells()[&CellId(1)].specific_h_j_kg, datum + 80.0);
    assert_eq!(next.cells()[&CellId(2)].specific_h_j_kg, datum + 80.0);
    assert_eq!(audit.residual_j, 0.0);
    assert_eq!(audit.external_heat_j, 0.0);
}

#[test]
fn g0_unrepresented_transfer_budget_accumulates_across_substeps() {
    let datum = 2.0_f64.powi(26);
    let power = 2.0_f64.powi(-28);
    let law = ReferencedPhase { datum };
    let initial = prepare(
        &law,
        vec![
            birth(1, 0.0, datum + 80.0, 0.0),
            birth(2, 0.0, datum + 80.0, 0.0),
        ],
    );
    // Each interval loses 2*power joules locally, while its net heat is zero.
    // One interval fits the tolerance; two must refuse the entire candidate.
    let result = initial.advance_observed(
        &Exchange(power),
        2.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 3.0 * power,
        },
        || Ok(()),
    );
    assert!(matches!(result, Err(DepositionError::Arithmetic)));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells()[&CellId(1)].specific_h_j_kg, datum + 80.0);
}
