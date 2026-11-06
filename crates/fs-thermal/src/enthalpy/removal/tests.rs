use super::*;
use crate::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use crate::enthalpy::{
    Boundary, ContactPatch, DepositionBudget, Event, StepConfig, ThermalLaw, ThermalLink,
    TimedEvent,
};

struct Caloric {
    datum: f64,
    basis: ThermalStorageBasis,
}
impl ThermalLaw for Caloric {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.basis
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        let mut identity = [41; 32];
        identity[..8].copy_from_slice(&self.datum.to_le_bytes());
        identity[8] = self.basis as u8;
        Some(identity)
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + (h - self.datum) / 2000.0;
        if (200.0..=600.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.0005
    }
}

fn config() -> StepConfig {
    StepConfig {
        max_dt_s: 1.0,
        energy_tolerance_j: 1.0e-8,
    }
}
fn tolerances() -> MaterialRemovalTolerances {
    MaterialRemovalTolerances {
        mass_kg: 1.0e-12,
        volume_m3: 1.0e-15,
    }
}
fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: 350.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    }
}

fn state<'a>(law: &'a Caloric, source_w: f64, extra: Vec<TimedEvent>) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for id in 1..=3 {
        let length = if id == 1 { 0.02 } else { 0.01 };
        let (start, end) = match id {
            1 => ([-0.02, 0.0, 0.0], [0.0; 3]),
            2 => ([0.0; 3], [0.01, 0.0, 0.0]),
            _ => ([1.0, 0.0, 0.0], [1.01, 0.0, 0.0]),
        };
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: id,
                material_id: MaterialId(1),
                mass_kg: 1000.0 * 0.1 * length,
                specific_h_j_kg: law.datum + 100_000.0,
                geometry: CellGeometry {
                    start_m: start,
                    end_m: end,
                    volume_m3: 0.1 * length,
                    surface_area_m2: 0.2 + 2.0 * length,
                },
                boundary: boundary(if id == 1 { source_w } else { 0.0 }),
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: CellId(1),
            b: CellId(2),
            area_m2: 0.1,
            conductance_w_k: 2.0,
        }),
    });
    events.extend(extra);
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 3,
            max_links: 4,
            max_events: 16,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0
}

fn removal(s: &DepositionSimulation<'_>, id: u64, volume: f64) -> MaterialRemoval {
    let expected = s.cells()[&CellId(id)].clone();
    let mut geometry = expected.geometry;
    geometry.volume_m3 -= volume;
    geometry.surface_area_m2 -= 2.0 * volume / 0.1;
    if id == 1 {
        geometry.start_m[0] += volume / 0.1;
    } else {
        geometry.end_m[0] -= volume / 0.1;
    }
    MaterialRemoval {
        caloric_identity: s
            .material_caloric_identity(expected.material_id)
            .unwrap()
            .unwrap(),
        expected,
        removed_volume_m3: volume,
        density_kg_m3: 1000.0,
        replacement_geometry: geometry,
    }
}
fn patch(s: &DepositionSimulation<'_>, conductance: f64) -> ContactPatchUpdate {
    let old = s.contact_patch(7).unwrap();
    ContactPatchUpdate {
        id: 7,
        expected: Some(old),
        replacement: Some(ContactPatch {
            link: ThermalLink {
                conductance_w_k: conductance,
                ..old.link
            },
            ..old
        }),
    }
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2.0e-12 * a.abs().max(b.abs()).max(1.0e-12),
        "{a:e} != {b:e}"
    );
}

#[test]
fn actual_mass_volume_and_outgoing_enthalpy_close_and_change_the_next_heat_capacity() {
    let law = Caloric {
        datum: 0.0,
        basis: ThermalStorageBasis::Enthalpy,
    };
    let s = state(&law, 100.0, vec![]);
    let requests = [removal(&s, 1, 0.0002), removal(&s, 2, 0.0001)];
    let update = patch(&s, 3.0);
    let (next, audit) = s
        .remove_material_candidate(&requests, &[update], tolerances(), || Ok(()))
        .unwrap();
    close(next.cells()[&CellId(1)].mass_kg, 1.8);
    close(next.cells()[&CellId(2)].mass_kg, 0.9);
    close(next.cells()[&CellId(1)].geometry.volume_m3, 0.0018);
    close(next.cells()[&CellId(2)].geometry.volume_m3, 0.0009);
    close(audit.removed_mass_kg, 0.3);
    close(audit.removed_volume_m3, 0.0003);
    close(audit.outgoing_enthalpy_j, 30_000.0);
    close(audit.retained_energy_change_j, -30_000.0);
    close(next.mass_kg() + audit.removed_mass_kg, s.mass_kg());
    close(next.energy_j() + audit.outgoing_enthalpy_j, s.energy_j());
    assert_eq!(audit.mass_balance_residual_kg, 0.0);
    assert_eq!(audit.energy_residual_j, 0.0);
    assert_eq!(next.cells()[&CellId(3)], s.cells()[&CellId(3)]);
    assert_eq!(next.contact_patch(7), update.replacement);
    assert_eq!(next.time_s(), s.time_s());
    assert_eq!(next.pending_events(), s.pending_events());
    for id in [CellId(1), CellId(2)] {
        assert_eq!(
            next.cells()[&id].specific_h_j_kg.to_bits(),
            s.cells()[&id].specific_h_j_kg.to_bits()
        );
        assert_eq!(next.temperature_k(id).unwrap(), 350.0);
    }
    for receipt in &audit.cells {
        assert_eq!(receipt.caloric_identity, law.caloric_identity().unwrap());
        assert!(receipt.mass_roundoff_kg.abs() <= tolerances().mass_kg);
        assert!(receipt.volume_roundoff_m3.abs() <= tolerances().volume_m3);
        assert!(receipt.maximum_density_residual_kg <= tolerances().mass_kg);
    }
    let (heated, heat) = next.advance_to(0.1, config(), || Ok(())).unwrap();
    assert_eq!(heat.substeps, 1);
    close(heat.external_heat_j, 10.0);
    close(
        heated.temperature_k(CellId(1)).unwrap(),
        350.0 + 10.0 / (1.8 * 2000.0),
    );
    let (unworn, _) = s.advance_to(0.1, config(), || Ok(())).unwrap();
    assert!(heated.temperature_k(CellId(1)).unwrap() > unworn.temperature_k(CellId(1)).unwrap());
    assert_eq!(s.cells()[&CellId(1)].mass_kg, 2.0);

    let (permuted, other) = s
        .remove_material_candidate(
            &[requests[1].clone(), requests[0].clone()],
            &[update],
            tolerances(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(other, audit);
    assert_eq!(permuted.cells(), next.cells());
}

#[test]
fn positive_and_negative_caloric_datums_change_only_the_carried_energy_coordinate() {
    let mut reference = None;
    for datum in [0.0, 36_028_797_018_963_968.0, -36_028_797_018_963_968.0] {
        let law = Caloric {
            datum,
            basis: ThermalStorageBasis::Enthalpy,
        };
        let s = state(&law, 0.0, vec![]);
        let request = removal(&s, 1, 0.0002);
        let (next, audit) = s
            .remove_material_candidate(&[request], &[patch(&s, 2.5)], tolerances(), || Ok(()))
            .unwrap();
        let receipt = &audit.cells[0];
        assert_eq!(
            next.temperature_k(CellId(1)).unwrap(),
            s.temperature_k(CellId(1)).unwrap()
        );
        assert_eq!(
            next.cells()[&CellId(1)].specific_h_j_kg.to_bits(),
            s.cells()[&CellId(1)].specific_h_j_kg.to_bits()
        );
        assert_eq!(audit.retained_energy_change_j, -audit.outgoing_enthalpy_j);
        assert_eq!(audit.energy_residual_j, 0.0);
        assert_eq!(audit.mass_balance_residual_kg, 0.0);
        if let Some((mass, volume, outgoing)) = reference {
            assert_eq!(receipt.removed_mass_kg, mass);
            assert_eq!(receipt.removed_volume_m3, volume);
            close(receipt.outgoing_enthalpy_j - outgoing, mass * datum);
        } else {
            reference = Some((
                receipt.removed_mass_kg,
                receipt.removed_volume_m3,
                receipt.outgoing_enthalpy_j,
            ));
        }
        if datum < 0.0 {
            assert!(receipt.outgoing_enthalpy_j < 0.0);
        }
    }
}

#[test]
fn complete_contact_updates_density_and_late_refusals_are_atomic() {
    let law = Caloric {
        datum: 0.0,
        basis: ThermalStorageBasis::Enthalpy,
    };
    let s = state(&law, 0.0, vec![]);
    let saved = s.cells().clone();
    let old_patch = s.contact_patch(7);
    let requests = [removal(&s, 1, 0.0002), removal(&s, 2, 0.0001)];
    let updates = [patch(&s, 3.0)];
    assert!(
        s.remove_material_candidate(&requests, &[], tolerances(), || Ok(()))
            .is_err()
    );
    let mut stale = requests.clone();
    stale[1].expected.specific_h_j_kg += 1.0;
    assert!(
        s.remove_material_candidate(&stale, &updates, tolerances(), || Ok(()))
            .is_err()
    );
    let mut wrong_density = requests.clone();
    wrong_density[1].density_kg_m3 *= 0.5;
    assert!(
        s.remove_material_candidate(&wrong_density, &updates, tolerances(), || Ok(()))
            .is_err()
    );
    let mut wrong_law = requests.clone();
    wrong_law[1].caloric_identity = [0; 32];
    assert!(
        s.remove_material_candidate(&wrong_law, &updates, tolerances(), || Ok(()))
            .is_err()
    );
    let mut coverage = requests.clone();
    coverage[1].replacement_geometry.surface_area_m2 = 0.05;
    assert!(
        s.remove_material_candidate(&coverage, &updates, tolerances(), || Ok(()))
            .is_err()
    );
    let retained_patch = patch(&s, old_patch.unwrap().link.conductance_w_k);
    assert!(
        s.remove_material_candidate(&coverage, &[retained_patch], tolerances(), || Ok(()))
            .is_err()
    );
    let mut count = 0;
    s.remove_material_candidate(&requests, &updates, tolerances(), || {
        count += 1;
        Ok(())
    })
    .unwrap();
    let mut calls = 0;
    assert_eq!(
        s.remove_material_candidate(&requests, &updates, tolerances(), || {
            calls += 1;
            if calls == count {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        })
        .err(),
        Some(DepositionError::Cancelled)
    );
    assert_eq!(s.cells(), &saved);
    assert_eq!(s.contact_patch(7), old_patch);
    assert_eq!(s.time_s(), 0.0);
    assert!(
        s.remove_material_candidate(&requests, &updates, tolerances(), || Ok(()))
            .is_ok()
    );
}

#[test]
fn unsupported_storage_exhaustion_and_unrepresented_removal_refuse() {
    for basis in [
        ThermalStorageBasis::InternalEnergy,
        ThermalStorageBasis::DeclaredSensibleStorage,
    ] {
        let law = Caloric { datum: 0.0, basis };
        let s = state(&law, 0.0, vec![]);
        assert!(
            s.remove_material_candidate(
                &[removal(&s, 1, 0.0002)],
                &[patch(&s, 3.0)],
                tolerances(),
                || Ok(())
            )
            .is_err()
        );
    }
    let law = Caloric {
        datum: 0.0,
        basis: ThermalStorageBasis::Enthalpy,
    };
    let s = state(&law, 0.0, vec![]);
    let mut exhausted = removal(&s, 1, 0.0002);
    exhausted.removed_volume_m3 = exhausted.expected.geometry.volume_m3;
    assert!(
        s.remove_material_candidate(&[exhausted], &[patch(&s, 3.0)], tolerances(), || Ok(()))
            .is_err()
    );
    let tiny = removal(&s, 1, f64::MIN_POSITIVE);
    assert_eq!(
        s.remove_material_candidate(&[tiny], &[patch(&s, 3.0)], tolerances(), || Ok(()))
            .err(),
        Some(DepositionError::Arithmetic)
    );
    let no_removal = removal(&s, 1, 0.0);
    let (same, audit) = s
        .remove_material_candidate(&[no_removal], &[], tolerances(), || Ok(()))
        .unwrap();
    assert_eq!(same.cells(), s.cells());
    assert_eq!(same.links(), s.links());
    assert_eq!(audit.removed_mass_kg, 0.0);
    assert_eq!(audit.outgoing_enthalpy_j, 0.0);
}

#[test]
fn explicit_left_limit_preserves_pending_events_and_their_original_order() {
    let law = Caloric {
        datum: 0.0,
        basis: ThermalStorageBasis::Enthalpy,
    };
    let event = TimedEvent {
        time_s: 1.0,
        event: Event::Boundary(CellId(1), boundary(99.0)),
    };
    let s = state(&law, 0.0, vec![event.clone()]);
    let step = BackwardEulerStep::new(&s, 1.0, || Ok(())).unwrap();
    let trial = step
        .evaluate(
            &[0.0; 3],
            |_| {
                Ok::<_, DepositionError>(StepCoefficients {
                    sources_w: BTreeMap::new(),
                    conductances_w_k: BTreeMap::new(),
                    payload: (),
                })
            },
            || Ok(()),
        )
        .unwrap();
    let (left, (), _) = trial.accept(0.0, 0.0, || Ok(())).unwrap();
    let request = removal(&left, 1, 0.0002);
    let updates = [patch(&left, 3.0)];
    assert!(
        left.remove_material_candidate(&[request.clone()], &updates, tolerances(), || Ok(()))
            .is_err()
    );
    let (removed, _) = left
        .remove_material_before_pending_events(&[request], &updates, tolerances(), || Ok(()))
        .unwrap();
    assert_eq!(removed.time_s(), 1.0);
    assert_eq!(removed.pending_events(), &[event]);
    assert_eq!(removed.cells()[&CellId(1)].boundary.source_w, 0.0);
    let (after_event, _) = removed.advance_to(1.0, config(), || Ok(())).unwrap();
    assert_eq!(after_event.cells()[&CellId(1)].boundary.source_w, 99.0);
    assert_eq!(
        after_event.cells()[&CellId(1)].mass_kg,
        removed.cells()[&CellId(1)].mass_kg
    );
    assert_eq!(
        after_event.cells()[&CellId(1)].specific_h_j_kg,
        removed.cells()[&CellId(1)].specific_h_j_kg
    );
    assert!(after_event.pending_events().is_empty());
}
