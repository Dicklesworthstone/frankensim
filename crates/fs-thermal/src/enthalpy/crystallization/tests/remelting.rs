//! Heating and later hot deposition remelt the existing kinetic crystal mass.
use super::*;

const MAT: MaterialId = MaterialId(3);
const OLD: CellId = CellId(10);
const NEW: CellId = CellId(11);

struct CycleBranch(f64);
impl CycleBranch {
    fn identity(&self) -> [u8; 32] {
        let mut id = [0x74; 32];
        id[..8].copy_from_slice(&self.0.to_bits().to_le_bytes());
        id
    }
}
impl CrystallizationCaloric for CycleBranch {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some(self.identity())
    }
    fn minimum_heat_capacity_j_kg_k(&self) -> Option<f64> {
        Some(1000.0)
    }
    fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, CrystallizationError> {
        if !(250.0..=650.0).contains(&t) {
            return Err(CrystallizationError::Invalid(
                "cycle source temperature domain",
            ));
        }
        Ok(1000.0 * (t - 300.0) + self.0)
    }
}
impl ThermalLaw for CycleBranch {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some(self.identity())
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.001
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + (h - self.0) / 1000.0;
        if (250.0..=650.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
}
static PARENT: CycleBranch = CycleBranch(0.0);
static CRYSTAL: CycleBranch = CycleBranch(-100_000.0);

fn cycle_law(upper: f64) -> CrystallizationLaw<'static> {
    CrystallizationLaw::new(
        &PARENT,
        &CRYSTAL,
        &RATE,
        CrystallizationConfig {
            identity: [0x75; 32],
            temperature_range_k: [250.0, upper],
            ..law().config()
        },
    )
    .unwrap()
}
fn birth(id: CellId, t: f64, mass: f64) -> Event {
    Event::Birth(DepositedCell {
        id,
        segment_id: id.0,
        body_id: id.0,
        material_id: MAT,
        mass_kg: mass,
        specific_h_j_kg: 1000.0 * (t - 300.0),
        geometry: CellGeometry {
            start_m: [id.0 as f64 * 0.01, 0.0, 0.0],
            end_m: [(id.0 + 1) as f64 * 0.01, 0.0, 0.0],
            volume_m3: mass / 1000.0,
            surface_area_m2: 0.01,
        },
        boundary: boundary(0.0),
    })
}
fn prepare<'a>(
    law: &'a CrystallizationLaw<'a>,
    mut later: Vec<TimedEvent>,
) -> CrystallizationSimulation<'a> {
    let mut events = vec![TimedEvent {
        time_s: 0.0,
        event: birth(OLD, 380.0, 0.01),
    }];
    events.append(&mut later);
    let thermal = DepositionSimulation::new(
        events,
        &[(MAT, &PARENT as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 4,
            max_substeps: 10_000,
        },
    )
    .unwrap();
    CrystallizationSimulation::new(
        thermal,
        &[CrystallizationMaterial {
            material_id: MAT,
            amorphous: &PARENT,
            law,
            maximum_temperature_slope: 0.001,
        }],
        || Ok(()),
    )
    .unwrap()
}
fn source(state: &mut CrystallizationSimulation<'_>, power: f64) {
    let mut value = boundary(0.0);
    value.source_w = power;
    state.set_boundary(OLD, value, || Ok(())).unwrap();
}

#[test]
fn heating_absorbs_latent_energy_then_superheats_and_recooling_recrystallizes() {
    let law = cycle_law(650.0);
    let remelting = law.equivalent_age_remelting().unwrap();
    let start = prepare(&law, vec![]);
    let selected = start
        .clone()
        .with_equilibrium_remelting(MAT, remelting)
        .unwrap();
    let (mut aged, _) = selected
        .advance_to(20.0, controls(0.02), || Ok(()))
        .unwrap();
    let original_crystals = aged.states()[&OLD].crystalline_mass_fraction();
    assert!(original_crystals > 0.3);
    assert!(aged.temperature_k(OLD).unwrap() < 500.0);
    close(aged.energy_j(), 800.0, 1e-10);
    source(&mut aged, 100.0);

    // The original path still refuses unsupported remelting at the same heat.
    let (mut irreversible, _) = start.advance_to(20.0, controls(0.02), || Ok(())).unwrap();
    source(&mut irreversible, 100.0);
    assert!(
        irreversible
            .advance_to(37.0, controls(0.02), || Ok(()))
            .is_err()
    );
    assert_eq!(irreversible.time_s(), 20.0);

    // M=0.01 kg, cp=1000 and L=100000. Nine seconds of 100 W
    // puts h=170000 J/kg: at Tm=500, exactly X=0.3 remains.
    let (partial, audit) = aged.advance_to(29.0, controls(0.02), || Ok(())).unwrap();
    close(partial.cells()[&OLD].specific_h_j_kg, 170_000.0, 2e-6);
    close(partial.temperature_k(OLD).unwrap(), 500.0, 1e-10);
    close(
        partial.states()[&OLD].crystalline_mass_fraction(),
        0.3,
        3e-11,
    );
    close(audit.thermal.external_heat_j, 900.0, 2e-8);
    close(audit.thermal.stored_energy_change_j, 900.0, 2e-8);
    close(
        audit.crystallized_mass_kg,
        0.01 * (0.3 - original_crystals),
        4e-13,
    );
    assert!(audit.crystallized_mass_kg < 0.0);
    assert_eq!(partial.mass_kg(), aged.mass_kg());
    assert_eq!(partial.remelting_model(MAT), Some(remelting));

    let mut held = partial.clone();
    source(&mut held, 0.0);
    let (held_next, held_audit) = held.advance_to(30.0, controls(0.02), || Ok(())).unwrap();
    assert_eq!(held_next.states(), held.states());
    assert_eq!(held_audit.crystallized_mass_kg, 0.0);
    assert_eq!(held_audit.thermal.external_heat_j, 0.0);

    let (mut melted, audit) = partial.advance_to(37.0, controls(0.02), || Ok(())).unwrap();
    let replay = partial.advance_to(37.0, controls(0.02), || Ok(())).unwrap();
    assert_eq!(melted.cells(), replay.0.cells());
    assert_eq!(melted.states(), replay.0.states());
    assert_eq!(audit, replay.1);
    close(melted.temperature_k(OLD).unwrap(), 550.0, 3e-9);
    assert_eq!(melted.states()[&OLD].crystalline_mass_fraction(), 0.0);
    assert_eq!(melted.states()[&OLD].exposure(), 0.0);
    close(audit.thermal.external_heat_j, 800.0, 2e-8);
    close(audit.crystallized_mass_kg, -0.003, 4e-13);

    source(&mut melted, -100.0);
    let (recooled, audit) = melted.advance_to(52.0, controls(0.02), || Ok(())).unwrap();
    assert!(recooled.temperature_k(OLD).unwrap() < 500.0);
    assert!(recooled.states()[&OLD].crystalline_mass_fraction() > 0.0);
    assert!(audit.crystallized_mass_kg > 0.0);
    close(recooled.energy_j(), 1000.0, 4e-8);
    close(audit.thermal.external_heat_j, -1500.0, 2e-8);
    assert_eq!(recooled.mass_kg(), aged.mass_kg());
}

#[test]
fn a_later_hot_deposit_remelts_older_crystals_through_original_finite_contact() {
    let law = cycle_law(650.0);
    let remelting = law.equivalent_age_remelting().unwrap();
    let events = vec![
        TimedEvent {
            time_s: 20.0,
            event: birth(NEW, 600.0, 0.05),
        },
        TimedEvent {
            time_s: 20.0,
            event: Event::Contact(ThermalLink {
                id: 1,
                a: OLD,
                b: NEW,
                area_m2: 0.001,
                conductance_w_k: 1.25,
            }),
        },
    ];
    let start = prepare(&law, events)
        .with_equilibrium_remelting(MAT, remelting)
        .unwrap();
    let (before, _) = start.advance_to(19.5, controls(0.02), || Ok(())).unwrap();
    assert!(before.states()[&OLD].crystalline_mass_fraction() > 0.3);
    let (born, birth_audit) = before.advance_to(20.0, controls(0.02), || Ok(())).unwrap();
    let old_fraction = born.states()[&OLD].crystalline_mass_fraction();
    assert!(old_fraction > 0.3);
    assert_eq!(born.states()[&NEW].exposure(), 0.0);
    assert_eq!(birth_audit.thermal.incoming_energy_j, 15_000.0);
    assert_eq!(birth_audit.thermal.external_heat_j, 0.0);
    close(born.energy_j(), 15_800.0, 1e-8);

    let (partial, audit) = born.advance_to(25.0, controls(0.02), || Ok(())).unwrap();
    close(partial.temperature_k(OLD).unwrap(), 500.0, 1e-10);
    let fraction = partial.states()[&OLD].crystalline_mass_fraction();
    assert!(fraction > 0.0 && fraction < old_fraction);
    assert!(audit.crystallized_mass_kg < 0.0);
    assert_eq!(audit.thermal.external_heat_j, 0.0);
    close(partial.energy_j(), born.energy_j(), 1e-7);

    let (melted, audit) = born.advance_to(100.0, controls(0.02), || Ok(())).unwrap();
    for id in [OLD, NEW] {
        assert_eq!(melted.states()[&id].crystalline_mass_fraction(), 0.0);
        assert_eq!(melted.states()[&id].exposure(), 0.0);
        close(
            melted.temperature_k(id).unwrap(),
            300.0 + 15_800.0 / 60.0,
            0.002,
        );
    }
    close(melted.energy_j(), born.energy_j(), 2e-7);
    close(audit.crystallized_mass_kg, -0.01 * old_fraction, 3e-12);
    close(melted.mass_kg(), 0.06, 1e-17);
    assert_eq!(audit.thermal.external_heat_j, 0.0);
    assert_eq!(audit.thermal.incoming_energy_j, 0.0);
    assert_eq!(melted.links(), born.links());
    assert_eq!(before.cells().len(), 1);
}

#[test]
fn remelting_binding_and_late_cancellation_preserve_the_original_history() {
    let law = cycle_law(650.0);
    let remelting = law.equivalent_age_remelting().unwrap();
    let start = prepare(&law, vec![]);
    assert!(
        start
            .clone()
            .with_equilibrium_remelting(MaterialId(99), remelting)
            .is_err()
    );
    let selected = start.with_equilibrium_remelting(MAT, remelting).unwrap();
    assert!(
        selected
            .clone()
            .with_equilibrium_remelting(MAT, remelting)
            .is_err()
    );
    let cold_only_law = cycle_law(490.0);
    assert!(
        prepare(&cold_only_law, vec![])
            .with_equilibrium_remelting(MAT, remelting)
            .is_err()
    );
    let (mut aged, _) = selected
        .advance_to(20.0, controls(0.02), || Ok(()))
        .unwrap();
    source(&mut aged, 100.0);
    let cells = aged.cells().clone();
    let states = aged.states().clone();
    let mut polls = 0;
    let accepted = aged
        .advance_to(37.0, controls(0.02), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut calls = 0;
    let result = aged.advance_to(37.0, controls(0.02), || {
        calls += 1;
        if calls == polls {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        result,
        Err(CrystallizationTimelineError::Thermal(
            DepositionError::Cancelled
        ))
    ));
    assert_eq!(calls, polls);
    assert_eq!(aged.cells(), &cells);
    assert_eq!(aged.states(), &states);
    assert_eq!(aged.time_s(), 20.0);
    assert_eq!(aged.remelting_model(MAT), Some(remelting));
    let replay = aged.advance_to(37.0, controls(0.02), || Ok(())).unwrap();
    assert_eq!(replay.0.cells(), accepted.0.cells());
    assert_eq!(replay.0.states(), accepted.0.states());
    assert_eq!(replay.1, accepted.1);
}
