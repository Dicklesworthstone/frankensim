//! Original simultaneous impact and explicit finite-cell caloric destinations.
use super::*;
use fs_couple::friction_heat::directional::impact::heat::*;
use std::cell::Cell;

fn controls() -> ImpactHeatControls {
    ImpactHeatControls {
        mechanical: impact_controls(),
        cell_energy_tolerance_j: 1e-9,
        energy_tolerance_j: 1e-9,
        combined_energy_tolerance_j: 3e-8,
    }
}
fn collision() -> (
    [RigidBodyState; 1],
    [MassProperties; 1],
    [RigidFrictionalContact<()>; 1],
) {
    (
        [body(Vec3::ZERO, Vec3::new(1.0, 0.0, -1.0))],
        [MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 2.0, 2.0)).unwrap()],
        [RigidFrictionalContact {
            id: 71,
            point: RigidImpactPoint {
                body: 0,
                arm_body: Vec3::ZERO,
            },
            other: RigidImpactCounterpart::PrescribedPlane {
                point_world: Vec3::ZERO,
                velocity_world: Vec3::ZERO,
            },
            normal_world: Vec3::new(0.0, 0.0, 1.0),
            tangent_world: Vec3::new(1.0, 0.0, 0.0),
            activation_tolerance_m: 0.0,
            friction_coefficient: (),
            normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
        }],
    )
}

#[test]
fn explicit_global_loss_partition_heats_original_cells_without_aging_or_redistributing_contacts() {
    let (_, law, query) = source_fixture(scalar(0.2));
    let curves = [material(0).1, material(1).1];
    let mut receiver = thermal(&curves, [300.0; 2]);
    let before = receiver.clone();
    let (mut bodies, properties, contacts) = collision();
    let mut histories = [law.initial_state()];
    let report = apply_source_rigid_elliptic_impacts_with_heat(
        &mut bodies,
        &properties,
        &contacts,
        &mut histories,
        &[binding(&law)],
        &mut receiver,
        &BTreeMap::from([(CELLS[0], 0.25), (CELLS[1], 0.5)]),
        controls(),
        |_, _| Ok::<_, &'static str>(conditions(&query)),
        &CancelGate::new(),
    )
    .unwrap();
    // Independent central impact: pn=2, pt=-.4, vx_after=.8, K_after=.64 J.
    close(report.impact.mechanical.contact.dissipated_energy, 1.36);
    close(receiver.temperature_k(CELLS[0]).unwrap(), 317.0);
    close(
        receiver.temperature_k(CELLS[1]).unwrap(),
        300.0 + 0.68 / 0.06,
    );
    close(report.heat.thermal.stored_energy_change_j, 1.02);
    close(report.heat.unallocated_energy_j, 0.34);
    close(report.combined_energy_residual_j, 0.0);
    assert_eq!(receiver.time_s(), before.time_s());
    assert_eq!(receiver.links(), before.links());
    assert_eq!(
        (report.heat.thermal.events, report.heat.thermal.substeps),
        (0, 0)
    );
    for id in CELLS {
        assert_eq!(receiver.cells()[&id].boundary, before.cells()[&id].boundary);
        assert_eq!(receiver.cells()[&id].mass_kg, before.cells()[&id].mass_kg);
    }
}

#[test]
fn empty_or_zero_allocations_reproduce_original_impact_and_histories_exactly() {
    let (_, law, query) = source_fixture(scalar(0.2));
    let curves = [material(0).1, material(1).1];
    let start = thermal(&curves, [300.0; 2]);
    let (initial, properties, contacts) = collision();
    let mut original_bodies = initial;
    let mut original_histories = [law.initial_state()];
    let original = apply_source_rigid_elliptic_impacts(
        &mut original_bodies,
        &properties,
        &contacts,
        &mut original_histories,
        &[binding(&law)],
        impact_controls(),
        |_, _| Ok::<_, &'static str>(conditions(&query)),
        &CancelGate::new(),
    )
    .unwrap();
    for fractions in [BTreeMap::new(), BTreeMap::from([(CELLS[0], 0.0)])] {
        let mut receiver = start.clone();
        let mut bodies = initial;
        let mut histories = [law.initial_state()];
        let report = apply_source_rigid_elliptic_impacts_with_heat(
            &mut bodies,
            &properties,
            &contacts,
            &mut histories,
            &[binding(&law)],
            &mut receiver,
            &fractions,
            controls(),
            |_, _| Ok::<_, &'static str>(conditions(&query)),
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(bodies, original_bodies);
        assert_eq!(histories, original_histories);
        assert_eq!(
            report.impact.mechanical.contact.impulses,
            original.mechanical.contact.impulses
        );
        assert_eq!(
            report.heat.unallocated_energy_j,
            original.mechanical.contact.dissipated_energy
        );
        assert_eq!(receiver.cells(), start.cells());
        assert_eq!(report.heat.thermal.external_heat_j, 0.0);
    }
}

struct CancellingLaw<'a> {
    original: &'a dyn ThermalLaw,
    gate: &'a CancelGate,
    armed: &'a Cell<bool>,
}
impl ThermalLaw for CancellingLaw<'_> {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let temperature = self.original.temperature_k(h)?;
        if self.armed.get() {
            self.gate.request();
        }
        Ok(temperature)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        self.original.maximum_temperature_slope()
    }
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.original.storage_basis()
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        self.original.caloric_identity()
    }
    fn validate_cell(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
        self.original.validate_cell(cell)
    }
}

#[test]
fn late_thermal_refusal_or_cancellation_preserves_bodies_histories_and_caloric_receiver() {
    let (_, law, query) = source_fixture(scalar(0.2));
    let curves = [material(0).1, material(1).1];
    let (initial, properties, contacts) = collision();
    for cancel in [false, true] {
        let gate = CancelGate::new();
        let armed = Cell::new(false);
        let guarded = CancellingLaw {
            original: curves[0].curve(),
            gate: &gate,
            armed: &armed,
        };
        let start = thermal(&curves, [if cancel { 300.0 } else { 399.0 }; 2]);
        let cell = start.cells()[&CELLS[0]].clone();
        let mut receiver = DepositionSimulation::new(
            vec![TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }],
            &[(MaterialId(0), &guarded as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: 1,
                max_links: 0,
                max_events: 1,
                max_substeps: 1,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-9,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
        let before = receiver.clone();
        let mut bodies = initial;
        let mut histories = [law.initial_state()];
        let old_histories = histories.clone();
        let result = apply_source_rigid_elliptic_impacts_with_heat(
            &mut bodies,
            &properties,
            &contacts,
            &mut histories,
            &[binding(&law)],
            &mut receiver,
            &BTreeMap::from([(CELLS[0], 0.5)]),
            controls(),
            |_, _| {
                armed.set(cancel);
                Ok::<_, &'static str>(conditions(&query))
            },
            &gate,
        );
        assert!(matches!(
            result,
            Err(ImpactHeatError::Thermal(
                DepositionError::MaterialDomain | DepositionError::Cancelled
            ))
        ));
        assert_eq!(bodies, initial);
        assert_eq!(histories, old_histories);
        assert_eq!(receiver.cells(), before.cells());
        assert_eq!(receiver.time_s(), before.time_s());
    }
}
