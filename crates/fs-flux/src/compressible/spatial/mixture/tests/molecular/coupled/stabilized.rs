//! Stabilization retains the original three-axis forces, wall efforts and inertia boundary.
use super::*;
use crate::compressible::spatial::mixture::molecular::{
    InertialWall3d, StabilizedMolecularStepControl, WallInertiaControl3d,
};
use fs_qty::Mass;

fn native() -> (MixtureEuler3d, MolecularTransport3d<ConstantTransport>) {
    let (original, _) = wall_case();
    let (flow, previous) = original.into_parts();
    let source = SpatialTransport::new(
        flow.gas().clone(),
        law(flow.gas(), 20.0, 0.03, 0.01, 0.04),
        || false,
    )
    .unwrap();
    let bytes = MolecularTransport3d::required_stabilized_workspace_bytes(&flow, &source).unwrap();
    let mut molecular = MolecularTransport3d::new(&flow, source, bytes, || false).unwrap();
    for side in BoundarySide::ALL {
        molecular = molecular.with_wall(side, previous.wall(side)).unwrap();
    }
    (flow, molecular)
}

fn owner(
    flow: MixtureEuler3d,
    molecular: MolecularTransport3d<ConstantTransport>,
) -> (Coupled, ForcingWorkspace3d, usize) {
    let scratch = flow
        .forcing_workspace(flow.required_forcing_workspace_bytes().unwrap())
        .unwrap();
    let bytes = Coupled::required_stabilized_workspace_bytes(&flow, &molecular).unwrap()
        + scratch.workspace_bytes();
    (
        Coupled::new(flow, molecular, bytes, || false).unwrap(),
        scratch,
        bytes,
    )
}

#[derive(Debug, PartialEq)]
struct Published {
    states: Vec<f64>,
    time: Time,
    split: Option<MolecularStepReport3d>,
    advection: Option<crate::compressible::spatial::mixture::StepReport>,
    molecular: Option<MolecularReport>,
    boundary: Vec<[f64; 7]>,
    molecular_boundary: Vec<[f64; 7]>,
    walls: [MolecularWall3d; 6],
}
fn published(flow: &Coupled) -> Published {
    Published {
        states: flow.flow().states().to_vec(),
        time: flow.flow().time(),
        split: flow.last_step().copied(),
        advection: flow.flow().last_step().copied(),
        molecular: flow.molecular().last_step().copied(),
        boundary: loads(flow),
        molecular_boundary: molecular_loads(flow),
        walls: BoundarySide::ALL.map(|side| flow.molecular().wall(side)),
    }
}

#[test]
fn g1_stabilized_spatial_split_retains_native_forcing_heat_and_species_efforts() {
    let controls = StabilizedMolecularStepControl::default();
    let force = [0.2, -0.1, 0.15];
    let source = BodyForce3d::UniformForceDensity(force);
    let dt = 0.002;
    let (mut manual, mut operator) = native();
    let mut manual_scratch = manual.forcing_workspace(usize::MAX).unwrap();
    let first = operator
        .advance_stabilized(&mut manual, Time::new(0.5 * dt), controls.molecular, || {
            false
        })
        .unwrap();
    let middle = manual
        .step_with_forcing(
            Time::new(dt),
            controls.advection,
            source,
            CAP,
            &mut manual_scratch,
            || false,
        )
        .unwrap();
    assert_eq!(middle.elapsed.value(), dt);
    let second = operator
        .advance_stabilized(&mut manual, Time::new(0.5 * dt), controls.molecular, || {
            false
        })
        .unwrap();

    let (flow, molecular) = native();
    let (mut flow, mut scratch, bytes) = owner(flow, molecular);
    let walls = BoundarySide::ALL.map(|side| flow.molecular().wall(side));
    let report = flow
        .step_stabilized_with_forcing(Time::new(dt), controls, source, CAP, &mut scratch, || false)
        .unwrap();
    let physical = report.transport;
    check(&physical);
    final_admission(&flow, &physical, CAP);
    assert_eq!(physical.coupling_retries, 0);
    assert_eq!(physical.first, first.transport);
    assert_eq!(physical.advection, middle);
    assert_eq!(physical.second, second.transport);
    assert_eq!(flow.flow().states(), manual.states());
    assert_eq!(flow.flow().time(), manual.time());
    assert_eq!(flow.last_step(), Some(&physical));
    assert_eq!(flow.molecular().last_step(), Some(&second.transport));
    assert_eq!(
        report.operator_evaluations,
        first.operator_evaluations + second.operator_evaluations
    );
    assert!(report.maximum_stages_used > 2);
    assert_eq!(
        report.maximum_stages_used,
        first.maximum_stages_used.max(second.maximum_stages_used)
    );
    assert!(physical.first.maximum_parabolic_number > controls.molecular.transport.parabolic_cfl);
    assert!(
        report.maximum_scaled_parabolic_number
            <= controls.molecular.transport.parabolic_cfl * (1.0 + 1e-13)
    );
    assert_eq!(
        BoundarySide::ALL.map(|side| flow.molecular().wall(side)),
        walls
    );
    assert_eq!(flow.workspace_bytes() + scratch.workspace_bytes(), bytes);
    assert!(physical.first.boundary_change.energy.abs() > 1e-8);
    close(physical.advection.boundary_change.energy, 0.0, 2e-12);

    // Independent extensive accounts include authored body work once and the
    // actual species-referenced wall energy, including formation references.
    let volume = flow.flow().grid().cell_volume() * flow.flow().cell_count() as f64;
    let mut outward = [0.0; 7];
    for row in loads(&flow) {
        for (total, value) in outward.iter_mut().zip(row) {
            *total += value;
        }
    }
    for axis in 0..3 {
        close(
            physical.source_change.impulse[axis],
            volume * dt * force[axis],
            3e-12,
        );
        close(
            physical.after.momentum[axis] - physical.before.momentum[axis] + outward[3 + axis]
                - physical.source_change.impulse[axis],
            0.0,
            3e-12,
        );
    }
    let references = flow.flow().gas().energy_references();
    let outward_energy = outward[6] + (0..3).map(|k| outward[k] * references[k]).sum::<f64>();
    close(
        physical.after.absolute_energy(flow.flow().gas()).unwrap()
            - physical.before.absolute_energy(flow.flow().gas()).unwrap()
            + outward_energy
            - physical.source_change.work,
        0.0,
        3e-11,
    );
    for k in 0..3 {
        close(
            physical.after.species_masses[k] - physical.before.species_masses[k] + outward[k],
            0.0,
            3e-12,
        );
    }
}

#[test]
fn g0_stabilized_split_cannot_freeze_an_original_finite_inertia_wall() {
    let flow = field(
        equal_calorics(),
        CartesianGrid::new(1, 1, 1, 1.0, 1.0, 1.0).unwrap(),
        Boundaries3d {
            y: WALLS,
            ..PERIODIC
        },
        |_| ([0.7, 0.3, 0.0], [0.0; 3], 1.0),
    );
    let source = SpatialTransport::new(
        flow.gas().clone(),
        law(flow.gas(), 0.001, 0.75, 0.0, 0.0),
        || false,
    )
    .unwrap();
    let bytes = MolecularTransport3d::required_stabilized_workspace_bytes(&flow, &source).unwrap()
        + WallInertiaControl3d::solver_workspace_bytes();
    let molecular = MolecularTransport3d::new(&flow, source, bytes, || false)
        .unwrap()
        .with_wall(
            BoundarySide::YMax,
            MolecularWall3d::NoSlip {
                velocity: [1.0, 0.0, -0.4],
                thermal: WallThermalCondition::Adiabatic,
            },
        )
        .unwrap()
        .with_wall_inertia(
            &[InertialWall3d {
                side: BoundarySide::YMax,
                mass: Mass::new(0.5),
            }],
            WallInertiaControl3d::default(),
        )
        .unwrap();
    let (mut flow, mut scratch, _) = owner(flow, molecular);
    flow.step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    let before = published(&flow);
    let controls = StabilizedMolecularStepControl::default();
    assert_eq!(
        flow.step_stabilized(Time::new(0.001), controls, || false),
        Err(DissipationError::Invalid(
            "stabilized molecular transport requires prescribed walls without finite inertia or acceleration"
        ))
    );
    assert_eq!(published(&flow), before);
    assert_eq!(
        flow.step_stabilized_with_forcing(
            Time::new(0.001),
            controls,
            BodyForce3d::UniformAcceleration([0.1, 0.0, 0.0]),
            CAP,
            &mut scratch,
            || false,
        ),
        Err(DissipationError::Invalid(
            "stabilized molecular transport requires prescribed walls without finite inertia or acceleration"
        ))
    );
    assert_eq!(published(&flow), before);
    // The original reciprocal carrier solve remains available after refusal.
    flow.step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    assert!(flow.flow().time().value() > before.time.value());
    assert_ne!(
        BoundarySide::ALL.map(|side| flow.molecular().wall(side)),
        before.walls
    );
}
