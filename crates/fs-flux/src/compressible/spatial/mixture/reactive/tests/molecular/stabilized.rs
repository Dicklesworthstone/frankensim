//! The complete reacting spatial split retains original molecular and wall physics.
use super::*;
use crate::compressible::spatial::mixture::molecular::{
    DissipationError, InertialWall3d, MolecularReport, MolecularWall3d,
    StabilizedDissipationControl, WallInertiaControl3d, WallThermalCondition,
};
use fs_qty::Mass;

fn prepared(
    flow: &MixtureEuler3d,
    diffusion: f64,
    extra: usize,
) -> MolecularTransport3d<ConstantTransport> {
    let ordinary = operator(flow, diffusion);
    let source =
        SpatialTransport::new(flow.gas().clone(), ordinary.law().clone(), || false).unwrap();
    let bytes =
        MolecularTransport3d::required_stabilized_workspace_bytes(flow, &source).unwrap() + extra;
    MolecularTransport3d::new(flow, source, bytes, || false).unwrap()
}

fn walls(flow: &MixtureEuler3d) -> MolecularTransport3d<ConstantTransport> {
    prepared(flow, 20.0, 0)
        .with_wall(
            BoundarySide::ZMin,
            MolecularWall3d::NoSlip {
                velocity: [0.2, -0.1, 0.0],
                thermal: WallThermalCondition::Isothermal(Temperature::new(1.2)),
            },
        )
        .unwrap()
        .with_wall(
            BoundarySide::ZMax,
            MolecularWall3d::NoSlip {
                velocity: [0.0, 0.2, 0.0],
                thermal: WallThermalCondition::Isothermal(Temperature::new(0.85)),
            },
        )
        .unwrap()
}

#[derive(Debug, PartialEq)]
struct Published {
    states: Vec<f64>,
    time: Time,
    coupled: Option<ReactiveStepReport3d>,
    advection: Option<crate::compressible::spatial::mixture::StepReport>,
    molecular: Option<MolecularReport>,
    boundary: Vec<[f64; 6]>,
    molecular_boundary: Vec<[f64; 6]>,
    walls: [MolecularWall3d; 6],
}
fn published(
    flow: &ReactingEuler3d,
    molecular: &MolecularTransport3d<ConstantTransport>,
) -> Published {
    Published {
        states: flow.flow().states().to_vec(),
        time: flow.flow().time(),
        coupled: flow.last_step().copied(),
        advection: flow.flow().last_step().copied(),
        molecular: molecular.last_step().copied(),
        boundary: loads(flow),
        molecular_boundary: molecular_loads(molecular),
        walls: BoundarySide::ALL.map(|side| molecular.wall(side)),
    }
}

#[test]
fn g1_stabilized_reacting_spatial_split_retains_native_wall_heat_species_and_body_efforts() {
    let (gas, mechanism) = mechanism(5.0, 8.0, 40.0);
    let make = || {
        field(
            gas.clone(),
            CartesianGrid::new(3, 2, 2, 0.2, 0.25, 0.3).unwrap(),
            Boundaries3d {
                x: WALLS,
                y: WALLS,
                z: WALLS,
            },
            |x| {
                let rho = 1.0 + 0.05 * x[2];
                let a = 0.65 + 0.1 * x[0] - 0.05 * x[1] + 0.02 * x[2];
                (
                    [rho * a, rho * (1.0 - a)],
                    [0.1 + 0.03 * x[1], -0.1 + 0.02 * x[2], 0.05],
                    rho,
                )
            },
        )
    };
    let control = ReactiveStepControl::default();
    let molecular = StabilizedDissipationControl::default();
    let cap = ForceControl3d { maximum_mach: 5.0 };
    let force = [0.2, -0.1, 0.15];
    let source = BodyForce3d::UniformForceDensity(force);
    let dt = 0.002;
    let mut manual = make();
    let mut native_operator = walls(&manual);
    let first = native_operator
        .advance_stabilized(&mut manual, Time::new(0.5 * dt), molecular, || false)
        .unwrap();
    let mut central = coupled(manual, mechanism.clone());
    let mut native_force = central.flow().forcing_workspace(usize::MAX).unwrap();
    let middle = central
        .step_with_forcing(
            Time::new(dt),
            control,
            source,
            cap,
            &mut native_force,
            || false,
        )
        .unwrap();
    assert_eq!(middle.elapsed.value(), dt);
    let mut manual = central.into_flow();
    let second = native_operator
        .advance_stabilized(&mut manual, Time::new(0.5 * dt), molecular, || false)
        .unwrap();

    let mut flow = coupled(make(), mechanism);
    let mut operator = walls(flow.flow());
    let mut work = flow.flow().forcing_workspace(usize::MAX).unwrap();
    let bytes =
        flow.required_stabilized_workspace_bytes(&operator).unwrap() + work.workspace_bytes();
    let original_walls = BoundarySide::ALL.map(|side| operator.wall(side));
    let report = flow
        .step_with_stabilized_dissipation_and_forcing(
            &mut operator,
            Time::new(dt),
            control,
            molecular,
            bytes,
            source,
            cap,
            &mut work,
            || false,
        )
        .unwrap();
    let physical = report.coupled;
    let split = physical.dissipation.unwrap();
    balance(&physical);
    assert_eq!(physical.coupling_retries, 0);
    assert_eq!(split.first, first.transport);
    assert_eq!(split.second, second.transport);
    assert_eq!(physical.transport, middle.transport);
    assert_eq!(physical.chemistry, middle.chemistry);
    assert_eq!(flow.flow().states(), manual.states());
    assert_eq!(flow.flow().time(), manual.time());
    assert_eq!(flow.last_step(), Some(&physical));
    assert_eq!(operator.last_step(), Some(&second.transport));
    assert_eq!(
        molecular_loads(&operator),
        molecular_loads(&native_operator)
    );
    assert_eq!(
        BoundarySide::ALL.map(|side| operator.wall(side)),
        original_walls
    );
    assert_eq!(
        report.operator_evaluations,
        first.operator_evaluations + second.operator_evaluations
    );
    assert_eq!(
        report.maximum_stages_used,
        first.maximum_stages_used.max(second.maximum_stages_used)
    );
    assert!(report.maximum_stages_used > 2);
    assert!(split.first.maximum_parabolic_number > molecular.transport.parabolic_cfl);
    assert!(
        report.maximum_scaled_parabolic_number <= molecular.transport.parabolic_cfl * (1.0 + 1e-13)
    );
    assert_eq!(budget(&flow, &operator) + work.workspace_bytes(), bytes);
    assert_eq!(physical.chemistry.cell_solves, 2 * flow.flow().cell_count());
    assert!(physical.chemistry.species_change_kg[0] < 0.0);
    assert!(physical.chemistry.referenced_energy_change_j > 0.0);
    assert!(split.first.boundary_change.energy.abs() > 1e-8);
    close(physical.transport.boundary_change.energy, 0.0, 2e-12);

    let mut outward = [0.0; 6];
    for row in loads(&flow) {
        for (sum, value) in outward.iter_mut().zip(row) {
            *sum += value;
        }
    }
    let volume = flow.flow().grid().cell_volume() * flow.flow().cell_count() as f64;
    for axis in 0..3 {
        close(
            physical.source_change.impulse[axis],
            volume * dt * force[axis],
            3e-12,
        );
        close(
            physical.after.momentum[axis] - physical.before.momentum[axis] + outward[2 + axis]
                - physical.source_change.impulse[axis],
            0.0,
            3e-12,
        );
    }
    for species in 0..2 {
        close(
            physical.after.species_masses[species] - physical.before.species_masses[species]
                + outward[species]
                - physical.chemistry.species_change_kg[species],
            0.0,
            3e-12,
        );
    }
    let references = flow.flow().gas().energy_references();
    let boundary_energy = outward[5] + (0..2).map(|k| outward[k] * references[k]).sum::<f64>();
    close(
        physical.after.absolute_energy(flow.flow().gas()).unwrap()
            - physical.before.absolute_energy(flow.flow().gas()).unwrap()
            + boundary_energy
            - physical.source_change.work,
        0.0,
        3e-10,
    );
}

#[test]
fn g0_reacting_stabilization_refuses_to_freeze_finite_wall_inertia_without_publishing() {
    let (gas, mechanism) = mechanism(1.0, 2.0, 37.0);
    let raw = field(
        gas,
        CartesianGrid::new(1, 1, 1, 1.0, 1.0, 1.0).unwrap(),
        Boundaries3d {
            y: WALLS,
            ..PERIODIC
        },
        |_| ([0.8, 0.2], [0.0; 3], 1.0),
    );
    let mut operator = prepared(&raw, 0.001, WallInertiaControl3d::solver_workspace_bytes())
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
    let mut flow = coupled(raw, mechanism);
    let mut work = flow.flow().forcing_workspace(usize::MAX).unwrap();
    let bytes =
        flow.required_stabilized_workspace_bytes(&operator).unwrap() + work.workspace_bytes();
    let controls = ReactiveStepControl::default();
    flow.step_with_dissipation(
        &mut operator,
        Time::new(0.001),
        controls,
        DissipationControl::default(),
        bytes,
        || false,
    )
    .unwrap();
    let before = published(&flow, &operator);
    let molecular = StabilizedDissipationControl::default();
    let refusal = ReactiveError::Dissipation(DissipationError::Invalid(
        "stabilized molecular transport requires prescribed walls without finite inertia",
    ));
    assert_eq!(
        flow.step_with_stabilized_dissipation(
            &mut operator,
            Time::new(0.001),
            controls,
            molecular,
            bytes,
            || false,
        ),
        Err(refusal.clone())
    );
    assert_eq!(published(&flow, &operator), before);
    assert_eq!(
        flow.step_with_stabilized_dissipation_and_forcing(
            &mut operator,
            Time::new(0.001),
            controls,
            molecular,
            bytes,
            BodyForce3d::UniformAcceleration([0.1, 0.0, 0.0]),
            ForceControl3d { maximum_mach: 5.0 },
            &mut work,
            || false,
        ),
        Err(refusal)
    );
    assert_eq!(published(&flow, &operator), before);
    flow.step_with_dissipation(
        &mut operator,
        Time::new(0.001),
        controls,
        DissipationControl::default(),
        bytes,
        || false,
    )
    .unwrap();
    assert!(flow.flow().time().value() > before.time.value());
    assert_ne!(
        BoundarySide::ALL.map(|side| operator.wall(side)),
        before.walls
    );
}
