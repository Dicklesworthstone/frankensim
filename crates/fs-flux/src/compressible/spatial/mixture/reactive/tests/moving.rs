use super::*;
use crate::compressible::spatial::MeshMotion;
use crate::compressible::spatial::mixture::molecular::{
    MolecularWall3d, SpatialTransport, WallThermalCondition,
};
use fs_gas::multicomponent::dissipation::{TransportCoefficients, TransportDomain};
use fs_qty::DynViscosity;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;

fn case(reference: f64) -> (ReactingEuler3d, MolecularTransport3d<ConstantTransport>) {
    let (gas, mechanism) = mechanism(5.0, 8.0, reference);
    let field = field(
        gas,
        CartesianGrid::new(2, 2, 2, 0.3, 0.25, 0.2).unwrap(),
        Boundaries3d {
            x: WALLS,
            y: WALLS,
            z: WALLS,
        },
        |x| {
            let a = 0.7 + 0.05 * x[0] - 0.03 * x[1];
            (
                [a, 1.0 - a],
                [0.09 + 0.02 * x[1], -0.06 + 0.03 * x[2], 0.04 - 0.01 * x[0]],
                1.0 + 0.1 * x[2],
            )
        },
    );
    let law = ConstantTransport::new(
        field.gas().source().clone(),
        "moving-reacting-fixture",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 0.001,
            maximum_pressure_pa: 100.0,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(0.03),
            bulk_viscosity: DynViscosity::new(0.01),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.02),
        },
        vec![0.0, 0.02, 0.02, 0.0],
    )
    .unwrap();
    let source = SpatialTransport::new(field.gas().clone(), law, || false).unwrap();
    let bytes = MolecularTransport3d::required_workspace_bytes(&field, &source).unwrap();
    let operator = MolecularTransport3d::new(&field, source, bytes, || false)
        .unwrap()
        .with_wall(
            BoundarySide::ZMin,
            MolecularWall3d::NoSlip {
                velocity: [0.12, -0.08, 0.0],
                thermal: WallThermalCondition::Isothermal(Temperature::new(1.4)),
            },
        )
        .unwrap()
        .with_wall(
            BoundarySide::ZMax,
            MolecularWall3d::NavierSlip {
                velocity: [-0.06, 0.1, 0.0],
                thermal: WallThermalCondition::Isothermal(Temperature::new(0.8)),
                slip_length_m: 0.05,
                fluid_heat_fraction: 0.4,
            },
        )
        .unwrap();
    (coupled(field, mechanism), operator)
}

fn step(
    flow: &mut ReactingEuler3d,
    operator: &mut MolecularTransport3d<ConstantTransport>,
    molecular: bool,
    duration: f64,
    motion: MeshMotion,
    cancelled: impl FnMut() -> bool,
) -> Result<MovingReactiveStepReport3d, ReactiveError> {
    if molecular {
        let bytes = flow.workspace_bytes() + operator.workspace_bytes();
        flow.step_moving_with_dissipation(
            operator,
            Time::new(duration),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            motion,
            cancelled,
        )
    } else {
        flow.step_moving(
            Time::new(duration),
            ReactiveStepControl::default(),
            motion,
            cancelled,
        )
    }
}

fn loads(flow: &ReactingEuler3d) -> Vec<[f64; 6]> {
    BoundarySide::ALL
        .into_iter()
        .flat_map(|side| {
            flow.boundary_impulses(side).map(|q| {
                [
                    q.species_masses[0],
                    q.species_masses[1],
                    q.momentum[0],
                    q.momentum[1],
                    q.momentum[2],
                    q.energy,
                ]
            })
        })
        .collect()
}
fn molecular_loads(operator: &MolecularTransport3d<ConstantTransport>) -> Vec<[f64; 6]> {
    BoundarySide::ALL
        .into_iter()
        .flat_map(|side| {
            operator.boundary_impulses(side).map(|q| {
                [
                    q.species_masses[0],
                    q.species_masses[1],
                    q.momentum[0],
                    q.momentum[1],
                    q.momentum[2],
                    q.energy,
                ]
            })
        })
        .collect()
}

fn moving_balance(flow: &ReactingEuler3d, report: MovingReactiveStepReport3d, motion: MeshMotion) {
    balance(&report.reactive);
    assert_eq!(
        report.reactive.chemistry.cell_solves,
        2 * flow.flow().cell_count()
    );
    let molecular_energy = report.reactive.dissipation.map_or(0.0, |d| {
        d.first.boundary_change.energy + d.second.boundary_change.energy
    });
    close(
        molecular_energy,
        report.molecular_wall_work_j.iter().sum::<f64>() + report.wall_heat_j.iter().sum::<f64>(),
        3e-12,
    );
    let mut rate: f64 = 0.0;
    let grid = flow.flow().grid();
    for i in 0..flow.flow().cell_count() {
        let p = flow
            .flow()
            .cell(i)
            .unwrap()
            .primitive(flow.flow().gas(), || false)
            .unwrap();
        let c = grid.coordinates(i);
        rate = rate.max(
            (0..3)
                .map(|axis| {
                    let lo = motion.velocity(grid, axis, c[axis]);
                    let hi = motion.velocity(grid, axis, c[axis] + 1);
                    ((p.velocity[axis] - lo)
                        .abs()
                        .max((p.velocity[axis] - hi).abs())
                        + p.sound_speed)
                        / grid.spacing()[axis]
                })
                .sum(),
        );
        assert!(
            (flow.chemistry_controls().minimum_temperature_k
                ..=flow.chemistry_controls().maximum_temperature_k)
                .contains(&p.temperature)
        );
    }
    assert!(
        rate * report.reactive.elapsed.value()
            <= ReactiveStepControl::default().transport.cfl * (1.0 + 1e-13)
    );
}

#[test]
fn g0_zero_motion_matches_original_reacting_splits_and_real_thermal_wall_transfer() {
    for molecular in [false, true] {
        let (mut original, mut aop) = case(0.0);
        let (mut moving, mut bop) = case(0.0);
        let dt = Time::new(0.001);
        let expected = if molecular {
            let bytes = original.workspace_bytes() + aop.workspace_bytes();
            original
                .step_with_dissipation(
                    &mut aop,
                    dt,
                    ReactiveStepControl::default(),
                    DissipationControl::default(),
                    bytes,
                    || false,
                )
                .unwrap()
        } else {
            original
                .step(dt, ReactiveStepControl::default(), || false)
                .unwrap()
        };
        let report = step(
            &mut moving,
            &mut bop,
            molecular,
            dt.value(),
            MeshMotion::default(),
            || false,
        )
        .unwrap();
        assert_eq!(report.reactive, expected);
        assert_eq!(moving.flow().states(), original.flow().states());
        assert_eq!(loads(&moving), loads(&original));
        assert_eq!(molecular_loads(&bop), molecular_loads(&aop));
        moving_balance(&moving, report, MeshMotion::default());
        if molecular {
            assert!(report.wall_heat_j[4] > 0.0 && report.wall_heat_j[5] < 0.0);
            assert!(report.molecular_wall_work_j[4].abs() > 1e-8);
        }
    }
}

#[test]
fn g1_translating_exothermic_gas_matches_native_reactor_and_final_source_cfl() {
    let mut flow = uniform();
    let before = flow.flow().totals();
    let bytes = flow.workspace_bytes();
    let volume = flow.flow().grid().cell_volume();
    let speed = [0.7, -0.4, 0.2];
    let motion = MeshMotion {
        lower_velocity: speed,
        upper_velocity: speed,
    };
    let source = flow.mechanism().clone();
    let report = flow
        .step_moving(
            Time::new(0.2),
            ReactiveStepControl::default(),
            motion,
            || false,
        )
        .unwrap();
    moving_balance(&flow, report, motion);
    assert!(report.reactive.coupling_retries > 0);
    assert_eq!(flow.workspace_bytes(), bytes);
    let reference = advance_adiabatic(
        &source,
        &[(SpeciesId::new("A").unwrap(), Amount::new(volume / R))],
        Volume::new(volume),
        Temperature::new(1.0),
        report.reactive.elapsed,
        chemistry(),
        || false,
    )
    .unwrap();
    for i in 0..flow.flow().cell_count() {
        let q = flow.flow().cell(i).unwrap();
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        close(p.temperature, reference.temperature().value(), 3e-5);
        close(
            q.partial_densities[1],
            reference.amounts()[1].value() * R / volume,
            3e-5,
        );
        for axis in 0..3 {
            close(p.velocity[axis], speed[axis], 3e-12);
        }
    }
    for axis in 0..3 {
        close(
            report.reactive.after.momentum[axis],
            before.momentum[axis],
            3e-12,
        );
        close(
            flow.flow().grid().origin()[axis],
            speed[axis] * report.reactive.elapsed.value(),
            3e-15,
        );
    }
    assert_eq!(report.molecular_wall_work_j, [0.0; 6]);
    assert_eq!(report.wall_heat_j, [0.0; 6]);
    // The accepted duration replays the entire source/ALE/source interval.
    let mut replay = uniform();
    let repeated = replay
        .step_moving(
            report.reactive.elapsed,
            ReactiveStepControl::default(),
            motion,
            || false,
        )
        .unwrap();
    assert_eq!(repeated.reactive.coupling_retries, 0);
    assert_eq!(repeated.reactive.chemistry, report.reactive.chemistry);
    assert_eq!(replay.flow().states(), flow.flow().states());
}

#[test]
fn g2_reacting_viscous_deformation_closes_all_momenta_work_heat_and_rebinds_source() {
    for axis in 0..3 {
        let (mut flow, mut op) = case(0.0);
        let bytes = flow.workspace_bytes() + op.workspace_bytes();
        let grid = *flow.flow().grid();
        let mut motion = MeshMotion::default();
        motion.lower_velocity[axis] = -0.04;
        motion.upper_velocity[axis] = 0.08;
        let report = step(&mut flow, &mut op, true, 0.001, motion, || false).unwrap();
        moving_balance(&flow, report, motion);
        assert_eq!(flow.workspace_bytes() + op.workspace_bytes(), bytes);
        let dt = report.reactive.elapsed.value();
        close(
            flow.flow().grid().spacing()[axis],
            grid.spacing()[axis] + 0.12 * dt / grid.shape()[axis] as f64,
            2e-14,
        );
        assert!(report.reactive.chemistry.species_change_kg[0] < 0.0);
        assert!(report.reactive.chemistry.referenced_energy_change_j > 0.0);
        for q in loads(&flow) {
            assert_eq!(&q[..2], &[0.0; 2]);
        }
        let d = report.reactive.dissipation.unwrap();
        assert_eq!(d.first.elapsed.value(), 0.5 * dt);
        assert_eq!(d.second.elapsed, d.first.elapsed);
        assert_eq!(op.last_step(), Some(&d.second));
        assert!(report.molecular_wall_work_j.iter().any(|w| w.abs() > 1e-8));
        assert!(report.wall_heat_j[4] > 0.0 && report.wall_heat_j[5] < 0.0);
        balance(
            &flow
                .step_with_dissipation(
                    &mut op,
                    Time::new(0.0005),
                    ReactiveStepControl::default(),
                    DissipationControl::default(),
                    bytes,
                    || false,
                )
                .unwrap(),
        );
    }
}

#[test]
fn g2_large_common_chemical_gauge_preserves_reacting_moving_solution() {
    let (mut a, mut ao) = case(0.0);
    let (mut b, mut bo) = case(1e6);
    let motion = MeshMotion {
        lower_velocity: [-0.04, 0.0, 0.0],
        upper_velocity: [0.08, 0.0, 0.0],
    };
    let ra = step(&mut a, &mut ao, true, 0.001, motion, || false).unwrap();
    let rb = step(&mut b, &mut bo, true, 0.001, motion, || false).unwrap();
    close(
        ra.reactive.elapsed.value(),
        rb.reactive.elapsed.value(),
        2e-15,
    );
    close(rb.reactive.absolute_energy_residual_j, 0.0, 2e-8);
    close(rb.reactive.mass_residual_kg, 0.0, 2e-12);
    for i in 0..a.flow().cell_count() {
        let qa = a.flow().cell(i).unwrap();
        let qb = b.flow().cell(i).unwrap();
        for k in 0..2 {
            close(qa.partial_densities[k], qb.partial_densities[k], 2e-8);
        }
        let pa = qa.primitive(a.flow().gas(), || false).unwrap();
        let pb = qb.primitive(b.flow().gas(), || false).unwrap();
        close(pa.temperature, pb.temperature, 2e-8);
        close(pa.pressure, pb.pressure, 2e-8);
        for axis in 0..3 {
            close(pa.velocity[axis], pb.velocity[axis], 2e-9);
        }
    }
    for side in 0..6 {
        close(ra.wall_heat_j[side], rb.wall_heat_j[side], 3e-10);
        close(
            ra.molecular_wall_work_j[side],
            rb.molecular_wall_work_j[side],
            3e-10,
        );
    }
}

#[test]
fn g4_moving_reactive_cancellation_preserves_every_owner_and_replays_whole_interval() {
    let motion = MeshMotion {
        lower_velocity: [-0.04, 0.0, 0.0],
        upper_velocity: [0.08, 0.0, 0.0],
    };
    for molecular in [false, true] {
        let make = || {
            let (mut flow, mut op) = case(0.0);
            step(&mut flow, &mut op, molecular, 0.0005, motion, || false).unwrap();
            (flow, op)
        };
        let (mut expected, mut expected_op) = make();
        let mut calls = 0;
        let accepted = step(
            &mut expected,
            &mut expected_op,
            molecular,
            0.001,
            motion,
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
        for cutoff in [1, calls / 3, 2 * calls / 3, calls] {
            let (mut flow, mut op) = make();
            let states = flow.flow().states().to_vec();
            let grid = *flow.flow().grid();
            let clock = flow.flow().time();
            let report = *flow.last_step().unwrap();
            let op_report = op.last_step().copied();
            let boundary = loads(&flow);
            let op_boundary = molecular_loads(&op);
            let bytes = flow.workspace_bytes() + op.workspace_bytes();
            let mut seen = 0;
            assert_eq!(
                step(&mut flow, &mut op, molecular, 0.001, motion, || {
                    seen += 1;
                    seen == cutoff
                }),
                Err(ReactiveError::Cancelled)
            );
            assert_eq!(flow.flow().states(), states);
            assert_eq!(*flow.flow().grid(), grid);
            assert_eq!(flow.flow().time(), clock);
            assert_eq!(flow.last_step(), Some(&report));
            assert_eq!(op.last_step(), op_report.as_ref());
            assert_eq!(loads(&flow), boundary);
            assert_eq!(molecular_loads(&op), op_boundary);
            assert_eq!(flow.workspace_bytes() + op.workspace_bytes(), bytes);
            let recovered = step(&mut flow, &mut op, molecular, 0.001, motion, || false).unwrap();
            assert_eq!(recovered, accepted);
            assert_eq!(flow.flow().states(), expected.flow().states());
            assert_eq!(loads(&flow), loads(&expected));
            assert_eq!(molecular_loads(&op), molecular_loads(&expected_op));
        }
    }
}
