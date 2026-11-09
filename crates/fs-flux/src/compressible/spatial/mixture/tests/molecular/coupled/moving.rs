use super::*;
use crate::compressible::spatial::MeshMotion;
use crate::compressible::spatial::mixture::molecular::MovingMolecularStepReport3d;

fn check_moving(report: &MovingMolecularStepReport3d) {
    check(&report.transport);
    let molecular_energy = report.transport.first.boundary_change.energy
        + report.transport.second.boundary_change.energy;
    close(
        molecular_energy,
        report.molecular_wall_work_j.iter().sum::<f64>() + report.wall_heat_j.iter().sum::<f64>(),
        4e-11,
    );
}

#[test]
fn g0_stationary_moving_split_retains_native_bits_and_real_belt_work_and_heat() {
    let (mut original, _) = wall_case();
    let (mut moving, _) = wall_case();
    let dt = Time::new(0.002);
    let expected = original
        .step(dt, MolecularStepControl::default(), || false)
        .unwrap();
    let actual = moving
        .step_moving(
            dt,
            MolecularStepControl::default(),
            MeshMotion::default(),
            || false,
        )
        .unwrap();
    check_moving(&actual);
    assert_eq!(actual.transport, expected);
    assert_eq!(moving.flow().states(), original.flow().states());
    assert_eq!(loads(&moving), loads(&original));
    assert_eq!(molecular_loads(&moving), molecular_loads(&original));
    assert!(actual.wall_heat_j[4] > 0.0 && actual.wall_heat_j[5] < 0.0);
    assert!(actual.molecular_wall_work_j[4].abs() > 1e-7);
    assert!(actual.molecular_wall_work_j[5].abs() > 1e-7);
}

fn translating_species(n: usize) -> f64 {
    let h = 1.0 / n as f64;
    let grid = CartesianGrid::new(n, n, n, h, h, h).unwrap();
    let speed = [0.3, -0.25, 0.18];
    let motion = MeshMotion {
        lower_velocity: speed,
        upper_velocity: speed,
    };
    let wave = |x: [f64; 3]| (2.0 * std::f64::consts::PI * x.iter().sum::<f64>()).sin();
    let z = std::f64::consts::PI * h;
    let amplitude = 0.15 * (z.sin() / z).powi(3);
    let flow = field(equal_mass_gas(), grid, PERIODIC, |x| {
        let first = 0.5 + amplitude * wave(x);
        ([first, 1.0 - first, 0.0], speed, 1.0)
    });
    let molecular = prepare(&flow, law(flow.gas(), 0.02, 0.03, 0.01, 0.04));
    let (mut flow, _) = coupled(flow, molecular);
    let bytes = flow.workspace_bytes();
    let end = 0.01;
    while flow.flow().time().value() < end {
        let report = flow
            .step_moving(
                Time::new((end - flow.flow().time().value()).min(0.001)),
                MolecularStepControl::default(),
                motion,
                || false,
            )
            .unwrap();
        check_moving(&report);
        assert_eq!(report.molecular_wall_work_j, [0.0; 6]);
        assert_eq!(report.wall_heat_j, [0.0; 6]);
    }
    assert_eq!(flow.workspace_bytes(), bytes);
    assert_eq!(flow.flow().grid().spacing(), grid.spacing());
    for axis in 0..3 {
        close(flow.flow().grid().origin()[axis], end * speed[axis], 2e-14);
    }
    let decay = (-0.02 * 12.0 * std::f64::consts::PI.powi(2) * end).exp();
    let mut error = 0.0;
    for index in 0..grid.len() {
        let q = flow.flow().cell(index).unwrap();
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        let expected = 0.5 + amplitude * decay * wave(grid.cell_center(index).unwrap());
        error += (q.partial_densities[0] - expected).abs();
        assert_eq!(q.partial_densities[2], 0.0);
        close(p.pressure, 1.0, 4e-10);
        close(p.temperature, 1.0, 4e-10);
        for axis in 0..3 {
            close(p.velocity[axis], speed[axis], 3e-11);
        }
    }
    error / grid.len() as f64
}

#[test]
fn g1_moving_three_axis_species_wave_retains_calorics_and_refines_diffusion() {
    let coarse = translating_species(4);
    let fine = translating_species(8);
    assert!(
        fine < 0.4 * coarse,
        "moving Maxwell-Stefan L1 {coarse} / {fine}"
    );
    assert!(fine < 2e-4, "moving Maxwell-Stefan L1 {fine}");
}

fn translated_column(boost: f64, kind: usize) -> (Coupled, MeshMotion) {
    let flow = field(
        equal_mass_gas(),
        CartesianGrid::new(2, 2, 1, 0.4, 0.3, 0.5).unwrap(),
        Boundaries3d {
            x: WALLS,
            ..PERIODIC
        },
        |x| {
            (
                [0.35 + 0.1 * x[0], 0.65 - 0.1 * x[0], 0.0],
                [0.05 + boost, 0.12, -0.08],
                1.0 + 0.03 * x[0],
            )
        },
    );
    let mut molecular = prepare(&flow, law(flow.gas(), 0.01, 0.08, 0.025, 0.03));
    for side in [BoundarySide::XMin, BoundarySide::XMax] {
        let wall = match kind {
            0 => MolecularWall3d::SlipAdiabatic,
            1 => MolecularWall3d::NoSlip {
                velocity: [0.0; 3],
                thermal: WallThermalCondition::Isothermal(Temperature::new(1.2)),
            },
            _ => MolecularWall3d::NavierSlip {
                velocity: [0.0; 3],
                thermal: WallThermalCondition::Adiabatic,
                slip_length_m: 0.15,
                fluid_heat_fraction: 0.4,
            },
        };
        molecular = molecular.with_wall(side, wall).unwrap();
    }
    let motion = MeshMotion {
        lower_velocity: [0.2 + boost, 0.0, 0.0],
        upper_velocity: [0.2 + boost, 0.0, 0.0],
    };
    (coupled(flow, molecular).0, motion)
}

#[test]
fn g1_normal_viscous_work_is_galilean_for_slip_no_slip_and_navier_walls() {
    let boost = 0.6;
    for kind in 0..3 {
        let (mut base, motion) = translated_column(0.0, kind);
        let (mut shifted, shifted_motion) = translated_column(boost, kind);
        let dt = Time::new(0.001);
        let a = base
            .step_moving(dt, MolecularStepControl::default(), motion, || false)
            .unwrap();
        let b = shifted
            .step_moving(dt, MolecularStepControl::default(), shifted_motion, || {
                false
            })
            .unwrap();
        check_moving(&a);
        check_moving(&b);
        assert!(a.molecular_wall_work_j[..2].iter().all(|w| *w > 0.0));
        for index in 0..base.flow().cell_count() {
            let qa = base.flow().cell(index).unwrap();
            let qb = shifted.flow().cell(index).unwrap();
            let pa = qa.primitive(base.flow().gas(), || false).unwrap();
            let pb = qb.primitive(shifted.flow().gas(), || false).unwrap();
            for species in 0..3 {
                close(
                    qa.partial_densities[species],
                    qb.partial_densities[species],
                    2e-12,
                );
            }
            close(pb.velocity[0], pa.velocity[0] + boost, 2e-12);
            for axis in 1..3 {
                close(pb.velocity[axis], pa.velocity[axis], 2e-12);
            }
            close(pb.temperature, pa.temperature, 2e-11);
            close(pb.pressure, pa.pressure, 2e-11);
        }
        for side in 0..6 {
            close(a.wall_heat_j[side], b.wall_heat_j[side], 2e-12);
        }
        let molecular_impulse = a.transport.first.boundary_change.momentum[0]
            + a.transport.second.boundary_change.momentum[0];
        close(
            b.molecular_wall_work_j.iter().sum::<f64>()
                - a.molecular_wall_work_j.iter().sum::<f64>(),
            boost * molecular_impulse,
            3e-12,
        );
        if kind == 2 {
            assert!(a.wall_heat_j[0] < 0.0 && a.wall_heat_j[1] < 0.0);
        }
    }
}

#[test]
fn g2_deforming_mixture_closes_all_impulses_heat_and_rebinds_stationary_owner() {
    for axis in 0..3 {
        let (mut flow, _) = wall_case();
        let before_grid = *flow.flow().grid();
        let bytes = flow.workspace_bytes();
        let mut motion = MeshMotion::default();
        motion.lower_velocity[axis] = -0.08;
        motion.upper_velocity[axis] = 0.13;
        let report = flow
            .step_moving(
                Time::new(0.002),
                MolecularStepControl::default(),
                motion,
                || false,
            )
            .unwrap();
        check_moving(&report);
        let dt = report.transport.elapsed.value();
        close(
            flow.flow().grid().spacing()[axis],
            before_grid.spacing()[axis] + dt * 0.21 / before_grid.shape()[axis] as f64,
            2e-14,
        );
        close(flow.flow().grid().origin()[axis], -0.08 * dt, 2e-14);
        assert_eq!(flow.workspace_bytes(), bytes);
        for row in loads(&flow) {
            assert_eq!(&row[..3], &[0.0; 3]);
        }
        assert_eq!(flow.molecular().last_step(), Some(&report.transport.second));
        // The old reusable frame and the molecular binding must follow the
        // accepted geometry when the caller next chooses the stationary API.
        let stationary = flow
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        check(&stationary);
        let (mut field, mut source) = flow.into_parts();
        super::super::check(
            source
                .advance(
                    &mut field,
                    Time::new(0.0005),
                    DissipationControl::default(),
                    || false,
                )
                .unwrap(),
        );
    }
}

#[test]
fn g3_moving_cancellation_keeps_geometry_source_and_both_published_loads_atomic() {
    let motion = MeshMotion {
        lower_velocity: [-0.04, 0.0, 0.0],
        upper_velocity: [0.08, 0.0, 0.0],
    };
    let make = || {
        let (mut flow, _) = wall_case();
        flow.step_moving(
            Time::new(0.001),
            MolecularStepControl::default(),
            motion,
            || false,
        )
        .unwrap();
        flow
    };
    let mut expected = make();
    let mut calls = 0;
    let accepted = expected
        .step_moving(
            Time::new(0.002),
            MolecularStepControl::default(),
            motion,
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, calls / 3, 2 * calls / 3, calls] {
        let mut flow = make();
        let states = flow.flow().states().to_vec();
        let grid = *flow.flow().grid();
        let time = flow.flow().time();
        let report = *flow.last_step().unwrap();
        let molecular = *flow.molecular().last_step().unwrap();
        let boundary = loads(&flow);
        let molecular_boundary = molecular_loads(&flow);
        let bytes = flow.workspace_bytes();
        let mut seen = 0;
        let failed = flow.step_moving(
            Time::new(0.002),
            MolecularStepControl::default(),
            motion,
            || {
                seen += 1;
                seen == stop
            },
        );
        assert!(
            matches!(failed, Err(DissipationError::Cancelled)),
            "poll {stop}/{calls}"
        );
        assert_eq!(flow.flow().states(), states);
        assert_eq!(*flow.flow().grid(), grid);
        assert_eq!(flow.flow().time(), time);
        assert_eq!(flow.last_step(), Some(&report));
        assert_eq!(flow.molecular().last_step(), Some(&molecular));
        assert_eq!(loads(&flow), boundary);
        assert_eq!(molecular_loads(&flow), molecular_boundary);
        assert_eq!(flow.workspace_bytes(), bytes);
        let recovered = flow
            .step_moving(
                Time::new(0.002),
                MolecularStepControl::default(),
                motion,
                || false,
            )
            .unwrap();
        assert_eq!(recovered, accepted);
        assert_eq!(flow.flow().states(), expected.flow().states());
        assert_eq!(loads(&flow), loads(&expected));
    }
}
