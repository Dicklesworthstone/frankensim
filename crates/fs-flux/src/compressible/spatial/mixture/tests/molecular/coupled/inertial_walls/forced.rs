use super::*;

fn work(flow: &Coupled) -> ForcingWorkspace3d {
    flow.flow()
        .forcing_workspace(flow.flow().required_forcing_workspace_bytes().unwrap())
        .unwrap()
}

fn check_forced(r: MolecularInertialWallsStep, gas: &MixtureGas) {
    let m = r.mechanics;
    let report = m.gas.transport;
    check(&report);
    assert!(m.combined_momentum_residual_n_s.abs() <= controls().momentum_tolerance_n_s);
    assert!(m.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
    for end in 0..2 {
        assert!(m.momentum_residual_n_s[end].abs() <= controls().momentum_tolerance_n_s);
        assert!(m.mechanical_energy_residual_j[end].abs() <= controls().energy_tolerance_j);
        for a in 0..3 {
            assert_eq!(
                r.guide_angular_impulse_n_m_s[end][a],
                -r.gas_angular_impulse_n_m_s[end][a]
            );
        }
    }
    let wall_change: f64 = (0..2)
        .map(|i| {
            m.kinetic_energy_change_j[i] + m.exterior_pressure_work_j[i]
                - m.external_force_work_j[i]
        })
        .sum();
    let heat = m.gas.wall_heat_j.iter().sum::<f64>();
    let referenced =
        report.after.energy - report.before.energy + wall_change - heat - report.source_change.work;
    assert!(referenced.abs() <= controls().energy_tolerance_j);
    let chemical: f64 = (0..gas.species_count())
        .map(|k| {
            (report.after.species_masses[k] - report.before.species_masses[k])
                * gas.energy_references()[k]
        })
        .sum();
    assert!((referenced + chemical).abs() <= controls().energy_tolerance_j);
    close(
        report.first.boundary_change.energy + report.second.boundary_change.energy,
        heat + m.gas.molecular_wall_work_j.iter().sum::<f64>(),
        4e-11,
    );
    for a in 0..3 {
        let piston = if a == r.axis {
            // The mechanical receipt reports impulse ON the pistons.
            -m.gas_impulse_n_s.iter().sum::<f64>()
        } else {
            r.guide_impulse_n_s[0][a] + r.guide_impulse_n_s[1][a]
        };
        close(
            report.after.momentum[a]
                - report.before.momentum[a]
                - report.source_change.impulse[a]
                - r.side_boundary_impulse_n_s[a]
                - piston,
            0.0,
            3e-11,
        );
    }
}

#[test]
fn g0_zero_body_force_retains_native_molecular_piston_reports_and_compensated_balances() {
    let mut a = case(0, 0.0);
    let mut b = case(0, 0.0);
    let mut wa = pistons();
    let mut wb = pistons();
    let mut scratch = work(&b);
    let expected = a
        .step_with_inertial_walls(
            0,
            &mut wa,
            Time::new(0.001),
            MolecularStepControl::default(),
            external_loads(),
            controls(),
            || false,
        )
        .unwrap();
    let actual = b
        .step_with_inertial_walls_and_forcing(
            0,
            &mut wb,
            Time::new(0.001),
            MolecularStepControl::default(),
            external_loads(),
            controls(),
            BodyForce3d::UniformAcceleration([0.0; 3]),
            CAP,
            &mut scratch,
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(wa, wb);
    assert_eq!(a.flow().states(), b.flow().states());
    assert_eq!(loads(&a), loads(&b));
    assert_eq!(molecular_loads(&a), molecular_loads(&b));
}

#[test]
fn g1_body_force_and_finite_molecular_pistons_close_work_heat_support_and_all_momenta() {
    let force_density = [0.2, -0.15, 0.1];
    for axis in 0..3 {
        let mut flow = case(axis, 0.0);
        let mut scratch = work(&flow);
        let bytes = flow.workspace_bytes() + scratch.workspace_bytes();
        let old_volume = flow.flow().grid().cell_volume() * flow.flow().cell_count() as f64;
        let mut walls = pistons();
        let r = flow
            .step_with_inertial_walls_and_forcing(
                axis,
                &mut walls,
                Time::new(0.001),
                MolecularStepControl::default(),
                external_loads(),
                controls(),
                BodyForce3d::UniformForceDensity(force_density),
                CAP,
                &mut scratch,
                || false,
            )
            .unwrap();
        check_forced(r, flow.flow().gas());
        let report = r.mechanics.gas.transport;
        let new_volume = flow.flow().grid().cell_volume() * flow.flow().cell_count() as f64;
        for a in 0..3 {
            close(
                report.source_change.impulse[a],
                force_density[a] * 0.5 * (old_volume + new_volume) * report.elapsed.value(),
                3e-12,
            );
        }
        assert!(report.source_change.work.abs() > 1e-8);
        assert!(r.side_boundary_impulse_n_s[axis].abs() > 1e-7);
        assert!(r.guide_impulse_n_s.iter().flatten().any(|x| x.abs() > 1e-7));
        assert!(r.mechanics.gas.wall_heat_j[2 * axis + 1] > 0.0);
        assert_eq!(flow.workspace_bytes() + scratch.workspace_bytes(), bytes);
        assert_eq!(flow.last_step(), Some(&report));
        assert_eq!(flow.molecular().last_step(), Some(&report.second));
        for side in BoundarySide::ALL {
            for q in flow.boundary_impulses(side) {
                assert_eq!(q.species_masses, &[0.0; 3]);
            }
        }
    }
}

#[test]
fn g1_species_energy_references_do_not_change_forced_finite_piston_mechanics() {
    let force = BodyForce3d::UniformAcceleration([0.17, -0.13, 0.09]);
    let mut a = case(0, 0.0);
    let mut b = case(0, 1024.0);
    let mut aw = work(&a);
    let mut bw = work(&b);
    let mut wa = pistons();
    let mut wb = pistons();
    for _ in 0..2 {
        let ra = a
            .step_with_inertial_walls_and_forcing(
                0,
                &mut wa,
                Time::new(0.001),
                MolecularStepControl::default(),
                external_loads(),
                controls(),
                force,
                CAP,
                &mut aw,
                || false,
            )
            .unwrap();
        let rb = b
            .step_with_inertial_walls_and_forcing(
                0,
                &mut wb,
                Time::new(0.001),
                MolecularStepControl::default(),
                external_loads(),
                controls(),
                force,
                CAP,
                &mut bw,
                || false,
            )
            .unwrap();
        check_forced(ra, a.flow().gas());
        check_forced(rb, b.flow().gas());
        assert_eq!(
            ra.mechanics.gas.transport.elapsed,
            rb.mechanics.gas.transport.elapsed
        );
        for end in 0..2 {
            close(wa[end].velocity_m_s(), wb[end].velocity_m_s(), 4e-10);
            close(
                ra.mechanics.gas_impulse_n_s[end],
                rb.mechanics.gas_impulse_n_s[end],
                4e-11,
            );
        }
    }
    for i in 0..a.flow().cell_count() {
        let qa = a.flow().cell(i).unwrap();
        let qb = b.flow().cell(i).unwrap();
        for k in 0..3 {
            close(qa.partial_densities[k], qb.partial_densities[k], 4e-10);
        }
        let pa = qa.primitive(a.flow().gas(), || false).unwrap();
        let pb = qb.primitive(b.flow().gas(), || false).unwrap();
        close(pa.temperature, pb.temperature, 4e-10);
        close(pa.pressure, pb.pressure, 4e-10);
        for d in 0..3 {
            close(pa.velocity[d], pb.velocity[d], 4e-10);
        }
    }
}

#[test]
fn g4_forced_piston_final_cancellation_and_exhaustion_publish_no_partial_owner() {
    let force = BodyForce3d::UniformAcceleration([0.17, -0.13, 0.09]);
    let make = || {
        let mut flow = case(0, 0.0);
        let mut scratch = work(&flow);
        let mut walls = pistons();
        flow.step_with_inertial_walls_and_forcing(
            0,
            &mut walls,
            Time::new(0.0005),
            MolecularStepControl::default(),
            external_loads(),
            controls(),
            force,
            CAP,
            &mut scratch,
            || false,
        )
        .unwrap();
        (flow, scratch, walls)
    };
    let (mut expected, mut expected_work, mut expected_walls) = make();
    let mut polls = 0;
    let accepted = expected
        .step_with_inertial_walls_and_forcing(
            0,
            &mut expected_walls,
            Time::new(0.001),
            MolecularStepControl::default(),
            external_loads(),
            controls(),
            force,
            CAP,
            &mut expected_work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [1, polls / 2, polls] {
        let (mut flow, mut scratch, mut walls) = make();
        let states = flow.flow().states().to_vec();
        let grid = *flow.flow().grid();
        let time = flow.flow().time();
        let old_walls = walls;
        let report = *flow.last_step().unwrap();
        let molecular = *flow.molecular().last_step().unwrap();
        let boundary = loads(&flow);
        let molecular_boundary = molecular_loads(&flow);
        let bytes = flow.workspace_bytes() + scratch.workspace_bytes();
        let mut seen = 0;
        assert_eq!(
            flow.step_with_inertial_walls_and_forcing(
                0,
                &mut walls,
                Time::new(0.001),
                MolecularStepControl::default(),
                external_loads(),
                controls(),
                force,
                CAP,
                &mut scratch,
                || {
                    seen += 1;
                    seen == cutoff
                },
            ),
            Err(MolecularInertialWallsError::Cancelled)
        );
        assert_eq!(flow.flow().states(), states);
        assert_eq!(*flow.flow().grid(), grid);
        assert_eq!(flow.flow().time(), time);
        assert_eq!(walls, old_walls);
        assert_eq!(flow.last_step(), Some(&report));
        assert_eq!(flow.molecular().last_step(), Some(&molecular));
        assert_eq!(loads(&flow), boundary);
        assert_eq!(molecular_loads(&flow), molecular_boundary);
        assert_eq!(flow.workspace_bytes() + scratch.workspace_bytes(), bytes);
        let mut exhausted = controls();
        exhausted.maximum_iterations = 1;
        exhausted.momentum_tolerance_n_s = 1e-30;
        assert!(
            flow.step_with_inertial_walls_and_forcing(
                0,
                &mut walls,
                Time::new(0.001),
                MolecularStepControl {
                    maximum_coupling_retries: 0,
                    ..MolecularStepControl::default()
                },
                external_loads(),
                exhausted,
                force,
                CAP,
                &mut scratch,
                || false,
            )
            .is_err()
        );
        assert_eq!(flow.flow().states(), states);
        assert_eq!(*flow.flow().grid(), grid);
        assert_eq!(walls, old_walls);
        assert_eq!(flow.last_step(), Some(&report));
        let replay = flow
            .step_with_inertial_walls_and_forcing(
                0,
                &mut walls,
                Time::new(0.001),
                MolecularStepControl::default(),
                external_loads(),
                controls(),
                force,
                CAP,
                &mut scratch,
                || false,
            )
            .unwrap();
        assert_eq!(replay, accepted);
        assert_eq!(walls, expected_walls);
        assert_eq!(flow.flow().states(), expected.flow().states());
        assert_eq!(loads(&flow), loads(&expected));
        assert_eq!(molecular_loads(&flow), molecular_loads(&expected));
    }
}
