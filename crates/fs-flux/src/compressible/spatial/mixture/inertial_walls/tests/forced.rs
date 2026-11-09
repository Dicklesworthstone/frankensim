use super::*;

const CAP: ForceControl3d = ForceControl3d { maximum_mach: 5.0 };

fn case() -> MixtureEuler3d {
    field(
        standard(0.4),
        CartesianGrid::new(2, 2, 2, 0.3, 0.4, 0.5).unwrap(),
        |x| {
            (
                [0.6 + 0.03 * x[0], 0.4 + 0.02 * x[1], 0.0],
                [0.04 + 0.02 * x[1], -0.03, 0.02],
                1.0 + 0.05 * x[1],
            )
        },
    )
}

fn work(flow: &MixtureEuler3d) -> ForcingWorkspace3d {
    flow.forcing_workspace(flow.required_forcing_workspace_bytes().unwrap())
        .unwrap()
}

fn loads() -> [InertialWallLoad; 2] {
    [
        InertialWallLoad {
            exterior_pressure_pa: 0.8,
            external_force_n: 0.02,
        },
        InertialWallLoad {
            exterior_pressure_pa: 0.7,
            external_force_n: -0.01,
        },
    ]
}

#[test]
fn g0_zero_source_preserves_original_inviscid_piston_transaction() {
    let mut a = case();
    let mut b = case();
    let mut wa = walls(0.0);
    let mut wb = walls(0.0);
    let mut scratch = work(&b);
    let ra = a
        .step_with_inertial_walls(
            0,
            &mut wa,
            Time::new(0.001),
            StepControl::default(),
            loads(),
            controls(0.0),
            || false,
        )
        .unwrap();
    let rb = b
        .step_with_inertial_walls_and_forcing(
            0,
            &mut wb,
            Time::new(0.001),
            StepControl::default(),
            loads(),
            controls(0.0),
            BodyForce3d::UniformAcceleration([0.0; 3]),
            CAP,
            &mut scratch,
            || false,
        )
        .unwrap();
    assert_eq!(ra, rb);
    assert_eq!(wa, wb);
    assert_eq!(a.states(), b.states());
    assert_eq!(a.boundary_impulses, b.boundary_impulses);
}

#[test]
fn g1_held_body_impulse_and_exact_kinetic_work_are_external_to_each_piston() {
    let value = [0.2, -0.15, 0.1];
    for axis in 0..3 {
        for density_force in [false, true] {
            let mut flow = case();
            let mut scratch = work(&flow);
            let bytes = flow.workspace_bytes() + scratch.workspace_bytes();
            let old_grid = *flow.grid();
            let before_cells = (0..flow.cell_count())
                .map(|i| {
                    flow.cell(i)
                        .unwrap()
                        .primitive(flow.gas(), || false)
                        .unwrap()
                })
                .collect::<Vec<_>>();
            let mut planes = walls(0.0);
            let original_planes = planes;
            let force = if density_force {
                BodyForce3d::UniformForceDensity(value)
            } else {
                BodyForce3d::UniformAcceleration(value)
            };
            let r = flow
                .step_with_inertial_walls_and_forcing(
                    axis,
                    &mut planes,
                    Time::new(0.001),
                    StepControl::default(),
                    loads(),
                    controls(0.0),
                    force,
                    CAP,
                    &mut scratch,
                    || false,
                )
                .unwrap();
            let m = r.mechanics;
            let gas = m.gas;
            let dt = gas.elapsed.value();
            let half = 0.5 * dt;
            assert!(m.combined_momentum_residual_n_s.abs() <= controls(0.0).momentum_tolerance_n_s);
            assert!(m.total_energy_residual_j.abs() <= controls(0.0).energy_tolerance_j);
            for residual in gas.conservation_residual.species_masses {
                close(residual, 0.0, 5e-12);
            }
            for residual in gas.conservation_residual.momentum {
                close(residual, 0.0, 5e-12);
            }
            close(gas.conservation_residual.energy, 0.0, 6e-12);
            close(gas.absolute_energy_residual_j, 0.0, 3e-10);
            let measure = if density_force {
                0.5 * (old_grid.cell_volume() + flow.grid().cell_volume())
                    * flow.cell_count() as f64
            } else {
                gas.before.mass()
            };
            for a in 0..3 {
                close(gas.source_change.impulse[a], value[a] * measure * dt, 3e-12);
            }
            let exact_work: f64 = if density_force {
                before_cells
                    .iter()
                    .map(|p| {
                        (0..3)
                            .map(|a| {
                                value[a]
                                    * half
                                    * old_grid.cell_volume()
                                    * (p.velocity[a] + 0.5 * value[a] * half / p.density)
                            })
                            .sum::<f64>()
                    })
                    .sum::<f64>()
                    + (0..flow.cell_count())
                        .map(|i| {
                            let p = flow
                                .cell(i)
                                .unwrap()
                                .primitive(flow.gas(), || false)
                                .unwrap();
                            (0..3)
                                .map(|a| {
                                    value[a]
                                        * half
                                        * flow.grid().cell_volume()
                                        * (p.velocity[a] - 0.5 * value[a] * half / p.density)
                                })
                                .sum::<f64>()
                        })
                        .sum::<f64>()
            } else {
                (0..3)
                    .map(|a| half * value[a] * (gas.before.momentum[a] + gas.after.momentum[a]))
                    .sum()
            };
            close(gas.source_change.work, exact_work, 3e-12);
            let sides = [BoundarySide::ALL[2 * axis], BoundarySide::ALL[2 * axis + 1]];
            let area = old_grid.face_area(sides[0]) * old_grid.face_count(sides[0]) as f64;
            for end in 0..2 {
                let sign = if end == 0 { 1.0 } else { -1.0 };
                let external_impulse = (loads()[end].external_force_n
                    + sign * loads()[end].exterior_pressure_pa * area)
                    * dt;
                close(
                    planes[end].mass_kg()
                        * (planes[end].velocity_m_s() - original_planes[end].velocity_m_s())
                        - m.gas_impulse_n_s[end]
                        - external_impulse,
                    0.0,
                    controls(0.0).momentum_tolerance_n_s,
                );
            }
            let wall_change: f64 = (0..2)
                .map(|i| {
                    m.kinetic_energy_change_j[i] + m.exterior_pressure_work_j[i]
                        - m.external_force_work_j[i]
                })
                .sum();
            let referenced =
                gas.after.energy - gas.before.energy + wall_change - gas.source_change.work;
            assert!(referenced.abs() <= controls(0.0).energy_tolerance_j);
            let reference_change: f64 = (0..flow.gas().species_count())
                .map(|k| {
                    (gas.after.species_masses[k] - gas.before.species_masses[k])
                        * flow.gas().energy_references()[k]
                })
                .sum();
            assert!((referenced + reference_change).abs() <= controls(0.0).energy_tolerance_j);
            assert_eq!(flow.workspace_bytes() + scratch.workspace_bytes(), bytes);
            assert_eq!(flow.last_step(), Some(&gas));
            for side in BoundarySide::ALL {
                for q in flow.boundary_impulses(side) {
                    assert_eq!(q.species_masses, &[0.0; 3]);
                }
            }
        }
    }
}

#[test]
fn g4_source_piston_cancellation_and_failed_roots_leave_all_accepted_owners_intact() {
    let force = BodyForce3d::UniformAcceleration([0.2, -0.15, 0.1]);
    let make = || {
        let mut flow = case();
        let mut scratch = work(&flow);
        let mut planes = walls(0.0);
        flow.step_with_inertial_walls_and_forcing(
            0,
            &mut planes,
            Time::new(0.0005),
            StepControl::default(),
            loads(),
            controls(0.0),
            force,
            CAP,
            &mut scratch,
            || false,
        )
        .unwrap();
        (flow, scratch, planes)
    };
    let (mut expected, mut expected_work, mut expected_planes) = make();
    let mut polls = 0;
    let accepted = expected
        .step_with_inertial_walls_and_forcing(
            0,
            &mut expected_planes,
            Time::new(0.001),
            StepControl::default(),
            loads(),
            controls(0.0),
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
        let (mut flow, mut scratch, mut planes) = make();
        let states = flow.states().to_vec();
        let grid = *flow.grid();
        let time = flow.time();
        let last = *flow.last_step().unwrap();
        let boundary = flow.boundary_impulses.clone();
        let old_planes = planes;
        let mut seen = 0;
        assert_eq!(
            flow.step_with_inertial_walls_and_forcing(
                0,
                &mut planes,
                Time::new(0.001),
                StepControl::default(),
                loads(),
                controls(0.0),
                force,
                CAP,
                &mut scratch,
                || {
                    seen += 1;
                    seen == cutoff
                },
            ),
            Err(PistonError::Cancelled)
        );
        assert_eq!(flow.states(), states);
        assert_eq!(*flow.grid(), grid);
        assert_eq!(flow.time(), time);
        assert_eq!(flow.last_step(), Some(&last));
        assert_eq!(flow.boundary_impulses, boundary);
        assert_eq!(planes, old_planes);
        let mut exhausted = controls(0.0);
        exhausted.maximum_iterations = 1;
        exhausted.momentum_tolerance_n_s = 1e-30;
        assert!(
            flow.step_with_inertial_walls_and_forcing(
                0,
                &mut planes,
                Time::new(0.001),
                StepControl {
                    maximum_halvings: 0,
                    ..StepControl::default()
                },
                loads(),
                exhausted,
                force,
                CAP,
                &mut scratch,
                || false,
            )
            .is_err()
        );
        assert_eq!(flow.states(), states);
        assert_eq!(*flow.grid(), grid);
        assert_eq!(planes, old_planes);
        assert_eq!(flow.last_step(), Some(&last));
        let replay = flow
            .step_with_inertial_walls_and_forcing(
                0,
                &mut planes,
                Time::new(0.001),
                StepControl::default(),
                loads(),
                controls(0.0),
                force,
                CAP,
                &mut scratch,
                || false,
            )
            .unwrap();
        assert_eq!(replay, accepted);
        assert_eq!(planes, expected_planes);
        assert_eq!(flow.states(), expected.states());
        assert_eq!(flow.boundary_impulses, expected.boundary_impulses);
    }
}
