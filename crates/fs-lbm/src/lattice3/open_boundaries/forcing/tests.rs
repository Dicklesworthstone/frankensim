use super::*;
use crate::d3q19::{CollisionModel3, FaceBoundary3};
use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27};
use crate::lattice3::{D3q19, Grid3, GridBudget3, Lattice3, StepData};
use core::marker::PhantomData;

mod acceleration;

const RHO: f64 = 1.2;
const VELOCITY: [f64; 3] = [0.012, -0.008, 0.009];

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

fn fixture<L: EquilibriumLattice3<Q>, const Q: usize>(
    dims: [usize; 3],
    collision: L::Collision,
    force: [f64; 3],
    axis: Option<usize>,
    stress: bool,
    solid: Option<[usize; 3]>,
) -> OpenBoundaryGrid3<L, Q> {
    let mut periodic = [true; 3];
    let mut faces = [FaceBoundary3::Periodic; 6];
    if let Some(axis) = axis {
        periodic[axis] = false;
        faces[2 * axis] = FaceBoundary3::Velocity { velocity: VELOCITY };
        faces[2 * axis + 1] = FaceBoundary3::Pressure { density: RHO };
    }
    let grid = Grid3::try_new(
        dims,
        collision,
        force,
        periodic,
        GridBudget3 {
            maximum_cells: 1024,
            maximum_mach: 0.2,
        },
        |p| solid == Some(p),
        |_| {
            let mut f = initial::equilibrium::<L, Q>(collision, RHO, VELOCITY);
            if stress {
                let positive = L::E.iter().position(|&e| e == [1, 0, 0]).unwrap();
                let negative = L::OPP[positive];
                f[positive] += 0.0002;
                f[negative] += 0.0002;
                f[0] -= 0.0004;
            }
            initial::attach_force::<L, Q>(&mut f, collision, VELOCITY, force);
            f
        },
    )
    .unwrap();
    OpenBoundaryGrid3::try_new_with_poll(grid, BoundarySpec3::new(faces), &mut || Ok(())).unwrap()
}

// Independent hydrodynamic identities, including the full symmetric stress
// correction; this does not reconstruct or collide any populations.
fn physical<L: Lattice3<Q>, const Q: usize>(
    f: &[f64; Q],
    force: [f64; 3],
) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let rho: f64 = f.iter().sum();
    let velocity: [f64; 3] = core::array::from_fn(|a| {
        (f.iter()
            .zip(L::E)
            .map(|(f, e)| f * f64::from(e[a]))
            .sum::<f64>()
            + 0.5 * force[a])
            / rho
    });
    let stress = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            f.iter()
                .zip(L::E)
                .map(|(f, e)| f * f64::from(e[a] * e[b]))
                .sum::<f64>()
                + 0.5 * (velocity[a] * force[b] + force[a] * velocity[b])
                - rho * (velocity[a] * velocity[b] + if a == b { L::CS2 } else { 0.0 })
        })
    });
    (rho, velocity, stress)
}

fn totals<L: EquilibriumLattice3<Q>, const Q: usize>(
    flow: &ForcedOpenBoundaryGrid3<L, Q>,
) -> (f64, [f64; 3]) {
    let mut mass = 0.0;
    let mut momentum = [0.0; 3];
    for cell in 0..flow.populations().len() {
        if flow.inner.grid.solid[cell] {
            continue;
        }
        let (rho, velocity) = flow.macroscopic(cell).unwrap();
        mass += rho;
        for axis in 0..3 {
            momentum[axis] += rho * velocity[axis];
        }
    }
    (mass, momentum)
}

fn body_impulse(start: &[[f64; 3]], end: &[[f64; 3]]) -> [f64; 3] {
    core::array::from_fn(|axis| {
        start
            .iter()
            .zip(end)
            .map(|(a, b)| 0.5 * a[axis] + 0.5 * b[axis])
            .sum()
    })
}

fn variable_force(cell: usize, epoch: usize) -> [f64; 3] {
    let sign = if epoch % 2 == 0 { 1.0 } else { -1.0 };
    [
        sign * 0.00013 * (1 + cell % 3) as f64,
        -sign * 0.00009 * (1 + cell % 5) as f64,
        sign * 0.00007 * (1 + cell % 7) as f64,
    ]
}

fn exercise<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let dims = [5, 4, 3];
    let cells: usize = dims.into_iter().product();
    let old_force = [0.0001, -0.00007, 0.00004];

    // Attach nonuniform force to a genuinely stressed state, preserving all
    // physical moments. A replacement by equilibrium fails this oracle.
    let original = fixture::<L, Q>(dims, collision, old_force, None, true, None);
    let before: Vec<_> = original
        .populations()
        .iter()
        .map(|f| physical::<L, Q>(f, old_force))
        .collect();
    let force: Vec<_> = (0..cells)
        .map(|cell| {
            if cell % 4 == 0 {
                [0.0; 3]
            } else {
                variable_force(cell, cell % 2)
            }
        })
        .collect();
    let mut spatial = ForcedOpenBoundaryGrid3::try_new_preserving_velocity_with_poll(
        original,
        force,
        &mut || Ok(()),
    )
    .unwrap();
    assert_eq!(spatial.completed_steps(), 0);
    for (cell, &(rho, velocity, stress)) in before.iter().enumerate() {
        let (actual_rho, actual_velocity, actual_stress) =
            physical::<L, Q>(&spatial.populations()[cell], spatial.forces()[cell]);
        near(actual_rho, rho, 3e-15);
        for a in 0..3 {
            near(actual_velocity[a], velocity[a], 3e-15);
            near(spatial.macroscopic(cell).unwrap().1[a], velocity[a], 3e-15);
            for b in 0..3 {
                near(actual_stress[a][b], stress[a][b], 3e-15);
            }
        }
    }

    // Actual periodic native collide/pull steps close the independently
    // measured trapezoidal impulse while the force reverses between epochs.
    for epoch in 1..=3 {
        let mut end: Vec<_> = (0..cells).map(|cell| variable_force(cell, epoch)).collect();
        let expected = body_impulse(spatial.forces(), &end);
        let (mass_before, momentum_before) = totals(&spatial);
        let (receipt, step) = spatial
            .stage_serial_with_poll(&end, None, [0.0; 3], &mut || Ok(()))
            .unwrap();
        spatial.publish_with_force_swap(step, &mut end);
        let (mass_after, momentum_after) = totals(&spatial);
        near(mass_after, mass_before, 3e-12);
        assert_eq!(receipt.reconstructed_cells, 0);
        for axis in 0..3 {
            near(
                momentum_after[axis] - momentum_before[axis],
                expected[axis],
                3e-13,
            );
        }
    }

    for axis in 0..3 {
        // A constant field must retain the complete original uniform-force
        // arithmetic, including native boundary populations and receipts.
        let mut uniform = fixture::<L, Q>(dims, collision, old_force, Some(axis), true, None);
        let original_bits = uniform.populations().to_vec();
        let same = fixture::<L, Q>(dims, collision, old_force, Some(axis), true, None);
        let mut held = ForcedOpenBoundaryGrid3::try_new_preserving_velocity_with_poll(
            same,
            vec![old_force; cells],
            &mut || Ok(()),
        )
        .unwrap();
        assert_eq!(held.populations(), original_bits);
        let expected = uniform
            .stage_serial_with_forces(
                uniform.plan,
                None,
                None,
                None,
                [0.5, -0.75, 1.25],
                &mut || Ok(()),
            )
            .unwrap();
        let end = held.forces.clone();
        let actual = held
            .stage_serial_with_poll(&end, None, [0.5, -0.75, 1.25], &mut || Ok(()))
            .unwrap();
        assert_eq!(held.inner.grid.next, uniform.grid.next);
        assert_eq!(actual, expected);

        // Distinct donor/target endpoint forces independently test BOTH half
        // source frames at real velocity and pressure reservoir faces.
        let original = fixture::<L, Q>(dims, collision, [0.0; 3], Some(axis), true, None);
        let start: Vec<_> = (0..cells).map(|cell| variable_force(cell, 0)).collect();
        let mut flow = ForcedOpenBoundaryGrid3::try_new_preserving_velocity_with_poll(
            original,
            start,
            &mut || Ok(()),
        )
        .unwrap();
        let mut end: Vec<_> = (0..cells).map(|cell| variable_force(cell + 1, 1)).collect();
        let expected_body = body_impulse(flow.forces(), &end);
        let (mass_before, momentum_before) = totals(&flow);
        let (receipt, step) = flow
            .stage_serial_with_poll(&end, None, [0.0; 3], &mut || Ok(()))
            .unwrap();
        let grid = &flow.inner.grid;
        let data = StepData::<L, Q> {
            dims: grid.dims,
            periodic: grid.periodic,
            solid: &grid.solid,
            collision: grid.collision,
            collisions: None,
            force: grid.force,
            forces: Some(&end),
            budget: grid.budget,
            marker: PhantomData,
        };
        for cell in flow.inner.plan.cells() {
            let (face, donor) = flow.inner.plan.locate(cell).unwrap();
            let mut source = [0.0; Q];
            data.pull_cell(&grid.post, &mut source, donor);
            let (rho_source, velocity_source, stress_source) =
                physical::<L, Q>(&source, end[donor]);
            let (rho, velocity, stress) = physical::<L, Q>(&grid.next[cell], end[cell]);
            let (target_rho, target_velocity) = match flow.boundaries().face(face) {
                FaceBoundary3::Velocity { velocity } => (rho_source, velocity),
                FaceBoundary3::Pressure { density } => (density, velocity_source),
                _ => unreachable!(),
            };
            near(rho, target_rho, 5e-13);
            for a in 0..3 {
                near(velocity[a], target_velocity[a], 5e-14);
                for b in 0..3 {
                    near(stress[a][b], stress_source[a][b], 5e-13);
                }
            }
        }
        flow.publish_with_force_swap(step, &mut end);
        let (mass_after, momentum_after) = totals(&flow);
        near(mass_after - mass_before, receipt.mass_change, 3e-12);
        for a in 0..3 {
            near(
                momentum_after[a] - momentum_before[a],
                expected_body[a] + receipt.impulse[a],
                3e-13,
            );
        }
    }
}

#[test]
fn d3q19_families_preserve_force_frames_trapezoidal_impulse_and_open_targets() {
    for collision in [
        CollisionModel3::Bgk { tau: 0.9 },
        CollisionModel3::Trt {
            tau: 0.9,
            magic: 3.0 / 16.0,
        },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.2,
            higher_order_rate: 1.6,
        },
        CollisionModel3::CentralMomentBulk {
            shear_rate: 1.2,
            bulk_rate: 0.9,
            higher_order_rate: 1.6,
        },
        CollisionModel3::ReducedCumulant {
            second_order_rate: 1.2,
            third_order_rate: 1.1,
            fourth_order_rate: 0.7,
        },
    ] {
        exercise::<D3q19, 19>(collision);
    }
}

#[test]
fn d3q27_families_preserve_force_frames_trapezoidal_impulse_and_open_targets() {
    for equilibrium in [Equilibrium27::Quadratic, Equilibrium27::Factorized] {
        for model in [
            CollisionModel27::Bgk {
                tau: 0.9,
                equilibrium,
            },
            CollisionModel27::Trt {
                tau: 0.9,
                magic: 3.0 / 16.0,
                equilibrium,
            },
        ] {
            exercise::<D3q27, 27>(CollisionConfig27 {
                model,
                maximum_mach: 0.2,
            });
        }
    }
    for model in [
        CollisionModel27::CentralMoment {
            rates: Rates27 { omega: [1.2; 10] },
        },
        CollisionModel27::Cumulant {
            rates: Rates27 { omega: [1.2; 10] },
        },
    ] {
        exercise::<D3q27, 27>(CollisionConfig27 {
            model,
            maximum_mach: 0.2,
        });
    }
}

#[test]
fn late_cancellation_and_endpoint_mach_refusal_preserve_source_epoch_and_retry_exactly() {
    let dims = [6, 5, 4];
    let cells: usize = dims.into_iter().product();
    let collision = CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    };
    let make = || {
        let inner = fixture::<D3q27, 27>(dims, collision, [0.0; 3], Some(0), true, None);
        ForcedOpenBoundaryGrid3::try_new_preserving_velocity_with_poll(
            inner,
            (0..cells).map(|cell| variable_force(cell, 0)).collect(),
            &mut || Ok(()),
        )
        .unwrap()
    };
    let mut actual = make();
    let mut control = make();
    let before_populations = actual.populations().to_vec();
    let before_forces = actual.forces().to_vec();
    let mut end: Vec<_> = (0..cells).map(|cell| variable_force(cell, 1)).collect();
    let mut count = 0;
    let (expected, expected_step) = control
        .stage_serial_with_poll(&end, None, [0.0; 3], &mut || {
            count += 1;
            Ok(())
        })
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        actual.stage_serial_with_poll(&end, None, [0.0; 3], &mut || {
            polls += 1;
            if polls == count {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(GridError3::Cancelled)
    ));
    assert_eq!(actual.populations(), before_populations);
    assert_eq!(actual.forces(), before_forces);
    assert_eq!(actual.completed_steps(), 0);
    let mut impossible = end.clone();
    impossible[3 + dims[0] * (2 + dims[1] * 2)] = [1.0, 0.0, 0.0];
    assert!(matches!(
        actual.stage_serial_with_poll(&impossible, None, [0.0; 3], &mut || Ok(())),
        Err(GridError3::Mach { .. })
    ));
    assert_eq!(actual.populations(), before_populations);
    assert_eq!(actual.forces(), before_forces);
    let (receipt, step) = actual
        .stage_serial_with_poll(&end, None, [0.0; 3], &mut || Ok(()))
        .unwrap();
    let mut control_end = end.clone();
    actual.publish_with_force_swap(step, &mut end);
    control.publish_with_force_swap(expected_step, &mut control_end);
    assert_eq!(receipt, expected);
    assert_eq!(actual.populations(), control.populations());
    assert_eq!(actual.forces(), control.forces());
    assert_eq!(
        end, before_forces,
        "accepted force storage must return as reusable scratch"
    );
    assert_eq!(actual.completed_steps(), 1);
}

#[test]
fn inactive_source_cells_and_invalid_force_fields_refuse_without_publishing() {
    let dims = [5, 4, 3];
    let cells: usize = dims.into_iter().product();
    let solid = [2, 2, 1];
    let solid_index = solid[0] + dims[0] * (solid[1] + dims[1] * solid[2]);
    let inner = fixture::<D3q19, 19>(
        dims,
        CollisionModel3::Bgk { tau: 0.9 },
        [0.0; 3],
        Some(0),
        false,
        Some(solid),
    );
    let solid_before = inner.populations()[solid_index];
    let mut force: Vec<_> = (0..cells).map(|cell| variable_force(cell, 0)).collect();
    force[solid_index] = [0.0; 3];
    let mut flow = ForcedOpenBoundaryGrid3::try_new_preserving_velocity_with_poll(
        inner,
        force,
        &mut || Ok(()),
    )
    .unwrap();
    assert_eq!(flow.populations()[solid_index], solid_before);
    let before = flow.populations().to_vec();
    let before_force = flow.forces().to_vec();
    let mut bad = before_force.clone();
    bad[solid_index][1] = 0.001;
    assert!(matches!(
        flow.stage_serial_with_poll(&bad, None, [0.0; 3], &mut || Ok(())),
        Err(GridError3::MaterialField { .. })
    ));
    bad[solid_index][1] = 0.0;
    bad[cells - 1][0] = f64::NAN;
    assert!(
        flow.stage_serial_with_poll(&bad, None, [0.0; 3], &mut || Ok(()))
            .is_err()
    );
    assert!(
        flow.stage_serial_with_poll(&before_force[..cells - 1], None, [0.0; 3], &mut || Ok(()))
            .is_err()
    );
    assert_eq!(flow.populations(), before);
    assert_eq!(flow.forces(), before_force);
    assert_eq!(flow.completed_steps(), 0);
}
