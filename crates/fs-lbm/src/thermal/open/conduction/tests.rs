use super::*;
use crate::thermal::boundary::conduction::WallCell;
use crate::thermal::boundary::{FiniteThermalWall, ThermalWall};
use crate::thermal::open::{ThermalOpenBalance3, ThermalOpenBoundaries3, ThermalOpenLimits3};

fn controls() -> WallConductionControls {
    WallConductionControls {
        max_nodes: 256,
        max_links: 512,
        max_factor_bytes: 4 * 1024 * 1024,
        balance_tolerance: 1e-9,
    }
}
fn edge(a: (usize, usize), b: (usize, usize), conductance: f64) -> WallConductionLink {
    WallConductionLink {
        a: WallCell {
            side: a.0,
            link: a.1,
        },
        b: WallCell {
            side: b.0,
            link: b.1,
        },
        conductance,
    }
}
fn make(
    dims: [usize; 3],
    normal_axes: &[usize],
    velocity: [f64; 3],
    solid: &[bool],
) -> ThermalOpenField3 {
    let cells: usize = dims.into_iter().product();
    let mut faces = [ThermalOpenFace3::Periodic; 6];
    let flow_axis = (0..3).find(|axis| !normal_axes.contains(axis)).unwrap();
    faces[2 * flow_axis] = ThermalOpenFace3::Reservoir { temperature: 0.4 };
    faces[2 * flow_axis + 1] = ThermalOpenFace3::Reservoir { temperature: 0.8 };
    for &axis in normal_axes {
        for side in 0..2 {
            faces[2 * axis + side] =
                ThermalOpenFace3::Wall(ThermalWall::FiniteSolid(FiniteThermalWall {
                    temperature: if side == 0 { 1.2 } else { 0.65 },
                    heat_capacity: if side == 0 { 0.05 } else { 0.07 },
                    resistance: if side == 0 { 0.3 } else { 0.6 },
                    heat_per_step: 0.0001,
                }));
        }
    }
    let initial: Vec<_> = (0..cells).map(|cell| 0.4 + 0.0003 * cell as f64).collect();
    ThermalOpenField3::try_new::<ThermalOpenError3>(
        dims,
        solid,
        &initial,
        0.85,
        ThermalOpenBoundaries3 {
            faces,
            internal_wall: ThermalWall::Adiabatic,
        },
        ThermalOpenLimits3 {
            maximum_cells: 4096,
            maximum_storage_bytes: 4 * 1024 * 1024,
            temperature_range: [-1.0, 4.0],
            balance_tolerance: 1e-9,
        },
        0,
        |_| Ok(velocity),
        || Ok(()),
    )
    .unwrap()
}
fn step(field: &mut ThermalOpenField3, heat: &[f64]) -> ThermalOpenBalance3 {
    let candidate = field
        .begin_step::<ThermalOpenError3>(heat, None, || Ok(()))
        .unwrap();
    let balance = field
        .finish_step::<ThermalOpenError3>(heat, candidate, || Ok(()))
        .unwrap();
    field.publish(field.completed_steps() + 1, balance);
    balance
}
fn solid_temperature(field: &ThermalOpenField3, axis: usize, side: usize, link: usize) -> f64 {
    let ThermalWall::FiniteSolid(solid) = field.wall_states(2 * axis + side).unwrap()[link] else {
        panic!("finite cell");
    };
    solid.temperature
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

#[test]
fn noncubic_rotated_face_pairs_sample_the_actual_outgoing_populations() {
    let dims = [5, 4, 3];
    let cells: usize = dims.into_iter().product();
    let velocity = [0.003, 0.006, 0.009];
    let heat: Vec<_> = (0..cells)
        .map(|cell| 0.0002 * (1 + cell % 5) as f64)
        .collect();
    // Explicit Cartesian endpoints independently check all surface inverses.
    let endpoints = [
        ([0, 1, 0], [4, 2, 2]),
        ([1, 0, 0], [3, 3, 2]),
        ([1, 0, 0], [3, 3, 2]),
    ];
    for axis in 0..3 {
        let plane = cells / dims[axis];
        let mut field = make(dims, &[axis], velocity, &vec![false; cells]);
        field
            .set_wall_conduction(axis, &[edge((0, 1), (1, plane - 2), 0.7)], controls())
            .unwrap();
        let address = |p: [usize; 3]| p[0] + dims[0] * (p[1] + dims[1] * p[2]);
        let a = address(endpoints[axis].0);
        let b = address(endpoints[axis].1);
        // Fresh history has no correction; the equilibrium first source kick
        // gives these exact postcollision normal traces before any streaming.
        let ta = (0.4 + 0.0003 * a as f64 + 0.5 * heat[a]) * (1.0 - 4.0 * velocity[axis]);
        let tb = (0.4 + 0.0003 * b as f64 + 0.5 * heat[b]) * (1.0 + 4.0 * velocity[axis]);
        let ka = 1.0 / (0.3 + 4.0);
        let kb = 1.0 / (0.6 + 4.0);
        let aa = 0.05 + ka + 0.7;
        let bb = 0.07 + kb + 0.7;
        let rhs_a = 0.05 * 1.2 + 0.0001 + ka * ta;
        let rhs_b = 0.07 * 0.65 + 0.0001 + kb * tb;
        let determinant = aa * bb - 0.7 * 0.7;
        let expected_a = (bb * rhs_a + 0.7 * rhs_b) / determinant;
        let expected_b = (aa * rhs_b + 0.7 * rhs_a) / determinant;
        let receipt = step(&mut field, &heat);
        near(solid_temperature(&field, axis, 0, 1), expected_a, 3e-14);
        near(
            solid_temperature(&field, axis, 1, plane - 2),
            expected_b,
            3e-14,
        );
        let lower = receipt.walls[2 * axis].unwrap().solid.unwrap();
        let upper = receipt.walls[2 * axis + 1].unwrap().solid.unwrap();
        near(lower.conducted_in, 0.7 * (expected_b - expected_a), 3e-14);
        near(lower.conducted_in + upper.conducted_in, 0.0, 1e-16);
        near(receipt.combined_defect, 0.0, 3e-13);
    }
}

fn lateral_links(plane: usize, width: usize, conductance: f64) -> Vec<WallConductionLink> {
    let mut links = Vec::new();
    for side in 0..2 {
        for link in 0..plane {
            if link % width + 1 < width {
                links.push(edge((side, link), (side, link + 1), conductance));
            }
            if link + width < plane {
                links.push(edge((side, link), (side, link + width), conductance));
            }
        }
    }
    links
}

#[test]
fn stiff_lateral_walls_share_heat_in_a_sourced_throughflow_and_close_energy() {
    let dims = [8, 4, 3];
    let cells: usize = dims.into_iter().product();
    let velocity = [0.025, 0.0, 0.0];
    let mut connected = make(dims, &[1, 2], velocity, &vec![false; cells]);
    let mut independent = make(dims, &[1, 2], velocity, &vec![false; cells]);
    for axis in [1, 2] {
        let plane = cells / dims[axis];
        // G/C exceeds 10,000: a same-step explicit solid update is unstable.
        connected
            .set_wall_conduction(axis, &lateral_links(plane, dims[0], 1000.0), controls())
            .unwrap();
    }
    let heat: Vec<_> = (0..cells)
        .map(|cell| 0.0002 * (1 + cell % dims[0]) as f64)
        .collect();
    for _ in 0..30 {
        let before: f64 = connected.populations().iter().flatten().sum();
        let receipt = step(&mut connected, &heat);
        step(&mut independent, &heat);
        let after: f64 = connected.populations().iter().flatten().sum();
        let mut stored = 0.0;
        let mut source = 0.0;
        let mut conducted = 0.0;
        for face in 2..6 {
            let solid = receipt.walls[face].unwrap().solid.unwrap();
            stored += solid.stored_change;
            source += solid.source;
            conducted += solid.conducted_in;
        }
        near(conducted, 0.0, 1e-12);
        near(
            after - before + stored,
            heat.iter().sum::<f64>()
                + source
                + receipt.reservoirs.iter().map(|r| r.heat_in).sum::<f64>(),
            3e-10,
        );
        near(receipt.combined_defect, 0.0, 3e-10);
    }
    for axis in [1, 2] {
        for side in 0..2 {
            let spread = |field: &ThermalOpenField3| {
                let values: Vec<_> = (0..cells / dims[axis])
                    .map(|link| solid_temperature(field, axis, side, link))
                    .collect();
                values.iter().copied().fold(f64::NEG_INFINITY, f64::max)
                    - values.iter().copied().fold(f64::INFINITY, f64::min)
            };
            assert!(
                spread(&connected) < 0.02 * spread(&independent),
                "stiff solid conduction did not redistribute wall heat"
            );
        }
    }
}

#[test]
fn failed_step_does_not_publish_conducted_walls_and_retry_rebuilds_graph_scratch() {
    let dims = [6, 4, 3];
    let cells: usize = dims.into_iter().product();
    let make_connected = || {
        let mut field = make(dims, &[1], [0.02, 0.0, 0.0], &vec![false; cells]);
        field
            .set_wall_conduction(
                1,
                &lateral_links(cells / dims[1], dims[0], 17.0),
                controls(),
            )
            .unwrap();
        field
    };
    let mut actual = make_connected();
    let mut expected = make_connected();
    let heat = vec![0.001; cells];
    step(&mut actual, &heat);
    step(&mut expected, &heat);
    let before_populations = actual.populations.clone();
    let before_walls = actual.walls.clone();
    let before_balance = actual.last_balance;
    let larger_heat = vec![0.03; cells];
    let candidate = actual
        .begin_step::<ThermalOpenError3>(&larger_heat, None, || Ok(()))
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        actual.finish_step::<ThermalOpenError3>(&larger_heat, candidate, || {
            polls += 1;
            if polls == 5 {
                Err(ThermalOpenError3::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(ThermalOpenError3::Cancelled)
    ));
    assert_eq!(actual.populations, before_populations);
    assert_eq!(actual.walls, before_walls);
    assert_eq!(actual.last_balance, before_balance);
    assert_eq!(actual.completed_steps(), 1);
    let retry = step(&mut actual, &heat);
    let control = step(&mut expected, &heat);
    assert_eq!(actual.populations, expected.populations);
    assert_eq!(actual.walls, expected.walls);
    assert_eq!(retry, control);
}

#[test]
fn inactive_endpoints_and_invalid_replacements_leave_prepared_plan_intact() {
    let dims = [5, 4, 3];
    let cells: usize = dims.into_iter().product();
    for axis in 0..3 {
        let mut solid = vec![false; cells];
        solid[boundary_cell(dims, axis, 0, 0)] = true;
        let mut field = make(dims, &[axis], [0.01, 0.0, 0.0], &solid);
        let valid = edge((0, 1), (0, 2), 0.2);
        field
            .set_wall_conduction(axis, &[valid], controls())
            .unwrap();
        let before = field.populations().to_vec();
        for conductance in [0.0, 0.3] {
            assert!(
                field
                    .set_wall_conduction(axis, &[edge((0, 0), (0, 1), conductance)], controls())
                    .is_err()
            );
            assert_eq!(field.wall_conduction(axis).unwrap().node_count(), 2);
            assert_eq!(field.wall_conduction(axis).unwrap().link_count(), 1);
        }
        let mut too_small = controls();
        too_small.max_nodes = 1;
        assert!(
            field
                .set_wall_conduction(axis, &[valid], too_small)
                .is_err()
        );
        assert_eq!(field.wall_conduction(axis).unwrap().node_count(), 2);
        assert_eq!(field.populations(), before);
        assert_eq!(field.completed_steps(), 0);
        assert!(field.set_wall_conduction(3, &[valid], controls()).is_err());
        field.clear_wall_conduction(axis).unwrap();
        assert!(field.wall_conduction(axis).is_none());
        assert_eq!(field.populations(), before);
    }
}
