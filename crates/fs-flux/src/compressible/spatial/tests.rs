use super::*;
use crate::compressible as planar;

// Reuse the original manufactured NASA9 source; no alternate test-only EOS.
use crate::compressible::tests::{close, gas};

const PERIODIC: Boundaries3d = Boundaries3d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
    z: AxisBoundaries::Periodic,
};
const WALL: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Reflecting,
    upper: BoundaryCondition::Reflecting,
};
const OPEN: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Transmissive,
    upper: BoundaryCondition::Transmissive,
};
const WALLS: Boundaries3d = Boundaries3d {
    x: WALL,
    y: WALL,
    z: WALL,
};

fn field(
    gas: FrozenGas,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    sample: impl Fn([f64; 3]) -> [f64; 5],
) -> Euler3d {
    let states = (0..grid.len())
        .map(|i| {
            let p = sample(grid.cell_center(i).unwrap());
            State::from_primitive(&gas, p[0], [p[1], p[2], p[3]], p[4]).unwrap()
        })
        .collect();
    Euler3d::new(gas, grid, states, boundaries, || false).unwrap()
}

fn balance(report: &StepReport) {
    assert!(report.maximum_cfl <= 0.40000000000001);
    for r in report.conservation_residual.array() {
        close(r, 0.0, 6e-12);
    }
}

fn run(flow: &mut Euler3d, end: f64) {
    let mut steps = 0;
    while flow.time().value() < end {
        let report = flow
            .step(
                Time::new(end - flow.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        balance(&report);
        steps += 1;
        assert!(steps <= 2000);
    }
}

#[test]
fn g0_spatial_uniform_state_and_physical_face_geometry() {
    let grid = CartesianGrid::new(5, 4, 3, 0.2, 0.3, 0.4).unwrap();
    let mut flow = field(gas(0.4), grid, PERIODIC, |_| [1.1, 0.5, -0.3, 0.7, 1.7]);
    let original = flow.states().to_vec();
    let total = flow.totals();
    close(total.mass, 1.1 * 5.0 * 4.0 * 3.0 * 0.2 * 0.3 * 0.4, 1e-14);
    assert_eq!(grid.face_count(BoundarySide::XMin), 12);
    assert_eq!(grid.face_count(BoundarySide::YMax), 15);
    assert_eq!(grid.face_count(BoundarySide::ZMin), 20);
    close(grid.face_area(BoundarySide::XMin), 0.12, 1e-15);
    close(grid.face_area(BoundarySide::YMax), 0.08, 1e-15);
    close(grid.face_area(BoundarySide::ZMin), 0.06, 1e-15);
    for (actual, expected) in grid
        .face_center(BoundarySide::XMax, 5)
        .unwrap()
        .into_iter()
        .zip([1.0, 0.45, 0.6])
    {
        close(actual, expected, 1e-15);
    }
    assert!(grid.cell_center(grid.len()).is_none());
    assert!(grid.face_center(BoundarySide::YMin, 15).is_none());
    run(&mut flow, 0.1);
    assert_eq!(flow.states(), original);
    assert_eq!(flow.totals(), total);
    for side in BoundarySide::ALL {
        assert_eq!(flow.boundary_impulses(side).len(), grid.face_count(side));
    }
}

#[test]
fn g1_spatial_oblique_contact_second_order_refinement() {
    let end = 0.08;
    let mut errors = Vec::new();
    let resolutions = [12, 24, 36];
    for n in resolutions {
        let h = 1.0 / n as f64;
        let z = std::f64::consts::PI * h;
        let average = (z.sin() / z).powi(3);
        let exact = |x: [f64; 3], time: f64| {
            1.0 + 0.12
                * average
                * (2.0 * std::f64::consts::PI * (x[0] + x[1] + x[2] - 0.35 * time)).sin()
        };
        let grid = CartesianGrid::new(n, n, n, h, h, h).unwrap();
        let mut flow = field(gas(0.0), grid, PERIODIC, |x| {
            [exact(x, 0.0), 0.3, -0.15, 0.2, 1.0]
        });
        run(&mut flow, end);
        let error = flow
            .states()
            .iter()
            .enumerate()
            .map(|(i, q)| (q.density - exact(grid.cell_center(i).unwrap(), end)).abs())
            .sum::<f64>()
            / grid.len() as f64;
        errors.push(error);
        for q in flow.states() {
            let p = q.primitive(flow.gas(), || false).unwrap();
            close(p.pressure, 1.0, 4e-11);
            for (u, v) in p.velocity.into_iter().zip([0.3, -0.15, 0.2]) {
                close(u, v, 3e-11);
            }
        }
    }
    eprintln!("spatial contact density L1: {errors:?}");
    for i in 0..2 {
        let order = (errors[i] / errors[i + 1]).ln()
            / (resolutions[i + 1] as f64 / resolutions[i] as f64).ln();
        eprintln!("spatial contact refinement order: {order:.6}");
        assert!(
            order > 1.8,
            "resolved contact refinement: {errors:?}, order={order}"
        );
    }
}

#[test]
fn g2_spatial_sod_reduces_to_original_one_dimensional_source_on_every_axis() {
    use fs_gas::dynamics::{Boundaries, Boundary, Euler1d, Mesh1d};
    use fs_qty::Length;
    let n = 32;
    let h = 1.0 / n as f64;
    for axis in 0..3 {
        let source = gas(0.0);
        let mut shape = [2; 3];
        shape[axis] = n;
        let mut spacing = [0.3; 3];
        spacing[axis] = h;
        let grid = CartesianGrid::new(
            shape[0], shape[1], shape[2], spacing[0], spacing[1], spacing[2],
        )
        .unwrap();
        let mut pairs = [AxisBoundaries::Periodic; 3];
        pairs[axis] = OPEN;
        let boundaries = Boundaries3d {
            x: pairs[0],
            y: pairs[1],
            z: pairs[2],
        };
        let mut flow = field(source.clone(), grid, boundaries, |x| {
            if x[axis] < 0.5 {
                [1.0, 0.0, 0.0, 0.0, 1.0]
            } else {
                [0.125, 0.0, 0.0, 0.0, 0.1]
            }
        });
        let mesh = Mesh1d::uniform(n, Length::new(h), 0.36).unwrap();
        let states = (0..n)
            .map(|i| {
                if i < n / 2 {
                    source.conserved(1.0, 0.0, 1.0).unwrap()
                } else {
                    source.conserved(0.125, 0.0, 0.1).unwrap()
                }
            })
            .collect();
        let mut line = Euler1d::new_on_mesh(
            source,
            mesh,
            states,
            Boundaries::Ends {
                left: Boundary::Transmissive,
                right: Boundary::Transmissive,
            },
            || false,
        )
        .unwrap();
        while flow.time().value() < 0.08 {
            let dt = (0.08 - flow.time().value()).min(0.001);
            let report = flow
                .step(Time::new(dt), StepControl::default(), || false)
                .unwrap();
            let line_report = line
                .step(report.elapsed, StepControl::default(), || false)
                .unwrap();
            close(line_report.elapsed.value(), report.elapsed.value(), 1e-15);
            balance(&report);
        }
        for (i, q) in flow.states().iter().enumerate() {
            let reference = line.states()[grid.coordinates(i)[axis]];
            close(q.density, reference.density, 5e-11);
            close(q.momentum[axis], reference.momentum, 5e-11);
            close(q.energy, reference.energy, 5e-11);
            for a in tangents(axis) {
                close(q.momentum[a], 0.0, 1e-13);
            }
        }
    }
}

#[test]
fn g2_spatial_extrusion_preserves_existing_planar_variable_caloric_solution() {
    let source = gas(0.3);
    let [nx, ny, nz] = [9, 7, 3];
    let [dx, dy, dz] = [0.12, 0.14, 0.2];
    let grid = CartesianGrid::new(nx, ny, nz, dx, dy, dz).unwrap();
    let planar_grid = planar::CartesianGrid::new(nx, ny, dx, dy, nz as f64 * dz).unwrap();
    let sample = |x: [f64; 3]| {
        let p = if x[0] + 0.4 * x[1] < 0.6 { 2.0 } else { 0.6 };
        [1.0, 0.2, -0.1, 0.0, p]
    };
    let mut flow = field(source.clone(), grid, PERIODIC, sample);
    let states = (0..nx * ny)
        .map(|i| {
            let x = planar_grid.cell_center(i).unwrap();
            let p = sample([x[0], x[1], 0.0]);
            planar::State::from_primitive(&source, p[0], [p[1], p[2]], p[4]).unwrap()
        })
        .collect();
    let mut two = planar::Euler2d::new(
        source,
        planar_grid,
        states,
        planar::Boundaries2d {
            x: AxisBoundaries::Periodic,
            y: AxisBoundaries::Periodic,
        },
        || false,
    )
    .unwrap();
    for _ in 0..12 {
        let r = flow
            .step(Time::new(0.002), StepControl::default(), || false)
            .unwrap();
        let p = two
            .step(r.elapsed, StepControl::default(), || false)
            .unwrap();
        close(r.elapsed.value(), p.elapsed.value(), 1e-15);
        balance(&r);
    }
    for (i, q) in flow.states().iter().enumerate() {
        let p = two.states()[i % (nx * ny)];
        close(q.density, p.density, 8e-11);
        close(q.energy, p.energy, 8e-11);
        close(q.momentum[0], p.momentum[0], 8e-11);
        close(q.momentum[1], p.momentum[1], 8e-11);
        close(q.momentum[2], 0.0, 1e-13);
    }
    close(flow.totals().mass, two.totals().mass, 4e-13);
    close(flow.totals().energy, two.totals().energy, 4e-13);
}

#[test]
fn g2_spatial_spherical_blast_conserves_and_preserves_axis_symmetry() {
    let n = 12;
    let h = 1.0 / n as f64;
    let grid = CartesianGrid::new(n, n, n, h, h, h).unwrap();
    let mut flow = field(gas(0.2), grid, WALLS, |x| {
        let r2 = x.into_iter().map(|s| (s - 0.5).powi(2)).sum::<f64>();
        [
            1.0,
            0.0,
            0.0,
            0.0,
            if r2 < 0.19_f64.powi(2) { 4.0 } else { 1.0 },
        ]
    });
    let before = flow.totals();
    run(&mut flow, 0.07);
    close(flow.totals().mass, before.mass, 3e-13);
    close(flow.totals().energy, before.energy, 3e-13);
    for (index, q) in flow.states().iter().enumerate() {
        let c = grid.coordinates(index);
        let rotated = flow.states()[grid.index([c[1], c[2], c[0]])];
        close(q.density, rotated.density, 4e-11);
        close(q.energy, rotated.energy, 4e-11);
        for a in 0..3 {
            close(q.momentum[(a + 1) % 3], rotated.momentum[a], 4e-11);
        }
        let mirrored = flow.states()[grid.index([n - 1 - c[0], c[1], c[2]])];
        close(q.density, mirrored.density, 4e-11);
        close(q.momentum[0], -mirrored.momentum[0], 4e-11);
    }
    for side in BoundarySide::ALL {
        for load in flow.boundary_impulses(side) {
            assert_eq!(load.mass, 0.0);
            assert_eq!(load.energy, 0.0);
            assert!(side.sign() * load.momentum[side.axis()] > 0.0);
            for a in tangents(side.axis()) {
                assert_eq!(load.momentum[a], 0.0);
            }
        }
    }
}

#[test]
fn g0_spatial_wall_impulse_and_cancelled_retries_preserve_accepted_state() {
    let grid = CartesianGrid::new(5, 4, 3, 0.2, 0.25, 0.3).unwrap();
    let make = || field(gas(0.2), grid, WALLS, |x| [1.0, 0.3, -0.2, 0.4, 1.0 + x[0]]);
    let mut flow = make();
    let mut replay = make();
    for f in [&mut flow, &mut replay] {
        let r = f
            .step(Time::new(0.004), StepControl::default(), || false)
            .unwrap();
        balance(&r);
        let mut support = [0.0; 3];
        for side in BoundarySide::ALL {
            for impulse in f.boundary_impulses(side) {
                for (s, v) in support.iter_mut().zip(impulse.momentum) {
                    *s += v;
                }
            }
        }
        for a in 0..3 {
            close(
                r.after.momentum[a] - r.before.momentum[a] + support[a],
                0.0,
                8e-13,
            );
        }
        close(r.after.mass, r.before.mass, 1e-14);
        close(r.after.energy, r.before.energy, 1e-14);
    }
    let state = flow.states().to_vec();
    let last = *flow.last_step().unwrap();
    let time = flow.time();
    let impulses: Vec<_> = BoundarySide::ALL
        .into_iter()
        .flat_map(|s| flow.boundary_impulses(s).to_vec())
        .collect();
    let mut polls = 0;
    let expected = replay
        .step(Time::new(0.004), StepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    assert_eq!(
        flow.step(Time::new(0.004), StepControl::default(), || {
            count += 1;
            count >= polls - 1
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(flow.states(), state);
    assert_eq!(flow.time(), time);
    assert_eq!(flow.last_step(), Some(&last));
    let after: Vec<_> = BoundarySide::ALL
        .into_iter()
        .flat_map(|s| flow.boundary_impulses(s).to_vec())
        .collect();
    assert_eq!(after, impulses);
    let actual = flow
        .step(Time::new(0.004), StepControl::default(), || false)
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(flow.states(), replay.states());
}

#[test]
fn g0_spatial_geometry_source_and_extensive_admission() {
    assert!(CartesianGrid::new(0, 1, 1, 1.0, 1.0, 1.0).is_err());
    assert!(CartesianGrid::new(usize::MAX, 2, 2, 1.0, 1.0, 1.0).is_err());
    assert!(CartesianGrid::new(MAX_EULER_CELLS + 1, 1, 1, 1.0, 1.0, 1.0).is_err());
    assert!(CartesianGrid::new(1, 1, 1, 1e-200, 1e-200, 1.0).is_err());
    assert!(CartesianGrid::new(1, 1, 1, 1e200, 1e200, 1.0).is_err());
    let source = gas(0.0);
    assert!(State::from_primitive(&source, 1.0, [0.0, 0.0, f64::NAN], 1.0).is_err());
    assert!(State::from_primitive(&source, 1.0, [0.0, 0.0, 0.0], 100.0).is_err());
    let grid = CartesianGrid::new(1, 1, 1, 1e100, 1e100, 1e100).unwrap();
    let q = State::from_primitive(&source, 1e100, [0.0; 3], 1e100).unwrap();
    assert!(matches!(
        Euler3d::new(source, grid, vec![q], PERIODIC, || false),
        Err(Error::Invalid(_))
    ));
    let small = CartesianGrid::new(1, 1, 1, 1.0, 1.0, 1.0).unwrap();
    let mut flow = field(gas(0.0), small, PERIODIC, |_| [1.0, 0.2, 0.3, -0.1, 1.0]);
    assert!(
        flow.step(
            Time::new(0.1),
            StepControl {
                cfl: 0.46,
                ..StepControl::default()
            },
            || false
        )
        .is_err()
    );
    assert_eq!(flow.time().value(), 0.0);
}
