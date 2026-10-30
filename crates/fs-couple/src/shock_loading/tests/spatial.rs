//! Volumetric shock footprints retain both transverse moments in the body frame.
use super::*;
use fs_flux::compressible::spatial::{Boundaries3d, BoundarySide, CartesianGrid, Euler3d, State};
use fs_flux::compressible::{AxisBoundaries, BoundaryCondition};

const AXES: [[f64; 3]; 3] = [
    [2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0],
    [2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0],
    [-1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
];
const REFLECTING: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Reflecting,
    upper: BoundaryCondition::Reflecting,
};
const CLOSED: Boundaries3d = Boundaries3d {
    x: REFLECTING,
    y: REFLECTING,
    z: REFLECTING,
};
const SIDES: [BoundarySide; 6] = [
    BoundarySide::XMin,
    BoundarySide::XMax,
    BoundarySide::YMin,
    BoundarySide::YMax,
    BoundarySide::ZMin,
    BoundarySide::ZMax,
];

fn rotate(v: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| (0..3).map(|j| AXES[j][i] * v[j]).sum())
}

fn translated(v: [f64; 3]) -> [f64; 3] {
    let v = rotate(v);
    let shift = [0.25, -0.5, 1.0];
    core::array::from_fn(|i| v[i] + shift[i])
}

#[test]
fn volumetric_shock_footprint_preserves_both_moments_and_affine_work() {
    let gas = gas();
    let grid = CartesianGrid::new(4, 2, 2, 0.5, 0.5, 0.5).unwrap();
    let states = (0..16)
        .map(|i| {
            let y = (i / 4) % 2;
            let z = i / 8;
            State::from_primitive(
                &gas,
                1.0,
                [1.0, 0.0, 0.0],
                1.0 + 2.0 * y as f64 + 3.0 * z as f64,
            )
            .unwrap()
        })
        .collect();
    let mut flow = Euler3d::new(gas, grid, states, CLOSED, || false).unwrap();
    assert!(WallPressureImpulse::from_euler3d(&flow, BoundarySide::XMax, 0, AXES, 0.0).is_err());
    let report = flow.step(Time::new(0.001), controls(), || false).unwrap();
    let exterior_pressure = 0.2;
    let impulses = (0..4)
        .map(|face| {
            WallPressureImpulse::from_euler3d(
                &flow,
                BoundarySide::XMax,
                face,
                AXES,
                exterior_pressure,
            )
            .unwrap()
        })
        .collect::<Vec<_>>();

    let nodes = (0..9)
        .map(|i| translated([2.0, 0.5 * (i % 3) as f64, 0.5 * (i / 3) as f64]))
        .collect();
    let mut triangles = Vec::new();
    let patches = (0..4)
        .map(|face| {
            let node = face % 2 + 3 * (face / 2);
            triangles.push([node, node + 1, node + 4]);
            triangles.push([node, node + 4, node + 3]);
            (2 * face..2 * face + 2)
                .map(|triangle| OrientedShellFacet {
                    triangle,
                    side: ShellFace::Positive,
                })
                .collect()
        })
        .collect();
    let mesh = ShellMesh::new(nodes, triangles).unwrap();
    let budget = SurfaceBudget {
        geometry: ProfileBudget {
            max_nodes: 9,
            max_triangles: 8,
            max_feature_evaluations: 0,
        },
        maximum_patches: 4,
        maximum_preparation_terms: 400,
    };
    let surface = ShockLoadSurface::new(&mesh, patches, budget, || false).unwrap();
    let load = surface
        .project(&impulses, translated([2.0, 0.5, 0.5]), || false)
        .unwrap();
    let dt = report.elapsed.value();
    let raw = flow.boundary_impulses(BoundarySide::XMax);
    assert!(raw[1].momentum[0] > raw[0].momentum[0]);
    assert!(raw[2].momentum[0] > raw[0].momentum[0]);
    let mut force_impulse = 0.0;
    let mut moment = [0.0; 3];
    let mut affine_work = 0.0;
    for (face, raw) in raw.iter().enumerate() {
        let yc = 0.25 + 0.5 * (face % 2) as f64;
        let zc = 0.25 + 0.5 * (face / 2) as f64;
        let impulse = raw.momentum[0] - exterior_pressure * 0.25 * dt;
        force_impulse += impulse;
        moment[1] += (zc - 0.5) * impulse;
        moment[2] -= (yc - 0.5) * impulse;
        affine_work += impulse / dt * (1.0 + 2.0 * yc + 3.0 * zc);
    }
    assert!(moment[1] > 0.0 && moment[2] < 0.0);
    for (actual, expected) in
        load.resultant_impulse_n_s
            .into_iter()
            .zip(rotate([force_impulse, 0.0, 0.0]))
    {
        close(actual, expected);
    }
    for (actual, expected) in load.angular_impulse_n_m_s.into_iter().zip(rotate(moment)) {
        close(actual, expected);
    }
    // P1 integrates an affine displacement exactly on each original triangle.
    let displacement = (0..9)
        .map(|i| rotate([1.0 + (i % 3) as f64 + 1.5 * (i / 3) as f64, 0.0, 0.0]))
        .collect::<Vec<_>>();
    close(load.virtual_work_j(&displacement).unwrap(), affine_work);
}

fn uniform_flow(boundaries: Boundaries3d) -> Euler3d {
    let gas = gas();
    let grid = CartesianGrid::new(2, 3, 4, 0.2, 0.3, 0.4).unwrap();
    let state = State::from_primitive(&gas, 1.0, [0.8, -0.3, 0.6], 1.0).unwrap();
    Euler3d::new(gas, grid, vec![state; 24], boundaries, || false).unwrap()
}

#[test]
fn every_rotated_wall_face_closes_momentum_and_cancelled_steps_retain_the_load() {
    let mut flow = uniform_flow(CLOSED);
    let report = flow.step(Time::new(0.001), controls(), || false).unwrap();
    let mut support = [0.0; 3];
    for side in SIDES {
        for face in 0..flow.grid().face_count(side) {
            let impulse = WallPressureImpulse::from_euler3d(&flow, side, face, AXES, 0.0).unwrap();
            let magnitude = impulse.gas_pressure_impulse_pa_s() * impulse.area_m2();
            for (sum, normal) in support.iter_mut().zip(impulse.outward_normal()) {
                *sum += normal * magnitude;
            }
        }
    }
    let fluid = rotate(core::array::from_fn(|axis| {
        report.after.momentum[axis] - report.before.momentum[axis]
    }));
    for axis in 0..3 {
        close(fluid[axis] + support[axis], 0.0);
    }

    let side = BoundarySide::ZMax;
    let accepted = WallPressureImpulse::from_euler3d(&flow, side, 0, AXES, 0.0).unwrap();
    let states = flow.states().to_vec();
    let time = flow.time();
    assert!(
        WallPressureImpulse::from_euler3d(&flow, side, flow.grid().face_count(side), AXES, 0.0,)
            .is_err()
    );
    assert!(WallPressureImpulse::from_euler3d(&flow, side, 0, [AXES[0]; 3], 0.0).is_err());
    assert!(WallPressureImpulse::from_euler3d(&flow, side, 0, AXES, -1.0).is_err());
    assert!(flow.step(Time::new(0.001), controls(), || true).is_err());
    assert_eq!(flow.states(), states);
    assert_eq!(flow.time(), time);
    assert_eq!(
        WallPressureImpulse::from_euler3d(&flow, side, 0, AXES, 0.0).unwrap(),
        accepted
    );
}

#[test]
fn transported_or_periodic_momentum_cannot_be_used_as_volumetric_wall_pressure() {
    for pair in [
        AxisBoundaries::Periodic,
        AxisBoundaries::Ends {
            lower: BoundaryCondition::Transmissive,
            upper: BoundaryCondition::Transmissive,
        },
    ] {
        let mut flow = uniform_flow(Boundaries3d { x: pair, ..CLOSED });
        flow.step(Time::new(0.001), controls(), || false).unwrap();
        for side in [BoundarySide::XMin, BoundarySide::XMax] {
            assert!(WallPressureImpulse::from_euler3d(&flow, side, 0, AXES, 0.0).is_err());
        }
        assert!(WallPressureImpulse::from_euler3d(&flow, BoundarySide::YMax, 0, AXES, 0.0).is_ok());
    }
}
