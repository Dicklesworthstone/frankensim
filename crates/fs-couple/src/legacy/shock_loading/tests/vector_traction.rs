//! Full resolved wall momentum reaches structural force, torque and P1 work.
use super::molecular::law;
use super::*;
use fs_flux::compressible::spatial::{
    Boundaries3d, BoundarySide, CartesianGrid, Euler3d, State,
    molecular::{MolecularEuler3d, MolecularWall3d, WallThermalCondition},
};
use fs_flux::compressible::{
    AxisBoundaries, Boundaries2d, BoundaryCondition, BoundarySide as PlanarSide,
    CartesianGrid as PlanarGrid, Euler2d, State as PlanarState,
    molecular::{MolecularEuler2d, MolecularStepControl, MolecularWall2d},
};
use fs_gas::multicomponent::dissipation::ConstantTransport;

const AXES: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const ROTATED: [[f64; 3]; 3] = [[0.6, 0.8, 0.0], [-0.8, 0.6, 0.0], [0.0, 0.0, 1.0]];
const STATIONARY: MolecularWall3d = MolecularWall3d::NoSlip {
    velocity: [0.0; 3],
    thermal: WallThermalCondition::Adiabatic,
};

fn channel(velocity: [f64; 3], wall: MolecularWall3d) -> MolecularEuler3d<ConstantTransport> {
    let gas = gas();
    let transport = law(&gas, 0.04);
    let state = State::from_primitive(&gas, 1.0, velocity, 1.0).unwrap();
    let flow = Euler3d::new(
        gas,
        CartesianGrid::new(4, 1, 1, 0.5, 1.0, 1.0).unwrap(),
        vec![state; 4],
        Boundaries3d {
            x: AxisBoundaries::Ends {
                lower: BoundaryCondition::Reflecting,
                upper: BoundaryCondition::Reflecting,
            },
            y: AxisBoundaries::Periodic,
            z: AxisBoundaries::Periodic,
        },
        || false,
    )
    .unwrap();
    MolecularEuler3d::new(flow, transport, || false)
        .unwrap()
        .with_wall(BoundarySide::XMin, wall)
        .unwrap()
        .with_wall(BoundarySide::XMax, wall)
        .unwrap()
}

fn planar_channel(
    velocity: [f64; 2],
    wall: MolecularWall2d,
) -> MolecularEuler2d<ConstantTransport> {
    let gas = gas();
    let transport = law(&gas, 0.04);
    let state = PlanarState::from_primitive(&gas, 1.0, velocity, 1.0).unwrap();
    let flow = Euler2d::new(
        gas,
        PlanarGrid::new(4, 1, 0.5, 1.0, 2.0).unwrap(),
        vec![state; 4],
        Boundaries2d {
            x: AxisBoundaries::Ends {
                lower: BoundaryCondition::Reflecting,
                upper: BoundaryCondition::Reflecting,
            },
            y: AxisBoundaries::Periodic,
        },
        || false,
    )
    .unwrap();
    MolecularEuler2d::new(flow, transport, || false)
        .unwrap()
        .with_wall(PlanarSide::XMin, wall)
        .unwrap()
        .with_wall(PlanarSide::XMax, wall)
        .unwrap()
}

fn rotate(v: [f64; 3]) -> [f64; 3] {
    [0.6 * v[0] - 0.8 * v[1], 0.8 * v[0] + 0.6 * v[1], v[2]]
}

fn affine_displacements(mesh: &ShellMesh) -> Vec<[f64; 3]> {
    mesh.nodes
        .iter()
        .map(|p| {
            [
                0.5 + 2.0 * p[1] - p[2],
                -1.0 + p[0] + 3.0 * p[2],
                2.0 - p[0] + p[1],
            ]
        })
        .collect()
}

#[test]
fn vector_p1_projection_preserves_independent_force_moment_and_affine_work() {
    let impulse =
        WallVectorTractionImpulse::new(0.2, 1.4, 1.0, [1.0, 0.0, 0.0], 7.0, [0.0, -3.0, 5.0], 2.0)
            .unwrap();
    let load = surface(ShellFace::Positive)
        .project_vector_tractions(&[impulse], [0.25, -0.5, 1.25], || false)
        .unwrap();
    // Unit square centroid (2, 1/2, 1/2), with exterior impulse (-0.4, 0, 0).
    // The moment is independently r x J with r = (1.75, 1, -0.75).
    for axis in 0..3 {
        close(load.gas_impulse_n_s[axis], [7.0, -3.0, 5.0][axis]);
        close(load.exterior_impulse_n_s[axis], [-0.4, 0.0, 0.0][axis]);
        close(load.resultant_impulse_n_s[axis], [6.6, -3.0, 5.0][axis]);
        close(
            load.angular_impulse_n_m_s[axis],
            [2.75, -13.7, -11.85][axis],
        );
    }
    // The exact area mean of this affine field is (1, 2.5, 0.5), giving 8 J.
    close(
        load.virtual_work_j(&affine_displacements(&mesh())).unwrap(),
        8.0,
    );
    let centroid = surface(ShellFace::Positive)
        .project_vector_tractions(&[impulse], [2.0, 0.5, 0.5], || false)
        .unwrap();
    for moment in centroid.angular_impulse_n_m_s {
        close(moment, 0.0);
    }
    assert_eq!(centroid.nodal_impulses_n_s, load.nodal_impulses_n_s);
}

#[test]
fn nonuniform_shear_retains_a_pure_couple_and_refuses_mixed_intervals() {
    let surface = ShockLoadSurface::new(
        &mesh(),
        vec![
            vec![OrientedShellFacet {
                triangle: 0,
                side: ShellFace::Positive,
            }],
            vec![OrientedShellFacet {
                triangle: 1,
                side: ShellFace::Positive,
            }],
        ],
        budget(),
        || false,
    )
    .unwrap();
    let impulses = [1.0, -1.0].map(|shear| {
        WallVectorTractionImpulse::new(1.0, 2.0, 0.5, AXES[0], 0.0, [0.0, shear, 0.0], 0.0).unwrap()
    });
    let load = surface
        .project_vector_tractions(&impulses, [0.0; 3], || false)
        .unwrap();
    for axis in 0..3 {
        close(load.resultant_impulse_n_s[axis], 0.0);
        close(
            load.angular_impulse_n_m_s[axis],
            [1.0 / 3.0, 0.0, 0.0][axis],
        );
    }
    // Equal opposite impulses act at z=1/3 and z=2/3. Their exact torque
    // does work 1/3 J on a unit infinitesimal rotation about the x axis.
    let rotation: Vec<_> = mesh().nodes.iter().map(|p| [0.0, -p[2], p[1]]).collect();
    close(load.virtual_work_j(&rotation).unwrap(), 1.0 / 3.0);
    let later =
        WallVectorTractionImpulse::new(1.0, 3.0, 0.5, AXES[0], 0.0, [0.0, -1.0, 0.0], 0.0).unwrap();
    assert!(
        surface
            .project_vector_tractions(&[impulses[0], later], [0.0; 3], || false)
            .is_err()
    );
    assert!(
        surface
            .project_vector_tractions(&impulses[..1], [0.0; 3], || false)
            .is_err()
    );
}

#[test]
fn balanced_large_normal_pressure_does_not_erase_small_rotated_wall_shear() {
    let rotated_mesh =
        ShellMesh::new(mesh().nodes.into_iter().map(rotate).collect(), mesh().tris).unwrap();
    let impulse =
        WallVectorTractionImpulse::new(1.0, 1.0, 1.0, ROTATED[0], 1e20, ROTATED[1], 1e20).unwrap();
    // Combining the two world-frame pressure vectors before subtraction
    // would lose this entire 1 N s shear impulse beside the 1e20 N s normal.
    assert_eq!(impulse.net_traction_impulse_pa_s(), ROTATED[1]);
    let load = ShockLoadSurface::new(
        &rotated_mesh,
        vec![selection(ShellFace::Positive)],
        budget(),
        || false,
    )
    .unwrap()
    .project_vector_tractions(&[impulse], [0.0; 3], || false)
    .unwrap();
    for axis in 0..3 {
        close(load.resultant_impulse_n_s[axis], ROTATED[1][axis]);
    }
}

#[test]
fn opposing_face_loads_preserve_small_shared_node_impulse_force_moment_and_work() {
    // The first and last triangles share the edge (0,1); their centroids
    // differ only along the direction of their equal opposite z impulses.
    // The middle triangle contributes a small impulse at their shared node.
    let mesh = ShellMesh::new(
        vec![
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [2.0, 0.0, 1.0],
            [2.0, 0.0, -1.0],
            [2.0, -1.0, 0.0],
            [2.0, -1.0, -1.0],
        ],
        vec![[0, 1, 2], [0, 4, 5], [0, 3, 1]],
    )
    .unwrap();
    let surface = ShockLoadSurface::new(
        &mesh,
        (0..3)
            .map(|triangle| {
                vec![OrientedShellFacet {
                    triangle,
                    side: ShellFace::Positive,
                }]
            })
            .collect(),
        SurfaceBudget {
            geometry: ProfileBudget {
                max_nodes: 6,
                max_triangles: 3,
                max_feature_evaluations: 0,
            },
            maximum_patches: 3,
            maximum_preparation_terms: 200,
        },
        || false,
    )
    .unwrap();
    let impulses = [1e16, 1.0, -1e16].map(|value| {
        WallVectorTractionImpulse::new(1.0, 1.0, 0.5, AXES[0], 0.0, [0.0, 0.0, value], 0.0).unwrap()
    });
    let load = surface
        .project_vector_tractions(&impulses, [0.0; 3], || false)
        .unwrap();
    close(load.nodal_impulses_n_s[0][2], 1.0 / 3.0);
    assert_eq!(load.nodal_impulses_n_s[1], [0.0; 3]);
    // Only the middle triangle's 1 N s impulse remains: its centroid is
    // (2, -2/3, -1/3), yielding torque (-2/3, -2, 0) N m s.
    for axis in 0..3 {
        close(load.gas_impulse_n_s[axis], [0.0, 0.0, 1.0][axis]);
        close(load.resultant_impulse_n_s[axis], [0.0, 0.0, 1.0][axis]);
        close(
            load.angular_impulse_n_m_s[axis],
            [-2.0 / 3.0, -2.0, 0.0][axis],
        );
    }
    close(
        load.virtual_work_j(&vec![[0.0, 0.0, 1.0]; mesh.nodes.len()])
            .unwrap(),
        1.0,
    );
}

#[test]
fn accepted_no_slip_shear_closes_gas_support_momentum_and_rotated_structural_work() {
    let mut flow = channel([0.0, 0.3, -0.2], STATIONARY);
    assert!(
        WallVectorTractionImpulse::from_molecular_euler3d(&flow, BoundarySide::XMax, 0, AXES, 0.0)
            .is_err()
    );
    let report = flow
        .step(Time::new(1e-6), MolecularStepControl::default(), || false)
        .unwrap();
    let right =
        WallVectorTractionImpulse::from_molecular_euler3d(&flow, BoundarySide::XMax, 0, AXES, 0.25)
            .unwrap();
    let left =
        WallVectorTractionImpulse::from_molecular_euler3d(&flow, BoundarySide::XMin, 0, AXES, 0.0)
            .unwrap();
    assert_eq!(right.duration_s(), report.elapsed.value());
    assert_eq!(right.end_time_s(), flow.flow().time().value());
    assert_eq!(right.area_m2(), 1.0);
    assert_eq!(right.outward_normal(), AXES[0]);
    close(
        right.exterior_pressure_impulse_pa_s(),
        0.25 * report.elapsed.value(),
    );
    // Independent half-cell Newtonian shear: mu U/(h/2). The small complete
    // interval approaches its initial value, while retaining both RK halves.
    for impulse in [left, right] {
        let mean = impulse.mean_gas_traction_pa();
        for axis in [1, 2] {
            let expected = 0.04 * [0.0, 0.3, -0.2][axis] / 0.25;
            assert!((mean[axis] - expected).abs() < 1e-6 * expected.abs());
        }
    }
    let raw = flow.boundary_impulses(BoundarySide::XMax)[0];
    assert_eq!(raw.mass, 0.0);
    assert_eq!(raw.energy, 0.0);
    assert_eq!(
        flow.flow().boundary_impulses(BoundarySide::XMax)[0].momentum[1],
        0.0
    );
    assert!(raw.momentum[1] > 0.0 && raw.momentum[2] < 0.0);
    // The earlier normal-only API explicitly refuses this real shear impulse.
    assert!(
        WallTractionImpulse::from_molecular_euler3d(&flow, BoundarySide::XMax, 0, AXES, 0.25)
            .is_err()
    );

    let load = surface(ShellFace::Positive)
        .project_vector_tractions(&[right], [0.0; 3], || false)
        .unwrap();
    let mut left_mesh = mesh();
    for node in &mut left_mesh.nodes {
        node[0] = 0.0;
    }
    let left_load = ShockLoadSurface::new(
        &left_mesh,
        vec![selection(ShellFace::Negative)],
        budget(),
        || false,
    )
    .unwrap()
    .project_vector_tractions(&[left], [0.0; 3], || false)
    .unwrap();
    let expected = [
        raw.momentum[0] - 0.25 * report.elapsed.value(),
        raw.momentum[1],
        raw.momentum[2],
    ];
    let moment = [
        0.5 * (expected[2] - expected[1]),
        0.5 * expected[0] - 2.0 * expected[2],
        2.0 * expected[1] - 0.5 * expected[0],
    ];
    for axis in 0..3 {
        close(load.gas_impulse_n_s[axis], raw.momentum[axis]);
        close(load.resultant_impulse_n_s[axis], expected[axis]);
        close(load.angular_impulse_n_m_s[axis], moment[axis]);
        close(
            report.after.momentum[axis] - report.before.momentum[axis]
                + load.gas_impulse_n_s[axis]
                + left_load.gas_impulse_n_s[axis],
            0.0,
        );
    }
    let displacement = affine_displacements(&mesh());
    close(
        load.virtual_work_j(&displacement).unwrap(),
        (expected[0] + 2.5 * expected[1] + 0.5 * expected[2]) / report.elapsed.value(),
    );

    let rotated_mesh =
        ShellMesh::new(mesh().nodes.into_iter().map(rotate).collect(), mesh().tris).unwrap();
    let rotated_impulse = WallVectorTractionImpulse::from_molecular_euler3d(
        &flow,
        BoundarySide::XMax,
        0,
        ROTATED,
        0.25,
    )
    .unwrap();
    let rotated = ShockLoadSurface::new(
        &rotated_mesh,
        vec![selection(ShellFace::Positive)],
        budget(),
        || false,
    )
    .unwrap()
    .project_vector_tractions(&[rotated_impulse], [0.0; 3], || false)
    .unwrap();
    for axis in 0..3 {
        close(rotated.resultant_impulse_n_s[axis], rotate(expected)[axis]);
        close(rotated.angular_impulse_n_m_s[axis], rotate(moment)[axis]);
    }
    close(
        rotated
            .virtual_work_j(&displacement.into_iter().map(rotate).collect::<Vec<_>>())
            .unwrap(),
        load.virtual_work_j(&affine_displacements(&mesh())).unwrap(),
    );
    assert!(
        WallVectorTractionImpulse::from_molecular_euler3d(&flow, BoundarySide::YMin, 0, AXES, 0.0)
            .is_err()
    );
    assert!(
        WallVectorTractionImpulse::from_molecular_euler3d(&flow, BoundarySide::XMax, 1, AXES, 0.0)
            .is_err()
    );
    assert!(
        WallVectorTractionImpulse::from_molecular_euler3d(
            &flow,
            BoundarySide::XMax,
            0,
            [AXES[0]; 3],
            0.0
        )
        .is_err()
    );
}

#[test]
fn moving_and_isothermal_wall_metadata_cannot_hide_behind_zero_energy_transfer() {
    let source = gas();
    let temperature = State::from_primitive(&source, 1.0, [0.0; 3], 1.0)
        .unwrap()
        .primitive(&source, || false)
        .unwrap()
        .temperature;
    for wall in [
        MolecularWall3d::NoSlip {
            velocity: [0.0, 0.3, -0.2],
            thermal: WallThermalCondition::Adiabatic,
        },
        MolecularWall3d::NoSlip {
            velocity: [0.0; 3],
            thermal: WallThermalCondition::Isothermal(Temperature::new(temperature)),
        },
    ] {
        let MolecularWall3d::NoSlip { velocity, .. } = wall else {
            unreachable!()
        };
        let mut flow = channel(velocity, wall);
        flow.step(Time::new(1e-6), MolecularStepControl::default(), || false)
            .unwrap();
        assert_eq!(flow.boundary_impulses(BoundarySide::XMax)[0].energy, 0.0);
        assert!(
            WallVectorTractionImpulse::from_molecular_euler3d(
                &flow,
                BoundarySide::XMax,
                0,
                AXES,
                0.0
            )
            .is_err()
        );
        // Reauthoring a supported wall cannot relabel the previously accepted
        // moving/isothermal impulse as a stationary adiabatic interval.
        let mut flow = flow.with_wall(BoundarySide::XMax, STATIONARY).unwrap();
        assert!(flow.last_step().is_none());
        assert!(
            WallVectorTractionImpulse::from_molecular_euler3d(
                &flow,
                BoundarySide::XMax,
                0,
                AXES,
                0.0
            )
            .is_err()
        );
        flow.step(Time::new(1e-6), MolecularStepControl::default(), || false)
            .unwrap();
        assert!(
            WallVectorTractionImpulse::from_molecular_euler3d(
                &flow,
                BoundarySide::XMax,
                0,
                AXES,
                0.0
            )
            .is_ok()
        );
    }
}

#[test]
fn identical_wall_configuration_and_cancelled_source_steps_retain_accepted_vector_loads() {
    let mut flow = channel([0.0, 0.3, -0.2], STATIONARY);
    let accepted = flow
        .step(Time::new(1e-6), MolecularStepControl::default(), || false)
        .unwrap();
    let before =
        WallVectorTractionImpulse::from_molecular_euler3d(&flow, BoundarySide::XMax, 0, AXES, 0.0)
            .unwrap();
    let mut flow = flow.with_wall(BoundarySide::XMax, STATIONARY).unwrap();
    assert_eq!(flow.last_step(), Some(&accepted));
    assert!(
        flow.step(Time::new(1e-6), MolecularStepControl::default(), || true)
            .is_err()
    );
    assert_eq!(
        WallVectorTractionImpulse::from_molecular_euler3d(&flow, BoundarySide::XMax, 0, AXES, 0.0)
            .unwrap(),
        before,
    );
    assert_eq!(
        surface(ShellFace::Positive)
            .project_vector_tractions(&[before], [0.0; 3], || true)
            .unwrap_err(),
        Error::Cancelled,
    );
}

#[test]
fn planar_no_slip_wall_drag_reaches_spatial_force_torque_and_work_with_physical_depth() {
    let mut flow = planar_channel(
        [0.0, 0.3],
        MolecularWall2d::NoSlip {
            velocity: [0.0; 2],
            thermal: WallThermalCondition::Adiabatic,
        },
    );
    let report = flow
        .step(Time::new(1e-6), MolecularStepControl::default(), || false)
        .unwrap();
    let axes = [AXES[0], AXES[1]];
    let right =
        WallVectorTractionImpulse::from_molecular_euler2d(&flow, PlanarSide::XMax, 0, axes, 0.25)
            .unwrap();
    let left =
        WallVectorTractionImpulse::from_molecular_euler2d(&flow, PlanarSide::XMin, 0, axes, 0.0)
            .unwrap();
    assert_eq!(right.area_m2(), 2.0);
    // Half-cell Newtonian shear is independent of the represented depth;
    // the extensive support momentum below includes that depth exactly once.
    let shear = 0.04 * 0.3 / 0.25;
    assert!((right.mean_gas_traction_pa()[1] - shear).abs() < 1e-6 * shear);
    assert!(
        WallTractionImpulse::from_molecular_euler2d(&flow, PlanarSide::XMax, 0, axes, 0.25)
            .is_err()
    );
    for axis in 0..2 {
        close(
            report.after.momentum[axis] - report.before.momentum[axis]
                + 2.0 * right.gas_traction_impulse_pa_s()[axis]
                + 2.0 * left.gas_traction_impulse_pa_s()[axis],
            0.0,
        );
    }

    let mut structural = mesh();
    for node in &mut structural.nodes {
        node[2] *= 2.0;
    }
    let raw = flow.boundary_impulses(PlanarSide::XMax)[0];
    let expected = [
        raw.momentum[0] - 0.25 * 2.0 * report.elapsed.value(),
        raw.momentum[1],
        0.0,
    ];
    // The physical patch centroid is (2, 1/2, 1), including its depth.
    let torque = [
        -expected[1],
        expected[0],
        2.0 * expected[1] - 0.5 * expected[0],
    ];
    let displacement = affine_displacements(&structural);
    for rotated_frame in [false, true] {
        let (positions, axes, expected, torque, displacement) = if rotated_frame {
            (
                structural.nodes.iter().copied().map(rotate).collect(),
                [ROTATED[0], ROTATED[1]],
                rotate(expected),
                rotate(torque),
                displacement.iter().copied().map(rotate).collect::<Vec<_>>(),
            )
        } else {
            (
                structural.nodes.clone(),
                axes,
                expected,
                torque,
                displacement.clone(),
            )
        };
        let mesh = ShellMesh::new(positions, structural.tris.clone()).unwrap();
        let surface = ShockLoadSurface::new(
            &mesh,
            vec![selection(ShellFace::Positive)],
            budget(),
            || false,
        )
        .unwrap();
        let impulse = WallVectorTractionImpulse::from_molecular_euler2d(
            &flow,
            PlanarSide::XMax,
            0,
            axes,
            0.25,
        )
        .unwrap();
        let load = surface
            .project_vector_tractions(&[impulse], [0.0; 3], || false)
            .unwrap();
        for axis in 0..3 {
            close(load.resultant_impulse_n_s[axis], expected[axis]);
            close(load.angular_impulse_n_m_s[axis], torque[axis]);
        }
        // Exact mean of the original affine field is (1/2, 4, 1/2).
        close(
            load.virtual_work_j(&displacement).unwrap(),
            (0.5 * (raw.momentum[0] - 0.5 * report.elapsed.value()) + 4.0 * raw.momentum[1])
                / report.elapsed.value(),
        );
    }
}

#[test]
fn planar_vector_loads_refuse_moving_or_isothermal_metadata_even_without_energy_transfer() {
    let source = gas();
    let temperature = PlanarState::from_primitive(&source, 1.0, [0.0; 2], 1.0)
        .unwrap()
        .primitive(&source, || false)
        .unwrap()
        .temperature;
    for wall in [
        MolecularWall2d::NoSlip {
            velocity: [0.0, 0.3],
            thermal: WallThermalCondition::Adiabatic,
        },
        MolecularWall2d::NoSlip {
            velocity: [0.0; 2],
            thermal: WallThermalCondition::Isothermal(Temperature::new(temperature)),
        },
    ] {
        let MolecularWall2d::NoSlip { velocity, .. } = wall else {
            unreachable!()
        };
        let mut flow = planar_channel(velocity, wall);
        flow.step(Time::new(1e-6), MolecularStepControl::default(), || false)
            .unwrap();
        assert_eq!(flow.boundary_impulses(PlanarSide::XMax)[0].energy, 0.0);
        assert!(
            WallVectorTractionImpulse::from_molecular_euler2d(
                &flow,
                PlanarSide::XMax,
                0,
                [AXES[0], AXES[1]],
                0.0,
            )
            .is_err()
        );
    }
}
