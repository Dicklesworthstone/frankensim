use super::*;
use crate::d3q19::equilibrium3;
use crate::d3q19::freesurface3::{
    ContactModel3, PlanarWall3, PlanarWallFieldControls3, PlanarWallVelocityField3,
};
use crate::kinematics::KinematicsError;
use crate::sgs3::SgsUse;

const WALLS: [PlanarWall3; 6] = [
    PlanarWall3::XMin,
    PlanarWall3::XMax,
    PlanarWall3::YMin,
    PlanarWall3::YMax,
    PlanarWall3::ZMin,
    PlanarWall3::ZMax,
];

fn controls() -> PlanarWallFieldControls3 {
    PlanarWallFieldControls3 {
        maximum_domain_cells: 512,
        maximum_wall_faces: 64,
        maximum_relative_mass_defect: 1e-10,
        maximum_relative_face_mass_defect: 1e-12,
    }
}

fn config(model: SgsModel) -> SgsConfig {
    SgsConfig {
        model,
        filter_width_lu: 1.3,
        maximum_eddy_to_molecular: 10.0,
        usage: SgsUse::Exploration,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

// Linear/quadratic wall-normal shear plus rigid rotation about the normal.
// The normal velocity is identically zero; the wall motion is tangent-only.
fn profile(
    wall: PlanarWall3,
    p: [usize; 3],
    linear: f64,
    quadratic: f64,
    spin: f64,
) -> ([f64; 3], [[f64; 3]; 3]) {
    let n = wall.axis();
    let a = (n + 1) % 3;
    let b = (n + 2) % 3;
    let sign = if wall.upper() { -1.0 } else { 1.0 };
    let s = if wall.upper() {
        7.5 - p[n] as f64
    } else {
        p[n] as f64 + 0.5
    };
    let mut u = [0.0; 3];
    u[a] = linear * s + quadratic * s * s - spin * (p[b] as f64 - 3.5);
    u[b] = spin * (p[a] as f64 - 3.5);
    let mut gradient = [[0.0; 3]; 3];
    gradient[a][n] = sign * (linear + 2.0 * quadratic * s);
    gradient[a][b] = -spin;
    gradient[b][a] = spin;
    (u, gradient)
}

fn fixture(
    wall: PlanarWall3,
    layers: usize,
    linear: f64,
    quadratic: f64,
    spin: f64,
    boost: [f64; 2],
) -> (FreeSurface3, Vec<[f64; 2]>) {
    let gate = CancelGate::new_clock_free();
    let n = wall.axis();
    let a = (n + 1) % 3;
    let b = (n + 2) % 3;
    let mut liquid = FreeSurface3::new(
        8,
        8,
        8,
        0.8,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        |x, y, z| {
            let normal = [x, y, z][n];
            if wall.upper() {
                normal >= 8 - layers
            } else {
                normal < layers
            }
        },
    )
    .unwrap();
    let mut forces = BTreeMap::new();
    for key in liquid.grid.active_keys().to_vec() {
        let slot = liquid.grid.slot_of(key).unwrap();
        for lane in 0..TILE_CELLS {
            if liquid.tiles[&key].cells[lane] == Cell3::Gas {
                continue;
            }
            let (x, y, z) = FreeSurface3::global_of(key, lane);
            let p = [x as usize, y as usize, z as usize];
            let (mut velocity, g) = profile(wall, p, linear, quadratic, spin);
            velocity[a] += boost[0];
            velocity[b] += boost[1];
            let mut populations = equilibrium3(1.0, velocity);
            crate::d3q19::freesurface3::forcing::preserve_density(&mut populations, 1.0).unwrap();
            liquid.grid.set_populations(slot, lane, populations);
            // A physical incoming viscous stress makes changed relaxation
            // observable in the first actual collision, not just its report.
            let strain = [
                g[0][0],
                g[1][1],
                g[2][2],
                0.5 * (g[0][1] + g[1][0]),
                0.5 * (g[0][2] + g[2][0]),
                0.5 * (g[1][2] + g[2][1]),
            ];
            liquid.initialize_strain_at(p, strain, 0.8).unwrap();
            forces.insert(
                (key, lane),
                [
                    0.0001 * (p[1] as f64 + 0.5),
                    -0.00007 * (p[2] as f64 + 0.5),
                    0.00004 * (p[0] as f64 + 0.5),
                ],
            );
        }
    }
    // This spatial Guo reframe must not alter any physical gradient.
    liquid.rebind_force_field(forces, &gate).unwrap();
    let wall_velocity = (0..64)
        .map(|i| {
            [
                boost[0] - spin * ((i / 8) as f64 - 3.5),
                boost[1] + spin * ((i % 8) as f64 - 3.5),
            ]
        })
        .collect();
    (liquid, wall_velocity)
}

fn motion(wall: PlanarWall3, velocity: &[[f64; 2]]) -> PlanarWallVelocityField3<'_> {
    PlanarWallVelocityField3 {
        wall,
        moment_origin: [0.0; 3],
        tangential_velocity: velocity,
    }
}

#[test]
fn g1_wall_distance_recovers_shear_rotation_and_two_layer_interface_gradients() {
    let gate = CancelGate::new_clock_free();
    for wall in WALLS {
        for layers in [2, 5] {
            for (linear, quadratic, spin) in [
                (0.002, 0.0, 0.0),
                (0.003, -0.00015, 0.0),
                (0.0, 0.0, 0.0015),
                (0.002, 0.0001, 0.0015),
            ] {
                for boost in [[0.0; 2], [0.007, -0.004]] {
                    let (liquid, walls) = fixture(wall, layers, linear, quadratic, spin, boost);
                    let mut workspace = FreeSurfaceKinematicsWorkspace3::new(
                        liquid.grid.active_tiles(),
                        usize::MAX,
                    )
                    .unwrap();
                    let field = liquid
                        .kinematics_with_wall(&gate, &mut workspace, Some(motion(wall, &walls)))
                        .unwrap();
                    for &key in liquid.grid.active_keys() {
                        for lane in 0..TILE_CELLS {
                            if liquid.tiles[&key].cells[lane] == Cell3::Gas {
                                continue;
                            }
                            let (x, y, z) = FreeSurface3::global_of(key, lane);
                            let p = [x as usize, y as usize, z as usize];
                            let sample = field.sample_at(p[0], p[1], p[2]).unwrap().unwrap();
                            let expected = profile(wall, p, linear, quadratic, spin).1;
                            for (actual, expected) in sample
                                .kinematics
                                .gradient()
                                .iter()
                                .flatten()
                                .zip(expected.iter().flatten())
                            {
                                near(*actual, *expected, 3e-14);
                            }
                            near(
                                sample.kinematics.vorticity()[wall.axis()],
                                2.0 * spin,
                                4e-14,
                            );
                            if layers == 2 {
                                assert_eq!(
                                    sample.stencils[wall.axis()],
                                    if wall.upper() {
                                        WetGradientStencil3::UpperWallSecondOrder
                                    } else {
                                        WetGradientStencil3::LowerWallSecondOrder
                                    }
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn g1_wall_gradient_reaches_native_sgs_collision_and_preserves_bulk_and_trt() {
    let gate = CancelGate::new_clock_free();
    for wall in WALLS {
        for family in 0..3 {
            let (mut liquid, walls) = fixture(wall, 2, 0.002, 0.0001, 0.0015, [0.0; 2]);
            if family == 1 {
                liquid.set_trt_magic(0.25).unwrap();
            }
            if family == 2 {
                liquid.set_central_moment_bulk(0.91, 1.1).unwrap();
                let bulk = liquid
                    .grid
                    .active_keys()
                    .iter()
                    .flat_map(|&key| (0..TILE_CELLS).map(move |lane| (key, lane)))
                    .filter(|&(key, lane)| liquid.tiles[&key].cells[lane] != Cell3::Gas)
                    .map(|key| (key, 0.85 + 0.01 * (key.1 % 4) as f64))
                    .collect();
                liquid.rebind_bulk_relaxations(bulk, &gate).unwrap();
            }
            let mut kinematics =
                FreeSurfaceKinematicsWorkspace3::new(liquid.grid.active_tiles(), usize::MAX)
                    .unwrap();
            let mut workspace = FreeSurfaceSgsWorkspace3::default();
            // The existing liquid-only route honestly cannot differentiate
            // this two-layer film. Its wall-driven route now has enough data.
            assert!(
                liquid
                    .prepare_gradient_sgs(
                        &gate,
                        config(SgsModel::Wale { cw: 0.4 }),
                        &mut kinematics,
                        &mut workspace
                    )
                    .is_err()
            );
            for model in [SgsModel::Wale { cw: 0.4 }, SgsModel::Vreman { cv: 0.1 }] {
                let controls = controls();
                let config = config(model);
                let motion = motion(wall, &walls);
                let mut exact_tau = BTreeMap::new();
                for &key in liquid.grid.active_keys() {
                    for lane in 0..TILE_CELLS {
                        if liquid.tiles[&key].cells[lane] == Cell3::Gas {
                            continue;
                        }
                        let (x, y, z) = FreeSurface3::global_of(key, lane);
                        let g = profile(
                            wall,
                            [x as usize, y as usize, z as usize],
                            0.002,
                            0.0001,
                            0.0015,
                        )
                        .1;
                        let eddy = eddy_viscosity(model, config.filter_width_lu, g).unwrap();
                        exact_tau.insert((key, lane), 0.8 + 3.0 * eddy);
                    }
                }
                // Independent analytic-gradient relaxation enters the same
                // original collision, reflection, VOF and conversion operator.
                let (reference, reference_wall) = liquid
                    .planar_wall_field_candidate_with_relaxation(
                        motion,
                        controls,
                        Some(&exact_tau),
                        &gate,
                    )
                    .unwrap();
                let (actual, actual_wall, report) = liquid
                    .planar_wall_velocity_field_gradient_sgs_candidate(
                        motion,
                        controls,
                        &gate,
                        config,
                        &mut kinematics,
                        &mut workspace,
                    )
                    .unwrap();
                assert!(report.maximum_eddy_viscosity_lu > 0.0);
                assert!(report.sgs_dissipation_lu > 0.0);
                assert_eq!(actual.collision_model(), liquid.collision_model());
                assert_eq!(actual.local_bulk_relaxation, liquid.local_bulk_relaxation);
                assert_eq!(actual.grid.active_keys(), reference.grid.active_keys());
                for (slot, &key) in actual.grid.active_keys().iter().enumerate() {
                    for lane in 0..TILE_CELLS {
                        for (a, b) in actual
                            .grid
                            .populations(slot, lane)
                            .iter()
                            .zip(reference.grid.populations(slot, lane))
                        {
                            near(*a, b, 8e-14);
                        }
                        if let Some(expected) = exact_tau.get(&(key, lane)) {
                            near(workspace.tau[&(key, lane)], *expected, 8e-14);
                        }
                    }
                }
                for axis in 0..3 {
                    near(
                        actual_wall.fluid_impulse[axis],
                        reference_wall.fluid_impulse[axis],
                        2e-11,
                    );
                    near(
                        actual_wall.fluid_angular_impulse[axis],
                        reference_wall.fluid_angular_impulse[axis],
                        2e-11,
                    );
                }
                near(actual_wall.fluid_work, reference_wall.fluid_work, 2e-12);
                near(actual.ledger_mass(), liquid.ledger_mass(), 2e-10);
                let (native, _) = liquid
                    .planar_wall_velocity_field_candidate(motion, controls, &gate)
                    .unwrap();
                assert_ne!(actual.grid.state_bits(), native.grid.state_bits());
                let disabled = SgsConfig {
                    model: match model {
                        SgsModel::Wale { .. } => SgsModel::Wale { cw: 0.0 },
                        _ => SgsModel::Vreman { cv: 0.0 },
                    },
                    ..config
                };
                let (zero, zero_wall, _) = liquid
                    .planar_wall_velocity_field_gradient_sgs_candidate(
                        motion,
                        controls,
                        &gate,
                        disabled,
                        &mut kinematics,
                        &mut workspace,
                    )
                    .unwrap();
                let (native, native_wall) = liquid
                    .planar_wall_velocity_field_candidate(motion, controls, &gate)
                    .unwrap();
                assert_eq!(zero.grid.state_bits(), native.grid.state_bits());
                assert_eq!(zero.mass_snapshot(), native.mass_snapshot());
                assert_eq!(zero_wall, native_wall);
            }
        }
    }
}

#[test]
fn g0_rigid_spin_and_laminar_shear_do_not_manufacture_modeled_shear_power() {
    let gate = CancelGate::new_clock_free();
    for (linear, quadratic, spin) in [
        (0.002, 0.0, 0.0),
        (0.003, -0.00015, 0.0),
        (0.0, 0.0, 0.0015),
    ] {
        let (liquid, walls) = fixture(
            PlanarWall3::YMin,
            2,
            linear,
            quadratic,
            spin,
            [0.006, -0.003],
        );
        let mut kinematics =
            FreeSurfaceKinematicsWorkspace3::new(liquid.grid.active_tiles(), usize::MAX).unwrap();
        let mut workspace = FreeSurfaceSgsWorkspace3::default();
        for model in [SgsModel::Wale { cw: 0.4 }, SgsModel::Vreman { cv: 0.1 }] {
            let report = liquid
                .prepare_gradient_sgs_with_wall(
                    &gate,
                    config(model),
                    Some(motion(PlanarWall3::YMin, &walls)),
                    &mut kinematics,
                    &mut workspace,
                )
                .unwrap();
            assert!(report.sgs_dissipation_lu < 1e-18);
            if spin == 0.0 {
                assert!(report.maximum_eddy_viscosity_lu < 1e-14);
            }
            if linear == 0.0 {
                assert!(report.resolved_dissipation_lu < 1e-22);
            }
        }
    }
}

#[test]
fn g4_unresolved_wet_layers_and_invalid_wall_epochs_refuse_without_publication() {
    let gate = CancelGate::new_clock_free();
    let (mut liquid, walls) = fixture(PlanarWall3::ZMin, 2, 0.002, 0.0001, 0.0015, [0.0; 2]);
    let mut kinematics =
        FreeSurfaceKinematicsWorkspace3::new(liquid.grid.active_tiles(), usize::MAX).unwrap();
    let mut workspace = FreeSurfaceSgsWorkspace3::default();
    let config = config(SgsModel::Wale { cw: 0.4 });
    let before = liquid.grid.state_bits();
    let before_mass = liquid.mass_snapshot();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        liquid
            .planar_wall_velocity_field_gradient_sgs_candidate(
                motion(PlanarWall3::ZMin, &walls),
                controls(),
                &cancelled,
                config,
                &mut kinematics,
                &mut workspace,
            )
            .is_err()
    );
    for invalid in [
        vec![[0.0; 2]; 63],
        vec![[f64::NAN, 0.0]; 64],
        vec![[0.3, 0.0]; 64],
    ] {
        assert!(
            liquid
                .planar_wall_velocity_field_gradient_sgs_candidate(
                    motion(PlanarWall3::ZMin, &invalid),
                    controls(),
                    &gate,
                    config,
                    &mut kinematics,
                    &mut workspace,
                )
                .is_err()
        );
    }
    assert_eq!(liquid.grid.state_bits(), before);
    assert_eq!(liquid.mass_snapshot(), before_mass);
    liquid
        .planar_wall_velocity_field_gradient_sgs_candidate(
            motion(PlanarWall3::ZMin, &walls),
            controls(),
            &gate,
            config,
            &mut kinematics,
            &mut workspace,
        )
        .unwrap();
    // Empty interface carriers are not usable velocity samples. A wall and
    // one occupied center alone cannot supply second-order SGS derivatives.
    for &key in liquid.grid.active_keys() {
        for lane in 0..TILE_CELLS {
            if liquid.tiles[&key].cells[lane] == Cell3::Interface {
                liquid.tiles.get_mut(&key).unwrap().mass[lane] = 0.0;
            }
        }
    }
    assert!(matches!(
        liquid.prepare_gradient_sgs_with_wall(
            &gate,
            config,
            Some(motion(PlanarWall3::ZMin, &walls)),
            &mut kinematics,
            &mut workspace
        ),
        Err(FreeSurfaceSgsError3::Kinematics(
            FreeSurfaceKinematicsError3::MissingDerivative { axis: 2, .. }
        ))
    ));
    let (one_layer, walls) = fixture(PlanarWall3::ZMin, 1, 0.0, 0.0, 0.0, [0.0; 2]);
    assert!(
        one_layer
            .prepare_gradient_sgs_with_wall(
                &gate,
                config,
                Some(motion(PlanarWall3::ZMin, &walls)),
                &mut kinematics,
                &mut workspace
            )
            .is_err()
    );
    let mut short = FreeSurfaceKinematicsWorkspace3::new(0, 0).unwrap();
    assert!(matches!(
        one_layer.prepare_gradient_sgs_with_wall(
            &gate,
            config,
            Some(motion(PlanarWall3::ZMin, &walls)),
            &mut short,
            &mut workspace
        ),
        Err(FreeSurfaceSgsError3::Kinematics(
            FreeSurfaceKinematicsError3::Kinematics(KinematicsError::Budget { .. })
        ))
    ));
}
