//! Complete original material tensors on native rigid work coordinates.
use fs_material::{
    Voigt,
    state_point::MaterialPropertySelection,
    visco_tensor::{
        TensorProjectionBudget, TensorProny, TensorPronyBranch,
        source::{
            ResolvedTensorMaterial, TensorBranchProperties, TensorElasticity, TensorPronyProperties,
        },
    },
};
use fs_mbd::{
    MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3,
    spatial_loads::{
        internal_impulse::{
            InternalWrenchImpulseControls,
            viscous::{
                BushingCoordinateScales, ViscousConnectorControls,
                network::{
                    driven::DrivenViscousNetworkLimits,
                    relaxation::{
                        PreparedRelaxationNetwork, RelaxationNetworkLimits, SpatialMemoryEvolution,
                    },
                },
            },
        },
        internal_motion::viscous::{
            driven::HeldInternalEffort,
            network::relaxation::{RelaxationNetworkMotionControls, step_relaxation_network},
        },
    },
};
use fs_mbd_material::spatial_relaxation::{
    TensorSpatialError, TensorSpatialGeometry, TensorSpatialRelaxation,
};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn budget() -> TensorProjectionBudget {
    TensorProjectionBudget {
        max_coordinates: 6,
        max_scalar_terms: 100_000,
    }
}
fn geometry() -> TensorSpatialGeometry {
    TensorSpatialGeometry {
        bodies: [0, 1],
        arms_body_m: [Vec3::new(0.02, -0.01, 0.03), Vec3::new(-0.01, 0.03, 0.02)],
        strain_columns: [
            [1.2, 0.1, 0.0, -0.2, 0.0, 0.03],
            [0.2, -0.8, 0.1, 0.25, -0.05, 0.0],
            [-0.1, 0.1, 0.6, 0.0, 0.2, -0.1],
            [0.02, 0.0, 0.03, 0.04, 0.0, -0.03],
            [0.0, -0.04, 0.01, 0.0, 0.03, 0.02],
            [0.03, 0.01, -0.02, -0.01, 0.02, 0.04],
        ],
        volume_m3: 0.02,
        coordinate_scales: BushingCoordinateScales {
            translation_m: 0.003,
            rotation_rad: 0.07,
        },
    }
}
fn gram(rows: [[f64; 6]; 6]) -> [[f64; 6]; 6] {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| (0..6).map(|r| rows[r][i] * rows[r][j]).sum())
    })
}
fn law() -> TensorProny {
    let cross = [0.4, -0.3, 0.2, 0.1, -0.2, 0.15];
    let equilibrium = core::array::from_fn(|i| {
        core::array::from_fn(|j| cross[i] * cross[j] + if i == j { 20.0 + i as f64 } else { 0.0 })
    });
    let mut rows = [[0.0; 6]; 6];
    rows[0] = [0.8, -0.2, 0.4, 0.15, -0.1, 0.3];
    rows[1] = rows[0].map(|x| 2.0 * x);
    rows[3] = [0.1, 0.3, -0.2, 0.4, 0.05, -0.1];
    let second = core::array::from_fn(|i| {
        core::array::from_fn(|j| if i == j { 0.3 + 0.1 * i as f64 } else { 0.0 })
    });
    let viscosity = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            0.04 * (cross[i] * cross[j]) + if i == j { 0.03 + 0.01 * i as f64 } else { 0.0 }
        })
    });
    TensorProny::new_mandel(
        equilibrium,
        vec![
            TensorPronyBranch::from_gram_rows(rows, 0.2).unwrap(),
            TensorPronyBranch::from_gram_rows(second, 3.0).unwrap(),
        ],
        0.2,
    )
    .unwrap()
    .with_kelvin_viscosity_mandel(viscosity)
    .unwrap()
    .with_band((0.0, 30.0))
    .unwrap()
}
fn mandel(v: Voigt) -> [f64; 6] {
    core::array::from_fn(|i| v[i] * if i < 3 { 1.0 } else { 2.0_f64.sqrt() })
}

#[test]
fn complete_full_tensor_spectrum_and_work_map_preserve_original_factor_rows() {
    let original = law();
    let geometry = geometry();
    let bound =
        TensorSpatialRelaxation::from_material(&original, geometry, budget(), || false).unwrap();
    assert_eq!(bound.projection().material(), &original);
    assert_eq!(bound.native_link().viscous.bodies, geometry.bodies);
    assert_eq!(
        bound.native_link().viscous.arms_body_m,
        geometry.arms_body_m
    );
    assert_eq!(bound.native_link().branches.len(), 3);
    assert_eq!(
        bound.native_link().branches[0].evolution(),
        SpatialMemoryEvolution::Equilibrium
    );
    assert_eq!(
        bound.native_link().branches[1].relaxation_time_s(),
        Some(0.2)
    );
    assert_eq!(
        bound.native_link().branches[2].relaxation_time_s(),
        Some(3.0)
    );
    for i in 0..6 {
        for j in 0..6 {
            assert_eq!(
                bound.native_link().branches[0].gram_rows()[i][j],
                bound.projection().equilibrium_gram_rows()[i][j]
            );
            for b in 0..2 {
                assert_eq!(
                    bound.native_link().branches[b + 1].gram_rows()[i][j],
                    bound.projection().branches()[b].rows()[i][j]
                );
            }
        }
    }
    for omega in [0.0, 0.7, 4.0, 29.0] {
        bound.admit_frequency_band([omega, omega]).unwrap();
        let (real, imag) = original.modulus_mandel(omega).unwrap();
        for i in 0..6 {
            for j in 0..6 {
                let left = mandel(geometry.strain_columns[i]);
                let right = mandel(geometry.strain_columns[j]);
                let project = |matrix: [[f64; 6]; 6]| -> f64 {
                    geometry.volume_m3
                        * (0..6)
                            .map(|r| {
                                (0..6)
                                    .map(|c| left[r] * matrix[r][c] * right[c])
                                    .sum::<f64>()
                            })
                            .sum::<f64>()
                };
                let mut storage = gram(*bound.native_link().branches[0].gram_rows())[i][j];
                let mut loss =
                    omega * bound.native_link().viscous.damping.native_coefficients()[i][j];
                for branch in &bound.native_link().branches[1..] {
                    let z = omega * branch.relaxation_time_s().unwrap();
                    let k = gram(*branch.gram_rows())[i][j];
                    storage += k * z * z / (1.0 + z * z);
                    loss += k * z / (1.0 + z * z);
                }
                close(storage, project(real), 5e-13);
                close(loss, project(imag), 5e-13);
            }
        }
    }
    let q = [0.01, -0.02, 0.015, 0.005, -0.004, 0.006];
    let strain = bound.strain(q).unwrap();
    let stress = [1.0, 2.0, -3.0, 0.5, -0.7, 0.2];
    let force = bound.projection().conjugate_force(&stress).unwrap();
    let work = (0..6).map(|i| force[i] * q[i]).sum::<f64>();
    let physical = geometry.volume_m3
        * (0..6)
            .map(|i| stress[i] * strain[i] * if i < 3 { 1.0 } else { 2.0 })
            .sum::<f64>();
    close(work, physical, 2e-17);
    let initial = bound
        .memory_from_elastic_coordinates(q, &[q, q.map(|x| -0.4 * x)])
        .unwrap();
    for b in 0..3 {
        let displacement = if b == 2 { q.map(|x| -0.4 * x) } else { q };
        assert_eq!(
            initial[b],
            bound.native_link().branches[b]
                .elastic_memory(displacement)
                .unwrap()
        );
    }
    let relaxed = bound.relaxed_memory(q).unwrap();
    assert_eq!(relaxed[0], initial[0]);
    assert_eq!(relaxed[1..], [[0.0; 6]; 2]);
}

fn controls() -> ViscousConnectorControls {
    ViscousConnectorControls {
        impulse: InternalWrenchImpulseControls {
            angular_origin_world_m: Vec3::ZERO,
            linear_tolerance_ns: 3e-10,
            angular_tolerance_nms: 3e-10,
            work_tolerance_j: 3e-10,
        },
        linear_speed_tolerance_m_s: 3e-10,
        angular_speed_tolerance_rad_s: 3e-10,
    }
}
fn limits() -> RelaxationNetworkLimits {
    RelaxationNetworkLimits {
        source: DrivenViscousNetworkLimits {
            maximum_bodies: 2,
            maximum_connectors: 1,
            maximum_live_bytes: 128_000_000,
            maximum_preparation_terms: 80_000_000_000,
        },
        maximum_branches: 2,
    }
}
fn axial_material() -> TensorProny {
    let diagonal =
        |x: f64| core::array::from_fn(|i| core::array::from_fn(|j| if i == j { x } else { 0.0 }));
    TensorProny::new_mandel(
        diagonal(200.0),
        vec![TensorPronyBranch::from_gram_rows(diagonal(50.0_f64.sqrt()), 0.2).unwrap()],
        0.05,
    )
    .unwrap()
    .with_kelvin_viscosity_mandel(diagonal(1.0))
    .unwrap()
}
fn axial_geometry() -> TensorSpatialGeometry {
    let mut geometry = geometry();
    geometry.strain_columns = [[0.0; 6]; 6];
    geometry.strain_columns[0][0] = 10.0;
    geometry.volume_m3 = 0.001;
    geometry.arms_body_m = [Vec3::ZERO; 2];
    geometry
}
fn bodies() -> ([RigidBodyState; 2], [MassProperties; 2]) {
    let properties = [MassProperties::new(2.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap(); 2];
    let states = core::array::from_fn(|i| {
        RigidBodyState::new(
            Pose::new(Vec3::new(i as f64, 0.0, 0.0), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(if i == 0 { -0.0006 } else { 0.0006 }, 0.0, 0.0),
            Vec3::ZERO,
        )
        .unwrap()
    });
    (states, properties)
}

#[test]
fn native_fixed_and_moving_receipts_continue_explicit_strain_and_complete_history() {
    let material = axial_material();
    let bound =
        TensorSpatialRelaxation::from_material(&material, axial_geometry(), budget(), || false)
            .unwrap();
    let q = [0.0002, 0.01, -0.02, 0.03, -0.01, 0.02];
    let memory = bound
        .memory_from_elastic_coordinates(q, &[[-0.0001, 0.0, 0.0, 0.0, 0.0, 0.0]])
        .unwrap();
    let bank = [memory.clone()];
    let (states, properties) = bodies();
    let links = [bound.native_link().clone()];
    let efforts = [[0.001, 0.0, 0.0, 0.0, 0.0, 0.0]];
    let mut source = PreparedRelaxationNetwork::new(
        &states.map(|s| s.pose()),
        &properties,
        &links,
        0.05,
        limits(),
        || false,
    )
    .unwrap();
    let out = source
        .advance(&states, &bank, &efforts, controls(), || false)
        .unwrap();
    let next = bound.coordinates_after_source(q, 0, &out).unwrap();
    let eq = links[0].branches[0].elastic_memory(next).unwrap();
    for i in 0..6 {
        close(out.memory_after[0][0][i], eq[i], 3e-15);
    }
    assert_eq!(&next[1..], &q[1..]);
    assert!(out.dissipated_energy_j > 0.0);
    assert_ne!(
        out.memory_after[0][1],
        links[0].branches[1].elastic_memory(next).unwrap()
    );
    let mut half = PreparedRelaxationNetwork::new(
        &states.map(|s| s.pose()),
        &properties,
        &links,
        0.025,
        limits(),
        || false,
    )
    .unwrap();
    let first = half
        .advance(&states, &bank, &efforts, controls(), || false)
        .unwrap();
    let middle = bound.coordinates_after_source(q, 0, &first).unwrap();
    let second = half
        .advance(
            &first.states_after,
            &first.memory_after,
            &efforts,
            controls(),
            || false,
        )
        .unwrap();
    let continued = bound.coordinates_after_source(middle, 0, &second).unwrap();
    for i in 0..6 {
        close(continued[i], next[i], 3e-15);
    }
    close(
        first.dissipated_energy_j + second.dissipated_energy_j,
        out.dissipated_energy_j,
        3e-15,
    );
    let moving = step_relaxation_network(
        &states,
        &properties,
        &links,
        &bank,
        &[HeldInternalEffort::MaterialFirst(efforts[0])],
        0.02,
        RelaxationNetworkMotionControls {
            source: controls(),
            limits: limits(),
        },
        || false,
    )
    .unwrap();
    let moved = bound.coordinates_after_motion(q, 0, &moving).unwrap();
    let eq = links[0].branches[0].elastic_memory(moved).unwrap();
    for i in 0..6 {
        close(moving.memory_after()[0][0][i], eq[i], 3e-15);
    }
    assert_eq!(bank[0], memory);
    assert_eq!(q, [0.0002, 0.01, -0.02, 0.03, -0.01, 0.02]);
}

#[test]
fn material_domain_geometry_budget_and_half_endpoint_refusals_are_atomic() {
    let material = axial_material();
    let original = material.clone();
    let geometry = axial_geometry();
    let bound =
        TensorSpatialRelaxation::from_material(&material, geometry, budget(), || false).unwrap();
    assert!(bound.strain([0.006, 0.0, 0.0, 0.0, 0.0, 0.0]).is_err());
    assert!(
        bound
            .memory_from_elastic_coordinates([0.0; 6], &[])
            .is_err()
    );
    assert!(
        bound
            .memory_from_elastic_coordinates([0.0; 6], &[[0.006; 6]])
            .is_err()
    );
    for bad in [
        TensorSpatialGeometry {
            volume_m3: 0.0,
            ..geometry
        },
        TensorSpatialGeometry {
            bodies: [1, 1],
            ..geometry
        },
        TensorSpatialGeometry {
            arms_body_m: [Vec3::new(f64::NAN, 0.0, 0.0), Vec3::ZERO],
            ..geometry
        },
    ] {
        assert!(
            TensorSpatialRelaxation::from_material(&material, bad, budget(), || false).is_err()
        );
    }
    assert!(
        TensorSpatialRelaxation::from_material(
            &material,
            geometry,
            TensorProjectionBudget {
                max_scalar_terms: 1,
                ..budget()
            },
            || false
        )
        .is_err()
    );
    assert!(matches!(
        TensorSpatialRelaxation::from_material(&material, geometry, budget(), || true),
        Err(TensorSpatialError::Cancelled)
    ));
    let (states, properties) = bodies();
    let links = [bound.native_link().clone()];
    let bank = [bound.relaxed_memory([0.0; 6]).unwrap()];
    let moving = step_relaxation_network(
        &states,
        &properties,
        &links,
        &bank,
        &[HeldInternalEffort::MaterialFirst([0.0; 6])],
        0.02,
        RelaxationNetworkMotionControls {
            source: controls(),
            limits: limits(),
        },
        || false,
    )
    .unwrap();
    let accepted = bound
        .coordinates_after_motion([0.0; 6], 0, &moving)
        .unwrap();
    let mut outside = moving.clone();
    // A complete receipt can have a valid final strain while a resolved source
    // half already exceeds the material domain. Both halves must be checked.
    outside.half_sources[0].average_relative_twist[0][0] = 1.0;
    outside.half_sources[1].average_relative_twist[0][0] = -1.0;
    assert!(
        bound
            .coordinates_after_motion([0.0; 6], 0, &outside)
            .is_err()
    );
    assert!(
        bound
            .coordinates_after_source([0.0; 6], 1, &moving.half_sources[0])
            .is_err()
    );
    assert_eq!(
        bound
            .coordinates_after_motion([0.0; 6], 0, &moving)
            .unwrap(),
        accepted
    );
    assert_eq!(material, original);
}

#[test]
fn calibrated_source_rotation_receipts_and_original_frequency_domain_survive() {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, QueryPoint, UncertaintyModel,
    };
    use fs_qty::{Density, Dims, Pressure};
    let scalar = |name: &str, dims: Dims, value: f64| PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained()
            .with("T", 280.0, 330.0)
            .with("omega", 0.5, 20.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured complete spatial material".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    };
    let branch = TensorBranchProperties {
        stiffness_upper: core::array::from_fn(|i| format!("arm_c{i}")),
        relaxation_time: "arm_tau".into(),
    };
    let properties = TensorPronyProperties {
        elasticity: TensorElasticity::Isotropic,
        strain_limit: "small_strain".into(),
        branches: vec![branch.clone()],
    };
    let mut claims = ClaimSet::new();
    for claim in [
        scalar("young_modulus", Pressure::DIMS, 30.0),
        scalar("poisson_ratio", Dims::NONE, 0.25),
        scalar("density", Density::DIMS, 900.0),
        scalar("small_strain", Dims::NONE, 0.1),
        scalar("arm_tau", Dims([0, 0, 1, 0, 0, 0]), 0.4),
    ] {
        claims.insert_claim(claim).unwrap();
    }
    let cross = [0.4, -0.3, 0.2, 0.1, -0.2, 0.15];
    let mut k = 0;
    for i in 0..6 {
        for j in i..6 {
            claims
                .insert_claim(scalar(
                    &branch.stiffness_upper[k],
                    Pressure::DIMS,
                    cross[i] * cross[j] + if i == j { 2.0 + i as f64 } else { 0.0 },
                ))
                .unwrap();
            k += 1;
        }
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured spatial law".into(),
            phase: "solid".into(),
            process: "source bridge test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with("omega", 2.0)
        .unwrap();
    let source = ResolvedTensorMaterial::resolve(
        &card,
        &query,
        properties.clone(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let before = source.clone();
    let rotation = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    let bound =
        TensorSpatialRelaxation::from_source(&source, rotation, geometry(), budget(), || false)
            .unwrap();
    assert!(core::ptr::eq(bound.source().unwrap(), &source));
    assert_eq!(bound.source().unwrap().properties(), &properties);
    assert_eq!(bound.material_to_first_body(), Some(rotation));
    assert_eq!(
        bound.projection().material(),
        &source.rotated_law(rotation).unwrap()
    );
    assert_eq!(
        bound.native_link().branches[1].relaxation_time_s(),
        Some(0.4)
    );
    bound.admit_frequency_band([0.5, 20.0]).unwrap();
    assert!(
        bound
            .relaxed_memory([0.001, 0.0, 0.0, 0.0, 0.0, 0.0])
            .is_err()
    );
    assert!(bound.relaxed_memory([0.0; 6]).is_ok());
    for band in [[0.0, 20.0], [0.5, 20.1], [2.0, 1.0], [f64::NAN, 10.0]] {
        assert!(bound.admit_frequency_band(band).is_err());
    }
    assert!(
        TensorSpatialRelaxation::from_source(&source, [[0.0; 3]; 3], geometry(), budget(), || {
            false
        })
        .is_err()
    );
    assert_eq!(source, before);
}
