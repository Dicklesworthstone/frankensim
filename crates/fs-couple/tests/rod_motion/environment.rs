//! Material-fluid substitutions reach the actual native time and balance owner.

use super::*;
use fs_couple::rod_motion::SpatialRodMotion;
use fs_material::gas::{GasSpec, GasState};
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, ResolvedLiquidState, resolve_liquid_state};
use fs_material::state_point::SPECIFIC_HEAT_CAPACITY_DIMS;
use fs_qty::Dims;
use fs_qty::semantic::{QuantityKind, QuantitySpec, SemanticType, ValueForm};
use fs_solid::rod::material_field::SpatialRodError;
use fs_solid::rod::material_field::gravity::ReferenceDisplacementVolume;

fn card(phase: &str, values: &[(&str, f64, Dims)]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for &(name, value, dims) in values {
        let key = if name == "thermal_conductivity" {
            PropertyKey::with_quantity(
                name,
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::ThermalConductivity,
                    ValueForm::Static,
                )),
            )
        } else {
            PropertyKey::new(name, dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured rod/fluid substitution".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: phase.into(),
            process: "analytic fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn liquid(density: f64) -> ResolvedLiquidState {
    let card = card(
        "liquid",
        &[
            ("density", density, Density::DIMS),
            (
                "specific_heat_capacity",
                4200.0,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            ("thermal_conductivity", 0.6, Dims([1, 1, -3, -1, 0, 0])),
            ("dynamic_viscosity", 0.001, DYNAMIC_VISCOSITY_DIMS),
        ],
    );
    resolve_liquid_state(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn other_section(density: f64) -> ConcentricShearSection {
    let values = [
        ("density", density, Density::DIMS),
        ("young_modulus", 1e4, Pressure::DIMS),
        ("shear_modulus", 4e3, Pressure::DIMS),
    ];
    let card = card("solid", &values);
    let requirements = values
        .into_iter()
        .map(|(name, _, dims)| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap()
        })
        .collect::<Vec<_>>();
    let material = resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ConcentricElasticSection::try_new(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: 0.05,
            material,
        }],
    )
    .unwrap()
    .with_shear_modulus("shear_modulus", 0.8)
    .unwrap()
}

fn submerged(density: f64) -> Arc<PreparedRodMotion> {
    let gate = CancelGate::new();
    Arc::new(
        SpatialMaterialRod::try_straight(0.8, vec![section(1.0); 2])
            .unwrap()
            .with_gravity([0.0, 0.0, -9.0], &gate)
            .unwrap()
            .with_submerged_fluid(
                &liquid(density),
                ReferenceDisplacementVolume::FloodedVoids,
                &gate,
            )
            .unwrap()
            .into_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    )
}

fn vector_near(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    for axis in 0..3 {
        near(actual[axis], expected[axis], tolerance);
    }
}

#[test]
fn liquid_card_density_controls_sinking_floating_and_native_momentum_balance() {
    let gate = CancelGate::new();
    let dry = SpatialMaterialRod::try_straight(0.8, vec![section(1.0); 2])
        .unwrap()
        .into_motion_model(RodRootBoundary::Free, &gate)
        .unwrap();
    let velocity = [0.1, 0.2, -0.15];
    for density in [500.0, 1000.0, 2000.0] {
        let model = submerged(density);
        assert_eq!(model.inertia(), dry.inertia());
        let initial_pose = model.initial_pose();
        let mut motion = RodMotion::new(
            model,
            initial_pose.clone(),
            &[velocity; 3],
            &[[2.0, 0.0, 0.0]; 3],
            0.0,
            &gate,
        )
        .unwrap();
        let acceleration = -9.0 * (1.0 - density / 1000.0);
        for _ in 0..40 {
            let (next, step) = motion
                .step_candidate(
                    0.00025,
                    &[RodNodalWrench::default(); 3],
                    controls(0.00025),
                    &gate,
                )
                .unwrap();
            vector_near(step.linear_balance_residual_n_s, [0.0; 3], 2e-12);
            vector_near(step.angular_balance_residual_nm_s, [0.0; 3], 2e-12);
            near(step.integration_defect_j, 0.0, 2e-12);
            assert_eq!(step.source_work_quadrature_j, 0.0);
            motion = next;
        }
        for node in 0..3 {
            for axis in 0..3 {
                let a = if axis == 2 { acceleration } else { 0.0 };
                near(
                    motion.pose().positions_m[node][axis],
                    initial_pose.positions_m[node][axis]
                        + 0.01 * velocity[axis]
                        + 0.5 * a * 0.01_f64.powi(2),
                    2e-12,
                );
                near(
                    motion.velocities_m_s().unwrap()[node][axis],
                    velocity[axis] + 0.01 * a,
                    2e-12,
                );
            }
        }
    }
}

#[test]
fn fluid_refresh_retains_history_and_accounts_exact_fixed_pose_parameter_work() {
    let gate = CancelGate::new();
    let model = submerged(500.0);
    let reference = model.initial_pose();
    let offset = [0.0, 0.1, 0.3];
    let velocity = [0.04, -0.02, -0.1];
    let mut pose = reference.clone();
    for position in &mut pose.positions_m {
        for axis in 0..3 {
            position[axis] += offset[axis];
        }
    }
    let initial = RodMotion::new(model, pose, &[velocity; 3], &[[0.0; 3]; 3], 0.0, &gate).unwrap();
    let before = evolve(initial, 0.001, 20);
    let (unchanged, noop) = before
        .rebind_submerged_fluid_candidate(&liquid(500.0), &gate)
        .unwrap();
    assert_eq!(noop.environmental_work_j, 0.0);
    assert_eq!(
        unchanged.observe([0.0; 3]).unwrap(),
        before.observe([0.0; 3]).unwrap()
    );
    let (refreshed, report) = before
        .rebind_submerged_fluid_candidate(&liquid(2000.0), &gate)
        .unwrap();
    assert_eq!(refreshed.pose(), before.pose());
    assert_eq!(
        refreshed.linear_momenta_kg_m_s(),
        before.linear_momenta_kg_m_s()
    );
    assert_eq!(refreshed.body_spins_kg_m2_s(), before.body_spins_kg_m2_s());
    assert_eq!(refreshed.model().inertia(), before.model().inertia());
    assert_eq!(
        refreshed.model().source().reference_mass(),
        before.model().source().reference_mass()
    );
    assert_eq!(refreshed.time_s(), before.time_s());
    assert_eq!(refreshed.accepted_steps(), before.accepted_steps());
    assert_eq!(report.time_s, before.time_s());
    let volumes = before
        .model()
        .source()
        .submerged_fluid_load()
        .unwrap()
        .nodal_volumes_m3();
    let expected_work = 1500.0
        * volumes
            .iter()
            .enumerate()
            .map(|(node, volume)| {
                -9.0 * volume
                    * (before.pose().positions_m[node][2] - reference.positions_m[node][2])
            })
            .sum::<f64>();
    near(report.environmental_work_j, expected_work, 2e-12);
    let old = before.observe([0.0; 3]).unwrap();
    let new = refreshed.observe([0.0; 3]).unwrap();
    assert_eq!(old.kinetic_energy_j, new.kinetic_energy_j);
    near(
        new.potential_energy_j - old.potential_energy_j,
        report.environmental_work_j,
        2e-12,
    );
    let (restored, undo) = refreshed
        .rebind_submerged_fluid_candidate(&liquid(500.0), &gate)
        .unwrap();
    near(
        undo.environmental_work_j + report.environmental_work_j,
        0.0,
        2e-12,
    );
    assert_eq!(restored.observe([0.0; 3]).unwrap(), old);
    let after = evolve(refreshed, 0.001, 20);
    // Independent piecewise constant acceleration, including the velocity at
    // the density switch. Restarting the trajectory would miss the cross term.
    for node in 0..3 {
        for axis in 0..3 {
            let a = if axis == 2 { -4.5 } else { 0.0 };
            let b = if axis == 2 { 9.0 } else { 0.0 };
            let exact = reference.positions_m[node][axis]
                + offset[axis]
                + 0.04 * velocity[axis]
                + (0.5 * a + a + 0.5 * b) * 0.02_f64.powi(2);
            near(after.pose().positions_m[node][axis], exact, 2e-12);
            near(
                after.velocities_m_s().unwrap()[node][axis],
                velocity[axis] + 0.02 * (a + b),
                2e-12,
            );
        }
    }
    assert_eq!(after.accepted_steps(), 40);
}

#[test]
fn eccentric_displacement_moment_uses_actual_frames_and_not_the_mass_centroid() {
    let gate = CancelGate::new();
    // Separate outer envelopes: the hydrostatic owner assumes nonoverlap.
    let offset = [0.0, 0.12, 0.08];
    let gravity = [0.7, -0.5, -2.0];
    let link_section = other_section(2000.0);
    let link_volume = link_section.elastic().area_m2() * 0.9;
    let link = MaterialTensionLink::try_new(
        link_section.elastic().clone(),
        [
            TensionLinkAttachment {
                node: 0,
                offset_body_m: offset,
            },
            TensionLinkAttachment {
                node: 1,
                offset_body_m: offset,
            },
        ],
        0.9,
        1,
        0.2,
    )
    .unwrap();
    let rod_section = section(1.0);
    let rod_volume = rod_section.elastic().area_m2() * 0.8;
    let source = SpatialMaterialRod::try_straight(0.8, vec![rod_section])
        .unwrap()
        .with_tension_links(vec![link], &gate)
        .unwrap()
        .with_gravity(gravity, &gate)
        .unwrap()
        .with_submerged_fluid(
            &liquid(500.0),
            ReferenceDisplacementVolume::FloodedVoids,
            &gate,
        )
        .unwrap();
    let model = Arc::new(
        source
            .into_spatial_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let mut pose = model.initial_pose();
    pose.frames
        .fill(So3::exp(So3Tangent::new(Vec3::new(0.2, 0.0, 0.0))).unwrap());
    let loads = model.body_loads(&pose, &gate).unwrap();
    let effective_mass = 0.5 * (rod_volume * 500.0 + link_volume * 1500.0);
    let endpoint_force = gravity.map(|g| 0.5 * link_volume * 1500.0 * g);
    for node in 0..2 {
        vector_near(
            loads.force_n[node],
            gravity.map(|g| effective_mass * g),
            2e-13,
        );
        let local = pose.frames[node]
            .inverse()
            .rotate(Vec3::new(
                endpoint_force[0],
                endpoint_force[1],
                endpoint_force[2],
            ))
            .unwrap();
        let torque = Vec3::new(offset[0], offset[1], offset[2]).cross(local);
        vector_near(
            loads.torque_body_nm[node],
            [torque.x, torque.y, torque.z],
            2e-13,
        );
        // The center of displaced volume differs from the center of material
        // mass because the two sections have different solid densities.
        let center = model.inertia()[node].center_of_mass_body_m();
        let naive = Vec3::new(center[0], center[1], center[2]).cross(
            pose.frames[node]
                .inverse()
                .rotate(Vec3::new(
                    loads.force_n[node][0],
                    loads.force_n[node][1],
                    loads.force_n[node][2],
                ))
                .unwrap(),
        );
        assert!((naive.z - torque.z).abs() > 1e-4);
    }
    let mut motion = SpatialRodMotion::new_spatial(
        model,
        pose,
        &[[0.05, -0.03, 0.02]; 2],
        &[[0.1, 0.2, -0.1]; 2],
        0.0,
        1e-12,
        &gate,
    )
    .unwrap();
    let initial = motion.observe([0.0; 3]).unwrap();
    let mut defects = 0.0;
    let mut parameter_work = 0.0;
    for interval in 0..40 {
        if interval == 20 {
            let (next, change) = motion
                .rebind_submerged_fluid_candidate(&liquid(1600.0), &gate)
                .unwrap();
            assert_eq!(next.model().inertia(), motion.model().inertia());
            assert_eq!(next.body_spins_kg_m2_s(), motion.body_spins_kg_m2_s());
            parameter_work += change.environmental_work_j;
            motion = next;
        }
        let (next, report) = motion
            .step_candidate(1e-5, &[RodNodalWrench::default(); 2], controls(1e-5), &gate)
            .unwrap();
        vector_near(report.linear_balance_residual_n_s, [0.0; 3], 3e-12);
        vector_near(report.angular_balance_residual_nm_s, [0.0; 3], 3e-12);
        defects += report.integration_defect_j;
        motion = next;
    }
    let final_state = motion.observe([0.0; 3]).unwrap();
    near(
        final_state.kinetic_energy_j + final_state.potential_energy_j
            - initial.kinetic_energy_j
            - initial.potential_energy_j,
        parameter_work + defects,
        3e-12,
    );
}

#[test]
fn gas_pressure_refresh_preserves_supports_and_rejected_refresh_preserves_motion() {
    let gate = CancelGate::new();
    let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 1e5).unwrap();
    let compressed = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 2e5).unwrap();
    let source = SpatialMaterialRod::try_straight(0.8, vec![section(1.0)])
        .unwrap()
        .with_supports(
            &[RodSupport {
                node: 1,
                position_m: [Some(0.8), Some(0.0), Some(0.0)],
                frame: Some(So3::identity()),
            }],
            &gate,
        )
        .unwrap()
        .with_gravity([0.0, 0.0, -9.0], &gate)
        .unwrap()
        .with_submerged_fluid(&air, ReferenceDisplacementVolume::SealedEnvelope, &gate)
        .unwrap();
    let model = Arc::new(
        source
            .into_motion_model(RodRootBoundary::Clamped, &gate)
            .unwrap(),
    );
    let pose = model.initial_pose();
    let state = RodMotion::new(
        model,
        pose.clone(),
        &[[0.0; 3]; 2],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap();
    let (refreshed, report) = state
        .rebind_submerged_fluid_candidate(&compressed, &gate)
        .unwrap();
    assert_eq!(report.environmental_work_j, 0.0);
    near(
        report.after.density_kg_m3(),
        2.0 * report.before.density_kg_m3(),
        1e-13,
    );
    let forces = refreshed
        .model()
        .body_loads(refreshed.pose(), &gate)
        .unwrap();
    let half_volume = 0.4 * std::f64::consts::PI * 0.05_f64.powi(2);
    for node in 0..2 {
        let weight_minus_displaced_mass =
            refreshed.model().inertia()[node].mass_kg() - compressed.density * half_volume;
        vector_near(
            forces.force_n[node],
            [0.0, 0.0, -9.0 * weight_minus_displaced_mass],
            2e-13,
        );
    }
    let (next, step) = refreshed
        .step_candidate(
            0.001,
            &[RodNodalWrench::default(); 2],
            controls(0.001),
            &gate,
        )
        .unwrap();
    assert_eq!(next.pose(), &pose);
    for node in 0..2 {
        vector_near(
            step.support_impulses[node].linear_n_s,
            forces.force_n[node].map(|f| -0.001 * f),
            2e-13,
        );
    }
    vector_near(step.linear_balance_residual_n_s, [0.0; 3], 2e-13);
    vector_near(step.angular_balance_residual_nm_s, [0.0; 3], 2e-13);
    let saved = next.observe([0.0; 3]).unwrap();
    let saved_time = next.time_s();
    let mut corrupt = air;
    corrupt.density *= 2.0;
    assert!(
        next.rebind_submerged_fluid_candidate(&corrupt, &gate)
            .is_err()
    );
    let dry_model = axial_model(false);
    let dry = RodMotion::new(
        dry_model.clone(),
        dry_model.initial_pose(),
        &[[0.0; 3]; 2],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap();
    assert!(dry.rebind_submerged_fluid_candidate(&air, &gate).is_err());
    gate.request();
    assert!(matches!(
        next.rebind_submerged_fluid_candidate(&compressed, &gate),
        Err(RodMotionError::Solid(SpatialRodError::Cancelled))
    ));
    assert_eq!(next.observe([0.0; 3]).unwrap(), saved);
    assert_eq!(next.pose(), &pose);
    assert_eq!(next.time_s(), saved_time);
    assert_eq!(next.accepted_steps(), 1);
    assert_eq!(
        next.model()
            .source()
            .submerged_fluid_load()
            .unwrap()
            .transport_state(),
        report.after
    );
}
