//! Native consistent momentum with finite plasticity, damage and physical work.
use super::*;
use fs_material::ductile::{JohnsonCookInitiation, JohnsonCookParameters};
use fs_solid::hyper3::plastic::ductile::dynamics::{
    FiniteDuctileMotionControls, FiniteDuctileTetMotion,
};

const RHO: f64 = 12.0;
const VOLUME: f64 = 1.0 / 6.0;

fn motion_controls() -> FiniteDuctileMotionControls {
    FiniteDuctileMotionControls {
        material: controls(),
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}
fn at(dt: f64) -> DuctileStepConditions<'static> {
    DuctileStepConditions {
        duration_s: dt,
        ..conditions()
    }
}
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn rate_model(width: f64) -> FiniteDuctileJ2<JohnsonCookInitiation> {
    let base = model(width, 1.0);
    let initiation = JohnsonCookInitiation::new(
        JohnsonCookParameters {
            d1: ONSET,
            d2: 0.0,
            d3: 0.0,
            d4: 0.05,
            d5: 0.0,
            reference_strain_rate_s_inv: 1.0,
            transition_temperature_k: 299.0,
            melting_temperature_k: 301.0,
        },
        DuctileDomain {
            triaxiality: [-20.0, 20.0],
            strain_rate_s_inv: [1e-8, 10.0],
            temperature_k: [300.0, 300.0],
        },
    )
    .unwrap();
    FiniteDuctileJ2::new(base.native().clone(), initiation, base.settings()).unwrap()
}

#[test]
fn free_body_uses_original_uniform_acceleration_and_retains_mass_through_continuation() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let material = model(bands[0].band_width_m(), 0.04);
        let materials = [&material];
        let solid = FiniteDuctileTetSolid::try_new(
            &geometry,
            &materials,
            &bands,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let v0 = [0.25, -0.125, 0.0625];
        let moving = FiniteDuctileTetMotion::try_new(solid, &[RHO], &[v0; 4], cx).unwrap();
        let dt = 0.125;
        let acceleration = [0.5, -0.25, 0.125];
        let body = acceleration.map(|a| RHO * a);
        let first = moving
            .advance(
                &[],
                &[],
                body,
                at(dt),
                IncrementSettings::default(),
                motion_controls(),
                cx,
            )
            .unwrap();
        let v1 = core::array::from_fn(|i| v0[i] + dt * acceleration[i]);
        for (u, v) in first
            .motion
            .solid()
            .displacement_m()
            .iter()
            .zip(first.motion.velocity_mps())
        {
            for axis in 0..3 {
                near(u[axis], dt * v1[axis]);
                near(v[axis], v1[axis]);
            }
        }
        let mass = RHO * VOLUME;
        near(first.inertia.kinetic_before_j, 0.5 * mass * dot3(v0, v0));
        near(first.inertia.kinetic_after_j, 0.5 * mass * dot3(v1, v1));
        near(
            first.inertia.numerical_loss_j,
            0.5 * mass * dt * dt * dot3(acceleration, acceleration),
        );
        let work = mass * dt * dot3(acceleration, v1);
        near(first.inertia.inertial_work_j, work);
        near(first.material_energy.dead_load_endpoint_work_j, work);
        near(first.material_energy.support_endpoint_work_j, 0.0);
        near(first.material_energy.plastic_dissipation_j, 0.0);
        near(first.material_energy.damage_release_j, 0.0);
        near(first.balance_error_j, 0.0);
        assert!(first.equilibrium.reactions_n.is_empty());
        assert_eq!(
            first.motion.solid().element_states()[0].native_state(),
            material.initial_state().native_state()
        );
        near(first.motion.solid().element_states()[0].integrity(), 1.0);
        for (actual, expected) in first.motion.solid().element_deformation()[0]
            .iter()
            .zip([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])
        {
            near(*actual, expected);
        }
        assert!(core::ptr::eq(
            moving.reference_density_kg_m3(),
            first.motion.reference_density_kg_m3()
        ));
        near(first.motion.reference_mass_kg(), mass);
        let second = first
            .motion
            .advance(
                &[],
                &[],
                body,
                at(dt),
                IncrementSettings::default(),
                motion_controls(),
                cx,
            )
            .unwrap();
        let v2 = core::array::from_fn::<_, 3, _>(|i| v1[i] + dt * acceleration[i]);
        for (u, v) in second
            .motion
            .solid()
            .displacement_m()
            .iter()
            .zip(second.motion.velocity_mps())
        {
            for axis in 0..3 {
                near(u[axis], dt * (v1[axis] + v2[axis]));
                near(v[axis], v2[axis]);
            }
        }
        near(second.motion.reference_mass_kg(), mass);
        assert_eq!(moving.velocity_mps(), &[v0; 4]);
        assert_eq!(moving.solid().displacement_m(), &[[0.0; 3]; 4]);
    });
}

#[test]
fn free_coordinate_solves_rate_sensitive_damage_and_inertia_in_one_physical_interval() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let material = rate_model(bands[0].band_width_m());
        let materials = [&material];
        let supports: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let lambda = 1.02_f64;
        let alpha = (1000.0 * lambda.ln() - 3.0) / 1500.0;
        let effective_mass = RHO * VOLUME / 10.0;
        let mut integrities = Vec::new();
        for dt in [0.0625, 0.125] {
            let solid = FiniteDuctileTetSolid::try_new(
                &geometry,
                &materials,
                &bands,
                TetAssemblyBudget::standard(),
                cx,
            )
            .unwrap();
            let moving =
                FiniteDuctileTetMotion::try_new(solid, &[RHO], &[[0.0; 3]; 4], cx).unwrap();
            // Independent confined logarithmic J2 plus authored JC onset.
            // Constant reference rho gives M_11=rho*V/10 for this free DOF.
            let onset = ONSET * (1.0 + 0.05 * (alpha / dt).ln());
            let g = (-bands[0].band_width_m() * 3.0 * (alpha - onset) / GP).exp();
            let piola = g * (1000.0 * lambda.ln() + 2.0) / lambda;
            let acceleration = (lambda - 1.0) / (dt * dt);
            let force = [(3, VOLUME * piola + effective_mass * acceleration)];
            let step = moving
                .advance(
                    &supports,
                    &force,
                    [0.0; 3],
                    at(dt),
                    IncrementSettings::default(),
                    motion_controls(),
                    cx,
                )
                .unwrap();
            near(step.motion.solid().displacement_m()[1][0], lambda - 1.0);
            near(step.motion.velocity_mps()[1][0], (lambda - 1.0) / dt);
            near(
                step.motion.solid().element_states()[0].native_state().alpha,
                alpha,
            );
            near(step.motion.solid().element_states()[0].integrity(), g);
            near(step.motion.solid().element_first_piola_pa()[0][0], piola);
            assert!(step.motion.solid().element_states()[0].damage() > 0.1);
            assert!(step.equilibrium.residual_history.len() > 1);
            let kinetic = 0.5 * effective_mass * ((lambda - 1.0) / dt).powi(2);
            near(step.inertia.kinetic_after_j, kinetic);
            near(step.inertia.numerical_loss_j, kinetic);
            near(step.inertia.inertial_work_j, 2.0 * kinetic);
            let stored = VOLUME * g * (0.5 * 1000.0 * lambda.ln().powi(2) + 9.0 / 3000.0);
            let plastic = VOLUME * g * 3.0 * alpha;
            let internal = VOLUME * piola * (lambda - 1.0);
            near(step.material_energy.stored_after_j, stored);
            near(step.material_energy.plastic_dissipation_j, plastic);
            near(step.material_energy.damage_release_j, 0.0);
            near(
                step.material_energy.endpoint_work_defect_j,
                internal - stored - plastic,
            );
            near(
                step.material_energy.dead_load_endpoint_work_j,
                force[0].1 * (lambda - 1.0),
            );
            near(step.material_energy.support_endpoint_work_j, 0.0);
            near(step.balance_error_j, 0.0);
            let reaction = supports
                .iter()
                .zip(&step.equilibrium.reactions_n)
                .find(|((dof, _), _)| *dof == 6)
                .unwrap()
                .1;
            near(*reaction, RHO * VOLUME / 20.0 * acceleration);
            integrities.push(g);

            let strict_material = FiniteDuctileMotionControls {
                material: FinitePlasticStepControls {
                    maximum_endpoint_work_defect_j: 0.5
                        * step.material_energy.absolute_endpoint_work_defect_j,
                    ..controls()
                },
                ..motion_controls()
            };
            assert!(matches!(
                moving.advance(
                    &supports,
                    &force,
                    [0.0; 3],
                    at(dt),
                    IncrementSettings::default(),
                    strict_material,
                    cx
                ),
                Err(FiniteDuctileTetError::Energy {
                    quantity: "sum of absolute element integration defects",
                    ..
                })
            ));
            assert_eq!(moving.solid().element_states()[0], material.initial_state());
            assert_eq!(moving.velocity_mps(), &[[0.0; 3]; 4]);

            let replay = moving
                .advance(
                    &supports,
                    &force,
                    [0.0; 3],
                    at(dt),
                    IncrementSettings::default(),
                    motion_controls(),
                    cx,
                )
                .unwrap();
            assert_eq!(
                replay.motion.solid().element_states(),
                step.motion.solid().element_states()
            );
            assert_eq!(replay.motion.velocity_mps(), step.motion.velocity_mps());
            // Hold the accepted configuration while momentum decays to rest.
            // Stopping kinetic loss is numerical, with no second plastic clock.
            let retained = step.motion.solid().element_states().to_vec();
            let static_force = VOLUME * step.motion.solid().element_first_piola_pa()[0][0];
            let stop_force =
                static_force - effective_mass * (step.motion.velocity_mps()[1][0] / dt);
            let stopped = step
                .motion
                .advance(
                    &supports,
                    &[(3, stop_force)],
                    [0.0; 3],
                    at(dt),
                    IncrementSettings::default(),
                    motion_controls(),
                    cx,
                )
                .unwrap();
            assert_eq!(stopped.motion.solid().element_states(), retained);
            near(stopped.motion.velocity_mps()[1][0], 0.0);
            near(stopped.material_energy.plastic_dissipation_j, 0.0);
            near(stopped.material_energy.damage_release_j, 0.0);
            near(stopped.inertia.kinetic_change_j, -kinetic);
            near(stopped.inertia.numerical_loss_j, kinetic);
            near(stopped.balance_error_j, 0.0);
        }
        assert!((integrities[0] - integrities[1]).abs() > 1e-4);
    });
}

#[test]
fn moving_supports_include_original_consistent_mass_and_nonzero_incoming_velocity() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let material = model(bands[0].band_width_m(), 0.04);
        let materials = [&material];
        let solid = FiniteDuctileTetSolid::try_new(
            &geometry,
            &materials,
            &bands,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let old_v = [0.25, -0.125, 0.0625];
        let new_v = [0.125, 0.25, -0.25];
        let dt = 0.125;
        let moving = FiniteDuctileTetMotion::try_new(solid, &[RHO], &[old_v; 4], cx).unwrap();
        let supports: Vec<_> = (0..12).map(|dof| (dof, dt * new_v[dof % 3])).collect();
        let step = moving
            .advance(
                &supports,
                &[],
                [0.0; 3],
                at(dt),
                IncrementSettings::default(),
                motion_controls(),
                cx,
            )
            .unwrap();
        let mass = RHO * VOLUME;
        let dv = core::array::from_fn(|i| new_v[i] - old_v[i]);
        for ((dof, _), reaction) in supports.iter().zip(&step.equilibrium.reactions_n) {
            near(*reaction, 0.25 * mass * dv[dof % 3] / dt);
        }
        for velocity in step.motion.velocity_mps() {
            for axis in 0..3 {
                near(velocity[axis], new_v[axis]);
            }
        }
        near(
            step.material_energy.support_endpoint_work_j,
            mass * dot3(dv, new_v),
        );
        near(step.material_energy.dead_load_endpoint_work_j, 0.0);
        near(
            step.inertia.kinetic_before_j,
            0.5 * mass * dot3(old_v, old_v),
        );
        near(
            step.inertia.kinetic_after_j,
            0.5 * mass * dot3(new_v, new_v),
        );
        near(step.inertia.numerical_loss_j, 0.5 * mass * dot3(dv, dv));
        near(step.material_energy.internal_endpoint_work_j, 0.0);
        near(step.material_energy.plastic_dissipation_j, 0.0);
        near(step.material_energy.damage_release_j, 0.0);
        near(step.balance_error_j, 0.0);
        assert_eq!(moving.velocity_mps(), &[old_v; 4]);
    });
}

#[test]
fn kinetic_budget_damage_and_cancellation_refuse_the_complete_motion_with_refined_retry() {
    with_cx(|gate, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let material = model(bands[0].band_width_m(), 0.04);
        let materials = [&material];
        let fresh = || {
            FiniteDuctileTetSolid::try_new(
                &geometry,
                &materials,
                &bands,
                TetAssemblyBudget::standard(),
                cx,
            )
            .unwrap()
        };
        assert!(FiniteDuctileTetMotion::try_new(fresh(), &[RHO], &[[0.0; 3]; 3], cx).is_err());
        assert!(FiniteDuctileTetMotion::try_new(fresh(), &[RHO], &[[f64::NAN; 3]; 4], cx).is_err());
        let moving = FiniteDuctileTetMotion::try_new(fresh(), &[RHO], &[[0.0; 3]; 4], cx).unwrap();
        let states = moving.solid().element_states().to_vec();
        let dt = 0.125;
        let body = [RHO, 0.0, 0.0];
        let limited = FiniteDuctileMotionControls {
            maximum_numerical_kinetic_loss_j: 0.25 * RHO * VOLUME * dt * dt,
            ..motion_controls()
        };
        assert!(matches!(
            moving.advance(
                &[],
                &[],
                body,
                at(dt),
                IncrementSettings::default(),
                limited,
                cx
            ),
            Err(FiniteDuctileTetError::Energy {
                quantity: "backward-Euler numerical kinetic loss",
                ..
            })
        ));
        let refined = moving
            .advance(
                &[],
                &[],
                body,
                at(0.5 * dt),
                IncrementSettings::default(),
                limited,
                cx,
            )
            .unwrap();
        near(refined.motion.velocity_mps()[0][0], 0.5 * dt);
        assert!(refined.inertia.numerical_loss_j < limited.maximum_numerical_kinetic_loss_j);
        let supports = affine(&NODES, stretch(0.08), [0.0; 3], 0);
        assert!(
            moving
                .advance(
                    &supports,
                    &[],
                    [0.0; 3],
                    at(dt),
                    IncrementSettings::default(),
                    motion_controls(),
                    cx
                )
                .is_err()
        );
        assert!(
            moving
                .advance(
                    &[],
                    &[],
                    body,
                    at(0.0),
                    IncrementSettings::default(),
                    motion_controls(),
                    cx
                )
                .is_err()
        );
        assert_eq!(moving.solid().element_states(), states);
        assert_eq!(moving.solid().displacement_m(), &[[0.0; 3]; 4]);
        assert_eq!(moving.velocity_mps(), &[[0.0; 3]; 4]);
        gate.request();
        assert!(matches!(
            moving.advance(
                &[],
                &[],
                body,
                at(dt),
                IncrementSettings::default(),
                motion_controls(),
                cx
            ),
            Err(FiniteDuctileTetError::Solve(
                fs_solid::plastic3::SmallStrainTetError::Cancelled
            ))
        ));
        assert_eq!(moving.solid().element_states(), states);
        assert_eq!(moving.velocity_mps(), &[[0.0; 3]; 4]);
    });
}
