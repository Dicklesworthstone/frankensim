//! Native rod Maxwell history: real material sources and independent dynamics.
#![cfg(feature = "rod-dynamics")]

use std::sync::Arc;

use fs_couple::rod_motion::{
    PreparedRodAxialRelaxation, RodMotion, RodMotionControls, RodMotionError, RodNodalWrench,
};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_ga::{So3, So3Tangent, Vec3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::relaxation::{
    LayerRelaxationProperties, axial::ConcentricPronyAxialSection,
};
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::{Density, Pressure, Time};
use fs_solid::rod::material_field::SpatialMaterialRod;
use fs_solid::rod::material_field::motion::{PreparedRodMotion, RodRootBoundary};
use fs_solid::rod::material_field::supports::RodSupport;

const LENGTH: f64 = 0.8;

fn source(
    arms: &[(f64, f64)],
    band: (f64, f64),
) -> (ConcentricShearSection, ConcentricPronyAxialSection) {
    let mut properties = vec![
        ("density".to_string(), 1000.0, Density::DIMS),
        ("young_modulus".to_string(), 1e4, Pressure::DIMS),
        ("shear_modulus".to_string(), 4e3, Pressure::DIMS),
    ];
    for (index, &(delta, tau)) in arms.iter().enumerate() {
        properties.push((format!("delta_{index}"), delta, Pressure::DIMS));
        properties.push((format!("tau_{index}"), tau, Time::DIMS));
    }
    let mut claims = ClaimSet::new();
    for (name, value, dims) in &properties {
        let mut validity = ValidityDomain::unconstrained();
        if name == "young_modulus" {
            validity = validity.with("omega", 0.0, f64::MAX);
        } else if name.starts_with("delta_") || name.starts_with("tau_") {
            validity = validity.with("omega", band.0, band.1);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, *dims),
                value: PropertyValue::Scalar {
                    value: *value,
                    dims: *dims,
                },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured axial Maxwell rod".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic rod".into(),
            phase: "solid".into(),
            process: "manufactured spectrum".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let resolve = |enriched: bool| {
        let requirements: Vec<_> = properties[..if enriched { properties.len() } else { 3 }]
            .iter()
            .map(|(name, _, dims)| {
                ScalarPropertyRequirement::try_new(name, *dims, ScalarAdmissibility::NonNegative)
                    .unwrap()
            })
            .collect();
        let material = resolve_material_state_point(
            &card,
            &QueryPoint::new().with("omega", 10.0).unwrap(),
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
    };
    let mechanical = resolve(false)
        .with_shear_modulus("shear_modulus", 0.8)
        .unwrap();
    let hereditary = resolve(true)
        .with_prony_axial(&[LayerRelaxationProperties {
            equilibrium_modulus: "young_modulus".into(),
            branches: (0..arms.len())
                .map(|i| (format!("delta_{i}"), format!("tau_{i}")))
                .collect(),
        }])
        .unwrap();
    (mechanical, hereditary)
}

fn specimen(
    arms: &[(f64, f64)],
    fixed_tip: bool,
    free: bool,
) -> (Arc<PreparedRodMotion>, Arc<PreparedRodAxialRelaxation>) {
    let gate = CancelGate::new();
    let (mechanical, hereditary) = source(arms, (0.0, f64::MAX));
    let mut rod = SpatialMaterialRod::try_straight(LENGTH, vec![mechanical]).unwrap();
    if !free {
        rod = rod
            .with_supports(
                &[RodSupport {
                    node: 1,
                    position_m: [fixed_tip.then_some(LENGTH), Some(0.0), Some(0.0)],
                    frame: Some(So3::identity()),
                }],
                &gate,
            )
            .unwrap();
    }
    let model = Arc::new(
        rod.into_motion_model(
            if free {
                RodRootBoundary::Free
            } else {
                RodRootBoundary::Clamped
            },
            &gate,
        )
        .unwrap(),
    );
    let law = Arc::new(PreparedRodAxialRelaxation::new(&model, vec![hereditary], &gate).unwrap());
    (model, law)
}

fn controls(dt: f64) -> RodMotionControls {
    RodMotionControls {
        maximum_step_s: dt,
        maximum_steps: 100_000,
        energy_defect_limit_j: 1e-4,
        maximum_rotation_rad: 0.25,
        balance_origin_m: [0.0; 3],
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn held_rod_relaxes_each_original_arm_with_exact_positive_heat() {
    let arms = [(5000.0, 0.08), (2500.0, 0.003)];
    let initial = [0.03, -0.02];
    let gate = CancelGate::new();
    let (model, law) = specimen(&arms, true, false);
    let mut state = RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[[0.0; 3]; 2],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap()
    .with_axial_relaxation(law.clone(), &initial, &gate)
    .unwrap();
    let storage0 = state.axial_relaxation_state().unwrap().stored_energy_j();
    let mut heat = [0.0; 2];
    let mut tip_reaction = 0.0;
    for step in 1..=20 {
        let (next, record) = state
            .step_candidate(0.01, &[RodNodalWrench::default(); 2], controls(0.01), &gate)
            .unwrap();
        tip_reaction += record.support_impulses[1].linear_n_s[0];
        let loss = record.axial_relaxation_heat.unwrap();
        assert!(loss.total_heat_j >= 0.0);
        near(loss.total_heat_j, loss.segment_heat_j.iter().sum(), 2e-16);
        for i in 0..2 {
            heat[i] += loss.branch_heat_j[i];
            let decay = fs_math::det::exp(-(step as f64 * 0.01) / arms[i].1);
            near(
                next.axial_relaxation_state().unwrap().elastic_strains()[i],
                initial[i] * decay,
                2e-15,
            );
            let scale = law.sections()[0].branches()[i].axial_stiffness_n * LENGTH;
            near(
                heat[i],
                0.5 * scale * initial[i] * initial[i] * (1.0 - decay * decay),
                2e-15,
            );
        }
        near(record.integration_defect_j, 0.0, 2e-15);
        assert_eq!(next.pose(), &model.initial_pose());
        state = next;
    }
    let memory = state.axial_relaxation_state().unwrap();
    near(
        memory.stored_energy_j() + memory.accumulated_heat_j(),
        storage0,
        3e-15,
    );
    let expected_reaction: f64 = arms
        .iter()
        .enumerate()
        .map(|(i, &(_, tau))| {
            law.sections()[0].branches()[i].axial_stiffness_n
                * initial[i]
                * tau
                * (1.0 - fs_math::det::exp(-0.2 / tau))
        })
        .sum();
    near(tip_reaction, expected_reaction, 2e-15);
    assert_eq!(state.linear_momenta_kg_m_s(), &[[0.0; 3]; 2]);
}

fn independent_bar_solution(
    mut x: [f64; 3],
    mass: f64,
    equilibrium_k: f64,
    excess_ea: f64,
    tau: f64,
    duration: f64,
) -> [f64; 3] {
    // Independent scalar continuum reduction: [extension, velocity, arm strain].
    // No rod/material stepping implementation participates in this RK4 oracle.
    let rhs = |x: [f64; 3]| {
        [
            x[1],
            -(equilibrium_k * x[0] + excess_ea * x[2]) / mass,
            x[1] / LENGTH - x[2] / tau,
        ]
    };
    let h = duration / 20_000.0;
    for _ in 0..20_000 {
        let k1 = rhs(x);
        let k2 = rhs(std::array::from_fn(|i| x[i] + 0.5 * h * k1[i]));
        let k3 = rhs(std::array::from_fn(|i| x[i] + 0.5 * h * k2[i]));
        let k4 = rhs(std::array::from_fn(|i| x[i] + h * k3[i]));
        x = std::array::from_fn(|i| x[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    x
}

#[test]
fn native_free_recovery_matches_standard_linear_solid_and_second_order_refinement() {
    let gate = CancelGate::new();
    let (model, law) = specimen(&[(5000.0, 0.08)], false, false);
    let k = model.source().sections()[0].elastic().axial_stiffness_n() / LENGTH;
    let mass = model.inertia()[1].mass_kg();
    let excess = law.sections()[0].branches()[0].axial_stiffness_n;
    let reference = independent_bar_solution([0.02, 0.0, 0.025], mass, k, excess, 0.08, 0.4);
    let mut errors = Vec::new();
    for steps in [100, 200, 400] {
        let dt = 0.4 / f64::from(steps);
        let mut pose = model.initial_pose();
        pose.positions_m[1][0] += 0.02;
        let mut state = RodMotion::new(
            model.clone(),
            pose,
            &[[0.0; 3]; 2],
            &[[0.0; 3]; 2],
            0.0,
            &gate,
        )
        .unwrap()
        .with_axial_relaxation(law.clone(), &[0.025], &gate)
        .unwrap();
        let energy0 = state.observe([0.0; 3]).unwrap().potential_energy_j;
        let mut numerical = 0.0;
        for _ in 0..steps {
            let (next, record) = state
                .step_candidate(dt, &[RodNodalWrench::default(); 2], controls(dt), &gate)
                .unwrap();
            assert!(record.axial_relaxation_heat.as_ref().unwrap().total_heat_j >= 0.0);
            numerical += record.integration_defect_j;
            state = next;
        }
        let x = [
            state.pose().positions_m[1][0] - LENGTH,
            state.velocities_m_s().unwrap()[1][0],
            state.axial_relaxation_state().unwrap().elastic_strains()[0],
        ];
        errors.push(
            x.iter()
                .zip(reference)
                .map(|(a, b)| (a - b) * (a - b))
                .sum::<f64>()
                .sqrt(),
        );
        let end = state.observe([0.0; 3]).unwrap();
        near(
            end.kinetic_energy_j
                + end.potential_energy_j
                + state.axial_relaxation_state().unwrap().accumulated_heat_j()
                - energy0,
            numerical,
            2e-14,
        );
    }
    assert!(
        errors[1] < 0.3 * errors[0] && errors[2] < 0.3 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 2e-6, "{errors:?}");
}

#[test]
fn finite_rotation_history_is_objective_and_keeps_internal_momentum_balance() {
    let gate = CancelGate::new();
    let (model, law) = specimen(&[(5000.0, 0.08)], false, true);
    let mut pose = model.initial_pose();
    pose.positions_m[1] = [0.81, 0.04, -0.02];
    pose.frames[1] = So3::exp(So3Tangent::new(Vec3::new(0.15, -0.2, 0.3))).unwrap();
    let velocity = [[0.02, -0.01, 0.03], [-0.01, 0.03, -0.02]];
    let omega = [[0.04, -0.02, 0.03], [-0.02, 0.01, 0.05]];
    let mut a = RodMotion::new(model.clone(), pose.clone(), &velocity, &omega, 0.0, &gate)
        .unwrap()
        .with_axial_relaxation(law.clone(), &[0.012], &gate)
        .unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(-0.4, 0.7, 0.2))).unwrap();
    let rotate = |x: [f64; 3]| {
        let q = rotation.rotate(Vec3::new(x[0], x[1], x[2])).unwrap();
        [q.x, q.y, q.z]
    };
    for node in 0..2 {
        let rotated = rotate(pose.positions_m[node]);
        pose.positions_m[node] = std::array::from_fn(|i| rotated[i] + [0.7, -0.3, 0.2][i]);
        pose.frames[node] = rotation.compose(pose.frames[node]).unwrap();
    }
    let mut b = RodMotion::new(model, pose, &velocity.map(rotate), &omega, 0.0, &gate)
        .unwrap()
        .with_axial_relaxation(law, &[0.012], &gate)
        .unwrap();
    for _ in 0..20 {
        let (na, ra) = a
            .step_candidate(1e-4, &[RodNodalWrench::default(); 2], controls(1e-4), &gate)
            .unwrap();
        let (nb, rb) = b
            .step_candidate(1e-4, &[RodNodalWrench::default(); 2], controls(1e-4), &gate)
            .unwrap();
        near(
            na.axial_relaxation_state().unwrap().elastic_strains()[0],
            nb.axial_relaxation_state().unwrap().elastic_strains()[0],
            2e-14,
        );
        near(
            ra.axial_relaxation_heat.unwrap().total_heat_j,
            rb.axial_relaxation_heat.unwrap().total_heat_j,
            2e-15,
        );
        for value in ra
            .linear_balance_residual_n_s
            .into_iter()
            .chain(ra.angular_balance_residual_nm_s)
        {
            assert!(value.abs() < 2e-12, "{value:e}");
        }
        a = na;
        b = nb;
    }
}

#[test]
fn relaxed_rigid_translation_and_empty_spectra_do_not_create_heat() {
    let gate = CancelGate::new();
    for arms in [vec![], vec![(5000.0, 0.08)]] {
        let (model, law) = specimen(&arms, false, true);
        let mut state = RodMotion::new(
            model.clone(),
            model.initial_pose(),
            &[[0.0, 0.2, 0.0]; 2],
            &[[0.0; 3]; 2],
            0.0,
            &gate,
        )
        .unwrap()
        .with_relaxed_axial_relaxation(law, &gate)
        .unwrap();
        for _ in 0..5 {
            let (next, record) = state
                .step_candidate(1e-3, &[RodNodalWrench::default(); 2], controls(1e-3), &gate)
                .unwrap();
            assert_eq!(record.axial_relaxation_heat.unwrap().total_heat_j, 0.0);
            assert_eq!(
                next.axial_relaxation_state().unwrap().stored_energy_j(),
                0.0
            );
            state = next;
        }
    }
}

#[test]
fn failed_steps_and_cloned_continuations_preserve_material_history() {
    let gate = CancelGate::new();
    let (model, law) = specimen(&[(5000.0, 0.08)], false, false);
    let mut pose = model.initial_pose();
    pose.positions_m[1][0] += 0.02;
    let state = RodMotion::new(model, pose, &[[0.0; 3]; 2], &[[0.0; 3]; 2], 0.0, &gate)
        .unwrap()
        .with_axial_relaxation(law.clone(), &[0.025], &gate)
        .unwrap();
    let mut tight = controls(0.004);
    tight.energy_defect_limit_j = 0.0;
    assert!(matches!(
        state.step_candidate(0.004, &[RodNodalWrench::default(); 2], tight, &gate),
        Err(RodMotionError::Energy { .. })
    ));
    assert_eq!(
        state.axial_relaxation_state().unwrap().elastic_strains(),
        &[0.025]
    );
    assert_eq!(
        state.axial_relaxation_state().unwrap().accumulated_heat_j(),
        0.0
    );
    let (a, ra) = state
        .step_candidate(
            0.001,
            &[RodNodalWrench::default(); 2],
            controls(0.001),
            &gate,
        )
        .unwrap();
    let (b, rb) = state
        .clone()
        .step_candidate(
            0.001,
            &[RodNodalWrench::default(); 2],
            controls(0.001),
            &gate,
        )
        .unwrap();
    assert_eq!(a.pose(), b.pose());
    assert_eq!(
        a.axial_relaxation_state().unwrap().elastic_strains(),
        b.axial_relaxation_state().unwrap().elastic_strains()
    );
    assert_eq!(ra.axial_relaxation_heat, rb.axial_relaxation_heat);
    assert!(a.with_relaxed_axial_relaxation(law, &gate).is_err());
}

#[test]
fn finite_band_spectra_cannot_authorize_arbitrary_native_transients() {
    let gate = CancelGate::new();
    let (section, law) = source(&[(5000.0, 0.08)], (1.0, 1000.0));
    let model = SpatialMaterialRod::try_straight(LENGTH, vec![section])
        .unwrap()
        .into_motion_model(RodRootBoundary::Free, &gate)
        .unwrap();
    assert!(PreparedRodAxialRelaxation::new(&model, vec![law], &gate).is_err());
}

#[test]
fn fast_fixed_stiffness_arms_do_not_leave_artificial_step_drag() {
    let gate = CancelGate::new();
    let (model, law) = specimen(&[(5000.0, 1e-8)], false, false);
    let initial = RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[[0.0; 3], [0.2, 0.0, 0.0]],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap();
    let mut elastic = initial.clone();
    let mut hereditary = initial.with_relaxed_axial_relaxation(law, &gate).unwrap();
    for _ in 0..20 {
        elastic = elastic
            .step_candidate(0.01, &[RodNodalWrench::default(); 2], controls(0.01), &gate)
            .unwrap()
            .0;
        hereditary = hereditary
            .step_candidate(0.01, &[RodNodalWrench::default(); 2], controls(0.01), &gate)
            .unwrap()
            .0;
    }
    // Here dt/tau=10^6. Physical excess viscosity is delta(EA)*tau/L,
    // so both drag and dissipated heat vanish with tau at fixed excess stiffness.
    // A split explicit branch-force kick leaves a finite O(k*dt) drag instead.
    near(
        hereditary.pose().positions_m[1][0],
        elastic.pose().positions_m[1][0],
        2e-8,
    );
    near(
        hereditary.velocities_m_s().unwrap()[1][0],
        elastic.velocities_m_s().unwrap()[1][0],
        4e-8,
    );
    let memory = hereditary.axial_relaxation_state().unwrap();
    assert!(memory.accumulated_heat_j() > 0.0 && memory.accumulated_heat_j() < 4e-9);
    assert!(memory.elastic_strains()[0].abs() < 3e-9);
}

fn independent_kelvin_solution(
    mut x: [f64; 2],
    mass: f64,
    stiffness: f64,
    resistance: f64,
    duration: f64,
) -> [f64; 2] {
    let rhs = |x: [f64; 2]| [x[1], -(stiffness * x[0] + resistance * x[1]) / mass];
    let h = duration / 20_000.0;
    for _ in 0..20_000 {
        let k1 = rhs(x);
        let k2 = rhs(std::array::from_fn(|i| x[i] + 0.5 * h * k1[i]));
        let k3 = rhs(std::array::from_fn(|i| x[i] + 0.5 * h * k2[i]));
        let k4 = rhs(std::array::from_fn(|i| x[i] + h * k3[i]));
        x = std::array::from_fn(|i| x[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    x
}

#[test]
fn unresolved_fast_arm_at_fixed_viscosity_reaches_the_kelvin_limit() {
    let gate = CancelGate::new();
    let tau = 1e-9;
    let (model, law) = specimen(&[(100.0 / tau, tau)], false, false);
    let mass = model.inertia()[1].mass_kg();
    let stiffness = model.source().sections()[0].elastic().axial_stiffness_n() / LENGTH;
    let resistance = law.sections()[0].branches()[0].axial_stiffness_n * tau / LENGTH;
    let reference = independent_kelvin_solution([0.02, 0.1], mass, stiffness, resistance, 0.1);
    let mut pose = model.initial_pose();
    pose.positions_m[1][0] += 0.02;
    let mut state = RodMotion::new(
        model,
        pose,
        &[[0.0; 3], [0.1, 0.0, 0.0]],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap()
    .with_axial_relaxation(law, &[tau * 0.1 / LENGTH], &gate)
    .unwrap();
    for _ in 0..200 {
        let (next, record) = state
            .step_candidate(
                0.0005,
                &[RodNodalWrench::default(); 2],
                controls(0.0005),
                &gate,
            )
            .unwrap();
        assert!(record.axial_relaxation_heat.as_ref().unwrap().total_heat_j >= 0.0);
        state = next;
    }
    near(state.pose().positions_m[1][0] - LENGTH, reference[0], 2e-7);
    near(state.velocities_m_s().unwrap()[1][0], reference[1], 2e-7);
    assert!(state.axial_relaxation_state().unwrap().accumulated_heat_j() > 0.0);
}

#[test]
fn strong_reciprocal_flow_preserves_small_native_momentum_tails() {
    let gate = CancelGate::new();
    let tau = 1e-3;
    let (model, law) = specimen(&[(3.2e7, tau)], false, false);
    let mass = model.inertia()[1].mass_kg();
    let rigidity = law.sections()[0].branches()[0].axial_stiffness_n;
    let omega_squared = rigidity / (mass * LENGTH);
    let rate = 1.0 / tau;
    // Independent stable roots for r''+r'/tau+mu*k*r=0, r'(0)=0.
    let slow =
        -2.0 * omega_squared / (rate + fs_math::det::sqrt(rate * rate - 4.0 * omega_squared));
    let fast = -rate - slow;
    let slow_weight = -fast / (slow - fast);
    let fast_weight = slow / (slow - fast);
    for sign in [-1.0, 1.0] {
        let mut state = RodMotion::new(
            model.clone(),
            model.initial_pose(),
            &[[0.0; 3], [sign, 0.0, 0.0]],
            &[[0.0; 3]; 2],
            0.0,
            &gate,
        )
        .unwrap()
        .with_relaxed_axial_relaxation(law.clone(), &gate)
        .unwrap();
        for step in 1..=2 {
            let (next, _) = state
                .step_candidate(1.0, &[RodNodalWrench::default(); 2], controls(1.0), &gate)
                .unwrap();
            // The first exact half-flow leaves motion below the position ulp.
            // The native elastic A stage therefore sees its unchanged reference
            // geometry: this fixture isolates repeated frozen-pose flow tails.
            assert_eq!(next.pose(), &model.initial_pose());
            let t = f64::from(step);
            let a = slow_weight * fs_math::det::exp(slow * t);
            let b = fast_weight * fs_math::det::exp(fast * t);
            let expected = sign * (a + b);
            let strain = -sign * (slow * a + fast * b) * mass / rigidity;
            let actual = next.velocities_m_s().unwrap()[1][0];
            assert!(
                actual * sign > 0.0,
                "surviving rate lost its physical sign: {actual:e}"
            );
            near(actual / expected, 1.0, 2e-8);
            near(
                next.axial_relaxation_state().unwrap().elastic_strains()[0] / strain,
                1.0,
                2e-8,
            );
            state = next;
        }
    }
}
