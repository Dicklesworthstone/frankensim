//! G0/G1/G3/G4: source-driven transverse shear loss on native rod motion.
use super::*;
use fs_couple::rod_motion::{
    PreparedRodLoads, PreparedRodTransverseLoss, RodAngularLosses, RodMaterialLosses,
};
use fs_qty::DynViscosity;

fn viscous_section(
    eta: f64,
    band: Option<(f64, f64)>,
    temperature: f64,
    correction: f64,
) -> ConcentricShearSection {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [
        ("density", 1000.0, Density::DIMS),
        ("young_modulus", 1e4, Pressure::DIMS),
        ("shear_modulus", 4e3, Pressure::DIMS),
        ("shear_eta", eta, DynViscosity::DIMS),
    ] {
        let mut validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
        if name == "shear_eta"
            && let Some((lo, hi)) = band
        {
            validity = validity.with("omega", lo, hi);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured transverse shear viscosity".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured shear rod".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new()
        .with("T", temperature)
        .unwrap()
        .with("omega", 10.0)
        .unwrap();
    let material = resolve_material_state_point(
        &card,
        &query,
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
    .with_shear_modulus("shear_modulus", correction)
    .unwrap()
}

fn specimen(
    eta: f64,
    free: bool,
    segments: usize,
) -> (Arc<PreparedRodMotion>, PreparedRodTransverseLoss) {
    let gate = CancelGate::new();
    let section = viscous_section(eta, Some((0.0, f64::MAX)), 300.0, 0.8);
    let mut rod = SpatialMaterialRod::try_straight(0.8, vec![section.clone(); segments]).unwrap();
    if !free {
        let supports: Vec<_> = (1..=segments)
            .map(|node| RodSupport {
                node,
                position_m: [Some(0.8 * node as f64 / segments as f64), None, None],
                frame: Some(So3::identity()),
            })
            .collect();
        rod = rod.with_supports(&supports, &gate).unwrap();
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
    let loss = PreparedRodTransverseLoss::new(
        &model,
        vec![section.with_transverse_viscosity("shear_eta").unwrap(); segments],
        &gate,
    )
    .unwrap();
    (model, loss)
}
fn losses(loss: &PreparedRodTransverseLoss) -> RodMaterialLosses<'_> {
    RodMaterialLosses {
        angular: RodAngularLosses {
            bending: None,
            torsional: None,
        },
        axial: None,
        transverse: Some(loss),
    }
}

#[test]
fn g1_nonuniform_reference_intervals_keep_their_actual_transverse_resistances() {
    let section = viscous_section(800.0, Some((0.0, f64::MAX)), 300.0, 0.8);
    let lengths = vec![0.125, 0.375, 0.5];
    let gate = CancelGate::new();
    let model = SpatialMaterialRod::try_nonuniform_straight(lengths.clone(), vec![section.clone(); 3]).unwrap()
        .into_motion_model(RodRootBoundary::Free, &gate).unwrap();
    let loss = PreparedRodTransverseLoss::new(&model,
        vec![section.with_transverse_viscosity("shear_eta").unwrap(); 3], &gate).unwrap();
    for (&actual, length) in loss.resistance_n_s_m().iter().zip(lengths) {
        let expected = 0.8 * 800.0 * core::f64::consts::PI * 0.05_f64.powi(2) / length;
        near(actual, expected, 3e-14 * expected);
    }
}

#[test]
fn g1_transverse_heat_and_supports_match_independent_split_oscillators() {
    let gate = CancelGate::new();
    for eta in [800.0, 8e8] {
        let (model, loss) = specimen(eta, false, 1);
        let mass = model.inertia()[1].mass_kg();
        let k = model.source().sections()[0].shear_stiffness_n() / 0.8;
        let c = 0.8 * eta * core::f64::consts::PI * 0.05_f64.powi(2) / 0.8;
        near(loss.resistance_n_s_m()[0], c, 2e-14 * c);
        let (dt, x0, v0, force) = (0.001, [0.008, -0.005], [0.1, -0.07], [0.002, -0.003]);
        let decay = (-0.5 * dt * c / mass).exp();
        let mut pose = model.initial_pose();
        pose.positions_m[1][1..].copy_from_slice(&x0);
        let state = RodMotion::new(
            model.clone(),
            pose,
            &[[0.0; 3], [0.0, v0[0], v0[1]]],
            &[[0.0; 3]; 2],
            0.0,
            &gate,
        )
        .unwrap();
        let loads = [
            RodNodalWrench::default(),
            RodNodalWrench {
                force_n: [0.0, force[0], force[1]],
                torque_nm: [0.0; 3],
            },
        ];
        let (next, report) = state
            .step_material_damped_candidate(dt, losses(&loss), &loads, controls(dt), &gate)
            .unwrap();
        let mut expected_heat = 0.0;
        for axis in 0..2 {
            let middle = decay * v0[axis] + 0.5 * dt * (force[axis] - k * x0[axis]) / mass;
            let x1 = x0[axis] + dt * middle;
            let before_loss = middle + 0.5 * dt * (force[axis] - k * x1) / mass;
            expected_heat += 0.5
                * mass
                * (1.0 - decay * decay)
                * (v0[axis] * v0[axis] + before_loss * before_loss);
            near(next.pose().positions_m[1][axis + 1], x1, 3e-14);
            near(
                next.velocities_m_s().unwrap()[1][axis + 1],
                decay * before_loss,
                3e-14,
            );
        }
        let heat = report.transverse_heat.as_ref().unwrap();
        near(heat.total_heat_j, expected_heat, 4e-15);
        assert_eq!(heat.segment_heat_j, vec![heat.total_heat_j]);
        assert!(
            report.axial_heat.is_none()
                && report.bending_heat.is_none()
                && report.torsional_heat.is_none()
        );
        for residual in report
            .linear_balance_residual_n_s
            .into_iter()
            .chain(report.angular_balance_residual_nm_s)
        {
            near(residual, 0.0, 4e-14);
        }
        let prepared = PreparedRodLoads::new(&loads, &[[0.0; 3]; 2], &[], &gate).unwrap();
        let (spatial, receipt) = state
            .step_spatial_material_damped_candidate(
                dt,
                losses(&loss),
                &prepared,
                controls(dt),
                &gate,
            )
            .unwrap();
        assert_eq!(spatial.pose(), next.pose());
        assert_eq!(receipt.transverse_heat, report.transverse_heat);
    }
}

#[test]
fn g1_transverse_ringdown_and_independent_heat_converge_at_second_order() {
    let gate = CancelGate::new();
    let (model, loss) = specimen(800.0, false, 1);
    let mass = model.inertia()[1].mass_kg();
    let k = model.source().sections()[0].shear_stiffness_n() / 0.8;
    let gamma = loss.resistance_n_s_m()[0] / (2.0 * mass);
    let wd = (k / mass - gamma * gamma).sqrt();
    let (x0, time) = (0.008, 0.4);
    let x = x0 * (-gamma * time).exp() * ((wd * time).cos() + gamma / wd * (wd * time).sin());
    let v = -x0 * (-gamma * time).exp() * (k / mass) / wd * (wd * time).sin();
    let exact_heat = 0.5 * k * (x0 * x0 - x * x) - 0.5 * mass * v * v;
    let mut errors = Vec::new();
    let mut heat_errors = Vec::new();
    for steps in [100, 200, 400] {
        let dt = time / steps as f64;
        let mut pose = model.initial_pose();
        pose.positions_m[1][1] = x0;
        let mut state = RodMotion::new(
            model.clone(),
            pose,
            &[[0.0; 3]; 2],
            &[[0.0; 3]; 2],
            0.0,
            &gate,
        )
        .unwrap();
        let mut heat = 0.0;
        for _ in 0..steps {
            let (next, report) = state
                .step_material_damped_candidate(
                    dt,
                    losses(&loss),
                    &[RodNodalWrench::default(); 2],
                    controls(dt),
                    &gate,
                )
                .unwrap();
            assert!(report.transverse_heat.as_ref().unwrap().total_heat_j >= 0.0);
            heat += report.transverse_heat.unwrap().total_heat_j;
            state = next;
        }
        errors.push((state.pose().positions_m[1][1] - x).abs());
        heat_errors.push((heat - exact_heat).abs());
    }
    for values in [errors, heat_errors] {
        for i in 0..2 {
            assert!(
                (3.8..4.2).contains(&(values[i] / values[i + 1])),
                "{values:?}"
            );
        }
    }
}

#[test]
fn g3_zero_source_keeps_native_motion_and_common_translation_is_lossless() {
    let gate = CancelGate::new();
    let (model, zero) = specimen(0.0, true, 2);
    let mut pose = model.initial_pose();
    pose.positions_m[1][1] = 0.003;
    let state = RodMotion::new(
        model,
        pose,
        &[[0.12, -0.03, 0.02]; 3],
        &[[0.0; 3]; 3],
        0.0,
        &gate,
    )
    .unwrap();
    let loads = [RodNodalWrench::default(); 3];
    let (plain, _) = state
        .step_candidate(1e-4, &loads, controls(1e-4), &gate)
        .unwrap();
    let (zeroed, report) = state
        .step_material_damped_candidate(1e-4, losses(&zero), &loads, controls(1e-4), &gate)
        .unwrap();
    assert_eq!(plain.pose(), zeroed.pose());
    assert_eq!(
        plain.velocities_m_s().unwrap(),
        zeroed.velocities_m_s().unwrap()
    );
    assert_eq!(report.transverse_heat.unwrap().total_heat_j, 0.0);
    let (model, loss) = specimen(800.0, true, 2);
    let state = RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[[0.12, -0.03, 0.02]; 3],
        &[[0.0; 3]; 3],
        0.0,
        &gate,
    )
    .unwrap();
    let (_, report) = state
        .step_material_damped_candidate(1e-4, losses(&loss), &loads, controls(1e-4), &gate)
        .unwrap();
    assert!(report.transverse_heat.unwrap().total_heat_j < 1e-29);
}

#[test]
fn g0_g4_shear_source_admission_and_cancelled_candidate_preserve_state() {
    let gate = CancelGate::new();
    let (model, loss) = specimen(800.0, true, 1);
    for section in [
        viscous_section(-1.0, Some((0.0, f64::MAX)), 300.0, 0.8),
        viscous_section(800.0, None, 300.0, 0.8),
    ] {
        assert!(section.with_transverse_viscosity("shear_eta").is_err());
    }
    let section = viscous_section(800.0, Some((0.0, f64::MAX)), 300.0, 0.8);
    for name in ["missing", "density"] {
        assert!(section.clone().with_transverse_viscosity(name).is_err());
    }
    for section in [
        viscous_section(800.0, Some((1.0, 100.0)), 300.0, 0.8),
        viscous_section(800.0, Some((0.0, f64::MAX)), 301.0, 0.8),
        viscous_section(800.0, Some((0.0, f64::MAX)), 300.0, 0.7),
    ] {
        assert!(
            PreparedRodTransverseLoss::new(
                &model,
                vec![section.with_transverse_viscosity("shear_eta").unwrap()],
                &gate
            )
            .is_err()
        );
    }
    let state = RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[[0.0; 3]; 2],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap();
    let before = state.pose().clone();
    gate.request();
    assert!(
        state
            .step_material_damped_candidate(
                1e-4,
                losses(&loss),
                &[RodNodalWrench::default(); 2],
                controls(1e-4),
                &gate
            )
            .is_err()
    );
    assert_eq!(*state.pose(), before);
    assert_eq!(state.time_s(), 0.0);
}
