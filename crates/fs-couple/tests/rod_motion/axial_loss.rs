//! G0/G1/G3/G4: axial material heat on the original finite-rotation trajectory.
use super::*;
use fs_couple::rod_motion::{PreparedRodAxialLoss, PreparedRodLoads, RodAngularLosses};
use fs_material::section::axial::ConcentricAxialViscosity;
use fs_qty::DynViscosity;

fn source(eta: f64, band: (f64, f64)) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density", 1000.0, Density::DIMS),
        ("young_modulus", 1e4, Pressure::DIMS),
        ("shear_modulus", 4e3, Pressure::DIMS),
        ("axial_eta", eta, DynViscosity::DIMS),
    ] {
        let mut validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
        if name == "axial_eta" {
            validity = validity.with("omega", band.0, band.1);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic axial rod specimen".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic axial rod".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn section(card: &MaterialCard, temperature: f64, with_loss: bool) -> ConcentricShearSection {
    let mut properties = vec![
        ("density", Density::DIMS),
        ("young_modulus", Pressure::DIMS),
        ("shear_modulus", Pressure::DIMS),
    ];
    if with_loss {
        properties.push(("axial_eta", DynViscosity::DIMS));
    }
    let requirements: Vec<_> = properties
        .into_iter()
        .map(|(name, dims)| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::NonNegative)
                .unwrap()
        })
        .collect();
    let point = QueryPoint::new()
        .with("T", temperature)
        .unwrap()
        .with("omega", 10.0)
        .unwrap();
    let material = resolve_material_state_point(
        card,
        &point,
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
fn viscous(card: &MaterialCard, temperature: f64) -> ConcentricAxialViscosity {
    section(card, temperature, true)
        .elastic()
        .clone()
        .with_axial_viscosity("axial_eta")
        .unwrap()
}
fn model(card: &MaterialCard, free: bool, segments: usize) -> Arc<PreparedRodMotion> {
    let gate = CancelGate::new();
    let mut rod =
        SpatialMaterialRod::try_straight(0.8, vec![section(card, 300.0, false); segments]).unwrap();
    if !free {
        let supports: Vec<_> = (1..=segments)
            .map(|node| RodSupport {
                node,
                position_m: [None, Some(0.0), Some(0.0)],
                frame: Some(So3::identity()),
            })
            .collect();
        rod = rod.with_supports(&supports, &gate).unwrap();
    }
    Arc::new(
        rod.into_motion_model(
            if free {
                RodRootBoundary::Free
            } else {
                RodRootBoundary::Clamped
            },
            &gate,
        )
        .unwrap(),
    )
}
fn angular() -> RodAngularLosses<'static> {
    RodAngularLosses {
        bending: None,
        torsional: None,
    }
}
fn vector(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn rotate(r: So3, a: [f64; 3]) -> [f64; 3] {
    vector(r.rotate(Vec3::new(a[0], a[1], a[2])).unwrap())
}

#[test]
fn g1_axial_velocity_heat_and_support_impulse_match_independent_split_solution() {
    let gate = CancelGate::new();
    for eta in [800.0, 8e8] {
        let card = source(eta, (0.0, f64::MAX));
        let model = model(&card, false, 1);
        let loss = PreparedRodAxialLoss::new(&model, vec![viscous(&card, 300.0)], &gate).unwrap();
        let mass = model.inertia()[1].mass_kg();
        let k = model.source().sections()[0].elastic().axial_stiffness_n() / 0.8;
        let c = eta * core::f64::consts::PI * 0.05_f64.powi(2) / 0.8;
        near(loss.resistance_n_s_m()[0], c, 2e-14 * c);
        let (x0, v0, dt, force) = (0.008, 0.1, 0.001, 0.002);
        let decay = (-0.5 * dt * c / mass).exp();
        let middle = decay * v0 + 0.5 * dt * (force - k * x0) / mass;
        let x1 = x0 + dt * middle;
        let before_loss = middle + 0.5 * dt * (force - k * x1) / mass;
        let expected_heat =
            0.5 * mass * (1.0 - decay * decay) * (v0 * v0 + before_loss * before_loss);
        let mut pose = model.initial_pose();
        pose.positions_m[1][0] += x0;
        let state = RodMotion::new(
            model.clone(),
            pose,
            &[[0.0; 3], [v0, 0.0, 0.0]],
            &[[0.0; 3]; 2],
            0.0,
            &gate,
        )
        .unwrap();
        let loads = [
            RodNodalWrench::default(),
            RodNodalWrench {
                force_n: [force, 0.0, 0.0],
                torque_nm: [0.0; 3],
            },
        ];
        let (next, report) = state
            .step_axially_damped_candidate(dt, &loss, angular(), &loads, controls(dt), &gate)
            .unwrap();
        near(next.pose().positions_m[1][0], 0.8 + x1, 2e-14);
        near(
            next.velocities_m_s().unwrap()[1][0],
            decay * before_loss,
            2e-14,
        );
        let heat = report.axial_heat.as_ref().unwrap();
        near(heat.total_heat_j, expected_heat, 3e-15);
        assert_eq!(heat.segment_heat_j, vec![heat.total_heat_j]);
        assert!(report.torsional_heat.is_none() && report.bending_heat.is_none());
        for residual in report
            .linear_balance_residual_n_s
            .into_iter()
            .chain(report.angular_balance_residual_nm_s)
        {
            near(residual, 0.0, 3e-14);
        }
        let prepared = PreparedRodLoads::new(&loads, &[[0.0; 3]; 2], &[], &gate).unwrap();
        let (spatial, spatial_report) = state
            .step_spatial_axially_damped_candidate(
                dt,
                &loss,
                angular(),
                &prepared,
                controls(dt),
                &gate,
            )
            .unwrap();
        assert_eq!(spatial.pose(), next.pose());
        assert_eq!(spatial_report.axial_heat, report.axial_heat);
    }
}

#[test]
fn g1_material_axial_ringdown_and_heat_converge_at_second_order() {
    let gate = CancelGate::new();
    let card = source(800.0, (0.0, f64::MAX));
    let model = model(&card, false, 1);
    let loss = PreparedRodAxialLoss::new(&model, vec![viscous(&card, 300.0)], &gate).unwrap();
    let mass = model.inertia()[1].mass_kg();
    let k = model.source().sections()[0].elastic().axial_stiffness_n() / 0.8;
    let gamma = loss.resistance_n_s_m()[0] / (2.0 * mass);
    let wd = (k / mass - gamma * gamma).sqrt();
    let (x0, time) = (0.008, 0.4);
    let x = x0 * (-gamma * time).exp() * ((wd * time).cos() + gamma / wd * (wd * time).sin());
    let v = -x0 * (-gamma * time).exp() * (k / mass) / wd * (wd * time).sin();
    let expected_heat = 0.5 * k * (x0 * x0 - x * x) - 0.5 * mass * v * v;
    let mut errors = Vec::new();
    let mut heat_errors = Vec::new();
    for steps in [100, 200, 400] {
        let dt = time / steps as f64;
        let mut pose = model.initial_pose();
        pose.positions_m[1][0] += x0;
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
                .step_axially_damped_candidate(
                    dt,
                    &loss,
                    angular(),
                    &[RodNodalWrench::default(); 2],
                    controls(dt),
                    &gate,
                )
                .unwrap();
            heat += report.axial_heat.unwrap().total_heat_j;
            state = next;
        }
        errors.push((state.pose().positions_m[1][0] - 0.8 - x).abs());
        heat_errors.push((heat - expected_heat).abs());
    }
    for errors in [errors, heat_errors] {
        for i in 0..2 {
            assert!(
                (3.8..4.2).contains(&(errors[i] / errors[i + 1])),
                "{errors:?}"
            );
        }
    }
}

#[test]
fn g3_axial_loss_keeps_rigid_translation_and_is_objective_under_world_rotation() {
    let gate = CancelGate::new();
    let card = source(800.0, (0.0, f64::MAX));
    let model = model(&card, true, 2);
    let loss = PreparedRodAxialLoss::new(&model, vec![viscous(&card, 300.0); 2], &gate).unwrap();
    let velocity = [0.1, -0.2, 0.3];
    let state = RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[velocity; 3],
        &[[0.0; 3]; 3],
        0.0,
        &gate,
    )
    .unwrap();
    let (next, report) = state
        .step_axially_damped_candidate(
            0.001,
            &loss,
            angular(),
            &[RodNodalWrench::default(); 3],
            controls(0.001),
            &gate,
        )
        .unwrap();
    assert!(report.axial_heat.unwrap().total_heat_j < 1e-26);
    for (a, b) in next
        .linear_momenta_kg_m_s()
        .iter()
        .flatten()
        .zip(state.linear_momenta_kg_m_s().iter().flatten())
    {
        near(*a, *b, 2e-13);
    }

    let world = So3::exp(So3Tangent::new(Vec3::new(0.4, -0.3, 0.2))).unwrap();
    let mut pose = model.initial_pose();
    pose.positions_m[1][1] = 0.01;
    pose.frames[2] = So3::exp(So3Tangent::new(Vec3::new(0.1, 0.2, -0.15))).unwrap();
    let velocities = [[0.1, 0.2, -0.1], [-0.1, 0.0, 0.1], [0.2, -0.2, 0.1]];
    let spins = [[0.2, -0.1, 0.3], [-0.2, 0.3, 0.1], [0.1, -0.2, -0.1]];
    let mut rotated = pose.clone();
    for node in 0..3 {
        rotated.positions_m[node] = rotate(world, pose.positions_m[node]);
        rotated.frames[node] = world.compose(pose.frames[node]).unwrap();
    }
    let a = RodMotion::new(model.clone(), pose, &velocities, &spins, 0.0, &gate).unwrap();
    let b = RodMotion::new(
        model.clone(),
        rotated,
        &velocities.map(|v| rotate(world, v)),
        &spins,
        0.0,
        &gate,
    )
    .unwrap();
    let dt = 1e-5;
    let (a, ar) = a
        .step_axially_damped_candidate(
            dt,
            &loss,
            angular(),
            &[RodNodalWrench::default(); 3],
            controls(dt),
            &gate,
        )
        .unwrap();
    let (b, br) = b
        .step_axially_damped_candidate(
            dt,
            &loss,
            angular(),
            &[RodNodalWrench::default(); 3],
            controls(dt),
            &gate,
        )
        .unwrap();
    near(
        ar.axial_heat.unwrap().total_heat_j,
        br.axial_heat.unwrap().total_heat_j,
        2e-16,
    );
    for node in 0..3 {
        for (x, y) in rotate(world, a.linear_momenta_kg_m_s()[node])
            .into_iter()
            .zip(b.linear_momenta_kg_m_s()[node])
        {
            near(x, y, 2e-13);
        }
    }
    for residual in ar
        .linear_balance_residual_n_s
        .into_iter()
        .chain(ar.angular_balance_residual_nm_s)
        .chain(br.linear_balance_residual_n_s)
        .chain(br.angular_balance_residual_nm_s)
    {
        near(residual, 0.0, 3e-13);
    }
}

#[test]
fn g0_extension_pullback_matches_finite_rotation_probes_and_annihilates_rigid_motion() {
    let gate = CancelGate::new();
    let card = source(800.0, (0.0, f64::MAX));
    let model = model(&card, true, 1);
    let mut pose = model.initial_pose();
    pose.positions_m = vec![[0.2, -0.1, 0.3], [0.95, 0.2, 0.1]];
    pose.frames[0] = So3::exp(So3Tangent::new(Vec3::new(0.2, -0.3, 0.1))).unwrap();
    pose.frames[1] = So3::exp(So3Tangent::new(Vec3::new(-0.1, 0.4, -0.2))).unwrap();
    let mut rows = [[[0.0; 6]; 2]; 1];
    model
        .axial_extension_gradients_into(&pose, &mut rows, &gate)
        .unwrap();
    let extension = |p: &fs_solid::rod::material_field::motion::RodMotionPose| {
        let half = So3::exp(p.frames[1].body_minus(p.frames[0]).unwrap().scale(0.5)).unwrap();
        let middle = p.frames[0].compose(half).unwrap();
        let delta = core::array::from_fn(|i| p.positions_m[1][i] - p.positions_m[0][i]);
        rotate(middle.inverse(), delta)[0] - 0.8
    };
    let h = 2e-6;
    for side in 0..2 {
        for component in 0..6 {
            let mut plus = pose.clone();
            let mut minus = pose.clone();
            if component < 3 {
                plus.positions_m[side][component] += h;
                minus.positions_m[side][component] -= h;
            } else {
                let mut step = [0.0; 3];
                step[component - 3] = h;
                let step = So3Tangent::new(Vec3::new(step[0], step[1], step[2]));
                plus.frames[side] = pose.frames[side].body_plus(step).unwrap();
                minus.frames[side] = pose.frames[side].body_plus(step.scale(-1.0)).unwrap();
            }
            near(
                (extension(&plus) - extension(&minus)) / (2.0 * h),
                rows[0][side][component],
                2e-9,
            );
        }
    }
    let omega = Vec3::new(0.3, -0.4, 0.2);
    let mut rate = 0.0;
    for side in 0..2 {
        let p = pose.positions_m[side];
        let v = vector(omega.cross(Vec3::new(p[0], p[1], p[2])));
        let w = rotate(pose.frames[side].inverse(), vector(omega));
        for c in 0..3 {
            rate += rows[0][side][c] * v[c] + rows[0][side][c + 3] * w[c];
        }
    }
    near(rate, 0.0, 2e-15);
}

#[test]
fn g0_g4_axial_zero_source_band_and_failed_step_preserve_native_state() {
    let gate = CancelGate::new();
    for eta in [0.0, 800.0] {
        let card = source(eta, (0.0, f64::MAX));
        let model = model(&card, false, 1);
        let loss = PreparedRodAxialLoss::new(&model, vec![viscous(&card, 300.0)], &gate).unwrap();
        assert!(PreparedRodAxialLoss::new(&model, vec![viscous(&card, 350.0)], &gate).is_err());
        assert!(PreparedRodAxialLoss::new(&model, vec![], &gate).is_err());
        let mut pose = model.initial_pose();
        pose.positions_m[1][0] += 0.008;
        let state = RodMotion::new(
            model,
            pose,
            &[[0.0; 3], [0.1, 0.0, 0.0]],
            &[[0.0; 3]; 2],
            0.0,
            &gate,
        )
        .unwrap();
        let dt = 0.001;
        let loads = [RodNodalWrench::default(); 2];
        let original_pose = state.pose().clone();
        let original_p = state.linear_momenta_kg_m_s().to_vec();
        let (a, ar) = state
            .step_axially_damped_candidate(dt, &loss, angular(), &loads, controls(dt), &gate)
            .unwrap();
        if eta == 0.0 {
            let (b, br) = state
                .step_candidate(dt, &loads, controls(dt), &gate)
                .unwrap();
            assert_eq!(a.pose(), b.pose());
            assert_eq!(a.linear_momenta_kg_m_s(), b.linear_momenta_kg_m_s());
            assert_eq!(ar.integration_defect_j, br.integration_defect_j);
            assert_eq!(ar.axial_heat.unwrap().total_heat_j, 0.0);
        }
        let mut strict = controls(dt);
        strict.energy_defect_limit_j = 0.0;
        assert!(
            state
                .step_axially_damped_candidate(dt, &loss, angular(), &loads, strict, &gate)
                .is_err()
        );
        let cancelled = CancelGate::new();
        cancelled.request();
        assert!(
            state
                .step_axially_damped_candidate(
                    dt,
                    &loss,
                    angular(),
                    &loads,
                    controls(dt),
                    &cancelled
                )
                .is_err()
        );
        assert_eq!(state.pose(), &original_pose);
        assert_eq!(state.linear_momenta_kg_m_s(), original_p);
    }
    let band = source(800.0, (1.0, 100.0));
    let model = model(&band, false, 1);
    assert!(PreparedRodAxialLoss::new(&model, vec![viscous(&band, 300.0)], &gate).is_err());
}
