//! Physical bending loss on the original finite-rotation rod trajectory.
use super::*;
use fs_couple::rod_motion::{
    PreparedRodBendingLoss, PreparedRodLoads, PreparedRodTorsionalLoss, RodAngularLosses,
    SpatialRodMotion,
};
use fs_material::section::ConcentricKelvinVoigtSection;
use fs_qty::DynViscosity;

fn card(eta: f64, band: (f64, f64)) -> MaterialCard {
    card_with_shear(eta, 0.3 * eta, band)
}

fn card_with_shear(eta: f64, shear_eta: f64, band: (f64, f64)) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density", 1000.0, Density::DIMS),
        ("young_modulus", 1e4, Pressure::DIMS),
        ("shear_modulus", 4e3, Pressure::DIMS),
        ("bending_eta", eta, DynViscosity::DIMS),
        ("shear_eta", shear_eta, DynViscosity::DIMS),
    ] {
        let mut validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
        if dims == DynViscosity::DIMS {
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
                    source: "synthetic independent bending and shear viscosities".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic bending specimen".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn sourced(card: &MaterialCard, temperature: f64, loss: bool) -> ConcentricShearSection {
    let mut properties = vec![
        (
            "density",
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        ),
        (
            "young_modulus",
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        ),
        (
            "shear_modulus",
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        ),
    ];
    if loss {
        properties.extend([
            (
                "bending_eta",
                DynViscosity::DIMS,
                ScalarAdmissibility::NonNegative,
            ),
            (
                "shear_eta",
                DynViscosity::DIMS,
                ScalarAdmissibility::NonNegative,
            ),
        ]);
    }
    let requirements: Vec<_> = properties
        .into_iter()
        .map(|(name, dims, domain)| ScalarPropertyRequirement::try_new(name, dims, domain).unwrap())
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

fn viscous(card: &MaterialCard, temperature: f64) -> ConcentricKelvinVoigtSection {
    sourced(card, temperature, true)
        .elastic()
        .clone()
        .with_kelvin_voigt_bending("bending_eta")
        .unwrap()
}

fn pinned_model(card: &MaterialCard, temperature: f64) -> Arc<PreparedRodMotion> {
    pinned_segments(card, temperature, 1)
}

fn pinned_segments(
    card: &MaterialCard,
    temperature: f64,
    segments: usize,
) -> Arc<PreparedRodMotion> {
    let gate = CancelGate::new();
    let supports: Vec<_> = (1..=segments)
        .map(|node| RodSupport {
            node,
            position_m: [
                Some(0.8 * node as f64 / segments as f64),
                Some(0.0),
                Some(0.0),
            ],
            frame: None,
        })
        .collect();
    Arc::new(
        SpatialMaterialRod::try_straight(0.8, vec![sourced(card, temperature, false); segments])
            .unwrap()
            .with_supports(&supports, &gate)
            .unwrap()
            .into_motion_model(RodRootBoundary::Clamped, &gate)
            .unwrap(),
    )
}

fn about(axis: usize, angle: f64) -> So3 {
    let mut vector = [0.0; 3];
    vector[axis] = angle;
    So3::exp(So3Tangent::new(Vec3::new(vector[0], vector[1], vector[2]))).unwrap()
}
fn vector(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
fn rotate(frame: So3, value: [f64; 3]) -> [f64; 3] {
    vector(
        frame
            .rotate(Vec3::new(value[0], value[1], value[2]))
            .unwrap(),
    )
}
fn vector_near(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    for c in 0..3 {
        near(a[c], b[c], tolerance);
    }
}

#[test]
fn g1_both_bending_planes_match_exact_pair_heat_and_nonlinear_elastic_motion() {
    let gate = CancelGate::new();
    for eta in [800.0, 8e8] {
        for axis in [1, 2] {
            let source = card(eta, (0.0, f64::MAX));
            let model = pinned_model(&source, 300.0);
            let loss =
                PreparedRodBendingLoss::new(&model, vec![viscous(&source, 300.0)], &gate).unwrap();
            let section = &model.source().sections()[0];
            let inertia = model.inertia()[1].principal_kg_m2()[axis];
            let ea = section.elastic().axial_stiffness_n();
            let ga = section.shear_stiffness_n();
            let kb = section.elastic().bending_stiffness_n_m2() / 0.8;
            let radius_squared = 0.05 * 0.05;
            let c = eta * core::f64::consts::PI * (radius_squared * radius_squared) / (4.0 * 0.8);
            near(loss.resistance_nm_s()[0], c, 2e-14 * c);
            // With both positions held, midpoint rotation contributes axial
            // and shear energy. Differentiate that energy independently.
            let restoring = |q: f64| {
                let (s, co) = (0.5 * q).sin_cos();
                kb * q + 0.4 * s * (ga * co - ea * (co - 1.0))
            };
            let (q0, w0, dt, torque) = (0.08, 0.3, 0.001, 0.002);
            let decay = (-0.5 * dt * c / inertia).exp();
            let middle = decay * w0 + 0.5 * dt * (torque - restoring(q0)) / inertia;
            let q1 = q0 + dt * middle;
            let wb = middle + 0.5 * dt * (torque - restoring(q1)) / inertia;
            let expected_heat = 0.5 * inertia * (1.0 - decay * decay) * (w0 * w0 + wb * wb);
            let mut pose = model.initial_pose();
            pose.frames[1] = about(axis, q0);
            let mut omega = [[0.0; 3]; 2];
            omega[1][axis] = w0;
            let state =
                RodMotion::new(model.clone(), pose, &[[0.0; 3]; 2], &omega, 0.0, &gate).unwrap();
            let mut loads = [RodNodalWrench::default(); 2];
            loads[1].torque_nm[axis] = torque;
            let selected = RodAngularLosses {
                bending: Some(&loss),
                torsional: None,
            };
            let (next, receipt) = state
                .step_damped_candidate(dt, selected, &loads, controls(dt), &gate)
                .unwrap();
            near(vector(next.pose().frames[1].log().angular)[axis], q1, 2e-13);
            near(
                next.angular_velocities_body_rad_s().unwrap()[1][axis],
                decay * wb,
                2e-12,
            );
            let heat = receipt.bending_heat.as_ref().unwrap();
            near(heat.total_heat_j, expected_heat, 2e-14);
            assert_eq!(heat.segment_heat_j, vec![heat.total_heat_j]);
            assert!(receipt.torsional_heat.is_none());
            vector_near(receipt.angular_balance_residual_nm_s, [0.0; 3], 2e-13);
            vector_near(receipt.linear_balance_residual_n_s, [0.0; 3], 2e-13);
            assert_eq!(next.pose().positions_m, model.initial_pose().positions_m);
            assert_eq!(next.pose().frames[0], So3::identity());
            let prepared = PreparedRodLoads::new(&loads, &[[0.0; 3]; 2], &[], &gate).unwrap();
            let (spatial, spatial_receipt) = state
                .step_spatial_damped_candidate(dt, selected, &prepared, controls(dt), &gate)
                .unwrap();
            assert_eq!(spatial.pose(), next.pose());
            assert_eq!(spatial_receipt.bending_heat, receipt.bending_heat);
        }
    }
}

#[test]
fn g1_small_bending_oscillator_and_dissipated_energy_refine_at_second_order() {
    let gate = CancelGate::new();
    let source = card(1600.0, (0.0, f64::MAX));
    let model = pinned_model(&source, 300.0);
    let loss = PreparedRodBendingLoss::new(&model, vec![viscous(&source, 300.0)], &gate).unwrap();
    let section = &model.source().sections()[0];
    let inertia = model.inertia()[1].principal_kg_m2()[2];
    // Linearization of the independently specified pinned-position energy.
    let k =
        section.elastic().bending_stiffness_n_m2() / 0.8 + section.shear_stiffness_n() * 0.8 / 4.0;
    let c = loss.resistance_nm_s()[0];
    let alpha = c / (2.0 * inertia);
    let frequency = (k / inertia - alpha * alpha).sqrt();
    let (q0, time) = (1e-6, 0.1);
    let (s, co) = (frequency * time).sin_cos();
    let exact_q = q0 * (-alpha * time).exp() * (co + alpha * s / frequency);
    let exact_w = -q0 * (-alpha * time).exp() * (k / inertia) * s / frequency;
    let exact_heat = 0.5 * k * (q0 * q0 - exact_q * exact_q) - 0.5 * inertia * exact_w * exact_w;
    let mut motion_errors = Vec::new();
    let mut heat_errors = Vec::new();
    for (dt, steps) in [(0.001, 100), (0.0005, 200), (0.00025, 400)] {
        let mut pose = model.initial_pose();
        pose.frames[1] = about(2, q0);
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
            let (next, receipt) = state
                .step_damped_candidate(
                    dt,
                    RodAngularLosses {
                        bending: Some(&loss),
                        torsional: None,
                    },
                    &[RodNodalWrench::default(); 2],
                    controls(dt),
                    &gate,
                )
                .unwrap();
            heat += receipt.bending_heat.unwrap().total_heat_j;
            state = next;
        }
        let q = state.pose().frames[1].log().angular.z;
        let w = state.angular_velocities_body_rad_s().unwrap()[1][2];
        motion_errors.push((q - exact_q).abs() + (w - exact_w).abs() / frequency);
        heat_errors.push((heat - exact_heat).abs());
    }
    for errors in [motion_errors, heat_errors] {
        for pair in errors.windows(2) {
            assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "{errors:?}");
        }
    }
}

#[test]
fn g0_zero_bending_keeps_native_motion_and_torsion_is_not_counted_as_bending() {
    let gate = CancelGate::new();
    for eta in [0.0, 800.0] {
        let source = card(eta, (0.0, f64::MAX));
        let model = pinned_model(&source, 300.0);
        let bending =
            PreparedRodBendingLoss::new(&model, vec![viscous(&source, 300.0)], &gate).unwrap();
        let mut pose = model.initial_pose();
        pose.frames[1] = about(if eta == 0.0 { 2 } else { 0 }, 0.05);
        let mut angular = [[0.0; 3]; 2];
        angular[1][if eta == 0.0 { 2 } else { 0 }] = 0.3;
        let state = RodMotion::new(model, pose, &[[0.0; 3]; 2], &angular, 0.0, &gate).unwrap();
        let loads = [RodNodalWrench::default(); 2];
        let (actual, receipt) = state
            .step_damped_candidate(
                0.001,
                RodAngularLosses {
                    bending: Some(&bending),
                    torsional: None,
                },
                &loads,
                controls(0.001),
                &gate,
            )
            .unwrap();
        let (native, _) = state
            .step_candidate(0.001, &loads, controls(0.001), &gate)
            .unwrap();
        assert_eq!(actual.pose(), native.pose());
        assert_eq!(actual.body_spins_kg_m2_s(), native.body_spins_kg_m2_s());
        assert_eq!(receipt.bending_heat.unwrap().total_heat_j, 0.0);
        assert!(receipt.torsional_heat.is_none());
    }
}

#[test]
fn g0_zero_bending_preserves_the_existing_multisegment_torsional_sweep() {
    let gate = CancelGate::new();
    let source = card_with_shear(0.0, 800.0, (0.0, f64::MAX));
    let model = pinned_segments(&source, 300.0, 2);
    let bending =
        PreparedRodBendingLoss::new(&model, vec![viscous(&source, 300.0); 2], &gate).unwrap();
    let torsion = PreparedRodTorsionalLoss::new(
        &model,
        vec![
            sourced(&source, 300.0, true)
                .with_torsional_viscosity("shear_eta")
                .unwrap();
            2
        ],
        &gate,
    )
    .unwrap();
    let mut pose = model.initial_pose();
    pose.frames[1] = about(0, 0.03);
    pose.frames[2] = about(0, -0.05);
    let state = RodMotion::new(
        model,
        pose,
        &[[0.0; 3]; 3],
        &[[0.0; 3], [0.3, 0.0, 0.0], [-0.1, 0.0, 0.0]],
        0.0,
        &gate,
    )
    .unwrap();
    let loads = [RodNodalWrench::default(); 3];
    let (original, expected) = state
        .step_torsional_candidate(0.001, &torsion, &loads, controls(0.001), &gate)
        .unwrap();
    let (next, receipt) = state
        .step_damped_candidate(
            0.001,
            RodAngularLosses {
                bending: Some(&bending),
                torsional: Some(&torsion),
            },
            &loads,
            controls(0.001),
            &gate,
        )
        .unwrap();
    assert_eq!(next.pose(), original.pose());
    assert_eq!(next.body_spins_kg_m2_s(), original.body_spins_kg_m2_s());
    assert_eq!(receipt.torsional_heat, expected.torsional_heat);
    assert_eq!(receipt.bending_heat.unwrap().total_heat_j, 0.0);
}

#[test]
fn g0_g3_simultaneous_bending_torsion_and_eccentric_inertia_are_objective() {
    let gate = CancelGate::new();
    let source = card(800.0, (0.0, f64::MAX));
    let mechanical = sourced(&source, 300.0, false);
    let offset = [0.02, 0.04, -0.03];
    let link = MaterialTensionLink::try_new(
        mechanical.elastic().clone(),
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
    let model = Arc::new(
        SpatialMaterialRod::try_straight(0.8, vec![mechanical])
            .unwrap()
            .with_tension_links(vec![link], &gate)
            .unwrap()
            .into_spatial_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    assert_ne!(model.inertia()[0].inertia_com_body_kg_m2()[0][1], 0.0);
    let bending =
        PreparedRodBendingLoss::new(&model, vec![viscous(&source, 300.0)], &gate).unwrap();
    let torsion = PreparedRodTorsionalLoss::new(
        &model,
        vec![
            sourced(&source, 300.0, true)
                .with_torsional_viscosity("shear_eta")
                .unwrap(),
        ],
        &gate,
    )
    .unwrap();
    let mut pose = model.initial_pose();
    pose.frames[1] = So3::exp(So3Tangent::new(Vec3::new(0.1, 0.2, -0.15))).unwrap();
    let velocity = [[0.1, -0.2, 0.03], [-0.03, 0.02, 0.04]];
    let angular = [[0.3, -0.2, 0.1], [-0.2, 0.1, -0.3]];
    let base = SpatialRodMotion::new_spatial(
        model.clone(),
        pose.clone(),
        &velocity,
        &angular,
        0.0,
        1e-12,
        &gate,
    )
    .unwrap();
    let global = So3::exp(So3Tangent::new(Vec3::new(0.4, -0.2, 0.3))).unwrap();
    for (point, frame) in pose.positions_m.iter_mut().zip(&mut pose.frames) {
        *point = rotate(global, *point);
        *frame = global.compose(*frame).unwrap();
    }
    let rotated = SpatialRodMotion::new_spatial(
        model,
        pose,
        &velocity.map(|v| rotate(global, v)),
        &angular,
        0.0,
        1e-12,
        &gate,
    )
    .unwrap();
    let mut accepted = Vec::new();
    for state in [base, rotated] {
        let (next, receipt) = state
            .step_damped_candidate(
                1e-4,
                RodAngularLosses {
                    bending: Some(&bending),
                    torsional: Some(&torsion),
                },
                &[RodNodalWrench::default(); 2],
                controls(1e-4),
                &gate,
            )
            .unwrap();
        let bending_heat = receipt.bending_heat.as_ref().unwrap().total_heat_j;
        let torsional_heat = receipt.torsional_heat.as_ref().unwrap().total_heat_j;
        assert!(bending_heat > 0.0 && torsional_heat > 0.0);
        near(
            receipt.integration_defect_j,
            receipt.after.kinetic_energy_j - receipt.before.kinetic_energy_j
                + receipt.after.potential_energy_j
                - receipt.before.potential_energy_j
                + bending_heat
                + torsional_heat,
            2e-14,
        );
        vector_near(receipt.linear_balance_residual_n_s, [0.0; 3], 2e-12);
        vector_near(receipt.angular_balance_residual_nm_s, [0.0; 3], 2e-12);
        accepted.push((next, receipt));
    }
    near(
        accepted[0].1.bending_heat.as_ref().unwrap().total_heat_j,
        accepted[1].1.bending_heat.as_ref().unwrap().total_heat_j,
        2e-13,
    );
    near(
        accepted[0].1.torsional_heat.as_ref().unwrap().total_heat_j,
        accepted[1].1.torsional_heat.as_ref().unwrap().total_heat_j,
        2e-13,
    );
    for node in 0..2 {
        vector_near(
            accepted[1].0.pose().positions_m[node],
            rotate(global, accepted[0].0.pose().positions_m[node]),
            2e-12,
        );
        vector_near(
            accepted[1].0.body_spins_kg_m2_s()[node],
            accepted[0].0.body_spins_kg_m2_s()[node],
            2e-12,
        );
    }
}

#[test]
fn g4_band_source_and_stale_state_refuse_before_motion_or_heat_is_published() {
    let gate = CancelGate::new();
    let source = card(800.0, (0.0, f64::MAX));
    let model = pinned_model(&source, 300.0);
    let loss = PreparedRodBendingLoss::new(&model, vec![viscous(&source, 300.0)], &gate).unwrap();
    let foreign = card(1600.0, (0.0, f64::MAX));
    assert!(PreparedRodBendingLoss::new(&model, vec![viscous(&foreign, 300.0)], &gate).is_err());
    let narrow = card(800.0, (1.0, 100.0));
    assert!(
        PreparedRodBendingLoss::new(
            &pinned_model(&narrow, 300.0),
            vec![viscous(&narrow, 300.0)],
            &gate
        )
        .is_err()
    );
    let hot_model = pinned_model(&source, 320.0);
    let mut pose = hot_model.initial_pose();
    pose.frames[1] = about(2, 0.05);
    let state = RodMotion::new(
        hot_model.clone(),
        pose.clone(),
        &[[0.0; 3]; 2],
        &[[0.0; 3], [0.0, 0.0, 0.3]],
        0.0,
        &gate,
    )
    .unwrap();
    let loads = [RodNodalWrench::default(); 2];
    assert!(
        state
            .step_damped_candidate(
                0.001,
                RodAngularLosses {
                    bending: Some(&loss),
                    torsional: None
                },
                &loads,
                controls(0.001),
                &gate
            )
            .is_err()
    );
    let fresh =
        PreparedRodBendingLoss::new(&hot_model, vec![viscous(&source, 320.0)], &gate).unwrap();
    let selected = RodAngularLosses {
        bending: Some(&fresh),
        torsional: None,
    };
    let mut strict = controls(0.001);
    strict.energy_defect_limit_j = 1e-30;
    assert!(
        state
            .step_damped_candidate(0.001, selected, &loads, strict, &gate)
            .is_err()
    );
    let (_, retry) = state
        .step_damped_candidate(0.001, selected, &loads, controls(0.001), &gate)
        .unwrap();
    assert!(retry.bending_heat.unwrap().total_heat_j > 0.0);
    gate.request();
    assert!(
        state
            .step_damped_candidate(0.001, selected, &loads, controls(0.001), &gate)
            .is_err()
    );
    assert_eq!(state.accepted_steps(), 0);
    assert_eq!(state.pose(), &pose);
}
