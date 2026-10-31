//! Source-derived torsional damping on the actual native rod trajectory.
use super::*;
use fs_couple::rod_motion::{PreparedRodLoads, PreparedRodTorsionalLoss, SpatialRodMotion};
use fs_material::section::torsion::ConcentricTorsionalSection;
use fs_material::state_point::INVERSE_TEMPERATURE_DIMS;
use fs_qty::DynViscosity;
use fs_solid::rod::material_field::thermal::source::{
    RodLayerSource, RodMaterialSource, RodSegmentSource,
};

fn card(eta: f64, band: (f64, f64)) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density", 1000.0, Density::DIMS),
        ("young_modulus", 1e4, Pressure::DIMS),
        ("shear_modulus", 4e3, Pressure::DIMS),
        ("shear_viscosity", eta, DynViscosity::DIMS),
        (
            "linear_thermal_expansion_coefficient",
            1e-5,
            INVERSE_TEMPERATURE_DIMS,
        ),
    ] {
        let mut validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
        if name == "shear_viscosity" {
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
                    source: "synthetic torsional viscosity".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic shaft".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", temperature)
        .unwrap()
        .with("omega", 10.0)
        .unwrap()
}
fn sourced(card: &MaterialCard, temperature: f64, loss: bool) -> ConcentricShearSection {
    let mut requirements = vec![
        ScalarPropertyRequirement::try_new(
            "density",
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
        ScalarPropertyRequirement::try_new(
            "young_modulus",
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
        ScalarPropertyRequirement::try_new(
            "shear_modulus",
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    ];
    if loss {
        requirements.push(
            ScalarPropertyRequirement::try_new(
                "shear_viscosity",
                DynViscosity::DIMS,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
        );
    }
    let material = resolve_material_state_point(
        card,
        &point(temperature),
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
fn viscous(card: &MaterialCard) -> ConcentricTorsionalSection {
    sourced(card, 300.0, true)
        .with_torsional_viscosity("shear_viscosity")
        .unwrap()
}
fn pinned_model(card: &MaterialCard, segments: usize) -> Arc<PreparedRodMotion> {
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
        SpatialMaterialRod::try_straight(0.8, vec![sourced(card, 300.0, false); segments])
            .unwrap()
            .with_supports(&supports, &gate)
            .unwrap()
            .into_motion_model(RodRootBoundary::Clamped, &gate)
            .unwrap(),
    )
}
fn rx(angle: f64) -> So3 {
    So3::exp(So3Tangent::new(Vec3::new(angle, 0.0, 0.0))).unwrap()
}
fn vector_near(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    for c in 0..3 {
        near(a[c], b[c], tolerance);
    }
}
fn rotate(frame: So3, value: [f64; 3]) -> [f64; 3] {
    let out = frame
        .rotate(Vec3::new(value[0], value[1], value[2]))
        .unwrap();
    [out.x, out.y, out.z]
}

#[test]
fn g1_exact_pair_relaxation_heat_and_native_elastic_split_agree() {
    let gate = CancelGate::new();
    for eta in [8000.0, 8e8] {
        let source = card(eta, (0.0, f64::MAX));
        let model = pinned_model(&source, 1);
        let loss = PreparedRodTorsionalLoss::new(&model, vec![viscous(&source)], &gate).unwrap();
        let k = model.source().sections()[0].torsional_stiffness_n_m2() / 0.8;
        let inertia = model.inertia()[1].principal_kg_m2()[0];
        let c = eta * (core::f64::consts::PI * 0.05_f64.powi(4) / 2.0) / 0.8;
        near(loss.resistance_nm_s()[0], c, 2e-14 * c);
        let theta0 = 0.08;
        let w0 = 0.3;
        let dt = 0.002;
        let torque = 0.002;
        let mut pose = model.initial_pose();
        pose.frames[1] = rx(theta0);
        let state = RodMotion::new(
            model,
            pose,
            &[[0.0; 3]; 2],
            &[[0.0; 3], [w0, 0.0, 0.0]],
            0.0,
            &gate,
        )
        .unwrap();
        let decay = (-0.5 * dt * c / inertia).exp();
        let wa = decay * w0;
        let middle = wa + 0.5 * dt * (torque - k * theta0) / inertia;
        let theta1 = theta0 + dt * middle;
        let wb = middle + 0.5 * dt * (torque - k * theta1) / inertia;
        let expected_heat = 0.5 * inertia * (1.0 - decay * decay) * (w0 * w0 + wb * wb);
        let loads = [
            RodNodalWrench::default(),
            RodNodalWrench {
                force_n: [0.0; 3],
                torque_nm: [torque, 0.0, 0.0],
            },
        ];
        let (next, report) = state
            .step_torsional_candidate(dt, &loss, &loads, controls(dt), &gate)
            .unwrap();
        near(
            next.pose().frames[1]
                .body_minus(So3::identity())
                .unwrap()
                .angular
                .x,
            theta1,
            2e-13,
        );
        near(
            next.angular_velocities_body_rad_s().unwrap()[1][0],
            decay * wb,
            2e-12,
        );
        let heat = report.torsional_heat.as_ref().unwrap();
        near(heat.total_heat_j, expected_heat, 2e-14);
        assert_eq!(heat.segment_heat_j, vec![heat.total_heat_j]);
        near(
            report.support_impulses[0].angular_nm_s[0],
            inertia * (decay * wb - w0) - dt * torque,
            2e-13,
        );
        vector_near(report.angular_balance_residual_nm_s, [0.0; 3], 2e-13);
        near(
            report.after.kinetic_energy_j + report.after.potential_energy_j
                - report.before.kinetic_energy_j
                - report.before.potential_energy_j
                + heat.total_heat_j
                - report.source_work_quadrature_j,
            report.integration_defect_j,
            2e-14,
        );
        let spatial = PreparedRodLoads::new(&loads, &[[0.0; 3]; 2], &[], &gate).unwrap();
        let (other, other_report) = state
            .step_spatial_torsional_candidate(dt, &loss, &spatial, controls(dt), &gate)
            .unwrap();
        assert_eq!(other.pose(), next.pose());
        assert_eq!(other.body_spins_kg_m2_s(), next.body_spins_kg_m2_s());
        assert_eq!(other_report.torsional_heat, report.torsional_heat);
        if eta == 8000.0 {
            let mut strict = controls(dt);
            strict.energy_defect_limit_j = 0.0;
            assert!(
                state
                    .step_torsional_candidate(dt, &loss, &loads, strict, &gate)
                    .is_err()
            );
            assert_eq!(state.accepted_steps(), 0);
        }
    }
}

#[test]
fn g1_two_segment_viscous_shaft_follows_analytic_mode_with_second_order_heat() {
    let gate = CancelGate::new();
    for eta in [800.0, 1600.0] {
        let source = card(eta, (0.0, f64::MAX));
        let model = pinned_model(&source, 2);
        let loss = PreparedRodTorsionalLoss::new(&model, vec![viscous(&source); 2], &gate).unwrap();
        let stiffness = model.source().sections()[0].torsional_stiffness_n_m2() / 0.4;
        let tip_inertia = model.inertia()[2].principal_kg_m2()[0];
        let shape = [0.0, 0.5_f64.sqrt(), 1.0];
        let omega2 = (stiffness / tip_inertia) * (1.0 - shape[1]);
        let alpha = 0.5 * (eta / 4000.0) * omega2;
        let frequency = (omega2 - alpha * alpha).sqrt();
        let horizon = 0.2;
        let a0 = 0.08;
        let v0 = 0.2;
        let sine = (frequency * horizon).sin();
        let cosine = (frequency * horizon).cos();
        let factor = (-alpha * horizon).exp();
        let exact_a = factor * (a0 * cosine + (v0 + alpha * a0) * sine / frequency);
        let exact_v =
            -alpha * exact_a + factor * (-a0 * frequency * sine + (v0 + alpha * a0) * cosine);
        let energy = |a: f64, w: f64| {
            0.5 * tip_inertia * (2.0 * shape[1] * shape[1] + 1.0) * w * w
                + 0.5 * stiffness * (shape[1].powi(2) + (1.0 - shape[1]).powi(2)) * a * a
        };
        let exact_heat = energy(a0, v0) - energy(exact_a, exact_v);
        let mut errors = Vec::new();
        let mut heat_errors = Vec::new();
        for count in [100, 200, 400] {
            let dt = horizon / count as f64;
            let mut pose = model.initial_pose();
            for node in 1..3 {
                pose.frames[node] = rx(a0 * shape[node]);
            }
            let mut state = RodMotion::new(
                model.clone(),
                pose,
                &[[0.0; 3]; 3],
                &shape.map(|s| [v0 * s, 0.0, 0.0]),
                0.0,
                &gate,
            )
            .unwrap();
            let mut total_heat = 0.0;
            for _ in 0..count {
                let (next, report) = state
                    .step_torsional_candidate(
                        dt,
                        &loss,
                        &[RodNodalWrench::default(); 3],
                        controls(dt),
                        &gate,
                    )
                    .unwrap();
                let heat = report.torsional_heat.unwrap();
                assert!(heat.segment_heat_j.iter().all(|q| *q >= 0.0));
                near(heat.segment_heat_j.iter().sum(), heat.total_heat_j, 1e-16);
                total_heat += heat.total_heat_j;
                state = next;
            }
            let a = state.pose().frames[2]
                .body_minus(So3::identity())
                .unwrap()
                .angular
                .x;
            let w = state.angular_velocities_body_rad_s().unwrap()[2][0];
            errors.push((a - exact_a).abs() + (w - exact_v).abs());
            heat_errors.push((total_heat - exact_heat).abs());
        }
        for j in 0..2 {
            assert!(
                (3.6..4.4).contains(&(errors[j] / errors[j + 1])),
                "{eta}: {errors:?}"
            );
            assert!(
                (3.6..4.4).contains(&(heat_errors[j] / heat_errors[j + 1])),
                "{eta}: {heat_errors:?}"
            );
        }
    }
}

#[test]
fn g0_collective_spin_has_no_material_heat_and_zero_viscosity_preserves_native_path() {
    let gate = CancelGate::new();
    for eta in [0.0, 8000.0] {
        let source = card(eta, (0.0, f64::MAX));
        let model = Arc::new(
            SpatialMaterialRod::try_straight(0.8, vec![sourced(&source, 300.0, false)])
                .unwrap()
                .into_motion_model(RodRootBoundary::Free, &gate)
                .unwrap(),
        );
        let loss = PreparedRodTorsionalLoss::new(&model, vec![viscous(&source)], &gate).unwrap();
        let state = RodMotion::new(
            model.clone(),
            model.initial_pose(),
            &[[0.0; 3]; 2],
            &[[2.0, 0.0, 0.0]; 2],
            0.0,
            &gate,
        )
        .unwrap();
        let loads =
            PreparedRodLoads::new(&[RodNodalWrench::default(); 2], &[[0.0; 3]; 2], &[], &gate)
                .unwrap();
        let (actual, report) = state
            .step_spatial_torsional_candidate(0.001, &loss, &loads, controls(0.001), &gate)
            .unwrap();
        let (native, plain) = state
            .step_spatial_candidate(0.001, &loads, controls(0.001), &gate)
            .unwrap();
        assert_eq!(actual.pose(), native.pose());
        assert_eq!(actual.body_spins_kg_m2_s(), native.body_spins_kg_m2_s());
        assert_eq!(report.torsional_heat.unwrap().total_heat_j, 0.0);
        assert!(plain.torsional_heat.is_none());
    }
}

#[test]
fn g0_g3_finite_bend_eccentric_inertia_and_global_frame_preserve_momentum_and_heat() {
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
    let loss = PreparedRodTorsionalLoss::new(&model, vec![viscous(&source)], &gate).unwrap();
    assert_ne!(model.inertia()[0].inertia_com_body_kg_m2()[0][1], 0.0);
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
    for (p, r) in pose.positions_m.iter_mut().zip(&mut pose.frames) {
        *p = rotate(global, *p);
        *r = global.compose(*r).unwrap();
    }
    let transformed = SpatialRodMotion::new_spatial(
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
    for state in [base, transformed] {
        let (next, report) = state
            .step_torsional_candidate(
                1e-4,
                &loss,
                &[RodNodalWrench::default(); 2],
                controls(1e-4),
                &gate,
            )
            .unwrap();
        vector_near(report.linear_balance_residual_n_s, [0.0; 3], 2e-12);
        vector_near(report.angular_balance_residual_nm_s, [0.0; 3], 2e-12);
        assert!(report.torsional_heat.as_ref().unwrap().total_heat_j > 0.0);
        accepted.push((next, report));
    }
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
fn g4_narrow_band_wrong_source_and_actual_temperature_refresh_refuse_atomically() {
    let gate = CancelGate::new();
    let source = card(800.0, (0.0, f64::MAX));
    let model = pinned_model(&source, 1);
    let foreign = card(1600.0, (0.0, f64::MAX));
    assert!(PreparedRodTorsionalLoss::new(&model, vec![viscous(&foreign)], &gate).is_err());
    let narrow = card(800.0, (1.0, 100.0));
    let narrow_model = pinned_model(&narrow, 1);
    assert!(PreparedRodTorsionalLoss::new(&narrow_model, vec![viscous(&narrow)], &gate).is_err());

    let thermal = RodMaterialSource::try_new(
        0.8,
        vec![RodSegmentSource {
            inner_radius_m: 0.0,
            shear_correction: 0.8,
            layers: vec![RodLayerSource {
                outer_radius_m: 0.05,
                card: &source,
                reference_point: point(300.0),
                mechanical_selection: MaterialPropertySelection::SingleClaimOnly,
                expansion_selection: MaterialPropertySelection::SingleClaimOnly,
            }],
        }],
        "shear_modulus",
        0.2,
        &gate,
    )
    .unwrap();
    let model = Arc::new(
        thermal
            .reference_rod()
            .clone()
            .into_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let loss = PreparedRodTorsionalLoss::new(&model, vec![viscous(&source)], &gate).unwrap();
    let state = RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[[0.0; 3]; 2],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap();
    let (hot, _) = state
        .rebind_material_temperatures_candidate(&thermal, &[350.0], &gate)
        .unwrap();
    assert!(
        hot.step_torsional_candidate(
            0.001,
            &loss,
            &[RodNodalWrench::default(); 2],
            controls(0.001),
            &gate
        )
        .is_err()
    );
    assert_eq!(hot.accepted_steps(), 0);
    assert_eq!(hot.body_spins_kg_m2_s(), state.body_spins_kg_m2_s());
    gate.request();
    assert!(
        state
            .step_torsional_candidate(
                0.001,
                &loss,
                &[RodNodalWrench::default(); 2],
                controls(0.001),
                &gate
            )
            .is_err()
    );
    assert_eq!(state.accepted_steps(), 0);
    assert_eq!(state.pose(), &model.initial_pose());
}
