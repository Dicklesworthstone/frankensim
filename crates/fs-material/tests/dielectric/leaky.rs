//! Maxwell relaxation and current/heat identities for a source-backed leaky dielectric.
use super::*;
use fs_material::{
    conductor::ELECTRICAL_CONDUCTIVITY_DIMS,
    field_storage::dielectric::leaky::ResolvedLeakyDielectric,
    state_point::{
        ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};

fn source(epsilon: f64, sigma: f64) -> MaterialCard {
    card(
        &[
            ("epsilon_inf", ELECTRIC_PERMITTIVITY_DIMS, epsilon),
            ("delta1", ELECTRIC_PERMITTIVITY_DIMS, 3e-6),
            ("tau1", Time::DIMS, 0.1),
            ("dc_leakage", ELECTRICAL_CONDUCTIVITY_DIMS, sigma),
            ("density", Dims([-3, 1, 0, 0, 0, 0]), 1234.0),
        ],
        ValidityDomain::unconstrained(),
    )
}

fn law(epsilon: f64, sigma: f64, pole_count: usize) -> ResolvedLeakyDielectric {
    ResolvedLeakyDielectric::resolve(
        &source(epsilon, sigma),
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &properties()[..pole_count],
        "dc_leakage",
    )
    .unwrap()
}

fn complete_bundle(source: &MaterialCard) -> ResolvedMaterialStatePoint {
    let requirements = [
        ("epsilon_inf", ELECTRIC_PERMITTIVITY_DIMS),
        ("delta1", ELECTRIC_PERMITTIVITY_DIMS),
        ("tau1", Time::DIMS),
        ("dc_leakage", ELECTRICAL_CONDUCTIVITY_DIMS),
        ("density", Dims([-3, 1, 0, 0, 0, 0])),
    ]
    .into_iter()
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap()
    })
    .collect::<Vec<_>>();
    resolve_material_state_point(
        source,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn close(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 5e-13 * scale, "{a} != {b}, scale={scale}");
}

#[test]
fn g1_one_original_bundle_retains_leakage_name_extra_properties_and_zero_limit() {
    let bundle = complete_bundle(&source(2e-6, 4e-6));
    let combined = ResolvedLeakyDielectric::try_from_material(
        bundle.clone(),
        "epsilon_inf",
        &properties()[..1],
        "dc_leakage",
    )
    .unwrap();
    assert_eq!(combined.material(), &bundle);
    assert_eq!(combined.material().properties().len(), 5);
    assert_eq!(combined.conductivity_property(), "dc_leakage");
    assert_eq!(combined.conductivity_s_m(), 4e-6);
    let (dielectric, name, sigma) = combined.into_parts();
    assert_eq!(dielectric.material(), &bundle);
    assert_eq!(name, "dc_leakage");
    assert_eq!(sigma, 4e-6);
    assert!(
        ResolvedLeakyDielectric::try_from_dielectric(dielectric, "electrical_conductivity",)
            .is_err()
    );

    let insulator = law(2e-6, 0.0, 1);
    let dt = 0.03;
    let e0 = [200.0, -80.0];
    let e1 = [-20.0, 130.0];
    let p0 = [[1e-4, -2e-4]];
    let pure = insulator.dielectric();
    assert_eq!(
        insulator.midpoint_permittivity_f_m(dt).unwrap(),
        pure.midpoint_permittivity_f_m(dt).unwrap(),
    );
    assert_eq!(
        insulator.midpoint_history_c_m2(dt, e0, &p0).unwrap(),
        pure.midpoint_history_c_m2(dt, e0, &p0).unwrap(),
    );
    let step = insulator.step_midpoint(dt, e0, e1, &p0).unwrap();
    assert_eq!(
        step.dielectric,
        pure.step_midpoint(dt, e0, e1, &p0).unwrap()
    );
    assert_eq!(step.joule_heat_j_m3, 0.0);
    assert_eq!(step.ohmic_current_density_a_m2, [0.0; 2]);
    assert_eq!(step.total_heat_j_m3, step.dielectric.heat_j_m3);
    assert_eq!(step.supplied_work_j_m3, step.dielectric.supplied_work_j_m3);
}

#[test]
fn g1_condensed_current_matches_recovery_and_separate_heat_channels() {
    let law = law(2e-6, 4e-6, 1);
    let dt = 0.037;
    let e0 = [230.0, -100.0, 40.0];
    let e1 = [-80.0, 310.0, 90.0];
    let p0 = [[1e-4, -2e-4, 3e-4]];
    let effective = law.midpoint_permittivity_f_m(dt).unwrap();
    let history = law.midpoint_history_c_m2(dt, e0, &p0).unwrap();
    let step = law.step_midpoint(dt, e0, e1, &p0).unwrap();
    near(
        effective,
        law.dielectric().midpoint_permittivity_f_m(dt).unwrap() + dt * law.conductivity_s_m() / 2.0,
    );
    let mut expected_joule = 0.0;
    let mut current_work = 0.0;
    for i in 0..3 {
        let midpoint = 0.5 * (e0[i] + e1[i]);
        let conduction = dt * 4e-6 * midpoint;
        near(step.integrated_ohmic_current_c_m2[i], conduction);
        close(
            step.dielectric.displacement_after_c_m2[i] + conduction,
            effective * e1[i] + history[i],
            1e-3,
        );
        close(
            step.total_current_increment_c_m2[i],
            step.dielectric.displacement_after_c_m2[i]
                - step.dielectric.displacement_before_c_m2[i]
                + conduction,
            1e-3,
        );
        expected_joule += dt * 4e-6 * midpoint * midpoint;
        current_work += midpoint * step.total_current_increment_c_m2[i];
    }
    near(step.joule_heat_j_m3, expected_joule);
    assert!(step.dielectric.heat_j_m3 > 0.0);
    near(
        step.total_heat_j_m3,
        step.dielectric.heat_j_m3 + expected_joule,
    );
    near(step.supplied_work_j_m3, current_work);
    close(
        step.dielectric.stored_change_j_m3 + step.total_heat_j_m3,
        step.supplied_work_j_m3,
        1.0,
    );
    close(step.energy_residual_j_m3, 0.0, 1.0);
}

#[test]
fn g1_open_circuit_maxwell_relaxation_is_second_order_and_passive() {
    let epsilon = 2e-6;
    let sigma = 4e-6;
    let law = law(epsilon, sigma, 0);
    let initial = 50.0;
    let duration = 0.3;
    let exact = initial * f64::exp(-duration * sigma / epsilon);
    let integrate = |steps: usize| {
        let dt = duration / steps as f64;
        let mut electric = initial;
        let mut heat = 0.0;
        for _ in 0..steps {
            let history = law.midpoint_history_c_m2(dt, [electric], &[]).unwrap()[0];
            let next = (epsilon * electric - history) / law.midpoint_permittivity_f_m(dt).unwrap();
            let step = law.step_midpoint(dt, [electric], [next], &[]).unwrap();
            close(step.total_current_increment_c_m2[0], 0.0, 1e-3);
            close(
                step.dielectric.stored_change_j_m3 + step.joule_heat_j_m3,
                0.0,
                1e-2,
            );
            assert!(step.joule_heat_j_m3 >= 0.0);
            heat += step.joule_heat_j_m3;
            electric = next;
        }
        close(
            heat,
            0.5 * epsilon * (initial * initial - electric * electric),
            1e-2,
        );
        electric
    };
    let error16 = (integrate(16) - exact).abs();
    let error32 = (integrate(32) - exact).abs();
    let error64 = (integrate(64) - exact).abs();
    assert!((3.99..4.01).contains(&(error16 / error32)));
    assert!((3.99..4.01).contains(&(error32 / error64)));

    // Midpoint may reverse sign for dt > 2 epsilon/sigma but still dissipates.
    let dt = 2.0;
    let h = law.midpoint_history_c_m2(dt, [initial], &[]).unwrap()[0];
    let end = (epsilon * initial - h) / law.midpoint_permittivity_f_m(dt).unwrap();
    let stiff = law.step_midpoint(dt, [initial], [end], &[]).unwrap();
    assert!(end < 0.0 && end.abs() < initial);
    assert!(stiff.joule_heat_j_m3 > 0.0);
    close(
        stiff.dielectric.stored_change_j_m3 + stiff.joule_heat_j_m3,
        0.0,
        1e-2,
    );
}

#[test]
fn g1_equilibrated_poles_do_not_double_count_dc_joule_heat() {
    let law = law(2e-6, 8e-6, 1);
    let e = [128.0, -64.0];
    let p = [[3e-6 * e[0], 3e-6 * e[1]]];
    let step = law.step_midpoint(0.25, e, e, &p).unwrap();
    assert_eq!(step.dielectric.heat_j_m3, 0.0);
    assert_eq!(step.dielectric.stored_change_j_m3, 0.0);
    assert_eq!(step.dielectric.supplied_work_j_m3, 0.0);
    near(
        step.joule_heat_j_m3,
        0.25 * 8e-6 * (128.0 * 128.0 + 64.0 * 64.0),
    );
    assert_eq!(step.total_heat_j_m3, step.joule_heat_j_m3);
    assert_eq!(step.supplied_work_j_m3, step.joule_heat_j_m3);
}

#[test]
fn g3_current_increments_survive_dc_bias_and_coefficient_overflow_reordering() {
    let biased = law(1e-12, 1e-30, 0);
    let e0 = 1e12_f64;
    let e1 = e0 + 1e-4;
    let step = biased.step_midpoint(1.0, [e0], [e1], &[]).unwrap();
    let expected = 1e-12 * (e1 - e0) + 1e-30 * f64::midpoint(e0, e1);
    near(step.total_current_increment_c_m2[0], expected);

    // dt*sigma overflows while dt*sigma/2 and the physical current moment fit.
    let scaled = law(1e-100, 2.0, 0);
    let dt = 1e308;
    near(scaled.midpoint_permittivity_f_m(dt).unwrap(), 1e308);
    let step = scaled.step_midpoint(dt, [1e-100], [1e-100], &[]).unwrap();
    near(step.integrated_ohmic_current_c_m2[0], 2e208);
    near(step.joule_heat_j_m3, 2e108);
    assert_eq!(step.energy_residual_j_m3, 0.0);

    // Round once at the half-product: (dt*sigma)/2 would first round the
    // full product to the minimum subnormal, then incorrectly return zero.
    let minimum = f64::from_bits(1);
    let subnormal = law(1e-320, 0.4, 0);
    assert_eq!(
        subnormal
            .midpoint_history_c_m2(3.0 * minimum, [1.0], &[])
            .unwrap(),
        [minimum],
    );
    assert_eq!(
        subnormal.midpoint_permittivity_f_m(3.0 * minimum).unwrap(),
        1e-320 + minimum,
    );
    // Conversely, halving dt first would erase a finite normal coefficient.
    let tiny_interval = law(1e-12, 2.0_f64.powi(1023), 0);
    assert_eq!(
        tiny_interval
            .midpoint_history_c_m2(minimum, [1.0], &[])
            .unwrap(),
        [2.0_f64.powi(-52)],
    );
}

#[test]
fn g3_negative_leakage_and_unowned_frequency_dependence_refuse() {
    let negative = complete_bundle(&source(2e-6, -1.0));
    assert!(
        ResolvedLeakyDielectric::try_from_material(
            negative,
            "epsilon_inf",
            &properties()[..1],
            "dc_leakage",
        )
        .is_err()
    );
    let frequency = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Frequency(FrequencyConvention::Cyclic),
        ValueForm::Static,
    ));
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("epsilon_inf", ELECTRIC_PERMITTIVITY_DIMS, 2e-6),
        ("dc_leakage", ELECTRICAL_CONDUCTIVITY_DIMS, 4e-6),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: if name == "dc_leakage" {
                    ValidityDomain::unconstrained().with_quantity(
                        "frequency",
                        frequency,
                        1.0,
                        100.0,
                    )
                } else {
                    ValidityDomain::unconstrained()
                },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic conductivity-domain adversary".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let source = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "dielectric with dispersive conductivity".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new()
        .with_quantity("frequency", frequency, 10.0)
        .unwrap();
    assert!(
        ResolvedLeakyDielectric::resolve(
            &source,
            &query,
            MaterialPropertySelection::SingleClaimOnly,
            "epsilon_inf",
            &[],
            "dc_leakage",
        )
        .is_err()
    );
    let law = law(2e-6, 4e-6, 0);
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(law.midpoint_permittivity_f_m(dt).is_err());
    }
    assert!(law.step_midpoint(0.1, [1.0], [2.0], &[[0.0]]).is_err());
}
