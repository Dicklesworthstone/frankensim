//! Source-pinned Hill saturation, reciprocal thermal work and native sections.
//! All coefficients are authored numerical fixtures, not alloy calibration.

use super::*;
use fs_material::plastic::hill::oriented::OrientedHill48;
use fs_material::plastic::hill::source::{
    HILL_KINEMATIC_MODULUS_PROPERTY as C_PROPERTY,
    HILL_VOCE_SATURATION_INCREMENT_PROPERTY as Q_PROPERTY,
    HILL_VOCE_SATURATION_RATE_PROPERTY as B_PROPERTY, ResolvedHill48, hill48_voce_requirements,
    resolve_hill48_voce,
};
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};

const Q: f64 = 2e6;
const B: f64 = 180.0;
const Q_SLOPE: f64 = -1e5;
const C: f64 = 4e8;
const C_SLOPE: f64 = -2e7;

fn voce_properties(
    expansion: f64,
    varying: bool,
    q: f64,
    slope: f64,
    b: f64,
) -> Vec<PropertyClaim> {
    let mut values = properties(expansion, varying);
    values.push(property(Q_PROPERTY, Pressure::DIMS, q, slope));
    values.push(property(B_PROPERTY, Dims([0; 6]), b, 0.0));
    values
}

fn voce_source(fixture: &Fixture, rotation: [[f64; 3]; 3]) -> ThermalHillSource<'_> {
    fixture
        .source(rotation)
        .with_voce_isotropic_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
}

// Independent primitive at ordinary b*alpha; the tiny-exponent case below
// uses its independently known leading limit instead of subtracting this form.
fn g(alpha: f64) -> f64 {
    alpha + (-B * alpha).exp_m1() / B
}

#[test]
fn fixed_source_retains_hill_specific_coefficients_axes_receipts_and_native_cycles() {
    let mut values = voce_properties(0.0, false, Q, 0.0, B);
    values.push(property(C_PROPERTY, Pressure::DIMS, C, 0.0));
    let fixture = Fixture::new(values);
    let mut requirements = hill48_voce_requirements().unwrap();
    for (name, dims) in [
        (C_PROPERTY, Pressure::DIMS),
        (DENSITY_PROPERTY, Density::DIMS),
    ] {
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::NonNegative)
                .unwrap(),
        );
    }
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    let bundle = resolve_material_state_point(
        &fixture.card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let rotation = fs_material::tensor::rotation([0.0, 0.0, 1.0], 0.37);
    let plain = ResolvedHill48::try_from_material(bundle.clone(), 0.04, rotation).unwrap();
    assert_eq!(plain.local_law().voce_saturation_pa(), 0.0);
    assert_eq!(plain.local_law().kinematic_modulus_pa(), 0.0);
    let selected = plain
        .clone()
        .with_voce_isotropic_hardening()
        .unwrap()
        .with_linear_kinematic_hardening()
        .unwrap();
    let reversed_order = plain
        .clone()
        .with_linear_kinematic_hardening()
        .unwrap()
        .with_voce_isotropic_hardening()
        .unwrap();
    assert_eq!(selected, reversed_order);
    assert_eq!(selected.material(), &bundle);
    for receipt in selected.material().properties() {
        fixture
            .card
            .claims()
            .verify_receipt(&receipt.answer().receipt)
            .unwrap();
    }
    let native = OrientedHill48::new(
        plain
            .local_law()
            .clone()
            .with_voce_isotropic_hardening(Q, B)
            .unwrap()
            .with_linear_kinematic_hardening(C)
            .unwrap(),
        rotation,
    )
    .unwrap();
    let mut state = HillState::default();
    for strain in [
        [0.008, -0.003, -0.005, 0.002, 0.0, 0.001],
        [-0.006, 0.002, 0.004, -0.003, 0.002, 0.001],
        [0.01, -0.004, -0.006, 0.005, -0.002, 0.0],
    ] {
        let result = selected.response(&strain, &state).unwrap();
        assert_eq!(result, native.response(&strain, &state).unwrap());
        state = result.state;
    }
    let direct = resolve_hill48_voce(
        &fixture.card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
        0.04,
        rotation,
    )
    .unwrap();
    assert_eq!(direct.local_law().voce_saturation_pa(), Q);
    assert_eq!(direct.local_law().voce_rate_per_strain(), B);
    assert_eq!(direct.local_law().kinematic_modulus_pa(), 0.0);
}

#[test]
fn affine_q_and_prager_softening_through_zero_close_the_rotated_tensor_first_law() {
    let expansion = 1e-5;
    let mut values = voce_properties(expansion, true, Q, Q_SLOPE, B);
    values.push(property(C_PROPERTY, Pressure::DIMS, C, C_SLOPE));
    let fixture = Fixture::new(values);
    let rotation = fs_material::tensor::rotation([0.0, 0.0, 1.0], 0.37);
    let source = voce_source(&fixture, rotation)
        .with_linear_kinematic_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap();
    let reversed_order = fixture
        .source(rotation)
        .with_linear_kinematic_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
        .with_voce_isotropic_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap();
    for temperature in [280.0, 300.0, 320.0] {
        let endpoint = source
            .resolve_section_step(0.0, h(temperature), 0.1)
            .unwrap();
        let other = reversed_order
            .resolve_section_step(0.0, h(temperature), 0.1)
            .unwrap();
        assert_eq!(endpoint.law.resolved().law(), other.law.resolved().law());
        for name in [Q_PROPERTY, B_PROPERTY, C_PROPERTY] {
            let receipt = endpoint.law.resolved().material().property(name).unwrap();
            assert_eq!(
                receipt.answer().receipt.selected,
                source
                    .reference()
                    .material()
                    .property(name)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
            );
            fixture
                .card
                .claims()
                .verify_receipt(&receipt.answer().receipt)
                .unwrap();
        }
    }
    let old_strain = [0.008, -0.001, -0.002, 0.001, 0.0005, -0.0003];
    let new_strain = [0.009, -0.0008, -0.0025, 0.0018, 0.0006, -0.0001];
    let old = source
        .section_increment(
            0.0,
            0.0,
            &[0.0; 6],
            &HillState::default(),
            &old_strain,
            tolerance(),
            0.1,
        )
        .unwrap();
    let next = source
        .section_increment(
            0.0,
            h(320.0),
            &old_strain,
            &old.state,
            &new_strain,
            tolerance(),
            0.2,
        )
        .unwrap();
    let (e0, s0, base0, entropy0) =
        endpoint_oracle(300.0, &old_strain, &old.state, &rotation, expansion);
    let (_, s1, base1, entropy1) =
        endpoint_oracle(320.0, &new_strain, &next.state, &rotation, expansion);
    let shape0 = 0.5 * dot(&old.state.plastic_strain, &old.state.plastic_strain);
    let shape1 = 0.5 * dot(&next.state.plastic_strain, &next.state.plastic_strain);
    let psi0 = base0 + C * shape0 + Q * g(old.state.alpha);
    let entropy0 = entropy0 - C_SLOPE * shape0 - Q_SLOPE * g(old.state.alpha);
    // Q=C=0 at the hot endpoint, while both free-energy derivatives persist.
    let entropy1 = entropy1 - C_SLOPE * shape1 - Q_SLOPE * g(next.state.alpha);
    let delta_ts = 320.0 * entropy1 - 300.0 * entropy0;
    let coefficient = 0.5 * dot(&e0, &s0) * E_SLOPE * 20.0 / E
        + 0.5 * H_SLOPE * 20.0 * old.state.alpha.powi(2)
        - C * shape0
        - Q * g(old.state.alpha);
    let reversible = expansion * (s1[0] + s1[1] + s1[2]) * 20.0 - coefficient - delta_ts;
    near(next.stored_change_j_m3, base1 - psi0, 1e-5);
    near(next.coefficient_storage_change_j_m3, coefficient, 1e-5);
    near(next.thermoelastic_storage_change_j_m3, delta_ts, 1e-5);
    near(next.reversible_heat_j_m3, reversible, 1e-5);
    near(
        next.physical_heat_j_m3,
        (Y + 20.0 * Y_SLOPE) * (next.state.alpha - old.state.alpha),
        1e-6,
    );
    assert!(Q_SLOPE * g(next.state.alpha) < -1.0);
    assert!(next.numerical_loss_j_m3 >= 0.0);
    near(
        next.mechanical_work_j_m3
            - next.stored_change_j_m3
            - next.thermoelastic_storage_change_j_m3
            - next.physical_heat_j_m3
            - next.reversible_heat_j_m3
            - next.numerical_loss_j_m3,
        0.0,
        1e-5,
    );
}

#[test]
fn heated_native_section_keeps_saturation_history_and_reciprocal_heat() {
    let fixture = Fixture::new(voce_properties(0.0, false, Q, Q_SLOPE, B));
    let layers = [PlaneStressLayer {
        bottom_z_m: -0.005,
        top_z_m: 0.005,
        subdivisions: 1,
        law: voce_source(&fixture, IDENTITY),
    }];
    let reduction = ReductionConfig {
        max_iterations: 48,
        stress_tolerance_pa: 1e-4,
        strain_limit: 0.04,
    };
    let budget = SectionBudget {
        maximum_layers: 1,
        maximum_quadrature_points: 3,
    };
    let cold = [0.0];
    let hot = [h(320.0)];
    let loading = ThermoplasticSectionStep::new(
        &layers,
        &cold,
        &cold,
        0.1,
        ThicknessQuadrature::Gauss3,
        reduction,
        budget,
    )
    .unwrap();
    let strain = [0.0, 0.0, 0.012, 0.0, 0.0, 0.0];
    let virgin = loading.section().initial_state();
    let loaded = loading.section().evaluate(strain, &virgin).unwrap();
    loading
        .audit_energy(&virgin, &loaded.state, tolerance(), 1e-6)
        .unwrap();
    let heating = ThermoplasticSectionStep::new(
        &layers,
        &cold,
        &hot,
        0.2,
        ThicknessQuadrature::Gauss3,
        reduction,
        budget,
    )
    .unwrap();
    let heated = heating.section().evaluate(strain, &loaded.state).unwrap();
    let audit = heating
        .audit_energy(&loaded.state, &heated.state, tolerance(), 1e-6)
        .unwrap();
    assert!(heated.resultants[2] < loaded.resultants[2]);
    let mut coefficient = 0.0;
    let mut physical = 0.0;
    let mut delta_ts = 0.0;
    for (index, point) in heating.section().points().iter().enumerate() {
        let a0 = loaded.state.points[index].history.alpha;
        let a1 = heated.state.points[index].history.alpha;
        assert!(a1 > a0);
        coefficient -= point.weight_m * Q * g(a0);
        delta_ts += point.weight_m * (-Q_SLOPE) * (320.0 * g(a1) - 300.0 * g(a0));
        physical += point.weight_m * (Y + 20.0 * Y_SLOPE) * (a1 - a0);
    }
    near(audit.total.generalized_work_j_m2, 0.0, 0.0);
    near(audit.total.physical_heat_j_m2, physical, 1e-7);
    near(
        audit.total.coefficient_storage_change_j_m2,
        coefficient,
        1e-7,
    );
    near(
        audit.total.thermoelastic_storage_change_j_m2,
        delta_ts,
        1e-7,
    );
    near(
        audit.total.reversible_heat_j_m2,
        -coefficient - delta_ts,
        1e-7,
    );
    near(audit.total.balance_residual_j_m2, 0.0, 1e-6);
}

#[test]
fn zero_q_reference_retains_nonzero_entropy_slope_during_cooling() {
    let fixture = Fixture::new(voce_properties(0.0, false, Q, Q_SLOPE, B));
    let mut query = fixture.query(IDENTITY);
    query.context = QueryPoint::new().with("T", 320.0).unwrap();
    let source = ThermalHillSource::try_new(query, h(320.0))
        .unwrap()
        .with_voce_isotropic_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap();
    assert_eq!(source.reference().local_law().voce_saturation_pa(), 0.0);
    assert_eq!(source.reference().local_law().voce_rate_per_strain(), B);
    let strain = [0.0, 0.0, 0.0, 0.006, 0.0, 0.0];
    let loaded = source
        .section_increment(
            h(320.0),
            h(320.0),
            &[0.0; 6],
            &HillState::default(),
            &strain,
            tolerance(),
            0.1,
        )
        .unwrap();
    let cooled = source
        .section_increment(
            h(320.0),
            0.0,
            &strain,
            &loaded.state,
            &strain,
            tolerance(),
            0.2,
        )
        .unwrap();
    assert_eq!(cooled.state, loaded.state);
    let coefficient = Q * g(loaded.state.alpha);
    assert!(coefficient > 1.0);
    near(cooled.coefficient_storage_change_j_m3, coefficient, 1e-6);
    near(cooled.thermoelastic_storage_change_j_m3, -coefficient, 1e-6);
    near(cooled.reversible_heat_j_m3, 0.0, 1e-6);
    near(cooled.physical_heat_j_m3, 0.0, 0.0);
}

#[test]
fn tiny_exponent_thermal_work_survives_an_underflowing_unit_q_shape() {
    let q = 2e300;
    let q_slope = -1e299;
    let b = 1e-320;
    let fixture = Fixture::new(voce_properties(0.0, false, q, q_slope, b));
    let source = voce_source(&fixture, IDENTITY);
    let strain = [0.0, 0.0, 0.0, 0.006, 0.0, 0.0];
    let loaded = source
        .section_increment(
            0.0,
            0.0,
            &[0.0; 6],
            &HillState::default(),
            &strain,
            tolerance(),
            0.1,
        )
        .unwrap();
    let next = source
        .section_increment(
            0.0,
            h(320.0),
            &strain,
            &loaded.state,
            &strain,
            tolerance(),
            0.2,
        )
        .unwrap();
    let entropy = |alpha: f64| 0.5 * (-q_slope * b) * alpha * alpha;
    let coefficient = -0.5 * (q * b) * loaded.state.alpha.powi(2);
    let delta_ts = 320.0 * entropy(next.state.alpha) - 300.0 * entropy(loaded.state.alpha);
    assert!(next.coefficient_storage_change_j_m3 < 0.0);
    assert!(next.thermoelastic_storage_change_j_m3 > 0.0);
    near(next.coefficient_storage_change_j_m3, coefficient, 1e-34);
    near(next.thermoelastic_storage_change_j_m3, delta_ts, 1e-33);
    near(next.reversible_heat_j_m3, -coefficient - delta_ts, 1e-33);
}

#[test]
fn zero_saturation_is_exactly_the_original_thermal_path() {
    let fixture = Fixture::new(voce_properties(1e-5, true, 0.0, 0.0, B));
    let plain = fixture.source(QUARTER_TURN);
    let selected = voce_source(&fixture, QUARTER_TURN);
    let mut state = HillState::default();
    let mut old_strain = [0.0; 6];
    let mut old_h = 0.0;
    for (temperature, strain) in [
        (300.0, [0.008, -0.002, -0.004, 0.002, 0.0, 0.001]),
        (320.0, [0.009, -0.001, -0.003, 0.003, 0.001, 0.0]),
        (280.0, [-0.006, 0.002, 0.004, -0.003, 0.0, 0.0]),
    ] {
        let new_h = h(temperature);
        let a = plain
            .section_increment(old_h, new_h, &old_strain, &state, &strain, tolerance(), 0.1)
            .unwrap();
        let b = selected
            .section_increment(old_h, new_h, &old_strain, &state, &strain, tolerance(), 0.1)
            .unwrap();
        assert_eq!(a.state, b.state);
        assert_eq!(a.stress_pa, b.stress_pa);
        let account =
            |r: &fs_material::plastic::section::thermal::source::ThermalSectionIncrement<
                HillState,
            >| {
                [
                    r.physical_heat_j_m3,
                    r.creep_heat_j_m3,
                    r.reversible_heat_j_m3,
                    r.stored_change_j_m3,
                    r.coefficient_storage_change_j_m3,
                    r.thermoelastic_storage_change_j_m3,
                    r.numerical_loss_j_m3,
                    r.mechanical_work_j_m3,
                ]
                .map(f64::to_bits)
            };
        assert_eq!(account(&a), account(&b));
        state = a.state;
        old_strain = strain;
        old_h = new_h;
    }
}

#[test]
fn source_selection_rejects_missing_negative_kinked_or_varying_rate_data() {
    let missing = Fixture::new(properties(0.0, false));
    assert!(
        missing
            .source(IDENTITY)
            .with_voce_isotropic_hardening(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
    let negative = Fixture::new(voce_properties(0.0, false, Q, 2.0 * Q_SLOPE, B));
    assert!(
        negative
            .source(IDENTITY)
            .with_voce_isotropic_hardening(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
    for (name, dims, samples) in [
        (Q_PROPERTY, Pressure::DIMS, [2.0 * Q, 1.2 * Q, 0.0]),
        (B_PROPERTY, Dims([0; 6]), [B, B + 1.0, B]),
    ] {
        let mut values = voce_properties(0.0, false, Q, 0.0, B);
        let claim = values.iter_mut().find(|p| p.key.name() == name).unwrap();
        claim.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![
                (280.0, samples[0]),
                (300.0, samples[1]),
                (320.0, samples[2]),
            ],
            dims,
        };
        claim.interpolation = InterpolationPolicy::LinearInside;
        let fixture = Fixture::new(values);
        assert!(
            fixture
                .source(IDENTITY)
                .with_voce_isotropic_hardening(MaterialPropertySelection::SingleClaimOnly)
                .is_err()
        );
    }
    let mut wrong_roles = properties(0.0, false);
    wrong_roles.push(property(
        "voce_saturation_increment",
        Pressure::DIMS,
        Q,
        0.0,
    ));
    wrong_roles.push(property("voce_saturation_rate", Dims([0; 6]), B, 0.0));
    let wrong = Fixture::new(wrong_roles);
    assert!(
        wrong
            .source(IDENTITY)
            .with_voce_isotropic_hardening(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
    let fixture = Fixture::new(voce_properties(0.0, false, Q, 0.0, B));
    assert!(
        voce_source(&fixture, IDENTITY)
            .with_voce_isotropic_hardening(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
}
