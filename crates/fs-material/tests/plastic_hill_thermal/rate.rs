//! G0/G1/G3: sourced anisotropic overstress on physical thermal intervals.
//! All coefficients are authored fixtures, not qualified material data.

use super::*;
use fs_material::plastic::hill::source::{
    HILL_KINEMATIC_MODULUS_PROPERTY, HILL_OVERSTRESS_VISCOSITY_PROPERTY,
    HILL_VOCE_SATURATION_INCREMENT_PROPERTY, HILL_VOCE_SATURATION_RATE_PROPERTY,
};

const VISCOSITY_DIMS: Dims = Dims([-1, 1, -1, 0, 0, 0]);

fn viscosity_property() -> PropertyClaim {
    let mut value = property(HILL_OVERSTRESS_VISCOSITY_PROPERTY, VISCOSITY_DIMS, 1e8, 0.0);
    value.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(280.0, 4e8), (300.0, 1e8), (310.0, 2e7), (320.0, 1e7)],
        dims: VISCOSITY_DIMS,
    };
    value.interpolation = InterpolationPolicy::LinearInside;
    value
}

fn rate_properties(expansion: f64, varying_moduli: bool) -> Vec<PropertyClaim> {
    let mut values = properties(expansion, varying_moduli);
    values.push(viscosity_property());
    values
}

fn source_with_rate(fixture: &Fixture, orientation: [[f64; 3]; 3]) -> ThermalHillSource<'_> {
    fixture
        .source(orientation)
        .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
}

#[test]
fn sourced_shear_relaxation_and_reversal_match_the_independent_scalar_update() {
    let mut data = rate_properties(0.0, false);
    *data
        .iter_mut()
        .find(|p| p.key.name() == HILL_REFERENCE_YIELD_STRESS_PROPERTY)
        .unwrap() = property(HILL_REFERENCE_YIELD_STRESS_PROPERTY, Pressure::DIMS, Y, 0.0);
    let fixture = Fixture::new(data);
    let source = source_with_rate(&fixture, IDENTITY);
    let twice_shear = E / (1.0 + NU);
    // For the authored N=2, q=2*|tau_xy| and dp_xy=sign(tau)*dalpha.
    // Therefore q_trial - Y - H*alpha_old =
    // (2*(2*mu) + H + eta/dt)*dalpha while flow is active.
    let a = 2.0 * twice_shear + H;
    let mut old_strain = [0.0; 6];
    let mut old = HillState::default();
    let mut old_temperature = 300.0;
    for (shear, temperature, dt, eta) in [
        (0.004, 300.0, 0.05, 1e8),
        (0.004, 310.0, 0.10, 2e7),
        (0.004, 300.0, 0.025, 1e8),
        (-0.004, 305.0, 0.05, 6e7),
        (-0.004, 300.0, 0.20, 1e8),
    ] {
        let mut strain = [0.0; 6];
        strain[3] = shear;
        let trial_shear = twice_shear * (shear - old.plastic_strain[3]);
        let excess = 2.0 * trial_shear.abs() - Y - H * old.alpha;
        let increment = excess.max(0.0) / (a + eta / dt);
        let plastic = old.plastic_strain[3] + trial_shear.signum() * increment;
        let alpha = old.alpha + increment;
        let shear_stress = twice_shear * (shear - plastic);
        let initial_storage = twice_shear * (old_strain[3] - old.plastic_strain[3]).powi(2)
            + 0.5 * H * old.alpha.powi(2);
        let final_storage = twice_shear * (shear - plastic).powi(2) + 0.5 * H * alpha.powi(2);
        let physical_heat = (Y + eta / dt * increment) * increment;
        let elastic_change = shear - plastic - old_strain[3] + old.plastic_strain[3];
        let numerical = twice_shear * elastic_change.powi(2) + 0.5 * H * increment.powi(2);
        let saved = old.clone();
        let endpoint = source
            .resolve_section_step(h(old_temperature), h(temperature), dt)
            .unwrap();
        let rate = endpoint.law.linear_overstress().unwrap();
        assert_eq!(rate.duration_s(), dt);
        near(rate.viscosity_pa_s(), eta, 1e-6);
        let next = source
            .section_increment(
                h(old_temperature),
                h(temperature),
                &old_strain,
                &old,
                &strain,
                tolerance(),
                dt,
            )
            .unwrap();
        near(next.state.alpha, alpha, 2e-15);
        near(next.state.plastic_strain[3], plastic, 2e-15);
        near(next.stress_pa[3], shear_stress, 1e-5);
        near(
            next.stored_change_j_m3,
            final_storage - initial_storage,
            1e-5,
        );
        near(next.physical_heat_j_m3, physical_heat, 1e-5);
        near(next.numerical_loss_j_m3, numerical, 1e-5);
        near(
            next.mechanical_work_j_m3,
            2.0 * shear_stress * (shear - old_strain[3]),
            1e-5,
        );
        assert_eq!(next.coefficient_storage_change_j_m3, 0.0);
        assert_eq!(next.thermoelastic_storage_change_j_m3, 0.0);
        assert_eq!(next.reversible_heat_j_m3, 0.0);
        assert_eq!(next.creep_heat_j_m3, 0.0);
        assert_eq!(next.state, endpoint.law.update_state(&strain, &old));
        assert_eq!(next.stress_pa, endpoint.law.stress(&strain, &old));
        source
            .validate_section_history(h(temperature), &strain, &next.state, 1e-3)
            .unwrap();
        assert!(
            endpoint
                .law
                .resolved()
                .law()
                .validate_committed_history(&strain, &next.state, 1e-3)
                .is_err(),
            "a finite-rate endpoint retains positive physical overstress"
        );
        assert_eq!(old, saved);
        old = next.state;
        old_strain = strain;
        old_temperature = temperature;
    }
}

#[test]
fn thermal_rate_with_rotated_voce_and_prager_history_closes_the_complete_first_law() {
    const C: f64 = 2e8;
    const C_SLOPE: f64 = -2e6;
    const Q: f64 = 2e6;
    const Q_SLOPE: f64 = -2e4;
    const B: f64 = 180.0;
    let expansion = 1e-5;
    let mut data = rate_properties(expansion, true);
    data.extend([
        property(HILL_KINEMATIC_MODULUS_PROPERTY, Pressure::DIMS, C, C_SLOPE),
        property(
            HILL_VOCE_SATURATION_INCREMENT_PROPERTY,
            Pressure::DIMS,
            Q,
            Q_SLOPE,
        ),
        property(HILL_VOCE_SATURATION_RATE_PROPERTY, Dims([0; 6]), B, 0.0),
    ]);
    let fixture = Fixture::new(data);
    let rotation = fs_material::tensor::rotation([0.0, 0.0, 1.0], 0.37);
    let rate_first = source_with_rate(&fixture, rotation)
        .with_linear_kinematic_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
        .with_voce_isotropic_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap();
    let rate_last = fixture
        .source(rotation)
        .with_voce_isotropic_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
        .with_linear_kinematic_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
        .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
        .unwrap();
    let thermodynamics = |temperature: f64, strain: &Voigt, state: &HillState| {
        let (_, stress, mut free, mut entropy) =
            endpoint_oracle(temperature, strain, state, &rotation, expansion);
        let norm_squared = dot(&state.plastic_strain, &state.plastic_strain);
        let shape = state.alpha + (-B * state.alpha).exp_m1() / B;
        free += 0.5 * (C + C_SLOPE * (temperature - 300.0)) * norm_squared
            + (Q + Q_SLOPE * (temperature - 300.0)) * shape;
        entropy -= 0.5 * C_SLOPE * norm_squared + Q_SLOPE * shape;
        (stress, free, entropy)
    };
    let mut old = HillState::default();
    let mut old_strain = [0.0; 6];
    let mut old_temperature = 300.0;
    for (temperature, dt, strain) in [
        (300.0, 0.05, [0.008, -0.002, -0.003, 0.001, 0.0005, -0.0003]),
        (315.0, 0.2, [0.008, -0.002, -0.003, 0.001, 0.0005, -0.0003]),
        (303.0, 0.1, [-0.006, 0.001, 0.002, -0.001, -0.0003, 0.0002]),
    ] {
        let next = rate_first
            .section_increment(
                h(old_temperature),
                h(temperature),
                &old_strain,
                &old,
                &strain,
                tolerance(),
                dt,
            )
            .unwrap();
        let reordered = rate_last
            .section_increment(
                h(old_temperature),
                h(temperature),
                &old_strain,
                &old,
                &strain,
                tolerance(),
                dt,
            )
            .unwrap();
        assert_eq!(next.state, reordered.state);
        assert_eq!(next.stress_pa, reordered.stress_pa);
        assert_eq!(next.physical_heat_j_m3, reordered.physical_heat_j_m3);
        let (_, psi0, entropy0) = thermodynamics(old_temperature, &old_strain, &old);
        let (stress, psi1, entropy1) = thermodynamics(temperature, &strain, &next.state);
        let delta_ts = temperature * entropy1 - old_temperature * entropy0;
        let delta_strain = core::array::from_fn(|i| strain[i] - old_strain[i]);
        let work = dot(&stress, &delta_strain);
        let endpoint = rate_first
            .resolve_section_step(h(old_temperature), h(temperature), dt)
            .unwrap();
        let eta = endpoint.law.linear_overstress().unwrap().viscosity_pa_s();
        let da = next.state.alpha - old.alpha;
        let expected_heat = (Y + Y_SLOPE * (temperature - 300.0) + eta / dt * da) * da;
        near(next.stored_change_j_m3, psi1 - psi0, 2e-5);
        near(next.thermoelastic_storage_change_j_m3, delta_ts, 2e-5);
        near(next.physical_heat_j_m3, expected_heat, 1e-5);
        near(next.mechanical_work_j_m3, work, 1e-5);
        near(
            work - (psi1 - psi0)
                - delta_ts
                - expected_heat
                - next.reversible_heat_j_m3
                - next.numerical_loss_j_m3,
            0.0,
            3e-5,
        );
        for property in endpoint.law.resolved().material().properties() {
            fixture
                .card
                .claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
        let pin = endpoint
            .law
            .resolved()
            .material()
            .property(HILL_OVERSTRESS_VISCOSITY_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .selected;
        assert_eq!(
            pin,
            rate_first
                .reference()
                .material()
                .property(HILL_OVERSTRESS_VISCOSITY_PROPERTY)
                .unwrap()
                .answer()
                .receipt
                .selected
        );
        rate_first
            .validate_section_history(h(temperature), &strain, &next.state, 1e-3)
            .unwrap();
        old = next.state;
        old_strain = strain;
        old_temperature = temperature;
    }
}

#[test]
fn existing_thermal_section_retains_rate_history_and_relaxes_during_a_hold() {
    let fixture = Fixture::new(rate_properties(0.0, false));
    let layers = [PlaneStressLayer {
        bottom_z_m: -0.005,
        top_z_m: 0.005,
        subdivisions: 1,
        law: source_with_rate(&fixture, QUARTER_TURN),
    }];
    let reduction = ReductionConfig {
        max_iterations: 64,
        stress_tolerance_pa: 1e-4,
        strain_limit: 0.04,
    };
    let budget = SectionBudget {
        maximum_layers: 2,
        maximum_quadrature_points: 8,
    };
    let loading = ThermoplasticSectionStep::new(
        &layers,
        &[0.0],
        &[0.0],
        0.05,
        ThicknessQuadrature::Gauss3,
        reduction,
        budget,
    )
    .unwrap();
    let initial = loading.section().initial_state();
    let strain = [0.008, 0.0, 0.0, 0.0, 0.0, 0.0];
    let loaded = loading.section().evaluate(strain, &initial).unwrap();
    loading
        .audit_energy(&initial, &loaded.state, tolerance(), 1e-6)
        .unwrap();
    let holding = ThermoplasticSectionStep::new(
        &layers,
        &[0.0],
        &[0.0],
        0.2,
        ThicknessQuadrature::Gauss3,
        reduction,
        budget,
    )
    .unwrap();
    let saved = loaded.state.clone();
    let held = holding.section().evaluate(strain, &loaded.state).unwrap();
    let account = holding
        .audit_energy(&loaded.state, &held.state, tolerance(), 1e-6)
        .unwrap();
    assert!(held.resultants[0] < loaded.resultants[0]);
    assert!(
        held.state
            .points
            .iter()
            .zip(&loaded.state.points)
            .all(|(new, old)| new.history.alpha > old.history.alpha)
    );
    assert_eq!(account.total.generalized_work_j_m2, 0.0);
    assert!(account.total.physical_heat_j_m2 > 0.0);
    near(account.total.balance_residual_j_m2, 0.0, 1e-6);
    assert_eq!(
        holding.section().evaluate(strain, &loaded.state).unwrap(),
        held
    );
    assert_eq!(loaded.state, saved);
}

#[test]
fn viscosity_requires_an_original_positive_source_path_and_a_real_interval() {
    let missing = Fixture::new(properties(0.0, false));
    assert!(
        missing
            .source(IDENTITY)
            .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
    for case in 0..3 {
        let mut data = rate_properties(0.0, false);
        let property = data
            .iter_mut()
            .find(|p| p.key.name() == HILL_OVERSTRESS_VISCOSITY_PROPERTY)
            .unwrap();
        match case {
            0 => {
                property.key = PropertyKey::new(HILL_OVERSTRESS_VISCOSITY_PROPERTY, Pressure::DIMS);
                property.value = PropertyValue::Scalar {
                    value: 1e8,
                    dims: Pressure::DIMS,
                };
                property.interpolation = InterpolationPolicy::ConstantWithinValidity;
            }
            1 => {
                property.value = PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, 4e8), (290.0, -1.0), (300.0, 1e8), (320.0, 1e7)],
                    dims: VISCOSITY_DIMS,
                }
            }
            _ => property.validity = ValidityDomain::unconstrained().with("T", 295.0, 305.0),
        }
        let fixture = Fixture::new(data);
        assert!(
            fixture
                .source(IDENTITY)
                .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
                .is_err(),
            "unsupported viscosity case {case}"
        );
    }
    let fixture = Fixture::new(rate_properties(0.0, false));
    let source = source_with_rate(&fixture, IDENTITY);
    assert!(
        source
            .clone()
            .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
    for dt in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(source.resolve_section_step(0.0, 0.0, dt).is_err());
    }
    let invalid = HillState {
        plastic_strain: [0.001, 0.0, 0.0, 0.0, 0.0, 0.0],
        alpha: 0.1,
    };
    assert!(
        source
            .validate_section_history(0.0, &[0.0; 6], &invalid, 1e-3)
            .is_err()
    );
    assert!(
        source
            .validate_section_history(0.0, &[0.0; 6], &HillState::default(), 0.0)
            .is_err()
    );
}
