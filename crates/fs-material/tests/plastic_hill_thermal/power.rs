//! G0/G1/G3: nonlinear sourced Hill kinetics at actual thermal intervals.
//! Authored coefficients specify mathematical fixtures, not calibrated metals.
use super::*;
use fs_material::plastic::hill::source::{
    HILL_KINEMATIC_MODULUS_PROPERTY, HILL_OVERSTRESS_EXPONENT_PROPERTY,
    HILL_OVERSTRESS_STRESS_SCALE_PROPERTY, HILL_OVERSTRESS_VISCOSITY_PROPERTY,
    HILL_VOCE_SATURATION_INCREMENT_PROPERTY, HILL_VOCE_SATURATION_RATE_PROPERTY,
};

const ETA_DIMS: Dims = Dims([-1, 1, -1, 0, 0, 0]);

fn kinetic_curve(name: &str, dims: Dims, knots: Vec<(f64, f64)>) -> PropertyClaim {
    let mut claim = property(name, dims, knots[0].1, 0.0);
    claim.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots,
        dims,
    };
    claim.interpolation = InterpolationPolicy::LinearInside;
    claim
}

fn power_properties(expansion: f64, varying_moduli: bool) -> Vec<PropertyClaim> {
    let mut data = properties(expansion, varying_moduli);
    data.extend([
        kinetic_curve(
            HILL_OVERSTRESS_VISCOSITY_PROPERTY,
            ETA_DIMS,
            vec![(280.0, 4e8), (300.0, 1e8), (310.0, 2e7), (320.0, 1e7)],
        ),
        kinetic_curve(
            HILL_OVERSTRESS_STRESS_SCALE_PROPERTY,
            Pressure::DIMS,
            vec![(280.0, 2e6), (300.0, 1e6), (310.0, 5e5), (320.0, 7.5e5)],
        ),
        kinetic_curve(
            HILL_OVERSTRESS_EXPONENT_PROPERTY,
            Dims([0; 6]),
            vec![(280.0, 1.5), (300.0, 2.0), (310.0, 4.0), (320.0, 3.0)],
        ),
    ]);
    data
}

fn power_source(fixture: &Fixture, axes: [[f64; 3]; 3]) -> ThermalHillSource<'_> {
    fixture
        .source(axes)
        .with_power_law_overstress(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
}

#[test]
fn sourced_power_law_matches_scalar_loading_relaxation_and_reversal() {
    let mut data = power_properties(0.0, false);
    *data
        .iter_mut()
        .find(|p| p.key.name() == HILL_REFERENCE_YIELD_STRESS_PROPERTY)
        .unwrap() = property(HILL_REFERENCE_YIELD_STRESS_PROPERTY, Pressure::DIMS, Y, 0.0);
    let fixture = Fixture::new(data);
    let source = power_source(&fixture, IDENTITY);
    let twice_shear = E / (1.0 + NU);
    let a = 2.0 * twice_shear + H;
    let (mut old, mut old_strain, mut old_temperature) = (HillState::default(), [0.0; 6], 300.0);
    for (shear, temperature, dt, eta, scale, exponent) in [
        (0.004, 300.0, 0.05, 1e8, 1e6, 2.0),
        (0.004, 310.0, 0.10, 2e7, 5e5, 4.0),
        (0.004, 300.0, 0.025, 1e8, 1e6, 2.0),
        (-0.004, 305.0, 0.05, 6e7, 7.5e5, 3.0),
        (-0.004, 300.0, 0.20, 1e8, 1e6, 2.0),
    ] {
        let strain = [0.0, 0.0, 0.0, shear, 0.0, 0.0];
        // N=2 gives q=2|tau| and dp_xy=sign(tau)*dalpha. The single
        // scalar balance below is independent of the five-mode return map.
        let trial = twice_shear * (shear - old.plastic_strain[3]);
        let excess = 2.0 * trial.abs() - Y - H * old.alpha;
        let resistance = |da: f64| scale * (eta * da / (dt * scale)).powf(1.0 / exponent);
        let (mut lo, mut hi) = (0.0, excess.max(0.0) / a);
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if a * mid + resistance(mid) < excess {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let da = 0.5 * (lo + hi);
        assert!(da > 0.0);
        let alpha = old.alpha + da;
        let plastic = old.plastic_strain[3] + trial.signum() * da;
        let stress = twice_shear * (shear - plastic);
        let heat = (Y + resistance(da)) * da;
        let initial_storage = twice_shear * (old_strain[3] - old.plastic_strain[3]).powi(2)
            + 0.5 * H * old.alpha.powi(2);
        let final_storage = twice_shear * (shear - plastic).powi(2) + 0.5 * H * alpha.powi(2);
        let de = shear - plastic - old_strain[3] + old.plastic_strain[3];
        let numerical = twice_shear * de.powi(2) + 0.5 * H * da.powi(2);
        let endpoint = source
            .resolve_section_step(h(old_temperature), h(temperature), dt)
            .unwrap();
        assert!(endpoint.law.linear_overstress().is_none());
        let kinetics = endpoint.law.power_law_overstress().unwrap();
        near(kinetics.viscosity_pa_s(), eta, 1e-6);
        near(kinetics.stress_scale_pa(), scale, 1e-8);
        near(kinetics.exponent(), exponent, 1e-14);
        assert_eq!(kinetics.duration_s(), dt);
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
        near(next.state.alpha, alpha, 5e-14);
        near(next.state.plastic_strain[3], plastic, 5e-14);
        near(next.stress_pa[3], stress, 1e-4);
        near(next.physical_heat_j_m3, heat, 1e-5);
        near(
            next.stored_change_j_m3,
            final_storage - initial_storage,
            1e-5,
        );
        near(next.numerical_loss_j_m3, numerical, 1e-5);
        near(
            next.mechanical_work_j_m3,
            2.0 * stress * (shear - old_strain[3]),
            1e-5,
        );
        assert_eq!(next.reversible_heat_j_m3, 0.0);
        assert_eq!(next.thermoelastic_storage_change_j_m3, 0.0);
        assert_eq!(next.coefficient_storage_change_j_m3, 0.0);
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
                .is_err()
        );
        old = next.state;
        old_strain = strain;
        old_temperature = temperature;
    }
}

#[test]
fn nonlinear_thermal_kinetics_keep_rotated_hardening_storage_and_source_pins() {
    let (c, c_slope, q, q_slope, b) = (2e8, -2e6, 2e6, -2e4, 180.0);
    let expansion = 1e-5;
    let mut data = power_properties(expansion, true);
    data.extend([
        property(HILL_KINEMATIC_MODULUS_PROPERTY, Pressure::DIMS, c, c_slope),
        property(
            HILL_VOCE_SATURATION_INCREMENT_PROPERTY,
            Pressure::DIMS,
            q,
            q_slope,
        ),
        property(HILL_VOCE_SATURATION_RATE_PROPERTY, Dims([0; 6]), b, 0.0),
    ]);
    let fixture = Fixture::new(data);
    let axes = fs_material::tensor::rotation([0.0, 0.0, 1.0], 0.37);
    let first = power_source(&fixture, axes)
        .with_linear_kinematic_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
        .with_voce_isotropic_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap();
    let last = fixture
        .source(axes)
        .with_voce_isotropic_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
        .with_linear_kinematic_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
        .with_power_law_overstress(MaterialPropertySelection::SingleClaimOnly)
        .unwrap();
    let thermodynamics = |temperature: f64, strain: &Voigt, state: &HillState| {
        let (_, stress, mut free, mut entropy) =
            endpoint_oracle(temperature, strain, state, &axes, expansion);
        let square = dot(&state.plastic_strain, &state.plastic_strain);
        let shape = state.alpha + (-b * state.alpha).exp_m1() / b;
        free += 0.5 * (c + c_slope * (temperature - 300.0)) * square
            + (q + q_slope * (temperature - 300.0)) * shape;
        entropy -= 0.5 * c_slope * square + q_slope * shape;
        (stress, free, entropy)
    };
    let (mut old, mut old_strain, mut old_t) = (HillState::default(), [0.0; 6], 300.0);
    for (t, dt, strain) in [
        (300.0, 0.05, [0.008, -0.002, -0.003, 0.001, 0.0005, -0.0003]),
        (315.0, 0.20, [0.008, -0.002, -0.003, 0.001, 0.0005, -0.0003]),
        (303.0, 0.10, [-0.006, 0.001, 0.002, -0.001, -0.0003, 0.0002]),
    ] {
        let next = first
            .section_increment(h(old_t), h(t), &old_strain, &old, &strain, tolerance(), dt)
            .unwrap();
        let other = last
            .section_increment(h(old_t), h(t), &old_strain, &old, &strain, tolerance(), dt)
            .unwrap();
        assert_eq!(next.state, other.state);
        assert_eq!(next.stress_pa, other.stress_pa);
        assert_eq!(next.physical_heat_j_m3, other.physical_heat_j_m3);
        let (_, psi0, s0) = thermodynamics(old_t, &old_strain, &old);
        let (stress, psi1, s1) = thermodynamics(t, &strain, &next.state);
        let delta_ts = t * s1 - old_t * s0;
        let work = dot(
            &stress,
            &core::array::from_fn(|i| strain[i] - old_strain[i]),
        );
        let endpoint = first.resolve_section_step(h(old_t), h(t), dt).unwrap();
        let law = endpoint.law.power_law_overstress().unwrap();
        let da = next.state.alpha - old.alpha;
        let resistance = law.stress_scale_pa()
            * (law.viscosity_pa_s() * da / (dt * law.stress_scale_pa())).powf(1.0 / law.exponent());
        let heat = (Y + Y_SLOPE * (t - 300.0) + resistance) * da;
        near(next.stored_change_j_m3, psi1 - psi0, 3e-5);
        near(next.thermoelastic_storage_change_j_m3, delta_ts, 3e-5);
        near(next.physical_heat_j_m3, heat, 2e-5);
        near(next.mechanical_work_j_m3, work, 1e-5);
        near(
            work - (psi1 - psi0)
                - delta_ts
                - heat
                - next.reversible_heat_j_m3
                - next.numerical_loss_j_m3,
            0.0,
            4e-5,
        );
        for name in [
            HILL_OVERSTRESS_VISCOSITY_PROPERTY,
            HILL_OVERSTRESS_STRESS_SCALE_PROPERTY,
            HILL_OVERSTRESS_EXPONENT_PROPERTY,
        ] {
            let property = endpoint.law.resolved().material().property(name).unwrap();
            fixture
                .card
                .claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                first
                    .reference()
                    .material()
                    .property(name)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
        old = next.state;
        old_strain = strain;
        old_t = t;
    }
}

#[test]
fn unit_exponent_reproduces_linear_thermal_history_and_heat_exactly() {
    let mut data = power_properties(0.0, false);
    *data
        .iter_mut()
        .find(|p| p.key.name() == HILL_OVERSTRESS_EXPONENT_PROPERTY)
        .unwrap() = property(HILL_OVERSTRESS_EXPONENT_PROPERTY, Dims([0; 6]), 1.0, 0.0);
    let fixture = Fixture::new(data);
    let linear = fixture
        .source(QUARTER_TURN)
        .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
        .unwrap();
    let power = power_source(&fixture, QUARTER_TURN);
    let (mut old, mut old_strain, mut old_t) = (HillState::default(), [0.0; 6], 300.0);
    for (t, dt, strain) in [
        (300.0, 0.05, [0.008, -0.002, -0.003, 0.001, 0.0, 0.0]),
        (310.0, 0.20, [0.008, -0.002, -0.003, 0.001, 0.0, 0.0]),
        (300.0, 0.10, [-0.006, 0.001, 0.002, -0.001, 0.0, 0.0]),
    ] {
        let a = linear
            .section_increment(h(old_t), h(t), &old_strain, &old, &strain, tolerance(), dt)
            .unwrap();
        let b = power
            .section_increment(h(old_t), h(t), &old_strain, &old, &strain, tolerance(), dt)
            .unwrap();
        assert_eq!(a.state, b.state);
        assert_eq!(a.stress_pa, b.stress_pa);
        assert_eq!(a.physical_heat_j_m3, b.physical_heat_j_m3);
        assert_eq!(a.stored_change_j_m3, b.stored_change_j_m3);
        assert_eq!(a.reversible_heat_j_m3, b.reversible_heat_j_m3);
        assert_eq!(a.numerical_loss_j_m3, b.numerical_loss_j_m3);
        old = a.state;
        old_strain = strain;
        old_t = t;
    }
}

#[test]
fn power_kinetics_require_complete_original_paths_and_one_exclusive_selection() {
    for case in 0..4 {
        let mut data = power_properties(0.0, false);
        match case {
            0 => data.retain(|p| p.key.name() != HILL_OVERSTRESS_STRESS_SCALE_PROPERTY),
            1 => {
                *data
                    .iter_mut()
                    .find(|p| p.key.name() == HILL_OVERSTRESS_STRESS_SCALE_PROPERTY)
                    .unwrap() = property(HILL_OVERSTRESS_STRESS_SCALE_PROPERTY, ETA_DIMS, 1e6, 0.0);
            }
            2 => {
                *data
                    .iter_mut()
                    .find(|p| p.key.name() == HILL_OVERSTRESS_EXPONENT_PROPERTY)
                    .unwrap() = kinetic_curve(
                    HILL_OVERSTRESS_EXPONENT_PROPERTY,
                    Dims([0; 6]),
                    vec![(280.0, 2.0), (290.0, 0.9), (300.0, 2.0), (320.0, 4.0)],
                );
            }
            _ => {
                data.iter_mut()
                    .find(|p| p.key.name() == HILL_OVERSTRESS_STRESS_SCALE_PROPERTY)
                    .unwrap()
                    .validity = ValidityDomain::unconstrained().with("T", 295.0, 305.0);
            }
        }
        let fixture = Fixture::new(data);
        assert!(
            fixture
                .source(IDENTITY)
                .with_power_law_overstress(MaterialPropertySelection::SingleClaimOnly)
                .is_err(),
            "unsupported kinetic source case {case}"
        );
    }
    let fixture = Fixture::new(power_properties(0.0, false));
    let selected = power_source(&fixture, IDENTITY);
    assert!(
        selected
            .clone()
            .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
    assert!(
        selected
            .with_power_law_overstress(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
    assert!(
        fixture
            .source(IDENTITY)
            .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
            .unwrap()
            .with_power_law_overstress(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
}
