//! G0/G1 thermal mixed hardening through original Hill sources and sections.

use super::*;
use fs_material::plastic::hill::source::HILL_KINEMATIC_MODULUS_PROPERTY;

const C: f64 = 4e8;
const C_SLOPE: f64 = -2e7;

fn mixed_properties(expansion: f64, varying_moduli: bool, slope: f64) -> Vec<PropertyClaim> {
    let mut values = properties(expansion, varying_moduli);
    values.push(property(
        HILL_KINEMATIC_MODULUS_PROPERTY,
        Pressure::DIMS,
        C,
        slope,
    ));
    values
}

fn mixed_source(fixture: &Fixture, rotation: [[f64; 3]; 3]) -> ThermalHillSource<'_> {
    fixture
        .source(rotation)
        .with_linear_kinematic_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
}

#[test]
fn cyclic_thermal_section_matches_shear_bauschinger_threshold_and_stored_backstress_work() {
    let fixture = Fixture::new(mixed_properties(0.0, false, 0.0));
    let thickness = 0.01;
    let layers = [PlaneStressLayer {
        bottom_z_m: -0.5 * thickness,
        top_z_m: 0.5 * thickness,
        subdivisions: 1,
        law: mixed_source(&fixture, IDENTITY),
    }];
    let step = ThermoplasticSectionStep::new(
        &layers,
        &[0.0],
        &[0.0],
        0.1,
        ThicknessQuadrature::Gauss3,
        ReductionConfig {
            max_iterations: 48,
            stress_tolerance_pa: 1e-4,
            strain_limit: 0.04,
        },
        SectionBudget {
            maximum_layers: 1,
            maximum_quadrature_points: 3,
        },
    )
    .unwrap();
    let mut accepted = step.section().initial_state();
    let twice_shear = E / (1.0 + NU);
    let n = 2.0_f64;
    let b = (2.0 * n).sqrt();
    let flow = 0.5 * b;
    let peak_strain = 0.006;
    let peak_alpha = (b * twice_shear * peak_strain - Y) / (n * (twice_shear + C) + H);
    let peak_plastic = flow * peak_alpha;
    let reverse_yield = C * peak_plastic - (Y + H * peak_alpha) / b;
    // The translated surface yields in reverse before reaching the magnitude
    // of the virgin reverse yield stress. The section must preserve that memory.
    assert!(reverse_yield < 0.0 && reverse_yield > -Y / b);
    let before_reverse = peak_plastic + (reverse_yield + 0.01 * Y / b) / twice_shear;
    let past_reverse = peak_plastic + (reverse_yield - 0.1 * Y / b) / twice_shear;
    let (mut old_strain, mut plastic, mut alpha, mut old_storage) = (0.0, 0.0, 0.0, 0.0);
    for (interval, strain) in [peak_strain, before_reverse, past_reverse, -0.006, 0.006]
        .into_iter()
        .enumerate()
    {
        // q=b|tau-Cp|; dp=sign(tau-Cp)*(b/2)*dalpha in tensor shear.
        let relative_trial = twice_shear * strain - (twice_shear + C) * plastic;
        let increment =
            (b * relative_trial.abs() - Y - H * alpha).max(0.0) / (n * (twice_shear + C) + H);
        let delta_plastic = relative_trial.signum() * flow * increment;
        plastic += delta_plastic;
        alpha += increment;
        let stress = twice_shear * (strain - plastic);
        let storage = twice_shear * (strain - plastic).powi(2)
            + 0.5 * H * alpha.powi(2)
            + C * plastic.powi(2);
        let numerical = twice_shear * (strain - old_strain - delta_plastic).powi(2)
            + 0.5 * H * increment.powi(2)
            + C * delta_plastic.powi(2);
        let next = step
            .section()
            .evaluate([0.0, 0.0, 2.0 * strain, 0.0, 0.0, 0.0], &accepted)
            .unwrap();
        let audit = step
            .audit_energy(&accepted, &next.state, tolerance(), 1e-6)
            .unwrap();
        near(next.resultants[2], thickness * stress, 1e-5);
        for point in &next.state.points {
            near(point.history.alpha, alpha, 2e-12);
            near(point.history.plastic_strain[3], plastic, 2e-12);
            for i in [0, 1, 2, 4, 5] {
                near(point.history.plastic_strain[i], 0.0, 1e-15);
            }
        }
        near(
            audit.total.physical_heat_j_m2,
            thickness * Y * increment,
            1e-7,
        );
        near(
            audit.total.stored_change_j_m2,
            thickness * (storage - old_storage),
            1e-7,
        );
        near(audit.total.numerical_loss_j_m2, thickness * numerical, 1e-7);
        near(
            audit.total.generalized_work_j_m2,
            thickness * 2.0 * stress * (strain - old_strain),
            1e-7,
        );
        near(audit.total.reversible_heat_j_m2, 0.0, 0.0);
        near(audit.total.balance_residual_j_m2, 0.0, 1e-6);
        if interval == 1 {
            assert_eq!(next.state.points[0].history, accepted.points[0].history);
            near(increment, 0.0, 0.0);
        } else if interval == 2 {
            assert!(next.state.points[0].history.alpha > accepted.points[0].history.alpha);
            assert!(
                next.state.points[0].history.plastic_strain[3]
                    < accepted.points[0].history.plastic_strain[3]
            );
        }
        old_storage = storage;
        old_strain = strain;
        accepted = next.state;
    }
}

fn mixed_oracle(
    temperature: f64,
    total: &Voigt,
    state: &HillState,
    rotation: &[[f64; 3]; 3],
    expansion: f64,
) -> (Voigt, Voigt, f64, f64) {
    let (elastic, stress, energy, entropy) =
        endpoint_oracle(temperature, total, state, rotation, expansion);
    let shape = 0.5 * dot(&state.plastic_strain, &state.plastic_strain);
    let modulus = C + C_SLOPE * (temperature - 300.0);
    (
        elastic,
        stress,
        energy + modulus * shape,
        entropy - C_SLOPE * shape,
    )
}

#[test]
fn affine_kinematic_softening_through_zero_keeps_tensor_storage_entropy_and_original_history() {
    let expansion = 1e-5;
    let fixture = Fixture::new(mixed_properties(expansion, true, C_SLOPE));
    let rotation = fs_material::tensor::rotation([0.0, 0.0, 1.0], 0.37);
    let source = mixed_source(&fixture, rotation);
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
    let saved = old.state.clone();
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
    let (elastic0, stress0, psi0, entropy0) =
        mixed_oracle(300.0, &old_strain, &old.state, &rotation, expansion);
    let (_, stress1, psi1, entropy1) =
        mixed_oracle(320.0, &new_strain, &next.state, &rotation, expansion);
    let shape = 0.5 * dot(&old.state.plastic_strain, &old.state.plastic_strain);
    let shear_shape: f64 = old.state.plastic_strain[3..].iter().map(|v| v * v).sum();
    assert!(
        C * shear_shape > 1.0,
        "the oracle must exercise tensor shear storage"
    );
    let coefficient = 0.5 * dot(&elastic0, &stress0) * E_SLOPE * 20.0 / E
        + 0.5 * H_SLOPE * 20.0 * old.state.alpha.powi(2)
        - C * shape;
    let delta_ts = 320.0 * entropy1 - 300.0 * entropy0;
    let reversible =
        expansion * (stress1[0] + stress1[1] + stress1[2]) * 20.0 - coefficient - delta_ts;
    let delta: Voigt = core::array::from_fn(|i| new_strain[i] - old_strain[i]);
    let work = dot(&stress1, &delta);
    for i in 0..6 {
        near(next.stress_pa[i], stress1[i], 1e-3);
    }
    near(next.stored_change_j_m3, psi1 - psi0, 1e-5);
    near(next.coefficient_storage_change_j_m3, coefficient, 1e-5);
    near(next.thermoelastic_storage_change_j_m3, delta_ts, 1e-5);
    near(next.reversible_heat_j_m3, reversible, 1e-5);
    near(next.mechanical_work_j_m3, work, 1e-5);
    near(
        next.physical_heat_j_m3,
        (Y + 20.0 * Y_SLOPE) * (next.state.alpha - old.state.alpha),
        1e-5,
    );
    near(
        work - (psi1 - psi0)
            - delta_ts
            - next.physical_heat_j_m3
            - reversible
            - next.numerical_loss_j_m3,
        0.0,
        1e-5,
    );
    let endpoint = source.resolve_section_step(0.0, h(320.0), 0.2).unwrap();
    assert_eq!(
        endpoint.law.resolved().local_law().kinematic_modulus_pa(),
        0.0
    );
    let mechanical: Voigt =
        core::array::from_fn(|i| new_strain[i] - if i < 3 { expansion * 20.0 } else { 0.0 });
    near(
        endpoint
            .law
            .resolved()
            .law()
            .stored_energy(&mechanical, &next.state)
            .unwrap()
            .kinematic_hardening_j_m3,
        0.0,
        0.0,
    );
    assert!(next.state.alpha > old.state.alpha && next.state.plastic_strain != [0.0; 6]);
    assert_eq!(
        next.state,
        endpoint.law.update_state(&new_strain, &old.state)
    );
    source
        .validate_section_history(h(320.0), &new_strain, &next.state, 1e-3)
        .unwrap();
    assert_eq!(old.state, saved);
}

#[test]
fn kinematic_sources_require_one_pinned_affine_role_and_reject_corrupted_history() {
    let missing = Fixture::new(properties(0.0, false));
    assert!(
        missing
            .source(IDENTITY)
            .with_linear_kinematic_hardening(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
    let negative = Fixture::new(mixed_properties(0.0, false, -3e7));
    assert!(
        negative
            .source(IDENTITY)
            .with_linear_kinematic_hardening(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
    let mut data = mixed_properties(0.0, false, 0.0);
    let claim = data
        .iter_mut()
        .find(|p| p.key.name() == HILL_KINEMATIC_MODULUS_PROPERTY)
        .unwrap();
    claim.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(280.0, 8e8), (300.0, C), (320.0, 1e8)],
        dims: Pressure::DIMS,
    };
    claim.interpolation = InterpolationPolicy::LinearInside;
    let kinked = Fixture::new(data);
    assert!(
        kinked
            .source(IDENTITY)
            .with_linear_kinematic_hardening(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );

    let fixture = Fixture::new(mixed_properties(0.0, false, 0.0));
    let source = mixed_source(&fixture, IDENTITY);
    assert!(
        source
            .clone()
            .with_linear_kinematic_hardening(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
    let endpoint = source.resolve_section_step(0.0, 0.0, 0.1).unwrap();
    let selected = endpoint
        .law
        .resolved()
        .material()
        .property(HILL_KINEMATIC_MODULUS_PROPERTY)
        .unwrap()
        .answer()
        .receipt
        .selected;
    let mut data = mixed_properties(0.0, false, 0.0);
    data.push(property(
        HILL_KINEMATIC_MODULUS_PROPERTY,
        Pressure::DIMS,
        0.5 * C,
        0.0,
    ));
    let ambiguous = Fixture::new(data);
    let candidate = ambiguous.source(IDENTITY);
    assert!(
        candidate
            .clone()
            .with_linear_kinematic_hardening(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
    let pinned = candidate
        .with_linear_kinematic_hardening(MaterialPropertySelection::PinnedByProperty(vec![(
            HILL_KINEMATIC_MODULUS_PROPERTY.into(),
            selected,
        )]))
        .unwrap();
    for temperature in [280.0, 300.0, 320.0] {
        let endpoint = pinned
            .resolve_section_step(0.0, h(temperature), 0.1)
            .unwrap();
        let property = endpoint
            .law
            .resolved()
            .material()
            .property(HILL_KINEMATIC_MODULUS_PROPERTY)
            .unwrap();
        assert_eq!(property.answer().receipt.selected, selected);
        near(
            endpoint.law.resolved().local_law().kinematic_modulus_pa(),
            C,
            0.0,
        );
        ambiguous
            .card
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    let strain = [0.0, 0.0, 0.0, 0.006, 0.0, 0.0];
    let loaded = pinned
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
    let saved = loaded.state.clone();
    let mut invalid = loaded.state.clone();
    invalid.plastic_strain[0] += 0.001;
    assert!(
        pinned
            .validate_section_history(0.0, &strain, &invalid, 1e-3)
            .is_err()
    );
    assert!(
        pinned
            .section_increment(0.0, h(310.0), &strain, &invalid, &strain, tolerance(), 0.1)
            .is_err()
    );
    assert_eq!(loaded.state, saved);
}
