//! Native coefficient curvature and source-qualified finite thermal derivatives.

use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::hyper::cards::ResolvedHyperelastic;
use fs_material::hyper::thermoelastic::source::ResolvedNeoHookeanThermalCoefficients;
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

const T: f64 = 300.0;
const YOUNGS: f64 = 3000.0;
const POISSON: f64 = 0.25;
const SLOPES: [f64; 2] = [-12.0, 0.004];

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            t,
        )
        .unwrap()
}
fn card(slopes: [f64; 2], kink: bool, model_scale: f64) -> MaterialCard {
    card_values([YOUNGS, POISSON], slopes, kink, model_scale, None)
}
fn card_values(
    values: [f64; 2],
    slopes: [f64; 2],
    kink: bool,
    model_scale: f64,
    lambda_override: Option<f64>,
) -> MaterialCard {
    let artifact = fs_blake3::hash_bytes(b"authored finite thermal coefficient fixture");
    let provenance = Provenance {
        source: "authored Neo-Hookean E/nu paths; not measurements".into(),
        license: "test-only".into(),
        artifact: Some(artifact),
    };
    let validity = ValidityDomain::unconstrained().with_quantity(
        "T",
        *point(T).axis_quantities().get("T").unwrap(),
        280.0,
        320.0,
    );
    let mut claims = ClaimSet::new();
    for (index, (name, dims, value)) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, values[0]),
        (POISSON_RATIO_PROPERTY, Dims::NONE, values[1]),
    ]
    .into_iter()
    .enumerate()
    {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (290.0, value - 10.0 * slopes[index]),
                        (T, value),
                        (
                            310.0,
                            value
                                + 10.0 * slopes[index] * if kink && index == 0 { 2.0 } else { 1.0 },
                        ),
                    ],
                    dims,
                },
                validity: validity.clone(),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let mu = values[0] / (2.0 * (1.0 + values[1]));
    let lambda = lambda_override.unwrap_or(2.0 * mu * values[1] / (1.0 - 2.0 * values[1]));
    let model = ConstitutiveModelCard {
        law: LawId("neo-hookean-log".into()),
        law_version: 1,
        parameters: [
            ("mu", mu * model_scale, Pressure::DIMS),
            ("lambda", lambda * model_scale, Pressure::DIMS),
            ("stretch_limit", 4.0, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity,
        sources: vec![artifact],
        provenance,
    };
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored thermal Neo-Hookean".into(),
            phase: "solid".into(),
            process: "analytic fixture".into(),
            revision: 0,
        },
        claims,
        vec![model],
    )
    .unwrap()
}
fn bind(
    card: &MaterialCard,
    range: [f64; 2],
) -> Result<
    ResolvedNeoHookeanThermalCoefficients,
    fs_material::hyper::thermoelastic::source::NeoHookeanThermalSourceError,
> {
    let at = point(T);
    let model = ResolvedHyperelastic::try_from_model_card_at(&card.models()[0], &at).unwrap();
    let coefficients = resolve_material_state_point(
        card,
        &at,
        &[
            ScalarPropertyRequirement::try_new(
                YOUNG_MODULUS_PROPERTY,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                POISSON_RATIO_PROPERTY,
                Dims::NONE,
                ScalarAdmissibility::OpenInterval {
                    lower: -1.0,
                    upper: 0.5,
                },
            )
            .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ResolvedNeoHookeanThermalCoefficients::try_new(card, &model, &coefficients, range)
}

// Independent invariant energy and analytic Piola, with explicit G inverses.
// Both the natural map and the original E/nu paths change with temperature.
fn independent(delta: f64, f: [f64; 9]) -> (f64, [f64; 9]) {
    let g = core::array::from_fn(|i| G[i] + delta * G1[i] + 0.5 * delta * delta * G2[i]);
    let jg = determinant(g);
    let cofactor = cof(g);
    let inverse = core::array::from_fn(|ij| cofactor[3 * (ij % 3) + ij / 3] / jg);
    let fe = multiply(f, inverse);
    let je = determinant(fe);
    let fe_cofactor = cof(fe);
    let e = YOUNGS + SLOPES[0] * delta;
    let nu = POISSON + SLOPES[1] * delta;
    let mu = e / (2.0 * (1.0 + nu));
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let log_j = je.ln();
    let energy = jg
        * (0.5 * mu * (fe.iter().map(|x| x * x).sum::<f64>() - 3.0) - mu * log_j
            + 0.5 * lambda * log_j * log_j);
    let pe: [f64; 9] =
        core::array::from_fn(|i| mu * fe[i] + (lambda * log_j - mu) * fe_cofactor[i] / je);
    let p = core::array::from_fn(|ij| {
        jg * (0..3)
            .map(|k| pe[3 * (ij / 3) + k] * inverse[3 * (ij % 3) + k])
            .sum::<f64>()
    });
    (energy, p)
}

#[test]
fn original_coefficients_and_noncommuting_expansion_share_full_thermal_chain_rule() {
    let card = card(SLOPES, false, 1.0);
    let source = bind(&card, [290.0, 310.0]).unwrap();
    let state = source.thermoelastic_state(&F, jet()).unwrap();
    let h = 0.01;
    let samples = [-2.0 * h, -h, h, 2.0 * h].map(|d| independent(d, F));
    let center = independent(0.0, F);
    let psi_t = derivative(samples.map(|s| s.0), h);
    let psi_tt = (-samples[3].0 + 16.0 * samples[2].0 - 30.0 * center.0 + 16.0 * samples[1].0
        - samples[0].0)
        / (12.0 * h * h);
    near(state.free_energy_j_m3(), center.0, 1e-10, 1e-12);
    near(state.entropy_j_m3_k(), -psi_t, 2e-8, 2e-8);
    near(
        state.internal_energy_j_m3(),
        center.0 - T * psi_t,
        1e-5,
        2e-8,
    );
    near(
        state.mechanical_heat_capacity_j_m3_k(),
        -T * psi_tt,
        2e-6,
        2e-7,
    );
    for component in 0..9 {
        let beta = derivative(samples.map(|s| s.1[component]), h);
        near(
            state.first_piola_pa()[component],
            center.1[component],
            2e-10,
            1e-12,
        );
        near(
            state.stress_temperature_derivative_pa_k()[component],
            beta,
            3e-8,
            2e-8,
        );
        let mut plus = F;
        plus[component] += 1e-6;
        let mut minus = F;
        minus[component] -= 1e-6;
        let entropy_f = (source
            .thermoelastic_state(&plus, jet())
            .unwrap()
            .entropy_j_m3_k()
            - source
                .thermoelastic_state(&minus, jet())
                .unwrap()
                .entropy_j_m3_k())
            / 2e-6;
        near(-entropy_f, beta, 2e-6, 2e-7);
    }
    for (name, line) in [YOUNG_MODULUS_PROPERTY, POISSON_RATIO_PROPERTY]
        .into_iter()
        .zip(source.source_lines())
    {
        let original = source
            .material()
            .property(name)
            .unwrap()
            .answer()
            .receipt
            .selected;
        assert_eq!(line.knots.len(), 3);
        assert!(
            line.knots
                .iter()
                .all(|(_, answer)| answer.receipt.selected == original)
        );
    }
    assert_eq!(source.material().card_identity(), card.content_hash());
}

#[test]
fn constant_limit_is_exact_and_variable_moduli_drive_fixed_map_capacity_and_waves() {
    let constant = bind(&card([0.0; 2], false, 1.0), [290.0, 310.0]).unwrap();
    assert_eq!(
        constant.thermoelastic_state(&F, jet()).unwrap(),
        FiniteThermoelasticState::try_new(constant.model().law(), &F, T, jet()).unwrap()
    );
    let source = bind(&card(SLOPES, false, 1.0), [290.0, 310.0]).unwrap();
    let held =
        ThermalReferenceJet::try_new(ReferenceDeformation::identity(), [0.0; 9], [0.0; 9]).unwrap();
    let state = source.thermoelastic_state(&F, held).unwrap();
    let fixed = FiniteThermoelasticState::try_new(source.model().law(), &F, T, held).unwrap();
    assert_eq!(state.first_piola_pa(), fixed.first_piola_pa());
    assert_eq!(state.isothermal_tangent_pa(), fixed.isothermal_tangent_pa());
    assert!(state.mechanical_heat_capacity_j_m3_k().abs() > 0.1);
    assert_eq!(fixed.mechanical_heat_capacity_j_m3_k(), 0.0);
    assert!(
        source
            .second_derivatives_pa_k2()
            .iter()
            .any(|d| d.abs() > 0.01)
    );
    let actual = state.adiabatic_acoustics(6.0, 1000.0).unwrap();
    let baseline = fixed.adiabatic_acoustics(6.0, 1000.0).unwrap();
    near(
        actual.fixed_deformation_heat_capacity_j_m3_k(),
        1000.0 + state.mechanical_heat_capacity_j_m3_k(),
        1e-12,
        1e-12,
    );
    let waves = actual
        .acoustic()
        .directional_waves([0.3, -0.4, 0.5])
        .unwrap();
    let old = baseline
        .acoustic()
        .directional_waves(waves.direction)
        .unwrap();
    assert!(
        (waves.modes[2].squared_phase_speed_m2_s2 - old.modes[2].squared_phase_speed_m2_s2).abs()
            > 1e-3
    );
    assert!(state.adiabatic_acoustics(6.0, 1e-8).is_err());
}

#[test]
fn coefficient_sources_refuse_kinks_wrong_model_context_and_unsupported_intervals() {
    // Finite E/nu can imply an infinite lambda near the incompressible limit.
    // A relative comparison must not accept inf <= inf against a finite model.
    let overflow = card_values(
        [1e300, 0.499_999_999_999_9],
        [0.0; 2],
        false,
        1.0,
        Some(1e299),
    );
    assert!(bind(&overflow, [290.0, 310.0]).is_err());
    for (slopes, kink, model_scale) in [
        (SLOPES, true, 1.0),
        (SLOPES, false, 1.01),
        ([1000.0, 0.004], false, 1.0),
        ([-12.0, 0.04], false, 1.0),
    ] {
        assert!(bind(&card(slopes, kink, model_scale), [290.0, 310.0]).is_err());
    }
    let original = card(SLOPES, false, 1.0);
    for range in [
        [300.0, 310.0],
        [290.0, 300.0],
        [280.0, 320.0],
        [0.0, 310.0],
        [290.0, f64::INFINITY],
    ] {
        assert!(bind(&original, range).is_err());
    }
    let source = bind(&original, [290.0, 310.0]).unwrap();
    let foreign = card([0.0; 2], false, 1.0);
    assert!(
        ResolvedNeoHookeanThermalCoefficients::try_new(
            &foreign,
            source.model(),
            source.material(),
            [290.0, 310.0]
        )
        .is_err()
    );
    let other_point = point(301.0);
    let other_model =
        ResolvedHyperelastic::try_from_model_card_at(&original.models()[0], &other_point).unwrap();
    assert!(
        ResolvedNeoHookeanThermalCoefficients::try_new(
            &original,
            &other_model,
            source.material(),
            [290.0, 310.0]
        )
        .is_err()
    );
    assert_eq!(source.material().query(), point(T));
}

#[test]
fn source_owned_rates_include_coefficient_changes_and_relaxed_thermal_motion() {
    use fs_material::hyper::mixed::VolumetricTreatment;
    use fs_material::hyper::plane_stress::{PlaneStressSettings, plane_stress_with_seed};
    use fs_material::hyper::thermoelastic::plane_stress::FiniteThermoelasticPlaneStressState;

    let source = bind(&card(SLOPES, false, 1.0), [290.0, 310.0]).unwrap();
    let thermal_at = |delta: f64| {
        ThermalReferenceJet::try_new(
            reference_at(delta),
            core::array::from_fn(|i| G1[i] + delta * G2[i]),
            G2,
        )
        .unwrap()
    };
    let tdot = 2.3;
    let fdot = [0.03, -0.04, 0.02, -0.01, 0.05, 0.03, -0.02, 0.01, -0.04];
    let h = 0.01;
    let state = source.thermoelastic_state(&F, jet()).unwrap();
    let rate = state.rate(fdot, tdot).unwrap();
    let samples = [-2.0 * h, -h, h, 2.0 * h]
        .map(|time| independent(tdot * time, core::array::from_fn(|i| F[i] + time * fdot[i])));
    near(
        rate.free_energy_rate_w_m3(),
        derivative(samples.map(|s| s.0), h),
        2e-7,
        2e-8,
    );
    for i in 0..9 {
        near(
            rate.first_piola_rate_pa_s()[i],
            derivative(samples.map(|s| s.1[i]), h),
            2e-7,
            2e-8,
        );
    }
    let settings = PlaneStressSettings {
        relative_tolerance: 2e-12,
        ..Default::default()
    };
    let a = [1.22, 0.16, -0.05, 0.91];
    let adot = [0.08, -0.03, 0.01, -0.05];
    let relaxed_at = |time: f64| {
        let delta = tdot * time;
        let a: [f64; 4] = core::array::from_fn(|i| a[i] + time * adot[i]);
        let thermal = thermal_at(delta);
        let reference = *thermal.reference();
        // This independent native relaxation uses the authored E/nu fixture
        // at the sample temperature. The thermal derivative owner below
        // keeps the original card and source interval.
        let young = YOUNGS + SLOPES[0] * delta;
        let poisson = POISSON + SLOPES[1] * delta;
        let mu = young / (2.0 * (1.0 + poisson));
        let lambda = 2.0 * mu * poisson / (1.0 - 2.0 * poisson);
        let law = Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 4.0).unwrap();
        let plane = plane_stress_with_seed(
            &ReferencedHyperelastic::new(&law, reference),
            &a,
            [
                0.0,
                0.0,
                reference.volume_ratio() / (a[0] * a[3] - a[1] * a[2]),
            ],
            settings,
        )
        .unwrap();
        let full = source
            .thermoelastic_state_at_temperature(T + delta, &plane.deformation_gradient, thermal)
            .unwrap();
        FiniteThermoelasticPlaneStressState::try_from_state(
            full,
            VolumetricTreatment::OriginalMaterial,
            settings,
        )
        .unwrap()
    };
    let plane = relaxed_at(0.0);
    let rate = plane.rate(adot, tdot).unwrap();
    let samples = [-2.0 * h, -h, h, 2.0 * h].map(relaxed_at);
    for i in 0..9 {
        near(
            rate.deformation_gradient_rate_per_s()[i],
            derivative(
                samples
                    .each_ref()
                    .map(|s| s.response().deformation_gradient[i]),
                h,
            ),
            5e-8,
            3e-7,
        );
        near(
            rate.first_piola_rate_pa_s()[i],
            derivative(
                samples
                    .each_ref()
                    .map(|s| s.full_state().first_piola_pa()[i]),
                h,
            ),
            3e-6,
            3e-7,
        );
    }
    for (actual, values) in [
        (
            rate.free_energy_rate_w_m3(),
            samples.each_ref().map(|s| s.free_energy_j_m3()),
        ),
        (
            rate.entropy_rate_w_m3_k(),
            samples.each_ref().map(|s| s.entropy_j_m3_k()),
        ),
        (
            rate.internal_energy_rate_w_m3(),
            samples.each_ref().map(|s| s.internal_energy_j_m3()),
        ),
    ] {
        near(actual, derivative(values, h), 3e-6, 3e-7);
    }
    for i in [2, 5, 8] {
        near(rate.first_piola_rate_pa_s()[i], 0.0, 2e-9, 0.0);
    }
}

// Closed constrained NH potential with independently specified exponential
// natural stretch and rational mu(T)=E(T)/(2*(1+nu(T))).
struct ExactThermalOracle {
    energy: f64,
    entropy: f64,
    internal_energy: f64,
    capacity: f64,
    pressure: f64,
    pressure_t: f64,
    piola: [f64; 4],
    beta: [f64; 4],
    deformation: [f64; 9],
    elastic_deformation: [f64; 9],
}
fn exact_thermal_jet(delta: f64) -> ThermalReferenceJet {
    let alpha = 0.002;
    let stretch = 1.1 * (alpha * delta).exp();
    ThermalReferenceJet::try_new(
        ReferenceDeformation::try_new(diagonal(stretch)).unwrap(),
        diagonal(alpha * stretch),
        diagonal(alpha * alpha * stretch),
    )
    .unwrap()
}
fn exact_thermal_oracle(delta: f64, a: [f64; 4]) -> ExactThermalOracle {
    let alpha = 0.002;
    let s = 1.1 * (alpha * delta).exp();
    let e = YOUNGS + SLOPES[0] * delta;
    let n = 1.0 + POISSON + SLOPES[1] * delta;
    let mu = e / (2.0 * n);
    let mu_t = SLOPES[0] / (2.0 * n) - e * SLOPES[1] / (2.0 * n * n);
    let mu_tt = -SLOPES[0] * SLOPES[1] / (n * n) + e * SLOPES[1] * SLOPES[1] / n.powi(3);
    let d = a[0] * a[3] - a[1] * a[2];
    let da = [a[3], -a[2], -a[1], a[0]];
    let ia: f64 = a.iter().map(|v| v * v).sum();
    let h = 0.5 * (s * ia + s.powi(7) / (d * d) - 3.0 * s.powi(3));
    let h_t = 0.5 * alpha * (s * ia + 7.0 * s.powi(7) / (d * d) - 9.0 * s.powi(3));
    let h_tt = 0.5 * alpha * alpha * (s * ia + 49.0 * s.powi(7) / (d * d) - 27.0 * s.powi(3));
    let energy = mu * h;
    let entropy = -mu_t * h - mu * h_t;
    let pressure_shape = s.powi(4) / (d * d) - 1.0;
    let deformation = [a[0], a[1], 0.0, a[2], a[3], 0.0, 0.0, 0.0, s.powi(3) / d];
    ExactThermalOracle {
        energy,
        entropy,
        internal_energy: energy + (T + delta) * entropy,
        capacity: -(T + delta) * (mu_tt * h + 2.0 * mu_t * h_t + mu * h_tt),
        pressure: mu * pressure_shape,
        pressure_t: mu_t * pressure_shape + 4.0 * mu * alpha * s.powi(4) / (d * d),
        piola: core::array::from_fn(|i| mu * (s * a[i] - s.powi(7) * da[i] / d.powi(3))),
        beta: core::array::from_fn(|i| {
            mu_t * (s * a[i] - s.powi(7) * da[i] / d.powi(3))
                + mu * alpha * (s * a[i] - 7.0 * s.powi(7) * da[i] / d.powi(3))
        }),
        deformation,
        elastic_deformation: deformation.map(|v| v / s),
    }
}

#[test]
fn exact_volume_sources_match_independent_modulus_and_natural_volume_product_rule() {
    use fs_material::hyper::plane_stress::PlaneStressSettings;
    let source = bind(&card(SLOPES, false, 1.0), [290.0, 310.0]).unwrap();
    let a = [1.22, 0.16, -0.05, 0.91];
    for delta in [0.0, 4.0] {
        let state = source
            .incompressible_plane_stress_state_at_temperature(
                T + delta,
                &a,
                exact_thermal_jet(delta),
                None,
                PlaneStressSettings::default(),
            )
            .unwrap();
        let expected = exact_thermal_oracle(delta, a);
        for (actual, expected) in [
            (state.free_energy_j_m3(), expected.energy),
            (state.entropy_j_m3_k(), expected.entropy),
            (state.internal_energy_j_m3(), expected.internal_energy),
            (state.mechanical_heat_capacity_j_m3_k(), expected.capacity),
            (state.response().pressure_pa, expected.pressure),
            (
                state.pressure_temperature_derivative_pa_k(),
                expected.pressure_t,
            ),
        ] {
            near(actual, expected, 3e-9, 3e-11);
        }
        for i in 0..4 {
            near(state.response().piola[i], expected.piola[i], 3e-9, 3e-11);
            near(
                state.stress_temperature_derivative_pa_k()[i],
                expected.beta[i],
                3e-10,
                3e-11,
            );
        }
    }
    for temperature in [289.9, 310.1, f64::NAN, f64::INFINITY] {
        assert!(
            source
                .incompressible_plane_stress_state_at_temperature(
                    temperature,
                    &a,
                    exact_thermal_jet(0.0),
                    None,
                    PlaneStressSettings::default(),
                )
                .is_err()
        );
    }
}

#[test]
fn exact_volume_source_rates_follow_the_same_closed_thermal_trajectory() {
    use fs_material::hyper::plane_stress::PlaneStressSettings;
    let source = bind(&card(SLOPES, false, 1.0), [290.0, 310.0]).unwrap();
    let a = [1.22, 0.16, -0.05, 0.91];
    let adot = [0.08, -0.03, 0.01, -0.05];
    let tdot = 2.3;
    let state = source
        .incompressible_plane_stress_state(
            &a,
            exact_thermal_jet(0.0),
            None,
            PlaneStressSettings::default(),
        )
        .unwrap();
    let rate = state.rate(adot, tdot).unwrap();
    let h = 0.01;
    let samples = [-2.0 * h, -h, h, 2.0 * h].map(|time| {
        exact_thermal_oracle(tdot * time, core::array::from_fn(|i| a[i] + time * adot[i]))
    });
    for i in 0..9 {
        near(
            rate.material().deformation_gradient_rate_per_s()[i],
            derivative(samples.each_ref().map(|s| s.deformation[i]), h),
            3e-10,
            3e-9,
        );
        near(
            rate.material().elastic_deformation_gradient_rate_per_s()[i],
            derivative(samples.each_ref().map(|s| s.elastic_deformation[i]), h),
            3e-10,
            3e-9,
        );
    }
    for (i, component) in [0, 1, 3, 4].into_iter().enumerate() {
        near(
            rate.material().first_piola_rate_pa_s()[component],
            derivative(samples.each_ref().map(|s| s.piola[i]), h),
            3e-7,
            3e-8,
        );
    }
    for (actual, samples) in [
        (
            rate.pressure_rate_pa_s(),
            samples.each_ref().map(|s| s.pressure),
        ),
        (
            rate.material().free_energy_rate_w_m3(),
            samples.each_ref().map(|s| s.energy),
        ),
        (
            rate.material().entropy_rate_w_m3_k(),
            samples.each_ref().map(|s| s.entropy),
        ),
        (
            rate.material().internal_energy_rate_w_m3(),
            samples.each_ref().map(|s| s.internal_energy),
        ),
    ] {
        near(actual, derivative(samples, h), 3e-7, 3e-8);
    }
    for component in [2, 5, 8] {
        near(
            rate.material().first_piola_rate_pa_s()[component],
            0.0,
            3e-10,
            0.0,
        );
    }
}

#[test]
fn exact_source_constant_limit_and_changing_bulk_modulus_keep_native_remainder_scope() {
    use fs_material::hyper::plane_stress::PlaneStressSettings;
    use fs_material::hyper::thermoelastic::incompressible_plane_stress::IncompressibleThermoelasticPlaneStressState;
    let a = [1.22, 0.16, -0.05, 0.91];
    let settings = PlaneStressSettings::default();
    let constant = bind(&card([0.0; 2], false, 1.0), [290.0, 310.0]).unwrap();
    let sourced = constant
        .incompressible_plane_stress_state(&a, jet(), None, settings)
        .unwrap();
    let native = IncompressibleThermoelasticPlaneStressState::try_new(
        constant.model().law(),
        &a,
        T,
        jet(),
        None,
        settings,
    )
    .unwrap();
    assert_eq!(sourced.response().full_piola, native.response().full_piola);
    assert_eq!(sourced.response().tangent, native.response().tangent);
    assert_eq!(sourced.entropy_j_m3_k(), native.entropy_j_m3_k());
    assert_eq!(
        sourced.mechanical_heat_capacity_j_m3_k(),
        native.mechanical_heat_capacity_j_m3_k()
    );
    assert_eq!(
        sourced.rate([0.08, -0.03, 0.01, -0.05], 2.3).unwrap(),
        native.rate([0.08, -0.03, 0.01, -0.05], 2.3).unwrap(),
    );
    // E=2*mu*(1+nu) changes while mu is constant. The original lambda
    // source has a nonzero derivative, which must not enter this exact mode.
    let variable = bind(&card([9.6, 0.004], false, 1.0), [290.0, 310.0]).unwrap();
    assert!(variable.first_derivatives_pa_k()[1].abs() > 1.0);
    let held =
        ThermalReferenceJet::try_new(ReferenceDeformation::identity(), [0.0; 9], [0.0; 9]).unwrap();
    let state = variable
        .incompressible_plane_stress_state(&a, held, None, settings)
        .unwrap();
    near(state.entropy_j_m3_k(), 0.0, 3e-11, 0.0);
    near(state.mechanical_heat_capacity_j_m3_k(), 0.0, 3e-10, 0.0);
    for beta in state.stress_temperature_derivative_pa_k() {
        near(*beta, 0.0, 3e-11, 0.0);
    }
    let rate = state.rate([0.0; 4], 2.3).unwrap();
    near(rate.pressure_rate_pa_s(), 0.0, 3e-10, 0.0);
    for value in rate.material().first_piola_rate_pa_s() {
        near(*value, 0.0, 3e-10, 0.0);
    }
    for value in [
        rate.material().free_energy_rate_w_m3(),
        rate.material().entropy_rate_w_m3_k(),
        rate.material().internal_energy_rate_w_m3(),
    ] {
        near(value, 0.0, 3e-10, 0.0);
    }
}

#[test]
fn exact_source_finite_increments_keep_changing_moduli_volume_and_constrained_entropy() {
    use fs_material::hyper::plane_stress::PlaneStressSettings;
    let source = bind(&card(SLOPES, false, 1.0), [290.0, 310.0]).unwrap();
    let a = [1.22, 0.16, -0.05, 0.91];
    let before = source
        .incompressible_plane_stress_state(
            &a,
            exact_thermal_jet(0.0),
            None,
            PlaneStressSettings::default(),
        )
        .unwrap();
    let after = source
        .incompressible_plane_stress_state_at_temperature(
            T + 4.0,
            &a,
            exact_thermal_jet(4.0),
            None,
            PlaneStressSettings::default(),
        )
        .unwrap();
    let first = exact_thermal_oracle(0.0, a);
    let second = exact_thermal_oracle(4.0, a);
    near(
        before.free_energy_increment_to(&after).unwrap(),
        second.energy - first.energy,
        3e-9,
        3e-11,
    );
    near(
        before.internal_energy_increment_to(&after).unwrap(),
        second.internal_energy - first.internal_energy,
        3e-9,
        3e-11,
    );
    near(
        after.internal_energy_increment_to(&before).unwrap(),
        first.internal_energy - second.internal_energy,
        3e-9,
        3e-11,
    );
    assert_eq!(before.free_energy_increment_to(&before).unwrap(), 0.0);
    assert_eq!(before.internal_energy_increment_to(&before).unwrap(), 0.0);
}
