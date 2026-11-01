//! Synthetic joint creep/yield references, not empirical material calibrations.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::plastic::finite::thermal::source::{
    ResolvedAffineThermalFiniteJ2, THERMAL_FINITE_J2_LAW_ID, THERMAL_FINITE_NORTON_J2_LAW_ID,
    THERMAL_FINITE_NORTON_LAW_ID, THERMAL_PERZYNA_FINITE_J2_LAW_ID,
};
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalError,
    FiniteJ2ThermalOverstress, FiniteJ2ThermalParameters, FiniteNortonThermalParameters,
};
use fs_material::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

const ID: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
const LIMITS: FiniteJ2Limits = FiniteJ2Limits {
    minimum_jacobian: 0.1,
    maximum_elastic_log_strain: 0.8,
    maximum_plastic_log_strain: 3.0,
    maximum_accumulated_plastic_strain: 3.0,
};
fn parameters() -> FiniteJ2ThermalParameters {
    FiniteJ2ThermalParameters {
        reference_temperature_k: 300.0,
        temperature_range_k: [280.0, 500.0],
        reference: FiniteJ2Parameters {
            shear_modulus_pa: 80e6,
            bulk_modulus_pa: 160e6,
            yield_stress_pa: 2e6,
            hardening_modulus_pa: 8e6,
        },
        slopes: FiniteJ2TemperatureSlopes {
            shear_modulus_pa_k: -20_000.0,
            bulk_modulus_pa_k: -30_000.0,
            yield_stress_pa_k: -4_000.0,
            hardening_modulus_pa_k: -10_000.0,
        },
        logarithmic_expansion_per_k: 2e-4,
        logarithmic_expansion_slope_per_k2: 1e-7,
    }
}
fn kinetics(exponent: f64) -> FiniteJ2ThermalOverstress {
    FiniteJ2ThermalOverstress {
        reference_viscosity_pa_s: 60e6,
        stress_scale_pa: 2e6,
        exponent,
        activation_temperature_k: 1800.0,
    }
}
fn family(exponent: f64) -> AffineThermalFiniteJ2 {
    AffineThermalFiniteJ2::new(parameters(), LIMITS)
        .unwrap()
        .with_norton_creep(kinetics(exponent), 3.0)
        .unwrap()
}
fn coaxial(mean: f64, axial: f64) -> [f64; 9] {
    [
        (mean + axial).exp(),
        0.0,
        0.0,
        0.0,
        (mean - 0.5 * axial).exp(),
        0.0,
        0.0,
        0.0,
        (mean - 0.5 * axial).exp(),
    ]
}
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

#[test]
fn activated_linear_creep_and_yield_match_an_independent_joint_return_and_thermal_volume() {
    let p = parameters();
    let k = kinetics(1.0);
    let family = family(1.0);
    let (mean, axial, dt) = (0.03, 0.15, 0.25);
    let mut previous_creep = 0.0;
    for temperature in [300.0, 420.0] {
        let delta = temperature - p.reference_temperature_k;
        let mu = p.reference.shear_modulus_pa + p.slopes.shear_modulus_pa_k * delta;
        let bulk = p.reference.bulk_modulus_pa + p.slopes.bulk_modulus_pa_k * delta;
        let yield_stress = p.reference.yield_stress_pa + p.slopes.yield_stress_pa_k * delta;
        let hardening = p.reference.hardening_modulus_pa + p.slopes.hardening_modulus_pa_k * delta;
        let theta = (p.logarithmic_expansion_per_k * delta
            + 0.5 * p.logarithmic_expansion_slope_per_k2 * delta * delta)
            .exp();
        let eta = k.reference_viscosity_pa_s
            * (k.activation_temperature_k * (1.0 / temperature - 1.0 / 300.0)).exp();
        let f = coaxial(mean, axial).map(|x| theta * x);

        // Solve the active linear complementarity equations directly. Both
        // mechanisms share q, while only plastic accumulation hardens yield.
        let a = 3.0 * mu;
        let r = eta / dt;
        let dp = (a * axial - yield_stress * (1.0 + a / r)) / (a + hardening * (1.0 + a / r));
        let q = yield_stress + hardening * dp;
        let dc = q / r;
        assert!(dp > 0.0 && dc > previous_creep);
        previous_creep = dc;
        let elastic_axial = axial - dp - dc;
        let stored = 1.5 * mu * elastic_axial.powi(2)
            + 0.5 * bulk * (3.0 * mean).powi(2)
            + 0.5 * hardening * dp.powi(2);
        let initial = family.initial_state();
        let response = family
            .at_temperature_for_step(temperature, dt)
            .unwrap()
            .response(&f, &initial)
            .unwrap();
        near(response.plastic_increment, dp, 3e-14, 3e-11);
        near(response.creep_increment, dc, 3e-14, 3e-11);
        near(response.state.alpha, dp, 3e-14, 3e-11);
        near(response.state.creep_alpha, dc, 3e-14, 3e-11);
        near(
            response.state.inverse_plastic_cauchy_green[0],
            (-2.0 * (dp + dc)).exp(),
            3e-14,
            3e-11,
        );
        for i in [1, 2] {
            near(
                response.state.inverse_plastic_cauchy_green[i],
                (dp + dc).exp(),
                3e-14,
                3e-11,
            );
        }
        near(response.equivalent_kirchhoff_pa, q, 3e-6, 3e-11);
        near(
            response.jacobian,
            theta.powi(3) * (3.0 * mean).exp(),
            3e-14,
            3e-12,
        );
        near(response.stored_energy_j_m3, stored, 3e-6, 3e-12);
        near(
            response.plastic_dissipation_j_m3,
            yield_stress * dp,
            3e-6,
            3e-11,
        );
        near(response.creep_dissipation_j_m3, q * dc, 3e-6, 3e-11);
        near(
            response.inelastic_dissipation_j_m3(),
            yield_stress * dp + q * dc,
            3e-6,
            3e-11,
        );
        near(
            response.incremental_potential_j_m3,
            stored + yield_stress * dp + 0.5 * q * dc,
            3e-6,
            3e-12,
        );
        for (i, tau) in [
            (0, 3.0 * bulk * mean + 2.0 * mu * elastic_axial),
            (4, 3.0 * bulk * mean - mu * elastic_axial),
            (8, 3.0 * bulk * mean - mu * elastic_axial),
        ] {
            near(response.first_piola_pa[i], tau / f[i], 3e-6, 3e-12);
        }
        family
            .validate_committed(&f, &response.state, temperature, 1e-3)
            .unwrap();
        assert_eq!(initial, family.initial_state());
    }
}

#[test]
fn held_deformation_creeps_without_further_yield_or_creep_hardening_storage() {
    let mut p = parameters();
    p.slopes = FiniteJ2TemperatureSlopes {
        shear_modulus_pa_k: 0.0,
        bulk_modulus_pa_k: 0.0,
        yield_stress_pa_k: 0.0,
        hardening_modulus_pa_k: 0.0,
    };
    p.logarithmic_expansion_per_k = 0.0;
    p.logarithmic_expansion_slope_per_k2 = 0.0;
    let constant = AffineThermalFiniteJ2::new(p, LIMITS)
        .unwrap()
        .with_norton_creep(kinetics(1.0), 3.0)
        .unwrap();
    let f = coaxial(0.0, 0.15);
    let load = constant
        .at_temperature_for_step(300.0, 0.05)
        .unwrap()
        .response(&f, &constant.initial_state())
        .unwrap();
    assert!(load.plastic_increment > 0.0 && load.creep_increment > 0.0);
    let old = load.state;
    let hold = constant
        .checked_increment_for_step(&f, &old, 300.0, &f, 300.0, 1e-3, 0.4)
        .unwrap();
    let q_trial = p.reference.yield_stress_pa + p.reference.hardening_modulus_pa * old.alpha;
    let r = kinetics(1.0).reference_viscosity_pa_s / 0.4;
    let dc = q_trial / (3.0 * p.reference.shear_modulus_pa + r);
    assert_eq!(hold.mechanical.plastic_increment, 0.0);
    assert_eq!(hold.mechanical.plastic_dissipation_j_m3, 0.0);
    assert_eq!(hold.mechanical.state.alpha, old.alpha);
    near(hold.mechanical.creep_increment, dc, 3e-14, 3e-11);
    near(
        hold.mechanical.state.creep_alpha,
        old.creep_alpha + dc,
        3e-14,
        3e-11,
    );
    near(hold.mechanical.equivalent_kirchhoff_pa, r * dc, 3e-6, 3e-11);
    near(hold.caloric_heat_j_m3(), r * dc * dc, 3e-6, 3e-11);
    near(hold.mechanical_work_j_m3, 0.0, 0.0, 0.0);
    near(hold.balance_error_j_m3, 0.0, 3e-6, 0.0);

    let family = family(1.0);
    let fixed = family.thermodynamics(&f, &old, 420.0).unwrap();
    let mut more_creep = old.clone();
    more_creep.creep_alpha += 0.2;
    assert_eq!(
        family.thermodynamics(&f, &more_creep, 420.0).unwrap(),
        fixed
    );
    let mut more_plastic = old.clone();
    more_plastic.alpha += 0.1;
    let changed = family.thermodynamics(&f, &more_plastic, 420.0).unwrap();
    let alpha_squared_change = more_plastic.alpha.powi(2) - old.alpha.powi(2);
    let p = parameters();
    let h = p.reference.hardening_modulus_pa + 120.0 * p.slopes.hardening_modulus_pa_k;
    near(
        changed.free_energy_j_m3 - fixed.free_energy_j_m3,
        0.5 * h * alpha_squared_change,
        3e-6,
        3e-12,
    );
    near(
        changed.entropy_j_m3_k - fixed.entropy_j_m3_k,
        -0.5 * p.slopes.hardening_modulus_pa_k * alpha_squared_change,
        3e-8,
        3e-12,
    );
    assert_eq!(changed.first_piola_pa, fixed.first_piola_pa);
    let mut forgotten_creep = old.clone();
    forgotten_creep.creep_alpha = 0.0;
    assert!(family.thermodynamics(&f, &forgotten_creep, 300.0).is_err());
}

#[test]
fn activated_nonlinear_joint_return_preserves_exact_tangent_and_two_temperature_heat() {
    let family = family(2.4);
    let old_f = coaxial(0.01, 0.06);
    let old = family
        .at_temperature_for_step(300.0, 0.02)
        .unwrap()
        .response(&old_f, &family.initial_state())
        .unwrap()
        .state;
    let f = [1.19, 0.08, 0.02, -0.02, 0.93, 0.04, 0.01, -0.01, 1.04];
    let (temperature, dt) = (410.0, 0.015);
    // The accepted old state belongs to its actual old temperature. Heating
    // softens the yield surface, and that new endpoint excess needs a return.
    family
        .validate_committed(&old_f, &old, 300.0, 1e-3)
        .unwrap();
    assert!(
        family
            .validate_committed(&old_f, &old, temperature, 1e-3)
            .is_err()
    );
    let law = family.at_temperature_for_step(temperature, dt).unwrap();
    let response = law.response(&f, &old).unwrap();
    assert!(response.plastic_increment > 0.0 && response.creep_increment > 0.0);
    let k = kinetics(2.4);
    let eta = k.reference_viscosity_pa_s * (1800.0 * (1.0 / temperature - 1.0 / 300.0)).exp();
    let q = k.stress_scale_pa
        * (eta * response.creep_increment / (dt * k.stress_scale_pa)).powf(1.0 / k.exponent);
    let p = parameters();
    let yield_stress = p.reference.yield_stress_pa + 110.0 * p.slopes.yield_stress_pa_k;
    let hardening = p.reference.hardening_modulus_pa + 110.0 * p.slopes.hardening_modulus_pa_k;
    near(response.equivalent_kirchhoff_pa, q, 3e-6, 3e-11);
    near(
        q,
        yield_stress + hardening * response.state.alpha,
        3e-6,
        3e-11,
    );
    near(
        response.plastic_dissipation_j_m3,
        yield_stress * response.plastic_increment,
        3e-6,
        3e-11,
    );
    near(
        response.creep_dissipation_j_m3,
        q * response.creep_increment,
        3e-6,
        3e-11,
    );
    near(
        response.incremental_potential_j_m3,
        response.stored_energy_j_m3
            + response.plastic_dissipation_j_m3
            + (2.4 / 3.4) * response.creep_dissipation_j_m3,
        3e-6,
        3e-12,
    );
    for col in 0..9 {
        let h = 1e-6;
        let mut plus = f;
        plus[col] += h;
        let mut minus = f;
        minus[col] -= h;
        let a = law.stress_and_potential(&plus, &old).unwrap();
        let b = law.stress_and_potential(&minus, &old).unwrap();
        near(
            (a.1 - b.1) / (2.0 * h),
            response.first_piola_pa[col],
            0.1,
            4e-7,
        );
        for row in 0..9 {
            near(
                (a.0[row] - b.0[row]) / (2.0 * h),
                response.tangent_pa[row][col],
                8.0,
                3e-6,
            );
            near(
                response.tangent_pa[row][col],
                response.tangent_pa[col][row],
                3e-6,
                3e-12,
            );
        }
    }
    let step = family
        .checked_increment_for_step(&old_f, &old, 300.0, &f, temperature, 1e-3, dt)
        .unwrap();
    assert_eq!(step.mechanical, response);
    assert_ne!(step.reversible_heat_j_m3, 0.0);
    assert_eq!(
        step.caloric_heat_j_m3(),
        response.inelastic_dissipation_j_m3() + step.reversible_heat_j_m3
    );
    near(
        step.mechanical_work_j_m3,
        step.internal_energy_change_j_m3()
            + step.caloric_heat_j_m3()
            + step.endpoint_work_defect_j_m3,
        3e-6,
        3e-12,
    );
    let fixed = family.thermodynamics(&f, &old, temperature).unwrap();
    let h = 0.01;
    let plus = family.thermodynamics(&f, &old, temperature + h).unwrap();
    let minus = family.thermodynamics(&f, &old, temperature - h).unwrap();
    near(
        -(plus.free_energy_j_m3 - minus.free_energy_j_m3) / (2.0 * h),
        fixed.entropy_j_m3_k,
        1e-3,
        1e-7,
    );
    near(
        (plus.internal_energy_j_m3 - minus.internal_energy_j_m3) / (2.0 * h),
        fixed.heat_capacity_j_m3_k,
        1e-3,
        1e-7,
    );
    for i in 0..9 {
        near(
            (plus.first_piola_pa[i] - minus.first_piola_pa[i]) / (2.0 * h),
            fixed.first_piola_temperature_derivative_fixed_history_pa_k[i],
            2e-3,
            2e-7,
        );
    }
}

#[test]
fn joint_selection_requires_time_and_keeps_equilibrium_yield_and_separate_history_bounds() {
    let plastic = AffineThermalFiniteJ2::new(parameters(), LIMITS).unwrap();
    let combined = family(1.0);
    assert_eq!(combined.maximum_accumulated_creep_strain(), Some(3.0));
    assert_eq!(combined.overstress(), Some(&kinetics(1.0)));
    assert!(!combined.is_norton_creep());
    assert!(matches!(
        combined.at_temperature(300.0),
        Err(FiniteJ2ThermalError::DurationRequired)
    ));
    let virgin = combined.initial_state();
    let inspector = combined.at_temperature_for_inspection(300.0).unwrap();
    assert!(matches!(
        inspector.response(&ID, &virgin),
        Err(FiniteJ2ThermalError::DurationRequired)
    ));
    assert!(matches!(
        inspector.stress_and_potential(&ID, &virgin),
        Err(FiniteJ2ThermalError::DurationRequired)
    ));
    assert!(matches!(
        combined.checked_increment(&ID, &virgin, 300.0, &ID, 300.0, 0.0),
        Err(FiniteJ2ThermalError::DurationRequired)
    ));
    let overstressed = coaxial(0.0, 0.02);
    assert!(
        inspector
            .validate_committed(&overstressed, &virgin, 0.0)
            .is_err()
    );
    let perzyna = plastic.clone().with_overstress(kinetics(1.0)).unwrap();
    perzyna
        .validate_committed(&overstressed, &virgin, 300.0, 0.0)
        .unwrap();
    let p = parameters();
    let norton = AffineThermalFiniteJ2::new_norton(
        FiniteNortonThermalParameters {
            reference_temperature_k: p.reference_temperature_k,
            temperature_range_k: p.temperature_range_k,
            shear_modulus_pa: p.reference.shear_modulus_pa,
            bulk_modulus_pa: p.reference.bulk_modulus_pa,
            shear_modulus_slope_pa_k: p.slopes.shear_modulus_pa_k,
            bulk_modulus_slope_pa_k: p.slopes.bulk_modulus_pa_k,
            logarithmic_expansion_per_k: p.logarithmic_expansion_per_k,
            logarithmic_expansion_slope_per_k2: p.logarithmic_expansion_slope_per_k2,
        },
        LIMITS,
        kinetics(1.0),
    )
    .unwrap();
    let mut extra_history = virgin.clone();
    extra_history.creep_alpha = 0.01;
    for legacy in [&plastic, &perzyna, &norton] {
        assert!(legacy.thermodynamics(&ID, &extra_history, 300.0).is_err());
    }
    assert!(perzyna.with_norton_creep(kinetics(1.0), 3.0).is_err());
    assert!(norton.with_norton_creep(kinetics(1.0), 3.0).is_err());
    assert!(combined.clone().with_overstress(kinetics(1.0)).is_err());
    for maximum in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            plastic
                .clone()
                .with_norton_creep(kinetics(1.0), maximum)
                .is_err()
        );
    }
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(combined.at_temperature_for_step(300.0, dt).is_err());
    }
    let limited = plastic.with_norton_creep(kinetics(1.0), 1e-5).unwrap();
    assert!(
        limited
            .checked_increment_for_step(&ID, &virgin, 300.0, &overstressed, 300.0, 1e-3, 1.0)
            .is_err()
    );
    assert_eq!(virgin, combined.initial_state());
}

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", temperature_quantity(), 300.0)
        .unwrap()
}
fn card() -> ConstitutiveModelCard {
    let p = parameters();
    let k = kinetics(2.4);
    let pa = Dims([-1, 1, -2, 0, 0, 0]);
    let pa_k = Dims([-1, 1, -2, -1, 0, 0]);
    let kelvin = Dims([0, 0, 0, 1, 0, 0]);
    let hash = fs_blake3::hash_bytes(b"authored simultaneous finite Norton/J2 thermal fixture");
    ConstitutiveModelCard {
        law: LawId(THERMAL_FINITE_NORTON_J2_LAW_ID.into()),
        law_version: 1,
        parameters: [
            ("reference_temperature_k", p.reference_temperature_k, kelvin),
            ("minimum_temperature_k", p.temperature_range_k[0], kelvin),
            ("maximum_temperature_k", p.temperature_range_k[1], kelvin),
            ("shear_modulus_pa", p.reference.shear_modulus_pa, pa),
            ("bulk_modulus_pa", p.reference.bulk_modulus_pa, pa),
            ("yield_stress_pa", p.reference.yield_stress_pa, pa),
            ("hardening_modulus_pa", p.reference.hardening_modulus_pa, pa),
            (
                "shear_modulus_slope_pa_k",
                p.slopes.shear_modulus_pa_k,
                pa_k,
            ),
            ("bulk_modulus_slope_pa_k", p.slopes.bulk_modulus_pa_k, pa_k),
            ("yield_stress_slope_pa_k", p.slopes.yield_stress_pa_k, pa_k),
            (
                "hardening_modulus_slope_pa_k",
                p.slopes.hardening_modulus_pa_k,
                pa_k,
            ),
            (
                "logarithmic_expansion_per_k",
                p.logarithmic_expansion_per_k,
                Dims([0, 0, 0, -1, 0, 0]),
            ),
            (
                "logarithmic_expansion_slope_per_k2",
                p.logarithmic_expansion_slope_per_k2,
                Dims([0, 0, 0, -2, 0, 0]),
            ),
            (
                "minimum_mechanical_jacobian",
                LIMITS.minimum_jacobian,
                Dims::NONE,
            ),
            (
                "maximum_elastic_log_strain",
                LIMITS.maximum_elastic_log_strain,
                Dims::NONE,
            ),
            (
                "maximum_plastic_log_strain",
                LIMITS.maximum_plastic_log_strain,
                Dims::NONE,
            ),
            (
                "maximum_accumulated_plastic_strain",
                LIMITS.maximum_accumulated_plastic_strain,
                Dims::NONE,
            ),
            (
                "creep_reference_viscosity_pa_s",
                k.reference_viscosity_pa_s,
                Dims([-1, 1, -1, 0, 0, 0]),
            ),
            ("creep_stress_scale_pa", k.stress_scale_pa, pa),
            ("creep_exponent", k.exponent, Dims::NONE),
            (
                "creep_activation_temperature_k",
                k.activation_temperature_k,
                kelvin,
            ),
            ("maximum_accumulated_creep_strain", 3.0, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 2,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with_quantity(
            "temperature",
            temperature_quantity(),
            280.0,
            500.0,
        ),
        sources: vec![hash],
        provenance: Provenance {
            source: "authored simultaneous finite Norton/J2 thermal fixture".into(),
            license: "test fixture".into(),
            artifact: Some(hash),
        },
    }
}

#[test]
fn explicit_twenty_two_parameter_source_preserves_identity_and_executes_both_mechanisms() {
    let card = card();
    assert_eq!(card.parameters.len(), 22);
    let context = context();
    let source =
        ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&card, &context, "temperature")
            .unwrap();
    assert_eq!(source.model_card(), &card);
    assert_eq!(source.model_card().content_hash(), card.content_hash());
    assert_eq!(source.context(), &context);
    assert_eq!(source.temperature_axis(), "temperature");
    let native = family(2.4);
    assert_eq!(source.law(), &native);
    let f = coaxial(0.02, 0.15);
    let result = source
        .law()
        .checked_increment_for_step(&ID, &native.initial_state(), 300.0, &f, 410.0, 1e-3, 0.02)
        .unwrap();
    assert!(result.mechanical.plastic_increment > 0.0);
    assert!(result.mechanical.creep_increment > 0.0);
    assert!(result.mechanical.plastic_dissipation_j_m3 > 0.0);
    assert!(result.mechanical.creep_dissipation_j_m3 > 0.0);
    assert_eq!(
        result,
        native
            .checked_increment_for_step(&ID, &native.initial_state(), 300.0, &f, 410.0, 1e-3, 0.02)
            .unwrap()
    );
}

#[test]
fn combined_source_requires_exact_keys_units_state_v2_and_complete_thermal_coverage() {
    let original = card();
    let resolve = |card: &ConstitutiveModelCard| {
        ResolvedAffineThermalFiniteJ2::try_from_model_card_at(card, &context(), "temperature")
    };
    for name in original.parameters.keys() {
        let mut missing = original.clone();
        missing.parameters.remove(name);
        assert!(resolve(&missing).is_err(), "missing {name} accepted");
        let mut wrong = original.clone();
        wrong.parameters.get_mut(name).unwrap().dims = Dims([1, 0, 0, 0, 0, 0]);
        assert!(
            resolve(&wrong).is_err(),
            "wrong dimensions for {name} accepted"
        );
    }
    for version in [0, 1, 3] {
        let mut invalid = original.clone();
        invalid.state_schema_version = version;
        assert!(resolve(&invalid).is_err());
    }
    let mut version = original.clone();
    version.law_version = 2;
    assert!(resolve(&version).is_err());
    for name in [
        "overstress_exponent",
        "creep_viscosity",
        "maximum_creep_strain",
    ] {
        let mut extra = original.clone();
        extra.parameters.insert(
            name.into(),
            LawParameter {
                value: 1.0,
                dims: Dims::NONE,
            },
        );
        assert!(resolve(&extra).is_err());
    }
    for (name, value) in [
        ("maximum_accumulated_creep_strain", 0.0),
        ("yield_stress_pa", 0.0),
        ("hardening_modulus_slope_pa_k", -100_000.0),
        ("creep_exponent", 0.9),
        ("creep_activation_temperature_k", 1e9),
    ] {
        let mut invalid = original.clone();
        invalid.parameters.get_mut(name).unwrap().value = value;
        assert!(
            resolve(&invalid).is_err(),
            "invalid {name}={value} accepted"
        );
    }
    let mut narrow = original.clone();
    narrow.validity = ValidityDomain::unconstrained().with_quantity(
        "temperature",
        temperature_quantity(),
        300.0,
        450.0,
    );
    assert!(resolve(&narrow).is_err());
    assert!(
        ResolvedAffineThermalFiniteJ2::try_from_model_card_at(
            &original,
            &context().with("temperature", 310.0).unwrap(),
            "temperature"
        )
        .is_err()
    );
    for name in [
        THERMAL_FINITE_J2_LAW_ID,
        THERMAL_PERZYNA_FINITE_J2_LAW_ID,
        THERMAL_FINITE_NORTON_LAW_ID,
    ] {
        let mut mismatched = original.clone();
        mismatched.law = LawId(name.into());
        assert!(resolve(&mismatched).is_err());
    }
    let mut legacy = original;
    legacy.law = LawId(THERMAL_FINITE_J2_LAW_ID.into());
    legacy.parameters.retain(|name, _| {
        !name.starts_with("creep_") && name != "maximum_accumulated_creep_strain"
    });
    legacy.state_schema_version = 1;
    assert!(resolve(&legacy).is_ok());
    legacy.state_schema_version = 2;
    assert!(resolve(&legacy).is_err());
}
