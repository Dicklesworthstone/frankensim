//! Explicit-time finite thermoplasticity: kinetic physics and strict source transport.
//! These authored fixtures verify the law, not experimental material calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::plastic::finite::thermal::source::{
    ResolvedAffineThermalFiniteJ2, THERMAL_FINITE_J2_LAW_ID, THERMAL_PERZYNA_FINITE_J2_LAW_ID,
};
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalError,
    FiniteJ2ThermalOverstress, FiniteJ2ThermalParameters,
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
fn constant() -> FiniteJ2ThermalParameters {
    FiniteJ2ThermalParameters {
        slopes: FiniteJ2TemperatureSlopes {
            shear_modulus_pa_k: 0.0,
            bulk_modulus_pa_k: 0.0,
            yield_stress_pa_k: 0.0,
            hardening_modulus_pa_k: 0.0,
        },
        logarithmic_expansion_per_k: 0.0,
        logarithmic_expansion_slope_per_k2: 0.0,
        ..parameters()
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
fn coaxial(volume_log: f64, axial: f64) -> [f64; 9] {
    [
        (volume_log + axial).exp(),
        0.0,
        0.0,
        0.0,
        (volume_log - 0.5 * axial).exp(),
        0.0,
        0.0,
        0.0,
        (volume_log - 0.5 * axial).exp(),
    ]
}
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

#[test]
fn explicit_duration_preserves_legacy_arithmetic_and_cannot_be_silently_invented() {
    let elastic = AffineThermalFiniteJ2::new(parameters(), LIMITS).unwrap();
    let f = coaxial(0.03, 0.16);
    let old = elastic.initial_state();
    let legacy = elastic
        .checked_increment(&ID, &old, 300.0, &f, 380.0, 1e-3)
        .unwrap();
    for duration in [1e-8, 0.01, 1.0, 1e6] {
        assert_eq!(
            elastic
                .checked_increment_for_step(&ID, &old, 300.0, &f, 380.0, 1e-3, duration)
                .unwrap(),
            legacy
        );
    }
    let kinetic = elastic.clone().with_overstress(kinetics(2.0)).unwrap();
    assert!(matches!(
        kinetic.at_temperature(300.0),
        Err(FiniteJ2ThermalError::DurationRequired)
    ));
    assert!(matches!(
        kinetic.checked_increment(&ID, &old, 300.0, &f, 380.0, 1e-3),
        Err(FiniteJ2ThermalError::DurationRequired)
    ));
    let inspector = kinetic.at_temperature_for_inspection(300.0).unwrap();
    assert!(matches!(
        inspector.response(&f, &old),
        Err(FiniteJ2ThermalError::DurationRequired)
    ));
    assert!(matches!(
        inspector.stress_and_potential(&f, &old),
        Err(FiniteJ2ThermalError::DurationRequired)
    ));
    // An unrelaxed but finite state is valid kinetic history, not equilibrium yield.
    inspector.validate_committed(&f, &old, 0.0).unwrap();
    assert!(elastic.validate_committed(&f, &old, 300.0, 0.0).is_err());
    for duration in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(kinetic.at_temperature_for_step(300.0, duration).is_err());
    }
    let mut invalid = old.clone();
    invalid.inverse_plastic_cauchy_green[0] = -1.0;
    assert!(inspector.validate_committed(&f, &invalid, 0.0).is_err());
    assert!(
        inspector
            .validate_committed(&ID.map(|x| -x), &old, 0.0)
            .is_err()
    );
    assert!(inspector.validate_committed(&f, &old, f64::NAN).is_err());
    assert!(kinetic.at_temperature_for_step(501.0, 1.0).is_err());
    assert_eq!(old, elastic.initial_state());
}

#[test]
fn arrhenius_activation_accelerates_the_exact_linear_implicit_flow_at_nonunit_volume() {
    let p = constant();
    let k = kinetics(1.0);
    let family = AffineThermalFiniteJ2::new(p, LIMITS)
        .unwrap()
        .with_overstress(k)
        .unwrap();
    let (axial, volume, duration) = (0.16, 0.04, 0.4);
    let f = coaxial(volume, axial);
    let old = family.initial_state();
    let mut previous_increment = 0.0;
    for temperature in [280.0, 300.0, 400.0, 500.0] {
        let eta = k.reference_viscosity_pa_s
            * (k.activation_temperature_k * (1.0 / temperature - 1.0 / 300.0)).exp();
        near(
            family
                .overstress_viscosity_pa_s(temperature)
                .unwrap()
                .unwrap(),
            eta,
            1e-6,
            4e-14,
        );
        let mu = p.reference.shear_modulus_pa;
        let hardening = p.reference.hardening_modulus_pa;
        let yield_stress = p.reference.yield_stress_pa;
        let gamma = (3.0 * mu * axial - yield_stress) / (3.0 * mu + hardening + eta / duration);
        let viscous = eta * gamma / duration;
        let law = family
            .at_temperature_for_step(temperature, duration)
            .unwrap();
        let r = law.response(&f, &old).unwrap();
        near(r.plastic_increment, gamma, 2e-14, 2e-12);
        assert!(r.plastic_increment > previous_increment);
        previous_increment = r.plastic_increment;
        near(
            r.equivalent_kirchhoff_pa,
            yield_stress + hardening * gamma + viscous,
            1e-6,
            2e-12,
        );
        near(r.jacobian, (3.0 * volume).exp(), 1e-14, 1e-14);
        near(
            r.plastic_dissipation_j_m3,
            (yield_stress + viscous) * gamma,
            1e-7,
            2e-12,
        );
        let stored = 1.5 * mu * (axial - gamma).powi(2)
            + 0.5 * p.reference.bulk_modulus_pa * (3.0 * volume).powi(2)
            + 0.5 * hardening * gamma * gamma;
        near(r.stored_energy_j_m3, stored, 1e-6, 2e-12);
        near(
            r.incremental_potential_j_m3,
            stored + yield_stress * gamma + 0.5 * viscous * gamma,
            1e-6,
            2e-12,
        );
        assert!(r.incremental_potential_j_m3 < stored + r.plastic_dissipation_j_m3);
    }
    assert_eq!(
        family.overstress_viscosity_pa_s(300.0).unwrap(),
        Some(k.reference_viscosity_pa_s)
    );
}

#[test]
fn fixed_deformation_relaxes_accepted_overstress_and_heats_physical_loss_once() {
    let p = constant();
    let k = kinetics(1.0);
    let family = AffineThermalFiniteJ2::new(p, LIMITS)
        .unwrap()
        .with_overstress(k)
        .unwrap();
    let f = coaxial(0.02, 0.18);
    let loaded = family
        .at_temperature_for_step(300.0, 0.03)
        .unwrap()
        .response(&f, &family.initial_state())
        .unwrap();
    let old = loaded.state.clone();
    family.validate_committed(&f, &old, 300.0, 0.0).unwrap();
    let duration = 0.6;
    let step = family
        .checked_increment_for_step(&f, &old, 300.0, &f, 300.0, 0.0, duration)
        .unwrap();
    let excess = loaded.equivalent_kirchhoff_pa
        - p.reference.yield_stress_pa
        - p.reference.hardening_modulus_pa * old.alpha;
    let gamma = excess
        / (3.0 * p.reference.shear_modulus_pa
            + p.reference.hardening_modulus_pa
            + k.reference_viscosity_pa_s / duration);
    let viscous = k.reference_viscosity_pa_s * gamma / duration;
    let physical = (p.reference.yield_stress_pa + viscous) * gamma;
    near(step.mechanical.plastic_increment, gamma, 1e-14, 2e-12);
    near(step.caloric_heat_j_m3(), physical, 1e-6, 2e-12);
    assert!(step.mechanical.state.alpha > old.alpha);
    assert!(step.mechanical.equivalent_kirchhoff_pa < loaded.equivalent_kirchhoff_pa);
    assert_eq!(step.mechanical_work_j_m3, 0.0);
    assert_eq!(step.reversible_heat_j_m3, 0.0);
    assert!(step.internal_energy_change_j_m3() < 0.0);
    assert!(step.endpoint_work_defect_j_m3 > 0.0);
    near(
        step.internal_energy_change_j_m3() + physical + step.endpoint_work_defect_j_m3,
        0.0,
        1e-6,
        0.0,
    );
    near(step.balance_error_j_m3, 0.0, 1e-6, 0.0);
    let reduced_loss = step.mechanical.incremental_potential_j_m3 - step.new.free_energy_j_m3;
    near(physical - reduced_loss, 0.5 * viscous * gamma, 1e-6, 2e-12);
    assert_eq!(old, loaded.state);
}

#[test]
fn viscosity_has_no_fake_thermodynamic_storage_and_two_temperature_work_closes() {
    let p = parameters();
    let base = AffineThermalFiniteJ2::new(p, LIMITS).unwrap();
    let family = base.clone().with_overstress(kinetics(2.5)).unwrap();
    let changed = base
        .clone()
        .with_overstress(FiniteJ2ThermalOverstress {
            reference_viscosity_pa_s: 5e3,
            activation_temperature_k: 0.0,
            ..kinetics(4.0)
        })
        .unwrap();
    let old_f = coaxial(0.01, 0.1);
    let history = family
        .at_temperature_for_step(300.0, 0.05)
        .unwrap()
        .response(&old_f, &family.initial_state())
        .unwrap()
        .state;
    let new_f = [1.22, 0.16, 0.01, -0.03, 0.94, 0.05, 0.03, -0.02, 1.02];
    for temperature in [280.0, 370.0, 500.0] {
        let expected = base.thermodynamics(&new_f, &history, temperature).unwrap();
        assert_eq!(
            family
                .thermodynamics(&new_f, &history, temperature)
                .unwrap(),
            expected
        );
        assert_eq!(
            changed
                .thermodynamics(&new_f, &history, temperature)
                .unwrap(),
            expected
        );
    }
    assert_eq!(
        family.minimum_mechanical_heat_capacity_j_m3_k(),
        base.minimum_mechanical_heat_capacity_j_m3_k()
    );
    let step = family
        .checked_increment_for_step(&old_f, &history, 300.0, &new_f, 410.0, 1e-3, 0.25)
        .unwrap();
    let total_work: f64 = step
        .mechanical
        .first_piola_pa
        .iter()
        .enumerate()
        .map(|(i, p)| p * (new_f[i] - old_f[i]))
        .sum();
    near(
        total_work,
        step.new.internal_energy_j_m3 - step.old.internal_energy_j_m3
            + step.caloric_heat_j_m3()
            + step.endpoint_work_defect_j_m3,
        2e-6,
        2e-13,
    );
    near(step.balance_error_j_m3, 0.0, 2e-6, 0.0);
    assert!(step.mechanical.plastic_increment > 0.0);
    let gamma = step.mechanical.plastic_increment;
    let yield_stress = p
        .slopes
        .yield_stress_pa_k
        .mul_add(110.0, p.reference.yield_stress_pa);
    let eta = family.overstress_viscosity_pa_s(410.0).unwrap().unwrap();
    let k = family.overstress().unwrap();
    let viscous =
        k.stress_scale_pa * (eta * gamma / (0.25 * k.stress_scale_pa)).powf(1.0 / k.exponent);
    near(
        step.mechanical.plastic_dissipation_j_m3,
        (yield_stress + viscous) * gamma,
        1e-6,
        3e-12,
    );
    near(
        step.mechanical.incremental_potential_j_m3,
        step.new.free_energy_j_m3
            + (yield_stress + k.exponent / (k.exponent + 1.0) * viscous) * gamma,
        2e-6,
        3e-12,
    );
}

#[test]
fn activated_noninteger_flow_preserves_total_f_potential_gradient_and_full_tangent() {
    let family = AffineThermalFiniteJ2::new(parameters(), LIMITS)
        .unwrap()
        .with_overstress(kinetics(2.4))
        .unwrap();
    let old = family
        .at_temperature_for_step(300.0, 0.1)
        .unwrap()
        .response(&coaxial(0.0, 0.1), &family.initial_state())
        .unwrap()
        .state;
    let law = family.at_temperature_for_step(410.0, 0.3).unwrap();
    let f = [1.2, 0.15, 0.03, -0.04, 0.94, 0.07, 0.02, -0.03, 1.01];
    let r = law.response(&f, &old).unwrap();
    assert!(r.plastic_increment > 0.0);
    assert_eq!(
        law.stress_and_potential(&f, &old).unwrap(),
        (r.first_piola_pa, r.incremental_potential_j_m3)
    );
    for col in 0..9 {
        let h = 1e-6;
        let mut plus = f;
        let mut minus = f;
        plus[col] += h;
        minus[col] -= h;
        let a = law.stress_and_potential(&plus, &old).unwrap();
        let b = law.stress_and_potential(&minus, &old).unwrap();
        near((a.1 - b.1) / (2.0 * h), r.first_piola_pa[col], 0.1, 4e-7);
        for row in 0..9 {
            near(
                (a.0[row] - b.0[row]) / (2.0 * h),
                r.tangent_pa[row][col],
                8.0,
                3e-6,
            );
            near(r.tangent_pa[row][col], r.tangent_pa[col][row], 2e-6, 3e-12);
        }
    }
}

fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", absolute_temperature(), 300.0)
        .unwrap()
}
fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn card() -> ConstitutiveModelCard {
    let p = parameters();
    let k = kinetics(2.5);
    let pa = Dims([-1, 1, -2, 0, 0, 0]);
    let pa_k = Dims([-1, 1, -2, -1, 0, 0]);
    let kelvin = Dims([0, 0, 0, 1, 0, 0]);
    let source =
        fs_blake3::hash_bytes(b"authored finite Perzyna thermal fixture; no measured calibration");
    ConstitutiveModelCard {
        law: LawId(THERMAL_PERZYNA_FINITE_J2_LAW_ID.into()),
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
                "overstress_reference_viscosity_pa_s",
                k.reference_viscosity_pa_s,
                Dims([-1, 1, -1, 0, 0, 0]),
            ),
            ("overstress_stress_scale_pa", k.stress_scale_pa, pa),
            ("overstress_exponent", k.exponent, Dims::NONE),
            (
                "overstress_activation_temperature_k",
                k.activation_temperature_k,
                kelvin,
            ),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with_quantity(
            "temperature",
            absolute_temperature(),
            p.temperature_range_k[0],
            p.temperature_range_k[1],
        ),
        sources: vec![source],
        provenance: Provenance {
            source: "authored synthetic finite Perzyna thermal fixture".into(),
            license: "test fixture".into(),
            artifact: Some(source),
        },
    }
}

#[test]
fn explicit_twenty_one_parameter_source_retains_identity_and_drives_native_kinetics() {
    let card = card();
    let context = context();
    let source =
        ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&card, &context, "temperature")
            .unwrap();
    assert_eq!(source.model_card(), &card);
    assert_eq!(source.model_card().content_hash(), card.content_hash());
    assert_eq!(source.context(), &context);
    assert_eq!(source.temperature_axis(), "temperature");
    assert_eq!(source.law().overstress(), Some(&kinetics(2.5)));
    assert_eq!(source.law().parameters(), &parameters());
    assert_eq!(source.law().limits(), &LIMITS);
    let native = AffineThermalFiniteJ2::new(parameters(), LIMITS)
        .unwrap()
        .with_overstress(kinetics(2.5))
        .unwrap();
    let f = coaxial(0.02, 0.15);
    assert_eq!(
        source
            .law()
            .checked_increment_for_step(&ID, &native.initial_state(), 300.0, &f, 410.0, 1e-3, 0.4)
            .unwrap(),
        native
            .checked_increment_for_step(&ID, &native.initial_state(), 300.0, &f, 410.0, 1e-3, 0.4)
            .unwrap()
    );
    assert!(source.law().at_temperature(410.0).is_err());
}

#[test]
fn kinetic_source_requires_explicit_units_parameters_full_span_and_original_reference() {
    let original = card();
    let resolve = |candidate: &ConstitutiveModelCard| {
        ResolvedAffineThermalFiniteJ2::try_from_model_card_at(candidate, &context(), "temperature")
    };
    for name in [
        "overstress_reference_viscosity_pa_s",
        "overstress_stress_scale_pa",
        "overstress_exponent",
        "overstress_activation_temperature_k",
    ] {
        let mut absent = original.clone();
        absent.parameters.remove(name);
        assert!(resolve(&absent).is_err());
        let mut wrong_units = original.clone();
        wrong_units.parameters.get_mut(name).unwrap().dims = Dims([1, 0, 0, 0, 0, 0]);
        assert!(resolve(&wrong_units).is_err());
    }
    for (name, value) in [
        ("overstress_reference_viscosity_pa_s", 0.0),
        ("overstress_stress_scale_pa", -1.0),
        ("overstress_exponent", 0.999),
        ("overstress_activation_temperature_k", -1.0),
        ("overstress_activation_temperature_k", 1e9),
    ] {
        let mut invalid = original.clone();
        invalid.parameters.get_mut(name).unwrap().value = value;
        assert!(resolve(&invalid).is_err(), "accepted {name}={value}");
    }
    let mut legacy = original.clone();
    legacy.law = LawId(THERMAL_FINITE_J2_LAW_ID.into());
    assert!(resolve(&legacy).is_err());
    let mut truncated = original.clone();
    truncated.validity = ValidityDomain::unconstrained().with_quantity(
        "temperature",
        absolute_temperature(),
        290.0,
        450.0,
    );
    assert!(resolve(&truncated).is_err());
    let mut version = original.clone();
    version.law_version = 2;
    assert!(resolve(&version).is_err());
    version.law_version = 1;
    version.state_schema_version = 2;
    assert!(resolve(&version).is_err());
    let shifted = context().with("temperature", 310.0).unwrap();
    assert!(
        ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&original, &shifted, "temperature")
            .is_err()
    );
    // The reference value can be finite while a distant endpoint underflows.
    let tiny = AffineThermalFiniteJ2::new(parameters(), LIMITS)
        .unwrap()
        .with_overstress(FiniteJ2ThermalOverstress {
            reference_viscosity_pa_s: 1e-280,
            activation_temperature_k: 100_000.0,
            ..kinetics(2.0)
        });
    assert!(tiny.is_err());
}
