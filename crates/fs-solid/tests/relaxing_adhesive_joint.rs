//! Original-source adhesive relaxation and creep in actual deformable members.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::cohesive::relaxation::{ViscoelasticBondLimits, source::*};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermovisco::source::{ResolvedThermoMaxwell, ThermoMaxwellSourceQuery};
use fs_material::thermovisco::{HistoryIntegration, Limits};
use fs_qty::{Dims, Length, Pressure, QuantitySpec, Temperature, Time};
use fs_solid::cohesive::thermal_joint::{relaxation::*, source::SourceMember};
use std::collections::BTreeMap;

const AREA: f64 = 1e-4;
const LENGTH: f64 = 0.005;
const THICKNESS: f64 = 0.001;
const E_INF: f64 = 1e6;
const E_BRANCH: f64 = 2e6;
const TAU: f64 = 1.0;

fn query(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn claim(name: &str, dims: Dims, values: [f64; 2]) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(280.0, values[0]), (520.0, values[1])],
            dims,
        },
        validity: ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            280.0,
            520.0,
        ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured series adhesive source; no measured validation claim".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn material(
    name: &str,
    e: f64,
    branch: bool,
    integration: HistoryIntegration,
) -> ResolvedThermoMaxwell {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        ("E", Pressure::DIMS, [e; 2]),
        ("branch", Pressure::DIMS, [2.0 * e; 2]),
        ("tau", Time::DIMS, [2.0; 2]),
        ("scale", Dims::NONE, [1.0; 2]),
        ("shift", Dims::NONE, [1.0; 2]),
        ("free", Dims::NONE, [-0.002, 0.022]),
    ] {
        claims.insert_claim(claim(name, dims, values)).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "manufactured relaxing joint".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let branches = if branch {
        vec![(
            PropertyKey::new("branch", Pressure::DIMS),
            PropertyKey::new("tau", Time::DIMS),
        )]
    } else {
        vec![]
    };
    ResolvedThermoMaxwell::try_new(
        &card,
        ThermoMaxwellSourceQuery {
            equilibrium_modulus: &PropertyKey::new("E", Pressure::DIMS),
            branches: &branches,
            modulus_scale: &PropertyKey::new("scale", Dims::NONE),
            relaxation_scale: &PropertyKey::new("shift", Dims::NONE),
            free_strain: &PropertyKey::new("free", Dims::NONE),
            temperature_axis: "T",
            lower: &query(280.0),
            upper: &query(520.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        Limits {
            max_strain: 0.04,
            max_strain_rate_s: 1.0,
        },
        || false,
    )
    .unwrap()
    .with_history_integration(integration)
}
fn bond(a: &ResolvedThermoMaxwell, b: &ResolvedThermoMaxwell) -> ResolvedViscoelasticBond {
    resolve_bond(
        &bond_card(a, b, [[E_INF; 2], [E_BRANCH; 2], [TAU; 2]]),
        &query(300.0),
    )
}
fn bond_card(
    a: &ResolvedThermoMaxwell,
    b: &ResolvedThermoMaxwell,
    coefficients: [[f64; 2]; 3],
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(claim(
            BONDLINE_THICKNESS_PROPERTY,
            Length::DIMS,
            [THICKNESS; 2],
        ))
        .unwrap();
    for prefix in [NORMAL_RELAXATION_PREFIX, SHEAR_RELAXATION_PREFIX] {
        for (suffix, dims, values) in [
            ("branch_count", Dims::NONE, [1.0; 2]),
            ("equilibrium_modulus", Pressure::DIMS, coefficients[0]),
            ("modulus_1", Pressure::DIMS, coefficients[1]),
            ("time_1", Time::DIMS, coefficients[2]),
        ] {
            claims
                .insert_claim(claim(&format!("{prefix}_{suffix}"), dims, values))
                .unwrap();
        }
    }
    let surface = |s: &ResolvedThermoMaxwell| SurfaceSpec {
        material: s.samples()[0].material().clone(),
        texture_frame: "fixed axial frame".into(),
    };
    InterfaceSystemCard::assemble(
        surface(a),
        surface(b),
        SystemContext {
            medium: "explicit intact adhesive".into(),
            third_body: None,
            environment: "manufactured prescribed bond temperature".into(),
            history: "stress-free bonded birth".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn resolve_bond(card: &InterfaceSystemCard, point: &QueryPoint) -> ResolvedViscoelasticBond {
    resolve_viscoelastic_bond(
        card,
        point,
        MaterialPropertySelection::SingleClaimOnly,
        ViscoelasticBondLimits {
            maximum_jump_m: 2e-5,
            maximum_rate_m_s: 0.05,
            maximum_elapsed_s: 100.0,
        },
    )
    .unwrap()
}
fn members<'a>(
    a: &'a ResolvedThermoMaxwell,
    b: &'a ResolvedThermoMaxwell,
) -> Vec<SourceMember<'a>> {
    [a, b]
        .into_iter()
        .enumerate()
        .map(|(i, material)| SourceMember {
            id: i as u64 + 1,
            reference_length_m: LENGTH,
            reference_area_m2: AREA,
            material,
            birth_temperature_k: 500.0,
        })
        .collect()
}
fn budget() -> Budget {
    Budget {
        maximum_members: 4,
        maximum_branch_states: 16,
    }
}
fn config() -> Config {
    Config {
        force_tolerance_n: 1e-9,
        extension_tolerance_m: 1e-14,
        local_bond_energy_tolerance_j_m2: 1e-9,
        energy_tolerance_j: 1e-12,
    }
}
fn initial<'a>(
    a: &'a ResolvedThermoMaxwell,
    b: &'a ResolvedThermoMaxwell,
    bond: &'a ResolvedViscoelasticBond,
) -> SourceRelaxingJoint<'a> {
    SourceRelaxingJoint::try_new(0.0, members(a, b), 7, AREA, bond, [1, 2], budget(), || {
        false
    })
    .unwrap()
}
fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "{actual:e} != {expected:e}"
    );
}
fn balance(old: &SourceRelaxingJoint<'_>, next: &SourceRelaxingJoint<'_>, report: &Report) {
    let bulk = &report.bulk_energy_j;
    near(
        next.stored_energy_j() - old.stored_energy_j(),
        bulk.stored_change + report.bond_stored_change_j,
        1e-10,
        1e-15,
    );
    let left = bulk.stored_change
        + report.bond_stored_change_j
        + bulk.dissipated
        + bulk.numerical_dissipation
        + report.bond_relaxation_heat_j
        + report.bond_work_correction_j
        + report.bulk_work_correction_j;
    let right = report.endpoint_work_j + bulk.thermal_strain_work + bulk.modulus_work;
    near(left, right, 1e-10, 1e-15);
    near(
        report.bond_work_j,
        report.bond_stored_change_j + report.bond_relaxation_heat_j,
        1e-11,
        1e-15,
    );
    assert!(report.energy_residual_j.abs() <= config().energy_tolerance_j);
    assert!(
        report.bond_relaxation_heat_j >= 0.0
            && bulk.dissipated >= 0.0
            && bulk.numerical_dissipation >= 0.0
    );
}

#[test]
fn restrained_sls_bond_relaxes_with_the_independent_series_time_and_heat() {
    let a = material("A", 1e7, false, HistoryIntegration::BackwardEuler);
    let b = material("B", 2e7, false, HistoryIntegration::BackwardEuler);
    let bond = bond(&a, &b);
    let temperatures = BTreeMap::from([(1, 500.0), (2, 500.0)]);
    let extension = 1e-6;
    let (loaded, first) = initial(&a, &b, &bond)
        .advance_controlled(
            0.1,
            Control::ExtensionM(extension),
            &temperatures,
            &query(300.0),
            config(),
            || false,
        )
        .unwrap();
    let ks = (LENGTH / (AREA * 1e7) + LENGTH / (AREA * 2e7)).recip();
    let kinf = AREA * E_INF / THICKNESS;
    let k0 = AREA * (E_INF + E_BRANCH) / THICKNESS;
    let time = TAU * (ks + k0) / (ks + kinf);
    let equilibrium = extension * ks * kinf / (ks + kinf);
    let duration = 2.5;
    let expected = equilibrium + (first.force_n - equilibrium) * (-duration / time).exp();
    let q0 = (first.force_n / AREA - E_INF * first.opening_m / THICKNESS) / E_BRANCH;
    let expected_heat = AREA * THICKNESS * E_BRANCH * q0 * q0 * time / (2.0 * TAU)
        * (1.0 - (-2.0 * duration / time).exp());
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let mut state = loaded.clone();
        let mut last = first.force_n;
        let mut heat = 0.0;
        for i in 1..=n {
            let (next, report) = state
                .advance_controlled(
                    0.1 + duration * f64::from(i) / f64::from(n),
                    Control::ExtensionM(extension),
                    &temperatures,
                    &query(300.0),
                    config(),
                    || false,
                )
                .unwrap();
            assert!(report.force_n < last && report.force_n > equilibrium);
            assert!(report.opening_m > state.bond().jump_m()[0]);
            assert_eq!(report.endpoint_work_j, 0.0);
            balance(&state, &next, &report);
            heat += report.bond_relaxation_heat_j;
            last = report.force_n;
            state = next;
        }
        errors.push((last - expected).abs());
        if n == 32 {
            near(heat, expected_heat, 3e-4, 1e-15);
        }
        assert!(core::ptr::eq(state.interface_source(), &bond));
        assert!(core::ptr::eq(state.member_sources()[&1], &a));
    }
    assert!(
        errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8,
        "{errors:?}"
    );
    near(errors[2], 0.0, 0.0, 3e-5 * expected.abs());
}

#[test]
fn constant_force_creep_uses_original_adhesive_memory_and_retardation_time() {
    let a = material("A", 1e7, false, HistoryIntegration::ExactStrainRamp);
    let b = material("B", 2e7, false, HistoryIntegration::ExactStrainRamp);
    let bond = bond(&a, &b);
    let temperatures = BTreeMap::from([(1, 500.0), (2, 500.0)]);
    let force = 0.1;
    let (loaded, first) = initial(&a, &b, &bond)
        .advance_controlled(
            0.1,
            Control::ForceN(force),
            &temperatures,
            &query(300.0),
            config(),
            || false,
        )
        .unwrap();
    assert_eq!(first.bulk_energy_j.numerical_dissipation, 0.0);
    assert!(first.bulk_work_correction_j > 0.0);
    let g_inf = force * THICKNESS / (AREA * E_INF);
    let time = TAU * (E_INF + E_BRANCH) / E_INF;
    let duration = 2.5;
    let expected = g_inf + (first.opening_m - g_inf) * (-duration / time).exp();
    let q0 = (force / AREA - E_INF * first.opening_m / THICKNESS) / E_BRANCH;
    let expected_heat = AREA * THICKNESS * E_BRANCH * q0 * q0 * time / (2.0 * TAU)
        * (1.0 - (-2.0 * duration / time).exp());
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let mut state = loaded.clone();
        let mut heat = 0.0;
        for i in 1..=n {
            let (next, report) = state
                .advance_controlled(
                    0.1 + duration * f64::from(i) / f64::from(n),
                    Control::ForceN(force),
                    &temperatures,
                    &query(300.0),
                    config(),
                    || false,
                )
                .unwrap();
            near(report.force_n, force, 1e-12, 1e-14);
            assert!(next.extension_m() > state.extension_m());
            near(
                next.members()[&1].point().total_strain(),
                force / (AREA * 1e7),
                1e-12,
                1e-15,
            );
            near(
                report.endpoint_work_j,
                force * (next.extension_m() - state.extension_m()),
                1e-12,
                1e-15,
            );
            balance(&state, &next, &report);
            heat += report.bond_relaxation_heat_j;
            state = next;
        }
        errors.push((state.bond().jump_m()[0] - expected).abs());
        if n == 32 {
            near(heat, expected_heat, 3e-4, 1e-15);
        }
    }
    assert!(
        errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8,
        "{errors:?}"
    );
    near(errors[2], 0.0, 0.0, 5e-5 * expected.abs());
}

#[test]
fn heated_unequal_members_retain_native_branch_updates_and_signed_work_accounts() {
    let a = material("A", 1e7, true, HistoryIntegration::ExactStrainRamp);
    let b = material("B", 2e7, true, HistoryIntegration::BackwardEuler);
    let bond = bond(&a, &b);
    let old = initial(&a, &b, &bond);
    let temperatures = BTreeMap::from([(1, 490.0), (2, 480.0)]);
    let dt = 0.2;
    let (next, report) = old
        .advance_controlled(
            dt,
            Control::ExtensionM(0.0),
            &temperatures,
            &query(300.0),
            config(),
            || false,
        )
        .unwrap();
    let psi = |x: f64| -(-x).exp_m1() / x;
    let ea = 1e7 * (1.0 + 2.0 * psi(dt / 2.0));
    let eb = 2e7 * (1.0 + 2.0 / (1.0 + dt / 2.0));
    let kb = AREA * (E_INF + E_BRANCH * psi(dt / TAU)) / THICKNESS;
    let free = LENGTH * 1e-4 * ((490.0 - 500.0) + (480.0 - 500.0));
    let force = -free / (LENGTH / (AREA * ea) + LENGTH / (AREA * eb) + 1.0 / kb);
    near(report.force_n, force, 2e-12, 1e-12);
    near(report.opening_m, force / kb, 2e-12, 1e-14);
    let strain_a = next.members()[&1].point().total_strain();
    let mechanical_a = strain_a - 1e-4 * (490.0 - 500.0);
    let mean_a = 1e7 * (0.5 + 2.0 * (1.0 - psi(dt / 2.0)) / (dt / 2.0)) * mechanical_a;
    near(
        report.bulk_work_correction_j,
        (force - AREA * mean_a) * LENGTH * strain_a,
        3e-11,
        1e-15,
    );
    near(
        next.members()[&1].point().branch_strains()[0],
        psi(dt / 2.0) * mechanical_a,
        2e-12,
        1e-15,
    );
    let strain_b = next.members()[&2].point().total_strain();
    near(
        next.members()[&2].point().branch_strains()[0],
        (strain_b - 1e-4 * (480.0 - 500.0)) / (1.0 + dt / 2.0),
        2e-12,
        1e-15,
    );
    assert!(
        report.bulk_energy_j.dissipated > 0.0 && report.bulk_energy_j.numerical_dissipation > 0.0
    );
    assert!(report.bond_relaxation_heat_j > 0.0 && report.bulk_work_correction_j < 0.0);
    balance(&old, &next, &report);
    let (held, after) = next
        .advance_controlled(
            0.4,
            Control::ExtensionM(0.0),
            &temperatures,
            &query(300.0),
            config(),
            || false,
        )
        .unwrap();
    assert!(after.force_n < report.force_n && after.bond_relaxation_heat_j > 0.0);
    assert!(core::ptr::eq(held.interface_source(), &bond));
    balance(&next, &held, &after);
}

#[test]
fn sources_domains_and_late_cancellation_leave_all_accepted_memory_retryable() {
    let a = material("A", 1e7, false, HistoryIntegration::BackwardEuler);
    let b = material("B", 2e7, false, HistoryIntegration::BackwardEuler);
    let bond = bond(&a, &b);
    let temperatures = BTreeMap::from([(1, 500.0), (2, 500.0)]);
    assert!(
        SourceRelaxingJoint::try_new(
            0.0,
            members(&a, &b),
            7,
            AREA,
            &bond,
            [2, 1],
            budget(),
            || false
        )
        .is_err()
    );
    assert!(
        SourceRelaxingJoint::try_new(
            0.0,
            members(&a, &b),
            7,
            AREA,
            &bond,
            [1, 2],
            Budget {
                maximum_branch_states: 2,
                ..budget()
            },
            || false
        )
        .is_err()
    );
    let (state, _) = initial(&a, &b, &bond)
        .advance_controlled(
            0.1,
            Control::ForceN(0.1),
            &temperatures,
            &query(300.0),
            config(),
            || false,
        )
        .unwrap();
    let before = (
        state.time_s(),
        state.extension_m(),
        state.bond().elapsed_s(),
        state.bond().jump_m(),
        state.bond().traction_pa(),
        state.bond().stored_energy_j_m2(),
    );
    let memory: Vec<_> = state
        .members()
        .values()
        .map(|m| m.point().clone())
        .collect();
    for (target, control, temperature) in [
        (0.1, Control::ForceN(0.1), 300.0),
        (101.0, Control::ForceN(0.1), 300.0),
        (0.2, Control::ForceN(f64::NAN), 300.0),
        (0.2, Control::ForceN(1e8), 300.0),
        (0.2, Control::ForceN(0.1), 310.0),
        (0.1 + 1e-10, Control::ForceN(1.0), 300.0),
    ] {
        assert!(
            state
                .advance_controlled(
                    target,
                    control,
                    &temperatures,
                    &query(temperature),
                    config(),
                    || false
                )
                .is_err()
        );
    }
    assert!(
        state
            .advance_controlled(
                0.2,
                Control::ForceN(0.1),
                &BTreeMap::from([(1, 270.0), (2, 500.0)]),
                &query(300.0),
                config(),
                || false
            )
            .is_err()
    );
    let mut polls = 0;
    let (expected, _) = state
        .advance_controlled(
            0.2,
            Control::ForceN(0.1),
            &temperatures,
            &query(300.0),
            config(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut calls = 0;
    assert!(
        state
            .advance_controlled(
                0.2,
                Control::ForceN(0.1),
                &temperatures,
                &query(300.0),
                config(),
                || {
                    calls += 1;
                    calls == polls
                }
            )
            .is_err()
    );
    assert_eq!(
        before,
        (
            state.time_s(),
            state.extension_m(),
            state.bond().elapsed_s(),
            state.bond().jump_m(),
            state.bond().traction_pa(),
            state.bond().stored_energy_j_m2()
        )
    );
    assert_eq!(
        memory,
        state
            .members()
            .values()
            .map(|m| m.point().clone())
            .collect::<Vec<_>>()
    );
    let (retry, _) = state
        .advance_controlled(
            0.2,
            Control::ForceN(0.1),
            &temperatures,
            &query(300.0),
            config(),
            || false,
        )
        .unwrap();
    assert_eq!(retry.bond().traction_pa(), expected.bond().traction_pa());
    assert_eq!(
        retry.bond().stored_energy_j_m2(),
        expected.bond().stored_energy_j_m2()
    );
}

fn temperature_balance(
    old: &SourceRelaxingJoint<'_>,
    next: &SourceRelaxingJoint<'_>,
    report: &TemperatureReport,
) {
    let r = &report.mechanical;
    let bulk = &r.bulk_energy_j;
    near(
        next.stored_energy_j() - old.stored_energy_j(),
        report.interface_stored_change_j + bulk.stored_change + r.bond_stored_change_j,
        2e-10,
        1e-15,
    );
    near(
        report.interface_parameter_work_j,
        report.interface_stored_change_j,
        1e-13,
        1e-20,
    );
    let stored = report.interface_stored_change_j + bulk.stored_change + r.bond_stored_change_j;
    let lost = bulk.dissipated + bulk.numerical_dissipation + r.bond_relaxation_heat_j;
    let work = r.endpoint_work_j - r.bond_work_correction_j - r.bulk_work_correction_j
        + bulk.thermal_strain_work
        + bulk.modulus_work
        + report.interface_parameter_work_j;
    near(stored + lost, work, 2e-10, 1e-15);
    assert!(report.energy_residual_j.abs() <= config().energy_tolerance_j);
    assert!(report.interface_energy_residual_j.abs() <= config().energy_tolerance_j);
}

#[test]
fn temperature_changes_real_series_response_without_resetting_adhesive_memory() {
    let a = material("A", 1e7, false, HistoryIntegration::ExactStrainRamp);
    let b = material("B", 2e7, false, HistoryIntegration::BackwardEuler);
    let card = bond_card(&a, &b, [[E_INF, 0.5e6], [E_BRANCH, 0.7e6], [TAU, 0.2]]);
    let cold = resolve_bond(&card, &query(280.0));
    let hot = cold.at_temperature(&card, "T", 520.0).unwrap();
    let temperatures = BTreeMap::from([(1, 500.0), (2, 500.0)]);
    let (loaded, _) = initial(&a, &b, &cold)
        .advance_controlled(
            0.1,
            Control::ExtensionM(1e-6),
            &temperatures,
            &query(280.0),
            config(),
            || false,
        )
        .unwrap();
    let compliance = LENGTH / (AREA * 1e7) + LENGTH / (AREA * 2e7);
    for control in [Control::ExtensionM(1e-6), Control::ForceN(0.1)] {
        let mut state = loaded.clone();
        // Heat, then cool the same retained spring strains. Both prescribed
        // extension and prescribed force must use the new source kinetics.
        for (source, dt) in [(&hot, 0.2), (&cold, 0.3)] {
            let old_spectrum = state.interface_source().law().normal_spectrum();
            let spectrum = source.law().normal_spectrum();
            let strain = state.bond().jump_m()[0] / THICKNESS;
            let q = (state.bond().traction_pa()[0] - old_spectrum.e_inf * strain)
                / old_spectrum.terms[0].0;
            let (modulus, tau) = spectrum.terms[0];
            let x = dt / tau;
            let decay = (-x).exp();
            let ramp = -(-x).exp_m1() / x;
            let held_force = AREA * (spectrum.e_inf * strain + modulus * decay * q);
            let stiffness = AREA * (spectrum.e_inf + modulus * ramp) / THICKNESS;
            let increment = match control {
                Control::ExtensionM(extension) => {
                    (extension - state.bond().jump_m()[0] - compliance * held_force)
                        / (1.0 + compliance * stiffness)
                }
                Control::ForceN(force) => (force - held_force) / stiffness,
            };
            let opening = state.bond().jump_m()[0] + increment;
            let expected_q = decay * q + ramp * increment / THICKNESS;
            let force = AREA * (spectrum.e_inf * opening / THICKNESS + modulus * expected_q);
            let parameter_work = 0.5
                * AREA
                * THICKNESS
                * ((spectrum.e_inf - old_spectrum.e_inf) * strain * strain
                    + (modulus - old_spectrum.terms[0].0) * q * q);
            let (next, report) = state
                .advance_with_interface_temperature_controlled(
                    state.time_s() + dt,
                    control,
                    &temperatures,
                    source,
                    "T",
                    config(),
                    || false,
                )
                .unwrap();
            near(report.mechanical.opening_m, opening, 1e-12, 1e-17);
            near(report.mechanical.force_n, force, 1e-12, 1e-14);
            near(
                report.interface_parameter_work_j,
                parameter_work,
                1e-12,
                1e-20,
            );
            let actual_q =
                (next.bond().traction_pa()[0] - spectrum.e_inf * opening / THICKNESS) / modulus;
            near(actual_q, expected_q, 2e-12, 1e-17);
            near(
                next.bond().elapsed_s(),
                state.bond().elapsed_s() + dt,
                1e-14,
                1e-16,
            );
            assert_eq!(next.interface_query(), &source.resolved().query());
            assert!(core::ptr::eq(next.interface_source(), source));
            assert!(core::ptr::eq(next.member_sources()[&1], &a));
            temperature_balance(&state, &next, &report);
            state = next;
        }
    }
}

#[test]
fn temperature_dependent_relaxation_time_matches_the_heated_series_decay() {
    let a = material("A", 1e7, false, HistoryIntegration::ExactStrainRamp);
    let b = material("B", 2e7, false, HistoryIntegration::ExactStrainRamp);
    let card = bond_card(&a, &b, [[E_INF; 2], [E_BRANCH; 2], [TAU, 0.2]]);
    let cold = resolve_bond(&card, &query(280.0));
    let hot = cold.at_temperature(&card, "T", 520.0).unwrap();
    let temperatures = BTreeMap::from([(1, 500.0), (2, 500.0)]);
    let extension = 1e-6;
    let (loaded, first) = initial(&a, &b, &cold)
        .advance_controlled(
            0.1,
            Control::ExtensionM(extension),
            &temperatures,
            &query(280.0),
            config(),
            || false,
        )
        .unwrap();
    let ks = (LENGTH / (AREA * 1e7) + LENGTH / (AREA * 2e7)).recip();
    let kinf = AREA * E_INF / THICKNESS;
    let k0 = AREA * (E_INF + E_BRANCH) / THICKNESS;
    let time = 0.2 * (ks + k0) / (ks + kinf);
    let equilibrium = extension * ks * kinf / (ks + kinf);
    let duration = 0.5;
    let expected = equilibrium + (first.force_n - equilibrium) * (-duration / time).exp();
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let mut state = loaded.clone();
        let mut last = first.force_n;
        for i in 1..=n {
            let (next, report) = state
                .advance_with_interface_temperature_controlled(
                    0.1 + duration * f64::from(i) / f64::from(n),
                    Control::ExtensionM(extension),
                    &temperatures,
                    &hot,
                    "T",
                    config(),
                    || false,
                )
                .unwrap();
            assert_eq!(report.interface_parameter_work_j, 0.0);
            assert_eq!(report.interface_stored_change_j, 0.0);
            assert!(report.mechanical.force_n < last);
            temperature_balance(&state, &next, &report);
            last = report.mechanical.force_n;
            state = next;
        }
        errors.push((last - expected).abs());
    }
    assert!(
        errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8,
        "{errors:?}"
    );
    near(errors[2], 0.0, 0.0, 3e-5 * expected.abs());
}

#[test]
fn interface_temperature_refusal_or_late_cancellation_cannot_publish_partial_history() {
    let a = material("A", 1e7, true, HistoryIntegration::ExactStrainRamp);
    let b = material("B", 2e7, true, HistoryIntegration::BackwardEuler);
    let card = bond_card(&a, &b, [[E_INF, 0.5e6], [E_BRANCH, 0.7e6], [TAU, 0.2]]);
    let query_cold = query(280.0).with("humidity", 0.3).unwrap();
    let cold = resolve_bond(&card, &query_cold);
    let hot = cold.at_temperature(&card, "T", 520.0).unwrap();
    let changed_condition = resolve_bond(&card, &query(520.0).with("humidity", 0.8).unwrap());
    let foreign = bond(&a, &b);
    let temperatures = BTreeMap::from([(1, 500.0), (2, 500.0)]);
    let (loaded, _) = initial(&a, &b, &cold)
        .advance_controlled(
            0.1,
            Control::ForceN(0.1),
            &temperatures,
            &query_cold,
            config(),
            || false,
        )
        .unwrap();
    let snapshot = format!("{loaded:?}");
    for (source, axis, target, force) in [
        (&foreign, "T", 0.2, 0.1),
        (&changed_condition, "T", 0.2, 0.1),
        (&hot, "humidity", 0.2, 0.1),
        (&hot, "T", 0.1, 0.1),
        (&hot, "T", 0.2, 1e8),
    ] {
        assert!(
            loaded
                .advance_with_interface_temperature_controlled(
                    target,
                    Control::ForceN(force),
                    &temperatures,
                    source,
                    axis,
                    config(),
                    || false,
                )
                .is_err()
        );
        assert_eq!(format!("{loaded:?}"), snapshot);
    }
    let mut polls = 0;
    let (expected, expected_report) = loaded
        .advance_with_interface_temperature_controlled(
            0.2,
            Control::ForceN(0.1),
            &temperatures,
            &hot,
            "T",
            config(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut calls = 0;
    assert!(
        loaded
            .advance_with_interface_temperature_controlled(
                0.2,
                Control::ForceN(0.1),
                &temperatures,
                &hot,
                "T",
                config(),
                || {
                    calls += 1;
                    calls == polls
                },
            )
            .is_err()
    );
    assert_eq!(format!("{loaded:?}"), snapshot);
    let (retry, retry_report) = loaded
        .advance_with_interface_temperature_controlled(
            0.2,
            Control::ForceN(0.1),
            &temperatures,
            &hot,
            "T",
            config(),
            || false,
        )
        .unwrap();
    assert_eq!(format!("{retry:?}"), format!("{expected:?}"));
    assert_eq!(format!("{retry_report:?}"), format!("{expected_report:?}"));
    temperature_balance(&loaded, &retry, &retry_report);
}
