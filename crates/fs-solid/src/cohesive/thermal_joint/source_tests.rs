//! Manufactured source data exercise the original coupled joint equations.
use super::source::{SourceJoint, SourceJointError, SourceMember};
use super::{Budget, Config};
use std::collections::BTreeMap;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::cohesive::source::*;
use fs_material::state_point::{MaterialPropertySelection, interface_stiffness_per_area_dims};
use fs_material::thermovisco::Limits;
use fs_material::thermovisco::source::{ResolvedThermoMaxwell, ThermoMaxwellSourceQuery};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature, Time};
use fs_solver::NewtonKrylovConfig;

const ONE: Dims = Dims([0; 6]);
fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new().with("T", temperature).unwrap()
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
        validity: ValidityDomain::unconstrained().with("T", 280.0, 520.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured thermal joint source; no measured specimen claim".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn material(name: &str, modulus: f64, branch: bool) -> ResolvedThermoMaxwell {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        ("E", Pressure::DIMS, [modulus; 2]),
        ("branch", Pressure::DIMS, [2.0 * modulus; 2]),
        ("tau", Time::DIMS, [2.0; 2]),
        ("scale", ONE, [1.0; 2]),
        ("shift", ONE, [1.0; 2]),
        ("free", ONE, [-0.002, 0.022]),
    ] {
        claims.insert_claim(claim(name, dims, values)).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "manufactured thermal joint fixture".into(),
            revision: 1,
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
            modulus_scale: &PropertyKey::new("scale", ONE),
            relaxation_scale: &PropertyKey::new("shift", ONE),
            free_strain: &PropertyKey::new("free", ONE),
            temperature_axis: "T",
            lower: &point(280.0),
            upper: &point(520.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        Limits {
            max_strain: 0.04,
            max_strain_rate_s: 1.0,
        },
        || false,
    )
    .unwrap()
}
fn interface(
    a: &ResolvedThermoMaxwell,
    b: &ResolvedThermoMaxwell,
    strength: f64,
    fracture_energy: f64,
) -> ResolvedBilinearCohesive {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (
            COHESIVE_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            1e10,
        ),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, strength),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            fracture_energy,
        ),
        (
            COHESIVE_COMPRESSION_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            2e10,
        ),
    ] {
        claims.insert_claim(claim(name, dims, [value; 2])).unwrap();
    }
    let surface = |source: &ResolvedThermoMaxwell| SurfaceSpec {
        material: source.samples()[0].material().clone(),
        texture_frame: "fixed axial reference frame".into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface(a),
        surface(b),
        SystemContext {
            medium: "explicit cohesive bond".into(),
            third_body: None,
            environment: "interface externally fixed at 300 K".into(),
            history: "initially bonded".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_bilinear_cohesive(
        &card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        0.01,
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
        .map(|(index, material)| SourceMember {
            id: index as u64 + 1,
            reference_length_m: 0.005,
            reference_area_m2: 1e-4,
            material,
            birth_temperature_k: 500.0,
        })
        .collect()
}
fn budget() -> Budget {
    Budget {
        max_members: 8,
        max_branch_states: 32,
        max_krylov_basis_entries: 64,
    }
}
fn config() -> Config {
    Config {
        equilibrium: super::super::Config {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-10,
                relative_tolerance: 1e-12,
                linear_restart: 4,
                forcing_maximum: 0.01,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 80,
            max_operator_evaluations: 2000,
            force_tolerance_n: 1e-7,
            local_energy_tolerance_j_m2: 1e-7,
            total_cohesive_energy_tolerance_j: 1e-9,
            maximum_displacement_m: 0.005,
        },
        extension_tolerance_m: 1e-11,
        energy_tolerance_j: 1e-7,
    }
}
fn joint<'a>(
    a: &'a ResolvedThermoMaxwell,
    b: &'a ResolvedThermoMaxwell,
    bond: &'a ResolvedBilinearCohesive,
) -> SourceJoint<'a> {
    SourceJoint::try_new(0.0, members(a, b), 7, 1e-4, bond, [1, 2], budget(), || {
        false
    })
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() < tolerance, "{a:e} != {b:e}");
}

#[test]
fn original_bulk_and_interface_sources_match_independent_series_and_thermal_work() {
    let a = material("A", 1e9, false);
    let b = material("B", 2e9, false);
    let bond = interface(&a, &b, 2e5, 10.0);
    let state = joint(&a, &b, &bond);
    let temperatures = BTreeMap::from([(1, 480.0), (2, 490.0)]);
    let (next, report) = state
        .advance_controlled(1.0, 0.0, &temperatures, &point(300.0), config(), || false)
        .unwrap();
    let compliance = 0.005 / (1e-4 * 1e9) + 0.005 / (1e-4 * 2e9);
    let bond_stiffness = 1e-4 * 1e10;
    let free_contraction = 0.005 * 1e-4 * (20.0 + 10.0);
    let force = free_contraction / (compliance + 1.0 / bond_stiffness);
    near(report.force_n, force, 1e-7);
    near(report.opening_m, force / bond_stiffness, 1e-13);
    near(
        next.joint().stored_energy_j(),
        0.5 * force * free_contraction,
        1e-11,
    );
    near(
        report.bulk_energy_j.thermal_strain_work,
        force * free_contraction,
        1e-11,
    );
    near(
        report.bulk_energy_j.numerical_dissipation,
        0.5 * force * force * compliance,
        1e-11,
    );
    near(report.endpoint_work_j, 0.0, 1e-15);
    near(report.bulk_energy_j.dissipated, 0.0, 1e-15);
    near(report.energy_residual_j, 0.0, config().energy_tolerance_j);
    assert!(core::ptr::eq(next.member_sources()[&1], &a));
    assert!(core::ptr::eq(next.member_sources()[&2], &b));
    assert!(core::ptr::eq(next.interface_source(), &bond));
    assert_eq!(next.surface_members(), [1, 2]);
}

#[test]
fn source_only_substitutions_change_response_and_keep_relaxation_and_damage_history() {
    let mut forces = Vec::new();
    for modulus in [1e9, 2e9] {
        let a = material("A", modulus, true);
        let b = material("B", modulus, true);
        let bond = interface(&a, &b, 2e5, 10.0);
        let state = joint(&a, &b, &bond);
        let temperatures = BTreeMap::from([(1, 500.0), (2, 500.0)]);
        let (loaded, report) = state
            .advance_controlled(
                0.1,
                1e-5,
                &temperatures,
                state.interface_query(),
                config(),
                || false,
            )
            .unwrap();
        forces.push(report.force_n);
        let q0 = loaded.joint().members()[&1].point().branch_strains()[0];
        let (relaxed, report) = loaded
            .advance_controlled(
                1.1,
                1e-5,
                &temperatures,
                loaded.interface_query(),
                config(),
                || false,
            )
            .unwrap();
        let gamma = 1e-4 / (1e6 * 0.01);
        let effective_tau = 2.0 * (1.0 + gamma * 3.0 * modulus) / (1.0 + gamma * modulus);
        near(
            relaxed.joint().members()[&1].point().branch_strains()[0],
            q0 / (1.0 + 1.0 / effective_tau),
            1e-13,
        );
        assert!(report.bulk_energy_j.dissipated > 0.0);
        assert!(report.force_n < forces[forces.len() - 1]);
        assert!(core::ptr::eq(relaxed.member_sources()[&1], &a));
    }
    assert!(forces[1] > forces[0]);

    let a = material("A", 1e9, false);
    let b = material("B", 1e9, false);
    let weak = interface(&a, &b, 2e5, 10.0);
    let strong = interface(&a, &b, 1e6, 2000.0);
    let temperatures = BTreeMap::from([(1, 300.0), (2, 300.0)]);
    let (broken, failure) = joint(&a, &b, &weak)
        .advance_controlled(1.0, 0.0, &temperatures, &point(300.0), config(), || false)
        .unwrap();
    let (surviving, _) = joint(&a, &b, &strong)
        .advance_controlled(1.0, 0.0, &temperatures, &point(300.0), config(), || false)
        .unwrap();
    assert!(broken.joint().bond().failed());
    assert!(!surviving.joint().bond().failed());
    near(failure.cohesive_energy_j.dissipated, 1e-4 * 10.0, 1e-10);
    let (reheated, report) = broken
        .advance_controlled(
            2.0,
            0.0,
            &BTreeMap::from([(1, 500.0), (2, 500.0)]),
            &point(300.0),
            config(),
            || false,
        )
        .unwrap();
    assert!(reheated.joint().bond().failed());
    near(report.cohesive_energy_j.dissipated, 0.0, 1e-15);
    assert!(core::ptr::eq(reheated.interface_source(), &weak));
}

#[test]
fn source_assignment_interface_state_domains_and_cancellation_refuse_atomically() {
    let a = material("A", 1e9, true);
    let b = material("B", 1e9, true);
    let bond = interface(&a, &b, 2e5, 10.0);
    assert!(matches!(
        SourceJoint::try_new(
            0.0,
            members(&a, &b),
            7,
            1e-4,
            &bond,
            [2, 1],
            budget(),
            || false
        ),
        Err(SourceJointError::Interface(_))
    ));
    assert!(
        SourceJoint::try_new(
            0.0,
            members(&a, &b),
            7,
            1e-4,
            &bond,
            [1, 99],
            budget(),
            || false
        )
        .is_err()
    );
    let mut invalid = members(&a, &b);
    invalid[0].birth_temperature_k = 600.0;
    assert!(
        SourceJoint::try_new(0.0, invalid, 7, 1e-4, &bond, [1, 2], budget(), || false).is_err()
    );
    assert!(matches!(
        SourceJoint::try_new(
            0.0,
            members(&a, &b),
            7,
            1e-4,
            &bond,
            [1, 2],
            budget(),
            || true
        ),
        Err(SourceJointError::Cancelled)
    ));
    let state = joint(&a, &b, &bond);
    let saved = state.joint().clone();
    let temperatures = BTreeMap::from([(1, 480.0), (2, 480.0)]);
    let typed = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap();
    for changed in [point(301.0), typed] {
        assert!(matches!(
            state.advance_controlled(1.0, 0.0, &temperatures, &changed, config(), || false),
            Err(SourceJointError::Invalid(_))
        ));
    }
    for outside in [
        BTreeMap::from([(1, 600.0), (2, 480.0)]),
        BTreeMap::from([(1, 480.0)]),
    ] {
        assert!(
            state
                .advance_controlled(
                    1.0,
                    0.0,
                    &outside,
                    state.interface_query(),
                    config(),
                    || false
                )
                .is_err()
        );
    }
    for allowed in [0, 2, 10] {
        let mut calls = 0;
        assert!(matches!(
            state.advance_controlled(
                1.0,
                0.0,
                &temperatures,
                state.interface_query(),
                config(),
                || {
                    calls += 1;
                    calls > allowed
                }
            ),
            Err(SourceJointError::Cancelled)
        ));
        assert_eq!(state.joint(), &saved);
    }
    let (a, _) = state
        .advance_controlled(
            1.0,
            0.0,
            &temperatures,
            state.interface_query(),
            config(),
            || false,
        )
        .unwrap();
    let (b, _) = state
        .clone()
        .advance_controlled(
            1.0,
            0.0,
            &temperatures,
            state.interface_query(),
            config(),
            || false,
        )
        .unwrap();
    assert_eq!(a.joint(), b.joint());
}
