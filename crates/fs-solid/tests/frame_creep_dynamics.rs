//! Manufactured momentum/work references for the native timed frame fibers.
//! These coefficients verify algorithms, not measured material qualification.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::creep::plastic::NortonJ2;
use fs_material::creep::plastic::source::NortonJ2Source;
use fs_material::creep::plastic::uniaxial::UniaxialNortonJ2;
use fs_material::creep::source::{NortonCreepSource, NortonSourceQuery, norton_requirements};
use fs_material::creep::thermal::expansion::ExpandingCreepSource;
use fs_material::creep::uniaxial::UniaxialNorton;
use fs_material::creep::{
    CreepTolerances, NortonCreep, NortonCreepParameters, StrainHardeningCreep,
};
use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::source::{HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY};
use fs_material::state_point::{
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, YIELD_STRESS_PROPERTY, resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec};
use fs_solid::fiber::{Fiber, FiberLaw, Section};
use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::frame2d::creep::{FrameCreepError, FrameCreepTolerances};
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel, SolidError};
use std::cell::Cell;

const E: f64 = 1e9;
const STRESS: f64 = 1e6;
const RATE: f64 = 1e-2;
const DT: f64 = 0.1;
const AREA: f64 = 5e-4;
const ARM: f64 = 0.05;
const MASS: f64 = 1e4;
const FORCE: f64 = 2.0 * AREA * STRESS;
const YIELD: f64 = 0.25e6;
const HARDENING: f64 = 2e8;
const ALPHA: f64 = 2e-5;

fn kinetics(duration_s: f64) -> NortonCreep {
    NortonCreep::new(
        IsotropicElastic::new(E, 0.25, 0.03).unwrap(),
        NortonCreepParameters {
            reference_stress_pa: STRESS,
            reference_rate_s: RATE,
            stress_exponent: 1.0,
            activation_energy_j_mol: 0.0,
            reference_temperature_k: 300.0,
            temperature_range_k: [300.0, 500.0],
        },
        300.0,
        duration_s,
    )
    .unwrap()
}

fn law(duration_s: f64) -> FiberLaw {
    FiberLaw::Norton(UniaxialNorton::new(kinetics(duration_s)))
}

fn frame(law: &FiberLaw, formulation: FrameMassFormulation, axial: bool) -> Frame2d {
    let mut result = Frame2d::with_mass_formulation(Geometry::Linear, formulation);
    result.add_node(0.0, 0.0);
    result.add_node(1.0, 0.0);
    result.fix(0, [true; 3]);
    if axial {
        result.fix(1, [false, true, true]);
    }
    result
        .add_element(0, 1, 3, MASS, &|| {
            SectionModel::Fiber(Section::new(
                [-ARM, ARM]
                    .into_iter()
                    .map(|y| Fiber {
                        y,
                        area: AREA,
                        law: law.clone(),
                    })
                    .collect(),
            ))
        })
        .unwrap();
    result
}

fn controls() -> FrameCreepTolerances {
    FrameCreepTolerances {
        material: CreepTolerances {
            stress_pa: 1e-3,
            energy_j_m3: 1e-5,
        },
        energy_j: 1e-7,
        force_n: 1e-5,
        moment_n_m: 1e-6,
        maximum_fiber_points: 6,
    }
}

fn axial_load(force: f64) -> [f64; 6] {
    [0.0, 0.0, 0.0, force, 0.0, 0.0]
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:.3e}",
    );
}

#[test]
fn axial_motion_matches_closed_form_on_loading_and_release_without_extra_aging() {
    for (formulation, mass) in [
        (FrameMassFormulation::Lumped, MASS / 2.0),
        (FrameMassFormulation::ConsistentReference, MASS / 3.0),
    ] {
        let mut accepted = frame(&law(DT), formulation, true);
        let compliance = 1.0 / E + DT * RATE / STRESS;
        let stiffness = 2.0 * AREA / compliance;
        let (mut u, mut v, mut creep, mut accumulated, mut stress) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for force in [FORCE, FORCE, FORCE, FORCE, 0.0, 0.0, 0.0, 0.0] {
            let old_state = accepted.creep_fiber(0, 0, 0).unwrap().1.clone();
            let old_u = accepted.displacements().to_vec();
            let old_v = accepted.velocities().to_vec();
            // The n=1 material compliance and scalar momentum equation are
            // linear. This is independent of frame/section Newton iterations.
            let acceleration =
                (force - stiffness * (u + DT * v - creep)) / (mass + DT * DT * stiffness);
            let next_v = v + DT * acceleration;
            let next_u = u + DT * next_v;
            let next_stress = (next_u - creep) / compliance;
            let creep_increment = DT * RATE * next_stress / STRESS;
            let step = accepted
                .advance_creep_backward_euler(&axial_load(force), DT, controls(), || false)
                .unwrap();
            assert_eq!(accepted.displacements(), old_u);
            assert_eq!(accepted.velocities(), old_v);
            assert_eq!(accepted.creep_fiber(0, 0, 0).unwrap().1, &old_state);
            close(step.frame.displacements()[3], next_u, 2e-11);
            close(step.frame.velocities()[3], next_v, 2e-10);
            close(step.frame.accelerations()[3], acceleration, 2e-9);
            close(
                step.frame.element_forces(0)[0],
                2.0 * AREA * next_stress,
                1e-5,
            );
            close(
                step.energy.material.creep_dissipation_j,
                2.0 * AREA * next_stress * creep_increment,
                1e-7,
            );
            close(
                step.energy.material.numerical_dissipation_j,
                AREA * (next_stress - stress).powi(2) / E,
                1e-7,
            );
            close(step.energy.initial_kinetic_j, 0.5 * mass * v * v, 1e-7);
            close(
                step.energy.final_kinetic_j,
                0.5 * mass * next_v * next_v,
                1e-7,
            );
            close(
                step.energy.inertial_numerical_dissipation_j,
                0.5 * mass * (next_v - v).powi(2),
                1e-7,
            );
            close(
                step.energy.caloric_heat_j(),
                step.energy.material.creep_dissipation_j,
                1e-12,
            );
            close(
                step.energy.inertial_work_residual_j,
                0.0,
                controls().energy_j,
            );
            close(
                step.energy.material.equilibrium_work_residual_j,
                0.0,
                controls().energy_j,
            );
            close(step.energy.balance_residual_j, 0.0, controls().energy_j);
            assert_eq!(step.duration_s, DT);
            assert_eq!(step.fibers.len(), 6);
            creep += creep_increment;
            accumulated += creep_increment.abs();
            for section in 0..3 {
                for fiber in 0..2 {
                    let state = step.frame.creep_fiber(0, section, fiber).unwrap().1;
                    close(state.creep_strain(), creep, 2e-11);
                    close(state.equivalent_creep_strain(), accumulated, 2e-11);
                }
            }
            (u, v, stress) = (next_u, next_v, next_stress);
            accepted = step.frame;
        }
        // The static route must not silently discard already accepted motion.
        assert!(matches!(
            accepted.solve_creep_step(&axial_load(0.0), controls(), || false),
            Err(FrameCreepError::Invalid(_))
        ));
    }
}

#[test]
fn bending_uses_original_translation_rotation_inertia_and_independent_moment_gate() {
    let effective_e = 1.0 / (1.0 / E + DT * RATE / STRESS);
    let ei = effective_e * 2.0 * AREA * ARM * ARM;
    for (formulation, inertia) in [
        (FrameMassFormulation::Lumped, [MASS / 2.0, 0.0, 0.0]),
        (
            FrameMassFormulation::ConsistentReference,
            [
                156.0 * MASS / 420.0,
                -22.0 * MASS / 420.0,
                4.0 * MASS / 420.0,
            ],
        ),
    ] {
        let accepted = frame(&law(DT), formulation, false);
        let a = inertia[0] + DT * DT * 12.0 * ei;
        let b = inertia[1] - DT * DT * 6.0 * ei;
        let c = inertia[2] + DT * DT * 4.0 * ei;
        let force = 10.0;
        let moment = 0.2;
        let det = a * c - b * b;
        let acceleration = [
            (c * force - b * moment) / det,
            (a * moment - b * force) / det,
        ];
        let step = accepted
            .advance_creep_backward_euler(
                &[0.0, 0.0, 0.0, 0.0, force, moment],
                DT,
                controls(),
                || false,
            )
            .unwrap();
        for (dof, expected) in [(4, acceleration[0]), (5, acceleration[1])] {
            close(step.frame.accelerations()[dof], expected, 2e-9);
            close(step.frame.velocities()[dof], DT * expected, 2e-10);
            close(step.frame.displacements()[dof], DT * DT * expected, 2e-11);
        }
        let velocity = [DT * acceleration[0], DT * acceleration[1]];
        let kinetic = 0.5
            * (inertia[0] * velocity[0].powi(2)
                + 2.0 * inertia[1] * velocity[0] * velocity[1]
                + inertia[2] * velocity[1].powi(2));
        close(step.energy.final_kinetic_j, kinetic, 1e-9);
        close(step.energy.inertial_numerical_dissipation_j, kinetic, 1e-9);
        assert!(step.energy.material.creep_dissipation_j > 0.0);
        close(step.energy.balance_residual_j, 0.0, controls().energy_j);
    }
}

#[test]
fn primary_creep_and_joint_plastic_histories_follow_the_same_momentum_interval() {
    let reference_strain = 4e-4;
    let creep_law = kinetics(DT)
        .with_strain_hardening(StrainHardeningCreep {
            reference_strain,
            exponent: 1.0,
            residual_rate_fraction: 0.0,
        })
        .unwrap();
    let plastic = J2Plasticity::new(creep_law.elastic().clone(), YIELD, HARDENING).unwrap();
    let law = FiberLaw::NortonJ2(UniaxialNortonJ2::new(
        NortonJ2::new(creep_law, plastic).unwrap(),
    ));
    let mut accepted = frame(&law, FrameMassFormulation::Lumped, true);
    let mass = MASS / 2.0;
    let (mut u, mut v, mut creep, mut plastic) = (0.0, 0.0, 0.0, 0.0);
    for _ in 0..4 {
        // For primary exponent one, delta_c has an independent quadratic
        // root at known tensile stress; linear hardening determines delta_p.
        let increments = |stress: f64| {
            let b = reference_strain + creep;
            let c = DT * RATE * stress / STRESS * reference_strain;
            (
                2.0 * c / (b + (b * b + 4.0 * c).sqrt()),
                ((stress - YIELD) / HARDENING - plastic).max(0.0),
            )
        };
        let residual = |stress: f64| {
            let (dc, dp) = increments(stress);
            let strain = stress / E + creep + dc + plastic + dp;
            2.0 * AREA * stress + mass * (strain - u - DT * v) / (DT * DT) - FORCE
        };
        let (mut lo, mut hi) = (0.0, 2.0 * STRESS);
        assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if residual(mid) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let stress = 0.5 * (lo + hi);
        let (dc, dp) = increments(stress);
        let next_u = stress / E + creep + dc + plastic + dp;
        let next_v = (next_u - u) / DT;
        let old_history = accepted.creep_plastic_fiber(0, 0, 0).unwrap().1.clone();
        let step = accepted
            .advance_creep_backward_euler(&axial_load(FORCE), DT, controls(), || false)
            .unwrap();
        assert_eq!(
            accepted.creep_plastic_fiber(0, 0, 0).unwrap().1,
            &old_history
        );
        close(step.frame.displacements()[3], next_u, 2e-11);
        close(step.frame.velocities()[3], next_v, 2e-10);
        close(
            step.energy.material.creep_dissipation_j,
            2.0 * AREA * stress * dc,
            1e-7,
        );
        close(
            step.energy.material.plastic_dissipation_j,
            2.0 * AREA * YIELD * dp,
            1e-7,
        );
        creep += dc;
        plastic += dp;
        for section in 0..3 {
            for fiber in 0..2 {
                let state = step.frame.creep_plastic_fiber(0, section, fiber).unwrap().1;
                close(state.creep_strain(), creep, 2e-11);
                close(state.equivalent_creep_strain(), creep, 2e-11);
                close(state.plastic_strain(), plastic, 2e-11);
                close(state.equivalent_plastic_strain(), plastic, 2e-11);
            }
        }
        close(step.energy.balance_residual_j, 0.0, controls().energy_j);
        (u, v) = (next_u, next_v);
        accepted = step.frame;
    }
    assert!(creep > 0.0 && plastic > 0.0);
}

fn temperature_spec() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_spec(), temperature)
        .unwrap()
}

fn claim(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with_quantity(
            "T",
            temperature_spec(),
            300.0,
            500.0,
        ),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured inertial thermocreep".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

#[test]
fn both_thermal_families_keep_reciprocal_heat_outside_numerical_inertia_loss() {
    let mut claims = ClaimSet::new();
    for (requirement, value) in norton_requirements()
        .unwrap()
        .iter()
        .zip([E, 0.25, STRESS, RATE, 1.0, 0.0, 300.0])
    {
        claims
            .insert_claim(claim(requirement.name(), requirement.dims(), value))
            .unwrap();
    }
    let expansion_dims = Dims([0, 0, 0, -1, 0, 0]);
    for (name, dims, value) in [
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            expansion_dims,
            ALPHA,
        ),
        (
            ISOTROPIC_HARDENING_MODULUS_PROPERTY,
            Pressure::DIMS,
            HARDENING,
        ),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, YIELD),
    ] {
        claims.insert_claim(claim(name, dims, value)).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic inertial thermocreep".into(),
            phase: "solid".into(),
            process: "manufactured".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let expansion = resolve_material_state_point(
        &card,
        &point(300.0),
        &[ScalarPropertyRequirement::try_new(
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            expansion_dims,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let lower = point(300.0);
    let upper = point(500.0);
    let query = || NortonSourceQuery {
        lower: &lower,
        upper: &upper,
        temperature_axis: "T",
        selection: MaterialPropertySelection::SingleClaimOnly,
        strain_limit: 0.03,
    };
    let plain = ExpandingCreepSource::new(
        &card,
        NortonCreepSource::new(&card, query(), || false).unwrap(),
        &expansion,
        || false,
    )
    .unwrap();
    let combined = ExpandingCreepSource::new(
        &card,
        NortonJ2Source::new(
            &card,
            query(),
            HardeningPropertyBasis::PlasticStrainModulus,
            || false,
        )
        .unwrap(),
        &expansion,
        || false,
    )
    .unwrap();
    for (cold, hot, plastic) in [
        (
            FiberLaw::ThermalNorton(plain.at_uniaxial_temperature(300.0, DT).unwrap()),
            FiberLaw::ThermalNorton(plain.at_uniaxial_temperature(350.0, DT).unwrap()),
            false,
        ),
        (
            FiberLaw::ThermalNortonJ2(combined.at_uniaxial_temperature(300.0, DT).unwrap()),
            FiberLaw::ThermalNortonJ2(combined.at_uniaxial_temperature(350.0, DT).unwrap()),
            true,
        ),
    ] {
        let accepted = frame(&cold, FrameMassFormulation::Lumped, true);
        let prepared = accepted
            .with_creep_laws(&vec![hot; 6], "T", || false)
            .unwrap();
        let step = prepared
            .advance_creep_backward_euler(&axial_load(1.5 * FORCE), DT, controls(), || false)
            .unwrap();
        let mass = MASS / 2.0;
        let compliance = 1.0 / E + DT * RATE / STRESS;
        let plastic_compliance = if plastic { 1.0 / HARDENING } else { 0.0 };
        let free_strain = ALPHA * 50.0;
        let stress = (DT * DT * 1.5 * FORCE / mass - free_strain + YIELD * plastic_compliance)
            / (compliance + DT * DT * 2.0 * AREA / mass + plastic_compliance);
        let displacement =
            free_strain + compliance * stress + (stress - YIELD) * plastic_compliance;
        close(step.frame.displacements()[3], displacement, 2e-11);
        close(step.frame.element_forces(0)[0], 2.0 * AREA * stress, 1e-5);
        assert!(step.energy.material.thermoelastic_storage_change_j > 0.0);
        assert!(step.energy.material.reversible_heat_j < 0.0);
        assert!(step.energy.inertial_numerical_dissipation_j > 0.0);
        close(
            step.energy.caloric_heat_j(),
            step.energy.material.physical_dissipation_j + step.energy.material.reversible_heat_j,
            0.0,
        );
        if plastic {
            assert!(step.energy.material.plastic_dissipation_j > 0.0);
            assert_eq!(
                step.frame
                    .thermal_creep_plastic_fiber(0, 0, 0)
                    .unwrap()
                    .1
                    .temperature_k(),
                350.0
            );
        } else {
            assert_eq!(
                step.frame
                    .thermal_creep_fiber(0, 0, 0)
                    .unwrap()
                    .1
                    .temperature_k(),
                350.0
            );
        }
        close(step.energy.balance_residual_j, 0.0, controls().energy_j);
    }
}

#[test]
fn interval_controls_and_cancellation_refuse_without_partial_motion_or_history() {
    let accepted = frame(&law(DT), FrameMassFormulation::ConsistentReference, true);
    let before = accepted.creep_fiber(0, 0, 0).unwrap().1.clone();
    for duration in [0.0, DT / 2.0, f64::NAN, f64::INFINITY, f64::MIN_POSITIVE] {
        assert!(matches!(
            accepted
                .advance_creep_backward_euler(&axial_load(FORCE), duration, controls(), || false),
            Err(FrameCreepError::Invalid(_))
        ));
    }
    let mut invalid_controls = controls();
    invalid_controls.maximum_fiber_points = 5;
    assert!(matches!(
        accepted.advance_creep_backward_euler(&axial_load(FORCE), DT, invalid_controls, || false),
        Err(FrameCreepError::Invalid(_))
    ));
    invalid_controls = controls();
    invalid_controls.energy_j = f64::NAN;
    assert!(matches!(
        accepted.advance_creep_backward_euler(&axial_load(FORCE), DT, invalid_controls, || false),
        Err(FrameCreepError::Invalid(_))
    ));
    let mut unequal = vec![law(DT); 6];
    unequal[0] = law(DT / 2.0);
    let prepared = accepted.with_creep_laws(&unequal, "T", || false).unwrap();
    assert!(matches!(
        prepared.advance_creep_backward_euler(&axial_load(FORCE), DT, controls(), || false),
        Err(FrameCreepError::Invalid(_))
    ));
    for at in [1, 9, 20] {
        let polls = Cell::new(0);
        let result =
            accepted.advance_creep_backward_euler(&axial_load(FORCE), DT, controls(), || {
                polls.set(polls.get() + 1);
                polls.get() >= at
            });
        assert!(matches!(result, Err(FrameCreepError::Cancelled)));
    }
    assert_eq!(accepted.displacements(), &[0.0; 6]);
    assert_eq!(accepted.velocities(), &[0.0; 6]);
    assert_eq!(accepted.creep_fiber(0, 0, 0).unwrap().1, &before);
    let mut untimed = accepted.clone();
    assert!(matches!(
        untimed.newmark(&[0.0, 0.0], DT, (1.0, 0.0), Rayleigh { a0: 0.0, a1: 0.0 }),
        Err(SolidError::MaterialRefused { .. })
    ));
    assert!(
        accepted
            .advance_creep_backward_euler(&axial_load(FORCE), DT, controls(), || false)
            .is_ok()
    );
}
