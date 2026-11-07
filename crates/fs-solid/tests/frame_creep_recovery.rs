//! Actual native members retain permanent creep while resistance recovers.
//! Manufactured coefficients verify history/clock/energy behavior, not service life.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::IsotropicElastic;
use fs_material::creep::recovery::source::{
    RecoveringNortonSource, recovering_norton_requirements,
};
use fs_material::creep::recovery::uniaxial::UniaxialRecoveringNorton;
use fs_material::creep::recovery::{CreepRecoveryParameters, RecoveringNortonCreep};
use fs_material::creep::source::NortonSourceQuery;
use fs_material::creep::uniaxial::UniaxialNorton;
use fs_material::creep::{
    CreepTolerances, NortonCreep, NortonCreepParameters, StrainHardeningCreep,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_solid::fiber::material::FiberMaterial;
use fs_solid::fiber::{Fiber, FiberLaw, Section};
use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::frame2d::creep::{FrameCreepError, FrameCreepTolerances};
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel, SolidError};
use std::cell::Cell;

const E: f64 = 1e9;
const STRESS: f64 = 1e6;
const RATE: f64 = 1e-5;
const DT: f64 = 10.0;
const PRIMARY: f64 = 4e-4;
const RECOVERY: f64 = 0.2;
const AREA: f64 = 5e-4;
const ARM: f64 = 0.05;
const FORCE: f64 = 2.0 * AREA * STRESS;

fn law(dt: f64, recovery: f64) -> UniaxialRecoveringNorton {
    let native = NortonCreep::new(
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
        dt,
    )
    .unwrap()
    .with_strain_hardening(StrainHardeningCreep {
        reference_strain: PRIMARY,
        exponent: 1.0,
        residual_rate_fraction: 0.0,
    })
    .unwrap();
    UniaxialRecoveringNorton::new(
        RecoveringNortonCreep::new(
            native,
            CreepRecoveryParameters {
                reference_rate_s: recovery,
                activation_energy_j_mol: 0.0,
                reference_temperature_k: 300.0,
            },
        )
        .unwrap(),
    )
}

fn frame(law: &UniaxialRecoveringNorton, mass: f64, formulation: FrameMassFormulation) -> Frame2d {
    // Exercise the complete fiber material dispatch while preserving authored
    // density status even when the law itself has original source receipts.
    let material =
        FiberMaterial::authored(FiberLaw::RecoveringNorton(law.clone()), mass / (2.0 * AREA))
            .unwrap();
    assert!(material.source().is_none());
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, formulation);
    frame.add_node(0.0, 0.0);
    frame.add_node(1.0, 0.0);
    frame.fix(0, [true; 3]);
    frame.fix(1, [false, true, true]);
    frame
        .add_element(0, 1, 3, mass, &|| {
            SectionModel::Fiber(Section::new(
                [-ARM, ARM]
                    .into_iter()
                    .map(|y| Fiber {
                        y,
                        area: AREA,
                        law: material.law().clone(),
                    })
                    .collect(),
            ))
        })
        .unwrap();
    frame
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

fn load(force: f64) -> [f64; 6] {
    [0.0, 0.0, 0.0, force, 0.0, 0.0]
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:.3e}"
    );
}

// Independent primary-exponent-one quadratic at a known endpoint stress.
fn delta(stress: f64, resistance: f64, dt: f64) -> f64 {
    let a = PRIMARY * (1.0 + RECOVERY * dt);
    let b = a + resistance;
    let c = dt * RATE * stress.abs() / STRESS * a;
    2.0 * c / (b + (b * b + 4.0 * c).sqrt())
}

#[test]
fn loaded_member_relaxes_when_held_and_dwell_renews_reload_without_erasing_permanent_set() {
    let mut accepted = frame(&law(DT, RECOVERY), 1.0, FrameMassFormulation::Lumped);
    let (mut creep, mut resistance) = (0.0, 0.0);
    for _ in 0..3 {
        let before = accepted.creep_recovery_fiber(0, 0, 0).unwrap().1.clone();
        let increment = delta(STRESS, resistance, DT);
        let step = accepted
            .solve_creep_step(&load(FORCE), controls(), || false)
            .unwrap();
        assert_eq!(accepted.creep_recovery_fiber(0, 0, 0).unwrap().1, &before);
        creep += increment;
        resistance = (resistance + increment) / (1.0 + RECOVERY * DT);
        close(step.frame.displacements()[3], STRESS / E + creep, 2e-11);
        close(step.frame.element_forces(0)[0], FORCE, 1e-6);
        close(step.energy.creep_dissipation_j, FORCE * increment, 1e-8);
        close(
            step.energy.caloric_heat_j(),
            step.energy.creep_dissipation_j,
            0.0,
        );
        close(step.energy.balance_residual_j, 0.0, controls().energy_j);
        for section in 0..3 {
            for fiber in 0..2 {
                let state = step
                    .frame
                    .creep_recovery_fiber(0, section, fiber)
                    .unwrap()
                    .1;
                close(state.creep_strain(), creep, 2e-11);
                close(state.equivalent_creep_strain(), creep, 2e-11);
                close(state.resistance_strain(), resistance, 2e-12);
            }
        }
        accepted = step.frame;
    }
    let mut held = accepted.clone();
    held.fix(1, [true; 3]);
    let relaxed = held
        .solve_creep_step(&[0.0; 6], controls(), || false)
        .unwrap();
    assert_eq!(relaxed.frame.displacements(), held.displacements());
    assert!(relaxed.frame.element_forces(0)[0] < FORCE);
    assert!(relaxed.frame.element_forces(0)[0] > 0.0);
    assert!(relaxed.energy.creep_dissipation_j > 0.0);
    close(
        relaxed.energy.applied_work_j + relaxed.energy.support_work_j,
        0.0,
        0.0,
    );
    close(relaxed.energy.balance_residual_j, 0.0, controls().energy_j);

    let released = accepted
        .solve_creep_step(&load(0.0), controls(), || false)
        .unwrap();
    close(released.frame.displacements()[3], creep, 2e-11);
    close(released.energy.creep_dissipation_j, 0.0, 1e-8);
    let released_history = released
        .frame
        .creep_recovery_fiber(0, 0, 0)
        .unwrap()
        .1
        .clone();
    let immediate = released
        .frame
        .solve_creep_step(&load(FORCE), controls(), || false)
        .unwrap();
    let mut rested = released.frame;
    for _ in 0..4 {
        let before = rested.creep_recovery_fiber(0, 0, 0).unwrap().1.clone();
        let step = rested
            .solve_creep_step(&load(0.0), controls(), || false)
            .unwrap();
        let after = step.frame.creep_recovery_fiber(0, 0, 0).unwrap().1;
        close(after.creep_strain(), released_history.creep_strain(), 2e-11);
        close(
            after.equivalent_creep_strain(),
            released_history.equivalent_creep_strain(),
            2e-11,
        );
        close(
            after.resistance_strain(),
            before.resistance_strain() / 3.0,
            2e-12,
        );
        close(step.energy.creep_dissipation_j, 0.0, 1e-8);
        close(step.energy.balance_residual_j, 0.0, controls().energy_j);
        rested = step.frame;
    }
    let reloaded = rested
        .solve_creep_step(&load(FORCE), controls(), || false)
        .unwrap();
    let after = reloaded.frame.creep_recovery_fiber(0, 0, 0).unwrap().1;
    let without_dwell = immediate.frame.creep_recovery_fiber(0, 0, 0).unwrap().1;
    assert!(after.creep_strain() > without_dwell.creep_strain());
    assert!(after.equivalent_creep_strain() > without_dwell.equivalent_creep_strain());
    assert!(after.equivalent_creep_strain() > released_history.equivalent_creep_strain());
    assert!(after.resistance_strain() < after.equivalent_creep_strain());
    close(reloaded.frame.element_forces(0)[0], FORCE, 1e-6);
    close(reloaded.energy.balance_residual_j, 0.0, controls().energy_j);
}

#[test]
fn variable_interval_motion_preserves_resistance_and_matches_independent_momentum() {
    let member_mass = 1e4;
    for (formulation, mass) in [
        (FrameMassFormulation::Lumped, member_mass / 2.0),
        (FrameMassFormulation::ConsistentReference, member_mass / 3.0),
    ] {
        let mut accepted = frame(&law(0.1, RECOVERY), member_mass, formulation);
        let (mut u, mut v, mut creep, mut alpha, mut resistance) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for (dt, force) in [
            (0.1, FORCE),
            (0.1, FORCE),
            (0.05, FORCE),
            (0.1, FORCE),
            (0.05, 0.0),
            (0.1, 0.0),
        ] {
            let history = accepted.creep_recovery_fiber(0, 0, 0).unwrap().1.clone();
            let prepared = accepted
                .with_creep_laws(&vec![law(dt, RECOVERY); 6], "T", || false)
                .unwrap();
            assert_eq!(prepared.displacements(), accepted.displacements());
            assert_eq!(prepared.velocities(), accepted.velocities());
            assert_eq!(prepared.creep_recovery_fiber(0, 0, 0).unwrap().1, &history);
            let strain =
                |stress: f64| stress / E + creep + stress.signum() * delta(stress, resistance, dt);
            let residual = |stress: f64| {
                2.0 * AREA * stress + mass * (strain(stress) - u - dt * v) / (dt * dt) - force
            };
            let (mut lo, mut hi) = (-4.0 * STRESS, 4.0 * STRESS);
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
            let next_u = strain(stress);
            let acceleration = (force - 2.0 * AREA * stress) / mass;
            let next_v = v + dt * acceleration;
            let dc = delta(stress, resistance, dt);
            let step = prepared
                .advance_creep_backward_euler(&load(force), dt, controls(), || false)
                .unwrap();
            assert_eq!(prepared.creep_recovery_fiber(0, 0, 0).unwrap().1, &history);
            close(step.frame.displacements()[3], next_u, 2e-11);
            close(step.frame.velocities()[3], next_v, 2e-10);
            close(step.frame.accelerations()[3], acceleration, 2e-9);
            close(
                step.energy.material.creep_dissipation_j,
                2.0 * AREA * stress.abs() * dc,
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
                0.0,
            );
            close(step.energy.balance_residual_j, 0.0, controls().energy_j);
            creep += stress.signum() * dc;
            alpha += dc;
            resistance = (resistance + dc) / (1.0 + RECOVERY * dt);
            for section in 0..3 {
                for fiber in 0..2 {
                    let (law, state) = step.frame.creep_recovery_fiber(0, section, fiber).unwrap();
                    close(state.creep_strain(), creep, 2e-11);
                    close(state.equivalent_creep_strain(), alpha, 2e-11);
                    close(state.resistance_strain(), resistance, 2e-12);
                    assert_eq!(law.kinetics().duration_s(), dt);
                    assert_eq!(
                        law.kinetics().strain_hardening().unwrap().reference_strain,
                        PRIMARY
                    );
                }
            }
            (u, v) = (next_u, next_v);
            accepted = step.frame;
        }
        assert!(resistance > 0.0 && resistance < alpha);
    }
}

fn temperature_spec() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn pressure_spec() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}

fn point(temperature: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_spec(), temperature)
        .unwrap()
        .with_quantity("pressure", pressure_spec(), pressure)
        .unwrap()
}

fn card(name: &str) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let values = [
        E, 0.25, STRESS, RATE, 1.0, 20_000.0, 300.0, PRIMARY, 1.0, 0.0, RECOVERY, 40_000.0, 300.0,
    ];
    for (requirement, value) in recovering_norton_requirements().unwrap().iter().zip(values) {
        let dims = requirement.dims();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(requirement.name(), dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature_spec(), 300.0, 500.0)
                    .with_quantity("pressure", pressure_spec(), 1e5, 3e5),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured recovering member".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "manufactured".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn source(card: &MaterialCard, pressure: f64) -> RecoveringNortonSource<'_> {
    RecoveringNortonSource::new(
        card,
        NortonSourceQuery {
            lower: &point(300.0, pressure),
            upper: &point(500.0, pressure),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.03,
        },
        || false,
    )
    .unwrap()
}

#[test]
fn source_refresh_retains_all_recovery_claims_history_and_fixed_context() {
    let material = card("original recovering member");
    let family = source(&material, 1e5);
    let cold = family.at_uniaxial_temperature(300.0, DT).unwrap();
    let accepted = frame(&cold, 1.0, FrameMassFormulation::Lumped)
        .solve_creep_step(&load(FORCE), controls(), || false)
        .unwrap()
        .frame;
    let history = accepted.creep_recovery_fiber(0, 0, 0).unwrap().1.clone();
    let warm = family.at_uniaxial_temperature(350.0, DT / 2.0).unwrap();
    let prepared = accepted
        .with_creep_laws(&vec![warm; 6], "T", || false)
        .unwrap();
    assert_eq!(prepared.displacements(), accepted.displacements());
    assert_eq!(prepared.creep_recovery_fiber(0, 0, 0).unwrap().1, &history);
    let next = prepared
        .solve_creep_step(&load(FORCE), controls(), || false)
        .unwrap();
    let (updated, updated_history) = next.frame.creep_recovery_fiber(0, 0, 0).unwrap();
    assert_eq!(updated.material().unwrap().properties().len(), 13);
    assert_eq!(
        updated.material().unwrap().state_coordinate("T"),
        Some(350.0)
    );
    assert_eq!(
        updated.material().unwrap().card_identity(),
        cold.material().unwrap().card_identity()
    );
    for (old, next) in cold
        .material()
        .unwrap()
        .properties()
        .iter()
        .zip(updated.material().unwrap().properties())
    {
        assert_eq!(old.requirement(), next.requirement());
        assert_eq!(
            old.answer().receipt.selected,
            next.answer().receipt.selected
        );
    }
    assert_eq!(
        updated
            .kinetics()
            .strain_hardening()
            .unwrap()
            .reference_strain,
        PRIMARY
    );
    assert_eq!(
        updated.backbone().recovery_parameters(),
        cold.backbone().recovery_parameters()
    );
    assert!(updated_history.equivalent_creep_strain() > history.equivalent_creep_strain());
    close(next.energy.balance_residual_j, 0.0, controls().energy_j);
    assert!(family.at_uniaxial_temperature(501.0, DT).is_err());
    let changed_context = source(&material, 2e5)
        .at_uniaxial_temperature(350.0, DT)
        .unwrap();
    assert!(
        accepted
            .with_creep_laws(&vec![changed_context; 6], "T", || false)
            .is_err()
    );
    let foreign = card("foreign recovering member");
    let foreign_law = source(&foreign, 1e5)
        .at_uniaxial_temperature(350.0, DT)
        .unwrap();
    assert!(
        accepted
            .with_creep_laws(&vec![foreign_law; 6], "T", || false)
            .is_err()
    );
    let warmed = family.at_uniaxial_temperature(350.0, DT).unwrap();
    assert!(
        accepted
            .with_creep_laws(&vec![warmed; 6], "pressure", || false)
            .is_err()
    );
    assert_eq!(accepted.creep_recovery_fiber(0, 0, 0).unwrap().1, &history);
}

#[test]
fn changing_recovery_family_clock_or_parameters_and_cancellation_cannot_publish_history() {
    let native = law(DT, RECOVERY);
    let accepted = frame(&native, 1.0, FrameMassFormulation::Lumped);
    let virgin = accepted.creep_recovery_fiber(0, 0, 0).unwrap().1.clone();
    assert!(
        accepted
            .with_creep_laws(&vec![law(DT, 2.0 * RECOVERY); 6], "T", || false)
            .is_err()
    );
    let ordinary = UniaxialNorton::new(native.kinetics().clone());
    assert!(
        accepted
            .with_creep_laws(&vec![ordinary; 6], "T", || false)
            .is_err()
    );
    assert!(matches!(
        accepted.advance_creep_backward_euler(&load(FORCE), DT / 2.0, controls(), || false),
        Err(FrameCreepError::Invalid(_))
    ));
    let mut mixed = vec![law(DT, RECOVERY); 6];
    mixed[1] = law(DT / 2.0, RECOVERY);
    let prepared = accepted.with_creep_laws(&mixed, "T", || false).unwrap();
    assert!(matches!(
        prepared.solve_creep_step(&load(FORCE), controls(), || false),
        Err(FrameCreepError::Invalid(_))
    ));
    for at in [1, 12] {
        let polls = Cell::new(0);
        assert!(matches!(
            accepted.solve_creep_step(&load(FORCE), controls(), || {
                polls.set(polls.get() + 1);
                polls.get() >= at
            }),
            Err(FrameCreepError::Cancelled)
        ));
    }
    assert_eq!(accepted.displacements(), &[0.0; 6]);
    assert_eq!(accepted.creep_recovery_fiber(0, 0, 0).unwrap().1, &virgin);
    let mut untimed = accepted.clone();
    assert!(matches!(
        untimed.static_load(&load(FORCE), 1),
        Err(SolidError::MaterialRefused { .. })
    ));
    assert!(matches!(
        untimed.newmark(&[0.0, 0.0], DT, (1.0, 0.0), Rayleigh { a0: 0.0, a1: 0.0 }),
        Err(SolidError::MaterialRefused { .. })
    ));
    assert!(
        accepted
            .solve_creep_step(&load(FORCE), controls(), || false)
            .is_ok()
    );
}
