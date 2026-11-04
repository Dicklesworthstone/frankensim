#![cfg(feature = "experimental-euler")]

//! G0/G1/G3/G4: original-source mechanical work, transport and atomic refusal.

use fs_evidence::ValidityDomain;
use fs_gas::body_force::{BodyForce, ForceChange, ForceControl};
use fs_gas::dynamics::{
    Boundaries, Boundary, CaloricDomain, Euler1d, FrozenGas, IsothermalAtmosphere, Mesh1d,
    StepControl,
};
use fs_gas::multicomponent::{MixtureEuler1d, MixtureGas};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Acceleration, Length, MolarMass, Pressure, Temperature, Time};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    ChargeVector, Composition, CompositionBasis, ElementId, ElementalMatrix,
    ElementalReferenceIdV1, Nasa9RegionV1, Nasa9StandardStateModelV1, ReferenceEquationOfStateV1,
    SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
};

const FORCE: ForceControl = ForceControl { maximum_mach: 80.0 };

fn source(
    name: &str,
    mass: f64,
    cp: f64,
    linear: f64,
    formation: f64,
) -> Nasa9StandardStateModelV1 {
    let coefficients = [0.0, 0.0, cp, linear, 0.0, 0.0, 0.0, formation, 0.0];
    let mut parameters = std::collections::BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(name.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 1.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()),
        law_version: 1,
        parameters,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 0.001, 50.0),
        sources: vec![],
        provenance: Provenance {
            source: "manufactured body-force caloric law".into(),
            license: "test fixture".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        SpeciesId::new(name).unwrap(),
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(1.0),
            ElementalReferenceIdV1::new("body-force-fixture").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}

fn frozen_gas() -> FrozenGas {
    FrozenGas::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![source("A", R, 3.5, 0.2, 1e18)],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}

fn mixture_gas() -> MixtureGas {
    let models = vec![
        source("A", R, 3.5, 0.2, 1e18),
        source("B", 2.0 * R, 4.0, 0.1, -2e18),
        source("C", R, 3.0, 0.3, 3e18),
    ];
    let ids: Vec<_> = models.iter().map(|model| model.species().clone()).collect();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1, 2, 1]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0, 0, 0]).unwrap(),
        models,
    )
    .unwrap();
    MixtureGas::new(
        bank,
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}

fn frozen(velocity: f64) -> Euler1d {
    let gas = frozen_gas();
    let states = vec![gas.conserved(1.0, velocity, 1.7).unwrap(); 8];
    Euler1d::new(
        gas,
        Length::new(0.125),
        0.7,
        states,
        Boundaries::Periodic,
        || false,
    )
    .unwrap()
}

fn mixture(velocity: f64) -> MixtureEuler1d {
    let gas = mixture_gas();
    let states = vec![gas.conserved(&[0.8, 0.2, 0.0], velocity, 2.7).unwrap(); 8];
    MixtureEuler1d::new(
        gas,
        Length::new(0.125),
        0.7,
        states,
        Boundaries::Periodic,
        || false,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_frozen_exact_work_includes_braking_reversal_and_galilean_boost() {
    for boost in [0.0, 25.0] {
        let mut flow = frozen(-1.2 + boost);
        let before = flow.totals();
        let exchange = flow
            .apply_body_force(
                Time::new(1.0),
                BodyForce::UniformAcceleration(2.0),
                FORCE,
                || false,
            )
            .unwrap();
        let after = flow.totals();
        close(exchange.impulse, 2.0 * before.mass, 2e-14);
        let exact_work = 0.5 * before.mass * ((0.8 + boost).powi(2) - (-1.2 + boost).powi(2));
        close(exchange.work, exact_work, 3e-12);
        close(
            after.energy - before.energy - exchange.work,
            exchange.energy_residual,
            2e-13,
        );
        assert_eq!(flow.time(), Time::new(0.0));
        for &q in flow.states() {
            let p = flow.gas().primitive(q, || false).unwrap();
            close(p.velocity, 0.8 + boost, 1e-13);
            close(p.temperature, 1.7, 2e-10);
            close(p.pressure, 1.7, 2e-10);
            assert_eq!(q.density, 1.0);
        }
    }
}

#[test]
fn g1_mixture_force_density_uses_each_mass_and_preserves_species_calorics() {
    let gas = mixture_gas();
    let partials = [[0.8, 0.2, 0.0], [0.4, 1.6, 0.0]];
    let velocities = [-0.6, 0.7];
    let states = partials
        .iter()
        .zip(velocities)
        .map(|(d, u)| gas.conserved(d, u, 3.0 * (d[0] + 0.5 * d[1])).unwrap())
        .collect();
    let mut flow = MixtureEuler1d::new(
        gas,
        Length::new(0.4),
        0.3,
        states,
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    let forces = [4.0, -3.0];
    let before = flow.totals();
    let source = flow
        .apply_body_force(
            Time::new(0.4),
            BodyForce::ForceDensity(&forces),
            FORCE,
            || false,
        )
        .unwrap();
    let mut expected_work = 0.0;
    for i in 0..2 {
        let q = flow.cell(i).unwrap();
        assert_eq!(q.partial_densities, &partials[i]);
        let p = flow.gas().primitive(q, || false).unwrap();
        let density: f64 = partials[i].iter().sum();
        let velocity = velocities[i] + forces[i] * 0.4 / density;
        close(p.velocity, velocity, 1e-13);
        close(p.temperature, 3.0, 1e-10);
        expected_work += 0.5 * density * (velocity * velocity - velocities[i].powi(2)) * 0.12;
    }
    close(source.impulse, 0.4 * (forces[0] + forces[1]) * 0.12, 1e-14);
    close(source.work, expected_work, 1e-13);
    let after = flow.totals();
    assert_eq!(before.species_masses, after.species_masses);
    close(
        after.energy - before.energy - source.work,
        source.energy_residual,
        4e-15,
    );
    let source = flow
        .apply_body_force(
            Time::new(0.1),
            BodyForce::Acceleration(&[0.5, -2.0]),
            FORCE,
            || false,
        )
        .unwrap();
    close(source.impulse, 0.12 * (0.05 - 2.0 * 0.2), 1e-14);
    for (i, expected_velocity) in [1.05, -0.1].into_iter().enumerate() {
        let state = flow.cell(i).unwrap();
        assert_eq!(state.partial_densities, &partials[i]);
        let p = flow.gas().primitive(state, || false).unwrap();
        close(p.velocity, expected_velocity, 1e-13);
        close(p.temperature, 3.0, 1e-10);
    }
}

#[test]
fn g0_unrepresented_impulse_is_not_a_heat_source() {
    let mut flow = frozen(1.0);
    let original = flow.states().to_vec();
    let source = flow
        .apply_body_force(
            Time::new(1.0),
            BodyForce::UniformForceDensity(1e-17),
            FORCE,
            || false,
        )
        .unwrap();
    assert_eq!(source, ForceChange::default());
    assert_eq!(flow.states(), original);
    let source = flow
        .apply_body_force(
            Time::new(1.0),
            BodyForce::UniformForceDensity(5e-16),
            FORCE,
            || false,
        )
        .unwrap();
    let impulse_density = flow.states()[0].momentum - original[0].momentum;
    assert!(impulse_density > 0.0);
    assert_ne!(impulse_density, 5e-16);
    close(source.impulse, 0.7 * impulse_density, 1e-30);
    let midpoint = 0.5 * (1.0 + flow.states()[0].momentum);
    close(source.work, source.impulse * midpoint, 1e-30);
    let mut mixture = mixture(1.0);
    let before: Vec<_> = (0..8).map(|i| mixture.cell(i).unwrap().momentum).collect();
    let source = mixture
        .apply_body_force(
            Time::new(1.0),
            BodyForce::UniformAcceleration(1e-17),
            FORCE,
            || false,
        )
        .unwrap();
    assert_eq!(source, ForceChange::default());
    for (i, momentum) in before.into_iter().enumerate() {
        assert_eq!(mixture.cell(i).unwrap().momentum, momentum);
    }
}

#[test]
fn g1_native_frozen_and_mixture_steps_integrate_actual_accepted_work() {
    let mut frozen = frozen(-0.4);
    let mut mixture = mixture(-0.4);
    let initial_frozen = frozen.totals();
    let initial_mixture = mixture.totals();
    let (mut frozen_work, mut mixture_work) = (0.0, 0.0);
    for _ in 0..20 {
        let a = frozen
            .step_with_body_force(
                Time::new(0.001),
                StepControl::default(),
                BodyForce::UniformAcceleration(5.0),
                FORCE,
                || false,
            )
            .unwrap();
        let b = mixture
            .step_with_body_force(
                Time::new(0.001),
                StepControl::default(),
                BodyForce::UniformAcceleration(5.0),
                FORCE,
                || false,
            )
            .unwrap();
        assert_eq!(a.elapsed, Time::new(0.001));
        assert_eq!(b.elapsed, a.elapsed);
        close(a.source_change.impulse, initial_frozen.mass * 0.005, 1e-14);
        close(
            b.source_change.impulse,
            initial_mixture.mass() * 0.005,
            1e-14,
        );
        close(a.conservation_residual.momentum, 0.0, 1e-13);
        close(a.conservation_residual.energy, 0.0, 1e-13);
        close(b.conservation_residual.momentum, 0.0, 1e-13);
        close(b.conservation_residual.energy, 0.0, 1e-13);
        frozen_work += a.source_change.work;
        mixture_work += b.source_change.work;
    }
    let work_per_mass = 0.5 * ((-0.3_f64).powi(2) - (-0.4_f64).powi(2));
    close(frozen_work, initial_frozen.mass * work_per_mass, 1e-12);
    close(mixture_work, initial_mixture.mass() * work_per_mass, 1e-12);
    close(
        frozen.totals().energy - initial_frozen.energy,
        frozen_work,
        1e-12,
    );
    close(
        mixture.totals().energy - initial_mixture.energy,
        mixture_work,
        1e-12,
    );
    for i in 0..8 {
        let p = frozen
            .gas()
            .primitive(frozen.states()[i], || false)
            .unwrap();
        close(p.temperature, 1.7, 1e-9);
        close(p.velocity, -0.3, 1e-12);
        let q = mixture.cell(i).unwrap();
        assert_eq!(q.partial_densities, &[0.8, 0.2, 0.0]);
        let p = mixture.gas().primitive(q, || false).unwrap();
        close(p.temperature, 3.0, 1e-9);
        close(p.velocity, -0.3, 1e-12);
    }
}

#[test]
fn g3_zero_force_retains_bit_identical_native_steps() {
    let mut ordinary = frozen(0.3);
    let mut forced = frozen(0.3);
    let a = ordinary
        .step(Time::new(0.01), StepControl::default(), || false)
        .unwrap();
    let b = forced
        .step_with_body_force(
            Time::new(0.01),
            StepControl::default(),
            BodyForce::Acceleration(&[0.0; 8]),
            FORCE,
            || false,
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(ordinary.states(), forced.states());
    let mut ordinary = mixture(0.3);
    let mut forced = mixture(0.3);
    let a = ordinary
        .step(Time::new(0.01), StepControl::default(), || false)
        .unwrap();
    let b = forced
        .step_with_body_force(
            Time::new(0.01),
            StepControl::default(),
            BodyForce::UniformForceDensity(0.0),
            FORCE,
            || false,
        )
        .unwrap();
    assert_eq!(a, b);
    for i in 0..8 {
        assert_eq!(ordinary.cell(i), forced.cell(i));
    }
}

#[test]
fn g4_force_cfl_retries_restart_both_source_halves() {
    let mut flow = frozen(0.0);
    let before = flow.totals();
    let report = flow
        .step_with_body_force(
            Time::new(0.1),
            StepControl::default(),
            BodyForce::UniformAcceleration(1000.0),
            FORCE,
            || false,
        )
        .unwrap();
    assert!(report.halvings > 0);
    assert!(!report.first_order_retry);
    assert!(report.maximum_cfl <= 0.4 * (1.0 + 64.0 * f64::EPSILON));
    let velocity = 1000.0 * report.elapsed.value();
    close(report.source_change.impulse, before.mass * velocity, 1e-12);
    close(
        report.source_change.work,
        0.5 * before.mass * velocity * velocity,
        1e-12,
    );
    close(report.conservation_residual.energy, 0.0, 1e-12);
    assert_eq!(flow.time(), report.elapsed);
    let mut flow = mixture(0.0);
    let before = flow.totals();
    let report = flow
        .step_with_body_force(
            Time::new(0.1),
            StepControl::default(),
            BodyForce::UniformForceDensity(1000.0),
            FORCE,
            || false,
        )
        .unwrap();
    assert!(report.halvings > 0);
    assert!(!report.first_order_retry);
    assert!(report.maximum_cfl <= 0.4 * (1.0 + 64.0 * f64::EPSILON));
    let velocity = 1000.0 * report.elapsed.value();
    close(
        report.source_change.impulse,
        before.mass() * velocity,
        1e-12,
    );
    close(
        report.source_change.work,
        0.5 * before.mass() * velocity * velocity,
        1e-12,
    );
    close(report.conservation_residual.energy, 0.0, 1e-12);
    assert_eq!(flow.totals().species_masses, before.species_masses);
    assert_eq!(flow.time(), report.elapsed);
}

#[test]
fn g4_invalid_and_late_cancelled_sources_publish_neither_cells_nor_time() {
    let mut reference = frozen(0.2);
    let mut polls = 0;
    let expected = reference
        .step_with_body_force(
            Time::new(0.001),
            StepControl::default(),
            BodyForce::UniformAcceleration(2.0),
            FORCE,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut flow = frozen(0.2);
    let initial = flow.states().to_vec();
    let mut count = 0;
    let error = flow
        .step_with_body_force(
            Time::new(0.001),
            StepControl::default(),
            BodyForce::UniformAcceleration(2.0),
            FORCE,
            || {
                count += 1;
                count == polls
            },
        )
        .unwrap_err();
    assert_eq!(error, fs_gas::dynamics::Error::Cancelled);
    assert_eq!(flow.states(), initial);
    assert_eq!(flow.time(), Time::new(0.0));
    let retried = flow
        .step_with_body_force(
            Time::new(0.001),
            StepControl::default(),
            BodyForce::UniformAcceleration(2.0),
            FORCE,
            || false,
        )
        .unwrap();
    assert_eq!(retried, expected);
    assert_eq!(flow.states(), reference.states());

    let mut flow = mixture(0.2);
    let initial = flow.totals();
    assert!(
        flow.apply_body_force(
            Time::new(0.01),
            BodyForce::ForceDensity(&[1.0; 7]),
            FORCE,
            || false
        )
        .is_err()
    );
    assert!(
        flow.apply_body_force(
            Time::new(0.01),
            BodyForce::UniformAcceleration(f64::NAN),
            FORCE,
            || false
        )
        .is_err()
    );
    assert!(
        flow.apply_body_force(
            Time::new(1.0),
            BodyForce::UniformAcceleration(200.0),
            ForceControl { maximum_mach: 1.0 },
            || false
        )
        .is_err()
    );
    assert_eq!(flow.totals(), initial);
    assert_eq!(flow.time(), Time::new(0.0));
    let mut flow = mixture(0.2);
    let mut reference = mixture(0.2);
    let mut polls = 0;
    let expected = reference
        .step_with_body_force(
            Time::new(0.001),
            StepControl::default(),
            BodyForce::UniformForceDensity(-3.0),
            FORCE,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut count = 0;
    let error = flow
        .step_with_body_force(
            Time::new(0.001),
            StepControl::default(),
            BodyForce::UniformForceDensity(-3.0),
            FORCE,
            || {
                count += 1;
                count == polls
            },
        )
        .unwrap_err();
    assert_eq!(error, fs_gas::multicomponent::Error::Cancelled);
    assert_eq!(flow.totals(), initial);
    assert_eq!(flow.time(), Time::new(0.0));
    let retried = flow
        .step_with_body_force(
            Time::new(0.001),
            StepControl::default(),
            BodyForce::UniformForceDensity(-3.0),
            FORCE,
            || false,
        )
        .unwrap();
    assert_eq!(retried, expected);
}

#[test]
fn g1_gravity_and_additional_body_work_are_separate_sources() {
    let mut flow = Euler1d::new_isothermal_atmosphere(
        frozen_gas(),
        Mesh1d::uniform(8, Length::new(0.125), 0.7).unwrap(),
        IsothermalAtmosphere {
            acceleration: Acceleration::new(-0.2),
            temperature: Temperature::new(1.7),
            reference_pressure: Pressure::new(1.7),
            reference_position: Length::new(0.0),
        },
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        || false,
    )
    .unwrap();
    let initial = flow.totals();
    let potential = flow.gravitational_potential_energy();
    let report = flow
        .step_with_body_force(
            Time::new(0.001),
            StepControl::default(),
            BodyForce::UniformAcceleration(0.3),
            FORCE,
            || false,
        )
        .unwrap();
    let dt = report.elapsed.value();
    close(report.source_change.impulse, 0.3 * initial.mass * dt, 1e-13);
    close(
        report.gravity_momentum_impulse,
        -0.2 * initial.mass * dt,
        1e-12,
    );
    close(report.conservation_residual.energy, 0.0, 1e-12);
    close(
        flow.totals().energy - initial.energy + flow.gravitational_potential_energy()
            - potential
            - report.source_change.work,
        0.0,
        2e-12,
    );
}

#[test]
fn g1_time_dependent_drive_refines_to_exact_velocity_and_work() {
    let error = |steps: usize| {
        let gas = frozen_gas();
        let q = gas.conserved(1.0, 0.0, 1.7).unwrap();
        let mut flow = Euler1d::new(
            gas,
            Length::new(2.0),
            0.7,
            vec![q; 2],
            Boundaries::Periodic,
            || false,
        )
        .unwrap();
        let dt = 1.0 / steps as f64;
        let initial = flow.totals();
        let mut work = 0.0;
        for i in 0..steps {
            let acceleration = ((i as f64 + 0.5) * dt).sin();
            let report = flow
                .step_with_body_force(
                    Time::new(dt),
                    StepControl::default(),
                    BodyForce::UniformAcceleration(acceleration),
                    FORCE,
                    || false,
                )
                .unwrap();
            assert_eq!(report.elapsed, Time::new(dt));
            work += report.source_change.work;
        }
        let exact_velocity = 1.0 - 1.0_f64.cos();
        let velocity = flow.states()[0].momentum / flow.states()[0].density;
        close(work, 0.5 * initial.mass * velocity * velocity, 1e-11);
        close(flow.totals().energy - initial.energy, work, 1e-11);
        (velocity - exact_velocity).abs()
    };
    let errors = [error(8), error(16), error(32)];
    for pair in errors.windows(2) {
        assert!((3.9..4.1).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}
