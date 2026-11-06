#![cfg(feature = "experimental-euler")]
//! Native evolving-species gravity: analytical rest, mechanical density,
//! conservative gravitational work, and atomic original Euler acceptance.

use fs_evidence::ValidityDomain;
use fs_gas::dynamics::IsothermalAtmosphere;
use fs_gas::multicomponent::{
    Boundaries, Boundary, CaloricDomain, Conserved, Error, IsothermalMixtureGravity, Mesh1d,
    MixtureEuler1d, MixtureGas, Reconstruction, StepControl,
};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Acceleration, Length, MolarMass, Pressure, Temperature, Time};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
};

const CLOSED: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting,
    right: Boundary::Reflecting,
};

// Original manufactured source construction shared in form with the native
// multicomponent tests; coefficients are numerical fixtures, not measured gas data.
fn source(
    name: &str,
    mass: f64,
    cp: f64,
    cp_linear: f64,
    formation: f64,
) -> Nasa9StandardStateModelV1 {
    let a = [0.0, 0.0, cp, cp_linear, 0.0, 0.0, 0.0, formation, 0.0];
    let mut parameters = std::collections::BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(a) {
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
            source: "manufactured multicomponent Euler caloric source".into(),
            license: "fixture".into(),
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
            ElementalReferenceIdV1::new("multicomponent-euler-fixture").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn bank(models: Vec<Nasa9StandardStateModelV1>) -> IdealGasEquilibriumV1 {
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1; ids.len()]],
        )
        .unwrap(),
        ChargeVector::new(ids.clone(), vec![0; ids.len()]).unwrap(),
        models,
    )
    .unwrap()
}
fn gas(models: Vec<Nasa9StandardStateModelV1>) -> MixtureGas {
    MixtureGas::new(
        bank(models),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}

fn material(formation: f64) -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.0, formation),
        source("B", 2.0 * R, 3.8, 0.15, -2.0 * formation),
        source("absent-C", 0.5 * R, 4.1, 0.1, 3.0 * formation),
    ])
}

fn mesh(n: usize, shift: f64) -> Mesh1d {
    let positions: Vec<_> = (0..=n)
        .map(|i| shift + (i as f64 / n as f64).powi(2))
        .collect();
    let volumes = positions.windows(2).map(|p| 0.7 * (p[1] - p[0])).collect();
    Mesh1d::new(positions, vec![0.7; n + 1], volumes).unwrap()
}

fn reference(g: f64, shift: f64) -> IsothermalAtmosphere {
    IsothermalAtmosphere {
        acceleration: Acceleration::new(g),
        temperature: Temperature::new(2.0),
        reference_pressure: Pressure::new(3.0),
        reference_position: Length::new(shift),
    }
}

fn gravity(n: usize, g: f64, formation: f64) -> IsothermalMixtureGravity {
    IsothermalMixtureGravity::try_new(
        material(formation),
        mesh(n, 0.0),
        reference(g, 0.0),
        &[0.4, 0.6, 0.0],
        || false,
    )
    .unwrap()
}

fn control() -> StepControl {
    StepControl {
        reconstruction: Reconstruction::FirstOrder,
        ..StepControl::default()
    }
}

fn snapshot(flow: &MixtureEuler1d) -> Vec<f64> {
    (0..flow.cell_count())
        .flat_map(|i| {
            let q = flow.cell(i).unwrap();
            q.partial_densities
                .iter()
                .copied()
                .chain([q.momentum, q.energy])
                .collect::<Vec<_>>()
        })
        .collect()
}

fn disturbed(reference: &IsothermalMixtureGravity, cold: bool) -> MixtureEuler1d {
    let n = reference.mesh().cell_count();
    let values: Vec<_> = (0..n)
        .map(|i| {
            let (q, p) = reference.reference_cell(i).unwrap();
            if cold {
                let pressure = (q.partial_densities[0] + 0.5 * q.partial_densities[1]) * 0.001001;
                reference
                    .gas()
                    .conserved(q.partial_densities, 0.0, pressure)
                    .unwrap()
            } else if i >= n / 3 && i < 2 * n / 3 {
                let fractions = [0.8, 0.2, 0.0];
                let rho = p.pressure / ((fractions[0] + 0.5 * fractions[1]) * 2.0);
                reference
                    .gas()
                    .conserved(&fractions.map(|y| y * rho), 0.0, p.pressure)
                    .unwrap()
            } else {
                Conserved {
                    partial_densities: q.partial_densities.to_vec(),
                    momentum: q.momentum,
                    energy: q.energy,
                }
            }
        })
        .collect();
    MixtureEuler1d::new_on_mesh(
        reference.gas().clone(),
        reference.mesh().clone(),
        values,
        CLOSED,
        || false,
    )
    .unwrap()
}

fn near(a: f64, b: f64, scale: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * scale.max(b.abs()).max(1e-300),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_integrated_original_composition_atmosphere_remains_at_exact_discrete_rest() {
    for g in [-4.0, 4.0, 0.0, 1e-9] {
        let reference = gravity(19, g, 1e12);
        assert_eq!(reference.reference_mass_fractions(), &[0.4, 0.6, 0.0]);
        let mut flow = reference.initial_flow(|| false).unwrap();
        let initial = snapshot(&flow);
        let beta = g / 1.4;
        let expected_mass = 0.7
            * (3.0 / 1.4)
            * if beta == 0.0 {
                1.0
            } else {
                beta.exp_m1() / beta
            };
        near(flow.totals().mass(), expected_mass, expected_mass, 5e-14);
        for i in 0..flow.cell_count() {
            let (q, p) = reference.reference_cell(i).unwrap();
            let left = flow.mesh().positions()[i];
            let dx = flow.mesh().positions()[i + 1] - left;
            let expected = (3.0 / 1.4)
                * (beta * left).exp()
                * if beta == 0.0 {
                    1.0
                } else {
                    (beta * dx).exp_m1() / (beta * dx)
                };
            near(q.density(), expected, expected, 5e-14);
            near(p.temperature, 2.0, 2.0, 2e-13);
            assert_eq!(q.partial_densities[2], 0.0);
        }
        for _ in 0..48 {
            let step = flow
                .step_with_isothermal_gravity(Time::new(0.02), control(), &reference, || false)
                .unwrap();
            assert_eq!(snapshot(&flow), initial);
            assert_eq!(step.transport.source_change.work, 0.0);
            assert_eq!(step.potential_energy_change_j, 0.0);
            assert_eq!(step.total_energy_residual_j, 0.0);
            near(
                step.transport.source_change.impulse,
                step.physical_gravity_impulse_kg_m_s
                    + step.reference_traction_impulse_defect_kg_m_s,
                step.transport.source_change.impulse.abs(),
                3e-14,
            );
        }
    }
}

#[test]
fn g1_composition_buoyancy_uses_actual_mass_and_conserves_gas_plus_potential_energy() {
    let mut original = None;
    for formation in [0.0, 1e12] {
        let reference = gravity(24, -4.0, formation);
        let mut flow = disturbed(&reference, false);
        let initial = snapshot(&flow);
        let before = flow.totals();
        let potential_before = reference.potential_energy(&flow, || false).unwrap();
        assert!(before.mass() < reference.initial_flow(|| false).unwrap().totals().mass());
        for iteration in 0..40 {
            let step = flow
                .step_with_isothermal_gravity(Time::new(0.001), control(), &reference, || false)
                .unwrap();
            let report = step.transport;
            near(
                step.physical_gravity_impulse_kg_m_s,
                -4.0 * before.mass() * report.elapsed.value(),
                report.source_change.impulse.abs(),
                3e-13,
            );
            near(
                report.source_change.impulse,
                step.physical_gravity_impulse_kg_m_s
                    + step.reference_traction_impulse_defect_kg_m_s,
                report.source_change.impulse.abs(),
                3e-14,
            );
            near(
                report.source_change.work + step.potential_energy_change_j,
                0.0,
                report.before.energy.abs(),
                2e-12,
            );
            near(
                step.total_energy_residual_j,
                0.0,
                report.before.energy.abs(),
                3e-12,
            );
            for residual in report.conservation_residual.species_masses {
                near(residual, 0.0, before.mass(), 2e-13);
            }
            if iteration == 0 {
                assert!(flow.cell(flow.cell_count() / 2).unwrap().momentum > 0.0);
            }
            for i in 0..flow.cell_count() {
                let q = flow.cell(i).unwrap();
                assert_eq!(q.partial_densities[2], 0.0);
                assert!(q.partial_densities.iter().all(|rho| *rho >= 0.0));
                flow.gas().primitive(q, || false).unwrap();
            }
        }
        assert_ne!(snapshot(&flow), initial);
        let after = flow.totals();
        for k in 0..3 {
            near(
                after.species_masses[k],
                before.species_masses[k],
                before.mass(),
                3e-13,
            );
        }
        let potential_after = reference.potential_energy(&flow, || false).unwrap();
        near(
            after.energy + potential_after,
            before.energy + potential_before,
            before.energy.abs() + potential_before.abs(),
            4e-12,
        );
        if let Some(expected) = &original {
            assert_eq!(&snapshot(&flow), expected);
        } else {
            original = Some(snapshot(&flow));
        }
    }
}

#[test]
fn g0_original_cfl_caloric_retries_and_late_cancellation_remain_atomic() {
    let reference = gravity(10, -4.0, 0.0);
    let mut cold = disturbed(&reference, true);
    let step = cold
        .step_with_isothermal_gravity(Time::new(0.05), control(), &reference, || false)
        .unwrap();
    assert!(step.transport.halvings > 0);
    assert!(step.transport.elapsed.value() < 0.05);
    assert!(step.transport.maximum_cfl <= control().cfl * (1.0 + 1e-13));
    for i in 0..cold.cell_count() {
        let p = cold
            .gas()
            .primitive(cold.cell(i).unwrap(), || false)
            .unwrap();
        assert!(p.temperature >= cold.gas().domain().minimum_temperature);
    }
    let mut trial = reference.initial_flow(|| false).unwrap();
    let mut total_polls = 0;
    trial
        .step_with_isothermal_gravity(Time::new(0.001), control(), &reference, || {
            total_polls += 1;
            false
        })
        .unwrap();
    let mut cancelled = reference.initial_flow(|| false).unwrap();
    let initial = snapshot(&cancelled);
    let initial_time = cancelled.time();
    let mut polls = 0;
    assert!(matches!(
        cancelled.step_with_isothermal_gravity(Time::new(0.001), control(), &reference, || {
            polls += 1;
            polls >= total_polls - 1
        }),
        Err(Error::Cancelled)
    ));
    assert_eq!(snapshot(&cancelled), initial);
    assert_eq!(cancelled.time(), initial_time);
    let bad = StepControl {
        reconstruction: Reconstruction::Muscl,
        ..control()
    };
    assert!(
        cancelled
            .step_with_isothermal_gravity(Time::new(0.001), bad, &reference, || false)
            .is_err()
    );
    let wrong_source = gravity(10, -4.0, 1e6);
    assert!(
        cancelled
            .step_with_isothermal_gravity(Time::new(0.001), control(), &wrong_source, || false)
            .is_err()
    );
    cancelled
        .set_boundaries(
            Boundaries::Ends {
                left: Boundary::Transmissive,
                right: Boundary::Reflecting,
            },
            || false,
        )
        .unwrap();
    assert!(
        cancelled
            .step_with_isothermal_gravity(Time::new(0.001), control(), &reference, || false)
            .is_err()
    );
    assert_eq!(snapshot(&cancelled), initial);
    assert_eq!(cancelled.time(), initial_time);
}

#[test]
fn g0_zero_gravity_preserves_native_unforced_bits_and_unsupported_geometry_refuses() {
    let reference = gravity(12, 0.0, 0.0);
    let mut ordinary = disturbed(&reference, false);
    let mut selected = disturbed(&reference, false);
    let a = ordinary
        .step(Time::new(0.001), control(), || false)
        .unwrap();
    let b = selected
        .step_with_isothermal_gravity(Time::new(0.001), control(), &reference, || false)
        .unwrap();
    assert_eq!(a, b.transport);
    assert_eq!(snapshot(&ordinary), snapshot(&selected));
    assert_eq!(b.potential_energy_before_j, 0.0);
    assert_eq!(b.potential_energy_after_j, 0.0);
    for fractions in [[0.4, 0.4, 0.0], [0.4, -0.1, 0.7], [f64::NAN, 0.6, 0.0]] {
        assert!(
            IsothermalMixtureGravity::try_new(
                material(0.0),
                mesh(10, 0.0),
                reference.atmosphere(),
                &fractions,
                || false,
            )
            .is_err()
        );
    }
    let tapered =
        Mesh1d::new(vec![0.0, 0.5, 1.0], vec![0.7, 0.8, 0.9], vec![0.375, 0.425]).unwrap();
    let wrong_volume = Mesh1d::new(vec![0.0, 0.5, 1.0], vec![0.7; 3], vec![0.3, 0.35]).unwrap();
    for unsupported in [tapered, wrong_volume] {
        assert!(
            IsothermalMixtureGravity::try_new(
                material(0.0),
                unsupported,
                reference.atmosphere(),
                &[0.4, 0.6, 0.0],
                || false,
            )
            .is_err()
        );
    }
}
