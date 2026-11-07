//! Original inertial polarization deposits only friction into caloric storage.

use super::*;
use fs_em::electrostatic::dispersion::lorentz::PlanarLorentzDielectric;
use fs_material::field_storage::dielectric::lorentz::{
    LorentzPoleProperties, LorentzState, ResolvedLorentzMaterial,
    leaky::ResolvedLeakyLorentzMaterial,
};
use fs_qty::{Dims, semantic::FrequencyConvention};

fn resonant_card(base: MaterialCard, gamma: f64, varying: bool) -> MaterialCard {
    let mut claims = base.claims().clone();
    let angular = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Frequency(FrequencyConvention::Angular),
        ValueForm::Static,
    ));
    for (name, quantity, coefficient) in [
        ("omega0", angular, 10.0),
        (
            "gamma",
            QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
            gamma,
        ),
    ] {
        let variable = varying && name == "gamma";
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: if variable {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(250.0, coefficient), (500.0, 1.2 * coefficient)],
                        dims: quantity.dims(),
                    }
                } else {
                    PropertyValue::Scalar {
                        value: coefficient,
                        dims: quantity.dims(),
                    }
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    250.0,
                    500.0,
                ),
                interpolation: if variable {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic inertial dielectric heat fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

fn inertial_field(
    card: &MaterialCard,
    geometry: &PlanarDebyeDielectric,
) -> PlanarLorentzDielectric {
    let material = ResolvedLorentzMaterial::resolve(
        card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &[LorentzPoleProperties {
            permittivity_increment: "delta".into(),
            angular_frequency: "omega0".into(),
            damping_rate: "gamma".into(),
        }],
    )
    .unwrap();
    PlanarLorentzDielectric::new(
        geometry.mesh().clone(),
        vec![material; 2],
        geometry.voltage_v().to_vec(),
        vec![
            vec![LorentzState {
                polarization_c_m2: [0.4, 0.0],
                current_a_m2: [3.0, 0.0]
            }];
            2
        ],
        &[0, 1, 2, 3],
        0.0,
        WeightedAssemblyLimits {
            max_cells: 10,
            max_dofs: 20,
            max_triplets: 200,
        },
        field_controls(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn original_lorentz_friction_heats_the_receiver_without_debiting_kinetic_storage() {
    for gamma in [0.0, 2.0] {
        let cards = [resonant_card(
            card(2.0, 3.0, 0.2, 10.0, 1.0, 500.0, false),
            gamma,
            false,
        )];
        with_case(
            &cards,
            &[300.0],
            true,
            false,
            |network, mut thermal, geometry| {
                let mut field = inertial_field(&cards[0], &geometry);
                let dt = 0.03;
                let prepared = field.prepare_step(dt, &mut || false).unwrap();
                let initial = field.stored_energy_j() + thermal.energy_j();
                let initial_thermal = thermal.energy_j();
                let (mut heat, mut work) = (0.0, 0.0);
                for _ in 0..8 {
                    let proposal = field
                        .step_candidate(&prepared, &ends(1), field_controls(), &mut || false)
                        .unwrap();
                    let before_thermal = thermal.energy_j();
                    let next = network
                        .step_dielectric_field_candidate(
                            &thermal,
                            &proposal,
                            controls(),
                            &[0.0],
                            || false,
                        )
                        .unwrap();
                    near(
                        next.thermal().energy_j() - before_thermal,
                        proposal.step().polarization_heat_j,
                        2e-8,
                    );
                    near(
                        next.report().coupled.deposited_joule_j,
                        proposal.step().polarization_heat_j,
                        1e-12,
                    );
                    assert_eq!(proposal.step().joule_heat_j, 0.0);
                    assert!(field.same_problem(next.field()));
                    assert_ne!(
                        next.field().polarization_states()[0][0].current_a_m2[0],
                        0.0
                    );
                    let snapshot = &next.report().coupled.source.dielectric[0];
                    near(
                        snapshot.query().axes()["T"],
                        next.thermal().temperature_k(CellId(7)).unwrap(),
                        1e-10,
                    );
                    assert!(snapshot.property("omega0").is_some());
                    assert!(snapshot.property("gamma").is_some());
                    heat += proposal.step().polarization_heat_j;
                    work += proposal.step().supplied_work_j;
                    (thermal, field, _) = next.into_parts();
                }
                near(thermal.energy_j() - initial_thermal, heat, 4e-8);
                near(
                    thermal.energy_j() + field.stored_energy_j() - initial,
                    work,
                    4e-8,
                );
                if gamma == 0.0 {
                    assert_eq!(heat, 0.0);
                    near(thermal.temperature_k(CellId(7)).unwrap(), 300.0, 1e-11);
                } else {
                    assert!(heat > 0.0);
                    assert!(thermal.temperature_k(CellId(7)).unwrap() > 300.0);
                }
            },
        );
    }
}

#[test]
fn changing_lorentz_coefficients_and_cancellation_cannot_publish_a_thermal_candidate() {
    for variable in [false, true] {
        let cards = [resonant_card(
            card(2.0, 3.0, 0.2, 10.0, 1.0, 500.0, false),
            2.0,
            variable,
        )];
        with_case(
            &cards,
            &[300.0],
            true,
            false,
            |network, thermal, geometry| {
                let field = inertial_field(&cards[0], &geometry);
                let prepared = field.prepare_step(0.03, &mut || false).unwrap();
                let proposal = field
                    .step_candidate(&prepared, &ends(1), field_controls(), &mut || false)
                    .unwrap();
                let initial_energy = thermal.energy_j();
                let history = field.polarization_states().to_vec();
                let result = network.step_dielectric_field_candidate(
                    &thermal,
                    &proposal,
                    controls(),
                    &[0.0],
                    || !variable,
                );
                assert!(result.is_err());
                assert_eq!(thermal.energy_j(), initial_energy);
                assert_eq!(thermal.time_s(), 0.0);
                assert_eq!(field.time_s(), 0.0);
                assert_eq!(field.polarization_states(), history.as_slice());
            },
        );
    }
}

fn conducting_inertial_field(
    card: &MaterialCard,
    geometry: &PlanarDebyeDielectric,
) -> PlanarLorentzDielectric {
    let source = ResolvedLeakyLorentzMaterial::resolve(
        card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &[LorentzPoleProperties {
            permittivity_increment: "delta".into(),
            angular_frequency: "omega0".into(),
            damping_rate: "gamma".into(),
        }],
        "sigma",
    )
    .unwrap();
    PlanarLorentzDielectric::new_leaky(
        geometry.mesh().clone(),
        vec![source; 2],
        geometry.voltage_v().to_vec(),
        vec![
            vec![LorentzState {
                polarization_c_m2: [0.4, 0.0],
                current_a_m2: [3.0, 0.0],
            }];
            2
        ],
        &[0, 1, 2, 3],
        0.0,
        WeightedAssemblyLimits {
            max_cells: 10,
            max_dofs: 20,
            max_triplets: 200,
        },
        field_controls(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn conducting_lorentz_history_deposits_each_physical_loss_once() {
    for sigma in [0.0, 4.0] {
        for gamma in [0.0, 2.0] {
            let cards = [with_leakage(
                resonant_card(card(2.0, 3.0, 0.2, 10.0, 1.0, 500.0, false), gamma, false),
                sigma,
                false,
            )];
            with_case(
                &cards,
                &[300.0],
                true,
                false,
                |network, mut thermal, geometry| {
                    let mut field = conducting_inertial_field(&cards[0], &geometry);
                    let dt = 0.03;
                    let prepared = field.prepare_step(dt, &mut || false).unwrap();
                    let total0 = field.stored_energy_j() + thermal.energy_j();
                    let heat0 = thermal.energy_j();
                    let (mut total_heat, mut total_work) = (0.0, 0.0);
                    for _ in 0..6 {
                        let proposal = field
                            .step_candidate(&prepared, &ends(1), field_controls(), &mut || false)
                            .unwrap();
                        let step = proposal.step();
                        let expected_friction: f64 = field
                            .polarization_states()
                            .iter()
                            .zip(proposal.field().polarization_states())
                            .zip(field.cell_volumes_m3())
                            .map(|((before, after), volume)| {
                                let j =
                                    0.5 * (before[0].current_a_m2[0] + after[0].current_a_m2[0]);
                                volume * dt * gamma * j * j / (3.0 * 100.0)
                            })
                            .sum();
                        near(step.polarization_heat_j, expected_friction, 2e-12);
                        near(step.joule_heat_j, dt * sigma, 2e-12);
                        near(step.heat_j, expected_friction + dt * sigma, 2e-12);
                        let next = network
                            .step_dielectric_field_candidate(
                                &thermal,
                                &proposal,
                                controls(),
                                &[0.0],
                                || false,
                            )
                            .unwrap();
                        near(
                            next.thermal().energy_j() - thermal.energy_j(),
                            step.heat_j,
                            2e-8,
                        );
                        near(next.report().coupled.deposited_joule_j, step.heat_j, 2e-12);
                        for snapshot in &next.report().coupled.source.dielectric {
                            for key in ["epsilon_inf", "delta", "omega0", "gamma", "sigma"] {
                                assert!(snapshot.property(key).is_some());
                            }
                            near(snapshot.property("sigma").unwrap().value_si(), sigma, 0.0);
                            near(
                                snapshot.query().axes()["T"],
                                next.thermal().temperature_k(CellId(7)).unwrap(),
                                1e-10,
                            );
                        }
                        assert!(field.same_problem(next.field()));
                        assert!(
                            next.field().polarization_states()[0][0].current_a_m2[0].abs() > 0.0
                        );
                        total_heat += step.heat_j;
                        total_work += step.supplied_work_j;
                        (thermal, field, _) = next.into_parts();
                    }
                    near(thermal.energy_j() - heat0, total_heat, 5e-8);
                    near(
                        thermal.energy_j() + field.stored_energy_j() - total0,
                        total_work,
                        5e-8,
                    );
                    if gamma == 0.0 && sigma == 0.0 {
                        assert_eq!(total_heat, 0.0);
                    } else {
                        assert!(total_heat > 0.0);
                    }
                },
            );
        }
    }
}

#[test]
fn varying_leakage_and_cancellation_refuse_without_consuming_resonant_heat() {
    for sigma in [0.0, 4.0] {
        for variable in [false, true] {
            let cards = [with_leakage(
                resonant_card(card(2.0, 3.0, 0.2, 10.0, 1.0, 500.0, false), 2.0, false),
                sigma,
                variable,
            )];
            with_case(
                &cards,
                &[300.0],
                true,
                false,
                |network, thermal, geometry| {
                    let field = conducting_inertial_field(&cards[0], &geometry);
                    let prepared = field.prepare_step(0.03, &mut || false).unwrap();
                    let proposal = field
                        .step_candidate(&prepared, &ends(1), field_controls(), &mut || false)
                        .unwrap();
                    let history = field.polarization_states().to_vec();
                    let energy = thermal.energy_j();
                    let rejected = network.step_dielectric_field_candidate(
                        &thermal,
                        &proposal,
                        controls(),
                        &[0.0],
                        || !variable,
                    );
                    assert!(rejected.is_err());
                    assert_eq!(field.time_s(), 0.0);
                    assert_eq!(thermal.time_s(), 0.0);
                    assert_eq!(field.polarization_states(), history.as_slice());
                    assert_eq!(thermal.energy_j(), energy);
                    if !variable {
                        let retried = network
                            .step_dielectric_field_candidate(
                                &thermal,
                                &proposal,
                                controls(),
                                &[0.0],
                                || false,
                            )
                            .unwrap();
                        near(
                            retried.thermal().energy_j() - energy,
                            proposal.step().heat_j,
                            2e-8,
                        );
                    }
                },
            );
        }
    }
}
