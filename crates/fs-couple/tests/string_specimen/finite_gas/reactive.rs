//! G1/G3/G4: actual chemistry, continuing sound and one thermal energy account.
//! Manufactured isomers exercise integration, not measured combustion fidelity.

use super::*;
use fs_couple::acoustic_realize::finite_gas::reactive::ReactiveFiniteGasStringSolve;
use fs_couple::cavity_phs::piston::reactive::{
    GasReaction, GasReactionSolve, implicit::coupled::ImplicitMechanismSolve,
};
use fs_qty::Temperature;
use fs_solver::NewtonKrylovConfig;
use fs_thermochem::equilibrium::*;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::*;

#[path = "../../material_exchange/gas.rs"]
#[allow(dead_code)]
mod source;

const RATE: f64 = 200_000.0;
const RELEASE: f64 = 3000.0;
const DT: f64 = 1.0 / 48000.0;
const R: f64 = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;

fn law(release: f64) -> GasReaction {
    let components =
        [("A", 0.0, 1.7e-5), ("B", -release, 4.0e-5)].map(|(name, formation, viscosity)| {
            let (thermo, mut transport) = source::component(name, 0.028, formation);
            let mut claims = ClaimSet::new();
            for (key, dims, value) in [
                ("dynamic_viscosity", DynViscosity::DIMS, viscosity),
                (
                    "thermal_conductivity",
                    fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS,
                    0.03,
                ),
            ] {
                claims
                    .insert_claim(PropertyClaim {
                        key: PropertyKey::new(key, dims),
                        value: PropertyValue::Curve {
                            abscissa: "temperature".into(),
                            abscissa_dims: Temperature::DIMS,
                            knots: vec![(200.0, value), (900.0, 2.0 * value)],
                            dims,
                        },
                        validity: ValidityDomain::unconstrained()
                            .with_quantity(
                                "temperature",
                                QuantitySpec::semantic(SemanticType::new(
                                    QuantityKind::AbsoluteTemperature,
                                    ValueForm::Static,
                                )),
                                200.0,
                                900.0,
                            )
                            .with_quantity(
                                "pressure",
                                QuantitySpec::semantic(SemanticType::new(
                                    QuantityKind::Pressure,
                                    ValueForm::Static,
                                )),
                                1000.0,
                                1e7,
                            ),
                        uncertainty: UncertaintyModel::Unstated,
                        interpolation: InterpolationPolicy::LinearInside,
                        observations: vec![],
                        provenance: Provenance {
                            source: "synthetic reacting string gas".into(),
                            license: "CC0-1.0".into(),
                            artifact: None,
                        },
                    })
                    .unwrap();
            }
            transport.material = MaterialCard::assemble(
                MaterialStateId {
                    chemistry: name.into(),
                    phase: "gas".into(),
                    process: "manufactured".into(),
                    revision: 0,
                },
                claims,
                vec![],
            )
            .unwrap();
            (thermo, transport)
        });
    let species = components
        .iter()
        .map(|(s, _)| s.species().clone())
        .collect::<Vec<_>>();
    let id = ReactionId::new("isomerization").unwrap();
    let matrix =
        StoichiometricMatrix::new(species.clone(), vec![id.clone()], vec![vec![-1], vec![1]])
            .unwrap();
    let elements = ElementalMatrix::new(
        vec![ElementId::new("X").unwrap()],
        species.clone(),
        vec![vec![1, 1]],
    )
    .unwrap();
    let charge = ChargeVector::new(species, vec![0, 0]).unwrap();
    let conserved = verify_conservation(&elements, &matrix, &charge).unwrap();
    let equilibrium = IdealGasReactionEquilibriumV1::new(
        matrix,
        conserved,
        id,
        components.iter().map(|(s, _)| s.clone()).collect(),
    )
    .unwrap();
    let catalog = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            components.iter().map(|(s, _)| s.clone()).collect(),
            Composition::new(CompositionBasis::MoleFraction, vec![0.9, 0.1]).unwrap(),
        )
        .unwrap(),
        components.into_iter().map(|(_, t)| t).collect(),
        None,
    )
    .unwrap();
    GasReaction::new(
        catalog,
        DeclaredStoichiometricMassActionV1::new(
            ThermodynamicReverseRateClosureV1::new(
                equilibrium,
                ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(RATE))
                    .unwrap(),
            ),
            MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients,
        ),
    )
    .unwrap()
}

fn policy() -> ImplicitMechanismSolve {
    ImplicitMechanismSolve {
        extent_brackets_mol: [(ReactionId::new("isomerization").unwrap(), [-0.04, 0.04])].into(),
        extent_tolerance_mol: 1e-12,
        interior_fraction: 0.05,
        jacobian_relative_step: 1e-5,
        maximum_iterations: 32,
        maximum_evaluations: 512,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-13,
            relative_tolerance: 1e-13,
            linear_restart: 8,
            max_linear_cycles: 8,
            forcing_maximum: 1e-6,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn controls() -> ReactiveFiniteGasStringSolve {
    ReactiveFiniteGasStringSolve {
        reaction: GasReactionSolve {
            gas: solve().gas,
            maximum_consumed_fraction: 0.75,
            species_tolerance_mol: 1e-12,
            entropy_tolerance_j_k: 1e-10,
        },
        maximum_gas_temperature_increment_k: 100.0,
    }
}
fn gas(law: &GasReaction) -> IsochoricGasParcel {
    let typed = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    IsochoricGasParcel::new(
        law.catalog().clone(),
        QueryPoint::new()
            .with_quantity(
                "temperature",
                typed(QuantityKind::AbsoluteTemperature),
                300.0,
            )
            .unwrap()
            .with_quantity("pressure", typed(QuantityKind::Pressure), 100000.0)
            .unwrap()
            .with("material-lot", 7.0)
            .unwrap(),
        0.001,
    )
    .unwrap()
}
fn prepared(
    cx: &fs_exec::Cx<'_>,
    law: &GasReaction,
) -> (ThermalMaterialStringRuntime, IsochoricGasParcel) {
    let mut gas = gas(law);
    let mut live = runtime(cx, gas.state().state());
    for epoch in 0..8 {
        live.step_with_isochoric_gas(
            cx,
            epoch,
            &[0.0; 2],
            &mut gas,
            FiniteGasStringHeat::default(),
            solve(),
        )
        .unwrap();
    }
    (live, gas)
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} (tolerance {tolerance})"
    );
}

#[test]
fn g1_heated_chemistry_closes_independent_extent_caloric_pressure_and_total_energy() {
    with_string_cx(|cx, _| {
        let reaction = law(RELEASE).with_implicit_mechanism(policy()).unwrap();
        let (mut live, mut gas) = prepared(cx, &reaction);
        let before = live.clone();
        let old = gas.clone();
        let heat = FiniteGasStringHeat {
            string_heat_j: 5e-8,
            gas_heat_j: 0.125,
        };
        let result = live
            .step_with_reacting_isochoric_gas(
                cx,
                8,
                &[0.0; 2],
                &mut gas,
                heat,
                &reaction,
                controls(),
            )
            .unwrap();
        let frame = &result.coupled;
        let chemical = &result.chemistry.reaction;
        let n = old.mass_kg() / 0.028;
        let xi = chemical.extent_mol;
        let q = heat.gas_heat_j + frame.drag_heat_j;
        let t = old.state().state().temperature + (q + RELEASE * xi) / (n * 2.5 * R);
        let equilibrium = (RELEASE / (R * t)).exp();
        let a = (0.9 * n - xi) * R * t / (old.volume_m3() * 100000.0);
        let b = (0.1 * n + xi) * R * t / (old.volume_m3() * 100000.0);
        near(
            xi,
            DT * old.volume_m3() * RATE * (a - b / equilibrium),
            2e-12,
        );
        near(gas.state().state().temperature, t, 2e-8);
        near(
            gas.state().state().pressure,
            old.state().state().pressure * t / old.state().state().temperature,
            2e-5,
        );
        assert!(xi > 0.0 && frame.drag_heat_j > 0.0);
        for (row, amount) in chemical.species.iter().zip([0.9 * n - xi, 0.1 * n + xi]) {
            near(row.final_mol, amount, 1e-12);
        }
        assert_eq!(gas.mass_kg(), old.mass_kg());
        assert_eq!(gas.volume_m3(), old.volume_m3());
        assert_eq!(gas.point().axes().get("material-lot"), Some(&7.0));
        assert_eq!(
            gas.model().thermodynamics().components(),
            old.model().thermodynamics().components()
        );
        for (new, old) in gas
            .state()
            .components()
            .iter()
            .zip(old.state().components())
        {
            assert_eq!(new.card_identity(), old.card_identity());
        }
        assert_eq!(frame.gas, result.chemistry.heat);
        assert_eq!(frame.gas.supplied_heat_j, q);
        assert_eq!(frame.external_gas_heat_j, heat.gas_heat_j);
        near(
            chemical.stored_energy_change_j,
            q,
            old.mass_kg() * controls().reaction.gas.energy_tolerance_j_kg,
        );
        assert!(
            result.chemistry.entropy_residual_j_k >= -controls().reaction.entropy_tolerance_j_k
        );
        let solid = live.mechanical().specimen().mass_kg()
            * (live.thermal().specific_enthalpy_j_kg() - before.thermal().specific_enthalpy_j_kg());
        let gas_storage = gas.mass_kg()
            * (gas.state().specific_internal_energy_j_kg()
                - old.state().specific_internal_energy_j_kg());
        let dissipation = &frame.string.vibration.dissipation;
        let defect = heat.string_heat_j + heat.gas_heat_j
            - (energy(&live) - energy(&before)
                + solid
                + gas_storage
                + dissipation.authored_loss_j
                + dissipation.roundoff_residual_j);
        assert!(defect.abs() <= frame.energy_tolerance_j);
        assert!(
            (defect + RELEASE * xi).abs() > 100.0 * frame.energy_tolerance_j,
            "crediting formation heat twice must fail"
        );
        assert!(frame.energy_balance_residual_j.abs() <= frame.energy_tolerance_j);
        assert_eq!(chemical.implicit_mechanism.as_ref().unwrap().duration_s, DT);
    });
}

#[test]
fn g3_final_mixture_rebuilds_drag_and_pressure_without_resetting_motion_or_history() {
    with_string_cx(|cx, _| {
        let reaction = law(RELEASE).with_implicit_mechanism(policy()).unwrap();
        let (mut live, mut gas) = prepared(cx, &reaction);
        let mut frozen = live.clone();
        let previous = accelerations(&live);
        let old = gas.clone();
        frozen.step(cx, 8, &[0.0; 2], 0.0).unwrap();
        let result = live
            .step_with_reacting_isochoric_gas(
                cx,
                8,
                &[0.0; 2],
                &mut gas,
                FiniteGasStringHeat::default(),
                &reaction,
                controls(),
            )
            .unwrap();
        assert_eq!(live.mechanical().states(), frozen.mechanical().states());
        assert_eq!(live.thermal(), frozen.thermal());
        assert_eq!(live.mechanical().epoch(), 9);
        assert_eq!(live.mechanical().accepted_samples(), 9);
        assert_eq!(live.mechanical().ambient(), gas.state().state());
        assert!(
            live.mechanical().modes()[0].damping_ratio
                > frozen.mechanical().modes()[0].damping_ratio
        );
        let string = live.mechanical().specimen().string();
        let pi = core::f64::consts::PI;
        let acceleration = accelerations(&live);
        let pressure = |medium: &fs_material::gas::GasState| {
            let factor = medium.density * string.width_m
                / (4.0 * pi * (string.lin_density_kg_m * string.length_m / 2.0).sqrt());
            factor * 2.0 * string.length_m / pi * acceleration[0]
                - factor * string.length_m.powi(2) / (2.0 * pi * medium.sound_speed)
                    * (acceleration[1] - previous[1])
                    / DT
        };
        near(
            result
                .coupled
                .string
                .vibration
                .acoustic
                .observer_pressure_pa,
            pressure(gas.state().state()),
            1e-11,
        );
        assert!((pressure(gas.state().state()) - pressure(old.state().state())).abs() > 1e-12);
        live.step_with_reacting_isochoric_gas(
            cx,
            9,
            &[0.0; 2],
            &mut gas,
            FiniteGasStringHeat::default(),
            &reaction,
            controls(),
        )
        .unwrap();
        frozen.step(cx, 9, &[0.0; 2], 0.0).unwrap();
        assert_ne!(live.mechanical().states(), frozen.mechanical().states());
    });
}

#[test]
fn g1_repeated_samples_retain_species_and_match_the_zero_formation_analytic_endpoint() {
    with_string_cx(|cx, _| {
        let reaction = law(0.0).with_implicit_mechanism(policy()).unwrap();
        let (mut live, mut gas) = prepared(cx, &reaction);
        let n = gas.mass_kg() / 0.028;
        let mut amounts = [0.9 * n, 0.1 * n];
        for epoch in 8..24 {
            let old_temperature = gas.state().state().temperature;
            let heat = FiniteGasStringHeat {
                string_heat_j: 0.0,
                gas_heat_j: 0.01,
            };
            let frame = live
                .step_with_reacting_isochoric_gas(
                    cx,
                    epoch,
                    &[0.0; 2],
                    &mut gas,
                    heat,
                    &reaction,
                    controls(),
                )
                .unwrap();
            let t = old_temperature + (heat.gas_heat_j + frame.coupled.drag_heat_j) / (n * 2.5 * R);
            let h = DT * RATE * R * t / 100000.0;
            let xi = h * (amounts[0] - amounts[1]) / (1.0 + 2.0 * h);
            near(frame.chemistry.reaction.extent_mol, xi, 2e-12);
            near(gas.state().state().temperature, t, 2e-8);
            for (row, expected) in frame.chemistry.reaction.species.iter().zip(amounts) {
                near(row.initial_mol, expected, 4e-12);
            }
            amounts[0] -= xi;
            amounts[1] += xi;
            assert_eq!(live.mechanical().epoch(), epoch + 1);
            assert_eq!(live.mechanical().accepted_samples(), epoch + 1);
            assert!(
                frame.coupled.energy_balance_residual_j.abs() <= frame.coupled.energy_tolerance_j
            );
        }
    });
}

#[test]
fn g4_chemical_and_late_publication_refusals_leave_both_owners_retryable() {
    with_string_cx(|cx, gate| {
        let raw = law(RELEASE);
        let reaction = raw.clone().with_implicit_mechanism(policy()).unwrap();
        let (before, old_gas) = prepared(cx, &reaction);
        let exhausted = raw
            .clone()
            .with_implicit_mechanism(ImplicitMechanismSolve {
                maximum_evaluations: 1,
                ..policy()
            })
            .unwrap();
        let wrong_source = law(RELEASE + 1.0)
            .with_implicit_mechanism(policy())
            .unwrap();
        let tiny_consumption = ReactiveFiniteGasStringSolve {
            reaction: GasReactionSolve {
                maximum_consumed_fraction: 1e-8,
                ..controls().reaction
            },
            ..controls()
        };
        let late_limit = ReactiveFiniteGasStringSolve {
            maximum_gas_temperature_increment_k: 1e-8,
            ..controls()
        };
        for (law, limits) in [
            (&raw, controls()),
            (&exhausted, controls()),
            (&wrong_source, controls()),
            (&reaction, tiny_consumption),
            (&reaction, late_limit),
        ] {
            let mut live = before.clone();
            let mut gas = old_gas.clone();
            assert!(
                live.step_with_reacting_isochoric_gas(
                    cx,
                    8,
                    &[0.0; 2],
                    &mut gas,
                    FiniteGasStringHeat::default(),
                    law,
                    limits
                )
                .is_err()
            );
            assert_eq!(gas, old_gas);
            assert_eq!(live.mechanical().states(), before.mechanical().states());
            assert_eq!(live.mechanical().ambient(), before.mechanical().ambient());
            assert_eq!(live.thermal(), before.thermal());
            assert_eq!(live.mechanical().epoch(), 8);
            assert_eq!(live.mechanical().accepted_samples(), 8);
            let mut control = before.clone();
            let mut control_gas = old_gas.clone();
            let retry = live
                .step_with_reacting_isochoric_gas(
                    cx,
                    8,
                    &[0.0; 2],
                    &mut gas,
                    FiniteGasStringHeat::default(),
                    &reaction,
                    controls(),
                )
                .unwrap();
            let reference = control
                .step_with_reacting_isochoric_gas(
                    cx,
                    8,
                    &[0.0; 2],
                    &mut control_gas,
                    FiniteGasStringHeat::default(),
                    &reaction,
                    controls(),
                )
                .unwrap();
            assert_eq!(retry.coupled.string, reference.coupled.string);
            assert_eq!(retry.chemistry, reference.chemistry);
            assert_eq!(gas, control_gas);
        }
        let mut live = before.clone();
        let mut gas = old_gas.clone();
        gate.request();
        assert!(
            live.step_with_reacting_isochoric_gas(
                cx,
                8,
                &[0.0; 2],
                &mut gas,
                FiniteGasStringHeat::default(),
                &reaction,
                controls()
            )
            .is_err()
        );
        assert_eq!(gas, old_gas);
        assert_eq!(live.mechanical().states(), before.mechanical().states());
        assert_eq!(live.thermal(), before.thermal());
        assert_eq!(live.mechanical().epoch(), 8);
    });
}
