//! Source temperature changes the original sliding mechanics on its next tick.

use super::*;
use fs_couple::render::plate::impact::material::heating::liquid_film::{
    LiquidFilmBinding, LiquidFilmImpactSystem, LiquidSlidingFilmLayerBinding,
};
use fs_couple::render::plate::impact::material::heating::transport::CoupledThermalStep;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::UniformEnthalpyStep;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};

const INTERFACE: f64 = 0.25;
const CP: f64 = 0.2;

fn source_card(hot_mu: f64, upper_t: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density", 1000.0, Dims([-3, 1, 0, 0, 0, 0])),
        ("specific_heat_capacity", CP, Dims([2, 0, -2, -1, 0, 0])),
        ("dynamic_viscosity", 0.01, Dims([-1, 1, -1, 0, 0, 0])),
        ("thermal_conductivity", 0.6, Dims([1, 1, -3, -1, 0, 0])),
    ] {
        let quantity = if name == "thermal_conductivity" {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        let (value, interpolation) = if name == "dynamic_viscosity" {
            (
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, value), (upper_t, hot_mu)],
                    dims,
                },
                InterpolationPolicy::LinearInside,
            )
        } else {
            (
                PropertyValue::Scalar { value, dims },
                InterpolationPolicy::ConstantWithinValidity,
            )
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value,
                interpolation,
                validity: ValidityDomain::unconstrained()
                    .with("T", 300.0, upper_t)
                    .with("p", 100000.0, 100000.0),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic sliding-film thermal feedback oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic lubricant".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn law(card: &MaterialCard) -> LiquidTemperatureLaw<'_> {
    LiquidTemperatureLaw::try_new(
        card,
        &QueryPoint::new()
            .with("T", 300.0)
            .unwrap()
            .with("p", 100000.0)
            .unwrap(),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn source_curve(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Liquid,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: CP * 100.0,
                temperature_k: 400.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
        ],
    )
    .unwrap()
}

fn source_inner(
    c: &EquilibriumEnthalpyPhaseCurve,
    layered: bool,
    route_count: usize,
) -> HeatedImpactSystem<'_> {
    let graph = graph().with_barus_viscosity(2e-5).unwrap();
    let rows = [1.0, -1.0]
        .into_iter()
        .map(|sign| WallVelocityPorts {
            lower: vec![sign / M1.sqrt(), 0.0],
            upper: vec![0.0, sign / M2.sqrt()],
        })
        .collect();
    let sliding = if layered {
        let profile = GapViscosityProfile::new(vec![
            ViscosityLayer {
                upper_fraction: INTERFACE,
                dynamic_viscosity_pa_s: 0.01,
            },
            ViscosityLayer {
                upper_fraction: 1.0,
                dynamic_viscosity_pa_s: 0.01,
            },
        ])
        .unwrap();
        ShearedFilm::from_stratified(StratifiedFilm::new(&graph, vec![profile; 2]).unwrap(), rows)
    } else {
        ShearedFilm::new(graph.clone(), rows)
    }
    .unwrap();
    // Normal and sliding islands deliberately share index zero. The fixed
    // normal graph has no swept volume and contributes no mechanical loss.
    let mechanical = bare()
        .with_squeeze_film(graph)
        .unwrap()
        .with_sliding_film(sliding)
        .unwrap();
    let routes = (0..route_count)
        .map(|channel| SlidingFilmLayerHeatRoute {
            film: 0,
            channel,
            layer_destinations: if layered {
                vec!["lower".into(), "upper".into()]
            } else {
                vec![if channel == 0 { "lower" } else { "upper" }.into()]
            },
        })
        .collect();
    let mut stores = stores(c);
    if !layered {
        stores[1] = CaloricReservoir::new(
            "upper".into(),
            c,
            c.state_at_specific_enthalpy(CP * 0.01).unwrap(),
            STORE_MASS,
            STORE_VOLUME,
            STORE_AREA,
        )
        .unwrap();
    }
    HeatedImpactSystem::new_with_sliding_losses(
        mechanical.prepare_analytic().unwrap(),
        stores,
        vec![],
        vec![],
        vec![],
        vec![],
        routes,
        vec![],
        1e-12,
    )
    .unwrap()
}

fn source_bindings(
    card: &MaterialCard,
    layered: bool,
    limit: f64,
) -> Vec<LiquidSlidingFilmLayerBinding<'_>> {
    let mut result = Vec::new();
    for channel in 0..2 {
        for layer in 0..if layered { 2 } else { 1 } {
            let lower = if layered { layer == 0 } else { channel == 0 };
            let reservoir = if lower { "lower" } else { "upper" };
            result.push(LiquidSlidingFilmLayerBinding {
                film: 0,
                channel,
                layer,
                reservoir: reservoir.into(),
                law: law(card),
                maximum_relative_viscosity_change: limit,
            });
        }
    }
    result.reverse(); // authored order must not alter physical layer identity
    result
}

fn source_system<'m>(
    card: &'m MaterialCard,
    c: &'m EquilibriumEnthalpyPhaseCurve,
    layered: bool,
    limit: f64,
) -> LiquidFilmImpactSystem<'m> {
    LiquidFilmImpactSystem::new_with_sliding(
        source_inner(c, layered, 2),
        vec![LiquidFilmBinding {
            film: 0,
            reservoir: "lower".into(),
            law: law(card),
            maximum_relative_viscosity_change: 0.1,
        }],
        vec![],
        vec![],
        source_bindings(card, layered, limit),
    )
    .unwrap()
}

#[test]
fn source_temperatures_update_actual_drag_and_preserve_the_completed_layers_work() {
    // G1/G3: solve the two free masses' midpoint drag independently. Unequal
    // layer thicknesses heat the equal stores differently, changing both the
    // harmonic drag and through-gap work fractions on subsequent ticks.
    let card = source_card(0.0025, 300.1);
    let c = source_curve(&card);
    let law = law(&card);
    for analytic in [false, true] {
        let mut live = source_system(&card, &c, true, 0.1);
        let mut frozen = source_inner(&c, true, 2);
        live.set_analytic_newton(analytic);
        frozen.set_analytic_newton(analytic);
        let initial = live.system().mechanics().stored_energy_j();
        let initial_momentum = momentum(live.system().mechanics().state());
        assert_eq!(
            live.sliding_layer_bindings()
                .iter()
                .map(|row| (row.channel, row.layer))
                .collect::<Vec<_>>(),
            [(0, 0), (0, 1), (1, 0), (1, 1)]
        );
        assert_eq!(live.sources().len(), 1);
        assert!(live.channel_sources().is_empty() && live.layer_sources().is_empty());
        for tick in 1..=64 {
            let before = relative(live.system().mechanics().state());
            let mu = [
                live.sliding_layer_sources()[0].dynamic_viscosity_pa_s(),
                live.sliding_layer_sources()[1].dynamic_viscosity_pa_s(),
            ];
            let fluidity = INTERFACE / mu[0] + (1.0 - INTERFACE) / mu[1];
            let drag = 2.0 * WIDTH * LENGTH / (GAP * fluidity);
            let rate = drag * (1.0 / M1 + 1.0 / M2);
            let expected = before * (1.0 - 0.5 * DT * rate) / (1.0 + 0.5 * DT * rate);
            let received: Vec<_> = live
                .system()
                .reservoirs()
                .iter()
                .map(CaloricReservoir::received_heat_j)
                .collect();
            let frame = live.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
            frozen
                .step(&[0.0; 2], &CancelGate::new_clock_free())
                .unwrap();
            let after = relative(live.system().mechanics().state());
            let heat = DT * drag * (0.5 * (before + after)).powi(2);
            near(after, expected, 2e-11);
            near(frame.thermal.impact.deposited_heat_j, heat, 2e-18);
            assert_eq!(frame.thermal.impact.mechanical.sample, tick);
            let fractions = [
                INTERFACE / mu[0] / fluidity,
                (1.0 - INTERFACE) / mu[1] / fluidity,
            ];
            for layer in 0..2 {
                let store = &live.system().reservoirs()[layer];
                near(
                    store.received_heat_j() - received[layer],
                    fractions[layer] * heat,
                    2e-18,
                );
                let expected = law
                    .at_temperature(Temperature::new(store.state().temperature_k()))
                    .unwrap();
                for channel in 0..2 {
                    assert_eq!(live.sliding_layer_sources()[2 * channel + layer], expected);
                    let film = live.system().mechanics().sliding_film(0).unwrap();
                    assert_eq!(
                        film.stratified_profiles().unwrap()[channel].layers()[layer]
                            .dynamic_viscosity_pa_s,
                        expected.dynamic_viscosity_pa_s()
                    );
                    assert_eq!(
                        film.pressure_operator()
                            .pressure_viscosity_coefficient_per_pa(),
                        2e-5
                    );
                }
            }
            assert_eq!(
                live.system()
                    .mechanics()
                    .squeeze_film(0)
                    .unwrap()
                    .viscosity_pa_s(),
                Some(live.sources()[0].dynamic_viscosity_pa_s())
            );
        }
        assert_ne!(
            live.sliding_layer_sources()[0].dynamic_viscosity_pa_s(),
            live.sliding_layer_sources()[1].dynamic_viscosity_pa_s()
        );
        assert_ne!(
            live.system().mechanics().state(),
            frozen.mechanics().state()
        );
        near(
            momentum(live.system().mechanics().state()),
            initial_momentum,
            2e-11,
        );
        near(
            live.system().mechanics().stored_energy_j()
                + live
                    .system()
                    .reservoirs()
                    .iter()
                    .map(CaloricReservoir::storage_change_j)
                    .sum::<f64>(),
            initial,
            2e-11,
        );
    }
}

#[test]
fn uniform_channels_keep_their_own_source_temperatures_and_physical_heat() {
    let card = source_card(0.0025, 300.1);
    let c = source_curve(&card);
    let law = law(&card);
    let mut live = source_system(&card, &c, false, 0.1);
    assert_eq!(live.sliding_layer_sources().len(), 2);
    assert_ne!(
        live.sliding_layer_sources()[0].dynamic_viscosity_pa_s(),
        live.sliding_layer_sources()[1].dynamic_viscosity_pa_s()
    );
    for _ in 0..32 {
        let before = relative(live.system().mechanics().state());
        let mu: Vec<_> = live
            .sliding_layer_sources()
            .iter()
            .map(|source| source.dynamic_viscosity_pa_s())
            .collect();
        let received: Vec<_> = live
            .system()
            .reservoirs()
            .iter()
            .map(CaloricReservoir::received_heat_j)
            .collect();
        let frame = live.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
        let speed = 0.5 * (before + relative(live.system().mechanics().state()));
        let film = live.system().mechanics().sliding_film(0).unwrap();
        assert!(film.stratified_profiles().is_none());
        let mut total = 0.0;
        for channel in 0..2 {
            let heat = DT * WIDTH * LENGTH * mu[channel] * speed * speed / GAP;
            total += heat;
            let store = &live.system().reservoirs()[channel];
            near(store.received_heat_j() - received[channel], heat, 2e-18);
            let source = law
                .at_temperature(Temperature::new(store.state().temperature_k()))
                .unwrap();
            assert_eq!(live.sliding_layer_sources()[channel], source);
            assert_eq!(
                film.pressure_operator().channel_viscosities_pa_s()[channel],
                source.dynamic_viscosity_pa_s()
            );
        }
        near(frame.thermal.impact.deposited_heat_j, total, 2e-18);
    }
}

#[test]
fn constant_source_keeps_uniform_and_retained_layer_trajectories_bitwise() {
    let card = source_card(0.01, 300.1);
    let c = source_curve(&card);
    for layered in [false, true] {
        let mut live = source_system(&card, &c, layered, 0.1);
        let mut frozen = source_inner(&c, layered, 2);
        for _ in 0..16 {
            let live_frame = live.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
            let frozen_frame = frozen
                .step(&[0.0; 2], &CancelGate::new_clock_free())
                .unwrap();
            assert_eq!(
                live.system().mechanics().state(),
                frozen.mechanics().state()
            );
            assert_eq!(live.system().reservoirs(), frozen.reservoirs());
            assert_eq!(
                live_frame.thermal.impact.deposited_heat_j,
                frozen_frame.deposited_heat_j
            );
        }
        let (inner, uniform, channels, layers, sliding) = live.into_sliding_parts();
        let rebound =
            LiquidFilmImpactSystem::new_with_sliding(inner, uniform, channels, layers, sliding)
                .unwrap();
        assert_eq!(
            rebound.system().mechanics().state(),
            frozen.mechanics().state()
        );
        assert_eq!(rebound.system().reservoirs(), frozen.reservoirs());
    }
}

#[test]
fn source_domain_step_change_and_late_cancellation_leave_all_accepted_owners_unchanged() {
    for (upper_t, limit, cancel) in [
        (300.0 + 1e-8, 0.1, false),
        (300.1, 1e-12, false),
        (300.1, 0.1, true),
    ] {
        let card = source_card(0.0025, upper_t);
        let c = source_curve(&card);
        let mut live = source_system(&card, &c, true, limit);
        let motion = live.system().mechanics().state().to_vec();
        let stores = live.system().reservoirs().to_vec();
        let sources = live.sliding_layer_sources().to_vec();
        let normal = live.sources().to_vec();
        let profiles = live
            .system()
            .mechanics()
            .sliding_film(0)
            .unwrap()
            .stratified_profiles()
            .unwrap()
            .to_vec();
        let gate = CancelGate::new_clock_free();
        let result = live.step_with_coupled_transport(
            &[0.0; 2],
            &gate,
            |_, _| Ok(None),
            |trial| {
                let transfers = trial
                    .reservoirs()
                    .iter()
                    .map(|input| {
                        let input = input.thermal;
                        UniformEnthalpyStep {
                            state: input
                                .curve
                                .advance_specific_energy(
                                    input.initial,
                                    input.internal_heat_j / input.mass_kg,
                                )
                                .unwrap(),
                            external_heat_j: 0.0,
                            energy_residual_tolerance_j: 1e-12,
                            report: (),
                        }
                    })
                    .collect();
                if cancel {
                    gate.request();
                }
                Ok(CoupledThermalStep {
                    transfers,
                    boundary_heat_j: 0.0,
                })
            },
        );
        assert!(result.is_err());
        assert_eq!(live.system().mechanics().state(), motion);
        assert_eq!(live.system().reservoirs(), stores);
        assert_eq!(live.sliding_layer_sources(), sources);
        assert_eq!(live.sources(), normal);
        assert_eq!(
            live.system()
                .mechanics()
                .sliding_film(0)
                .unwrap()
                .stratified_profiles()
                .unwrap(),
            profiles
        );
        if cancel {
            let retry = live.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
            assert_eq!(retry.thermal.impact.mechanical.sample, 1);
        }
    }
}

#[test]
fn sliding_bindings_require_complete_original_layers_and_their_existing_heat_destinations() {
    let card = source_card(0.0025, 300.1);
    let c = source_curve(&card);
    for case in 0..5 {
        let mut bindings = source_bindings(&card, true, 0.1);
        match case {
            0 => {
                bindings.pop();
            }
            1 => {
                bindings[1] = bindings[0].clone();
            }
            2 => {
                bindings[0].reservoir = "lower".into();
            }
            3 => {
                bindings[0].film = 7;
            }
            _ => {}
        }
        let inner = source_inner(&c, true, if case == 4 { 1 } else { 2 });
        assert!(
            LiquidFilmImpactSystem::new_with_sliding(inner, vec![], vec![], vec![], bindings)
                .is_err()
        );
    }
}
