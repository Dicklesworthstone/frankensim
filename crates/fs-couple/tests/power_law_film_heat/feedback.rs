//! Actual finite-temperature feedback into the original nonlinear pressure law.

use super::*;
use fs_couple::render::plate::impact::material::heating::liquid_film::{
    LiquidFilmBinding, LiquidFilmImpactSystem, LiquidPowerLawFilmChannelBinding,
};
use fs_couple::render::plate::impact::material::heating::network::ReservoirCell;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};

const REFERENCE_MU: f64 = 0.01;

fn source_card(upper_t: f64, hot_mu: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density", 1000.0, Dims([-3, 1, 0, 0, 0, 0])),
        ("specific_heat_capacity", CP, Dims([2, 0, -2, -1, 0, 0])),
        (
            "dynamic_viscosity",
            REFERENCE_MU,
            Dims([-1, 1, -1, 0, 0, 0]),
        ),
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
                    source:
                        "manufactured power-law consistency temperature shift; not lubricant data"
                            .into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn source_law(card: &MaterialCard) -> LiquidTemperatureLaw<'_> {
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
                specific_enthalpy_j_kg: 20.0,
                temperature_k: 400.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
        ],
    )
    .unwrap()
}

fn bindings(
    card: &MaterialCard,
    n: f64,
    local: bool,
    limit: f64,
) -> Vec<LiquidPowerLawFilmChannelBinding<'_>> {
    [(1, if local { "b" } else { "a" }), (0, "a")]
        .into_iter()
        .map(|(channel, reservoir)| LiquidPowerLawFilmChannelBinding {
            film: 0,
            channel,
            reservoir: reservoir.into(),
            law: source_law(card),
            reference_consistency_pa_s_n: law(n).consistency_pa_s_n,
            maximum_relative_viscosity_change: limit,
        })
        .collect()
}

fn bind<'m>(
    inner: HeatedImpactSystem<'m>,
    rows: Vec<LiquidPowerLawFilmChannelBinding<'m>>,
) -> Result<LiquidFilmImpactSystem<'m>, ImpactError> {
    LiquidFilmImpactSystem::new_with_power_law(inner, vec![], vec![], vec![], vec![], rows)
}

fn shifted_laws<'s>(system: &'s LiquidFilmImpactSystem<'_>) -> &'s [RegularizedPowerLaw] {
    system
        .system()
        .mechanics()
        .power_law_film(0)
        .unwrap()
        .laws()
}

#[test]
fn thermal_consistency_changes_actual_thinning_newtonian_and_thickening_mechanics() {
    let card = source_card(300.2, 0.001);
    let c = source_curve(&card);
    let source = source_law(&card);
    for n in [0.5, 1.0, 2.0] {
        for analytic in [false, true] {
            let mut live = bind(
                heated(&c, n, vec![route(None, 1.0)], false).unwrap(),
                bindings(&card, n, false, 0.2),
            )
            .unwrap();
            let mut frozen = heated(&c, n, vec![route(None, 1.0)], false).unwrap();
            live.set_analytic_newton(analytic);
            frozen.set_analytic_newton(analytic);
            let initial = live.system().mechanics().stored_energy_j();
            let initial_momentum = momentum(live.system().mechanics().state());
            for tick in 1..=16 {
                let before = live.system().mechanics().state().to_vec();
                let old_law = shifted_laws(&live)[0];
                let frame = live.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
                frozen
                    .step(&[0.0; 2], &CancelGate::new_clock_free())
                    .unwrap();
                let after = live.system().mechanics().state();
                if tick == 1 {
                    assert_eq!(after, frozen.mechanics().state());
                }
                let pressure =
                    (relative(&before) - relative(after)) / (DT * (1.0 / M1 + 1.0 / M2) * AREA);
                let gap = GAP - 0.5 * (closure(&before) + closure(after));
                let (flow, power) = profile_work(pressure, gap, old_law);
                near(
                    flow,
                    AREA * 0.5 * (relative(&before) + relative(after)),
                    1e-12,
                );
                // Independent stress-profile integral uses the accepted OLD K.
                near(
                    frame.thermal.impact.deposited_heat_j,
                    DT * power,
                    2e-15 + 1e-7 * DT * power,
                );
                near(frame.thermal.impact.remaining_mechanical_loss_j, 0.0, 2e-15);
                let temperature = live.system().reservoirs()[0].state().temperature_k();
                let expected = source
                    .at_temperature(Temperature::new(temperature))
                    .unwrap();
                assert_eq!(
                    live.power_law_channel_sources(),
                    &[expected.clone(), expected.clone()]
                );
                let k =
                    law(n).consistency_pa_s_n * (expected.dynamic_viscosity_pa_s() / REFERENCE_MU);
                for shifted in shifted_laws(&live) {
                    assert_eq!(shifted.consistency_pa_s_n, k);
                    assert_eq!(shifted.flow_index, n);
                    assert_eq!(
                        shifted.minimum_shear_rate_per_s,
                        law(n).minimum_shear_rate_per_s
                    );
                }
                assert_eq!(live.system().mechanics().samples(), tick);
                near(momentum(after), initial_momentum, 1e-12);
            }
            assert!(shifted_laws(&live)[0].consistency_pa_s_n < law(n).consistency_pa_s_n);
            assert!(live.system().mechanics().state()[3] < frozen.mechanics().state()[3]);
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
}

#[test]
fn channel_temperatures_keep_original_flow_split_indices_and_reference_after_rebinding() {
    let card = source_card(300.2, 0.001);
    let c = source_curve(&card);
    let make = |reverse: bool| {
        let inner = heated(
            &c,
            0.5,
            vec![route(Some(0), 1.0), route(Some(1), 0.0)],
            true,
        )
        .unwrap();
        let mut rows = bindings(&card, 0.5, true, 0.2);
        if reverse {
            rows.reverse();
        }
        LiquidFilmImpactSystem::new_with_power_law(
            inner,
            vec![LiquidFilmBinding {
                film: 0,
                reservoir: "b".into(),
                law: source_law(&card),
                maximum_relative_viscosity_change: 0.2,
            }],
            vec![],
            vec![],
            vec![],
            rows,
        )
        .unwrap()
    };
    let mut live = make(false);
    let mut canonical = make(true);
    for _ in 0..12 {
        let before = live.system().mechanics().state().to_vec();
        let old_laws = shifted_laws(&live).to_vec();
        let old_heat: Vec<_> = live
            .system()
            .reservoirs()
            .iter()
            .map(CaloricReservoir::received_heat_j)
            .collect();
        live.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
        canonical
            .step(&[0.0; 2], &CancelGate::new_clock_free())
            .unwrap();
        let after = live.system().mechanics().state();
        let pressure = (relative(&before) - relative(after)) / (DT * (1.0 / M1 + 1.0 / M2) * AREA);
        let gap = GAP - 0.5 * (closure(&before) + closure(after));
        let mut total_flow = 0.0;
        for (channel, old_law) in old_laws.into_iter().enumerate() {
            let (flow, power) = profile_work(pressure, gap, old_law);
            let width_fraction = (channel + 1) as f64 / 3.0;
            total_flow += width_fraction * flow;
            let store = &live.system().reservoirs()[channel];
            near(
                store.received_heat_j() - old_heat[channel],
                DT * width_fraction * power,
                2e-15,
            );
            let expected = source_law(&card)
                .at_temperature(Temperature::new(store.state().temperature_k()))
                .unwrap();
            assert_eq!(live.power_law_channel_sources()[channel], expected);
        }
        near(
            total_flow,
            AREA * 0.5 * (relative(&before) + relative(after)),
            1e-12,
        );
        assert_eq!(
            live.system().mechanics().state(),
            canonical.system().mechanics().state()
        );
        assert_eq!(live.system().reservoirs(), canonical.system().reservoirs());
    }
    assert!(shifted_laws(&live)[1].consistency_pa_s_n < shifted_laws(&live)[0].consistency_pa_s_n);
    assert_eq!(live.sources().len(), 1);
    assert!(live.sliding_layer_sources().is_empty());
    assert_eq!(live.power_law_channel_sources().len(), 2);
    assert_eq!(live.system().mechanics().squeeze_film_count(), 1);
    let old_laws = shifted_laws(&live).to_vec();
    let (inner, uniform, channels, layers, sliding, power) = live.into_power_law_parts();
    let mut rebound = LiquidFilmImpactSystem::new_with_power_law(
        inner, uniform, channels, layers, sliding, power,
    )
    .unwrap();
    assert_eq!(shifted_laws(&rebound), old_laws);
    assert_eq!(rebound.system().mechanics().samples(), 12);
    rebound
        .step(&[0.0; 2], &CancelGate::new_clock_free())
        .unwrap();
    canonical
        .step(&[0.0; 2], &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(
        rebound.system().mechanics().state(),
        canonical.system().mechanics().state()
    );
    assert_eq!(
        rebound.system().reservoirs(),
        canonical.system().reservoirs()
    );
    assert_eq!(
        rebound.power_law_channel_sources(),
        canonical.power_law_channel_sources()
    );
}

#[test]
fn rejected_thermal_proposals_and_late_cancellation_preserve_sources_laws_and_retry() {
    let card = source_card(300.2, 0.001);
    let c = source_curve(&card);
    for cancel in [false, true] {
        let make = || {
            bind(
                heated(&c, 0.5, vec![route(None, 1.0)], false).unwrap(),
                bindings(&card, 0.5, false, 0.2),
            )
            .unwrap()
        };
        let mut live = make();
        let mut fresh = make();
        let motion = live.system().mechanics().state().to_vec();
        let thermal = live.system().reservoirs().to_vec();
        let sources = live.power_law_channel_sources().to_vec();
        let laws = shifted_laws(&live).to_vec();
        let gate = CancelGate::new_clock_free();
        let result = live.step_with_coupled_transport(
            &[0.0; 2],
            &gate,
            |_, _| Ok(None),
            |trial| {
                let transfers = trial
                    .reservoirs()
                    .iter()
                    .map(|row| {
                        let input = row.thermal;
                        // In the noncancelled case all new source laws compile, then
                        // the original enthalpy transaction rejects this invented heat.
                        let factor = if cancel { 1.0 } else { 1.5 };
                        UniformEnthalpyStep {
                            state: input
                                .curve
                                .advance_specific_energy(
                                    input.initial,
                                    factor * input.internal_heat_j / input.mass_kg,
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
        if cancel {
            assert!(matches!(result, Err(ImpactError::Cancelled)));
        }
        assert_eq!(live.system().mechanics().state(), motion);
        assert_eq!(live.system().mechanics().samples(), 0);
        assert_eq!(live.system().reservoirs(), thermal);
        assert_eq!(live.power_law_channel_sources(), sources);
        assert_eq!(shifted_laws(&live), laws);
        live.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
        fresh
            .step(&[0.0; 2], &CancelGate::new_clock_free())
            .unwrap();
        assert_eq!(
            live.system().mechanics().state(),
            fresh.system().mechanics().state()
        );
        assert_eq!(live.system().reservoirs(), fresh.system().reservoirs());
        assert_eq!(shifted_laws(&live), shifted_laws(&fresh));
        assert_eq!(
            live.power_law_channel_sources(),
            fresh.power_law_channel_sources()
        );
    }
}

#[test]
fn source_support_and_step_change_refuse_before_any_mechanical_or_thermal_publication() {
    for (upper_t, limit) in [(300.0 + 1e-8, 0.2), (300.2, 1e-14)] {
        let card = source_card(upper_t, 0.001);
        let c = source_curve(&card);
        let mut live = bind(
            heated(&c, 0.5, vec![route(None, 1.0)], false).unwrap(),
            bindings(&card, 0.5, false, limit),
        )
        .unwrap();
        let motion = live.system().mechanics().state().to_vec();
        let stores = live.system().reservoirs().to_vec();
        let sources = live.power_law_channel_sources().to_vec();
        let laws = shifted_laws(&live).to_vec();
        assert!(live.step(&[0.0; 2], &CancelGate::new_clock_free()).is_err());
        assert_eq!(live.system().mechanics().state(), motion);
        assert_eq!(live.system().mechanics().samples(), 0);
        assert_eq!(live.system().reservoirs(), stores);
        assert_eq!(live.power_law_channel_sources(), sources);
        assert_eq!(shifted_laws(&live), laws);
    }
}

#[test]
fn incomplete_duplicate_unknown_or_unjustified_reference_bindings_refuse() {
    let card = source_card(300.2, 0.001);
    let c = source_curve(&card);
    for invalid_case in 0..7 {
        let mut rows = bindings(&card, 0.5, true, 0.2);
        match invalid_case {
            0 => {
                rows.pop();
            }
            1 => {
                rows[0].channel = 0;
            }
            2 => {
                rows[0].film = 1;
            }
            3 => {
                rows[0].reservoir = "absent".into();
            }
            4 => {
                rows[0].reference_consistency_pa_s_n = 0.0;
            }
            5 => {
                rows[0].reference_consistency_pa_s_n = f64::INFINITY;
            }
            _ => {
                rows[0].maximum_relative_viscosity_change = 0.0;
            }
        }
        assert!(
            bind(
                heated(&c, 0.5, vec![route(None, 1.0)], false).unwrap(),
                rows
            )
            .is_err()
        );
    }
    // Equal numeric h(T) is insufficient: both stores and viscosity must retain
    // the same original material card identity.
    let foreign = source_card(300.2, 0.002);
    assert!(
        bind(
            heated(&c, 0.5, vec![route(None, 1.0)], false).unwrap(),
            bindings(&foreign, 0.5, false, 0.2)
        )
        .is_err()
    );
}

fn network(c: &EquilibriumEnthalpyPhaseCurve) -> DepositionSimulation<'_> {
    let mut events: Vec<_> = [1, 2]
        .into_iter()
        .map(|id| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: id,
                material_id: MaterialId(1),
                mass_kg: MASS,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [id as f64, 0.0, 0.0],
                    end_m: [id as f64 + 0.01, 0.0, 0.0],
                    volume_m3: 6e-7,
                    surface_area_m2: 0.0012,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: CellId(1),
            b: CellId(2),
            area_m2: 1e-6,
            conductance_w_k: 0.1,
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), c as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 8,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT / 4.0,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

#[test]
fn conduction_and_original_viscous_work_update_both_local_laws_on_one_clock() {
    let card = source_card(300.2, 0.001);
    let c = source_curve(&card);
    let mut live = bind(
        heated(&c, 0.5, vec![route(None, 1.0)], false).unwrap(),
        bindings(&card, 0.5, true, 0.2),
    )
    .unwrap();
    let mut thermal = network(&c);
    let reservoir_cells = [
        ReservoirCell {
            reservoir: "a".into(),
            cell: CellId(1),
        },
        ReservoirCell {
            reservoir: "b".into(),
            cell: CellId(2),
        },
    ];
    let initial = live.system().mechanics().stored_energy_j();
    for tick in 1..=8 {
        let frame = live
            .step_with_thermal_network(
                &mut thermal,
                &mut (),
                &reservoir_cells,
                StepConfig {
                    max_dt_s: DT / 4.0,
                    energy_tolerance_j: 1e-12,
                },
                &[0.0; 2],
                &CancelGate::new_clock_free(),
                |_, _| Ok(None),
            )
            .unwrap();
        assert!(frame.network.substeps >= 4);
        near(frame.external_heat_j, 0.0, 1e-18);
        near(
            frame
                .impact
                .thermal
                .transfers
                .iter()
                .map(|x| x.report.deposited_heat_j)
                .sum(),
            frame.impact.thermal.impact.deposited_heat_j,
            1e-18,
        );
        for (channel, binding) in reservoir_cells.iter().enumerate() {
            let store = &live.system().reservoirs()[channel];
            assert_eq!(
                store.state().specific_enthalpy_j_kg(),
                thermal.cells()[&binding.cell].specific_h_j_kg
            );
            let expected = source_law(&card)
                .at_temperature(Temperature::new(store.state().temperature_k()))
                .unwrap();
            assert_eq!(live.power_law_channel_sources()[channel], expected);
        }
        assert_eq!(live.system().mechanics().samples(), tick);
    }
    assert_eq!(live.system().reservoirs()[1].received_heat_j(), 0.0);
    assert!(live.system().reservoirs()[1].state().temperature_k() > 300.0);
    assert!(shifted_laws(&live)[1].consistency_pa_s_n < law(0.5).consistency_pa_s_n);
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
