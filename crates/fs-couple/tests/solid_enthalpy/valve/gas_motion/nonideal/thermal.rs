//! G1/G3/G4: finite source heat participates in the same compression endpoint.
//! Manufactured caloric/conductive data and original six-tet/PR fixtures.
use super::*;
use core::convert::Infallible;
use fs_couple::cavity_phs::modal::thermal::ModalGasThermalError;
use fs_couple::cavity_phs::piston::PistonGasState;
use fs_couple::cavity_phs::piston::thermal::layers::MaterialConductionLayers;
use fs_couple::cavity_phs::piston::thermal::{ConductanceSample, HeatExchangeSolve};
use fs_material::fluid_transport::FluidTransport;
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStepInput};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;
use fs_material::thermochemical_gas::nonideal::ResolvedThermochemicalRealGasState;

const BODY_MASS: f64 = 0.001;

fn exchange() -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-2.0, 2.0],
        source_tolerance_j: 1e-10,
        energy_tolerance_j: 1e-6,
        maximum_iterations: 96,
    }
}

fn wall(k: f64, cp: f64) -> (MaterialCard, PropertyKey, ResolvedHeatCapacityCurve) {
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let capacity = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let mut claims = ClaimSet::new();
    for (key, low, high) in [
        (conductivity.clone(), k, 2.0 * k),
        (capacity.clone(), cp, cp),
        (density.clone(), 1000.0, 1000.0),
    ] {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(250.0, low), (550.0, high)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        250.0,
                        550.0,
                    )
                    .with("lot", 7.0, 7.0),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured finite compression heat store".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal store".into(),
            phase: "authored caloric chart".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &capacity,
            density: &density,
            temperature_axis: "T",
            lower: &point(250.0),
            upper: &point(550.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-7,
    )
    .unwrap();
    (card, conductivity, caloric)
}

fn body(
    curve: &EquilibriumEnthalpyPhaseCurve,
    h: f64,
    heat_j: f64,
) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve,
        initial: curve.state_at_specific_enthalpy(h).unwrap(),
        mass_kg: BODY_MASS,
        volume_m3: BODY_MASS / 1000.0,
        surface_area_m2: 0.2,
        internal_heat_j: heat_j,
        duration_s: 1.0 / RATE as f64,
    }
}

fn layers<'a>(
    card: &'a MaterialCard,
    key: &'a PropertyKey,
    context: &'a QueryPoint,
) -> MaterialConductionLayers<'a> {
    MaterialConductionLayers {
        material: card,
        context,
        temperature_axis: "T",
        liquid_fraction_axis: None,
        conductivity: key,
        selection: MaterialPropertySelection::SingleClaimOnly,
        area_m2: 0.2,
        gas_distance_m: 0.001,
        body_distance_m: 0.02,
    }
}

#[test]
fn ideal_compression_matches_independent_two_storage_midpoint_elimination() {
    with_cx(|cx| {
        for (g, cp) in [(50.0, 50.0), (200.0, 50.0), (200.0, 100.0)] {
            let wall_motion = motion(1000.0, -0.04, 0.0, AREA, RATE, cx);
            let a = wall_motion.valve().valve().mass_normalized_pressure_load()[0];
            let state = wall_motion.model().states()[0];
            let (qq, qv, _, _, bq, _) = wall_motion.model().probe_linear_map(0);
            let free_sweep = a
                * ((qq - 1.0) * state.displacement_m_sqrt_kg + qv * state.velocity_m_sqrt_kg_per_s);
            let compliance = a * a * bq;
            let run = GasValveMotion::new(gas(0.75), wall_motion, cx).unwrap();
            let (_, _, caloric) = wall(1.0, cp);
            let input = body(caloric.curve(), cp * 100.0, 0.2);
            let exterior = EXTERIOR + 500.0;
            let next = run
                .thermal_candidate(
                    0.1,
                    exterior,
                    input,
                    controls(),
                    exchange(),
                    |_, _| {
                        Ok::<_, Infallible>(ConductanceSample {
                            conductance_w_k: g,
                            source: (),
                        })
                    },
                    cx,
                )
                .unwrap();

            // Constant NASA capacities give nR=P0 V0/T0 and Cv=(1.5+x) nR.
            // Eliminate body and gas midpoint temperatures analytically, then
            // solve only p_mid=P_EOS at their independently reconstructed state.
            let nr = EXTERIOR * VOLUME / 400.0;
            let gas_capacity = 2.25 * nr;
            let body_capacity = BODY_MASS * cp;
            let dg = g / RATE as f64;
            let sample = |p: f64| {
                let dv = free_sweep + compliance * (p - exterior);
                let heat_without_exchange = 0.1 - p * dv;
                let q = dg
                    * (350.0 - 400.0 + 0.2 / (2.0 * body_capacity)
                        - heat_without_exchange / (2.0 * gas_capacity))
                    / (1.0 + dg * (0.5 / body_capacity + 0.5 / gas_capacity));
                let tg = 400.0 + (heat_without_exchange + q) / gas_capacity;
                let tb = 350.0 + (0.2 - q) / body_capacity;
                let eos = nr * (400.0 + tg) / (2.0 * VOLUME + dv);
                (p - eos, q, tg, tb, VOLUME + dv)
            };
            let (mut lo, mut hi) = (0.5 * EXTERIOR, 1.5 * EXTERIOR);
            assert!(sample(lo).0 < 0.0 && sample(hi).0 > 0.0);
            for _ in 0..100 {
                let mid = 0.5 * lo + 0.5 * hi;
                if sample(mid).0 < 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let pressure = 0.5 * lo + 0.5 * hi;
            let (_, q, tg, tb, volume) = sample(pressure);
            let report = next.step();
            near(report.exchanged_heat_j, q, 2e-8);
            near(
                report.coupled.midpoint.temperature_k(),
                0.5 * (400.0 + tg),
                2e-7,
            );
            near(report.coupled.midpoint.pressure_pa(), pressure, 3e-4);
            near(next.motion().gas().state().temperature_k(), tg, 2e-7);
            near(next.thermal().state.temperature_k(), tb, 2e-7);
            near(
                report.thermal_midpoint.temperature_k(),
                0.5 * (350.0 + tb),
                2e-7,
            );
            near(next.motion().gas().volume_m3(), volume, 1e-14);
            assert!(report.coupled.motion.swept_volume_change_m3 < 0.0);
            assert!(report.exchanged_heat_j < 0.0); // Compressed gas heats the cooler finite wall.
            assert_eq!(next.thermal().external_heat_j, -report.exchanged_heat_j);
            assert!(report.source_residual_j.abs() <= exchange().source_tolerance_j);
            assert!(report.total_energy_residual_j.abs() <= exchange().energy_tolerance_j);
            assert!(report.coupled_trials > report.coupled.trials);
            assert_eq!(run.motion().accepted_steps(), 0);
            assert_eq!(next.motion().motion().accepted_steps(), 1);
        }
    });
}

#[test]
fn nonideal_source_layers_share_heat_pressure_work_and_continuing_modal_history() {
    with_cx(|cx| {
        let initial_gas = parcel();
        let initial_u = initial_gas.state().specific_internal_energy_j_kg();
        let original_mass = initial_gas.mass_kg();
        let wall_motion = motion(1000.0, -0.05, 0.1, AREA, RATE, cx);
        let initial_energy = energy(&wall_motion);
        let identity = wall_motion.valve().identity();
        let mut run = GasValveMotion::new(initial_gas.clone(), wall_motion, cx).unwrap();
        let (card, key, caloric) = wall(1.0, 100.0);
        let context = point(500.0);
        let mut input = body(caloric.curve(), 25_000.0, 0.0);
        let initial_h = input.initial.specific_enthalpy_j_kg();
        let (mut exterior_work, mut loss) = (0.0, 0.0);
        for tick in 1..=6 {
            let old_volume = run.gas().volume_m3();
            let next = run
                .thermal_layers_candidate(
                    25.0,
                    PRESSURE,
                    input,
                    solve(),
                    exchange(),
                    layers(&card, &key, &context),
                    cx,
                )
                .unwrap();
            let report = next.step();
            assert!(report.exchanged_heat_j > 0.0);
            assert!(next.thermal().state.temperature_k() < input.initial.temperature_k());
            let midpoint = &report.coupled.midpoint;
            let actual_k = midpoint
                .transport_state()
                .unwrap()
                .thermal_conductivity_w_m_k();
            let tb = report.thermal_midpoint.temperature_k();
            let body_k = 1.0 + (tb - 250.0) / 300.0;
            near(
                report.conductance_w_k,
                0.2 / (0.001 / actual_k + 0.02 / body_k),
                1e-12,
            );
            near(
                report.exchanged_heat_j,
                report.conductance_w_k * (tb - midpoint.temperature_k()) / RATE as f64,
                1e-10,
            );
            near(
                report.coupled.gas.work_by_gas_j,
                midpoint.pressure_pa() * (next.motion().gas().volume_m3() - old_volume),
                1e-10,
            );
            assert!(
                (midpoint.density_kg_m3() / midpoint.source().state().density - 1.0).abs() > 0.05
            );
            assert_eq!(
                midpoint.constitutive_law(),
                initial_gas.state().constitutive_law()
            );
            assert_eq!(next.thermal().report.card_identity(), card.content_hash());
            assert_eq!(next.thermal().report.query().axes().get("lot"), Some(&7.0));
            near(
                *next.thermal().report.query().axes().get("T").unwrap(),
                tb,
                0.0,
            );
            let original_source = layers(&card, &key, &context)
                .resolve(midpoint, report.thermal_midpoint)
                .unwrap();
            assert_eq!(next.thermal().report, original_source.source);
            assert_eq!(
                next.motion().gas().point().axes().get("specimen-lot"),
                Some(&42.0)
            );
            assert_eq!(next.motion().gas().model(), initial_gas.model());
            assert_eq!(next.motion().gas().mass_kg(), original_mass);
            assert_eq!(next.motion().motion().valve().identity(), identity);
            assert_eq!(next.motion().motion().accepted_steps(), tick);
            assert_eq!(run.motion().accepted_steps(), tick - 1);
            if tick == 1 {
                assert!(next.motion().gas().volume_m3() < old_volume);
            }
            exterior_work += report.coupled.exterior_pressure_work_j;
            loss += report.coupled.motion.modal.viscous_dissipation_j;
            let measured = original_mass
                * (next.motion().gas().state().specific_internal_energy_j_kg() - initial_u)
                + energy(next.motion().motion())
                - initial_energy
                + BODY_MASS * (next.thermal().state.specific_enthalpy_j_kg() - initial_h)
                + exterior_work
                + loss
                - 25.0 * tick as f64;
            near(measured, 0.0, 6.0 * exchange().energy_tolerance_j);
            assert!(report.body_energy_residual_j.abs() <= exchange().energy_tolerance_j);
            assert!(report.total_energy_residual_j.abs() <= exchange().energy_tolerance_j);
            let (accepted, thermal, _) = next.into_parts();
            run = accepted;
            input.initial = thermal.state;
        }
        assert!(loss > 0.0);
    });
}

#[test]
fn latent_storage_uses_enthalpy_midpoints_and_retains_the_original_phase_chart() {
    with_cx(|cx| {
        let run =
            GasValveMotion::new(parcel(), motion(1000.0, -0.01, 0.0, AREA, RATE, cx), cx).unwrap();
        let (card, _, _) = wall(1.0, 100.0);
        let knot = |h, t, f| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: 1000.0,
        };
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(
            card.content_hash(),
            vec![
                knot(0.0, 250.0, 0.0),
                knot(25_000.0, 500.0, 0.0),
                knot(35_000.0, 500.0, 1.0),
                knot(40_000.0, 550.0, 1.0),
            ],
        )
        .unwrap();
        let input = body(&curve, 30_000.0, 0.0);
        let next = run
            .thermal_candidate(
                0.0,
                PRESSURE,
                input,
                solve(),
                exchange(),
                |_, phase| {
                    Ok::<_, Infallible>(ConductanceSample {
                        conductance_w_k: 20.0,
                        source: phase,
                    })
                },
                cx,
            )
            .unwrap();
        assert_eq!(next.thermal().state.temperature_k(), 500.0);
        assert_eq!(next.step().thermal_midpoint.temperature_k(), 500.0);
        assert!(next.thermal().state.liquid_mass_fraction() < 0.5);
        near(
            next.step().thermal_midpoint.liquid_mass_fraction(),
            0.5 * (input.initial.liquid_mass_fraction()
                + next.thermal().state.liquid_mass_fraction()),
            1e-14,
        );
        near(
            next.thermal().state.specific_enthalpy_j_kg(),
            30_000.0 - next.step().exchanged_heat_j / BODY_MASS,
            1e-10,
        );
        assert_eq!(next.thermal().report, next.step().thermal_midpoint);
        assert_eq!(input.initial.liquid_mass_fraction(), 0.5);
        assert!(next.step().total_energy_residual_j.abs() <= exchange().energy_tolerance_j);
    });
}

#[test]
fn explicit_insulation_reduces_exactly_to_the_native_nonideal_motion_step() {
    with_cx(|cx| {
        let run =
            GasValveMotion::new(parcel(), motion(1000.0, 0.01, 0.1, AREA, RATE, cx), cx).unwrap();
        let (_, _, caloric) = wall(1.0, 100.0);
        let input = body(caloric.curve(), 10_000.0, 0.0);
        let mut reference = run.clone();
        let native = reference.step(10.0, PRESSURE, solve(), cx).unwrap();
        let next = run
            .thermal_candidate(
                10.0,
                PRESSURE,
                input,
                solve(),
                exchange(),
                |_, _| {
                    Ok::<_, Infallible>(ConductanceSample {
                        conductance_w_k: 0.0,
                        source: (),
                    })
                },
                cx,
            )
            .unwrap();
        assert_eq!(next.motion().gas(), reference.gas());
        assert_eq!(
            next.motion().motion().model().states(),
            reference.motion().model().states()
        );
        assert_eq!(next.step().coupled.gas, native.gas);
        assert_eq!(next.step().coupled.midpoint, native.midpoint);
        assert_eq!(next.thermal().state, input.initial);
        assert_eq!(next.step().exchanged_heat_j, 0.0);
        assert_eq!(next.step().exchange_trials, 1);
        assert_eq!(next.step().coupled_trials, native.trials);
        assert_eq!(run.motion().accepted_steps(), 0);
    });
}

fn constant_conductance(
    _: &ResolvedThermochemicalRealGasState,
    _: EquilibriumPhaseState,
) -> Result<ConductanceSample<()>, &'static str> {
    Ok(ConductanceSample {
        conductance_w_k: 5.0,
        source: (),
    })
}

#[test]
fn late_interface_domain_clock_and_budget_failures_preserve_all_accepted_owners() {
    with_cx(|cx| {
        let mut run =
            GasValveMotion::new(parcel(), motion(1000.0, 0.01, 0.1, AREA, RATE, cx), cx).unwrap();
        run.step(10.0, PRESSURE, solve(), cx).unwrap();
        let before = run.clone();
        let (_, _, caloric) = wall(1.0, 100.0);
        let input = body(caloric.curve(), 10_000.0, 0.0);
        let resolve = constant_conductance;
        assert!(matches!(
            run.thermal_candidate(
                0.0,
                PRESSURE,
                UniformEnthalpyStepInput {
                    duration_s: input.duration_s * 2.0,
                    ..input
                },
                solve(),
                exchange(),
                resolve,
                cx
            ),
            Err(ModalGasThermalError::Invalid(_))
        ));
        let mut calls = 0;
        assert!(matches!(
            run.thermal_candidate(
                0.0,
                PRESSURE,
                input,
                solve(),
                exchange(),
                |_, _| {
                    calls += 1;
                    if calls == 2 {
                        Err("late source refusal")
                    } else {
                        Ok(ConductanceSample {
                            conductance_w_k: 5.0,
                            source: (),
                        })
                    }
                },
                cx
            ),
            Err(ModalGasThermalError::Interface("late source refusal"))
        ));
        assert_eq!(calls, 2);
        assert!(matches!(
            run.thermal_candidate(
                0.0,
                PRESSURE,
                input,
                solve(),
                HeatExchangeSolve {
                    heat_bracket_j: [-1e-6, 1e-6],
                    ..exchange()
                },
                resolve,
                cx
            ),
            Err(ModalGasThermalError::OutsideBracket { .. })
        ));
        assert!(matches!(
            run.thermal_candidate(
                0.0,
                PRESSURE,
                input,
                ModalGasSolve {
                    maximum_trials: 1,
                    ..solve()
                },
                exchange(),
                resolve,
                cx
            ),
            Err(ModalGasThermalError::Motion(ModalGasError::TrialBudget))
        ));
        assert!(matches!(
            run.thermal_candidate(
                0.0,
                PRESSURE,
                input,
                solve(),
                HeatExchangeSolve {
                    heat_bracket_j: [-100.0, 100.0],
                    ..exchange()
                },
                resolve,
                cx
            ),
            Err(ModalGasThermalError::Thermal(_))
        ));
        let next = run
            .thermal_candidate(0.0, PRESSURE, input, solve(), exchange(), resolve, cx)
            .unwrap();
        let repeated = before
            .thermal_candidate(0.0, PRESSURE, input, solve(), exchange(), resolve, cx)
            .unwrap();
        assert_eq!(next.motion().gas(), repeated.motion().gas());
        assert_eq!(
            next.motion().motion().model().states(),
            repeated.motion().motion().model().states()
        );
        assert_eq!(next.thermal(), repeated.thermal());
        assert_eq!(run.gas(), before.gas());
        assert_eq!(
            run.motion().model().states(),
            before.motion().model().states()
        );
        assert_eq!(
            run.motion().accepted_steps(),
            before.motion().accepted_steps()
        );
        near(input.initial.temperature_k(), 350.0, 1e-12);
    });
}

#[test]
fn cancellation_after_actual_interface_evaluation_publishes_no_partial_heat_or_motion() {
    with_cx(|construction_cx| {
        let run = GasValveMotion::new(
            parcel(),
            motion(1000.0, -0.01, 0.0, AREA, RATE, construction_cx),
            construction_cx,
        )
        .unwrap();
        let before = run.clone();
        let (_, _, caloric) = wall(1.0, 100.0);
        let input = body(caloric.curve(), 10_000.0, 0.0);
        let gate = CancelGate::new_clock_free();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            let mut calls = 0;
            let result = run.thermal_candidate(
                0.0,
                PRESSURE,
                input,
                solve(),
                exchange(),
                |gas, _| {
                    assert!(gas.pressure_pa().is_finite());
                    calls += 1;
                    gate.request();
                    Ok::<_, Infallible>(ConductanceSample {
                        conductance_w_k: 5.0,
                        source: (),
                    })
                },
                &cx,
            );
            assert_eq!(calls, 1);
            assert!(matches!(result, Err(ModalGasThermalError::Cancelled)));
        });
        assert_eq!(run.gas(), before.gas());
        assert_eq!(
            run.motion().model().states(),
            before.motion().model().states()
        );
        assert_eq!(run.motion().accepted_steps(), 0);
        let clean = run
            .thermal_candidate(
                0.0,
                PRESSURE,
                input,
                solve(),
                exchange(),
                |_, _| {
                    Ok::<_, Infallible>(ConductanceSample {
                        conductance_w_k: 5.0,
                        source: (),
                    })
                },
                construction_cx,
            )
            .unwrap();
        assert_eq!(clean.motion().motion().accepted_steps(), 1);
        assert!(clean.step().exchanged_heat_j > 0.0);
    });
}
