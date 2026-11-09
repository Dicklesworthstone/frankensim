//! Shared finite-gas energy with distinct native sensible/latent bodies.
use super::*;
use fs_couple::finite_gas_thermal::{ConvectionCoefficient,
    bodies::{FiniteGasBodiesEnvironment, FiniteGasBodiesError, FiniteGasBody}};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};

fn input(source: &Sources, temperature: f64, heat: f64) -> UniformEnthalpyStepInput<'_> {
    let old = Oscillator::try_new(&source.mechanical, &source.caloric,
        geometry(0.0), temperature, 0.0, 0.0).unwrap();
    UniformEnthalpyStepInput { curve: source.caloric.curve(), initial: *old.point().caloric_state(),
        mass_kg: old.point().mass_kg(), volume_m3: old.point().volume_m3(),
        surface_area_m2: old.geometry().thermal_surface_m2, internal_heat_j: heat, duration_s: 0.01 }
}
fn shared(gas: &IsochoricGasParcel, heat: f64) -> FiniteGasBodiesEnvironment<'_> {
    FiniteGasBodiesEnvironment { gas, gas_heat_j: heat, gas_solve: gas_solve(),
        exchange: exchange(0.02), maximum_gas_temperature_increment_k: 20.0, maximum_bodies: 8 }
}

#[test]
fn distinct_capacities_share_one_analytic_gas_endpoint_independent_of_input_order() {
    scope(|cx, _| {
        let sources = [sources(2.0, 1000.0, false, false), sources(4.0, 1000.0, false, false),
            sources(6.0, 1000.0, false, false)];
        let boundaries = [boundary(&sources[0], 100.0), boundary(&sources[1], 50.0), boundary(&sources[2], 75.0)];
        let inputs = [input(&sources[0], 320.0, 2e-4), input(&sources[1], 290.0, 1e-4), input(&sources[2], 310.0, -2e-4)];
        let bodies: Vec<_> = (0..3).map(|i| FiniteGasBody {
            id: [30, 10, 20][i], environment: &boundaries[i], input: inputs[i],
        }).collect();
        let gas = gas(300.0, 1e-6);
        let saved = gas.clone();
        let environment = shared(&gas, 1e-3);
        let result = environment.advance_with_convection(cx, &bodies, |body, gas| {
            Ok::<_, core::convert::Infallible>(ConvectionCoefficient {
                coefficient_w_per_m2_k: body.environment.convection_w_per_m2_k,
                source: gas.state().temperature,
            })
        }).unwrap();
        let capacities = [2.0, 4.0, 6.0].map(|cp| cp * inputs[0].mass_kg);
        let beta: Vec<_> = (0..3).map(|i| {
            let a = inputs[i].duration_s * inputs[i].surface_area_m2 * boundaries[i].convection_w_per_m2_k;
            a * capacities[i] / (a + capacities[i])
        }).collect();
        let isolated = [320.0, 290.0, 310.0];
        let cg = 2.5 * 100000.0 * gas.volume_m3() / 300.0;
        let expected_gas = (cg * 300.0 + environment.gas_heat_j
            + (0..3).map(|i| beta[i] * (isolated[i] + inputs[i].internal_heat_j / capacities[i])).sum::<f64>())
            / (cg + beta.iter().sum::<f64>());
        close(result.gas().state().state().temperature, expected_gas, 2e-8);
        for next in result.bodies() {
            let i = [30, 10, 20].iter().position(|id| *id == next.report.id).unwrap();
            let source_t = isolated[i] + inputs[i].internal_heat_j / capacities[i];
            let q = beta[i] * (source_t - expected_gas);
            close(next.report.solid_to_gas_heat_j, q, 2e-11);
            close(next.state.temperature_k(), source_t - q / capacities[i], 2e-8);
            assert_eq!(next.report.convection.source, result.gas().state().state().temperature);
            assert_eq!(next.external_heat_j, -next.report.solid_to_gas_heat_j);
        }
        assert_eq!(result.bodies().iter().map(|b| b.report.id).collect::<Vec<_>>(), vec![10, 20, 30]);
        assert_eq!(result.report().body_trials, 3 * result.report().exchange_trials);
        assert!(result.report().exchange_trials <= environment.exchange.maximum_iterations + 2);
        assert!(result.report().absolute_energy_defect_j <= environment.exchange.energy_tolerance_j);
        assert!(result.report().entropy.balance_residual_j_k.abs() < 1e-15);
        let mut reversed = bodies.clone(); reversed.reverse();
        let a = environment.advance(cx, &bodies).unwrap();
        let b = environment.advance(cx, &reversed).unwrap();
        assert_eq!(a, b);
        // Native sequential reservoir updates impose a different physical split.
        let sequential = |order: &[FiniteGasBody<'_>]| {
            let mut current = gas.clone();
            for body in order {
                let next = FiniteGasThermalEnvironment { gas: &current, body: body.environment,
                    gas_heat_j: 0.0, gas_solve: gas_solve(), exchange: exchange(0.02),
                    maximum_gas_temperature_increment_k: 20.0 }.advance(cx, body.input).unwrap();
                current = next.report.next_gas;
            }
            current.state().state().temperature
        };
        assert!((sequential(&bodies) - sequential(&reversed)).abs() > 0.1);
        assert_eq!(gas, saved);
    });
}

#[test]
fn zero_net_gas_heat_still_freezes_one_latent_body_and_melts_the_other() {
    scope(|cx, _| {
        let source = sources(2.0, 1000.0, false, false);
        let curves: Vec<_> = [310.0, 290.0].map(|t| EquilibriumEnthalpyPhaseCurve::try_new(
            source.caloric.curve().material_card_identity(),
            [(2.0 * (280.0 - t), 280.0, 0.0), (0.0, t, 0.0),
                (100.0, t, 1.0), (100.0 + 2.0 * (340.0 - t), 340.0, 1.0)]
                .map(|(h, t, liquid)| EnthalpyPhaseKnot { specific_enthalpy_j_kg: h,
                    temperature_k: t, liquid_mass_fraction: liquid, bulk_density_kg_m3: 1000.0 }).to_vec(),
        ).unwrap()).to_vec();
        let boundary = boundary(&source, 100.0);
        let bodies: Vec<_> = curves.iter().enumerate().map(|(i, curve)| FiniteGasBody {
            id: i as u64, environment: &boundary,
            input: UniformEnthalpyStepInput { curve, initial: curve.state_at_specific_enthalpy(50.0).unwrap(),
                mass_kg: 0.001, volume_m3: 1e-6, surface_area_m2: 0.0013,
                internal_heat_j: 0.0, duration_s: 0.01 },
        }).collect();
        let gas = gas(300.0, 1e-6);
        // Native exact idle must not query this unused unsupported endpoint.
        let environment = FiniteGasBodiesEnvironment { gas_solve: GasEnergySolve {
            maximum_temperature_k: 1000.0, ..gas_solve()
        }, maximum_gas_temperature_increment_k: 1000.0, ..shared(&gas, 0.0) };
        let result = environment.advance(cx, &bodies).unwrap();
        assert_eq!(result.gas(), &gas);
        assert_eq!(result.report().exchange_trials, 1);
        assert_eq!(result.report().body_trials, 2);
        assert_eq!(result.report().heat_into_gas_j, 0.0);
        for (i, q) in [0.013, -0.013].into_iter().enumerate() {
            let body = &result.bodies()[i];
            close(body.report.solid_to_gas_heat_j, q, 1e-12);
            assert_eq!(body.state.temperature_k(), [310.0, 290.0][i]);
            close(body.state.liquid_mass_fraction(), 0.5 - q / 0.1, 1e-10);
            close(body.report.entropy.storage_change_j_k, -q / body.state.temperature_k(), 1e-14);
        }
        let entropy = result.report().entropy;
        close(entropy.storage_change_j_k[0], 0.013 * (1.0 / 290.0 - 1.0 / 310.0), 1e-14);
        assert_eq!(entropy.storage_change_j_k[1], 0.0);
        assert!(entropy.convection_production_j_k > 0.0);
        close(entropy.temporal_defect_j_k[0], 0.0, 1e-14);
        close(entropy.balance_residual_j_k, 0.0, 1e-15);
    });
}

#[test]
fn gas_sources_are_counted_once_and_each_radiation_enclosure_stays_external() {
    scope(|cx, _| {
        let a = sources(2.0, 1000.0, false, false);
        let b = sources(4.0, 1000.0, false, false);
        let boundaries = [LumpedThermalEnvironment { radiation_temperature_k: 400.0,
            transport: LumpedThermalTransport::try_declared(10.0, 0.4).unwrap(), ..boundary(&a, 0.0) },
            LumpedThermalEnvironment { radiation_temperature_k: 250.0,
                transport: LumpedThermalTransport::try_declared(10.0, 0.6).unwrap(), ..boundary(&b, 0.0) }];
        let bodies = [FiniteGasBody { id: 1, environment: &boundaries[0], input: input(&a, 300.0, 0.001) },
            FiniteGasBody { id: 2, environment: &boundaries[1], input: input(&b, 300.0, -0.002) }];
        let gas = gas(300.0, 1e-6);
        let result = shared(&gas, 0.002).advance(cx, &bodies).unwrap();
        close(result.gas().state().state().temperature, 302.4, 1e-8);
        close(result.report().gas_internal_energy_change_j, 0.002, 1e-11);
        assert_eq!(result.report().heat_from_bodies_j, 0.0);
        assert!(result.bodies()[0].report.radiation_heat_j > 0.0);
        assert!(result.bodies()[1].report.radiation_heat_j < 0.0);
        for (original, next) in bodies.iter().zip(result.bodies()) {
            close(next.report.enthalpy_change_j,
                original.input.internal_heat_j + next.report.radiation_heat_j, 1e-10);
            assert_eq!(next.external_heat_j, next.report.radiation_heat_j);
        }
        let reservoir_entropy: f64 = result.bodies().iter().zip(&boundaries)
            .map(|(body, env)| body.report.radiation_heat_j / env.radiation_temperature_k).sum();
        close(result.report().entropy.radiation_reservoir_input_j_k, reservoir_entropy, 1e-15);
        close(result.report().entropy.balance_residual_j_k, 0.0, 1e-15);
    });
}

#[test]
fn one_body_matches_existing_owner_and_late_failure_keeps_every_input() {
    scope(|cx, gate| {
        let source = sources(2.0, 1000.0, false, false);
        let boundary = boundary(&source, 100.0);
        let gas = gas(300.0, 1e-6);
        let saved = gas.clone();
        let original = input(&source, 320.0, 2e-4);
        let descriptors = [FiniteGasBody { id: 1, environment: &boundary, input: original },
            FiniteGasBody { id: 2, environment: &boundary, input: input(&source, 290.0, 0.0) }];
        let shared = shared(&gas, 0.0);
        let one = shared.advance(cx, &descriptors[..1]).unwrap();
        let existing = FiniteGasThermalEnvironment { gas: &gas, body: &boundary, gas_heat_j: 0.0,
            gas_solve: gas_solve(), exchange: exchange(0.02), maximum_gas_temperature_increment_k: 20.0 }
            .advance(cx, original).unwrap();
        close(one.gas().state().state().temperature, existing.report.next_gas.state().state().temperature, 1e-8);
        close(one.bodies()[0].state.temperature_k(), existing.state.temperature_k(), 1e-8);
        close(one.bodies()[0].report.solid_to_gas_heat_j, existing.report.solid_to_gas_heat_j, 1e-11);
        let duplicate = [descriptors[0], descriptors[0]];
        assert!(shared.advance(cx, &duplicate).is_err());
        assert!(FiniteGasBodiesEnvironment { maximum_bodies: 1, ..shared.clone() }
            .advance(cx, &descriptors).is_err());
        let mut mismatched = descriptors;
        mismatched[1].input.duration_s *= 2.0;
        assert!(shared.advance(cx, &mismatched).is_err());
        let mut calls = 0;
        let refused = shared.advance_with_convection(cx, &descriptors, |body, _| {
            calls += 1;
            if body.id == 2 { Err("second material refused") } else {
                Ok(ConvectionCoefficient { coefficient_w_per_m2_k: 100.0, source: () })
            }
        });
        assert!(matches!(refused, Err(FiniteGasBodiesError::Body { id: 2, .. })));
        assert_eq!(calls, 2);
        let cancelled = shared.advance_with_convection(cx, &descriptors, |body, _| {
            if body.id == 2 { gate.request(); }
            Ok::<_, core::convert::Infallible>(ConvectionCoefficient { coefficient_w_per_m2_k: 100.0, source: () })
        });
        assert!(matches!(cancelled, Err(FiniteGasBodiesError::Cancelled)));
        assert_eq!(gas, saved);
        assert_eq!(original.initial.temperature_k(), 320.0);
    });
}

#[test]
fn dense_gas_keeps_original_fixed_density_departure_energy_with_multiple_bodies() {
    use fs_material::thermochemical_gas::nonideal::{PengRobinsonSpeciesBinding,
        ThermochemicalPengRobinsonModel, energy::parcel::IsochoricRealGasParcel};
    use fs_qty::{Pressure, Temperature};
    use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
        eos::peng_robinson::PengRobinsonPhaseV1};
    scope(|cx, _| {
        let ideal = gas(300.0, 1e-6);
        // Manufactured critical data, not a measured material identification.
        let model = ThermochemicalPengRobinsonModel::new(ideal.model().clone(),
            vec![PengRobinsonSpeciesBinding { species: SpeciesId::new("manufactured-thermal-gas").unwrap(),
                critical_temperature: Temperature::new(190.564), critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142 }], vec![]).unwrap();
        let typed = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
        let point = QueryPoint::new()
            .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), 300.0).unwrap()
            .with_quantity("pressure", typed(QuantityKind::Pressure), 5e6).unwrap();
        let gas = IsochoricRealGasParcel::new(model, point, 1e-6, PengRobinsonPhaseV1::Vapor).unwrap();
        let sources = [sources(2.0, 1000.0, false, false), sources(4.0, 1000.0, false, false)];
        let boundaries = [boundary(&sources[0], 100.0), boundary(&sources[1], 100.0)];
        let inputs = [input(&sources[0], 320.0, 0.0), input(&sources[1], 290.0, 0.0)];
        let bodies = [FiniteGasBody { id: 1, environment: &boundaries[0], input: inputs[0] },
            FiniteGasBody { id: 2, environment: &boundaries[1], input: inputs[1] }];
        let environment = FiniteGasBodiesEnvironment { gas: &gas, gas_heat_j: 1e-3,
            gas_solve: GasEnergySolve { minimum_temperature_k: 280.0,
                maximum_temperature_k: 340.0, ..gas_solve() }, exchange: exchange(0.02),
            maximum_gas_temperature_increment_k: 20.0, maximum_bodies: 2 };
        let result = environment.advance(cx, &bodies).unwrap();
        // Independent pure-fluid PR internal-energy departure at constant molar
        // volume. No native heat inversion or material energy query is an oracle.
        let v = gas.volume_m3() * 0.028 / gas.mass_kg();
        let n = gas.mass_kg() / 0.028;
        let b = 0.07780 * R * 190.564 / 4.5992e6;
        let integral = ((v + (1.0 + 2.0_f64.sqrt()) * b)
            / (v + (1.0 - 2.0_f64.sqrt()) * b)).ln() / (2.0 * 2.0_f64.sqrt() * b);
        let departure = |t: f64| {
            let k = 0.37464 + 1.54226 * 0.01142 - 0.26992 * 0.01142_f64.powi(2);
            let a0 = 0.45724 * R * R * 190.564_f64.powi(2) / 4.5992e6;
            let root = 1.0 + k * (1.0 - (t / 190.564).sqrt());
            let a = a0 * root * root;
            let da = -a0 * k * root / (190.564 * t).sqrt();
            (t * da - a) * integral
        };
        let energy = |t| n * (2.5 * R * (t - 300.0) + departure(t) - departure(300.0));
        let beta = [0.002, 0.004].map(|c| 0.0013 * c / (0.0013 + c));
        let (mut lo, mut hi) = (280.0, 340.0);
        for _ in 0..100 {
            let t = f64::midpoint(lo, hi);
            let residual = energy(t) + beta[0] * (t - 320.0) + beta[1] * (t - 290.0) - environment.gas_heat_j;
            if residual < 0.0 { lo = t; } else { hi = t; }
        }
        let expected = f64::midpoint(lo, hi);
        close(result.gas().state().temperature_k(), expected, 1e-7);
        close(result.report().gas_internal_energy_change_j, energy(expected), 2e-10);
        assert!((energy(expected) - n * 2.5 * R * (expected - 300.0)).abs() > 1e-6);
        assert_eq!(result.gas().mass_kg(), gas.mass_kg());
        assert_eq!(result.gas().volume_m3(), gas.volume_m3());
        close(result.report().entropy.balance_residual_j_k, 0.0, 1e-15);
    });
}
