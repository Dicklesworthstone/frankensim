//! G1/G3/G4: original material contact, mechanical feedback and atomic candidates.
use super::*;
use crate::cavity_phs::membrane::{
    GasMembrane, GasMembraneMotion, MembraneGas, SealedMembraneExterior,
    tests::{
        motion::{motion_controls, real_gas},
        near, pyramid, source_gas,
    },
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use std::{cell::Cell, convert::Infallible};

fn controls() -> GasMembraneThermalControls {
    let old = motion_controls();
    GasMembraneThermalControls {
        gas: GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1800.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
        exterior_gas: None,
        time: old.time,
        maximum_evaluations: 2000,
        maximum_steps: old.maximum_steps,
        maximum_krylov_restart: old.maximum_krylov_restart,
        energy_tolerance_j: 1e-5,
        pressure_work_tolerance_j: 1e-8,
    }
}
fn exchange() -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-4.0, 4.0],
        source_tolerance_j: 1e-10,
        energy_tolerance_j: 1e-5,
        maximum_iterations: 96,
    }
}
fn sensible() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        fs_blake3::hash_bytes(b"synthetic membrane contact cp=500"),
        SolidLiquidPhase::Solid,
        [200.0, 1000.0]
            .into_iter()
            .map(|temperature_k| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 500.0 * (temperature_k - 200.0),
                temperature_k,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 4000.0,
            })
            .collect(),
    )
    .unwrap()
}
fn latent(temperature: f64) -> EquilibriumEnthalpyPhaseCurve {
    let knot = |h, t, f| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: f,
        bulk_density_kg_m3: 4000.0,
    };
    EquilibriumEnthalpyPhaseCurve::try_new(
        fs_blake3::hash_bytes(
            format!("synthetic membrane latent contact T={temperature}").as_bytes(),
        ),
        vec![
            knot(0.0, 200.0, 0.0),
            knot(100000.0, temperature, 0.0),
            knot(160000.0, temperature, 1.0),
            knot(400000.0, 900.0, 1.0),
        ],
    )
    .unwrap()
}
fn body(curve: &EquilibriumEnthalpyPhaseCurve, h: f64) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve,
        initial: curve.state_at_specific_enthalpy(h).unwrap(),
        mass_kg: 0.001,
        volume_m3: 0.001 / 4000.0,
        surface_area_m2: 0.02,
        internal_heat_j: 0.0,
        duration_s: 2e-5,
    }
}
fn link(
    g: f64,
) -> impl FnMut(
    &MembraneGasSample,
    EquilibriumPhaseState,
) -> Result<ConductanceSample<(MembraneGasSample, EquilibriumPhaseState)>, Infallible> {
    move |gas, phase| {
        Ok(ConductanceSample {
            conductance_w_k: g,
            source: (gas.clone(), phase),
        })
    }
}
fn balance<S>(candidate: &GasMembraneContactCandidate<'_, '_, S>) {
    let step = candidate.step();
    assert!(step.source_residual_j.abs() <= exchange().source_tolerance_j);
    assert!(step.thermal_energy_residual_j.abs() <= exchange().energy_tolerance_j);
    assert!(step.total_energy_residual_j.abs() <= exchange().energy_tolerance_j);
    assert!(step.contact_entropy_production_j_k >= 0.0);
    assert!(step.entropy_residual_j_k.abs() <= step.entropy_tolerance_j_k);
    assert_eq!(candidate.thermal().external_heat_j, -step.exchanged_heat_j);
    assert_eq!(
        step.motion.work_pressure_pa,
        step.motion.gas.work_pressure_pa
    );
    assert_eq!(
        step.motion.gas_sampling_temperature_k,
        step.gas_thermodynamic_mean.temperature_k()
    );
    assert_eq!(
        step.motion.gas_sampling_pressure_pa,
        step.gas_thermodynamic_mean.pressure_pa()
    );
    assert!(step.motion_evaluations >= step.motion.evaluations);
    assert!(step.gas_trials >= step.motion.gas_trials);
}
// Independent NASA9 integrals for the original 75:25 synthetic source at slope=.002.
fn gas_changes(old: &MembraneGas, new: &MembraneGas) -> (f64, f64) {
    let (MembraneGas::Ideal(a), MembraneGas::Ideal(b)) = (old, new) else {
        panic!("ideal fixture")
    };
    let (t0, t1) = (a.state().state().temperature, b.state().state().temperature);
    let dt = t1 - t0;
    let mr = a.mass_kg() * a.state().state().specific_gas_constant;
    (
        mr * dt * (2.25 + 0.00075 * (t0 + t1)),
        mr * (2.25 * (dt / t0).ln_1p()
            + 0.0015 * dt
            + ((b.volume_m3() - a.volume_m3()) / a.volume_m3()).ln_1p()),
    )
}

#[test]
fn g1_sensible_contact_drives_pressure_motion_and_retains_native_history() {
    let (patch, reference, fixed) = pyramid();
    let model = GasMembrane::new(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        100000.0,
        [0.0007, 0.0015],
        || false,
    )
    .unwrap();
    let loads = [[0.0; 3]; 5];
    let mut prefix = GasMembraneMotion::new(
        &model,
        &reference,
        &loads,
        &loads,
        motion_controls(),
        || false,
    )
    .unwrap()
    .into_thermal();
    prefix
        .step(1e-5, &loads, 0.01, None, controls(), || false)
        .unwrap();
    let curve = sensible();
    for temperature in [300.0, 650.0] {
        let mut motion = prefix.clone();
        let mut thermal = body(&curve, 500.0 * (temperature - 200.0));
        for _ in 0..3 {
            let old_state = motion.state().clone();
            let old = motion.accepted().clone();
            let old_kinetic = motion.kinetic_energy_j().unwrap();
            let candidate = motion
                .thermal_discrete_gradient_candidate(
                    thermal.duration_s,
                    &loads,
                    0.0,
                    None,
                    thermal,
                    controls(),
                    exchange(),
                    link(300.0),
                    || false,
                )
                .unwrap();
            balance(&candidate);
            assert_eq!(motion.state(), &old_state);
            assert_eq!(motion.accepted().gas, old.gas);
            let report = candidate.step();
            assert!(report.exchanged_heat_j * (temperature - 400.0) > 0.0);
            assert!(report.contact_entropy_production_j_k > 0.0);
            assert_eq!(candidate.thermal().report.0, report.gas_thermodynamic_mean);
            assert_eq!(candidate.thermal().report.1, report.thermal_midpoint);
            let t0 = thermal.initial.temperature_k();
            let t1 = candidate.thermal().state.temperature_k();
            let mean_body = (t1 - t0) / ((t1 - t0) / t0).ln_1p();
            near(report.thermal_mean_temperature_k, mean_body, 2e-9, 1e-12);
            let after = candidate.motion().accepted();
            let (gas_energy, gas_entropy) = gas_changes(&old.gas, &after.gas);
            let body_energy = 0.5 * (t1 - t0);
            let body_entropy = 0.5 * ((t1 - t0) / t0).ln_1p();
            near(
                report.total_entropy_change_j_k,
                gas_entropy + body_entropy,
                2e-11,
                1e-8,
            );
            near(
                report.total_entropy_change_j_k,
                report.contact_entropy_production_j_k,
                1e-9,
                1e-8,
            );
            let total = candidate.motion().kinetic_energy_j().unwrap() - old_kinetic
                + after.strain_energy_j
                - old.strain_energy_j
                + gas_energy
                + body_energy
                + 100000.0 * (after.gas.volume_m3() - old.gas.volume_m3());
            near(report.total_energy_residual_j, total, 2e-9, 1e-7);
            let mut direct = motion.clone();
            let direct_report = direct
                .step_discrete_gradient(
                    thermal.duration_s,
                    &loads,
                    report.exchanged_heat_j,
                    None,
                    controls(),
                    || false,
                )
                .unwrap();
            assert_eq!(candidate.motion().state(), direct.state());
            assert_eq!(after.gas, direct.accepted().gas);
            assert_eq!(report.motion.gas, direct_report.gas);
            assert_eq!(
                &candidate.motion().state().history[..old_state.history.len()],
                old_state.history.as_slice()
            );
            assert_eq!(candidate.motion().state().steps, old_state.steps + 1);
            assert!(candidate.motion().state().q != old_state.q);
            let (next, accepted, _) = candidate.into_parts();
            thermal.initial = accepted.state;
            motion = next;
        }
    }
    // Exact insulating limit uses the unchanged prescribed-heat Gibbs path.
    let thermal = body(&curve, 225000.0);
    let candidate = prefix
        .thermal_discrete_gradient_candidate(
            thermal.duration_s,
            &loads,
            0.02,
            None,
            thermal,
            controls(),
            exchange(),
            link(0.0),
            || false,
        )
        .unwrap();
    let mut direct = prefix.clone();
    direct
        .step_discrete_gradient(thermal.duration_s, &loads, 0.02, None, controls(), || false)
        .unwrap();
    assert_eq!(candidate.motion().state(), direct.state());
    assert_eq!(candidate.thermal().state, thermal.initial);
    assert_eq!(candidate.step().exchange_trials, 1);
    assert_eq!(candidate.step().exchanged_heat_j, 0.0);
    assert_eq!(candidate.step().contact_entropy_production_j_k, 0.0);
}

#[test]
fn g1_latent_contact_and_opposing_gas_keep_separate_heat_entropy_and_work() {
    let (patch, reference, fixed) = pyramid();
    let model = GasMembrane::with_sealed_exterior(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        [0.0007, 0.0015],
        SealedMembraneExterior {
            gas: source_gas(0.75, 0.002).into(),
            volume_interval_m3: [0.0007, 0.0015],
            entropy: motion_controls().gas,
        },
        || false,
    )
    .unwrap();
    let loads = [[0.0; 3]; 5];
    let mut velocity = loads;
    velocity[0][2] = 0.03;
    let motion = GasMembraneMotion::new(
        &model,
        &reference,
        &velocity,
        &loads,
        motion_controls(),
        || false,
    )
    .unwrap()
    .into_thermal();
    let settings = GasMembraneThermalControls {
        exterior_gas: Some(controls().gas),
        ..controls()
    };
    for transition in [350.0, 500.0] {
        let curve = latent(transition);
        let thermal = UniformEnthalpyStepInput {
            internal_heat_j: 0.01,
            ..body(&curve, 130000.0)
        };
        let candidate = motion
            .thermal_discrete_gradient_candidate(
                thermal.duration_s,
                &loads,
                0.02,
                Some(-0.02),
                thermal,
                settings,
                exchange(),
                |gas, phase| {
                    Ok::<_, Infallible>(ConductanceSample {
                        conductance_w_k: 300.0 * (1.0 + phase.liquid_mass_fraction()),
                        source: (gas.clone(), phase),
                    })
                },
                || false,
            )
            .unwrap();
        balance(&candidate);
        let report = candidate.step();
        assert_eq!(candidate.thermal().state.temperature_k(), transition);
        near(report.thermal_mean_temperature_k, transition, 1e-10, 0.0);
        near(
            report.thermal_entropy_change_j_k,
            report.thermal_energy_change_j / transition,
            1e-12,
            1e-12,
        );
        assert!(
            (candidate.thermal().state.liquid_mass_fraction()
                - thermal.initial.liquid_mass_fraction())
                * (400.0 - transition)
                > 0.0
        );
        assert_eq!(
            report.conductance_w_k,
            300.0 * (1.0 + report.thermal_midpoint.liquid_mass_fraction())
        );
        assert_eq!(
            report.motion.gas.supplied_heat_j,
            0.02 + report.exchanged_heat_j
        );
        assert_eq!(report.motion.exterior_gas.unwrap().supplied_heat_j, -0.02);
        let after = candidate.motion().accepted();
        let (_, first_entropy) = gas_changes(&motion.accepted().gas, &after.gas);
        let (_, other_entropy) = gas_changes(
            motion.accepted().exterior_gas.as_ref().unwrap(),
            after.exterior_gas.as_ref().unwrap(),
        );
        near(
            report.total_entropy_change_j_k,
            first_entropy + other_entropy + report.thermal_energy_change_j / transition,
            2e-11,
            1e-8,
        );
        let prescribed = 0.02 / report.gas_thermodynamic_mean.temperature_k()
            - 0.02 / report.motion.exterior_sampling_temperature_k.unwrap()
            + 0.01 / report.thermal_mean_temperature_k;
        near(
            report.prescribed_source_entropy_j_k,
            prescribed,
            1e-17,
            1e-13,
        );
        assert_eq!(report.motion.exterior_pressure_work_j, 0.0);
        near(
            after.gas.volume_m3() + after.exterior_gas.as_ref().unwrap().volume_m3(),
            0.002,
            1e-18,
            0.0,
        );
    }
}

#[test]
fn g3_real_gas_contact_retains_the_original_sampling_state_and_work_pressure() {
    let (patch, reference, fixed) = pyramid();
    let source = real_gas();
    let exterior = source.pressure_response().unwrap().pressure_pa;
    let model = GasMembrane::new(
        &patch,
        source.clone(),
        &reference,
        &fixed,
        exterior,
        [0.0007, 0.0015],
        || false,
    )
    .unwrap();
    let loads = [[0.0; 3]; 5];
    let mut velocity = loads;
    velocity[0][2] = 0.3;
    let motion = GasMembraneMotion::new(
        &model,
        &reference,
        &velocity,
        &loads,
        motion_controls(),
        || false,
    )
    .unwrap()
    .into_thermal();
    let settings = GasMembraneThermalControls {
        gas: GasEnergySolve {
            minimum_temperature_k: 280.0,
            maximum_temperature_k: 330.0,
            ..controls().gas
        },
        ..controls()
    };
    let curve = sensible();
    let thermal = body(&curve, 225000.0);
    let candidate = motion
        .thermal_discrete_gradient_candidate(
            thermal.duration_s,
            &loads,
            0.0,
            None,
            thermal,
            settings,
            exchange(),
            |gas, phase| {
                let MembraneGasSample::PengRobinson(state) = gas else {
                    panic!("native PR sample")
                };
                Ok::<_, Infallible>(ConductanceSample {
                    conductance_w_k: state.temperature_k(),
                    source: (state.clone(), phase),
                })
            },
            || false,
        )
        .unwrap();
    balance(&candidate);
    let report = candidate.step();
    let native = source
        .into_isochoric()
        .discrete_gradient_heat_volume(
            candidate.motion().accepted().gas.volume_m3(),
            report.exchanged_heat_j,
            settings.gas,
            || false,
        )
        .unwrap();
    assert_eq!(report.motion.gas, *native.step());
    assert_eq!(candidate.thermal().report.0, *native.thermodynamic_mean());
    assert_eq!(
        report.gas_thermodynamic_mean,
        MembraneGasSample::PengRobinson(native.thermodynamic_mean().clone())
    );
    assert_eq!(
        report.conductance_w_k,
        native.thermodynamic_mean().temperature_k()
    );
    assert!(
        (report.motion.work_pressure_pa - report.gas_thermodynamic_mean.pressure_pa()).abs() > 1e-8
    );
    let MembraneGas::PengRobinson(endpoint) = &candidate.motion().accepted().gas else {
        panic!("native PR endpoint")
    };
    assert_eq!(&endpoint.clone().into_isochoric(), native.endpoint());
    assert_eq!(
        candidate.motion().accepted().pressure.pressure_pa,
        endpoint.pressure_response().unwrap().pressure_pa
    );
}

#[test]
fn g4_contact_refusal_late_cancellation_and_retry_preserve_both_gases_and_history() {
    let (patch, reference, fixed) = pyramid();
    let model = GasMembrane::with_sealed_exterior(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        [0.0007, 0.0015],
        SealedMembraneExterior {
            gas: source_gas(0.75, 0.002).into(),
            volume_interval_m3: [0.0007, 0.0015],
            entropy: motion_controls().gas,
        },
        || false,
    )
    .unwrap();
    let loads = [[0.0; 3]; 5];
    let settings = GasMembraneThermalControls {
        exterior_gas: Some(controls().gas),
        ..controls()
    };
    let mut motion = GasMembraneMotion::new(
        &model,
        &reference,
        &loads,
        &loads,
        motion_controls(),
        || false,
    )
    .unwrap()
    .into_thermal();
    motion
        .step_discrete_gradient(1e-5, &loads, 0.01, Some(0.0), settings, || false)
        .unwrap();
    let curve = sensible();
    let thermal = body(&curve, 225000.0);
    let state = motion.state().clone();
    let accepted = motion.accepted().clone();
    let unchanged = || {
        assert_eq!(motion.state(), &state);
        assert_eq!(motion.accepted().gas, accepted.gas);
        assert_eq!(motion.accepted().exterior_gas, accepted.exterior_gas);
        assert_eq!(motion.accepted().displacement_m, accepted.displacement_m);
    };
    let polls = Cell::new(0usize);
    let expected = motion
        .thermal_discrete_gradient_candidate(
            thermal.duration_s,
            &loads,
            0.0,
            Some(0.0),
            thermal,
            settings,
            exchange(),
            link(300.0),
            || {
                polls.set(polls.get() + 1);
                false
            },
        )
        .unwrap();
    let final_poll = polls.get();
    let late = Cell::new(0usize);
    assert!(matches!(
        motion.thermal_discrete_gradient_candidate(
            thermal.duration_s,
            &loads,
            0.0,
            Some(0.0),
            thermal,
            settings,
            exchange(),
            link(300.0),
            || {
                late.set(late.get() + 1);
                late.get() == final_poll
            }
        ),
        Err(GasMembraneContactError::Cancelled)
    ));
    unchanged();
    let calls = Cell::new(0usize);
    assert!(matches!(
        motion.thermal_discrete_gradient_candidate(
            thermal.duration_s,
            &loads,
            0.0,
            Some(0.0),
            thermal,
            settings,
            exchange(),
            |_, _| {
                calls.set(calls.get() + 1);
                if calls.get() > 1 {
                    Err("late interface refusal")
                } else {
                    Ok(ConductanceSample {
                        conductance_w_k: 300.0,
                        source: (),
                    })
                }
            },
            || false
        ),
        Err(GasMembraneContactError::Interface("late interface refusal"))
    ));
    unchanged();
    let invalid = UniformEnthalpyStepInput {
        duration_s: thermal.duration_s * 2.0,
        ..thermal
    };
    assert!(matches!(
        motion.thermal_discrete_gradient_candidate(
            thermal.duration_s,
            &loads,
            0.0,
            Some(0.0),
            invalid,
            settings,
            exchange(),
            link(300.0),
            || false
        ),
        Err(GasMembraneContactError::Invalid(_))
    ));
    let limited = GasMembraneThermalControls {
        maximum_evaluations: 1,
        ..settings
    };
    assert!(
        motion
            .thermal_discrete_gradient_candidate(
                thermal.duration_s,
                &loads,
                0.0,
                Some(0.0),
                thermal,
                limited,
                exchange(),
                link(300.0),
                || false
            )
            .is_err()
    );
    assert!(
        motion
            .thermal_discrete_gradient_candidate(
                thermal.duration_s,
                &loads,
                0.0,
                Some(0.0),
                thermal,
                settings,
                exchange(),
                link(-1.0),
                || false
            )
            .is_err()
    );
    unchanged();
    let retry = motion
        .thermal_discrete_gradient_candidate(
            thermal.duration_s,
            &loads,
            0.0,
            Some(0.0),
            thermal,
            settings,
            exchange(),
            link(300.0),
            || false,
        )
        .unwrap();
    assert_eq!(retry.motion().state(), expected.motion().state());
    assert_eq!(
        retry.motion().accepted().gas,
        expected.motion().accepted().gas
    );
    assert_eq!(
        retry.motion().accepted().exterior_gas,
        expected.motion().accepted().exterior_gas
    );
    assert_eq!(retry.thermal(), expected.thermal());
    assert_eq!(
        retry.step().exchanged_heat_j,
        expected.step().exchanged_heat_j
    );
    unchanged();
}
