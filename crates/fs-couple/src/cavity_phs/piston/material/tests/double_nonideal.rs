//! Original solid inertia retained while both dense-gas inventories evolve.
use super::body;
use crate::cavity_phs::piston::double_acting::{
    DoubleActingGasPiston, DoublePistonInput, DoublePistonSolve,
};
use crate::cavity_phs::piston::material::{MaterialDoubleActingGasPiston, MaterialGasPiston};
use crate::cavity_phs::piston::tests::real_gas::{real_controls, real_parcel};
use crate::cavity_phs::piston::{GasPiston, PistonError, PistonInput};
use fs_ports::PortedModule;

fn input() -> DoublePistonInput {
    DoublePistonInput {
        duration_s: 2e-6,
        heat_j: [20.0, -5.0],
        external_force_n: 0.0,
    }
}
fn controls() -> DoublePistonSolve {
    DoublePistonSolve {
        piston: real_controls(),
        opposing_gas: real_controls().gas,
        volume_tolerance_m3: 1e-18,
    }
}

#[test]
fn source_density_changes_shared_motion_and_both_densities_with_one_fixed_solid_mass() {
    let primary = real_parcel(300.0, 5e6);
    let opposing = real_parcel(330.0, 5e6);
    let run = |density| {
        let specimen = body(density);
        let mut native = DoubleActingGasPiston::new(
            GasPiston::new(primary.clone(), specimen.inertia().mass(), 0.01, 0.0).unwrap(),
            opposing.clone(),
            0.01,
        )
        .unwrap();
        let mut bound = MaterialDoubleActingGasPiston::new(
            primary.clone(),
            opposing.clone(),
            specimen.clone(),
            [0.01, 0.01],
            0.0,
        )
        .unwrap();
        let mut heat = 0.0;
        for _ in 0..3 {
            let expected = native.step(input(), controls(), || false).unwrap();
            let actual = bound.step(input(), controls(), || false).unwrap();
            assert_eq!(actual, expected);
            assert_eq!(bound.piston(), &native);
            assert_eq!(bound.body(), &specimen);
            assert_eq!(
                bound.piston().primary().mass_kg(),
                specimen.inertia().mass()
            );
            assert_eq!(
                bound.piston().last_window_transfers(),
                native.last_window_transfers()
            );
            assert_eq!(
                bound.piston().audit_last_window().unwrap(),
                native.audit_last_window().unwrap()
            );
            heat += input().heat_j.into_iter().sum::<f64>();
            let mut energy =
                0.5 * specimen.inertia().mass() * bound.piston().primary().velocity_m_s().powi(2)
                    - heat;
            for (gas, initial) in bound.gases().into_iter().zip([&primary, &opposing]) {
                assert_eq!(gas.mass_kg(), initial.mass_kg());
                assert_eq!(gas.model(), initial.model());
                assert_eq!(gas.point().axes().get("specimen-lot"), Some(&42.0));
                assert_eq!(gas.density_kg_m3(), gas.mass_kg() / gas.volume_m3());
                assert_ne!(gas.density_kg_m3(), gas.state().source().state().density);
                energy += gas.mass_kg()
                    * gas
                        .model()
                        .specific_internal_energy_change_between_densities(
                            initial.density_kg_m3(),
                            initial.state().temperature_k(),
                            gas.density_kg_m3(),
                            gas.state().temperature_k(),
                        )
                        .unwrap();
            }
            assert!(energy.abs() <= 4.0 * controls().piston.energy_tolerance_j);
            assert!(actual.total_energy_residual_j.abs() <= controls().piston.energy_tolerance_j);
        }
        bound
    };
    let light = run(1000.0);
    let dense = run(8000.0);
    assert!(light.piston().primary().velocity_m_s() > dense.piston().primary().velocity_m_s());
    assert!(dense.piston().primary().velocity_m_s() > 0.0);
    assert!(light.gases()[0].density_kg_m3() < dense.gases()[0].density_kg_m3());
    assert!(light.gases()[1].density_kg_m3() > dense.gases()[1].density_kg_m3());
}

#[test]
fn adding_a_real_opposing_volume_retains_history_and_refused_steps_publish_neither_side() {
    let mut single =
        MaterialGasPiston::new(real_parcel(300.0, 5e6), body(8000.0), 0.01, 0.0).unwrap();
    single
        .step(
            PistonInput {
                duration_s: input().duration_s,
                heat_j: 1.0,
                exterior_pressure_pa: 5e6,
                external_force_n: 0.0,
            },
            real_controls(),
            || false,
        )
        .unwrap();
    let accepted = single.piston().clone();
    let specimen = single.body().clone();
    let opposing = real_parcel(330.0, 5e6);
    let mut bound = single.with_opposing_gas(opposing.clone(), 0.01).unwrap();
    assert_eq!(bound.body(), &specimen);
    assert_eq!(bound.piston().primary(), &accepted);
    assert_eq!(bound.piston().opposing(), &opposing);
    assert_eq!(
        bound.piston().primary().last_window_transfers(),
        accepted.last_window_transfers()
    );
    bound.step(input(), controls(), || false).unwrap();
    let saved = bound.clone();
    assert!(
        bound
            .step(
                DoublePistonInput {
                    heat_j: [0.0, 1e9],
                    ..input()
                },
                controls(),
                || false
            )
            .is_err()
    );
    assert_eq!(bound, saved);
    let mut probe = saved.clone();
    let mut polls = 0;
    probe
        .step(input(), controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        bound.step(input(), controls(), || {
            seen += 1;
            seen == polls
        }),
        Err(PistonError::Cancelled)
    ));
    assert_eq!(bound, saved);
    assert_eq!(
        bound.piston().last_window_transfers(),
        saved.piston().last_window_transfers()
    );
    assert_eq!(
        bound.piston().audit_last_window().unwrap(),
        saved.piston().audit_last_window().unwrap()
    );
}

#[test]
fn one_sourced_solid_mass_controls_both_inertia_and_finite_heat_capacity() {
    use crate::cavity_phs::piston::thermal::double_acting::DoubleHeatExchangeSolve;
    use crate::cavity_phs::piston::thermal::{
        ConductanceSample, HeatExchangeSolve, ThermalPistonError,
    };
    use fs_material::phase::{
        EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, SolidLiquidPhase,
        UniformEnthalpyStepInput,
    };
    use fs_material::thermochemical_gas::nonideal::ResolvedThermochemicalRealGasState;
    use std::convert::Infallible;

    let primary = real_parcel(300.0, 5e6);
    let opposing = real_parcel(330.0, 5e6);
    let heat = DoubleHeatExchangeSolve {
        faces: [HeatExchangeSolve {
            heat_bracket_j: [-0.001, 0.001],
            source_tolerance_j: 1e-9,
            energy_tolerance_j: 2e-6,
            maximum_iterations: 128,
        }; 2],
        maximum_joint_trials: 256,
    };
    let resolve = |gases: &[ResolvedThermochemicalRealGasState; 2], _: EquilibriumPhaseState| {
        Ok::<_, Infallible>([
            ConductanceSample {
                conductance_w_k: 1.0,
                source: gases[0].clone(),
            },
            ConductanceSample {
                conductance_w_k: 0.8,
                source: gases[1].clone(),
            },
        ])
    };
    let run = |density| {
        let specimen = body(density);
        let bound = MaterialDoubleActingGasPiston::new(
            primary.clone(),
            opposing.clone(),
            specimen.clone(),
            [0.01; 2],
            0.0,
        )
        .unwrap();
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            super::card(density).content_hash(),
            SolidLiquidPhase::Solid,
            [200.0, 400.0]
                .map(|temperature_k| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 500.0 * (temperature_k - 200.0),
                    temperature_k,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: density,
                })
                .to_vec(),
        )
        .unwrap();
        let thermal = UniformEnthalpyStepInput {
            curve: &curve,
            initial: curve
                .state_at_specific_enthalpy(500.0 * (325.0 - 200.0))
                .unwrap(),
            mass_kg: specimen.inertia().mass(),
            volume_m3: 0.0001,
            surface_area_m2: 0.02,
            duration_s: input().duration_s,
            internal_heat_j: 0.0,
        };
        let saved = bound.clone();
        let mut polls = 0;
        let candidate = bound
            .thermal_candidate(input(), thermal, controls(), heat, resolve, || {
                polls += 1;
                false
            })
            .unwrap();
        assert_eq!(bound, saved);
        assert_eq!(candidate.piston().body(), &specimen);
        assert_eq!(candidate.thermal().report, candidate.step().piston.midpoint);
        assert!(candidate.step().exchanged_heat_j[0] > 0.0);
        assert!(candidate.step().exchanged_heat_j[1] < 0.0);
        let body_energy = thermal.mass_kg
            * (candidate.thermal().state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        let mut stored = body_energy
            + 0.5
                * specimen.inertia().mass()
                * candidate.piston().piston().primary().velocity_m_s().powi(2);
        for (gas, initial) in candidate
            .piston()
            .gases()
            .into_iter()
            .zip([&primary, &opposing])
        {
            assert_eq!(gas.mass_kg(), initial.mass_kg());
            assert_eq!(gas.model(), initial.model());
            stored += gas.mass_kg()
                * gas
                    .model()
                    .specific_internal_energy_change_between_densities(
                        initial.density_kg_m3(),
                        initial.state().temperature_k(),
                        gas.density_kg_m3(),
                        gas.state().temperature_k(),
                    )
                    .unwrap();
        }
        assert!(
            (stored - input().heat_j.into_iter().sum::<f64>()).abs()
                <= heat.faces[0].energy_tolerance_j
        );
        assert!(
            (body_energy + candidate.step().exchanged_heat_j.into_iter().sum::<f64>()).abs()
                <= heat.faces[0].energy_tolerance_j
        );
        let mut seen = 0;
        assert!(matches!(
            bound.thermal_candidate(input(), thermal, controls(), heat, resolve, || {
                seen += 1;
                seen == polls
            }),
            Err(ThermalPistonError::Cancelled)
        ));
        assert_eq!(bound, saved);
        assert_eq!(
            bound.piston().last_window_transfers(),
            saved.piston().last_window_transfers()
        );
        (
            candidate.piston().piston().primary().velocity_m_s(),
            candidate.thermal().state.temperature_k(),
        )
    };
    let (light_speed, light_temperature) = run(1000.0);
    let (dense_speed, dense_temperature) = run(8000.0);
    assert!(light_speed > dense_speed && dense_speed > 0.0);
    assert!(light_temperature < dense_temperature && dense_temperature < 325.0);
}
