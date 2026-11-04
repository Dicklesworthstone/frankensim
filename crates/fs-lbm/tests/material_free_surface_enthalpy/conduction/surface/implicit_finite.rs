//! Actual coupled gas/liquid backward Euler, original sources and atomic refusal.
use super::*;
use fs_convection::{CorrelationId, CorrelationInputs, ThermalConductivity};
use fs_evidence::Evidence;
use fs_lbm::material::free_surface::enthalpy::conduction::surface::convection::{
    FreeSurfaceConvectionQuery, FreeSurfaceGasConvection, FreeSurfaceHeatTransferCoefficient,
};
use fs_lbm::material::free_surface::enthalpy::conduction::surface::{
    FreeSurfaceGasHeatTransfer, FreeSurfaceImplicitSolve,
};
use fs_material::fluid_transport::FluidTransport;
use fs_qty::Length;
use std::cell::RefCell;

const R_GAS: f64 = 8.314_462_618_153_24 / 0.028;
const LIQUID_CAPACITY: f64 = 0.002;

fn implicit() -> FreeSurfaceInterfaceSolve {
    FreeSurfaceInterfaceSolve {
        implicit: Some(FreeSurfaceImplicitSolve::default()),
        ..FreeSurfaceInterfaceSolve::default()
    }
}

fn isolated_face() -> FreeSurfaceInterfaceFace {
    FreeSurfaceInterfaceFace {
        emissivity: 0.0.into(),
        ..face(400.0)
    }
}

// Independent manufactured NASA9 integral: cv/R = 2 + 0.001 T. A formation
// energy reference cancels before arithmetic, independently of the gas solver.
fn delta_u(t0: f64, t1: f64) -> f64 {
    R_GAS * (t1 - t0) * (2.0 + 0.0005 * (t1 + t0))
}

fn conductance(tg: f64) -> f64 {
    let pressure = 100000.0 * tg / 450.0;
    let kg = 2.6 * (0.5 + (pressure - 1000.0) / 249000.0);
    let gl = 0.5 * 1e-6 / 0.0005;
    let gg = kg * 1e-6 / 0.0005;
    gl * gg / (gl + gg)
}

fn reference_be(mass: f64, dt: f64) -> (f64, f64) {
    let (mut lo, mut hi) = (400.0, 450.0);
    for _ in 0..80 {
        let tg = 0.5 * (lo + hi);
        let stored = mass * delta_u(450.0, tg);
        let tl = 400.0 - stored / LIQUID_CAPACITY;
        if stored + dt * conductance(tg) * (tg - tl) < 0.0 {
            lo = tg;
        } else {
            hi = tg;
        }
    }
    let tg = 0.5 * (lo + hi);
    (tg, 400.0 - mass * delta_u(450.0, tg) / LIQUID_CAPACITY)
}

#[test]
fn stiff_finite_atmosphere_uses_one_joint_interval_and_original_energy_sources() {
    let liquid = card(0.02, 1000.0);
    let run = |implicit_step: bool, formation_reference: f64| {
        with_exposed_liquid(&liquid, false, |thermal, gate| {
            let mut gas = IsochoricGasParcel::new(
                gas_source::model(formation_reference, 250000.0, 100.0),
                gas_source::point(450.0, 100000.0),
                1e-10,
            )
            .unwrap();
            let before = thermal.total_enthalpy_j().unwrap();
            let report = thermal
                .step_with_finite_gas_surface_exchange(
                    &[],
                    &[],
                    FreeSurfaceConductionControls {
                        maximum_dt_s: 0.001,
                        ..controls()
                    },
                    &mut [],
                    &[],
                    &[],
                    &mut gas,
                    &[isolated_face()],
                    gas_controls(),
                    if implicit_step {
                        implicit()
                    } else {
                        FreeSurfaceInterfaceSolve::default()
                    },
                    false,
                    gate,
                )
                .unwrap()
                .conduction;
            (report, gas, thermal.total_enthalpy_j().unwrap() - before)
        })
    };
    let (report, gas, liquid_change) = run(true, 0.0);
    let surface = report.surface.as_ref().unwrap();
    let finite = report.finite_gas.as_ref().unwrap();
    let work = surface.implicit.as_ref().unwrap();
    assert_eq!(work.intervals, 1);
    assert_eq!(finite.substeps, 1);
    let (expected_gas, expected_liquid) = reference_be(gas.mass_kg(), 0.001);
    assert!((gas.state().state().temperature - expected_gas).abs() < 2e-7);
    assert!((surface.thermal_endpoint[0].liquid_temperature_k - expected_liquid).abs() < 2e-9);
    let gas_change = gas.mass_kg() * delta_u(450.0, gas.state().state().temperature);
    assert!((liquid_change + gas_change).abs() < 1e-9);
    assert!((gas_change + surface.heat_from_gas_j).abs() <= finite.gas_energy_tolerance_j);
    assert!(
        work.maximum_gas_residual_j <= FreeSurfaceImplicitSolve::default().cell_energy_tolerance_j
    );
    assert_eq!(
        finite.implicit_energy_tolerance_j,
        FreeSurfaceImplicitSolve::default().cell_energy_tolerance_j
    );
    near(
        gas.state().state().pressure,
        100000.0 * gas.state().state().temperature / 450.0,
        2e-12,
    );
    let (explicit, _, _) = run(false, 0.0);
    assert!(
        explicit.thermal.substeps >= 10 * report.thermal.substeps,
        "stiff gas must force substantially more original explicit CFL intervals"
    );
    let (_, shifted, shifted_liquid_change) = run(true, 1e12);
    assert_eq!(
        gas.state().state().temperature,
        shifted.state().state().temperature
    );
    assert_eq!(liquid_change, shifted_liquid_change);
}

#[derive(Default)]
struct FilmPlate<'a> {
    calls: RefCell<Vec<[f64; 4]>>,
    cancel: Option<(&'a CancelGate, f64)>,
}
impl FreeSurfaceGasConvection for FilmPlate<'_> {
    fn coefficient(
        &self,
        q: FreeSurfaceConvectionQuery<'_>,
    ) -> Result<Option<Evidence<FreeSurfaceHeatTransferCoefficient>>, FreeSurfaceEnthalpyError>
    {
        self.calls.borrow_mut().push([
            q.surface_temperature_k,
            q.bulk.state().temperature,
            q.film.state().temperature,
            q.bulk.state().pressure,
        ]);
        if let Some((gate, threshold)) = self.cancel {
            if q.bulk.state().temperature > threshold {
                gate.request();
            }
        }
        let state = q
            .film
            .transport_state()
            .map_err(|e| FreeSurfaceEnthalpyError::Convection(e.to_string()))?;
        let re = state
            .reynolds(10.0, 0.01)
            .map_err(|e| FreeSurfaceEnthalpyError::Convection(e.to_string()))?;
        fs_convection::evaluate(
            CorrelationId::FlatPlateLaminarAverage,
            CorrelationInputs::forced(re, state.prandtl()),
        )
        .and_then(|nu| {
            nu.heat_transfer_coefficient(
                ThermalConductivity::new(state.thermal_conductivity_w_m_k()),
                Length::new(0.01),
            )
        })
        .map(Some)
        .map_err(|e| FreeSurfaceEnthalpyError::Convection(e.to_string()))
    }
}

#[test]
fn joint_convection_radiation_uses_evolving_original_gas_and_independent_union_heat() {
    let liquid = card(0.02, 1000.0);
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let mut gas = IsochoricGasParcel::new(
            gas_source::model(1e12, 250000.0, 1.0),
            gas_source::point(400.0, 100000.0),
            1e-9,
        )
        .unwrap();
        let before = thermal.total_enthalpy_j().unwrap();
        let plate = FilmPlate::default();
        let authored = 1e-5;
        let report = thermal
            .step_with_finite_gas_surface_convection(
                &[],
                &[FreeSurfaceCellHeat {
                    coordinates: EXPOSED,
                    energy_j: authored,
                }],
                FreeSurfaceConductionControls {
                    maximum_dt_s: 0.001,
                    ..controls()
                },
                &mut [],
                &[],
                &[],
                &mut gas,
                &[face(1000.0)],
                gas_controls(),
                implicit(),
                &plate,
                None,
                false,
                gate,
            )
            .unwrap()
            .conduction;
        let exchange = report.surface.as_ref().unwrap();
        let finite = report.finite_gas.as_ref().unwrap();
        let gas_change = gas.mass_kg() * delta_u(400.0, gas.state().state().temperature);
        assert!(gas_change > 0.0 && exchange.heat_from_gas_j < 0.0);
        assert!(
            (thermal.total_enthalpy_j().unwrap() - before + gas_change
                - authored
                - exchange.heat_from_radiation_j)
                .abs()
                < 1e-9
        );
        assert!((gas_change + exchange.heat_from_gas_j).abs() <= finite.gas_energy_tolerance_j);
        assert!(plate.calls.borrow().iter().any(|q| q[1] > 401.0));
        for &[ts, tg, film, pressure] in plate.calls.borrow().iter() {
            assert_eq!(film, 0.5 * ts + 0.5 * tg);
            near(pressure, 100000.0 * tg / 400.0, 2e-12);
        }
        let endpoint = &exchange.thermal_endpoint[0];
        assert!(matches!(
            endpoint.gas_heat_transfer,
            FreeSurfaceGasHeatTransfer::Convection(Some(_))
        ));
        assert_eq!(
            endpoint.gas,
            gas.model()
                .evaluate(&gas_source::point(
                    0.5 * endpoint.surface_temperature_k + 0.5 * endpoint.gas_temperature_k,
                    gas.state().state().pressure,
                ))
                .unwrap()
        );
        assert!(finite.combined_thermal_balance_defect_j.abs() < 1e-9);
    });
}

#[test]
fn implicit_gas_late_domain_and_cancellation_preserve_the_entire_accepted_prefix() {
    let liquid = card(0.02, 1000.0);
    let run = |refuse_first: bool| {
        with_exposed_liquid(&liquid, true, |thermal, gate| {
            let mut gas = IsochoricGasParcel::new(
                gas_source::model(0.0, 250000.0, 1.0),
                gas_source::point(400.0, 100000.0),
                1e-9,
            )
            .unwrap();
            if refuse_first {
                let grid = thermal.flow().lattice().grid().state_bits();
                let energy = thermal.total_enthalpy_j().unwrap();
                let source = thermal.material_at(EXPOSED).unwrap().clone();
                let original_gas = gas.clone();
                let plate = FilmPlate::default();
                let restricted = GasEnergySolve {
                    maximum_temperature_k: 403.0,
                    ..gas_controls()
                };
                let domain_failure = thermal
                    .step_with_finite_gas_surface_convection(
                        &[],
                        &[],
                        controls(),
                        &mut [],
                        &[],
                        &[],
                        &mut gas,
                        &[face(1000.0)],
                        restricted,
                        implicit(),
                        &plate,
                        None,
                        false,
                        gate,
                    )
                    .unwrap_err();
                assert!(
                    matches!(domain_failure, FreeSurfaceEnthalpyError::Gas(_)),
                    "expected the original bounded gas-energy owner: {domain_failure:?}"
                );
                assert!(
                    plate.calls.borrow().iter().any(|q| q[1] > 401.0),
                    "the original gas must reach a later supported trial before refusal"
                );
                assert_eq!(thermal.flow().lattice().grid().state_bits(), grid);
                assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
                assert_eq!(thermal.material_at(EXPOSED).unwrap(), &source);
                assert_eq!(gas, original_gas);
                let cancellation = CancelGate::new_clock_free();
                let cancelling = FilmPlate {
                    cancel: Some((&cancellation, 401.0)),
                    ..FilmPlate::default()
                };
                let failure = thermal
                    .step_with_finite_gas_surface_convection(
                        &[],
                        &[],
                        controls(),
                        &mut [],
                        &[],
                        &[],
                        &mut gas,
                        &[face(1000.0)],
                        gas_controls(),
                        implicit(),
                        &cancelling,
                        None,
                        false,
                        &cancellation,
                    )
                    .unwrap_err();
                assert!(cancellation.is_requested());
                assert_eq!(
                    failure,
                    FreeSurfaceEnthalpyError::from(
                        fs_conduction::enthalpy::DepositionError::Cancelled
                    )
                );
                assert_eq!(thermal.flow().lattice().grid().state_bits(), grid);
                assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
                assert_eq!(thermal.material_at(EXPOSED).unwrap(), &source);
                assert_eq!(gas, original_gas);
                let exhausted = FreeSurfaceInterfaceSolve {
                    implicit: Some(FreeSurfaceImplicitSolve {
                        maximum_residual_evaluations: 1,
                        ..FreeSurfaceImplicitSolve::default()
                    }),
                    ..implicit()
                };
                assert!(
                    thermal
                        .step_with_finite_gas_surface_exchange(
                            &[],
                            &[],
                            controls(),
                            &mut [],
                            &[],
                            &[],
                            &mut gas,
                            &[face(1000.0)],
                            gas_controls(),
                            exhausted,
                            false,
                            gate,
                        )
                        .is_err()
                );
                assert_eq!(gas, original_gas);
                assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
                // The fixture has 144 liquid stores. The appended gas is a
                // real 145th unknown and must consume the same explicit budget.
                let dimension = FreeSurfaceInterfaceSolve {
                    implicit: Some(FreeSurfaceImplicitSolve {
                        maximum_unknowns: 144,
                        ..FreeSurfaceImplicitSolve::default()
                    }),
                    ..implicit()
                };
                assert!(
                    thermal
                        .step_with_finite_gas_surface_exchange(
                            &[],
                            &[],
                            controls(),
                            &mut [],
                            &[],
                            &[],
                            &mut gas,
                            &[face(1000.0)],
                            gas_controls(),
                            dimension,
                            false,
                            gate,
                        )
                        .is_err()
                );
                assert_eq!(gas, original_gas);
                assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
            }
            let report = thermal
                .step_with_finite_gas_surface_convection(
                    &[],
                    &[],
                    controls(),
                    &mut [],
                    &[],
                    &[],
                    &mut gas,
                    &[face(1000.0)],
                    gas_controls(),
                    implicit(),
                    &FilmPlate::default(),
                    None,
                    false,
                    gate,
                )
                .unwrap()
                .conduction;
            (
                thermal.flow().lattice().grid().state_bits(),
                thermal.total_enthalpy_j().unwrap(),
                gas,
                report.surface.unwrap().heat_from_gas_j,
            )
        })
    };
    assert_eq!(run(true), run(false));
}
