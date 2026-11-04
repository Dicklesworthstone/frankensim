//! Actual surface/thermal/flow owner checks using the original plate correlation.
use super::*;
use fs_convection::{CorrelationId, CorrelationInputs, ThermalConductivity};
use fs_evidence::Evidence;
use fs_lbm::material::free_surface::enthalpy::conduction::surface::convection::{
    FreeSurfaceConvectionQuery, FreeSurfaceGasConvection, FreeSurfaceGasHeatTransfer,
    FreeSurfaceHeatTransferCoefficient,
};
use fs_material::fluid_transport::FluidTransport;
use fs_qty::Length;
use std::cell::RefCell;

const SPEED: f64 = 10.0;
const LENGTH: f64 = 0.01;

// The higher production adapter lives in fs-airflow. This lower-owner test
// binds the same original correlation directly, avoiding a crate cycle.
#[derive(Default)]
struct Plate {
    minimum_bulk_k: Option<f64>,
    calls: RefCell<Vec<[f64; 3]>>,
}
impl FreeSurfaceGasConvection for Plate {
    fn coefficient(
        &self,
        q: FreeSurfaceConvectionQuery<'_>,
    ) -> Result<Option<Evidence<FreeSurfaceHeatTransferCoefficient>>, FreeSurfaceEnthalpyError>
    {
        self.calls.borrow_mut().push([
            q.surface_temperature_k,
            q.bulk.state().temperature,
            q.film.state().temperature,
        ]);
        assert_eq!(q.outward_normal, [1.0, 0.0, 0.0]);
        if self
            .minimum_bulk_k
            .is_some_and(|limit| q.bulk.state().temperature < limit)
        {
            return Err(FreeSurfaceEnthalpyError::Convection(
                "manufactured late source refusal".into(),
            ));
        }
        let transport = q
            .film
            .transport_state()
            .map_err(|e| FreeSurfaceEnthalpyError::Convection(e.to_string()))?;
        let reynolds = transport
            .reynolds(SPEED, LENGTH)
            .map_err(|e| FreeSurfaceEnthalpyError::Convection(e.to_string()))?;
        let result = fs_convection::evaluate(
            CorrelationId::FlatPlateLaminarAverage,
            CorrelationInputs::forced(reynolds, transport.prandtl()),
        )
        .and_then(|value| {
            value.heat_transfer_coefficient(
                ThermalConductivity::new(transport.thermal_conductivity_w_m_k()),
                Length::new(LENGTH),
            )
        })
        .map_err(|e| FreeSurfaceEnthalpyError::Convection(e.to_string()))?;
        Ok(Some(result))
    }
}

fn plate_gas() -> ThermochemicalGasModel {
    // Manufactured source k(T)=0.00005 T, with the existing NASA9 caloric
    // owner and pressure-dependent original viscosity; Pr remains admitted.
    gas_with_conductivity(PropertyValue::Curve {
        abscissa: "temperature".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(200.0, 0.01), (700.0, 0.035)],
        dims: Dims([1, 1, -3, -1, 0, 0]),
    })
}

fn convective_face() -> FreeSurfaceInterfaceFace {
    FreeSurfaceInterfaceFace {
        emissivity: 0.0.into(),
        ..face(400.0)
    }
}

#[test]
fn original_plate_law_uses_trial_surface_film_and_liquid_series_resistance() {
    let liquid = card(0.02, 1000.0);
    let model = plate_gas();
    let point = gas_source::point(450.0, 100000.0);
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let plate = Plate::default();
        let mut surface = convective_face();
        surface.face.distances_m[1] = 0.05;
        let before = thermal.total_enthalpy_j().unwrap();
        let report = thermal
            .step_with_surface_convection(
                &[],
                &[],
                controls(),
                &mut [],
                &[],
                &[],
                &[FreeSurfaceInterfaceContact {
                    surface,
                    model: &model,
                    context: &point,
                }],
                FreeSurfaceInterfaceSolve::default(),
                &plate,
                false,
                gate,
            )
            .unwrap()
            .conduction;
        let exchange = report.surface.as_ref().unwrap();
        assert!(plate.calls.borrow().len() > report.thermal.substeps);
        assert!(
            plate
                .calls
                .borrow()
                .iter()
                .all(|[ts, tg, film]| *film == 0.5 * ts + 0.5 * tg)
        );
        for sample in exchange
            .last_applied
            .iter()
            .chain(&exchange.thermal_endpoint)
        {
            let gas = sample.gas.state();
            // Independent average-plate equation and group definitions.
            let re = gas.density * SPEED * LENGTH / gas.dynamic_viscosity;
            let pr = gas.specific_heat_cp * gas.dynamic_viscosity / gas.thermal_conductivity;
            let h = 0.664 * re.sqrt() * pr.cbrt() * gas.thermal_conductivity / LENGTH;
            let FreeSurfaceGasHeatTransfer::Convection(Some(evidence)) = &sample.gas_heat_transfer
            else {
                panic!("lost original correlation evidence");
            };
            near(evidence.value.value(), h, 2e-12);
            near(sample.gas_conductance_w_k, h * surface.face.area_m2, 2e-12);
            assert!(
                sample.gas_conductance_w_k
                    > 10.0 * gas.thermal_conductivity * surface.face.area_m2
                        / surface.face.distances_m[1]
            );
            let gl = 0.5 * surface.face.area_m2 / surface.face.distances_m[0];
            let gg = h * surface.face.area_m2;
            let ts = (gl * sample.liquid_temperature_k + gg * 450.0) / (gl + gg);
            near(sample.surface_temperature_k, ts, 2e-11);
            assert!(sample.surface_temperature_k > sample.liquid_temperature_k + 0.1);
            near(sample.heat_from_gas_w, gg * (450.0 - ts), 2e-10);
            near(
                sample.liquid_decay_conductance_w_k,
                gl * gg / (gl + gg),
                2e-12,
            );
            assert_eq!(
                sample.gas,
                model
                    .evaluate(&gas_source::point(
                        0.5 * sample.surface_temperature_k + 225.0,
                        100000.0
                    ))
                    .unwrap()
            );
            assert_eq!(
                sample.liquid.material().card_identity(),
                liquid.content_hash()
            );
        }
        near(
            thermal.total_enthalpy_j().unwrap() - before,
            exchange.heat_from_gas_j,
            2e-8,
        );
        near(report.flow.external_heat_j, exchange.heat_from_gas_j, 2e-12);
        assert_eq!(exchange.heat_from_radiation_j, 0.0);
    });
}

#[test]
fn finite_gas_convection_is_reciprocal_and_uses_the_original_capacity_step_bound() {
    let liquid = card(0.02, 1000.0);
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let mut gas =
            IsochoricGasParcel::new(plate_gas(), gas_source::point(450.0, 100000.0), 2e-11)
                .unwrap();
        let gas_before = gas.state().specific_internal_energy_j_kg();
        let liquid_before = thermal.total_enthalpy_j().unwrap();
        let plate = Plate::default();
        let report = thermal
            .step_with_finite_gas_surface_convection(
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
                &[convective_face()],
                gas_controls(),
                FreeSurfaceInterfaceSolve::default(),
                &plate,
                None,
                false,
                gate,
            )
            .unwrap()
            .conduction;
        let exchange = report.surface.as_ref().unwrap();
        let finite = report.finite_gas.as_ref().unwrap();
        assert!(
            report.thermal.substeps > 1,
            "finite gas must shorten the otherwise single thermal interval"
        );
        assert_eq!(report.thermal.substeps, finite.substeps);
        let gas_change = gas.mass_kg() * (gas.state().specific_internal_energy_j_kg() - gas_before);
        assert!(gas_change < 0.0);
        assert!((gas_change + exchange.heat_from_gas_j).abs() <= finite.gas_energy_tolerance_j);
        assert!((thermal.total_enthalpy_j().unwrap() - liquid_before + gas_change).abs() < 1e-9);
        assert!(finite.combined_thermal_balance_defect_j.abs() < 1e-9);
        assert!((400.0..=450.0).contains(&gas.state().state().temperature));
        near(
            gas.state().state().pressure,
            100000.0 * gas.state().state().temperature / 450.0,
            2e-12,
        );
        assert!(plate.calls.borrow().iter().any(|q| q[1] < 449.0));
    });
}

#[test]
fn late_convective_refusal_rolls_back_liquid_wall_gas_and_the_next_accepted_step() {
    let liquid = card(0.02, 1000.0);
    let solid = body::wall_card(0.5, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    let run = |reject_first: bool| {
        with_exposed_liquid(&liquid, true, |thermal, gate| {
            let mut gas =
                IsochoricGasParcel::new(plate_gas(), gas_source::point(450.0, 100000.0), 1e-8)
                    .unwrap();
            let mut bodies = [body::wall(&solid, &chart, &context, 450.0)];
            if reject_first {
                let grid = thermal.flow().lattice().grid().state_bits();
                let before = thermal.total_enthalpy_j().unwrap();
                let liquid_source = thermal.material_at(EXPOSED).unwrap().clone();
                let wall_before = bodies[0].specific_enthalpy_j_kg();
                let wall_source = bodies[0].material().clone();
                let gas_before = gas.state().clone();
                let gas_point = gas.point().clone();
                let refusing = Plate {
                    minimum_bulk_k: Some(449.99),
                    ..Plate::default()
                };
                let failure = thermal
                    .step_with_finite_gas_surface_convection(
                        &[],
                        &[],
                        controls(),
                        &mut bodies,
                        &[body::contact()],
                        &[],
                        &mut gas,
                        &[convective_face()],
                        gas_controls(),
                        FreeSurfaceInterfaceSolve::default(),
                        &refusing,
                        None,
                        true,
                        gate,
                    )
                    .unwrap_err();
                assert!(format!("{failure}").contains("manufactured late source refusal"));
                assert!(refusing.calls.borrow().iter().any(|q| q[1] == 450.0));
                assert!(refusing.calls.borrow().iter().any(|q| q[1] < 449.99));
                assert_eq!(thermal.flow().lattice().grid().state_bits(), grid);
                assert_eq!(thermal.total_enthalpy_j().unwrap(), before);
                assert_eq!(thermal.material_at(EXPOSED).unwrap(), &liquid_source);
                assert_eq!(bodies[0].specific_enthalpy_j_kg(), wall_before);
                assert_eq!(bodies[0].material(), &wall_source);
                assert_eq!(gas.state(), &gas_before);
                assert_eq!(gas.point(), &gas_point);
            }
            let report = thermal
                .step_with_finite_gas_surface_convection(
                    &[],
                    &[],
                    controls(),
                    &mut bodies,
                    &[body::contact()],
                    &[],
                    &mut gas,
                    &[convective_face()],
                    gas_controls(),
                    FreeSurfaceInterfaceSolve::default(),
                    &Plate::default(),
                    None,
                    true,
                    gate,
                )
                .unwrap()
                .conduction;
            (
                thermal.flow().lattice().grid().state_bits(),
                thermal.total_enthalpy_j().unwrap(),
                bodies[0].specific_enthalpy_j_kg(),
                gas.state().clone(),
                report.surface.unwrap().heat_from_gas_j,
            )
        })
    };
    assert_eq!(run(true), run(false));
}

struct EquilibriumOnly;
impl FreeSurfaceGasConvection for EquilibriumOnly {
    fn coefficient(
        &self,
        _: FreeSurfaceConvectionQuery<'_>,
    ) -> Result<Option<Evidence<FreeSurfaceHeatTransferCoefficient>>, FreeSurfaceEnthalpyError>
    {
        Ok(None)
    }
}

#[test]
fn undefined_coefficient_is_accepted_only_at_exact_surface_gas_equilibrium() {
    let liquid = card(0.02, 1000.0);
    let model = plate_gas();
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let before = thermal.total_enthalpy_j().unwrap();
        let warm = gas_source::point(450.0, 100000.0);
        assert!(
            thermal
                .step_with_surface_convection(
                    &[],
                    &[],
                    controls(),
                    &mut [],
                    &[],
                    &[],
                    &[FreeSurfaceInterfaceContact {
                        surface: convective_face(),
                        model: &model,
                        context: &warm
                    }],
                    FreeSurfaceInterfaceSolve::default(),
                    &EquilibriumOnly,
                    false,
                    gate,
                )
                .is_err()
        );
        assert_eq!(thermal.total_enthalpy_j().unwrap(), before);
        assert!(
            thermal
                .step_with_surface_convection(
                    &[],
                    &[],
                    controls(),
                    &mut [],
                    &[],
                    &[],
                    &[],
                    FreeSurfaceInterfaceSolve::default(),
                    &EquilibriumOnly,
                    false,
                    gate,
                )
                .is_err(),
            "a provider cannot silently vanish with an empty boundary"
        );
        let equal = gas_source::point(400.0, 100000.0);
        let report = thermal
            .step_with_surface_convection(
                &[],
                &[],
                controls(),
                &mut [],
                &[],
                &[],
                &[FreeSurfaceInterfaceContact {
                    surface: convective_face(),
                    model: &model,
                    context: &equal,
                }],
                FreeSurfaceInterfaceSolve::default(),
                &EquilibriumOnly,
                false,
                gate,
            )
            .unwrap()
            .conduction;
        let exchange = report.surface.unwrap();
        assert_eq!(exchange.heat_from_gas_j, 0.0);
        assert_eq!(
            exchange.thermal_endpoint[0].gas_heat_transfer,
            FreeSurfaceGasHeatTransfer::Convection(None)
        );
        assert_eq!(thermal.total_enthalpy_j().unwrap(), before);
    });
}

#[cfg(feature = "implicit-enthalpy-coupling")]
#[test]
fn implicit_surface_trials_reuse_the_same_convection_law_and_endpoint_heat() {
    use fs_lbm::material::free_surface::enthalpy::conduction::surface::FreeSurfaceImplicitSolve;
    let liquid = card(0.02, 1000.0);
    let model = plate_gas();
    let context = gas_source::point(450.0, 100000.0);
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let surface = [FreeSurfaceInterfaceContact {
            surface: convective_face(),
            model: &model,
            context: &context,
        }];
        let before = thermal.total_enthalpy_j().unwrap();
        let grid = thermal.flow().lattice().grid().state_bits();
        let refused = FreeSurfaceInterfaceSolve {
            implicit: Some(FreeSurfaceImplicitSolve {
                maximum_residual_evaluations: 1,
                ..FreeSurfaceImplicitSolve::default()
            }),
            ..FreeSurfaceInterfaceSolve::default()
        };
        assert!(
            thermal
                .step_with_surface_convection(
                    &[],
                    &[],
                    controls(),
                    &mut [],
                    &[],
                    &[],
                    &surface,
                    refused,
                    &Plate::default(),
                    false,
                    gate,
                )
                .is_err()
        );
        assert_eq!(thermal.total_enthalpy_j().unwrap(), before);
        assert_eq!(thermal.flow().lattice().grid().state_bits(), grid);
        let plate = Plate::default();
        let solve = FreeSurfaceInterfaceSolve {
            implicit: Some(FreeSurfaceImplicitSolve::default()),
            ..FreeSurfaceInterfaceSolve::default()
        };
        let report = thermal
            .step_with_surface_convection(
                &[],
                &[],
                FreeSurfaceConductionControls {
                    maximum_dt_s: 0.001,
                    ..controls()
                },
                &mut [],
                &[],
                &[],
                &surface,
                solve,
                &plate,
                false,
                gate,
            )
            .unwrap()
            .conduction;
        let exchange = report.surface.unwrap();
        let implicit = exchange.implicit.unwrap();
        assert_eq!(implicit.intervals, 1);
        assert!(implicit.residual_evaluations > 1);
        let endpoint = &exchange.thermal_endpoint[0];
        assert!(matches!(
            endpoint.gas_heat_transfer,
            FreeSurfaceGasHeatTransfer::Convection(Some(_))
        ));
        // Independent backward-Euler cell balance. The endpoint correlation,
        // rather than an initial-temperature h, supplies the accepted watts.
        assert!(
            (0.002 * (endpoint.liquid_temperature_k - 400.0) - 0.001 * endpoint.heat_from_gas_w)
                .abs()
                < 1e-10
        );
        assert!((exchange.heat_from_gas_j - 0.001 * endpoint.heat_from_gas_w).abs() < 1e-10);
        assert!(plate.calls.borrow().len() > implicit.residual_evaluations);
        assert!(
            (thermal.total_enthalpy_j().unwrap() - before - exchange.heat_from_gas_j).abs() < 1e-9
        );
    });
}
