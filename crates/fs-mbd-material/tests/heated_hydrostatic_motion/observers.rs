//! Same native buoyancy step with generated heat and existing spatial transport.
use super::*;
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveReservoir, ConductivityLaw};
use fs_thermal::observers::{CoupledObservers, ObserverCouplingError};
use std::collections::BTreeMap;

const HEAT_DT: f64 = DT / 2.0;

/// Synthetic source fixture with independent history and a tighter time bound.
#[derive(Clone, Debug, PartialEq)]
struct GeneratedHeat {
    watts: f64,
    time_s: f64,
    supplied_j: f64,
    refuse_after_s: Option<f64>,
}
impl GeneratedHeat {
    fn new(watts: f64) -> Self {
        Self {
            watts,
            time_s: 0.0,
            supplied_j: 0.0,
            refuse_after_s: None,
        }
    }
}
impl DepositionObserver for GeneratedHeat {
    type Error = DepositionError;

    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }

    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        Ok(BTreeMap::from([(CellId(0), self.watts)]))
    }

    fn limit_timestep_s(
        &mut self,
        _: &DepositionSimulation<'_>,
        proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        checkpoint()?;
        Ok(proposed_s.min(HEAT_DT))
    }

    fn step(
        &mut self,
        start_s: f64,
        dt_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        assert_eq!(self.time_s, start_s);
        assert!(dt_s <= HEAT_DT);
        self.time_s = thermal.time_s();
        self.supplied_j += self.watts * dt_s;
        if self.refuse_after_s.is_some_and(|time| self.time_s > time) {
            return Err(DepositionError::InvalidInput(
                "generated heat endpoint refused",
            ));
        }
        Ok(())
    }
}

fn layered<'a>(
    surface: &'a FluidTransportState,
    lower: &'a FluidTransportState,
) -> LayeredHydrostaticMedium<'a> {
    LayeredHydrostaticMedium::new(
        None,
        &[
            HydrostaticLayer {
                depth_m: 0.0,
                source: surface,
            },
            HydrostaticLayer {
                depth_m: 0.5,
                source: lower,
            },
        ],
        gravity(),
        Vec3::default(),
        || false,
    )
    .unwrap()
}

fn velocity(body: &HeatedRigidBody<'_>) -> f64 {
    body.state()
        .center_of_mass_velocity_world(body.principal_frame().properties())
        .unwrap()
        .z
}

fn second_step_velocity(first_heat_j: f64) -> f64 {
    // Independent rectangular-column integrals at both native Verlet endpoints
    // after the first heat update. The initial cube is neutral at z=0.
    let stretch = 1.0 + 0.01 * first_heat_j / 20.0;
    let depth = stretch / 2.0;
    let displaced_mass = stretch.powi(2) * (4.0 * depth + 4.0 * (depth - 0.5));
    let acceleration_before = 10.0 * displaced_mass / 2.0 - 10.0;
    let endpoint_height = 0.5 * HEAT_DT * HEAT_DT * acceleration_before;
    let endpoint_depth = depth - endpoint_height;
    let endpoint_mass =
        stretch.powi(2) * (4.0 * endpoint_depth + 4.0 * (endpoint_depth - 0.5).max(0.0));
    let acceleration_after = 10.0 * endpoint_mass / 2.0 - 10.0;
    0.5 * HEAT_DT * (acceleration_before + acceleration_after)
}

#[test]
fn generated_heat_changes_layered_lift_on_its_own_substeps_and_restarts_exactly() {
    let material = card(10.0, 0.01, 500.0);
    let source = reference(&material, 1);
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let body = moving(
        heated(&source, &material, &chart, 0.0, vec![]),
        0.0,
        controls(),
    );
    let parts = [envelope("a", 0.0)];
    let surface = fluid(4.0);
    let lower = fluid(8.0);
    let medium = layered(&surface, &lower);
    let input = GeneratedHeat::new(40.0);
    let (half, half_heat, half_audit) = body
        .advance_hydrostatic_with(&input, &medium, &parts, DT, config(), || Ok(()))
        .unwrap();

    close(velocity(&half), second_step_velocity(40.0 * HEAT_DT));
    close(half.thermal().energy_j(), 40.0 * DT);
    close(half_heat.supplied_j, half_audit.thermal.external_heat_j);
    close(half_heat.time_s, half.thermal().time_s());
    close(
        half.configuration().inertia().mass(),
        source.inertia().mass(),
    );
    close(half_audit.combined_residual_j, 0.0);
    assert!(half_audit.displacement_parameter_work_j > 0.0);

    let (whole, whole_heat, whole_audit) = body
        .advance_hydrostatic_with(&input, &medium, &parts, 2.0 * DT, config(), || Ok(()))
        .unwrap();
    let (split, split_heat, split_audit) = half
        .advance_hydrostatic_with(&half_heat, &medium, &parts, 2.0 * DT, config(), || Ok(()))
        .unwrap();
    assert_eq!(whole.state(), split.state());
    assert_eq!(whole.thermal().cells(), split.thermal().cells());
    assert_eq!(whole_heat, split_heat);
    close(
        whole_audit.parameter_work_j,
        half_audit.parameter_work_j + split_audit.parameter_work_j,
    );
    close(
        whole_audit.thermal.external_heat_j,
        half_audit.thermal.external_heat_j + split_audit.thermal.external_heat_j,
    );
    close(whole_audit.combined_residual_j, 0.0);
    assert_eq!(input, GeneratedHeat::new(40.0));
    close(body.thermal().energy_j(), 0.0);
    close(velocity(&body), 0.0);
}

struct ConstantConductivity;
impl ConductivityLaw for ConstantConductivity {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(1.0)
    }
}

#[test]
fn native_conduction_and_generated_heat_are_integrated_once_before_buoyant_motion() {
    let material = card(10.0, 0.01, 500.0);
    let source = reference(&material, 1);
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let mut thermal_region = region(&material, &chart, "a", 0, 0.0);
    thermal_region.boundary.reservoir = Some(Reservoir {
        area_m2: 1.0,
        conductance_w_k: 0.0,
        temperature_k: 320.0,
    });
    let body = moving(
        HeatedBody::new(
            &source,
            vec![thermal_region],
            vec![],
            budget(),
            config(),
            || Ok(()),
        )
        .unwrap(),
        0.0,
        controls(),
    );
    let conductivity = ConstantConductivity;
    let conduction = ConductionObserver::new(
        body.thermal(),
        &[(MaterialId(10), &conductivity)],
        &[],
        &[ConductiveReservoir {
            cell_id: CellId(0),
            half_distance_m: 0.5,
            surface_resistance_m2_k_w: 0.0,
        }],
        || Ok(()),
    )
    .unwrap();
    let input = CoupledObservers::new(GeneratedHeat::new(20.0), conduction);
    let parts = [envelope("a", 0.0)];
    let surface = fluid(4.0);
    let lower = fluid(8.0);
    let medium = layered(&surface, &lower);
    let (next, accepted, audit) = body
        .advance_hydrostatic_with(&input, &medium, &parts, DT, config(), || Ok(()))
        .unwrap();

    // The production conduction owner supplies G=A*k/d=2 W/K.
    // Its second power sample sees the first interval's actual temperature.
    let first_heat = (20.0 + 2.0 * (320.0 - 300.0)) * HEAT_DT;
    let first_temperature = 300.0 + first_heat / 20.0;
    let second_heat = (20.0 + 2.0 * (320.0 - first_temperature)) * HEAT_DT;
    close(next.thermal().energy_j(), first_heat + second_heat);
    close(audit.thermal.external_heat_j, first_heat + second_heat);
    close(accepted.first.supplied_j, 20.0 * DT);
    close(velocity(&next), second_step_velocity(first_heat));
    close(audit.combined_residual_j, 0.0);
    close(input.first.supplied_j, 0.0);
    close(body.thermal().energy_j(), 0.0);
}

#[test]
fn late_heat_motion_and_cancellation_refusals_leave_all_original_histories_intact() {
    let material = card(10.0, 0.01, 500.0);
    let source = reference(&material, 1);
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let thermal = heated(&source, &material, &chart, 0.0, vec![]);
    let body = moving(thermal.clone(), 0.0, controls());
    let parts = [envelope("a", 0.0)];
    let surface = fluid(4.0);
    let lower = fluid(8.0);
    let medium = layered(&surface, &lower);
    let input = GeneratedHeat::new(40.0);
    let refusing = GeneratedHeat {
        refuse_after_s: Some(HEAT_DT),
        ..input.clone()
    };
    assert!(matches!(
        body.advance_hydrostatic_with(&refusing, &medium, &parts, DT, config(), || Ok(())),
        Err(ObserverCouplingError::First(DepositionError::InvalidInput(
            "generated heat endpoint refused"
        )))
    ));

    let mut tight = controls();
    tight.transfer.maximum_parameter_work_j = 1e-8;
    let limited = moving(thermal, 0.0, tight);
    assert!(matches!(
        limited.advance_hydrostatic_with(&input, &medium, &parts, DT, config(), || Ok(())),
        Err(ObserverCouplingError::Second(
            fs_mbd_material::thermal::heating::motion::HeatedMotionError::Budget(_)
        ))
    ));

    let mut calls = 0;
    body.advance_hydrostatic_with(&input, &medium, &parts, DT, config(), || {
        calls += 1;
        Ok(())
    })
    .unwrap();
    let mut visited = 0;
    let refused = body.advance_hydrostatic_with(&input, &medium, &parts, DT, config(), || {
        visited += 1;
        if visited == calls {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(refused.is_err());
    assert_eq!(visited, calls);
    assert_eq!(input, GeneratedHeat::new(40.0));
    close(refusing.time_s, 0.0);
    close(refusing.supplied_j, 0.0);
    close(body.thermal().time_s(), 0.0);
    close(limited.thermal().energy_j(), 0.0);
    assert_eq!(limited.state(), body.state());
}
