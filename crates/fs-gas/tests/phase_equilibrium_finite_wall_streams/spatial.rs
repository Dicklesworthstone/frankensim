//! The native spatial wall remains the only wall inventory and time integrator.

use super::*;
use fs_gas::material::phase_heat_exchanger::spatial_wall::{
    MixtureStreamWallBudget, MixtureStreamWallError, MixtureStreamWallPort, MixtureStreamWallSource,
};
use fs_material::thermochemical_gas::nonideal::flash::equilibrium::acoustics::PhaseEquilibriumAcousticControls;
use fs_thermal::enthalpy::{
    DepositionBudget, DepositionSimulation, Event, StepConfig, ThermalLink, TimedEvent,
    conduction::{
        ConductionObserver, ConductiveLink, ConductivityLaw,
        implicit::{
            ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionLaw,
            ImplicitConductionState, ImplicitConductionWorkspace, source::ImplicitHeatSourceLaw,
        },
    },
};

impl ConductivityLaw for AnalyticWall {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.temperature_k(h)?;
        Ok(50.0)
    }
}
impl ImplicitConductionLaw for AnalyticWall {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.range()
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: self.temperature_k(h)?,
            temperature_slope_kg_k_j: if self.latent && (0.0..=100_000.0).contains(&h) {
                0.0
            } else {
                self.maximum_temperature_slope()
            },
            conductivity_w_m_k: self.conductivity_w_m_k(h)?,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
}

fn spatial_port(
    model: &ThermochemicalPengRobinsonModel,
    tin: f64,
    g: f64,
) -> MixtureStreamWallPort<'_> {
    let response = PhaseEquilibriumAcousticControls::default();
    let mut stream = port(model, tin, g);
    stream.controls.equilibrium.stability = response.phase.stability;
    stream.controls.equilibrium.enthalpy.flash_options = response.phase.flash;
    MixtureStreamWallPort { stream, response }
}
fn source_budget() -> MixtureStreamWallBudget {
    MixtureStreamWallBudget {
        maximum_duty_evaluations: 10_000,
        maximum_pt_evaluations: 1_000_000,
    }
}
fn native_controls(dt: f64) -> ImplicitConductionControls {
    ImplicitConductionControls {
        maximum_dt_s: dt,
        maximum_steps: 1,
        maximum_halvings: 0,
        cell_tolerance_j: 1e-3,
        total_tolerance_j: 2e-3,
        ..ImplicitConductionControls::default()
    }
}
fn slab(law: &AnalyticWall, h: [f64; 2]) -> (DepositionSimulation<'_>, ConductionObserver<'_>) {
    let mut events: Vec<_> = h
        .into_iter()
        .enumerate()
        .map(|(i, specific_h_j_kg)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(0),
                mass_kg: 2.0,
                specific_h_j_kg,
                geometry: CellGeometry {
                    start_m: [0.1 * i as f64, 0.0, 0.0],
                    end_m: [0.1 * (i + 1) as f64, 0.0, 0.0],
                    volume_m3: 0.002,
                    surface_area_m2: 2.0,
                },
                boundary: ThermalBoundary {
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
            id: 0,
            a: CellId(0),
            b: CellId(1),
            area_m2: 1.0,
            conductance_w_k: 0.0,
        }),
    });
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(0), law)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 8,
        },
    )
    .unwrap();
    let (thermal, _) = thermal
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 0.0,
            },
            || Ok(()),
        )
        .unwrap();
    let observer = ConductionObserver::new(
        &thermal,
        &[(MaterialId(0), law)],
        &[ConductiveLink {
            link_id: 0,
            half_distances_m: [0.05; 2],
            contact_resistance_m2_k_w: 0.0,
        }],
        &[],
        || Ok(()),
    )
    .unwrap();
    (thermal, observer)
}

#[test]
fn g1_original_phase_flow_cools_a_spatial_wall_through_native_conduction_and_storage() {
    let gas = source();
    let law = AnalyticWall::constant(ThermalStorageBasis::InternalEnergy);
    let initial_h = [500_000.0, 0.0]; // 400 K and 300 K at fixed reference masses.
    let (thermal, observer) = slab(&law, initial_h);
    let flowing = MixtureStreamWallSource::new(
        &law,
        &[spatial_port(&gas, 300.0, 500.0)],
        source_budget(),
        || Ok(()),
    )
    .unwrap();
    assert!(std::ptr::eq(
        ImplicitHeatSourceLaw::caloric_law(&flowing),
        ImplicitConductionLaw::caloric_law(&law)
    ));
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_heat_sources(&thermal, &[(CellId(0), &flowing)], || Ok(()))
            .unwrap();
    let dt = 20.0;
    let (next, report) = chain
        .advance_to(
            &thermal,
            dt,
            native_controls(dt),
            &mut ImplicitConductionWorkspace::new(2).unwrap(),
            || Ok(()),
        )
        .unwrap();

    // Independent 2x2 backward-Euler equations in the dilute ideal-gas limit.
    // The actual code retains finite-pressure PR departures and original PH.
    let c = 10_000.0;
    let k = 500.0; // A / (0.05 / 50 + 0.05 / 50).
    let f = effective_conductance();
    let a = c + dt * (k + f);
    let d = c + dt * k;
    let b0 = c * 400.0 + dt * f * 300.0;
    let b1 = c * 300.0;
    let determinant = a * d - (dt * k).powi(2);
    let expected = [
        (d * b0 + dt * k * b1) / determinant,
        (a * b1 + dt * k * b0) / determinant,
    ];
    let temperature = [
        next.temperature_k(CellId(0)).unwrap(),
        next.temperature_k(CellId(1)).unwrap(),
    ];
    for i in 0..2 {
        close(temperature[i], expected[i], 2e-4);
    }
    assert!(temperature[0] > temperature[1]);
    assert!(temperature[0] < 400.0 && temperature[1] > 300.0);

    let receipt = flowing
        .evaluate(next.cells()[&CellId(0)].specific_h_j_kg, || Ok(()))
        .unwrap();
    let gas_heat = dt * receipt.streams[0].enthalpy_flow_change_w;
    let wall_change: f64 = (0..2)
        .map(|i| {
            let cell = &next.cells()[&CellId(i as u64)];
            assert_eq!(cell.mass_kg, 2.0);
            assert_eq!(cell.boundary.source_w, 0.0);
            cell.mass_kg * (cell.specific_h_j_kg - initial_h[i])
        })
        .sum();
    assert_eq!(report.energy.substeps, 1);
    close(wall_change + gas_heat, 0.0, 2e-3);
    close(report.energy.external_heat_j, -gas_heat, 2e-3);
    close(report.energy.stored_energy_change_j, wall_change, 1e-8);
    assert!(report.entropy.has_unclassified_sources);
    assert!(receipt.streams[0].entropy_flow_change_w_k > 0.0);
    assert!(receipt.responses[0].is_some());
    // Native source power has not been authored into the input timeline.
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, initial_h[0]);
    assert_eq!(thermal.cells()[&CellId(1)].specific_h_j_kg, initial_h[1]);
}

#[test]
fn g1_two_mixed_ports_change_native_latent_inventory_without_a_temperature_surrogate() {
    let gas = source();
    let law = AnalyticWall {
        latent: true,
        ..AnalyticWall::constant(ThermalStorageBasis::Enthalpy)
    };
    let (thermal, observer) = slab(&law, [50_000.0; 2]);
    let flowing = MixtureStreamWallSource::new(
        &law,
        &[
            spatial_port(&gas, 300.0, 500.0),
            spatial_port(&gas, 390.0, 700.0),
        ],
        source_budget(),
        || Ok(()),
    )
    .unwrap();
    let initial = flowing.evaluate(50_000.0, || Ok(())).unwrap();
    assert!(initial.responses.iter().all(Option::is_none));
    assert_eq!(initial.source.power_enthalpy_derivative_w_kg_j, 0.0);
    assert!(initial.streams[0].enthalpy_flow_change_w > 0.0);
    assert!(initial.streams[1].enthalpy_flow_change_w < 0.0);
    close(
        initial.source.power_w,
        -initial
            .streams
            .iter()
            .map(|s| s.enthalpy_flow_change_w)
            .sum::<f64>(),
        0.0,
    );
    assert_eq!(
        flowing.work().charged_pt_evaluations,
        2 + initial
            .streams
            .iter()
            .map(|s| s.pt_evaluations)
            .sum::<usize>()
    );
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_heat_sources(&thermal, &[(CellId(0), &flowing)], || Ok(()))
            .unwrap();
    let (next, report) = chain
        .advance_to(
            &thermal,
            1.0,
            native_controls(1.0),
            &mut ImplicitConductionWorkspace::new(2).unwrap(),
            || Ok(()),
        )
        .unwrap();
    let final_h = next.cells()[&CellId(0)].specific_h_j_kg;
    assert!((0.0..100_000.0).contains(&final_h));
    assert_ne!(final_h, 50_000.0);
    assert_eq!(next.temperature_k(CellId(0)).unwrap(), 340.0);
    assert_eq!(next.temperature_k(CellId(1)).unwrap(), 340.0);
    close(2.0 * (final_h - 50_000.0), initial.source.power_w, 2e-3);
    close(next.cells()[&CellId(1)].specific_h_j_kg, 50_000.0, 1e-8);
    close(report.energy.residual_j, 0.0, 2e-3);
}

#[test]
fn g1_source_energy_tangent_matches_original_stream_trials_and_charges_all_queries() {
    let gas = source();
    let law = AnalyticWall::constant(ThermalStorageBasis::InternalEnergy);
    let declaration = spatial_port(&gas, 300.0, 500.0);
    let flowing =
        MixtureStreamWallSource::new(&law, &[declaration], source_budget(), || Ok(())).unwrap();
    let center = flowing.evaluate(250_000.0, || Ok(())).unwrap();
    let step = 50.0;
    let lower = flowing.evaluate(250_000.0 - step, || Ok(())).unwrap();
    let upper = flowing.evaluate(250_000.0 + step, || Ok(())).unwrap();
    close(
        center.source.power_enthalpy_derivative_w_kg_j,
        (upper.source.power_w - lower.source.power_w) / (2.0 * step),
        1e-6,
    );
    assert!(center.source.power_enthalpy_derivative_w_kg_j < 0.0);
    let results = [&center, &lower, &upper];
    assert_eq!(
        flowing.work().charged_duty_evaluations,
        results
            .iter()
            .map(|r| r.streams[0].duty_evaluations)
            .sum::<usize>()
    );
    assert_eq!(
        flowing.work().charged_pt_evaluations,
        1
            + results
                .iter()
                .map(|r| r.streams[0].pt_evaluations
                    + r.responses[0].as_ref().unwrap().pt_evaluations())
                .sum::<usize>()
    );

    // A whole-source budget survives repeated equal-state queries: it is not
    // recreated each time Newton or adaptive error estimation revisits h.
    let duty = center.streams[0].duty_evaluations;
    let pt = center.streams[0].pt_evaluations + 9;
    let limited = MixtureStreamWallSource::new(
        &law,
        &[declaration],
        MixtureStreamWallBudget {
            maximum_duty_evaluations: 2 * duty,
            maximum_pt_evaluations: 1 + 2 * pt,
        },
        || Ok(()),
    )
    .unwrap();
    for _ in 0..2 {
        assert_eq!(limited.evaluate(250_000.0, || Ok(())).unwrap(), center);
    }
    assert_eq!(limited.work().charged_duty_evaluations, 2 * duty);
    assert_eq!(limited.work().charged_pt_evaluations, 1 + 2 * pt);
    assert!(matches!(
        limited.evaluate(250_000.0, || Ok(())),
        Err(MixtureStreamWallError::Thermal(DepositionError::Budget(_)))
    ));
}

#[test]
fn g0_source_preserves_checkpoint_refusals_and_bounds_failed_native_trials() {
    let gas = source();
    let law = AnalyticWall::constant(ThermalStorageBasis::InternalEnergy);
    let declaration = spatial_port(&gas, 300.0, 500.0);
    let budget = MixtureStreamWallBudget {
        maximum_duty_evaluations: declaration.stream.controls.maximum_duty_evaluations,
        maximum_pt_evaluations: 1 + declaration.stream.controls.maximum_pt_evaluations,
    };
    let flowing = MixtureStreamWallSource::new(&law, &[declaration], budget, || Ok(())).unwrap();
    let expected = DepositionError::Budget("caller checkpoint");
    let mut polls = 0;
    let error = flowing
        .state_at_specific_enthalpy(250_000.0, &mut || {
            polls += 1;
            if polls == 3 {
                Err(expected.clone())
            } else {
                Ok(())
            }
        })
        .unwrap_err();
    assert_eq!(error, expected);
    assert_eq!(
        flowing.last_failure(),
        Some(MixtureStreamWallError::Thermal(expected))
    );
    assert_eq!(flowing.work().refused_stream_queries, 1);
    assert_eq!(
        flowing.work().charged_duty_evaluations,
        budget.maximum_duty_evaluations
    );
    assert_eq!(
        flowing.work().charged_pt_evaluations,
        budget.maximum_pt_evaluations
    );
    assert!(matches!(
        flowing.state_at_specific_enthalpy(250_000.0, &mut || Ok(())),
        Err(DepositionError::Budget(_))
    ));
}
