//! Independent surface, storage and reservoir-entropy checks for two actual baths.
use super::super::boundary::{
    GrayAmbientBoundary, GrayAmbientControls, GrayTwoBathBoundary, ImplicitBoundaryGeometry,
};
use super::*;

const SIGMA: f64 = 5.670_374_419e-8;

fn geometry(distance: f64, resistance: f64) -> ImplicitBoundaryGeometry {
    ImplicitBoundaryGeometry {
        area_m2: 1.0,
        half_distance_m: distance,
        surface_resistance_m2_k_w: resistance,
    }
}

// The oracle uses absolute Kelvin, direct fourth powers and bisection. It does
// not use the production surface-offset coordinate or its analytic Newton step.
fn surface(t: f64, gas: f64, radiation: f64, resistance: f64, eps: f64, htc: f64) -> f64 {
    let (mut lo, mut hi) = (t.min(gas).min(radiation), t.max(gas).max(radiation));
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let flux = htc * (mid - gas) + eps * SIGMA * (mid.powi(4) - radiation.powi(4));
        if (mid - t) / resistance + flux > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

fn powers(t: f64, gas: f64, radiation: f64, eps: f64, htc: f64) -> [f64; 2] {
    let ts = surface(t, gas, radiation, 1.0, eps, htc);
    [
        htc * (ts - gas),
        eps * SIGMA * (ts.powi(4) - radiation.powi(4)),
    ]
}

fn endpoint(law: &Law, initial: f64, dt: f64, gas: f64, radiation: f64, eps: f64, htc: f64) -> f64 {
    let (mut lo, mut hi) = (100.0, 1000.0);
    for _ in 0..80 {
        let h = 0.5 * (lo + hi);
        let t = law.temperature_k(h).unwrap();
        // The independently resolved half-cell has k=A=d=1.
        let outward = t - surface(t, gas, radiation, 1.0, eps, htc);
        if h - initial + dt * outward > 0.0 {
            hi = h;
        } else {
            lo = h;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn distinct_baths_share_one_surface_and_the_original_series_drop() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (eps, htc) = (0.8, 9.0);
    for (cell, ts, gas, radiation) in [
        (600.0_f64, 500.0_f64, 300.0_f64, 400.0_f64),
        (400.0, 500.0, 300.0, 700.0),
    ] {
        let exchange =
            GrayTwoBathBoundary::new(&law, eps, htc, radiation, GrayAmbientControls::default())
                .unwrap();
        let expected = [
            htc * (ts - gas),
            eps * SIGMA * (ts.powi(4) - radiation.powi(4)),
        ];
        let resistance = (cell - ts) / (expected[0] + expected[1]);
        let mut face = geometry(0.75 * resistance, 0.25 * resistance);
        face.area_m2 = 2.0;
        let result = exchange.evaluate(cell, gas, face, &mut || Ok(())).unwrap();
        near(result.surface_temperature_k, ts, 2e-10);
        near(result.convection_power_w, 2.0 * expected[0], 1e-8);
        near(result.radiation_power_w, 2.0 * expected[1], 1e-8);
        near(
            result.outward_power_w,
            2.0 * (expected[0] + expected[1]),
            1e-8,
        );
        near(result.surface_residual_w, 0.0, 1e-8);
        assert_eq!(result.boundary.baths[0].temperature_k, gas);
        assert_eq!(result.boundary.baths[1].temperature_k, radiation);
        let slope = htc + 4.0 * eps * SIGMA * ts.powi(3);
        near(
            result.boundary.power_enthalpy_derivative_w_kg_j,
            2.0 * slope / (1.0 + resistance * slope),
            1e-10,
        );
        assert!(result.iterations < 16);
    }
}

struct AffineConductivity {
    law: Law,
    intercept: f64,
    slope: f64,
}
impl ThermalLaw for AffineConductivity {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.law.temperature_k(h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.law.specific_entropy_change_j_kg_k(a, b)
    }
}
impl ConductivityLaw for AffineConductivity {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.temperature_k(h)?;
        Ok(self.intercept + self.slope * h)
    }
}
impl ImplicitConductionLaw for AffineConductivity {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [100.0, 1000.0]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let mut state = self.law.state_at_specific_enthalpy(h)?;
        state.conductivity_w_m_k = self.conductivity_w_m_k(h)?;
        state.conductivity_slope_kg_w_m_k_j = self.slope;
        Ok(state)
    }
}

#[test]
fn latent_conductivity_feedback_and_one_bath_limit_retain_actual_enthalpy_tangent() {
    let phase = AffineConductivity {
        law: Law {
            nonlinear_k: false,
            plateau: true,
        },
        intercept: 0.0,
        slope: 1.0 / 400.0,
    };
    let sensible = Law {
        nonlinear_k: true,
        plateau: false,
    };
    for (law, h, gas, radiation) in [
        (&phase as &dyn ImplicitConductionLaw, 400.0, 200.0, 400.0),
        (&sensible, 600.0, 300.0, 450.0),
    ] {
        let exchange =
            GrayTwoBathBoundary::new(law, 0.9, 8.0, radiation, GrayAmbientControls::default())
                .unwrap();
        let evaluate = |h| {
            exchange
                .evaluate(h, gas, geometry(0.1, 0.02), &mut || Ok(()))
                .unwrap()
        };
        let current = evaluate(h);
        let dh = 1e-3;
        let difference =
            (evaluate(h + dh).outward_power_w - evaluate(h - dh).outward_power_w) / (2.0 * dh);
        near(
            current.boundary.power_enthalpy_derivative_w_kg_j,
            difference,
            2e-7,
        );
        assert!(current.boundary.power_enthalpy_derivative_w_kg_j.abs() > 1e-5);

        let single =
            GrayAmbientBoundary::new(law, 0.9, 8.0, GrayAmbientControls::default()).unwrap();
        let same =
            GrayTwoBathBoundary::new(law, 0.9, 8.0, gas, GrayAmbientControls::default()).unwrap();
        let expected = single
            .evaluate(h, gas, geometry(0.1, 0.02), &mut || Ok(()))
            .unwrap();
        let actual = same
            .evaluate(h, gas, geometry(0.1, 0.02), &mut || Ok(()))
            .unwrap();
        near(
            actual.surface_temperature_k,
            expected.surface_temperature_k,
            1e-9,
        );
        near(actual.outward_power_w, expected.outward_power_w, 1e-8);
        near(
            actual.boundary.power_enthalpy_derivative_w_kg_j,
            expected.boundary.power_enthalpy_derivative_w_kg_j,
            1e-9,
        );
    }
    assert_eq!(
        phase
            .state_at_specific_enthalpy(400.0)
            .unwrap()
            .temperature_slope_kg_k_j,
        0.0
    );
}

#[test]
fn an_insulating_bulk_retains_exterior_bath_transfer_and_zero_body_tangent() {
    let law = AffineConductivity {
        law: Law {
            nonlinear_k: false,
            plateau: false,
        },
        intercept: 0.0,
        slope: 0.0,
    };
    let (gas, radiation, eps, htc) = (300.0, 700.0, 0.8, 10.0);
    let exchange =
        GrayTwoBathBoundary::new(&law, eps, htc, radiation, GrayAmbientControls::default())
            .unwrap();
    let expected = surface(600.0, gas, radiation, f64::INFINITY, eps, htc);
    for h in [200.0, 600.0, 900.0] {
        let response = exchange
            .evaluate(h, gas, geometry(0.1, 0.02), &mut || Ok(()))
            .unwrap();
        near(response.surface_temperature_k, expected, 1e-10);
        near(response.outward_power_w, 0.0, 1e-10);
        assert!(response.convection_power_w > 0.0);
        assert!(response.radiation_power_w < 0.0);
        assert_eq!(response.boundary.power_enthalpy_derivative_w_kg_j, 0.0);
        assert!(response.convection_power_w / gas + response.radiation_power_w / radiation > 0.0);
    }
}

#[test]
fn zero_body_heat_keeps_separate_bath_entropy_in_native_be_and_sdirk() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (initial, gas, radiation, eps) = (400.0_f64, 300.0_f64, 600.0_f64, 0.8);
    let power = eps * SIGMA * (radiation.powi(4) - initial.powi(4));
    let htc = power / (initial - gas);
    let exchange =
        GrayTwoBathBoundary::new(&law, eps, htc, radiation, GrayAmbientControls::default())
            .unwrap();
    let (thermal, conduction) = setup(&law, &[initial], Some(gas), vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_boundary_exchanges(&thermal, &[(CellId(0), &exchange)], &[], || Ok(()))
            .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &conduction,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_boundary_exchanges(&thermal, &[(CellId(0), &exchange)], &[], || Ok(()))
    .unwrap();
    let dt = 1.0;
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    for result in [
        chain.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
        graph.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
        chain.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
        graph.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
    ] {
        let (next, report) = result.unwrap();
        near(next.cells()[&CellId(0)].specific_h_j_kg, initial, 1e-10);
        near(report.energy.external_heat_j, 0.0, 1e-9);
        near(report.energy.stored_energy_change_j, 0.0, 1e-9);
        near(report.energy.residual_j, 0.0, 1e-9);
        let expected_entropy = dt * power * (1.0 / radiation - 1.0 / gas);
        near(
            report.entropy.reservoir_entropy_input_j_k,
            expected_entropy,
            1e-10,
        );
        near(
            report.entropy.conduction_generation_j_k,
            -expected_entropy,
            1e-10,
        );
        near(report.entropy.temporal_defect_j_k, 0.0, 1e-10);
        assert!(!report.entropy.has_unclassified_sources);
        assert_eq!(report.energy.substeps, 1);
    }
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, initial);
}

#[test]
fn native_sensible_and_latent_storage_agree_with_independent_two_bath_backward_euler() {
    for (plateau, initial, gas, radiation, dt) in [
        (false, 600.0, 300.0, 450.0, 1.0),
        (true, 450.0, 250.0, 400.0, 0.5),
    ] {
        let law = Law {
            nonlinear_k: false,
            plateau,
        };
        let (eps, htc) = (0.8, 20.0);
        let exchange =
            GrayTwoBathBoundary::new(&law, eps, htc, radiation, GrayAmbientControls::default())
                .unwrap();
        let (thermal, conduction) = setup(&law, &[initial], Some(gas), vec![]);
        let chain = ImplicitConductionChain::new(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            || Ok(()),
        )
        .unwrap()
        .with_boundary_exchanges(&thermal, &[(CellId(0), &exchange)], &[], || Ok(()))
        .unwrap();
        let graph = ImplicitConductionGraph::new(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .unwrap()
        .with_boundary_exchanges(&thermal, &[(CellId(0), &exchange)], &[], || Ok(()))
        .unwrap();
        let expected = endpoint(&law, initial, dt, gas, radiation, eps, htc);
        let t = law.temperature_k(expected).unwrap();
        let [convection, emitted] = powers(t, gas, radiation, eps, htc);
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        for result in [
            chain.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
            graph.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
        ] {
            let (next, report) = result.unwrap();
            near(next.cells()[&CellId(0)].specific_h_j_kg, expected, 1e-9);
            near(report.energy.external_heat_j, expected - initial, 1e-9);
            near(report.energy.residual_j, 0.0, 1e-9);
            near(
                report.entropy.reservoir_entropy_input_j_k,
                -dt * (convection / gas + emitted / radiation),
                1e-10,
            );
            near(
                report.entropy.conduction_generation_j_k,
                dt * (convection * (1.0 / gas - 1.0 / t) + emitted * (1.0 / radiation - 1.0 / t)),
                1e-10,
            );
            assert!(report.entropy.conduction_generation_j_k > 0.0);
            assert!(!report.entropy.has_unclassified_sources);
            if plateau {
                assert_eq!(next.temperature_k(CellId(0)).unwrap(), 300.0);
            }
        }
    }
}

#[test]
fn both_sdirk_stages_retain_actual_bath_powers_and_entropy_weights() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (eps, htc, gas, radiation, initial, dt) = (0.6, 2.0, 300.0, 500.0, 600.0, 0.4);
    let exchange =
        GrayTwoBathBoundary::new(&law, eps, htc, radiation, GrayAmbientControls::default())
            .unwrap();
    let (thermal, conduction) = setup(&law, &[initial], Some(gas), vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_boundary_exchanges(&thermal, &[(CellId(0), &exchange)], &[], || Ok(()))
            .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &conduction,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_boundary_exchanges(&thermal, &[(CellId(0), &exchange)], &[], || Ok(()))
    .unwrap();
    let gamma = 1.0 - 1.0 / 2.0_f64.sqrt();
    let first = endpoint(&law, initial, gamma * dt, gas, radiation, eps, htc);
    let first_outward = first - surface(first, gas, radiation, 1.0, eps, htc);
    // Stage 2 still starts from the original Hn; earlier heat is an explicit
    // quadrature term, not a different physical old caloric state.
    let expected = endpoint(
        &law,
        initial - (1.0 - gamma) * dt * first_outward,
        gamma * dt,
        gas,
        radiation,
        eps,
        htc,
    );
    let mut expected_entropy = 0.0;
    let mut expected_generation = 0.0;
    for (t, weight) in [(first, 1.0 - gamma), (expected, gamma)] {
        let [convection, emitted] = powers(t, gas, radiation, eps, htc);
        expected_entropy -= dt * weight * (convection / gas + emitted / radiation);
        expected_generation += dt
            * weight
            * (convection * (1.0 / gas - 1.0 / t) + emitted * (1.0 / radiation - 1.0 / t));
    }
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    for result in [
        chain.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
        graph.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
    ] {
        let (next, report) = result.unwrap();
        near(next.cells()[&CellId(0)].specific_h_j_kg, expected, 1e-9);
        near(report.energy.external_heat_j, expected - initial, 1e-9);
        near(report.energy.residual_j, 0.0, 1e-9);
        near(
            report.entropy.reservoir_entropy_input_j_k,
            expected_entropy,
            1e-10,
        );
        near(
            report.entropy.conduction_generation_j_k,
            expected_generation,
            1e-10,
        );
        assert!(!report.entropy.has_unclassified_sources);
        assert_eq!(report.energy.substeps, 1);
    }
}

#[test]
fn original_gas_bath_events_split_exchange_and_entropy_at_the_real_event_time() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (initial, old_gas, new_gas, radiation, eps, htc) = (600.0, 300.0, 350.0, 500.0, 0.7, 5.0);
    let event_time = 0.25;
    let (thermal, conduction) = setup(
        &law,
        &[initial],
        Some(old_gas),
        vec![TimedEvent {
            time_s: event_time,
            event: Event::Boundary(CellId(0), boundary(Some(new_gas))),
        }],
    );
    let exchange =
        GrayTwoBathBoundary::new(&law, eps, htc, radiation, GrayAmbientControls::default())
            .unwrap();
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_boundary_exchanges(&thermal, &[(CellId(0), &exchange)], &[], || Ok(()))
            .unwrap();
    let first = endpoint(&law, initial, event_time, old_gas, radiation, eps, htc);
    let expected = endpoint(&law, first, event_time, new_gas, radiation, eps, htc);
    let mut expected_entropy = 0.0;
    for (t, gas) in [(first, old_gas), (expected, new_gas)] {
        let [convection, emitted] = powers(t, gas, radiation, eps, htc);
        expected_entropy -= event_time * (convection / gas + emitted / radiation);
    }
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, report) = chain
        .advance_to(&thermal, 0.5, controls(event_time), &mut workspace, || {
            Ok(())
        })
        .unwrap();
    near(next.temperature_k(CellId(0)).unwrap(), expected, 1e-9);
    near(report.energy.external_heat_j, expected - initial, 1e-9);
    near(
        report.entropy.reservoir_entropy_input_j_k,
        expected_entropy,
        1e-10,
    );
    assert_eq!(report.energy.events, 1);
    assert_eq!(report.energy.substeps, 2);
    assert!(!report.entropy.has_unclassified_sources);
    // SDIRK also splits at the real event. Both stages of each interval must
    // retain that interval's gas temperature, including separate bath entropy.
    let gamma = 1.0 - 1.0 / 2.0_f64.sqrt();
    let mut expected_sdirk = initial;
    let mut sdirk_entropy = 0.0;
    let mut sdirk_generation = 0.0;
    for gas in [old_gas, new_gas] {
        let first = endpoint(
            &law,
            expected_sdirk,
            gamma * event_time,
            gas,
            radiation,
            eps,
            htc,
        );
        let first_outward = first - surface(first, gas, radiation, 1.0, eps, htc);
        let second = endpoint(
            &law,
            expected_sdirk - (1.0 - gamma) * event_time * first_outward,
            gamma * event_time,
            gas,
            radiation,
            eps,
            htc,
        );
        for (t, weight) in [(first, 1.0 - gamma), (second, gamma)] {
            let [convection, emitted] = powers(t, gas, radiation, eps, htc);
            sdirk_entropy -= event_time * weight * (convection / gas + emitted / radiation);
            sdirk_generation += event_time
                * weight
                * (convection * (1.0 / gas - 1.0 / t) + emitted * (1.0 / radiation - 1.0 / t));
        }
        expected_sdirk = second;
    }
    let (next, report) = chain
        .advance_sdirk2_to(&thermal, 0.5, controls(event_time), &mut workspace, || {
            Ok(())
        })
        .unwrap();
    near(next.temperature_k(CellId(0)).unwrap(), expected_sdirk, 1e-9);
    near(
        report.energy.external_heat_j,
        expected_sdirk - initial,
        1e-9,
    );
    near(report.energy.residual_j, 0.0, 1e-9);
    near(
        report.entropy.reservoir_entropy_input_j_k,
        sdirk_entropy,
        1e-10,
    );
    near(
        report.entropy.conduction_generation_j_k,
        sdirk_generation,
        1e-10,
    );
    assert_eq!(report.energy.events, 1);
    assert_eq!(report.energy.substeps, 2);
    assert!(!report.entropy.has_unclassified_sources);
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, initial);
}

#[test]
fn exposed_ambient_exchanges_use_uncovered_area_and_replace_authored_htc() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let expected = [500.0_f64, 400.0_f64];
    let masses = [2.0, 3.0];
    let gas = [300.0, 350.0];
    let total_area = [3.0, 2.5];
    // The original internal contact covers 1 m² on each cell. These are the
    // remaining exterior areas, not the complete authored surface areas.
    let exposed_area = [2.0, 1.5];
    let half_distances = [0.3, 0.7];
    let coatings = [0.2, 0.1];
    let (eps, htc, radiation, dt) = (0.8, 10.0, 600.0_f64, 0.2);
    let mut bath_powers = [[0.0; 2]; 2];
    for i in 0..2 {
        let ts = surface(
            expected[i],
            gas[i],
            radiation,
            half_distances[i] + coatings[i],
            eps,
            htc,
        );
        bath_powers[i] = [
            exposed_area[i] * htc * (ts - gas[i]),
            exposed_area[i] * eps * SIGMA * (ts.powi(4) - radiation.powi(4)),
        ];
    }
    // Manufacture the initial caloric states from independent BE balances at
    // the prescribed final temperatures. k=1, contact A=1 and half-lengths
    // 0.5+0.5 give an internal conductance of exactly 1 W/K.
    let link_power = expected[0] - expected[1];
    let outward = [
        bath_powers[0][0] + bath_powers[0][1],
        bath_powers[1][0] + bath_powers[1][1],
    ];
    let initial = [
        expected[0] + dt * (outward[0] + link_power) / masses[0],
        expected[1] + dt * (outward[1] - link_power) / masses[1],
    ];
    let boundaries = gas.map(|ambient_k| Boundary {
        ambient_k,
        htc_w_m2_k: 999.0,
        source_w: 0.0,
        reservoir: None,
    });
    let mut events = Vec::new();
    for i in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                material_id: MaterialId(0),
                body_id: 0,
                segment_id: 0,
                mass_kg: masses[i],
                specific_h_j_kg: initial[i],
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: total_area[i],
                },
                boundary: boundaries[i],
            }),
        });
    }
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
        &[(MaterialId(0), &law)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 4,
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
    let conduction = ConductionObserver::new(
        &thermal,
        &[(MaterialId(0), &law)],
        &[ConductiveLink {
            link_id: 0,
            half_distances_m: [0.5; 2],
            contact_resistance_m2_k_w: 0.0,
        }],
        &[],
        || Ok(()),
    )
    .unwrap()
    .with_ambient_faces(
        &thermal,
        &[
            ConductiveAmbient {
                cell_id: CellId(0),
                half_distance_m: half_distances[0],
                surface_resistance_m2_k_w: coatings[0],
            },
            ConductiveAmbient {
                cell_id: CellId(1),
                half_distance_m: half_distances[1],
                surface_resistance_m2_k_w: coatings[1],
            },
        ],
        || Ok(()),
    )
    .unwrap();
    let exchange =
        GrayTwoBathBoundary::new(&law, eps, htc, radiation, GrayAmbientControls::default())
            .unwrap();
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_boundary_exchanges(
                &thermal,
                &[],
                &[(CellId(0), &exchange), (CellId(1), &exchange)],
                || Ok(()),
            )
            .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &conduction,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_boundary_exchanges(
        &thermal,
        &[],
        &[(CellId(0), &exchange), (CellId(1), &exchange)],
        || Ok(()),
    )
    .unwrap();
    let expected_heat = -dt * (outward[0] + outward[1]);
    let mut expected_entropy = 0.0;
    let mut expected_storage_entropy = 0.0;
    let mut expected_generation = dt * link_power * (1.0 / expected[1] - 1.0 / expected[0]);
    for i in 0..2 {
        let [convection, emitted] = bath_powers[i];
        expected_entropy -= dt * (convection / gas[i] + emitted / radiation);
        expected_storage_entropy += masses[i] * (expected[i] / initial[i]).ln();
        expected_generation += dt
            * (convection * (1.0 / gas[i] - 1.0 / expected[i])
                + emitted * (1.0 / radiation - 1.0 / expected[i]));
    }
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    for result in [
        chain.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
        graph.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
    ] {
        let (next, report) = result.unwrap();
        for i in 0..2 {
            near(
                next.temperature_k(CellId(i as u64)).unwrap(),
                expected[i],
                1e-9,
            );
            assert_eq!(next.cells()[&CellId(i as u64)].boundary, boundaries[i]);
        }
        near(report.energy.external_heat_j, expected_heat, 1e-9);
        near(report.energy.stored_energy_change_j, expected_heat, 1e-9);
        near(report.energy.residual_j, 0.0, 1e-9);
        near(
            report.entropy.reservoir_entropy_input_j_k,
            expected_entropy,
            1e-10,
        );
        near(
            report.entropy.thermal_storage_change_j_k,
            expected_storage_entropy,
            1e-10,
        );
        near(
            report.entropy.conduction_generation_j_k,
            expected_generation,
            1e-10,
        );
        near(
            report.entropy.temporal_defect_j_k,
            expected_storage_entropy - expected_entropy - expected_generation,
            1e-10,
        );
        assert!(!report.entropy.has_unclassified_sources);
        assert_eq!(report.energy.substeps, 1);
    }
}

#[test]
fn exchange_refusals_preserve_original_calorics_and_sole_channel_ownership() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let other = Law {
        nonlinear_k: false,
        plateau: false,
    };
    for radiation in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            GrayTwoBathBoundary::new(&law, 0.8, 10.0, radiation, GrayAmbientControls::default())
                .is_err()
        );
    }
    let exchange =
        GrayTwoBathBoundary::new(&law, 0.8, 10.0, 600.0, GrayAmbientControls::default()).unwrap();
    assert_eq!(
        exchange.evaluate(400.0, 300.0, geometry(0.1, 0.0), &mut || Err(
            DepositionError::Cancelled
        )),
        Err(DepositionError::Cancelled)
    );
    assert!(
        exchange
            .evaluate(99.0, 300.0, geometry(0.1, 0.0), &mut || Ok(()))
            .is_err()
    );
    let exhausted = GrayTwoBathBoundary::new(
        &law,
        0.8,
        10.0,
        600.0,
        GrayAmbientControls {
            maximum_iterations: 1,
            absolute_power_tolerance_w: 0.0,
            relative_power_tolerance: 0.0,
        },
    )
    .unwrap();
    assert!(
        exhausted
            .evaluate(400.0, 300.0, geometry(0.1, 0.0), &mut || Ok(()))
            .is_err()
    );
    let (thermal, conduction) = setup(&law, &[400.0], Some(300.0), vec![]);
    let prepare = || {
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
    };
    let legacy = GrayAmbientBoundary::new(&law, 0.8, 10.0, GrayAmbientControls::default()).unwrap();
    assert!(
        prepare()
            .with_boundary_laws(&thermal, &[(CellId(0), &legacy)], &[], || Ok(()))
            .unwrap()
            .with_boundary_exchanges(&thermal, &[(CellId(0), &exchange)], &[], || Ok(()))
            .is_err()
    );
    assert!(
        prepare()
            .with_boundary_exchanges(&thermal, &[(CellId(0), &exchange)], &[], || Ok(()))
            .unwrap()
            .with_boundary_laws(&thermal, &[(CellId(0), &legacy)], &[], || Ok(()))
            .is_err()
    );
    let foreign =
        GrayTwoBathBoundary::new(&other, 0.8, 10.0, 600.0, GrayAmbientControls::default()).unwrap();
    assert!(
        prepare()
            .with_boundary_exchanges(&thermal, &[(CellId(0), &foreign)], &[], || Ok(()))
            .is_err()
    );
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 400.0);
}
