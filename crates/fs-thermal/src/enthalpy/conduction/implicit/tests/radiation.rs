//! G1/G3: original enthalpy storage with a resolved nonlinear exterior face.
use super::super::adaptive::AdaptiveImplicitConductionControls;
use super::super::boundary::{GrayAmbientBoundary, GrayAmbientControls, ImplicitBoundaryGeometry};
use super::*;

const SIGMA: f64 = 5.670_374_419e-8;

fn geometry(distance: f64, resistance: f64) -> ImplicitBoundaryGeometry {
    ImplicitBoundaryGeometry {
        area_m2: 1.0,
        half_distance_m: distance,
        surface_resistance_m2_k_w: resistance,
    }
}

// Independent oracle uses absolute fourth powers and bisection on temperature,
// rather than the production secant/fraction/Newton surface coordinates.
fn surface_temperature(t: f64, bath: f64, resistance: f64, eps: f64, htc: f64) -> f64 {
    let (mut lo, mut hi) = (t.min(bath), t.max(bath));
    for _ in 0..70 {
        let mid = 0.5 * (lo + hi);
        let flux = htc * (mid - bath) + eps * SIGMA * (mid.powi(4) - bath.powi(4));
        if (mid - t) / resistance + flux > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn gray_surface_resolves_radiation_convection_series_drop_and_analytic_tangent() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (eps, htc) = (0.8, 9.0);
    let boundary =
        GrayAmbientBoundary::new(&law, eps, htc, GrayAmbientControls::default()).unwrap();
    for (cell, surface, bath) in [(600.0_f64, 500.0_f64, 300.0_f64), (400.0, 500.0, 700.0)] {
        let convection = htc * (surface - bath);
        let radiation = eps * SIGMA * (surface.powi(4) - bath.powi(4));
        let total = convection + radiation;
        let r = (cell - surface) / total;
        let mut face = geometry(0.75 * r, 0.25 * r);
        face.area_m2 = 2.0;
        let response = boundary.evaluate(cell, bath, face, &mut || Ok(())).unwrap();
        near(response.surface_temperature_k, surface, 2e-10);
        near(response.convection_power_w, 2.0 * convection, 1e-8);
        near(response.radiation_power_w, 2.0 * radiation, 1e-8);
        near(response.outward_power_w, 2.0 * total, 1e-8);
        near(response.surface_residual_w, 0.0, 1e-8);
        let slope = htc + 4.0 * eps * SIGMA * surface.powi(3);
        near(
            response.boundary.power_enthalpy_derivative_w_kg_j,
            2.0 * slope / (1.0 + r * slope),
            1e-10,
        );
        assert!(response.iterations < 16);
    }
    let at_equilibrium = boundary
        .evaluate(400.0, 400.0, geometry(0.1, 0.02), &mut || Ok(()))
        .unwrap();
    let slope = htc + 4.0 * eps * SIGMA * 400.0_f64.powi(3);
    assert_eq!(at_equilibrium.outward_power_w, 0.0);
    near(
        at_equilibrium.boundary.conductance_w_k,
        slope / (1.0 + 0.12 * slope),
        1e-13,
    );
    near(
        at_equilibrium.boundary.power_enthalpy_derivative_w_kg_j,
        at_equilibrium.boundary.conductance_w_k,
        1e-13,
    );
}

struct PhaseConductivity(Law);
impl ThermalLaw for PhaseConductivity {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.0.temperature_k(h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.0.specific_entropy_change_j_kg_k(a, b)
    }
}
impl ConductivityLaw for PhaseConductivity {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.temperature_k(h)?;
        Ok(h / 400.0)
    }
}
impl ImplicitConductionLaw for PhaseConductivity {
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
        let mut state = self.0.state_at_specific_enthalpy(h)?;
        state.conductivity_w_m_k = self.conductivity_w_m_k(h)?;
        state.conductivity_slope_kg_w_m_k_j = 1.0 / 400.0;
        Ok(state)
    }
}

#[test]
fn latent_plateau_keeps_conductivity_feedback_in_the_surface_enthalpy_tangent() {
    let phase = PhaseConductivity(Law {
        nonlinear_k: false,
        plateau: true,
    });
    let sensible = Law {
        nonlinear_k: true,
        plateau: false,
    };
    for (law, h, bath) in [
        (&phase as &dyn ImplicitConductionLaw, 400.0, 200.0),
        (&sensible, 600.0, 300.0),
    ] {
        let boundary =
            GrayAmbientBoundary::new(law, 0.9, 8.0, GrayAmbientControls::default()).unwrap();
        let face = geometry(0.1, 0.02);
        let evaluate = |h| boundary.evaluate(h, bath, face, &mut || Ok(())).unwrap();
        let current = evaluate(h);
        let step = 1e-3;
        let derivative = (evaluate(h + step).outward_power_w - evaluate(h - step).outward_power_w)
            / (2.0 * step);
        near(
            current.boundary.power_enthalpy_derivative_w_kg_j,
            derivative,
            2e-8,
        );
        assert!(current.boundary.power_enthalpy_derivative_w_kg_j > 0.0);
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
fn native_chain_and_graph_radiative_cooling_store_sensible_and_latent_heat_once() {
    for (plateau, initial, bath, dt) in [(false, 600.0, 300.0, 2.0), (true, 450.0, 250.0, 0.5)] {
        let law = Law {
            nonlinear_k: false,
            plateau,
        };
        let boundary =
            GrayAmbientBoundary::new(&law, 0.8, 20.0, GrayAmbientControls::default()).unwrap();
        let (thermal, conduction) = setup(&law, &[initial], Some(bath), vec![]);
        let chain = ImplicitConductionChain::new(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            || Ok(()),
        )
        .unwrap()
        .with_boundary_laws(&thermal, &[(CellId(0), &boundary)], &[], || Ok(()))
        .unwrap();
        let graph = ImplicitConductionGraph::new(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .unwrap()
        .with_boundary_laws(&thermal, &[(CellId(0), &boundary)], &[], || Ok(()))
        .unwrap();
        let (mut lo, mut hi) = (bath, initial);
        for _ in 0..80 {
            let h = 0.5 * (lo + hi);
            let t = law.temperature_k(h).unwrap();
            let surface = surface_temperature(t, bath, 1.0, 0.8, 20.0);
            let q = t - surface;
            if h - initial + dt * q > 0.0 {
                hi = h;
            } else {
                lo = h;
            }
        }
        let expected = 0.5 * (lo + hi);
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        for result in [
            chain.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
            graph.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
        ] {
            let (next, report) = result.unwrap();
            near(next.cells()[&CellId(0)].specific_h_j_kg, expected, 1e-9);
            near(report.energy.external_heat_j, expected - initial, 1e-9);
            near(
                report.energy.stored_energy_change_j,
                expected - initial,
                1e-9,
            );
            near(report.energy.residual_j, 0.0, 1e-9);
            near(
                report.entropy.reservoir_entropy_input_j_k,
                (expected - initial) / bath,
                1e-11,
            );
            assert!(report.entropy.conduction_generation_j_k > 0.0);
            assert_eq!(report.energy.substeps, 1);
            if plateau {
                assert_eq!(next.temperature_k(CellId(0)).unwrap(), 300.0);
            }
        }
        assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, initial);
    }
}

#[test]
fn adaptive_radiative_cooling_refines_to_independent_continuous_calorimetry() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (eps, htc, bath, initial, final_temperature) = (0.6, 2.0, 300.0, 600.0, 450.0);
    let boundary =
        GrayAmbientBoundary::new(&law, eps, htc, GrayAmbientControls::default()).unwrap();
    let (thermal, conduction) = setup(&law, &[initial], Some(bath), vec![]);
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &conduction,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_boundary_laws(&thermal, &[(CellId(0), &boundary)], &[], || Ok(()))
    .unwrap();
    // For constant C=A=R=1, Tc=Ts+f(Ts), so elapsed time is
    // integral (1+f'(Ts))/f(Ts) dTs. Simpson quadrature is independent of
    // both the native time stepper and its surface-temperature solve.
    let from = surface_temperature(final_temperature, bath, 1.0, eps, htc);
    let to = surface_temperature(initial, bath, 1.0, eps, htc);
    let integrand = |t: f64| {
        (1.0 + htc + 4.0 * eps * SIGMA * t.powi(3))
            / (htc * (t - bath) + eps * SIGMA * (t.powi(4) - bath.powi(4)))
    };
    let intervals = 2048;
    let spacing = (to - from) / f64::from(intervals);
    let mut sum = integrand(from) + integrand(to);
    for i in 1..intervals {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * integrand(from + spacing * f64::from(i));
    }
    let elapsed = spacing * sum / 3.0;
    let mut errors = Vec::new();
    for tolerance in [0.1, 1e-3] {
        let config = AdaptiveImplicitConductionControls {
            step: controls(elapsed),
            initial_dt_s: elapsed,
            maximum_trials: 10_000,
            enthalpy_absolute_tolerance_j_kg: tolerance,
            temperature_absolute_tolerance_k: tolerance,
            phase_fraction_absolute_tolerance: None,
            relative_tolerance: 0.0,
        };
        let (next, report) = graph
            .advance_adaptive_to(
                &thermal,
                elapsed,
                config,
                &mut ImplicitConductionWorkspace::new(1).unwrap(),
                || Ok(()),
            )
            .unwrap();
        errors.push((next.temperature_k(CellId(0)).unwrap() - final_temperature).abs());
        assert_eq!(next.time_s(), elapsed);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        assert!(report.accuracy_rejections > 0);
        near(
            report.solve.energy.external_heat_j,
            next.cells()[&CellId(0)].specific_h_j_kg - initial,
            1e-8,
        );
        near(report.solve.energy.residual_j, 0.0, 1e-8);
    }
    assert!(errors[1] < 0.2 * errors[0], "{errors:?}");
    assert!(errors[1] < 0.3, "{errors:?}");
}

#[test]
fn gray_boundary_refuses_invalid_data_failed_surface_solves_and_cancellation() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    for eps in [-0.1, 1.1, f64::NAN] {
        assert!(GrayAmbientBoundary::new(&law, eps, 0.0, GrayAmbientControls::default()).is_err());
    }
    let boundary =
        GrayAmbientBoundary::new(&law, 0.8, 10.0, GrayAmbientControls::default()).unwrap();
    assert_eq!(
        boundary.evaluate(600.0, 300.0, geometry(0.1, 0.0), &mut || Err(
            DepositionError::Cancelled
        )),
        Err(DepositionError::Cancelled)
    );
    assert!(
        boundary
            .evaluate(99.0, 300.0, geometry(0.1, 0.0), &mut || Ok(()))
            .is_err()
    );
    assert!(
        boundary
            .evaluate(600.0, 0.0, geometry(0.1, 0.0), &mut || Ok(()))
            .is_err()
    );
    let exhausted = GrayAmbientBoundary::new(
        &law,
        0.8,
        10.0,
        GrayAmbientControls {
            maximum_iterations: 1,
            absolute_power_tolerance_w: 0.0,
            relative_power_tolerance: 0.0,
        },
    )
    .unwrap();
    assert!(
        exhausted
            .evaluate(600.0, 300.0, geometry(0.1, 0.0), &mut || Ok(()))
            .is_err()
    );
    let insulated =
        GrayAmbientBoundary::new(&law, 0.0, 0.0, GrayAmbientControls::default()).unwrap();
    let response = insulated
        .evaluate(600.0, 300.0, geometry(0.1, 0.0), &mut || Ok(()))
        .unwrap();
    assert_eq!(response.outward_power_w, 0.0);
    assert_eq!(response.surface_temperature_k, 600.0);
    assert_eq!(response.boundary.power_enthalpy_derivative_w_kg_j, 0.0);
}
