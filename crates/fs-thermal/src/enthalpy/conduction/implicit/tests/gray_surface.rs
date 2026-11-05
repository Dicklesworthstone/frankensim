//! Original enthalpy surface sources, independent balances and native accounts.
use super::super::boundary::{
    GrayAmbientBoundary, GrayAmbientControls, GrayTwoBathBoundary, ImplicitBoundaryGeometry,
    ImplicitGraySurfaceLaw, ImplicitGraySurfaceState,
};
use super::*;
use std::cell::Cell;

const SIGMA: f64 = 5.670_374_419e-8;

struct Surface<'a> {
    caloric: &'a dyn ThermalLaw,
    bounds: [f64; 2],
    anchor_h: f64,
    values: [f64; 2],
    slopes: [f64; 2],
    queries: Cell<usize>,
    last_h: Cell<f64>,
}
impl<'a> Surface<'a> {
    fn new(caloric: &'a dyn ThermalLaw) -> Self {
        Self {
            caloric,
            bounds: [100.0, 1000.0],
            anchor_h: 400.0,
            values: [0.45, 7.0],
            slopes: [0.0004, 0.004],
            queries: Cell::new(0),
            last_h: Cell::new(0.0),
        }
    }
    // Declared original source law. The independent physical oracle uses these
    // values directly and never calls the production face/source adapter.
    fn coefficients(&self, h: f64) -> [f64; 2] {
        [
            self.values[0] + self.slopes[0] * (h - self.anchor_h),
            self.values[1] + self.slopes[1] * (h - self.anchor_h),
        ]
    }
}
impl ImplicitGraySurfaceLaw for Surface<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.bounds
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitGraySurfaceState, DepositionError> {
        self.queries.set(self.queries.get() + 1);
        self.last_h.set(h);
        if !(self.bounds[0]..=self.bounds[1]).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        let [emissivity, convection_w_m2_k] = self.coefficients(h);
        Ok(ImplicitGraySurfaceState {
            emissivity,
            convection_w_m2_k,
            emissivity_enthalpy_derivative_kg_j: self.slopes[0],
            convection_enthalpy_derivative_kg_w_m2_k_j: self.slopes[1],
        })
    }
}

fn face() -> ImplicitBoundaryGeometry {
    ImplicitBoundaryGeometry {
        area_m2: 1.3,
        half_distance_m: 0.2,
        surface_resistance_m2_k_w: 0.01,
    }
}
fn native_face() -> ImplicitBoundaryGeometry {
    ImplicitBoundaryGeometry {
        area_m2: 1.0,
        half_distance_m: 1.0,
        surface_resistance_m2_k_w: 0.0,
    }
}

// Independent physical solve in absolute Kelvin with direct fourth powers and
// bisection. No production surface coordinate, secant or Newton tangent.
fn surface(t: f64, gas: f64, radiation: f64, bulk: f64, coefficients: [f64; 2]) -> f64 {
    let [eps, htc] = coefficients;
    let (mut lo, mut hi) = (t.min(gas).min(radiation), t.max(gas).max(radiation));
    for _ in 0..90 {
        let ts = 0.5 * (lo + hi);
        let outward = htc * (ts - gas) + eps * SIGMA * (ts.powi(4) - radiation.powi(4));
        if bulk * (t - ts) > outward {
            lo = ts;
        } else {
            hi = ts;
        }
    }
    0.5 * (lo + hi)
}
fn reference(
    material: &dyn ImplicitConductionLaw,
    source: &Surface<'_>,
    h: f64,
    gas: f64,
    radiation: f64,
    geometry: ImplicitBoundaryGeometry,
) -> (f64, f64, [f64; 2]) {
    let state = material.state_at_specific_enthalpy(h).unwrap();
    let bulk = if state.conductivity_w_m_k == 0.0 {
        0.0
    } else {
        1.0 / (geometry.half_distance_m / state.conductivity_w_m_k
            + geometry.surface_resistance_m2_k_w)
    };
    let coefficients = source.coefficients(h);
    let ts = surface(state.temperature_k, gas, radiation, bulk, coefficients);
    let powers = [
        geometry.area_m2 * coefficients[1] * (ts - gas),
        geometry.area_m2 * coefficients[0] * SIGMA * (ts.powi(4) - radiation.powi(4)),
    ];
    (
        ts,
        geometry.area_m2 * bulk * (state.temperature_k - ts),
        powers,
    )
}
fn stage(
    material: &dyn ImplicitConductionLaw,
    source: &Surface<'_>,
    base: f64,
    dt: f64,
    gas: f64,
    radiation: f64,
) -> f64 {
    let (mut lo, mut hi) = (100.0, 1000.0);
    for _ in 0..90 {
        let h = 0.5 * (lo + hi);
        let q = reference(material, source, h, gas, radiation, native_face()).1;
        if h - base + dt * q > 0.0 {
            hi = h;
        } else {
            lo = h;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn original_surface_law_is_queried_once_and_retains_sensible_and_latent_heat_tangents() {
    for plateau in [false, true] {
        let law = Law {
            nonlinear_k: true,
            plateau,
        };
        let source = Surface::new(&law);
        let h = if plateau { 400.0 } else { 600.0 };
        let gas = 200.0;
        let ambient = GrayAmbientBoundary::new(&law, 0.0, 0.0, GrayAmbientControls::default())
            .unwrap()
            .with_surface_law(&source)
            .unwrap();
        let count = source.queries.get();
        let result = ambient.evaluate(h, gas, face(), &mut || Ok(())).unwrap();
        assert_eq!(source.queries.get(), count + 1);
        assert_eq!(source.last_h.get(), h);
        assert!(result.iterations > 1);
        let (ts, q, powers) = reference(&law, &source, h, gas, gas, face());
        near(result.surface_temperature_k, ts, 2e-10);
        near(result.outward_power_w, q, 2e-8);
        near(result.convection_power_w, powers[0], 2e-8);
        near(result.radiation_power_w, powers[1], 2e-8);
        let dh = 0.001;
        let derivative = (reference(&law, &source, h + dh, gas, gas, face()).1
            - reference(&law, &source, h - dh, gas, gas, face()).1)
            / (2.0 * dh);
        near(
            result.boundary.power_enthalpy_derivative_w_kg_j,
            derivative,
            3e-7,
        );
        for radiation in [gas, 500.0] {
            let boundary =
                GrayTwoBathBoundary::new(&law, 0.0, 0.0, radiation, GrayAmbientControls::default())
                    .unwrap()
                    .with_surface_law(&source)
                    .unwrap();
            let count = source.queries.get();
            let actual = boundary.evaluate(h, gas, face(), &mut || Ok(())).unwrap();
            assert_eq!(source.queries.get(), count + 1);
            assert_eq!(source.last_h.get(), h);
            let (ts, q, powers) = reference(&law, &source, h, gas, radiation, face());
            near(actual.surface_temperature_k, ts, 2e-10);
            near(actual.outward_power_w, q, 2e-8);
            near(actual.convection_power_w, powers[0], 2e-8);
            near(actual.radiation_power_w, powers[1], 2e-8);
            let derivative = (reference(&law, &source, h + dh, gas, radiation, face()).1
                - reference(&law, &source, h - dh, gas, radiation, face()).1)
                / (2.0 * dh);
            near(
                actual.boundary.power_enthalpy_derivative_w_kg_j,
                derivative,
                3e-7,
            );
            if radiation == gas {
                near(actual.outward_power_w, result.outward_power_w, 2e-8);
                near(
                    actual.boundary.power_enthalpy_derivative_w_kg_j,
                    result.boundary.power_enthalpy_derivative_w_kg_j,
                    2e-9,
                );
            }
            if plateau {
                let state = law.state_at_specific_enthalpy(h).unwrap();
                assert_eq!(state.temperature_slope_kg_k_j, 0.0);
                assert_eq!(state.conductivity_slope_kg_w_m_k_j, 0.0);
                // Both ordinary caloric/k terms vanish. This is surface-source
                // phase response, not apparent capacity or a numerical slope.
                assert!(actual.boundary.power_enthalpy_derivative_w_kg_j.abs() > 1e-4);
            }
        }
    }
}

#[test]
fn zero_emissivity_and_film_onset_keep_the_original_inward_enthalpy_derivative() {
    let law = Law {
        nonlinear_k: false,
        plateau: true,
    };
    let source = Surface {
        bounds: [400.0, 500.0],
        values: [0.0; 2],
        slopes: [0.003, 0.02],
        ..Surface::new(&law)
    };
    let h = 400.0;
    let t = 300.0_f64;
    let gas = 200.0_f64;
    let ambient = GrayAmbientBoundary::new(&law, 0.0, 0.0, GrayAmbientControls::default())
        .unwrap()
        .with_surface_law(&source)
        .unwrap();
    let result = ambient.evaluate(h, gas, face(), &mut || Ok(())).unwrap();
    assert_eq!(result.outward_power_w, 0.0);
    assert_eq!(result.surface_temperature_k, t);
    let expected = face().area_m2
        * (source.slopes[1] * (t - gas) + source.slopes[0] * SIGMA * (t.powi(4) - gas.powi(4)));
    near(
        result.boundary.power_enthalpy_derivative_w_kg_j,
        expected,
        1e-12,
    );
    let dh = 1e-4;
    near(
        ambient
            .evaluate(h + dh, gas, face(), &mut || Ok(()))
            .unwrap()
            .outward_power_w
            / dh,
        expected,
        3e-5,
    );

    let radiation = 500.0_f64;
    let separate =
        GrayTwoBathBoundary::new(&law, 0.0, 0.0, radiation, GrayAmbientControls::default())
            .unwrap()
            .with_surface_law(&source)
            .unwrap();
    let result = separate.evaluate(h, gas, face(), &mut || Ok(())).unwrap();
    assert_eq!(result.outward_power_w, 0.0);
    assert_eq!(result.surface_temperature_k, t);
    let expected = face().area_m2
        * (source.slopes[1] * (t - gas)
            + source.slopes[0] * SIGMA * (t.powi(4) - radiation.powi(4)));
    near(
        result.boundary.power_enthalpy_derivative_w_kg_j,
        expected,
        1e-12,
    );
    near(
        separate
            .evaluate(h + dh, gas, face(), &mut || Ok(()))
            .unwrap()
            .outward_power_w
            / dh,
        expected,
        4e-5,
    );
}

#[test]
fn native_be_and_sdirk_use_current_surface_sources_and_actual_bath_entropy() {
    for (two_baths, plateau) in [(false, false), (false, true), (true, false), (true, true)] {
        let law = Law {
            nonlinear_k: true,
            plateau,
        };
        let source = Surface::new(&law);
        let (initial, gas, dt) = (if plateau { 400.0 } else { 600.0 }, 200.0, 0.1);
        let radiation = if two_baths { 500.0 } else { gas };
        let (thermal, conduction) = setup(&law, &[initial], Some(gas), vec![]);
        let common = GrayAmbientBoundary::new(&law, 0.0, 0.0, GrayAmbientControls::default())
            .unwrap()
            .with_surface_law(&source)
            .unwrap();
        let distinct =
            GrayTwoBathBoundary::new(&law, 0.0, 0.0, radiation, GrayAmbientControls::default())
                .unwrap()
                .with_surface_law(&source)
                .unwrap();
        let chain = ImplicitConductionChain::new(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            || Ok(()),
        )
        .unwrap();
        let chain = if two_baths {
            chain.with_boundary_exchanges(&thermal, &[(CellId(0), &distinct)], &[], || Ok(()))
        } else {
            chain.with_boundary_laws(&thermal, &[(CellId(0), &common)], &[], || Ok(()))
        }
        .unwrap();
        let graph = ImplicitConductionGraph::new(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .unwrap();
        let graph = if two_baths {
            graph.with_boundary_exchanges(&thermal, &[(CellId(0), &distinct)], &[], || Ok(()))
        } else {
            graph.with_boundary_laws(&thermal, &[(CellId(0), &common)], &[], || Ok(()))
        }
        .unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        for sdirk in [false, true] {
            let (expected, stages) = if sdirk {
                let gamma = 1.0 - 1.0 / 2.0_f64.sqrt();
                let first = stage(&law, &source, initial, gamma * dt, gas, radiation);
                let first_q = reference(&law, &source, first, gas, radiation, native_face()).1;
                let final_h = stage(
                    &law,
                    &source,
                    initial - (1.0 - gamma) * dt * first_q,
                    gamma * dt,
                    gas,
                    radiation,
                );
                (final_h, [(first, 1.0 - gamma), (final_h, gamma)])
            } else {
                let final_h = stage(&law, &source, initial, dt, gas, radiation);
                (final_h, [(final_h, 1.0), (final_h, 0.0)])
            };
            let (mut entropy, mut generation) = (0.0, 0.0);
            for (h, weight) in stages {
                let (_, _, [conv, rad]) =
                    reference(&law, &source, h, gas, radiation, native_face());
                let t = law.temperature_k(h).unwrap();
                entropy -= dt * weight * (conv / gas + rad / radiation);
                generation += dt
                    * weight
                    * (conv * (1.0 / gas - 1.0 / t) + rad * (1.0 / radiation - 1.0 / t));
            }
            let results = if sdirk {
                [
                    chain.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
                    graph.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
                ]
            } else {
                [
                    chain.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
                    graph.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
                ]
            };
            for result in results {
                let (next, report) = result.unwrap();
                near(next.cells()[&CellId(0)].specific_h_j_kg, expected, 2e-8);
                near(report.energy.external_heat_j, expected - initial, 2e-8);
                near(report.energy.residual_j, 0.0, 1e-9);
                near(report.entropy.reservoir_entropy_input_j_k, entropy, 2e-10);
                near(report.entropy.conduction_generation_j_k, generation, 2e-10);
                assert!(!report.entropy.has_unclassified_sources);
                assert_eq!(report.energy.substeps, 1);
                if plateau {
                    assert_eq!(next.temperature_k(CellId(0)).unwrap(), 300.0);
                }
            }
        }
        assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, initial);
    }
}

struct Insulator<'a>(&'a Law);
impl ConductivityLaw for Insulator<'_> {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.0.temperature_k(h)?;
        Ok(0.0)
    }
}
impl ImplicitConductionLaw for Insulator<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.0
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
        state.conductivity_w_m_k = 0.0;
        state.conductivity_slope_kg_w_m_k_j = 0.0;
        Ok(state)
    }
}

#[test]
fn sourced_ambient_counterflow_keeps_actual_free_area_and_reservoir_entropy() {
    let law = Law {
        nonlinear_k: true,
        plateau: true,
    };
    let source = Surface::new(&law);
    let (h, t, gas, dt) = (400.0, 300.0_f64, 200.0_f64, 0.25);
    let [eps, htc] = source.coefficients(h);
    let radiation = (t.powi(4) + htc * (t - gas) / (eps * SIGMA)).sqrt().sqrt();
    let bath = Boundary {
        ambient_k: gas,
        htc_w_m2_k: 999.0,
        source_w: 0.0,
        reservoir: None,
    };
    let (thermal, conduction) = setup(
        &law,
        &[h, h],
        None,
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Boundary(CellId(0), bath),
        }],
    );
    // Each original cell has total area 2 and its contact covers 1. The first
    // cell's actual free area is exactly 1; authored htc=999 must be replaced.
    let conduction = conduction
        .with_ambient_faces(
            &thermal,
            &[ConductiveAmbient {
                cell_id: CellId(0),
                half_distance_m: 1.0,
                surface_resistance_m2_k_w: 0.0,
            }],
            || Ok(()),
        )
        .unwrap();
    let boundary =
        GrayTwoBathBoundary::new(&law, 0.0, 0.0, radiation, GrayAmbientControls::default())
            .unwrap()
            .with_surface_law(&source)
            .unwrap();
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_boundary_exchanges(&thermal, &[], &[(CellId(0), &boundary)], || Ok(()))
            .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &conduction,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_boundary_exchanges(&thermal, &[], &[(CellId(0), &boundary)], || Ok(()))
    .unwrap();
    let expected_generation = dt * htc * (t - gas) * (1.0 / gas - 1.0 / radiation);
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    for result in [
        chain.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
        graph.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
    ] {
        let (next, report) = result.unwrap();
        for id in [CellId(0), CellId(1)] {
            near(next.cells()[&id].specific_h_j_kg, h, 1e-10);
        }
        near(report.energy.external_heat_j, 0.0, 1e-9);
        near(report.entropy.thermal_storage_change_j_k, 0.0, 1e-12);
        near(
            report.entropy.reservoir_entropy_input_j_k,
            -expected_generation,
            1e-10,
        );
        near(
            report.entropy.conduction_generation_j_k,
            expected_generation,
            1e-10,
        );
        assert!(report.entropy.conduction_generation_j_k > 0.0);
        assert_eq!(next.cells()[&CellId(0)].boundary.htc_w_m2_k, 999.0);
    }
    // With k=0 the exterior still transfers heat between the two baths; the
    // original body derivative remains zero despite nonzero source derivatives.
    let insulating = Insulator(&law);
    let isolated =
        GrayTwoBathBoundary::new(&insulating, 0.0, 0.0, 500.0, GrayAmbientControls::default())
            .unwrap()
            .with_surface_law(&source)
            .unwrap();
    let result = isolated.evaluate(h, gas, face(), &mut || Ok(())).unwrap();
    let (ts, _, powers) = reference(&insulating, &source, h, gas, 500.0, face());
    near(result.surface_temperature_k, ts, 2e-10);
    near(result.convection_power_w, powers[0], 2e-8);
    near(result.radiation_power_w, powers[1], 2e-8);
    near(result.outward_power_w, 0.0, 1e-9);
    assert_eq!(result.boundary.power_enthalpy_derivative_w_kg_j, 0.0);
    assert!(result.convection_power_w > 0.0 && result.radiation_power_w < 0.0);
}

#[test]
fn surface_sources_preserve_fixed_limits_and_refuse_identity_domain_ownership_and_cancellation() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let constant = Surface {
        slopes: [0.0; 2],
        ..Surface::new(&law)
    };
    let sourced = GrayAmbientBoundary::new(&law, 0.0, 0.0, GrayAmbientControls::default())
        .unwrap()
        .with_surface_law(&constant)
        .unwrap();
    let fixed = GrayAmbientBoundary::new(
        &law,
        constant.values[0],
        constant.values[1],
        GrayAmbientControls::default(),
    )
    .unwrap();
    assert_eq!(
        sourced
            .evaluate(600.0, 200.0, face(), &mut || Ok(()))
            .unwrap(),
        fixed
            .evaluate(600.0, 200.0, face(), &mut || Ok(()))
            .unwrap()
    );
    let separate = GrayTwoBathBoundary::new(&law, 0.0, 0.0, 500.0, GrayAmbientControls::default())
        .unwrap()
        .with_surface_law(&constant)
        .unwrap();
    let fixed = GrayTwoBathBoundary::new(
        &law,
        constant.values[0],
        constant.values[1],
        500.0,
        GrayAmbientControls::default(),
    )
    .unwrap();
    assert_eq!(
        separate
            .evaluate(600.0, 200.0, face(), &mut || Ok(()))
            .unwrap(),
        fixed
            .evaluate(600.0, 200.0, face(), &mut || Ok(()))
            .unwrap()
    );
    assert!(sourced.with_surface_law(&constant).is_err());
    assert!(separate.with_surface_law(&constant).is_err());
    assert!(fixed.with_surface_law(&constant).is_err());
    let alien = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let wrong = Surface::new(&alien);
    assert!(
        GrayAmbientBoundary::new(&law, 0.0, 0.0, GrayAmbientControls::default())
            .unwrap()
            .with_surface_law(&wrong)
            .is_err()
    );
    let invalid = Surface {
        bounds: [500.0, 400.0],
        ..Surface::new(&law)
    };
    assert!(
        GrayTwoBathBoundary::new(&law, 0.0, 0.0, 500.0, GrayAmbientControls::default())
            .unwrap()
            .with_surface_law(&invalid)
            .is_err()
    );
    for (values, slopes) in [
        ([1.1, 7.0], [0.0; 2]),
        ([0.4, -1.0], [0.0; 2]),
        ([0.4, 7.0], [f64::NAN, 0.0]),
    ] {
        let bad = Surface {
            values,
            slopes,
            ..Surface::new(&law)
        };
        let boundary = GrayAmbientBoundary::new(&law, 0.0, 0.0, GrayAmbientControls::default())
            .unwrap()
            .with_surface_law(&bad)
            .unwrap();
        assert_eq!(
            boundary.evaluate(400.0, 200.0, face(), &mut || Ok(())),
            Err(DepositionError::MaterialDomain)
        );
    }
    let narrow = Surface {
        bounds: [400.0, 500.0],
        ..Surface::new(&law)
    };
    let insulating = Insulator(&law);
    let boundary =
        GrayTwoBathBoundary::new(&insulating, 0.0, 0.0, 500.0, GrayAmbientControls::default())
            .unwrap()
            .with_surface_law(&narrow)
            .unwrap();
    assert_eq!(
        boundary.evaluate(600.0, 200.0, face(), &mut || Ok(())),
        Err(DepositionError::MaterialDomain)
    );
    assert_eq!(narrow.queries.get(), 0);
    let source = Surface::new(&law);
    let boundary = GrayTwoBathBoundary::new(&law, 0.0, 0.0, 500.0, GrayAmbientControls::default())
        .unwrap()
        .with_surface_law(&source)
        .unwrap();
    let mut calls = 0;
    assert_eq!(
        boundary.evaluate(600.0, 200.0, face(), &mut || {
            calls += 1;
            if calls == 3 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    );
    assert_eq!(source.queries.get(), 1);
    let (thermal, conduction) = setup(&law, &[600.0], Some(200.0), vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_boundary_exchanges(&thermal, &[(CellId(0), &boundary)], &[], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    assert!(
        chain
            .advance_to(&thermal, 0.1, controls(0.1), &mut workspace, || Err(
                DepositionError::Cancelled
            ))
            .is_err()
    );
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 600.0);
}
