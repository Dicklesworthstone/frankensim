//! G1/G4: actual exterior-state coefficients on the original Cp/Cv stores.

use super::*;
use fs_thermal::enthalpy::conduction::implicit::boundary::{
    ImplicitGrayExteriorFluxState, ImplicitGrayExteriorLaw, ImplicitGrayExteriorState,
};

#[derive(Clone)]
struct Exterior<'a> {
    caloric: &'a dyn ThermalLaw,
    energies: [f64; 2],
    surfaces: [f64; 2],
    gases: [f64; 2],
    anchor: [f64; 3],
    values: [f64; 2],
    energy_slopes: [f64; 2],
    surface_slopes: [f64; 2],
    gas_slopes: [f64; 2],
    queries: Cell<usize>,
    last: Cell<[f64; 3]>,
    outside_surface: Cell<bool>,
}

impl<'a> Exterior<'a> {
    fn new(caloric: &'a dyn ThermalLaw, energies: [f64; 2], anchor: f64) -> Self {
        Self {
            caloric,
            energies,
            surfaces: [300.0, 700.0],
            gases: [100.0, 1100.0],
            anchor: [anchor, 500.0, 400.0],
            values: [0.55, 25.0],
            energy_slopes: [8e-7, 3e-5],
            surface_slopes: [0.0007, 0.025],
            gas_slopes: [0.0001, 0.01],
            queries: Cell::new(0),
            last: Cell::new([f64::NAN; 3]),
            outside_surface: Cell::new(false),
        }
    }

    fn coefficients(&self, e: f64, surface: f64, gas: f64) -> [f64; 2] {
        std::array::from_fn(|i| {
            self.values[i]
                + self.energy_slopes[i] * (e - self.anchor[0])
                + self.surface_slopes[i] * (surface - self.anchor[1])
                + self.gas_slopes[i] * (gas - self.anchor[2])
        })
    }
}

impl ImplicitGrayExteriorLaw for Exterior<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }

    fn energy_bounds_j_kg(&self) -> [f64; 2] {
        self.energies
    }

    fn surface_temperature_bounds_k(&self) -> [f64; 2] {
        self.surfaces
    }

    fn state_at_surface(
        &self,
        e: f64,
        surface: f64,
        gas: f64,
    ) -> Result<ImplicitGrayExteriorState, DepositionError> {
        self.queries.set(self.queries.get() + 1);
        self.last.set([e, surface, gas]);
        if !(self.surfaces[0]..=self.surfaces[1]).contains(&surface) {
            self.outside_surface.set(true);
            return Err(DepositionError::MaterialDomain);
        }
        if !(self.energies[0]..=self.energies[1]).contains(&e)
            || !(self.gases[0]..=self.gases[1]).contains(&gas)
        {
            return Err(DepositionError::MaterialDomain);
        }
        let [emissivity, convection_w_m2_k] = self.coefficients(e, surface, gas);
        Ok(ImplicitGrayExteriorState {
            coefficients: ImplicitGraySurfaceState {
                emissivity,
                convection_w_m2_k,
                emissivity_enthalpy_derivative_kg_j: self.energy_slopes[0],
                convection_enthalpy_derivative_kg_w_m2_k_j: self.energy_slopes[1],
            },
            emissivity_temperature_derivative_per_k: self.surface_slopes[0],
            convection_temperature_derivative_w_m2_k2: self.surface_slopes[1],
        })
    }
}

// Independent quadratic material/coating elimination and absolute-temperature
// bisection. The prescribed affine surface law is evaluated at EVERY candidate
// exterior temperature. Neither the production integral/inverse nor its power
// coordinate, surface root or tangent participates in this reference.
fn reference_exterior(
    side: usize,
    e: f64,
    gas: f64,
    radiation: f64,
    source: &Exterior<'_>,
    face: ImplicitBoundaryGeometry,
) -> (f64, f64, [f64; 2]) {
    let cell = temperature(side, e);
    let a = [0.04, 0.02][side] / (2.0 * face.half_distance_m);
    let p = face.surface_resistance_m2_k_w * a;
    let material = |surface: f64| {
        let rhs = surface + p * cell * cell;
        2.0 * rhs / (1.0 + (1.0 + 4.0 * p * rhs).sqrt())
    };
    let powers = |surface: f64| {
        let [eps, htc] = source.coefficients(e, surface, gas);
        [
            face.area_m2 * htc * (surface - gas),
            face.area_m2 * eps * SIGMA * (surface.powi(4) - radiation.powi(4)),
        ]
    };
    let (mut lo, mut hi) = (cell.min(gas).min(radiation), cell.max(gas).max(radiation));
    for _ in 0..100 {
        let surface = 0.5 * (lo + hi);
        let inner = material(surface);
        let q = powers(surface);
        if face.area_m2 * a * (cell * cell - inner * inner) > q[0] + q[1] {
            lo = surface;
        } else {
            hi = surface;
        }
    }
    let surface = 0.5 * (lo + hi);
    (material(surface), surface, powers(surface))
}

#[test]
fn g1_actual_exterior_sources_close_both_caloric_chains_and_coated_surface_tangents() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    for (side, layer) in layers.iter().enumerate() {
        for resistance in [0.0, 0.015] {
            let face = ImplicitBoundaryGeometry {
                surface_resistance_m2_k_w: resistance,
                ..geometry(0.08)
            };
            for [cell, gas, radiation] in [
                [620.0, 320.0, 360.0],
                [340.0, 650.0, 650.0],
                [450.0, 300.0, 650.0],
            ] {
                let source = Exterior::new(
                    layer.caloric_law(),
                    layer.energy_bounds_j_kg(),
                    energy(side, 500.0),
                );
                let boundary = IntegratedCaloricGrayBoundary::new(
                    layer,
                    face,
                    0.0,
                    0.0,
                    radiation,
                    layer_controls(),
                )
                .unwrap()
                .with_exterior_law(&source)
                .unwrap();
                assert!(std::ptr::addr_eq(
                    boundary.caloric_law(),
                    layer.caloric_law()
                ));
                let e = energy(side, cell);
                let result = boundary.evaluate(e, gas, &mut || Ok(())).unwrap();
                let (inner, surface, powers) =
                    reference_exterior(side, e, gas, radiation, &source, face);
                near(result.surface_temperature_k, surface, 2e-9);
                near(result.convection_power_w, powers[0], 2e-8);
                near(result.radiation_power_w, powers[1], 2e-8);
                near(result.outward_power_w, powers[0] + powers[1], 2e-8);
                near(result.surface_residual_w, 0.0, 2e-8);
                assert_eq!(source.last.get(), [e, result.surface_temperature_k, gas]);
                assert!(source.queries.get() > 1);
                assert!(source.queries.get() <= layer_controls().max_iterations + 2);
                assert!(!source.outside_surface.get());
                assert!((surface - cell).abs() > 1.0);

                // Differentiate the independent material/coating/surface
                // balance. Source derivatives are partial at fixed Ts and
                // gas, while 1/Cv or 1/Cp belongs to the original body only.
                let [eps, htc] = source.coefficients(e, surface, gas);
                let fourth = SIGMA * (surface.powi(4) - radiation.powi(4));
                let exterior = AREA
                    * (htc
                        + source.surface_slopes[1] * (surface - gas)
                        + 4.0 * eps * SIGMA * surface.powi(3)
                        + source.surface_slopes[0] * fourth);
                let partial = AREA
                    * (source.energy_slopes[1] * (surface - gas)
                        + source.energy_slopes[0] * fourth);
                let alpha = [0.04, 0.02][side];
                let bulk = AREA * alpha * cell / face.half_distance_m;
                let material_face = AREA * alpha * inner / face.half_distance_m;
                let capacity = [200.0 + 0.5 * cell, 400.0 + cell][side];
                let expected = (exterior * bulk / capacity + material_face * partial)
                    / (material_face + exterior + resistance / AREA * material_face * exterior);
                near(
                    result.boundary.power_enthalpy_derivative_w_kg_j,
                    expected,
                    3e-10,
                );
                let q = |energy| {
                    let (_, _, p) = reference_exterior(side, energy, gas, radiation, &source, face);
                    p[0] + p[1]
                };
                near(
                    result.boundary.power_enthalpy_derivative_w_kg_j,
                    (q(e + 0.25) - q(e - 0.25)) / 0.5,
                    3e-8,
                );
            }
        }
    }
}

#[test]
fn g4_exterior_support_is_distinct_from_bulk_and_baths_and_source_ownership_is_unique() {
    let fixture = Fixture::new(true);
    let cv = fixture.cv();
    let k = IsochoricConductionLaw::try_new(&fixture.card, &cv, &fixture.conductivity, SELECTION)
        .unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&k).unwrap();
    let e = energy(0, 500.0);
    let source = Exterior {
        surfaces: [750.0, 900.0],
        gases: [800.0, 1100.0],
        values: [0.6, 20.0],
        energy_slopes: [0.0; 2],
        surface_slopes: [0.0005, 0.01],
        gas_slopes: [0.0; 2],
        ..Exterior::new(&cv, layer.energy_bounds_j_kg(), e)
    };
    let face = ImplicitBoundaryGeometry {
        surface_resistance_m2_k_w: 0.01,
        ..geometry(0.02)
    };
    let blank = || {
        IntegratedCaloricGrayBoundary::new(&layer, face, 0.0, 0.0, 1000.0, layer_controls())
            .unwrap()
    };
    let boundary = blank().with_exterior_law(&source).unwrap();
    let response = boundary.evaluate(e, 1000.0, &mut || Ok(())).unwrap();
    let (inner, surface, powers) = reference_exterior(0, e, 1000.0, 1000.0, &source, face);
    assert!((300.0..=700.0).contains(&inner));
    assert!((750.0..=900.0).contains(&surface));
    near(response.surface_temperature_k, surface, 2e-9);
    near(response.outward_power_w, powers[0] + powers[1], 2e-8);
    assert!(!source.outside_surface.get());
    // Neither the material chart nor the exterior chart must contain a bath.
    assert!(surface > 700.0 && inner < 750.0 && 1000.0 > source.surfaces[1]);
    assert!(matches!(
        boundary.evaluate(e, 700.0, &mut || Ok(())),
        Err(DepositionError::MaterialDomain)
    ));
    let unsupported = Exterior {
        surfaces: [900.0, 950.0],
        ..source.clone()
    };
    assert!(
        blank()
            .with_exterior_law(&unsupported)
            .unwrap()
            .evaluate(e, 1000.0, &mut || Ok(()))
            .is_err()
    );
    assert!(!unsupported.outside_surface.get());

    let other = fixture.cv();
    let wrong_object = Exterior {
        caloric: &other,
        ..source.clone()
    };
    assert!(blank().with_exterior_law(&wrong_object).is_err());
    assert!(
        IntegratedCaloricGrayBoundary::new(&layer, face, 0.1, 0.0, 1000.0, layer_controls())
            .unwrap()
            .with_exterior_law(&source)
            .is_err()
    );
    let bulk = Surface {
        caloric: &cv,
        bounds: layer.energy_bounds_j_kg(),
        anchor: e,
        values: [0.6, 20.0],
        slopes: [0.0; 2],
        queries: Cell::new(0),
        last_energy: Cell::new(f64::NAN),
    };
    assert!(
        blank()
            .with_surface_law(&bulk)
            .unwrap()
            .with_exterior_law(&source)
            .is_err()
    );
    assert!(
        blank()
            .with_exterior_law(&source)
            .unwrap()
            .with_surface_law(&bulk)
            .is_err()
    );
}

#[test]
fn g1_native_exterior_evolution_keeps_original_energy_entropy_and_late_retry() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let materials: [&dyn ImplicitConductionLaw; 2] = [&a, &b];
    let face = ImplicitBoundaryGeometry {
        surface_resistance_m2_k_w: 0.025,
        ..geometry(0.04)
    };
    let (gas, radiation, dt) = (300.0, 680.0, 150.0);
    for side in 0..2 {
        let source = Exterior::new(
            layers[side].caloric_law(),
            layers[side].energy_bounds_j_kg(),
            energy(side, 570.0),
        );
        let boundary = IntegratedCaloricGrayBoundary::new(
            &layers[side],
            face,
            0.0,
            0.0,
            radiation,
            layer_controls(),
        )
        .unwrap()
        .with_exterior_law(&source)
        .unwrap();
        let thermal = single_cell(layers[side].caloric_law(), source.anchor[0], gas);
        let solver = chain(&thermal, materials[side], &boundary, face);
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        let (accepted, _) = solver
            .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        let initial = accepted.cells()[&CELL].specific_h_j_kg;
        assert!(initial < source.anchor[0]);
        let mut calls = 0;
        let (next, report) = solver
            .advance_to(&accepted, 2.0 * dt, controls(dt), &mut workspace, || {
                calls += 1;
                Ok(())
            })
            .unwrap();
        let (mut lo, mut hi) = (300.0, 700.0);
        for _ in 0..80 {
            let t = 0.5 * (lo + hi);
            let (_, _, powers) =
                reference_exterior(side, energy(side, t), gas, radiation, &source, face);
            if MASS * (energy(side, t) - initial) + dt * (powers[0] + powers[1]) > 0.0 {
                hi = t;
            } else {
                lo = t;
            }
        }
        let t = 0.5 * (lo + hi);
        let (_, _, powers) =
            reference_exterior(side, energy(side, t), gas, radiation, &source, face);
        assert!(powers[0] > 0.0 && powers[1] < 0.0);
        near(next.temperature_k(CELL).unwrap(), t, 2e-8);
        near(
            report.energy.stored_energy_change_j,
            MASS * (energy(side, t) - initial),
            2e-6,
        );
        near(
            report.energy.external_heat_j,
            -dt * (powers[0] + powers[1]),
            2e-6,
        );
        near(report.energy.residual_j, 0.0, 1e-6);
        near(
            report.entropy.reservoir_entropy_input_j_k,
            -dt * (powers[0] / gas + powers[1] / radiation),
            2e-8,
        );
        let old_t = temperature(side, initial);
        let entropy =
            MASS * ([200.0, 400.0][side] * (t / old_t).ln() + [0.5, 1.0][side] * (t - old_t));
        near(report.entropy.thermal_storage_change_j_k, entropy, 2e-8);
        near(
            report.entropy.conduction_generation_j_k,
            dt * (powers[0] * (1.0 / gas - 1.0 / t) + powers[1] * (1.0 / radiation - 1.0 / t)),
            2e-8,
        );
        assert_eq!(
            next.material_storage_basis(MATERIAL).unwrap(),
            [
                ThermalStorageBasis::InternalEnergy,
                ThermalStorageBasis::Enthalpy,
            ][side]
        );
        assert_eq!(report.energy.substeps, 1);
        assert!(!source.outside_surface.get());
        let mut cancelled = 0;
        assert!(matches!(
            solver.advance_to(&accepted, 2.0 * dt, controls(dt), &mut workspace, || {
                cancelled += 1;
                if cancelled == calls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(ImplicitConductionError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(accepted.time_s(), dt);
        assert_eq!(accepted.cells()[&CELL].specific_h_j_kg, initial);
        let (retry, retry_report) = solver
            .advance_to(&accepted, 2.0 * dt, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        assert_eq!(retry.cells(), next.cells());
        assert_eq!(retry_report, report);
    }
}

#[test]
fn g1_exterior_feedback_retains_signed_stable_tangents_and_zero_coefficient_onset() {
    let fixture = Fixture::new(true);
    let cv = fixture.cv();
    let k = IsochoricConductionLaw::try_new(&fixture.card, &cv, &fixture.conductivity, SELECTION)
        .unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&k).unwrap();
    let face = ImplicitBoundaryGeometry {
        surface_resistance_m2_k_w: 0.005,
        ..geometry(0.02)
    };
    let source = Exterior {
        surfaces: [300.0, 650.0],
        values: [0.0, 20.0],
        energy_slopes: [0.0; 2],
        surface_slopes: [0.0, -0.1],
        gas_slopes: [0.0; 2],
        ..Exterior::new(&cv, layer.energy_bounds_j_kg(), energy(0, 600.0))
    };
    let boundary =
        IntegratedCaloricGrayBoundary::new(&layer, face, 0.0, 0.0, 300.0, layer_controls())
            .unwrap()
            .with_exterior_law(&source)
            .unwrap();
    let e = energy(0, 600.0);
    let response = boundary.evaluate(e, 300.0, &mut || Ok(())).unwrap();
    let (_, surface, powers) = reference_exterior(0, e, 300.0, 300.0, &source, face);
    assert!(source.coefficients(e, surface, 300.0)[1] > 0.0);
    assert!(response.outward_power_w > 0.0);
    assert!(response.boundary.power_enthalpy_derivative_w_kg_j < 0.0);
    near(response.outward_power_w, powers[0], 2e-8);
    let q = |energy| reference_exterior(0, energy, 300.0, 300.0, &source, face).2[0];
    near(
        response.boundary.power_enthalpy_derivative_w_kg_j,
        (q(e + 0.25) - q(e - 0.25)) / 0.5,
        3e-8,
    );

    // A positive film is passive even if its feedback makes a particular
    // zero-storage surface root unstable. Here Ts=500 exactly balances the
    // k(T)=0.04T leg from Tc=600, but dQfilm/dTs=-29 W/K exceeds the
    // material restoring slope 20 W/K in magnitude. That branch must refuse.
    let unstable = Exterior {
        surfaces: [450.0, 500.0],
        values: [0.0, 550.0],
        surface_slopes: [0.0, -10.0],
        ..source.clone()
    };
    assert!(matches!(
        IntegratedCaloricGrayBoundary::new(
            &layer,
            geometry(0.02),
            0.0,
            0.0,
            300.0,
            layer_controls(),
        )
        .unwrap()
        .with_exterior_law(&unstable)
        .unwrap()
        .evaluate(e, 300.0, &mut || Ok(())),
        Err(DepositionError::InvalidInput(_))
    ));

    let e = energy(0, 500.0);
    let t = layer.state_at_specific_energy(e).unwrap().temperature_k;
    let onset = Exterior {
        energies: [e, e + 500.0],
        surfaces: [t, 600.0],
        anchor: [e, t, 400.0],
        values: [0.0; 2],
        energy_slopes: [1e-6, 0.001],
        surface_slopes: [0.001, 0.05],
        gas_slopes: [0.0; 2],
        ..source
    };
    let face = geometry(0.02);
    let boundary =
        IntegratedCaloricGrayBoundary::new(&layer, face, 0.0, 0.0, 450.0, layer_controls())
            .unwrap()
            .with_exterior_law(&onset)
            .unwrap();
    let response = boundary.evaluate(e, 400.0, &mut || Ok(())).unwrap();
    assert_eq!(response.outward_power_w, 0.0);
    near(response.surface_temperature_k, t, 1e-12);
    let fourth = SIGMA * (t.powi(4) - 450.0_f64.powi(4));
    let f_t = AREA * (onset.surface_slopes[1] * (t - 400.0) + onset.surface_slopes[0] * fourth);
    let f_e = AREA * (onset.energy_slopes[1] * (t - 400.0) + onset.energy_slopes[0] * fourth);
    let material = AREA / face.half_distance_m * 0.04 * t;
    let expected = material * (f_t / (200.0 + 0.5 * t) + f_e) / (material + f_t);
    near(
        response.boundary.power_enthalpy_derivative_w_kg_j,
        expected,
        2e-12,
    );
    near(
        boundary
            .evaluate(e + 0.1, 400.0, &mut || Ok(()))
            .unwrap()
            .outward_power_w
            / 0.1,
        expected,
        2e-7,
    );
}


// This manufactured film has a finite heat-flux derivative at equilibrium,
// although its coefficient derivative does not exist there. It tests the
// constitutive interface without introducing a new convection correlation.
struct EquilibriumFlux<'a> {
    caloric: &'a dyn ThermalLaw,
    energies: [f64; 2],
    equilibrium_queries: Cell<usize>,
}

impl ImplicitGrayExteriorLaw for EquilibriumFlux<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }

    fn energy_bounds_j_kg(&self) -> [f64; 2] {
        self.energies
    }

    fn surface_temperature_bounds_k(&self) -> [f64; 2] {
        [300.0, 700.0]
    }

    fn state_at_surface(
        &self,
        e: f64,
        surface: f64,
        gas: f64,
    ) -> Result<ImplicitGrayExteriorState, DepositionError> {
        let flux = self.flux_state_at_surface(e, surface, gas)?;
        let delta = surface - gas;
        if delta == 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(ImplicitGrayExteriorState {
            coefficients: flux.coefficients,
            emissivity_temperature_derivative_per_k: 0.0,
            convection_temperature_derivative_w_m2_k2:
                0.5 * delta.signum() * delta.abs().powf(-0.75),
        })
    }

    fn flux_state_at_surface(
        &self,
        e: f64,
        surface: f64,
        gas: f64,
    ) -> Result<ImplicitGrayExteriorFluxState, DepositionError> {
        if !(self.energies[0]..=self.energies[1]).contains(&e)
            || !(300.0..=700.0).contains(&surface)
            || !(300.0..=700.0).contains(&gas)
        {
            return Err(DepositionError::MaterialDomain);
        }
        let contrast = (surface - gas).abs().powf(0.25);
        if surface == gas {
            self.equilibrium_queries.set(self.equilibrium_queries.get() + 1);
        }
        Ok(ImplicitGrayExteriorFluxState {
            coefficients: ImplicitGraySurfaceState {
                emissivity: 0.4,
                convection_w_m2_k: 10.0 + 2.0 * contrast,
                emissivity_enthalpy_derivative_kg_j: 0.0,
                convection_enthalpy_derivative_kg_w_m2_k_j: 0.0,
            },
            emissivity_temperature_derivative_per_k: 0.0,
            convection_temperature_tangent_w_m2_k: 10.0 + 2.5 * contrast,
        })
    }
}

#[test]
fn g1_finite_flux_tangent_crosses_equilibrium_without_a_coefficient_derivative() {
    let fixture = Fixture::new(true);
    let cv = fixture.cv();
    let k = IsochoricConductionLaw::try_new(&fixture.card, &cv, &fixture.conductivity, SELECTION)
        .unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&k).unwrap();
    let source = EquilibriumFlux {
        caloric: layer.caloric_law(),
        energies: layer.energy_bounds_j_kg(),
        equilibrium_queries: Cell::new(0),
    };
    let gas = 500.0_f64;
    assert!(source.state_at_surface(energy(0, gas), gas, gas).is_err());
    let flux = source.flux_state_at_surface(energy(0, gas), gas, gas).unwrap();
    assert_eq!(flux.convection_temperature_tangent_w_m2_k, 10.0);

    for resistance in [0.0, 0.015] {
        let face = ImplicitBoundaryGeometry {
            surface_resistance_m2_k_w: resistance,
            ..geometry(0.08)
        };
        let boundary =
            IntegratedCaloricGrayBoundary::new(&layer, face, 0.0, 0.0, gas, layer_controls())
                .unwrap()
                .with_exterior_law(&source)
                .unwrap();

        // Independent temperature-coordinate reference: eliminate the
        // original quadratic k(T) material leg and coating, then bisect its
        // balance against the manufactured film and gray radiation.
        let reference = |e: f64| {
            let cell = temperature(0, e);
            let a = 0.04 / (2.0 * face.half_distance_m);
            let p = resistance * a;
            let material_face = |surface: f64| {
                let rhs = surface + p * cell * cell;
                2.0 * rhs / (1.0 + (1.0 + 4.0 * p * rhs).sqrt())
            };
            let powers = |surface: f64| {
                let delta = surface - gas;
                [
                    face.area_m2 * (10.0 + 2.0 * delta.abs().powf(0.25)) * delta,
                    face.area_m2 * 0.4 * SIGMA * (surface.powi(4) - gas.powi(4)),
                ]
            };
            let (mut lo, mut hi) = (cell.min(gas), cell.max(gas));
            for _ in 0..100 {
                let surface = 0.5 * (lo + hi);
                let inner = material_face(surface);
                let q = powers(surface);
                if face.area_m2 * a * (cell * cell - inner * inner) > q[0] + q[1] {
                    lo = surface;
                } else {
                    hi = surface;
                }
            }
            let surface = 0.5 * (lo + hi);
            (surface, powers(surface))
        };

        for cell in [480.0, 500.0, 520.0] {
            let e = energy(0, cell);
            let result = boundary.evaluate(e, gas, &mut || Ok(())).unwrap();
            let (surface, powers) = reference(e);
            near(result.surface_temperature_k, surface, 2e-9);
            near(result.convection_power_w, powers[0], 2e-8);
            near(result.radiation_power_w, powers[1], 2e-8);
            near(result.outward_power_w, powers[0] + powers[1], 2e-8);
            near(result.surface_residual_w, 0.0, 2e-8);
            if cell == gas {
                assert_eq!(result.outward_power_w, 0.0);
                // At zero contrast use the exact flux limit, not a finite
                // difference through a singular coefficient derivative.
                let material = face.area_m2 * 0.04 * cell / face.half_distance_m;
                let exterior = face.area_m2 * (10.0 + 4.0 * 0.4 * SIGMA * cell.powi(3));
                let expected = material * exterior
                    / ((200.0 + 0.5 * cell)
                        * (material + exterior
                            + resistance / face.area_m2 * material * exterior));
                near(result.boundary.power_enthalpy_derivative_w_kg_j, expected, 2e-12);
            } else {
                let (_, plus) = reference(e + 0.25);
                let (_, minus) = reference(e - 0.25);
                let expected = (plus[0] + plus[1] - minus[0] - minus[1]) / 0.5;
                near(result.boundary.power_enthalpy_derivative_w_kg_j, expected, 3e-8);
                assert_eq!(result.outward_power_w.signum(), (cell - gas).signum());
            }
        }
    }
    assert!(source.equilibrium_queries.get() > 2);
}
