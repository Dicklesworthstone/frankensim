//! G1/G3/G4: source-integrated Cp/Cv surfaces and the original two-bath accounts.

use super::*;
use fs_thermal::enthalpy::Reservoir;
use fs_thermal::enthalpy::conduction::ConductiveReservoir;
use fs_thermal::enthalpy::conduction::implicit::boundary::{
    ImplicitBoundaryExchangeLaw, ImplicitBoundaryGeometry, ImplicitGraySurfaceLaw,
    ImplicitGraySurfaceState,
};
use std::cell::Cell;

const SIGMA: f64 = 5.670_374_419e-8;
const CELL: CellId = CellId(0);
const MATERIAL: MaterialId = MaterialId(11);
const MASS: f64 = 0.5;

fn geometry(distance: f64) -> ImplicitBoundaryGeometry {
    ImplicitBoundaryGeometry {
        area_m2: AREA,
        half_distance_m: distance,
        surface_resistance_m2_k_w: 0.0,
    }
}

// Independent physical oracle: the quadratic primitive for k(T)=a*T, direct
// fourth powers and absolute-Kelvin bisection. No production source, integral,
// surface coordinate or Newton derivative is used.
fn reference(
    side: usize,
    cell: f64,
    gas: f64,
    radiation: f64,
    coefficients: [f64; 2],
    face: ImplicitBoundaryGeometry,
) -> (f64, [f64; 2]) {
    let [eps, htc] = coefficients;
    let primitive = [0.04, 0.02][side] / (2.0 * face.half_distance_m);
    let (mut lo, mut hi) = (cell.min(gas).min(radiation), cell.max(gas).max(radiation));
    for _ in 0..90 {
        let surface = 0.5 * (lo + hi);
        let outward = htc * (surface - gas) + eps * SIGMA * (surface.powi(4) - radiation.powi(4));
        if primitive * (cell * cell - surface * surface) > outward {
            lo = surface;
        } else {
            hi = surface;
        }
    }
    let surface = 0.5 * (lo + hi);
    (
        surface,
        [
            face.area_m2 * htc * (surface - gas),
            face.area_m2 * eps * SIGMA * (surface.powi(4) - radiation.powi(4)),
        ],
    )
}

fn controls(dt: f64) -> ImplicitConductionControls {
    ImplicitConductionControls {
        maximum_dt_s: dt,
        maximum_halvings: 0,
        cell_tolerance_j: 1e-6,
        total_tolerance_j: 1e-6,
        ..ImplicitConductionControls::default()
    }
}

fn single_cell(law: &dyn ThermalLaw, initial: f64, gas: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELL,
                segment_id: 1,
                body_id: 1,
                material_id: MATERIAL,
                mass_kg: MASS,
                specific_h_j_kg: initial,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [MASS / (RHO * AREA), 0.0, 0.0],
                    volume_m3: MASS / RHO,
                    surface_area_m2: 3.0 * AREA,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: Some(Reservoir {
                        area_m2: AREA,
                        conductance_w_k: 999.0,
                        temperature_k: gas,
                    }),
                },
            }),
        }],
        &[(MATERIAL, law)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn chain<'a>(
    thermal: &DepositionSimulation<'_>,
    material: &'a dyn ImplicitConductionLaw,
    exchange: &'a dyn ImplicitBoundaryExchangeLaw,
    face: ImplicitBoundaryGeometry,
) -> ImplicitConductionChain<'a> {
    let observer = ConductionObserver::new(
        thermal,
        &[(MATERIAL, material.conductivity_law())],
        &[],
        &[ConductiveReservoir {
            cell_id: CELL,
            half_distance_m: face.half_distance_m,
            surface_resistance_m2_k_w: face.surface_resistance_m2_k_w,
        }],
        || Ok(()),
    )
    .unwrap();
    ImplicitConductionChain::new(thermal, &observer, &[(MATERIAL, material)], || Ok(()))
        .unwrap()
        .with_boundary_exchanges(thermal, &[(CELL, exchange)], &[], || Ok(()))
        .unwrap()
}

struct Surface<'a> {
    caloric: &'a dyn ThermalLaw,
    bounds: [f64; 2],
    anchor: f64,
    values: [f64; 2],
    slopes: [f64; 2],
    queries: Cell<usize>,
    last_energy: Cell<f64>,
}

impl Surface<'_> {
    fn coefficients(&self, e: f64) -> [f64; 2] {
        [
            self.values[0] + self.slopes[0] * (e - self.anchor),
            self.values[1] + self.slopes[1] * (e - self.anchor),
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
        e: f64,
    ) -> Result<ImplicitGraySurfaceState, DepositionError> {
        self.queries.set(self.queries.get() + 1);
        self.last_energy.set(e);
        if !(self.bounds[0]..=self.bounds[1]).contains(&e) {
            return Err(DepositionError::MaterialDomain);
        }
        let [emissivity, convection_w_m2_k] = self.coefficients(e);
        Ok(ImplicitGraySurfaceState {
            emissivity,
            convection_w_m2_k,
            emissivity_enthalpy_derivative_kg_j: self.slopes[0],
            convection_enthalpy_derivative_kg_w_m2_k_j: self.slopes[1],
        })
    }
}

#[test]
fn g1_original_cp_cv_gray_surface_matches_quartic_balance_and_complete_tangents() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let face = geometry(0.12);
    let coefficients = [0.8, 20.0];
    for (side, layer) in layers.iter().enumerate() {
        for [cell, gas, radiation] in [
            [620.0, 320.0, 360.0],
            [340.0, 650.0, 650.0],
            [450.0, 300.0, 650.0],
            [500.0, 500.0, 500.0],
        ] {
            let boundary = IntegratedCaloricGrayBoundary::new(
                layer,
                face,
                coefficients[0],
                coefficients[1],
                radiation,
                layer_controls(),
            )
            .unwrap();
            assert!(std::ptr::addr_eq(
                boundary.caloric_law(),
                layer.caloric_law()
            ));
            let e = energy(side, cell);
            let result = boundary.evaluate(e, gas, &mut || Ok(())).unwrap();
            let (surface, powers) = reference(side, cell, gas, radiation, coefficients, face);
            near(result.surface_temperature_k, surface, 2e-9);
            near(result.convection_power_w, powers[0], 2e-8);
            near(result.radiation_power_w, powers[1], 2e-8);
            near(result.outward_power_w, powers[0] + powers[1], 2e-8);
            near(result.surface_residual_w, 0.0, 2e-8);
            assert_eq!(result.boundary.baths[0].temperature_k, gas);
            assert_eq!(result.boundary.baths[1].temperature_k, radiation);
            let q = |e| {
                let (_, p) = reference(
                    side,
                    temperature(side, e),
                    gas,
                    radiation,
                    coefficients,
                    face,
                );
                p[0] + p[1]
            };
            near(
                result.boundary.power_enthalpy_derivative_w_kg_j,
                (q(e + 0.25) - q(e - 0.25)) / 0.5,
                2e-8,
            );
            // Differentiate the independent integral equation directly. The
            // interior derivative uses k(Tcell); the surface derivative uses
            // k(Tsurface), rather than freezing either value across the layer.
            let alpha = [0.04, 0.02][side];
            let capacity = [200.0 + 0.5 * cell, 400.0 + cell][side];
            let exterior = coefficients[1] + 4.0 * coefficients[0] * SIGMA * surface.powi(3);
            let expected = face.area_m2 * exterior * alpha * cell
                / (capacity * (alpha * surface + face.half_distance_m * exterior));
            near(
                result.boundary.power_enthalpy_derivative_w_kg_j,
                expected,
                2e-10,
            );
        }
    }
}

#[test]
fn g1_bulk_source_is_queried_once_and_retains_signed_and_zero_onset_derivatives() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let (face, gas, radiation) = (geometry(0.12), 320.0, 400.0);
    for (side, layer) in layers.iter().enumerate() {
        let e = energy(side, 620.0);
        let source = Surface {
            caloric: layer.caloric_law(),
            bounds: [e - 500.0, e + 500.0],
            anchor: e,
            values: [0.5, 20.0],
            slopes: [-0.0001, -0.01],
            queries: Cell::new(0),
            last_energy: Cell::new(f64::NAN),
        };
        let boundary =
            IntegratedCaloricGrayBoundary::new(layer, face, 0.0, 0.0, radiation, layer_controls())
                .unwrap()
                .with_surface_law(&source)
                .unwrap();
        let before = source.queries.get();
        let result = boundary.evaluate(e, gas, &mut || Ok(())).unwrap();
        assert_eq!(source.queries.get(), before + 1);
        assert_eq!(source.last_energy.get(), e);
        let q = |e| {
            let (_, p) = reference(
                side,
                temperature(side, e),
                gas,
                radiation,
                source.coefficients(e),
                face,
            );
            p[0] + p[1]
        };
        near(result.outward_power_w, q(e), 2e-8);
        assert!(result.outward_power_w > 0.0);
        assert!(result.boundary.power_enthalpy_derivative_w_kg_j < 0.0);
        near(
            result.boundary.power_enthalpy_derivative_w_kg_j,
            (q(e + 0.25) - q(e - 0.25)) / 0.5,
            3e-8,
        );
        assert!(
            boundary
                .evaluate(source.bounds[1] + 1.0, gas, &mut || Ok(()))
                .is_err()
        );

        let e = energy(side, 500.0);
        let onset = Surface {
            bounds: [e, e + 1000.0],
            anchor: e,
            values: [0.0; 2],
            slopes: [0.00001, 0.02],
            ..source
        };
        let boundary =
            IntegratedCaloricGrayBoundary::new(layer, face, 0.0, 0.0, radiation, layer_controls())
                .unwrap()
                .with_surface_law(&onset)
                .unwrap();
        let result = boundary.evaluate(e, gas, &mut || Ok(())).unwrap();
        assert_eq!(result.outward_power_w, 0.0);
        near(result.surface_temperature_k, 500.0, 1e-12);
        let expected = AREA
            * (onset.slopes[1] * (500.0 - gas)
                + onset.slopes[0] * SIGMA * (500.0_f64.powi(4) - radiation.powi(4)));
        near(
            result.boundary.power_enthalpy_derivative_w_kg_j,
            expected,
            2e-12,
        );
        near(
            boundary
                .evaluate(e + 0.01, gas, &mut || Ok(()))
                .unwrap()
                .outward_power_w
                / 0.01,
            expected,
            2e-6,
        );
    }
}

#[test]
fn g1_native_gray_evolution_and_neutral_body_keep_original_storage_and_bath_entropy() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let material: [&dyn ImplicitConductionLaw; 2] = [&a, &b];
    let (face, gas, radiation, dt) = (geometry(0.04), 300.0, 680.0, 150.0);
    let coefficients = [0.8, 40.0];
    for side in 0..2 {
        let thermal = single_cell(layers[side].caloric_law(), energy(side, 570.0), gas);
        let boundary = IntegratedCaloricGrayBoundary::new(
            &layers[side],
            face,
            coefficients[0],
            coefficients[1],
            radiation,
            layer_controls(),
        )
        .unwrap();
        let solver = chain(&thermal, material[side], &boundary, face);
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        let (accepted, _) = solver
            .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        assert!(accepted.temperature_k(CELL).unwrap() < 570.0);
        let initial = accepted.cells()[&CELL].specific_h_j_kg;
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
            let (_, p) = reference(side, t, gas, radiation, coefficients, face);
            if MASS * (energy(side, t) - initial) + dt * (p[0] + p[1]) > 0.0 {
                hi = t;
            } else {
                lo = t;
            }
        }
        let t = 0.5 * (lo + hi);
        let (_, p) = reference(side, t, gas, radiation, coefficients, face);
        assert!(p[0] > 0.0 && p[1] < 0.0);
        near(next.temperature_k(CELL).unwrap(), t, 2e-8);
        near(report.energy.external_heat_j, -dt * (p[0] + p[1]), 2e-6);
        near(
            report.energy.stored_energy_change_j,
            MASS * (energy(side, t) - initial),
            2e-6,
        );
        near(report.energy.residual_j, 0.0, 1e-6);
        near(
            report.entropy.reservoir_entropy_input_j_k,
            -dt * (p[0] / gas + p[1] / radiation),
            2e-8,
        );
        let old_t = temperature(side, initial);
        let storage =
            MASS * ([200.0, 400.0][side] * (t / old_t).ln() + [0.5, 1.0][side] * (t - old_t));
        near(report.entropy.thermal_storage_change_j_k, storage, 2e-8);
        near(
            report.entropy.conduction_generation_j_k,
            dt * (p[0] * (1.0 / gas - 1.0 / t) + p[1] * (1.0 / radiation - 1.0 / t)),
            2e-8,
        );
        assert!(!report.entropy.has_unclassified_sources);
        assert_eq!(
            next.cells()[&CELL].boundary,
            thermal.cells()[&CELL].boundary
        );
        assert_eq!(
            next.material_storage_basis(MATERIAL).unwrap(),
            [
                ThermalStorageBasis::InternalEnergy,
                ThermalStorageBasis::Enthalpy
            ][side]
        );
        assert_eq!(report.energy.substeps, 1);
        let mut cancelled_calls = 0;
        assert!(matches!(
            solver.advance_to(&accepted, 2.0 * dt, controls(dt), &mut workspace, || {
                cancelled_calls += 1;
                if cancelled_calls == calls {
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

        // Exact opposing bath powers leave the body unchanged, but still
        // produce bath entropy. Check the same fact through native SDIRK2.
        let (t, radiation) = (400.0_f64, 600.0_f64);
        let radiation_flux = coefficients[0] * SIGMA * (radiation.powi(4) - t.powi(4));
        let neutral = single_cell(layers[side].caloric_law(), energy(side, t), gas);
        let boundary = IntegratedCaloricGrayBoundary::new(
            &layers[side],
            face,
            coefficients[0],
            radiation_flux / (t - gas),
            radiation,
            layer_controls(),
        )
        .unwrap();
        let solver = chain(&neutral, material[side], &boundary, face);
        let (next, report) = solver
            .advance_sdirk2_to(&neutral, 10.0, controls(10.0), &mut workspace, || Ok(()))
            .unwrap();
        near(next.temperature_k(CELL).unwrap(), t, 2e-9);
        near(report.energy.external_heat_j, 0.0, 1e-7);
        near(report.energy.stored_energy_change_j, 0.0, 1e-7);
        let bath_entropy = 10.0 * AREA * radiation_flux * (1.0 / radiation - 1.0 / gas);
        near(
            report.entropy.reservoir_entropy_input_j_k,
            bath_entropy,
            2e-9,
        );
        near(
            report.entropy.conduction_generation_j_k,
            -bath_entropy,
            2e-9,
        );
    }
}

#[test]
fn g4_external_baths_need_only_supported_material_faces_and_failed_steps_preserve_history() {
    let fixture = Fixture::new(true);
    let cv = fixture.cv();
    let k = IsochoricConductionLaw::try_new(&fixture.card, &cv, &fixture.conductivity, SELECTION)
        .unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&k).unwrap();
    let face = geometry(0.02);
    let boundary =
        IntegratedCaloricGrayBoundary::new(&layer, face, 0.1, 20.0, 1000.0, layer_controls())
            .unwrap();
    let response = boundary
        .evaluate(energy(0, 500.0), 1000.0, &mut || Ok(()))
        .unwrap();
    let expected = reference(0, 500.0, 1000.0, 1000.0, [0.1, 20.0], face);
    near(response.surface_temperature_k, expected.0, 2e-9);
    assert!(response.outward_power_w < 0.0 && response.surface_temperature_k < 700.0);
    assert!(
        boundary
            .evaluate(energy(0, 695.0), 1000.0, &mut || Ok(()))
            .is_err()
    );
    assert!(
        boundary
            .state_at_specific_enthalpy(energy(0, 500.0), 1000.0, geometry(0.01), &mut || Ok(()))
            .is_err()
    );
    assert!(
        IntegratedCaloricGrayBoundary::new(
            &layer,
            ImplicitBoundaryGeometry {
                surface_resistance_m2_k_w: 0.01,
                ..face
            },
            0.1,
            20.0,
            1000.0,
            layer_controls(),
        )
        .is_err()
    );
    let other = fixture.cv();
    let source = Surface {
        caloric: &other,
        bounds: layer.energy_bounds_j_kg(),
        anchor: energy(0, 500.0),
        values: [0.1, 20.0],
        slopes: [0.0; 2],
        queries: Cell::new(0),
        last_energy: Cell::new(0.0),
    };
    assert!(
        IntegratedCaloricGrayBoundary::new(&layer, face, 0.0, 0.0, 1000.0, layer_controls())
            .unwrap()
            .with_surface_law(&source)
            .is_err()
    );
    let unsupported = Surface {
        caloric: &cv,
        bounds: [
            layer.energy_bounds_j_kg()[0],
            layer.energy_bounds_j_kg()[1] + 1.0,
        ],
        ..source
    };
    assert!(
        IntegratedCaloricGrayBoundary::new(&layer, face, 0.0, 0.0, 1000.0, layer_controls())
            .unwrap()
            .with_surface_law(&unsupported)
            .is_err()
    );

    let thermal = single_cell(&cv, energy(0, 680.0), 1000.0);
    let solver = chain(&thermal, &k, &boundary, face);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (accepted, _) = solver
        .advance_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    let initial = accepted.cells()[&CELL].specific_h_j_kg;
    assert!(accepted.temperature_k(CELL).unwrap() > 680.0);
    assert!(
        solver
            .advance_to(&accepted, 1001.0, controls(1000.0), &mut workspace, || Ok(
                ()
            ))
            .is_err()
    );
    assert_eq!(accepted.time_s(), 1.0);
    assert_eq!(accepted.cells()[&CELL].specific_h_j_kg, initial);
    let (retry, _) = solver
        .advance_to(&accepted, 2.0, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    assert!(retry.temperature_k(CELL).unwrap() > accepted.temperature_k(CELL).unwrap());
    assert!(retry.temperature_k(CELL).unwrap() < 700.0);
}
