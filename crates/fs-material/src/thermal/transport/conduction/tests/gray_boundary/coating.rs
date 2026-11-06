//! G1/G4: original Cp/Cv storage behind a finite zero-storage surface coating.

use super::*;
use crate::thermal::conduction::{ConductivityCurve, ReservoirLayerGeometry};

// Independently eliminate the quadratic material/coating leg at a prescribed
// exterior temperature, then bisect its quartic bath balance. The material
// face obeys p*Tm²+Tm=Ts+p*Tc², p=R''*alpha/(2d), for k(T)=alpha*T.
// No production integral, inverse, surface coordinate or flux solver is used.
fn reference_coated(
    side: usize,
    cell: f64,
    gas: f64,
    radiation: f64,
    coefficients: [f64; 2],
    face: ImplicitBoundaryGeometry,
) -> (f64, f64, [f64; 2]) {
    let [eps, htc] = coefficients;
    let a = [0.04, 0.02][side] / (2.0 * face.half_distance_m);
    let p = face.surface_resistance_m2_k_w * a;
    let material = |surface: f64| {
        let rhs = surface + p * cell * cell;
        2.0 * rhs / (1.0 + (1.0 + 4.0 * p * rhs).sqrt())
    };
    let (mut lo, mut hi) = (cell.min(gas).min(radiation), cell.max(gas).max(radiation));
    for _ in 0..90 {
        let surface = 0.5 * (lo + hi);
        let inner = material(surface);
        let outward = htc * (surface - gas) + eps * SIGMA * (surface.powi(4) - radiation.powi(4));
        if a * (cell * cell - inner * inner) > outward {
            lo = surface;
        } else {
            hi = surface;
        }
    }
    let surface = 0.5 * (lo + hi);
    (
        material(surface),
        surface,
        [
            face.area_m2 * htc * (surface - gas),
            face.area_m2 * eps * SIGMA * (surface.powi(4) - radiation.powi(4)),
        ],
    )
}

#[test]
fn g1_coated_faces_preserve_original_caloric_tangents_and_series_resistor_limits() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let face = ImplicitBoundaryGeometry {
        surface_resistance_m2_k_w: 0.015,
        ..geometry(0.08)
    };
    for (side, layer) in layers.iter().enumerate() {
        for (temperatures, values, slopes) in [
            ([620.0, 320.0, 360.0], [0.8, 20.0], [0.0; 2]),
            ([340.0, 650.0, 650.0], [0.8, 20.0], [0.0; 2]),
            ([450.0, 300.0, 650.0], [0.8, 20.0], [0.0; 2]),
            ([500.0, 500.0, 500.0], [0.8, 20.0], [0.0; 2]),
            ([620.0, 320.0, 400.0], [0.5, 20.0], [-0.0001, -0.01]),
        ] {
            let [cell, gas, radiation] = temperatures;
            let e = energy(side, cell);
            let source = Surface {
                caloric: layer.caloric_law(),
                bounds: [e - 500.0, e + 500.0],
                anchor: e,
                values,
                slopes,
                queries: Cell::new(0),
                last_energy: Cell::new(f64::NAN),
            };
            let boundary = IntegratedCaloricGrayBoundary::new(
                layer,
                face,
                0.0,
                0.0,
                radiation,
                layer_controls(),
            )
            .unwrap()
            .with_surface_law(&source)
            .unwrap();
            let result = boundary.evaluate(e, gas, &mut || Ok(())).unwrap();
            assert_eq!(source.queries.get(), 1);
            let (inner, surface, powers) =
                reference_coated(side, cell, gas, radiation, values, face);
            near(result.surface_temperature_k, surface, 2e-9);
            near(result.convection_power_w, powers[0], 2e-8);
            near(result.radiation_power_w, powers[1], 2e-8);
            near(result.outward_power_w, powers[0] + powers[1], 2e-8);
            near(result.surface_residual_w, 0.0, 2e-8);
            let resistance = face.surface_resistance_m2_k_w / face.area_m2;
            near(
                result.surface_temperature_k + resistance * result.outward_power_w,
                inner,
                2e-8,
            );
            let alpha = [0.04, 0.02][side];
            let capacity = [200.0 + 0.5 * cell, 400.0 + cell][side];
            let bulk = AREA / face.half_distance_m * alpha * cell;
            let material_face = AREA / face.half_distance_m * alpha * inner;
            let exterior = AREA * (values[1] + 4.0 * values[0] * SIGMA * surface.powi(3));
            let source_response = AREA
                * (slopes[1] * (surface - gas)
                    + slopes[0] * SIGMA * (surface.powi(4) - radiation.powi(4)));
            let expected = (exterior * bulk / capacity + material_face * source_response)
                / (material_face + exterior + resistance * material_face * exterior);
            near(
                result.boundary.power_enthalpy_derivative_w_kg_j,
                expected,
                3e-10,
            );
            let q = |e| {
                let (_, _, p) = reference_coated(
                    side,
                    temperature(side, e),
                    gas,
                    radiation,
                    source.coefficients(e),
                    face,
                );
                p[0] + p[1]
            };
            near(
                result.boundary.power_enthalpy_derivative_w_kg_j,
                (q(e + 0.25) - q(e - 0.25)) / 0.5,
                3e-8,
            );
            if slopes == [0.0; 2] {
                let fixed = IntegratedCaloricGrayBoundary::new(
                    layer,
                    face,
                    values[0],
                    values[1],
                    radiation,
                    layer_controls(),
                )
                .unwrap();
                assert_eq!(fixed.evaluate(e, gas, &mut || Ok(())).unwrap(), result);
            } else {
                assert!(result.outward_power_w > 0.0);
                assert!(result.boundary.power_enthalpy_derivative_w_kg_j < 0.0);
            }
        }
    }

    // Independent three-node resistor network checks the generic kernel's
    // secants, which the native gray adapter does not itself consume.
    let curve = ConductivityCurve::new(vec![[300.0, 10.0], [700.0, 10.0]]).unwrap();
    let leg = ReservoirLayerGeometry::planar(0.02, 0.02).unwrap();
    let film = |surface: f64| {
        Ok::<_, DepositionError>((
            [3.0 * (surface - 300.0), 5.0 * (surface - 700.0)],
            [3.0, 5.0],
            surface.to_bits(),
        ))
    };
    let result = leg
        .evaluate_material_to_two_baths_with_resistance(
            &curve,
            500.0,
            [300.0, 700.0],
            [300.0, 700.0],
            0.4,
            layer_controls(),
            film,
            || Ok(()),
        )
        .unwrap();
    near(result.material_face_temperature_k, 508.0, 2e-9);
    near(result.surface_temperature_k, 540.0, 2e-9);
    near(result.solid_power_to_surface_w, -80.0, 2e-8);
    for side in 0..2 {
        near(result.powers_to_baths_w[side], [720.0, -800.0][side], 2e-8);
        near(
            result.material_bath_conductances_w_k[side],
            [0.6, 1.0][side],
            2e-12,
        );
    }
    near(result.bath_pair_conductance_w_k, 1.5, 2e-12);
    let original = leg
        .evaluate_material_to_two_baths(
            &curve,
            500.0,
            [300.0, 700.0],
            [300.0, 700.0],
            layer_controls(),
            film,
            || Ok(()),
        )
        .unwrap();
    let zero = leg
        .evaluate_material_to_two_baths_with_resistance(
            &curve,
            500.0,
            [300.0, 700.0],
            [300.0, 700.0],
            0.0,
            layer_controls(),
            film,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(zero, original);
}

#[test]
fn g1_native_coated_source_boundary_keeps_original_storage_bath_entropy_and_late_retry() {
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
        let source = Surface {
            caloric: layers[side].caloric_law(),
            bounds: layers[side].energy_bounds_j_kg(),
            anchor: energy(side, 570.0),
            values: [0.55, 35.0],
            slopes: [0.0000005, 0.00002],
            queries: Cell::new(0),
            last_energy: Cell::new(f64::NAN),
        };
        let boundary = IntegratedCaloricGrayBoundary::new(
            &layers[side],
            face,
            0.0,
            0.0,
            radiation,
            layer_controls(),
        )
        .unwrap()
        .with_surface_law(&source)
        .unwrap();
        let thermal = single_cell(layers[side].caloric_law(), source.anchor, gas);
        let solver = chain(&thermal, materials[side], &boundary, face);
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        let (accepted, _) = solver
            .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        let initial = accepted.cells()[&CELL].specific_h_j_kg;
        assert!(initial < source.anchor);
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
            let (_, _, powers) = reference_coated(
                side,
                t,
                gas,
                radiation,
                source.coefficients(energy(side, t)),
                face,
            );
            if MASS * (energy(side, t) - initial) + dt * (powers[0] + powers[1]) > 0.0 {
                hi = t;
            } else {
                lo = t;
            }
        }
        let t = 0.5 * (lo + hi);
        let (_, _, powers) = reference_coated(
            side,
            t,
            gas,
            radiation,
            source.coefficients(energy(side, t)),
            face,
        );
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
                ThermalStorageBasis::Enthalpy
            ][side]
        );
        assert_eq!(
            next.cells()[&CELL].boundary,
            thermal.cells()[&CELL].boundary
        );
        assert_eq!(report.energy.substeps, 1);
        let mut cancelled = 0;
        assert!(matches!(
            solver.advance_to(&accepted, 2.0 * dt, controls(dt), &mut workspace, || {
                cancelled += 1;
                if cancelled == calls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },),
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
fn g4_coating_admits_external_surface_outside_material_chart_without_extrapolating_material() {
    let fixture = Fixture::new(true);
    let cv = fixture.cv();
    let k = IsochoricConductionLaw::try_new(&fixture.card, &cv, &fixture.conductivity, SELECTION)
        .unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&k).unwrap();
    for (bath, resistance) in [(1000.0, 0.01), (100.0, 0.1)] {
        let face = ImplicitBoundaryGeometry {
            surface_resistance_m2_k_w: resistance,
            ..geometry(0.02)
        };
        let boundary =
            IntegratedCaloricGrayBoundary::new(&layer, face, 0.8, 20.0, bath, layer_controls())
                .unwrap();
        let result = boundary
            .evaluate(energy(0, 500.0), bath, &mut || Ok(()))
            .unwrap();
        let (inner, surface, powers) = reference_coated(0, 500.0, bath, bath, [0.8, 20.0], face);
        assert!((300.0..=700.0).contains(&inner));
        assert!(!(300.0..=700.0).contains(&surface));
        near(result.surface_temperature_k, surface, 2e-9);
        near(result.outward_power_w, powers[0] + powers[1], 2e-8);
        near(
            surface + resistance / AREA * result.outward_power_w,
            inner,
            2e-8,
        );
        let unsupported = if bath > 700.0 { 695.0 } else { 300.0 };
        assert!(
            boundary
                .evaluate(energy(0, unsupported), bath, &mut || Ok(()))
                .is_err()
        );
        assert!(
            boundary
                .state_at_specific_enthalpy(energy(0, 500.0), bath, geometry(0.02), &mut || Ok(()),)
                .is_err()
        );
    }
    let face = ImplicitBoundaryGeometry {
        surface_resistance_m2_k_w: 0.01,
        ..geometry(0.02)
    };
    let boundary =
        IntegratedCaloricGrayBoundary::new(&layer, face, 0.8, 20.0, 1000.0, layer_controls())
            .unwrap();
    let thermal = single_cell(&cv, energy(0, 680.0), 1000.0);
    let solver = chain(&thermal, &k, &boundary, face);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (accepted, _) = solver
        .advance_to(&thermal, 0.1, controls(0.1), &mut workspace, || Ok(()))
        .unwrap();
    let initial = accepted.cells()[&CELL].specific_h_j_kg;
    assert!(
        solver
            .advance_to(&accepted, 100.1, controls(100.0), &mut workspace, || Ok(()),)
            .is_err()
    );
    assert_eq!(accepted.time_s(), 0.1);
    assert_eq!(accepted.cells()[&CELL].specific_h_j_kg, initial);
    let (retry, _) = solver
        .advance_to(&accepted, 0.2, controls(0.1), &mut workspace, || Ok(()))
        .unwrap();
    assert!(retry.temperature_k(CELL).unwrap() > accepted.temperature_k(CELL).unwrap());
}
