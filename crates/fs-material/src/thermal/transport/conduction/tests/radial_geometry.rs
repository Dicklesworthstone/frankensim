//! G1/G3/G4: explicit radial metrics, original caloric faces and native heat.
use super::*;
use fs_thermal::enthalpy::Reservoir;
use fs_thermal::enthalpy::conduction::ConductiveReservoir;
use fs_thermal::enthalpy::conduction::implicit::boundary::{
    ImplicitBoundaryGeometry, ImplicitBoundaryLaw,
};

const CELL: CellId = CellId(0);
const MATERIAL: MaterialId = MaterialId(11);

fn radial_geometry(
    cylinder: bool,
    radii_and_length: [f64; 3],
    fraction: f64,
) -> RadialCaloricBoundaryGeometry {
    let [surface, trace, length] = radii_and_length;
    if cylinder {
        RadialCaloricBoundaryGeometry::cylindrical(surface, trace, length, fraction)
    } else {
        RadialCaloricBoundaryGeometry::spherical(surface, trace, fraction)
    }
    .unwrap()
}

// Independent physical radial integrals, never the production A/d lowering.
fn radial_metric(cylinder: bool, r: [f64; 3], fraction: f64) -> (f64, f64) {
    let [surface, trace, length] = r;
    if cylinder {
        let perimeter_length = core::f64::consts::TAU * length * fraction;
        (
            perimeter_length / (surface / trace).ln().abs(),
            perimeter_length * surface,
        )
    } else {
        let angle = 4.0 * core::f64::consts::PI * fraction;
        (
            angle / (1.0 / surface - 1.0 / trace).abs(),
            angle * surface * surface,
        )
    }
}

// Exact quadratic face root for k(T)=coefficient*T and one finite film.
fn film_response(c: f64, area: f64, t: f64, bath: f64, k: f64, film: f64) -> (f64, f64) {
    let resistance = film / area;
    let p = 0.5 * c * k * resistance;
    let rhs = bath + p * t * t;
    let face = 2.0 * rhs / (1.0 + (1.0 + 4.0 * p * rhs).sqrt());
    (face, (face - bath) / resistance)
}

fn controls(dt: f64) -> ImplicitConductionControls {
    ImplicitConductionControls {
        maximum_dt_s: dt,
        maximum_halvings: 0,
        cell_tolerance_j: 1e-6,
        total_tolerance_j: 1e-6,
        ..Default::default()
    }
}

#[test]
fn radial_caloric_faces_match_thick_shell_integrals_and_original_energy_tangents() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let film = 0.02;
    for cylinder in [false, true] {
        for r in [[0.04, 0.1, 0.13], [0.1, 0.04, 0.13]] {
            for fraction in [1.0, 0.25] {
                let geometry = radial_geometry(cylinder, r, fraction);
                let (c, area) = radial_metric(cylinder, r, fraction);
                near(geometry.conductance_multiplier_m(), c, 2e-15);
                near(geometry.surface_area_m2(), area, 2e-16);
                assert_eq!(geometry.surface_radius_m(), r[0]);
                assert_eq!(geometry.trace_radius_m(), r[1]);
                assert_eq!(geometry.cylinder_length_m(), cylinder.then_some(r[2]));
                assert_eq!(geometry.area_fraction(), fraction);
                // A thick radial path is not a slab with thickness |r0-r1|.
                let planar_c = area / (r[0] - r[1]).abs();
                assert!((planar_c / c - 1.0).abs() > 0.1);
                for (side, layer) in layers.iter().enumerate() {
                    let boundary = IntegratedRadialCaloricBoundary::new(
                        layer,
                        geometry,
                        film,
                        layer_controls(),
                    )
                    .unwrap();
                    assert!(core::ptr::eq(boundary.caloric_law(), layer.caloric_law()));
                    assert_eq!(boundary.geometry(), geometry);
                    assert!(!boundary.supports_current_geometry());
                    let native = boundary.native_geometry();
                    near(native.half_distance_m, area / c, 3e-17);
                    assert_eq!(native.surface_resistance_m2_k_w, film);
                    let k = [0.04, 0.02][side];
                    for [t, bath] in [
                        [620.0, 340.0],
                        [340.0, 620.0],
                        [500.0, 1000.0],
                        [500.0, 500.0],
                    ] {
                        let (face, expected) = film_response(c, area, t, bath, k, film);
                        let response = boundary
                            .response_at_specific_energy(energy(side, t), bath, &mut || Ok(()))
                            .unwrap();
                        near(response.material_face_temperature_k, face, 2e-10);
                        near(response.power_to_reservoir_w, expected, 3e-8);
                        let state = boundary
                            .state_at_specific_enthalpy(energy(side, t), bath, native, &mut || {
                                Ok(())
                            })
                            .unwrap();
                        let capacity = if side == 0 {
                            200.0 + 0.5 * t
                        } else {
                            400.0 + t
                        };
                        let derivative =
                            c * k * t / ((1.0 + c * k * face * film / area) * capacity);
                        near(state.power_enthalpy_derivative_w_kg_j, derivative, 2e-11);
                        near(state.conductance_w_k * (t - bath), expected, 3e-8);
                        let lo = film_response(
                            c,
                            area,
                            temperature(side, energy(side, t) - 0.25),
                            bath,
                            k,
                            film,
                        )
                        .1;
                        let hi = film_response(
                            c,
                            area,
                            temperature(side, energy(side, t) + 0.25),
                            bath,
                            k,
                            film,
                        )
                        .1;
                        near(
                            state.power_enthalpy_derivative_w_kg_j,
                            (hi - lo) / 0.5,
                            3e-8,
                        );
                    }
                    let ideal = IntegratedRadialCaloricBoundary::new(
                        layer,
                        geometry,
                        0.0,
                        layer_controls(),
                    )
                    .unwrap();
                    let q = ideal
                        .response_at_specific_energy(energy(side, 620.0), 340.0, &mut || Ok(()))
                        .unwrap()
                        .power_to_reservoir_w;
                    near(
                        q,
                        0.5 * c * k * (620.0_f64.powi(2) - 340.0_f64.powi(2)),
                        3e-8,
                    );
                    assert!(
                        (q / (0.5 * planar_c * k * (620.0_f64.powi(2) - 340.0_f64.powi(2))) - 1.0)
                            .abs()
                            > 0.1
                    );
                    assert!(
                        ideal
                            .response_at_specific_energy(
                                energy(side, 500.0),
                                1000.0,
                                &mut || Ok(()),
                            )
                            .is_err()
                    );
                    assert!(
                        boundary
                            .state_at_specific_enthalpy(
                                energy(side, 500.0),
                                350.0,
                                ImplicitBoundaryGeometry {
                                    half_distance_m: (r[0] - r[1]).abs(),
                                    ..native
                                },
                                &mut || Ok(()),
                            )
                            .is_err()
                    );
                    assert!(
                        boundary
                            .state_at_specific_enthalpy_with_geometry_tangent(
                                energy(side, 500.0),
                                350.0,
                                native,
                                0.0,
                                &mut || Ok(()),
                            )
                            .unwrap()
                            .is_none()
                    );
                    assert!(matches!(
                        boundary.response_at_specific_energy(
                            energy(side, 500.0),
                            350.0,
                            &mut || Err(DepositionError::Cancelled),
                        ),
                        Err(DepositionError::Cancelled)
                    ));
                }
            }
        }
    }
    for bad in [
        RadialCaloricBoundaryGeometry::cylindrical(0.0, 0.1, 0.13, 1.0),
        RadialCaloricBoundaryGeometry::cylindrical(0.04, 0.1, 0.0, 1.0),
        RadialCaloricBoundaryGeometry::spherical(0.04, 0.04, 1.0),
        RadialCaloricBoundaryGeometry::spherical(0.04, 0.1, 0.0),
        RadialCaloricBoundaryGeometry::spherical(0.04, 0.1, 1.01),
        RadialCaloricBoundaryGeometry::spherical(f64::NAN, 0.1, 1.0),
    ] {
        assert!(bad.is_err());
    }
}

// Actual reference shell volume and a volume-average radial temperature trace.
// No resistance coordinate is used to construct physical inventory.
fn shell(cylinder: bool, inner: f64, outer: f64, length: f64) -> (f64, f64, f64) {
    if cylinder {
        (
            core::f64::consts::PI * length * (outer * outer - inner * inner),
            core::f64::consts::TAU * (length * (inner + outer) + outer * outer - inner * inner),
            (2.0 / 3.0) * (outer.powi(3) - inner.powi(3)) / (outer * outer - inner * inner),
        )
    } else {
        (
            4.0 * core::f64::consts::PI / 3.0 * (outer.powi(3) - inner.powi(3)),
            4.0 * core::f64::consts::PI * (inner * inner + outer * outer),
            0.75 * (outer.powi(4) - inner.powi(4)) / (outer.powi(3) - inner.powi(3)),
        )
    }
}

fn shell_cell<'a>(
    law: &'a dyn ThermalLaw,
    initial_energy: f64,
    cylinder: bool,
    radii_and_length: [f64; 3],
    bath: f64,
) -> DepositionSimulation<'a> {
    let [inner, outer, length] = radii_and_length;
    let (volume, area, trace) = shell(cylinder, inner, outer, length);
    let (_, footprint) = radial_metric(cylinder, [outer, trace, length], 1.0);
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELL,
                segment_id: 1,
                body_id: 1,
                material_id: MATERIAL,
                mass_kg: RHO * volume,
                specific_h_j_kg: initial_energy,
                geometry: CellGeometry {
                    start_m: [inner, 0.0, 0.0],
                    end_m: [outer, 0.0, 0.0],
                    volume_m3: volume,
                    surface_area_m2: area,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: Some(Reservoir {
                        area_m2: footprint,
                        temperature_k: bath,
                        conductance_w_k: 999.0,
                    }),
                },
            }),
        }],
        &[(MATERIAL, law)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 8,
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

#[test]
fn native_radial_boundary_keeps_physical_shell_mass_and_exact_bath_heat() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let conductivity: [&dyn ConductivityLaw; 2] = [&a, &b];
    let implicit: [&dyn ImplicitConductionLaw; 2] = [&a, &b];
    let (inner, outer, length, initial, bath, dt) = (0.03, 0.08, 0.1, 320.0, 650.0, 20.0);
    for cylinder in [false, true] {
        let (volume, _, trace) = shell(cylinder, inner, outer, length);
        let (c, _) = radial_metric(cylinder, [outer, trace, length], 1.0);
        let mass = RHO * volume;
        for side in 0..2 {
            let thermal = shell_cell(
                layers[side].caloric_law(),
                energy(side, initial),
                cylinder,
                [inner, outer, length],
                bath,
            );
            let boundary = IntegratedRadialCaloricBoundary::new(
                &layers[side],
                radial_geometry(cylinder, [outer, trace, length], 1.0),
                0.0,
                layer_controls(),
            )
            .unwrap();
            let native = boundary.native_geometry();
            let observer = ConductionObserver::new(
                &thermal,
                &[(MATERIAL, conductivity[side])],
                &[],
                &[ConductiveReservoir {
                    cell_id: CELL,
                    half_distance_m: native.half_distance_m,
                    surface_resistance_m2_k_w: 0.0,
                }],
                || Ok(()),
            )
            .unwrap();
            let chain = ImplicitConductionChain::new(
                &thermal,
                &observer,
                &[(MATERIAL, implicit[side])],
                || Ok(()),
            )
            .unwrap()
            .with_boundary_laws(&thermal, &[(CELL, &boundary)], &[], || Ok(()))
            .unwrap();
            let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
            let (next, report) = chain
                .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
                .unwrap();
            let (linear, quadratic, k) = if side == 0 {
                (200.0, 0.25, 0.04)
            } else {
                (400.0, 0.5, 0.02)
            };
            let g = 0.5 * c * k;
            let coefficient = mass * quadratic + dt * g;
            let rhs =
                mass * (linear * initial + quadratic * initial * initial) + dt * g * bath * bath;
            let expected = 2.0 * rhs
                / (mass * linear + ((mass * linear).powi(2) + 4.0 * coefficient * rhs).sqrt());
            let heat = dt * g * (bath * bath - expected * expected);
            near(next.temperature_k(CELL).unwrap(), expected, 2e-8);
            near(report.energy.external_heat_j, heat, 2e-6);
            near(
                report.energy.stored_energy_change_j,
                mass * (energy(side, expected) - energy(side, initial)),
                2e-6,
            );
            near(
                report.entropy.reservoir_entropy_input_j_k,
                heat / bath,
                2e-8,
            );
            assert_eq!(next.cells()[&CELL].mass_kg, mass);
            assert_eq!(
                next.cells()[&CELL].geometry,
                thermal.cells()[&CELL].geometry
            );
            assert_eq!(
                next.cells()[&CELL].boundary,
                thermal.cells()[&CELL].boundary
            );
            assert_eq!(report.energy.substeps, 1);
            assert_eq!(thermal.time_s(), 0.0);
        }
    }
}
