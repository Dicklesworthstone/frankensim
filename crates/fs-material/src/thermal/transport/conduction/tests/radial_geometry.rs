//! G1/G3/G4: explicit radial metrics, original caloric faces and native heat.
use super::*;
use fs_thermal::enthalpy::Reservoir;
use fs_thermal::enthalpy::conduction::implicit::boundary::{
    ImplicitBoundaryGeometry, ImplicitBoundaryLaw,
};
use fs_thermal::enthalpy::conduction::implicit::contact::{
    ImplicitContactLaw, ImplicitContactState,
};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionGraph, ImplicitConductionLinearControls,
};
use fs_thermal::enthalpy::conduction::{ConductiveReservoir, CurrentConductionGeometry};

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

fn radial_link_geometry(cylinder: bool, traces: [f64; 2], face: f64, length: f64) -> LayerGeometry {
    if cylinder {
        LayerGeometry::cylindrical(traces, face, length)
    } else {
        LayerGeometry::spherical(traces, face)
    }
    .unwrap()
}

// Exact shell volumes set inventory independently of either heat metric.
fn radial_cells<'a>(
    laws: [&'a dyn ThermalLaw; 2],
    cylinder: bool,
    edges: [f64; 3],
    length: f64,
) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for side in 0..2 {
        let (volume, surface, _) = shell(cylinder, edges[side], edges[side + 1], length);
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(side as u64),
                segment_id: side as u64,
                body_id: 1,
                material_id: [CV_MATERIAL, CP_MATERIAL][side],
                mass_kg: RHO * volume,
                specific_h_j_kg: energy(side, [620.0, 360.0][side]),
                geometry: CellGeometry {
                    start_m: [edges[side], 0.0, 0.0],
                    end_m: [edges[side + 1], 0.0, 0.0],
                    volume_m3: volume,
                    surface_area_m2: surface,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: None,
                },
            }),
        });
    }
    let (_, area) = radial_metric(cylinder, [edges[1], edges[0], length], 1.0);
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: LINK,
            a: CellId(0),
            b: CellId(1),
            area_m2: area,
            conductance_w_k: 0.0,
        }),
    });
    DepositionSimulation::new(
        events,
        &[(CV_MATERIAL, laws[0]), (CP_MATERIAL, laws[1])],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
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

struct RadialContact<'a> {
    laws: [&'a dyn ThermalLaw; 2],
    bounds: [[f64; 2]; 2],
}
fn radial_resistance(e: [f64; 2]) -> f64 {
    0.002 * (1.0 + e[0] / 200_000.0 + e[1] / 400_000.0)
}
impl ImplicitContactLaw for RadialContact<'_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        self.laws
    }
    fn enthalpy_bounds_j_kg(&self) -> [[f64; 2]; 2] {
        self.bounds
    }
    fn state_at_specific_enthalpies(
        &self,
        e: [f64; 2],
    ) -> Result<ImplicitContactState, DepositionError> {
        Ok(ImplicitContactState {
            resistance_m2_k_w: radial_resistance(e),
            resistance_enthalpy_derivatives_m2_kg_k_w_j: [1e-8, 5e-9],
        })
    }
}

// Independent total-power root using physical radial integrals and one
// contact jump. Production equivalent distances and face roots are not used.
fn radial_link_power(e: [f64; 2], integrals: [f64; 2], area: f64) -> f64 {
    let t = [temperature(0, e[0]), temperature(1, e[1])];
    let sign = if t[0] >= t[1] { 1.0 } else { -1.0 };
    let (mut lo, mut hi) = (
        0.0,
        (t[0] * t[0] - t[1] * t[1]).abs() / (2.0 * (integrals[0] / 0.04 + integrals[1] / 0.02)),
    );
    for _ in 0..80 {
        let q = 0.5 * (lo + hi);
        let face_a = (t[0] * t[0] - 2.0 * sign * q * integrals[0] / 0.04).sqrt();
        let face_b = (t[1] * t[1] + 2.0 * sign * q * integrals[1] / 0.02).sqrt();
        if sign * (face_a - face_b) > q * radial_resistance(e) / area {
            lo = q;
        } else {
            hi = q;
        }
    }
    sign * 0.5 * (lo + hi)
}

#[test]
fn current_radial_contact_retains_source_tangents_and_native_shell_inventory() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let contact = RadialContact {
        laws: [&cv, &cp],
        bounds: [
            cv.specific_internal_energy_span_j_kg(),
            cp.specific_enthalpy_span_j_kg(),
        ],
    };
    let edges: [f64; 3] = [0.02, 0.06, 0.1];
    let (reference_length, current_length, new_inner): (f64, f64, f64) = (0.1, 0.13, 0.035);
    let initial_e = [energy(0, 620.0), energy(1, 360.0)];
    for cylinder in [false, true] {
        let current_edges = edges.map(|r| {
            if cylinder {
                (new_inner * new_inner
                    + (r * r - edges[0] * edges[0]) * reference_length / current_length)
                    .sqrt()
            } else {
                (new_inner.powi(3) + r.powi(3) - edges[0].powi(3)).cbrt()
            }
        });
        let original_shells =
            [0, 1].map(|i| shell(cylinder, edges[i], edges[i + 1], reference_length));
        let current_shells = [0, 1].map(|i| {
            shell(
                cylinder,
                current_edges[i],
                current_edges[i + 1],
                current_length,
            )
        });
        let masses = original_shells.map(|s| RHO * s.0);
        for side in 0..2 {
            near(current_shells[side].0, original_shells[side].0, 2e-18);
        }
        let original_geometry = radial_link_geometry(
            cylinder,
            original_shells.map(|s| s.2),
            edges[1],
            reference_length,
        );
        let geometry = radial_link_geometry(
            cylinder,
            current_shells.map(|s| s.2),
            current_edges[1],
            current_length,
        );
        let thermal = radial_cells([&cv, &cp], cylinder, edges, reference_length);
        let original = IntegratedCaloricLink::new_with_contact(
            &thermal,
            LINK,
            [&layers[0], &layers[1]],
            original_geometry,
            &contact,
            layer_controls(),
            || Ok(()),
        )
        .unwrap();
        let current = original
            .with_current_radial_geometry(&thermal, geometry, || Ok(()))
            .unwrap();
        assert_eq!(current.reference_geometry(), original_geometry);
        assert_eq!(current.geometry(), geometry);
        assert_eq!(
            current.geometry().cylinder_length_m(),
            cylinder.then_some(current_length)
        );
        assert!(current.supports_current_geometry());
        assert_eq!(current.native_geometry().contact_resistance_m2_k_w, 0.0);
        assert!(core::ptr::addr_eq(
            current.caloric_laws()[0],
            &cv as &dyn ThermalLaw
        ));
        assert!(core::ptr::addr_eq(
            current.caloric_laws()[1],
            &cp as &dyn ThermalLaw
        ));
        let independent = current_shells
            .map(|s| radial_metric(cylinder, [current_edges[1], s.2, current_length], 1.0));
        let integrals = independent.map(|(c, _)| 1.0 / c);
        let area = independent[0].1;
        let mut patch = thermal.contact_patch(LINK).unwrap();
        patch.link.area_m2 = current.current_geometry().area_m2;
        near(patch.link.area_m2, area, 1e-16);
        for t in [[620.0, 360.0], [350.0, 610.0], [500.0, 500.0]] {
            let e = [energy(0, t[0]), energy(1, t[1])];
            let state = current
                .state_at_specific_enthalpies(e, current.native_geometry(), patch, &mut || Ok(()))
                .unwrap();
            near(
                state.conductance_w_k * (t[0] - t[1]),
                radial_link_power(e, integrals, area),
                3e-8,
            );
            for side in 0..2 {
                let mut lo = e;
                let mut hi = e;
                lo[side] -= 0.25;
                hi[side] += 0.25;
                near(
                    state.power_enthalpy_derivatives_w_kg_j[side],
                    (radial_link_power(hi, integrals, area)
                        - radial_link_power(lo, integrals, area))
                        / 0.5,
                    3e-8,
                );
            }
        }
        assert!(
            current
                .state_at_specific_enthalpies(
                    initial_e,
                    original.native_geometry(),
                    patch,
                    &mut || Ok(()),
                )
                .is_err()
        );
        let mut changed_coverage = patch;
        changed_coverage.covered_area_m2[0] *= 0.9;
        assert!(
            current
                .state_at_specific_enthalpies(
                    initial_e,
                    current.native_geometry(),
                    changed_coverage,
                    &mut || Ok(()),
                )
                .is_err()
        );
        assert!(
            current
                .state_at_specific_enthalpies_with_geometry_tangent(
                    initial_e,
                    current.native_geometry(),
                    patch,
                    [0.0; 2],
                    &mut || Ok(()),
                )
                .unwrap()
                .is_none()
        );
        let (traces, face) = geometry.radial_coordinates_m().unwrap();
        for invalid in [
            LayerGeometry::planar(geometry.equivalent_distances_m(), area).unwrap(),
            radial_link_geometry(!cylinder, traces, face, current_length),
            radial_link_geometry(cylinder, [traces[1], traces[0]], face, current_length),
        ] {
            assert!(
                original
                    .with_current_radial_geometry(&thermal, invalid, || Ok(()))
                    .is_err()
            );
        }
        assert!(matches!(
            original.with_current_radial_geometry(&thermal, geometry, || Err(
                DepositionError::Cancelled
            ),),
            Err(DepositionError::Cancelled)
        ));

        let observer = ConductionObserver::new(
            &thermal,
            &[
                (CV_MATERIAL, &a as &dyn ConductivityLaw),
                (CP_MATERIAL, &b as &dyn ConductivityLaw),
            ],
            &[original.native_geometry()],
            &[],
            || Ok(()),
        )
        .unwrap()
        .with_current_geometry(
            &thermal,
            CurrentConductionGeometry {
                links: [(LINK, current.current_geometry())].into(),
                ..Default::default()
            },
            || Ok(()),
        )
        .unwrap();
        let graph = ImplicitConductionGraph::new(
            &thermal,
            &observer,
            &[
                (CV_MATERIAL, &a as &dyn ImplicitConductionLaw),
                (CP_MATERIAL, &b as &dyn ImplicitConductionLaw),
            ],
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .unwrap()
        .with_link_laws(&thermal, &[(LINK, &current)], || Ok(()))
        .unwrap();
        let dt = 20.0;
        let energies = |q| [initial_e[0] - q / masses[0], initial_e[1] + q / masses[1]];
        let residual = |q| q - dt * radial_link_power(energies(q), integrals, area);
        let (mut lo, mut hi) = (0.0, masses[0] * (initial_e[0] - energy(0, 360.0)));
        assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
        for _ in 0..80 {
            let q = 0.5 * (lo + hi);
            if residual(q) < 0.0 {
                lo = q;
            } else {
                hi = q;
            }
        }
        let expected_e = energies(0.5 * (lo + hi));
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        let (next, report) = graph
            .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        for side in 0..2 {
            let id = CellId(side as u64);
            near(next.cells()[&id].specific_h_j_kg, expected_e[side], 1e-5);
            assert_eq!(next.cells()[&id].mass_kg, masses[side]);
            assert_eq!(next.cells()[&id].geometry, thermal.cells()[&id].geometry);
        }
        near(report.energy.stored_energy_change_j, 0.0, 1e-6);
        near(report.energy.external_heat_j, 0.0, 0.0);
        assert!(report.entropy.conduction_generation_j_k > 0.0);
        assert_eq!(next.links(), thermal.links());
        assert_eq!(next.contact_patch(LINK), thermal.contact_patch(LINK));
        assert_eq!(report.energy.substeps, 1);
        assert_eq!(thermal.time_s(), 0.0);
    }
}

#[test]
fn current_radial_bath_uses_current_film_area_and_original_shell_mass() {
    let fixture = Fixture::new(true);
    let caloric = fixture.cv();
    let conductivity =
        IsochoricConductionLaw::try_new(&fixture.card, &caloric, &fixture.conductivity, SELECTION)
            .unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conductivity).unwrap();
    let (inner, outer, length, new_inner, new_length): (f64, f64, f64, f64, f64) =
        (0.03, 0.08, 0.1, 0.05, 0.13);
    let (initial, bath, film, dt) = (350.0, 650.0, 0.02, 20.0);
    for cylinder in [false, true] {
        let new_outer = if cylinder {
            (new_inner * new_inner + (outer * outer - inner * inner) * length / new_length).sqrt()
        } else {
            (new_inner.powi(3) + outer.powi(3) - inner.powi(3)).cbrt()
        };
        let (volume, _, trace) = shell(cylinder, inner, outer, length);
        let (new_volume, _, new_trace) = shell(cylinder, new_inner, new_outer, new_length);
        near(new_volume, volume, 2e-18);
        let mass = RHO * volume;
        let reference_geometry = radial_geometry(cylinder, [outer, trace, length], 1.0);
        let geometry = radial_geometry(cylinder, [new_outer, new_trace, new_length], 1.0);
        let original = IntegratedRadialCaloricBoundary::new(
            &layer,
            reference_geometry,
            film,
            layer_controls(),
        )
        .unwrap();
        let current = original.with_current_radial_geometry(geometry).unwrap();
        assert_eq!(current.reference_geometry(), reference_geometry);
        assert_eq!(current.geometry(), geometry);
        assert!(current.supports_current_geometry());
        assert!(core::ptr::addr_eq(
            current.caloric_law(),
            &caloric as &dyn ThermalLaw
        ));
        let (c, area) = radial_metric(cylinder, [new_outer, new_trace, new_length], 1.0);
        let (_, reference_area) = radial_metric(cylinder, [outer, trace, length], 1.0);
        let expected = film_response(c, area, initial, bath, 0.04, film);
        // Reusing R''/A_reference after motion would produce a different heat.
        let stale = film_response(c, reference_area, initial, bath, 0.04, film);
        assert!((expected.1 / stale.1 - 1.0).abs() > 0.05);
        let response = current
            .response_at_specific_energy(energy(0, initial), bath, &mut || Ok(()))
            .unwrap();
        near(response.material_face_temperature_k, expected.0, 2e-10);
        near(response.power_to_reservoir_w, expected.1, 3e-8);
        let state = current
            .state_at_specific_enthalpy(
                energy(0, initial),
                bath,
                current.native_geometry(),
                &mut || Ok(()),
            )
            .unwrap();
        near(
            state.power_enthalpy_derivative_w_kg_j,
            c * 0.04 * initial
                / ((1.0 + c * 0.04 * expected.0 * film / area) * (200.0 + 0.5 * initial)),
            2e-11,
        );
        assert!(
            current
                .state_at_specific_enthalpy(
                    energy(0, initial),
                    bath,
                    original.native_geometry(),
                    &mut || Ok(())
                )
                .is_err()
        );
        assert!(
            current
                .state_at_specific_enthalpy_with_geometry_tangent(
                    energy(0, initial),
                    bath,
                    current.native_geometry(),
                    0.0,
                    &mut || Ok(()),
                )
                .unwrap()
                .is_none()
        );
        for invalid in [
            radial_geometry(!cylinder, [new_outer, new_trace, new_length], 1.0),
            radial_geometry(cylinder, [new_trace, new_outer, new_length], 1.0),
            radial_geometry(cylinder, [new_outer, new_trace, new_length], 0.5),
        ] {
            assert!(original.with_current_radial_geometry(invalid).is_err());
        }
        let thermal = shell_cell(
            &caloric,
            energy(0, initial),
            cylinder,
            [inner, outer, length],
            bath,
        );
        let observer = ConductionObserver::new(
            &thermal,
            &[(MATERIAL, &conductivity as &dyn ConductivityLaw)],
            &[],
            &[ConductiveReservoir {
                cell_id: CELL,
                half_distance_m: original.native_geometry().half_distance_m,
                surface_resistance_m2_k_w: film,
            }],
            || Ok(()),
        )
        .unwrap()
        .with_current_geometry(
            &thermal,
            CurrentConductionGeometry {
                reservoirs: [(CELL, current.current_geometry())].into(),
                ..Default::default()
            },
            || Ok(()),
        )
        .unwrap();
        let chain =
            ImplicitConductionChain::new(&thermal, &observer, &[(MATERIAL, &conductivity)], || {
                Ok(())
            })
            .unwrap()
            .with_boundary_laws(&thermal, &[(CELL, &current)], &[], || Ok(()))
            .unwrap();
        let (mut lo, mut hi) = (initial, bath);
        for _ in 0..80 {
            let t = 0.5 * (lo + hi);
            if mass * (energy(0, t) - energy(0, initial))
                + dt * film_response(c, area, t, bath, 0.04, film).1
                < 0.0
            {
                lo = t;
            } else {
                hi = t;
            }
        }
        let expected_t = 0.5 * (lo + hi);
        let heat = -dt * film_response(c, area, expected_t, bath, 0.04, film).1;
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        let (next, report) = chain
            .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        near(next.temperature_k(CELL).unwrap(), expected_t, 2e-8);
        near(report.energy.external_heat_j, heat, 2e-6);
        near(report.energy.stored_energy_change_j, heat, 2e-6);
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
