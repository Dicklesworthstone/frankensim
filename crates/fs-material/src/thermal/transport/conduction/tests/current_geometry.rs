//! Current planar source geometry with original Cp/Cv storage and exact fluxes.
use super::*;
use fs_thermal::enthalpy::conduction::implicit::boundary::{
    ImplicitBoundaryGeometry, ImplicitBoundaryLaw,
};
use fs_thermal::enthalpy::conduction::implicit::contact::{
    ImplicitContactLaw, ImplicitContactState,
};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionGraph, ImplicitConductionLinearControls,
};
use fs_thermal::enthalpy::conduction::{
    ConductiveAmbient, ConductiveReservoir, CurrentConductionGeometry, CurrentConductiveBoundary,
    CurrentConductiveLink,
};
use fs_thermal::enthalpy::{ContactPatchUpdate, Reservoir};
use std::cell::Cell;

fn controls(dt: f64) -> ImplicitConductionControls {
    ImplicitConductionControls {
        maximum_dt_s: dt,
        maximum_halvings: 0,
        cell_tolerance_j: 1e-6,
        total_tolerance_j: 1e-6,
        ..Default::default()
    }
}

fn resistance(e: [f64; 2]) -> f64 {
    0.002 * (1.0 + e[0] / 200_000.0 + e[1] / 400_000.0)
}

struct Contact<'a> {
    laws: [&'a dyn ThermalLaw; 2],
    bounds: [[f64; 2]; 2],
    calls: Cell<usize>,
}
impl ImplicitContactLaw for Contact<'_> {
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
        self.calls.set(self.calls.get() + 1);
        Ok(ImplicitContactState {
            resistance_m2_k_w: resistance(e),
            resistance_enthalpy_derivatives_m2_kg_k_w_j: [1e-8, 5e-9],
        })
    }
}

// Independent quadratic Kirchhoff primitives, with a scalar contact jump.
// No production layer root, caloric inverse or derivative enters this oracle.
fn link_power(e: [f64; 2], metric: CurrentConductiveLink, r: f64) -> f64 {
    let t = [temperature(0, e[0]), temperature(1, e[1])];
    let sign = if t[0] >= t[1] { 1.0 } else { -1.0 };
    let d = metric.half_distances_m;
    let (mut lo, mut hi) = (
        0.0,
        (t[0] * t[0] - t[1] * t[1]).abs() / (2.0 * (d[0] / 0.04 + d[1] / 0.02)),
    );
    for _ in 0..80 {
        let j = 0.5 * (lo + hi);
        let a = (t[0] * t[0] - 2.0 * sign * j * d[0] / 0.04).sqrt();
        let b = (t[1] * t[1] + 2.0 * sign * j * d[1] / 0.02).sqrt();
        if sign * (a - b) > r * j {
            lo = j;
        } else {
            hi = j;
        }
    }
    sign * metric.area_m2 * 0.5 * (lo + hi)
}

fn current_patch(thermal: &DepositionSimulation<'_>, area: f64) -> ContactPatch {
    let mut patch = thermal.contact_patch(LINK).unwrap();
    patch.link.area_m2 = area;
    patch
}

#[test]
fn current_link_reuses_original_integrals_contact_and_reference_identity() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let thermal = timeline([&cv, &cp], [620.0, 360.0], [0.4, 0.32], AREA);
    let geometry = LayerGeometry::planar([0.01, 0.008], AREA).unwrap();
    let contact = Contact {
        laws: [&cv, &cp],
        bounds: [
            cv.specific_internal_energy_span_j_kg(),
            cp.specific_enthalpy_span_j_kg(),
        ],
        calls: Cell::new(0),
    };
    let original = IntegratedCaloricLink::new_with_contact(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        geometry,
        &contact,
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    assert!(!original.supports_current_geometry());
    // Both A*d products stay fixed, compatible with the original Cv density.
    for scale in [2.0, 0.5] {
        let metric = CurrentConductiveLink {
            area_m2: scale * AREA,
            half_distances_m: [0.01 / scale, 0.008 / scale],
        };
        contact.calls.set(0);
        let current = original
            .with_current_planar_geometry(&thermal, metric, || Ok(()))
            .unwrap();
        assert_eq!(contact.calls.get(), 1);
        assert!(current.supports_current_geometry());
        assert_eq!(current.reference_geometry(), geometry);
        assert_eq!(current.current_geometry(), metric);
        assert_eq!(current.native_geometry().contact_resistance_m2_k_w, 0.0);
        assert!(core::ptr::addr_eq(
            current.caloric_laws()[0],
            &cv as &dyn ThermalLaw
        ));
        for t in [[620.0, 360.0], [350.0, 610.0], [500.0, 500.0]] {
            let e = [energy(0, t[0]), energy(1, t[1])];
            contact.calls.set(0);
            let state = current
                .state_at_specific_enthalpies(
                    e,
                    current.native_geometry(),
                    current_patch(&thermal, metric.area_m2),
                    &mut || Ok(()),
                )
                .unwrap();
            assert_eq!(contact.calls.get(), 1);
            near(
                state.conductance_w_k * (t[0] - t[1]),
                link_power(e, metric, resistance(e)),
                3e-8,
            );
            for side in 0..2 {
                let mut lo = e;
                let mut hi = e;
                lo[side] -= 0.25;
                hi[side] += 0.25;
                let derivative = (link_power(hi, metric, resistance(hi))
                    - link_power(lo, metric, resistance(lo)))
                    / 0.5;
                near(
                    state.power_enthalpy_derivatives_w_kg_j[side],
                    derivative,
                    3e-8,
                );
            }
        }
        assert!(
            current
                .state_at_specific_enthalpies(
                    [energy(0, 620.0), energy(1, 360.0)],
                    original.native_geometry(),
                    thermal.contact_patch(LINK).unwrap(),
                    &mut || Ok(()),
                )
                .is_err()
        );
        assert!(
            current
                .state_at_specific_enthalpies_with_geometry_tangent(
                    [energy(0, 620.0), energy(1, 360.0)],
                    current.native_geometry(),
                    current_patch(&thermal, metric.area_m2),
                    [0.0; 2],
                    &mut || Ok(()),
                )
                .unwrap()
                .is_none()
        );
    }
    let metric = CurrentConductiveLink {
        area_m2: 2.0 * AREA,
        half_distances_m: [0.005, 0.004],
    };
    let old = thermal.contact_patch(LINK).unwrap();
    let changed = thermal
        .replace_contact_patches(
            &[ContactPatchUpdate {
                id: LINK,
                expected: Some(old),
                replacement: Some(ContactPatch {
                    covered_area_m2: [0.9 * AREA, 0.8 * AREA],
                    ..old
                }),
            }],
            || Ok(()),
        )
        .unwrap();
    assert!(
        original
            .with_current_planar_geometry(&changed, metric, || Ok(()))
            .is_err()
    );
    assert!(
        original
            .with_current_planar_geometry(&thermal, metric, || Err(DepositionError::Cancelled),)
            .is_err()
    );
    for invalid in [
        CurrentConductiveLink {
            area_m2: 0.0,
            ..metric
        },
        CurrentConductiveLink {
            half_distances_m: [f64::NAN, 0.004],
            ..metric
        },
    ] {
        assert!(
            original
                .with_current_planar_geometry(&thermal, invalid, || Ok(()))
                .is_err()
        );
    }
    // A radial equivalent resistance cannot silently become a planar thickness.
    let radial = LayerGeometry::spherical([0.08, 0.12], 0.1).unwrap();
    let radial_thermal = timeline(
        [&cv, &cp],
        [620.0, 360.0],
        [0.4, 0.32],
        radial.contact_area_m2(),
    );
    let radial_link = IntegratedCaloricLink::new(
        &radial_thermal,
        LINK,
        [&layers[0], &layers[1]],
        radial,
        0.002,
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    assert!(
        radial_link
            .with_current_planar_geometry(&radial_thermal, metric, || Ok(()))
            .is_err()
    );
    assert_eq!(thermal.contact_patch(LINK), Some(old));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(original.geometry(), geometry);
}

#[test]
fn current_integrated_link_drives_native_be_and_sdirk_without_remapping_storage() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let initial_t = [650.0, 350.0];
    let masses = [0.4, 0.32];
    let thermal = timeline([&cv, &cp], initial_t, masses, AREA);
    let original = IntegratedCaloricLink::new(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        LayerGeometry::planar([0.01, 0.008], AREA).unwrap(),
        0.0,
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    let metric = CurrentConductiveLink {
        area_m2: 2.0 * AREA,
        half_distances_m: [0.005, 0.004],
    };
    let current = original
        .with_current_planar_geometry(&thermal, metric, || Ok(()))
        .unwrap();
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
    let laws = [
        (CV_MATERIAL, &a as &dyn ImplicitConductionLaw),
        (CP_MATERIAL, &b as &dyn ImplicitConductionLaw),
    ];
    assert!(
        ImplicitConductionChain::new(&thermal, &observer, &laws, || Ok(()))
            .unwrap()
            .with_link_laws(&thermal, &[(LINK, &original)], || Ok(()))
            .is_err()
    );
    let chain = ImplicitConductionChain::new(&thermal, &observer, &laws, || Ok(()))
        .unwrap()
        .with_link_laws(&thermal, &[(LINK, &current)], || Ok(()))
        .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &observer,
        &laws,
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_link_laws(&thermal, &[(LINK, &current)], || Ok(()))
    .unwrap();
    let initial_e = [energy(0, initial_t[0]), energy(1, initial_t[1])];
    let energies = |q| [initial_e[0] - q / masses[0], initial_e[1] + q / masses[1]];
    let power = |q| {
        let e = energies(q);
        let t = [temperature(0, e[0]), temperature(1, e[1])];
        metric.area_m2 * (t[0] * t[0] - t[1] * t[1])
            / (2.0 * (metric.half_distances_m[0] / 0.04 + metric.half_distances_m[1] / 0.02))
    };
    let solve_heat = |earlier: f64, duration: f64| {
        let (mut lo, mut hi) = (0.0, masses[0] * (initial_e[0] - energy(0, initial_t[1])));
        assert!(lo - earlier - duration * power(lo) <= 0.0);
        assert!(hi - earlier - duration * power(hi) > 0.0);
        for _ in 0..80 {
            let q = 0.5 * (lo + hi);
            if q - earlier - duration * power(q) < 0.0 {
                lo = q;
            } else {
                hi = q;
            }
        }
        0.5 * (lo + hi)
    };
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    for second_order in [false, true] {
        let dt = if second_order { 2.0 } else { 50.0 };
        let expected_heat = if second_order {
            let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
            let first = solve_heat(0.0, gamma * dt);
            solve_heat((1.0 - gamma) * dt * power(first), gamma * dt)
        } else {
            solve_heat(0.0, dt)
        };
        for sparse in [false, true] {
            let (next, report) = match (sparse, second_order) {
                (false, false) => {
                    chain.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
                }
                (true, false) => {
                    graph.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
                }
                (false, true) => {
                    chain.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
                }
                (true, true) => {
                    graph.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
                }
            }
            .unwrap();
            for (side, expected) in energies(expected_heat).into_iter().enumerate() {
                let id = CellId(side as u64);
                near(next.cells()[&id].specific_h_j_kg, expected, 1e-5);
                assert_eq!(next.cells()[&id].mass_kg, masses[side]);
                assert_eq!(next.cells()[&id].geometry, thermal.cells()[&id].geometry);
            }
            near(report.energy.external_heat_j, 0.0, 0.0);
            near(report.energy.stored_energy_change_j, 0.0, 1e-6);
            assert!(report.entropy.conduction_generation_j_k > 0.0);
            assert_eq!(report.energy.substeps, 1);
            assert_eq!(next.links(), thermal.links());
            assert_eq!(
                next.material_storage_basis(CV_MATERIAL).unwrap(),
                ThermalStorageBasis::InternalEnergy
            );
            assert_eq!(
                next.material_storage_basis(CP_MATERIAL).unwrap(),
                ThermalStorageBasis::Enthalpy
            );
        }
    }
    assert_eq!(thermal.time_s(), 0.0);
}

// Exact material-face solution for k(T)=coefficient*T and a finite areal film.
fn bath_response(t: f64, bath: f64, k: f64, geometry: ImplicitBoundaryGeometry) -> (f64, f64) {
    let p = geometry.surface_resistance_m2_k_w * k / (2.0 * geometry.half_distance_m);
    let rhs = bath + p * t * t;
    let face = 2.0 * rhs / (1.0 + (1.0 + 4.0 * p * rhs).sqrt());
    (
        face,
        geometry.area_m2 / geometry.surface_resistance_m2_k_w * (face - bath),
    )
}

#[test]
fn current_boundary_reprepares_film_and_exact_caloric_tangent_without_source_widening() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let reference = ImplicitBoundaryGeometry {
        area_m2: AREA,
        half_distance_m: 0.015,
        surface_resistance_m2_k_w: 0.004,
    };
    for (side, layer) in layers.iter().enumerate() {
        let original = IntegratedCaloricBoundary::new(layer, reference, layer_controls()).unwrap();
        let current = original
            .with_current_planar_geometry(CurrentConductiveBoundary {
                area_m2: 3.0 * AREA,
                half_distance_m: 0.005,
            })
            .unwrap();
        assert!(!original.supports_current_geometry());
        assert!(current.supports_current_geometry());
        assert!(core::ptr::eq(current.caloric_law(), layer.caloric_law()));
        assert_eq!(current.reference_geometry(), reference);
        let geometry = current.native_geometry();
        let k = [0.04, 0.02][side];
        for [t, bath] in [
            [620.0, 340.0],
            [340.0, 620.0],
            [500.0, 1000.0],
            [500.0, 500.0],
        ] {
            let e = energy(side, t);
            let (face, power) = bath_response(t, bath, k, geometry);
            let response = current
                .response_at_specific_energy(e, bath, &mut || Ok(()))
                .unwrap();
            near(response.material_face_temperature_k, face, 2e-10);
            near(response.power_to_reservoir_w, power, 3e-8);
            let state = current
                .state_at_specific_enthalpy(e, bath, geometry, &mut || Ok(()))
                .unwrap();
            let capacity = if side == 0 {
                200.0 + 0.5 * t
            } else {
                400.0 + t
            };
            let derivative = geometry.area_m2 * k * t
                / ((geometry.half_distance_m + geometry.surface_resistance_m2_k_w * k * face)
                    * capacity);
            near(state.power_enthalpy_derivative_w_kg_j, derivative, 2e-11);
            let lo = bath_response(temperature(side, e - 0.25), bath, k, geometry).1;
            let hi = bath_response(temperature(side, e + 0.25), bath, k, geometry).1;
            near(
                state.power_enthalpy_derivative_w_kg_j,
                (hi - lo) / 0.5,
                3e-8,
            );
        }
        let e = energy(side, 500.0);
        for stale in [
            reference,
            ImplicitBoundaryGeometry {
                surface_resistance_m2_k_w: 2.0 * reference.surface_resistance_m2_k_w,
                ..geometry
            },
        ] {
            assert!(
                current
                    .state_at_specific_enthalpy(e, 350.0, stale, &mut || Ok(()))
                    .is_err()
            );
        }
        assert!(
            current
                .state_at_specific_enthalpy(energy(side, 710.0), 350.0, geometry, &mut || Ok(()),)
                .is_err()
        );
        assert!(
            current
                .state_at_specific_enthalpy_with_geometry_tangent(
                    e,
                    350.0,
                    geometry,
                    0.0,
                    &mut || Ok(()),
                )
                .unwrap()
                .is_none()
        );
        assert!(
            original
                .with_current_planar_geometry(CurrentConductiveBoundary {
                    area_m2: 0.0,
                    half_distance_m: 0.005,
                })
                .is_err()
        );
        let same = original
            .with_current_planar_geometry(CurrentConductiveBoundary {
                area_m2: reference.area_m2,
                half_distance_m: reference.half_distance_m,
            })
            .unwrap();
        assert_eq!(
            same.state_at_specific_enthalpy(e, 350.0, reference, &mut || Ok(()))
                .unwrap(),
            original
                .state_at_specific_enthalpy(e, 350.0, reference, &mut || Ok(()))
                .unwrap(),
        );
    }
}

#[test]
fn current_integrated_baths_drive_original_native_heat_and_separate_entropy() {
    let source = Fixture::new(true);
    let cv = source.cv();
    let conductivity =
        IsochoricConductionLaw::try_new(&source.card, &cv, &source.conductivity, SELECTION)
            .unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conductivity).unwrap();
    let (initial, mass, dt) = (320.0, 0.5, 40.0);
    let baths = Boundary {
        ambient_k: 650.0,
        htc_w_m2_k: 999.0,
        source_w: 0.0,
        reservoir: Some(Reservoir {
            area_m2: AREA,
            conductance_w_k: 999.0,
            temperature_k: 350.0,
        }),
    };
    let thermal = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(0),
                material_id: CV_MATERIAL,
                body_id: 1,
                segment_id: 1,
                mass_kg: mass,
                specific_h_j_kg: energy(0, initial),
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [mass / (RHO * AREA), 0.0, 0.0],
                    volume_m3: mass / RHO,
                    surface_area_m2: 3.0 * AREA,
                },
                boundary: baths,
            }),
        }],
        &[(CV_MATERIAL, &cv)],
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
    .0;
    let reference = [
        ImplicitBoundaryGeometry {
            area_m2: AREA,
            half_distance_m: 0.0125,
            surface_resistance_m2_k_w: 0.004,
        },
        ImplicitBoundaryGeometry {
            area_m2: 2.0 * AREA,
            half_distance_m: 0.02,
            surface_resistance_m2_k_w: 0.006,
        },
    ];
    let reservoir = IntegratedCaloricBoundary::new(&layer, reference[0], layer_controls()).unwrap();
    let ambient = IntegratedCaloricBoundary::new(&layer, reference[1], layer_controls()).unwrap();
    let current_r = reservoir
        .with_current_planar_geometry(CurrentConductiveBoundary {
            area_m2: 2.0 * AREA,
            half_distance_m: 0.00625,
        })
        .unwrap();
    let current_a = ambient
        .with_current_planar_geometry(CurrentConductiveBoundary {
            area_m2: 4.0 * AREA,
            half_distance_m: 0.01,
        })
        .unwrap();
    let observer = ConductionObserver::new(
        &thermal,
        &[(CV_MATERIAL, &conductivity)],
        &[],
        &[ConductiveReservoir {
            cell_id: CellId(0),
            half_distance_m: reference[0].half_distance_m,
            surface_resistance_m2_k_w: reference[0].surface_resistance_m2_k_w,
        }],
        || Ok(()),
    )
    .unwrap()
    .with_ambient_faces(
        &thermal,
        &[ConductiveAmbient {
            cell_id: CellId(0),
            half_distance_m: reference[1].half_distance_m,
            surface_resistance_m2_k_w: reference[1].surface_resistance_m2_k_w,
        }],
        || Ok(()),
    )
    .unwrap()
    .with_current_geometry(
        &thermal,
        CurrentConductionGeometry {
            reservoirs: [(CellId(0), current_r.current_geometry())].into(),
            ambient_faces: [(CellId(0), current_a.current_geometry())].into(),
            ..Default::default()
        },
        || Ok(()),
    )
    .unwrap();
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(CV_MATERIAL, &conductivity)], || {
            Ok(())
        })
        .unwrap()
        .with_boundary_laws(
            &thermal,
            &[(CellId(0), &current_r)],
            &[(CellId(0), &current_a)],
            || Ok(()),
        )
        .unwrap();
    let heat = |t| {
        [
            -dt * bath_response(t, 350.0, 0.04, current_r.native_geometry()).1,
            -dt * bath_response(t, 650.0, 0.04, current_a.native_geometry()).1,
        ]
    };
    let (mut lo, mut hi) = (initial, 650.0);
    for _ in 0..80 {
        let t = 0.5 * (lo + hi);
        if mass * (energy(0, t) - energy(0, initial)) - heat(t).iter().sum::<f64>() < 0.0 {
            lo = t;
        } else {
            hi = t;
        }
    }
    let expected = 0.5 * (lo + hi);
    let expected_heat = heat(expected);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let mut calls = 0;
    let (next, report) = chain
        .advance_to(&thermal, dt, controls(dt), &mut workspace, || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    near(next.temperature_k(CellId(0)).unwrap(), expected, 2e-8);
    near(
        report.energy.external_heat_j,
        expected_heat.iter().sum(),
        1e-6,
    );
    near(
        report.energy.stored_energy_change_j,
        mass * (energy(0, expected) - energy(0, initial)),
        2e-6,
    );
    near(
        report.entropy.reservoir_entropy_input_j_k,
        expected_heat[0] / 350.0 + expected_heat[1] / 650.0,
        2e-8,
    );
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(next.cells()[&CellId(0)].mass_kg, mass);
    assert_eq!(
        next.cells()[&CellId(0)].geometry,
        thermal.cells()[&CellId(0)].geometry
    );
    assert_eq!(next.cells()[&CellId(0)].boundary, baths);
    let mut cancelled = 0;
    assert!(
        chain
            .advance_to(&thermal, dt, controls(dt), &mut workspace, || {
                cancelled += 1;
                if cancelled == calls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    let (retry, retry_report) = chain
        .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
        .unwrap();
    assert_eq!(retry.cells(), next.cells());
    assert_eq!(retry_report, report);
    assert_eq!(thermal.time_s(), 0.0);
}
