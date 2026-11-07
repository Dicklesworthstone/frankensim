//! Source-exact density, integrated flux and conserved native Cp evolution.
use super::*;
use crate::thermal::conduction::{LayerConductionConfig, LayerGeometry};
use crate::thermal::transport::{
    CONDUCTIVITY,
    conduction::{CaloricConductionLayer, IntegratedCaloricBoundary, IntegratedCaloricLink},
    isobaric::{
        conduction::IsobaricConductionLaw, expansion::conduction::IsobaricExpansionConductionLaw,
    },
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, Reservoir, StepConfig, ThermalLaw, ThermalLink, ThermalStorageBasis, TimedEvent,
    conduction::{
        ConductionObserver, ConductiveReservoir, PlanarDensityLayer,
        implicit::{
            ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionGraph,
            ImplicitConductionLaw, ImplicitConductionLinearControls, ImplicitConductionWorkspace,
            boundary::{ImplicitBoundaryGeometry, ImplicitBoundaryLaw},
            contact::{ImplicitContactLaw, ImplicitContactState},
            link::ImplicitLinkLaw,
        },
    },
};
use std::cell::Cell;

const AREA: f64 = 0.02;
const RHO: f64 = 1000.0;
const K: [(f64, f64); 3] = [(290.0, 1.0), (303.0, 4.0), (310.0, 2.0)];

fn specimen() -> (
    MaterialCard,
    IsobaricHeatCapacityLaw,
    PreparedFreeExpansionDensity,
    PropertyKey,
) {
    let keys = keys();
    let conductivity = PropertyKey::new("normal-conductivity", CONDUCTIVITY);
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim(&keys[0], &CP, true)).unwrap();
    claims.insert_claim(claim(&conductivity, &K, true)).unwrap();
    let mut density = claim(&keys[1], &[(300.0, RHO), (300.0, RHO)], true);
    density.value = PropertyValue::Scalar {
        value: RHO,
        dims: keys[1].quantity().dims(),
    };
    density.interpolation = InterpolationPolicy::ConstantWithinValidity;
    claims.insert_claim(density).unwrap();
    // Original free expansion is normal only: the native fixed footprint is
    // physically compatible, without imposing lateral constraints or stress.
    for (key, knots) in [
        (
            &keys[2],
            &[(290.0, -0.01), (302.0, 0.002), (310.0, 0.018)][..],
        ),
        (&keys[3], &[(290.0, 0.0), (310.0, 0.0)][..]),
        (&keys[4], &[(290.0, 0.0), (310.0, 0.0)][..]),
    ] {
        claims.insert_claim(claim(key, knots, true)).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic expansion/conduction specimen".into(),
            phase: "single phase".into(),
            process: "maintained pressure, normal free expansion".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let (cp, density) = prepare(&card, 3, 290.0, 310.0, 300.0, PRESSURE);
    (card, cp, density, conductivity)
}

fn flux_controls() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 96,
        absolute_heat_tolerance_w: 1e-11,
        relative_heat_tolerance: 1e-12,
    }
}

fn step_controls(dt: f64) -> ImplicitConductionControls {
    ImplicitConductionControls {
        maximum_dt_s: dt,
        maximum_halvings: 0,
        cell_tolerance_j: 1e-9,
        total_tolerance_j: 1e-8,
        ..ImplicitConductionControls::default()
    }
}

fn distance(expansion: &IsobaricFreeExpansion<'_>, h: f64, mass: f64) -> (f64, f64) {
    let rho = expansion.bulk_density_kg_m3(h).unwrap().unwrap();
    (
        mass / (2.0 * AREA * rho),
        rho * expansion
            .specific_volume_derivative_at_specific_enthalpy(h)
            .unwrap(),
    )
}

fn timeline<'a>(
    expansion: &'a IsobaricFreeExpansion<'a>,
    temperatures: &[f64],
    masses: &[f64],
    bath: Option<f64>,
) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for (i, (&temperature, &mass)) in temperatures.iter().zip(masses).enumerate() {
        let h = expansion
            .caloric()
            .specific_enthalpy_j_kg(temperature)
            .unwrap();
        let volume = mass / expansion.bulk_density_kg_m3(h).unwrap().unwrap();
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                material_id: MaterialId(0),
                body_id: 0,
                segment_id: i as u64,
                mass_kg: mass,
                specific_h_j_kg: h,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + volume / AREA, 0.0, 0.0],
                    volume_m3: volume,
                    surface_area_m2: 2.0 * AREA,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: if i == 0 {
                        bath.map(|temperature_k| Reservoir {
                            temperature_k,
                            area_m2: AREA,
                            conductance_w_k: 0.0,
                        })
                    } else {
                        None
                    },
                },
            }),
        });
    }
    if temperatures.len() == 2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 7,
                a: CellId(0),
                b: CellId(1),
                area_m2: AREA,
                conductance_w_k: 0.0,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), expansion)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 128,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

#[test]
fn original_density_tangent_uses_current_cp_and_source_breakpoint_conventions() {
    let (card, cp, density, key) = specimen();
    let expansion = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    let original = IsobaricConductionLaw::try_new(&card, &cp, &key, POLICY).unwrap();
    let law = IsobaricExpansionConductionLaw::try_new(&expansion, &original).unwrap();
    let layer = CaloricConductionLayer::from_expanding_isobaric(&law).unwrap();
    assert!(core::ptr::addr_eq(
        layer.caloric_law(),
        &expansion as &dyn ThermalLaw
    ));
    assert!(core::ptr::eq(layer.source_line(), original.source_line()));
    assert!(layer.supports_material_card(card.content_hash()));
    assert_eq!(expansion.storage_basis(), ThermalStorageBasis::Enthalpy);
    assert_ne!(expansion.caloric_identity(), cp.caloric_identity());
    assert_eq!(cp.bulk_density_kg_m3(0.0).unwrap(), None);
    for t in [290.0, 295.0, 300.0, 302.0, 305.0, 310.0] {
        let h = cp.specific_enthalpy_j_kg(t).unwrap();
        let cp_current = if t <= 300.0 {
            2.0
        } else {
            2.0 + 0.2 * (t - 300.0)
        };
        let expansion_slope = if t < 302.0 { 0.001 } else { 0.002 };
        let dv = law
            .specific_volume_derivative_at_specific_enthalpy(h)
            .unwrap()
            .unwrap();
        close(dv, expansion_slope / (RHO * cp_current), 1e-19);
        assert_eq!(
            ImplicitConductionLaw::state_at_specific_enthalpy(&law, h).unwrap(),
            ImplicitConductionLaw::state_at_specific_enthalpy(&original, h).unwrap(),
        );
        assert_eq!(
            expansion.specific_entropy_change_j_kg_k(0.0, h).unwrap(),
            cp.specific_entropy_change_j_kg_k(0.0, h).unwrap(),
        );
        if t != 290.0 && t != 300.0 && t != 302.0 && t != 310.0 {
            let epsilon = 1e-3;
            let volume = |energy| {
                expansion
                    .bulk_density_kg_m3(energy)
                    .unwrap()
                    .unwrap()
                    .recip()
            };
            close(
                dv,
                (volume(h + epsilon) - volume(h - epsilon)) / (2.0 * epsilon),
                1e-14,
            );
        }
    }
    // A byte-equivalent clone is still a different original Cp object.
    let other_cp = cp.clone();
    let other_k = IsobaricConductionLaw::try_new(&card, &other_cp, &key, POLICY).unwrap();
    assert!(IsobaricExpansionConductionLaw::try_new(&expansion, &other_k).is_err());
    // The density bridge also retains the actual cubic free-volume response,
    // including zero and negative expansion; no positive-expansion clipping.
    for slope in [0.0, -0.001] {
        let lengths = [(290.0, -10.0 * slope), (310.0, 10.0 * slope)];
        let (card, cp, density) = fixture(&CP, &[&lengths], true, 300.0);
        let path = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
        let h = cp.specific_enthalpy_j_kg(305.0).unwrap();
        let s = 1.0 + 5.0 * slope;
        close(
            path.specific_volume_derivative_at_specific_enthalpy(h)
                .unwrap(),
            (3.0 * s * s * slope) / (RHO * 3.0),
            1e-19,
        );
    }
}

#[test]
fn integrated_boundary_density_tangent_matches_full_original_profile_differences() {
    let (card, cp, density, key) = specimen();
    let expansion = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    let original = IsobaricConductionLaw::try_new(&card, &cp, &key, POLICY).unwrap();
    let law = IsobaricExpansionConductionLaw::try_new(&expansion, &original).unwrap();
    let layer = CaloricConductionLayer::from_expanding_isobaric(&law).unwrap();
    for t in [295.0, 305.0] {
        let h = cp.specific_enthalpy_j_kg(t).unwrap();
        for bath in [292.0, 308.0] {
            for film in [0.0, 0.1] {
                let geometry = |energy| ImplicitBoundaryGeometry {
                    area_m2: AREA,
                    half_distance_m: distance(&expansion, energy, 1.0).0,
                    surface_resistance_m2_k_w: film,
                };
                let boundary =
                    IntegratedCaloricBoundary::new(&layer, geometry(h), flux_controls()).unwrap();
                let state = boundary
                    .state_at_specific_enthalpy_with_geometry_tangent(
                        h,
                        bath,
                        geometry(h),
                        distance(&expansion, h, 1.0).1,
                        &mut || Ok(()),
                    )
                    .unwrap()
                    .unwrap();
                let power = |energy| {
                    IntegratedCaloricBoundary::new(&layer, geometry(energy), flux_controls())
                        .unwrap()
                        .response_at_specific_energy(energy, bath, &mut || Ok(()))
                        .unwrap()
                        .power_to_reservoir_w
                };
                close(
                    state.power_enthalpy_derivative_w_kg_j,
                    (power(h + 1e-3) - power(h - 1e-3)) / 2e-3,
                    2e-6,
                );
                assert_eq!(
                    boundary
                        .state_at_specific_enthalpy(h, bath, geometry(h), &mut || Ok(()))
                        .unwrap(),
                    boundary
                        .state_at_specific_enthalpy_with_geometry_tangent(
                            h,
                            bath,
                            geometry(h),
                            0.0,
                            &mut || Ok(()),
                        )
                        .unwrap()
                        .unwrap(),
                );
                assert!(
                    boundary
                        .state_at_specific_enthalpy(h, bath, geometry(h + 0.1), &mut || Ok(()))
                        .is_err()
                );
            }
        }
    }

    let geometry = ImplicitBoundaryGeometry {
        area_m2: AREA,
        half_distance_m: distance(&expansion, 0.0, 1.0).0,
        surface_resistance_m2_k_w: 0.1,
    };
    let boundary = IntegratedCaloricBoundary::new(&layer, geometry, flux_controls()).unwrap();
    let mut calls = 0;
    let cancelled = boundary.state_at_specific_enthalpy_with_geometry_tangent(
        0.0,
        305.0,
        geometry,
        distance(&expansion, 0.0, 1.0).1,
        &mut || {
            calls += 1;
            if calls == 3 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(cancelled, Err(DepositionError::Cancelled)));
}

struct Contact<'a> {
    caloric: &'a dyn ThermalLaw,
    queries: Cell<usize>,
}
fn contact_resistance(h: [f64; 2]) -> f64 {
    0.1 + 1e-4 * (h[0] + h[1])
}
impl ImplicitContactLaw for Contact<'_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        [self.caloric; 2]
    }
    fn enthalpy_bounds_j_kg(&self) -> [[f64; 2]; 2] {
        [[-20.0, 30.0]; 2]
    }
    fn state_at_specific_enthalpies(
        &self,
        h: [f64; 2],
    ) -> Result<ImplicitContactState, DepositionError> {
        self.queries.set(self.queries.get() + 1);
        Ok(ImplicitContactState {
            resistance_m2_k_w: contact_resistance(h),
            resistance_enthalpy_derivatives_m2_kg_k_w_j: [1e-4; 2],
        })
    }
}

#[test]
fn integrated_link_composes_density_and_contact_derivatives_in_one_original_solve() {
    let (card, cp, density, key) = specimen();
    let expansion = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    let original = IsobaricConductionLaw::try_new(&card, &cp, &key, POLICY).unwrap();
    let law = IsobaricExpansionConductionLaw::try_new(&expansion, &original).unwrap();
    let layer = CaloricConductionLayer::from_expanding_isobaric(&law).unwrap();
    let masses = [0.4, 0.6];
    let initial = timeline(&expansion, &[308.0, 294.0], &masses, None);
    let h = [
        initial.cells()[&CellId(0)].specific_h_j_kg,
        initial.cells()[&CellId(1)].specific_h_j_kg,
    ];
    let distances = |energies: [f64; 2]| {
        std::array::from_fn(|i| distance(&expansion, energies[i], masses[i]).0)
    };
    let slopes = std::array::from_fn(|i| distance(&expansion, h[i], masses[i]).1);
    let contact = Contact {
        caloric: &expansion,
        queries: Cell::new(0),
    };
    let link = IntegratedCaloricLink::new_with_contact(
        &initial,
        7,
        [&layer; 2],
        LayerGeometry::planar(distances(h), AREA).unwrap(),
        &contact,
        flux_controls(),
        || Ok(()),
    )
    .unwrap();
    contact.queries.set(0);
    let state = link
        .state_at_specific_enthalpies_with_geometry_tangent(
            h,
            link.native_geometry(),
            initial.contact_patch(7).unwrap(),
            slopes,
            &mut || Ok(()),
        )
        .unwrap()
        .unwrap();
    assert_eq!(contact.queries.get(), 1);
    let power = |energy: [f64; 2]| {
        LayerGeometry::planar(distances(energy), AREA)
            .unwrap()
            .conduction([layer.conductivity_curve(); 2], contact_resistance(energy))
            .unwrap()
            .evaluate(
                energy.map(|e| cp.temperature_k(e).unwrap()),
                flux_controls(),
                || Ok(()),
            )
            .unwrap()
            .power_a_to_b_w
    };
    for side in 0..2 {
        let mut plus = h;
        let mut minus = h;
        plus[side] += 1e-3;
        minus[side] -= 1e-3;
        close(
            state.power_enthalpy_derivatives_w_kg_j[side],
            (power(plus) - power(minus)) / 2e-3,
            2e-6,
        );
    }
}

// Independent source primitives: no production caloric inverse or face solve.
fn enthalpy(t: f64) -> f64 {
    let x = t - 300.0;
    2.0 * x + if x > 0.0 { 0.1 * x * x } else { 0.0 }
}
fn stretch(t: f64) -> f64 {
    if t <= 302.0 {
        1.0 + 0.001 * (t - 300.0)
    } else {
        1.002 + 0.002 * (t - 302.0)
    }
}
fn kirchhoff(t: f64) -> f64 {
    let x = t.min(303.0) - 290.0;
    let y = (t - 303.0).max(0.0);
    x + (3.0 / 26.0) * x * x + 4.0 * y - y * y / 7.0
}

#[test]
fn native_expanding_cp_boundary_matches_energy_oracle_without_double_pressure_work() {
    let (card, cp, density, key) = specimen();
    let expansion = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    let original = IsobaricConductionLaw::try_new(&card, &cp, &key, POLICY).unwrap();
    let law = IsobaricExpansionConductionLaw::try_new(&expansion, &original).unwrap();
    let layer = CaloricConductionLayer::from_expanding_isobaric(&law).unwrap();
    for (temperature, bath) in [(305.0, 290.0), (295.0, 310.0)] {
        let initial = timeline(&expansion, &[temperature], &[1.0], Some(bath));
        let h = initial.cells()[&CellId(0)].specific_h_j_kg;
        let geometry = ImplicitBoundaryGeometry {
            area_m2: AREA,
            half_distance_m: distance(&expansion, h, 1.0).0,
            surface_resistance_m2_k_w: 0.0,
        };
        let boundary = IntegratedCaloricBoundary::new(&layer, geometry, flux_controls()).unwrap();
        let observer = ConductionObserver::new(
            &initial,
            &[(MaterialId(0), &law)],
            &[],
            &[ConductiveReservoir {
                cell_id: CellId(0),
                half_distance_m: geometry.half_distance_m,
                surface_resistance_m2_k_w: 0.0,
            }],
            || Ok(()),
        )
        .unwrap()
        .with_material_density_layers(
            &initial,
            &[PlanarDensityLayer {
                cell_id: CellId(0),
                footprint_area_m2: AREA,
            }],
            || Ok(()),
        )
        .unwrap();
        let chain =
            ImplicitConductionChain::new(&initial, &observer, &[(MaterialId(0), &law)], || Ok(()))
                .unwrap()
                .with_boundary_laws(&initial, &[(CellId(0), &boundary)], &[], || Ok(()))
                .unwrap();
        let dt = 0.5;
        let (next, report) = chain
            .advance_to(
                &initial,
                dt,
                step_controls(dt),
                &mut ImplicitConductionWorkspace::new(1).unwrap(),
                || Ok(()),
            )
            .unwrap();
        let mut lo = temperature.min(bath);
        let mut hi = temperature.max(bath);
        for _ in 0..96 {
            let t = 0.5 * (lo + hi);
            let d = stretch(t) / (2.0 * AREA * RHO);
            let residual = enthalpy(t) - enthalpy(temperature)
                + dt * (AREA / d) * (kirchhoff(t) - kirchhoff(bath));
            if residual > 0.0 {
                hi = t;
            } else {
                lo = t;
            }
        }
        let expected = 0.5 * (lo + hi);
        close(next.temperature_k(CellId(0)).unwrap(), expected, 2e-8);
        close(
            report.energy.external_heat_j,
            enthalpy(expected) - enthalpy(temperature),
            5e-8,
        );
        close(report.energy.residual_j, 0.0, 1e-8);
        let first_law = expansion
            .advance_heat(temperature, 1.0, report.energy.external_heat_j)
            .unwrap();
        close(first_law.state.temperature_k(), expected, 2e-8);
        assert_ne!(first_law.boundary_work_j, 0.0);
        close(
            first_law.internal_energy_change_j + first_law.boundary_work_j,
            report.energy.external_heat_j,
            1e-8,
        );
        assert_eq!(next.cells()[&CellId(0)].mass_kg, 1.0);
        assert_eq!(
            next.cells()[&CellId(0)].geometry,
            initial.cells()[&CellId(0)].geometry
        );
        assert_eq!(report.energy.substeps, 1);
    }
}

#[test]
fn native_integrated_density_link_preserves_closed_inventory_in_graph_sdirk() {
    let (card, cp, density, key) = specimen();
    let expansion = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    let original = IsobaricConductionLaw::try_new(&card, &cp, &key, POLICY).unwrap();
    let law = IsobaricExpansionConductionLaw::try_new(&expansion, &original).unwrap();
    let layer = CaloricConductionLayer::from_expanding_isobaric(&law).unwrap();
    let masses = [0.4, 0.6];
    let initial = timeline(&expansion, &[308.0, 294.0], &masses, None);
    let distances = std::array::from_fn(|i| {
        distance(
            &expansion,
            initial.cells()[&CellId(i as u64)].specific_h_j_kg,
            masses[i],
        )
        .0
    });
    let link = IntegratedCaloricLink::new(
        &initial,
        7,
        [&layer; 2],
        LayerGeometry::planar(distances, AREA).unwrap(),
        0.02,
        flux_controls(),
        || Ok(()),
    )
    .unwrap();
    let observer = ConductionObserver::new(
        &initial,
        &[(MaterialId(0), &law)],
        &[link.native_geometry()],
        &[],
        || Ok(()),
    )
    .unwrap()
    .with_material_density_layers(
        &initial,
        &[
            PlanarDensityLayer {
                cell_id: CellId(0),
                footprint_area_m2: AREA,
            },
            PlanarDensityLayer {
                cell_id: CellId(1),
                footprint_area_m2: AREA,
            },
        ],
        || Ok(()),
    )
    .unwrap();
    let graph = ImplicitConductionGraph::new(
        &initial,
        &observer,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_link_laws(&initial, &[(7, &link)], || Ok(()))
    .unwrap();
    let (next, report) = graph
        .advance_sdirk2_to(
            &initial,
            0.05,
            step_controls(0.05),
            &mut ImplicitConductionWorkspace::new(2).unwrap(),
            || Ok(()),
        )
        .unwrap();
    close(next.energy_j(), initial.energy_j(), 1e-8);
    close(report.energy.external_heat_j, 0.0, 0.0);
    close(report.energy.residual_j, 0.0, 1e-8);
    assert!(next.temperature_k(CellId(0)).unwrap() < 308.0);
    assert!(next.temperature_k(CellId(1)).unwrap() > 294.0);
    assert!(report.entropy.conduction_generation_j_k > 0.0);
    for (&id, cell) in next.cells() {
        assert_eq!(cell.mass_kg, initial.cells()[&id].mass_kg);
        assert_eq!(cell.geometry, initial.cells()[&id].geometry);
    }
}
