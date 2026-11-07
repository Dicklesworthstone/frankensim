//! Density-aware radiation and actual surface sources on original Cp storage.
use super::*;
use crate::thermal::conduction::GrayLayerGap;
use crate::thermal::transport::conduction::{
    IntegratedCaloricGrayBoundary, IntegratedCaloricGrayLink, IntegratedCaloricSurfaceLink,
};
use fs_thermal::enthalpy::conduction::{
    ConductiveLink,
    implicit::{
        boundary::{
            ImplicitBoundaryExchangeLaw, ImplicitGrayExteriorLaw, ImplicitGrayExteriorState,
            ImplicitGraySurfaceLaw, ImplicitGraySurfaceState,
        },
        link::{ImplicitSurfaceExchangeLaw, ImplicitSurfaceExchangeState},
    },
};

fn face(expansion: &IsobaricFreeExpansion<'_>, h: f64, coating: f64) -> ImplicitBoundaryGeometry {
    ImplicitBoundaryGeometry {
        area_m2: AREA,
        half_distance_m: distance(expansion, h, 1.0).0,
        surface_resistance_m2_k_w: coating,
    }
}

fn pair_geometry(expansion: &IsobaricFreeExpansion<'_>, h: [f64; 2]) -> LayerGeometry {
    LayerGeometry::planar(
        std::array::from_fn(|i| distance(expansion, h[i], [0.4, 0.6][i]).0),
        AREA,
    )
    .unwrap()
}

fn pair_slopes(expansion: &IsobaricFreeExpansion<'_>, h: [f64; 2]) -> [f64; 2] {
    std::array::from_fn(|i| distance(expansion, h[i], [0.4, 0.6][i]).1)
}

#[test]
fn integrated_gray_gap_density_and_original_resistance_tangents_match_full_flux_differences() {
    let (card, cp, density, key) = specimen();
    let expansion = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    let original = IsobaricConductionLaw::try_new(&card, &cp, &key, POLICY).unwrap();
    let law = IsobaricExpansionConductionLaw::try_new(&expansion, &original).unwrap();
    let layer = CaloricConductionLayer::from_expanding_isobaric(&law).unwrap();
    let gap = GrayLayerGap::equal_area_full_view([0.7, 0.8], 0.0).unwrap();
    let source = Contact {
        caloric: &expansion,
        queries: Cell::new(0),
    };
    for temperatures in [[308.0, 294.0], [294.0, 308.0]] {
        let initial = timeline(&expansion, &temperatures, &[0.4, 0.6], None);
        let h = temperatures.map(|t| cp.specific_enthalpy_j_kg(t).unwrap());
        let link = IntegratedCaloricGrayLink::new_with_gap_resistance(
            &initial,
            7,
            [&layer; 2],
            pair_geometry(&expansion, h),
            gap,
            &source,
            flux_controls(),
            || Ok(()),
        )
        .unwrap();
        source.queries.set(0);
        let state = link
            .state_at_specific_enthalpies_with_geometry_tangent(
                h,
                link.native_geometry(),
                initial.contact_patch(7).unwrap(),
                pair_slopes(&expansion, h),
                &mut || Ok(()),
            )
            .unwrap()
            .unwrap();
        assert_eq!(source.queries.get(), 1);
        let power = |energy: [f64; 2]| {
            pair_geometry(&expansion, energy)
                .conduction([layer.conductivity_curve(); 2], 0.0)
                .unwrap()
                .evaluate_gray_gap(
                    energy.map(|e| cp.temperature_k(e).unwrap()),
                    gap.with_nonradiative_resistance(contact_resistance(energy))
                        .unwrap(),
                    flux_controls(),
                    || Ok(()),
                )
                .unwrap()
                .power_a_to_b_w
        };
        for side in 0..2 {
            let mut plus = h;
            let mut minus = h;
            plus[side] += 0.01;
            minus[side] -= 0.01;
            close(
                state.power_enthalpy_derivatives_w_kg_j[side],
                (power(plus) - power(minus)) / 0.02,
                2e-7,
            );
        }
        assert_eq!(
            link.state_at_specific_enthalpies(
                h,
                link.native_geometry(),
                initial.contact_patch(7).unwrap(),
                &mut || Ok(())
            )
            .unwrap(),
            link.state_at_specific_enthalpies_with_geometry_tangent(
                h,
                link.native_geometry(),
                initial.contact_patch(7).unwrap(),
                [0.0; 2],
                &mut || Ok(()),
            )
            .unwrap()
            .unwrap(),
        );
    }
}

struct Surface<'a> {
    caloric: &'a dyn ThermalLaw,
    queries: Cell<usize>,
}
impl ImplicitSurfaceExchangeLaw for Surface<'_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        [self.caloric; 2]
    }
    fn energy_bounds_j_kg(&self) -> [[f64; 2]; 2] {
        [[-20.0, 30.0]; 2]
    }
    fn surface_temperature_bounds_k(&self) -> [[f64; 2]; 2] {
        [[290.0, 310.0]; 2]
    }
    fn state_at_surfaces(
        &self,
        energies: [f64; 2],
        surfaces: [f64; 2],
    ) -> Result<ImplicitSurfaceExchangeState, DepositionError> {
        self.queries.set(self.queries.get() + 1);
        if energies.iter().any(|e| !(-20.0..=30.0).contains(e))
            || surfaces.iter().any(|t| !(290.0..=310.0).contains(t))
        {
            return Err(DepositionError::MaterialDomain);
        }
        let g = 10.0
            + 0.02 * (surfaces[0] - 300.0)
            + 0.01 * (surfaces[1] - 300.0)
            + 0.005 * energies[0]
            + 0.003 * energies[1];
        ImplicitSurfaceExchangeState::from_conductance(surfaces, g, [0.02, 0.01], [0.005, 0.003])
    }
}

#[test]
fn integrated_surface_source_keeps_accepted_sample_with_complete_density_chain_rule() {
    let (card, cp, density, key) = specimen();
    let expansion = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    let original = IsobaricConductionLaw::try_new(&card, &cp, &key, POLICY).unwrap();
    let law = IsobaricExpansionConductionLaw::try_new(&expansion, &original).unwrap();
    let layer = CaloricConductionLayer::from_expanding_isobaric(&law).unwrap();
    let source = Surface {
        caloric: &expansion,
        queries: Cell::new(0),
    };
    for temperatures in [[308.0, 294.0], [294.0, 308.0]] {
        let initial = timeline(&expansion, &temperatures, &[0.4, 0.6], None);
        let h = temperatures.map(|t| cp.specific_enthalpy_j_kg(t).unwrap());
        let link = IntegratedCaloricSurfaceLink::new(
            &initial,
            7,
            [&layer; 2],
            pair_geometry(&expansion, h),
            &source,
            flux_controls(),
            || Ok(()),
        )
        .unwrap();
        source.queries.set(0);
        let baseline = link.evaluate(h, || Ok(())).unwrap();
        let original_queries = source.queries.replace(0);
        assert!(original_queries > 0);
        let state = link
            .state_at_specific_enthalpies_with_geometry_tangent(
                h,
                link.native_geometry(),
                initial.contact_patch(7).unwrap(),
                pair_slopes(&expansion, h),
                &mut || Ok(()),
            )
            .unwrap()
            .unwrap();
        assert_eq!(
            source.queries.get(),
            original_queries,
            "the geometry tangent adds no source query or second surface solve"
        );
        assert_eq!(state.conductance_w_k, baseline.conductance_w_k);
        let power = |energy: [f64; 2]| {
            IntegratedCaloricSurfaceLink::new(
                &initial,
                7,
                [&layer; 2],
                pair_geometry(&expansion, energy),
                &source,
                flux_controls(),
                || Ok(()),
            )
            .unwrap()
            .evaluate(energy, || Ok(()))
            .unwrap()
            .power_a_to_b_w
        };
        for side in 0..2 {
            let mut plus = h;
            let mut minus = h;
            plus[side] += 0.01;
            minus[side] -= 0.01;
            close(
                state.power_enthalpy_derivatives_w_kg_j[side],
                (power(plus) - power(minus)) / 0.02,
                2e-7,
            );
        }
    }
}

struct BulkSurface<'a> {
    caloric: &'a dyn ThermalLaw,
    queries: Cell<usize>,
}
impl ImplicitGraySurfaceLaw for BulkSurface<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [-20.0, 30.0]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitGraySurfaceState, DepositionError> {
        self.queries.set(self.queries.get() + 1);
        Ok(ImplicitGraySurfaceState {
            emissivity: 0.4 + 0.001 * h,
            emissivity_enthalpy_derivative_kg_j: 0.001,
            convection_w_m2_k: 3.0 + 0.01 * h,
            convection_enthalpy_derivative_kg_w_m2_k_j: 0.01,
        })
    }
}

#[test]
fn coated_gray_boundary_combines_density_with_one_original_bulk_surface_sample() {
    let (card, cp, density, key) = specimen();
    let expansion = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    let original = IsobaricConductionLaw::try_new(&card, &cp, &key, POLICY).unwrap();
    let law = IsobaricExpansionConductionLaw::try_new(&expansion, &original).unwrap();
    let layer = CaloricConductionLayer::from_expanding_isobaric(&law).unwrap();
    let source = BulkSurface {
        caloric: &expansion,
        queries: Cell::new(0),
    };
    for temperature in [295.0, 305.0] {
        let h = cp.specific_enthalpy_j_kg(temperature).unwrap();
        for coating in [0.0, 0.02, 0.3] {
            let geometry = face(&expansion, h, coating);
            let boundary = IntegratedCaloricGrayBoundary::new(
                &layer,
                geometry,
                0.0,
                0.0,
                315.0,
                flux_controls(),
            )
            .unwrap()
            .with_surface_law(&source)
            .unwrap();
            source.queries.set(0);
            let state = boundary
                .state_at_specific_enthalpy_with_geometry_tangent(
                    h,
                    285.0,
                    geometry,
                    distance(&expansion, h, 1.0).1,
                    &mut || Ok(()),
                )
                .unwrap()
                .unwrap();
            assert_eq!(source.queries.get(), 1);
            assert!(state.baths[0].outward_power_w > 0.0);
            assert!(state.baths[1].outward_power_w < 0.0);
            assert_eq!(state.baths.map(|bath| bath.temperature_k), [285.0, 315.0]);
            let power = |energy| {
                IntegratedCaloricGrayBoundary::new(
                    &layer,
                    face(&expansion, energy, coating),
                    0.0,
                    0.0,
                    315.0,
                    flux_controls(),
                )
                .unwrap()
                .with_surface_law(&source)
                .unwrap()
                .evaluate(energy, 285.0, &mut || Ok(()))
                .unwrap()
                .outward_power_w
            };
            close(
                state.power_enthalpy_derivative_w_kg_j,
                (power(h + 0.01) - power(h - 0.01)) / 0.02,
                5e-8,
            );
            assert_eq!(
                boundary
                    .state_at_specific_enthalpy(h, 285.0, geometry, &mut || Ok(()))
                    .unwrap(),
                boundary
                    .state_at_specific_enthalpy_with_geometry_tangent(
                        h,
                        285.0,
                        geometry,
                        0.0,
                        &mut || Ok(()),
                    )
                    .unwrap()
                    .unwrap(),
            );
        }
    }
}

struct Exterior<'a> {
    caloric: &'a dyn ThermalLaw,
    negative_feedback: bool,
    queries: Cell<usize>,
    last: Cell<[f64; 3]>,
}
impl ImplicitGrayExteriorLaw for Exterior<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }
    fn energy_bounds_j_kg(&self) -> [f64; 2] {
        [-20.0, 30.0]
    }
    fn surface_temperature_bounds_k(&self) -> [f64; 2] {
        [290.0, 310.0]
    }
    fn state_at_surface(
        &self,
        h: f64,
        surface: f64,
        gas: f64,
    ) -> Result<ImplicitGrayExteriorState, DepositionError> {
        self.queries.set(self.queries.get() + 1);
        self.last.set([h, surface, gas]);
        if !(-20.0..=30.0).contains(&h)
            || !(290.0..=310.0).contains(&surface)
            || !(280.0..=320.0).contains(&gas)
        {
            return Err(DepositionError::MaterialDomain);
        }
        let (base, decay) = if self.negative_feedback {
            (0.1, 0.2)
        } else {
            (3.0, -0.02)
        };
        let htc = base * (0.001 * h - decay * (surface - 300.0)).exp();
        let epsilon = if self.negative_feedback {
            0.0
        } else {
            0.4 + 0.001 * h + 0.002 * (surface - 300.0)
        };
        Ok(ImplicitGrayExteriorState {
            coefficients: ImplicitGraySurfaceState {
                emissivity: epsilon,
                emissivity_enthalpy_derivative_kg_j: if self.negative_feedback {
                    0.0
                } else {
                    0.001
                },
                convection_w_m2_k: htc,
                convection_enthalpy_derivative_kg_w_m2_k_j: 0.001 * htc,
            },
            emissivity_temperature_derivative_per_k: if self.negative_feedback {
                0.0
            } else {
                0.002
            },
            convection_temperature_derivative_w_m2_k2: -decay * htc,
        })
    }
}

#[test]
fn density_coating_and_actual_exterior_feedback_share_the_accepted_nonlinear_balance() {
    let (card, cp, density, key) = specimen();
    let expansion = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    let original = IsobaricConductionLaw::try_new(&card, &cp, &key, POLICY).unwrap();
    let law = IsobaricExpansionConductionLaw::try_new(&expansion, &original).unwrap();
    let layer = CaloricConductionLayer::from_expanding_isobaric(&law).unwrap();
    let h = cp.specific_enthalpy_j_kg(305.0).unwrap();
    for negative_feedback in [false, true] {
        let source = Exterior {
            caloric: &expansion,
            negative_feedback,
            queries: Cell::new(0),
            last: Cell::new([0.0; 3]),
        };
        for coating in [0.0, 0.02, 1.0] {
            let geometry = face(&expansion, h, coating);
            let boundary = IntegratedCaloricGrayBoundary::new(
                &layer,
                geometry,
                0.0,
                0.0,
                315.0,
                flux_controls(),
            )
            .unwrap()
            .with_exterior_law(&source)
            .unwrap();
            source.queries.set(0);
            let baseline = boundary.evaluate(h, 285.0, &mut || Ok(())).unwrap();
            let original_queries = source.queries.replace(0);
            let state = boundary
                .state_at_specific_enthalpy_with_geometry_tangent(
                    h,
                    285.0,
                    geometry,
                    distance(&expansion, h, 1.0).1,
                    &mut || Ok(()),
                )
                .unwrap()
                .unwrap();
            assert_eq!(source.queries.get(), original_queries);
            assert_eq!(
                source.last.get(),
                [h, baseline.surface_temperature_k, 285.0]
            );
            if negative_feedback {
                assert!(
                    state.power_enthalpy_derivative_w_kg_j < 0.0,
                    "passive heat can have a negative differential response"
                );
            }
            let power = |energy| {
                IntegratedCaloricGrayBoundary::new(
                    &layer,
                    face(&expansion, energy, coating),
                    0.0,
                    0.0,
                    315.0,
                    flux_controls(),
                )
                .unwrap()
                .with_exterior_law(&source)
                .unwrap()
                .evaluate(energy, 285.0, &mut || Ok(()))
                .unwrap()
                .outward_power_w
            };
            close(
                state.power_enthalpy_derivative_w_kg_j,
                (power(h + 0.01) - power(h - 0.01)) / 0.02,
                5e-8,
            );
        }
    }
}

#[test]
fn density_integrated_gray_boundary_keeps_actual_bath_entropy_and_enthalpy_heat() {
    let (card, cp, density, key) = specimen();
    let expansion = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    let original = IsobaricConductionLaw::try_new(&card, &cp, &key, POLICY).unwrap();
    let law = IsobaricExpansionConductionLaw::try_new(&expansion, &original).unwrap();
    let layer = CaloricConductionLayer::from_expanding_isobaric(&law).unwrap();
    let initial = timeline(&expansion, &[305.0], &[1.0], Some(285.0));
    let h = initial.cells()[&CellId(0)].specific_h_j_kg;
    let geometry = face(&expansion, h, 0.02);
    let boundary =
        IntegratedCaloricGrayBoundary::new(&layer, geometry, 0.6, 3.0, 315.0, flux_controls())
            .unwrap();
    let observer = ConductionObserver::new(
        &initial,
        &[(MaterialId(0), &law)],
        &[],
        &[ConductiveReservoir {
            cell_id: CellId(0),
            half_distance_m: geometry.half_distance_m,
            surface_resistance_m2_k_w: geometry.surface_resistance_m2_k_w,
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
    let graph = ImplicitConductionGraph::new(
        &initial,
        &observer,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_boundary_exchanges(&initial, &[(CellId(0), &boundary)], &[], || Ok(()))
    .unwrap();
    let dt = 0.5;
    let (next, report) = graph
        .advance_to(
            &initial,
            dt,
            step_controls(dt),
            &mut ImplicitConductionWorkspace::new(1).unwrap(),
            || Ok(()),
        )
        .unwrap();
    let final_h = next.cells()[&CellId(0)].specific_h_j_kg;
    let exchange = IntegratedCaloricGrayBoundary::new(
        &layer,
        face(&expansion, final_h, 0.02),
        0.6,
        3.0,
        315.0,
        flux_controls(),
    )
    .unwrap()
    .evaluate(final_h, 285.0, &mut || Ok(()))
    .unwrap();
    close(final_h - h, -dt * exchange.outward_power_w, 1e-8);
    close(report.energy.external_heat_j, final_h - h, 1e-8);
    close(
        report.entropy.reservoir_entropy_input_j_k,
        -dt * exchange
            .boundary
            .baths
            .iter()
            .map(|bath| bath.outward_power_w / bath.temperature_k)
            .sum::<f64>(),
        1e-10,
    );
    close(report.energy.residual_j, 0.0, 1e-8);
    assert_eq!(next.cells()[&CellId(0)].mass_kg, 1.0);
    assert_eq!(
        next.cells()[&CellId(0)].geometry,
        initial.cells()[&CellId(0)].geometry
    );
    let work = expansion
        .advance_heat(305.0, 1.0, report.energy.external_heat_j)
        .unwrap();
    close(
        work.internal_energy_change_j + work.boundary_work_j,
        report.energy.external_heat_j,
        1e-8,
    );
    assert_ne!(work.boundary_work_j, 0.0);
}

#[test]
fn density_radiative_and_surface_links_preserve_closed_graph_energy_and_inventory() {
    let (card, cp, density, key) = specimen();
    let expansion = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    let original = IsobaricConductionLaw::try_new(&card, &cp, &key, POLICY).unwrap();
    let law = IsobaricExpansionConductionLaw::try_new(&expansion, &original).unwrap();
    let layer = CaloricConductionLayer::from_expanding_isobaric(&law).unwrap();
    let initial = timeline(&expansion, &[308.0, 294.0], &[0.4, 0.6], None);
    let h = [
        initial.cells()[&CellId(0)].specific_h_j_kg,
        initial.cells()[&CellId(1)].specific_h_j_kg,
    ];
    let geometry = pair_geometry(&expansion, h);
    let source = Surface {
        caloric: &expansion,
        queries: Cell::new(0),
    };
    let gray = IntegratedCaloricGrayLink::new(
        &initial,
        7,
        [&layer; 2],
        geometry,
        GrayLayerGap::equal_area_full_view([0.7, 0.8], 10.0).unwrap(),
        flux_controls(),
        || Ok(()),
    )
    .unwrap();
    let surface = IntegratedCaloricSurfaceLink::new(
        &initial,
        7,
        [&layer; 2],
        geometry,
        &source,
        flux_controls(),
        || Ok(()),
    )
    .unwrap();
    let native = ConductiveLink {
        link_id: 7,
        half_distances_m: geometry.equivalent_distances_m(),
        contact_resistance_m2_k_w: 0.0,
    };
    let observer = ConductionObserver::new(
        &initial,
        &[(MaterialId(0), &law)],
        &[native],
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
    for link in [
        &gray as &dyn ImplicitLinkLaw,
        &surface as &dyn ImplicitLinkLaw,
    ] {
        let graph = ImplicitConductionGraph::new(
            &initial,
            &observer,
            &[(MaterialId(0), &law)],
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .unwrap()
        .with_link_laws(&initial, &[(7, link)], || Ok(()))
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
        assert!(report.entropy.conduction_generation_j_k > 0.0);
        assert!(next.temperature_k(CellId(0)).unwrap() < 308.0);
        assert!(next.temperature_k(CellId(1)).unwrap() > 294.0);
        for (&id, cell) in next.cells() {
            assert_eq!(cell.mass_kg, initial.cells()[&id].mass_kg);
            assert_eq!(cell.geometry, initial.cells()[&id].geometry);
        }
    }
}
