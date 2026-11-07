//! Density-driven thickness in the original gray surface and gap balances.
use super::*;
use crate::enthalpy::conduction::implicit::boundary::{
    GrayAmbientBoundary, GrayAmbientControls, GrayTwoBathBoundary, ImplicitGraySurfaceLaw,
    ImplicitGraySurfaceState,
};
use crate::enthalpy::conduction::implicit::link::{GrayRadiativeControls, GrayRadiativeLink};
use std::cell::Cell;

struct Plateau;
static PLATEAU: Plateau = Plateau;
impl ThermalLaw for Plateau {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        EXPANDING.temperature_k(h)?;
        Ok(if h < 300.0 {
            h
        } else if h <= 500.0 {
            300.0
        } else {
            h - 200.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn bulk_density_kg_m3(&self, h: f64) -> Result<Option<f64>, DepositionError> {
        EXPANDING.bulk_density_kg_m3(h)
    }
}
impl ImplicitConductionLaw for Plateau {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        &EXPANDING
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [100.0, 1000.0]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: self.temperature_k(h)?,
            temperature_slope_kg_k_j: if (300.0..500.0).contains(&h) {
                0.0
            } else {
                1.0
            },
            conductivity_w_m_k: 1.0,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
    fn specific_volume_derivative_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.temperature_k(h)?;
        Ok(Some(1.0 / 300.0))
    }
}

fn face(h: f64) -> ImplicitBoundaryGeometry {
    ImplicitBoundaryGeometry {
        area_m2: 0.3,
        half_distance_m: 0.1 * h / 400.0,
        surface_resistance_m2_k_w: 0.02,
    }
}

#[test]
fn latent_gray_boundary_geometry_matches_analytic_convection_and_radiation_differences() {
    let h = 400.0;
    for emissivity in [0.0, 0.8] {
        let boundary =
            GrayAmbientBoundary::new(&PLATEAU, emissivity, 5.0, GrayAmbientControls::default())
                .unwrap();
        for bath in [200.0, 400.0] {
            let state = boundary
                .state_at_specific_enthalpy_with_geometry_tangent(
                    h,
                    bath,
                    face(h),
                    1.0 / h,
                    &mut || Ok(()),
                )
                .unwrap()
                .unwrap();
            let fixed = boundary
                .state_at_specific_enthalpy(h, bath, face(h), &mut || Ok(()))
                .unwrap();
            // Temperature and conductivity are both constant on this plateau.
            // All of this derivative comes from the original volume change.
            near(fixed.power_enthalpy_derivative_w_kg_j, 0.0, 0.0);
            assert!(state.power_enthalpy_derivative_w_kg_j.abs() > 1e-3);
            let power = |energy| {
                boundary
                    .evaluate(energy, bath, face(energy), &mut || Ok(()))
                    .unwrap()
                    .outward_power_w
            };
            near(
                state.power_enthalpy_derivative_w_kg_j,
                (power(h + 0.01) - power(h - 0.01)) / 0.02,
                2e-7,
            );
            if emissivity == 0.0 {
                let resistance = face(h).half_distance_m + 0.02 + 1.0 / 5.0;
                near(state.conductance_w_k, 0.3 / resistance, 1e-12);
                near(
                    state.power_enthalpy_derivative_w_kg_j,
                    -0.3 * (300.0 - bath) * (0.1 / 400.0) / resistance.powi(2),
                    1e-12,
                );
            }
            assert_eq!(
                fixed,
                boundary
                    .state_at_specific_enthalpy_with_geometry_tangent(
                        h,
                        bath,
                        face(h),
                        0.0,
                        &mut || Ok(()),
                    )
                    .unwrap()
                    .unwrap()
            );
        }
    }
}

struct Surface {
    queries: Cell<usize>,
}
impl ImplicitGraySurfaceLaw for Surface {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        &PLATEAU
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [100.0, 1000.0]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitGraySurfaceState, DepositionError> {
        self.queries.set(self.queries.get() + 1);
        Ok(ImplicitGraySurfaceState {
            emissivity: 0.4 + 0.0004 * (h - 400.0),
            emissivity_enthalpy_derivative_kg_j: 0.0004,
            convection_w_m2_k: 4.0 + 0.004 * (h - 400.0),
            convection_enthalpy_derivative_kg_w_m2_k_j: 0.004,
        })
    }
}

#[test]
fn geometry_and_original_surface_response_share_one_sample_and_full_tangent() {
    let surface = Surface {
        queries: Cell::new(0),
    };
    let common = GrayAmbientBoundary::new(&PLATEAU, 0.0, 0.0, GrayAmbientControls::default())
        .unwrap()
        .with_surface_law(&surface)
        .unwrap();
    let distinct =
        GrayTwoBathBoundary::new(&PLATEAU, 0.0, 0.0, 450.0, GrayAmbientControls::default())
            .unwrap()
            .with_surface_law(&surface)
            .unwrap();
    let h = 400.0;
    for two_baths in [false, true] {
        surface.queries.set(0);
        let derivative = if two_baths {
            distinct
                .state_at_specific_enthalpy_with_geometry_tangent(
                    h,
                    200.0,
                    face(h),
                    1.0 / h,
                    &mut || Ok(()),
                )
                .unwrap()
                .unwrap()
                .power_enthalpy_derivative_w_kg_j
        } else {
            common
                .state_at_specific_enthalpy_with_geometry_tangent(
                    h,
                    200.0,
                    face(h),
                    1.0 / h,
                    &mut || Ok(()),
                )
                .unwrap()
                .unwrap()
                .power_enthalpy_derivative_w_kg_j
        };
        assert_eq!(surface.queries.get(), 1);
        let power = |energy| {
            if two_baths {
                distinct
                    .evaluate(energy, 200.0, face(energy), &mut || Ok(()))
                    .unwrap()
                    .outward_power_w
            } else {
                common
                    .evaluate(energy, 200.0, face(energy), &mut || Ok(()))
                    .unwrap()
                    .outward_power_w
            }
        };
        near(derivative, (power(h + 0.01) - power(h - 0.01)) / 0.02, 2e-7);
    }
}

#[test]
fn two_bath_geometry_preserves_each_signed_exchange_and_actual_bath_entropy() {
    let boundary =
        GrayTwoBathBoundary::new(&PLATEAU, 0.8, 6.0, 450.0, GrayAmbientControls::default())
            .unwrap();
    let h = 400.0;
    let state = boundary
        .state_at_specific_enthalpy_with_geometry_tangent(
            h,
            250.0,
            face(h),
            1.0 / h,
            &mut || Ok(()),
        )
        .unwrap()
        .unwrap();
    assert!(state.baths[0].outward_power_w > 0.0);
    assert!(state.baths[1].outward_power_w < 0.0);
    assert_eq!(state.baths.map(|bath| bath.temperature_k), [250.0, 450.0]);
    let power = |energy| {
        boundary
            .evaluate(energy, 250.0, face(energy), &mut || Ok(()))
            .unwrap()
            .outward_power_w
    };
    near(
        state.power_enthalpy_derivative_w_kg_j,
        (power(h + 0.01) - power(h - 0.01)) / 0.02,
        2e-7,
    );

    let (initial, observer) = setup(&EXPANDING, &[400.0], Some((300.0, None)));
    let observer = layers(&initial, &observer);
    let exchange =
        GrayTwoBathBoundary::new(&EXPANDING, 0.8, 3.0, 600.0, GrayAmbientControls::default())
            .unwrap();
    let graph = ImplicitConductionGraph::new(
        &initial,
        &observer,
        &[(MaterialId(0), &EXPANDING)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_boundary_exchanges(&initial, &[(CellId(0), &exchange)], &[], || Ok(()))
    .unwrap();
    let dt = 0.25;
    let (next, report) = graph
        .advance_to(
            &initial,
            dt,
            controls(dt),
            &mut ImplicitConductionWorkspace::new(1).unwrap(),
            || Ok(()),
        )
        .unwrap();
    let final_h = next.cells()[&CellId(0)].specific_h_j_kg;
    let final_exchange = exchange
        .evaluate(
            final_h,
            300.0,
            ImplicitBoundaryGeometry {
                area_m2: 1.0,
                half_distance_m: final_h / 600.0,
                surface_resistance_m2_k_w: 0.0,
            },
            &mut || Ok(()),
        )
        .unwrap();
    near(final_h - 400.0, -dt * final_exchange.outward_power_w, 1e-9);
    near(
        report.energy.external_heat_j,
        -dt * final_exchange.outward_power_w,
        1e-9,
    );
    let entropy_input = -dt
        * final_exchange
            .boundary
            .baths
            .iter()
            .map(|bath| bath.outward_power_w / bath.temperature_k)
            .sum::<f64>();
    near(
        report.entropy.reservoir_entropy_input_j_k,
        entropy_input,
        1e-10,
    );
    near(report.energy.residual_j, 0.0, 1e-9);
    assert_eq!(
        next.cells()[&CellId(0)].mass_kg,
        initial.cells()[&CellId(0)].mass_kg
    );
    assert_eq!(
        next.cells()[&CellId(0)].geometry,
        initial.cells()[&CellId(0)].geometry
    );
    assert_eq!(report.energy.substeps, 1);
}

struct Gap<'a> {
    materials: [&'a dyn ImplicitConductionLaw; 2],
    queries: Cell<usize>,
}
impl ImplicitContactLaw for Gap<'_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        self.materials.map(|material| material.caloric_law())
    }
    fn enthalpy_bounds_j_kg(&self) -> [[f64; 2]; 2] {
        [[100.0, 1000.0]; 2]
    }
    fn state_at_specific_enthalpies(
        &self,
        h: [f64; 2],
    ) -> Result<ImplicitContactState, DepositionError> {
        self.queries.set(self.queries.get() + 1);
        Ok(ImplicitContactState {
            resistance_m2_k_w: 0.25 + 1e-4 * (h[0] - 400.0) + 2e-4 * (h[1] - 500.0),
            resistance_enthalpy_derivatives_m2_kg_k_w_j: [1e-4, 2e-4],
        })
    }
}

#[test]
fn gray_gap_differentiates_both_density_layers_and_samples_original_contact_once() {
    let sensible = Expanding {
        nonlinear_k: true,
        ..EXPANDING
    };
    let gap = Gap {
        materials: [&PLATEAU, &sensible],
        queries: Cell::new(0),
    };
    let gray = GrayRadiativeLink::new(
        gap.materials,
        [0.6, 0.8],
        0.0,
        [0.02, 0.03],
        GrayRadiativeControls::default(),
    )
    .unwrap()
    .with_gap_resistance(&gap)
    .unwrap();
    let patch = ContactPatch::symmetric(ThermalLink {
        id: 7,
        a: CellId(0),
        b: CellId(1),
        area_m2: 0.5,
        conductance_w_k: 0.0,
    });
    let geometry = |h: [f64; 2]| ConductiveLink {
        link_id: 7,
        half_distances_m: [0.1 * h[0] / 400.0, 0.2 * h[1] / 500.0],
        contact_resistance_m2_k_w: 0.0,
    };
    let h = [400.0, 500.0];
    let state = gray
        .state_at_specific_enthalpies_with_geometry_tangent(
            h,
            geometry(h),
            patch,
            [1.0 / h[0], 1.0 / h[1]],
            &mut || Ok(()),
        )
        .unwrap()
        .unwrap();
    assert_eq!(
        gap.queries.get(),
        1,
        "geometry and source tangents share one original face solve"
    );
    for side in 0..2 {
        let mut plus = h;
        let mut minus = h;
        plus[side] += 0.01;
        minus[side] -= 0.01;
        let high = gray
            .evaluate(plus, geometry(plus), patch, &mut || Ok(()))
            .unwrap()
            .power_w;
        let low = gray
            .evaluate(minus, geometry(minus), patch, &mut || Ok(()))
            .unwrap()
            .power_w;
        near(
            state.power_enthalpy_derivatives_w_kg_j[side],
            (high - low) / 0.02,
            2e-7,
        );
    }
}

#[test]
fn both_density_boundary_channels_match_the_implicit_energy_oracle_and_cancellation() {
    let (initial, observer) = setup(&EXPANDING, &[300.0], Some((600.0, Some(600.0))));
    let observer = layers(&initial, &observer);
    // Reservoir film R=1/2; ambient coating+film R=0.4+1/10.
    // Both share d=T/600 and the same physical bath, on distinct unit faces.
    let reservoir =
        GrayAmbientBoundary::new(&EXPANDING, 0.0, 2.0, GrayAmbientControls::default()).unwrap();
    let ambient =
        GrayAmbientBoundary::new(&EXPANDING, 0.0, 10.0, GrayAmbientControls::default()).unwrap();
    let chain =
        ImplicitConductionChain::new(&initial, &observer, &[(MaterialId(0), &EXPANDING)], || {
            Ok(())
        })
        .unwrap()
        .with_boundary_laws(
            &initial,
            &[(CellId(0), &reservoir)],
            &[(CellId(0), &ambient)],
            || Ok(()),
        )
        .unwrap();
    let (next, report) = chain
        .advance_to(
            &initial,
            1.0,
            controls(1.0),
            &mut ImplicitConductionWorkspace::new(1).unwrap(),
            || Ok(()),
        )
        .unwrap();
    // T-300 = 1200*(600-T)/(T+300), with the conserved unit mass.
    let expected = 0.5 * (4_680_000.0_f64.sqrt() - 1200.0);
    near(next.temperature_k(CellId(0)).unwrap(), expected, 1e-9);
    near(report.energy.external_heat_j, expected - 300.0, 1e-9);
    near(report.energy.residual_j, 0.0, 1e-9);
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(next.cells()[&CellId(0)].mass_kg, 1.0);
    assert_eq!(
        next.cells()[&CellId(0)].geometry,
        initial.cells()[&CellId(0)].geometry
    );

    let mut calls = 0;
    let cancelled = reservoir.state_at_specific_enthalpy_with_geometry_tangent(
        400.0,
        600.0,
        face(400.0),
        1.0 / 400.0,
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
    assert_eq!(initial.temperature_k(CellId(0)).unwrap(), 300.0);
}

#[test]
fn density_gray_link_keeps_closed_inventory_and_energy_through_graph_sdirk() {
    let (initial, observer) = setup(&EXPANDING, &[500.0, 200.0], None);
    let observer = layers(&initial, &observer);
    let gray = GrayRadiativeLink::new(
        [&EXPANDING, &EXPANDING],
        [0.7, 0.9],
        2.0,
        [0.02, 0.03],
        GrayRadiativeControls::default(),
    )
    .unwrap();
    let graph = ImplicitConductionGraph::new(
        &initial,
        &observer,
        &[(MaterialId(0), &EXPANDING)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_link_laws(&initial, &[(0, &gray)], || Ok(()))
    .unwrap();
    let (next, report) = graph
        .advance_sdirk2_to(
            &initial,
            0.5,
            controls(0.5),
            &mut ImplicitConductionWorkspace::new(2).unwrap(),
            || Ok(()),
        )
        .unwrap();
    assert!(next.temperature_k(CellId(0)).unwrap() < 500.0);
    assert!(next.temperature_k(CellId(1)).unwrap() > 200.0);
    near(next.energy_j(), initial.energy_j(), 1e-9);
    near(report.energy.external_heat_j, 0.0, 0.0);
    near(report.energy.residual_j, 0.0, 1e-9);
    assert!(report.entropy.conduction_generation_j_k > 0.0);
    for (&id, cell) in next.cells() {
        assert_eq!(cell.mass_kg, initial.cells()[&id].mass_kg);
        assert_eq!(cell.geometry, initial.cells()[&id].geometry);
    }
}
