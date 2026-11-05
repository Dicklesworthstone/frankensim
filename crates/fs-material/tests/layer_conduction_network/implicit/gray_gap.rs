//! Original source profiles, independent surface balances and native energy.
use super::*;
use fs_material::thermal::conduction::{ConductivityCurve, GrayLayerGap, TwoLayerConduction};
use fs_thermal::enthalpy::conduction::implicit::link::{GrayRadiativeControls, GrayRadiativeLink};
use fs_thermal::enthalpy::conduction::implicit::{ImplicitConductionLaw, ImplicitConductionState};

const SIGMA: f64 = 5.670_374_419e-8;
const EMISSIVITIES: [f64; 2] = [0.8, 0.6];

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}

fn radiative_flux(a: f64, b: f64) -> f64 {
    SIGMA * (a.powi(4) - b.powi(4)) / (1.0 / EMISSIVITIES[0] + 1.0 / EMISSIVITIES[1] - 1.0)
}

#[test]
fn nonlinear_profiles_and_gray_surfaces_match_three_independent_fluxes() {
    // Independent manufactured state: hot k=2+.01T, cold k=8-.005T.
    // Their own nonoverlapping source spans cover only their actual profiles.
    let hot = ConductivityCurve::new(vec![[580.0, 7.8], [900.0, 11.0]]).unwrap();
    let cold = ConductivityCurve::new(vec![[250.0, 6.75], [420.0, 5.9]]).unwrap();
    let radiation = radiative_flux(600.0, 400.0);
    let flux = 20.0 * 200.0 + radiation;
    let distances = [850.0 / flux, 625.0 / flux];
    let model = TwoLayerConduction::new([&hot, &cold], distances, 0.004, 0.0).unwrap();
    let gap = GrayLayerGap::equal_area_full_view(EMISSIVITIES, 20.0).unwrap();
    let result = model
        .evaluate_gray_gap([700.0, 300.0], gap, config(), || Ok(()))
        .unwrap();
    near(result.power_a_to_b_w, 0.004 * flux, 1e-9);
    near(result.interface_temperature_k[0], 600.0, 1e-9);
    near(result.interface_temperature_k[1], 400.0, 1e-9);
    near(result.gap_conduction_power_w, 16.0, 1e-9);
    near(result.radiation_power_w, 0.004 * radiation, 1e-9);
    near(result.layer_heat_residual_w, 0.0, 1e-9);
    near(result.surface_heat_residual_w, 0.0, 1e-9);
    let [u, v] = result.interface_temperature_k;
    // Reconstruct using closed-form primitives, independently of chart inversion.
    near(
        0.004 * (2.0 * (700.0 - u) + 0.005 * (700.0_f64.powi(2) - u * u)) / distances[0],
        result.power_a_to_b_w,
        1e-9,
    );
    near(
        0.004 * (8.0 * (v - 300.0) - 0.0025 * (v * v - 300.0_f64.powi(2))) / distances[1],
        result.power_a_to_b_w,
        1e-9,
    );
    let swapped = TwoLayerConduction::new([&cold, &hot], [distances[1], distances[0]], 0.004, 0.0)
        .unwrap()
        .evaluate_gray_gap([300.0, 700.0], gap, config(), || Ok(()))
        .unwrap();
    assert_eq!(result.power_a_to_b_w, -swapped.power_a_to_b_w);
    assert_eq!(
        result.interface_temperature_k,
        [
            swapped.interface_temperature_k[1],
            swapped.interface_temperature_k[0]
        ]
    );
    for side in 0..2 {
        let mut plus = [700.0, 300.0];
        let mut minus = plus;
        plus[side] += 0.001;
        minus[side] -= 0.001;
        let derivative = (model
            .evaluate_gray_gap(plus, gap, config(), || Ok(()))
            .unwrap()
            .power_a_to_b_w
            - model
                .evaluate_gray_gap(minus, gap, config(), || Ok(()))
                .unwrap()
                .power_a_to_b_w)
            / 0.002;
        near(result.endpoint_derivative_w_k[side], derivative, 1e-7);
    }
}

struct ConstantMaterial(f64);
impl ThermalLaw for ConstantMaterial {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        Ok(h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
impl ConductivityLaw for ConstantMaterial {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(self.0)
    }
}
impl ImplicitConductionLaw for ConstantMaterial {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [200.0, 1000.0]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: h,
            temperature_slope_kg_k_j: 1.0,
            conductivity_w_m_k: self.0,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
}

#[test]
fn constant_profiles_recover_the_native_gray_link_and_conduction_only_limits() {
    let curves = [8.0, 2.0].map(|k| ConductivityCurve::new(vec![[200.0, k], [1000.0, k]]).unwrap());
    let materials = [ConstantMaterial(8.0), ConstantMaterial(2.0)];
    let geometry = ConductiveLink {
        link_id: 7,
        half_distances_m: [0.02, 0.03],
        contact_resistance_m2_k_w: 0.0,
    };
    let patch = ContactPatch {
        link: ThermalLink {
            id: 7,
            a: CellId(1),
            b: CellId(2),
            area_m2: 0.004,
            conductance_w_k: 0.0,
        },
        covered_area_m2: [0.004; 2],
    };
    let native = GrayRadiativeLink::new(
        [&materials[0], &materials[1]],
        EMISSIVITIES,
        20.0,
        [0.0; 2],
        GrayRadiativeControls::default(),
    )
    .unwrap();
    let model =
        TwoLayerConduction::new([&curves[0], &curves[1]], [0.02, 0.03], 0.004, 0.0).unwrap();
    let gap = GrayLayerGap::equal_area_full_view(EMISSIVITIES, 20.0).unwrap();
    for t in [[600.0, 300.0], [300.0, 600.0], [400.0; 2]] {
        let result = model
            .evaluate_gray_gap(t, gap, config(), || Ok(()))
            .unwrap();
        let expected = native.evaluate(t, geometry, patch, &mut || Ok(())).unwrap();
        near(result.power_a_to_b_w, expected.power_w, 1e-8);
        for side in 0..2 {
            near(
                result.interface_temperature_k[side],
                expected.surface_temperatures_k.unwrap()[side],
                1e-8,
            );
            near(
                result.endpoint_derivative_w_k[side],
                expected.link.power_enthalpy_derivatives_w_kg_j[side],
                1e-9,
            );
        }
        let conductive = model
            .evaluate_gray_gap(
                t,
                GrayLayerGap::equal_area_full_view([0.0, 0.9], 20.0).unwrap(),
                config(),
                || Ok(()),
            )
            .unwrap();
        let resistance =
            TwoLayerConduction::new([&curves[0], &curves[1]], [0.02, 0.03], 0.004, 1.0 / 20.0)
                .unwrap()
                .evaluate(t, config(), || Ok(()))
                .unwrap();
        assert_eq!(conductive.power_a_to_b_w, resistance.power_a_to_b_w);
        assert_eq!(conductive.conductance_w_k, resistance.conductance_w_k);
        assert_eq!(
            conductive.interface_temperature_k,
            resistance.interface_temperature_k
        );
        assert_eq!(
            conductive.endpoint_derivative_w_k,
            resistance.endpoint_derivative_w_k
        );
        assert_eq!(conductive.radiation_power_w, 0.0);
        let sealed = model
            .evaluate_gray_gap(
                t,
                GrayLayerGap::equal_area_full_view([0.0; 2], 0.0).unwrap(),
                config(),
                || Ok(()),
            )
            .unwrap();
        assert_eq!(sealed.power_a_to_b_w, 0.0);
        assert_eq!(sealed.endpoint_derivative_w_k, [0.0; 2]);
    }
}

#[test]
fn sourced_gray_profiles_drive_the_native_chain_and_graph_once_with_full_tangents() {
    let card = source();
    let caloric = curve(&card);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let layer = layer(&card, &caloric, &context);
    let old = initial(&caloric);
    let law = transport(&card, &caloric, &context);
    let geometry = ConductiveLink {
        link_id: 7,
        half_distances_m: [0.02, 0.01665],
        contact_resistance_m2_k_w: 0.0,
    };
    // Independent target: bulk [590,310], surfaces [510,400] K. Both original
    // affine k(T) integrals give j=30000 W/m². Set only the gap HTC to match it.
    let htc = (30_000.0 - radiative_flux(510.0, 400.0)) / 110.0;
    let gap = GrayLayerGap::equal_area_full_view(EMISSIVITIES, htc).unwrap();
    let prepared = PreparedLayerConduction::new(
        &old,
        vec![LayerLinkBinding {
            distances_m: geometry.half_distances_m,
            ..binding(&layer)
        }],
        1,
        || Ok(()),
    )
    .unwrap();
    let link = prepared
        .implicit_link_laws(config())
        .next()
        .unwrap()
        .with_gray_gap(gap)
        .unwrap();
    let patch = old.contact_patch(7).unwrap();
    let at = |h: [f64; 2]| {
        let state = link
            .state_at_specific_enthalpies(h, geometry, patch, &mut || Ok(()))
            .unwrap();
        let t = h.map(|v| caloric.temperature_at_specific_enthalpy(v).unwrap());
        (state, state.conductance_w_k * (t[0] - t[1]))
    };
    for h in [[39_000.0, 11_000.0], [11_000.0, 39_000.0], [25_000.0; 2]] {
        let (state, _) = at(h);
        for side in 0..2 {
            let mut plus = h;
            let mut minus = h;
            plus[side] += 0.25;
            minus[side] -= 0.25;
            near(
                state.power_enthalpy_derivatives_w_kg_j[side],
                (at(plus).1 - at(minus).1) / 0.5,
                1e-8,
            );
        }
    }
    let observer = native_observer(&old, &law, &prepared);
    let chain = ImplicitConductionChain::new(&old, &observer, &[(MaterialId(1), &law)], || Ok(()))
        .unwrap()
        .with_link_laws(&old, &[(7, &link)], || Ok(()))
        .unwrap();
    let graph = ImplicitConductionGraph::new(
        &old,
        &observer,
        &[(MaterialId(1), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_link_laws(&old, &[(7, &link)], || Ok(()))
    .unwrap();
    let dt = 1.0 / 12.0; // 120 W * dt = 10 J; each original cell has C=1 J/K.
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let a = chain
        .advance_to(&old, dt, controls(dt), &mut workspace, || Ok(()))
        .unwrap();
    let b = graph
        .advance_to(&old, dt, controls(dt), &mut workspace, || Ok(()))
        .unwrap();
    for (next, report) in [a, b] {
        near(next.temperature_k(CellId(1)).unwrap(), 590.0, 1e-7);
        near(next.temperature_k(CellId(2)).unwrap(), 310.0, 1e-7);
        assert_eq!(report.energy.substeps, 1);
        assert_eq!(report.energy.external_heat_j, 0.0);
        near(report.energy.residual_j, 0.0, 1e-7);
        assert!(report.entropy.conduction_generation_j_k > 0.0);
        assert_eq!(next.links()[&7], old.links()[&7]);
    }
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.temperature_k(CellId(1)).unwrap(), 600.0);
}

#[test]
fn gray_profiles_enforce_source_support_interface_ownership_and_cancellation() {
    let card = source();
    let caloric = curve(&card);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let layer = layer(&card, &caloric, &context);
    let old = initial(&caloric);
    let gap = GrayLayerGap::equal_area_full_view(EMISSIVITIES, 20.0).unwrap();
    let prepared = PreparedLayerConduction::new(&old, vec![binding(&layer)], 1, || Ok(())).unwrap();
    let make = || prepared.implicit_link_laws(config()).next().unwrap();
    let contact = AffineLayerContact {
        caloric: [&caloric; 2],
        bounds: [[0.0, 80_000.0]; 2],
        reference_h: [35_000.0, 15_000.0],
        resistance: 0.02,
        derivatives: [0.0; 2],
    };
    assert!(
        make()
            .with_contact(&contact)
            .unwrap()
            .with_gray_gap(gap)
            .is_err()
    );
    assert!(
        make()
            .with_gray_gap(gap)
            .unwrap()
            .with_contact(&contact)
            .is_err()
    );
    assert!(
        make()
            .with_gray_gap(gap)
            .unwrap()
            .with_gray_gap(gap)
            .is_err()
    );
    let fixed = PreparedLayerConduction::new(
        &old,
        vec![LayerLinkBinding {
            contact_area_resistance_m2_k_w: 0.1,
            ..binding(&layer)
        }],
        1,
        || Ok(()),
    )
    .unwrap();
    assert!(
        fixed
            .implicit_link_laws(config())
            .next()
            .unwrap()
            .with_gray_gap(gap)
            .is_err()
    );
    let link = make().with_gray_gap(gap).unwrap();
    let geometry = ConductiveLink {
        link_id: 7,
        half_distances_m: [0.02, 0.03],
        contact_resistance_m2_k_w: 0.0,
    };
    assert!(
        link.state_at_specific_enthalpies(
            [90_000.0, 10_000.0],
            geometry,
            old.contact_patch(7).unwrap(),
            &mut || Ok(())
        )
        .is_err()
    );
    assert!(matches!(
        link.state_at_specific_enthalpies(
            [40_000.0, 10_000.0],
            geometry,
            old.contact_patch(7).unwrap(),
            &mut || Err(DepositionError::Cancelled)
        ),
        Err(DepositionError::Cancelled)
    ));
    let narrow = ConductivityCurve::new(vec![[590.0, 7.9], [650.0, 8.5]]).unwrap();
    let model = TwoLayerConduction::new(
        [&narrow, layer.conductivity_curve()],
        [0.02, 0.03],
        0.004,
        0.0,
    )
    .unwrap();
    assert!(matches!(
        model.evaluate_gray_gap([600.0, 300.0], gap, config(), || Ok(())),
        Err(DepositionError::MaterialDomain)
    ));
    let model =
        TwoLayerConduction::new([layer.conductivity_curve(); 2], [0.02, 0.03], 0.004, 0.0).unwrap();
    assert!(matches!(
        model.evaluate_gray_gap(
            [600.0, 300.0],
            gap,
            LayerConductionConfig {
                max_iterations: 1,
                ..config()
            },
            || Ok(())
        ),
        Err(DepositionError::Budget(_))
    ));
    for (emissivities, htc) in [([-0.1, 0.5], 1.0), ([0.5, 1.1], 1.0), ([0.5; 2], -1.0)] {
        assert!(GrayLayerGap::equal_area_full_view(emissivities, htc).is_err());
    }
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.temperature_k(CellId(1)).unwrap(), 600.0);
}
