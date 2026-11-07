use super::*;
#[path = "boundary.rs"]
mod boundary;
#[path = "implicit/gray_gap.rs"]
mod gray_gap;
use fs_material::thermal::conduction::phase::PreparedPhaseConductivity;
use fs_thermal::enthalpy::conduction::{ConductiveLink, ConductivityLaw, ConductionObserver};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionGraph,
    ImplicitConductionLinearControls, ImplicitConductionWorkspace,
    link::ImplicitLinkLaw,
};

fn transport<'a>(card: &'a MaterialCard, caloric: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint) -> PreparedPhaseConductivity<'a> {
    PreparedPhaseConductivity::new(caloric, ThermalMaterialPoint {
        material: card, context, temperature_axis: "T", liquid_fraction_axis: None,
    }, &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly), [0.0, 80_000.0], || Ok(())).unwrap()
}

fn native_observer<'a>(thermal: &DepositionSimulation<'_>, law: &'a dyn ConductivityLaw,
    prepared: &PreparedLayerConduction<'_, '_>) -> ConductionObserver<'a> {
    let geometry: Vec<_> = prepared.bindings().map(|b| ConductiveLink {
        link_id: b.link, half_distances_m: b.distances_m,
        contact_resistance_m2_k_w: b.contact_area_resistance_m2_k_w,
    }).collect();
    ConductionObserver::new(thermal, &[(MaterialId(1), law)], &geometry, &[], || Ok(())).unwrap()
}

fn controls(dt: f64) -> ImplicitConductionControls {
    ImplicitConductionControls { maximum_dt_s: dt, cell_tolerance_j: 1e-8,
        total_tolerance_j: 1e-7, ..ImplicitConductionControls::default() }
}

#[test]
fn kirchhoff_layers_drive_native_chain_and_graph_with_analytic_heat_and_tangents() {
    let card = source(); let caloric = curve(&card);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let layer = layer(&card, &caloric, &context);
    let law = transport(&card, &caloric, &context);
    let old = initial(&caloric);
    let prepared = PreparedLayerConduction::new(&old, vec![binding(&layer)], 1, || Ok(())).unwrap();
    let links: Vec<_> = prepared.implicit_link_laws(config()).collect();
    let bindings: Vec<_> = links.iter().map(|link| (link.link_id(), link as &dyn ImplicitLinkLaw)).collect();
    let observer = native_observer(&old, &law, &prepared);
    let geometry = ConductiveLink { link_id: 7, half_distances_m: [0.02, 0.03],
        contact_resistance_m2_k_w: 0.0 };
    let sample = links[0].state_at_specific_enthalpies([40_000.0, 10_000.0], geometry,
        old.contact_patch(7).unwrap(), &mut || Ok(())).unwrap();
    assert!((sample.conductance_w_k - 0.52).abs() < 1e-12);
    // Independent dQ/dh=(A/L)*k(T)*(dT/dh), including the receiving-side sign.
    assert!((sample.power_enthalpy_derivatives_w_kg_j[0] - 0.0064).abs() < 1e-12);
    assert!((sample.power_enthalpy_derivatives_w_kg_j[1] + 0.004).abs() < 1e-12);
    let chain = ImplicitConductionChain::new(&old, &observer, &[(MaterialId(1), &law)], || Ok(()))
        .unwrap().with_link_laws(&old, &bindings, || Ok(())).unwrap();
    let graph = ImplicitConductionGraph::new(&old, &observer, &[(MaterialId(1), &law)],
        ImplicitConductionLinearControls::default(), || Ok(()))
        .unwrap().with_link_laws(&old, &bindings, || Ok(())).unwrap();
    let dt = 0.25;
    let difference = 300.0 / (1.0 + 2.0 * 0.52 * dt);
    let mut work = ImplicitConductionWorkspace::new(2).unwrap();
    let (a, report) = chain.advance_to(&old, dt, controls(dt), &mut work, || Ok(())).unwrap();
    let (b, graph_report) = graph.advance_to(&old, dt, controls(dt), &mut work, || Ok(())).unwrap();
    for (id, expected) in [(CellId(1), 450.0 + 0.5 * difference), (CellId(2), 450.0 - 0.5 * difference)] {
        assert!((a.temperature_k(id).unwrap() - expected).abs() < 1e-8);
        assert!((b.temperature_k(id).unwrap() - expected).abs() < 1e-8);
    }
    for report in [report, graph_report] {
        assert_eq!(report.energy.substeps, 1);
        assert_eq!(report.energy.external_heat_j, 0.0);
        assert!(report.energy.residual_j.abs() < 1e-7);
        assert!(report.entropy.conduction_generation_j_k > 0.0);
    }
    let sampled = prepared.resolve(&a, config(), || Ok(())).unwrap();
    assert!(Arc::ptr_eq(&sampled.payload[0].source_lines[0], layer.source_line()));
    assert_eq!(a.links()[&7], old.links()[&7]);
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn disjoint_source_spans_and_contact_use_the_native_nonlinear_solve() {
    let card = source(); let caloric = curve(&card);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let prepare = |range| ResolvedConductionLayer::new(&caloric, ThermalMaterialPoint {
        material: &card, context: &context, temperature_axis: "T", liquid_fraction_axis: None,
    }, &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly), SolidLiquidPhase::Solid,
        range, || Ok(())).unwrap();
    let hot = prepare([500.0, 800.0]); let cold = prepare([250.0, 450.0]);
    let old = initial(&caloric);
    let law = transport(&card, &caloric, &context);
    let prepared = PreparedLayerConduction::new(&old, vec![LayerLinkBinding {
        link: 7, layers: [&hot, &cold], distances_m: [0.02, 0.01665],
        contact_area_resistance_m2_k_w: 110.0 / 30000.0,
    }], 1, || Ok(())).unwrap();
    let links: Vec<_> = prepared.implicit_link_laws(config()).collect();
    let observer = native_observer(&old, &law, &prepared);
    let chain = ImplicitConductionChain::new(&old, &observer, &[(MaterialId(1), &law)], || Ok(()))
        .unwrap().with_link_laws(&old, &[(7, &links[0])], || Ok(())).unwrap();
    let mut work = ImplicitConductionWorkspace::new(2).unwrap();
    let dt = 1.0 / 12.0;
    let (next, report) = chain.advance_to(&old, dt, controls(dt), &mut work, || Ok(())).unwrap();
    // Independent affine integrals: endpoint [590,310] K, surfaces [510,400] K,
    // 30000 W/m² * .004 m² * (1/12) s = 10 J transferred internally.
    assert!((next.temperature_k(CellId(1)).unwrap() - 590.0).abs() < 1e-7);
    assert!((next.temperature_k(CellId(2)).unwrap() - 310.0).abs() < 1e-7);
    let response = prepared.resolve(&next, config(), || Ok(())).unwrap().payload[0].response;
    assert!((response.interface_temperature_k[0] - 510.0).abs() < 1e-6);
    assert!((response.interface_temperature_k[1] - 400.0).abs() < 1e-6);
    assert_eq!(report.energy.external_heat_j, 0.0);
    assert!(report.energy.residual_j.abs() < 1e-7);
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn integrated_source_binding_refuses_geometry_and_domain_changes() {
    let card = source(); let caloric = curve(&card);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let layer = layer(&card, &caloric, &context);
    let old = initial(&caloric);
    let prepared = PreparedLayerConduction::new(&old, vec![binding(&layer)], 1, || Ok(())).unwrap();
    let links: Vec<_> = prepared.implicit_link_laws(config()).collect();
    let geometry = ConductiveLink { link_id: 7, half_distances_m: [0.02, 0.03],
        contact_resistance_m2_k_w: 0.0 };
    let patch = old.contact_patch(7).unwrap();
    assert!(links[0].state_at_specific_enthalpies([90_000.0, 10_000.0], geometry,
        patch, &mut || Ok(())).is_err());
    assert!(links[0].state_at_specific_enthalpies([40_000.0, 10_000.0],
        ConductiveLink { half_distances_m: [0.03, 0.02], ..geometry }, patch, &mut || Ok(())).is_err());
    assert!(links[0].state_at_specific_enthalpies([40_000.0, 10_000.0], geometry,
        ContactPatch { covered_area_m2: [0.005, 0.006], ..patch }, &mut || Ok(())).is_err());
    assert!(matches!(links[0].state_at_specific_enthalpies([40_000.0, 10_000.0], geometry,
        patch, &mut || Err(DepositionError::Cancelled)), Err(DepositionError::Cancelled)));
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn solid_source_endpoint_uses_its_admitted_left_caloric_derivative() {
    let card = source(); let context = QueryPoint::new().with("T", 300.0).unwrap();
    let caloric = EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(),
        [(0.0, 200.0, 0.0), (40_000.0, 600.0, 0.0), (60_000.0, 600.0, 1.0), (100_000.0, 1000.0, 1.0)]
        .into_iter().map(|(h, t, f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg: h,
            temperature_k: t, liquid_mass_fraction: f, bulk_density_kg_m3: 1000.0 }).collect()).unwrap();
    let layer = ResolvedConductionLayer::new(&caloric, ThermalMaterialPoint {
        material: &card, context: &context, temperature_axis: "T", liquid_fraction_axis: None,
    }, &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly), SolidLiquidPhase::Solid,
        [200.0, 600.0], || Ok(())).unwrap();
    let old = initial(&caloric);
    let prepared = PreparedLayerConduction::new(&old, vec![binding(&layer)], 1, || Ok(())).unwrap();
    let links: Vec<_> = prepared.implicit_link_laws(config()).collect();
    let geometry = ConductiveLink { link_id: 7, half_distances_m: [0.02, 0.03],
        contact_resistance_m2_k_w: 0.0 };
    let sample = links[0].state_at_specific_enthalpies([40_000.0, 10_000.0], geometry,
        old.contact_patch(7).unwrap(), &mut || Ok(())).unwrap();
    assert_eq!(caloric.temperature_derivative_at_specific_enthalpy(40_000.0).unwrap(), 0.0);
    assert!((sample.power_enthalpy_derivatives_w_kg_j[0] - 0.0064).abs() < 1e-12);
    assert!(links[0].state_at_specific_enthalpies([40_000.0_f64.next_up(), 10_000.0], geometry,
        old.contact_patch(7).unwrap(), &mut || Ok(())).is_err());
}

fn affine_capacity_curve(card: &MaterialCard, cp: [f64; 2]) -> EquilibriumEnthalpyPhaseCurve {
    use fs_material::phase::HeatCapacityKnot;
    EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity_exact(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        0.0,
        &[(200.0, cp[0]), (1000.0, cp[1])].map(|(temperature_k, specific_heat_capacity_j_kg_k)| {
            HeatCapacityKnot {
                temperature_k,
                specific_heat_capacity_j_kg_k,
                bulk_density_kg_m3: 1000.0,
            }
        }),
    )
    .unwrap()
}

#[test]
fn source_endpoint_link_and_bath_tangents_retain_original_affine_heat_capacity() {
    use fs_material::thermal::conduction::network::ImplicitLayerBoundary;
    use fs_thermal::enthalpy::conduction::implicit::boundary::{
        ImplicitBoundaryGeometry, ImplicitBoundaryLaw,
    };

    for cp in [[50.0, 150.0], [150.0, 50.0]] {
        let card = source();
        let caloric = affine_capacity_curve(&card, cp);
        let context = QueryPoint::new().with("T", 300.0).unwrap();
        let old = initial(&caloric);
        // Independent integral of Cp(T)=cp[0]+(cp[1]-cp[0])*(T-200)/800.
        let enthalpy = |t: f64| {
            let delta = t - 200.0;
            delta * (cp[0] + 0.5 * (cp[1] - cp[0]) * delta / 800.0)
        };
        // Exercise both a bound inside the original Cp interval and its exact
        // upper knot. The enthalpy chord has Cp=100 in both source directions.
        for upper in [900.0, 1000.0] {
            let layer = ResolvedConductionLayer::new(
                &caloric,
                ThermalMaterialPoint {
                    material: &card,
                    context: &context,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
                ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
                SolidLiquidPhase::Solid,
                [200.0, upper],
                || Ok(()),
            )
            .unwrap();
            let prepared =
                PreparedLayerConduction::new(&old, vec![binding(&layer)], 1, || Ok(())).unwrap();
            let link = prepared.implicit_link_laws(config()).next().unwrap();
            let geometry = ConductiveLink {
                link_id: 7,
                half_distances_m: [0.02, 0.03],
                contact_resistance_m2_k_w: 0.0,
            };
            let h = enthalpy(upper);
            assert_eq!(layer.temperature_at_enthalpy(h).unwrap(), upper);
            let sample = link
                .state_at_specific_enthalpies(
                    [h, enthalpy(300.0)],
                    geometry,
                    old.contact_patch(7).unwrap(),
                    &mut || Ok(()),
                )
                .unwrap();
            let capacity = cp[0] + (cp[1] - cp[0]) * (upper - 200.0) / 800.0;
            // For the same material on both sides, Q=A/L*integral(k dT).
            // Its exact hot-side derivative is A/L*k(T)/Cp(T), not k/Cp_mean.
            let expected = 0.004 / 0.05 * (2.0 + 0.01 * upper) / capacity;
            assert!((sample.power_enthalpy_derivatives_w_kg_j[0] - expected).abs() < 1e-11);

            let face = ImplicitBoundaryGeometry {
                area_m2: 0.004,
                half_distance_m: 0.02,
                surface_resistance_m2_k_w: 0.0,
            };
            let boundary = ImplicitLayerBoundary::new(&layer, face, config()).unwrap();
            let sample = boundary
                .state_at_specific_enthalpy(h, 300.0, face, &mut || Ok(()))
                .unwrap();
            let expected = 0.004 / 0.02 * (2.0 + 0.01 * upper) / capacity;
            assert!((sample.power_enthalpy_derivative_w_kg_j - expected).abs() < 1e-11);
        }
    }
}

#[derive(Clone, Copy)]
struct AffineLayerContact<'a> {
    caloric: [&'a EquilibriumEnthalpyPhaseCurve; 2],
    bounds: [[f64; 2]; 2],
    reference_h: [f64; 2],
    resistance: f64,
    derivatives: [f64; 2],
}

impl fs_thermal::enthalpy::conduction::implicit::contact::ImplicitContactLaw
    for AffineLayerContact<'_>
{
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        self.caloric.map(|law| law as &dyn ThermalLaw)
    }

    fn enthalpy_bounds_j_kg(&self) -> [[f64; 2]; 2] {
        self.bounds
    }

    fn state_at_specific_enthalpies(
        &self,
        h: [f64; 2],
    ) -> Result<
        fs_thermal::enthalpy::conduction::implicit::contact::ImplicitContactState,
        DepositionError,
    > {
        use fs_thermal::enthalpy::conduction::implicit::contact::ImplicitContactState;
        Ok(ImplicitContactState {
            resistance_m2_k_w: self.derivatives[0].mul_add(
                h[0] - self.reference_h[0],
                self.derivatives[1].mul_add(h[1] - self.reference_h[1], self.resistance),
            ),
            resistance_enthalpy_derivatives_m2_kg_k_w_j: self.derivatives,
        })
    }
}

#[test]
fn contact_bounds_inside_affine_heat_capacity_spans_keep_exact_bulk_tangents() {
    for cp in [[50.0, 150.0], [150.0, 50.0]] {
        let card = source();
        let caloric = affine_capacity_curve(&card, cp);
        let context = QueryPoint::new().with("T", 300.0).unwrap();
        let layer = layer(&card, &caloric, &context);
        let old = initial(&caloric);
        let prepared =
            PreparedLayerConduction::new(&old, vec![binding(&layer)], 1, || Ok(())).unwrap();
        let upper = caloric
            .specific_enthalpy_on_sensible_interval(0, 700.0)
            .unwrap();
        let cold = caloric
            .specific_enthalpy_on_sensible_interval(0, 300.0)
            .unwrap();
        let contact = AffineLayerContact {
            caloric: [&caloric; 2],
            bounds: [[0.0, upper]; 2],
            reference_h: [upper, cold],
            resistance: 0.02,
            derivatives: [2e-7, -1e-7],
        };
        let link = prepared
            .implicit_link_laws(config())
            .next()
            .unwrap()
            .with_contact(&contact)
            .unwrap();
        let geometry = ConductiveLink {
            link_id: 7,
            half_distances_m: [0.02, 0.03],
            contact_resistance_m2_k_w: 0.0,
        };
        // A contact can truncate either endpoint inside the same original
        // sensible interval. Check the full source/contact flux derivative
        // from within that bound, including both signs of heat transfer.
        for side in 0..2 {
            let at = |h| {
                let mut enthalpies = [cold; 2];
                enthalpies[side] = h;
                let sample = link
                    .state_at_specific_enthalpies(
                        enthalpies,
                        geometry,
                        old.contact_patch(7).unwrap(),
                        &mut || Ok(()),
                    )
                    .unwrap();
                let temperatures =
                    enthalpies.map(|h| caloric.temperature_at_specific_enthalpy(h).unwrap());
                (
                    sample,
                    sample.conductance_w_k * (temperatures[0] - temperatures[1]),
                )
            };
            let (sample, power) = at(upper);
            let inward = (3.0 * power - 4.0 * at(upper - 1.0).1 + at(upper - 2.0).1) / 2.0;
            assert!((sample.power_enthalpy_derivatives_w_kg_j[side] - inward).abs() < 1e-8);
        }
    }
}

#[test]
fn integrated_bulk_and_contact_tangents_differentiate_the_actual_combined_heat() {
    let card = source();
    let caloric = curve(&card);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let layer = layer(&card, &caloric, &context);
    let old = initial(&caloric);
    let prepared = PreparedLayerConduction::new(&old, vec![binding(&layer)], 1, || Ok(())).unwrap();
    let contact = AffineLayerContact {
        caloric: [&caloric; 2],
        bounds: [[0.0, 80_000.0]; 2],
        reference_h: [35_000.0, 15_000.0],
        resistance: 0.02,
        derivatives: [2e-7, -1e-7],
    };
    let link = prepared
        .implicit_link_laws(config())
        .next()
        .unwrap()
        .with_contact(&contact)
        .unwrap();
    let geometry = ConductiveLink {
        link_id: 7,
        half_distances_m: [0.02, 0.03],
        contact_resistance_m2_k_w: 0.0,
    };
    let patch = old.contact_patch(7).unwrap();
    let at = |h: [f64; 2]| {
        let state = link
            .state_at_specific_enthalpies(h, geometry, patch, &mut || Ok(()))
            .unwrap();
        let temperatures = h.map(|value| caloric.temperature_at_specific_enthalpy(value).unwrap());
        let power = state.conductance_w_k * (temperatures[0] - temperatures[1]);
        (state, power)
    };
    // These finite differences vary BOTH the source contact and integrated
    // bulk conductivity, including opposite heat flow and thermal equilibrium.
    for h in [[35_000.0, 15_000.0], [15_000.0, 35_000.0], [25_000.0; 2]] {
        let (state, _) = at(h);
        for side in 0..2 {
            let mut plus = h;
            let mut minus = h;
            plus[side] += 0.25;
            minus[side] -= 0.25;
            let difference = (at(plus).1 - at(minus).1) / 0.5;
            assert!(
                (difference - state.power_enthalpy_derivatives_w_kg_j[side]).abs() < 1e-8,
                "h={h:?}, side={side}: finite difference {difference}, tangent {:?}",
                state.power_enthalpy_derivatives_w_kg_j,
            );
        }
    }
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn integrated_contact_domains_and_ownership_preserve_transactional_heat_history() {
    let card = source();
    let caloric = curve(&card);
    let foreign = caloric.clone();
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let layer = layer(&card, &caloric, &context);
    let old = initial(&caloric);
    let law = transport(&card, &caloric, &context);
    let prepared = PreparedLayerConduction::new(&old, vec![binding(&layer)], 1, || Ok(())).unwrap();
    let make = || prepared.implicit_link_laws(config()).next().unwrap();
    let contact = AffineLayerContact {
        caloric: [&caloric; 2],
        bounds: [[0.0, 80_000.0]; 2],
        reference_h: [40_000.0, 10_000.0],
        resistance: 0.02,
        derivatives: [2e-7, -1e-7],
    };
    let other = AffineLayerContact {
        caloric: [&caloric, &foreign],
        ..contact
    };
    assert!(make().with_contact(&other).is_err());
    let invalid = AffineLayerContact {
        bounds: [[f64::NAN, 80_000.0]; 2],
        ..contact
    };
    assert!(make().with_contact(&invalid).is_err());
    assert!(
        make()
            .with_contact(&contact)
            .unwrap()
            .with_contact(&contact)
            .is_err()
    );
    let mut fixed_binding = binding(&layer);
    fixed_binding.contact_area_resistance_m2_k_w = 0.01;
    let fixed = PreparedLayerConduction::new(&old, vec![fixed_binding], 1, || Ok(())).unwrap();
    assert!(
        fixed
            .implicit_link_laws(config())
            .next()
            .unwrap()
            .with_contact(&contact)
            .is_err()
    );

    let geometry = ConductiveLink {
        link_id: 7,
        half_distances_m: [0.02, 0.03],
        contact_resistance_m2_k_w: 0.0,
    };
    let patch = old.contact_patch(7).unwrap();
    let limited = AffineLayerContact {
        bounds: [[0.0, 30_000.0]; 2],
        ..contact
    };
    let limited_link = make().with_contact(&limited).unwrap();
    assert!(matches!(
        limited_link.state_at_specific_enthalpies(
            [40_000.0, 10_000.0],
            geometry,
            patch,
            &mut || Ok(()),
        ),
        Err(DepositionError::MaterialDomain),
    ));
    // Preserve the original contact admission, including the zero-resistance
    // boundary: a contact cannot report a negative value or a nonzero tangent
    // at a zero value. Neither is clipped into an admissible heat coefficient.
    for resistance in [-1.0, 0.0] {
        let invalid = AffineLayerContact {
            resistance,
            ..contact
        };
        let link = make().with_contact(&invalid).unwrap();
        assert!(matches!(
            link.state_at_specific_enthalpies(
                [40_000.0, 10_000.0],
                geometry,
                patch,
                &mut || Ok(()),
            ),
            Err(DepositionError::MaterialDomain),
        ));
    }

    let link = make().with_contact(&contact).unwrap();
    let observer = native_observer(&old, &law, &prepared);
    let chain = ImplicitConductionChain::new(&old, &observer, &[(MaterialId(1), &law)], || Ok(()))
        .unwrap()
        .with_link_laws(&old, &[(7, &link)], || Ok(()))
        .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let mut polls = 0;
    let (expected, report) = chain
        .advance_to(&old, 0.1, controls(0.1), &mut workspace, || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    assert!(expected.temperature_k(CellId(1)).unwrap() < 600.0);
    assert_eq!(report.energy.external_heat_j, 0.0);
    assert!(report.energy.residual_j.abs() < 1e-7);
    let mut calls = 0;
    assert!(matches!(
        chain.advance_to(&old, 0.1, controls(0.1), &mut workspace, || {
            calls += 1;
            if calls == polls - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(
            fs_thermal::enthalpy::conduction::implicit::ImplicitConductionError::Thermal(
                DepositionError::Cancelled,
            )
        ),
    ));
    let (retry, _) = chain
        .advance_to(&old, 0.1, controls(0.1), &mut workspace, || Ok(()))
        .unwrap();
    assert_eq!(retry.cells(), expected.cells());
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.temperature_k(CellId(1)).unwrap(), 600.0);
}

#[test]
fn integrated_contacts_share_the_bulk_branch_and_its_admitted_endpoint_tangent() {
    use fs_thermal::enthalpy::conduction::implicit::contact::{
        ImplicitContactLaw, ImplicitContactState,
    };
    struct TemperatureContact<'a>(&'a EquilibriumEnthalpyPhaseCurve, f64);
    impl ImplicitContactLaw for TemperatureContact<'_> {
        fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
            [self.0; 2]
        }
        fn enthalpy_bounds_j_kg(&self) -> [[f64; 2]; 2] {
            [[0.0, self.1]; 2]
        }
        fn state_at_specific_enthalpies(
            &self,
            h: [f64; 2],
        ) -> Result<ImplicitContactState, DepositionError> {
            let t = self
                .0
                .temperature_at_specific_enthalpy(h[0])
                .map_err(|_| DepositionError::MaterialDomain)?;
            // Same endpoint convention as the original prepared source owner.
            let at = if h[0] == self.1 {
                h[0].next_down()
            } else {
                h[0]
            };
            let dt = self
                .0
                .temperature_derivative_at_specific_enthalpy(at)
                .map_err(|_| DepositionError::MaterialDomain)?;
            Ok(ImplicitContactState {
                resistance_m2_k_w: 0.01 + 0.0001 * (t - 200.0),
                resistance_enthalpy_derivatives_m2_kg_k_w_j: [0.0001 * dt, 0.0],
            })
        }
    }
    for sensible_knot in [false, true] {
        let card = source();
        let knots = if sensible_knot {
            vec![
                (0.0, 200.0, 0.0),
                (40_000.0, 600.0, 0.0),
                (80_000.0, 800.0, 0.0),
            ]
        } else {
            vec![
                (0.0, 200.0, 0.0),
                (40_000.0, 600.0, 0.0),
                (60_000.0, 600.0, 1.0),
                (100_000.0, 1000.0, 1.0),
            ]
        };
        let knots = knots
            .into_iter()
            .map(|(h, t, f)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: f,
                bulk_density_kg_m3: 1000.0,
            })
            .collect();
        let caloric = if sensible_knot {
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                card.content_hash(),
                SolidLiquidPhase::Solid,
                knots,
            )
        } else {
            EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), knots)
        }
        .unwrap();
        let context = QueryPoint::new().with("T", 300.0).unwrap();
        let layer = ResolvedConductionLayer::new(
            &caloric,
            ThermalMaterialPoint {
                material: &card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            SolidLiquidPhase::Solid,
            [200.0, if sensible_knot { 800.0 } else { 600.0 }],
            || Ok(()),
        )
        .unwrap();
        let old = initial(&caloric);
        let prepared =
            PreparedLayerConduction::new(&old, vec![binding(&layer)], 1, || Ok(())).unwrap();
        let make = || prepared.implicit_link_laws(config()).next().unwrap();
        let wider = TemperatureContact(&caloric, if sensible_knot { 80_000.0 } else { 100_000.0 });
        // A phase-changing contact domain cannot be attached to a pure-solid
        // layer. The caller prepares the original contact on its shared branch.
        assert_eq!(make().with_contact(&wider).is_ok(), sensible_knot);
        let contact = TemperatureContact(&caloric, 40_000.0);
        let link = make().with_contact(&contact).unwrap();
        let geometry = ConductiveLink {
            link_id: 7,
            half_distances_m: [0.02, 0.03],
            contact_resistance_m2_k_w: 0.0,
        };
        let at = |h| {
            link.state_at_specific_enthalpies(
                [h, 10_000.0],
                geometry,
                old.contact_patch(7).unwrap(),
                &mut || Ok(()),
            )
            .unwrap()
        };
        let power = |h| {
            at(h).conductance_w_k * (caloric.temperature_at_specific_enthalpy(h).unwrap() - 300.0)
        };
        // The contact may end at an interior caloric knot of a wider bulk
        // chart. BOTH derivatives use the shared left branch there, whether
        // the excluded right branch is latent or simply has another Cp.
        let inward = (3.0 * power(40_000.0) - 4.0 * power(39_999.0) + power(39_998.0)) / 2.0;
        assert!((at(40_000.0).power_enthalpy_derivatives_w_kg_j[0] - inward).abs() < 1e-8);
        assert!(
            link.state_at_specific_enthalpies(
                [40_000.0_f64.next_up(), 10_000.0],
                geometry,
                old.contact_patch(7).unwrap(),
                &mut || Ok(()),
            )
            .is_err()
        );
    }
}


#[test]
fn integrated_contacts_preserve_radial_geometry_and_complete_tangents() {
    use fs_material::thermal::conduction::{
        LayerGeometry,
        network::GeometricLayerLinkBinding,
    };

    let card = source();
    let caloric = curve(&card);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let layer = layer(&card, &caloric, &context);
    let cylinder_length = 0.004 / (2.0 * core::f64::consts::PI * 0.02);
    let geometries = [
        LayerGeometry::cylindrical([0.01, 0.03], 0.02, cylinder_length).unwrap(),
        LayerGeometry::spherical([0.01, 0.03], 0.02).unwrap(),
    ];

    for geometry in geometries {
        let area = geometry.contact_area_m2();
        let mut events: Vec<_> = [(1, 600.0), (2, 300.0)]
            .into_iter()
            .map(|(id, t)| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(id),
                    segment_id: id,
                    body_id: 1,
                    material_id: MaterialId(1),
                    mass_kg: 0.01,
                    specific_h_j_kg: 100.0 * (t - 200.0),
                    geometry: CellGeometry {
                        start_m: [id as f64, 0.0, 0.0],
                        end_m: [id as f64, 0.1, 0.0],
                        volume_m3: 1e-5,
                        surface_area_m2: 0.1,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            })
            .collect();
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 7,
                a: CellId(1),
                b: CellId(2),
                area_m2: area,
                conductance_w_k: 0.0,
            }),
        });
        let old = DepositionSimulation::new(
            events,
            &[(MaterialId(1), &caloric)],
            DepositionBudget {
                max_cells: 2,
                max_links: 1,
                max_events: 3,
                max_substeps: 100,
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

        let prepared = PreparedLayerConduction::new_geometric(
            &old,
            vec![GeometricLayerLinkBinding {
                link: 7,
                layers: [&layer; 2],
                geometry,
                contact_area_resistance_m2_k_w: 0.0,
            }],
            1,
            || Ok(()),
        )
        .unwrap();
        assert_eq!(prepared.geometries().next().unwrap().1, geometry);

        let contact = AffineLayerContact {
            caloric: [&caloric; 2],
            bounds: [[0.0, 80_000.0]; 2],
            reference_h: [35_000.0, 15_000.0],
            resistance: 0.02,
            derivatives: [2e-7, -1e-7],
        };
        let link = prepared
            .implicit_link_laws(config())
            .next()
            .unwrap()
            .with_contact(&contact)
            .unwrap();
        let native = ConductiveLink {
            link_id: 7,
            half_distances_m: geometry.equivalent_distances_m(),
            contact_resistance_m2_k_w: 0.0,
        };
        let patch = old.contact_patch(7).unwrap();
        let power = |h: [f64; 2]| {
            let state = link
                .state_at_specific_enthalpies(h, native, patch, &mut || Ok(()))
                .unwrap();
            let t = h.map(|value| caloric.temperature_at_specific_enthalpy(value).unwrap());
            (state, state.conductance_w_k * (t[0] - t[1]))
        };
        let h = [35_000.0, 15_000.0];
        let (state, _) = power(h);
        for side in 0..2 {
            let mut plus = h;
            let mut minus = h;
            plus[side] += 0.25;
            minus[side] -= 0.25;
            let finite_difference = (power(plus).1 - power(minus).1) / 0.5;
            assert!(
                (finite_difference - state.power_enthalpy_derivatives_w_kg_j[side]).abs() < 1e-8,
                "radial geometry {geometry:?}, side {side}: finite difference {finite_difference}, tangent {:?}",
                state.power_enthalpy_derivatives_w_kg_j,
            );
        }
        assert_eq!(old.time_s(), 0.0);
    }
}
