use super::*;
#[path = "boundary.rs"]
mod boundary;
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
