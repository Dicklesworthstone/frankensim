//! Original selected-stream calorics, signed separation and one energy debit.
use super::*;
use crate::enthalpy::advection::SelectiveOutletLaw;

fn laws() -> [MixtureLaw; 3] {
    [MixtureLaw::new(vec![0.75, 0.25]), MixtureLaw::new(vec![0.375, 0.625]),
        MixtureLaw::new(vec![0.0, 1.0])]
}

fn species(flow: Advected<()>) -> Advected<()> {
    flow.with_stationary_species(vec!["A".into(), "B".into()], vec![
        InletSpecies { boundary_id: 10, mass_fractions: vec![0.375, 0.625] },
    ], 1, 2).unwrap()
}

fn stream(laws: &[MixtureLaw], inlet_t: f64, selected_t: f64, bind: bool) -> Advected<()> {
    let flow = Advected::new((), vec![], 0).unwrap().with_boundaries(vec![
        MassBoundary::Inlet { id: 10, cell: CellId(0), kg_s: 1.0,
            material_id: MaterialId(1), specific_h_j_kg: laws[1].h(inlet_t) },
        MassBoundary::Outlet { id: 11, cell: CellId(0), kg_s: 0.5 },
        MassBoundary::Outlet { id: 12, cell: CellId(0), kg_s: 0.5 },
    ], 3).unwrap().with_selective_outlets(vec![SelectiveOutlet {
        boundary_id: 12, mass_fractions: vec![0.0, 1.0],
        specific_h_j_kg: laws[2].h(selected_t),
        enthalpy_reference: "analytic-species-reference".into(),
    }], 2).unwrap();
    let flow = species(flow);
    if bind { flow.with_selective_outlet_laws(vec![SelectiveOutletLaw {
        boundary_id: 12, material_id: MaterialId(2),
    }]).unwrap() } else { flow }
}

#[test]
fn g1_selected_isothermal_separation_keeps_its_entropy_deficit_under_reference_shifts() {
    for formation in [[1000.0, -400.0], [2.0_f64.powi(40), -2.0_f64.powi(41)],
        [-2.0_f64.powi(40), 2.0_f64.powi(41)]] {
        let mut laws = laws();
        for law in &mut laws { law.formation = formation; }
        let initial = prepared(&laws, 300.0);
        let flow = stream(&laws, 300.0, 300.0, true);
        let expected = laws[1].mixing_entropy() - 0.5 * laws[0].mixing_entropy();
        assert!(expected > 0.0);
        let (_, endpoint) = flow.endpoint_transport_with_entropy(&initial, || Ok(())).unwrap();
        assert!((endpoint.boundary_entropy_input_w_k - expected).abs() < 1e-14);
        assert!((endpoint.upwind_mixing_w_k + expected).abs() < 1e-14);
        let (next, observed, energy, entropy) = initial.advance_advected_with_entropy(
            &flow, 0.25, cfg(), || Ok(()),
        ).unwrap();
        assert_eq!(next.cells(), initial.cells());
        assert_eq!(energy.external_heat_j, 0.0);
        assert_eq!(entropy.thermal_storage_change_j_k, 0.0);
        let carried = entropy.advection.unwrap();
        assert!((carried.boundary_entropy_input_j_k - 0.25 * expected).abs() < 1e-14);
        // A prescribed isothermal separator needs external thermodynamic
        // closure. The missing work/entropy must remain visible, not clipped.
        assert!((carried.upwind_mixing_j_k + 0.25 * expected).abs() < 1e-14);
        assert!(carried.balance_residual_j_k.abs() < 1e-14);
        assert!(entropy.temporal_defect_j_k.abs() < 1e-14);
        assert_eq!(observed.stationary_species(), flow.stationary_species());
        assert_eq!(flow.audit().substeps, 0);
    }
}

#[test]
fn g1_selected_stream_uses_its_own_phase_entropy_and_never_adds_latent_heat_twice() {
    let mut laws = laws();
    // Explicit source-owned coexistence offset: h_liquid=h_gas-L and
    // s_liquid=s_gas-L/Tm. No phase equilibrium is inferred by transport.
    laws[2].phase_enthalpy_j_kg = -3000.0;
    laws[2].phase_entropy_j_kg_k = -10.0;
    let flow = stream(&laws, 300.0, 300.0, true);
    let initial = prepared(&laws, 300.0);
    assert!(laws[0].temperature_k(laws[2].h(300.0)).is_err());
    let (plain, _, plain_energy) = initial.advance_observed(&flow, 0.1, cfg(), || Ok(())).unwrap();
    let (next, observed, energy, entropy) = initial.advance_advected_with_entropy(
        &flow, 0.1, cfg(), || Ok(()),
    ).unwrap();
    assert_eq!(next.cells(), plain.cells());
    assert_eq!(energy, plain_energy);
    assert!((energy.external_heat_j - 150.0).abs() < 1e-12);
    assert!((next.temperature_k(CellId(0)).unwrap() - 325.0).abs() < 1e-13);
    let separation = laws[1].mixing_entropy() - 0.5 * laws[0].mixing_entropy();
    let carried = entropy.advection.unwrap();
    assert!((carried.boundary_entropy_input_j_k - 0.1 * (5.0 + separation)).abs() < 1e-13);
    assert!((carried.upwind_mixing_j_k + 0.1 * separation).abs() < 1e-13);
    assert!(carried.balance_residual_j_k.abs() < 1e-13);
    let storage = 6.0 * (325.0_f64 / 300.0).ln();
    assert!((entropy.thermal_storage_change_j_k - storage).abs() < 1e-13);
    assert!((entropy.temporal_defect_j_k - (storage - 0.5)).abs() < 1e-13);
    assert!(observed.entropy_audit().unwrap().other_energy_at_cell_temperature_j_k.abs() < 1e-13);
    assert!(!entropy.has_unclassified_sources);
}

#[test]
fn g1_hot_feed_and_selective_removal_share_the_original_explicit_quadrature() {
    let laws = laws();
    let flow = stream(&laws, 600.0, 300.0, true);
    let initial = prepared(&laws, 300.0);
    let (plain, _, plain_energy) = initial.advance_observed(&flow, 0.1, cfg(), || Ok(())).unwrap();
    let (next, _, energy, entropy) = initial.advance_advected_with_entropy(
        &flow, 0.1, cfg(), || Ok(()),
    ).unwrap();
    assert_eq!(next.cells(), plain.cells());
    assert_eq!(energy, plain_energy);
    assert!((next.temperature_k(CellId(0)).unwrap() - 322.5).abs() < 1e-13);
    assert!((energy.external_heat_j - 135.0).abs() < 1e-12);
    let boundary = 4.5 * 2.0_f64.ln() + laws[1].mixing_entropy()
        - 0.5 * laws[0].mixing_entropy();
    let production = 4.5 - boundary;
    assert!(production > 0.0);
    let carried = entropy.advection.unwrap();
    assert!((carried.boundary_entropy_input_j_k - 0.1 * boundary).abs() < 1e-13);
    assert!((carried.upwind_mixing_j_k - 0.1 * production).abs() < 1e-13);
    assert!(carried.balance_residual_j_k.abs() < 1e-13);
    assert!((entropy.temporal_defect_j_k
        - (6.0 * (322.5_f64 / 300.0).ln() - 0.45)).abs() < 1e-13);
}

#[test]
fn g4_selected_entropy_requires_exact_sources_and_invalidates_replaced_bindings() {
    let laws = laws();
    let initial = prepared(&laws, 300.0);
    let unbound = stream(&laws, 400.0, 200.0, false);
    assert!(initial.advance_observed(&unbound, 0.1, cfg(), || Ok(())).is_ok());
    assert!(matches!(unbound.endpoint_transport_with_entropy(&initial, || Ok(())),
        Err(DepositionError::InvalidInput("selective outlet entropy requires every original stream law"))));
    for ids in [vec![], vec![10], vec![11], vec![99], vec![12, 12]] {
        assert!(unbound.clone().with_selective_outlet_laws(ids.into_iter().map(|boundary_id|
            SelectiveOutletLaw { boundary_id, material_id: MaterialId(2) }).collect()).is_err());
    }
    let bound = stream(&laws, 400.0, 200.0, true);
    let descriptor = SelectiveOutlet { boundary_id: 12, mass_fractions: vec![0.0, 1.0],
        specific_h_j_kg: laws[2].h(200.0), enthalpy_reference: laws[2].reference.into() };
    let replaced = species(bound.clone().with_selective_outlets(vec![descriptor.clone()], 2).unwrap());
    assert!(replaced.endpoint_transport_with_entropy(&initial, || Ok(())).is_err());
    let replaced = species(bound.clone().with_boundaries(bound.boundaries().to_vec(), 3).unwrap()
        .with_selective_outlets(vec![descriptor], 2).unwrap());
    assert!(replaced.endpoint_transport_with_entropy(&initial, || Ok(())).is_err());
    let missing = bound.clone().with_selective_outlet_laws(vec![SelectiveOutletLaw {
        boundary_id: 12, material_id: MaterialId(99),
    }]).unwrap();
    assert!(missing.endpoint_transport(&initial, || Ok(())).is_err());
    let bad_domain = stream(&laws, 400.0, -20.0, true);
    assert!(bad_domain.endpoint_transport(&initial, || Ok(())).is_err());

    let mut polls = 0;
    initial.advance_advected_with_entropy(&bound, 0.5, cfg(), || { polls += 1; Ok(()) }).unwrap();
    let mut visited = 0;
    assert!(matches!(initial.advance_advected_with_entropy(&bound, 0.5, cfg(), || {
        visited += 1;
        if visited + 1 >= polls { Err(DepositionError::Cancelled) } else { Ok(()) }
    }), Err(DepositionError::Cancelled)));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(bound.audit().substeps, 0);
    assert!(bound.entropy_audit().is_none());

    for invalid_source in 0..5 {
        let mut changed = self::laws();
        match invalid_source {
            0 => changed[2].fractions = vec![1.0, 0.0],
            1 => changed[2].reference = "different-phase-energy-reference",
            2 => changed[2].entropy_reference = None,
            3 => changed[2].entropy_reference = Some("different-phase-entropy-reference"),
            _ => changed[2].pressure_ratio = 0.0,
        }
        let state = prepared(&changed, 300.0);
        let flow = stream(&changed, 400.0, 200.0, true);
        assert!(flow.endpoint_transport_with_entropy(&state, || Ok(())).is_err());
        assert!(state.advance_advected_with_entropy(&flow, 0.1, cfg(), || Ok(())).is_err());
        assert_eq!(state.time_s(), 0.0);
        assert_eq!(flow.audit().substeps, 0);
    }
}
