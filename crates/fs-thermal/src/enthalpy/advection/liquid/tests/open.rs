//! Open liquid streams use the same fixed-mass enthalpy and boundary accounts.
use super::*;

fn boundaries(inlet: CellId, outlet: CellId, rate: f64, h: f64) -> Vec<MassBoundary> {
    vec![
        MassBoundary::Inlet {
            id: 10, cell: inlet, kg_s: rate, material_id: MAT, specific_h_j_kg: h,
        },
        MassBoundary::Outlet { id: 11, cell: outlet, kg_s: rate },
    ]
}

#[test]
fn liquid_at_melting_temperature_flushes_without_transporting_stationary_solid() {
    for datum in [0.0, 1_099_511_627_776.0] {
        let law = PhaseLaw { datum, identity: 7, melting_range_k: 0.0 };
        for h in [10.0, 80.0] {
            let initial = prepare(&law, &[h], 2.0);
            let flow = Advected::new((), vec![], 2).unwrap()
                .with_liquid_carrier(carrier(&law)).unwrap()
                .with_boundaries(boundaries(CellId(0), CellId(0), 0.5, datum + 100.0), 2).unwrap();
            let (next, observed, audit) = initial.advance_observed(&flow, 2.0, config(0.1), || Ok(())).unwrap();
            assert_eq!(next.cells(), initial.cells());
            assert_eq!(audit.external_heat_j, 0.0);
            assert_eq!(audit.stored_energy_change_j, 0.0);
            assert_eq!(observed.audit().external_enthalpy_j, 0.0);
            close(observed.audit().inlet_mass_kg, 1.0, 1e-14);
            close(observed.audit().outlet_mass_kg, 1.0, 1e-14);
            assert_eq!(next.mass_kg(), initial.mass_kg());
        }
    }
}

#[test]
fn superheated_inlet_advances_the_melt_front_and_discharge_carries_liquid_enthalpy() {
    for datum in [0.0, 1_099_511_627_776.0] {
        let law = PhaseLaw { datum, identity: 7, melting_range_k: 0.0 };
        for (rate, start, increments, inlet, outlet) in [
            (1.0, [20.0, 40.0], [90.0, 10.0], CellId(0), CellId(1)),
            (-1.0, [40.0, 20.0], [10.0, 90.0], CellId(1), CellId(0)),
        ] {
            let initial = prepare(&law, &start, 2.0);
            let flow = Advected::new((), vec![MassFlow { id: 0, a: CellId(0), b: CellId(1), kg_s: rate }], 3).unwrap()
                .with_boundaries(boundaries(inlet, outlet, 1.0, datum + 200.0), 3).unwrap()
                .with_liquid_carrier(carrier(&law)).unwrap();
            // BE at dt=1: h_interior,1=(20+200)/2=110, then
            // h_interior,2=40+(110-100)=50. Only 100 J enters overall:
            // 90 J melts solid and the remaining 10 J is liquid superheat.
            let step = BackwardEulerStep::new(&initial, 1.0, || Ok(())).unwrap();
            let trial = step.evaluate_with_transport(&increments,
                |state| flow.endpoint_transport(state, || Ok(())), || Ok(())).unwrap();
            let (next, (), audit) = trial.accept(1e-10, 1e-10, || Ok(())).unwrap();
            assert_eq!(next.cells()[&inlet].specific_h_j_kg, datum + 110.0);
            assert_eq!(next.cells()[&outlet].specific_h_j_kg, datum + 50.0);
            assert_eq!(next.temperature_k(inlet).unwrap(), 505.0);
            assert_eq!(next.temperature_k(outlet).unwrap(), 500.0);
            assert_eq!(audit.external_heat_j, 100.0);
            assert_eq!(audit.stored_energy_change_j, 100.0);
            assert_eq!(audit.residual_j, 0.0);
            assert_eq!(next.mass_kg(), initial.mass_kg());
        }
    }
}

#[test]
fn hot_liquid_fills_the_latent_budget_before_sensible_heating_and_refines_in_time() {
    let law = PhaseLaw { datum: 0.0, identity: 7, melting_range_k: 0.0 };
    let initial = prepare(&law, &[20.0], 2.0);
    let flow = Advected::new((), vec![], 2).unwrap()
        .with_liquid_carrier(carrier(&law)).unwrap()
        .with_boundaries(boundaries(CellId(0), CellId(0), 0.5, 150.0), 2).unwrap();
    let exact_final_h = 150.0 - 50.0 * (-0.9_f64).exp();
    let mut previous_error = f64::INFINITY;
    for dt in [0.2, 0.1, 0.05] {
        // M=1 and q=0.5 give 25 W while h<=100. Filling the 80 J
        // latent deficit takes exactly 3.2 s; subsequent superheat obeys
        // dh/dt=0.5*(150-h), with no second latent or mass-energy credit.
        let (melted, first, first_audit) = initial.advance_observed(&flow, 3.2, config(dt), || Ok(())).unwrap();
        close(melted.cells()[&CellId(0)].specific_h_j_kg, 100.0, 2e-12);
        close(first_audit.external_heat_j, 80.0, 2e-12);
        let (next, observed, audit) = melted.advance_observed(&first, 5.0, config(dt), || Ok(())).unwrap();
        let h = next.cells()[&CellId(0)].specific_h_j_kg;
        let error = (h - exact_final_h).abs();
        assert!(error < 0.6 * previous_error, "{error} did not refine from {previous_error}");
        previous_error = error;
        close(audit.external_heat_j, h - 100.0, 2e-12);
        close(observed.audit().external_enthalpy_j, audit.external_heat_j, 2e-12);
        close(observed.audit().inlet_mass_kg, 0.9, 2e-14);
        close(observed.audit().outlet_mass_kg, 0.9, 2e-14);
        assert!(observed.audit().maximum_courant <= 0.8);
        assert_eq!(next.mass_kg(), initial.mass_kg());
    }
}

#[test]
fn fully_liquid_open_transport_matches_original_in_both_attachment_and_flow_directions() {
    let law = PhaseLaw { datum: 0.0, identity: 7, melting_range_k: 0.0 };
    let initial = prepare(&law, &[130.0, 170.0], 2.0);
    for (rate, inlet, outlet) in [(0.2, CellId(0), CellId(1)), (-0.2, CellId(1), CellId(0))] {
        let graph = Advected::new((), vec![MassFlow { id: 0, a: CellId(0), b: CellId(1), kg_s: rate }], 3).unwrap();
        let streams = boundaries(inlet, outlet, 0.2, 220.0);
        let original = graph.clone().with_boundaries(streams.clone(), 3).unwrap();
        let baseline = initial.advance_observed(&original, 1.0, config(0.1), || Ok(())).unwrap();
        for selected in [
            original.clone().with_liquid_carrier(carrier(&law)).unwrap(),
            graph.with_liquid_carrier(carrier(&law)).unwrap().with_boundaries(streams, 3).unwrap(),
        ] {
            let actual = initial.advance_observed(&selected, 1.0, config(0.1), || Ok(())).unwrap();
            assert_eq!(actual.0.cells(), baseline.0.cells());
            assert_eq!(actual.1.audit(), baseline.1.audit());
            assert_eq!(actual.2, baseline.2);
        }
    }
}

#[test]
fn open_liquid_carriage_refuses_unsupported_streams_inventory_and_unbalanced_mass() {
    let law = PhaseLaw { datum: 0.0, identity: 7, melting_range_k: 0.0 };
    let original = Advected::new((), vec![], 2).unwrap();
    let bound = carrier(&law);
    for h in [-20.0, 0.0, 75.0, 99.0, 301.0] {
        let streams = boundaries(CellId(0), CellId(0), 0.5, h);
        assert!(original.clone().with_boundaries(streams.clone(), 2).unwrap().with_liquid_carrier(bound).is_err());
        assert!(original.clone().with_liquid_carrier(bound).unwrap().with_boundaries(streams, 2).is_err());
    }
    let mut wrong_material = boundaries(CellId(0), CellId(0), 0.5, 150.0);
    if let MassBoundary::Inlet { material_id, .. } = &mut wrong_material[0] { *material_id = MaterialId(2); }
    assert!(original.clone().with_liquid_carrier(bound).unwrap().with_boundaries(wrong_material, 2).is_err());
    let valid = original.with_liquid_carrier(bound).unwrap()
        .with_boundaries(boundaries(CellId(0), CellId(0), 0.5, 150.0), 2).unwrap();
    for initial in [prepare(&law, &[-20.0], 2.0), prepare(&law, &[25.0], 4.0)] {
        assert!(initial.advance_observed(&valid, 0.1, config(0.1), || Ok(())).is_err());
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(valid.audit().substeps, 0);
    }
    let mut streams = boundaries(CellId(0), CellId(0), 0.5, 150.0);
    if let MassBoundary::Outlet { kg_s, .. } = &mut streams[1] { *kg_s = 0.25; }
    let unbalanced = valid.with_boundaries(streams, 2).unwrap();
    let initial = prepare(&law, &[25.0], 2.0);
    assert!(matches!(unbalanced.endpoint_transport(&initial, || Ok(())),
        Err(DepositionError::InvalidInput("fixed-mass advection requires balanced incoming and outgoing mass at every cell"))));
}
