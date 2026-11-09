//! Original prepared chemistry shared by finite reservoirs and nozzle probes.
use super::*;
use super::super::nozzle::{
    EquilibriumNozzleControl, publish_on_path, resolve_on_isentrope,
};

#[test]
fn g3_weak_subsonic_flow_after_depletion_uses_the_direct_current_caloric_anchor() {
    let chemistry = owner(vec![model("A", 0.01, 4.0, 1e11, 1e8)], vec![1]);
    let feed = [(id("A"), Amount::new(1.0))];
    let c = EquilibriumNozzleControl::new(controls().isentropic, Pressure::new(100.0));
    let mut path = EquilibriumIsentropeV1::new(&chemistry, Temperature::new(3000.0),
        Pressure::new(P0), &feed, c.isentropic, || false).unwrap();
    let birth = path.reference_state().clone();
    let stagnation = path.sample(Pressure::new(1000.0), || false).unwrap();
    let inlet = path.publish_current(|| false).unwrap();
    let anchor = path.caloric_anchor(|| false).unwrap();
    let back = stagnation.pressure.value() * (1.0 - 1e-10);
    let solved = resolve_on_isentrope(&mut path, stagnation, &anchor, back, 0.01, c, &mut || false).unwrap();
    let flow = publish_on_path(&mut path, inlet, stagnation, solved, back, 0.01, &mut || false).unwrap();
    let expected = (4.0 * R / 0.01)
        * (flow.sample.temperature.value() - stagnation.temperature.value());
    assert!(expected < 0.0);
    close(flow.specific_enthalpy_change_j_kg, expected, expected.abs() * 2e-15);
    let rounded_birth_difference = flow.sample.specific_enthalpy_change_j_kg
        - stagnation.specific_enthalpy_change_j_kg;
    assert!((rounded_birth_difference - expected).abs() > 1e-11,
        "the fixture must expose cancellation: {rounded_birth_difference:e}, {expected:e}");
    close(flow.velocity_m_s, (-2.0 * expected).sqrt(), 1e-14);
    close(flow.mass_flow.value(), 0.01 * flow.sample.response.density_kg_m3
        * (-2.0 * expected).sqrt(), 1e-16);
    assert!(!flow.choked);
    assert_eq!(path.reference_state(), &birth);
    assert!(flow.sample.entropy_change_j_k.abs() <= c.isentropic.entropy_tolerance_j_k);
}

#[test]
fn g1_prepared_reacting_flow_retains_birth_entropy_and_actual_local_energy() {
    let chemistry = dissociation();
    let feed = [(id("A2"), Amount::new(1.0))];
    let c = EquilibriumNozzleControl::new(controls().isentropic, Pressure::new(100.0));
    let mut path = EquilibriumIsentropeV1::new(&chemistry, Temperature::new(T0),
        Pressure::new(P0), &feed, c.isentropic, || false).unwrap();
    let birth = path.reference_state().clone();
    for (chamber_pressure, back) in [(0.8 * P0, 0.7 * P0), (0.6 * P0, 0.0)] {
        let stagnation = path.sample(Pressure::new(chamber_pressure), || false).unwrap();
        let inlet = path.publish_current(|| false).unwrap();
        let anchor = path.caloric_anchor(|| false).unwrap();
        let resolved = resolve_on_isentrope(&mut path, stagnation, &anchor, back, 0.003, c, &mut || false).unwrap();
        let flow = publish_on_path(&mut path, inlet, stagnation, resolved, back, 0.003, &mut || false).unwrap();
        let actual = (dissociation_h(flow.sample.temperature.value(), flow.sample.pressure.value())
            - dissociation_h(stagnation.temperature.value(), chamber_pressure)) / 0.002;
        close(flow.specific_enthalpy_change_j_kg, actual, 2e-5);
        close(flow.specific_enthalpy_change_j_kg + 0.5 * flow.velocity_m_s.powi(2), 0.0, 2e-8);
        check_products(&flow.outlet, &chemistry);
        assert_eq!(flow.choked, back == 0.0);
        assert_eq!(path.reference_state(), &birth);
        assert!(flow.sample.entropy_change_j_k.abs() <= c.isentropic.entropy_tolerance_j_k);
    }
}

#[test]
fn g4_caloric_anchors_reject_changed_sources_normalization_and_stale_queries() {
    let chemistry = single();
    let changed = owner(vec![model("A", 0.01, 4.0, 1.0, 0.0)], vec![1]);
    let mut path = EquilibriumIsentropeV1::new(&chemistry, Temperature::new(T0),
        Pressure::new(P0), &[(id("A"), Amount::new(1.0))], controls().isentropic, || false).unwrap();
    let anchor = path.caloric_anchor(|| false).unwrap();
    for (owner, amount) in [(&changed, 1.0), (&chemistry, 2.0)] {
        let other = EquilibriumIsentropeV1::new(owner, Temperature::new(T0),
            Pressure::new(P0), &[(id("A"), Amount::new(amount))], controls().isentropic, || false).unwrap();
        assert!(matches!(other.specific_enthalpy_change_from(&anchor, || false),
            Err(EquilibriumExpansionErrorV1::Invalid(_))));
    }
    assert!(path.sample(Pressure::new(f64::NAN), || false).is_err());
    assert!(path.caloric_anchor(|| false).is_err());
    assert!(path.specific_enthalpy_change_from(&anchor, || false).is_err());
    assert_eq!(anchor.sample().pressure, Pressure::new(P0));
}

#[test]
fn g4_anchor_capture_and_direct_energy_refuse_final_cancellation_without_rebasing() {
    let chemistry = dissociation();
    let mut path = EquilibriumIsentropeV1::new(&chemistry, Temperature::new(T0),
        Pressure::new(P0), &[(id("A2"), Amount::new(1.0))], controls().isentropic, || false).unwrap();
    let birth = path.reference_state().clone();
    let mut polls = 0;
    let anchor = path.caloric_anchor(|| { polls += 1; false }).unwrap();
    let mut count = 0;
    assert!(matches!(path.caloric_anchor(|| { count += 1; count == polls }),
        Err(EquilibriumExpansionErrorV1::Equilibrium(CompositionEquilibriumErrorV1::Cancelled))));
    assert_eq!(path.caloric_anchor(|| false).unwrap(), anchor);
    path.sample(Pressure::new(0.8 * P0), || false).unwrap();
    polls = 0;
    let expected = path.specific_enthalpy_change_from(&anchor, || { polls += 1; false }).unwrap();
    count = 0;
    assert!(matches!(path.specific_enthalpy_change_from(&anchor, || { count += 1; count == polls }),
        Err(EquilibriumExpansionErrorV1::Equilibrium(CompositionEquilibriumErrorV1::Cancelled))));
    assert_eq!(path.specific_enthalpy_change_from(&anchor, || false).unwrap(), expected);
    assert_eq!(path.reference_state(), &birth);
}
