//! A cylinder gas loading per string: `cylinder-gas-loading R N M strings COUNT INDEX...` loads only the named strings (zero-based,
//! strictly increasing). An unnamed string keeps its unloaded modes bit for bit and carries no fluid coordinates; an axial string
//! without gas loading has no axial fluid memory; a selection that names no string, a wrong one or a non-increasing one refuses.
use super::*;

const GAS: &str = "cylinder-gas-loading 0.001 30 4";

/// The canonical two-bow card (two bonded strings, gas loading, a cavity), with `tail` appended to the gas record; with
/// `longitudinal`, two retained axial modes per string replace the geometric stretching.
fn text(tail: &str, longitudinal: bool) -> String {
    let t = include_str!("../nonlinear-double-bow.material-strings")
        .replace("48000 4096 0.01 128", "48000 129 0.01 128")
        .replace(GAS, &format!("{GAS}{tail}"));
    if longitudinal {
        t.replace("geometric-stretching 0.2 64 100000 1e-11 1e-10", "longitudinal-waves 0.2 0.01 64 100000 1e-11 1e-10 2 1 2")
    } else {
        t
    }
}
fn build(cx: &Cx<'_>, gate: &CancelGate, text: &str) -> Ready {
    compile(cx, gate, input::Spec::parse(text.as_bytes()).unwrap(), 7800.0)
}
fn sizes(ready: &Ready) -> Vec<usize> {
    ready.performance.bank().system().components().iter().map(|c| c.modes().len()).collect()
}
fn modes(ready: &Ready, component: usize) -> Vec<(u64, u64)> {
    ready.performance.bank().system().components()[component].modes().iter()
        .map(|m| (m.angular_frequency_rad_s.to_bits(), m.damping_ratio.to_bits())).collect()
}

#[test]
fn a_selection_loads_only_the_named_string_and_leaves_the_other_exactly_unloaded() {
    scoped(|cx, gate| {
        let all = build(cx, gate, &text("", false));
        let one = build(cx, gate, &text(" strings 1 1", false));
        let dry = build(cx, gate, &text("", false).replace(&format!("{GAS}\n"), ""));
        let loading = one.performance.bank().source().bank().cylinder_gas_loading().unwrap();
        assert_eq!(loading.strings().len(), 1);
        assert!(loading.string(0).is_none());
        assert_eq!(loading.string(1).unwrap().string_index(), 1);
        assert_eq!(all.performance.bank().source().bank().cylinder_gas_loading().unwrap().strings().len(), 2);
        assert_eq!(sizes(&all), vec![3, 3, 2, 1, 12, 12]);
        assert_eq!(sizes(&one), vec![3, 3, 2, 1, 12], "one fluid component, for the loaded string only");
        assert_eq!(sizes(&dry), vec![3, 3, 2, 1]);
        assert_eq!(modes(&one, 0), modes(&dry, 0), "the unloaded string keeps its modes bit for bit");
        assert_ne!(modes(&one, 1), modes(&dry, 1), "the loaded string gains the gas's mass and loss");
        assert_eq!(modes(&one, 1), modes(&all, 1), "the loaded string is the all-string binding's");
    });
}

#[test]
fn the_default_loads_every_string_and_a_full_explicit_selection_is_the_same_bank() {
    scoped(|cx, gate| {
        let default = build(cx, gate, &text("", false));
        let explicit = build(cx, gate, &text(" strings 2 0 1", false));
        assert_eq!(sizes(&default), sizes(&explicit));
        for component in 0..sizes(&default).len() {
            assert_eq!(modes(&default, component), modes(&explicit, component));
        }
        let mut a = default;
        let mut b = explicit;
        assert_eq!(pressure(&mut a, 37), pressure(&mut b, 37), "naming every string changes no byte");
    });
}

#[test]
fn a_gas_selection_that_names_no_string_or_a_wrong_one_refuses_with_its_own_text() {
    for (tail, expect) in [
        (" strings 0", "names no string"),
        (" strings 1 2", "exceeds"),
        (" strings 2 1 0", "strictly increasing"),
        (" strings 1 1 1", "extra field"),
    ] {
        let err = input::Spec::parse(text(tail, false).as_bytes()).err().expect("the selection must refuse");
        assert!(err.to_string().contains(expect), "{tail}: {err}");
    }
}

#[test]
fn the_partially_loaded_bank_renders_bit_stably_across_block_sizes_and_closes_its_energy_balance() {
    scoped(|cx, gate| {
        let mut a = build(cx, gate, &text(" strings 1 1", false));
        let first = pressure(&mut a, 37);
        assert!(first.iter().any(|p| p.abs() > 1e-15), "the partially loaded bank must sound");
        for block in [1, 128] {
            assert_eq!(pressure(&mut build(cx, gate, &text(" strings 1 1", false)), block), first);
        }
        let frame = a.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

#[test]
fn an_axial_string_without_gas_has_no_axial_fluid_and_the_loaded_one_follows_the_component_order() {
    scoped(|cx, gate| {
        let mut ready = build(cx, gate, &text(" strings 1 1", true));
        let axial = ready.performance.bank().source().bank().longitudinal_waves().unwrap();
        assert_eq!(axial.strings().len(), 2, "both bonded strings carry the axial law");
        // strings 3 + 3, body 2, cavity 1, string 1's fluid 12, axial 2 + 2, string 1's axial fluid 2 x 4.
        assert_eq!(sizes(&ready), vec![3, 3, 2, 1, 12, 2, 2, 8]);
        assert_eq!(axial.string(0).unwrap().axial_component(), 5);
        assert_eq!(axial.string(0).unwrap().axial_fluid_component(), None);
        assert_eq!(axial.string(1).unwrap().axial_component(), 6);
        assert_eq!(axial.string(1).unwrap().axial_fluid_component(), Some(7));
        let first = pressure(&mut ready, 37);
        assert!(first.iter().any(|p| p.abs() > 1e-15), "the mixed axial and gas bank must sound");
        assert_eq!(pressure(&mut build(cx, gate, &text(" strings 1 1", true)), 128), first);
        let frame = ready.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}
