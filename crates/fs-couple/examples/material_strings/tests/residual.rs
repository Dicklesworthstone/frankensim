//! `residual-flexibility exact-static-green`: the omitted bending modes in series with the bridge and every fixed
//! stop (pexfp.299, MagentaBadger 7592). The exact-string oracle is the library's; this is the input-to-PCM path.
use super::*;
use crate::stops::Profile;
use fs_couple::bowed_string::runtime::reciprocal::bank::residual;

const STOP: &str = "frankensim-string-stops-v1\nstop support 0 0.43 positive fixed 0 5000 1 0.2 10 0.001\n";

fn flexible(text: &str) -> input::Spec {
    let text = text.replace("cavity neck", "residual-flexibility exact-static-green\ncavity neck");
    input::Spec::parse(text.as_bytes()).unwrap()
}
fn pluck() -> String { include_str!("../pluck.material-strings").replace("48000 4096 0.01 128", "48000 256 0.01 128") }
// From rest, the 0.002 N port load held for 64 samples and released: no static preload.
fn struck() -> String { pluck().replace("forces preload", "forces retain").replace("force 0 0 0", "force 64 0 0") }
fn stopped(cx: &Cx<'_>, gate: &CancelGate, spec: input::Spec) -> Result<Ready> {
    spec.compile_with_profiles(cx, gate, None, Some(Profile::parse(STOP.as_bytes()).unwrap()),
        |name| Ok(pack(name, 7800.0, 200_000.0, None)))
}

#[test]
fn the_record_softens_the_bridge_through_the_omitted_modes_and_keeps_the_basis() {
    scoped(|cx, gate| {
        let flexible = flexible(&pluck());
        assert!(flexible.residual_flexibility && !spec().residual_flexibility);
        let mut stiff = compile(cx, gate, spec(), 7800.0);
        let mut soft = compile(cx, gate, flexible, 7800.0);
        let bank = soft.performance.bank();
        assert!(bank.source().bank().residual_flexibility());
        for (a, b) in stiff.performance.bank().system().components().iter().zip(bank.system().components()) {
            assert_eq!(a.modes(), b.modes(), "the residual must not retune the retained basis");
        }
        let (member, card) = (&bank.source().bank().strings()[0], &bank.cards()[0]);
        let connection = &member.connections[0];
        let condensed = residual::condense_bridge(&member.source, card, connection).unwrap();
        let g = residual::residual_compliance_m_n(card, connection.station_fraction, connection.station_fraction).unwrap();
        assert_eq!(condensed.residual_m_n, g);
        assert!(g > 0.0 && condensed.stiffness_n_m < connection.stiffness_n_m);
        // The bridge term of the stiffness bound shrinks with it.
        assert!(bank.system().linear_frequency_upper_bound_rad_s()
            < stiff.performance.bank().system().linear_frequency_upper_bound_rad_s());
        let (x, y) = (pressure(&mut stiff, 128), pressure(&mut soft, 128));
        assert!(x.iter().zip(&y).any(|(x, y)| (x-y).abs() > 1e-12));
    });
}

#[test]
fn a_fixed_stop_gains_the_bridge_cross_term_and_its_own_omitted_compliance() {
    scoped(|cx, gate| {
        let mut plain = stopped(cx, gate, input::Spec::parse(struck().as_bytes()).unwrap()).unwrap();
        let mut soft = stopped(cx, gate, flexible(&struck())).unwrap();
        assert!(plain.performance.bank().system().obstacle_residual_compliance().is_none());
        let bank = soft.performance.bank();
        let member = &bank.source().bank().strings()[0];
        let bridges = residual::residual_bridges(bank.source().bank().strings(), bank.cards()).unwrap();
        let (cross, compliance) = residual::stop_residual(&bridges, &[residual::ResidualStop {
            string: 0, station_fraction: 0.43, direction: [1.0, 0.0] }]).unwrap();
        assert_eq!(bank.system().obstacle_residual_compliance(), Some(compliance.as_slice()));
        assert!(compliance[0] > 0.0 && cross[0] != 0.0);
        assert!(bank.system().obstacle_pairs().is_empty());
        // Closure minus the cross term times the bridge link: its body part is the receiver side.
        let normal = bank.system().obstacle(0).unwrap();
        let bridge = &member.connections[0].body;
        assert_eq!(normal.right.component, bank.cards().len()+bridge.component);
        assert_eq!(normal.right.shapes, bridge.shapes.iter().map(|b| -cross[0]*b).collect::<Vec<_>>());
        let fixed = plain.performance.bank().system().obstacle(0).unwrap();
        assert!(fixed.right.shapes.iter().all(|x| *x == 0.0));
        assert_eq!(normal.left.shapes.len(), fixed.left.shapes.len());
        assert!(normal.left.shapes.iter().zip(&fixed.left.shapes).any(|(a, b)| a != b));
        let (x, y) = (pressure(&mut plain, 128), pressure(&mut soft, 128));
        assert!(y.iter().any(|p| p.abs() > 1e-12));
        assert!(x.iter().zip(&y).any(|(x, y)| (x-y).abs() > 1e-12));
    });
}

#[test]
fn what_the_residual_does_not_own_refuses_by_name() {
    scoped(|cx, gate| {
        // The released static preload's equilibrium owner has no residual.
        let error = stopped(cx, gate, flexible(&pluck())).err().unwrap().to_string();
        assert!(error.contains("residual compliance"), "{error}");
        // A stop on the bridge's own receiver is not wired.
        let receiver = flexible(&struck()).compile_with_profiles(cx, gate, None,
            Some(Profile::parse(STOP.replace("fixed", "bridge:0").as_bytes()).unwrap()),
            |name| Ok(pack(name, 7800.0, 200_000.0, None)));
        assert!(receiver.err().unwrap().to_string().contains("receiver stop"));
    });
    for record in ["residual-flexibility", "residual-flexibility exact-static-green extra", "residual-flexibility modal"] {
        let text = pluck().replace("cavity neck", &format!("{record}\ncavity neck"));
        assert!(input::Spec::parse(text.as_bytes()).is_err(), "accepted {record:?}");
    }
}
