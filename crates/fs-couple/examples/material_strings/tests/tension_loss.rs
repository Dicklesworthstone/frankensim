//! The string construction's tension loss factor eta_F, selected by row id from `tension-loss.tsv`.
use super::*;

const LAMPIS: &str = "lampis-2024-dominant-cello-g";

fn selected(record: &str) -> String {
    include_str!("../pluck.material-strings").replace("48000 4096 0.01 128","48000 256 0.01 128")
        .replace("cavity neck",&format!("{record}\ncavity neck"))
}

#[test]
fn the_table_holds_unique_cited_rows_and_the_lampis_label_as_ruled() {
    let rows = input::tension_loss_rows().unwrap();
    assert_eq!(rows.iter().map(|row| row.id.as_str()).collect::<Vec<_>>(), [LAMPIS,
        "mansour-2016-kaplan-cello-a3","mansour-2016-kaplan-cello-d3","mansour-2016-kaplan-cello-g2",
        "mansour-2016-kaplan-cello-c2","galluzzo-2017-steel-violin-e"]);
    assert_eq!(rows.iter().map(|row| row.eta_f).collect::<Vec<_>>(), [9.0e-4,2.2e-4,2.3e-4,2.0e-4,1.2e-4,5e-5]);
    assert_eq!(rows[0].source, "Lampis 2024 Dominant cello G (chrome-wound): mean Q 1064 over its first 30 modes, \
        inverted on that string's own T, mu, EI and diameter with the same air law; bending loss taken as 0, so an upper \
        bound; the violin G itself is unmeasured.");
}

#[test]
fn a_selected_row_adds_eta_f_on_the_tension_share_of_every_string_mode_and_reaches_the_pressure() {
    scoped(|cx,gate| {
        let row = input::tension_loss_rows().unwrap().into_iter().find(|row| row.id == LAMPIS).unwrap();
        let lossy = input::Spec::parse(selected(&format!("tension-loss {LAMPIS}")).as_bytes()).unwrap();
        assert_eq!(lossy.strings[0].tension_loss.as_ref(), Some(&row));
        assert_eq!(spec().strings[0].tension_loss, None);
        let mut plain = compile(cx,gate,spec(),7800.0);
        let mut lossy = compile(cx,gate,lossy,7800.0);
        let (a,b) = (&plain.performance.bank().cards()[0],&lossy.performance.bank().cards()[0]);
        assert_eq!(a.mode_count,b.mode_count);
        for k in 1..=b.mode_count {
            // The same modes: eta_F is a loss, not a stiffness.
            let omega = b.mode_omega_rad_s(k);
            assert_eq!(omega.to_bits(),a.mode_omega_rad_s(k).to_bits());
            // zeta_F = eta_F T k^2 / (2 mu omega^2), beside the pack's bending and air losses.
            let wave = k as f64*core::f64::consts::PI/b.length_m;
            let expected = 0.5*row.eta_f*b.tension_n*wave*wave/(b.linear_density_kg_m*omega*omega);
            let added = b.zetas[k-1]-a.zetas[k-1];
            assert!((added/expected-1.0).abs() < 1e-9,"mode {k}: added {added:e}, expected {expected:e}");
        }
        let x = pressure(&mut plain,128);
        let y = pressure(&mut lossy,128);
        assert!(x.iter().zip(&y).any(|(x,y)| (x-y).abs() > 1e-12));
    });
}

#[test]
fn an_absent_malformed_or_repeated_selection_refuses() {
    for record in ["tension-loss no-such-row","tension-loss",&format!("tension-loss {LAMPIS} extra"),
        &format!("tension-loss {LAMPIS}\ntension-loss {LAMPIS}")] {
        assert!(input::Spec::parse(selected(record).as_bytes()).is_err(),"accepted {record:?}");
    }
}
