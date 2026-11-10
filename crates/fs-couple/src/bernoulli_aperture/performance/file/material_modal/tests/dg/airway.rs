//! A player's airway as the DG valve card's upstream owner (`player_tract`), in the upstream tract's place.
use super::{CancelGate, MaterialModalDgValvePerformance, build, input, pack, render};
use crate::render::player_tract::TractPressure;
use crate::render::schedule::reed::PlayerTractInput;

const AIRWAY: &str = "player_tract lung\nplayer_tract_section 0.035 3e-4\nplayer_tract_section 0.02 1e-4\n";
const TRACT: &str = "upstream_tract\ndg_tube 0.04 0.001 closed 1 2 8 50 2000 0.01 1e-11 256\n";

fn card(rows: &str) -> String {
    input().replace("slit_length", &format!("{rows}slit_length"))
}

fn refusal(rows: &str) -> String {
    match MaterialModalDgValvePerformance::from_bytes_with_sources(
        card(rows).as_bytes(), 37, &CancelGate::new_clock_free(), None, |_| Ok(pack(900.0, false)),
    ) {
        Ok(_) => panic!("admitted {rows:?}"),
        Err(error) => format!("{error:?}"),
    }
}

#[test]
fn a_player_tract_block_makes_the_airway_the_upstream_owner_on_the_cards_clock_and_gas() {
    let airway_card = build(&card(AIRWAY), 37);
    let system = airway_card.renderer().system();
    assert!(system.has_upstream() && system.upstream_load().is_none());
    let tract = system.upstream_airway().expect("the card's airway");
    assert_eq!(tract.gas(), system.load().duct().inlet_gas());
    assert_eq!(tract.airway().sample_rate_hz(), 96_000);
    assert_eq!(tract.airway().receipt().sections, 2);
    assert!(matches!(tract.pressure(), TractPressure::Lung));
    assert!(tract.is_cold());
    // It plays: the valve sees the airway, so the trace differs from the card without it, and every accepted
    // sample carries the airway's ledger frame and no DG tract frame.
    let (plain, _) = render(build(&input(), 37));
    let (trace, played) = render(airway_card);
    assert_ne!(plain, trace, "the airway must change the accepted bore pressure");
    assert!(trace.iter().all(|x| x.is_finite()));
    let frame = played.last_frame().unwrap();
    assert!(frame.upstream.is_none());
    let ledger = frame.airway.expect("the airway's ledger frame");
    assert!(ledger.cumulative_absorbed_work_j >= -1e-6 * ledger.wave_energy_j);
    assert_eq!(played.system().upstream_airway().unwrap().cumulative_absorbed_work_j(),
        ledger.cumulative_absorbed_work_j);
    // Without the block the card is the original: no upstream owner at all.
    assert!(!build(&input(), 37).renderer().system().has_upstream());
}

#[test]
fn the_block_refuses_what_the_reed_performance_file_refuses_and_both_upstream_owners_at_once() {
    for (rows, what) in [
        ("player_tract unexpected\nplayer_tract_section 0.035 3e-4\n", "player_tract takes"),
        ("player_tract lung\n", "at least one player_tract_section"),
        ("player_tract mouth\nplayer_tract_section 0.035 3e-4\n", "missing field"),
        ("player_tract lung\nplayer_tract_section 0.035 0\n", "must be positive"),
        ("player_tract lung\nplayer_tract_section 0.035 3e-4\nplayer_tract_wall 15 -1 1\n", "nonnegative"),
        ("player_tract lung\nplayer_tract_section 0.035 3e-4\nplayer_airway weibel-1963 saab-2021-table-1-B 0\n",
            "1 to 24 generations"),
        ("player_tract lung\nplayer_tract_section 0.035 3e-4\nplayer_airway weibel-1963 other 3\n",
            "player_airway takes"),
    ] {
        let text = refusal(rows);
        assert!(text.contains(what), "{rows:?}: {text}");
    }
    let many = format!("player_tract lung\n{}", "player_tract_section 0.001 3e-4\n".repeat(257));
    assert!(refusal(&many).contains("at most 256 sections"));
    // A DG upstream tract and a player's airway are exclusive: the second block is an unexpected record.
    assert!(refusal(&format!("{TRACT}{AIRWAY}")).contains("unexpected record"));
    assert!(refusal(&format!("{AIRWAY}{TRACT}")).contains("unexpected record"));
}

#[test]
fn the_shared_admission_builds_one_chain_and_refuses_in_the_order_the_voice_file_reads() {
    let mut input = PlayerTractInput::new(TractPressure::Lung);
    assert_eq!(input.require_sections(), Err("player_tract needs at least one player_tract_section row"));
    for _ in 0..256 {
        input.push_section(0.001, 3e-4).unwrap();
    }
    assert_eq!(input.push_section(0.001, 3e-4), Err("player_tract takes at most 256 sections"));
    let wall = fs_phs::WallPin { surface_density: 15.0, stiffness_per_area: -1.0, resistance: 1.0 };
    assert_eq!(input.set_wall(wall), Err("player_tract_wall needs a nonnegative stiffness and resistance"));
    assert_eq!(input.set_generations(0), Err("player_airway takes 1 to 24 generations"));
    input.set_generations(3).unwrap();
    let (chain, pressure) = input.chain().unwrap();
    assert!(matches!(pressure, TractPressure::Lung));
    assert!(chain.sections().len() > 256, "the duct's 256 sections and the tree ahead of which it sits");
}
