//! Native finite-source ownership and its full escape/upstream composition.
use super::*;

const DECLARATION: &str = "thermochemical 300 100000 mole 1 synthetic-alpha 1";
const FINITE: &str = "finite_reservoir 0 0.0001 0.1 0.1 1e-16 1e-9 1e-9";

fn build(joint: bool, block: usize, edit: impl FnOnce(String) -> String)
    -> Result<ResolvedDgFilePerformance, ResolvedDgInputError>
{
    let packs = support::packs("synthetic-alpha", 0.028, 3.5, 200.0);
    let bindings = [packs.binding()];
    let sources = ThermochemicalSources {
        species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    };
    let record = read_gas_record(&mut DECLARATION.split_ascii_whitespace(),
        None, Some(&sources)).unwrap();
    let density = 100002.0 / (record.state().specific_gas_constant * 300.15);
    let opening = format!("reservoir_aperture 1 low {density:.17e} 0 0 0 300.15 1");
    let text = if joint {
        include_str!("../../../../../examples/resolved-dg-separate-openings.performance")
            .replace("gas dry-air 300 100000", &format!("gas {DECLARATION}"))
            .replace("dg_exterior 300 100000 0", &format!("dg_exterior {DECLARATION}"))
            .replace("bulk_viscosity ", "energy_solve 1e-9 96\nbulk_viscosity ")
            .replace("reservoir_aperture 1 low reference 1", &opening)
    } else {
        fixture(*record.state(), DECLARATION).replace("initial rest",
            &format!("{opening}\nreservoir_cell 5 0 1\ninitial rest"))
    };
    let text = text.replace("initial rest", &format!("{FINITE}\ninitial rest"));
    let (head, _) = text.split_once("\nboundary ").unwrap();
    let samples = if joint { 64 } else { SAMPLES };
    let text = edit(format!("{head}\nboundary finite-reservoir 0 1\nhold {samples}\nend\n"));
    ResolvedDgFilePerformance::from_bytes_with_sources(text.as_bytes(), block,
        Arc::new(CancelGate::new_clock_free()), Some(&sources),
        |_| Err("this fixture has no material wall attachment".into()))
}

#[test]
fn finite_holds_reach_the_retained_source_and_replay_across_callback_partitions() {
    let file = build(false, 8, |text| text).unwrap();
    let initial_mass = file.renderer().system().finite_reservoir().unwrap().mass_kg();
    assert_eq!(file.renderer().system().finite_reservoir_aperture(), Some(0));
    assert!(file.renderer().drive_history().iter().all(|drive|
        *drive == NsfDgDrive::FiniteReservoir { aperture: 0 }));
    let (whole, a) = render(file, 8);
    let (parts, b) = render(build(false, 3, |text| text).unwrap(), 3);
    assert_eq!(whole, parts);
    assert_eq!(a.system().zone().state_hash(), b.system().zone().state_hash());
    assert_eq!(a.system().finite_reservoir(), b.system().finite_reservoir());
    assert_ne!(a.system().finite_reservoir().unwrap().mass_kg(), initial_mass);
    let frame = a.last_frame().unwrap();
    let account = frame.finite_reservoir.as_ref().unwrap();
    assert_eq!(account.source_exchange, frame.reservoir_exchanges[0]);
    assert_eq!(account.final_mass_kg, a.system().finite_reservoir().unwrap().mass_kg());
    match a.system().zone().medium() {
        NsfGasMedium::Thermochemical(source) =>
            assert_eq!(a.system().finite_reservoir().unwrap().model(), source.model()),
        _ => panic!("native finite source lost its original thermochemical model"),
    }
}

#[test]
fn native_source_requires_owned_inventory_and_exclusive_controls() {
    assert!(build(false, 4, |text| text.replace(format!("{FINITE}\n").as_str(), "")).is_err());
    for (from, to) in [
        ("finite_reservoir 0 ", "finite_reservoir 1 "),
        ("0.0001 0.1 0.1", "0 0.1 0.1"),
        ("0.0001 0.1 0.1", "0.0001 0.1 1"),
        ("boundary finite-reservoir 0 1", "boundary finite-reservoir 1 1"),
        ("boundary finite-reservoir 0 1", "boundary aperture 0 1"),
        ("hold 32\n", "hold 32 100002 300.15\n"),
        ("hold 32\n", "ramp 32\n"),
    ] {
        assert!(build(false, 4, |text| text.replace(from, to)).is_err(), "{from}");
    }
    assert!(build(false, 4, |text| text.replace(DECLARATION, "dry-air 300 100000")
        .replace("energy_solve 1e-9 96\n", "")).is_err());
    assert!(build(true, 4, |text| text
        .replace("finite_reservoir 0 ", "finite_reservoir 1 ")
        .replace("boundary finite-reservoir 0 ", "boundary finite-reservoir 1 ")).is_err());
}

#[test]
fn native_finite_breath_advances_with_passive_escape_and_both_axial_ducts() {
    let mut renderer = build(true, 4, |text| text).unwrap().into_renderer();
    assert_eq!(renderer.system().reservoir_radiation_aperture(), Some(1));
    assert_eq!(renderer.reservoir_observation_aperture(), Some(1));
    assert_eq!(renderer.system().zone().reservoir_aperture_count(), 2);
    renderer.validate_sample_rate(8000).unwrap();
    let mut pressure = [0.0];
    renderer.block(&mut pressure).unwrap();
    assert!(pressure[0].is_finite());
    let frame = renderer.last_frame().unwrap();
    assert_eq!(frame.drive, NsfDgDrive::FiniteReservoir { aperture: 0 });
    assert_eq!(frame.finite_reservoir.as_ref().unwrap().source_exchange.aperture, 0);
    assert_eq!(frame.reservoir_radiation.as_ref().unwrap().resolved_exchange.aperture, 1);
    assert!(frame.upstream.is_some());
    assert!(renderer.system().accepted_steps() > 0);
}
