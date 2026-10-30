//! G0/G3 owner composition. Synthetic materials are numerical fixtures, not
//! measured nozzle or instrument validation, and inlet pressure is diagnostic.
use super::*;
use crate::pcm_wav::observation::PressureRenderer;
use fs_material::gas_species::input::{ThermochemicalPhase, read_gas_record};

mod pressure;
mod temperature;
mod upstream;
mod molecular;
mod patch;
mod reservoir;
mod flue_geometry;
mod aperture;
mod sgs;

#[path = "../../../render/schedule/reed/field_input/tests/support.rs"]
mod support;

const RATE: u32 = 96_000;
const SAMPLES: usize = 32;

fn reference(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, pressure).unwrap()
}

fn fixture(gas: GasState, declaration: &str) -> String {
    // This fixture explicitly chooses one physical mesh at its declared clock.
    // Tests that replace material properties keep that mesh text unchanged.
    let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
    let radius = (9.0 * dx * dx / core::f64::consts::PI).sqrt();
    let energy = if declaration.starts_with("thermochemical ") {
        "energy_solve 1e-9 96\n"
    } else {
        ""
    };
    format!(
        "{RESOLVED_DG_PERFORMANCE_SCHEMA}\n\
clock {RATE} 1\n\
output {SAMPLES} 20\n\
gas {declaration}\n\
grid 12 3 3 {dx:.17e} 0 0 0\n\
fluid_limits 108 250 600 0.2 1 17 1\n\
{energy}bulk_viscosity 0.73 authored-source-for-numerical-coupling-test\n\
wall adiabatic\n\
characteristic 0 high 0.05\n\
initial rest\n\
solids 0\n\
states 0\n\
coupling 64 6912 -0.5 0.5 1e-7 1e-9 0 1\n\
dg_tube 0.1 {radius:.17e} matched 1 2 32 50 2000 0.01 1e-11 512\n\
observation inlet-pressure\n\
boundary 1\n\
hold {SAMPLES} -0.1\n\
end\n"
    )
}

fn read(text: &str, block: usize) -> ResolvedDgFilePerformance {
    ResolvedDgFilePerformance::from_bytes(
        text.as_bytes(),
        block,
        Arc::new(CancelGate::new_clock_free()),
    )
    .unwrap()
}

fn render(
    performance: ResolvedDgFilePerformance,
    block: usize,
) -> (Vec<f64>, NsfDgPerformance<D3q27>) {
    let samples = performance.info().samples as usize;
    let mut renderer = performance.into_renderer();
    renderer.validate_sample_rate(RATE).unwrap();
    let mut pressure = vec![0.0; samples];
    for chunk in pressure.chunks_mut(block) {
        renderer.block(chunk).unwrap();
    }
    assert!(pressure.iter().all(|p| p.is_finite()));
    (pressure, renderer)
}

#[test]
fn physical_initial_state_geometry_and_owned_source_reach_the_original_owners() {
    let gas = reference(100000.0);
    let file = {
        let source = fixture(gas, "dry-air 300 100000").replace(
            "states 0",
            &format!(
                "states 1\nstate 3 1 1 {:.17e} 0.002 -0.001 0.0003 300.02",
                gas.density * 1.00001
            ),
        );
        let file = read(&source, 8);
        assert_eq!(
            file.info().input_hash,
            hash_domain(
                "org.frankensim.fs-couple.resolved-dg-performance.v1",
                source.as_bytes()
            )
        );
        file
    }; // The text backing every source record has gone away.
    assert_eq!(file.info().dimensions, [12, 3, 3]);
    assert_eq!(file.info().initial_state_cells, 1);
    assert_eq!(file.info().solid_cells, 0);
    let system = file.renderer().system();
    let zone = system.zone();
    assert_eq!(zone.cell_count(), 108);
    assert_eq!(
        zone.bulk_ratio_source(),
        "authored-source-for-numerical-coupling-test"
    );
    assert_eq!(system.reference_gas(), &gas);
    assert_eq!(system.load().duct().inlet_gas(), &gas);
    let (rho, velocity, temperature) = zone.cell([3, 1, 1]);
    assert!((rho - 1.00001).abs() < 2e-12);
    let unit = zone.law().dx_m / zone.law().dt_s;
    for (observed, expected) in velocity.into_iter().zip([0.002, -0.001, 0.0003]) {
        assert!((observed * unit - expected).abs() < 1e-10);
    }
    assert!((temperature - 300.02).abs() < 1e-9);
    assert!(
        (system.area_m2() - 9.0 * zone.law().dx_m.powi(2)).abs()
            < 8.0 * f64::EPSILON * system.area_m2()
    );
    assert_eq!(zone.steps(), 0);
    assert_eq!(file.renderer().remaining_samples(), SAMPLES as u64);
}

#[test]
fn boundary_playback_is_partition_invariant_and_pressure_and_geometry_change_physics() {
    let text = fixture(reference(100000.0), "dry-air 300 100000");
    let (baseline, complete) = render(read(&text, 8), 8);
    let (partitioned, _) = render(read(&text, 3), 3);
    assert_eq!(baseline, partitioned);
    assert_eq!(complete.system().zone().steps(), SAMPLES);
    assert_eq!(
        complete.system().load().state().step_count(),
        SAMPLES as u64
    );
    assert_eq!(complete.remaining_samples(), 0);
    assert!(complete.validate_sample_count(1).is_err());
    let peak = baseline.iter().copied().map(f64::abs).fold(0.0, f64::max);
    assert!(peak > 1e-6, "actual physical pressure must be live");

    // Pressure alone changes density, impedance and transport in the SAME mesh.
    let higher_pressure = text.replace("gas dry-air 300 100000", "gas dry-air 300 200000");
    let (dense, dense_system) = render(read(&higher_pressure, 8), 8);
    assert_eq!(
        complete.system().zone().law().dx_m,
        dense_system.system().zone().law().dx_m
    );
    assert_eq!(
        dense_system.system().reference_gas().density,
        2.0 * complete.system().reference_gas().density
    );
    assert!(dense.iter().copied().map(f64::abs).fold(0.0, f64::max) > 1.5 * peak);

    let constricted = text.replace(
        "solids 0",
        "solids 6\nsolid 5 0 0\nsolid 5 0 1\nsolid 5 0 2\nsolid 5 2 0\nsolid 5 2 1\nsolid 5 2 2",
    );
    let file = read(&constricted, 8);
    assert_eq!(file.info().solid_cells, 6);
    assert!(
        file.renderer()
            .system()
            .zone()
            .thermochemical_state([5, 0, 0])
            .is_err()
    );
    let (restricted, _) = render(file, 8);
    assert!(
        baseline
            .iter()
            .zip(&restricted)
            .any(|(a, b)| (a - b).abs() > 1e-6),
        "solid geometry must change the solved bore pressure, not only the file hash"
    );
}

#[test]
fn temperature_changes_the_physical_clock_without_resizing_geometry_or_output_time() {
    let text = fixture(reference(100000.0), "dry-air 300 100000");
    let baseline = read(&text, 8);
    let changed = read(&text.replace("gas dry-air 300 100000", "gas dry-air 310 100000"), 8);
    assert_eq!(baseline.renderer().system().zone().law().dx_m,
        changed.renderer().system().zone().law().dx_m);
    assert_ne!(baseline.renderer().system().zone().law().dt_s,
        changed.renderer().system().zone().law().dt_s);
    assert_eq!(baseline.info().samples, changed.info().samples);
    assert_eq!(baseline.info().sample_rate_hz, changed.info().sample_rate_hz);
    let expected = super::super::source_intervals_for(
        changed.renderer().system().load().propagator().dt(), RATE, SAMPLES as u64).unwrap();
    let (a, _) = render(baseline, 8);
    let (b, changed) = render(changed, 3);
    assert_ne!(a, b);
    assert_eq!(changed.samples_rendered(), SAMPLES as u64);
    assert_eq!(changed.system().accepted_steps(), expected as u64);
    assert_eq!(changed.remaining_samples(), 0);
}

#[test]
fn ambiguous_cells_unowned_sources_and_incomplete_boundary_histories_refuse() {
    let text = fixture(reference(100000.0), "dry-air 300 100000");
    let density = reference(100000.0).density;
    for invalid in [
        text.replace("clock 96000 1", "clock 0 1"),
        text.replace("clock 96000 1", "clock 96000 0"),
        text.replace("solids 0", "solids 2\nsolid 5 1 1\nsolid 5 1 1"),
        text.replace("solids 0", "solids 1\nsolid 0 1 1"),
        text.replace(
            "states 0",
            &format!("states 2\nstate 3 1 1 {density} 0 0 0 300\nstate 3 1 1 {density} 0 0 0 300"),
        ),
        text.replace("solids 0", "solids 1\nsolid 5 1 1").replace(
            "states 0",
            &format!("states 1\nstate 5 1 1 {density} 0 0 0 300"),
        ),
        text.replace("hold 32 -0.1", "hold 31 -0.1"),
        text.replace("boundary 1", "boundary 0"),
        text.replace("gas dry-air 300 100000", "gas dry-air 300 100000 extra"),
        text.replace(
            "gas dry-air 300 100000",
            "gas mixture 300 100000 mole 1 N2 1",
        ),
        text.replace(
            "gas dry-air 300 100000",
            "gas thermochemical 300 100000 mole 1 missing-source 1",
        ),
    ] {
        assert!(
            ResolvedDgFilePerformance::from_bytes(
                invalid.as_bytes(),
                8,
                Arc::new(CancelGate::new_clock_free())
            )
            .is_err()
        );
    }
    let gate = Arc::new(CancelGate::new_clock_free());
    gate.request();
    assert!(ResolvedDgFilePerformance::from_bytes(text.as_bytes(), 8, gate).is_err());
}

#[test]
fn retained_cancellation_gate_does_not_advance_either_physical_owner_on_refusal() {
    let text = fixture(reference(100000.0), "dry-air 300 100000");
    let gate = Arc::new(CancelGate::new_clock_free());
    let file =
        ResolvedDgFilePerformance::from_bytes(text.as_bytes(), 8, Arc::clone(&gate)).unwrap();
    let mut renderer = file.into_renderer();
    renderer.block(&mut [0.0]).unwrap();
    let before = (
        renderer.system().zone().state_hash(),
        renderer.system().load().state().clone(),
        renderer.samples_rendered(),
    );
    gate.request();
    assert!(renderer.block(&mut [0.0]).is_err());
    assert_eq!(
        before,
        (
            renderer.system().zone().state_hash(),
            renderer.system().load().state().clone(),
            renderer.samples_rendered()
        )
    );
    assert!(matches!(
        renderer.block(&mut [0.0]),
        Err(RenderError::Poisoned)
    ));
}

fn thermochemical_file(cp: f64) -> ResolvedDgFilePerformance {
    let packs = support::packs("synthetic-alpha", 0.028, cp, 200.0);
    let bindings = [packs.binding()];
    let sources = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let declaration = "thermochemical 300 100000 mole 1 synthetic-alpha 1";
    let original = read_gas_record(
        &mut declaration.split_ascii_whitespace(),
        None,
        Some(&sources),
    )
    .unwrap();
    // This mesh is fixed independently of the selected source's gas constant.
    let text = fixture(reference(100000.0), declaration);
    let file = ResolvedDgFilePerformance::from_bytes_with_sources(
        text.as_bytes(),
        8,
        Arc::new(CancelGate::new_clock_free()),
        Some(&sources),
        |_| Err("this fixture has no material wall attachment".into()),
    )
    .unwrap();
    assert_eq!(file.gas_records(), &[original]);
    let material = file
        .renderer()
        .system()
        .zone()
        .thermochemical_state([3, 1, 1])
        .unwrap()
        .unwrap();
    assert_eq!(
        material.source_packs()[0].species_pack,
        packs.species.content_hash()
    );
    assert_eq!(
        material.source_packs()[0].nasa9_pack,
        packs.nasa9.content_hash()
    );
    assert_eq!(
        material.source_packs()[0].transport_pack,
        packs.transport.content_hash()
    );
    let phase_checked = ThermochemicalSources {
        phase: ThermochemicalPhase::LiquidWaterLimit,
        ..sources
    };
    assert!(
        ResolvedDgFilePerformance::from_bytes_with_sources(
            text.as_bytes(),
            8,
            Arc::new(CancelGate::new_clock_free()),
            Some(&phase_checked),
            |_| Err("no material wall attachment".into())
        )
        .is_err()
    );
    file // All borrowed pack data and registry lifetimes end here.
}

#[test]
fn original_thermochemical_sources_survive_ingress_and_caloric_data_changes_sound() {
    let baseline = thermochemical_file(3.5);
    let higher_cp = thermochemical_file(4.5);
    assert_eq!(
        baseline.renderer().system().zone().law().dx_m,
        higher_cp.renderer().system().zone().law().dx_m
    );
    assert_eq!(
        baseline.renderer().system().reference_gas().density,
        higher_cp.renderer().system().reference_gas().density
    );
    assert!(
        higher_cp.renderer().system().reference_gas().sound_speed
            < baseline.renderer().system().reference_gas().sound_speed
    );
    let (a, solved) = render(baseline, 8);
    let (b, _) = render(higher_cp, 8);
    assert_ne!(
        a, b,
        "changing only source calorics must reach the original resolved and DG solves"
    );
    assert!(
        solved
            .system()
            .zone()
            .thermochemical_state([3, 1, 1])
            .unwrap()
            .is_some()
    );
}

#[test]
fn composition_changes_gas_constant_and_clock_without_changing_the_physical_mesh() {
    let build = |beta_fraction: f64| {
        let a = support::packs("synthetic-alpha", 0.028, 3.5, 200.0);
        let b = support::packs("synthetic-beta", 0.044, 4.0, 200.0);
        let bindings = [a.binding(), b.binding()];
        let sources = ThermochemicalSources { species: &bindings, water_species: None,
            phase: ThermochemicalPhase::FrozenGas };
        let declaration = format!("thermochemical 300 100000 mole 2 synthetic-alpha {} synthetic-beta {beta_fraction}",
            1.0 - beta_fraction);
        let text = fixture(reference(100000.0), &declaration);
        ResolvedDgFilePerformance::from_bytes_with_sources(text.as_bytes(), 8,
            Arc::new(CancelGate::new_clock_free()), Some(&sources),
            |_| Err("no material wall attachment".into())).unwrap()
    };
    let a = build(0.2);
    let b = build(0.8);
    assert_eq!(a.renderer().system().zone().law().dx_m, b.renderer().system().zone().law().dx_m);
    assert_ne!(a.renderer().system().zone().law().dt_s, b.renderer().system().zone().law().dt_s);
    assert_ne!(a.renderer().system().reference_gas().specific_gas_constant,
        b.renderer().system().reference_gas().specific_gas_constant);
    let (a, a_state) = render(a, 8);
    let (b, b_state) = render(b, 3);
    assert_ne!(a, b);
    assert_eq!(a.len(), SAMPLES);
    assert_eq!(b.len(), SAMPLES);
    assert_eq!(a_state.samples_rendered(), b_state.samples_rendered());
    for system in [a_state.system(), b_state.system()] {
        assert_eq!(system.zone().thermochemical_state([3, 1, 1]).unwrap().unwrap().source_packs().len(), 2);
    }
}

#[test]
fn shipped_resolved_example_uses_the_original_exterior_pressure_and_gas_owner() {
    let file = read(
        include_str!("../../../../examples/resolved-dg.performance"),
        8,
    );
    assert_eq!(file.info().dimensions, [16, 3, 3]);
    assert_eq!(file.info().solid_cells, 4);
    assert!(matches!(
        file.renderer().observation(),
        ApertureObservation::TubeBaffled(_)
    ));
    assert!(file.renderer().baffled_receiver().is_some());
    assert!(file.renderer().observer_gas().is_some());
    assert_eq!(file.gas_records().len(), 2);
}
