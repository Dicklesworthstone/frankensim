//! Resolved thermal/molecular fields between two finite DG ducts and a microphone.
//!
//! Run with a fresh destination:
//! cargo run --release -p fs-couple --example thermoviscous_bore -- output.wav
//!
//! This is a manufactured small-domain pressure-transient composition, not an
//! instrument, a resolved jet or a converged mesh. All lengths, forcing and
//! numerical choices below are explicit. Pressure is never peak-normalized.

use fs_couple::bernoulli_aperture::performance::ApertureObservation;
use fs_couple::bernoulli_aperture::radiation::environment::multiport::source::ComputedBaffledInput;
use fs_couple::cavity_modal::interior::thermoviscous::{
    ThermoviscousInterior, source_medium,
    dg::{
        ThermoviscousDg, ThermoviscousDgLimits,
        performance::{
            ThermoviscousDgControl, ThermoviscousDgPerformance,
            ThermoviscousDgPerformanceConfig,
        },
    },
};
use fs_couple::dg_port::DgPortLoad;
use fs_couple::pcm_wav::baffled::CircularOutletReceiver;
use fs_couple::pcm_wav::stream::{Pcm16WavStream, render_pressure_pcm16};
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, MAX_DG_PROPAGATOR_DOFS, PortId};
use fs_duct::{Duct, Segment};
use fs_exec::CancelGate;
use fs_feec::acoustics::helmholtz::multiport::BoundaryMode;
use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_feec::acoustics::thermoviscous::{
    ThermoviscousBoundary, ThermoviscousFe, ThermoviscousLimits,
    time::ThermoviscousTimeLimits,
};
use fs_feec::fixtures::kuhn_cube;
use fs_material::gas_species::input::read_gas_record;
use fs_material::relaxation::{RelaxationAbsorption, input::medium_from_record};
use std::fs::OpenOptions;
use std::sync::Arc;

const RATE: u32 = 48_000;
const SAMPLES: usize = 960;
const BLOCK: usize = 32;
const LENGTH_M: f64 = 0.015;
const HEIGHT_M: f64 = 0.002;
const WIDTH_M: f64 = 0.004;
const BAND_HZ: [f64; 2] = [20.0, 2_000.0];

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut arguments = std::env::args_os().skip(1);
    let output = arguments.next().ok_or("usage: thermoviscous_bore OUTPUT.wav")?;
    if arguments.next().is_some() {
        return Err("usage: thermoviscous_bore OUTPUT.wav".into());
    }
    let gate = Arc::new(CancelGate::new_clock_free());
    // Original explicit named-species source; the same composition supplies the
    // existing composition-limited terrestrial-air molecular law in the field
    // and BOTH original ducts.
    let record = read_gas_record(
        &mut "mixture 293.15 101325 mole 4 N2 0.78084 O2 0.209476 Ar 0.009365 CO2 0.000319"
            .split_ascii_whitespace(),
        None,
        None,
    )?;
    let gas = *record.state();
    let molecular = RelaxationAbsorption::iso9613_terrestrial_air()?
        .evaluate(&medium_from_record(&record)?)?;

    let (mesh, unit_points) = kuhn_cube(1);
    let points: Vec<_> = unit_points
        .iter()
        .map(|p| [p[0] * LENGTH_M, p[1] * HEIGHT_M, p[2] * WIDTH_M])
        .collect();
    let geometry = HelmholtzFe::new(
        &mesh,
        &points,
        2,
        Medium { density_kg_m3: gas.density, sound_speed_m_s: gas.sound_speed },
    )?;
    let faces = geometry.boundary_faces();
    let end = |x: f64| {
        faces.iter().copied().filter(|face| {
            face.iter().all(|&v| (points[v as usize][0] - x).abs() < 1e-14)
        }).collect::<Vec<_>>()
    };
    let inlet = end(0.0);
    let outlet = end(LENGTH_M);
    let drive: Vec<_> = faces.iter().copied().filter(|face| {
        face.iter().all(|&v| (points[v as usize][1] - HEIGHT_M).abs() < 1e-14)
    }).collect();
    let walls: Vec<_> = faces.iter().copied()
        .filter(|face| !inlet.contains(face) && !outlet.contains(face) && !drive.contains(face)).collect();
    let field = ThermoviscousFe::new(
        &mesh,
        &points,
        source_medium(&gas, 0.0)?,
        &ThermoviscousBoundary {
            no_slip_faces: walls.clone(),
            isothermal_faces: walls,
            ports: vec![
                BoundaryMode::uniform(&inlet), BoundaryMode::uniform(&outlet),
                BoundaryMode::uniform(&drive),
            ],
        },
        ThermoviscousLimits::default(),
        &mut || gate.is_requested(),
    )?;
    let source = ThermoviscousInterior::from_relaxation(field, molecular.clone(), BAND_HZ)?;
    let prepared = source.prepare_time(
        1.0 / f64::from(RATE),
        ThermoviscousTimeLimits::default(),
        &mut || gate.is_requested(),
    )?;

    // This connection retains only the uniform pressure/total-flow coordinate.
    // Equal area is its projection condition; it does not establish acoustic
    // equivalence of a rectangular opening and a circular pipe transition.
    let radius = (HEIGHT_M * WIDTH_M / std::f64::consts::PI).sqrt();
    let duct = DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius, length: 0.05 }] },
        &gas,
        DgMesh { elements: 2, order: 2 },
        DgBoundary::Port(PortId(0)),
        DgBoundary::Port(PortId(1)),
    )?
    .with_viscothermal(12, BAND_HZ, 0.01)?
    .with_molecular_relaxation(molecular.clone())?;
    let load = DgPortLoad::new_with_computed_baffled(
        duct,
        RATE,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        &record,
        ComputedBaffledInput::parse(
            "20 65 4 8 12 20 0.1 67108864 2 100 1e-8 64 1048576 257",
        )?,
        &gate,
    )?;
    // A second original finite duct provides a closed cork-side return. Its
    // propagation, molecular reservoirs and wall memory continue independently;
    // the full FE cross-port map transmits its reaction to the radiating bore.
    let cork_duct = DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius, length: 0.02 }] },
        &gas,
        DgMesh { elements: 2, order: 2 },
        DgBoundary::Port(PortId(0)),
        DgBoundary::Closed,
    )?
    .with_viscothermal(12, BAND_HZ, 0.01)?
    .with_molecular_relaxation(molecular)?;
    let cork_load = DgPortLoad::new(
        cork_duct, 1.0 / f64::from(RATE), 1e-11, MAX_DG_PROPAGATOR_DOFS,
    )?;
    let system = ThermoviscousDg::new(
        &source,
        prepared,
        load,
        1,
        ThermoviscousDgLimits {
            maximum_steps: SAMPLES as u64,
            maximum_abs_pressure_pa: 100.0,
            maximum_total_energy_j: 1e-4,
            ..ThermoviscousDgLimits::default()
        },
        &gate,
    )?.with_additional_duct(0, cork_load, &gate)?;
    // An independent top opening supplies an authored 2 ms, one-pascal
    // triangular pressure pulse. Both internal duct slots stay zero in the
    // controls: their reactions are solved together, including returning waves.
    // The pulse's broad spectrum and this coarse mesh require refinement;
    // BAND_HZ is an admission range, not an assertion that input is band-limited.
    let controls = (0..SAMPLES).map(|sample| {
        let pressure = if sample < 96 {
            1.0 - ((sample as f64 + 0.5) / 48.0 - 1.0).abs()
        } else { 0.0 };
        ThermoviscousDgControl {
            traction_pa: vec![0.0, 0.0, pressure],
            wall_forces_n_per_sqrt_kg: Vec::new(),
        }
    }).collect();
    let mut renderer = ThermoviscousDgPerformance::new(
        system,
        controls,
        ApertureObservation::TubeBaffled(CircularOutletReceiver {
            position_m: [0.0, 0.0, 0.1],
            radial_rings: 4,
            angular_points: 16,
            maximum_frequency_hz: BAND_HZ[1],
        }),
        ThermoviscousDgPerformanceConfig {
            sample_rate_hz: RATE,
            maximum_block: BLOCK,
            maximum_control_values: 3 * SAMPLES,
        },
        gate.clone(),
    )?;
    // Cold source and observer admission precedes filesystem publication.
    // create_new refuses an existing output instead of replacing it.
    let file = OpenOptions::new().write(true).create_new(true).open(output)?;
    let mut stream = Pcm16WavStream::new(file, RATE, 1.0, BLOCK)?;
    let progress = render_pressure_pcm16(
        &mut renderer, &mut stream, &gate, &mut [0.0; BLOCK], SAMPLES as u64,
    )?;
    let (_, summary) = stream.finish()?;
    let frame = renderer.last_frame().ok_or("no accepted physical interval")?;
    println!("{progress:?}; {} samples, {} clipped, {} Pa full scale",
        summary.samples, summary.clipped_samples, summary.full_scale_pa);
    println!("last physical interval: stored={} J, balance={} J, interface={} J",
        frame.stored_energy_j, frame.balance_residual_j, frame.interface_work_residual_j);
    for cork in &frame.additional_ducts {
        println!("finite cork passage at FE port {}: pressure={} Pa, flow={} m^3/s, interface={} J",
            cork.port, cork.mean_pressure_pa, cork.mean_flow_m3_s, cork.interface_work_residual_j);
    }
    Ok(())
}
