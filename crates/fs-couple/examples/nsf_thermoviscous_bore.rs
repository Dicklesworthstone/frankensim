//! Original nonlinear flow -> complete thermoviscous field -> finite ducts -> microphone.
//!
//! cargo run --release -p fs-couple --example nsf_thermoviscous_bore -- output.wav
//!
//! This manufactured channel/constriction is a composition example, not measured
//! anatomy or a qualified flute jet. The pressure history acts on the original
//! NSF source boundary; pressure and flow at the acoustic field are solved.

use fs_couple::bernoulli_aperture::performance::ApertureObservation;
use fs_couple::bernoulli_aperture::radiation::environment::multiport::source::ComputedBaffledInput;
use fs_couple::cavity_modal::interior::thermoviscous::{
    ThermoviscousInterior, source_medium,
    dg::{ThermoviscousDg, ThermoviscousDgLimits, load::ThermoviscousDgLoad},
};
use fs_couple::dg_port::DgPortLoad;
use fs_couple::nsf_dg::{
    NsfDgDrive, NsfDgLimits, performance::source_intervals_for,
    thermoviscous::{NsfThermoviscousDg, performance::NsfThermoviscousDgPerformance},
};
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
use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::{
    NamedMixtureNsfGas, NsfBoundaryDrive, NsfCharacteristicAxis,
    NsfCharacteristicBoundary, NsfGrid3, NsfSpec, WallThermal,
};
use fs_material::gas_species::input::{GasRecord, read_gas_record};
use fs_material::relaxation::{RelaxationAbsorption, input::medium_from_record};
use fs_qty::semantic::{Composition, CompositionBasis};
use std::fs::OpenOptions;
use std::sync::Arc;

const RATE: u32 = 48_000;
const SAMPLES: u64 = 480;
const BLOCK: usize = 32;
const SUBSTEPS: usize = 8;
const DIMS: [usize; 3] = [8, 8, 2];
const CELLS: usize = DIMS[0] * DIMS[1] * DIMS[2];
const EVALUATIONS: usize = 64;
const BAND_HZ: [f64; 2] = [20.0, 2_000.0];

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut arguments = std::env::args_os().skip(1);
    let output = arguments.next().ok_or("usage: nsf_thermoviscous_bore OUTPUT.wav")?;
    if arguments.next().is_some() {
        return Err("usage: nsf_thermoviscous_bore OUTPUT.wav".into());
    }
    let gate = Arc::new(CancelGate::new_clock_free());
    let record = read_gas_record(
        &mut "mixture 293.15 101325 mole 4 N2 0.78084 O2 0.209476 Ar 0.009365 CO2 0.000319"
            .split_ascii_whitespace(), None, None,
    )?;
    let GasRecord::Mixture(reference, parts) = record else {
        return Err("the example requires its original named mixture".into());
    };
    // Resolve the shared reader's retained mole fractions once for every owner,
    // exactly as the native NSF performance reader does. The NSF volume keeps
    // this source law and evaluates its actual cell temperatures and pressures.
    let composition = Composition::new(CompositionBasis::MoleFraction,
        parts.iter().map(|(_, fraction)| *fraction).collect())?;
    let named = NamedMixtureNsfGas::new(
        parts.iter().map(|(species, _)| *species).collect(), composition,
        reference.temperature, reference.pressure,
    )?;
    let gas = *named.reference_state();
    let record = GasRecord::Mixture(gas, parts);
    let molecular = RelaxationAbsorption::iso9613_terrestrial_air()?
        .evaluate(&medium_from_record(&record)?)?;
    let dt = 1.0 / f64::from(RATE);
    let dx = fs_math::det::sqrt(3.0 * gas.specific_gas_constant * gas.temperature)
        * dt / SUBSTEPS as f64;
    let face = NsfCharacteristicBoundary {
        reference_pressure_pa: gas.pressure,
        reference_temperature_k: gas.temperature,
        reference_velocity_m_s: [0.0; 3],
        incoming_characteristic_m_s: 0.0,
    };
    // Physical NSF box: [0,8dx] x [-8dx,0] x [0,2dx]. The high-y
    // characteristic face adjoins the FE low-y opening at y=0. These are
    // distinct volumes. A one-cell-thick interior constriction is manufactured.
    let zone = NsfGrid3::<D3q27>::new(
        NsfSpec {
            dims: DIMS, periodic: [false; 3], gas: named,
            temperature_k: gas.temperature, pressure_pa: gas.pressure, dx_m: dx,
            temperature_window_k: [273.15, 323.15], maximum_mach: 0.2,
            maximum_temperature_ratio: 0.1, bulk_to_shear_ratio: 0.73,
            bulk_ratio_source: "authored numerical composition; Cramer2012-N2-293K value, not mixture validation",
            maximum_cells: CELLS, wall: WallThermal::Adiabatic,
            workers: 1, transport_points: 17,
        },
        |cell| cell[1] == 3 && (cell[0] < 2 || cell[0] >= 6),
        |_| (1.0, [0.0; 3], gas.temperature),
    )?.with_characteristic_axis(NsfCharacteristicAxis {
        axis: 1, faces: [face; 2], maximum_relative_perturbation: 0.05,
    })?;
    let physical_dt = zone.law().dt_s * SUBSTEPS as f64;
    let intervals = source_intervals_for(physical_dt, RATE, SAMPLES)?;
    // The actual derived NSF clock is checked against the original acoustic
    // clock by the coupling constructor below; no clock is changed to fit it.

    let length = DIMS[0] as f64 * dx;
    let height = 2.0 * dx;
    let width = DIMS[2] as f64 * dx;
    let (mesh, unit_points) = kuhn_cube(1);
    let points: Vec<_> = unit_points.iter()
        .map(|p| [p[0] * length, p[1] * height, p[2] * width]).collect();
    let geometry = HelmholtzFe::new(&mesh, &points, 2,
        Medium { density_kg_m3: gas.density, sound_speed_m_s: gas.sound_speed })?;
    let faces = geometry.boundary_faces();
    let opening = |axis: usize, coordinate: f64| faces.iter().copied().filter(|face|
        face.iter().all(|&v| (points[v as usize][axis] - coordinate).abs() < 1e-14))
        .collect::<Vec<_>>();
    let cork_face = opening(0, 0.0);
    let bore_face = opening(0, length);
    let source_face = opening(1, 0.0);
    let walls: Vec<_> = faces.iter().copied().filter(|face|
        !cork_face.contains(face) && !bore_face.contains(face) && !source_face.contains(face)).collect();
    let field = ThermoviscousFe::new(&mesh, &points, source_medium(&gas, 0.0)?,
        &ThermoviscousBoundary { no_slip_faces: walls.clone(), isothermal_faces: walls,
            ports: vec![BoundaryMode::uniform(&cork_face), BoundaryMode::uniform(&bore_face),
                BoundaryMode::uniform(&source_face)] },
        ThermoviscousLimits::default(), &mut || gate.is_requested())?;
    let source = ThermoviscousInterior::from_relaxation(field, molecular.clone(), BAND_HZ)?;
    let prepared = source.prepare_time(dt, ThermoviscousTimeLimits::default(),
        &mut || gate.is_requested())?;
    // The circular connections retain one uniform pressure/total-flow mode.
    // Equal area does not prove rectangular-to-circular shape equivalence.
    let radius = fs_math::det::sqrt(height * width / std::f64::consts::PI);
    let bore = DgDuct::new(&Duct {
        segments: vec![Segment::Cylinder { radius, length: 0.05 }],
    }, &gas, DgMesh { elements: 2, order: 2 },
        DgBoundary::Port(PortId(0)), DgBoundary::Port(PortId(1)))?
        .with_viscothermal(12, BAND_HZ, 0.01)?
        .with_molecular_relaxation(molecular.clone())?;
    let bore = DgPortLoad::new_with_computed_baffled(bore, RATE, 1e-11,
        MAX_DG_PROPAGATOR_DOFS, &record, ComputedBaffledInput::parse(
            "20 65 4 8 12 20 0.1 67108864 2 100 1e-8 64 1048576 257")?, &gate)?;
    let cork = DgDuct::new(&Duct {
        segments: vec![Segment::Cylinder { radius, length: 0.02 }],
    }, &gas, DgMesh { elements: 2, order: 2 },
        DgBoundary::Port(PortId(0)), DgBoundary::Closed)?
        .with_viscothermal(12, BAND_HZ, 0.01)?.with_molecular_relaxation(molecular)?;
    let cork = DgPortLoad::new(cork, dt, 1e-11, MAX_DG_PROPAGATOR_DOFS)?;
    let acoustic = ThermoviscousDg::new(&source, prepared, bore, 1,
        ThermoviscousDgLimits { maximum_steps: intervals as u64,
            maximum_abs_pressure_pa: 100.0, maximum_total_energy_j: 1e-4,
            ..ThermoviscousDgLimits::default() }, &gate)?
        .with_additional_duct(0, cork, &gate)?;
    let load = ThermoviscousDgLoad::new(acoustic, 2, &gate)?;
    let system = NsfThermoviscousDg::new(zone, load, true, NsfDgLimits {
        substeps: SUBSTEPS, maximum_cells: CELLS,
        maximum_cell_steps: CELLS * SUBSTEPS * EVALUATIONS,
        maximum_evaluations: EVALUATIONS, maximum_steps: intervals as u64,
        incoming_bracket_m_s: [-0.05, 0.05],
        pressure_absolute_tolerance_pa: 1e-7, pressure_relative_tolerance: 1e-9,
        // An explicit exploratory projection budget, not a transparency claim.
        // A refusal retains every physical state; no covariance is called heat.
        projection_absolute_tolerance_j: 1e-16, projection_relative_tolerance: 0.05,
    }, &gate)?;
    let controls = (0..intervals).map(|sample| {
        let time = (sample as f64 + 0.5) * physical_dt;
        let pressure = if time < 0.002 { time / 0.002 }
            else if time < 0.006 { 1.0 }
            else if time < 0.008 { (0.008 - time) / 0.002 }
            else { 0.0 };
        NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(pressure))
    }).collect();
    let mut renderer = NsfThermoviscousDgPerformance::new_with_drives(system, 1,
        ApertureObservation::TubeBaffled(CircularOutletReceiver {
            position_m: [0.0, 0.0, 0.1], radial_rings: 4, angular_points: 16,
            maximum_frequency_hz: BAND_HZ[1],
        }), RATE, controls, BLOCK, Arc::clone(&gate))?;
    let file = OpenOptions::new().write(true).create_new(true).open(output)?;
    let mut stream = Pcm16WavStream::new(file, RATE, 1.0, BLOCK)?;
    let progress = render_pressure_pcm16(&mut renderer, &mut stream, &gate,
        &mut [0.0; BLOCK], SAMPLES)?;
    let (_, summary) = stream.finish()?;
    let frame = renderer.last_frame().ok_or("no accepted coupled interval")?;
    println!("{progress:?}; {} samples, {} clipped, {} Pa full scale",
        summary.samples, summary.clipped_samples, summary.full_scale_pa);
    println!("NSF dx={dx} m, coupled dt={physical_dt} s; pressure={} Pa, flow={} m^3/s",
        frame.inlet_pressure_pa, frame.inlet_flow_m3_s);
    println!("acoustic storage={} J, balance={} J, projection covariance={} J, work residual={} J",
        frame.acoustic_stored_energy_j, frame.acoustic_balance_residual_j,
        frame.projection_covariance_j, frame.acoustic_work_residual_j);
    println!("NSF mass change={} kg, total-energy change={} J, mass projection difference={} kg",
        frame.flow_mass_change_kg, frame.flow_energy_change_j, frame.mass_projection_difference_kg);
    Ok(())
}
