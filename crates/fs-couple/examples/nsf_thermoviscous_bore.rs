//! Original nonlinear flow -> complete thermoviscous field -> finite ducts -> microphone.
//!
//! cargo run --release -p fs-couple --example nsf_thermoviscous_bore -- output.wav
//!
//! A finite material-owned chamber feeds a manufactured windway, free gap and
//! labium through an original transverse NSF opening. Only chamber volume and
//! heat are authored; its pressure, temperature and emitted flow are solved.
//! A separate escape shares reciprocal radiation with the original bore; its
//! signed tile flows and the bore outlets reach the same microphone.
//! A second complete acoustic cavity and closed passage react on the upstream
//! face. All three physical volumes and the finite supply retain their owners.
//! This unexecuted coarse composition is not measured anatomy or a qualified jet.

use fs_couple::bernoulli_aperture::performance::ApertureObservation;
use fs_couple::bernoulli_aperture::radiation::environment::multiport::source::ComputedBaffledInput;
use fs_couple::cavity_modal::interior::thermoviscous::{
    ThermoviscousInterior, source_medium,
    dg::{ThermoviscousDg, ThermoviscousDgLimits, load::ThermoviscousDgLoad},
};
use fs_couple::dg_port::DgPortLoad;
use fs_couple::nsf_dg::{
    NsfDgDrive, NsfDgLimits, NsfDgRadiationSolve, NsfDgReservoirRadiation,
    NsfFiniteReservoirLimits, performance::{ReservoirObservation, source_intervals_for},
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
    NsfCharacteristicBoundary, NsfCharacteristicPatch, NsfGrid3, NsfOpenState,
    NsfReservoirAperture, NsfSpec, WallThermal,
};
use fs_material::gas_species::input::{GasRecord, read_gas_record};
use fs_material::gas::GasState;
use fs_material::relaxation::{RelaxationAbsorption, input::medium_from_record};
use fs_qty::semantic::{Composition, CompositionBasis};
use fs_rep_mesh::TetComplex;
use std::fs::OpenOptions;
use std::sync::Arc;

const RATE: u32 = 48_000;
const SAMPLES: u64 = 960;
const BLOCK: usize = 32;
const SUBSTEPS: usize = 8;
const DIMS: [usize; 3] = [12, 32, 5];
const CELLS: usize = DIMS[0] * DIMS[1] * DIMS[2];
const EVALUATIONS: usize = 64;
const BAND_HZ: [f64; 2] = [20.0, 2_000.0];
const CHAMBER_VOLUME_M3: f64 = 1e-5;
const PEAK_VOLUME_RATE_M3_S: f64 = -1e-6;
const PEAK_HEAT_POWER_W: f64 = 5e-4;

// Both acoustic regions use the original complete mixed field, with explicit
// planar ports in the supplied order. This helper only selects mesh faces; it
// introduces no volume reduction or separate propagation law.
fn acoustic_field<'c>(mesh: &'c TetComplex, points: &[[f64; 3]], gas: &GasState,
    planes: &[(usize, f64)], gate: &CancelGate)
    -> Result<ThermoviscousFe<'c>, Box<dyn std::error::Error>>
{
    let geometry = HelmholtzFe::new(mesh, points, 2,
        Medium { density_kg_m3: gas.density, sound_speed_m_s: gas.sound_speed })?;
    let faces = geometry.boundary_faces();
    let openings: Vec<Vec<_>> = planes.iter().map(|&(axis, coordinate)|
        faces.iter().copied().filter(|face| face.iter().all(|&v|
            (points[v as usize][axis] - coordinate).abs() < 1e-14)).collect()).collect();
    let walls: Vec<_> = faces.into_iter().filter(|face|
        !openings.iter().any(|opening| opening.contains(face))).collect();
    Ok(ThermoviscousFe::new(mesh, points, source_medium(gas, 0.0)?,
        &ThermoviscousBoundary { no_slip_faces: walls.clone(), isothermal_faces: walls,
            ports: openings.iter().map(|opening| BoundaryMode::uniform(opening)).collect() },
        ThermoviscousLimits::default(), &mut || gate.is_requested())?)
}

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
    // Physical NSF box: [0,12dx] x [-32dx,0] x [0,5dx]. The high-y
    // lower half adjoins FE port 2 at y=0; its upper half keeps the original
    // zero incoming characteristic. The complete low-y face adjoins a separate
    // upstream FE cavity ending at y=-32dx. All three volumes are disjoint.
    // Cell geometry below is a manufactured halfway-bounce-back staircase.
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
        |[x, y, _]| {
            // Both characteristic faces require two all-fluid collar layers.
            // The low-y collar now connects to a closed acoustic assembly;
            // no solid back wall may isolate that feedback from the plenum.
            let windway = (8..13).contains(&y) && !(5..7).contains(&x);
            // Seven fluid layers separate the windway exit and edge tip.
            // The wedge thickens above x=6; y=30,31 retain the fluid collar.
            let edge = (20..30).contains(&y) && x >= 6 && x < 7 + (y - 20) / 3;
            windway || edge
        },
        |_| (1.0, [0.0; 3], gas.temperature),
    )?.with_characteristic_axis(NsfCharacteristicAxis {
        axis: 1, faces: [face; 2], maximum_relative_perturbation: 0.05,
    })?.with_characteristic_patch(true, NsfCharacteristicPatch {
        start: [0, 0], end: [6, DIMS[2]],
    }, NsfBoundaryDrive::IncomingCharacteristic(0.0))?
        .with_reservoir_apertures(vec![NsfReservoirAperture {
            axis: 0, high: false,
            // Aperture 0 is the exclusive finite supply owner.
            // The original ghost boundary requires a 3x3 tangential and
            // two-layer inward fluid collar. All nine tiles satisfy it.
            cells: (4..7).flat_map(|y| (1..4).map(move |z| [0, y, z])).collect(),
            state: NsfOpenState { density_kg_m3: gas.density,
                velocity_m_s: [0.0; 3], temperature_k: gas.temperature },
        }, NsfReservoirAperture {
            axis: 0, high: false,
            // Aperture 1 is a manufactured passive escape in the free gap.
            // Its own complete fluid collar is downstream of the windway;
            // the high-y remainder still has its explicit ambient condition.
            cells: (16..19).flat_map(|y| (1..4).map(move |z| [0, y, z])).collect(),
            state: NsfOpenState { density_kg_m3: gas.density,
                velocity_m_s: [0.0; 3], temperature_k: gas.temperature },
        }])?;
    let physical_dt = zone.law().dt_s * SUBSTEPS as f64;
    let intervals = source_intervals_for(physical_dt, RATE, SAMPLES)?;
    // The actual derived NSF clock is checked against the original acoustic
    // clock by the coupling constructor below; no clock is changed to fit it.

    let length = 6.0 * dx;
    let height = 2.0 * dx;
    let width = DIMS[2] as f64 * dx;
    let (mesh, unit_points) = kuhn_cube(1);
    let points: Vec<_> = unit_points.iter()
        .map(|p| [p[0] * length, p[1] * height, p[2] * width]).collect();
    let field = acoustic_field(&mesh, &points, &gas,
        &[(0, 0.0), (0, length), (1, 0.0)], &gate)?;
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
        .with_viscothermal(12, BAND_HZ, 0.01)?.with_molecular_relaxation(molecular.clone())?;
    let cork = DgPortLoad::new(cork, dt, 1e-11, MAX_DG_PROPAGATOR_DOFS)?;
    let acoustic = ThermoviscousDg::new(&source, prepared, bore, 1,
        ThermoviscousDgLimits { maximum_steps: intervals as u64,
            maximum_abs_pressure_pa: 100.0, maximum_total_energy_j: 1e-4,
            ..ThermoviscousDgLimits::default() }, &gate)?
        .with_additional_duct(0, cork, &gate)?;
    let load = ThermoviscousDgLoad::new(acoustic, 2, &gate)?;
    // The upstream region occupies [0,12dx] x [-36dx,-32dx] x [0,5dx].
    // Its port 1 faces the low-y NSF interface; port 0 leads farther outward
    // into a finite closed passage. These are manufactured geometry data,
    // not acquired player anatomy or an acoustically equivalent tract claim.
    let upstream_width = DIMS[0] as f64 * dx;
    let upstream_end = -(DIMS[1] as f64) * dx;
    let upstream_start = upstream_end - 4.0 * dx;
    let (upstream_mesh, upstream_unit) = kuhn_cube(1);
    let upstream_points: Vec<_> = upstream_unit.iter().map(|p|
        [p[0] * upstream_width, upstream_start + p[1] * 4.0 * dx, p[2] * width]).collect();
    let upstream_field = acoustic_field(&upstream_mesh, &upstream_points, &gas,
        &[(1, upstream_start), (1, upstream_end)], &gate)?;
    let upstream_source = ThermoviscousInterior::from_relaxation(
        upstream_field, molecular.clone(), BAND_HZ)?;
    let upstream_prepared = upstream_source.prepare_time(dt, ThermoviscousTimeLimits::default(),
        &mut || gate.is_requested())?;
    let upstream_radius = fs_math::det::sqrt(upstream_width * width / std::f64::consts::PI);
    let upstream_duct = DgDuct::new(&Duct {
        segments: vec![Segment::Cylinder { radius: upstream_radius, length: 0.04 }],
    }, &gas, DgMesh { elements: 2, order: 2 },
        DgBoundary::Port(PortId(0)), DgBoundary::Closed)?
        .with_viscothermal(12, BAND_HZ, 0.01)?.with_molecular_relaxation(molecular)?;
    let upstream_duct = DgPortLoad::new(upstream_duct, dt, 1e-11, MAX_DG_PROPAGATOR_DOFS)?;
    let upstream_acoustic = ThermoviscousDg::new(&upstream_source, upstream_prepared,
        upstream_duct, 0, ThermoviscousDgLimits {
            maximum_steps: intervals as u64, maximum_abs_pressure_pa: 100.0,
            maximum_total_energy_j: 1e-4, ..ThermoviscousDgLimits::default()
        }, &gate)?;
    let upstream_load = ThermoviscousDgLoad::new(upstream_acoustic, 1, &gate)?;
    let system = NsfThermoviscousDg::new_on_patch(zone, load, true, NsfDgLimits {
        substeps: SUBSTEPS, maximum_cells: CELLS,
        maximum_cell_steps: CELLS * SUBSTEPS * EVALUATIONS,
        maximum_evaluations: EVALUATIONS, maximum_steps: intervals as u64,
        incoming_bracket_m_s: [-0.5, 0.5],
        pressure_absolute_tolerance_pa: 1e-7, pressure_relative_tolerance: 1e-9,
        // An explicit exploratory projection budget, not a transparency claim.
        // A refusal retains every physical state; no covariance is called heat.
        projection_absolute_tolerance_j: 1e-16, projection_relative_tolerance: 0.05,
    }, &gate)?.with_finite_reservoir(0, CHAMBER_VOLUME_M3, NsfFiniteReservoirLimits {
        maximum_relative_mass_exchange: 0.01,
        maximum_relative_energy_exchange: 0.01,
        mass_absolute_tolerance_kg: 1e-18, energy_absolute_tolerance_j: 1e-10,
        relative_tolerance: 1e-9,
    }, &gate)?;
    // One declared common rigid baffle holds the original bore outlets and
    // this distinct rectangular escape. This manufactured placement is not
    // inferred from the three-dimensional source/bore body or a point source.
    let mouth_radiation = NsfDgReservoirRadiation {
        origin_m: [-0.08, 0.0], cells: [16, 16],
        terminal_radial_rings: 8, terminal_angular_points: 32,
        pressure_bracket_pa: [-20.0, 20.0],
    };
    // One explicit policy solves downstream, passive escape and upstream
    // efforts together. Both attachments require this SAME policy.
    let joint_solve = NsfDgRadiationSolve {
            difference_step: 1e-5, maximum_iterations: 12,
            newton: fs_solver::NewtonKrylovConfig {
                absolute_tolerance: 1e-8, relative_tolerance: 0.0,
                linear_restart: 3, max_linear_cycles: 8,
                forcing_minimum: 1e-5, forcing_maximum: 0.5,
                forcing_gamma: 0.9, forcing_exponent: 1.5,
                globalization: fs_solver::Globalization::LineSearch(fs_solver::LineSearchConfig {
                    armijo: 1e-4, contraction: 0.5, minimum_step: 1e-6,
                }),
            },
        };
    let system = system.with_upstream_load(upstream_load, joint_solve, &gate)?
        .with_reservoir_radiation_at(1, mouth_radiation, joint_solve, &gate)?;
    // One gentle squeeze and explicit heating of the SAME finite gas owner.
    // Its held p/T comes from its accepted inventory, then actual NSF flux
    // and original midpoint heat/volume work update that inventory together.
    // There is no prescribed jet velocity, oscillating pressure or jet card.
    let controls = (0..intervals).map(|sample| {
        let time = (sample as f64 + 0.5) * physical_dt;
        let envelope = if time < 0.004 { time / 0.004 }
            else if time < 0.012 { 1.0 }
            else if time < 0.016 { (0.016 - time) / 0.004 }
            else { 0.0 };
        NsfDgDrive::FiniteReservoirActuation { aperture: 0,
            volume_rate_m3_s: PEAK_VOLUME_RATE_M3_S * envelope,
            heat_power_w: PEAK_HEAT_POWER_W * envelope,
            maximum_relative_volume_change: 0.01 }
    }).collect();
    let mut renderer = NsfThermoviscousDgPerformance::new_with_drives(system, 1,
        ApertureObservation::TubeBaffled(CircularOutletReceiver {
            position_m: [0.0, 0.0, 0.1],
            radial_rings: mouth_radiation.terminal_radial_rings,
            angular_points: mouth_radiation.terminal_angular_points,
            maximum_frequency_hz: BAND_HZ[1],
        }), RATE, controls, BLOCK, Arc::clone(&gate))?
        .with_reservoir_observation_at(1, ReservoirObservation {
            origin_m: mouth_radiation.origin_m, angle_rad: 0.0,
        })?;
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
    println!("downstream acoustic storage={} J, balance={} J, projection covariance={} J, work residual={} J",
        frame.acoustic_stored_energy_j, frame.acoustic_balance_residual_j,
        frame.projection_covariance_j, frame.acoustic_work_residual_j);
    println!("NSF mass change={} kg, total-energy change={} J, mass projection difference={} kg",
        frame.flow_mass_change_kg, frame.flow_energy_change_j, frame.mass_projection_difference_kg);
    let upstream = frame.upstream.as_ref().ok_or("the accepted upstream field receipt is missing")?;
    println!("upstream field: pressure={} Pa, outward source flow={} m^3/s, inlet work={} J, series source work={} J",
        upstream.inlet_pressure_pa, upstream.inlet_flow_m3_s,
        upstream.inlet_work_j, upstream.source_work_j);
    println!("upstream field: acoustic storage={} J, balance={} J, pressure residual={} Pa, projection covariance={} J, work residual={} J",
        upstream.acoustic_stored_energy_j, upstream.acoustic_balance_residual_j,
        upstream.pressure_residual_pa, upstream.projection_covariance_j,
        upstream.acoustic_work_residual_j);
    let escape = frame.reservoir_radiations.first()
        .filter(|mouth| mouth.resolved_exchange.aperture == 1)
        .ok_or("the accepted passive escape receipt is missing")?;
    println!("passive escape 1: outward flow={} m^3/s, radiation pressure={} Pa, radiation work={} J",
        escape.outward_flow_m3_s, escape.radiation_pressure_pa, escape.radiation_work_j);
    println!("passive escape 1: pressure residual={} Pa, projection covariance={} J, work residual={} J",
        escape.pressure_residual_pa, escape.projection_covariance_j, escape.acoustic_work_residual_j);
    println!("passive escape 1: streamed mass={} kg, streamed total energy={} J, mass projection difference={} kg",
        escape.resolved_exchange.mass_out_kg, escape.resolved_exchange.total_energy_out_j,
        escape.mass_projection_difference_kg);
    let chamber = renderer.system().finite_reservoir_inventory()
        .ok_or("the finite source inventory is missing")?;
    let source_frame = frame.finite_reservoir.as_ref()
        .ok_or("the accepted finite source receipt is missing")?;
    println!("finite chamber: volume={} m^3, mass={} kg, pressure={} Pa, temperature={} K, sensible energy={} J",
        chamber.volume_m3, chamber.mass_kg, chamber.gas.pressure,
        chamber.gas.temperature, chamber.sensible_energy_j);
    println!("last source interval: streamed mass={} kg, streamed total energy={} J, mass residual={} kg, energy residual={} J",
        source_frame.source_exchange.mass_out_kg, source_frame.source_exchange.total_energy_out_j,
        source_frame.mass_residual_kg, source_frame.energy_residual_j);
    if let Some(actuation) = source_frame.actuation {
        println!("last source actuation: midpoint pressure={} Pa, heat={} J, work by gas={} J",
            actuation.midpoint_pressure_pa, actuation.step.supplied_heat_j,
            actuation.step.work_by_gas_j);
    }
    let fluid = renderer.system().zone();
    let windway_velocity = fluid.cell([5, 12, 2]).1[1]
        * (fluid.law().dx_m / fluid.law().dt_s);
    println!("final resolved windway axial velocity={windway_velocity} m/s");
    Ok(())
}
