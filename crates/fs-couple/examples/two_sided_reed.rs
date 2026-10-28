//! two_sided_reed NEW_OUTPUT.wav UPSTREAM_LENGTH_M SECOND_BLADE_DENSITY_KG_M3 [EXTERIOR_K]
//! Optional trailing --supply-slit LENGTH_M WIDTH_M HEIGHT_M TRANSVERSE_CELLS.
//! Authored plate sections and tract termination, not identified cane/anatomy.
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::paired::{BladeState, PairedAperture, PairedApertureLimits};
use fs_couple::bernoulli_aperture::plate::{PlateApertureOptions, PlateApertureReduction};
use fs_couple::bernoulli_aperture::two_sided::{ApertureDuctSpec,
    TwoSidedLimits, ValveSide, combined_inlet_state, supply::combined_inlet_state_with_supply};
use fs_couple::bernoulli_aperture::slit_channel::{SlitChannel, SlitChannelSpec, SlitExit};
use fs_couple::bernoulli_aperture::two_sided::observer::TwoSidedPressure;
use fs_couple::bernoulli_aperture::two_sided::observer::schedule::ScheduledTwoSidedPressure;
use fs_couple::bernoulli_aperture::performance::AperturePerformanceConfig;
use fs_scenario::gesture::GestureSchedule;
use fs_couple::pcm_wav::{encode_pcm16_wav, baffled::CircularOutletReceiver,
    observation::PressureRenderer};
use fs_couple::pcm_wav::baffled::openings::duct::{DuctOpening,DuctRadiationSpec,DuctReconstruction};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_plate::{AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateSection, SliceOptions};
use std::io::Write;

fn plate(density: f64, gate: &CancelGate) -> Result<PlateApertureReduction, String> {
    let section = PlateSection::isotropic(4e9, 0.3, 0.0003, density).map_err(|e| e.to_string())?;
    let mesh = PlateMesh::rectangle(0.025, 0.01, 6, 2);
    let supports = mesh.nodes.iter().enumerate().filter(|(_, p)| p.0 == 0.0).map(|(i, _)| i).collect();
    let slit_edges = mesh.boundary_edges().into_iter().filter(|&(a, b)|
        (mesh.nodes[a].0-0.025).abs() < 1e-12 && (mesh.nodes[b].0-0.025).abs() < 1e-12)
        .map(|(a, b)| [a, b]).collect();
    let chart = PlateChart::with_boundary_and_regions(mesh, section, supports, vec![])
        .map_err(|e| e.to_string())?;
    PlateApertureReduction::from_chart(chart, PlateApertureOptions {
        assembly: AssemblyOptions { pretension: 0.0, support: EdgeSupport::Clamped },
        eigenvalue_window: (1.0, (std::f64::consts::TAU*450.0).powi(2)), mode_index: 0,
        eigensolver: SliceOptions::default(), max_nodes: 100, max_triangles: 100,
        slit_edges, rest_opening_m: 0.0002, damping_ratio: 0.02,
        max_slit_mode_variation: 0.25, max_slope: 0.1,
    }, gate).map_err(|e| e.to_string())
}
fn run() -> Result<(), String> {
    let mut args: Vec<_> = std::env::args().skip(1).collect();
    let supply_input = if args.len() >= 5 && args[args.len()-5] == "--supply-slit" {
        let start = args.len()-5;
        let length = args[start+1].parse::<f64>().map_err(|_| "supply length must be metres")?;
        let width = args[start+2].parse::<f64>().map_err(|_| "supply width must be metres")?;
        let height = args[start+3].parse::<f64>().map_err(|_| "supply gap must be metres")?;
        let cells = args[start+4].parse::<usize>().map_err(|_| "supply transverse cells must be an integer")?;
        args.truncate(start);
        Some((length, width, height, cells))
    } else { None };
    let (destination, length, density, exterior_k) = match args.as_slice() {
        [destination,length,density] => (destination,length,density,293.15),
        [destination,length,density,temperature] => (destination,length,density,
            temperature.parse::<f64>().map_err(|_| "exterior temperature must be kelvin".to_string())?),
        _ => return Err("usage: two_sided_reed NEW_OUTPUT.wav UPSTREAM_LENGTH_M SECOND_BLADE_DENSITY_KG_M3 [EXTERIOR_K] [--supply-slit LENGTH_M WIDTH_M HEIGHT_M TRANSVERSE_CELLS]".into()),
    };
    let length: f64 = length.parse().map_err(|_| "upstream length must be metres".to_string())?;
    let density: f64 = density.parse().map_err(|_| "blade density must be kg/m3".to_string())?;
    if !length.is_finite() || !(0.03..=1.0).contains(&length)
        || !density.is_finite() || density <= 0.0 || !exterior_k.is_finite() || exterior_k <= 0.0 {
        return Err("require upstream length in [0.03,1] m, positive finite density and exterior kelvin".into());
    }
    let gate = CancelGate::new_clock_free();
    let rate = 48_000;
    let dt = 1.0/f64::from(rate);
    let gas = GasState::try_new_moist_air(293.15, 101_325.0, 0.0).map_err(|e| format!("{e:?}"))?;
    let exterior = GasState::try_new_moist_air(exterior_k, 101_325.0, 0.0).map_err(|e| format!("{e:?}"))?;
    let make_side = |length, terminal| ApertureDuctSpec {
        graph: TubeNetworkSpec {
            nodes: vec![NetworkNode::Inlet, terminal],
            sections: vec![TubeSection { nodes: [0, 1], length_m: length,
                radius_m: 0.007, max_length_error_m: 0.005 }],
            sound_speed_m_s: gas.sound_speed, max_wave_memory_bytes: 1 << 20,
        }, gases: vec![gas], viscothermal: vec![],
    };
    let sides = [make_side(length, NetworkNode::Termination { reflection: 0.7 }),
        // Explicit cold endpoint: the common environment installs radiation here.
        make_side(0.25, NetworkNode::Termination { reflection: 0.0 })];
    let supply = supply_input.map(|(length_m, width_m, height_m, transverse_cells)|
        SlitChannel::from_gas(&gas, SlitChannelSpec { length_m, width_m, height_m,
            transverse_cells, sample_period_s: dt, exit: SlitExit::StaticPressure }))
        .transpose().map_err(|e| e.to_string())?;
    let (inlet, load) = match &supply {
        Some(channel) => combined_inlet_state_with_supply(&sides, channel),
        None => combined_inlet_state(&sides),
    }.map_err(|e| e.to_string())?;
    let plates = [plate(900.0, &gate)?, plate(density, &gate)?];
    for (i, p) in plates.iter().enumerate() {
        eprintln!("blade={i}; mass_kg={:.10e}; stiffness_n_m={:.10e}; pressure_area_m2={:.10e}",
            p.mass_kg(), p.stiffness_n_m(), p.pressure_area_m2());
    }
    let contact = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "authored effective opposed-blade contact, not measured cane".into())
        .map_err(|e| e.to_string())?;
    let valve = PairedAperture::from_plates(plates, inlet, load, dt, 4800,
        PairedApertureLimits { maximum_abs_displacement_m: 0.01,
            maximum_abs_velocity_m_s: 100.0, maximum_energy_j: 1.0, max_evaluations: 256,
            velocity_tolerance_m_s: 1e-13, relative_tolerance: 1e-11 },
        [BladeState::default(); 2], contact).map_err(|e| e.to_string())?;
    let limits = TwoSidedLimits {
        maximum_energy_j: 1.0, energy_absolute_tolerance_j: 1e-16, energy_relative_tolerance: 1e-7,
        flow_absolute_tolerance_m3_s: 1e-16, flow_relative_tolerance: 1e-8,
    };
    let radiation = DuctRadiationSpec {
        openings: vec![DuctOpening {network:1,node:1,center_m:[0.0;2]}], exterior:Some(exterior),
        receiver:CircularOutletReceiver {position_m:[0.0,0.0,0.25],radial_rings:4,
            angular_points:16,maximum_frequency_hz:1000.0},
        reconstruction:DuctReconstruction::FiniteDifference {sample_rate_hz:rate},
    };
    let (source,boundaries) = match supply {
        Some(channel) => TwoSidedPressure::from_baffled_environment_with_supply(
            valve.into(), sides, [None,None], channel, limits, None, radiation, 240, &gate),
        None => TwoSidedPressure::from_baffled_environment(
            valve.into(), sides, [None,None], limits, None, radiation, 240, &gate),
    }.map_err(|e| e.to_string())?;
    if let Some(channel) = source.system().supply_channel() {
        let channel = channel.slit().expect("this example declares a planar supply passage");
        eprintln!("supply_length_m={}; supply_width_m={}; supply_height_m={}; transverse_cells={}; gas_density={}; gas_viscosity={}",
            channel.spec().length_m, channel.spec().width_m, channel.spec().height_m,
            channel.spec().transverse_cells, channel.gas().density, channel.gas().dynamic_viscosity);
    }
    eprintln!("represented_upstream_m={}; represented_downstream_m={}",
        source.system().represented_sections(ValveSide::Upstream)[0].represented_length_m,
        source.system().represented_sections(ValveSide::Downstream)[0].represented_length_m);
    eprintln!("interior_k=293.15; exterior_k={exterior_k}; exterior_density={}; exterior_speed={}; radiation_loads={}",
        exterior.density,exterior.sound_speed,boundaries.len());
    // The common gesture compiler, not a host-block envelope, owns the source.
    let schedule = GestureSchedule::from_canonical_bytes(include_bytes!("two_sided_reed.gesture"))
        .map_err(|e| e.to_string())?;
    let mut renderer = ScheduledTwoSidedPressure::new(source, schedule, AperturePerformanceConfig {
        sample_rate_hz: rate, samples: 4800, max_block: 240,
        max_compile_work: 100_000, max_controls: 128,
    }).map_err(|e| e.to_string())?;
    let mut pressure = vec![0.0; 4800];
    for block in pressure.chunks_mut(240) {
        renderer.block(block).map_err(|e| e.to_string())?;
    }
    // Fixed scale across every material/geometry experiment; never normalize.
    let (wav, clips) = encode_pcm16_wav(&pressure, rate, 0.005).map_err(|e| e.to_string())?;
    let mut file = std::fs::OpenOptions::new().write(true).create_new(true).open(destination)
        .map_err(|e| e.to_string())?;
    file.write_all(&wav).and_then(|()| file.flush()).map_err(|e| e.to_string())?;
    eprintln!("full_scale_pa=0.005; clips={clips}; final_energy_j={:.10e}; source released without resets",
        renderer.source().system().stored_energy_j());
    Ok(())
}
fn main() {
    if let Err(error) = run() { eprintln!("two-sided reed refused: {error}"); std::process::exit(1); }
}
