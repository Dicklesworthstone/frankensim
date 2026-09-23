//! Geometry/material-driven valve and reciprocal tube, emitting physical CSV.
//! Usage: plate_aperture [YOUNG_PA DENSITY_KG_M3 THICKNESS_M]
//! The numeric sections and lay here are authored estimates, not measured cane.
use fs_couple::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture};
use fs_couple::bernoulli_aperture::plate::{PlateApertureOptions, PlateApertureReduction};
use fs_couple::bernoulli_aperture::tube::{ApertureTube, TubeDrive, UniformTubeSpec};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_plate::{AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateSection, SliceOptions};
use std::io::Write;

fn run() -> Result<(), String> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let physical = match args.as_slice() {
        [] => [4e9, 900.0, 0.0003],
        [e, rho, h] => {
            let mut values = [0.0; 3];
            for (slot, value) in values.iter_mut().zip([e, rho, h]) {
                *slot = value.parse::<f64>().map_err(|_| "expected numeric E, density and thickness".to_string())?;
            }
            values
        }
        _ => return Err("usage: plate_aperture [YOUNG_PA DENSITY_KG_M3 THICKNESS_M]".into()),
    };
    let section = PlateSection::isotropic(physical[0], 0.3, physical[2], physical[1]).map_err(|e| e.to_string())?;
    let mesh = PlateMesh::rectangle(0.025, 0.01, 6, 2);
    let root = mesh.nodes.iter().enumerate().filter(|(_, p)| p.0 == 0.0).map(|(i, _)| i).collect();
    let slit_edges: Vec<_> = mesh.boundary_edges().into_iter().filter(|&(a,b)|
        (mesh.nodes[a].0-0.025).abs() < 1e-12 && (mesh.nodes[b].0-0.025).abs() < 1e-12)
        .map(|(a,b)| [a,b]).collect();
    let tip = slit_edges[0][0];
    let chart = PlateChart::with_boundary_and_regions(mesh, section, root, vec![]).map_err(|e| e.to_string())?;
    let plate = PlateApertureReduction::from_chart(chart, PlateApertureOptions {
        assembly: AssemblyOptions { pretension: 0.0, support: EdgeSupport::Clamped },
        eigenvalue_window: (1.0, (core::f64::consts::TAU*450.0).powi(2)), mode_index: 0,
        eigensolver: SliceOptions::default(), max_nodes: 100, max_triangles: 100,
        slit_edges, rest_opening_m: 0.0002, damping_ratio: 0.02,
        max_slit_mode_variation: 0.25, max_slope: 0.1,
    }, &CancelGate::new()).map_err(|e| e.to_string())?;
    eprintln!("authored single-mode plate; mass_kg={:.17e}; stiffness_n_m={:.17e}; pressure_area_m2={:.17e}; closing_pressure_pa={:.17e}; slit_width_m={:.17e}",
        plate.mass_kg(), plate.stiffness_n_m(), plate.pressure_area_m2(), plate.closing_pressure_pa(), plate.width_m());
    // Ambient gas does not silently change the separately supplied solid card.
    let gas = GasState::try_new_moist_air(293.15, 101325.0, 0.0).map_err(|e| format!("{e:?}"))?;
    let dt = 1e-5;
    let tube = UniformTubeSpec {
        length_m: 0.25, radius_m: 0.007, sound_speed_m_s: gas.sound_speed,
        terminal_reflection: -0.8, max_length_error_m: 0.002, max_wave_memory_bytes: 1 << 20,
    };
    let z = tube.characteristic_impedance(gas.density).map_err(|e| e.to_string())?;
    let initial = ApertureState { opening_m: plate.options().rest_opening_m, opening_velocity_m_s: 0.0 };
    let lay = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "authored generalized plate-slit lay; no measured material claim".into())
        .and_then(|c| c.with_internal_loss(5.0)).map_err(|e| e.to_string())?;
    let valve = DynamicAperture::from_plate(plate, gas.density, z, dt, 4096, initial, lay)
        .map_err(|e| e.to_string())?;
    let mut model = ApertureTube::new(valve, tube).map_err(|e| e.to_string())?;
    eprintln!("requested_length_m={}; represented_length_m={}; pressure is internal, not an exterior microphone",
        tube.length_m, model.represented_length_m());
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    writeln!(out, "time_s,midpoint_bore_pressure_pa,opening_m,opening_velocity_m_s,tip_displacement_m,tip_velocity_m_s,max_slope,swept_flow_m3_s,jet_flow_m3_s,total_energy_j,upstream_work_j,loss_j,balance_residual_j")
        .map_err(|e| e.to_string())?;
    for n in 0..4096 {
        // Fixed physical drive for every material substitution: no compliance-
        // dependent pressure retuning or normalization hides a source change.
        let frame = model.step(TubeDrive { upstream_pressure_pa: if n < 1024 { 5.0 } else { 0.0 }, body_flow_m3_s: 0.0 })
            .map_err(|e| e.to_string())?;
        let plate = model.aperture().plate_reduction().expect("retained source");
        let a = frame.aperture;
        let motion = plate.nodal_motion(tip, a.state.opening_m, a.state.opening_velocity_m_s).map_err(|e| e.to_string())?;
        writeln!(out, "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            a.time_s, a.bore_pressure_pa, a.state.opening_m, a.state.opening_velocity_m_s,
            motion.displacement_rotation[0], motion.velocity_rotation_rate[0], plate.max_slope_at(a.state.opening_m),
            a.swept_flow_m3_s, a.jet_flow_m3_s, frame.stored_energy_j, frame.upstream_work_j,
            frame.dissipated_energy_j, frame.balance_residual_j()).map_err(|e| e.to_string())?;
    }
    out.flush().map_err(|e| e.to_string())
}
fn main() {
    if let Err(error) = run() { eprintln!("plate aperture refused: {error}"); std::process::exit(1); }
}
