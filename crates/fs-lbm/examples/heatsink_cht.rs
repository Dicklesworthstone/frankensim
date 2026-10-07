//! Ducted plate-fin heatsink, end to end through the conjugate pipeline:
//! voxelize aluminium fins on a base with a chip footprint, solve the steady
//! airflow (finite-volume SIMPLEC by default, or D3Q19 LBM on every core),
//! solve the conjugate energy equation, and report the junction
//! temperature, thermal resistance, effective film coefficient, pressure
//! drop and energy closure.
//!
//! ```text
//! cargo run --release -p fs-lbm --example heatsink_cht -- fv [inlet_velocity_m_s] [voxel_mm]
//! cargo run --release -p fs-lbm --example heatsink_cht -- lbm \
//!     [inlet_velocity_m_s] [lattice_inlet_velocity] [auto|bgk|central] [max_steps]
//! ```
//!
//! Geometry (metres, flow along +x): a 60 x 20 x 14 mm duct; 10 mm of
//! inlet run, a 30 mm long, 20 mm wide, 2 mm thick base carrying five 1 mm
//! fins 10 mm tall at a 4 mm pitch (2 mm bypass above the fin tips), and a
//! 20 mm outlet run. A 10 x 10 mm, 2 W chip is a uniform source in the
//! lowest base layer under the fin centre. Air at 300 K enters at 0.25 m/s;
//! the duct walls are adiabatic. Voxel edge 0.5 mm (the `fv` mode accepts
//! any edge that divides the 0.5 mm features, e.g. 1.0 or 0.25 mm).
//!
//! The output is one JSON object per line. Every value is Estimated
//! numerical evidence at this single resolution: no mesh-convergence,
//! turbulence, radiation, or physical-validation claim.

use fs_exec::CancelGate;
use fs_lbm::Face3;
use fs_lbm::conjugate::{
    EnergyConfig, FlowField, FluidProperties, FvBoundary, LbmCollisionChoice, LbmFlowConfig,
    SimpleConfig, SolidMaterial, ThermalFace, ThermalSetup, Voxel, VoxelDomain, lbm_duct_flow,
    simple_flow, solve_energy,
};

const MM: f64 = 1e-3;

fn heatsink(p: [f64; 3]) -> Voxel {
    let [x, y, z] = p;
    let in_sink_x = (10.0 * MM..40.0 * MM).contains(&x);
    let base = in_sink_x && z < 2.0 * MM;
    // Fin k occupies y in [1.5 + 4k, 2.5 + 4k] mm.
    let fin = in_sink_x
        && (2.0 * MM..12.0 * MM).contains(&z)
        && (0..5).any(|k| {
            let y0 = (1.5 + 4.0 * f64::from(k)) * MM;
            (y0..y0 + MM).contains(&y)
        });
    if base || fin {
        Voxel::Solid(0)
    } else {
        Voxel::Fluid
    }
}

/// Flow stage outputs shared by both solvers.
struct Flow {
    field: FlowField,
    inflow_m3_s: f64,
    pressure_drop_pa: f64,
}

fn lbm_flow(domain: &VoxelDomain, air: &FluidProperties, args: &[String]) -> Flow {
    let gate = CancelGate::new();
    let number = |index: usize, default: f64| {
        args.get(index)
            .map_or(default, |a| a.parse().expect("numeric argument"))
    };
    let collision = match args.get(2).map(String::as_str) {
        None | Some("auto") => LbmCollisionChoice::Auto,
        Some("bgk") => LbmCollisionChoice::Bgk,
        Some("central") => LbmCollisionChoice::CentralMoment,
        Some(other) => panic!("unknown collision {other}; use auto, bgk, or central"),
    };
    let config = LbmFlowConfig {
        inlet_velocity_m_s: number(0, 0.25),
        lattice_inlet_velocity: number(1, 0.05),
        collision,
        steady_tolerance: 1e-6,
        max_steps: args
            .get(3)
            .map_or(400_000, |a| a.parse().expect("integer max_steps")),
        ..LbmFlowConfig::default()
    };
    let started = std::time::Instant::now();
    let flow = lbm_duct_flow(domain, air, &config, &gate).expect("steady duct flow");
    let flow_s = started.elapsed().as_secs_f64();
    let r = &flow.report;
    println!(
        "{{\"stage\":\"flow\",\"solver\":\"lbm\",\"cells\":{},\"fluid_cells\":{},\"steps\":{},\"tau\":{:.5},\"collision\":\"{:?}\",\"inlet_mach\":{:.4},\"max_lattice_speed\":{:.4},\"realized_inflow_m3_s\":{:.6e},\"nominal_inflow_m3_s\":{:.6e},\"pressure_drop_pa\":{:.5},\"projection_max_correction_m3_s\":{:.3e},\"wall_s\":{flow_s:.1}}}",
        domain.cell_count(),
        domain.fluid_count(),
        r.steps,
        r.tau,
        r.collision,
        r.inlet_mach,
        r.max_lattice_speed,
        r.realized_inflow_m3_s,
        r.nominal_inflow_m3_s,
        r.pressure_drop_pa,
        r.projection.max_correction_m3_s,
    );
    Flow {
        inflow_m3_s: r.realized_inflow_m3_s,
        pressure_drop_pa: r.pressure_drop_pa,
        field: flow.field,
    }
}

fn fv_flow(domain: &VoxelDomain, air: &FluidProperties, args: &[String]) -> Flow {
    let gate = CancelGate::new();
    let inlet = args
        .first()
        .map_or(0.25, |a| a.parse().expect("numeric inlet velocity"));
    let mut config = SimpleConfig::new([
        FvBoundary::Inlet {
            velocity: [inlet, 0.0, 0.0],
        },
        FvBoundary::Outlet,
        FvBoundary::wall(),
        FvBoundary::wall(),
        FvBoundary::wall(),
        FvBoundary::wall(),
    ]);
    config.tolerance = 1e-6;
    let started = std::time::Instant::now();
    let flow = simple_flow(domain, air, &config, &gate).expect("steady duct flow");
    let flow_s = started.elapsed().as_secs_f64();
    let r = &flow.report;
    let [nx, _, _] = domain.dims();
    let pressure_drop_pa = flow.mean_pressure(domain, 0, 0).expect("fluid inlet layer")
        - flow
            .mean_pressure(domain, 0, nx - 1)
            .expect("fluid outlet layer");
    println!(
        "{{\"stage\":\"flow\",\"solver\":\"fv-simplec\",\"cells\":{},\"fluid_cells\":{},\"iterations\":{},\"mass_residual\":{:.2e},\"momentum_residual\":{:.2e},\"max_cell_reynolds\":{:.2},\"inflow_m3_s\":{:.6e},\"outflow_m3_s\":{:.6e},\"max_divergence_m3_s\":{:.2e},\"pressure_drop_pa\":{pressure_drop_pa:.5},\"wall_s\":{flow_s:.1}}}",
        domain.cell_count(),
        domain.fluid_count(),
        r.iterations,
        r.mass_residual,
        r.momentum_residual,
        r.max_cell_reynolds,
        r.inflow_m3_s,
        r.outflow_m3_s,
        r.max_divergence_m3_s,
    );
    Flow {
        inflow_m3_s: r.inflow_m3_s,
        pressure_drop_pa,
        field: flow.field,
    }
}

#[allow(clippy::too_many_lines)] // one linear report: flow, energy, result
fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let gate = CancelGate::new();
    let dx = match args.first().map(String::as_str) {
        Some("lbm") => 0.5 * MM,
        _ => {
            args.get(2)
                .map_or(0.5, |a| a.parse::<f64>().expect("numeric voxel_mm"))
                * MM
        }
    };
    let cells = |length_mm: f64| {
        let n = (length_mm * MM / dx).round();
        assert!(
            (n * dx - length_mm * MM).abs() < 1e-9,
            "voxel edge must divide {length_mm} mm"
        );
        n as usize
    };
    let (nx, ny, nz) = (cells(60.0), cells(20.0), cells(14.0));
    let domain = VoxelDomain::from_fn(nx, ny, nz, dx, heatsink).expect("admitted domain");
    let air = FluidProperties::dry_air_300k();
    let aluminium = [SolidMaterial::new(
        "AA6061-T6 (fixture card, 167 W/m/K)",
        167.0,
    )];
    let flow = match args.first().map(String::as_str) {
        None | Some("fv") => fv_flow(&domain, &air, args.get(1..).unwrap_or(&[])),
        Some("lbm") => lbm_flow(&domain, &air, &args[1..]),
        Some(other) => panic!("unknown flow solver {other}; use fv or lbm"),
    };

    let inlet_k = 300.0;
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[Face3::XMin as usize] = ThermalFace::Inflow {
        temperature: inlet_k,
    };
    faces[Face3::XMax as usize] = ThermalFace::Outflow {
        backflow_temperature: inlet_k,
    };
    let mut setup = ThermalSetup::new(faces);
    let chip_w = 2.0;
    let chip = setup.add_uniform_power(&domain, chip_w, |p| {
        p[2] < dx && (20.0 * MM..30.0 * MM).contains(&p[0]) && (5.0 * MM..15.0 * MM).contains(&p[1])
    });
    let started = std::time::Instant::now();
    let energy = solve_energy(
        &domain,
        &air,
        &aluminium,
        &flow.field,
        &setup,
        &EnergyConfig::default(),
        &gate,
    )
    .expect("steady conjugate energy");
    let energy_s = started.elapsed().as_secs_f64();
    let (_, junction_k) = energy
        .max_where(|c| !domain.is_fluid(c))
        .expect("solid cells exist");
    let solid_cells: Vec<usize> = (0..domain.cell_count())
        .filter(|&c| !domain.is_fluid(c))
        .collect();
    let mean_solid_k = solid_cells
        .iter()
        .map(|&c| energy.temperature[c])
        .sum::<f64>()
        / solid_cells.len() as f64;
    // Net upwind advective outflow is referenced to 0 K, which is exact
    // because the net boundary mass flow vanishes: it equals
    // rho c_p Q (T_out,bulk - T_in).
    let outlet_k = inlet_k
        + energy.report.balance.advective_outflow_w
            / (air.volumetric_heat_capacity() * flow.inflow_m3_s);
    // Wetted area: fluid/solid faces.
    let mut wetted = 0usize;
    for &c in &solid_cells {
        for f in 0..6 {
            let [x, y, z] = domain.coords(c);
            let neighbour = match f {
                0 => (x > 0).then(|| domain.index(x - 1, y, z)),
                1 => (x + 1 < nx).then(|| domain.index(x + 1, y, z)),
                2 => (y > 0).then(|| domain.index(x, y - 1, z)),
                3 => (y + 1 < ny).then(|| domain.index(x, y + 1, z)),
                4 => (z > 0).then(|| domain.index(x, y, z - 1)),
                _ => (z + 1 < nz).then(|| domain.index(x, y, z + 1)),
            };
            wetted += usize::from(neighbour.is_some_and(|n| domain.is_fluid(n)));
        }
    }
    let wetted_m2 = wetted as f64 * dx * dx;
    let interface_w = energy.solid_to_fluid_heat_w(&domain);
    let mean_air_k = 0.5 * (inlet_k + outlet_k);
    let h_eff = interface_w / (wetted_m2 * (mean_solid_k - mean_air_k));
    let b = energy.report.balance;
    println!(
        "{{\"stage\":\"energy\",\"chip_cells\":{chip},\"chip_w\":{chip_w},\"iterations\":{},\"relative_residual\":{:.2e},\"balance_relative_residual\":{:.2e},\"advective_outflow_w\":{:.6},\"interface_heat_w\":{:.6},\"max_cell_peclet\":{:.2},\"wall_s\":{energy_s:.1}}}",
        energy.report.iterations,
        energy.report.relative_residual,
        b.relative_residual,
        b.advective_outflow_w,
        interface_w,
        energy.report.max_cell_peclet,
    );
    println!(
        "{{\"stage\":\"result\",\"junction_k\":{junction_k:.4},\"mean_solid_k\":{mean_solid_k:.4},\"outlet_bulk_k\":{outlet_k:.4},\"thermal_resistance_k_per_w\":{:.4},\"wetted_area_m2\":{wetted_m2:.4e},\"effective_h_w_m2_k\":{h_eff:.3},\"pressure_drop_pa\":{:.5},\"evidence\":\"Estimated\"}}",
        (junction_k - inlet_k) / chip_w,
        flow.pressure_drop_pa,
    );
}
