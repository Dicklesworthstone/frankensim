//! Conjugate heat-transfer pipeline (bead frankensim-rc-root-q61wp.34).
//!
//! G1: an exact two-material conduction slab (harmonic face conductance is
//! the exact series resistance), and second-order convergence of the
//! fully developed parallel-plate Nusselt number toward 7.541.
//! G2: canonical laminar duct results — Nu_T = 7.541 (parallel plates,
//! isothermal walls), Nu_H = 8.235 (parallel plates, uniform flux through
//! conducting walls), Nu_T = 3.391 for a rectangular duct of aspect 0.5
//! (Shah & London 1978, Tables 41 and 43) — and an LBM-driven duct whose
//! downstream Nusselt number matches the analytic-profile solve.
//! G3: metamorphic wall-conductivity response.
//! Release-scale (ignored): simultaneously developing LBM duct at Pr = 0.72
//! against the fs-convection Shah-London Table 52 card, on two resolutions:
//! `cargo test -p fs-lbm --release --test conjugate -- --ignored --nocapture`.

use fs_exec::CancelGate;
use fs_lbm::Face3;
use fs_lbm::conjugate::{
    BuoyancyConfig, ChtError, ConvectionScheme, EnergyConfig, FlowFace, FlowField, FluidProperties,
    FvBoundary, FvBuoyancyConfig, LbmCollisionChoice, LbmFlowConfig, SimpleConfig, SolidMaterial,
    ThermalFace, ThermalSetup, TransientConfig, Voxel, VoxelDomain, fv_natural_convection,
    lbm_duct_flow, march_energy, natural_convection, simple_flow, solve_energy,
};

const OPEN_X: [FlowFace; 6] = [
    FlowFace::Fixed,
    FlowFace::Free,
    FlowFace::Wall,
    FlowFace::Wall,
    FlowFace::Wall,
    FlowFace::Wall,
];

fn unit_fluid() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0,
        kinematic_viscosity_m2_s: 1.0,
    }
}

#[test]
fn composite_slab_conduction_is_exact() {
    let gate = CancelGate::new();
    let dx = 0.01;
    let domain =
        VoxelDomain::from_fn(20, 1, 1, dx, |p| Voxel::Solid(u16::from(p[0] > 0.1))).unwrap();
    let solids = [
        SolidMaterial::new("k1", 1.0),
        SolidMaterial::new("k10", 10.0),
    ];
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(400.0);
    faces[1] = ThermalFace::Temperature(300.0);
    let solution = solve_energy(
        &domain,
        &unit_fluid(),
        &solids,
        &FlowField::quiescent(&domain),
        &ThermalSetup::new(faces),
        &EnergyConfig::default(),
        &gate,
    )
    .unwrap();
    let q = 100.0 / (0.1 / 1.0 + 0.1 / 10.0); // W/m^2
    for x in 0..20 {
        let xc = (x as f64 + 0.5) * dx;
        let exact = if xc < 0.1 {
            400.0 - q * xc
        } else {
            400.0 - q * 0.1 - q * (xc - 0.1) / 10.0
        };
        let t = solution.temperature[domain.index(x, 0, 0)];
        assert!((t - exact).abs() < 1e-9, "x={x}: {t} vs {exact}");
    }
    let balance = solution.report.balance;
    assert!(balance.boundary_outflow_w.abs() < 1e-12, "{balance:?}");
    assert!(balance.relative_residual < 1e-12, "{balance:?}");
}

/// Parallel plates of gap `h` resolved by `n` cells, analytic Poiseuille
/// profile with mean velocity giving `pe_dh = U D_h / alpha`.
fn plates(n: usize, pe_dh: f64, length_gaps: usize) -> (VoxelDomain, FlowField, f64) {
    let h = 1.0;
    let dx = h / n as f64;
    let domain = VoxelDomain::new(length_gaps * n, 1, n, dx).unwrap();
    let mean = pe_dh / (2.0 * h); // alpha = 1
    let flow = FlowField::from_face_velocity(&domain, OPEN_X, |axis, p| {
        if axis == 0 {
            let zeta = p[2] / h;
            6.0 * mean * zeta * (1.0 - zeta)
        } else {
            0.0
        }
    })
    .unwrap();
    (domain, flow, h)
}

fn isothermal_plate_nusselt(n: usize) -> (f64, f64) {
    let gate = CancelGate::new();
    let (domain, flow, h) = plates(n, 100.0, 12);
    let tw = 1.0;
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Inflow { temperature: 0.0 };
    faces[1] = ThermalFace::Outflow {
        backflow_temperature: 0.0,
    };
    faces[4] = ThermalFace::Temperature(tw);
    faces[5] = ThermalFace::Temperature(tw);
    let solution = solve_energy(
        &domain,
        &unit_fluid(),
        &[],
        &flow,
        &ThermalSetup::new(faces),
        &EnergyConfig::default(),
        &gate,
    )
    .unwrap();
    assert!(
        solution.report.balance.relative_residual < 1e-9,
        "{:?}",
        solution.report.balance
    );
    // Fully developed station: 8 gaps downstream (x* = x/(D_h Pe) = 0.04).
    let x = 8 * n;
    let dx = domain.dx();
    let bottom = solution.temperature[domain.index(x, 0, 0)];
    let top = solution.temperature[domain.index(x, 0, n - 1)];
    let q = 0.5 * ((tw - bottom) + (tw - top)) / (0.5 * dx);
    let tb = solution.bulk_temperature_x(&domain, &flow, x).unwrap();
    let nu = q * 2.0 * h / (tw - tb);
    (nu, solution.report.max_cell_peclet)
}

#[test]
fn parallel_plate_isothermal_nusselt_converges_to_7_541() {
    let (coarse, _) = isothermal_plate_nusselt(8);
    let (fine, peclet) = isothermal_plate_nusselt(16);
    let exact = 7.541;
    let (e8, e16) = ((coarse - exact).abs(), (fine - exact).abs());
    eprintln!("plates Nu_T: n=8 {coarse:.5} n=16 {fine:.5} (cell Pe {peclet:.2})");
    assert!(e16 / exact < 0.01, "n=16 Nu {fine} vs {exact}");
    assert!(
        e8 / e16 > 3.0,
        "observed order too low: errors {e8} -> {e16}"
    );
}

#[test]
fn conjugate_heated_walls_reach_uniform_flux_nusselt() {
    let gate = CancelGate::new();
    let n = 16;
    let wall = 2;
    let h = 1.0;
    let dx = h / n as f64;
    let nz = n + 2 * wall;
    let nx = 12 * n;
    let domain = VoxelDomain::from_fn(nx, 1, nz, dx, |p| {
        if p[2] < wall as f64 * dx || p[2] > (n + wall) as f64 * dx {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let mean = 100.0 / (2.0 * h);
    let flow = FlowField::from_face_velocity(&domain, OPEN_X, |axis, p| {
        let zeta = (p[2] - wall as f64 * dx) / h;
        if axis == 0 && (0.0..=1.0).contains(&zeta) {
            6.0 * mean * zeta * (1.0 - zeta)
        } else {
            0.0
        }
    })
    .unwrap();
    let flux = 2.0; // W/m^2 into each outer face
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Inflow { temperature: 0.0 };
    faces[1] = ThermalFace::Outflow {
        backflow_temperature: 0.0,
    };
    faces[4] = ThermalFace::HeatFlux(flux);
    faces[5] = ThermalFace::HeatFlux(flux);
    let solids = [SolidMaterial::new("wall", 10.0)];
    let solution = solve_energy(
        &domain,
        &unit_fluid(),
        &solids,
        &flow,
        &ThermalSetup::new(faces),
        &EnergyConfig::default(),
        &gate,
    )
    .unwrap();
    let balance = solution.report.balance;
    let input = 2.0 * flux * (nx as f64 * dx) * dx; // two faces, unit-cell depth
    // Flux faces count as negative outflow, so the net boundary flux closes.
    assert!(
        balance.boundary_outflow_w.abs() < 1e-9 * input,
        "{balance:?}"
    );
    assert!(balance.relative_residual < 1e-12, "{balance:?}");
    // Most heat leaves as enthalpy; the inlet face is a Dirichlet boundary
    // and conducts the remainder upstream.
    assert!(
        (balance.advective_outflow_w - input).abs() < 0.02 * input,
        "{balance:?}"
    );
    assert!((solution.solid_to_fluid_heat_w(&domain) - input).abs() < 1e-9 * input);
    // Fully developed station, interface temperature from the series
    // half-cell conductances on each side.
    let x = 8 * n;
    let mut nu = 0.0;
    for (solid_z, fluid_z) in [(wall - 1, wall), (wall + n, wall + n - 1)] {
        let ts = solution.temperature[domain.index(x, 0, solid_z)];
        let tf = solution.temperature[domain.index(x, 0, fluid_z)];
        let q = (ts - tf) / (0.5 * dx / 10.0 + 0.5 * dx / 1.0);
        let ti = tf + q * 0.5 * dx;
        let tb = solution.bulk_temperature_x(&domain, &flow, x).unwrap();
        nu += 0.5 * q * 2.0 * h / (ti - tb);
    }
    eprintln!("conjugate plates Nu_H = {nu:.5}");
    assert!((nu - 8.235).abs() / 8.235 < 0.01, "Nu_H {nu}");
}

#[test]
fn wall_conductivity_spreads_a_hot_spot() {
    // Metamorphic: raising only the wall conductivity must lower the peak
    // wall temperature of a localized heat source, while every watt still
    // leaves through the air.
    let run = |k_wall: f64| {
        let gate = CancelGate::new();
        let n = 8;
        let dx = 1.0 / n as f64;
        let domain = VoxelDomain::from_fn(96, 1, n + 2, dx, |p| {
            if p[2] < 2.0 * dx {
                Voxel::Solid(0)
            } else {
                Voxel::Fluid
            }
        })
        .unwrap();
        let flow = FlowField::from_face_velocity(&domain, OPEN_X, |axis, p| {
            let zeta = p[2] - 2.0 * dx;
            if axis == 0 && (0.0..=1.0).contains(&zeta) {
                300.0 * zeta * (1.0 - zeta)
            } else {
                0.0
            }
        })
        .unwrap();
        let mut faces = [ThermalFace::Adiabatic; 6];
        faces[0] = ThermalFace::Inflow { temperature: 300.0 };
        faces[1] = ThermalFace::Outflow {
            backflow_temperature: 300.0,
        };
        let mut setup = ThermalSetup::new(faces);
        let heated =
            setup.add_uniform_power(&domain, 5.0, |p| p[2] < dx && (3.0..4.0).contains(&p[0]));
        assert_eq!(heated, 8);
        let solution = solve_energy(
            &domain,
            &unit_fluid(),
            &[SolidMaterial::new("wall", k_wall)],
            &flow,
            &setup,
            &EnergyConfig::default(),
            &gate,
        )
        .unwrap();
        let balance = solution.report.balance;
        assert!((balance.source_w - 5.0).abs() < 1e-12);
        assert!(
            (balance.boundary_outflow_w - 5.0).abs() < 1e-7,
            "{balance:?}"
        );
        assert!(balance.relative_residual < 1e-10, "{balance:?}");
        // Discrete maximum principle: nothing is colder than the inflow.
        assert!(solution.temperature.iter().all(|&t| t >= 300.0 - 1e-9));
        solution.max_where(|c| !domain.is_fluid(c)).unwrap().1
    };
    let (soft, stiff) = (run(2.0), run(20.0));
    eprintln!("peak wall temperature: k=2 {soft:.4} K, k=20 {stiff:.4} K");
    assert!(
        stiff < soft,
        "k x10 must lower the hot spot: {soft} -> {stiff}"
    );
}

/// Fully developed rectangular-duct velocity shape from the same FV
/// operator: -lap(w) = 1, w = 0 on the walls, on an ny x nz cross-section.
fn duct_profile(ny: usize, nz: usize, dx: f64) -> Vec<f64> {
    let gate = CancelGate::new();
    let section = VoxelDomain::new(1, ny, nz, dx).unwrap();
    let mut faces = [ThermalFace::Temperature(0.0); 6];
    faces[0] = ThermalFace::Adiabatic;
    faces[1] = ThermalFace::Adiabatic;
    let mut setup = ThermalSetup::new(faces);
    setup.add_uniform_power(&section, (ny * nz) as f64 * dx * dx * dx, |_| true);
    let solution = solve_energy(
        &section,
        &unit_fluid(),
        &[],
        &FlowField::quiescent(&section),
        &setup,
        &EnergyConfig::default(),
        &gate,
    )
    .unwrap();
    let mean = solution.temperature.iter().sum::<f64>() / (ny * nz) as f64;
    solution.temperature.iter().map(|w| w / mean).collect()
}

fn rect_duct_nusselt(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    flow: &FlowField,
    station: usize,
    tw: f64,
) -> (f64, ChtReport) {
    let gate = CancelGate::new();
    let k = fluid.conductivity_w_m_k;
    let mut faces = [ThermalFace::Temperature(tw); 6];
    faces[0] = ThermalFace::Inflow { temperature: 0.0 };
    faces[1] = ThermalFace::Outflow {
        backflow_temperature: 0.0,
    };
    let solution = solve_energy(
        domain,
        fluid,
        &[],
        flow,
        &ThermalSetup::new(faces),
        &EnergyConfig::default(),
        &gate,
    )
    .unwrap();
    let [_, ny, nz] = domain.dims();
    let dx = domain.dx();
    let mut wall_heat = 0.0;
    for z in 0..nz {
        for y in 0..ny {
            let walls = usize::from(y == 0)
                + usize::from(y == ny - 1)
                + usize::from(z == 0)
                + usize::from(z == nz - 1);
            let t = solution.temperature[domain.index(station, y, z)];
            wall_heat += walls as f64 * 2.0 * dx * k * (tw - t); // W per layer
        }
    }
    let perimeter = 2.0 * (ny + nz) as f64 * dx;
    let q = wall_heat / (perimeter * dx);
    let dh = 4.0 * (ny * nz) as f64 * dx * dx / perimeter;
    let tb = solution.bulk_temperature_x(domain, flow, station).unwrap();
    (
        q * dh / (k * (tw - tb)),
        ChtReport {
            balance: solution.report.balance.relative_residual,
        },
    )
}

#[derive(Debug)]
struct ChtReport {
    balance: f64,
}

#[test]
fn rectangular_duct_aspect_half_reaches_shah_london_nu_t() {
    let (ny, nz, dx) = (12, 24, 1.0);
    let profile = duct_profile(ny, nz, dx);
    let perimeter = 2.0 * (ny + nz) as f64 * dx;
    let dh = 4.0 * (ny * nz) as f64 * dx * dx / perimeter;
    let mean = 50.0 / dh; // Pe_Dh = 50, alpha = 1; station x* = x/(D_h Pe) = 0.125
    let nx = 120;
    let domain = VoxelDomain::new(nx, ny, nz, dx).unwrap();
    let flow = FlowField::from_face_velocity(&domain, OPEN_X, |axis, p| {
        if axis == 0 {
            let (y, z) = (p[1] as usize, p[2] as usize);
            mean * profile[z * ny + y]
        } else {
            0.0
        }
    })
    .unwrap();
    assert!(flow.max_divergence(&domain) < 1e-9 * flow.max_flux());
    let (nu, report) = rect_duct_nusselt(&domain, &unit_fluid(), &flow, 100, 1.0);
    eprintln!("rectangular duct a=0.5: Nu_T = {nu:.5} (Shah-London 3.391)");
    assert!(report.balance < 1e-9, "{report:?}");
    assert!((nu - 3.391).abs() / 3.391 < 0.03, "Nu_T {nu}");
}

#[test]
fn lbm_duct_flow_drives_the_same_heat_transfer_as_the_analytic_profile() {
    let gate = CancelGate::new();
    // Aspect-0.5 duct, 8 x 16 cells, Re_Dh = 10, Pr = 0.7.
    let (nx, ny, nz) = (24, 8, 16);
    let dx = 1e-3;
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0 / 0.7,
        kinematic_viscosity_m2_s: 1.0,
    };
    let dh = 4.0 * (ny * nz) as f64 * dx * dx / (2.0 * (ny + nz) as f64 * dx);
    let mean = 10.0 * fluid.kinematic_viscosity_m2_s / dh;
    let domain = VoxelDomain::new(nx, ny, nz, dx).unwrap();
    let config = LbmFlowConfig {
        inlet_velocity_m_s: mean,
        lattice_inlet_velocity: 0.05,
        steady_tolerance: 1e-6,
        ..LbmFlowConfig::default()
    };
    let lbm = lbm_duct_flow(&domain, &fluid, &config, &gate).unwrap();
    let (lbm_flow, velocities, report) = (&lbm.field, &lbm.velocity_m_s, &lbm.report);
    eprintln!("{report:?}");
    assert!(report.tau > 0.6 && report.inlet_mach < 0.1);
    assert!(
        report.projection.max_divergence_after_m3_s < 1e-9 * report.projection.max_flux_m3_s,
        "{:?}",
        report.projection
    );
    // Mass: the projected field carries the inlet flow out unchanged. The
    // inlet rim cells share bounce-back links with the walls, so the admitted
    // flow sits slightly below the nominal plug.
    let inflow = -lbm_flow.boundary_outflow(&domain, Face3::XMin);
    let outflow = lbm_flow.boundary_outflow(&domain, Face3::XMax);
    let nominal = mean * (ny * nz) as f64 * dx * dx;
    eprintln!("inflow / nominal = {:.5}", inflow / nominal);
    for x in [1usize, 12, 22] {
        let mut sum = 0.0;
        for z in 0..nz {
            for y in 0..ny {
                sum += velocities[domain.index(x, y, z)][0];
            }
        }
        // Mass flux per layer is uniform to steady tolerance.
        let layer = sum * dx * dx;
        assert!(
            (layer - inflow).abs() < 2e-3 * inflow,
            "layer {x}: {layer} vs {inflow}"
        );
    }
    assert!((report.realized_inflow_m3_s - inflow).abs() < 1e-15 * inflow.max(1.0));
    assert!((inflow - nominal).abs() < 0.05 * nominal);
    assert!(
        (inflow - outflow).abs() < 1e-9 * inflow,
        "{inflow} vs {outflow}"
    );
    let mean = inflow / (ny * nz) as f64 / (dx * dx);
    // Developed pressure gradient against Shah & London's Darcy
    // f Re = 62.19 for aspect 0.5: dp/dx = f Re mu U / (2 D_h^2).
    let (p_a, p_b) = (
        lbm.mean_pressure_x(&domain, 10).unwrap(),
        lbm.mean_pressure_x(&domain, 20).unwrap(),
    );
    let gradient = (p_a - p_b) / (10.0 * dx);
    let darcy =
        62.19 * fluid.density_kg_m3 * fluid.kinematic_viscosity_m2_s * mean / (2.0 * dh * dh);
    eprintln!("developed dp/dx: lbm {gradient:.6e} Pa/m, Darcy fRe=62.19 {darcy:.6e} Pa/m");
    assert!(
        (gradient - darcy).abs() / darcy < 0.05,
        "{gradient} vs {darcy}"
    );
    assert!(report.pressure_drop_pa > 0.0);
    // Downstream, the LBM profile is the developed duct profile.
    let profile = duct_profile(ny, nz, dx);
    let station = 18;
    let mut worst = 0.0f64;
    for z in 0..nz {
        for y in 0..ny {
            let u = velocities[domain.index(station, y, z)][0];
            worst = worst.max((u / mean - profile[z * ny + y]).abs());
        }
    }
    eprintln!("max |u_lbm/U - w| at x={station}: {worst:.4}");
    assert!(worst < 0.06, "developed LBM profile deviates by {worst}");
    let analytic = FlowField::from_face_velocity(&domain, OPEN_X, |axis, p| {
        if axis == 0 {
            let (y, z) = ((p[1] / dx) as usize, (p[2] / dx) as usize);
            mean * profile[z * ny + y]
        } else {
            0.0
        }
    })
    .unwrap();
    // Same thermal problem on both flux fields: the downstream local
    // Nusselt numbers agree (both include the same axial conduction at
    // Pe = 14; the LBM field additionally carries its developing entrance).
    let station = 18;
    let (nu_lbm, lbm_report) = rect_duct_nusselt(&domain, &fluid, lbm_flow, station, 1.0);
    let (nu_analytic, _) = rect_duct_nusselt(&domain, &fluid, &analytic, station, 1.0);
    eprintln!("local Nu_T at x={station}: lbm {nu_lbm:.4} analytic {nu_analytic:.4}");
    assert!(lbm_report.balance < 1e-9, "{lbm_report:?}");
    assert!(
        (nu_lbm - nu_analytic).abs() / nu_analytic < 0.05,
        "{nu_lbm} vs {nu_analytic}"
    );
}

#[test]
fn refusals_are_structured() {
    let gate = CancelGate::new();
    let domain = VoxelDomain::new(8, 4, 4, 1.0).unwrap();
    let flow =
        FlowField::from_face_velocity(&domain, OPEN_X, |axis, _| if axis == 0 { 1.0 } else { 0.0 })
            .unwrap();
    // Flow leaves through a face declared adiabatic.
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Inflow { temperature: 0.0 };
    let err = solve_energy(
        &domain,
        &unit_fluid(),
        &[],
        &flow,
        &ThermalSetup::new(faces),
        &EnergyConfig::default(),
        &gate,
    )
    .unwrap_err();
    assert!(
        matches!(
            err,
            ChtError::FlowThroughClosedFace {
                face: Face3::XMax,
                ..
            }
        ),
        "{err}"
    );
    // Pure Neumann conduction has no anchor.
    let err = solve_energy(
        &domain,
        &unit_fluid(),
        &[],
        &FlowField::quiescent(&domain),
        &ThermalSetup::new([ThermalFace::HeatFlux(1.0); 6]),
        &EnergyConfig::default(),
        &gate,
    )
    .unwrap_err();
    assert!(
        matches!(
            err,
            ChtError::InvalidInput {
                field: "thermal.faces",
                ..
            }
        ),
        "{err}"
    );
    // A solid voxel naming an undeclared material.
    let mut bad = domain.clone();
    bad.set(1, 1, 1, Voxel::Solid(3));
    let err = solve_energy(
        &bad,
        &unit_fluid(),
        &[],
        &FlowField::quiescent(&bad),
        &ThermalSetup::new([ThermalFace::Temperature(1.0); 6]),
        &EnergyConfig::default(),
        &gate,
    )
    .unwrap_err();
    assert!(
        matches!(err, ChtError::UnknownMaterial { material: 3, .. }),
        "{err}"
    );
    // Air at 3 m/s on 1 mm voxels is under-resolved for the lattice.
    let err = lbm_duct_flow(
        &VoxelDomain::new(16, 8, 8, 1e-3).unwrap(),
        &FluidProperties::dry_air_300k(),
        &LbmFlowConfig {
            inlet_velocity_m_s: 3.0,
            ..LbmFlowConfig::default()
        },
        &gate,
    )
    .unwrap_err();
    assert!(matches!(err, ChtError::LatticeResolution { .. }), "{err}");
    // Non-tile dimensions refuse instead of panicking.
    let err = lbm_duct_flow(
        &VoxelDomain::new(10, 8, 8, 1e-3).unwrap(),
        &unit_fluid(),
        &LbmFlowConfig::default(),
        &gate,
    )
    .unwrap_err();
    assert!(matches!(err, ChtError::InvalidDomain { .. }), "{err}");
    // A tripped gate publishes nothing.
    let tripped = CancelGate::new();
    tripped.request();
    let err = solve_energy(
        &domain,
        &unit_fluid(),
        &[],
        &FlowField::quiescent(&domain),
        &ThermalSetup::new([ThermalFace::Temperature(1.0); 6]),
        &EnergyConfig {
            scheme: ConvectionScheme::Upwind,
            ..EnergyConfig::default()
        },
        &tripped,
    )
    .unwrap_err();
    assert_eq!(err, ChtError::Cancelled);
}

/// One simultaneously developing rung: uniform-inlet LBM flow through an
/// isothermal aspect-0.5 duct of `b x 2b` cells and length `nx`, Pr = 0.72.
/// Returns (Graetz number, channel-mean Nu_m,T, Reynolds number).
fn developing_duct_rung(b: usize, nx: usize, re_nominal: f64) -> (f64, f64, f64) {
    let gate = CancelGate::new();
    let (ny, nz) = (b, 2 * b);
    let dx = 1e-3;
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0 / 0.72,
        kinematic_viscosity_m2_s: 1.0,
    };
    let area = (ny * nz) as f64 * dx * dx;
    let dh = 4.0 * area / (2.0 * (ny + nz) as f64 * dx);
    let domain = VoxelDomain::new(nx, ny, nz, dx).unwrap();
    let config = LbmFlowConfig {
        inlet_velocity_m_s: re_nominal / dh,
        lattice_inlet_velocity: 0.05,
        steady_tolerance: 1e-7,
        collision: LbmCollisionChoice::Bgk,
        ..LbmFlowConfig::default()
    };
    let lbm = lbm_duct_flow(&domain, &fluid, &config, &gate).unwrap();
    let (flow, report) = (lbm.field, lbm.report);
    eprintln!(
        "rung b={b} nx={nx}: steps {} tau {:.4} {:?} realized/nominal {:.4}",
        report.steps,
        report.tau,
        report.collision,
        report.realized_inflow_m3_s / report.nominal_inflow_m3_s
    );
    let mut faces = [ThermalFace::Temperature(1.0); 6];
    faces[0] = ThermalFace::Inflow { temperature: 0.0 };
    faces[1] = ThermalFace::Outflow {
        backflow_temperature: 0.0,
    };
    let solution = solve_energy(
        &domain,
        &fluid,
        &[],
        &flow,
        &ThermalSetup::new(faces),
        &EnergyConfig::default(),
        &gate,
    )
    .unwrap();
    let q = report.realized_inflow_m3_s;
    let t_out = solution.report.balance.advective_outflow_w / q; // rho c = 1, T_in = 0
    let re = q / area * dh / fluid.kinematic_viscosity_m2_s;
    let gz = re * 0.72 * dh / (nx as f64 * dx);
    let nu_m = gz / 4.0 * fs_math::det::ln(1.0 / (1.0 - t_out));
    (gz, nu_m, re)
}

#[test]
#[ignore = "release-scale G2/Level-B: two LBM rungs vs Shah-London Table 52 (minutes)"]
fn developing_duct_matches_shah_london_table_52_card() {
    use fs_convection::{CorrelationId, CorrelationInputs, evaluate};
    let mut rows = Vec::new();
    for (b, nx) in [(12usize, 40usize), (16, 56)] {
        let (gz, nu, re) = developing_duct_rung(b, nx, 50.0);
        let length_ratio = re * 0.72 / gz;
        let card = evaluate(
            CorrelationId::RectangularDuctLaminarCwtDevelopingPr072,
            CorrelationInputs::forced(re, 0.72)
                .with_length_ratio(length_ratio)
                .with_aspect_ratio(0.5),
        )
        .unwrap();
        let reference = card.evidence().value;
        let discrepancy = (nu - reference) / reference;
        eprintln!(
            "Level-B rung b={b}: Re {re:.2} Gz {gz:.3} Nu_m lbm+fv {nu:.4} card {reference:.4} discrepancy {:+.2}%",
            100.0 * discrepancy
        );
        rows.push((nu, reference, discrepancy));
    }
    // The measured band, not a declared one: both rungs within 10% of the
    // card, and refinement does not move the answer away from it.
    for &(_, _, d) in &rows {
        assert!(d.abs() < 0.10, "discrepancy {d}");
    }
    assert!(rows[1].2.abs() <= rows[0].2.abs() + 0.01, "{rows:?}");
}

/// de Vahl Davis differentially heated square cavity (Int. J. Numer. Meth.
/// Fluids 3, 1983): hot wall x = 0, cold wall x = H, adiabatic floor and
/// ceiling, gravity along -y, Pr = 0.71, extruded along a periodic z.
/// Returns the mean hot-wall Nusselt number and the run report.
fn de_vahl_davis(
    n: usize,
    rayleigh: f64,
    tau: f64,
    tolerance: f64,
) -> (f64, fs_lbm::conjugate::NaturalConvectionReport) {
    let gate = CancelGate::new();
    let dx = 1.0;
    let h = n as f64 * dx;
    let alpha = 1e-3;
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: alpha,
        kinematic_viscosity_m2_s: 0.71 * alpha,
    };
    let g = 9.81;
    let beta = rayleigh * fluid.kinematic_viscosity_m2_s * alpha / (g * h * h * h);
    let domain = VoxelDomain::new(n, n, 4, dx).unwrap();
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(1.0);
    faces[1] = ThermalFace::Temperature(0.0);
    let config = BuoyancyConfig {
        gravity_m_s2: [0.0, -g, 0.0],
        expansion_per_k: beta,
        reference_temperature_k: 0.5,
        tau,
        periodic: [false, false, true],
        steady_tolerance: tolerance,
        ..BuoyancyConfig::default()
    };
    let run = natural_convection(
        &domain,
        &fluid,
        &[],
        &ThermalSetup::new(faces),
        &config,
        &gate,
    )
    .unwrap();
    let mut heat = 0.0;
    for z in 0..4 {
        for y in 0..n {
            let t = run.energy.temperature[domain.index(0, y, z)];
            heat += 2.0 * dx * alpha * (1.0 - t);
        }
    }
    let nu = heat / (alpha * 1.0 / h * h * 4.0 * dx);
    assert!(
        run.energy.report.balance.relative_residual < 1e-9,
        "{:?}",
        run.energy.report.balance
    );
    (nu, run.report)
}

#[test]
fn natural_convection_cavity_matches_de_vahl_davis_at_ra_1e3() {
    let (nu, report) = de_vahl_davis(16, 1e3, 0.8, 1e-6);
    eprintln!("de Vahl Davis Ra=1e3 n=16: Nu = {nu:.4} (reference 1.118); {report:?}");
    assert!(report.coupling_residual_k < 1e-3, "{report:?}");
    assert!((nu - 1.118).abs() / 1.118 < 0.03, "Nu {nu}");
}

#[test]
fn heated_block_enclosure_closes_energy_and_responds_to_buoyancy() {
    // A conducting block on the floor of a cold-walled enclosure: every watt
    // leaves through the walls, and stronger buoyancy (larger beta) must
    // cool the block (metamorphic).
    let run = |beta: f64| {
        let gate = CancelGate::new();
        let dx = 1.0;
        let domain = VoxelDomain::from_fn(16, 16, 4, dx, |p| {
            if (6.0..10.0).contains(&p[0]) && p[1] < 4.0 {
                Voxel::Solid(0)
            } else {
                Voxel::Fluid
            }
        })
        .unwrap();
        let alpha = 1e-3;
        let fluid = FluidProperties {
            density_kg_m3: 1.0,
            specific_heat_j_kg_k: 1.0,
            conductivity_w_m_k: alpha,
            kinematic_viscosity_m2_s: 0.71 * alpha,
        };
        let mut faces = [ThermalFace::Adiabatic; 6];
        faces[0] = ThermalFace::Temperature(0.0);
        faces[1] = ThermalFace::Temperature(0.0);
        faces[3] = ThermalFace::Temperature(0.0);
        let mut setup = ThermalSetup::new(faces);
        setup.add_uniform_power(&domain, 4e-3, |p| (6.0..10.0).contains(&p[0]) && p[1] < 1.0);
        let config = BuoyancyConfig {
            gravity_m_s2: [0.0, -9.81, 0.0],
            expansion_per_k: beta,
            reference_temperature_k: 0.0,
            tau: 0.8,
            periodic: [false, false, true],
            steady_tolerance: 1e-6,
            ..BuoyancyConfig::default()
        };
        let result = natural_convection(
            &domain,
            &fluid,
            &[SolidMaterial::new("block", 100.0 * alpha)],
            &setup,
            &config,
            &gate,
        )
        .unwrap();
        let balance = result.energy.report.balance;
        assert!(
            (balance.boundary_outflow_w - 4e-3).abs() < 1e-9,
            "{balance:?}"
        );
        assert!(
            result.energy.temperature.iter().all(|&t| t >= -1e-9),
            "maximum principle"
        );
        result.energy.max_where(|c| !domain.is_fluid(c)).unwrap().1
    };
    // beta for Rayleigh numbers of roughly 6e2 and 3e3 on the enclosure
    // height at the block's ~0.36 K temperature rise.
    let (weak, strong) = (run(3e-8), run(1.5e-7));
    eprintln!("block peak rise: weak buoyancy {weak:.5} K, strong {strong:.5} K");
    assert!(
        strong < weak,
        "stronger buoyancy must cool the block: {weak} -> {strong}"
    );
}

#[test]
#[ignore = "release-scale G2: de Vahl Davis Ra = 1e4 and 1e5 on 32 x 32 (minutes)"]
fn natural_convection_cavity_matches_de_vahl_davis_release() {
    // Lower viscosity at Ra = 1e5 keeps the buoyant lattice speed small.
    for (ra, reference, tau) in [(1e4, 2.243, 0.8), (1e5, 4.519, 0.56)] {
        let (nu, report) = de_vahl_davis(32, ra, tau, 1e-7);
        eprintln!(
            "de Vahl Davis Ra={ra:e} n=32: Nu = {nu:.4} (reference {reference}); steps {}",
            report.steps
        );
        assert!(
            (nu - reference).abs() / reference < 0.05,
            "Ra {ra}: Nu {nu}"
        );
    }
}

#[test]
fn transient_lumped_cube_follows_backward_euler_and_the_exponential() {
    // A near-isothermal solid (Biot ~ 4e-7) cooled through convective faces:
    // the discrete answer is the backward-Euler recursion of the lumped ODE
    // exactly, and the continuous exponential to O(dt).
    let gate = CancelGate::new();
    let dx = 0.01;
    let domain = VoxelDomain::from_fn(4, 4, 4, dx, |_| Voxel::Solid(0)).unwrap();
    let (h, k, rho_c) = (10.0, 1e6, 1e3);
    let solids = [SolidMaterial::new("block", k).with_heat_capacity(rho_c)];
    let setup = ThermalSetup::new([ThermalFace::Convective { h, ambient: 0.0 }; 6]);
    let face = dx * dx / (1.0 / h + 0.5 * dx / k);
    let conductance = 6.0 * 16.0 * face;
    let capacity = rho_c * 64.0 * dx * dx * dx;
    let tau = capacity / conductance;
    let dt = tau / 50.0;
    let steps = 100;
    let run = |dt: f64, steps: usize| {
        march_energy(
            &domain,
            &unit_fluid(),
            &solids,
            &FlowField::quiescent(&domain),
            &setup,
            &vec![1.0; domain.cell_count()],
            |_| 1.0,
            &TransientConfig {
                time_step_s: dt,
                steps,
                energy: EnergyConfig {
                    // k = 1e6 against a 10 W/m^2K film is conditioned near
                    // the double-precision floor (~5e-11 residual).
                    tolerance: 1e-10,
                    ..EnergyConfig::default()
                },
            },
            &gate,
        )
        .unwrap()
    };
    let coarse = run(dt, steps);
    let recursion = (1.0 / (1.0 + dt / tau)).powi(steps as i32);
    let t_end = coarse.temperature[0];
    println!(
        "lumped cube: T(2 tau) = {t_end:.10}, BE recursion {recursion:.10}, exp {:.10}",
        (-2.0f64).exp()
    );
    assert!((t_end - recursion).abs() < 1e-6 * recursion);
    for record in &coarse.records {
        assert!(record.closure_j.abs() < 1e-7 * record.stored_energy_change_j.abs());
    }
    // Halving dt halves the distance to the continuous solution.
    let fine = run(dt / 2.0, 2 * steps);
    let (e1, e2) = (
        (t_end - (-2.0f64).exp()).abs(),
        (fine.temperature[0] - (-2.0f64).exp()).abs(),
    );
    println!("time error: dt {e1:.3e}, dt/2 {e2:.3e}");
    assert!((e1 / e2 - 2.0).abs() < 0.1, "first-order ratio {}", e1 / e2);
}

#[test]
#[allow(clippy::too_many_lines)] // fixture, march, steady reference, refusal
fn transient_heated_channel_closes_energy_and_settles_to_the_steady_solution() {
    let gate = CancelGate::new();
    let n = 8;
    let dx = 1.0 / n as f64;
    let domain = VoxelDomain::from_fn(48, 1, n + 2, dx, |p| {
        if p[2] < 2.0 * dx {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let flow = FlowField::from_face_velocity(&domain, OPEN_X, |axis, p| {
        let zeta = p[2] - 2.0 * dx;
        if axis == 0 && (0.0..=1.0).contains(&zeta) {
            60.0 * zeta * (1.0 - zeta)
        } else {
            0.0
        }
    })
    .unwrap();
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Inflow { temperature: 0.0 };
    faces[1] = ThermalFace::Outflow {
        backflow_temperature: 0.0,
    };
    let mut setup = ThermalSetup::new(faces);
    setup.add_uniform_power(&domain, 2.0, |p| p[2] < dx && (1.0..3.0).contains(&p[0]));
    let solids = [SolidMaterial::new("plate", 20.0).with_heat_capacity(5.0)];
    let steady = solve_energy(
        &domain,
        &unit_fluid(),
        &solids,
        &flow,
        &setup,
        &EnergyConfig::default(),
        &gate,
    )
    .unwrap();
    let march = march_energy(
        &domain,
        &unit_fluid(),
        &solids,
        &flow,
        &setup,
        &vec![0.0; domain.cell_count()],
        |_| 1.0,
        &TransientConfig {
            time_step_s: 0.05,
            steps: 400,
            energy: EnergyConfig::default(),
        },
        &gate,
    )
    .unwrap();
    let stored: f64 = march.records.iter().map(|r| r.stored_energy_change_j).sum();
    let net: f64 = march
        .records
        .iter()
        .map(|r| r.source_j - r.boundary_outflow_j)
        .sum();
    let worst = march
        .records
        .iter()
        .map(|r| r.closure_j.abs())
        .fold(0.0, f64::max);
    let gap = march
        .temperature
        .iter()
        .zip(&steady.temperature)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    let peak = steady.temperature.iter().copied().fold(0.0, f64::max);
    println!(
        "channel: stored {stored:.9} J vs net {net:.9} J, worst step closure {worst:.2e} J, |T(20 s) - T_steady| = {gap:.2e} K of {peak:.4}"
    );
    assert!((stored - net).abs() < 1e-9 * stored.abs());
    assert!(gap < 1e-3 * peak, "not settled: {gap}");
    // Heating from cold is monotone in the peak solid temperature.
    assert!(
        march
            .records
            .windows(2)
            .all(|w| w[1].max_solid_temperature_k >= w[0].max_solid_temperature_k - 1e-12)
    );
    // A solid without heat capacity refuses.
    let err = march_energy(
        &domain,
        &unit_fluid(),
        &[SolidMaterial::new("plate", 20.0)],
        &flow,
        &setup,
        &vec![0.0; domain.cell_count()],
        |_| 1.0,
        &TransientConfig {
            time_step_s: 0.05,
            steps: 1,
            energy: EnergyConfig::default(),
        },
        &gate,
    )
    .unwrap_err();
    assert!(
        matches!(
            err,
            ChtError::InvalidInput {
                field: "solid.volumetric_heat_capacity_j_m3_k",
                ..
            }
        ),
        "{err}"
    );
}

/// Ghia, Ghia & Shin (1982), Table I, Re = 100: u on the vertical centreline.
const GHIA_RE100: [(f64, f64); 9] = [
    (0.0547, -0.03717),
    (0.1719, -0.10150),
    (0.2813, -0.15662),
    (0.4531, -0.21090),
    (0.5, -0.20581),
    (0.6172, -0.13641),
    (0.7344, 0.00332),
    (0.8516, 0.23151),
    (0.9531, 0.68717),
];

/// Lid-driven cavity at Re = 100 on n x n cells (one cell deep, symmetry in
/// z); returns the worst deviation from Ghia's centreline u.
fn lid_cavity(n: usize) -> (f64, fs_lbm::conjugate::SimpleReport) {
    let gate = CancelGate::new();
    let dx = 1.0 / n as f64;
    let domain = VoxelDomain::new(n, n, 1, dx).unwrap();
    let fluid = FluidProperties {
        kinematic_viscosity_m2_s: 0.01,
        ..unit_fluid()
    };
    let mut config = SimpleConfig::new([
        FvBoundary::wall(),
        FvBoundary::wall(),
        FvBoundary::wall(),
        FvBoundary::Wall {
            velocity: [1.0, 0.0, 0.0],
        },
        FvBoundary::Symmetry,
        FvBoundary::Symmetry,
    ]);
    config.tolerance = 1e-6;
    let flow = simple_flow(&domain, &fluid, &config, &gate).unwrap();
    // u on x = 1/2 (between columns n/2 - 1 and n/2), extended to the
    // walls (u = 0 at y = 0, u = 1 at the lid) for linear interpolation.
    let mut ys = vec![0.0];
    let mut us = vec![0.0];
    for y in 0..n {
        ys.push((y as f64 + 0.5) * dx);
        us.push(
            0.5 * (flow.velocity_m_s[domain.index(n / 2 - 1, y, 0)][0]
                + flow.velocity_m_s[domain.index(n / 2, y, 0)][0]),
        );
    }
    ys.push(1.0);
    us.push(1.0);
    let worst = GHIA_RE100
        .iter()
        .map(|&(y, u)| {
            let i = ys.iter().position(|&v| v >= y).unwrap().max(1);
            let t = (y - ys[i - 1]) / (ys[i] - ys[i - 1]);
            (us[i - 1] + t * (us[i] - us[i - 1]) - u).abs()
        })
        .fold(0.0, f64::max);
    (worst, flow.report)
}

#[test]
fn simplec_lid_driven_cavity_matches_ghia_at_re_100() {
    // G2: measured worst centreline deviations 0.0268 (16^2) and 0.0059
    // (32^2) against Ghia's 129^2 multigrid solution.
    let (coarse, _) = lid_cavity(16);
    let (fine, report) = lid_cavity(32);
    eprintln!("cavity worst |u - ghia|: 16^2 {coarse:.4} 32^2 {fine:.4} {report:?}");
    assert!(fine < 0.01, "32^2 deviates from Ghia by {fine}");
    assert!(
        fine < 0.5 * coarse,
        "no refinement gain: {coarse} -> {fine}"
    );
    // Closed cavity: the pinned pressure correction still conserves mass.
    assert!(report.max_divergence_m3_s < 1e-12 * report.max_cell_reynolds.max(1.0));
    assert_eq!(report.inflow_m3_s, 0.0);
}

/// Developed plane Poiseuille flow between walls `n` cells apart (Re_H = 1):
/// the developed pressure gradient.
fn poiseuille(n: usize) -> (f64, fs_lbm::conjugate::SimpleReport) {
    let gate = CancelGate::new();
    let dx = 1.0 / n as f64;
    let domain = VoxelDomain::new(6 * n, n, 1, dx).unwrap();
    let mut config = SimpleConfig::new([
        FvBoundary::Inlet {
            velocity: [1.0, 0.0, 0.0],
        },
        FvBoundary::Outlet,
        FvBoundary::wall(),
        FvBoundary::wall(),
        FvBoundary::Symmetry,
        FvBoundary::Symmetry,
    ]);
    config.tolerance = 1e-9;
    let flow = simple_flow(&domain, &unit_fluid(), &config, &gate).unwrap();
    let (a, b) = (3 * n, 5 * n);
    let gradient = (flow.mean_pressure(&domain, 0, a).unwrap()
        - flow.mean_pressure(&domain, 0, b).unwrap())
        / ((b - a) as f64 * dx);
    (gradient, flow.report)
}

#[test]
fn simplec_plane_poiseuille_matches_the_exact_discrete_solution() {
    // G1. The staggered stencil with half-cell wall diffusion is solved
    // exactly by u_j = (G / 2 mu) (y_j (H - y_j) + dx^2 / 4); its cell mean
    // (midpoint rule: H^2 / 6 + dx^2 / 12) gives the discrete law
    // G = 12 mu U / H^2 * n^2 / (n^2 + 2): second-order convergence to the
    // continuum 12 mu U / H^2 with a known constant.
    for n in [4usize, 8, 16] {
        let (gradient, report) = poiseuille(n);
        let n2 = (n * n) as f64;
        let discrete = 12.0 * n2 / (n2 + 2.0);
        eprintln!("poiseuille n={n}: dp/dx {gradient:.9} discrete {discrete:.9} {report:?}");
        assert!(
            (gradient - discrete).abs() < 1e-6 * discrete,
            "n={n}: {gradient} vs {discrete}"
        );
        assert!((report.outflow_m3_s - report.inflow_m3_s).abs() < 1e-12);
    }
}

/// Quarter of a square duct (symmetry on y- and z-min) of half-side `m`
/// cells at Re_Dh = 2: the Darcy friction constant f Re.
fn quarter_square_duct(m: usize) -> f64 {
    let gate = CancelGate::new();
    let dx = 1.0 / m as f64;
    let domain = VoxelDomain::new(8 * m, m, m, dx).unwrap();
    let mut config = SimpleConfig::new([
        FvBoundary::Inlet {
            velocity: [1.0, 0.0, 0.0],
        },
        FvBoundary::Outlet,
        FvBoundary::Symmetry,
        FvBoundary::wall(),
        FvBoundary::Symmetry,
        FvBoundary::wall(),
    ]);
    config.tolerance = 1e-8;
    let flow = simple_flow(&domain, &unit_fluid(), &config, &gate).unwrap();
    let (a, b) = (4 * m, 6 * m);
    let gradient = (flow.mean_pressure(&domain, 0, a).unwrap()
        - flow.mean_pressure(&domain, 0, b).unwrap())
        / ((b - a) as f64 * dx);
    // Full side 2 = D_h: f Re = 2 (dp/dx) D_h^2 / (mu U).
    8.0 * gradient
}

#[test]
fn simplec_square_duct_friction_converges_to_shah_london() {
    // G1/G2: Shah & London (1978) Darcy f Re = 56.91 for the square duct.
    // Measured 53.749 (m = 4) and 56.069 (m = 8): error ratio 3.8, and the
    // Richardson value 56.842 is within 0.12 %.
    let coarse = quarter_square_duct(4);
    let fine = quarter_square_duct(8);
    let richardson = fine + (fine - coarse) / 3.0;
    eprintln!("square duct fRe: {coarse:.4} {fine:.4} richardson {richardson:.4}");
    let reference = 56.91;
    let ratio = (reference - coarse) / (reference - fine);
    assert!((3.0..5.0).contains(&ratio), "order ratio {ratio}");
    assert!((fine - reference).abs() < 0.02 * reference, "{fine}");
    assert!(
        (richardson - reference).abs() < 0.005 * reference,
        "{richardson}"
    );
}

#[test]
fn simplec_flow_drives_the_conjugate_duct_like_the_analytic_profile() {
    // The same aspect-0.5 duct and thermal problem as the LBM rung above:
    // the FV flux field is exactly divergence-free, carries the inlet flow,
    // develops the Shah-London profile and gives the same downstream Nusselt
    // number as the analytic developed profile.
    let gate = CancelGate::new();
    let (nx, ny, nz) = (24, 8, 16);
    let dx = 1e-3;
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0 / 0.7,
        kinematic_viscosity_m2_s: 1.0,
    };
    let dh = 4.0 * (ny * nz) as f64 * dx * dx / (2.0 * (ny + nz) as f64 * dx);
    let mean = 10.0 * fluid.kinematic_viscosity_m2_s / dh;
    let domain = VoxelDomain::new(nx, ny, nz, dx).unwrap();
    let mut config = SimpleConfig::new([
        FvBoundary::Inlet {
            velocity: [mean, 0.0, 0.0],
        },
        FvBoundary::Outlet,
        FvBoundary::wall(),
        FvBoundary::wall(),
        FvBoundary::wall(),
        FvBoundary::wall(),
    ]);
    config.tolerance = 1e-8;
    let flow = simple_flow(&domain, &fluid, &config, &gate).unwrap();
    let report = &flow.report;
    eprintln!("{report:?}");
    let nominal = mean * (ny * nz) as f64 * dx * dx;
    assert!((report.inflow_m3_s - nominal).abs() < 1e-12 * nominal);
    assert!((report.outflow_m3_s - nominal).abs() < 1e-7 * nominal);
    assert!(report.max_divergence_m3_s < 1e-7 * nominal);
    let (p_a, p_b) = (
        flow.mean_pressure(&domain, 0, 10).unwrap(),
        flow.mean_pressure(&domain, 0, 20).unwrap(),
    );
    let gradient = (p_a - p_b) / (10.0 * dx);
    let darcy =
        62.19 * fluid.density_kg_m3 * fluid.kinematic_viscosity_m2_s * mean / (2.0 * dh * dh);
    eprintln!("developed dp/dx: fv {gradient:.6e} Pa/m, Darcy fRe=62.19 {darcy:.6e} Pa/m");
    // Measured -2.98 %: the second-order discretization error with 8 cells
    // across the short side (cf. the square-duct rungs above).
    assert!(
        (gradient - darcy).abs() / darcy < 0.04,
        "{gradient} vs {darcy}"
    );
    let profile = duct_profile(ny, nz, dx);
    let station = 18;
    let mut worst = 0.0f64;
    for z in 0..nz {
        for y in 0..ny {
            let u = flow.velocity_m_s[domain.index(station, y, z)][0];
            worst = worst.max((u / mean - profile[z * ny + y]).abs());
        }
    }
    eprintln!("max |u_fv/U - w| at x={station}: {worst:.4}");
    assert!(worst < 0.03, "developed FV profile deviates by {worst}");
    let analytic = FlowField::from_face_velocity(&domain, OPEN_X, |axis, p| {
        if axis == 0 {
            let (y, z) = ((p[1] / dx) as usize, (p[2] / dx) as usize);
            mean * profile[z * ny + y]
        } else {
            0.0
        }
    })
    .unwrap();
    let (nu_fv, fv_report) = rect_duct_nusselt(&domain, &fluid, &flow.field, station, 1.0);
    let (nu_analytic, _) = rect_duct_nusselt(&domain, &fluid, &analytic, station, 1.0);
    eprintln!("local Nu_T at x={station}: fv {nu_fv:.4} analytic {nu_analytic:.4}");
    assert!(fv_report.balance < 1e-9, "{fv_report:?}");
    assert!(
        (nu_fv - nu_analytic).abs() / nu_analytic < 0.03,
        "{nu_fv} vs {nu_analytic}"
    );
}

#[test]
fn simplec_refusals_are_structured() {
    let gate = CancelGate::new();
    let domain = VoxelDomain::new(6, 3, 3, 1.0).unwrap();
    let open = [
        FvBoundary::Inlet {
            velocity: [1.0, 0.0, 0.0],
        },
        FvBoundary::Outlet,
        FvBoundary::wall(),
        FvBoundary::wall(),
        FvBoundary::wall(),
        FvBoundary::wall(),
    ];
    let field_of = |config: &SimpleConfig| match simple_flow(&domain, &unit_fluid(), config, &gate)
    {
        Err(ChtError::InvalidInput { field, .. }) => field,
        other => panic!("expected an input refusal, got {other:?}"),
    };
    let mut outward = open;
    outward[0] = FvBoundary::Inlet {
        velocity: [-1.0, 0.0, 0.0],
    };
    assert_eq!(
        field_of(&SimpleConfig::new(outward)),
        "simple.inlet_velocity"
    );
    let mut piercing = open;
    piercing[2] = FvBoundary::Wall {
        velocity: [0.0, 0.5, 0.0],
    };
    assert_eq!(
        field_of(&SimpleConfig::new(piercing)),
        "simple.wall_velocity"
    );
    let mut relaxed = SimpleConfig::new(open);
    relaxed.velocity_relaxation = 1.0;
    assert_eq!(field_of(&relaxed), "simple.velocity_relaxation");
    // An exhausted budget is a refusal, not a silently unconverged field.
    let mut short = SimpleConfig::new(open);
    short.max_iterations = 1;
    assert!(matches!(
        simple_flow(&domain, &unit_fluid(), &short, &gate),
        Err(ChtError::FlowNotSteady { steps: 1, .. })
    ));
    let solid = VoxelDomain::from_fn(4, 4, 4, 1.0, |_| Voxel::Solid(0)).unwrap();
    assert!(matches!(
        simple_flow(&solid, &unit_fluid(), &SimpleConfig::new(open), &gate),
        Err(ChtError::InvalidDomain { .. })
    ));
    let tripped = CancelGate::new();
    tripped.request();
    assert_eq!(
        simple_flow(&domain, &unit_fluid(), &SimpleConfig::new(open), &tripped),
        Err(ChtError::Cancelled)
    );
}

#[test]
fn simplec_outlet_admits_undeveloped_outflow_at_one_fixed_point() {
    // G3 metamorphic: a block two cells upstream of the outlet sends a
    // recirculating, undeveloped stream through the outlet plane. The
    // ghost-pressure outlet keeps the discrete problem square, so SIMPLEC
    // converges (an extrapolated outlet stalls here) and the converged field
    // does not depend on the under-relaxation factor.
    let gate = CancelGate::new();
    let dx = 1e-3;
    let domain = VoxelDomain::from_fn(24, 10, 1, dx, |p| {
        let (x, y) = (p[0] / dx, p[1] / dx);
        if (18.0..20.0).contains(&x) && y < 6.0 {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let fluid = FluidProperties {
        kinematic_viscosity_m2_s: 1e-4,
        ..unit_fluid()
    };
    let solve = |alpha: f64| {
        let mut config = SimpleConfig::new([
            FvBoundary::Inlet {
                velocity: [0.5, 0.0, 0.0],
            },
            FvBoundary::Outlet,
            FvBoundary::wall(),
            FvBoundary::wall(),
            FvBoundary::Symmetry,
            FvBoundary::Symmetry,
        ]);
        config.velocity_relaxation = alpha;
        config.tolerance = 1e-10;
        simple_flow(&domain, &fluid, &config, &gate).unwrap()
    };
    let slow = solve(0.5);
    let fast = solve(0.8);
    eprintln!("{:?}\n{:?}", slow.report, fast.report);
    for flow in [&slow, &fast] {
        let r = &flow.report;
        assert!((r.outflow_m3_s - r.inflow_m3_s).abs() < 1e-12 * r.inflow_m3_s);
        assert!(r.max_divergence_m3_s < 1e-14 * r.inflow_m3_s);
    }
    // Reversed flow crosses the outlet plane: the outflow is undeveloped.
    let backflow = (0..10)
        .map(|y| slow.velocity_m_s[domain.index(23, y, 0)][0])
        .fold(f64::INFINITY, f64::min);
    eprintln!("min outlet-layer u {backflow:.4e}");
    let worst = (0..domain.cell_count())
        .map(|c| (slow.velocity_m_s[c][0] - fast.velocity_m_s[c][0]).abs())
        .fold(0.0, f64::max);
    eprintln!("relaxation 0.5 vs 0.8: max |du| {worst:.3e}");
    assert!(worst < 1e-7, "fixed point depends on relaxation: {worst}");
}

/// De Vahl Davis square cavity (hot x = 0, cold x = H, adiabatic y walls)
/// on n x n cells, one cell deep with symmetry in z, by FV-SIMPLEC
/// natural convection: the hot-wall mean Nusselt number.
fn fv_de_vahl_davis(
    n: usize,
    rayleigh: f64,
) -> (f64, fs_lbm::conjugate::FvNaturalConvectionReport) {
    let gate = CancelGate::new();
    let dx = 1.0 / n as f64;
    let alpha = 1e-2;
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: alpha,
        kinematic_viscosity_m2_s: 0.71 * alpha,
    };
    let g = 9.81;
    let beta = rayleigh * fluid.kinematic_viscosity_m2_s * alpha / g;
    let domain = VoxelDomain::new(n, n, 1, dx).unwrap();
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(1.0);
    faces[1] = ThermalFace::Temperature(0.0);
    let mut flow = SimpleConfig::new([
        FvBoundary::wall(),
        FvBoundary::wall(),
        FvBoundary::wall(),
        FvBoundary::wall(),
        FvBoundary::Symmetry,
        FvBoundary::Symmetry,
    ]);
    flow.tolerance = 1e-7;
    let config = FvBuoyancyConfig::new([0.0, -g, 0.0], beta, 0.5, flow);
    let run = fv_natural_convection(
        &domain,
        &fluid,
        &[],
        &ThermalSetup::new(faces),
        &config,
        &gate,
    )
    .unwrap();
    assert!(run.energy.report.balance.relative_residual < 1e-9);
    // Hot-wall heat over the conduction heat k dT / H through the same wall.
    let heat: f64 = (0..n)
        .map(|y| 2.0 * dx * alpha * (1.0 - run.energy.temperature[domain.index(0, y, 0)]))
        .sum();
    (heat / (alpha * dx), run.report)
}

#[test]
fn fv_natural_convection_cavity_matches_de_vahl_davis() {
    // G2: de Vahl Davis (1983) benchmark Nu = 1.118 (Ra 1e3), 2.243 (Ra 1e4).
    for (rayleigh, reference) in [(1e3, 1.118), (1e4, 2.243)] {
        let (nu, report) = fv_de_vahl_davis(32, rayleigh);
        eprintln!("FV de Vahl Davis Ra {rayleigh:e}: Nu {nu:.4} (ref {reference}) {report:?}");
        assert!(
            (nu - reference).abs() < 0.03 * reference,
            "Ra {rayleigh}: {nu} vs {reference}"
        );
    }
}

/// Open vertical channel (chimney) between isothermal plates `n` cells
/// apart and `n * aspect` tall, open at the bottom and top: the induced
/// volume flow per unit depth over the Elenbaas fully developed limit
/// g beta dT b^3 / (12 nu) of the discrete stencil.
fn fv_chimney(n: usize, aspect: usize) -> (f64, f64, fs_lbm::conjugate::FvNaturalConvectionReport) {
    let gate = CancelGate::new();
    let dx = 1e-3;
    let b = n as f64 * dx;
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1e-5,
        kinematic_viscosity_m2_s: 0.71e-5,
    };
    let alpha = fluid.conductivity_w_m_k;
    // Channel Rayleigh number on the gap: g beta dT b^3 / (nu alpha) = 20.
    let g = 9.81;
    let beta = 20.0 * fluid.kinematic_viscosity_m2_s * alpha / (g * b * b * b);
    let domain = VoxelDomain::new(n, n * aspect, 1, dx).unwrap();
    let mut thermal = [ThermalFace::Adiabatic; 6];
    thermal[0] = ThermalFace::Temperature(1.0);
    thermal[1] = ThermalFace::Temperature(1.0);
    thermal[2] = ThermalFace::Outflow {
        backflow_temperature: 0.0,
    };
    thermal[3] = ThermalFace::Outflow {
        backflow_temperature: 0.0,
    };
    let mut flow = SimpleConfig::new([
        FvBoundary::wall(),
        FvBoundary::wall(),
        FvBoundary::Outlet,
        FvBoundary::Outlet,
        FvBoundary::Symmetry,
        FvBoundary::Symmetry,
    ]);
    flow.tolerance = 1e-8;
    let config = FvBuoyancyConfig::new([0.0, -g, 0.0], beta, 0.0, flow);
    let run = fv_natural_convection(
        &domain,
        &fluid,
        &[],
        &ThermalSetup::new(thermal),
        &config,
        &gate,
    )
    .unwrap();
    let flow_report = &run.flow.report;
    // The exact discrete Poiseuille law of this stencil (see the plane
    // Poiseuille test) carries (n^2 + 2) / n^2 more flow at a given forcing.
    let n2 = (n * n) as f64;
    let limit =
        g * beta * 1.0 * b * b * b / (12.0 * fluid.kinematic_viscosity_m2_s) * dx * (n2 + 2.0) / n2;
    let q = flow_report.outflow_m3_s;
    // Energy: wall heat leaves as advection through the top opening.
    let balance = &run.energy.report.balance;
    assert!(balance.relative_residual < 1e-9, "{balance:?}");
    assert!((flow_report.outflow_m3_s - flow_report.inflow_m3_s).abs() < 1e-12 * q.max(1e-30));
    (q / limit, flow_report.max_divergence_m3_s / q, run.report)
}

#[test]
fn fv_open_chimney_approaches_the_elenbaas_developed_limit() {
    // G1/G2 for open-boundary buoyancy: at Ra_b = 20 the fluid reaches the
    // wall temperature within a few gaps, so the induced flow tends to the
    // fully developed limit g beta dT b^3 / (12 nu) from below as the
    // channel lengthens (the cool entrance region carries less buoyancy).
    // Measured 0.9795 (L/b = 10) and 0.9931 (L/b = 30) of the discrete
    // limit, which itself is 1.031 x the continuum limit at 8 cells.
    let (short, _, short_report) = fv_chimney(8, 10);
    let (long, divergence, long_report) = fv_chimney(8, 30);
    eprintln!(
        "chimney Q/Q_limit: L/b=10 {short:.4} L/b=30 {long:.4}\n{short_report:?}\n{long_report:?}"
    );
    assert!(divergence < 1e-10);
    assert!(short < long && long < 1.0, "{short} {long}");
    assert!(
        long > 0.98,
        "long channel far from the developed limit: {long}"
    );
}
