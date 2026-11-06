//! Manufactured multidimensional transport through the native temperature API.
//! Spatial refinement and physical energy closure; no multidimensional DMP claim.

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, ConvectionScheme, EnergyConfig, FlowFace, FlowField, FluidProperties,
    SolidMaterial, ThermalFace, ThermalSetup, TransientConfig, Voxel, VoxelDomain,
    march_energy, solve_energy,
};

fn properties() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 2.0,
        conductivity_w_m_k: 0.02,
        kinematic_viscosity_m2_s: 0.1,
    }
}

fn energy(scheme: ConvectionScheme) -> EnergyConfig {
    EnergyConfig { scheme, tolerance: 1e-12, ..EnergyConfig::default() }
}

fn polynomial(x: f64) -> f64 {
    x * x * (1.0 - x) * (1.0 - x)
}

fn derivative(x: f64) -> f64 {
    2.0 * x * (1.0 - x) * (1.0 - 2.0 * x)
}

fn primitive(x: f64) -> f64 {
    x * x * x * (1.0 / 3.0 - 0.5 * x + x * x / 5.0)
}

fn manufactured(n: usize, direction: f64, scheme: ConvectionScheme) -> (Vec<f64>, f64) {
    let width = 1.0 / n as f64;
    let depth = 0.25;
    let domain = VoxelDomain::graded([
        vec![width; n], vec![width; n], vec![depth],
    ]).unwrap();
    let fluid = properties();
    let rho_cp = fluid.density_kg_m3 * fluid.specific_heat_j_kg_k;
    let velocity = [direction, 0.6 * direction, 0.0];
    let flow = FlowField::from_face_velocity(
        &domain,
        [FlowFace::Fixed, FlowFace::Fixed, FlowFace::Fixed, FlowFace::Fixed,
            FlowFace::Wall, FlowFace::Wall],
        |axis, _| velocity[axis],
    ).unwrap();
    assert!(flow.max_divergence(&domain) < 1e-15);
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[..4].fill(ThermalFace::Inflow { temperature: 300.0 });
    let mut thermal = ThermalSetup::new(faces);
    let amplitude = 2560.0; // The exact interior peak is 310 K.
    // T = 300 + amplitude * p(x) * p(y), p = x^2 (1-x)^2.
    // Integrate rho Cp u.grad(T) - k laplacian(T) analytically over each cell.
    // Both p and p' vanish at the sides, avoiding a first-order boundary
    // stencil dominating the interior reconstruction's refinement evidence.
    thermal.power_w = (0..domain.cell_count()).map(|cell| {
        let [i, j, _] = domain.coords(cell);
        let (xl, xr) = (domain.face_coord(0, i), domain.face_coord(0, i + 1));
        let (yl, yr) = (domain.face_coord(1, j), domain.face_coord(1, j + 1));
        let ix = primitive(xr) - primitive(xl);
        let iy = primitive(yr) - primitive(yl);
        amplitude * depth * (
            rho_cp * (velocity[0] * (polynomial(xr) - polynomial(xl)) * iy
                + velocity[1] * ix * (polynomial(yr) - polynomial(yl)))
            - fluid.conductivity_w_m_k * (
                (derivative(xr) - derivative(xl)) * iy
                    + ix * (derivative(yr) - derivative(yl)))
        )
    }).collect();
    let solution = solve_energy(
        &domain, &fluid, &[], &flow, &thermal, &energy(scheme), &CancelGate::new(),
    ).unwrap();
    let scale = thermal.power_w.iter().map(|q| q.abs()).sum::<f64>().max(1.0);
    assert!(solution.report.balance.residual_w.abs() < 2e-8 * scale);
    let error = solution.temperature.iter().enumerate().map(|(cell, &t)| {
        let [i, j, k] = domain.coords(cell);
        let [x, y, _] = domain.center(i, j, k);
        (t - (300.0 + amplitude * polynomial(x) * polynomial(y))).abs()
    }).sum::<f64>() / domain.cell_count() as f64;
    (solution.temperature, error)
}

#[test]
fn g1_g3_oblique_temperature_convection_refines_and_reverses() {
    let (coarse, coarse_error) = manufactured(12, 1.0, ConvectionScheme::EnthalpyMinmod);
    let (_, fine_error) = manufactured(24, 1.0, ConvectionScheme::EnthalpyMinmod);
    let (_, donor_error) = manufactured(24, 1.0, ConvectionScheme::Upwind);
    assert!(fine_error < 0.45 * coarse_error,
        "limited spatial refinement: coarse={coarse_error}, fine={fine_error}");
    assert!(fine_error < 0.4 * donor_error,
        "limited={fine_error}, donor={donor_error}");
    let (reverse, _) = manufactured(12, -1.0, ConvectionScheme::EnthalpyMinmod);
    // A half-turn of the square reverses both velocity components and the
    // manufactured advective source while preserving the thermal boundaries.
    for j in 0..12 {
        for i in 0..12 {
            let reflected = reverse[(11 - j) * 12 + 11 - i];
            assert!((coarse[j * 12 + i] - reflected).abs() < 2e-7);
        }
    }
}

fn circulation() -> (VoxelDomain, FlowField) {
    let domain = VoxelDomain::graded([
        vec![0.125; 8], vec![0.125; 8], vec![0.25],
    ]).unwrap();
    // Curl of psi = 4*x*(1-x)*y*(1-y): exact zero face divergence on this
    // binary grid, with no normal velocity at any enclosing wall.
    let flow = FlowField::from_face_velocity(
        &domain, [FlowFace::Wall; 6], |axis, [x, y, _]| match axis {
            0 => 4.0 * x * (1.0 - x) * (1.0 - 2.0 * y),
            1 => -4.0 * (1.0 - 2.0 * x) * y * (1.0 - y),
            _ => 0.0,
        },
    ).unwrap();
    assert!(flow.max_divergence(&domain) < 1e-16);
    (domain, flow)
}

#[test]
fn g0_g3_closed_circulation_preserves_equilibrium_and_physical_energy() {
    let (domain, flow) = circulation();
    let fluid = properties();
    let thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let run = |initial: &[f64], scheme| march_energy(
        &domain, &fluid, &[], &flow, &thermal, initial, |_| 1.0,
        &TransientConfig { time_step_s: 0.05, steps: 4, energy: energy(scheme) },
        None, &CancelGate::new(),
    ).unwrap();
    let uniform = run(&vec![317.0; domain.cell_count()], ConvectionScheme::EnthalpyMinmod);
    assert!(uniform.temperature.iter().all(|t| (t - 317.0).abs() < 2e-8));
    let initial: Vec<_> = (0..domain.cell_count()).map(|cell| {
        let [i, j, k] = domain.coords(cell);
        let [x, y, _] = domain.center(i, j, k);
        300.0 + 15.0 * x + 3.0 * y + 4.0 * x * y
    }).collect();
    let limited = run(&initial, ConvectionScheme::EnthalpyMinmod);
    let donor = run(&initial, ConvectionScheme::Upwind);
    assert!(limited.temperature.iter().zip(&donor.temperature)
        .any(|(a, b)| (a - b).abs() > 1e-6));
    let stored = limited.temperature.iter().zip(&initial).enumerate()
        .map(|(cell, (&t, &old))| fluid.density_kg_m3 * fluid.specific_heat_j_kg_k
            * domain.volume(cell) * (t - old))
        .sum::<f64>();
    assert!(stored.abs() < 2e-8, "closed-domain stored heat={stored}");
    let recorded = limited.records.iter()
        .map(|record| record.stored_energy_change_j).sum::<f64>();
    assert!((stored - recorded).abs() < 2e-9);
    for record in &limited.records {
        assert_eq!(record.source_j, 0.0);
        assert_eq!(record.boundary_outflow_j, 0.0);
        assert!(record.closure_j.abs() < 2e-8);
    }
}

#[test]
fn g4_native_limited_temperature_solve_and_march_obey_cancellation() {
    let (domain, flow) = circulation();
    let fluid = properties();
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(solve_energy(
        &domain, &fluid, &[], &flow,
        &ThermalSetup::new([ThermalFace::Temperature(317.0); 6]),
        &energy(ConvectionScheme::EnthalpyMinmod), &gate,
    ), Err(ChtError::Cancelled)));
    // Trip after transient setup, inside the authored power schedule. The
    // subsequent endpoint solve must refuse instead of publishing a step.
    let gate = CancelGate::new();
    assert!(matches!(march_energy(
        &domain, &fluid, &[], &flow,
        &ThermalSetup::new([ThermalFace::Adiabatic; 6]),
        &vec![317.0; domain.cell_count()], |_| { gate.request(); 1.0 },
        &TransientConfig {
            time_step_s: 0.05, steps: 1, energy: energy(ConvectionScheme::EnthalpyMinmod),
        },
        None, &gate,
    ), Err(ChtError::Cancelled)));
}

#[test]
fn g1_limited_reconstruction_preserves_resolved_wall_conductivity() {
    // The reference solid has k=1, but the supplied physical field owns k=10.
    // A reconstruction warm start must not switch the boundary back to k=1.
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let flow = FlowField::from_face_velocity(
        &domain, [FlowFace::Wall; 6], |_, _| 0.0,
    ).unwrap();
    let solids = [SolidMaterial::new("resolved-solid", 1.0)];
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(300.0);
    let mut thermal = ThermalSetup::new(faces);
    thermal.power_w = vec![10.0];
    for scheme in [ConvectionScheme::Upwind, ConvectionScheme::EnthalpyMinmod] {
        let solution = energy(scheme).solve_with_conductivity(
            &domain, &properties(), &solids, &flow, &thermal, &[[10.0; 3]],
            &CancelGate::new(),
        ).unwrap();
        // Q=10 W, A=1 m2, half-width=0.5 m: delta(T)=Q*0.5/(A*k).
        assert!((solution.temperature[0] - 300.5).abs() < 1e-9);
        assert!(solution.report.balance.residual_w.abs() < 1e-9);
    }
}
