//! Vreman integration through native MAC gradients, momentum and energy.
use super::*;
use crate::conjugate::{
    EnergyConfig, FluidProperties, SimpleConfig, ThermalFace, ThermalSetup, TimeScheme,
    TransientConfig, UnsteadyConfig, VoxelDomain, march_conjugate_with_sgs, march_energy,
    simple_unsteady, simple_unsteady_with_sgs,
};

fn fluid() -> FluidProperties {
    FluidProperties { density_kg_m3: 1.0, specific_heat_j_kg_k: 2.0,
        conductivity_w_m_k: 0.002, kinematic_viscosity_m2_s: 0.01 }
}

fn config(cv: f64) -> FvSgsConfig {
    FvSgsConfig { model: SgsModel::Vreman { cv }, filter_width_factor: 1.0,
        maximum_eddy_to_molecular: 100.0, usage: SgsUse::Exploration }
}

fn fill_faces(solver: &mut Solver<'_>, field: impl Fn([f64; 3]) -> [f64; 3]) {
    for a in 0..3 {
        for f in 0..solver.vel[a].len() {
            let at = solver.comps[a].coords(f);
            let point = [0, 1, 2].map(|b| {
                let low = solver.domain.face_coord(b, at[b]);
                if a == b { low } else { low + 0.5 * solver.domain.width(b, at[b]) }
            });
            solver.vel[a][f] = field(point)[a];
        }
    }
}

#[test]
fn graded_mac_strain_uses_each_width_and_reuses_velocity_storage() {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::graded([
        vec![0.125, 0.25, 0.5], vec![0.0625, 0.125, 0.25], vec![0.25, 0.5, 1.0],
    ]).unwrap();
    let controls = SimpleConfig::new([FvBoundary::Outlet; 6]);
    let fluid = fluid();
    let mut solver = Solver::new(&domain, &fluid, &controls);
    solver.start_time(0.01, TimeScheme::Bdf2);
    solver.enable_sgs(config(0.07), &gate).unwrap();
    fill_faces(&mut solver, |p| [2.0 * p[0], -3.0 * p[1], p[2]]);
    let allocation = solver.sgs.as_ref().unwrap().velocity.as_ptr();
    assert!(solver.sgs.as_ref().unwrap().filter_width.is_empty());
    for _ in 0..3 {
        solver.update_sgs(&gate).unwrap();
        assert_eq!(solver.sgs.as_ref().unwrap().velocity.as_ptr(), allocation);
        for c in 0..domain.cell_count() {
            let d = domain.widths(c);
            let beta = [4.0 * d[0] * d[0], 9.0 * d[1] * d[1], d[2] * d[2]];
            let exact = 0.07 * ((beta[0] * beta[1] + beta[0] * beta[2]
                + beta[1] * beta[2]) / 14.0).sqrt();
            assert!((solver.eddy_viscosity[c] / exact - 1.0).abs() < 1e-13);
        }
    }
    let before = solver.eddy_viscosity.clone();
    solver.enable_sgs(FvSgsConfig { filter_width_factor: 0.5, ..config(0.07) }, &gate).unwrap();
    solver.update_sgs(&gate).unwrap();
    for (coarse, fine) in before.iter().zip(&solver.eddy_viscosity) {
        assert!((fine / coarse - 0.25).abs() < 1e-14);
    }
    assert_eq!(solver.heat_eddy_viscosity(), solver.eddy_viscosity);
    let flow = solver.finish_projected(&fluid, 0, (0.0, 0.0));
    for (k, nu) in flow.eddy_conductivity(&fluid).iter().zip(&flow.eddy_viscosity_m2_s) {
        assert!((k / (2.0 * nu / 0.9) - 1.0).abs() < 1e-14);
    }
}

#[test]
fn laminar_shear_vanishes_but_multidirectional_transport_respects_caps() {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::graded([
        vec![0.1, 0.2, 0.3], vec![0.05, 0.1, 0.15], vec![0.2, 0.3, 0.4],
    ]).unwrap();
    let controls = SimpleConfig::new([FvBoundary::Outlet; 6]);
    let mut solver = Solver::new(&domain, &fluid(), &controls);
    solver.start_time(0.01, TimeScheme::Bdf2);
    solver.enable_sgs(config(0.07), &gate).unwrap();
    fill_faces(&mut solver, |p| [p[1], 0.0, 0.0]);
    solver.update_sgs(&gate).unwrap();
    assert!(solver.eddy_viscosity.iter().all(|&v| v == 0.0));
    fill_faces(&mut solver, |p| [2.0 * p[0], -3.0 * p[1], p[2]]);
    solver.enable_sgs(FvSgsConfig { maximum_eddy_to_molecular: 0.0, ..config(0.07) }, &gate).unwrap();
    assert!(solver.update_sgs(&gate).is_err(), "ratio gate must refuse, not clip");
    assert!(solver.enable_sgs(FvSgsConfig { usage: SgsUse::AcceptanceRender, ..config(0.07) }, &gate).is_err());
    gate.request();
    assert!(matches!(solver.enable_sgs(config(0.07), &gate), Err(ChtError::Cancelled)));
}

#[test]
fn live_vreman_momentum_and_heat_share_the_same_transport_and_zero_is_exact() {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::new(5, 5, 1, 0.1).unwrap();
    let fluid = fluid();
    let mut faces = [FvBoundary::wall(); 6];
    faces[3] = FvBoundary::Wall { velocity: [0.1, 0.0, 0.0] };
    faces[4] = FvBoundary::Symmetry;
    faces[5] = FvBoundary::Symmetry;
    let mut controls = SimpleConfig::new(faces);
    controls.pressure_tolerance = 1e-11;
    let mut unsteady = UnsteadyConfig::new(0.01, 1);
    unsteady.inner_iterations = 400;
    unsteady.inner_tolerance = 1e-9;
    let laminar = simple_unsteady(&domain, &fluid, &controls, &unsteady, |_| 1.0, &gate).unwrap();
    let disabled = simple_unsteady_with_sgs(&domain, &fluid, &controls, &unsteady,
        config(0.0), |_| 1.0, &gate).unwrap();
    assert_eq!(disabled, laminar);
    let initial: Vec<_> = (0..domain.cell_count())
        .map(|c| 300.0 + 10.0 * domain.coords(c)[1] as f64).collect();
    let setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let coupled = march_conjugate_with_sgs(&domain, &fluid, &[], &setup, &controls,
        &unsteady, config(0.07), None, &initial, |_| 1.0, |_| 1.0,
        &EnergyConfig::default(), None, &gate).unwrap();
    assert!(coupled.flow.flow.eddy_viscosity_m2_s.iter().any(|&nu| nu > 0.0));
    assert_ne!(coupled.flow.flow.field, laminar.flow.field, "closure must reach momentum");
    let mut explicit = setup.clone();
    explicit.eddy_conductivity_w_m_k = coupled.flow.flow.eddy_conductivity(&fluid);
    let thermal = TransientConfig { time_step_s: unsteady.time_step_s,
        steps: 1, energy: EnergyConfig::default() };
    let frozen = march_energy(&domain, &fluid, &[], &coupled.flow.flow.field,
        &explicit, &initial, |_| 1.0, &thermal, None, &gate).unwrap();
    let molecular = march_energy(&domain, &fluid, &[], &coupled.flow.flow.field,
        &setup, &initial, |_| 1.0, &thermal, None, &gate).unwrap();
    let difference = |other: &[f64]| coupled.temperature.iter().zip(other)
        .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
    assert!(difference(&frozen.temperature) < 1e-8);
    assert!(difference(&molecular.temperature) > 1e-8, "closure must reach heat");
    assert!(coupled.energy_records[0].closure_j.abs() < 1e-8);
}
