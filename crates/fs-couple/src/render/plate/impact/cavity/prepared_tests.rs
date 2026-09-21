use super::*;
use crate::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use crate::render::plate::impact::{ImpactConfig, VolumeSpring};
use crate::render::schedule::force::coupled::{ModalCouplingConfig,
    contact::{ModalContactConfig, multiple::MultiContactConfig}};
use crate::vibroacoustic::CavityModes;

fn config(rate: u32) -> LinearImpactConfig {
    LinearImpactConfig { sample_rate_hz: rate, max_steps: 20_000,
        maximum_generalized_force: 10_000.0,
        component: ModalAcousticTimeBudget { maximum_total_energy_j: 20.0,
            ..ModalAcousticTimeBudget::audible_reference() },
        coupling: ModalCouplingConfig { max_modes: 256, max_connections: 8,
            max_setup_terms: 100_000, nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 20.0, maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e5, solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-10, energy_relative_tolerance: 1e-7 },
        contact: ModalContactConfig { max_iterations: 100, maximum_force_n: 1e4,
            maximum_penetration_m: 0.01, force_absolute_tolerance_n: 1e-9,
            force_relative_tolerance: 1e-9 },
        multiple: MultiContactConfig { max_contacts: 512, max_sweeps: 100,
            max_setup_terms: 50_000_000 } }
}
fn body(w: f64, q: f64, v: f64) -> ImpactBody {
    ImpactBody { potential: BodyPotential::Linear(vec![w]), damping_per_s: vec![0.0],
        initial: vec![ModalAcousticState { displacement_m_sqrt_kg: q,
            velocity_m_sqrt_kg_per_s: v }] }
}
fn basis() -> CavityModes {
    CavityModes { omegas: vec![0.0, 1000.0], lambdas: vec![0.01; 2],
        interface: vec![vec![1.0]; 2], loss_factor: 0.0, rho0: 1.2, c0: 343.0 }
}

#[test]
fn uniform_prepared_cavity_is_the_original_volume_without_extra_coordinates() {
    let mut air = basis(); air.omegas.truncate(1); air.lambdas.truncate(1); air.interface.truncate(1);
    let gate = CancelGate::new_clock_free(); let cfg = config(500_000);
    let parts = || vec![body(800.0, 1e-5, 0.0), body(900.0, 0.0, 0.0)];
    let volume = VolumeConnection { spring: VolumeSpring {
        bulk_modulus_pa: air.rho0*air.c0*air.c0, volume_m3: 0.01,
        areas: vec![0.1, -0.1] }, reference_area_m2: 0.1 };
    let mut original = LinearImpactSystem::new(parts(), vec![], vec![volume], cfg, &gate).unwrap();
    let (mut actual, probe) = CavityCoupling::new(&air, 2, &[0.1, -0.1], &[0.0]).unwrap()
        .build_linear(parts(), vec![], 0.1, cfg, &gate).unwrap();
    assert_eq!(actual.mode_count(), 2); assert_eq!(probe.total_modes(), 2);
    for _ in 0..400 {
        let a = original.step(&[0.0; 2], &gate).unwrap();
        let b = actual.step(&[0.0; 2], &gate).unwrap();
        for (x, y) in original.state().iter().zip(actual.state()) { assert!((x-y).abs() < 1e-12); }
        assert!((a.stored_energy_j-b.stored_energy_j).abs() < 1e-12);
        let expected = -(air.rho0*air.c0*air.c0/0.01)
            * (0.1*actual.state()[0]-0.1*actual.state()[2]);
        assert!((probe.pressure_at(actual.state(), &[1.0]).unwrap()-expected).abs() < 1e-10);
    }
}

#[test]
fn three_body_air_port_converges_to_independent_coupled_eigen_dynamics() {
    let air = basis(); let overlaps = [0.05, 0.03, -0.04, 0.02];
    let a = air.rho0*air.c0*air.c0/0.01;
    // Independently assemble the full coupled Hessian, NOT a pairwise spring
    // approximation. The air cross term couples both heads to one inertia.
    let mut k = vec![0.0; 9]; k[0] = 800.0_f64.powi(2); k[4] = 900.0_f64.powi(2);
    for column in [[0.05, -0.04, 0.0], [0.03, 0.02, 1000.0/a.sqrt()]] {
        for i in 0..3 { for j in 0..3 { k[i*3+j] += a*column[i]*column[j]; } }
    }
    let eigen = fs_modal::eigh_gen_dense(&k, &[1.,0.,0.,0.,1.,0.,0.,0.,1.], 3).unwrap();
    let time = 0.004; let mut exact = [0.0; 6];
    for pair in eigen {
        let w = pair.lambda.sqrt(); let amplitude = pair.phi[0]*1e-5;
        for i in 0..3 {
            exact[2*i] += amplitude*pair.phi[i]*(w*time).cos();
            exact[2*i+1] -= amplitude*pair.phi[i]*w*(w*time).sin();
        }
    }
    let run = |rate| {
        let gate = CancelGate::new_clock_free();
        let (mut system, _) = CavityCoupling::new(&air, 2, &overlaps, &[0.0; 2]).unwrap()
            .build_linear(vec![body(800.0, 1e-5, 0.0), body(900.0, 0.0, 0.0)],
                vec![], 0.1, config(rate), &gate).unwrap();
        let energy = system.frame().stored_energy_j;
        for _ in 0..(time*f64::from(rate)).round() as usize { system.step(&[0.0; 3], &gate).unwrap(); }
        assert!(system.state()[3].abs() > 1e-6 && system.state()[5].abs() > 1e-6);
        assert!((system.frame().stored_energy_j-energy).abs() < 1e-10);
        system.state().iter().zip(exact).enumerate().map(|(i, (x, y))|
            (x-y).abs()*if i%2 == 0 {1000.0} else {1.0}).fold(0.0, f64::max)
    };
    let coarse = run(20_000); let fine = run(40_000);
    assert!(coarse > 1e-10 && fine < 0.35*coarse, "second-order refinement: {coarse:e} -> {fine:e}");
}

#[test]
fn rescaled_pressure_basis_keeps_the_same_prepared_motion_and_pressure() {
    let air = basis(); let mut scaled = air.clone();
    for norm in &mut scaled.lambdas { *norm *= 9.0; }
    for row in &mut scaled.interface { for value in row { *value *= 3.0; } }
    let gate = CancelGate::new_clock_free();
    let build = |air: &CavityModes, overlap: &[f64]| {
        CavityCoupling::new(air, 2, overlap, &[0.0; 2]).unwrap()
            .build_linear(vec![body(800.0, 1e-5, 0.02), body(900.0, 0.0, 0.0)],
                vec![], 0.1, config(500_000), &gate).unwrap()
    };
    let (mut a, pa) = build(&air, &[0.1, 0.07, -0.1, 0.04]);
    let (mut b, pb) = build(&scaled, &[0.3, 0.21, -0.3, 0.12]);
    for _ in 0..200 {
        a.step(&[0.0; 3], &gate).unwrap(); b.step(&[0.0; 3], &gate).unwrap();
        for (x, y) in a.state().iter().zip(b.state()) { assert!((x-y).abs() < 1e-10); }
        assert!((pa.pressure_at(a.state(), &[1.0, 0.4]).unwrap()
            -pb.pressure_at(b.state(), &[3.0, 1.2]).unwrap()).abs() < 1e-6);
    }
}

#[test]
fn many_distinct_wire_coordinates_contact_and_cavity_share_one_transaction() {
    let air = basis(); let gate = CancelGate::new_clock_free(); let n = 82;
    let mut c = vec![0.0; n*2]; c[..4].copy_from_slice(&[-0.03, -0.02, 0.03, -0.02]);
    assert!(CavityCoupling::new(&air, n, &c, &[0.0; 2]).is_err());
    let build = || {
        let mut bodies = vec![body(800.0, 0.0, 0.1), body(900.0, 0.0, 0.0)];
        for i in 0..80 { bodies.push(body(600.0+i as f64, 0.0, 0.0)); }
        let mut rows = vec![0.0; 2*n]; rows[0] = 1.0; rows[2] = -1.0;
        rows[n+1] = 1.0; rows[n+3] = -1.0;
        let contact = Obstacle::new(rows, 2, n, vec![0.0; 2], vec![1.0; 2],
            1e6, 1.5, "two synthetic independent strand/head contacts".into()).unwrap()
            .with_internal_loss(0.05).unwrap();
        CavityCoupling::new_with_mode_budget(&air, n, &c, &[0.0; 2], 256).unwrap()
            .build_linear(bodies, vec![contact], 0.1, config(500_000), &gate).unwrap()
    };
    let (mut a, pressure) = build(); let (mut b, _) = build();
    assert_eq!(a.mode_count(), n+1); assert_eq!(a.contact_count(), 2);
    assert_eq!(a.state()[1], 0.1); assert!(a.state()[2..].iter().all(|x| *x == 0.0));
    let forces = vec![0.0; n+1]; let energy = a.frame().stored_energy_j; let mut loss = 0.0;
    for _ in 0..128 {
        let frame = a.step(&forces, &gate).unwrap(); b.step(&forces, &gate).unwrap();
        loss += frame.dissipated_energy_j; assert_eq!(a.state(), b.state());
    }
    assert!(a.state()[5] != 0.0, "contact must drive the first wire");
    assert!(a.state()[2*n+1] != 0.0, "head motion must drive air inertia");
    assert!((a.frame().stored_energy_j+loss-energy).abs() < 1e-9);
    assert!(pressure.pressure_at(a.state(), &[1.0, 0.3]).unwrap().is_finite());
    let before = a.state().to_vec(); let report = *a.frame();
    let stopped = CancelGate::new_clock_free(); stopped.request();
    assert!(matches!(a.step(&forces, &stopped), Err(ImpactError::Cancelled)));
    assert_eq!(a.state(), before); assert_eq!(*a.frame(), report);
    let mut invalid_force = forces.clone(); invalid_force[n] = f64::NAN;
    assert!(a.step(&invalid_force, &gate).is_err()); assert_eq!(a.state(), before);
    a.step(&forces, &gate).unwrap(); b.step(&forces, &gate).unwrap();
    assert_eq!(a.state(), b.state()); assert_eq!(a.frame(), b.frame());
}

#[test]
fn broader_prepared_budget_never_weakens_reference_or_loss_admission() {
    let air = basis(); let gate = CancelGate::new_clock_free();
    for maximum in [0, 1, 4097, usize::MAX] {
        assert!(CavityCoupling::new_with_mode_budget(&air, 2, &[0.1; 4], &[0.0; 2], maximum).is_err());
    }
    let lossy = CavityCoupling::new(&air, 1, &[0.1; 2], &[0.0, 5.0]).unwrap();
    assert!(lossy.build_linear(vec![body(800.0, 0.0, 0.0)], vec![], 0.1,
        config(500_000), &gate).is_err());
    let large = CavityCoupling::new_with_mode_budget(&air, 80, &vec![0.0; 160], &[0.0; 2], 256).unwrap();
    let parts = || (0..80).map(|_| body(800.0, 0.0, 0.0)).collect();
    let cfg = ImpactConfig { dt_s: 2e-6, max_steps: 10, maximum_energy_j: 20.0,
        energy_absolute_tolerance_j: 1e-10, energy_relative_tolerance: 1e-7,
        maximum_generalized_force: 10_000.0 };
    assert!(large.clone().build(parts(), vec![], vec![], cfg, &gate).is_err());
    let mut small = config(500_000); small.coupling.max_modes = 64;
    assert!(large.clone().build_linear(parts(), vec![], 0.1, small, &gate).is_err());
    gate.request();
    assert!(matches!(large.build_linear(parts(), vec![], 0.1, config(500_000), &gate),
        Err(ImpactError::Cancelled)));
}
