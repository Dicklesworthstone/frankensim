use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState,
    ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::plate::impact::cavity::CavityCoupling;
use fs_couple::vibroacoustic::CavityModes;
use fs_math::c64::C64;
use fs_phs::PortExchangeBudget;

fn budget() -> PortExchangeBudget {
    PortExchangeBudget { max_left:128, max_right:8, max_setup_terms:1_000_000,
        maximum_dt_coupling:0.125 }
}

fn basis(omegas: Vec<f64>) -> CavityModes {
    let count = omegas.len();
    CavityModes { omegas, lambdas:vec![0.1;count], interface:vec![vec![1.0];count],
        loss_factor:0.0, rho0:1.2, c0:343.0 }
}

#[test]
fn uniform_compression_has_no_momentum_and_converges_to_the_coupled_frequency() {
    let air = basis(vec![0.0]);
    let root_a = (air.rho0*air.c0*air.c0/air.lambdas[0]).sqrt();
    let coupling = 300.0;
    let model = CavityCoupling::new(&air,1,&[coupling/root_a],&[0.0]).unwrap();
    assert_eq!(model.total_modes(),1);
    let structural_omega = 800.0_f64;
    let full_omega = (structural_omega.powi(2)+coupling.powi(2)).sqrt();
    let initial_velocity = 0.04;
    let time = 0.02;
    let exact_q = initial_velocity/full_omega*(full_omega*time).sin();
    let exact_v = initial_velocity*(full_omega*time).cos();
    let run = |rate| {
        let mut cavity = model.prepare_exchange(rate,budget()).unwrap();
        let mut mechanical = ModalAcousticTimeModel::try_new(rate,vec![ModalAcousticMode {
            angular_frequency_rad_s:structural_omega, damping_ratio:0.0,
            pressure_per_modal_velocity:C64::ZERO,
        }],ModalAcousticTimeBudget::audible_reference()).unwrap();
        let mut q = 0.0; let mut v = [initial_velocity]; let mut loss = 0.0;
        for _ in 0..(time*f64::from(rate)).round() as usize {
            loss += cavity.before(&mut v).unwrap();
            mechanical.restore_states(&[ModalAcousticState {
                displacement_m_sqrt_kg:q, velocity_m_sqrt_kg_per_s:v[0],
            }]).unwrap();
            mechanical.step(&[0.0]).unwrap();
            q = mechanical.states()[0].displacement_m_sqrt_kg;
            v[0] = mechanical.states()[0].velocity_m_sqrt_kg_per_s;
            loss += cavity.after(&mut v).unwrap();
        }
        assert_eq!(loss,0.0,"uniform air has no dissipative momentum");
        let energy = 0.5*((structural_omega*q).powi(2)+v[0]*v[0])+cavity.energy();
        assert!((energy-0.5*initial_velocity.powi(2)).abs() < 2e-13);
        let mut pressure = [0.0]; cavity.pressures_into(&mut pressure).unwrap();
        let pressure_energy = -pressure[0]/root_a;
        assert!((cavity.energy()-0.5*pressure_energy*pressure_energy).abs() < 1e-16,
            "uniform air must contain compression only, with no added momentum");
        [(q-exact_q).abs()*full_omega,(v[0]-exact_v).abs(),
            (pressure_energy-coupling*exact_q).abs()].into_iter().fold(0.0,f64::max)
    };
    let coarse = run(8_000); let fine = run(16_000);
    assert!(coarse > 1e-9 && fine < 0.3*coarse,
        "uniform cavity must approach sqrt(omega_board²+A*C²): {coarse:e} -> {fine:e}");
}

#[test]
fn distributed_air_exchanges_work_and_drag_removes_only_acoustic_momentum_energy() {
    let air = basis(vec![0.0,700.0,1800.0]);
    let overlaps = [0.08,0.04,-0.03,-0.02,0.06,0.07];
    for drag in [0.0,30.0] {
        let model = CavityCoupling::new(&air,2,&overlaps,&[0.0,drag,2.0*drag]).unwrap();
        let mut cavity = model.prepare_exchange(48_000,budget()).unwrap();
        assert_eq!(model.cavity_modes(),3); assert_eq!(model.total_modes(),4);
        let mut velocity = [0.04,-0.03];
        let initial = 0.5*(velocity[0]*velocity[0]+velocity[1]*velocity[1]);
        let mut loss = 0.0; let mut peak_pressure = 0.0_f64; let mut peak_inertia = 0.0_f64;
        for _ in 0..4000 {
            loss += cavity.before(&mut velocity).unwrap();
            loss += cavity.after(&mut velocity).unwrap();
            let mut pressure = [0.0;3]; cavity.pressures_into(&mut pressure).unwrap();
            peak_pressure = peak_pressure.max(pressure.iter().map(|p|p.abs()).fold(0.0,f64::max));
            let compression = pressure.iter().zip(&air.lambdas).map(|(p,lambda)|
                0.5*p*p*lambda/(air.rho0*air.c0*air.c0)).sum::<f64>();
            peak_inertia = peak_inertia.max(cavity.energy()-compression);
            let mechanical_energy = 0.5*(velocity[0]*velocity[0]+velocity[1]*velocity[1]);
            assert!((cavity.energy()+mechanical_energy+loss-initial).abs() < 2e-12);
        }
        assert!(peak_pressure > 0.1 && cavity.energy() > 1e-8);
        assert!(peak_inertia > 1e-8,"standing-wave inertia must move");
        if drag == 0.0 { assert!(loss.abs() < 1e-12); }
        else { assert!(loss > 1e-6,"declared momentum drag must remove energy"); }
    }
}

#[test]
fn impedance_retains_uniform_compliance_and_the_declared_momentum_drag_numerator() {
    let uniform = basis(vec![0.0]);
    let a = uniform.rho0*uniform.c0*uniform.c0/uniform.lambdas[0];
    let model = CavityCoupling::new(&uniform,2,&[0.08,-0.02],&[0.0]).unwrap();
    let z = model.impedance(500.0).unwrap();
    assert_eq!(z[0].re,0.0); assert_eq!(z[1],z[2]);
    assert!((z[0].im-a*0.08*0.08/500.0).abs() < 1e-13);
    assert!((z[1].im+a*0.08*0.02/500.0).abs() < 1e-13);

    let standing = basis(vec![700.0]); let drag = 40.0; let query = 100.0;
    let model = CavityCoupling::new(&standing,2,&[0.08,-0.02],&[drag]).unwrap();
    let z = model.impedance(query).unwrap();
    let denominator = (700.0_f64.powi(2)-query*query).powi(2)+(drag*query).powi(2);
    let expected_real = a*0.08*0.08*drag*700.0_f64.powi(2)/denominator;
    assert!((z[0].re-expected_real).abs() < 1e-13*expected_real);
    for velocity in [[1.0,0.0],[0.0,1.0],[1.0,-2.0],[1.0,4.0]] {
        let power = (0..2).map(|i|(0..2).map(|j|
            velocity[i]*z[2*i+j].re*velocity[j]).sum::<f64>()).sum::<f64>();
        assert!(power >= -1e-13);
    }
    let lossless = CavityCoupling::new(&standing,2,&[0.08,-0.02],&[0.0]).unwrap();
    assert!(lossless.impedance(700.0).is_err());
}

#[test]
fn failed_half_steps_and_frame_restore_preserve_all_pressure_and_momentum_history() {
    let model = CavityCoupling::new(&basis(vec![0.0,700.0]),2,
        &[0.08,0.04,-0.02,0.06],&[0.0,30.0]).unwrap();
    let mut cavity = model.prepare_exchange(48_000,budget()).unwrap();
    let mut velocity = [0.04,-0.03];
    cavity.before(&mut velocity).unwrap(); cavity.after(&mut velocity).unwrap();
    cavity.checkpoint(); let mut old_pressure = [0.0;2];
    cavity.pressures_into(&mut old_pressure).unwrap(); let old_energy = cavity.energy();
    let old_velocity = velocity;
    let loss = cavity.before(&mut velocity).unwrap()+cavity.after(&mut velocity).unwrap();
    let mut next_pressure = [0.0;2]; cavity.pressures_into(&mut next_pressure).unwrap();
    let next_energy = cavity.energy(); let next_velocity = velocity;
    cavity.restore(); let mut restored_pressure = [0.0;2];
    cavity.pressures_into(&mut restored_pressure).unwrap();
    assert_eq!(restored_pressure,old_pressure); assert_eq!(cavity.energy(),old_energy);
    // after() first advances its staged damped history, then the existing
    // exchange owner refuses this finite input's overflowing kinetic energy.
    let mut overflow = [f64::MAX,0.0];
    assert!(cavity.after(&mut overflow).is_err());
    assert_eq!(overflow,[f64::MAX,0.0]);
    cavity.pressures_into(&mut restored_pressure).unwrap();
    assert_eq!(restored_pressure,old_pressure); assert_eq!(cavity.energy(),old_energy);
    let mut wrong_shape = [0.0]; assert!(cavity.before(&mut wrong_shape).is_err());
    let mut wrong_pressure = [42.0]; assert!(cavity.pressures_into(&mut wrong_pressure).is_err());
    assert_eq!(wrong_pressure,[42.0]);
    velocity = old_velocity;
    let replay_loss = cavity.before(&mut velocity).unwrap()+cavity.after(&mut velocity).unwrap();
    assert_eq!(loss.to_bits(),replay_loss.to_bits());
    cavity.pressures_into(&mut restored_pressure).unwrap();
    assert_eq!(restored_pressure,next_pressure); assert_eq!(cavity.energy(),next_energy);
    assert_eq!(velocity,next_velocity);
    assert!(model.prepare_exchange(0,budget()).is_err());
    assert!(model.prepare_exchange(100,budget()).is_err());
    let mut narrow = budget(); narrow.max_left = 1;
    assert!(model.prepare_exchange(48_000,narrow).is_err());
}

#[test]
fn bordered_dynamic_stiffness_reuses_the_same_springs_and_survives_fixed_wall_poles() {
    let air = basis(vec![0.0,700.0]);
    let a = air.rho0*air.c0*air.c0/air.lambdas[0];
    for drag in [0.0,30.0] {
        let model = CavityCoupling::new(&air,2,&[0.08,0.04,-0.02,0.06],&[0.0,drag]).unwrap();
        let n = model.total_modes(); assert_eq!(n,3);
        for query in [100.0,700.0] {
            let actual = model.dynamic_stiffness(query).unwrap();
            let mut expected = [C64::ZERO;9];
            for column in [[0.08,-0.02,0.0],[0.04,0.06,700.0/a.sqrt()]] {
                for r in 0..3 { for c in 0..3 { expected[3*r+c].re += a*column[r]*column[c]; } }
            }
            expected[8] = expected[8]+C64::new(-query*query,-query*drag);
            for (&actual,expected) in actual.iter().zip(expected) {
                assert!((actual-expected).abs() < 1e-9);
            }
            // Away from a fixed-wall pole, eliminating the retained acoustic
            // inertia must give s*Z in the original structural velocity basis.
            if query != 700.0 || drag != 0.0 {
                let z = model.impedance(query).unwrap();
                for r in 0..2 { for c in 0..2 {
                    let schur = actual[3*r+c]-actual[3*r+2]*actual[6+c]/actual[8];
                    let eliminated = C64::new(0.0,-query)*z[2*r+c];
                    assert!((schur-eliminated).abs() < 1e-8);
                } }
            }
        }
    }
}
