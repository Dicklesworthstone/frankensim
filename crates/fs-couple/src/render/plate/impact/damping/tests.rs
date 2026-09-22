use super::*;
use super::super::{BodyPotential, ImpactBody, ImpactConfig, ImpactSystem};
use super::super::linear::{LinearImpactConfig, LinearImpactSystem};
use crate::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use crate::render::schedule::force::coupled::{ModalCouplingConfig,
    contact::{ModalContactConfig, multiple::MultiContactConfig}};
use fs_exec::CancelGate;

fn config() -> ImpactConfig {
    ImpactConfig { dt_s: 1.0/20_000.0, max_steps: 200, maximum_energy_j: 10.0,
        energy_absolute_tolerance_j: 1e-10, energy_relative_tolerance: 1e-8,
        maximum_generalized_force: 1e6 }
}
fn linear_config() -> LinearImpactConfig {
    LinearImpactConfig { sample_rate_hz: 20_000, max_steps: 200, maximum_generalized_force: 1e6,
        component: ModalAcousticTimeBudget { maximum_total_energy_j: 10.0,
            ..ModalAcousticTimeBudget::audible_reference() },
        coupling: ModalCouplingConfig { max_modes: 64, max_connections: 8, max_setup_terms: 100_000,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 10.0, maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e5, solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-10, energy_relative_tolerance: 1e-8 },
        contact: ModalContactConfig { max_iterations: 100, maximum_force_n: 1e4,
            maximum_penetration_m: 0.01, force_absolute_tolerance_n: 1e-10, force_relative_tolerance: 1e-9 },
        multiple: MultiContactConfig { max_contacts: 32, max_sweeps: 100, max_setup_terms: 100_000 } }
}
fn parts() -> (Vec<ImpactBody>, ViscousDamper) {
    let (a, ba) = ImpactBody::free_mass(0.04, 0.0, 0.03).unwrap();
    let (b, bb) = ImpactBody::free_mass(0.06, 0.0, -0.01).unwrap();
    (vec![a,b], ViscousDamper { weights: vec![ba,-bb], damping_n_s_m: 0.8 })
}
fn reference(bodies: Vec<ImpactBody>, dampers: Vec<ViscousDamper>) -> ImpactSystem {
    ImpactSystem::new_with_dampers(bodies, vec![], vec![], vec![], dampers, config()).unwrap()
}

#[test]
fn spatial_resistance_keeps_cross_terms_and_nonnegative_port_power() {
    let dampers = vec![ViscousDamper { weights: vec![2.0,-3.0,0.0], damping_n_s_m: 0.7 }];
    let mut r = vec![0.0; 49]; // Includes one untouched Kelvin coordinate.
    add_resistance(&dampers, 3, 7, &mut r).unwrap();
    for v in [[1.0,0.0,0.0],[0.0,1.0,2.0],[3.0,2.0,7.0],[-0.4,0.6,-2.0]] {
        let power = (0..3).map(|i| (0..3).map(|j| v[i]*r[(2*i+1)*7+2*j+1]*v[j]).sum::<f64>()).sum::<f64>();
        let expected = 0.7*(2.0*v[0]-3.0*v[1]).powi(2);
        assert!((power-expected).abs() < 1e-12 && power > -1e-12);
    }
    assert!(r[1*7+3] < 0.0);
    for i in [0,2,4,6] { assert!(r[i*7..(i+1)*7].iter().all(|x| *x == 0.0)); }
}

#[test]
fn reciprocal_drag_matches_two_mass_solution_and_accounts_for_loss() {
    let gate = CancelGate::new_clock_free();
    let (bodies, damper) = parts();
    let mut model = reference(bodies.clone(), vec![damper.clone()]);
    let mut modal = LinearImpactSystem::new_with_dampers(bodies, vec![], vec![],
        vec![damper], linear_config(), &gate).unwrap();
    let initial = model.stored_energy_j(); let mut lost = 0.0;
    let lambda = 0.8*(1.0/0.04+1.0/0.06);
    let decay: f64 = (1.0-0.5*config().dt_s*lambda)/(1.0+0.5*config().dt_s*lambda);
    for step in 1..=80 {
        let frame = model.step(&[0.0;2], &gate).unwrap();
        let fast = modal.step(&[0.0;2], &gate).unwrap();
        lost += frame.dissipated_energy_j;
        let va = model.state()[1]/0.04_f64.sqrt();
        let vb = model.state()[3]/0.06_f64.sqrt();
        // Momentum and decay are independent force/trajectory checks, not just H.
        assert!((0.04*va+0.06*vb-0.0006).abs() < 2e-9);
        assert!((va-vb-0.04*decay.powi(step)).abs() < 2e-8);
        assert!(frame.dissipated_energy_j >= 0.0 && fast.dissipated_energy_j >= 0.0);
        for (a,b) in model.state().iter().zip(modal.state()) { assert!((a-b).abs() < 2e-8); }
    }
    assert!(lost > 0.0 && (model.stored_energy_j()+lost-initial).abs() < 1e-9);
    assert!(model.stored_energy_j() < initial*0.9);
}

#[test]
fn attachment_node_is_undamped_and_prepared_reference_retries_exactly() {
    let body = ImpactBody { potential: BodyPotential::Linear(vec![0.0;2]), damping_per_s: vec![0.0;2],
        initial: vec![ModalAcousticState { displacement_m_sqrt_kg: 0.0,
            velocity_m_sqrt_kg_per_s: 0.01 };2] };
    let damper = ViscousDamper { weights: vec![1.0,0.0], damping_n_s_m: 80.0 };
    let mut a = reference(vec![body.clone()], vec![damper.clone()]);
    let mut b = reference(vec![body], vec![damper]).prepare().unwrap();
    let gate = CancelGate::new_clock_free(); let cancel = CancelGate::new_clock_free(); cancel.request();
    let before = b.state().to_vec();
    assert!(b.step(&[0.0;2], &cancel).is_err());
    assert_eq!(b.state(), before); assert_eq!(b.samples(), 0);
    for _ in 0..40 {
        a.step(&[0.0;2], &gate).unwrap(); b.step(&[0.0;2], &gate).unwrap();
        assert!((a.state()[3]-0.01).abs() < 1e-12);
        for (x,y) in a.state().iter().zip(b.state()) { assert!((x-y).abs() < 1e-11); }
    }
    assert!(a.state()[1] < 0.009);
}

#[test]
fn zero_resistance_preserves_both_original_execution_paths() {
    let gate = CancelGate::new_clock_free(); let (bodies, mut damper) = parts(); damper.damping_n_s_m = 0.0;
    let mut a = ImpactSystem::new(bodies.clone(),vec![],vec![],vec![],config()).unwrap();
    let mut b = reference(bodies.clone(),vec![damper.clone()]);
    let mut c = LinearImpactSystem::new(bodies.clone(),vec![],vec![],linear_config(),&gate).unwrap();
    let mut d = LinearImpactSystem::new_with_dampers(bodies,vec![],vec![],vec![damper],linear_config(),&gate).unwrap();
    for _ in 0..20 {
        a.step(&[0.0;2],&gate).unwrap(); b.step(&[0.0;2],&gate).unwrap();
        let fc = c.step(&[0.0;2],&gate).unwrap(); let fd = d.step(&[0.0;2],&gate).unwrap();
        assert_eq!(a.state(),b.state()); assert_eq!(c.state(),d.state()); assert_eq!(fc,fd);
    }
}

#[test]
fn invalid_or_unrepresentable_dampers_refuse_without_silent_fallback() {
    for (weights, resistance) in [(vec![1.0],-1.0),(vec![1.0],f64::NAN),
        (vec![f64::INFINITY],1.0),(vec![0.0],1.0),(vec![],1.0),
        (vec![1e200],1.0),(vec![1e-200],1.0)] {
        assert!(validate(&[ViscousDamper {weights,damping_n_s_m:resistance}],1).is_err());
    }
    let damper = ViscousDamper { weights:vec![1.0],damping_n_s_m:1.0 };
    assert!(validate(&vec![damper;MAX_VISCOUS_DAMPERS+1],1).is_err());
    let (bodies, damper) = parts(); let mut limits = linear_config(); limits.coupling.max_connections = 0;
    assert!(LinearImpactSystem::new_with_dampers(bodies,vec![],vec![],vec![damper],limits,
        &CancelGate::new_clock_free()).is_err());
}

#[test]
fn acoustic_extension_never_attaches_a_solid_muffler_to_gas_inertia() {
    let (_, damper) = parts(); let weights = damper.weights.clone();
    let extended = extend(vec![damper],2,5).unwrap();
    assert_eq!(&extended[0].weights[..2],weights);
    assert_eq!(&extended[0].weights[2..],&[0.0;3]);
    let mut resistance = vec![0.0;100];
    add_resistance(&extended,5,10,&mut resistance).unwrap();
    for row in 4..10 { assert!(resistance[row*10..(row+1)*10].iter().all(|v| *v == 0.0)); }
}
