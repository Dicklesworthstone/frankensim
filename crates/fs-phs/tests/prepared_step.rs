use fs_phs::{PortHamiltonian, QuadraticStorage, StepWorkspace, Storage, step};

fn oscillator(damping: f64) -> PortHamiltonian {
    PortHamiltonian::new(2, 1, vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, damping], vec![0.0, 1.0],
        Box::new(QuadraticStorage::new(vec![9.0, 0.0, 0.0, 1.0], 2).unwrap())).unwrap()
}

#[derive(Debug)]
struct CoupledQuartic;
impl Storage for CoupledQuartic {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        let stretch = x[0] * x[0] + 0.3 * x[2] * x[2];
        0.5 * (7.0 * x[0] * x[0] + x[1] * x[1] + 11.0 * x[2] * x[2] + x[3] * x[3])
            + 20.0 * stretch * stretch
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        let stretch = x[0] * x[0] + 0.3 * x[2] * x[2];
        out[0] = 7.0 * x[0] + 80.0 * stretch * x[0];
        out[1] = x[1];
        out[2] = 11.0 * x[2] + 24.0 * stretch * x[2];
        out[3] = x[3];
    }
}
fn nonlinear() -> PortHamiltonian {
    let mut j = vec![0.0; 16];
    j[1] = 1.0; j[4] = -1.0; j[11] = 1.0; j[14] = -1.0;
    PortHamiltonian::new(4, 1, j, vec![0.0; 16], vec![0.0, 1.0, 0.0, 0.4], Box::new(CoupledQuartic)).unwrap()
}

#[test]
fn prepared_step_matches_analytic_implicit_midpoint() {
    let sys = oscillator(0.0);
    let mut work = StepWorkspace::new(&sys).unwrap();
    let (mut x, mut y) = ([0.0; 2], [0.0]);
    let dt = 0.01;
    let record = work.step_into(&sys, &[0.3, 0.7], &[0.0], dt, &mut x, &mut y).unwrap();
    let denominator = 1.0 + 0.25 * dt * dt * 9.0;
    let q = ((1.0 - 0.25 * dt * dt * 9.0) * 0.3 + dt * 0.7) / denominator;
    let p = (-dt * 9.0 * 0.3 + (1.0 - 0.25 * dt * dt * 9.0) * 0.7) / denominator;
    assert!((x[0] - q).abs() < 1e-11);
    assert!((x[1] - p).abs() < 1e-11);
    assert!((y[0] - 0.5 * (0.7 + p)).abs() < 1e-10);
    assert!(record.balance_residual().abs() < 1e-11);
}

#[test]
fn prepared_step_retains_nonlinear_mode_exchange_and_energy() {
    let sys = nonlinear();
    let mut work = StepWorkspace::new(&sys).unwrap();
    let mut actual = [0.25, 0.4, -0.13, 0.1];
    let mut reference = actual.to_vec();
    let initial_energy = sys.hamiltonian(&actual);
    let mut next = [0.0; 4];
    let mut y = [0.0];
    for _ in 0..400 {
        let expected = step(&sys, &reference, &[0.0], 0.002).unwrap();
        let record = work.step_into(&sys, &actual, &[0.0], 0.002, &mut next, &mut y).unwrap();
        for (a, b) in next.iter().zip(&expected.x) { assert!((a - b).abs() < 2e-8); }
        assert!((y[0] - expected.y[0]).abs() < 2e-8);
        assert!(record.balance_residual().abs() < 1e-9);
        assert!((sys.hamiltonian(&next) - initial_energy).abs() < 1e-8);
        actual = next;
        reference = expected.x;
    }
}

#[test]
fn prepared_step_accounts_for_external_work_and_dissipation() {
    let sys = oscillator(0.7);
    let mut work = StepWorkspace::new(&sys).unwrap();
    let mut state = [0.2, -0.4];
    let mut next = [0.0; 2]; let mut y = [0.0];
    for i in 0..100 {
        let force = if i < 30 { 0.8 } else { 0.0 };
        let record = work.step_into(&sys, &state, &[force], 0.003, &mut next, &mut y).unwrap();
        assert!((record.supplied - 0.003 * force * y[0]).abs() < 1e-14);
        assert!(record.dissipated >= 0.0);
        assert!(record.supply_defect() <= 1e-10);
        assert!(record.balance_residual().abs() < 1e-10);
        state = next;
    }
}

#[test]
fn prepared_step_failures_preserve_both_outputs_and_do_not_poison_retry() {
    let sys = oscillator(0.0);
    let mut work = StepWorkspace::new(&sys).unwrap();
    let (mut next, mut y) = ([123.0, 456.0], [789.0]);
    for (state, input, dt) in [([f64::NAN, 0.0], [0.0], 0.01), ([0.0, 0.0], [f64::INFINITY], 0.01), ([0.0, 0.0], [0.0], -0.01)] {
        assert!(work.step_into(&sys, &state, &input, dt, &mut next, &mut y).is_err());
        assert_eq!(next, [123.0, 456.0]); assert_eq!(y, [789.0]);
    }
    assert!(work.step_into(&sys, &[0.0], &[0.0], 0.01, &mut next, &mut y).is_err());
    work.set_iteration_limit(0).unwrap();
    assert!(work.step_into(&sys, &[0.3, 0.7], &[0.0], 0.01, &mut next, &mut y).is_err());
    assert_eq!(next, [123.0, 456.0]); assert_eq!(y, [789.0]);
    work.set_iteration_limit(50).unwrap();
    work.step_into(&sys, &[0.3, 0.7], &[0.0], 0.01, &mut next, &mut y).unwrap();
    let expected = step(&sys, &[0.3, 0.7], &[0.0], 0.01).unwrap();
    for (a, b) in next.iter().zip(expected.x) { assert!((a - b).abs() < 1e-10); }
}

#[test]
fn zero_duration_is_a_noop_with_the_correct_port_output() {
    let sys = oscillator(0.5);
    let mut work = StepWorkspace::new(&sys).unwrap();
    let (mut next, mut y) = ([0.0; 2], [0.0]);
    let record = work.step_into(&sys, &[0.3, 0.7], &[2.0], 0.0, &mut next, &mut y).unwrap();
    assert_eq!(next, [0.3, 0.7]); assert_eq!(y, [0.7]);
    assert_eq!(record.delta_h, 0.0); assert_eq!(record.supplied, 0.0); assert_eq!(record.dissipated, 0.0);
}

#[test]
fn supply_audit_still_detects_a_nonpassive_raw_system() {
    let sys = PortHamiltonian::from_raw_parts(1, 0, vec![0.2], vec![0.0], vec![],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()));
    let mut work = StepWorkspace::new(&sys).unwrap();
    let mut next = [0.0];
    let record = work.step_into(&sys, &[1.0], &[], 0.01, &mut next, &mut []).unwrap();
    assert!(record.supply_defect() > 0.001);
}

#[test]
fn cancellation_inside_jacobian_preserves_outputs_and_exact_retry() {
    let sys = nonlinear();
    let x = [0.25, 0.4, -0.13, 0.1];
    let mut work = StepWorkspace::new(&sys).unwrap();
    let (mut next, mut y) = ([123.0; 4], [456.0]);
    let mut polls = 0;
    let error = work.step_into_controlled(&sys, &x, &[0.2], 0.002,
        &mut next, &mut y, || { polls += 1; polls == 5 }).unwrap_err();
    assert_eq!(error, fs_phs::PreparedStepError::Cancelled);
    assert_eq!(polls, 5);
    assert_eq!(next, [123.0; 4]);
    assert_eq!(y, [456.0]);
    work.step_into(&sys, &x, &[0.2], 0.002, &mut next, &mut y).unwrap();
    let mut clean = StepWorkspace::new(&sys).unwrap();
    let (mut expected, mut expected_y) = ([0.0; 4], [0.0]);
    clean.step_into(&sys, &x, &[0.2], 0.002, &mut expected, &mut expected_y).unwrap();
    assert_eq!(next.map(f64::to_bits), expected.map(f64::to_bits));
    assert_eq!(y.map(f64::to_bits), expected_y.map(f64::to_bits));
}

#[test]
fn cancellation_at_publication_does_not_leak_a_converged_candidate() {
    let sys = oscillator(0.7);
    let mut work = StepWorkspace::new(&sys).unwrap();
    let (mut next, mut y) = ([0.0; 2], [0.0]);
    let mut total_polls = 0;
    work.step_into_controlled(&sys, &[0.2, -0.4], &[0.8], 0.003,
        &mut next, &mut y, || { total_polls += 1; false }).unwrap();
    assert!(total_polls > 5);
    next = [123.0; 2]; y = [456.0];
    let mut polls = 0;
    let error = work.step_into_controlled(&sys, &[0.2, -0.4], &[0.8], 0.003,
        &mut next, &mut y, || { polls += 1; polls == total_polls }).unwrap_err();
    assert_eq!(error, fs_phs::PreparedStepError::Cancelled);
    assert!(work.last_step_discrete_effort().is_none());
    assert_eq!(next, [123.0; 2]);
    assert_eq!(y, [456.0]);
}

#[test]
fn prepared_workspace_reads_the_current_operators_not_a_cached_model() {
    let undamped = oscillator(0.0);
    let damped = oscillator(1.5);
    let mut work = StepWorkspace::new(&undamped).unwrap();
    let mut clean = StepWorkspace::new(&damped).unwrap();
    let (mut actual, mut actual_y) = ([0.0; 2], [0.0]);
    let (mut expected, mut expected_y) = ([0.0; 2], [0.0]);
    work.step_into(&undamped, &[0.2, 0.3], &[0.1], 0.002, &mut actual, &mut actual_y).unwrap();
    work.step_into(&damped, &[0.2, 0.3], &[0.1], 0.002, &mut actual, &mut actual_y).unwrap();
    clean.step_into(&damped, &[0.2, 0.3], &[0.1], 0.002, &mut expected, &mut expected_y).unwrap();
    assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
    assert_eq!(actual_y.map(f64::to_bits), expected_y.map(f64::to_bits));
}

#[test]
fn prepared_maxwell_history_retains_reference_storage_and_work() {
    let sys = oscillator(0.1).with_relaxation_branches(vec![fs_phs::RelaxationBranch {
        projection: vec![1.0, 0.0], stiffness: 3.0, relaxation_time_s: 0.012,
    }]).unwrap();
    let mut work = StepWorkspace::new(&sys).unwrap();
    let mut x = [0.02, -0.01, 0.0];
    let mut reference = x.to_vec();
    let (mut next, mut y) = ([0.0; 3], [0.0]);
    for _ in 0..100 {
        let expected = step(&sys, &reference, &[0.02], 0.0005).unwrap();
        let record = work.step_into(&sys, &x, &[0.02], 0.0005, &mut next, &mut y).unwrap();
        assert!(record.dissipated >= 0.0);
        assert!(record.balance_residual().abs() < 1e-11);
        for (a, b) in next.iter().zip(&expected.x) { assert!((a - b).abs() < 1e-9); }
        x = next;
        reference = expected.x;
    }
    assert_ne!(x[2], 0.0);
}

struct DatumQuartic(f64);
impl Storage for DatumQuartic {
    fn hamiltonian(&self, x: &[f64]) -> f64 { self.0 + 0.25*x[0].powi(4) }
    fn gradient(&self, x: &[f64], out: &mut [f64]) { out[0] = x[0].powi(3); }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        let (left, right, m) = (a[0]-mid[0], b[0]-mid[0], mid[0]);
        Some(1.5*m*m*(right*right-left*left) + m*(right.powi(3)-left.powi(3))
            + 0.25*(right.powi(4)-left.powi(4)))
    }
}

#[test]
fn retained_effort_preserves_nonlinear_work_and_large_energy_datum() {
    let hessian = |x: &[f64], d: &[f64], out: &mut [f64]| {
        out[0] = 3.0*x[0]*x[0]*d[0]; true
    };
    for datum in [0.0, 1e30] { for analytic in [false, true] {
        let sys = PortHamiltonian::new(1,1,vec![0.0],vec![0.7],vec![1.0],
            Box::new(DatumQuartic(datum))).unwrap();
        let mut work = StepWorkspace::new(&sys).unwrap();
        assert!(work.last_step_discrete_effort().is_none());
        let initial = [1.0]; let input = [0.2]; let dt = 0.1;
        let (mut next, mut y) = ([0.0], [0.0]);
        let record = if analytic {
            work.step_into_analytic(&sys,&initial,&input,dt,&mut next,&mut y,&hessian)
        } else { work.step_into(&sys,&initial,&input,dt,&mut next,&mut y) }.unwrap();
        let effort = work.last_step_discrete_effort().unwrap()[0];
        // Independent secant of x^4/4, without subtracting either datum or
        // neighboring endpoint energies. In one dimension this is the DG.
        let exact = 0.25*(initial[0]+next[0])
            *(initial[0]*initial[0]+next[0]*next[0]);
        assert!((effort-exact).abs() < 3e-14);
        assert!((effort-f64::midpoint(initial[0],next[0]).powi(3)).abs() > 1e-5,
            "exercise the correction overwritten by storage-energy scratch");
        assert_eq!(effort.to_bits(),y[0].to_bits());
        assert!((record.dissipated-dt*0.7*effort*effort).abs() < 3e-15);
        assert!((record.supplied-dt*input[0]*effort).abs() < 3e-15);
        assert!((record.delta_h-exact*(next[0]-initial[0])).abs() < 3e-15);
    }}
}

#[test]
fn failed_attempts_hide_old_effort_and_retry_publishes_current_constitutive_work() {
    let sys = oscillator(0.7);
    let mut work = StepWorkspace::new(&sys).unwrap();
    let initial = [0.2,-0.4];
    let (mut next, mut y) = ([0.0;2],[0.0]);
    for failure in 0..4 {
        work.set_iteration_limit(50).unwrap();
        work.step_into(&sys,&initial,&[0.8],0.003,&mut next,&mut y).unwrap();
        assert!(work.last_step_discrete_effort().is_some());
        let accepted = (next,y);
        let refused = match failure {
            0 => work.step_into_controlled(&sys,&initial,&[0.8],0.003,
                &mut next,&mut y,||true),
            1 => work.step_into_controlled(&sys,&initial[..1],&[0.8],0.003,
                &mut next,&mut y,||false),
            2 => work.step_into_controlled(&sys,&initial,&[0.8],f64::NAN,
                &mut next,&mut y,||false),
            _ => {
                work.set_iteration_limit(0).unwrap();
                work.step_into_controlled(&sys,&initial,&[0.8],0.003,
                    &mut next,&mut y,||false)
            }
        };
        assert!(refused.is_err());
        assert!(work.last_step_discrete_effort().is_none());
        assert_eq!((next,y),accepted);
    }
    work.set_iteration_limit(50).unwrap();
    let passive = |_: &[f64], e: &[f64], out: &mut [f64]| {
        out[0]=0.0; out[1]=0.3*e[1]; true
    };
    let record = work.step_into_dissipative_controlled(&sys,&initial,&[0.8],0.003,
        &mut next,&mut y,&passive,None,None,||false).unwrap();
    let effort = work.last_step_discrete_effort().unwrap();
    assert_eq!(effort[1].to_bits(),y[0].to_bits());
    assert!((record.dissipated-0.003*effort[1]*effort[1]).abs()<1e-15);
}

// Synthetic, mass-normalized striker/two-receiver fixture. It deliberately
// crosses the Hertz surface from free flight with a large launch momentum.
// The old all-time-best stagnation counter refused at first contact even while
// Newton was recovering. This is not a measured drum or a new production law.
struct ContactRecovery {
    stiffness: [[f64; 3]; 3],
    contact: [f64; 3],
    hertz: f64,
}
impl Storage for ContactRecovery {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        let mut energy = (0..3).map(|i| 0.5*x[2*i+1]*x[2*i+1]).sum::<f64>();
        for i in 0..3 { for j in 0..3 {
            energy += 0.5*x[2*i]*self.stiffness[i][j]*x[2*j];
        }}
        let penetration = self.contact.iter().enumerate().map(|(i,b)|b*x[2*i]).sum::<f64>().max(0.0);
        energy + self.hertz/2.5*fs_math::det::pow(penetration,2.5)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        let penetration = self.contact.iter().enumerate().map(|(i,b)|b*x[2*i]).sum::<f64>().max(0.0);
        let force = self.hertz*fs_math::det::pow(penetration,1.5);
        for i in 0..3 {
            out[2*i] = (0..3).map(|j|self.stiffness[i][j]*x[2*j]).sum::<f64>()+self.contact[i]*force;
            out[2*i+1] = x[2*i+1];
        }
    }
}

#[test]
fn recovering_hertz_newton_is_not_stalled_until_it_beats_the_free_flight_residual() {
    let mass = 0.00724314_f64;
    let roots = [0.0, core::f64::consts::TAU*204.6189730766819,
        core::f64::consts::TAU*264.1200719782793];
    let areas = [0.0,0.22,-0.3];
    let bulk_over_volume = 1.2*343.0*343.0/(core::f64::consts::PI*0.1703*0.1703*0.1651);
    let stiffness: [[f64;3];3] = core::array::from_fn(|i|core::array::from_fn(|j|
        (if i==j {roots[i]*roots[i]} else {0.0}) + bulk_over_volume*areas[i]*areas[j]));
    let mut j = vec![0.0;36]; let mut r = vec![0.0;36];
    for i in 0..3 {
        j[(2*i)*6+2*i+1] = 1.0; j[(2*i+1)*6+2*i] = -1.0;
        if i>0 {r[(2*i+1)*6+2*i+1] = 0.002*stiffness[i][i].sqrt();}
    }
    let system = PortHamiltonian::new(6,0,j,r,vec![],Box::new(ContactRecovery {
        stiffness,contact:[1.0/mass.sqrt(),-16.0,0.0],hertz:4.0/3.0*0.8e9*0.003_f64.sqrt(),
    })).unwrap();
    let mut workspace = StepWorkspace::new(&system).unwrap();
    let mut state = [-0.0002*mass.sqrt(),4.0*mass.sqrt(),0.0,0.0,0.0,0.0];
    let mut reference = state.to_vec(); let mut candidate = [0.0;6];
    let initial_energy = system.hamiltonian(&state);
    let mut loss = 0.0; let mut iterations = 0; let mut receiver_peak = 0.0_f64;
    for _ in 0..512 {
        let expected = step(&system,&reference,&[],2e-6).unwrap();
        let actual = workspace.step_into(&system,&state,&[],2e-6,&mut candidate,&mut []).unwrap();
        for (a,b) in candidate.iter().zip(&expected.x) {assert!((a-b).abs()<1e-8);}
        assert!(actual.balance_residual().abs()<1e-9);
        assert!(actual.supply_defect()<=1e-9);
        assert!(actual.newton_iters<=50);
        iterations = iterations.max(actual.newton_iters);
        receiver_peak = receiver_peak.max(candidate[2].abs());
        loss += actual.dissipated; state = candidate; reference = expected.x;
    }
    assert!(iterations>3,"exercise recovery past the former premature stagnation stop");
    assert!(receiver_peak>1e-6,"contact must actually excite the receiver");
    assert!((system.hamiltonian(&state)+loss-initial_energy).abs()<1e-8);
}
