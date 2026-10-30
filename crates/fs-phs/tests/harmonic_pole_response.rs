//! Physical participation in original mechanical, circuit and material models.
use fs_la::pole_response::PoleResponseOptions;
use fs_math::c64::C64;
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_phs::harmonic::{HarmonicBudget, HarmonicModel, HarmonicPole, PoleOptions};
use fs_phs::harmonic::response::HarmonicPoleResponseOptions;
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch};
use std::cell::Cell;

fn budget() -> HarmonicBudget {
    HarmonicBudget { max_dimension: 32, max_dense_bytes: 2_000_000,
        max_dense_work: 10_000_000, relative_residual_tolerance: 1e-11 }
}
fn options(scale: f64) -> HarmonicPoleResponseOptions {
    HarmonicPoleResponseOptions {
        numerical: PoleResponseOptions {
            spectral_scale: scale, pole_residual_tolerance: 1e-12,
            pole_step_tolerance: 1e-12, linear_residual_tolerance: 1e-11,
            max_dimension: 32, max_vectors: 8, max_dense_bytes: 2_000_000,
            max_work: 10_000_000,
        }, max_physical_bytes: 2_000_000, max_physical_work: 10_000_000,
    }
}
fn pole(model: &HarmonicModel, s: C64, seed: &[C64], scale: f64) -> HarmonicPole {
    model.resonance_pole(s, seed, &PoleOptions {
        spectral_scale: scale, step_tolerance: 1e-13, residual_tolerance: 1e-13,
        max_iterations: 40, max_dimension: 32, max_dense_bytes: 2_000_000,
        max_work: 10_000_000,
    }, || false).unwrap()
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance*b.abs().max(1e-14), "{a:?} != {b:?}");
}
fn mechanical() -> (HarmonicModel, HarmonicPole) {
    let system = fs_phs::mass_spring_damper(2.5, 1600.0, 0.6).unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let s = C64::new(-0.12, -(640.0_f64-0.12*0.12).sqrt());
    let p = pole(&model, s, &[C64::ONE, s.scale(2.5)], 26.0);
    (model, p)
}

#[test]
fn g1_mechanical_force_residue_keeps_units_phase_and_both_frequency_charts() {
    let (model, pole) = mechanical();
    let prepared = model.prepare_pole_response(&pole, options(26.0), || false).unwrap();
    let force = C64::new(1.7, -0.3);
    let response = prepared.port_residues(&[&[force]], || false).unwrap().remove(0);
    let s = pole.laplace_pole_per_s;
    let q = force/(s.scale(5.0)+C64::from_re(0.6));
    close(response.state[0], q, 1e-10);
    close(response.state[1], s.scale(2.5)*q, 1e-10);
    close(response.effort[0], q.scale(1600.0), 1e-10);
    close(response.output[0], s*q, 1e-10);
    // Independent resolvent, not the prepared factor: H_q(s)=1/(ms^2+ds+k).
    let mut errors = Vec::new();
    for h in [0.2, 0.1, 0.05] {
        let delta = C64::new(h, 0.4*h);
        let z = s+delta;
        let x = force/((z*z).scale(2.5)+z.scale(0.6)+C64::from_re(1600.0));
        errors.push((delta*x-q).abs());
        let omega = C64::new(-s.im, s.re)+C64::new(-delta.im, delta.re);
        let iw = C64::new(omega.im, -omega.re);
        let actual = force/((iw*iw).scale(2.5)+iw.scale(0.6)+C64::from_re(1600.0));
        close(actual, x, 1e-12);
    }
    assert!(errors[0]/errors[1] > 1.9 && errors[1]/errors[2] > 1.9);
    // Every original component, not just the output, shares one participation.
    for (&x,&v) in response.state.iter().zip(prepared.right_state()) {
        close(x, v*response.participation, 1e-12);
    }
}

#[test]
fn g1_kirchhoff_descriptor_retains_constraint_reactions_in_the_residue() {
    let (r,l,c) = (0.2_f64,0.1_f64,2.2e-7_f64);
    let circuit = assemble_circuit(&CircuitGraph {
        node_count: 4, branches: vec![
            (1,0,Branch::VoltageSource { port:0 }),
            (1,2,Branch::Inductor { henries:l }),
            (2,3,Branch::Capacitor { farads:c }),
            (3,0,Branch::Resistor { ohms:r }),
        ], transformers: vec![],
    }).unwrap();
    let model = HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap();
    let omega = 1.0/(l*c).sqrt();
    let seed = model.factor(omega,||false).unwrap().response(&[C64::ONE],||false).unwrap().state;
    let s = C64::new(-r/(2.0*l), -(omega*omega-(r/(2.0*l)).powi(2)).sqrt());
    let p = pole(&model,s,&seed,omega);
    let prepared = model.prepare_pole_response(&p,options(omega),||false).unwrap();
    let voltage = C64::new(0.7,0.2);
    let out = prepared.port_residues(&[&[voltage]],||false).unwrap().remove(0);
    let current = s*voltage/(s.scale(2.0*l)+C64::from_re(r));
    close(out.state[circuit.flux_index[0]],current.scale(l),1e-9);
    close(out.state[circuit.node_potential_index[2]],current.scale(r),1e-9);
    close(out.output[0],current,1e-9);
    assert!(out.state[circuit.node_potential_index[0]].abs()<1e-10);
    assert_eq!(out.state.len(),seed.len());
    assert!(prepared.physical_residual()<1e-11);
}

#[test]
fn g1_maxwell_oscillation_and_pure_relaxation_keep_their_full_material_residues() {
    let system = fs_phs::mass_spring_damper(1.0,60_000.0,0.0).unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection:vec![1.0,0.0],stiffness:50_000.0,relaxation_time_s:1.0/300.0,
        }]).unwrap();
    let model = HarmonicModel::from_ode(&system,budget(),||false).unwrap();
    for s in [C64::new(-50.0,-87_500.0_f64.sqrt()),C64::from_re(-200.0)] {
        let memory = C64::from_re(50_000.0_f64.sqrt())/(C64::ONE+s.scale(1.0/300.0));
        let p = pole(&model,s,&[C64::ONE,s,memory],400.0);
        let prepared = model.prepare_pole_response(&p,options(400.0),||false).unwrap();
        let actual = prepared.port_residues(&[&[C64::ONE]],||false).unwrap().remove(0);
        let d = C64::ONE+s.scale(1.0/300.0);
        let expected = (s.scale(2.0)+C64::from_re(50_000.0/300.0)/(d*d)).recip();
        close(actual.state[0],expected,1e-10);
        close(actual.state[1],s*expected,1e-10);
        close(actual.state[2],memory*expected,1e-10);
        close(actual.output[0],s*expected,1e-10);
    }
}

#[test]
fn g3_an_unobservable_internal_mode_is_retained_with_zero_port_participation() {
    let j = vec![0.,0.,1.,0., 0.,0.,0.,1., -1.,0.,0.,0., 0.,-1.,0.,0.];
    let mut r = vec![0.;16]; r[10]=0.2; r[15]=0.4;
    let mut q = vec![0.;16]; q[0]=100.;q[5]=225.;q[10]=1.;q[15]=1.;
    let system = PortHamiltonian::new(4,1,j,r,vec![0.,0.,1.,0.],
        Box::new(QuadraticStorage::new(q,4).unwrap())).unwrap();
    let model = HarmonicModel::from_ode(&system,budget(),||false).unwrap();
    let s = C64::new(-0.2,-(225.0_f64-0.04).sqrt());
    let p = pole(&model,s,&[C64::ZERO,C64::ONE,C64::ZERO,s],15.0);
    let prepared = model.prepare_pole_response(&p,options(15.0),||false).unwrap();
    let out = prepared.port_residues(&[&[C64::ONE]],||false).unwrap().remove(0);
    assert_eq!(prepared.right_state().len(),4);
    assert!(prepared.right_state()[1].abs()>0.0);
    assert!(out.state.iter().chain(&out.output).all(|z|z.abs()<1e-14));
    // A simple zero pole is an integrator, not an invalid resonance frequency.
    let integrator = PortHamiltonian::new(1,1,vec![0.],vec![0.],vec![1.],
        Box::new(QuadraticStorage::new(vec![2.],1).unwrap())).unwrap();
    let model = HarmonicModel::from_ode(&integrator,budget(),||false).unwrap();
    let p = pole(&model,C64::ZERO,&[C64::ONE],1.0);
    let prepared = model.prepare_pole_response(&p,options(1.0),||false).unwrap();
    let out = prepared.port_residues(&[&[C64::ONE]],||false).unwrap().remove(0);
    close(out.state[0],C64::ONE,1e-12);
    close(out.output[0],C64::from_re(2.),1e-12);
}

#[test]
fn g3_g5_complex_batches_gauge_changes_and_zero_inputs_preserve_the_same_physics() {
    let (model, p) = mechanical();
    let prepared = model.prepare_pole_response(&p,options(26.),||false).unwrap();
    let a=[C64::new(1.,0.2)]; let b=[C64::new(-0.3,0.9)];
    let sum=[a[0]+b[0]];let zero=[C64::ZERO];
    let batch=prepared.port_residues(&[&a,&b,&sum,&zero],||false).unwrap();
    for i in 0..2 { close(batch[2].state[i],batch[0].state[i]+batch[1].state[i],1e-12); }
    assert!(batch[3].state.iter().all(|z|z.abs()==0.));
    let replay=prepared.port_residues(&[&a],||false).unwrap();
    assert_eq!(batch[0].state,replay[0].state);
    let mut scaled=p.clone();
    for x in &mut scaled.state { *x=*x*C64::new(-7.,3.); }
    let same=model.prepare_pole_response(&scaled,options(26.),||false).unwrap();
    let other=same.port_residues(&[&a],||false).unwrap();
    for i in 0..2 { close(other[0].state[i],batch[0].state[i],1e-11); }
}

#[test]
fn g0_g4_wrong_model_budgets_and_late_cancellation_cannot_publish_a_partial_batch() {
    let (model,p)=mechanical();let opt=options(26.);
    let ac=model.factor(26.,||false).unwrap();
    let before=ac.response(&[C64::ONE],||false).unwrap();
    let prepared=model.prepare_pole_response(&p,opt,||false).unwrap();
    let f=[C64::ONE];let count=Cell::new(0);
    let baseline=prepared.port_residues(&[&f,&f],||{count.set(count.get()+1);false}).unwrap();
    let total=count.get();
    for cut in [1,total/2,total] {
        let hit=Cell::new(0);
        assert!(prepared.port_residues(&[&f,&f],||{hit.set(hit.get()+1);hit.get()>=cut}).is_err());
    }
    assert_eq!(baseline[0].state,prepared.port_residues(&[&f],||false).unwrap()[0].state);
    assert_eq!(before.state,ac.response(&[&C64::ONE].map(|x|*x),||false).unwrap().state);
    let other=fs_phs::mass_spring_damper(2.5,1800.,0.6).unwrap();
    let other=HarmonicModel::from_ode(&other,budget(),||false).unwrap();
    assert!(other.prepare_pole_response(&p,opt,||false).is_err());
    let mut small=opt;small.max_physical_bytes=0;
    assert!(model.prepare_pole_response(&p,small,||false).is_err());
    small=opt;small.max_physical_work=0;
    assert!(model.prepare_pole_response(&p,small,||false).is_err());
    assert!(model.prepare_pole_response(&p,opt,||true).is_err());
    assert!(prepared.port_residues(&[],||false).is_err());
    assert!(prepared.port_residues(&[&[]],||false).is_err());
    assert!(prepared.port_residues(&[&[C64::new(f64::NAN,0.)]],||false).is_err());
}
