//! Independent circuit, scattering and work checks, not material identification.
use fs_math::c64::C64;
use fs_vfit::impedance::{ImpedanceState, SeriesImpedanceSpec};
use fs_vfit::relaxation::{ComplianceRelaxationTerm as RcTerm, ParallelCell, ParallelCellLoad,
    ParallelCellSpec, RelaxationImpedance, RelaxationImpedanceSpec, RelaxationTerm, MAX_PARALLEL_CELLS,
    MAX_RELAXATION_BRANCHES};
use fs_vfit::waveguide::network::{NetworkNode, NetworkSegment, WaveguideNetwork};
const DT: f64 = 0.02;
fn base() -> SeriesImpedanceSpec {
    SeriesImpedanceSpec { resistance_pa_s_m3: 0.4, inertance_pa_s2_m3: 0.3, compliance_m3_pa: Some(0.7) }
}
fn rl() -> [RelaxationTerm; 1] { [RelaxationTerm { resistance_pa_s_m3: 0.6, rate_per_s: 2.0 }] }
fn rc() -> [RcTerm; 2] { [RcTerm { resistance_pa_s_m3: 0.8, rate_per_s: 3.0 },
    RcTerm { resistance_pa_s_m3: 0.2, rate_per_s: 40.0 }] }
fn spec() -> RelaxationImpedanceSpec {
    RelaxationImpedanceSpec::new(base(), &rl()).unwrap().with_compliance_terms(&rc()).unwrap()
}
fn model() -> RelaxationImpedance { RelaxationImpedance::new(spec(), 1.0, DT).unwrap() }
fn close(a: f64, b: f64) { assert!((a-b).abs() < 3e-11 * (1.0+a.abs()+b.abs()), "{a} != {b}"); }

#[test]
fn mixed_memories_share_the_actual_flow_and_close_independent_storage() {
    let mut m = model();
    let mut old_base = ImpedanceState::default();
    let mut released = false;
    for n in 0..2000 {
        let old_q = m.branch_flows_m3_s().to_vec();
        let old_p = m.branch_pressures_pa().to_vec();
        let before = m.stored_energy_j();
        let f = m.step(if n < 300 { f64::from(n % 13) - 6.0 } else { 0.0 }).unwrap();
        let s = f.port.state;
        let u = f.port.flow_m3_s;
        let mut pressure = 0.4*u + 0.3*(s.inertive_flow_m3_s-old_base.inertive_flow_m3_s)/DT
            + 0.5*(s.compliance_pressure_pa+old_base.compliance_pressure_pa);
        let mut energy = 0.15*s.inertive_flow_m3_s.powi(2)+0.35*s.compliance_pressure_pa.powi(2);
        let mut loss = 0.4*u*u*DT;
        for (i,t) in rl().iter().enumerate() {
            let q = f.branch_flows_m3_s()[i]; let mid = 0.5*(q+old_q[i]);
            close((q-old_q[i])/DT,t.rate_per_s*(u-mid));
            pressure += t.resistance_pa_s_m3*(u-mid);
            energy += 0.5*t.resistance_pa_s_m3/t.rate_per_s*q*q;
            loss += t.resistance_pa_s_m3*(u-mid).powi(2)*DT;
        }
        for (i,t) in rc().iter().enumerate() {
            let p = f.branch_pressures_pa()[i]; let mid = 0.5*(p+old_p[i]);
            close((p-old_p[i])/DT,t.rate_per_s*(t.resistance_pa_s_m3*u-mid));
            pressure += mid;
            energy += 0.5*p*p/t.resistance_pa_s_m3/t.rate_per_s;
            loss += mid*mid/t.resistance_pa_s_m3*DT;
        }
        close(pressure,f.port.pressure_pa); close(energy,f.port.stored_energy_j);
        close(loss,f.port.dissipated_energy_j);
        let scale = (before+energy+loss+f.port.supplied_work_j.abs()).max(f64::MIN_POSITIVE);
        assert!(f.port.balance_residual_j().abs()<2e-12*scale);
        if n>=300 { assert!(energy<=before+2e-12*scale); released |= f.port.supplied_work_j<0.0; }
        old_base=s;
    }
    assert!(released);
}

#[test]
fn mixed_reflection_matches_independent_complex_impedance_and_warp() {
    for bin in [3,13,53,151,401] {
        let theta = core::f64::consts::TAU*f64::from(bin)/2048.0;
        let mut m = model(); let mut measured = C64::ZERO;
        for n in 0..12288 {
            let phase = theta*f64::from(n);
            let b = m.step(phase.cos()).unwrap().port.reflected_pressure_pa;
            if n>=10240 { measured=measured+C64::new(phase.cos(),-phase.sin()).scale(b/1024.0); }
        }
        let s=C64::new(0.0,2.0/DT*(theta*0.5).tan());
        let mut z=C64::from_re(0.4)+s.scale(0.3)+s.scale(0.7).recip();
        for t in rl() { z=z+s.scale(t.resistance_pa_s_m3)*(s+C64::from_re(t.rate_per_s)).recip(); }
        for t in rc() { z=z+(s+C64::from_re(t.rate_per_s)).recip().scale(t.resistance_pa_s_m3*t.rate_per_s); }
        let expected=(z-C64::ONE)*(z+C64::ONE).recip();
        assert!((measured-expected).abs()<1e-10);
    }
}

#[test]
fn pure_compliance_arm_matches_independent_parallel_rc_recurrence() {
    let zero=SeriesImpedanceSpec{resistance_pa_s_m3:0.0,inertance_pa_s2_m3:0.0,compliance_m3_pa:None};
    let t=rc()[0]; let c=1.0/(t.resistance_pa_s_m3*t.rate_per_s);
    let mut m=RelaxationImpedance::new(RelaxationImpedanceSpec::new(zero,&[]).unwrap()
        .with_compliance_terms(&[t]).unwrap(),1.0,DT).unwrap();
    let mut p=0.0;
    for n in 0..1000 {
        let a=if n<300 {1.0} else {0.0};
        // C p' = 2a - p - p/R, with the unit characteristic impedance.
        let k=0.5*DT/c*(1.0+1.0/t.resistance_pa_s_m3);
        let next=((1.0-k)*p+DT/c*2.0*a)/(1.0+k);
        let f=m.step(a).unwrap();
        close(f.branch_pressures_pa()[0],next);
        close(f.port.pressure_pa,0.5*(p+next));
        p=next;
    }
}

#[test]
fn empty_attachment_and_refused_previews_preserve_all_histories() {
    let original=RelaxationImpedanceSpec::new(base(),&rl()).unwrap();
    let mut a=RelaxationImpedance::new(original,1.0,DT).unwrap();
    let mut b=RelaxationImpedance::new(original.with_compliance_terms(&[]).unwrap(),1.0,DT).unwrap();
    for n in 0..300 { assert_eq!(a.step(f64::from(n).sin()).unwrap(),b.step(f64::from(n).sin()).unwrap()); }
    let mut a=model(); a.step(20.0).unwrap(); let mut b=a;
    let q=a.branch_flows_m3_s().to_vec(); let p=a.branch_pressures_pa().to_vec();
    let expected=a.preview_step(7.0).unwrap();
    for invalid in [f64::NAN,f64::INFINITY,f64::MAX] { assert!(a.step(invalid).is_err()); }
    assert_eq!(a.branch_flows_m3_s(),q); assert_eq!(a.branch_pressures_pa(),p);
    assert_eq!(a.step(7.0).unwrap(),expected); assert_eq!(expected,b.step(7.0).unwrap());
}

#[test]
fn combined_budget_and_nonrepresentable_storage_refuse_without_truncation() {
    let original=RelaxationImpedanceSpec::new(base(),&rl()).unwrap();
    assert!(original.with_compliance_terms(&[rc()[0];MAX_RELAXATION_BRANCHES-1]).is_ok());
    assert!(original.with_compliance_terms(&[rc()[0];MAX_RELAXATION_BRANCHES]).is_err());
    assert!(spec().with_compliance_terms(&rc()).is_err());
    for bad in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        assert!(original.with_compliance_terms(&[RcTerm{resistance_pa_s_m3:bad,..rc()[0]}]).is_err());
        assert!(original.with_compliance_terms(&[RcTerm{rate_per_s:bad,..rc()[0]}]).is_err());
    }
    let tiny=RcTerm{resistance_pa_s_m3:1e-300,rate_per_s:1e-300};
    let tiny=original.with_compliance_terms(&[tiny]).unwrap();
    assert!(RelaxationImpedance::new(tiny,1.0,DT).is_err());
}

#[test]
fn existing_series_and_shunt_networks_publish_rc_memory_with_wave_energy() {
    for kind in [NetworkNode::Series{load:spec()},NetworkNode::Shunt{load:spec()}] {
        let make=||WaveguideNetwork::new(&[NetworkNode::Inlet,kind,NetworkNode::Termination{reflection:0.0}],
            &[NetworkSegment{nodes:[0,1],one_way_samples:3,impedance_pa_s_m3:1.0},
              NetworkSegment{nodes:[1,2],one_way_samples:5,impedance_pa_s_m3:2.0}],DT,1<<20).unwrap();
        let mut a=make(); let mut b=make(); let mut heard=false;
        for n in 0..1200 {
            let before=a.stored_energy_j();
            let old=a.load_compliance_pressures(1).unwrap().to_vec();
            assert!(a.step(f64::NAN).is_err());
            assert_eq!(old,a.load_compliance_pressures(1).unwrap());
            let input=if n<300 {f64::from(n%17)-8.0} else {0.0};
            let f=a.step(input).unwrap(); assert_eq!(f,b.step(input).unwrap());
            let node=a.node_frame(1).unwrap(); let pressures=a.load_compliance_pressures(1).unwrap();
            for (i,t) in rc().iter().enumerate() {
                close((pressures[i]-old[i])/DT,t.rate_per_s*(t.resistance_pa_s_m3*node.load_flow_m3_s
                    -0.5*(pressures[i]+old[i])));
            }
            let scale=(before+f.stored_energy_j+f.inlet_work_j.abs()+f.interior_loss_j+f.terminal_loss_j)
                .max(f64::MIN_POSITIVE);
            assert!(f.balance_residual_j().abs()<2e-12*scale);
            heard |= pressures.iter().any(|p| p.abs()>1e-8);
        }
        assert!(heard);
    }
}

// One RL cell (C = 0) and resonances at 7.1 and 71 rad/s, below the 157 rad/s Nyquist of DT.
fn cells() -> [ParallelCell; 3] {
    [ParallelCell { resistance_pa_s_m3: 0.6, inertance_pa_s2_m3: 0.3, compliance_m3_pa: 0.0 },
     ParallelCell { resistance_pa_s_m3: 0.9, inertance_pa_s2_m3: 0.05, compliance_m3_pa: 0.4 },
     ParallelCell { resistance_pa_s_m3: 2.5, inertance_pa_s2_m3: 0.01, compliance_m3_pa: 0.02 }]
}
fn chain() -> ParallelCellLoad { ParallelCellLoad::new(ParallelCellSpec::new(&cells()).unwrap(), DT).unwrap() }

#[test]
fn parallel_cells_close_their_own_recurrence_storage_and_loss_on_the_driven_flow() {
    let mut m=chain(); let mut released=false;
    for n in 0..2000 {
        let old_i=m.inductor_flows_m3_s().to_vec(); let old_v=m.pressures_pa().to_vec();
        let before=m.stored_energy_j();
        let u=if n<300 {f64::from(n%11)-5.0} else {0.0};
        let (direct,history)=m.affine_pressure();
        let f=m.step(u).unwrap();
        close(f.pressure_pa,direct*u+history);
        let (mut pressure,mut energy,mut loss)=(0.0,0.0,0.0);
        for (k,t) in cells().iter().enumerate() {
            let (i,v)=(f.inductor_flows_m3_s()[k],f.pressures_pa()[k]);
            let mid=if t.compliance_m3_pa>0.0 {0.5*(v+old_v[k])} else {v};
            // L i' = v and C v' = U - v/R - i, both at the step midpoint.
            close(t.inertance_pa_s2_m3*(i-old_i[k])/DT,mid);
            close(t.compliance_m3_pa*(v-old_v[k])/DT,u-mid/t.resistance_pa_s_m3-0.5*(i+old_i[k]));
            pressure+=mid;
            energy+=0.5*t.inertance_pa_s2_m3*i*i+0.5*t.compliance_m3_pa*v*v;
            loss+=mid*mid/t.resistance_pa_s_m3*DT;
        }
        close(pressure,f.pressure_pa); close(energy,f.stored_energy_j); close(loss,f.dissipated_energy_j);
        let scale=(before+energy+loss+f.supplied_work_j.abs()).max(f64::MIN_POSITIVE);
        assert!(f.balance_residual_j().abs()<2e-12*scale);
        // A prescribed zero flow takes no work, so the cells return stored
        // energy to the driver only while it drives; afterwards they decay.
        released |= f.supplied_work_j<0.0;
        if n>=300 { assert!(energy<=before+2e-12*scale); }
    }
    assert!(released && m.stored_energy_j()<1e-12);
}

#[test]
fn a_cell_without_compliance_is_the_existing_high_pass_term_on_the_same_flow() {
    let zero=SeriesImpedanceSpec{resistance_pa_s_m3:0.0,inertance_pa_s2_m3:0.0,compliance_m3_pa:None};
    let term=rl()[0];
    let mut wave=RelaxationImpedance::new(RelaxationImpedanceSpec::new(zero,&[term]).unwrap(),1.0,DT).unwrap();
    let cell=ParallelCell{resistance_pa_s_m3:term.resistance_pa_s_m3,
        inertance_pa_s2_m3:term.resistance_pa_s_m3/term.rate_per_s,compliance_m3_pa:0.0};
    let mut flow=ParallelCellLoad::new(ParallelCellSpec::new(&[cell]).unwrap(),DT).unwrap();
    for n in 0..1500 {
        let f=wave.step(if n<200 {f64::from(n%7)-3.0} else {0.0}).unwrap();
        let g=flow.step(f.port.flow_m3_s).unwrap();
        close(g.pressure_pa,f.port.pressure_pa);
        close(g.inductor_flows_m3_s()[0],f.branch_flows_m3_s()[0]);
        close(g.stored_energy_j,f.port.stored_energy_j); close(g.dissipated_energy_j,f.port.dissipated_energy_j);
    }
}

#[test]
fn parallel_cells_answer_a_driven_flow_with_their_warped_impedance() {
    let spec=ParallelCellSpec::new(&cells()).unwrap();
    for bin in [3,13,53,151,401] {
        let theta=core::f64::consts::TAU*f64::from(bin)/2048.0;
        let mut m=ParallelCellLoad::new(spec,DT).unwrap(); let mut measured=C64::ZERO;
        for n in 0..12288 {
            // The step-mean flow of a sampled cosine is its sample-pair average.
            let mean=0.5*((theta*f64::from(n)).cos()+(theta*f64::from(n+1)).cos());
            let p=m.step(mean).unwrap().pressure_pa;
            if n>=10240 {
                let phase=theta*(f64::from(n)+0.5);
                measured=measured+C64::new(phase.cos(),-phase.sin()).scale(p/1024.0);
            }
        }
        let s=C64::new(0.0,2.0/DT*(theta*0.5).tan());
        // Midpoint pressure over midpoint flow is Z at the warped s; the drive's
        // midpoint phasor at the half-sample phase has amplitude cos(theta/2).
        let expected=spec.impedance(s).scale((theta*0.5).cos());
        assert!((measured-expected).abs()<1e-10*expected.abs().max(1.0), "{bin}: {measured:?} {expected:?}");
    }
}

#[test]
fn parallel_cells_refuse_bad_data_budgets_and_keep_history_on_refusal() {
    assert!(ParallelCellSpec::new(&[]).is_err());
    assert!(ParallelCellSpec::new(&[cells()[1];MAX_PARALLEL_CELLS]).is_ok());
    assert!(ParallelCellSpec::new(&[cells()[1];MAX_PARALLEL_CELLS+1]).is_err());
    for bad in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        assert!(ParallelCellSpec::new(&[ParallelCell{resistance_pa_s_m3:bad,..cells()[1]}]).is_err());
        assert!(ParallelCellSpec::new(&[ParallelCell{inertance_pa_s2_m3:bad,..cells()[1]}]).is_err());
    }
    for bad in [-1.0,f64::NAN,f64::INFINITY] {
        assert!(ParallelCellSpec::new(&[ParallelCell{compliance_m3_pa:bad,..cells()[1]}]).is_err());
    }
    let tiny=ParallelCellSpec::new(&[ParallelCell{inertance_pa_s2_m3:1e-320,..cells()[1]}]).unwrap();
    assert!(ParallelCellLoad::new(tiny,DT).is_err());
    assert!(ParallelCellLoad::new(ParallelCellSpec::new(&cells()).unwrap(),0.0).is_err());
    let mut m=chain(); m.step(3.0).unwrap();
    let (i,v,e)=(m.inductor_flows_m3_s().to_vec(),m.pressures_pa().to_vec(),m.stored_energy_j());
    let expected=m.preview_step(2.0).unwrap();
    for invalid in [f64::NAN,f64::INFINITY] { assert!(m.step(invalid).is_err()); }
    assert_eq!((m.inductor_flows_m3_s().to_vec(),m.pressures_pa().to_vec(),m.stored_energy_j()),(i,v,e));
    assert_eq!(m.step(2.0).unwrap(),expected);
    close(ParallelCellSpec::new(&cells()).unwrap().static_inertance_pa_s2_m3(),0.36);
}
