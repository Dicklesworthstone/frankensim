//! Actual shared-support storage and time owner; no fitted sound or mock dynamics.
use super::*;
use super::super::{ImpactConfig, ImpactSystem};
use fs_exec::CancelGate;

fn rail() -> TranslatingSupport {
    TranslatingSupport { mass_kg: 0.01, stiffness_n_m: 12.0, damping_n_s_m: 0.02,
        initial_position_m: 0.001, initial_velocity_m_s: 0.07,
        maximum_travel_m: 0.02, maximum_slope: 0.2 }
}
fn wires() -> [WireSpan; 2] {
    [3, 2].map(|n| WireSpan { endpoints_m: [[0.0, 0.0], [0.2, 0.0]],
        linear_density_kg_m: 0.01, tension_n: 0.8, bending_stiffness_n_m2: 1e-7,
        damping_per_s: vec![2.0; n] })
}
fn law(ea: f64) -> Option<StringStretching> {
    Some(StringStretching { axial_rigidity_n: ea, maximum_slope: 0.15 })
}
fn close(a: f64, b: f64) {
    assert!((a-b).abs() <= 2e-11 * (a.abs()+b.abs()).max(1e-20), "{a} != {b}");
}

#[test]
fn each_wire_keeps_its_own_stretching_and_uniform_inputs_keep_legacy_arithmetic() {
    let wires = wires();
    let old = SupportedStrings::new(&wires, law(100.0), rail()).unwrap();
    let same = SupportedStrings::with_individual_stretching(&wires, &[law(100.0); 2], rail()).unwrap();
    assert_eq!(old.initial_state(), same.initial_state());
    assert_eq!(old.clock_bounds(), same.clock_bounds());
    let mut z = old.body().initial.iter().map(|s| s.displacement_m_sqrt_kg).collect::<Vec<_>>();
    z[0] += 1e-5; z[3] += 2e-5;
    assert_eq!(old.potential(&z), same.potential(&z));
    let changed = SupportedStrings::with_individual_stretching(&wires,
        &[law(100.0), law(700.0)], rail()).unwrap();
    assert_eq!(changed.0.translation, old.0.translation);
    assert_eq!(changed.total_mass_kg(), old.total_mass_kg());
    assert_eq!(changed.clock_bounds(), old.clock_bounds()); // EA is only the nonlinear channel.
    assert_eq!(changed.string_observation(0, &z).unwrap(), old.string_observation(0, &z).unwrap());
    let a = old.string_observation(1, &z).unwrap();
    let b = changed.string_observation(1, &z).unwrap();
    close(b.stretching_energy_j, 7.0*a.stretching_energy_j);
    close(b.tension_n-wires[1].tension_n, 7.0*(a.tension_n-wires[1].tension_n));
    assert_eq!(a.additional_strain, b.additional_strain);
    let linear = SupportedStrings::with_individual_stretching(&wires, &[law(100.0), None], rail()).unwrap();
    assert_eq!(linear.string_observation(1, &z).unwrap().stretching_energy_j, 0.0);
    assert!(SupportedStrings::with_individual_stretching(&wires, &[law(100.0)], rail()).is_err());
    assert!(SupportedStrings::with_individual_stretching(&wires, &[law(100.0), law(f64::NAN)], rail()).is_err());
}

#[test]
fn initial_relative_motion_keeps_physical_kinetic_cross_terms_and_point_work() {
    let wires = wires(); let rail = rail();
    let s = SupportedStrings::with_individual_stretching(&wires, &[law(100.0),law(700.0)], rail).unwrap();
    let relative: Vec<_> = (0..5).map(|i| ModalAcousticState {
        displacement_m_sqrt_kg: 1e-6*(i as f64-1.2),
        velocity_m_sqrt_kg_per_s: 0.001*(i as f64+0.2),
    }).collect();
    let body = s.body_with_relative_states(&relative).unwrap();
    let kinetic = body.initial.iter().map(|q| 0.5*q.velocity_m_sqrt_kg_per_s.powi(2)).sum::<f64>();
    let physical = 0.5*s.total_mass_kg()*rail.initial_velocity_m_s.powi(2)
        + s.0.translation.iter().zip(&relative).map(|(c,q)|
            c*rail.initial_velocity_m_s*q.velocity_m_sqrt_kg_per_s
            + 0.5*q.velocity_m_sqrt_kg_per_s.powi(2)).sum::<f64>();
    close(kinetic, physical);
    for wire in 0..2 { for station in [0.0, 0.37, 1.0] {
        let b = s.absolute_shape(wire, station).unwrap();
        let position = b.iter().zip(&body.initial).map(|(b,q)| b*q.displacement_m_sqrt_kg).sum::<f64>();
        let velocity = b.iter().zip(&body.initial).map(|(b,q)| b*q.velocity_m_sqrt_kg_per_s).sum::<f64>();
        let (mut x, mut v) = (rail.initial_position_m, rail.initial_velocity_m_s);
        if station != 0.0 && station != 1.0 {
            let shape = fs_dcontact::string_collocation(0.2,0.01,&[station*0.2],wires[wire].damping_per_s.len()).unwrap();
            for (phi,q) in shape.iter().zip(&relative[s.0.ranges[wire].clone()]) {
                x += phi*q.displacement_m_sqrt_kg; v += phi*q.velocity_m_sqrt_kg_per_s;
            }
        }
        close(position,x); close(velocity,v);
        let force = -0.03;
        close(b.iter().zip(&body.initial).map(|(b,q)| force*b*q.velocity_m_sqrt_kg_per_s).sum(), force*v);
    }}
    assert!(s.absolute_shape(2,0.5).is_err()); assert!(s.absolute_shape(0,1.1).is_err());
    assert!(s.body_with_relative_states(&relative[..4]).is_err());
    let mut invalid = relative.clone();invalid[0].displacement_m_sqrt_kg=1.0;
    assert!(s.body_with_relative_states(&invalid).is_err());
}

#[test]
fn heterogeneous_nonlinear_support_gradient_and_damping_keep_the_same_energy() {
    let s = SupportedStrings::with_individual_stretching(&wires(), &[law(100.0),law(700.0)], rail()).unwrap();
    let n=s.mode_count();
    let z:Vec<_>=(0..n).map(|i|1e-5*(i as f64-1.2)).collect();
    let d:Vec<_>=(0..n).map(|i|0.003*(i as f64+0.7)).collect();
    let (mut g,mut hv)=(vec![0.0;n],vec![0.0;n]);s.gradient(&z,&mut g);s.hessian_vector(&z,&d,&mut hv);
    for i in 0..n {
        let h=1e-8;let mut a=z.clone();let mut b=z.clone();a[i]+=h;b[i]-=h;
        let fd=(s.potential(&a)-s.potential(&b))/(2.0*h);
        assert!((g[i]-fd).abs()<1e-6*g[i].abs().max(1e-5));
    }
    let h=1e-6;let a:Vec<_>=z.iter().zip(&d).map(|(z,d)|z+h*d).collect();
    let b:Vec<_>=z.iter().zip(&d).map(|(z,d)|z-h*d).collect();
    let (mut ga,mut gb)=(vec![0.0;n],vec![0.0;n]);s.gradient(&a,&mut ga);s.gradient(&b,&mut gb);
    for i in 0..n {assert!((hv[i]-(ga[i]-gb[i])/(2.0*h)).abs()<1e-7*hv[i].abs().max(1.0));}
    let dim=2*n;let mut r=vec![0.0;dim*dim];s.add_resistance(0,dim,&mut r).unwrap();
    let mut power=0.0;for i in 0..n {for j in 0..n {power+=d[i]*r[(2*i+1)*dim+2*j+1]*d[j];}}
    let velocity=d[n-1]*s.force_weight();
    let expected=rail().damping_n_s_m*velocity*velocity+(0..n-1)
        .map(|i|s.0.damping[i]*(d[i]-s.0.translation[i]*velocity).powi(2)).sum::<f64>();
    close(power,expected);
}

#[test]
fn a_material_only_stretching_change_reacts_on_the_shared_support_and_other_string() {
    let mut ends = Vec::new();
    for ea in [100.0,700.0] {
        let s = SupportedStrings::with_individual_stretching(&wires(), &[law(100.0),law(ea)],
            TranslatingSupport { initial_position_m:0.0,initial_velocity_m_s:0.0,..rail() }).unwrap();
        let mut initial=vec![ModalAcousticState::default();5];initial[3].displacement_m_sqrt_kg=2e-5;
        let body=s.body_with_relative_states(&initial).unwrap();
        let config=ImpactConfig {dt_s:2e-5,max_steps:400,maximum_energy_j:1.0,
            energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-7,maximum_generalized_force:1e4};
        let system=ImpactSystem::new(vec![body],vec![],vec![],vec![],config).unwrap();
        let start=system.stored_energy_j();let mut net=0.0;
        let mut runtime=system.prepare_analytic().unwrap();let gate=CancelGate::new_clock_free();
        for _ in 0..400 {
            let f=runtime.step(&vec![0.0;s.mode_count()],&gate).unwrap();
            net+=f.supplied_work_j-f.dissipated_energy_j;
            assert!((f.stored_energy_j-start-net).abs()<1e-7);
        }
        assert!(runtime.support_observation(0).unwrap().position_m.abs()>1e-10);
        ends.push(runtime.state().to_vec());
    }
    assert_ne!(ends[0][0],ends[1][0]); // Other string, not an output gain.
    assert_ne!(ends[0][10],ends[1][10]); // Shared support coordinate.
}
