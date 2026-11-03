//! Original constrained impedance composition: G1 analytic dynamics, G3 parity.
use fs_math::c64::C64;
use fs_phs::harmonic::{HarmonicBudget, HarmonicModel, PoleOptions};
use fs_phs::{
    ConstitutiveDiscreteGradient, DescriptorPortHamiltonian, PortHamiltonian,
    QuadraticStorage, ResistanceWitness, Storage, mass_spring_damper, step,
    step_descriptor,
};

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 32, max_dense_bytes: 1 << 20, max_dense_work: 1_000_000,
        relative_residual_tolerance: 2e-12,
    }
}
fn close(actual: C64, expected: C64) {
    assert!((actual - expected).abs() < 2e-10 * (1.0 + expected.abs()),
        "{actual:?} != {expected:?}");
}
fn inertance(mass: f64, drag: f64) -> DescriptorPortHamiltonian {
    // p'=-drag*p/mass+lambda; 0=-p/mass+v; y=lambda.
    DescriptorPortHamiltonian::new(
        2, 1, 1, vec![0.0,1.0,-1.0,0.0], vec![drag,0.0,0.0,0.0],
        vec![0.0,1.0], Box::new(QuadraticStorage::new(vec![1.0/mass],1).unwrap()),
    ).unwrap()
}

#[test]
fn g1_original_constraint_reaction_has_output_but_no_stored_energy() {
    let source = inertance(3.0, 0.7);
    assert_eq!(source.output(&[0.6, -2.5]), vec![-2.5]);
    assert_eq!(source.output(&[0.6, 4.0]), vec![4.0]);
    assert_eq!(source.hamiltonian(&[0.6, -2.5]), source.hamiltonian(&[0.6, 4.0]));
    let mut gradient = vec![f64::NAN; 2];
    source.storage().gradient(&[0.6, 4.0], &mut gradient);
    assert_eq!(gradient, vec![0.6/3.0, 0.0]);
}

#[test]
fn g1_constrained_inertance_matches_analytic_forcing_poles_and_power() {
    let (mass, stiffness, drag, load_mass, load_drag) = (2.0,18.0,0.4,3.0,0.7);
    for b in [0.5, -0.5] {
        let system = mass_spring_damper(mass, stiffness, drag).unwrap()
            .with_parallel_descriptor_load_columns(inertance(load_mass,load_drag), &[0.0,b])
            .unwrap();
        assert_eq!(system.state_dim(),4);
        assert_eq!(system.differential_dim(),3);
        assert_eq!(system.port_dim(),1);
        let model = HarmonicModel::from_descriptor(&system,budget(),||false).unwrap();
        let effective_mass = mass+b*b*load_mass;
        let effective_drag = drag+b*b*load_drag;
        for omega in [0.2, 1.0, 3.0, 6.0] {
            let s = C64::new(0.0,-omega);
            let force = C64::new(0.3,-0.2);
            let q = force/(C64::from_re(stiffness)+(s*s).scale(effective_mass)+s.scale(effective_drag));
            let v = s*q;
            let load_flow = v.scale(b);
            let reaction = (s.scale(load_mass)+C64::from_re(load_drag))*load_flow;
            let response = model.factor(omega,||false).unwrap().response(&[force],||false).unwrap();
            for (actual,expected) in response.state.iter().zip([
                q,v.scale(mass),load_flow.scale(load_mass),reaction,
            ]) { close(*actual,expected); }
            close(response.output[0],v);
            close(response.effort[3],reaction);
            close(C64::from_re(response.mean_stored_energy_j),
                C64::from_re(0.25*(stiffness*q.norm_sq()+effective_mass*v.norm_sq())));
            close(C64::from_re(response.dissipated_power_w),
                C64::from_re(0.5*effective_drag*v.norm_sq()));
            assert!(response.power_residual_w.abs()<1e-12);
        }
        let gamma=effective_drag/(2.0*effective_mass);
        let expected=C64::new(-gamma,-(stiffness/effective_mass-gamma*gamma).sqrt());
        let seed=expected+C64::new(0.02,0.03);
        let seed_flow=seed.scale(b);
        let pole=model.resonance_pole(seed,&[
            C64::ONE,seed.scale(mass),seed_flow.scale(load_mass),
            (seed.scale(load_mass)+C64::from_re(load_drag))*seed_flow,
        ],&PoleOptions {
            spectral_scale:3.0,step_tolerance:1e-13,residual_tolerance:1e-13,
            max_iterations:40,max_dimension:32,max_dense_bytes:1<<20,max_work:1_000_000,
        },||false).unwrap();
        close(pole.laplace_pole_per_s,expected);
        assert!(pole.decay_power_residual_w.abs()<1e-10*pole.dissipated_power_w);
    }
}

#[test]
fn g1_pure_algebraic_mutual_resistance_retains_signed_nonorthogonal_ports() {
    let host=PortHamiltonian::new(
        2,2,vec![0.0;4],vec![0.0;4],vec![1.0,0.0,0.0,1.0],
        Box::new(QuadraticStorage::new(vec![1.0,0.0,0.0,1.0],2).unwrap()),
    ).unwrap();
    let load=DescriptorPortHamiltonian::new(
        2,0,2,vec![0.0;4],vec![2.0,0.5,0.5,1.0],vec![1.0,0.0,0.0,1.0],
        Box::new(QuadraticStorage::new(vec![],0).unwrap()),
    ).unwrap();
    let b=[1.0,-0.5,0.4,0.8];
    let system=host.with_parallel_descriptor_load_columns(load,&b).unwrap();
    assert_eq!(system.differential_dim(),2);
    let omega=3.2;
    let s=C64::new(0.0,-omega);
    let force=[C64::new(1.0,0.2),C64::new(-0.4,0.7)];
    let inverse=[1.0/1.75,-0.5/1.75,-0.5/1.75,2.0/1.75];
    let mut z=[0.0;4];
    for i in 0..2 { for j in 0..2 { for a in 0..2 { for c in 0..2 {
        z[i*2+j]+=b[i*2+a]*inverse[a*2+c]*b[j*2+c];
    }}}}
    let a=s+C64::from_re(z[0]); let d=s+C64::from_re(z[3]);
    let determinant=a*d-C64::from_re(z[1]*z[2]);
    let velocity=[
        (d*force[0]-force[1].scale(z[1]))/determinant,
        (a*force[1]-force[0].scale(z[2]))/determinant,
    ];
    let flow=[
        velocity[0].scale(b[0])+velocity[1].scale(b[2]),
        velocity[0].scale(b[1])+velocity[1].scale(b[3]),
    ];
    let reactions=[
        flow[0].scale(inverse[0])+flow[1].scale(inverse[1]),
        flow[0].scale(inverse[2])+flow[1].scale(inverse[3]),
    ];
    let model=HarmonicModel::from_descriptor(&system,budget(),||false).unwrap();
    let response=model.factor(omega,||false).unwrap().response(&force,||false).unwrap();
    for i in 0..2 { close(response.state[i],velocity[i]); close(response.state[2+i],reactions[i]); }
    close(C64::from_re(response.mean_stored_energy_j),
        C64::from_re(0.25*(velocity[0].norm_sq()+velocity[1].norm_sq())));
    let heat=0.5*((reactions[0].conj()*flow[0]).re+(reactions[1].conj()*flow[1]).re);
    close(C64::from_re(response.dissipated_power_w),C64::from_re(heat));
    assert!(response.power_residual_w.abs()<1e-12);
}

struct Quartic { mass:f64, k:f64 }
impl Storage for Quartic {
    fn hamiltonian(&self,x:&[f64])->f64 {0.25*self.k*x[0].powi(4)+0.5*x[1]*x[1]/self.mass}
    fn gradient(&self,x:&[f64],out:&mut[f64]) {out[0]=self.k*x[0].powi(3);out[1]=x[1]/self.mass;}
    fn coordinate_scale(&self,coordinate:usize)->Option<f64> {
        match coordinate {0=>Some(0.3),1=>Some(self.mass),_=>None}
    }
    fn discrete_gradient_rule(&self)->Option<&dyn ConstitutiveDiscreteGradient> {Some(self)}
}
impl ConstitutiveDiscreteGradient for Quartic {
    fn evaluate(&self,a:&[f64],b:&[f64],_:&[f64],out:&mut[f64]) {
        out[0]=0.25*self.k*(a[0]+b[0])*(a[0]*a[0]+b[0]*b[0]);
        out[1]=f64::midpoint(a[1],b[1])/self.mass;
    }
    fn end_jacobian_action(&self,a:&[f64],b:&[f64],d:&[f64],out:&mut[f64])->bool {
        out[0]=0.25*self.k*(a[0]*a[0]+2.0*a[0]*b[0]+3.0*b[0]*b[0])*d[0];
        out[1]=0.5*d[1]/self.mass;
        true
    }
}
fn quartic(mass:f64,k:f64,drag:f64)->PortHamiltonian {
    PortHamiltonian::new(2,1,vec![0.0,1.0,-1.0,0.0],vec![0.0,0.0,0.0,drag],
        vec![0.0,1.0],Box::new(Quartic{mass,k})).unwrap()
}

#[test]
fn g3_nonlinear_constitutive_callbacks_survive_constrained_composition() {
    let (mass,load_mass,drag,load_drag,k,b)=(2.0,3.0,0.2,0.7,15.0,-0.5);
    let system=quartic(mass,k,drag)
        .with_parallel_descriptor_load_columns(inertance(load_mass,load_drag),&[0.0,b]).unwrap();
    let equivalent=quartic(mass+b*b*load_mass,k,drag+b*b*load_drag);
    assert!(system.storage().discrete_gradient_rule().is_some());
    assert_eq!(system.storage().coordinate_scale(0),Some(0.3));
    assert_eq!(system.storage().coordinate_scale(1),Some(mass));
    assert_eq!(system.storage().coordinate_scale(3),None);
    assert!(HarmonicModel::from_descriptor(&system,budget(),||false).is_err());
    let mut x=vec![0.3,mass*0.2,load_mass*b*0.2,0.0];
    let mut reference=vec![0.3,(mass+b*b*load_mass)*0.2];
    let mut changed=x.clone(); changed[3]=1e12;
    assert_eq!(system.hamiltonian(&x),system.hamiltonian(&changed));
    let initial=system.hamiltonian(&x);
    let (mut supply,mut heat)=(0.0,0.0);
    for tick in 0..40 {
        let input=[if tick<20 {0.15}else{0.0}];
        let record=step_descriptor(&system,&x,&input,0.01).unwrap();
        let oracle=step(&equivalent,&reference,&input,0.01).unwrap();
        assert!((record.x[0]-oracle.x[0]).abs()<1e-10);
        assert!((record.x[1]/mass-oracle.x[1]/(mass+b*b*load_mass)).abs()<1e-10);
        assert!((record.x[2]/load_mass-b*record.x[1]/mass).abs()<1e-10);
        supply+=record.supplied;heat+=record.dissipated;
        x=record.x;reference=oracle.x;
    }
    assert!(heat>0.0);
    assert!((system.hamiltonian(&x)+heat-initial-supply).abs()<1e-10);
}

#[test]
fn g3_all_differential_load_preserves_ode_response_and_source_resistance_witness() {
    let make=|| {
        let host=mass_spring_damper(2.0,18.0,0.25).unwrap();
        let r=host.structure().1.to_vec();
        host.with_resistance_witness(ResistanceWitness::try_from_source_rows(
            2,r,vec![0.0,0.5],
        ).unwrap()).unwrap()
    };
    let ode=mass_spring_damper(3.0,7.0,0.5).unwrap();
    let desc=DescriptorPortHamiltonian::new(
        2,2,1,ode.structure().0.to_vec(),ode.structure().1.to_vec(),ode.structure().2.to_vec(),
        Box::new(QuadraticStorage::new(vec![7.0,0.0,0.0,1.0/3.0],2).unwrap()),
    ).unwrap();
    let ordinary=make().with_parallel_load_columns(ode,&[0.0,-0.5]).unwrap();
    assert!(ordinary.resistance_witness().is_some());
    let descriptor=make().with_parallel_descriptor_load_columns(desc,&[0.0,-0.5]).unwrap();
    let a=HarmonicModel::from_ode(&ordinary,budget(),||false).unwrap();
    let b=HarmonicModel::from_descriptor(&descriptor,budget(),||false).unwrap();
    for omega in [0.1,2.0,12.0] {
        let input=[C64::new(1.0,-0.2)];
        let ra=a.factor(omega,||false).unwrap().response(&input,||false).unwrap();
        let rb=b.factor(omega,||false).unwrap().response(&input,||false).unwrap();
        assert_eq!(ra,rb);
    }
    for columns in [vec![],vec![0.0],vec![0.0,f64::NAN]] {
        assert!(make().with_parallel_descriptor_load_columns(inertance(3.0,0.5),&columns).is_err());
    }
}

