use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_couple::render::plate::impact::{BodyPotential,ImpactBody,ImpactConfig,ImpactError,ImpactSystem,VolumeSpring};
use fs_couple::render::plate::impact::cavity::CavityCoupling;
use fs_couple::render::plate::impact::felt::{FeltPad,KelvinBranch};
use fs_couple::vibroacoustic::{CavityModes,StructuralModes,VibroacousticModel};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_material::fiber::WoolFelt;
use fs_math::c64::C64;
use fs_phs::{PortHamiltonian,QuadraticStorage,StepWorkspace};

fn config() -> ImpactConfig { ImpactConfig {dt_s:2e-6,max_steps:2000,maximum_energy_j:20.0,
    energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-7,maximum_generalized_force:10000.0} }
fn body(w:f64,q:f64,v:f64) -> ImpactBody { ImpactBody {potential:BodyPotential::Linear(vec![w]),
    damping_per_s:vec![0.0],initial:vec![ModalAcousticState {displacement_m_sqrt_kg:q,velocity_m_sqrt_kg_per_s:v}]} }
fn cavity(omegas:Vec<f64>) -> CavityModes {let n=omegas.len(); CavityModes {omegas,
    lambdas:vec![0.01;n],interface:vec![vec![1.0];n],loss_factor:0.0,rho0:1.2,c0:343.0} }

#[test]
fn uniform_mode_is_the_original_sealed_volume_without_extra_inertia() {
    let basis=cavity(vec![0.0]); let gate=CancelGate::new_clock_free();
    let make=|| vec![body(800.0,1e-5,0.0),body(900.0,0.0,0.0)];
    let spring=VolumeSpring {bulk_modulus_pa:basis.rho0*basis.c0*basis.c0,volume_m3:0.01,areas:vec![0.1,-0.1]};
    let mut old=ImpactSystem::new(make(),vec![],vec![],vec![spring],config()).unwrap().prepare().unwrap();
    let compiled=CavityCoupling::new(&basis,2,&[0.1,-0.1],&[0.0]).unwrap();
    assert_eq!(compiled.total_modes(),2);
    let (system,probe)=compiled.build(make(),vec![],vec![],config(),&gate).unwrap();
    let mut new=system.prepare().unwrap();
    assert!(probe.pressure_at(new.state(),&[1.0]).unwrap()<0.0,"outward motion rarefies the gas");
    for _ in 0..200 {
        let a=old.step(&[0.0;2],&gate).unwrap(); let b=new.step(&[0.0;2],&gate).unwrap();
        assert_eq!(old.state(),new.state()); assert_eq!(a.stored_energy_j.to_bits(),b.stored_energy_j.to_bits());
    }
}

#[test]
fn standing_wave_matches_the_existing_frequency_domain_coupling_and_midpoint_dynamics() {
    let (ws,wa,c)=(800.0_f64,1000.0_f64,0.1_f64);
    let basis=cavity(vec![wa]); let scale=basis.rho0*basis.c0*basis.c0/basis.lambdas[0];
    let cross=wa*scale.sqrt()*c;
    let k00=ws*ws+scale*c*c;
    let model=VibroacousticModel::try_new(&StructuralModes {omegas:vec![ws],shapes:vec![vec![1.0]],loss_factor:0.0},
        &basis,vec![c],None).unwrap();
    let roots=model.undamped_natural_frequencies().unwrap();
    let trace=k00+wa*wa; let det=ws*ws*wa*wa;
    let expected=[((trace-(trace*trace-4.0*det).sqrt())/2.0).sqrt(),
        ((trace+(trace*trace-4.0*det).sqrt())/2.0).sqrt()];
    for (&a,&b) in roots.iter().zip(&expected) {assert!((a-b).abs()<1e-6*b);}
    for w in [300.0_f64,900.0,1400.0] {
        let response=model.frf(w,&[C64::ONE]).unwrap();
        let q=1.0/(k00-w*w-cross*cross/(wa*wa-w*w));
        let z=-cross*q/(wa*wa-w*w);
        let pressure=-scale*(c*q+wa/scale.sqrt()*z);
        assert!((response.b[0].re-q).abs()<1e-9*q.abs());
        assert!((response.a[0].re-pressure).abs()<1e-9*pressure.abs());
    }
    // Independently assembled quadratic Hamiltonian, not the volume adapter.
    let mut q=vec![0.0;16];q[0]=k00;q[2]=cross;q[8]=cross;q[10]=wa*wa;q[5]=1.0;q[15]=1.0;
    let mut j=vec![0.0;16];j[1]=1.0;j[4]=-1.0;j[11]=1.0;j[14]=-1.0;
    let reference=PortHamiltonian::new(4,0,j,vec![0.0;16],vec![],Box::new(QuadraticStorage::new(q,4).unwrap())).unwrap();
    let gate=CancelGate::new_clock_free();
    let (system,probe)=CavityCoupling::new(&basis,1,&[c],&[0.0]).unwrap()
        .build(vec![body(ws,1e-5,0.0)],vec![],vec![],config(),&gate).unwrap();
    let mut actual=system.prepare().unwrap();let h0=actual.stored_energy_j();
    let mut state=[1e-5,0.0,0.0,0.0];let mut next=[0.0;4];
    let mut workspace=StepWorkspace::new(&reference).unwrap();let mut acoustic_peak=0.0_f64;
    for _ in 0..400 {
        workspace.step_into(&reference,&state,&[],config().dt_s,&mut next,&mut []).unwrap();
        actual.step(&[0.0;2],&gate).unwrap();
        for (&a,&b) in actual.state().iter().zip(&next) {assert!((a-b).abs()<1e-9);}
        acoustic_peak=acoustic_peak.max(actual.state()[3].abs());
        assert!(probe.pressure_at(actual.state(),&[1.0]).unwrap().is_finite());state=next;
    }
    assert!(acoustic_peak>1e-6);
    assert!((actual.stored_energy_j()-h0).abs()<1e-9);
}

#[test]
fn basis_rescaling_changes_coefficients_but_not_physical_motion_or_pressure() {
    let gate=CancelGate::new_clock_free();let a=cavity(vec![0.0,1000.0]);let mut b=a.clone();
    for l in &mut b.lambdas {*l*=9.0;}for row in &mut b.interface {for x in row {*x*=3.0;}}
    let build=|basis:&CavityModes,c:&[f64]| CavityCoupling::new(basis,1,c,&[0.0,5.0]).unwrap()
        .build(vec![body(800.0,1e-5,0.02)],vec![],vec![],config(),&gate).unwrap();
    let (sa,pa)=build(&a,&[0.1,0.07]);let(sb,pb)=build(&b,&[0.3,0.21]);
    let mut sa=sa.prepare().unwrap();let mut sb=sb.prepare().unwrap();
    for _ in 0..200 {
        sa.step(&[0.0;2],&gate).unwrap();sb.step(&[0.0;2],&gate).unwrap();
        for (&x,&y) in sa.state().iter().zip(sb.state()) {assert!((x-y).abs()<1e-9);}
        let x=pa.pressure_at(sa.state(),&[1.0,0.4]).unwrap();let y=pb.pressure_at(sb.state(),&[3.0,1.2]).unwrap();
        assert!((x-y).abs()<1e-5);
    }
}

#[test]
fn physical_contact_and_felt_continue_to_address_solids_not_appended_gas_modes() {
    let gate=CancelGate::new_clock_free();let basis=cavity(vec![0.0,1800.0]);
    let (stick,w)=ImpactBody::free_mass(0.02,-0.0002,0.8).unwrap();
    let contact=Obstacle::new(vec![w,-2.0],1,2,vec![0.0],vec![1.0],2e7,1.5,"synthetic cavity impact".into()).unwrap();
    let pad=FeltPad {area_m2:0.001,thickness_m:0.006,precompression_m:0.0006,weights:vec![0.0,2.0],
        law:WoolFelt::new(30000.0,0.2,2.2,3.0,0.15,0.7).unwrap(),prior_maximum_strain:0.1,
        creep:vec![KelvinBranch {stiffness_n_m:1500.0,viscosity_n_s_m:8.0}]};
    let (system,probe)=CavityCoupling::new(&basis,2,&[0.0,0.0,0.03,0.02],&[0.0,5.0]).unwrap()
        .build(vec![stick,body(1000.0,0.0,0.0)],vec![contact],vec![pad],config(),&gate).unwrap();
    assert_eq!(system.state().len(),7,"q/p solid + q/p gas + Kelvin history");
    let mut system=system.prepare().unwrap();let h0=system.stored_energy_j();let mut loss=0.0;
    for _ in 0..1200 {let f=system.step(&[0.0;3],&gate).unwrap();loss+=f.dissipated_energy_j;}
    assert!(system.state()[1]<0.0);assert_ne!(system.state()[3],0.0);assert_ne!(system.state()[5],0.0);
    assert_ne!(system.state()[6],0.0);assert!((system.stored_energy_j()+loss-h0).abs()<1e-7);
    assert!(probe.pressure_at(system.state(),&[1.0,0.2]).unwrap().is_finite());
}

#[test]
fn invalid_frequency_loss_layout_and_cancelled_builds_are_refused() {
    let mut basis=cavity(vec![1000.0]);basis.loss_factor=0.1;
    assert!(CavityCoupling::new(&basis,1,&[0.1],&[0.0]).is_err());basis.loss_factor=0.0;
    assert!(CavityCoupling::new(&basis,1,&[],&[0.0]).is_err());
    let compiled=CavityCoupling::new(&basis,1,&[0.1],&[0.0]).unwrap();
    let gate=CancelGate::new_clock_free();let mut cfg=config();cfg.dt_s=0.01;
    assert!(compiled.clone().build(vec![body(1.0,0.0,0.0)],vec![],vec![],cfg,&gate).is_err());
    gate.request();assert!(matches!(compiled.clone().build(vec![body(1.0,0.0,0.0)],vec![],vec![],config(),&gate),Err(ImpactError::Cancelled)));
    let mut output=[123.0];assert!(compiled.pressures_into(&[f64::NAN;4],&mut output).is_err());assert_eq!(output,[123.0]);
    let zero=cavity(vec![0.0]);assert!(CavityCoupling::new(&zero,1,&[0.1],&[1.0]).is_err());
}

fn neck(averages:Vec<f64>,resistance:f64) -> fs_couple::render::plate::impact::cavity::neck::CavityNeck {
    fs_couple::render::plate::impact::cavity::neck::CavityNeck {
        area_m2:1e-4,effective_length_m:0.01,resistance_pa_s_m3:resistance,
        pressure_shape_averages:averages,initial_volume_m3:1e-6,initial_flow_m3_s:0.0,
    }
}

#[test]
fn g1_compact_neck_has_helmholtz_frequency_not_an_authored_oscillator() {
    let basis=cavity(vec![0.0]);let gate=CancelGate::new_clock_free();
    let opening=neck(vec![1.0],0.0);
    let omega=basis.c0*(opening.area_m2/(basis.lambdas[0]*opening.effective_length_m)).sqrt();
    let initial=opening.initial_volume_m3;
    let compiled=CavityCoupling::new(&basis,1,&[0.0],&[0.0]).unwrap()
        .with_necks(vec![opening],&gate).unwrap();
    assert_eq!(compiled.structural_modes(),1);assert_eq!(compiled.total_modes(),2);
    assert_eq!(compiled.neck_count(),1);
    let (system,probe)=compiled.build(vec![body(800.0,0.0,0.0)],vec![],vec![],config(),&gate).unwrap();
    let mut system=system.prepare().unwrap();let energy=system.stored_energy_j();
    // Quadratic Gonzalez is implicit midpoint; use its known phase, not the
    // continuous phase, so this also catches a changed discrete time owner.
    let phase=2.0*(0.5*omega*config().dt_s).atan();
    for step in 1..=1000 {
        let frame=system.step(&[0.0;2],&gate).unwrap();
        let n=probe.neck_observation(system.state(),0).unwrap();
        let angle=step as f64*phase;
        assert!((n.displaced_volume_m3-initial*angle.cos()).abs()<1e-12);
        assert!((n.volume_flow_m3_s+omega*initial*angle.sin()).abs()<1e-9);
        assert_eq!(n.coordinate,1);assert_eq!(n.dissipated_power_w,0.0);
        assert_eq!(&system.state()[..2],&[0.0,0.0]);
        assert!(frame.balance_residual_j.abs()<1e-10);
    }
    assert!((system.stored_energy_j()-energy).abs()<1e-10);
}

#[test]
fn g1_neck_resistance_obeys_physical_flow_equation_and_discrete_power_balance() {
    let basis=cavity(vec![0.0]);let gate=CancelGate::new_clock_free();
    let mut opening=neck(vec![1.0],12000.0);opening.initial_volume_m3=0.0;opening.initial_flow_m3_s=0.001;
    let inertance=basis.rho0*opening.effective_length_m/opening.area_m2;
    let compliance=basis.lambdas[0]/(basis.rho0*basis.c0*basis.c0);
    let w2=1.0/(inertance*compliance);let drag=opening.resistance_pa_s_m3/inertance;
    let resistance=opening.resistance_pa_s_m3;
    let (mut volume,mut flow)=(opening.initial_volume_m3,opening.initial_flow_m3_s);
    let (system,probe)=CavityCoupling::new(&basis,1,&[0.0],&[0.0]).unwrap()
        .with_necks(vec![opening],&gate).unwrap()
        .build(vec![body(800.0,0.0,0.0)],vec![],vec![],config(),&gate).unwrap();
    let mut system=system.prepare().unwrap();let initial=system.stored_energy_j();let mut loss=0.0;
    for _ in 0..500 {
        let h=0.5*config().dt_s;
        let next_flow=(flow*(1.0-h*drag-h*h*w2)-2.0*h*w2*volume)/(1.0+h*drag+h*h*w2);
        let next_volume=volume+h*(flow+next_flow);
        let expected_loss=config().dt_s*resistance*(0.5*(flow+next_flow)).powi(2);
        let frame=system.step(&[0.0;2],&gate).unwrap();loss+=frame.dissipated_energy_j;
        let n=probe.neck_observation(system.state(),0).unwrap();
        assert!((n.displaced_volume_m3-next_volume).abs()<1e-12);
        assert!((n.volume_flow_m3_s-next_flow).abs()<1e-9);
        assert!((n.driving_pressure_pa+next_volume/compliance).abs()<1e-5);
        assert!((n.resistive_pressure_drop_pa-resistance*next_flow).abs()<1e-5);
        assert!((frame.dissipated_energy_j-expected_loss).abs()<1e-11);
        assert!(n.dissipated_power_w>=0.0);
        (volume,flow)=(next_volume,next_flow);
    }
    assert!(loss>1e-6);assert!(system.stored_energy_j()<initial);
    assert!((system.stored_energy_j()+loss-initial).abs()<1e-9);
}

#[test]
fn g3_multiple_necks_preserve_basis_scaling_and_empty_addition_preserves_sealed_path() {
    let gate=CancelGate::new_clock_free();let basis=cavity(vec![0.0,1000.0]);let mut scaled=basis.clone();
    for norm in &mut scaled.lambdas {*norm*=9.0;}
    for row in &mut scaled.interface {for p in row {*p*=3.0;}}
    let build=|b:&CavityModes,scale:f64| {
        let first=neck(vec![scale,0.2*scale],1000.0);let second=neck(vec![scale,-0.4*scale],2000.0);
        CavityCoupling::new(b,1,&[0.1*scale,0.07*scale],&[0.0,5.0]).unwrap()
            .with_necks(vec![first],&gate).unwrap().with_necks(vec![second],&gate).unwrap()
            .build(vec![body(800.0,1e-5,0.02)],vec![],vec![],config(),&gate).unwrap()
    };
    let (a,pa)=build(&basis,1.0);let (b,pb)=build(&scaled,3.0);
    let mut a=a.prepare().unwrap();let mut b=b.prepare().unwrap();
    assert_eq!(pa.total_modes(),4);assert_eq!(pa.neck_count(),2);
    for _ in 0..200 {
        a.step(&[0.0;4],&gate).unwrap();b.step(&[0.0;4],&gate).unwrap();
        for (&x,&y) in a.state().iter().zip(b.state()) {assert!((x-y).abs()<1e-9);}
        for index in 0..2 {
            let x=pa.neck_observation(a.state(),index).unwrap();let y=pb.neck_observation(b.state(),index).unwrap();
            assert_eq!(x.coordinate,2+index);
            assert!((x.driving_pressure_pa-y.driving_pressure_pa).abs()<1e-5);
            assert!((x.volume_flow_m3_s-y.volume_flow_m3_s).abs()<1e-9);
        }
    }
    let sealed=CavityCoupling::new(&basis,1,&[0.1,0.07],&[0.0,5.0]).unwrap();
    let make=|c:CavityCoupling| c.build(vec![body(800.0,1e-5,0.02)],vec![],vec![],config(),&gate)
        .unwrap().0.prepare().unwrap();
    let mut old=make(sealed.clone());let mut empty=make(sealed.with_necks(vec![],&gate).unwrap());
    for _ in 0..20 {
        old.step(&[0.0;2],&gate).unwrap();empty.step(&[0.0;2],&gate).unwrap();
        assert_eq!(old.state(),empty.state());
    }
}

#[test]
fn g0_g4_neck_admission_and_failed_steps_do_not_publish_motion() {
    let basis=cavity(vec![0.0,1000.0]);let gate=CancelGate::new_clock_free();
    let compiled=CavityCoupling::new(&basis,1,&[0.1,0.07],&[0.0,0.0]).unwrap();
    let valid=neck(vec![1.0,0.2],0.0);
    let mut bad=valid.clone();bad.resistance_pa_s_m3=-1.0;
    assert!(compiled.clone().with_necks(vec![bad],&gate).is_err());
    let mut bad=valid.clone();bad.area_m2=0.0;
    assert!(compiled.clone().with_necks(vec![bad],&gate).is_err());
    let mut bad=valid.clone();bad.effective_length_m=1.0;
    assert!(compiled.clone().with_necks(vec![bad],&gate).is_err());
    for averages in [vec![],vec![1.0],vec![1.0,f64::NAN],vec![0.0,0.0]] {
        let mut bad=valid.clone();bad.pressure_shape_averages=averages;
        assert!(compiled.clone().with_necks(vec![bad],&gate).is_err());
    }
    assert!(compiled.clone().with_necks(vec![valid.clone();9],&gate).is_err());
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(matches!(compiled.clone().with_necks(vec![valid.clone()],&cancelled),Err(ImpactError::Cancelled)));
    let compiled=compiled.with_necks(vec![valid],&gate).unwrap();
    let build=||compiled.clone().build(vec![body(800.0,1e-5,0.0)],vec![],vec![],config(),&gate).unwrap();
    let (system,probe)=build();let mut system=system.prepare().unwrap();let before=system.state().to_vec();
    assert!(matches!(system.step(&[0.0;3],&cancelled),Err(ImpactError::Cancelled)));
    assert_eq!(system.state(),before);assert_eq!(system.samples(),0);
    system.set_iteration_limit(0).unwrap();
    assert!(system.step(&[0.0;3],&gate).is_err());
    assert_eq!(system.state(),before);assert_eq!(system.samples(),0);
    system.set_iteration_limit(50).unwrap();
    assert!(probe.neck_observation(&[],0).is_err());
    assert!(probe.neck_observation(system.state(),1).is_err());
    let mut invalid=before.clone();invalid[5]=f64::INFINITY;
    assert!(probe.neck_observation(&invalid,0).is_err());
    let mut reference=build().0.prepare().unwrap();
    for _ in 0..20 {
        system.step(&[0.0;3],&gate).unwrap();reference.step(&[0.0;3],&gate).unwrap();
        assert_eq!(system.state(),reference.state());
    }
}
