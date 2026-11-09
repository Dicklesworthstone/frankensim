//! Independent original differential/algebraic equations for singular Kelvin
//! plane stress. No production projection, modes, or viscosity inverse is used.
use super::*;

fn startup()->[[C64;3];2] {
    STRAINS.map(|tone|tone.map(|v|C64::new(0.0,v.im)))
}

pub(super) fn dae_rate(ode:&OriginalOde,y:&[f64;20],plane:MaterialPlane,imposed:Voigt)->Voigt {
    let mut matrix=ode.d;let mut rhs=ode.spring(y).map(|v|-v);
    for axis in axes(plane).1 {
        if ode.d[axis].iter().all(|v|*v==0.0) {
            // A zero original dashpot row imposes zero spring traction.
            // Its derivative is C0*epsilon' = sum Cj*qj/tauj.
            matrix[axis]=core::array::from_fn(|j|
                ode.eq[axis][j]+ode.branches[0][axis][j]+ode.branches[1][axis][j]);
            rhs[axis]=(0..2).map(|arm|(0..6).map(|j|
                ode.branches[arm][axis][j]*y[6+6*arm+j]/TAUS[arm]).sum::<f64>()).sum();
        }
    }
    for axis in axes(plane).0 {
        matrix[axis]=[0.0;6];matrix[axis][axis]=1.0;rhs[axis]=imposed[axis];
    }
    solve(matrix,rhs)
}

fn integrate_dae(ode:&OriginalOde,mut y:[f64;20],plane:MaterialPlane,
    start:f64,duration:f64,frequencies:&[f64],strains:&[[C64;3]],bias:[f64;3])->[f64;20]
{
    y[18]=0.0;y[19]=0.0;let steps=4096;let h=duration/steps as f64;
    let derivative=|state:&[f64;20],t:f64| {
        let (value,rate)=waveform(frequencies,strains,bias,t);let mut state=*state;
        let e=mandel(prescribed(plane,value));for i in axes(plane).0 {state[i]=e[i];}
        let rate=dae_rate(ode,&state,plane,mandel(prescribed(plane,rate)));
        let mut result=[0.0;20];result[..6].copy_from_slice(&rate);
        for arm in 0..2 {for i in 0..6 {
            result[6+6*arm+i]=rate[i]-state[6+6*arm+i]/TAUS[arm];
        }}
        result[18]=dot(ode.stress(&state,rate),rate);
        result[19]=ode.heat_rate(&state,rate);result
    };
    for step in 0..steps {
        let t=start+step as f64*h;let k1=derivative(&y,t);
        let k2=derivative(&core::array::from_fn(|i|y[i]+0.5*h*k1[i]),t+0.5*h);
        let k3=derivative(&core::array::from_fn(|i|y[i]+0.5*h*k2[i]),t+0.5*h);
        let k4=derivative(&core::array::from_fn(|i|y[i]+h*k3[i]),t+h);
        for i in 0..20 {y[i]+=h*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])/6.0;}
    }
    let end=mandel(prescribed(plane,waveform(frequencies,strains,bias,start+duration).0));
    for i in axes(plane).0 {y[i]=end[i];}y
}

fn compare_account(point:&fs_material::visco_tensor::kelvin_control::TensorKelvinPlaneStressHarmonicPoint,
    ode:&OriginalOde,reference:&[f64;20],plane:MaterialPlane,rate:[f64;3])
{
    let rate=dae_rate(ode,reference,plane,mandel(prescribed(plane,rate)));
    ode.compare(&point.history,reference,rate);
    for i in 0..6 {near(point.strain_rate[i],rate[i]/W[i],2e-10);}
    near(point.account.stress_work_j_m3,reference[18],2e-13);
    near(point.account.material_heat_j_m3,reference[19],2e-13);
    assert!(point.account.material_heat_j_m3>0.0);
    near(point.account.energy_balance_residual_j_m3,0.0,2e-13);
    for axis in axes(plane).1 {assert_eq!(point.history.stress()[axis],0.0);}
}

#[test]
fn directional_shear_viscosity_rank_two_matches_original_dae_in_every_plane() {
    let (law,ode)=fixture(diagonal([0.0,0.0,0.0,0.6,0.9,1.3]));
    let tones=startup();let duration=0.29;
    for plane in [MaterialPlane::Axes12,MaterialPlane::Axes21,MaterialPlane::Axes13,MaterialPlane::Axes23] {
        let history=law.history();let original=history.clone();
        let control=Control::new(&history,plane,TensorKelvinStressControlBudget::default(),||false).unwrap();
        assert_eq!(control.transverse_viscosity_rank(),2);
        assert_eq!(control.state_coordinates(),15);
        let point=control.prepare_harmonic_window(&FREQUENCIES,Dd::ZERO,duration,budget(),||false).unwrap()
            .evaluate(&tones,&[0.0;3],&history,||false).unwrap();
        let reference=integrate_dae(&ode,[0.0;20],plane,0.0,duration,&FREQUENCIES,&tones,[0.0;3]);
        compare_account(&point,&ode,&reference,plane,waveform(&FREQUENCIES,&tones,[0.0;3],duration).1);
        near(point.account.stored_energy_change_j_m3,ode.stored(&reference),2e-13);
        assert_eq!(history,original);
    }
}

#[test]
fn ranks_zero_and_one_preserve_memory_across_subdivision_and_a_changed_bank() {
    let plane=MaterialPlane::Axes12;let tones=startup();let duration=0.375;
    for (rank,diagonal_values) in [
        (0,[0.7,0.4,0.0,0.6,0.0,0.0]),
        (1,[0.7,0.4,1.1,0.6,0.0,0.0]),
    ] {
        let (law,ode)=fixture(diagonal(diagonal_values));let history=law.history();
        let control=Control::new(&history,plane,TensorKelvinStressControlBudget::default(),||false).unwrap();
        assert_eq!(control.transverse_viscosity_rank(),rank);
        assert_eq!(control.state_coordinates(),15);
        // The original no-Kelvin mixed controller must retain its refusal.
        assert!(PreparedTensorMixedControl::new(&history,[false,false,true,false,true,true],
            TensorStressControlBudget::default(),||false).is_err());
        let whole=control.prepare_harmonic_window(&FREQUENCIES,Dd::ZERO,duration,budget(),||false).unwrap()
            .evaluate(&tones,&[0.0;3],&history,||false).unwrap();
        let left=control.prepare_harmonic_window(&FREQUENCIES,Dd::ZERO,0.125,budget(),||false).unwrap()
            .evaluate(&tones,&[0.0;3],&history,||false).unwrap();
        let right=control.prepare_harmonic_window(&FREQUENCIES,Dd::from_f64(0.125),0.25,budget(),||false).unwrap()
            .evaluate(&tones,&[0.0;3],&left.history,||false).unwrap();
        let reference=integrate_dae(&ode,[0.0;20],plane,0.0,duration,&FREQUENCIES,&tones,[0.0;3]);
        compare_account(&whole,&ode,&reference,plane,waveform(&FREQUENCIES,&tones,[0.0;3],duration).1);
        for i in 0..6 {
            near(right.history.strain()[i],whole.history.strain()[i],3e-14);
            near(right.history.stress()[i],whole.history.stress()[i],3e-12);
        }
        near(whole.account.material_heat_j_m3,left.account.material_heat_j_m3+right.account.material_heat_j_m3,3e-13);
        near(whole.account.stress_work_j_m3,left.account.stress_work_j_m3+right.account.stress_work_j_m3,3e-13);
        let strain=right.history.strain();let held=[strain[0],strain[1],2.0*strain[3]];
        let hold=control.prepare_harmonic_window(&[],Dd::from_f64(duration),0.2,budget(),||false).unwrap()
            .evaluate(&[],&held,&right.history,||false).unwrap();
        let held_reference=integrate_dae(&ode,reference,plane,duration,0.2,&[],&[],held);
        compare_account(&hold,&ode,&held_reference,plane,[0.0;3]);
        assert_eq!(hold.account.stress_work_j_m3,0.0);
        near(hold.account.stored_energy_change_j_m3,-hold.account.material_heat_j_m3,2e-13);
        assert!(hold.history.branch_stored_energy()>0.0);
    }
}

#[test]
fn zero_dynamic_modes_keep_original_elastic_closure_and_direct_kelvin_heat() {
    let mut eq=diagonal([2.0,3.0,4.0,5.0,6.0,7.0]);eq[0][2]=1.0;eq[2][0]=1.0;
    let eta=0.7;
    let law=TensorProny::new_mandel(eq,vec![],0.1).unwrap()
        .with_kelvin_viscosity_mandel(diagonal([eta,0.0,0.0,0.0,0.0,0.0])).unwrap();
    let history=law.history();let original=history.clone();
    let control=Control::new(&history,MaterialPlane::Axes12,TensorKelvinStressControlBudget::default(),||false).unwrap();
    assert_eq!(control.transverse_viscosity_rank(),0);assert_eq!(control.state_coordinates(),3);
    let (amplitude,omega,duration)=(0.001_f64,3.0_f64,0.2_f64);
    let point=control.prepare_harmonic_window(&[omega],Dd::ZERO,duration,budget(),||false).unwrap()
        .evaluate(&[[C64::new(0.0,amplitude),C64::ZERO,C64::ZERO]],&[0.0;3],&history,||false).unwrap();
    let value=amplitude*(omega*duration).sin();let rate=amplitude*omega*(omega*duration).cos();
    let elastic=0.5*(2.0-1.0/4.0)*value*value;
    let heat=eta*amplitude*amplitude*omega*omega*(duration/2.0+(2.0*omega*duration).sin()/(4.0*omega));
    near(point.history.strain()[2],-value/4.0,2e-15);
    near(point.strain_rate[2],-rate/4.0,2e-15);
    near(point.account.end_stored_energy_j_m3,elastic,2e-18);
    near(point.account.material_heat_j_m3,heat,2e-18);
    near(point.account.stress_work_j_m3,elastic+heat,2e-18);
    assert_eq!(point.history.branch_stored_energy(),0.0);
    for axis in [2,4,5] {assert_eq!(point.history.stress()[axis],0.0);}
    assert_eq!(history,original);
}

#[test]
fn incompatible_or_unresolved_algebraic_states_refuse_without_projecting_history() {
    let plane=MaterialPlane::Axes12;
    let (law,_)=fixture(diagonal([0.0,0.0,0.0,0.6,0.9,1.3]));
    let (history,_)=seed(&law,plane,BIAS);let original=history.clone();
    let control=Control::new(&history,plane,TensorKelvinStressControlBudget::default(),||false).unwrap();
    let window=control.prepare_harmonic_window(&[],Dd::ZERO,0.2,budget(),||false).unwrap();
    // A full-strain ramp may retain nonzero algebraic spring traction. The
    // causal plane-stress path must refuse instead of making an initial jump.
    assert!(window.evaluate(&[],&BIAS,&history,||false).is_err());
    assert_eq!(window.evaluate(&[],&BIAS,&history,||true).unwrap_err(),TensorPronyError::Cancelled);
    assert_eq!(history,original);

    // Two independent original full-space loss channels have collinear
    // transverse projections. This bounded lane refuses that ambiguous rank;
    // it cannot discard a small singular value or add a damping floor.
    // Original passive factoring pivots yy (16), then xx (9), leaving
    // two collinear nonzero zz Gram projections (2 and 1).
    let a=[3.0,0.0,1.0,0.0,0.0,0.0];let b=[0.0,4.0,2.0,0.0,0.0,0.0];
    let redundant=core::array::from_fn(|i|core::array::from_fn(|j|a[i]*a[j]+b[i]*b[j]));
    let (law,_)=fixture(redundant);let history=law.history();let original=history.clone();
    assert!(Control::new(&history,plane,TensorKelvinStressControlBudget::default(),||false).is_err());
    assert_eq!(history,original);

    let (law,_)=fixture(diagonal([0.7,0.4,0.0,0.6,0.0,0.0]));let history=law.history();
    let control=Control::new(&history,plane,TensorKelvinStressControlBudget::default(),||false).unwrap();
    let huge=[C64::new(1e100,-2e100);3];let opposite=huge.map(|v|-v);
    let point=control.prepare_harmonic_window(&[3.0,3.0],Dd::ZERO,0.2,budget(),||false).unwrap()
        .evaluate(&[huge,opposite],&[0.0;3],&history,||false).unwrap();
    assert_eq!(point.history.strain(),history.strain());assert_eq!(point.account.material_heat_j_m3,0.0);
    assert_eq!(point.account.stress_work_j_m3,0.0);
}
