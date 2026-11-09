//! Original Kelvin/Maxwell equations, independent of the production modes.
use fs_material::plane_stress_viscosity::MaterialPlane;
use fs_material::visco_tensor::kelvin_control::{PreparedTensorKelvinPlaneStress as Control,
    TensorKelvinStressControlBudget};
use fs_material::visco_tensor::stress_control::{PreparedTensorMixedControl,
    TensorPlaneStressHarmonicBudget,TensorStressControlBudget};
use fs_material::visco_tensor::{TensorProny,TensorPronyBranch,TensorPronyHistory,TensorPronyError};
use fs_material::{Tangent6,Voigt};
use fs_math::{c64::C64,dd::Dd};

const TAUS:[f64;2]=[0.23,1.1];
const FREQUENCIES:[f64;2]=[2.3,3.9];
const BIAS:[f64;3]=[0.001,-0.0005,0.0007];
const STRAINS:[[C64;3];2]=[
    [C64::new(0.0012,-0.0007),C64::new(-0.0005,0.0009),C64::new(0.0006,0.0011)],
    [C64::new(-0.0008,0.0004),C64::new(0.0007,-0.0003),C64::new(-0.0009,0.0002)],
];
const W:Voigt=[1.0,1.0,1.0,std::f64::consts::SQRT_2,std::f64::consts::SQRT_2,std::f64::consts::SQRT_2];
fn budget()->TensorPlaneStressHarmonicBudget {
    TensorPlaneStressHarmonicBudget {maximum_payload_bytes:4_000_000,maximum_scalar_work:5_000_000}
}
fn axes(plane:MaterialPlane)->([usize;3],[usize;3]) {
    match plane {
        MaterialPlane::Axes12=>([0,1,3],[2,4,5]),MaterialPlane::Axes21=>([1,0,3],[2,4,5]),
        MaterialPlane::Axes13=>([0,2,5],[1,3,4]),MaterialPlane::Axes23=>([1,2,4],[0,3,5]),
    }
}
fn prescribed(plane:MaterialPlane,e:[f64;3])->Voigt {
    let mut v=[0.0;6];let k=axes(plane).0;
    v[k[0]]=e[0];v[k[1]]=e[1];v[k[2]]=e[2]/2.0;v
}
fn mandel(v:Voigt)->Voigt {core::array::from_fn(|i|v[i]*W[i])}
fn diagonal(v:Voigt)->Tangent6 {core::array::from_fn(|i|core::array::from_fn(|j|if i==j {v[i]}else{0.0}))}
fn apply(a:&Tangent6,x:Voigt)->Voigt {core::array::from_fn(|i|(0..6).map(|j|a[i][j]*x[j]).sum())}
fn dot(a:Voigt,b:Voigt)->f64 {(0..6).map(|i|a[i]*b[i]).sum()}
fn near(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<=tol,"actual={a:e}, expected={b:e}, tolerance={tol:e}");}
fn viscosity()->Tangent6 {
    let a=[0.5,-0.4,0.3,0.2,-0.1,0.6];let d=[0.7,0.4,1.1,0.6,0.9,1.3];
    core::array::from_fn(|i|core::array::from_fn(|j|a[i]*a[j]+if i==j {d[i]}else{0.0}))
}
fn fixture(d:Tangent6)->(TensorProny,OriginalOde) {
    let mut eq=diagonal([12.0,17.0,10.0,9.0,7.0,13.0]);
    eq[0][2]=-1.3;eq[2][0]=-1.3;eq[1][4]=0.9;eq[4][1]=0.9;
    let mut factors=[[[0.0;6];6];2];
    factors[0][0]=[1.7,-0.3,1.4,0.7,0.1,-0.2];
    factors[1][0]=[0.1,1.2,-0.2,0.6,0.9,-0.4];
    factors[1][1]=[-0.3,0.2,0.8,-0.4,0.0,1.1];
    let branches=factors.map(|b|core::array::from_fn(|i|core::array::from_fn(|j|(0..6).map(|k|b[k][i]*b[k][j]).sum())));
    let arms=factors.into_iter().zip(TAUS).map(|(b,tau)|TensorPronyBranch::from_gram_rows(b,tau).unwrap()).collect();
    let law=TensorProny::new_mandel(eq,arms,0.1).unwrap().with_kelvin_viscosity_mandel(d).unwrap().with_band((0.0,100.0)).unwrap();
    (law,OriginalOde {eq,d,branches})
}
fn waveform(frequencies:&[f64],strains:&[[C64;3]],bias:[f64;3],t:f64)->([f64;3],[f64;3]) {
    let mut e=bias;let mut rate=[0.0;3];
    for (&w,strain) in frequencies.iter().zip(strains) {
        let (s,c)=(w*t).sin_cos();for i in 0..3 {
            e[i]+=strain[i].re*c+strain[i].im*s;rate[i]+=w*(-strain[i].re*s+strain[i].im*c);
        }
    }
    (e,rate)
}
fn seed(law:&TensorProny,plane:MaterialPlane,initial:[f64;3])->(TensorPronyHistory,[f64;20]) {
    let mut e=[0.002,-0.001,0.0015,0.0007,-0.0004,0.0005];
    let imposed=prescribed(plane,initial);for i in axes(plane).0 {e[i]=imposed[i];}
    let mut history=law.history();history.advance(e,0.17,||false).unwrap();
    let e=mandel(history.strain());let mut state=[0.0;20];state[..6].copy_from_slice(&e);
    for (j,tau) in TAUS.into_iter().enumerate() {
        let r=0.17/tau;let phi=-(-r).exp_m1()/r;
        for i in 0..6 {state[6+6*j+i]=phi*e[i];}
    }
    (history,state)
}
fn solve(mut a:Tangent6,mut rhs:Voigt)->Voigt {
    for j in 0..6 {
        let pivot=(j..6).max_by(|&i,&k|a[i][j].abs().total_cmp(&a[k][j].abs())).unwrap();
        a.swap(j,pivot);rhs.swap(j,pivot);let row=a[j];
        for i in j+1..6 {let f=a[i][j]/row[j];for k in j..6 {a[i][k]-=f*row[k];}rhs[i]-=f*rhs[j];}
    }
    let mut x=[0.0;6];for i in (0..6).rev() {x[i]=(rhs[i]-(i+1..6).map(|j|a[i][j]*x[j]).sum::<f64>())/a[i][i];}x
}

// State: full Mandel strain, two raw original Maxwell spring strains, work,
// heat. Physical transverse rate is solved from D_FF, never a modal transform.
struct OriginalOde {eq:Tangent6,d:Tangent6,branches:[Tangent6;2]}
impl OriginalOde {
    fn spring(&self,y:&[f64;20])->Voigt {
        let mut s=apply(&self.eq,core::array::from_fn(|i|y[i]));
        for j in 0..2 {let arm=apply(&self.branches[j],core::array::from_fn(|i|y[6+6*j+i]));for i in 0..6 {s[i]+=arm[i];}}
        s
    }
    fn rate(&self,y:&[f64;20],plane:MaterialPlane,imposed:Voigt)->Voigt {
        let mut a=self.d;let mut rhs=self.spring(y).map(|v|-v);
        for i in axes(plane).0 {a[i]=[0.0;6];a[i][i]=1.0;rhs[i]=imposed[i];}
        solve(a,rhs)
    }
    fn stress(&self,y:&[f64;20],rate:Voigt)->Voigt {
        let s=self.spring(y);let viscous=apply(&self.d,rate);core::array::from_fn(|i|s[i]+viscous[i])
    }
    fn branch_energy(&self,y:&[f64;20])->f64 {(0..2).map(|j| {
        let q=core::array::from_fn(|i|y[6+6*j+i]);0.5*dot(q,apply(&self.branches[j],q))
    }).sum()}
    fn stored(&self,y:&[f64;20])->f64 {let e=core::array::from_fn(|i|y[i]);0.5*dot(e,apply(&self.eq,e))+self.branch_energy(y)}
    fn heat_rate(&self,y:&[f64;20],rate:Voigt)->f64 {
        dot(rate,apply(&self.d,rate))+(0..2).map(|j| {
            let q=core::array::from_fn(|i|y[6+6*j+i]);dot(q,apply(&self.branches[j],q))/TAUS[j]
        }).sum::<f64>()
    }
    fn integrate(&self,mut y:[f64;20],plane:MaterialPlane,start:f64,duration:f64,
        frequencies:&[f64],strains:&[[C64;3]],bias:[f64;3])->[f64;20]
    {
        y[18]=0.0;y[19]=0.0;let steps=4096;let h=duration/steps as f64;
        let derivative=|state:&[f64;20],t:f64| {
            let (value,rate)=waveform(frequencies,strains,bias,t);let mut state=*state;
            let e=mandel(prescribed(plane,value));for i in axes(plane).0 {state[i]=e[i];}
            let rate=self.rate(&state,plane,mandel(prescribed(plane,rate)));
            let mut result=[0.0;20];result[..6].copy_from_slice(&rate);
            for j in 0..2 {for i in 0..6 {result[6+6*j+i]=rate[i]-state[6+6*j+i]/TAUS[j];}}
            result[18]=dot(self.stress(&state,rate),rate);result[19]=self.heat_rate(&state,rate);result
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
    fn compare(&self,h:&TensorPronyHistory,y:&[f64;20],rate:Voigt) {
        let stress=self.stress(y,rate);
        for i in 0..6 {near(h.strain()[i],y[i]/W[i],3e-13);near(h.stress()[i],stress[i]/W[i],2e-11);}
        near(h.stored_energy(),self.stored(y),2e-13);near(h.branch_stored_energy(),self.branch_energy(y),2e-13);
    }
}

#[test]
fn full_anisotropic_kelvin_maxwell_matches_original_ode_in_all_ordered_planes() {
    let (law,ode)=fixture(viscosity());
    for plane in [MaterialPlane::Axes12,MaterialPlane::Axes21,MaterialPlane::Axes13,MaterialPlane::Axes23] {
        let (t,h)=(0.37,0.29);let (history,y)=seed(&law,plane,waveform(&FREQUENCIES,&STRAINS,BIAS,t).0);
        let original=history.clone();
        let control=Control::new(&history,plane,TensorKelvinStressControlBudget::default(),||false).unwrap();
        assert_eq!(control.state_coordinates(),15);
        let point=control.prepare_harmonic_window(&FREQUENCIES,Dd::from_f64(t),h,budget(),||false).unwrap()
            .evaluate(&STRAINS,&BIAS,&history,||false).unwrap();
        let reference=ode.integrate(y,plane,t,h,&FREQUENCIES,&STRAINS,BIAS);
        let rate=ode.rate(&reference,plane,mandel(prescribed(plane,waveform(&FREQUENCIES,&STRAINS,BIAS,t+h).1)));
        ode.compare(&point.history,&reference,rate);
        for i in 0..6 {near(point.strain_rate[i],rate[i]/W[i],2e-10);}
        near(point.account.stress_work_j_m3,reference[18],2e-13);
        near(point.account.material_heat_j_m3,reference[19],2e-13);
        near(point.account.stored_energy_change_j_m3,ode.stored(&reference)-ode.stored(&y),2e-13);
        assert!(point.account.material_heat_j_m3>0.0);near(point.account.energy_balance_residual_j_m3,0.0,2e-13);
        assert_eq!(history,original);
        // Handoff into the unchanged original full-strain ramp. Its endpoint
        // spring strains have an independent exact original-branch expression.
        let old=mandel(point.history.strain());let target=old.map(|e|0.8*e);let dt=0.13;
        let mut resumed=point.history.clone();resumed.advance(core::array::from_fn(|i|target[i]/W[i]),dt,||false).unwrap();
        let mut expected=reference;expected[..6].copy_from_slice(&target);
        for (j,tau) in TAUS.into_iter().enumerate() {
            let r=dt/tau;let decay=(-r).exp();let phi=-(-r).exp_m1()/r;
            for i in 0..6 {expected[6+6*j+i]=decay*reference[6+6*j+i]+phi*(target[i]-old[i]);}
        }
        ode.compare(&resumed,&expected,core::array::from_fn(|i|(target[i]-old[i])/dt));
    }
}

#[test]
fn semidefinite_full_viscosity_with_spd_transverse_block_retains_physical_loss() {
    let (law,ode)=fixture(diagonal([0.0,0.0,1.1,0.0,0.9,1.3]));let plane=MaterialPlane::Axes12;
    let (history,y)=seed(&law,plane,BIAS);
    assert!(PreparedTensorMixedControl::new(&history,[false,false,true,false,true,true],TensorStressControlBudget::default(),||false).is_err());
    let control=Control::new(&history,plane,TensorKelvinStressControlBudget::default(),||false).unwrap();
    let h=0.31;let point=control.prepare_harmonic_window(&[],Dd::ZERO,h,budget(),||false).unwrap().evaluate(&[],&BIAS,&history,||false).unwrap();
    let reference=ode.integrate(y,plane,0.0,h,&[],&[],BIAS);let rate=ode.rate(&reference,plane,[0.0;6]);
    ode.compare(&point.history,&reference,rate);near(point.account.material_heat_j_m3,reference[19],2e-13);
    assert!(point.account.material_heat_j_m3>0.0);assert_eq!(point.account.stress_work_j_m3,0.0);
    near(point.account.stored_energy_change_j_m3,-point.account.material_heat_j_m3,2e-13);
}

#[test]
fn tiny_retained_interval_keeps_positive_heat_and_direct_storage_loss() {
    let (law,ode)=fixture(viscosity());let plane=MaterialPlane::Axes12;let bias=[0.001,-0.0004,0.0];
    let (history,y)=seed(&law,plane,bias);let rate=ode.rate(&y,plane,[0.0;6]);
    let h=1e-20;let expected=ode.heat_rate(&y,rate)*h;
    let control=Control::new(&history,plane,TensorKelvinStressControlBudget::default(),||false).unwrap();
    let point=control.prepare_harmonic_window(&[],Dd::from_pair(1e16,0.3),h,budget(),||false).unwrap()
        .evaluate(&[],&bias,&history,||false).unwrap();
    assert!(expected>0.0);near(point.account.material_heat_j_m3,expected,5e-11*expected);
    near(point.account.stored_energy_change_j_m3,-expected,5e-11*expected);
    assert_eq!(point.account.stress_work_j_m3,0.0);
    assert_eq!(point.account.start_stored_energy_j_m3,point.account.end_stored_energy_j_m3);
    assert!(point.strain_increment.iter().any(|v|*v!=0.0));
}

#[test]
fn zero_maxwell_startup_recovers_cubic_kelvin_heat_without_a_loss_floor() {
    let mut eq=diagonal([2.0,3.0,4.0,5.0,6.0,7.0]);eq[0][2]=1.0;eq[2][0]=1.0;
    let law=TensorProny::new_mandel(eq,vec![],0.1).unwrap()
        .with_kelvin_viscosity_mandel(diagonal([0.0,0.0,2.0,0.0,1.0,1.0])).unwrap();
    let history=law.history();let control=Control::new(&history,MaterialPlane::Axes12,TensorKelvinStressControlBudget::default(),||false).unwrap();
    let (a,w,h)=(0.001_f64,3.0_f64,1e-10_f64);
    let point=control.prepare_harmonic_window(&[w],Dd::ZERO,h,budget(),||false).unwrap()
        .evaluate(&[[C64::new(0.0,a),C64::ZERO,C64::ZERO]],&[0.0;3],&history,||false).unwrap();
    // eta*ezz'=-Cxz*ex-Czz*ezz; ezz'=-Cxz*a*w*t/eta+O(t²).
    let expected=a*a*w*w*h*h*h/6.0;
    assert!(point.account.material_heat_j_m3>0.0);
    near(point.account.material_heat_j_m3,expected,2e-9*expected);
    assert_eq!(point.history.branch_stored_energy(),0.0);
    for i in [2,4,5] {assert_eq!(point.history.stress()[i],0.0);}
}

#[test]
fn dd_subdivision_and_coherent_cancellation_keep_the_original_kelvin_history() {
    let (law,_)=fixture(viscosity());let plane=MaterialPlane::Axes13;
    let start=Dd::from_pair(1e16,0.3);let mut initial=BIAS.map(Dd::from_f64);
    for k in 0..2 {let p=fs_fft::harmonic_window::phase_angular_interval_from_clock_dd(FREQUENCIES[k],start,0.125).unwrap();
        for j in 0..3 {initial[j]=initial[j]+Dd::from_f64(STRAINS[k][j].re)*p.start[0]-Dd::from_f64(STRAINS[k][j].im)*p.start[1];}}
    let (history,_)=seed(&law,plane,initial.map(Dd::to_f64));
    let control=Control::new(&history,plane,TensorKelvinStressControlBudget::default(),||false).unwrap();
    let whole=control.prepare_harmonic_window(&FREQUENCIES,start,0.375,budget(),||false).unwrap().evaluate(&STRAINS,&BIAS,&history,||false).unwrap();
    let left=control.prepare_harmonic_window(&FREQUENCIES,start,0.125,budget(),||false).unwrap().evaluate(&STRAINS,&BIAS,&history,||false).unwrap();
    let right=control.prepare_harmonic_window(&FREQUENCIES,start+Dd::from_f64(0.125),0.25,budget(),||false).unwrap().evaluate(&STRAINS,&BIAS,&left.history,||false).unwrap();
    for i in 0..6 {near(whole.history.strain()[i],right.history.strain()[i],3e-14);near(whole.history.stress()[i],right.history.stress()[i],3e-12);}
    near(whole.account.material_heat_j_m3,left.account.material_heat_j_m3+right.account.material_heat_j_m3,3e-13);
    near(whole.account.stress_work_j_m3,left.account.stress_work_j_m3+right.account.stress_work_j_m3,3e-13);
    let held=[right.history.strain()[0],right.history.strain()[2],2.0*right.history.strain()[5]];
    let huge=[C64::new(1e100,-2e100);3];let opposite=huge.map(|v|-v);
    let time=start+Dd::from_f64(0.375);
    let cancelled=control.prepare_harmonic_window(&[3.0,3.0],time,0.2,budget(),||false).unwrap().evaluate(&[huge,opposite],&held,&right.history,||false).unwrap();
    let empty=control.prepare_harmonic_window(&[],time,0.2,budget(),||false).unwrap().evaluate(&[],&held,&right.history,||false).unwrap();
    assert_eq!(cancelled.history,empty.history);assert_eq!(cancelled.account.material_heat_j_m3,empty.account.material_heat_j_m3);
    assert!(empty.account.material_heat_j_m3>0.0);assert_eq!(empty.account.stress_work_j_m3,0.0);
}

#[test]
fn source_domain_positive_transverse_admission_and_cancellation_are_atomic() {
    let (law,_)=fixture(viscosity());let plane=MaterialPlane::Axes12;let (history,_)=seed(&law,plane,BIAS);
    let original=history.clone();let control=Control::new(&history,plane,TensorKelvinStressControlBudget::default(),||false).unwrap();
    let window=control.prepare_harmonic_window(&[],Dd::ZERO,0.2,budget(),||false).unwrap();
    assert!(window.evaluate(&[],&BIAS,&law.history(),||false).is_err());
    assert!(window.evaluate(&[],&[0.02,0.0,0.0],&history,||false).is_err());
    assert!(control.prepare_harmonic_window(&[101.0],Dd::ZERO,0.2,budget(),||false).is_err());
    assert!(control.prepare_harmonic_window(&[],Dd::ZERO,0.2,TensorPlaneStressHarmonicBudget{maximum_payload_bytes:0,..budget()},||false).is_err());
    let mut calls=0;assert_eq!(window.evaluate(&[],&BIAS,&history,||{calls+=1;calls>3}).unwrap_err(),TensorPronyError::Cancelled);
    assert!(matches!(Control::new(&history,plane,TensorKelvinStressControlBudget::default(),||true),Err(TensorPronyError::Cancelled)));
    assert_eq!(history,original);
    let (singular,_)=fixture(diagonal([1.0,1.0,0.0,1.0,1.0,1.0]));
    let resolved=Control::new(&singular.history(),plane,TensorKelvinStressControlBudget::default(),||false).unwrap();
    assert_eq!(resolved.transverse_viscosity_rank(),2);
    let restricted=law.with_band((1.0,10.0)).unwrap();let (h,_)=seed(&restricted,plane,BIAS);
    let c=Control::new(&h,plane,TensorKelvinStressControlBudget::default(),||false).unwrap();
    assert!(c.prepare_harmonic_window(&[],Dd::ZERO,0.2,budget(),||false).unwrap().evaluate(&[],&BIAS,&h,||false).is_err());
}

#[path = "visco_kelvin_plane_stress_harmonic/semidefinite.rs"]
mod semidefinite;

#[path = "visco_kelvin_plane_stress_harmonic/ramp.rs"]
mod ramp;
