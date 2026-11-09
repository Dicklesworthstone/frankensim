//! G1/G3/G4: original nonlinear device Q(V), exact physical stored energy,
//! work-conjugate electrode-charge chords, and finite commutation dynamics.
//! The synthetic source has Cdiff=3c below |V|=1 and c above it. Its energy
//! above the knee is c*(1+V^2/2), rather than Cdiff(V)*V^2/2.

use fs_phs::circuit::capacitance::CapacitorChargeCurve;
use fs_phs::circuit::converter::commutation::capacitance::charge::{
    ChargeCommutationParameters, ChargeCommutationStorage,
};
use fs_phs::circuit::converter::commutation::{
    BodyDiode, CommutationParameters, DeviceOperatingBounds, SwitchDevice,
};
use fs_phs::{
    ConstitutiveDiscreteGradient, ContinuousStorageTangent, Storage,
};

const DEVICE: SwitchDevice = SwitchDevice {
    on_resistance_ohm: 0.2,
    off_resistance_ohm: 40.0,
    body_diode: BodyDiode { forward_voltage_v: 0.5, differential_resistance_ohm: 0.1 },
    bounds: DeviceOperatingBounds {
        maximum_reverse_voltage_v: 100.0,
        maximum_channel_current_a: 100.0,
        maximum_diode_current_a: 100.0,
    },
};
fn close(actual:f64, expected:f64, tolerance:f64) {
    assert!((actual-expected).abs()<=tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance {tolerance:e}");
}
fn curve(c:f64)->CapacitorChargeCurve {
    CapacitorChargeCurve::try_new(vec![
        (-22.0*c,-20.0),(-3.0*c,-1.0),(0.0,0.0),(3.0*c,1.0),(22.0*c,20.0),
    ]).unwrap()
}
fn parameters(bus:[f64;2],cap:[f64;4])->ChargeCommutationParameters {
    ChargeCommutationParameters {
        converter:CommutationParameters {
            inductance_h:0.4,bus_capacitance_f:bus,winding_resistance_ohm:0.1,
            devices:[DEVICE;4],
        },
        output_charge_curves:cap.map(curve),
        maximum_cells_per_chord:64,
    }
}
// Independent source integrals in PHYSICAL voltage coordinates.
fn charge(c:f64,v:f64)->f64 {
    if v.abs()<=1.0 {3.0*c*v} else {c*(v+2.0*v.signum())}
}
fn energy(c:f64,v:f64)->f64 {
    if v.abs()<=1.0 {1.5*c*v*v} else {c*(1.0+0.5*v*v)}
}
fn expected_state(p:&ChargeCommutationParameters,c:[f64;4],i:f64,b:[f64;2],n:[f64;2])
    ->[f64;5]
{
    let h=[charge(c[0],b[0]-n[0]),charge(c[2],b[1]-n[1])];
    [p.converter.inductance_h*i,
     p.converter.bus_capacitance_f[0]*b[0]+h[0],
     p.converter.bus_capacitance_f[1]*b[1]+h[1],
     charge(c[1],n[0])-h[0],charge(c[3],n[1])-h[1]]
}
#[derive(Clone,Copy)]
struct Segment {low:f64,high:f64,c:f64,offset:f64}
fn segments(c:f64)->[Segment;3] {
    [
        Segment{low:-20.0,high:-1.0,c,offset:-2.0*c},
        Segment{low:-1.0,high:1.0,c:3.0*c,offset:0.0},
        Segment{low:1.0,high:20.0,c,offset:2.0*c},
    ]
}
fn inverse(cb:f64,h:Segment,l:Segment)->[f64;3] {
    let det=cb*(h.c+l.c)+h.c*l.c;
    [(h.c+l.c)/det,h.c/det,(cb+h.c)/det]
}
fn in_cell(q:[f64;2],h:Segment,l:Segment,inverse:[f64;3])->[f64;2] {
    let b=q[0]-h.offset;
    let n=q[1]-l.offset+h.offset;
    [inverse[0]*b+inverse[1]*n,inverse[1]*b+inverse[2]*n]
}
fn leg_voltage(cb:f64,ch:f64,cl:f64,q:[f64;2])->[f64;2] {
    // Enumerate the nine analytical voltage cells of this FIXED fixture.
    // Production instead traces original source knots in device charge.
    for high in segments(ch) {for low in segments(cl) {
        let v=in_cell(q,high,low,inverse(cb,high,low));
        let h=v[0]-v[1];
        if h>=high.low-1e-12&&h<=high.high+1e-12
            &&v[1]>=low.low-1e-12&&v[1]<=low.high+1e-12 {return v;}
    }}
    panic!("fixture electrode charge is outside its explicit voltage support");
}
fn physical(p:&ChargeCommutationParameters,c:[f64;4],x:&[f64;5])->[f64;5] {
    let mut out=[0.0;5];out[0]=x[0]/p.converter.inductance_h;
    for k in 0..2 {
        let v=leg_voltage(p.converter.bus_capacitance_f[k],c[2*k],c[2*k+1],[x[1+k],x[3+k]]);
        out[1+k]=v[0];out[3+k]=v[1];
    }
    out
}
fn parts(p:&ChargeCommutationParameters,c:[f64;4],x:&[f64;5])->[f64;7] {
    let v=physical(p,c,x);
    [0.5*p.converter.inductance_h*v[0]*v[0],
     0.5*p.converter.bus_capacitance_f[0]*v[1]*v[1],
     0.5*p.converter.bus_capacitance_f[1]*v[2]*v[2],
     energy(c[0],v[1]-v[3]),energy(c[1],v[3]),
     energy(c[2],v[2]-v[4]),energy(c[3],v[4])]
}
fn clip_interval(value:f64,derivative:f64,low:f64,high:f64,t:&mut[f64;2])->bool {
    if derivative==0.0 {return value>=low&&value<=high;}
    let a=(low-value)/derivative;let b=(high-value)/derivative;
    t[0]=t[0].max(a.min(b));t[1]=t[1].min(a.max(b));t[1]>t[0]
}
fn analytic_chord(p:&ChargeCommutationParameters,c:[f64;4],a:&[f64;5],b:&[f64;5],d:&[f64;5])
    ->([f64;5],[f64;5])
{
    let mut mean=[0.0;5];let mut jac=[0.0;5];
    mean[0]=f64::midpoint(a[0],b[0])/p.converter.inductance_h;
    jac[0]=0.5*d[0]/p.converter.inductance_h;
    for k in 0..2 {
        let mut covered=0.0;
        for high in segments(c[2*k]) {for low in segments(c[2*k+1]) {
            let h=inverse(p.converter.bus_capacitance_f[k],high,low);
            let v0=in_cell([a[1+k],a[3+k]],high,low,h);
            let dq=[b[1+k]-a[1+k],b[3+k]-a[3+k]];
            let dv=[h[0]*dq[0]+h[1]*dq[1],h[1]*dq[0]+h[2]*dq[1]];
            let mut t=[0.0,1.0];
            if !clip_interval(v0[0]-v0[1],dv[0]-dv[1],high.low,high.high,&mut t)
                ||!clip_interval(v0[1],dv[1],low.low,low.high,&mut t) {continue;}
            if t[1]-t[0]<1e-14 {continue;}
            let width=t[1]-t[0];let weight=0.5*(t[1]*t[1]-t[0]*t[0]);
            covered+=width;
            mean[1+k]+=width*v0[0]+weight*dv[0];
            mean[3+k]+=width*v0[1]+weight*dv[1];
            jac[1+k]+=weight*(h[0]*d[1+k]+h[1]*d[3+k]);
            jac[3+k]+=weight*(h[1]*d[1+k]+h[2]*d[3+k]);
        }}
        close(covered,1.0,3e-12);
    }
    (mean,jac)
}
fn mean(storage:&ChargeCommutationStorage,a:&[f64;5],b:&[f64;5])->[f64;5] {
    let mid=core::array::from_fn::<_,5,_>(|i|f64::midpoint(a[i],b[i]));
    let mut out=[0.0;5];
    ConstitutiveDiscreteGradient::evaluate(storage,a,b,&mid,&mut out);
    out
}
#[test]
fn charge_coordinates_and_true_energy_follow_four_original_nonconstant_device_laws() {
    let c=[0.1,0.2,0.3,0.4];let p=parameters([0.5,0.7],c);
    let x=p.state(1.2,[4.0,-2.0],[-0.6,3.0]).unwrap();
    let expected=expected_state(&p,c,1.2,[4.0,-2.0],[-0.6,3.0]);
    for (a,b) in x.into_iter().zip(expected) {close(a,b,3e-14);}
    let energy: f64=parts(&p,c,&x).iter().sum();
    close(energy,13.204,3e-12);
    close(p.stored_energy_j(&x).unwrap(),energy,3e-12);
    let storage=ChargeCommutationStorage::new(&p).unwrap();let mut effort=[0.0;5];
    storage.gradient(&x,&mut effort);
    for (actual,expected) in effort.into_iter().zip([1.2,4.0,-2.0,-0.6,3.0]) {close(actual,expected,3e-13);}
    // This fixture includes two reverse-polarity capacitor voltages. Neither
    // their sign nor the nonlinear charge/energy datum may be discarded.
    assert!(expected[2]<0.0&&expected[3]<0.0&&expected[4]>0.0);
}

#[test]
fn exact_avf_and_end_jacobian_cross_original_cells_without_global_energy_correction() {
    let c=[0.1,0.2,0.3,0.4];let p=parameters([0.5,0.7],c);
    let a=p.state(1.2,[4.0,-2.0],[-0.6,3.0]).unwrap();
    let b=p.state(-0.8,[-3.0,5.0],[2.0,-2.5]).unwrap();
    let direction=[0.3,-0.7,0.2,0.4,-0.6];
    let storage=ChargeCommutationStorage::new(&p).unwrap();
    let (expected,jacobian)=analytic_chord(&p,c,&a,&b,&direction);
    let actual=mean(&storage,&a,&b);
    let reversed=mean(&storage,&b,&a);
    for i in 0..5 {close(actual[i],expected[i],2e-11);close(actual[i],reversed[i],2e-11);}
    let before: f64=parts(&p,c,&a).iter().sum();let after: f64=parts(&p,c,&b).iter().sum();
    let work: f64=(0..5).map(|i|actual[i]*(b[i]-a[i])).sum();
    close(work,after-before,3e-11);
    let mut tangent=[0.0;5];
    assert!(ConstitutiveDiscreteGradient::end_jacobian_action(&storage,&a,&b,&direction,&mut tangent));
    let epsilon=1e-6;
    let plus=mean(&storage,&a,&core::array::from_fn(|i|b[i]+epsilon*direction[i]));
    let minus=mean(&storage,&a,&core::array::from_fn(|i|b[i]-epsilon*direction[i]));
    for i in 0..5 {
        close(tangent[i],jacobian[i],2e-11);
        close(tangent[i],(plus[i]-minus[i])/(2.0*epsilon),2e-8);
    }
    let mid: [f64; 5] = core::array::from_fn(|i|f64::midpoint(a[i],b[i]));
    let physical_mid=physical(&p,c,&mid);
    let midpoint_work: f64=(0..5).map(|i|physical_mid[i]*(b[i]-a[i])).sum();
    close(storage.midpoint_energy_remainder(&a,&b,&mid).unwrap(),after-before-midpoint_work,3e-11);
    let mut continuous=[0.0;5];
    assert!(ContinuousStorageTangent::apply(&storage,&b,&direction,&mut continuous));
    let plus=physical(&p,c,&core::array::from_fn(|i|b[i]+epsilon*direction[i]));
    let minus=physical(&p,c,&core::array::from_fn(|i|b[i]-epsilon*direction[i]));
    for i in 0..5 {close(continuous[i],(plus[i]-minus[i])/(2.0*epsilon),2e-8);}
    let mut limited=p.clone();limited.maximum_cells_per_chord=1;
    let limited=ChargeCommutationStorage::new(&limited).unwrap();
    assert!(!mean(&limited,&a,&b).iter().all(|v|v.is_finite()));
    assert!(!ConstitutiveDiscreteGradient::end_jacobian_action(&limited,&a,&b,&direction,&mut tangent));
}

#[test]
fn supported_endpoints_do_not_admit_an_unsupported_charge_chord_or_a_nonsmooth_tangent() {
    let p=parameters([1.0;2],[1.0;4]);
    let a=p.state(0.0,[20.4,4.0],[0.5,2.0]).unwrap();
    let b=p.state(0.0,[21.9,4.0],[2.0,2.0]).unwrap();
    let storage=ChargeCommutationStorage::new(&p).unwrap();
    assert!(storage.hamiltonian(&a).is_finite()&&storage.hamiltonian(&b).is_finite());
    // Both high-device endpoint voltages are19.9V<20V. The electrode-charge
    // midpoint instead implies Vlow=27/28,Vhigh=1403/70>20V.
    let mid: [f64; 5] = core::array::from_fn(|i|f64::midpoint(a[i],b[i]));
    assert!(!storage.hamiltonian(&mid).is_finite());
    assert!(!mean(&storage,&a,&b).iter().all(|v|v.is_finite()));
    let mut tangent=[0.0;5];
    assert!(!ConstitutiveDiscreteGradient::end_jacobian_action(&storage,&a,&b,&[1.0;5],&mut tangent));
    let p=parameters([0.5,0.75],[0.25,0.5,0.125,0.25]);
    let kink=p.state(0.0,[3.0,4.0],[2.0,2.0]).unwrap();
    let storage=ChargeCommutationStorage::new(&p).unwrap();
    assert!(!ContinuousStorageTangent::apply(&storage,&kink,&[1.0;5],&mut tangent));
}

#[test]
fn final_chord_span_after_a_knee_keeps_the_original_avf_and_end_jacobian() {
    let cap = [0.02, 0.03, 0.025, 0.015];
    let mut parameters = parameters([0.04, 0.05], cap);
    parameters.maximum_cells_per_chord = 2;
    // This supported Newton chord crosses exactly one input-high source knee.
    // Its DD crossing coordinate has a nonzero low word. Adding the remaining
    // span to it need not reproduce an exact DD one, although the span itself
    // is the original, exact terminal remainder.
    let a = [0.0, 0.533, 0.43, 0.32149999999999995, -0.02];
    let b = [0.0, 0.5285254119813318, 0.43, 0.3291129534722174, -0.02];
    let direction = [0.0, 0.7, 0.0, -0.2, 0.0];
    let before = physical(&parameters, cap, &a);
    let after = physical(&parameters, cap, &b);
    assert!(before[1] - before[3] > 1.0);
    assert!(after[1] - after[3] < 1.0);
    let storage = ChargeCommutationStorage::new(&parameters).unwrap();
    let (expected, expected_tangent) =
        analytic_chord(&parameters, cap, &a, &b, &direction);
    let actual = mean(&storage, &a, &b);
    let reversed = mean(&storage, &b, &a);
    let mut tangent = [0.0; 5];
    assert!(ConstitutiveDiscreteGradient::end_jacobian_action(
        &storage, &a, &b, &direction, &mut tangent,
    ));
    for i in 0..5 {
        close(actual[i], expected[i], 2e-11);
        close(reversed[i], expected[i], 2e-11);
        close(tangent[i], expected_tangent[i], 2e-11);
    }
    let energy_change: f64 = parts(&parameters, cap, &b)
        .into_iter()
        .zip(parts(&parameters, cap, &a))
        .map(|(after, before)| after - before)
        .sum();
    let work: f64 = (0..5).map(|i| actual[i] * (b[i] - a[i])).sum();
    close(work, energy_change, 3e-12);
}
