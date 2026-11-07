use super::*;

const MU:f64=1.256_637_062_12e-6;
pub(in super::super) fn fixture(ni:f64,nt:f64,rows:&[(f64,f64,f64,[f64;3])])->Frequency {
    let mut entrance=0.0;
    let layers=rows.iter().map(|&(no,ne,thickness,axis)| {
        let exit=entrance+thickness;
        let layer=ResolvedLayer {no,ne,mu:MU,axis,thickness,entrance,exit};
        entrance=exit; layer
    }).collect();
    Frequency {omega:C0,ni,nt,mui:MU,mut_:MU,layers}
}
fn close(a:f64,b:f64,t:f64) { assert!((a-b).abs()<=t,"{a:e} != {b:e} (tol {t:e})"); }
fn zclose(a:C64,b:C64,t:f64) {assert!((a-b).abs()<=t,"{a:?} != {b:?} (tol {t:e})");}
fn answer(data:&Frequency,p:f64,e:[C64;2])->UniaxialStackResponse {
    solve(data,Incidence::TangentialIndex(p),e,[1.0,0.0],&mut ||false).unwrap()
}

// Independent scalar characteristic matrices for a SHORT, propagating stack.
// These intentionally do not use the production Schur recurrence.
fn characteristic(data:&Frequency,p:f64,tm:bool)->(C64,C64) {
    let adm=|n:f64,mu:f64| {let q=(n*n-p*p).sqrt(); if tm {n*n/(q*(mu/MU))} else {q/(mu/MU)}};
    let mut m=[[C64::ONE,C64::ZERO],[C64::ZERO,C64::ONE]];
    for layer in &data.layers {
        let q=(layer.no*layer.no-p*p).sqrt();
        let y=adm(layer.no,layer.mu);
        let (s,c)=(q*layer.thickness).sin_cos();
        let next=[[C64::from_re(c),C64::new(0.0,s/y)],[C64::new(0.0,s*y),C64::from_re(c)]];
        m=core::array::from_fn(|j| core::array::from_fn(|k| m[j][0]*next[0][k]+m[j][1]*next[1][k]));
    }
    let yi=adm(data.ni,data.mui); let yt=adm(data.nt,data.mut_);
    let a=m[0][0]+m[0][1].scale(yt); let b=m[1][0]+m[1][1].scale(yt);
    ((a.scale(yi)-b)/(a.scale(yi)+b),C64::from_re(2.0*yi)/(a.scale(yi)+b))
}

#[test]
fn isotropic_multilayers_match_independent_characteristic_te_tm_with_magnetic_contrast() {
    let mut data=fixture(1.3,1.2,&[
        (1.6,1.6,0.7,[0.2,0.5,0.8]),(2.1,2.1,0.43,[1.0,0.0,0.0]),
        (1.4,1.4,0.29,[0.0,1.0,0.0]),
    ]);
    data.layers[1].mu=1.4*MU; data.mut_=0.9*MU;
    for p in [-0.8,0.0,0.65] { for tm in [false,true] {
        let j=if tm {0} else {1};
        let mut e=[C64::ZERO;2];e[j]=C64::new(0.6,-0.3);
        let result=answer(&data,p,e);
        let (r,t)=characteristic(&data,p,tm);
        zclose(result.reflected_electric_v_m[j],e[j]*r,2e-14);
        zclose(result.transmitted_electric_v_m[j],e[j]*t,2e-14);
        close(result.reflected_power_fraction,r.norm_sq(),2e-14);
        close(result.reflected_power_fraction+result.transmitted_power_fraction,1.0,3e-14);
    }}
}

#[test]
fn rotated_resonant_crystals_form_noncommuting_coherent_waveplates() {
    let angle=core::f64::consts::PI/8.0;
    let x=[1.0,0.0,0.0]; let rotated=[angle.cos(),angle.sin(),0.0];
    let d=core::f64::consts::PI;
    let e=[C64::ONE,C64::new(0.0,0.4)];
    for reverse in [false,true] {
        let axes=if reverse {[rotated,x]} else {[x,rotated]};
        let data=fixture(1.0,1.0,&[(1.0,2.0,d,axes[0]),(1.0,2.0,d,axes[1])]);
        let result=answer(&data,0.0,e);
        let (s,c)=(2.0*angle).sin_cos(); let s=if reverse {-s} else {s};
        zclose(result.transmitted_electric_v_m[0],e[0].scale(c)-e[1].scale(s),3e-14);
        zclose(result.transmitted_electric_v_m[1],e[0].scale(s)+e[1].scale(c),3e-14);
        assert!(result.reflected_power_fraction<1e-26);
        close(result.transmitted_power_fraction,1.0,2e-14);
    }
}

#[test]
fn tilted_crystal_layer_subdivision_preserves_fields_and_all_interface_fluxes() {
    let axis=[0.3,0.7,0.5];
    let original=fixture(1.8,1.4,&[(1.3,2.1,0.8,axis),(1.7,2.3,0.6,[0.6,-0.2,0.5])]);
    let split=fixture(1.8,1.4,&[(1.3,2.1,0.3,axis),(1.3,2.1,0.5,axis),(1.7,2.3,0.6,[0.6,-0.2,0.5])]);
    let e=[C64::new(0.7,0.3),C64::new(-0.2,0.6)];
    let a=answer(&original,0.9,e);let b=answer(&split,0.9,e);
    for (x,y) in a.transmitted_electric_v_m.into_iter().zip(b.transmitted_electric_v_m) {zclose(x,y,2e-14);}
    for (x,y) in a.reflected_electric_v_m.into_iter().zip(b.reflected_electric_v_m) {zclose(x,y,2e-14);}
    assert!(a.reflected_electric_v_m[1].abs()>1e-3);
    for (za,ia,ib) in [(0.15,0,0),(0.6,0,1),(1.1,1,2)] {
        let fa=a.field_at(StackRegion::Layer(ia),[0.2,0.4,za]).unwrap();
        let fb=b.field_at(StackRegion::Layer(ib),[0.2,0.4,za]).unwrap();
        for (x,y) in fa.electric_v_m.into_iter().zip(fb.electric_v_m) {zclose(x,y,3e-14);}
        for (x,y) in fa.magnetic_a_m.into_iter().zip(fb.magnetic_a_m) {zclose(x,y,3e-16);}
        close(fa.poynting_w_m2[2]/a.incident_normal_power_w_m2,a.transmitted_power_fraction,2e-13);
    }
    let flipped=fixture(1.8,1.4,&[(1.3,2.1,0.8,axis.map(|x|-x)),(1.7,2.3,0.6,[-0.6,0.2,-0.5])]);
    let c=answer(&flipped,0.9,e);
    for (x,y) in a.transmitted_electric_v_m.into_iter().zip(c.transmitted_electric_v_m) {zclose(x,y,2e-14);}
}

#[test]
fn thick_evanescent_barriers_retain_positive_tunneling_without_growing_factors() {
    for depth in [0.2,10.0,100.0] {
        let data=fixture(2.0,2.1,&[
            (1.0,1.2,depth,[0.6,0.4,0.3]),(1.1,1.25,depth,[0.2,-0.3,0.8]),
        ]);
        let result=answer(&data,1.5,[C64::new(0.6,0.3),C64::ONE]);
        assert!(result.transmitted_power_fraction>0.0);
        close(result.reflected_power_fraction+result.transmitted_power_fraction,1.0,2e-13);
        for layer in &result.layers {
            for wave in &layer.forward {assert!(wave.normal_wavenumber_rad_m.im<0.0);}
            for wave in &layer.backward {assert!(wave.normal_wavenumber_rad_m.im>0.0);}
        }
        let field=result.field_at(StackRegion::Layer(1),[0.0,0.0,1.5*depth]).unwrap();
        assert!(field.electric_v_m.iter().all(|z| finite(*z)));
    }
}

#[test]
fn many_layer_bragg_reflection_matches_closed_form_without_cubic_storage() {
    for pairs in [1,8,32] {
        let mut rows=Vec::new();
        for _ in 0..pairs {
            rows.push((2.0,2.0,core::f64::consts::PI/4.0,[0.0,1.0,0.0]));
            rows.push((1.5,1.5,core::f64::consts::PI/3.0,[1.0,0.0,0.0]));
        }
        let data=fixture(1.0,1.0,&rows);
        let result=answer(&data,0.0,[C64::ONE,C64::new(0.0,0.7)]);
        let y=(4.0_f64/3.0).powi(2*pairs);
        let expected=4.0*y/((1.0+y)*(1.0+y));
        close(result.transmitted_power_fraction/expected,1.0,5e-13);
        assert_eq!(result.layers.len(),(2*pairs) as usize);
    }
    let a=workspace_bytes(10,3).unwrap();let b=workspace_bytes(20,3).unwrap();
    assert!(b>a && b<2*a);
    assert!(workspace_bytes(usize::MAX,2).is_err());
}

#[test]
fn critical_exit_and_internal_critical_refusal_and_cancellation_are_explicit() {
    let data=fixture(2.0,1.0,&[(1.5,1.8,0.7,[0.2,0.8,0.3]),(1.4,1.6,0.4,[0.3,0.1,0.8])]);
    let e=[C64::ONE,C64::new(0.0,0.4)];
    for p in [1.0,1.1] {
        let r=answer(&data,p,e);
        assert_eq!(r.transmitted_power_fraction,0.0);
        close(r.reflected_power_fraction,1.0,2e-13);
        assert!(r.field_at(StackRegion::Transmitted,[0.0,0.0,r.thickness_m+0.5]).is_ok());
    }
    let internal=fixture(2.0,1.5,&[(1.0,1.0,0.7,[1.0,0.0,0.0])]);
    assert!(solve(&internal,Incidence::TangentialIndex(1.0),e,[1.0,0.0],&mut ||false).is_err());
    let mut calls=0;
    assert!(matches!(solve(&data,Incidence::TangentialIndex(0.5),e,[1.0,0.0],&mut ||{calls+=1;calls==5}),
        Err(crate::plane_wave::angular_spectrum::AngularSpectrumError::Cancelled)));
    let r=answer(&data,0.5,e);
    assert!(r.field_at(StackRegion::Layer(2),[0.0,0.0,0.3]).is_err());
    assert!(r.field_at(StackRegion::Layer(0),[0.0,0.0,0.8]).is_err());
    assert!(r.field_at(StackRegion::Incident,[0.0,0.0,0.1]).is_err());
}
