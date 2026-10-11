//! Full-coordinate oracles for the explicitly reduced acoustic port chart.
use super::*;
use super::super::engine::radiation::{Model,Pole};

fn grid()->Vec<f64> {
    (0..65).map(|i|std::f64::consts::TAU*(40.+360.*i as f64/64.)).collect()
}
fn directions(n:usize)->(Vec<f64>,Vec<f64>) {
    let scale=1./(n as f64).sqrt();
    ((0..n).map(|i|if i%2==0 {scale}else{-scale}).collect(),
     (0..n).map(|i|if i%4<2 {scale}else{-scale}).collect())
}
fn exact_rank_two(n:usize)->(Vec<Vec<f64>>,Vec<f64>) {
    let (a,b)=directions(n);
    (a.iter().zip(&b).map(|(a,b)|vec![4.*a,*b,0.,0.]).collect(),vec![1.,1.,2.,4.])
}
fn blank_fit(model:Model)->Fit {
    Fit {model,peak_error:0.,rms_error:0.,resistance_peak_error:0.,resistance_rms_error:0.}
}

#[test]
fn geometry_selects_deterministic_signed_ports_without_removing_board_coordinates() {
    for n in [36,68,128] {
        let (weights,areas)=exact_rank_two(n);
        let (basis,error)=geometry_basis(&weights,&areas,8).unwrap();
        let (again,repeated)=geometry_basis(&weights,&areas,8).unwrap();
        assert_eq!(basis.board_ports(),n);assert_eq!(basis.ports(),2);
        assert_eq!(basis.vectors(),again.vectors());assert_eq!(error.to_bits(),repeated.to_bits());
        assert!(error<1e-12,"n={n}: {error:e}");
        for row in basis.vectors() {
            assert!(row.iter().any(|x|*x<0.) && row.iter().any(|x|*x>0.));
        }
        let (one,omitted)=geometry_basis(&weights,&areas,1).unwrap();
        assert_eq!(one.ports(),1);
        assert!((omitted-1./17.0_f64.sqrt()).abs()<1e-12);
    }
}

#[test]
fn signed_rank_two_complete_load_fits_above_the_old_structural_limit() {
    let n=36;let omega=grid();let (a,b)=directions(n);
    let (weights,areas)=exact_rank_two(n);let (basis,_)=geometry_basis(&weights,&areas,2).unwrap();
    // Independent analytic relaxation loads on signed physical directions.
    // No call to project/lift or to the fitter constructs this full oracle.
    let data:Vec<Vec<_>>=omega.iter().map(|&w| {
        let s=C64::new(0.,-w);
        let h1=s/(s+C64::new(6800.,0.));let h2=s/(s+C64::new(10000.,0.));
        (0..n*n).map(|k|h1.scale(a[k/n]*a[k%n])+h2.scale(0.7*b[k/n]*b[k%n])).collect()
    }).collect();
    let fitted=fit_projected(&omega,&data,&basis).unwrap();
    assert_eq!(fitted.model.ports,2);assert!(fitted.rms_error<0.02);
    assert!(fitted.resistance_rms_error<0.05);
    let index=31;let full=basis.lift_impedance(&fitted.model.impedance(omega[index]).unwrap()).unwrap();
    let velocity:Vec<_>=(0..n).map(|i|C64::new(a[i]+0.3*b[i],b[i]-0.2*a[i])).collect();
    let power=|z:&[C64]|(0..n).fold(C64::ZERO,|sum,i|sum+(0..n)
        .fold(C64::ZERO,|row,j|row+velocity[i].conj()*z[i*n+j]*velocity[j])).re*0.5;
    let exact=power(&data[index]);let actual=power(&full);
    assert!(actual>0. && exact>0.);
    assert!((actual-exact).abs()<0.05*exact,"{actual:e} versus {exact:e}");
    assert!(full.iter().any(|z|z.re<0.),"signed cross-port loads must survive");
}

#[test]
fn omitted_weak_resistance_cannot_hide_under_dominant_reactive_storage() {
    let n=36;let omega=grid();let (a,b)=directions(n);
    let basis=PortBasis::new(n,vec![a.clone()]).unwrap();
    let pole=8.*omega[omega.len()-1];
    // The retained image is exactly lossless. Its imaginary load dominates the
    // norm, while the small omitted direction owns ALL radiation resistance.
    let mut fitted=blank_fit(Model {ports:1,poles:vec![Pole {
        omega:pole,zeta:0.,coupling:vec![pole],
    }]});
    let data:Vec<Vec<_>>=omega.iter().map(|&w| {
        let retained=C64::new(0.,-pole*pole*w/(pole*pole-w*w));
        let s=C64::new(0.,-w);let omitted=s/(s+C64::new(6800.,0.));
        (0..n*n).map(|k|retained.scale(a[k/n]*a[k%n])+omitted.scale(b[k/n]*b[k%n])).collect()
    }).collect();
    let mut relative_complex=0.0_f64;
    for (f,&w) in omega.iter().enumerate() {
        let restored=basis.lift_impedance(&fitted.model.impedance(w).unwrap()).unwrap();
        let delta=restored.iter().zip(&data[f]).fold(0.0_f64,|sum,(x,y)|sum.hypot((*x-*y).abs()));
        let reference=data[f].iter().fold(0.0_f64,|sum,z|sum.hypot(z.abs()));
        relative_complex=relative_complex.max(delta/reference);
    }
    assert!(relative_complex<0.001);
    let error=validate_complete_fit(&omega,&data,&basis,&mut fitted).unwrap_err();
    assert!(error.contains("power peak/RMS=1.00000/1.00000"),"{error}");
}

#[test]
fn held_out_omitted_direction_is_checked_in_the_complete_board_basis() {
    let n=36;let omega=grid();let (a,b)=directions(n);
    let basis=PortBasis::new(n,vec![a.clone()]).unwrap();
    let pole=8.*omega[omega.len()-1];
    let model=Model {ports:1,poles:vec![Pole {omega:pole,zeta:0.2,coupling:vec![200.]}]};
    let mut fitted=blank_fit(model);
    let data:Vec<Vec<_>>=omega.iter().enumerate().map(|(f,&w)| {
        let h=C64::new(0.,-40000.*w)/C64::new(pole*pole-w*w,-0.4*pole*w);
        (0..n*n).map(|k|h.scale(a[k/n]*a[k%n]+if f%2==1 {b[k/n]*b[k%n]}else{0.})).collect()
    }).collect();
    assert!(validate_complete_fit(&omega,&data,&basis,&mut fitted).unwrap_err()
        .contains("complete board"));
}

#[test]
fn zero_surface_is_exact_and_incomplete_or_unbounded_geometry_refuses() {
    let (weights,areas)=exact_rank_two(36);
    for ports in [0,33,usize::MAX] {assert!(geometry_basis(&weights,&areas,ports).is_err());}
    assert!(geometry_basis(&[],&areas,2).is_err());
    assert!(geometry_basis(&vec![vec![0.;4];129],&areas,2).is_err());
    assert!(geometry_basis(&weights,&[1.,1.,2.],2).is_err());
    assert!(geometry_basis(&weights,&[1.,1.,f64::NAN,4.],2).is_err());
    let mut bad=weights.clone();bad[35][2]=f64::INFINITY;
    assert!(geometry_basis(&bad,&areas,2).is_err());
    let (basis,error)=geometry_basis(&vec![vec![0.;4];36],&areas,32).unwrap();
    assert_eq!(basis.board_ports(),36);assert_eq!(basis.ports(),1);assert_eq!(error,0.);
    let omega=grid();let zero=vec![vec![C64::ZERO;36*36];omega.len()];
    let fitted=fit_projected(&omega,&zero,&basis).unwrap();
    assert_eq!(fitted.peak_error,0.);assert_eq!(fitted.resistance_rms_error,0.);
}
