use super::*;
fn budget()->Budget {Budget {maximum_unknowns:1024,maximum_cache_bytes:64*1024*1024,
    maximum_work:20_000_000_000,angular_refinement:0}}
fn medium()->Medium {Medium {density:1.0,sound_speed:1.0}}
fn fields(s:&SpherePanels)->Vec<Vec<C64>> {
    let mut fields=vec![vec![C64::ONE;s.areas().len()]];
    for d in 0..3 {fields.push(s.normals().iter().map(|n|C64::from_re(n[d])).collect());}
    fields
}
fn impedance(s:&SpherePanels,v:&[Vec<C64>],p:&[Solution<'_,'_>])->Vec<C64> {
    let mut out=vec![C64::ZERO;v.len()*v.len()];
    for i in 0..v.len() {for j in 0..v.len() {
        out[i*v.len()+j]=v[i].iter().zip(&p[j].pressure).zip(s.areas())
            .fold(C64::ZERO,|a,((v,p),s)|a+v.conj()*p.scale(*s));
    }}out
}

#[test]
fn g1_sphere_impedance_and_point_pressure_converge_without_authored_resistance() {
    let mut previous=f64::INFINITY;
    for level in [0,1] {
        let surface=SpherePanels::icosphere(1.0,level).unwrap();
        let geometry=Prepared::new(&surface,1.0,budget(),||false).unwrap();
        let velocities=vec![C64::ONE;surface.areas().len()];
        let mut worst=0.0_f64;
        for k in [0.1,0.5,1.0] {
            let p=geometry.solve_batch(k,medium(),&[&velocities],||false).unwrap();
            let area: f64=surface.areas().iter().sum();
            let actual=p[0].pressure.iter().zip(surface.areas())
                .fold(C64::ZERO,|s,(p,a)|s+p.scale(*a)).scale(1.0/area);
            let exact=C64::new(0.0,k)*C64::new(-1.0,k).recip();
            worst=worst.max((actual-exact).abs()/exact.abs());
            assert!(p[0].backward_residual<1e-10 && p[0].flux_balance<1e-8);
            assert!(p[0].surface_power.re>0.0 && p[0].outgoing_power_w>0.0);
        }
        assert!(worst<if level==0 {0.13}else{0.04});
        assert!(worst<0.4*previous);previous=worst;
        if level==1 {
            let p=geometry.solve_batch(0.5,medium(),&[&velocities],||false).unwrap();
            let exact=C64::new(0.0,0.5)*C64::new(-1.0,0.5).recip()
                *C64::new(det::cos(1.0),det::sin(1.0)).scale(1.0/3.0);
            let actual=p[0].pressure_at([0.0,0.0,3.0]).unwrap();
            assert!((actual-exact).abs()<0.08*exact.abs());
        }
    }
}

#[test]
fn g0_full_complex_multiport_is_reciprocal_and_passive_including_weak_radiators() {
    let surface=SpherePanels::icosphere(1.0,0).unwrap();
    let geometry=Prepared::new(&surface,1.0,budget(),||false).unwrap();
    let mut v=fields(&surface);
    // Signed nonuniform sources and an exact linearly dependent source. No
    // modal rank truncation or positive diagonal loading may remove these.
    v.push(surface.normals().iter().map(|n|C64::from_re(n[0]*n[1])).collect());
    v.push(v[1].iter().map(|x|x.scale(-0.7)).collect());
    let refs:Vec<_>=v.iter().map(Vec::as_slice).collect();
    for k in [0.02,0.1,0.5,1.0] {
        let p=geometry.solve_batch(k,medium(),&refs,||false).unwrap();
        let z=impedance(&surface,&v,&p);let n=v.len();
        let scale=z.iter().map(|z|z.abs()).fold(0.0_f64,f64::max);
        let mut h=vec![0.0;4*n*n];
        for i in 0..n {for j in 0..n {
            assert!((z[i*n+j]-z[j*n+i]).abs()<1e-10*scale);
            let value=(z[i*n+j]+z[j*n+i].conj()).scale(0.5);
            h[i*2*n+j]=value.re;h[(i+n)*2*n+j+n]=value.re;
            h[i*2*n+j+n]=-value.im;h[(i+n)*2*n+j]=value.im;
        }}
        let (values,_)=fs_la::eigen::jacobi_eigh(&h,2*n);
        assert!(values.iter().all(|v|*v>=-1e-11*scale));
        assert!(z[n+5].re<0.0,"signed cross resistance must survive");
    }
}

#[test]
fn g3_phase_mixtures_close_surface_work_against_independent_far_field_quadrature() {
    let surface=SpherePanels::icosphere(1.0,0).unwrap();
    let geometry=Prepared::new(&surface,1.0,budget(),||false).unwrap();
    let velocity:Vec<_>=surface.normals().iter().map(|n|C64::new(n[0]+0.2,n[2]-0.4*n[1])).collect();
    let p=geometry.solve_batch(0.8,medium(),&[&velocity],||false).unwrap();
    let mut integral=0.0;
    // Twice the assembly's azimuthal resolution and different sample phase.
    for (mu,w) in X8.into_iter().zip(W8) {for i in 0..32 {
        let mu=2.0*mu-1.0;let phi=TAU*i as f64/32.0;
        let d=[(1.0-mu*mu).sqrt()*det::cos(phi),(1.0-mu*mu).sqrt()*det::sin(phi),mu];
        let field=p[0].far_field(d).unwrap();integral+=2.0*w*TAU/32.0*field.abs().powi(2)/2.0;
    }}
    assert!((integral-p[0].surface_power.re).abs()<1e-7*integral);
    assert!(p[0].pressure_at([0.0;3]).is_err());
    assert!(p[0].far_field([0.0;3]).is_err());
}

#[test]
fn g0_linear_static_moments_match_independent_gauss_integrals() {
    let surface=SpherePanels::icosphere(1.0,0).unwrap();
    let (facets,_,_,_)=facets(&surface).unwrap();let facet=&facets[0];
    let x=[2.3,-1.5,3.2];let (s,sm,dm)=static_moments(x,facet,false);
    let (mut single,mut moments,mut double)=(0.0,[0.0;3],[0.0;3]);
    for (u,wu) in X8.into_iter().zip(W8) {for (v,wv) in X8.into_iter().zip(W8) {
        let shape=[(1.0-u)*(1.0-v),u,(1.0-u)*v];
        let point=core::array::from_fn(|d|(0..3).map(|a|shape[a]*facet.triangle[a][d]).sum());
        let delta=sub(x,point);let r=norm(delta);let weight=2.0*facet.area*(1.0-u)*wu*wv;
        let g=weight/(4.0*PI*r);let d=g*dot(facet.normal,delta)/(r*r);
        single+=g;for a in 0..3 {moments[a]+=g*shape[a];double[a]+=d*shape[a];}
    }}
    assert!((s-single).abs()<1e-11);
    for a in 0..3 {assert!((sm[a]-moments[a]).abs()<1e-11);assert!((dm[a]-double[a]).abs()<1e-11);}
    assert!((sm.iter().sum::<f64>()-s).abs()<1e-13);
}

#[test]
fn g4_refused_preparation_and_batches_do_not_poison_retry_or_publish_partial_sources() {
    let surface=SpherePanels::icosphere(1.0,0).unwrap();
    assert!(Prepared::new(&surface,1.0,Budget{maximum_cache_bytes:1,..budget()},||false).is_err());
    assert!(Prepared::new(&surface,1.0,Budget{maximum_work:1,..budget()},||false).is_err());
    assert!(Prepared::new(&surface,1.0,Budget{maximum_unknowns:20,..budget()},||false).is_err());
    assert!(Prepared::new(&surface,2.0,budget(),||false).is_err());
    let mut calls=0;
    assert!(matches!(Prepared::new(&surface,1.0,budget(),||{calls+=1;calls==3}),Err(HelmholtzError::Cancelled)));
    let geometry=Prepared::new(&surface,1.0,budget(),||false).unwrap();
    let v=vec![C64::ONE;surface.areas().len()];let first=geometry.solve_batch(0.5,medium(),&[&v],||false).unwrap();
    let mut invalid=v.clone();invalid[3]=C64::new(f64::NAN,0.0);
    assert!(geometry.solve_batch(0.5,medium(),&[&v,&invalid],||false).is_err());
    assert!(geometry.solve_batch(1.1,medium(),&[&v],||false).is_err());
    assert!(matches!(geometry.solve_batch(0.5,medium(),&[&v],||true),Err(HelmholtzError::Cancelled)));
    let retry=geometry.solve_batch(0.5,medium(),&[&v],||false).unwrap();
    assert_eq!(first[0].pressure,retry[0].pressure);
    let zero=vec![C64::ZERO;surface.areas().len()];let rest=geometry.solve_batch(0.5,medium(),&[&zero],||false).unwrap();
    assert_eq!(rest[0].outgoing_power_w,0.0);assert!(rest[0].pressure.iter().all(|p|p.abs()==0.0));
}

#[test]
fn tiled_gram_term_is_bit_identical_to_the_direction_by_direction_sweep() {
    let surface=SpherePanels::icosphere(1.0,2).unwrap();
    let geometry=Prepared::new(&surface,1.0,budget(),||false).unwrap();
    let n=geometry.unknowns();
    assert!(n>2*GRAM_TILE && n%GRAM_TILE!=0,"several tiles and a partial one");
    let k=0.7;
    let waves=geometry.waves(k,&mut ||false).unwrap();
    // A symmetric start whose imaginary parts are already nonzero, so every
    // accumulation step shows in the bits.
    let mut tiled:Vec<C64>=(0..n*n).map(|e|{
        let (a,b)=((e/n).min(e%n),(e/n).max(e%n));
        C64::new(((7*a+3*b)%101) as f64*0.01-0.5,((5*a+11*b)%97) as f64*1e-3)
    }).collect();
    let mut sweep=tiled.clone();
    // The direction-by-direction sweep add_gram replaced, verbatim.
    for wave in &waves {
        let scale=k*wave.weight/(8.0*PI*PI);
        for i in 0..n {for j in i..n {
            sweep[i*n+j].im+=scale*(wave.coefficients[i].conj()*wave.coefficients[j]).re;
            sweep[j*n+i]=sweep[i*n+j];
        }}
    }
    let mut polls=0;
    add_gram(&mut tiled,n,k,&waves,&mut ||{polls+=1;false}).unwrap();
    assert_eq!(tiled.iter().map(|c|(c.re.to_bits(),c.im.to_bits())).collect::<Vec<_>>(),
        sweep.iter().map(|c|(c.re.to_bits(),c.im.to_bits())).collect::<Vec<_>>());
    let tiles=n.div_ceil(GRAM_TILE);
    assert_eq!(polls,tiles*(tiles+1)/2+n);
    let mut calls=0;
    assert!(matches!(add_gram(&mut tiled,n,k,&waves,&mut ||{calls+=1;calls==2}),Err(HelmholtzError::Cancelled)));
}
