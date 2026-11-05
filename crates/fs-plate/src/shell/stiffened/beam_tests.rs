use super::*;
fn spec() -> RoundBeamSpec {
    RoundBeamSpec {young_pa:12e9,density_kg_m3:800.,pivot_m:0.12,contact_m:0.39,
        hand_m:0.2,subdivisions:12,maximum_hz:1500.,maximum_modes:12}
}
fn shaft() -> [RoundStation;2] { [RoundStation{x_m:0.,radius_m:0.005},RoundStation{x_m:0.4,radius_m:0.005}] }
fn quad(a:&[f64;16],x:[f64;4]) -> f64 { (0..4).map(|i|x[i]*(0..4).map(|j|a[4*i+j]*x[j]).sum::<f64>()).sum() }

#[test]
fn tapered_element_keeps_rigid_kinetic_energy_and_exact_bending_work() {
    let l=0.4; let radius=0.005; let e=12e9; let rho=800.;
    let (k,m)=round_element(l,[radius;2],e,rho).unwrap();
    let area=std::f64::consts::PI*radius*radius; let inertia=area*radius*radius/4.;
    // Constant curvature w=x^2/2, and rigid translation/rotation respectively.
    assert!((quad(&k,[0.,0.,l*l/2.,l])/(e*inertia*l)-1.).abs()<1e-12);
    assert!((quad(&m,[1.,0.,1.,0.])/(rho*area*l)-1.).abs()<1e-13);
    assert!((quad(&m,[0.,1.,l,1.])/(rho*(area*l*l*l/3.+inertia*l))-1.).abs()<1e-13);
    assert!(quad(&k,[0.,1.,l,1.]).abs()<1e-10);
    let (kt,mt)=round_element(l,[radius,radius/2.],e,rho).unwrap();
    assert!(quad(&mt,[1.,0.,1.,0.])<quad(&m,[1.,0.,1.,0.]));
    assert!(quad(&kt,[0.,0.,l*l/2.,l])<quad(&k,[0.,0.,l*l/2.,l]));
}
#[test]
fn original_profile_mass_and_pinned_force_ports_are_not_tip_mass_proxies() {
    let s=spec(); let b=RoundBeamModes::new(&shaft(),s).unwrap();
    let l=0.4; let radius=0.005; let area=std::f64::consts::PI*radius*radius;
    let inertia=s.density_kg_m3*area*((l-s.pivot_m).powi(3)+s.pivot_m.powi(3))/3.
        +s.density_kg_m3*area*radius*radius*l/4.;
    assert!((b.pivot_inertia_kg_m2/inertia-1.).abs()<1e-12);
    assert_eq!(b.omega[0],0.); assert!(b.omega.len()>2);
    assert!((b.tip[0]-(s.contact_m-s.pivot_m)/inertia.sqrt()).abs()<1e-12);
    assert_eq!(b.point(s.pivot_m).unwrap(),vec![0.;b.omega.len()]);
    assert!(b.tip.iter().zip(&b.hand).skip(1).any(|(t,h)|(t-h).abs()>1e-3));
    for i in 0..b.omega.len() { assert!((b.point(s.contact_m).unwrap()[i]-b.tip[i]).abs()<1e-12); }
    let qdot:Vec<_>=(0..b.omega.len()).map(|i|(i as f64+1.)*0.03).collect();
    let force=1.7; let v=b.tip.iter().zip(&qdot).map(|(b,v)|b*v).sum::<f64>();
    assert!((force*v-b.tip.iter().zip(&qdot).map(|(b,v)|b*force*v).sum::<f64>()).abs()<1e-13);
}
#[test]
fn pin_free_uniform_beam_converges_to_its_independent_elastic_root() {
    // No-shear slender limit: tan(beta L)=tanh(beta L), beta L=3.926602312...
    // Rotary inertia is retained; choose r/L=0.00025 so its difference is small.
    let shaft=[RoundStation{x_m:0.,radius_m:0.0001},RoundStation{x_m:0.4,radius_m:0.0001}];
    let mut s=spec();s.pivot_m=0.;s.maximum_hz=4.;s.maximum_modes=2;
    let expected=3.926602312047919_f64.powi(2)/0.4_f64.powi(2)*(s.young_pa*0.0001_f64.powi(2)/(4.*s.density_kg_m3)).sqrt();
    s.subdivisions=4; let a=RoundBeamModes::new(&shaft,s).unwrap();
    s.subdivisions=8; let b=RoundBeamModes::new(&shaft,s).unwrap();
    assert!((b.omega[1]/expected-1.).abs()<5e-5);
    assert!((b.omega[1]-expected).abs()<0.1*(a.omega[1]-expected).abs());
}
#[test]
fn material_scaling_moves_elastic_frequencies_without_relabeling_the_rigid_mode() {
    let s=spec();let a=RoundBeamModes::new(&shaft(),s).unwrap();
    let mut t=s;t.young_pa*=4.;t.maximum_hz*=2.;let b=RoundBeamModes::new(&shaft(),t).unwrap();
    assert_eq!(a.pivot_inertia_kg_m2,b.pivot_inertia_kg_m2);
    assert_eq!(a.tip[0],b.tip[0]);assert_eq!(a.omega.len(),b.omega.len());
    for (a,b) in a.omega.iter().zip(&b.omega).skip(1) {assert!((b/(2.*a)-1.).abs()<1e-6);}
    let mut t=s;t.density_kg_m3*=4.;t.maximum_hz/=2.;let b=RoundBeamModes::new(&shaft(),t).unwrap();
    assert!((b.tip[0]/a.tip[0]-0.5).abs()<1e-12);
    for (a,b) in a.omega.iter().zip(&b.omega).skip(1) {assert!((b/(0.5*a)-1.).abs()<1e-6);}
}
#[test]
fn missing_or_excess_elastic_modes_and_invalid_geometry_refuse_without_truncation() {
    for change in 0..6 {
        let mut s=spec();match change {0=>s.pivot_m=0.4,1=>s.maximum_hz=1.,2=>s.maximum_modes=2,
            3=>s.subdivisions=0,4=>s.young_pa=f64::NAN,_=>s.hand_m=s.pivot_m}
        assert!(RoundBeamModes::new(&shaft(),s).is_err());
    }
    assert!(round_element(1.,[0.,0.],1.,1.).is_err());
    assert!(RoundBeamModes::new(&[RoundStation{x_m:0.,radius_m:0.},RoundStation{x_m:0.4,radius_m:0.}],spec()).is_err());
}

#[test]
fn attached_head_changes_the_source_spectrum_and_keeps_loaded_mass_orthonormality() {
    let s=spec(); let a=RoundInertia{x_m:s.contact_m,mass_kg:0.02,rotary_kg_m2:1.2e-6};
    let bare=RoundBeamModes::new(&shaft(),s).unwrap();
    let b=RoundBeamModes::with_inertias(&shaft(),s,&[a]).unwrap();
    let expected=bare.pivot_inertia_kg_m2+a.mass_kg*(a.x_m-s.pivot_m).powi(2)+a.rotary_kg_m2;
    assert!((b.pivot_inertia_kg_m2/expected-1.).abs()<1e-12);
    assert!(b.omega.len()>=bare.omega.len());
    assert!(b.omega[1]<0.98*bare.omega[1],"mass must enter the eigensolve, not an output scale");
    let tip=b.point(a.x_m).unwrap();let slope=b.slope(a.x_m).unwrap();
    assert!((slope[0]*expected.sqrt()-1.).abs()<1e-12);
    for i in 0..b.omega.len() {for j in 0..b.omega.len() {
        let mut product=a.mass_kg*tip[i]*tip[j]+a.rotary_kg_m2*slope[i]*slope[j];
        // Reintegrate PHYSICAL shapes in the unscaled source mesh, including
        // the whole shaft behind the pin, separately from the assembly code.
        for (cell,x) in b.nodes_m.windows(2).enumerate() {
            let (_,m)=round_element(x[1]-x[0],[0.005;2],s.young_pa,s.density_kg_m3).unwrap();
            for row in 0..4 {for col in 0..4 {
                product+=b.shapes[i][2*cell+row]*m[4*row+col]*b.shapes[j][2*cell+col];
            }}
        }
        assert!((product-if i==j{1.}else{0.}).abs()<1e-6,"loaded M[{i},{j}]={product}");
    }}
}

#[test]
fn empty_inertia_is_identical_and_invalid_physical_attachments_refuse() {
    let s=spec();let a=RoundBeamModes::new(&shaft(),s).unwrap();
    let b=RoundBeamModes::with_inertias(&shaft(),s,&[]).unwrap();
    assert_eq!(a.omega,b.omega);assert_eq!(a.shapes,b.shapes);assert_eq!(a.tip,b.tip);
    let good=RoundInertia{x_m:0.39,mass_kg:0.02,rotary_kg_m2:1e-6};
    for bad in [RoundInertia{x_m:0.5,..good},RoundInertia{mass_kg:-1.,..good},
        RoundInertia{rotary_kg_m2:f64::NAN,..good},RoundInertia{mass_kg:0.,rotary_kg_m2:0.,..good}] {
        assert!(RoundBeamModes::with_inertias(&shaft(),s,&[bad]).is_err());
    }
    assert!(RoundBeamModes::with_inertias(&shaft(),s,&[good;9]).is_err());
    assert!(b.slope(f64::NAN).is_err());
}

#[test]
fn graded_tapered_element_matches_independent_polynomial_kinetic_integrals() {
    // Coefficient integration is independent of the production Gauss rule.
    fn product(a: &[f64], b: &[f64]) -> Vec<f64> {
        let mut c=vec![0.;a.len()+b.len()-1];
        for (i,a) in a.iter().enumerate() {for (j,b) in b.iter().enumerate() {c[i+j]+=a*b;}}
        c
    }
    fn integral(a: &[f64]) -> f64 {a.iter().enumerate().map(|(i,a)|a/(i+1) as f64).sum()}
    let l=0.7; let e=12e9; let radii=[0.003,0.009]; let rho=[430.,1800.];
    let (k,m)=round_element_with_reference_density(l,radii,e,rho).unwrap();
    assert_eq!(k,round_element(l,radii,e,800.).unwrap().0);
    let r2=product(&[radii[0],radii[1]-radii[0]],&[radii[0],radii[1]-radii[0]]);
    let r4=product(&r2,&r2); let density=[rho[0],rho[1]-rho[0]];
    let area=product(&density,&r2); let inertia=product(&density,&r4);
    let n=[[1.,0.,-3.,2.],[0.,l,-2.*l,l],[0.,0.,3.,-2.],[0.,0.,-l,l]];
    let d=[[0.,-6./l,6./l],[1.,-4.,3.],[0.,6./l,-6./l],[0.,-2.,3.]];
    for i in 0..4 {for j in 0..4 {
        let expected=std::f64::consts::PI*l*(integral(&product(&area,&product(&n[i],&n[j])))
            +0.25*integral(&product(&inertia,&product(&d[i],&d[j]))));
        assert!((m[4*i+j]-expected).abs()<2e-12*expected.abs().max(1e-8),"M[{i},{j}]: {} != {expected}",m[4*i+j]);
    }}
}

#[test]
fn uniform_reference_profile_keeps_original_modes_and_attached_inertia_bits() {
    let s=spec(); let attached=[RoundInertia{x_m:0.35,mass_kg:0.01,rotary_kg_m2:2e-7}];
    let original=RoundBeamModes::with_inertias(&shaft(),s,&attached).unwrap();
    let explicit=RoundBeamModes::with_reference_density(&shaft(),s,&[[s.density_kg_m3;2]],&attached).unwrap();
    assert_eq!(original.omega,explicit.omega); assert_eq!(original.shapes,explicit.shapes);
    assert_eq!(original.pivot_inertia_kg_m2,explicit.pivot_inertia_kg_m2);
    assert_eq!(original.tip,explicit.tip); assert_eq!(original.hand,explicit.hand);
    assert_eq!(round_element(0.4,[0.005,0.002],12e9,800.).unwrap(),
        round_element_with_reference_density(0.4,[0.005,0.002],12e9,[800.;2]).unwrap());
}

#[test]
fn density_interfaces_and_subdivision_preserve_the_physical_pivot_inertia() {
    let stations=[RoundStation{x_m:0.,radius_m:0.005},RoundStation{x_m:0.17,radius_m:0.005},
        RoundStation{x_m:0.4,radius_m:0.005}];
    let density=[[400.,1200.],[1900.,900.]];
    let attached=RoundInertia{x_m:0.39,mass_kg:0.012,rotary_kg_m2:1.7e-7};
    let s=spec(); let area=std::f64::consts::PI*0.005_f64.powi(2);
    let second=area*0.005_f64.powi(2)/4.;
    let mut expected=attached.mass_kg*(attached.x_m-s.pivot_m).powi(2)+attached.rotary_kg_m2;
    for (x,rho) in stations.windows(2).zip(density) {
        let l=x[1].x_m-x[0].x_m; let a=x[0].x_m-s.pivot_m; let dr=rho[1]-rho[0];
        expected+=area*l*(rho[0]*(a*a+a*l+l*l/3.)+dr*(a*a/2.+2.*a*l/3.+l*l/4.))
            +second*l*(rho[0]+rho[1])/2.;
    }
    for subdivisions in [4,8,12] {
        let mut s=s; s.subdivisions=subdivisions;
        // The resolved profile replaces the unused homogeneous density.
        s.density_kg_m3=f64::NAN;
        let modes=RoundBeamModes::with_reference_density(&stations,s,&density,&[attached]).unwrap();
        assert!((modes.pivot_inertia_kg_m2/expected-1.).abs()<2e-13);
        assert!((modes.tip[0]*expected.sqrt()-(s.contact_m-s.pivot_m)).abs()<2e-14);
        assert!((modes.slope(s.hand_m).unwrap()[0]*expected.sqrt()-1.).abs()<2e-14);
    }
}

#[test]
fn equal_total_mass_with_different_density_distribution_changes_the_spectrum() {
    let mut s=spec(); s.pivot_m=0.;
    // Both constant-area shafts have rho_mean=800 kg/m³ and identical stiffness.
    let light_tip=RoundBeamModes::with_reference_density(&shaft(),s,&[[1200.,400.]],&[]).unwrap();
    let heavy_tip=RoundBeamModes::with_reference_density(&shaft(),s,&[[400.,1200.]],&[]).unwrap();
    assert!(heavy_tip.pivot_inertia_kg_m2>1.5*light_tip.pivot_inertia_kg_m2);
    assert!((light_tip.omega[1]/heavy_tip.omega[1]-1.).abs()>0.01,
        "equal mean density must not replace the spatial inertia");
    let mut scaled=s; scaled.maximum_hz/=2.;
    let heavier=RoundBeamModes::with_reference_density(&shaft(),scaled,&[[1600.,4800.]],&[]).unwrap();
    assert_eq!(heavy_tip.omega.len(),heavier.omega.len());
    for (a,b) in heavy_tip.omega.iter().zip(&heavier.omega).skip(1) {assert!((b/(0.5*a)-1.).abs()<2e-7);}
    assert!((heavier.tip[0]/heavy_tip.tip[0]-0.5).abs()<1e-13);
}

#[test]
fn graded_density_rejects_bad_counts_values_and_overflow() {
    for density in [vec![],vec![[800.;2];2],vec![[0.,800.]],vec![[-1.,800.]],
        vec![[800.,f64::NAN]],vec![[f64::INFINITY,800.]]] {
        assert!(RoundBeamModes::with_reference_density(&shaft(),spec(),&density,&[]).is_err());
    }
    for rho in [[800.,0.],[-1.,800.],[800.,f64::NAN],[f64::INFINITY,800.]] {
        assert!(round_element_with_reference_density(0.4,[0.005;2],12e9,rho).is_err());
    }
    assert!(round_element_with_reference_density(4.,[1.;2],1.,[f64::MAX/2.,f64::MAX]).is_err());
}
