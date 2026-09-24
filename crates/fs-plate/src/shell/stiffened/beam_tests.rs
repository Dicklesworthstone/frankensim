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
