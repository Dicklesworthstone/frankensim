//! Synthetic constitutive checks; no measured wood, cane or metal authority.
use fs_material::plane_stress_viscosity::{MaterialPlane, PlaneStressViscosity,
    PLANE_STRESS_VISCOSITY_PROPERTIES as KEYS};
use fs_material::state_point::{ResolvedMaterialStatePoint, MaterialPropertySelection,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_matdb::{ClaimSet, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, InterpolationPolicy, UncertaintyModel};
use fs_evidence::ValidityDomain;
use fs_qty::{Dims, DynViscosity, Pressure};

fn source(values: [f64; 6], omit: Option<usize>, wrong: bool, curve: bool) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new(); let mut requirements = Vec::new();
    for (i, (key, value)) in KEYS.into_iter().zip(values).enumerate() {
        if omit == Some(i) { continue; }
        let dims = if wrong && i == 0 { Pressure::DIMS } else { DynViscosity::DIMS };
        let varying = curve && i == 0;
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims),
            value: if varying { PropertyValue::Curve { abscissa: "omega".into(),
                abscissa_dims: Dims([0,0,-1,0,0,0]), knots: vec![(10.0,value),(10000.0,2.0*value)], dims }
            } else { PropertyValue::Scalar { value, dims } },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 310.0).with("omega", 10.0, 10000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: if varying { InterpolationPolicy::LinearInside } else { InterpolationPolicy::ConstantWithinValidity },
            observations: vec![], provenance: Provenance { source: "synthetic reduced-viscosity tensor".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(key,dims,ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(),
        phase: "solid".into(), process: "constitutive test".into(), revision: 0 }, claims, vec![]).unwrap();
    resolve_material_state_point(&card,&QueryPoint::new().with("T",293.15).unwrap()
        .with("omega",1000.0).unwrap(),&requirements,MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn law(values: [f64; 6]) -> PlaneStressViscosity {
    PlaneStressViscosity::from_material(source(values,None,false,false), MaterialPlane::Axes12, KEYS).unwrap()
}
fn close(a:f64,b:f64) { assert!((a-b).abs() <= 3e-12*(a.abs()+b.abs()).max(1e-25),"{a} != {b}"); }

#[test]
fn source_tensor_drives_consistent_stress_power_and_engineering_shear() {
    let material = source([4.0,9.0,16.0,1.0,0.5,-0.75],None,false,false);
    let law = PlaneStressViscosity::from_material(material.clone(),MaterialPlane::Axes12,KEYS).unwrap();
    assert_eq!(law.source(),&material); assert_eq!(law.plane(),MaterialPlane::Axes12);
    assert_eq!(law.properties(),&KEYS.map(str::to_string));
    assert_eq!(law.matrix_pa_s(),&[4.0,1.0,0.5,1.0,9.0,-0.75,0.5,-0.75,16.0]);
    let (stress,power) = law.stress_and_power([2.0,-1.0,3.0]).unwrap();
    assert_eq!(stress,[8.5,-9.25,49.75]); assert_eq!(power,175.5);
    // Three independent finite differences of stress recover the original V.
    for j in 0..3 {
        let mut plus = [2.0,-1.0,3.0]; let mut minus = plus;
        plus[j] += 1e-3; minus[j] -= 1e-3;
        let (a,_) = law.stress_and_power(plus).unwrap(); let (b,_) = law.stress_and_power(minus).unwrap();
        for i in 0..3 { assert!(((a[i]-b[i])/2e-3-law.matrix_pa_s()[3*i+j]).abs() < 2e-11); }
    }
    let shear = self::law([4.0,9.0,16.0,0.0,0.0,0.0]);
    assert_eq!(shear.stress_and_power([0.0,0.0,2.0]).unwrap(),([0.0,0.0,32.0],64.0));
}

#[test]
fn arbitrary_orientation_preserves_tensor_work_and_thickness_cubed() {
    let law = law([4.0,9.0,16.0,1.0,0.5,-0.75]);
    assert_eq!(law.rotated_matrix_pa_s(0.0).unwrap(),*law.matrix_pa_s());
    for i in 0..64 {
        let angle = i as f64*0.137; let (s,c) = angle.sin_cos();
        let e = [0.3,-0.7,1.1];
        let local = [c*c*e[0]+s*s*e[1]+c*s*e[2],
            s*s*e[0]+c*c*e[1]-c*s*e[2],-2.0*c*s*e[0]+2.0*c*s*e[1]+(c*c-s*s)*e[2]];
        let (_, expected) = law.stress_and_power(local).unwrap();
        let v = law.rotated_matrix_pa_s(angle).unwrap();
        assert_eq!(v[1],v[3]); assert_eq!(v[2],v[6]); assert_eq!(v[5],v[7]);
        let actual = (0..3).map(|r| (0..3).map(|k|e[r]*v[3*r+k]*e[k]).sum::<f64>()).sum();
        close(actual,expected);
        let a = law.bending_matrix_n_m_s(0.002,angle).unwrap();
        let b = law.bending_matrix_n_m_s(0.004,angle).unwrap();
        for j in 0..9 { close(b[j],8.0*a[j]); close(a[j],v[j]*(0.002_f64.powi(3)/12.0)); }
    }
}

#[test]
fn absent_indefinite_unresolved_or_frequency_sampled_data_never_becomes_a_default() {
    let good = [4.0,9.0,16.0,1.0,0.5,-0.75];
    for values in [[-1.0,9.0,16.0,0.0,0.0,0.0],[1.0,1.0,1.0,2.0,0.0,0.0],
        [1.0,1.0,1.0,0.0,0.0,2.0],[1.0,1.0,0.0,0.0,0.5,0.0]] {
        assert!(PlaneStressViscosity::from_material(source(values,None,false,false),MaterialPlane::Axes12,KEYS).is_err());
    }
    for material in [source(good,Some(5),false,false),source(good,None,true,false),source(good,None,false,true)] {
        assert!(PlaneStressViscosity::from_material(material,MaterialPlane::Axes12,KEYS).is_err());
    }
    let mut duplicate = KEYS; duplicate[5] = duplicate[0];
    assert!(PlaneStressViscosity::from_material(source(good,None,false,false),MaterialPlane::Axes12,duplicate).is_err());
    let zero = law([0.0;6]); assert_eq!(zero.stress_and_power([2.0,-3.0,4.0]).unwrap(),([0.0;3],0.0));
    for omega in [0.0,10001.0,f64::NAN] { assert!(zero.validate_frequency(omega).is_err()); }
    zero.validate_frequency(10.0).unwrap(); zero.validate_frequency(10000.0).unwrap();
    let law = law(good);
    assert!(law.rotated_matrix_pa_s(f64::NAN).is_err());
    assert!(law.bending_matrix_n_m_s(1e-200,0.0).is_err());
    assert!(law.stress_and_power([f64::MAX;3]).is_err());
}

#[test]
fn directional_source_retains_nullspaces_cross_terms_and_rotated_work() {
    let axial=law([4.0,0.0,0.0,0.0,0.0,0.0]);
    assert_eq!(axial.stress_and_power([0.0,3.0,-7.0]).unwrap(),([0.0;3],0.0));
    let coupled=law([4.0,9.0,1.0,6.0,-2.0,-3.0]); // [2,3,-1]^T [2,3,-1]
    for rate in [[1.0,2.0,-4.0],[-3.0,2.0,0.0],[0.0,1.0,3.0]] {
        let extension=2.0*rate[0]+3.0*rate[1]-rate[2];
        let (stress,power)=coupled.stress_and_power(rate).unwrap();
        assert!((power-extension*extension).abs()<2e-12);
        for (actual,row) in stress.into_iter().zip([2.0,3.0,-1.0]) {
            assert!((actual-row*extension).abs()<2e-12);
        }
    }
    for angle in [0.17_f64,0.43,1.11] {
        let (s,c)=angle.sin_cos(); let rate=[0.3,-0.7,1.1];
        let local=[c*c*rate[0]+s*s*rate[1]+c*s*rate[2],
            s*s*rate[0]+c*c*rate[1]-c*s*rate[2],
            -2.0*c*s*rate[0]+2.0*c*s*rate[1]+(c*c-s*s)*rate[2]];
        let expected=(2.0*local[0]+3.0*local[1]-local[2]).powi(2);
        let matrix=coupled.rotated_matrix_pa_s(angle).unwrap();
        let actual:f64=(0..3).flat_map(|i|(0..3).map(move |j|
            rate[i]*matrix[3*i+j]*rate[j])).sum();
        close(actual,expected);
        let bending=coupled.bending_matrix_n_m_s(0.003,angle).unwrap();
        for i in 0..9 {close(bending[i],matrix[i]*0.003_f64.powi(3)/12.0);}
    }
}

#[test]
fn passive_resultant_is_scale_independent_and_rejects_active_cross_coupling() {
    use fs_material::plane_stress_viscosity::factor_passive_resultant;
    for scale in [1e-200,1.0,1e200] {
        let v=[4.0,6.0,-2.0,6.0,9.0,-3.0,-2.0,-3.0,1.0].map(|x|x*scale);
        let factor=factor_passive_resultant(&v).unwrap();
        for i in 0..3 {for j in 0..3 {
            let actual:f64=(0..6).map(|k|factor[i][k]*factor[j][k]/scale).sum();
            assert!((actual-v[3*i+j]/scale).abs()<2e-12);
        }}
        let bad=[0.0,scale,0.0,scale,0.0,0.0,0.0,0.0,0.0];
        assert!(factor_passive_resultant(&bad).is_err());
        let active=[1.0,1.01,0.0,1.01,1.0,0.0,0.0,0.0,0.0].map(|x|x*scale);
        assert!(factor_passive_resultant(&active).is_err());
    }
}
