//! Effective sheet properties remain independent and source-bound.
use fs_evidence::ValidityDomain;
use fs_material::sheet::*;
use fs_material::state_point::{MaterialPropertySelection as Pick, ScalarPropertyRequirement,
    ScalarAdmissibility, resolve_material_state_point};
use fs_matdb::{ClaimSet, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, InterpolationPolicy, Provenance, QueryPoint, UncertaintyModel};
use fs_qty::Dims;

fn entries() -> Vec<(&'static str, f64, Dims)> {
    vec![
        (AREAL_DENSITY, 0.2, Dims([-2,1,0,0,0,0])), (THICKNESS, 0.001, Dims([1,0,0,0,0,0])),
        ("warp_stiffness", 8000.0, MEMBRANE_DIMS), ("weft_stiffness", 2000.0, MEMBRANE_DIMS),
        ("poisson_warp_weft", 0.2, Dims::NONE), ("shear_stiffness", 400.0, MEMBRANE_DIMS),
        ("bending_warp", 3e-5, BENDING_DIMS), ("bending_weft", 8e-5, BENDING_DIMS),
        ("bending_coupling", -1e-5, BENDING_DIMS), ("bending_twist", 2e-5, BENDING_DIMS),
        ("specific_heat_capacity", 1400.0, Dims([2,0,-2,-1,0,0])),
    ]
}
fn try_card(rows: &[(&str,f64,Dims)]) -> Result<MaterialCard, fs_matdb::MatDbError> {
    let mut claims = ClaimSet::new();
    for &(name,value,dims) in rows {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name,dims), value: PropertyValue::Scalar {value,dims},
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            validity: ValidityDomain::unconstrained().with("T",290.0,310.0),
            uncertainty: UncertaintyModel::Unstated, observations:vec![],
            provenance:Provenance {source:"synthetic effective sheet, not a measured fabric".into(),
                license:"test-only".into(),artifact:None},
        })?;
    }
    MaterialCard::assemble(MaterialStateId {chemistry:"synthetic sheet".into(),phase:"solid".into(),
        process:"characterized construction".into(),revision:1},claims,vec![])
}
fn card(rows: &[(&str,f64,Dims)]) -> MaterialCard { try_card(rows).unwrap() }
fn adopt(rows: &[(&str,f64,Dims)]) -> ResolvedSheet {
    let requirements:Vec<_> = rows.iter().map(|&(key,_,dims)|
        ScalarPropertyRequirement::try_new(key,dims,ScalarAdmissibility::Finite).unwrap()).collect();
    let bundle = resolve_material_state_point(&card(rows), &QueryPoint::new().with("T",300.0).unwrap(),
        &requirements, Pick::SingleClaimOnly).unwrap();
    let copy = bundle.clone();
    let sheet = ResolvedSheet::try_from_material(bundle).unwrap();
    assert_eq!(sheet.material(), &copy);
    sheet
}
fn close(a:f64,b:f64) { assert!((a-b).abs() < 1e-8*a.abs().max(b.abs()).max(1e-8),"{a} != {b}"); }

#[test]
fn source_substitution_changes_only_independent_channels() {
    let rows=entries();let base=adopt(&rows);
    assert!(base.material().property("specific_heat_capacity").is_some());
    let mut altered=rows.clone();altered[0].1 *= 2.0;altered[1].1 *= 3.0;
    let heavy=adopt(&altered);
    assert_eq!(base.membrane().unwrap(),heavy.membrane().unwrap());
    assert_eq!(base.bending().unwrap(),heavy.bending().unwrap());
    assert_eq!(heavy.areal_density_kg_m2(),2.0*base.areal_density_kg_m2());
    let mut altered=rows;altered[6].1 *= 4.0;
    let bend=adopt(&altered);
    assert_eq!(base.membrane().unwrap(),bend.membrane().unwrap());
    assert_ne!(base.bending().unwrap(),bend.bending().unwrap());
    assert_ne!(base.material().identity(),bend.material().identity());
}
#[test]
fn missing_sheet_data_never_becomes_bulk_plate_or_default_fabric() {
    let mut rows=entries(); rows.retain(|r|r.0!="bending_twist");
    let partial=adopt(&rows);
    assert!(partial.membrane().is_ok());assert!(partial.bending().is_err());
    let mut wrong=entries();wrong[2].2=Dims([-1,1,-2,0,0,0]);
    // The common source boundary already rejects bulk Pa under this sheet key.
    assert!(try_card(&wrong).is_err());
    let mut wrong=entries();wrong[6].2=Dims([3,1,-2,0,0,0]);
    // A consumer also refuses wrong units when a source key has no registered
    // dimension contract; either boundary must reject, never reinterpret.
    if try_card(&wrong).is_ok() { assert!(adopt(&wrong).bending().is_err()); }
    for invalid in [0.0,-1.0] {
        let mut wrong=entries();wrong[2].1=invalid;assert!(adopt(&wrong).membrane().is_err());
    }
}
#[test]
fn energy_stress_and_tangent_are_work_conjugate() {
    let law=adopt(&entries()).membrane().unwrap();let e=[0.021,-0.008,0.037];
    let r=law.response(e).unwrap();let delta=1e-6;
    for i in 0..3 {
        let mut lo=e;let mut hi=e;lo[i]-=delta;hi[i]+=delta;
        let (a,b)=(law.response(lo).unwrap(),law.response(hi).unwrap());
        close((b.energy_j_m2-a.energy_j_m2)/(2.0*delta),r.resultants_n_m[i]);
        for j in 0..3 {close((b.resultants_n_m[j]-a.resultants_n_m[j])/(2.0*delta),r.tangent_n_m[3*j+i]);}
    }
    close(r.resultants_n_m[2],400.0*e[2]); // no hidden factor two in engineering shear
    assert_eq!(law.response([0.0;3]).unwrap().energy_j_m2,0.0);
}
#[test]
fn grain_rotation_preserves_energy_and_exchanges_principal_axes() {
    let d=adopt(&entries()).bending().unwrap();let r=d.rotated(std::f64::consts::FRAC_PI_2).unwrap();
    close(r.coefficients()[0],d.coefficients()[4]);close(r.coefficients()[4],d.coefficients()[0]);
    let angle:f64=0.37;let(s,c)=angle.sin_cos();let e=[0.8,-0.3,0.5];
    let material=[c*c*e[0]+s*s*e[1]+c*s*e[2],s*s*e[0]+c*c*e[1]-c*s*e[2],
        -2.0*c*s*e[0]+2.0*c*s*e[1]+(c*c-s*s)*e[2]];
    let energy=|m:&[f64;9],x:[f64;3]|(0..3).flat_map(|i|(0..3).map(move |j|x[i]*m[3*i+j]*x[j])).sum::<f64>();
    close(energy(d.rotated(angle).unwrap().coefficients(),e),energy(d.coefficients(),material));
    assert_eq!(d,d.rotated(0.0).unwrap());
}
#[test]
fn query_domain_and_constitutive_stability_are_independent_gates() {
    let c=card(&entries());let caps=SheetCapabilities {membrane:true,bending:true};
    assert!(resolve_sheet(&c,&QueryPoint::new().with("T",300.0).unwrap(),caps,Pick::SingleClaimOnly).is_ok());
    assert!(resolve_sheet(&c,&QueryPoint::new().with("T",350.0).unwrap(),caps,Pick::SingleClaimOnly).is_err());
    assert!(SheetMatrix::orthotropic(1.0,1.0,1.0,1.0).is_err());
    assert!(SheetMatrix::orthotropic(1.0,1.0,-1.0,1.0).is_err());
    assert!(QuadraticMembrane::new(1.0,1.0,1.0,1.0).is_err());
    assert!(adopt(&entries()).membrane().unwrap().response([f64::NAN;3]).is_err());
}
