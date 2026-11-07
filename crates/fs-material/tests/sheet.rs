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
        process:"characterized construction".into(),revision:0},claims,vec![])
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
    // An open material key may retain bulk Pa; the membrane owner still
    // requires effective sheet stiffness in N/m and must not reinterpret it.
    if try_card(&wrong).is_ok() { assert!(adopt(&wrong).membrane().is_err()); }
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

// G3: one original source law reaches sheet inertia, membrane and independent
// bending at different held environments. These are independently initialized
// specimens; no live moisture transfer or pressure work is inferred.
fn environmental_sheet(
    axis: &str,
    quantity: Option<fs_qty::QuantitySpec>,
    bounds: (f64, f64),
    curve: bool,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let factors = [2.0, 1.25, 0.5, 0.5, 1.0, 0.5, 3.0, 3.0, 3.0, 3.0, 1.0];
    for ((name, value, dims), factor) in entries().into_iter().zip(factors) {
        let validity = ValidityDomain::unconstrained().with("T", 290.0, 310.0);
        let validity = match quantity {
            Some(quantity) => validity.with_quantity(axis, quantity, bounds.0, bounds.1),
            None => validity.with(axis, bounds.0, bounds.1),
        };
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: if curve {
                PropertyValue::Curve {
                    abscissa: axis.into(),
                    abscissa_dims: quantity.map_or(Dims::NONE, |q| q.dims()),
                    knots: vec![(bounds.0, value), (bounds.1, factor * value)],
                    dims,
                }
            } else {
                PropertyValue::Scalar { value, dims }
            },
            validity,
            interpolation: if curve {
                InterpolationPolicy::LinearInside
            } else {
                InterpolationPolicy::ConstantWithinValidity
            },
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured held-environment sheet".into(),
                license: "test-only".into(),
                artifact: None,
            },
        }).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic sheet".into(),
            phase: "solid".into(),
            process: "held environment".into(),
            revision: 0,
        },
        claims,
        vec![],
    ).unwrap()
}

fn environmental_query(
    axis: &str,
    quantity: Option<fs_qty::QuantitySpec>,
    value: f64,
) -> QueryPoint {
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    match quantity {
        Some(quantity) => point.with_quantity(axis, quantity, value).unwrap(),
        None => point.with(axis, value).unwrap(),
    }
}

#[test]
fn g3_held_environment_reaches_independent_sheet_channels() {
    use fs_qty::{QuantitySpec, semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm}};
    for (axis, kind, bounds) in [
        ("moisture_content", QuantityKind::Moisture(MoistureBasis::DryMass), (0.0, 0.2)),
        ("ambient_pressure", QuantityKind::Pressure, (100_000.0, 200_000.0)),
    ] {
        let quantity = Some(QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static)));
        let card = environmental_sheet(axis, quantity, bounds, true);
        let caps = SheetCapabilities { membrane: true, bending: true };
        let low_query = environmental_query(axis, quantity, bounds.0);
        let high_query = environmental_query(axis, quantity, bounds.1);
        let low = resolve_sheet(&card, &low_query, caps, Pick::SingleClaimOnly).unwrap();
        let high = resolve_sheet(&card, &high_query, caps, Pick::SingleClaimOnly).unwrap();
        assert_eq!(low.material().card_identity(), card.content_hash());
        assert_eq!(high.material().card_identity(), card.content_hash());
        assert_eq!(low.material().query(), low_query);
        assert_eq!(high.material().query(), high_query);
        assert_ne!(low.material().identity(), high.material().identity());
        close(high.areal_density_kg_m2(), 2.0 * low.areal_density_kg_m2());
        close(high.thickness_m(), 1.25 * low.thickness_m());

        let strain = [0.02, -0.01, 0.03];
        let a = low.membrane().unwrap().response(strain).unwrap();
        let b = high.membrane().unwrap().response(strain).unwrap();
        close(b.energy_j_m2, 0.5 * a.energy_j_m2);
        for i in 0..3 {
            close(b.resultants_n_m[i], 0.5 * a.resultants_n_m[i]);
        }
        let a = low.bending().unwrap();
        let b = high.bending().unwrap();
        for i in [0, 1, 3, 4, 8] {
            close(b.coefficients()[i], 3.0 * a.coefficients()[i]);
        }

        let mid_query = environmental_query(axis, quantity, 0.5 * (bounds.0 + bounds.1));
        let requirements = entries().into_iter().map(|(name, _, dims)|
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap()
        ).collect::<Vec<_>>();
        let bundle = resolve_material_state_point(
            &card, &mid_query, &requirements, Pick::SingleClaimOnly,
        ).unwrap();
        let original = bundle.clone();
        let mid = ResolvedSheet::try_from_material(bundle).unwrap();
        assert_eq!(mid.material(), &original);
        assert!(mid.material().property("specific_heat_capacity").is_some());
        close(mid.areal_density_kg_m2(), 1.5 * low.areal_density_kg_m2());
        let property = mid.material().property(AREAL_DENSITY).unwrap();
        assert_eq!(property.source_axes(), &[(axis.into(), quantity.unwrap().dims())]);
        assert!(matches!(property.answer().receipt.decision,
            fs_matdb::EvaluationDecision::LinearInside { .. }));
        assert!(resolve_sheet(
            &card, &environmental_query(axis, quantity, bounds.1 * 1.1),
            caps, Pick::SingleClaimOnly,
        ).is_err());
    }
}

#[test]
fn g0_sheet_never_freezes_unowned_response_axes() {
    use fs_qty::{QuantitySpec, semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm}};
    for (axis, quantity, bounds) in [
        ("moisture_content", None, (0.0, 0.2)),
        ("normal_stress", Some(QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Stress, ValueForm::Static))), (1.0, 2.0)),
        ("omega", Some(QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Frequency(FrequencyConvention::Angular), ValueForm::Static))), (0.0, 20_000.0)),
    ] {
        let card = environmental_sheet(axis, quantity, bounds, true);
        let query = environmental_query(axis, quantity, 0.5 * (bounds.0 + bounds.1));
        // The original data query succeeds. The constant sheet model cannot
        // freeze these response axes into reference inventory or stiffness.
        let requirements = entries().into_iter().map(|(name, _, dims)|
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap()
        ).collect::<Vec<_>>();
        let bundle = resolve_material_state_point(
            &card, &query, &requirements, Pick::SingleClaimOnly,
        ).unwrap();
        assert!(ResolvedSheet::try_from_material(bundle).is_err());
    }
}

#[test]
fn legacy_characterization_and_spectral_applicability_stay_source_bound() {
    use fs_qty::{QuantitySpec, semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm}};
    let angular = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Frequency(FrequencyConvention::Angular), ValueForm::Static));
    for (axis, quantity, bounds) in [
        ("specimen", None, (3.0, 3.0)),
        ("ftx_characterization", None, (0.0, 0.0)),
        ("omega", Some(angular), (0.0, 1000.0)),
    ] {
        let card = environmental_sheet(axis, quantity, bounds, false);
        let query = environmental_query(axis, quantity, bounds.0);
        let sheet = resolve_sheet(&card, &query,
            SheetCapabilities { membrane: true, bending: true }, Pick::SingleClaimOnly,
        ).unwrap();
        assert_eq!(sheet.material().card_identity(), card.content_hash());
        assert_eq!(sheet.material().query(), query);
        close(sheet.areal_density_kg_m2(), 0.2);
        close(sheet.thickness_m(), 0.001);
        close(sheet.bending().unwrap().coefficients()[0], 3e-5);
        for property in sheet.material().properties() {
            assert_eq!(property.answer().evidence.model.validity.bound(axis), Some(bounds));
            assert_eq!(property.angular_frequency_band_rad_s().unwrap(),
                if axis == "omega" { Some(bounds) } else { None });
        }
        // Generic adoption does not erase source limits. The eventual
        // constitutive owner remains responsible for its spectral contract.
        assert!(resolve_sheet(&card, &environmental_query(axis, quantity, bounds.1 + 1.0),
            SheetCapabilities { membrane: true, bending: true }, Pick::SingleClaimOnly,
        ).is_err());
    }
}

#[path="sheet/thermal.rs"]
mod thermal;
