use super::*;
use crate::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::DesignControl;

const MODEL: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/examples/equilibrium-design.model"));
const FIT: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/examples/equilibrium-response-limits.fit"));

#[test]
fn optional_response_constraints_preserve_old_objective_and_exact_source_identity() {
    let gate = CancelGate::new();
    let legacy = FIT.split_once("constraint_limits").unwrap().0;
    let old = EquilibriumDesignFile::from_bytes(MODEL, legacy.as_bytes(), &gate).unwrap();
    let new = EquilibriumDesignFile::from_bytes(MODEL, FIT.as_bytes(), &gate).unwrap();
    let a = old.problem().evaluate(&[0.0], &mut DesignControl::new(1,1), &gate).unwrap();
    let b = new.problem().evaluate(&[0.0], &mut DesignControl::new(1,1), &gate).unwrap();
    assert_eq!(a.value,b.value); assert_eq!(a.gradient,b.gradient); assert_eq!(a.cases,b.cases);
    assert!(old.problem().constraints().is_empty()); assert_eq!(new.problem().constraints().len(),4);
    assert_eq!(new.design_hash(),hash_domain(EQUILIBRIUM_DESIGN_HASH_DOMAIN,FIT.as_bytes()));
    assert_eq!(old.design_hash(),hash_domain(EQUILIBRIUM_DESIGN_HASH_DOMAIN,legacy.as_bytes()));
    assert_ne!(new.design_hash(),old.design_hash()); assert_eq!(new.model_info(),old.model_info());
    assert!(b.constraints[0].value>0.8 && b.constraints[0].residual>0.0);
    assert!(b.constraints[2].value<0.0 && b.constraints[2].residual<0.0, "retain signed spring reaction and lower-bound slack");
    assert!(matches!(&new.problem().constraints()[3].quantity,ResponseQuantity::Displacement(_)));
}

#[test]
fn incomplete_unknown_or_oversized_constraint_sections_cannot_be_ignored() {
    let gate=CancelGate::new();
    for text in [
        FIT.replace("constraint_limits 4","constraint_limits 3"),
        FIT.replace("constraint_limits 4","constraint_limits 65"),
        FIT.replace("constraints 4","constraints 3"),
        FIT.replace("contact-force 0","contact-force 99"),
        FIT.replace("normal-cap 0","normal-cap 99"),
        FIT.replace("normal-cap","travel-cap"),
        FIT.replace("displacement 1 0","displacement 1 99"),
        FIT.replace("at-most 0.8 1","at-most NaN 1"),
        FIT.replace("at-most 0.8 1","at-most 0.8 0"),
        FIT.replace("at-most 0.8 1","at-most 0.8 -1"),
        FIT.replace("at-most 0.8 1","approximately 0.8 1"),
        FIT.replace("contact-force 0","unknown-response 0"),
        format!("{FIT}ignored\n"),
    ] { assert!(EquilibriumDesignFile::from_bytes(MODEL,text.as_bytes(),&gate).is_err(),"{text}"); }
    let start=FIT.find("constraint_limits").unwrap();
    for (end,_) in FIT.match_indices('\n').filter(|(i,_)|*i>=start && *i+1<FIT.len()) {
        assert!(EquilibriumDesignFile::from_bytes(MODEL,&FIT.as_bytes()[..end+1],&gate).is_err());
    }
}
