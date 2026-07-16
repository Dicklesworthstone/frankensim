//! fs-matdb PR-3 conformance: interface systems are ORDERED
//! system+history identities, never unordered pair constants.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MatDbError, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_qty::Dims;

fn provenance() -> Provenance {
    Provenance {
        source: "pin-on-disk campaign POD-7".to_string(),
        license: "internal-use".to_string(),
        artifact: None,
    }
}

fn steel() -> SurfaceSpec {
    SurfaceSpec {
        material: MaterialStateId {
            chemistry: "AISI52100".to_string(),
            phase: "bearing-steel".to_string(),
            process: "hardened-60HRC".to_string(),
            revision: 0,
        },
        texture_frame: "ground-Ra0.2-frame-17".to_string(),
    }
}

fn ptfe() -> SurfaceSpec {
    SurfaceSpec {
        material: MaterialStateId {
            chemistry: "PTFE".to_string(),
            phase: "sintered".to_string(),
            process: "as-molded".to_string(),
            revision: 0,
        },
        texture_frame: "molded-frame-3".to_string(),
    }
}

fn friction_claims() -> ClaimSet {
    let mut set = ClaimSet::new();
    set.insert_claim(PropertyClaim {
        key: PropertyKey::new("kinetic-friction-coefficient", Dims::NONE),
        value: PropertyValue::Scalar {
            value: 0.08,
            dims: Dims::NONE,
        },
        validity: ValidityDomain::unconstrained()
            .with("T", 273.15, 353.15)
            .with("normal_pressure", 1.0e5, 5.0e6),
        uncertainty: UncertaintyModel::HalfWidth {
            half_width: 0.02,
            confidence: 0.9,
        },
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: provenance(),
    })
    .expect("friction claim inserts");
    set
}

fn dry_air(history: &str) -> SystemContext {
    SystemContext {
        medium: "dry".to_string(),
        third_body: None,
        environment: "air".to_string(),
        history: history.to_string(),
    }
}

fn dry_system(a: SurfaceSpec, b: SurfaceSpec, history: &str) -> InterfaceSystemCard {
    InterfaceSystemCard::assemble(a, b, dry_air(history), friction_claims(), Vec::new())
        .expect("system assembles")
}

#[test]
fn interface_systems_are_ordered_and_history_bearing() {
    let steel_on_ptfe = dry_system(steel(), ptfe(), "run-in-1000-cycles");
    let ptfe_on_steel = dry_system(ptfe(), steel(), "run-in-1000-cycles");
    assert_ne!(
        steel_on_ptfe.content_hash(),
        ptfe_on_steel.content_hash(),
        "surface order is semantic: (a,b) and (b,a) are different systems"
    );

    let virgin = dry_system(steel(), ptfe(), "virgin");
    assert_ne!(
        steel_on_ptfe.content_hash(),
        virgin.content_hash(),
        "history is identity-bearing, not a footnote"
    );

    assert_eq!(steel_on_ptfe.medium(), "dry");
    assert_eq!(steel_on_ptfe.environment(), "air");
    assert_eq!(steel_on_ptfe.history(), "run-in-1000-cycles");
    assert_eq!(
        steel_on_ptfe
            .claims_for("kinetic-friction-coefficient")
            .len(),
        1
    );
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"interface-order-history\",\"verdict\":\"pass\",\
         \"detail\":\"swapping surfaces or history moves the system identity\"}}"
    );
}

#[test]
fn wetting_is_a_three_phase_system_with_hysteresis_claims_coexisting() {
    let mut claims = ClaimSet::new();
    for (name, angle) in [
        ("advancing-contact-angle", 115.0_f64.to_radians()),
        ("receding-contact-angle", 95.0_f64.to_radians()),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, Dims::NONE),
                value: PropertyValue::Scalar {
                    value: angle,
                    dims: Dims::NONE,
                },
                validity: ValidityDomain::unconstrained().with("T", 288.15, 308.15),
                uncertainty: UncertaintyModel::HalfWidth {
                    half_width: 2.0_f64.to_radians(),
                    confidence: 0.9,
                },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: Vec::new(),
                provenance: provenance(),
            })
            .expect("angle claim inserts");
    }
    let sessile = SurfaceSpec {
        material: MaterialStateId {
            chemistry: "PTFE".to_string(),
            phase: "sintered".to_string(),
            process: "as-molded".to_string(),
            revision: 0,
        },
        texture_frame: "polished-frame-9".to_string(),
    };
    // Wetting: the liquid IS the medium, the gas IS the environment;
    // the second "surface" is the liquid's free surface, named by its
    // material state at revision 0.
    let liquid = SurfaceSpec {
        material: MaterialStateId {
            chemistry: "H2O-deionized".to_string(),
            phase: "liquid".to_string(),
            process: "as-supplied".to_string(),
            revision: 0,
        },
        texture_frame: "free-surface".to_string(),
    };
    let system = InterfaceSystemCard::assemble(
        sessile,
        liquid,
        SystemContext {
            medium: "water-deionized".to_string(),
            third_body: None,
            environment: "air-50pct-RH".to_string(),
            history: "virgin".to_string(),
        },
        claims,
        Vec::new(),
    )
    .expect("wetting system assembles");
    assert_eq!(system.claims_for("advancing-contact-angle").len(), 1);
    assert_eq!(system.claims_for("receding-contact-angle").len(), 1);
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"wetting-three-phase\",\"verdict\":\"pass\",\
         \"detail\":\"solid-liquid-gas system holds advancing/receding hysteresis as separate claims\"}}"
    );
}

#[test]
fn incomplete_system_identities_refuse() {
    let mut unnamed = steel();
    unnamed.texture_frame = "  ".to_string();
    assert!(matches!(
        InterfaceSystemCard::assemble(
            unnamed,
            ptfe(),
            dry_air("virgin"),
            ClaimSet::new(),
            Vec::new()
        ),
        Err(MatDbError::MissingTextureFrame { .. })
    ));

    let mut blank_medium = dry_air("virgin");
    blank_medium.medium = " ".to_string();
    assert!(matches!(
        InterfaceSystemCard::assemble(steel(), ptfe(), blank_medium, ClaimSet::new(), Vec::new()),
        Err(MatDbError::MissingSystemField { field: "medium" })
    ));

    let mut blank_history = dry_air("virgin");
    blank_history.history = String::new();
    assert!(matches!(
        InterfaceSystemCard::assemble(steel(), ptfe(), blank_history, ClaimSet::new(), Vec::new()),
        Err(MatDbError::MissingSystemField { field: "history" })
    ));

    let mut blank_third_body = dry_air("virgin");
    blank_third_body.medium = "oil-SAE30".to_string();
    blank_third_body.third_body = Some(String::new());
    assert!(matches!(
        InterfaceSystemCard::assemble(
            steel(),
            ptfe(),
            blank_third_body,
            ClaimSet::new(),
            Vec::new()
        ),
        Err(MatDbError::MissingSystemField {
            field: "third_body"
        })
    ));
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"interface-gates\",\"verdict\":\"pass\",\
         \"detail\":\"unnamed texture frame and blank medium/history/third-body refuse typed\"}}"
    );
}
