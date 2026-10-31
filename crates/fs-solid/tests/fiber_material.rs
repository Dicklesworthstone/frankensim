//! G0/G1/G3 for the native material-driven RC fiber path.
//! Source fixtures are authored numerical data, not measured RC calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fiber::{ManderConcrete, MenegottoPintoSteel};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Density, Dims, Pressure};
use fs_solid::fiber::material::{
    ConcreteFiberProperties, FiberMaterial, RcSection, RcSectionGeometry, SteelFiberProperties,
};

fn properties() -> ConcreteFiberProperties {
    ConcreteFiberProperties {
        density: "density".into(),
        compressive_strength: "fcc".into(),
        peak_strain: "eps_cc".into(),
        initial_modulus: "ec".into(),
        ultimate_strain: "eps_cu".into(),
    }
}
fn steel_properties() -> SteelFiberProperties {
    SteelFiberProperties {
        density: "density".into(),
        young_modulus: "e0".into(),
        yield_stress: "fy".into(),
        hardening_ratio: "b".into(),
        initial_curvature: "r0".into(),
        curvature_degradation: "a1".into(),
        curvature_scale: "a2".into(),
    }
}
fn card(values: &[(&str, f64, Dims)], strain_dependent: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for &(name, value, dims) in values {
        let (property, interpolation, validity) = if strain_dependent && name == "ec" {
            (
                PropertyValue::Curve {
                    abscissa: "strain".into(),
                    abscissa_dims: Dims::NONE,
                    knots: vec![(0.0, value), (0.01, 0.5 * value)],
                    dims,
                },
                InterpolationPolicy::LinearInside,
                ValidityDomain::unconstrained().with("strain", 0.0, 0.01),
            )
        } else {
            (
                PropertyValue::Scalar { value, dims },
                InterpolationPolicy::ConstantWithinValidity,
                ValidityDomain::unconstrained(),
            )
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: property,
                validity,
                uncertainty: UncertaintyModel::RelativeHalfWidth {
                    fraction: 0.01,
                    confidence: 0.95,
                },
                interpolation,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic RC coefficient fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "numerical RC constituent".into(),
            phase: "solid".into(),
            process: "authored constitutive fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn concrete_card(density: f64, modulus: f64, strain_dependent: bool) -> MaterialCard {
    card(
        &[
            ("density", density, Density::DIMS),
            ("fcc", 30e6, Pressure::DIMS),
            ("eps_cc", 0.002, Dims::NONE),
            ("ec", modulus, Pressure::DIMS),
            ("eps_cu", 0.005, Dims::NONE),
        ],
        strain_dependent,
    )
}
fn resolved_concrete(card: &MaterialCard) -> FiberMaterial {
    FiberMaterial::resolve_concrete(
        card,
        &QueryPoint::new(),
        &properties(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn steel() -> FiberMaterial {
    FiberMaterial::steel(
        MenegottoPintoSteel::new(200e9, 450e6, 0.01).unwrap(),
        7850.0,
    )
    .unwrap()
}
fn geometry(layers: usize) -> RcSectionGeometry {
    RcSectionGeometry {
        depth_m: 0.5,
        width_m: 0.4,
        cover_m: 0.04,
        concrete_layers: layers,
        steel_area_per_face_m2: 0.0015,
        steel_centroid_y_m: 0.18,
    }
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn material_areas_replace_concrete_and_mass_is_refinement_invariant() {
    let core = FiberMaterial::concrete(
        ManderConcrete::new(42e6, 0.004, 32e9, 0.015).unwrap(),
        2500.0,
    )
    .unwrap();
    let cover = FiberMaterial::concrete(
        ManderConcrete::new(30e6, 0.002, 30e9, 0.005).unwrap(),
        2300.0,
    )
    .unwrap();
    let g = geometry(16);
    let gross = g.depth_m * g.width_m;
    let core_area = (g.depth_m - 2.0 * g.cover_m) * (g.width_m - 2.0 * g.cover_m);
    let reinforcement = 2.0 * g.steel_area_per_face_m2;
    let expected_mass = 2500.0 * (core_area - reinforcement)
        + 2300.0 * (gross - core_area)
        + 7850.0 * reinforcement;
    for layers in [1, 7, 16, 25, 63] {
        let section =
            RcSection::new(geometry(layers), core.clone(), cover.clone(), steel()).unwrap();
        close(section.line_mass_kg_m(), expected_mass);
        close(section.section().fibers.iter().map(|f| f.area).sum(), gross);
        assert!(section.section().fibers.iter().all(|f| f.area > 0.0));
    }
}

#[test]
fn one_source_change_reaches_stiffness_mass_and_retained_receipts() {
    let first = resolved_concrete(&concrete_card(2300.0, 25e9, false));
    let second = resolved_concrete(&concrete_card(2800.0, 45e9, false));
    assert_eq!(first.source().unwrap().properties().len(), 5);
    assert_ne!(
        first.source().unwrap().identity(),
        second.source().unwrap().identity()
    );
    let a = RcSection::new(geometry(16), first.clone(), first, steel()).unwrap();
    let b = RcSection::new(geometry(16), second.clone(), second, steel()).unwrap();
    let ra = a.section().respond(-1e-6, 0.0);
    let rb = b.section().respond(-1e-6, 0.0);
    assert!(rb.tangent[0][0] > 1.5 * ra.tangent[0][0]);
    assert!(rb.n.abs() > 1.5 * ra.n.abs());
    let concrete_area = 0.5 * 0.4 - 2.0 * 0.0015;
    close(
        b.line_mass_kg_m() - a.line_mass_kg_m(),
        (2800.0 - 2300.0) * concrete_area,
    );
    let scaled = b.scaled(1.7).unwrap();
    let rs = scaled.section().respond(-1e-6, 0.0);
    close(rs.n, 1.7 * rb.n);
    close(rs.tangent[0][0], 1.7 * rb.tangent[0][0]);
    close(scaled.line_mass_kg_m(), 1.7 * b.line_mass_kg_m());
    assert_eq!(scaled.geometry(), b.geometry());
    assert_eq!(
        scaled.materials()[0].source().unwrap().identity(),
        b.materials()[0].source().unwrap().identity()
    );
}

#[test]
fn steel_source_requires_every_cyclic_parameter() {
    let source = card(
        &[
            ("density", 7850.0, Density::DIMS),
            ("e0", 190e9, Pressure::DIMS),
            ("fy", 400e6, Pressure::DIMS),
            ("b", 0.02, Dims::NONE),
            ("r0", 18.0, Dims::NONE),
            ("a1", 15.0, Dims::NONE),
            ("a2", 0.2, Dims::NONE),
        ],
        false,
    );
    let material = FiberMaterial::resolve_steel(
        &source,
        &QueryPoint::new(),
        &steel_properties(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(material.source().unwrap().properties().len(), 7);
    match material.law() {
        fs_solid::fiber::FiberLaw::Steel(law) => {
            assert_eq!(law.e0, 190e9);
            assert_eq!(law.b, 0.02);
            assert_eq!(law.r0, 18.0);
            assert_eq!(law.a1, 15.0);
            assert_eq!(law.a2, 0.2);
        }
        _ => panic!("steel source did not produce the steel law"),
    }
    let mut missing = steel_properties();
    missing.curvature_scale = "missing-a2".into();
    assert!(
        FiberMaterial::resolve_steel(
            &source,
            &QueryPoint::new(),
            &missing,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .is_err()
    );
}

#[test]
fn invalid_geometry_and_unowned_state_dependence_refuse() {
    let concrete = resolved_concrete(&concrete_card(2400.0, 30e9, false));
    for g in [
        RcSectionGeometry {
            cover_m: 0.25,
            ..geometry(16)
        },
        RcSectionGeometry {
            steel_centroid_y_m: 0.001,
            ..geometry(16)
        },
        RcSectionGeometry {
            steel_area_per_face_m2: 0.1,
            ..geometry(16)
        },
        RcSectionGeometry {
            steel_area_per_face_m2: f64::MIN_POSITIVE,
            ..geometry(16)
        },
        RcSectionGeometry {
            concrete_layers: 0,
            ..geometry(16)
        },
        RcSectionGeometry {
            depth_m: f64::INFINITY,
            ..geometry(16)
        },
    ] {
        assert!(RcSection::new(g, concrete.clone(), concrete.clone(), steel()).is_err());
    }
    let source = concrete_card(2400.0, 30e9, true);
    let point = QueryPoint::new().with("strain", 0.0).unwrap();
    assert!(
        FiberMaterial::resolve_concrete(
            &source,
            &point,
            &properties(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .is_err()
    );
    let mut invalid_steel = MenegottoPintoSteel::new(200e9, 450e6, 0.01).unwrap();
    invalid_steel.a1 = invalid_steel.r0;
    assert!(FiberMaterial::steel(invalid_steel, 7850.0).is_err());
    let section = RcSection::new(geometry(16), concrete.clone(), concrete, steel()).unwrap();
    assert!(section.scaled(f64::INFINITY).is_err());
    assert!(section.scaled(0.0).is_err());
}
