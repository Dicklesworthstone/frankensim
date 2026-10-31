//! Native multi-backstress fibers, immutable source bindings and structural use.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::IsotropicElastic;
use fs_material::creep::uniaxial::UniaxialNorton;
use fs_material::creep::{NortonCreep, NortonCreepParameters};
use fs_material::fiber::{ManderConcrete, Uniaxial};
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::chaboche::axial::UniaxialChaboche;
use fs_material::plastic::chaboche::source::BackstressProperties;
use fs_material::plastic::chaboche::{BackstressBranch, ChabocheJ2};
use fs_material::plastic::source::HardeningPropertyBasis;
use fs_material::plastic::voce::VoceJ2;
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY, VOCE_SATURATION_RATE_PROPERTY,
};
use fs_material::state_point::{MaterialPropertySelection, MaterialStatePointError};
use fs_qty::{Density, Dims, Pressure};
use fs_solid::fiber::chaboche::{ChabocheFiberLaw, ChabocheFiberState};
use fs_solid::fiber::material::{FiberMaterial, RcSection, RcSectionError, RcSectionGeometry};
use fs_solid::fiber::{Fiber, FiberLaw, FiberState, Section};
use fs_solid::{ForceBasedElement, Frame2d, Geometry, SectionModel, SolidError};

const E: f64 = 200e9;
const SY: f64 = 200e6;
const HI: f64 = 2e9;
const HK: f64 = 18e9;
const Q: f64 = 80e6;
const RATE: f64 = 500.0;
const C: [f64; 3] = [12e9, 30e9, 8e9];
const GAMMA: [f64; 3] = [50.0, 500.0, 2500.0];
const AREA: f64 = 1e-3;
const ARM: f64 = 0.1;

fn backbone() -> VoceJ2 {
    VoceJ2::new(
        J2Plasticity::new_mixed(IsotropicElastic::new(E, 0.3, 0.04).unwrap(), SY, HI, HK).unwrap(),
        Q,
        RATE,
    )
    .unwrap()
}
fn axial<const N: usize>() -> UniaxialChaboche<N> {
    UniaxialChaboche::new(
        ChabocheJ2::new(
            backbone(),
            std::array::from_fn(|i| {
                BackstressBranch::new(6e9 + 2e9 * i as f64, 50.0 * (i + 1).pow(2) as f64).unwrap()
            }),
        )
        .unwrap(),
    )
    .unwrap()
}
fn three() -> UniaxialChaboche<3> {
    UniaxialChaboche::new(
        ChabocheJ2::new(
            backbone(),
            std::array::from_fn(|i| BackstressBranch::new(C[i], GAMMA[i]).unwrap()),
        )
        .unwrap(),
    )
    .unwrap()
}
fn section(law: &FiberLaw) -> Section {
    Section::new(
        [-ARM, ARM]
            .map(|y| Fiber {
                y,
                area: AREA,
                law: law.clone(),
            })
            .to_vec(),
    )
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "{actual:.15e} != {expected:.15e}",
    );
}

fn check_dispatch<const N: usize>() {
    let law = axial::<N>();
    let material = FiberMaterial::chaboche_branches(law.clone(), 7850.0).unwrap();
    assert!(material.source().is_none());
    assert!(material.chaboche_branch_properties().is_none());
    match material.law() {
        FiberLaw::Chaboche(_) => assert_eq!(N, 1),
        FiberLaw::ChabocheBranches(branches) => assert_eq!(branches.branch_count(), N),
        _ => panic!("native Chaboche dispatch required"),
    }
    let mut element = section(material.law());
    let mut old = std::array::from_fn::<_, 2, _>(|_| law.initial_state());
    for (eps0, kappa) in [
        (0.0, 0.004),
        (0.0, 0.04),
        (0.0002, -0.03),
        (-0.0004, 0.055),
        (0.0, 0.0),
    ] {
        let response = element.try_respond(eps0, kappa).unwrap();
        let (mut n, mut m, mut k) = (0.0, 0.0, [[0.0; 2]; 2]);
        for (i, y) in [-ARM, ARM].into_iter().enumerate() {
            let next = law.response(eps0 - y * kappa, &old[i]).unwrap();
            n += AREA * next.stress_pa;
            m -= AREA * y * next.stress_pa;
            let ea = AREA * next.tangent_pa;
            k[0][0] += ea;
            k[0][1] -= ea * y;
            k[1][0] -= ea * y;
            k[1][1] += ea * y * y;
            old[i] = next.state;
        }
        close(response.n, n, 2e-12);
        close(response.m, m, 2e-12);
        for i in 0..2 {
            for j in 0..2 {
                close(response.tangent[i][j], k[i][j], 2e-12);
            }
        }
        element.commit(eps0, kappa);
    }
}

#[test]
fn all_supported_branch_counts_preserve_the_exact_native_cyclic_response() {
    check_dispatch::<1>();
    check_dispatch::<2>();
    check_dispatch::<3>();
    check_dispatch::<4>();
    check_dispatch::<5>();
    check_dispatch::<6>();
    check_dispatch::<7>();
    check_dispatch::<8>();
    let empty = UniaxialChaboche::new(ChabocheJ2::new(backbone(), []).unwrap()).unwrap();
    assert!(FiberMaterial::chaboche_branches(empty, 7850.0).is_err());
    assert!(FiberMaterial::chaboche_branches(three(), f64::NAN).is_err());
}

fn stress_at_plastic(p: f64) -> f64 {
    SY + (HI + HK) * p
        + Q * (1.0 - (-RATE * p).exp())
        + (0..3).map(|i| C[i] * p / (1.0 + GAMMA[i] * p)).sum::<f64>()
}

#[test]
fn three_distinct_recovery_scales_drive_bending_tangent_and_each_stored_backstress() {
    let material = FiberMaterial::chaboche_branches(three(), 7850.0).unwrap();
    let mut element = section(material.law());
    let p = 0.0015;
    let stress = stress_at_plastic(p);
    let curvature = (p + stress / E) / ARM;
    let hardening = HI
        + HK
        + Q * RATE * (-RATE * p).exp()
        + (0..3)
            .map(|i| C[i] / (1.0 + GAMMA[i] * p).powi(2))
            .sum::<f64>();
    let response = element.try_respond(0.0, curvature).unwrap();
    close(response.m, 2.0 * AREA * ARM * stress, 3e-12);
    close(
        response.tangent[1][1],
        2.0 * AREA * ARM * ARM * E * hardening / (E + hardening),
        3e-12,
    );
    element.commit(0.0, curvature);
    for (i, state) in element.states.iter().enumerate() {
        let FiberState::ChabocheBranches(ChabocheFiberState::Three(state)) = state else {
            panic!("three-branch state required")
        };
        close(state.equivalent_plastic_strain(), p, 2e-14);
        let sign = if i == 0 { 1.0 } else { -1.0 };
        for branch in 0..3 {
            close(
                state.backstress_pa()[branch],
                sign * C[branch] * p / (1.0 + GAMMA[branch] * p),
                3e-12,
            );
        }
    }
    close(element.try_respond(0.0, p / ARM).unwrap().m, 0.0, 3e-7);
    let reverse = -0.04;
    let h = 1e-7;
    let tangent = element.try_respond(0.0, reverse).unwrap().tangent[1][1];
    let fd = (element.try_respond(0.0, reverse + h).unwrap().m
        - element.try_respond(0.0, reverse - h).unwrap().m)
        / (2.0 * h);
    close(tangent, fd, 2e-8);
}

#[test]
fn domain_and_branch_count_refusals_preserve_history_and_raw_nan_behavior() {
    let mut element = section(&FiberLaw::chaboche_branches(three()).unwrap());
    element.commit(0.0, 0.03);
    let old: Vec<_> = element
        .states
        .iter()
        .map(|s| {
            let FiberState::ChabocheBranches(s) = s else {
                unreachable!()
            };
            s.clone()
        })
        .collect();
    assert!(matches!(
        element.try_respond(0.039, 0.0),
        Err(SolidError::MaterialRefused { .. })
    ));
    assert!(element.respond(0.039, 0.0).n.is_nan());
    for (actual, expected) in element.states.iter().zip(old) {
        let FiberState::ChabocheBranches(actual) = actual else {
            unreachable!()
        };
        assert_eq!(*actual, expected);
    }
    let mut wrong_multi = element.clone();
    wrong_multi.fibers[0].law = FiberLaw::chaboche_branches(axial::<2>()).unwrap();
    let mut wrong_single = element.clone();
    wrong_single.fibers[0].law = FiberLaw::chaboche_branches(axial::<1>()).unwrap();
    let mut wrong_old_family = section(&FiberLaw::chaboche_branches(axial::<1>()).unwrap());
    wrong_old_family.fibers[0].law = FiberLaw::chaboche_branches(three()).unwrap();
    for bad in [wrong_multi, wrong_single, wrong_old_family] {
        assert!(matches!(
            bad.try_respond(0.001, 0.0),
            Err(SolidError::MaterialRefused { .. })
        ));
        let raw = bad.respond(0.001, 0.0);
        assert!(raw.n.is_nan() && raw.tangent[0][0].is_nan());
    }
}

fn schema() -> [BackstressProperties; 3] {
    // Deliberately neither lexical nor magnitude order.
    [
        BackstressProperties::new("z_modulus", "z_recovery"),
        BackstressProperties::new("a_modulus", "a_recovery"),
        BackstressProperties::new("m_modulus", "m_recovery"),
    ]
}
fn insert(claims: &mut ClaimSet, name: &str, value: f64, dims: Dims) {
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "authored independent recovery branches".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
}
fn source() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("rho", 7850.0, Density::DIMS),
        ("young_modulus", E, Pressure::DIMS),
        ("poisson_ratio", 0.3, Dims::NONE),
        ("yield_stress", SY, Pressure::DIMS),
        ("isotropic_hardening_modulus", HI, Pressure::DIMS),
        ("kinematic_hardening_modulus", HK, Pressure::DIMS),
        (VOCE_SATURATION_INCREMENT_PROPERTY, Q, Pressure::DIMS),
        (VOCE_SATURATION_RATE_PROPERTY, RATE, Dims::NONE),
    ] {
        insert(&mut claims, name, value, dims);
    }
    for (i, role) in schema().into_iter().enumerate() {
        insert(&mut claims, &role.modulus, C[i], Pressure::DIMS);
        insert(&mut claims, &role.recovery, GAMMA[i], Dims::NONE);
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "numerical metal".into(),
            phase: "solid".into(),
            process: "authored multibranch fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn resolved(card: &MaterialCard) -> FiberMaterial {
    FiberMaterial::resolve_chaboche_branches(
        card,
        &QueryPoint::new(),
        "rho",
        HardeningPropertyBasis::MixedPlasticStrainModuli,
        0.04,
        MaterialPropertySelection::SingleClaimOnly,
        schema(),
    )
    .unwrap()
}

#[test]
fn named_source_branches_keep_order_complete_receipts_and_rc_reference_mass() {
    let card = source();
    let material = resolved(&card);
    let bundle = material.source().unwrap();
    assert_eq!(bundle.card_identity(), card.content_hash());
    assert_eq!(bundle.properties().len(), 14);
    assert_eq!(material.density_kg_m3(), 7850.0);
    assert_eq!(
        material.chaboche_branch_properties().unwrap(),
        schema().as_slice()
    );
    let FiberLaw::ChabocheBranches(ChabocheFiberLaw::Three(law)) = material.law() else {
        panic!("three-branch native axial law required")
    };
    assert_eq!(law, &three());
    for (i, role) in material
        .chaboche_branch_properties()
        .unwrap()
        .iter()
        .enumerate()
    {
        assert_eq!(bundle.property(&role.modulus).unwrap().value_si(), C[i]);
        assert_eq!(
            bundle.property(&role.recovery).unwrap().value_si(),
            GAMMA[i]
        );
    }
    let concrete = FiberMaterial::concrete(
        ManderConcrete::new(30e6, 0.002, 30e9, 0.005).unwrap(),
        2400.0,
    )
    .unwrap();
    let rc = RcSection::new(
        RcSectionGeometry {
            depth_m: 0.5,
            width_m: 0.4,
            cover_m: 0.04,
            concrete_layers: 16,
            steel_area_per_face_m2: 0.0015,
            steel_centroid_y_m: 0.18,
        },
        concrete.clone(),
        concrete,
        material.clone(),
    )
    .unwrap();
    close(
        rc.line_mass_kg_m(),
        (0.5 * 0.4 - 0.003) * 2400.0 + 0.003 * 7850.0,
        2e-13,
    );
    assert_eq!(rc.materials()[2].source(), material.source());
    assert_eq!(
        rc.materials()[2].chaboche_branch_properties(),
        material.chaboche_branch_properties()
    );
    assert!(rc.section().try_respond(0.0, 0.008).unwrap().m > 0.0);
    let single = FiberMaterial::resolve_chaboche_branches(
        &card,
        &QueryPoint::new(),
        "rho",
        HardeningPropertyBasis::MixedPlasticStrainModuli,
        0.04,
        MaterialPropertySelection::SingleClaimOnly,
        [schema()[0].clone()],
    )
    .unwrap();
    assert!(matches!(single.law(), FiberLaw::Chaboche(_)));
    assert_eq!(
        single.chaboche_branch_properties().unwrap(),
        &[schema()[0].clone()]
    );
    assert_eq!(single.source().unwrap().properties().len(), 10);
}

#[test]
fn density_name_cannot_alias_a_branch_and_conflicting_quantity_claims_are_rejected() {
    let card = source();
    let mut claims = card.claims().clone();
    // The card owner already forbids a same-name Pa/density ambiguity. The
    // fiber boundary must also reject the role alias before querying sources.
    let incompatible = claims.insert_claim(PropertyClaim {
        key: PropertyKey::new("z_modulus", Density::DIMS),
        value: PropertyValue::Scalar {
            value: 7850.0,
            dims: Density::DIMS,
        },
        validity: ValidityDomain::unconstrained(),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "deliberate density-role alias".into(),
            license: "test-only".into(),
            artifact: None,
        },
    });
    assert!(matches!(
        incompatible,
        Err(fs_matdb::MatDbError::DimsMismatch { expected, found, .. })
            if expected == Pressure::DIMS && found == Density::DIMS
    ));
    let result = FiberMaterial::resolve_chaboche_branches(
        &card,
        &QueryPoint::new(),
        "z_modulus",
        HardeningPropertyBasis::MixedPlasticStrainModuli,
        0.04,
        MaterialPropertySelection::SingleClaimOnly,
        schema(),
    );
    assert!(matches!(
        result,
        Err(RcSectionError::Source(
            MaterialStatePointError::InvalidRequirement { property, .. }
        )) if property == "z_modulus"
    ));
}

fn reference_deflection(element: &ForceBasedElement, force: f64, residual: bool) -> f64 {
    element
        .xi
        .iter()
        .zip(&element.wi)
        .map(|(&xi, &weight)| {
            let arm = element.length * (1.0 - xi);
            let stress = force * arm / (2.0 * AREA * ARM);
            let p = if stress <= SY {
                0.0
            } else {
                let (mut lower, mut upper) = (0.0, 0.02);
                assert!(stress_at_plastic(upper) > stress);
                for _ in 0..80 {
                    let middle = 0.5 * (lower + upper);
                    if stress_at_plastic(middle) < stress {
                        lower = middle;
                    } else {
                        upper = middle;
                    }
                }
                0.5 * (lower + upper)
            };
            let strain = p + if residual { 0.0 } else { stress / E };
            weight * element.length * strain / ARM * arm
        })
        .sum()
}

#[test]
fn source_bound_multi_branch_beam_and_frame_retain_permanent_and_reversed_response() {
    let material = resolved(&source());
    let make_section = || section(material.law());
    let mut beam = ForceBasedElement::new(2.0, &make_section);
    let mut frame = Frame2d::new(Geometry::Linear);
    let base = frame.add_node(0.0, 0.0);
    let tip = frame.add_node(2.0, 0.0);
    frame.fix(base, [true; 3]);
    frame
        .add_element(base, tip, 5, 0.0, &|| SectionModel::Fiber(make_section()))
        .unwrap();
    let force = 30_000.0;
    let expected_loaded = reference_deflection(&beam, force, false);
    let expected_permanent = reference_deflection(&beam, force, true);
    let mut load = vec![0.0; frame.ndof()];
    load[3 * tip + 1] = force;
    // One physical constitutive increment in both owners: subdividing recovery
    // changes backward-Euler integration and is not an exact comparison.
    frame.static_load(&load, 1).unwrap();
    let loaded = beam.tip_deflection_under_shear(force).unwrap();
    close(loaded, expected_loaded, 2e-9);
    close(frame.displacements()[3 * tip + 1], expected_loaded, 2e-8);
    let permanent = beam.tip_deflection_under_shear(0.0).unwrap();
    frame.static_load(&vec![0.0; frame.ndof()], 1).unwrap();
    close(permanent, expected_permanent, 2e-9);
    close(frame.displacements()[3 * tip + 1], expected_permanent, 2e-8);
    assert!(permanent > 0.0 && permanent < loaded);
    load[3 * tip + 1] = -force;
    let reversed = beam.tip_deflection_under_shear(-force).unwrap();
    frame.static_load(&load, 1).unwrap();
    assert!(reversed < 0.0);
    close(frame.displacements()[3 * tip + 1], reversed, 2e-8);
}

#[test]
fn preserved_norton_fiber_dispatch_keeps_one_physical_step_and_reference_mass() {
    let law = UniaxialNorton::new(
        NortonCreep::new(
            IsotropicElastic::new(E, 0.3, 0.04).unwrap(),
            NortonCreepParameters {
                reference_stress_pa: 100e6,
                reference_rate_s: 1e-5,
                stress_exponent: 1.0,
                activation_energy_j_mol: 0.0,
                reference_temperature_k: 300.0,
                temperature_range_k: [290.0, 310.0],
            },
            300.0,
            10.0,
        )
        .unwrap(),
    );
    let material = FiberMaterial::authored(FiberLaw::Norton(law), 7850.0).unwrap();
    close(material.density_kg_m3(), 7850.0, 0.0);
    assert!(material.source().is_none());
    assert!(material.chaboche_branch_properties().is_none());
    let mut element = section(material.law());
    let curvature = 0.02;
    let relaxation = 1.0 + E * 10.0 * 1e-5 / 100e6;
    let stress = E * ARM * curvature / relaxation;
    let response = element.try_respond(0.0, curvature).unwrap();
    close(response.m, 2.0 * AREA * ARM * stress, 2e-12);
    close(
        response.tangent[1][1],
        2.0 * AREA * ARM * ARM * E / relaxation,
        2e-12,
    );
    close(
        element.try_respond(0.0, curvature).unwrap().m,
        response.m,
        0.0,
    );
    for state in &element.states {
        let FiberState::Norton(state) = state else {
            panic!("native creep state required")
        };
        assert_eq!(state.creep_strain(), 0.0);
    }
    element.commit(0.0, curvature);
    for (i, state) in element.states.iter().enumerate() {
        let FiberState::Norton(state) = state else {
            unreachable!()
        };
        let sign = if i == 0 { 1.0 } else { -1.0 };
        close(
            state.creep_strain(),
            sign * (ARM * curvature - stress / E),
            2e-15,
        );
    }
    close(
        element.try_respond(0.0, curvature).unwrap().m,
        response.m / relaxation,
        2e-12,
    );
}
