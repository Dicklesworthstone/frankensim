//! G0/G1/G3: the shared J2 law drives existing metal fibers, beam-columns and
//! frame sections. Source fixtures are authored coefficients, not calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::elastic::IsotropicElastic;
use fs_material::fiber::ManderConcrete;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::chaboche::axial::UniaxialChaboche;
use fs_material::plastic::chaboche::source::{
    CHABOCHE_BACKSTRESS_MODULUS_PROPERTY, CHABOCHE_DYNAMIC_RECOVERY_PROPERTY,
};
use fs_material::plastic::chaboche::{BackstressBranch, ChabocheJ2};
use fs_material::plastic::reduction::{UniaxialJ2, UniaxialVoce};
use fs_material::plastic::source::HardeningPropertyBasis;
use fs_material::plastic::voce::VoceJ2;
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY, VOCE_SATURATION_RATE_PROPERTY,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Density, Dims, Pressure};
use fs_solid::fiber::material::{FiberMaterial, RcSection, RcSectionGeometry};
use fs_solid::fiber::{Fiber, FiberLaw, FiberState, Section};
use fs_solid::{ForceBasedElement, Frame2d, Geometry, SectionModel};

const E: f64 = 200e9;
const YIELD: f64 = 200e6;
const H_ISO: f64 = 2e9;
const H_KIN: f64 = 18e9;
const AREA: f64 = 1e-3;
const ARM: f64 = 0.1;

fn law(limit: f64) -> UniaxialJ2 {
    UniaxialJ2::new(
        J2Plasticity::new_mixed(
            IsotropicElastic::new(E, 0.3, limit).unwrap(),
            YIELD,
            H_ISO,
            H_KIN,
        )
        .unwrap(),
    )
    .unwrap()
}

fn section() -> Section {
    Section::new(
        [-ARM, ARM]
            .map(|y| Fiber {
                y,
                area: AREA,
                law: FiberLaw::J2(law(0.04)),
            })
            .to_vec(),
    )
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn source(include_kinematic: bool, tangent_basis: bool) -> MaterialCard {
    let h = H_ISO + H_KIN;
    let mut values = vec![
        ("density", 7850.0, Density::DIMS),
        ("young_modulus", E, Pressure::DIMS),
        ("poisson_ratio", 0.3, Dims::NONE),
        ("yield_stress", YIELD, Pressure::DIMS),
    ];
    if tangent_basis {
        values.push((
            "post_yield_tangent_modulus",
            E * (h / (E + h)),
            Pressure::DIMS,
        ));
    } else {
        values.push(("isotropic_hardening_modulus", H_ISO, Pressure::DIMS));
        if include_kinematic {
            values.push(("kinematic_hardening_modulus", H_KIN, Pressure::DIMS));
        }
    }
    let mut claims = ClaimSet::new();
    for (name, value, dims) in values {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::RelativeHalfWidth {
                    fraction: 0.01,
                    confidence: 0.95,
                },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored mixed-J2 metal fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "numerical metal".into(),
            phase: "solid".into(),
            process: "authored coefficient fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn resolved(card: &MaterialCard, basis: HardeningPropertyBasis) -> FiberMaterial {
    FiberMaterial::resolve_j2(
        card,
        &QueryPoint::new(),
        "density",
        basis,
        0.04,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn sourced_j2_reinforcement_keeps_native_coefficients_receipts_and_reference_mass() {
    let card = source(true, false);
    let metal = resolved(&card, HardeningPropertyBasis::MixedPlasticStrainModuli);
    assert_eq!(metal.source().unwrap().properties().len(), 6);
    assert_eq!(metal.density_kg_m3(), 7850.0);
    let FiberLaw::J2(j2) = metal.law() else {
        panic!("expected the native J2 law")
    };
    assert_eq!(j2.backbone(), law(0.04).backbone());
    let concrete = FiberMaterial::concrete(
        ManderConcrete::new(30e6, 0.002, 30e9, 0.005).unwrap(),
        2400.0,
    )
    .unwrap();
    let geometry = RcSectionGeometry {
        depth_m: 0.5,
        width_m: 0.4,
        cover_m: 0.04,
        concrete_layers: 16,
        steel_area_per_face_m2: 0.0015,
        steel_centroid_y_m: 0.18,
    };
    let rc = RcSection::new(geometry, concrete.clone(), concrete, metal.clone()).unwrap();
    let steel_area = 2.0 * geometry.steel_area_per_face_m2;
    close(
        rc.line_mass_kg_m(),
        2400.0 * (geometry.depth_m * geometry.width_m - steel_area) + 7850.0 * steel_area,
        2e-13,
    );
    assert_eq!(rc.materials()[2].source(), metal.source());
    assert_eq!(
        rc.section()
            .states
            .iter()
            .filter(|s| matches!(s, FiberState::J2(_)))
            .count(),
        2
    );
    assert!(rc.section().try_respond(0.0, 0.008).unwrap().m > 0.0);
}

#[test]
fn hardening_sources_are_explicit_and_monotonic_tangent_does_not_invent_a_split() {
    assert!(
        FiberMaterial::resolve_j2(
            &source(false, false),
            &QueryPoint::new(),
            "density",
            HardeningPropertyBasis::MixedPlasticStrainModuli,
            0.04,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let metal = resolved(
        &source(false, true),
        HardeningPropertyBasis::UniaxialTangentModulus,
    );
    let FiberLaw::J2(j2) = metal.law() else {
        panic!("expected J2")
    };
    close(j2.backbone().hardening, H_ISO + H_KIN, 2e-14);
    assert_eq!(j2.backbone().kinematic_hardening, 0.0);
    assert!(FiberMaterial::j2(law(0.04), f64::NAN).is_err());
    assert!(
        FiberMaterial::resolve_j2(
            &source(true, false),
            &QueryPoint::new(),
            "missing-density",
            HardeningPropertyBasis::MixedPlasticStrainModuli,
            0.04,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}

#[test]
fn metal_section_has_exact_bilinear_bending_and_residual_curvature() {
    let mut section = section();
    let curvature = 0.025;
    let h = H_ISO + H_KIN;
    let et = E * (h / (E + h));
    let expected_stress = YIELD + et * (ARM * curvature - YIELD / E);
    let r = section.try_respond(0.0, curvature).unwrap();
    close(r.n, 0.0, 1e-9);
    close(r.m, 2.0 * AREA * ARM * expected_stress, 2e-13);
    close(r.tangent[1][1], 2.0 * et * AREA * ARM * ARM, 2e-13);
    let delta = 1e-7;
    let fd = (section.try_respond(0.0, curvature + delta).unwrap().m
        - section.try_respond(0.0, curvature - delta).unwrap().m)
        / (2.0 * delta);
    close(fd, r.tangent[1][1], 1e-9);
    section.commit(0.0, curvature);
    let residual = curvature - r.m / (2.0 * E * AREA * ARM * ARM);
    assert!(residual > 0.0);
    let unloaded = section.try_respond(0.0, residual).unwrap();
    close(unloaded.m, 0.0, 2e-8);
    close(unloaded.tangent[1][1], 2.0 * E * AREA * ARM * ARM, 2e-13);
    let mut full_history = Vec::new();
    for state in &section.states {
        let FiberState::J2(s) = state else {
            panic!("expected J2 state")
        };
        full_history.push(s.tensor_history());
    }
    assert!(full_history.iter().all(|s| s.alpha > 0.0));
    close(
        full_history[0].plastic_strain.iter().take(3).sum(),
        0.0,
        1e-15,
    );
}

fn reference_deflection(element: &ForceBasedElement, load: f64, residual: bool) -> f64 {
    let et = E * ((H_ISO + H_KIN) / (E + H_ISO + H_KIN));
    element
        .xi
        .iter()
        .zip(&element.wi)
        .map(|(&xi, &weight)| {
            let arm = element.length * (1.0 - xi);
            let sigma = load * arm / (2.0 * AREA * ARM);
            let axial = if sigma <= YIELD {
                sigma / E
            } else {
                YIELD / E + (sigma - YIELD) / et
            };
            let curvature = (axial - if residual { sigma / E } else { 0.0 }) / ARM;
            weight * element.length * curvature * arm
        })
        .sum()
}

#[test]
fn existing_force_based_beam_loads_yields_and_keeps_permanent_deflection() {
    let mut beam = ForceBasedElement::new(2.0, &section);
    let load = 30_000.0;
    let expected = reference_deflection(&beam, load, false);
    for fraction in [0.25, 0.5, 0.75] {
        beam.tip_deflection_under_shear(fraction * load).unwrap();
    }
    let deflection = beam.tip_deflection_under_shear(load).unwrap();
    close(deflection, expected, 2e-10);
    assert!(deflection > load * beam.length.powi(3) / (6.0 * E * AREA * ARM * ARM));
    let permanent = beam.tip_deflection_under_shear(0.0).unwrap();
    close(permanent, reference_deflection(&beam, load, true), 2e-10);
    assert!(permanent > 0.0);
}

#[test]
fn existing_frame_uses_the_same_metal_section_history() {
    let mut frame = Frame2d::new(Geometry::Linear);
    let base = frame.add_node(0.0, 0.0);
    let tip = frame.add_node(2.0, 0.0);
    frame.fix(base, [true; 3]);
    frame
        .add_element(base, tip, 5, 0.0, &|| SectionModel::Fiber(section()))
        .unwrap();
    let mut load = vec![0.0; frame.ndof()];
    load[3 * tip + 1] = 30_000.0;
    frame.static_load(&load, 8).unwrap();
    let beam = ForceBasedElement::new(2.0, &section);
    close(
        frame.displacements()[3 * tip + 1],
        reference_deflection(&beam, 30_000.0, false),
        1e-8,
    );
    frame.static_load(&vec![0.0; frame.ndof()], 8).unwrap();
    close(
        frame.displacements()[3 * tip + 1],
        reference_deflection(&beam, 30_000.0, true),
        1e-8,
    );
}

const C_RECOVERY: f64 = 8e9;
const GAMMA_RECOVERY: f64 = 500.0;

fn recovery_section() -> Section {
    let axial = UniaxialChaboche::new(
        ChabocheJ2::new(
            VoceJ2::new(law(0.04).backbone().clone(), 80e6, 500.0).unwrap(),
            [BackstressBranch::new(C_RECOVERY, GAMMA_RECOVERY).unwrap()],
        )
        .unwrap(),
    )
    .unwrap();
    Section::new(
        [-ARM, ARM]
            .map(|y| Fiber {
                y,
                area: AREA,
                law: FiberLaw::Chaboche(axial.clone()),
            })
            .to_vec(),
    )
}

#[test]
fn recovering_backstress_drives_beam_deflection_permanent_set_and_reversed_section_tangent() {
    let mut beam = ForceBasedElement::new(2.0, &recovery_section);
    let force = 30_000.0;
    // Independent stress-parametrized monotonic backward-Euler law. Inverting
    // this scalar reference does not use the tensor or fiber return map.
    let plastic_at_stress = |stress: f64| {
        if stress <= YIELD {
            return 0.0;
        }
        let (mut lo, mut hi) = (0.0_f64, 0.03_f64);
        for _ in 0..70 {
            let p = 0.5 * (lo + hi);
            let s = YIELD
                + (H_ISO + H_KIN) * p
                + 80e6 * (1.0 - (-500.0 * p).exp())
                + C_RECOVERY * p / (1.0 + GAMMA_RECOVERY * p);
            if s < stress {
                lo = p;
            } else {
                hi = p;
            }
        }
        0.5 * (lo + hi)
    };
    let mut loaded_reference = 0.0;
    let mut permanent_reference = 0.0;
    for (&xi, &weight) in beam.xi.iter().zip(&beam.wi) {
        let arm = beam.length * (1.0 - xi);
        let stress = force * arm / (2.0 * AREA * ARM);
        let plastic = plastic_at_stress(stress);
        let scale = weight * beam.length * arm / ARM;
        loaded_reference += scale * (plastic + stress / E);
        permanent_reference += scale * plastic;
    }
    let loaded = beam.tip_deflection_under_shear(force).unwrap();
    close(loaded, loaded_reference, 3e-9);
    let permanent = beam.tip_deflection_under_shear(0.0).unwrap();
    close(permanent, permanent_reference, 3e-9);
    assert!(permanent > 0.0);
    let mut section = recovery_section();
    section.try_respond(0.0, 0.04).unwrap();
    section.commit(0.0, 0.04);
    let accepted: Vec<_> = section
        .states
        .iter()
        .map(|s| {
            let FiberState::Chaboche(s) = s else {
                panic!("recovery history required")
            };
            s.tensor_history().clone()
        })
        .collect();
    let curvature = -0.03;
    let h = 1e-7;
    let response = section.try_respond(0.0, curvature).unwrap();
    let fd = (section.try_respond(0.0, curvature + h).unwrap().m
        - section.try_respond(0.0, curvature - h).unwrap().m)
        / (2.0 * h);
    close(response.tangent[1][1], fd, 2e-8);
    assert!(section.try_respond(0.04, 0.0).is_err());
    assert!(section.respond(0.04, 0.0).n.is_nan());
    for (state, old) in section.states.iter().zip(&accepted) {
        let FiberState::Chaboche(state) = state else {
            panic!("recovery history required")
        };
        assert_eq!(state.tensor_history(), old);
    }
    section.commit(0.0, curvature);
    for (state, old) in section.states.iter().zip(&accepted) {
        let FiberState::Chaboche(state) = state else {
            panic!("recovery history required")
        };
        assert!(state.equivalent_plastic_strain() > old.plastic.alpha);
        assert!(state.tensor_history().backstress_pa[0][0] * old.backstress_pa[0][0] < 0.0);
    }
}

#[test]
fn sourced_recovery_keeps_complete_receipts_mass_and_reaches_native_fibers() {
    let original = source(true, false);
    let mut claims = original.claims().clone();
    for (name, value, dims) in [
        (VOCE_SATURATION_INCREMENT_PROPERTY, 80e6, Pressure::DIMS),
        (VOCE_SATURATION_RATE_PROPERTY, 500.0, Dims::NONE),
        (
            CHABOCHE_BACKSTRESS_MODULUS_PROPERTY,
            C_RECOVERY,
            Pressure::DIMS,
        ),
        (
            CHABOCHE_DYNAMIC_RECOVERY_PROPERTY,
            GAMMA_RECOVERY,
            Dims::NONE,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored nonlinear-recovery fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let material = FiberMaterial::resolve_chaboche(
        &card,
        &QueryPoint::new(),
        "density",
        HardeningPropertyBasis::MixedPlasticStrainModuli,
        0.04,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(
        material.source().unwrap().card_identity(),
        card.content_hash()
    );
    assert_eq!(material.source().unwrap().properties().len(), 10);
    assert_eq!(material.density_kg_m3(), 7850.0);
    assert!(matches!(material.law(), FiberLaw::Chaboche(_)));
    let sourced = Section::new(
        [-ARM, ARM]
            .map(|y| Fiber {
                y,
                area: AREA,
                law: material.law().clone(),
            })
            .to_vec(),
    );
    let expected = recovery_section().try_respond(0.0, 0.03).unwrap();
    let actual = sourced.try_respond(0.0, 0.03).unwrap();
    assert_eq!(actual.m, expected.m);
    assert_eq!(actual.tangent, expected.tangent);
    assert!(
        FiberMaterial::resolve_chaboche(
            &original,
            &QueryPoint::new(),
            "density",
            HardeningPropertyBasis::MixedPlasticStrainModuli,
            0.04,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}

#[test]
fn checked_fiber_response_admits_full_tensor_strain_and_preserves_state_on_refusal() {
    let mut section = section();
    section.commit(0.0, 0.02);
    let snapshot: Vec<_> = section
        .states
        .iter()
        .map(|s| {
            let FiberState::J2(state) = s else {
                panic!("expected J2")
            };
            state.tensor_history()
        })
        .collect();
    // Axial strain alone is below 0.04, but the compatible lateral plastic
    // contraction makes the complete tensor exceed the admitted norm.
    assert!(section.try_respond(0.039, 0.0).is_err());
    assert!(section.try_respond(f64::NAN, 0.0).is_err());
    for (state, expected) in section.states.iter().zip(snapshot) {
        let FiberState::J2(state) = state else {
            panic!("expected J2")
        };
        assert_eq!(state.tensor_history(), expected);
    }
}

fn saturation_section() -> Section {
    let voce =
        UniaxialVoce::new(VoceJ2::new(law(0.04).backbone().clone(), 80e6, 500.0).unwrap()).unwrap();
    Section::new(
        [-ARM, ARM]
            .map(|y| Fiber {
                y,
                area: AREA,
                law: FiberLaw::Voce(voce.clone()),
            })
            .to_vec(),
    )
}

#[test]
fn saturating_section_uses_the_nonlinear_tangent_and_keeps_residual_curvature() {
    let mut section = saturation_section();
    let p: f64 = 0.002;
    let stress = YIELD + (H_ISO + H_KIN) * p + 80e6 * (1.0 - (-500.0 * p).exp());
    let curvature = (p + stress / E) / ARM;
    let response = section.try_respond(0.0, curvature).unwrap();
    let h = H_ISO + H_KIN + 80e6 * 500.0 * (-500.0 * p).exp();
    close(response.m, 2.0 * AREA * ARM * stress, 2e-12);
    close(
        response.tangent[1][1],
        2.0 * AREA * ARM * ARM * E * (h / (E + h)),
        2e-12,
    );
    section.commit(0.0, curvature);
    let residual = section.try_respond(0.0, p / ARM).unwrap();
    close(residual.m, 0.0, 2e-7);
    for state in &section.states {
        let FiberState::Voce(state) = state else {
            panic!("Voce history required")
        };
        close(state.equivalent_plastic_strain(), p, 2e-14);
    }
    let accepted: Vec<_> = section
        .states
        .iter()
        .map(|s| {
            let FiberState::Voce(s) = s else {
                panic!("Voce history required")
            };
            s.tensor_history()
        })
        .collect();
    assert!(section.try_respond(0.04, 0.0).is_err());
    for (state, expected) in section.states.iter().zip(accepted) {
        let FiberState::Voce(state) = state else {
            panic!("Voce history required")
        };
        assert_eq!(state.tensor_history(), expected);
    }
}

#[test]
fn source_bound_saturation_reaches_metal_fibers_with_all_receipts_and_mass() {
    let original = source(true, false);
    let mut claims = original.claims().clone();
    for (name, value, dims) in [
        (VOCE_SATURATION_INCREMENT_PROPERTY, 80e6, Pressure::DIMS),
        (VOCE_SATURATION_RATE_PROPERTY, 500.0, Dims::NONE),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored saturation fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let material = FiberMaterial::resolve_voce(
        &card,
        &QueryPoint::new(),
        "density",
        HardeningPropertyBasis::MixedPlasticStrainModuli,
        0.04,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(
        material.source().unwrap().card_identity(),
        card.content_hash()
    );
    assert_eq!(material.source().unwrap().properties().len(), 8);
    assert_eq!(material.density_kg_m3(), 7850.0);
    let section = Section::new(
        [-ARM, ARM]
            .map(|y| Fiber {
                y,
                area: AREA,
                law: material.law().clone(),
            })
            .to_vec(),
    );
    let expected = saturation_section().try_respond(0.0, 0.03).unwrap();
    let actual = section.try_respond(0.0, 0.03).unwrap();
    assert_eq!(actual.n, expected.n);
    assert_eq!(actual.m, expected.m);
    assert_eq!(actual.tangent, expected.tangent);
    assert!(
        FiberMaterial::resolve_voce(
            &original,
            &QueryPoint::new(),
            "density",
            HardeningPropertyBasis::MixedPlasticStrainModuli,
            0.04,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}

#[test]
fn existing_beam_and_frame_share_saturating_history_and_permanent_deflection() {
    let force = 36_000.0;
    let mut beam = ForceBasedElement::new(2.0, &saturation_section);
    let mut frame = Frame2d::new(Geometry::Linear);
    let base = frame.add_node(0.0, 0.0);
    let tip = frame.add_node(2.0, 0.0);
    frame.fix(base, [true; 3]);
    frame
        .add_element(base, tip, 5, 0.0, &|| {
            SectionModel::Fiber(saturation_section())
        })
        .unwrap();
    let mut load = vec![0.0; frame.ndof()];
    load[3 * tip + 1] = force;
    frame.static_load(&load, 10).unwrap();
    let loaded = beam.tip_deflection_under_shear(force).unwrap();
    close(frame.displacements()[3 * tip + 1], loaded, 2e-8);
    let elastic = force * 2.0_f64.powi(3) / (6.0 * E * AREA * ARM * ARM);
    assert!(loaded > elastic);
    let permanent = beam.tip_deflection_under_shear(0.0).unwrap();
    frame.static_load(&vec![0.0; frame.ndof()], 10).unwrap();
    close(permanent, loaded - elastic, 2e-9);
    close(frame.displacements()[3 * tip + 1], permanent, 2e-8);
    assert!(permanent > 0.0);
}
