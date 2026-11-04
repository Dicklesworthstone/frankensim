//! Native finite-strain sheet laws; all sources below are synthetic fixtures.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::sheet::{
    AREAL_DENSITY, MEMBRANE_DIMS, MembraneLaw, ResolvedSheet, SheetCapabilities, SheetError,
    SheetLaw, THICKNESS, resolve_sheet,
};
use fs_material::state_point::{
    MaterialPropertySelection as Pick, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::Dims;

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const MU: f64 = 3.0e5;
const LAMBDA: f64 = 8.0e5;

fn provenance() -> Provenance {
    Provenance {
        source: "synthetic hyperelastic sheet verification; not empirical material data".into(),
        license: "test-only".into(),
        artifact: None,
    }
}

fn model(kind: &str) -> ConstitutiveModelCard {
    let mut values = vec![("stretch_limit", 5.0, Dims::NONE)];
    match kind {
        "neo-hookean-log" => values.extend([("mu", MU, PA), ("lambda", LAMBDA, PA)]),
        "mooney-rivlin-decoupled" => {
            values.extend([("c10", 1.1e5, PA), ("c01", 0.4e5, PA), ("kappa", 1.2e6, PA)]);
        }
        "ogden-decoupled" => values.extend([
            ("volumetric_kind", 1.0, Dims::NONE),
            ("kappa", 1.2e6, PA),
            ("term_count", 2.0, Dims::NONE),
            ("mu_0", 2.0e5, PA),
            ("alpha_0", 1.3, Dims::NONE),
            ("mu_1", -0.8e5, PA),
            ("alpha_1", -2.0, Dims::NONE),
        ]),
        "hgo-dispersed" => values.extend([
            ("volumetric_kind", 0.0, Dims::NONE),
            ("kappa", 1.2e6, PA),
            ("c10", 0.7e5, PA),
            ("k1", 1.8e5, PA),
            ("k2", 3.2, Dims::NONE),
            ("dispersion", 0.08, Dims::NONE),
            ("fibre_invariants", 1.0, Dims::NONE),
            ("fibre_count", 1.0, Dims::NONE),
            ("fibre_0_x", 0.8, Dims::NONE),
            ("fibre_0_y", 0.0, Dims::NONE),
            ("fibre_0_z", 0.6, Dims::NONE),
        ]),
        _ => panic!("unknown test law"),
    }
    ConstitutiveModelCard {
        law: LawId(kind.into()),
        law_version: 1,
        parameters: values
            .into_iter()
            .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
            .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
        sources: vec![],
        provenance: provenance(),
    }
}

fn source(
    models: Vec<ConstitutiveModelCard>,
    thickness: f64,
    mass: f64,
    reference_domain: ValidityDomain,
    extra: &[(&str, f64, Dims)],
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut rows = vec![
        (AREAL_DENSITY, mass, Dims([-2, 1, 0, 0, 0, 0])),
        (THICKNESS, thickness, Dims([1, 0, 0, 0, 0, 0])),
    ];
    rows.extend_from_slice(extra);
    for (name, value, dims) in rows {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                validity: reference_domain.clone(),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic hyperelastic sheet".into(),
            phase: "solid".into(),
            process: "constitutive integration test".into(),
            revision: 0,
        },
        claims,
        models,
    )
    .unwrap()
}

fn card(kind: &str, thickness: f64, mass: f64) -> MaterialCard {
    source(
        vec![model(kind)],
        thickness,
        mass,
        ValidityDomain::unconstrained(),
        &[],
    )
}

fn controls() -> PlaneStressSettings {
    PlaneStressSettings {
        relative_tolerance: 2.0e-13,
        ..PlaneStressSettings::default()
    }
}

fn point() -> QueryPoint {
    QueryPoint::new().with("T", 300.0).unwrap()
}

fn reference(c: &MaterialCard) -> ResolvedSheet {
    resolve_sheet(
        c,
        &point(),
        SheetCapabilities {
            membrane: false,
            bending: false,
        },
        Pick::SingleClaimOnly,
    )
    .unwrap()
}

fn bind(c: &MaterialCard) -> ResolvedSheet {
    reference(c)
        .with_hyperelastic_membrane(c, 0, controls())
        .unwrap()
}

fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn small_strain_plane_stress_and_independent_reference_channels() {
    let c = card("neo-hookean-log", 0.0013, 0.24);
    let sheet = bind(&c);
    let law = sheet.membrane_law().unwrap();
    let r = law.response([0.0; 3]).unwrap();
    let lambda_ps = 2.0 * MU * LAMBDA / (LAMBDA + 2.0 * MU);
    let expected = [
        2.0 * MU + lambda_ps,
        lambda_ps,
        0.0,
        lambda_ps,
        2.0 * MU + lambda_ps,
        0.0,
        0.0,
        0.0,
        MU,
    ];
    close(r.energy_j_m2, 0.0, 1e-12, 0.0);
    for s in r.resultants_n_m {
        close(s, 0.0, 1e-10, 0.0);
    }
    for (a, b) in r.tangent_n_m.into_iter().zip(expected) {
        close(a, sheet.thickness_m() * b, 1e-9, 2e-13);
    }
    law.admit_linearization([0.0; 3]).unwrap();
    assert!(sheet.membrane().is_err());
    assert!(sheet.bending().is_err());
    let MembraneLaw::Hyperelastic(bound) = &law else {
        panic!("wrong selected law")
    };
    assert_eq!(bound.material(), sheet.material());
    assert_eq!(bound.model().model_card(), &c.models()[0]);
    assert_eq!(bound.model().context(), &sheet.material().query());
    assert_eq!(bound.reference_thickness_m(), sheet.thickness_m());

    let e = [0.17, -0.04, 0.13];
    let original = law.response(e).unwrap();
    let thick = bind(&card("neo-hookean-log", 0.0039, 0.24));
    let heavy = bind(&card("neo-hookean-log", 0.0013, 0.72));
    let doubled = thick.membrane_law().unwrap().response(e).unwrap();
    let mass_changed = heavy.membrane_law().unwrap().response(e).unwrap();
    close(
        doubled.energy_j_m2,
        3.0 * original.energy_j_m2,
        1e-11,
        1e-13,
    );
    assert_eq!(mass_changed.energy_j_m2, original.energy_j_m2);
    assert_eq!(mass_changed.resultants_n_m, original.resultants_n_m);
    assert_eq!(mass_changed.tangent_n_m, original.tangent_n_m);
    for i in 0..3 {
        close(
            doubled.resultants_n_m[i],
            3.0 * original.resultants_n_m[i],
            1e-10,
            1e-13,
        );
    }
    for i in 0..9 {
        close(
            doubled.tangent_n_m[i],
            3.0 * original.tangent_n_m[i],
            1e-9,
            1e-13,
        );
    }
    assert_eq!(thick.areal_density_kg_m2(), sheet.areal_density_kg_m2());
    assert_eq!(heavy.areal_density_kg_m2(), 0.72);
}

// Independent compressible Neo-Hookean reduction in C, with a scalar bisection
// for thickness. No native plane-stress, Cholesky map or AD enters this oracle.
fn neo_hookean_oracle(e: [f64; 3], thickness: f64) -> (f64, [f64; 3], [f64; 9]) {
    let c = [[1.0 + 2.0 * e[0], e[2]], [e[2], 1.0 + 2.0 * e[1]]];
    let determinant = c[0][0] * c[1][1] - c[0][1] * c[1][0];
    let area = determinant.sqrt();
    let mut lo: f64 = 1e-8;
    let mut hi: f64 = 10.0;
    for _ in 0..100 {
        let z = 0.5 * (lo + hi);
        let traction = MU * (z * z - 1.0) + LAMBDA * (area * z).ln();
        if traction > 0.0 {
            hi = z;
        } else {
            lo = z;
        }
    }
    let z = 0.5 * (lo + hi);
    let log_j = (area * z).ln();
    let energy = thickness
        * (0.5 * MU * (c[0][0] + c[1][1] + z * z - 3.0) - MU * log_j
            + 0.5 * LAMBDA * log_j * log_j);
    let inv = [
        [c[1][1] / determinant, -c[0][1] / determinant],
        [-c[1][0] / determinant, c[0][0] / determinant],
    ];
    let slots = [(0, 0), (1, 1), (0, 1)];
    let resultants =
        slots.map(|(i, j)| thickness * MU * ((if i == j { 1.0 } else { 0.0 }) - z * z * inv[i][j]));
    let mut tangent = [0.0; 9];
    for col in 0..3 {
        let dc = match col {
            0 => [[2.0, 0.0], [0.0, 0.0]],
            1 => [[0.0, 0.0], [0.0, 2.0]],
            _ => [[0.0, 1.0], [1.0, 0.0]],
        };
        let (i, j) = slots[col];
        let dz2 = -2.0 * LAMBDA * z * z * inv[i][j] / (2.0 * MU * z * z + LAMBDA);
        for (row, (i, j)) in slots.into_iter().enumerate() {
            let mut inverse_variation = 0.0;
            for a in 0..2 {
                for b in 0..2 {
                    inverse_variation -= inv[i][a] * dc[a][b] * inv[b][j];
                }
            }
            tangent[3 * row + col] =
                -thickness * MU * (dz2 * inv[i][j] + z * z * inverse_variation);
        }
    }
    (energy, resultants, tangent)
}

#[test]
fn finite_compression_biaxial_stretch_and_shear_match_independent_reduction() {
    let sheet = bind(&card("neo-hookean-log", 0.0013, 0.24));
    let law = sheet.membrane_law().unwrap();
    for strain in [[0.16, 0.07, 0.11], [-0.12, 0.03, 0.06], [1.1, 0.17, 0.35]] {
        let r = law.response(strain).unwrap();
        let (energy, stress, tangent) = neo_hookean_oracle(strain, sheet.thickness_m());
        close(r.energy_j_m2, energy, 2e-10, 2e-12);
        for i in 0..3 {
            close(r.resultants_n_m[i], stress[i], 2e-9, 2e-11);
        }
        for i in 0..9 {
            close(r.tangent_n_m[i], tangent[i], 2e-8, 2e-11);
        }
    }
}

#[test]
fn native_models_have_work_conjugate_green_strain_tangents() {
    for kind in [
        "neo-hookean-log",
        "mooney-rivlin-decoupled",
        "ogden-decoupled",
        "hgo-dispersed",
    ] {
        let sheet = bind(&card(kind, 0.0013, 0.24));
        let law = sheet.membrane_law().unwrap();
        for e in [[0.16, 0.07, 0.11], [-0.08, 0.02, -0.05]] {
            let r = law.response(e).unwrap();
            let h = 2e-6;
            for col in 0..3 {
                let mut lo = e;
                let mut hi = e;
                lo[col] -= h;
                hi[col] += h;
                let a = law.response(lo).unwrap();
                let b = law.response(hi).unwrap();
                close(
                    (b.energy_j_m2 - a.energy_j_m2) / (2.0 * h),
                    r.resultants_n_m[col],
                    2e-6,
                    2e-7,
                );
                for row in 0..3 {
                    close(
                        (b.resultants_n_m[row] - a.resultants_n_m[row]) / (2.0 * h),
                        r.tangent_n_m[3 * row + col],
                        2e-5,
                        3e-7,
                    );
                    close(
                        r.tangent_n_m[3 * row + col],
                        r.tangent_n_m[3 * col + row],
                        2e-10,
                        2e-12,
                    );
                }
            }
        }
    }
}

#[test]
fn tilted_fibres_relax_transverse_shear_and_keep_activation_explicit() {
    for invariant_kind in [0.0, 1.0] {
        let mut m = model("hgo-dispersed");
        m.parameters.get_mut("fibre_invariants").unwrap().value = invariant_kind;
        let c = source(vec![m], 0.0013, 0.24, ValidityDomain::unconstrained(), &[]);
        let law = bind(&c).membrane_law().unwrap();
        let MembraneLaw::Hyperelastic(bound) = &law else {
            panic!("wrong law")
        };
        let e = [0.22, 0.0, 0.05];
        let r = (1.0_f64 + 2.0 * e[0]).sqrt();
        let s = e[2] / r;
        let t = (1.0 + 2.0 * e[1] - s * s).sqrt();
        let native = bound
            .model()
            .law()
            .plane_stress(&[r, s, 0.0, t], controls())
            .unwrap();
        assert!(
            native.deformation_gradient[2].abs() > 1e-3,
            "tilted fibres must relax transverse shear"
        );
        let piola = bound
            .model()
            .law()
            .piola(&native.deformation_gradient)
            .unwrap();
        for i in [2, 5, 8] {
            close(piola[i], 0.0, 1e-6, 0.0);
        }
        close(
            law.response(e).unwrap().energy_j_m2,
            0.0013 * native.energy_density_pa,
            1e-11,
            1e-12,
        );
        law.admit_linearization(e).unwrap();
        assert!(law.admit_linearization([0.0; 3]).is_err());
        assert!(law.response([0.0; 3]).is_ok()); // original one-sided nonlinear convention
        let x = law.response([0.15, 0.0, 0.0]).unwrap();
        let y = law.response([0.0, 0.15, 0.0]).unwrap();
        assert!(x.energy_j_m2 > 1.03 * y.energy_j_m2);
    }
    let mut m = model("hgo-dispersed");
    m.parameters.get_mut("k1").unwrap().value = 0.0;
    let c = source(vec![m], 0.0013, 0.24, ValidityDomain::unconstrained(), &[]);
    bind(&c)
        .membrane_law()
        .unwrap()
        .admit_linearization([0.0; 3])
        .unwrap();
}

#[test]
fn full_source_bundle_model_selection_and_fixed_context_survive_lowering() {
    let cp = ("specific_heat_capacity", 1400.0, Dims([2, 0, -2, -1, 0, 0]));
    let c = source(
        vec![model("neo-hookean-log"), model("mooney-rivlin-decoupled")],
        0.0013,
        0.24,
        ValidityDomain::unconstrained(),
        &[cp],
    );
    let requirements = [
        (AREAL_DENSITY, Dims([-2, 1, 0, 0, 0, 0])),
        (THICKNESS, Dims([1, 0, 0, 0, 0, 0])),
        (cp.0, cp.2),
    ]
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    });
    let original =
        resolve_material_state_point(&c, &point(), &requirements, Pick::SingleClaimOnly).unwrap();
    let sheet = ResolvedSheet::try_from_material(original.clone())
        .unwrap()
        .with_hyperelastic_membrane(&c, 1, controls())
        .unwrap();
    let law = sheet.membrane_law().unwrap();
    let MembraneLaw::Hyperelastic(bound) = &law else {
        panic!("wrong law")
    };
    assert_eq!(sheet.material(), &original);
    assert_eq!(bound.material(), &original);
    assert_eq!(bound.model().model_card(), &c.models()[1]);
    assert_eq!(bound.model().context(), &point());
    assert!(bound.material().property(cp.0).is_some());
    assert!(sheet.with_hyperelastic_membrane(&c, 0, controls()).is_err());

    let different = card("neo-hookean-log", 0.002, 0.24);
    assert!(
        reference(&c)
            .with_hyperelastic_membrane(&different, 0, controls())
            .is_err()
    );
    assert!(
        reference(&c)
            .with_hyperelastic_membrane(&c, 2, controls())
            .is_err()
    );
    let mut wrong_context = model("neo-hookean-log");
    wrong_context.validity = ValidityDomain::unconstrained().with("T", 320.0, 330.0);
    let wrong = source(
        vec![wrong_context],
        0.0013,
        0.24,
        ValidityDomain::unconstrained(),
        &[],
    );
    assert!(
        reference(&wrong)
            .with_hyperelastic_membrane(&wrong, 0, controls())
            .is_err()
    );
}

#[test]
fn reference_mechanical_conditioning_ambiguous_capabilities_and_bad_states_refuse() {
    let c = card("neo-hookean-log", 0.0013, 0.24);
    let law = bind(&c).membrane_law().unwrap();
    for e in [
        [f64::NAN, 0.0, 0.0],
        [-0.5, 0.0, 0.0],
        [0.0, -0.51, 0.0],
        [0.0, 0.0, 1.1],
        [20.0, 0.0, 0.0],
    ] {
        assert!(law.response(e).is_err());
    }
    let invalid = PlaneStressSettings {
        relative_tolerance: 0.0,
        ..controls()
    };
    assert!(
        reference(&c)
            .with_hyperelastic_membrane(&c, 0, invalid)
            .is_err()
    );
    let mut point_only = model("neo-hookean-log");
    point_only
        .parameters
        .get_mut("stretch_limit")
        .unwrap()
        .value = 1.0;
    let point_only = source(
        vec![point_only],
        0.0013,
        0.24,
        ValidityDomain::unconstrained(),
        &[],
    );
    let point_only = bind(&point_only).membrane_law().unwrap();
    assert!(point_only.response([0.0; 3]).is_ok());
    assert!(point_only.admit_linearization([0.0; 3]).is_err());
    let mut incomplete = model("neo-hookean-log");
    incomplete.parameters.remove("mu");
    let c = source(
        vec![incomplete],
        0.0013,
        0.24,
        ValidityDomain::unconstrained(),
        &[],
    );
    assert!(matches!(
        reference(&c).with_hyperelastic_membrane(&c, 0, controls()),
        Err(SheetError::Hyperelastic(_))
    ));

    let c = source(
        vec![model("neo-hookean-log")],
        0.0013,
        0.24,
        ValidityDomain::unconstrained().with("strain", -0.3, 0.3),
        &[],
    );
    let q = point().with("strain", 0.0).unwrap();
    let mechanical_reference = resolve_sheet(
        &c,
        &q,
        SheetCapabilities {
            membrane: false,
            bending: false,
        },
        Pick::SingleClaimOnly,
    )
    .unwrap();
    assert!(
        mechanical_reference
            .with_hyperelastic_membrane(&c, 0, controls())
            .is_err()
    );

    let extra = [("warp_stiffness", 5000.0, MEMBRANE_DIMS)];
    let c = source(
        vec![model("neo-hookean-log")],
        0.0013,
        0.24,
        ValidityDomain::unconstrained(),
        &extra,
    );
    let requirements = [
        (AREAL_DENSITY, Dims([-2, 1, 0, 0, 0, 0])),
        (THICKNESS, Dims([1, 0, 0, 0, 0, 0])),
        (extra[0].0, extra[0].2),
    ]
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    });
    let original =
        resolve_material_state_point(&c, &point(), &requirements, Pick::SingleClaimOnly).unwrap();
    assert!(
        ResolvedSheet::try_from_material(original)
            .unwrap()
            .with_hyperelastic_membrane(&c, 0, controls())
            .is_err()
    );
}
