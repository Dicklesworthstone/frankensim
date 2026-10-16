//! G0/G3 numerical adapter checks. The authored Sitka card transcribes the
//! committed FPL-GTR-282 seed, not a new specimen or measured uncertainty.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{plane_stress_viscosity::MaterialPlane, state_point::*};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

fn rows() -> Vec<(String, Dims, f64)> {
    [
        ("young_modulus_longitudinal", Pressure::DIMS, 1.188e10),
        ("er_over_el", Dims::NONE, 0.078),
        ("et_over_el", Dims::NONE, 0.043),
        ("glr_over_el", Dims::NONE, 0.064),
        ("grt_over_el", Dims::NONE, 0.003),
        ("glt_over_el", Dims::NONE, 0.061),
        ("nu_lr", Dims::NONE, 0.372),
        ("nu_lt", Dims::NONE, 0.467),
        ("nu_rt", Dims::NONE, 0.435),
        ("nu_rl", Dims::NONE, 0.04),
        ("nu_tl", Dims::NONE, 0.025),
        ("nu_tr", Dims::NONE, 0.245),
        ("density", Density::DIMS, 448.0),
    ]
    .map(|(n, d, v)| (n.to_owned(), d, v))
    .to_vec()
}
fn card(rows: Vec<(String, Dims, f64)>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in rows {
        claims.insert_claim(PropertyClaim {
            key:PropertyKey::new(name,dims),value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained().with("moisture_content",12.0,12.0),
            uncertainty:UncertaintyModel::Unstated,
            interpolation:InterpolationPolicy::ConstantWithinValidity,observations:vec![],
            provenance:Provenance {
                source:"test-fixture: Sitka FPL-GTR-282 ch5 Tables5-1/5-2/5-3a transcribed seed rows".into(),
                license:"Work-of-the-US-Government-Public-Use-Permitted".into(),artifact:None,
            },
        }).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "test-fixture: Sitka".into(),
            phase: "clear wood".into(),
            process: "authored source-table transcription".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn point(mc: f64) -> QueryPoint {
    QueryPoint::new().with("moisture_content", mc).unwrap()
}
fn resolve(
    card: &MaterialCard,
    set: OrthotropicRequirementSet,
) -> Result<ReferenceRatioOrthotropicStatePoint, MaterialStatePointError> {
    resolve_orthotropic_elastic_state_point_from_reference_ratios(
        card,
        &point(12.0),
        MaterialPropertySelection::SingleClaimOnly,
        0.01,
        set,
        PoissonReciprocityRule::MajorRatios,
    )
}
fn full(card: &MaterialCard) -> ReferenceRatioOrthotropicStatePoint {
    resolve(card, OrthotropicRequirementSet::Full3d).unwrap()
}
fn relative(a: f64, b: f64) -> f64 {
    (a / b - 1.0).abs()
}

// Independent pivoted Gauss-Jordan inversion; production stiffness uses
// three-dimensional compliance cofactors and tensor-shear convention.
fn inverse(mut a: [[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut b = [[0.0; 6]; 6];
    for (i, row) in b.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for c in 0..6 {
        let pivot = (c..6)
            .max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs()))
            .unwrap();
        a.swap(c, pivot);
        b.swap(c, pivot);
        let d = a[c][c];
        assert!(d.is_finite() && d != 0.0);
        for j in 0..6 {
            a[c][j] /= d;
            b[c][j] /= d;
        }
        let pivot_a = a[c];
        let pivot_b = b[c];
        for i in 0..6 {
            if i == c {
                continue;
            }
            let f = a[i][c];
            for j in 0..6 {
                a[i][j] -= f * pivot_a[j];
                b[i][j] -= f * pivot_b[j];
            }
        }
    }
    b
}
#[test]
fn g0_round_trip_compliance_recovers_pack_ratios() {
    let state = full(&card(rows()));
    let law = state.law().unwrap();
    let c = law.stiffness();
    let s = inverse(c);
    let recovered = [
        s[0][0] / s[1][1],
        s[0][0] / s[2][2],
        s[0][0] / (2.0 * s[3][3]),
        s[0][0] / (2.0 * s[4][4]),
        s[0][0] / (2.0 * s[5][5]),
        -s[0][1] / s[0][0],
        -s[0][2] / s[0][0],
        -s[1][2] / s[1][1],
    ];
    for (got, want) in recovered
        .into_iter()
        .zip([0.078, 0.043, 0.064, 0.003, 0.061, 0.372, 0.467, 0.435])
    {
        assert!(relative(got, want) < 2e-15, "{got} vs {want}");
        assert!(
            got.to_bits().abs_diff(want.to_bits()) <= 8,
            "{got} vs {want}"
        );
    }
    for i in 0..6 {
        for j in 0..6 {
            let product = (0..6).map(|k| c[i][k] * s[k][j]).sum::<f64>();
            assert!((product - if i == j { 1.0 } else { 0.0 }).abs() < 1e-14);
        }
    }
}
#[test]
fn g0_products_are_bitwise_receipted_values() {
    let source = card(rows());
    let state = full(&source);
    let law = state.law().unwrap();
    let v = |n| state.resolved().property(n).unwrap().value_si();
    let el = v(REFERENCE_YOUNG_MODULUS_PROPERTY);
    assert_eq!(
        law.e.map(f64::to_bits),
        [el, el * v("er_over_el"), el * v("et_over_el")].map(f64::to_bits)
    );
    assert_eq!(
        law.g.map(f64::to_bits),
        REFERENCE_SHEAR_RATIO_PROPERTIES.map(|p| (el * v(p)).to_bits())
    );
    assert_eq!(
        law.nu.map(f64::to_bits),
        REFERENCE_MAJOR_POISSON_PROPERTIES.map(|p| v(p).to_bits())
    );
    assert_eq!(state.resolved().properties().len(), 13);
    for p in state.resolved().properties() {
        source.claims().verify_receipt(&p.answer().receipt).unwrap();
    }
    let again = full(&source);
    assert_eq!(state, again);
    assert_eq!(
        state.derivation().unwrap().identity,
        again.derivation().unwrap().identity
    );
    let mut reversed = rows();
    reversed.reverse();
    assert_eq!(
        state.derivation().unwrap().identity,
        full(&card(reversed)).derivation().unwrap().identity
    );
}
#[test]
fn g0_missing_shear_ratio_refuses_full3d_but_plane12_resolves() {
    let source = card(
        rows()
            .into_iter()
            .filter(|(n, _, _)| n != "grt_over_el")
            .collect(),
    );
    assert!(matches!(resolve(&source,OrthotropicRequirementSet::Full3d),
        Err(MaterialStatePointError::Query {property,..}) if property=="grt_over_el"));
    let plane = resolve(
        &source,
        OrthotropicRequirementSet::Plane(MaterialPlane::Axes12),
    )
    .unwrap();
    assert!(plane.law().is_none());
    let ReferenceRatioOrthotropicStatePoint::Plane(state) = plane else {
        panic!("plane law expected")
    };
    assert_eq!(state.law().e, [1.188e10, 1.188e10 * 0.078]);
    assert_eq!(state.law().g, 1.188e10 * 0.064);
    let minimal = card(
        rows()
            .into_iter()
            .filter(|(n, _, _)| {
                [
                    "young_modulus_longitudinal",
                    "er_over_el",
                    "nu_lr",
                    "glr_over_el",
                    "density",
                ]
                .contains(&n.as_str())
            })
            .collect(),
    );
    assert!(
        resolve(
            &minimal,
            OrthotropicRequirementSet::Plane(MaterialPlane::Axes12)
        )
        .is_ok()
    );
}
#[test]
fn g0_missing_density_refuses_by_name() {
    let source = card(
        rows()
            .into_iter()
            .filter(|(n, _, _)| n != "density")
            .collect(),
    );
    for set in [
        OrthotropicRequirementSet::Full3d,
        OrthotropicRequirementSet::Plane(MaterialPlane::Axes12),
    ] {
        assert!(
            matches!(resolve(&source,set),Err(MaterialStatePointError::Query {property,..}) if property=="density")
        );
    }
}
#[test]
fn g0_moisture_outside_validity_refuses() {
    let source = card(rows());
    for mc in [11.0, 13.0] {
        assert!(matches!(
            resolve_orthotropic_elastic_state_point_from_reference_ratios(
                &source,
                &point(mc),
                MaterialPropertySelection::SingleClaimOnly,
                0.01,
                OrthotropicRequirementSet::Full3d,
                PoissonReciprocityRule::MajorRatios
            ),
            Err(MaterialStatePointError::Query { .. })
        ));
    }
    assert!(resolve(&source, OrthotropicRequirementSet::Full3d).is_ok());
    assert!(
        resolve_orthotropic_elastic_state_point_from_reference_ratios(
            &source,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            0.01,
            OrthotropicRequirementSet::Full3d,
            PoissonReciprocityRule::MajorRatios
        )
        .is_err()
    );
}
#[test]
fn g0_reciprocity_residuals_recorded_not_enforced() {
    let state = full(&card(rows()));
    // Exact hand expressions from the printed rows, not the rounded
    // descriptions 0.379,0.245,0.022 (which cannot be 1e-12 oracles).
    let expected = [
        0.04 / (0.372 * 0.078) - 1.0,
        0.025 / (0.467 * 0.043) - 1.0,
        0.245 / (0.435 * 0.043 / 0.078) - 1.0,
    ];
    for (got, want) in state
        .derivation()
        .unwrap()
        .reciprocity_residuals
        .into_iter()
        .zip(expected)
    {
        assert!((got.unwrap() - want).abs() < 1e-12);
    }
    let absent = card(
        rows()
            .into_iter()
            .filter(|(n, _, _)| !REFERENCE_MINOR_POISSON_PROPERTIES.contains(&n.as_str()))
            .collect(),
    );
    assert_eq!(
        full(&absent).derivation().unwrap().reciprocity_residuals,
        [None; 3]
    );
    assert_eq!(state.law().unwrap().nu, [0.372, 0.467, 0.435]);
}
#[test]
fn g0_indefinite_compliance_refuses() {
    let mut data = rows();
    data.iter_mut().find(|(n, _, _)| n == "nu_rt").unwrap().2 = 8.0;
    assert!(matches!(
        resolve(&card(data), OrthotropicRequirementSet::Full3d),
        Err(MaterialStatePointError::InvalidDerived {
            quantity: "orthotropic_compliance"
        })
    ));
}
#[test]
fn g0_axis_convention_file_matches_adapter() {
    let text = include_str!("../../../data/matdb/seed-v1/instrument-axis-convention.tsv");
    let axis_rows: Vec<_> = text
        .lines()
        .filter(|l| l.starts_with("axis\t"))
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .collect();
    assert_eq!(
        REFERENCE_RATIO_AXES,
        [WoodAxis::L, WoodAxis::R, WoodAxis::T]
    );
    for (row, (label, axis)) in axis_rows.iter().zip([("L", 1), ("R", 2), ("T", 3)]) {
        assert_eq!(row[1], label);
        assert!(row[3].contains(&format!("axis {axis}")));
    }
    assert_eq!(axis_rows.len(), 3);
    assert_eq!(
        REFERENCE_SHEAR_RATIO_PROPERTIES,
        ["glr_over_el", "grt_over_el", "glt_over_el"]
    );
    for name in REFERENCE_SHEAR_RATIO_PROPERTIES {
        assert!(
            text.lines()
                .any(|l| l.starts_with(&format!("ratio_property\t{name}\t")))
        );
    }
}
#[test]
fn g0_requirement_caps_and_invalid_inputs() {
    for length in [128, 129] {
        assert_eq!(
            ScalarPropertyRequirement::try_new(
                "x".repeat(length),
                Dims::NONE,
                ScalarAdmissibility::Finite
            )
            .is_ok(),
            length == 128
        );
    }
    let many = card(
        (0..64)
            .map(|i| (format!("property_{i:02}"), Dims::NONE, 1.0))
            .collect(),
    );
    let requirements: Vec<_> = (0..65)
        .map(|i| {
            ScalarPropertyRequirement::try_new(
                format!("property_{i:02}"),
                Dims::NONE,
                ScalarAdmissibility::Finite,
            )
            .unwrap()
        })
        .collect();
    assert!(
        resolve_material_state_point(
            &many,
            &point(12.0),
            &requirements[..64],
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_ok()
    );
    assert!(matches!(
        resolve_material_state_point(
            &many,
            &point(12.0),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly
        ),
        Err(MaterialStatePointError::RequirementCount {
            observed: 65,
            maximum: 64
        })
    ));
    for strain in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            resolve_orthotropic_elastic_state_point_from_reference_ratios(
                &card(rows()),
                &point(12.0),
                MaterialPropertySelection::SingleClaimOnly,
                strain,
                OrthotropicRequirementSet::Full3d,
                PoissonReciprocityRule::MajorRatios
            )
            .is_err()
        );
    }
    let mut wrong = rows();
    wrong
        .iter_mut()
        .find(|(n, _, _)| n == "er_over_el")
        .unwrap()
        .1 = Pressure::DIMS;
    assert!(matches!(
        resolve(&card(wrong), OrthotropicRequirementSet::Full3d),
        Err(MaterialStatePointError::DimensionMismatch { .. })
    ));
}
#[test]
fn g0_ordered_planes_retain_only_their_constants() {
    let source = card(rows());
    let mut identities = vec![];
    for (plane, e, nu, g) in [
        (
            MaterialPlane::Axes12,
            [1.188e10, 1.188e10 * 0.078],
            0.372,
            1.188e10 * 0.064,
        ),
        (
            MaterialPlane::Axes21,
            [1.188e10 * 0.078, 1.188e10],
            0.372 * (1.188e10 * 0.078) / 1.188e10,
            1.188e10 * 0.064,
        ),
        (
            MaterialPlane::Axes13,
            [1.188e10, 1.188e10 * 0.043],
            0.467,
            1.188e10 * 0.061,
        ),
        (
            MaterialPlane::Axes23,
            [1.188e10 * 0.078, 1.188e10 * 0.043],
            0.435,
            1.188e10 * 0.003,
        ),
    ] {
        let state = resolve(&source, OrthotropicRequirementSet::Plane(plane)).unwrap();
        identities.push(state.derivation().unwrap().identity);
        let ReferenceRatioOrthotropicStatePoint::Plane(state) = state else {
            panic!("plane")
        };
        assert_eq!(state.law().e, e);
        assert_eq!(state.law().nu, nu);
        assert_eq!(state.law().g, g);
        let q = state.law().stiffness();
        assert_eq!(q[1], q[3]);
        assert_eq!(q[8], g);
    }
    identities.sort();
    identities.dedup();
    assert_eq!(identities.len(), 4);
}
