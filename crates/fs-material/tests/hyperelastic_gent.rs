//! G0/G1/G3 checks of the Gent constitutive model, not empirical calibration.
//! Gent (1996), doi:10.5254/1.3538357; Horgan & Murphy (2022),
//! doi:10.1007/s10659-022-09908-1, Eq. (6.11). Jm is I1_max - 3.
//! The independent shear and uniaxial formulas below use incompressibility;
//! the finite-bulk checks retain the explicitly declared decoupled extension.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::MaterialError;
use fs_material::graph::LawRegistry;
use fs_material::hyper::cards::{ResolvedHyperelastic, register_hyperelastic_laws};
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::{Hyperelastic, HyperelasticModel, Volumetric};
use fs_material::sheet::{
    AREAL_DENSITY, MembraneLaw, SheetCapabilities, SheetLaw, THICKNESS, resolve_sheet,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::tensor::{matmul3, rotation};
use fs_qty::Dims;

const MU: f64 = 2.4e5;
const JM: f64 = 6.0;
const KAPPA: f64 = 4.0e6;
const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const IDENTITY: [f64; 9] = [1., 0., 0., 0., 1., 0., 0., 0., 1.];

fn gent(jm: f64, volumetric: Volumetric) -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::Gent {
            mu: MU,
            locking_invariant: jm,
            volumetric,
        },
        8.0,
    )
    .unwrap()
}
fn quadratic(jm: f64) -> Hyperelastic {
    gent(jm, Volumetric::Quadratic { kappa: KAPPA })
}
fn shear(amount: f64) -> [f64; 9] {
    let mut f = IDENTITY;
    f[1] = amount;
    f
}
fn determinant(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}
fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn g1_simple_shear_resolves_small_energy_and_chain_stiffening() {
    let law = quadratic(JM);
    for gamma in [1e-12_f64, 1e-7, 0.001, 0.8, 1.8, 2.35] {
        let x = gamma * gamma / JM;
        let energy = -0.5 * MU * JM * (-x).ln_1p();
        let stress = MU * gamma / (1.0 - x);
        let stiffness = MU * (1.0 + x) / (1.0 - x).powi(2);
        let f = shear(gamma);
        close(law.energy(&f), energy, energy, 2e-12);
        close(law.piola(&f).unwrap()[1], stress, stress, 3e-12);
        close(law.tangent(&f).unwrap()[1][1], stiffness, stiffness, 4e-12);
    }
    assert_eq!(law.energy(&IDENTITY), 0.0);
    assert_eq!(law.piola(&IDENTITY).unwrap(), [0.0; 9]);
    assert!(law.tangent(&shear(2.35)).unwrap()[1][1] > 100.0 * MU);
}

#[test]
fn g1_exactly_incompressible_uniaxial_and_large_chain_limit() {
    let law = quadratic(JM);
    for stretch in [0.5_f64, 1.0, 1.3, 2.0, 2.7] {
        let lateral = 1.0 / stretch.sqrt();
        let response = law
            .incompressible_plane_stress(
                &[stretch, 0.0, 0.0, lateral],
                PlaneStressSettings::default(),
            )
            .unwrap();
        let shift = stretch * stretch + 2.0 / stretch - 3.0;
        let expected = MU * (stretch - 1.0 / (stretch * stretch)) / (1.0 - shift / JM);
        close(response.piola[0], expected, MU.max(expected.abs()), 2e-11);
        close(response.piola[3], 0.0, MU, 2e-11);
        close(determinant(&response.deformation_gradient), 1.0, 1.0, 2e-15);
        close(
            response.energy_density_pa,
            -0.5 * MU * JM * (-shift / JM).ln_1p(),
            MU,
            2e-11,
        );
    }
    // mu*Jm would overflow, but the limiting energy and its derivatives do not.
    let large_chain = quadratic(1e308);
    let neo_hookean = Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 0.5 * MU,
            c01: 0.0,
            kappa: KAPPA,
        },
        8.0,
    )
    .unwrap();
    for f in [
        shear(1e-12),
        shear(0.8),
        [1.2, 0.1, 0., 0., 0.9, 0.04, 0., 0., 1.1],
    ] {
        close(large_chain.energy(&f), neo_hookean.energy(&f), MU, 2e-12);
        for (a, b) in large_chain
            .piola(&f)
            .unwrap()
            .into_iter()
            .zip(neo_hookean.piola(&f).unwrap())
        {
            close(a, b, KAPPA, 2e-12);
        }
        for (a, b) in large_chain
            .tangent(&f)
            .unwrap()
            .into_iter()
            .flatten()
            .zip(neo_hookean.tangent(&f).unwrap().into_iter().flatten())
        {
            close(a, b, KAPPA, 2e-12);
        }
    }
}

#[test]
fn g3_objectivity_and_consistent_tangent_for_both_bulk_laws() {
    let axis = [
        1.0 / 14.0_f64.sqrt(),
        2.0 / 14.0_f64.sqrt(),
        3.0 / 14.0_f64.sqrt(),
    ];
    let q = rotation(axis, 0.62);
    let f = [1.24, 0.18, -0.03, 0.01, 0.94, 0.05, 0.04, -0.02, 1.04];
    let matrix =
        |values: [f64; 9]| core::array::from_fn(|i| core::array::from_fn(|j| values[3 * i + j]));
    let qf = matmul3(&q, &matrix(f));
    let rotated = core::array::from_fn(|i| qf[i / 3][i % 3]);
    for volume in [
        Volumetric::Quadratic { kappa: KAPPA },
        Volumetric::LogSquared { kappa: KAPPA },
    ] {
        let law = gent(JM, volume);
        close(law.energy(&rotated), law.energy(&f), MU, 2e-12);
        let qp = matmul3(&q, &matrix(law.piola(&f).unwrap()));
        for (i, value) in law.piola(&rotated).unwrap().iter().enumerate() {
            close(*value, qp[i / 3][i % 3], KAPPA, 2e-12);
        }
        let tangent = law.tangent(&f).unwrap();
        for j in 0..9 {
            let (mut plus, mut minus) = (f, f);
            plus[j] += 1e-6;
            minus[j] -= 1e-6;
            let (p, m) = (law.piola(&plus).unwrap(), law.piola(&minus).unwrap());
            for i in 0..9 {
                close(tangent[i][j], (p[i] - m[i]) / 2e-6, KAPPA, 3e-8);
                close(tangent[i][j], tangent[j][i], KAPPA, 3e-14);
            }
        }
    }
}

#[test]
fn g0_strict_chain_limit_rejects_derivatives_and_native_reductions() {
    let law = quadratic(4.0);
    for f in [shear(2.0), shear(2.1)] {
        assert!(
            matches!(law.check_state(&f), Err(MaterialError::State { what }) if what.contains("Gent chain limit"))
        );
        assert!(law.piola(&f).is_err());
        assert!(law.tangent(&f).is_err());
        assert!(
            law.mixed_state(&f, 0.0, VolumetricTreatment::ExactIncompressible)
                .is_err()
        );
        assert!(
            law.incompressible_plane_stress(&[1.0, f[1], 0.0, 1.0], PlaneStressSettings::default())
                .is_err()
        );
    }
    assert!(law.piola(&shear(1.99)).is_ok());
    assert!(law.card().discrepancy_rel.is_infinite());
    assert!(law.card().calibration.is_none());
    for limit in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        let mut invalid = law.clone();
        if let HyperelasticModel::Gent {
            locking_invariant, ..
        } = &mut invalid.model
        {
            *locking_invariant = limit;
        }
        assert!(Hyperelastic::new(invalid.model.clone(), 8.0).is_err());
        assert!(matches!(
            invalid.check_state(&IDENTITY),
            Err(MaterialError::Parameters { .. })
        ));
    }
}

#[test]
fn g1_native_mixed_condensation_and_compressible_plane_stress() {
    let f = [1.2, 0.1, 0., 0., 0.9, 0.04, 0., 0., 1.1];
    let j = determinant(&f);
    for volume in [
        Volumetric::Quadratic { kappa: KAPPA },
        Volumetric::LogSquared { kappa: KAPPA },
    ] {
        let law = gent(JM, volume);
        let pressure = match volume {
            Volumetric::Quadratic { .. } => KAPPA * (j - 1.0),
            Volumetric::LogSquared { .. } => KAPPA * j.ln() / j,
        };
        law.validate_mixed().unwrap();
        let response = law
            .mixed_response(&f, pressure, VolumetricTreatment::OriginalMaterial)
            .unwrap();
        close(response.state.pressure_residual, 0.0, 1.0, 1e-14);
        close(response.state.potential_pa, law.energy(&f), MU, 2e-12);
        let tangent = law.tangent(&f).unwrap();
        let piola = law.piola(&f).unwrap();
        for a in 0..9 {
            close(response.state.first_piola_pa[a], piola[a], KAPPA, 2e-12);
            for b in 0..9 {
                let condensed = response.tangent_pa[a][b]
                    - response.pressure_coupling[a] * response.pressure_coupling[b]
                        / response.pressure_tangent_per_pa;
                close(condensed, tangent[a][b], KAPPA, 2e-12);
            }
        }
        let a = [1.45, 0.08, 0.0, 0.92];
        let controls = PlaneStressSettings::default();
        let plane = law.plane_stress(&a, controls).unwrap();
        let stress = law.piola(&plane.deformation_gradient).unwrap();
        for i in [2, 5, 8] {
            close(stress[i], 0.0, MU, 3e-10);
        }
        for column in 0..4 {
            let (mut plus, mut minus) = (a, a);
            plus[column] += 1e-6;
            minus[column] -= 1e-6;
            let p = law.plane_stress_state(&plus, controls).unwrap();
            let m = law.plane_stress_state(&minus, controls).unwrap();
            for row in 0..4 {
                close(
                    plane.tangent[row][column],
                    (p.piola[row] - m.piola[row]) / 2e-6,
                    KAPPA,
                    1e-7,
                );
            }
        }
    }
}

fn provenance() -> Provenance {
    Provenance {
        source: "synthetic Gent verification, not measured rubber".into(),
        license: "test-only".into(),
        artifact: None,
    }
}
fn model_card() -> ConstitutiveModelCard {
    ConstitutiveModelCard {
        law: LawId("gent-decoupled".into()),
        law_version: 1,
        parameters: [
            ("mu", MU, PA),
            ("locking_invariant", JM, Dims::NONE),
            ("kappa", KAPPA, PA),
            ("volumetric_kind", 0.0, Dims::NONE),
            ("stretch_limit", 8.0, Dims::NONE),
        ]
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

#[test]
fn g1_original_source_card_reaches_graph_and_finite_strain_sheet() {
    let model = model_card();
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    let resolved = ResolvedHyperelastic::try_from_model_card_at(&model, &point).unwrap();
    assert_eq!(resolved.model_card(), &model);
    assert_eq!(resolved.context(), &point);
    assert!(ResolvedHyperelastic::try_from_model_card(&model).is_err());
    for required in ["locking_invariant", "volumetric_kind"] {
        let mut missing = model.clone();
        missing.parameters.remove(required);
        assert!(Hyperelastic::from_model_card(&missing).is_err());
    }
    let mut mistyped = model.clone();
    mistyped
        .parameters
        .get_mut("locking_invariant")
        .unwrap()
        .dims = PA;
    assert!(Hyperelastic::from_model_card(&mistyped).is_err());
    let mut registry = LawRegistry::new();
    register_hyperelastic_laws(&mut registry);
    let node = registry.instantiate("solid/gent", &model).unwrap();
    let mut h = [0.0; 9];
    h[1] = 0.8;
    close(
        node.evaluate(&[], &h).unwrap().outputs[1],
        MU * 0.8 / (1.0 - 0.64 / JM),
        MU,
        2e-12,
    );
    h[1] = 3.0;
    assert!(node.evaluate(&[], &h).is_err());
    assert!(node.tangent(&[], &h).is_none());

    let thickness = 0.002;
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        (AREAL_DENSITY, 2.2, Dims([-2, 1, 0, 0, 0, 0])),
        (THICKNESS, thickness, Dims([1, 0, 0, 0, 0, 0])),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic Gent".into(),
            phase: "solid".into(),
            process: "constitutive test".into(),
            revision: 0,
        },
        claims,
        vec![model.clone()],
    )
    .unwrap();
    let sheet = resolve_sheet(
        &card,
        &point,
        SheetCapabilities {
            membrane: false,
            bending: false,
        },
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
    .with_incompressible_hyperelastic_membrane(&card, 0, PlaneStressSettings::default())
    .unwrap();
    let MembraneLaw::Hyperelastic(law) = sheet.membrane_law().unwrap() else {
        panic!("native hyperelastic law required");
    };
    assert_eq!(law.model().model_card(), &model);
    assert_eq!(law.material().card_identity(), card.content_hash());
    let stretch = 2.0;
    let shift = stretch * stretch + 2.0 / stretch - 3.0;
    let response = law
        .response([
            0.5 * (stretch * stretch - 1.0),
            0.5 * (1.0 / stretch - 1.0),
            0.0,
        ])
        .unwrap();
    let nominal = MU * (stretch - 1.0 / (stretch * stretch)) / (1.0 - shift / JM);
    close(
        response.resultants_n_m[0],
        thickness * nominal / stretch,
        thickness * MU,
        2e-11,
    );
    close(response.resultants_n_m[1], 0.0, thickness * MU, 2e-11);
    close(
        response.energy_j_m2,
        -0.5 * thickness * MU * JM * (-shift / JM).ln_1p(),
        thickness * MU,
        2e-11,
    );
    assert!(law.response([4.0, 0.0, 0.0]).is_err());
}
