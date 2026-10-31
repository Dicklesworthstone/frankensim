//! Directional finite expansion, original-mass mechanics and derivative oracles.
use super::*;
use fs_blake3::ContentHash;
use fs_couple::finite_solid::thermal::{
    ExpansionThermalFrame, ExpansionThermalReference, ThermalReference,
};
use fs_matdb::{ClaimSelection, SelectionPolicy};
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_material::thermal::transport::expansion::{
    FreeExpansionGeometryQuery, FreeExpansionStrainSources, PreparedFreeExpansion,
};

const PRINCIPAL: ContentHash = ContentHash([1; 32]);
const MESH: ContentHash = ContentHash([2; 32]);
const Q: [[f64; 3]; 3] = [[0.6, -0.8, 0.0], [0.8, 0.6, 0.0], [0.0, 0.0, 1.0]];
const STRETCHES: [f64; 3] = [1.10, 0.96, 1.04];

fn key(axis: usize) -> PropertyKey {
    PropertyKey::new(format!("engineering_length_change_{axis}"), Dims::NONE)
}

fn card(stretches: [f64; 3], source_axis: &str, stress_axis: bool) -> MaterialCard {
    let base = source_card(6.0, 1.0, false, false);
    let mut claims = base.claims().clone();
    // Nonzero fitted datum offsets ensure rebasing is a ratio of lengths,
    // not subtraction of engineering strains or reuse of a secant CTE.
    for (i, datum) in [0.2, -0.1, 0.05].into_iter().enumerate() {
        let change = (1.0 + datum) * (stretches[i] - 1.0);
        let mut validity = ValidityDomain::unconstrained().with_quantity(
            source_axis,
            QuantitySpec::dimensional(Temperature::DIMS),
            290.0,
            310.0,
        );
        if stress_axis {
            validity = validity.with_quantity(
                "stress",
                QuantitySpec::dimensional(Pressure::DIMS),
                0.0,
                100.0,
            );
        }
        claims
            .insert_claim(PropertyClaim {
                key: key(i),
                value: PropertyValue::Curve {
                    abscissa: source_axis.into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (290.0, datum - change),
                        (300.0, datum),
                        (310.0, datum + change),
                    ],
                    dims: Dims::NONE,
                },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: base.models()[0].provenance.clone(),
            })
            .unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn path(
    card: &MaterialCard,
    isotropic: bool,
    axis: &str,
    reference: QueryPoint,
) -> PreparedFreeExpansion {
    let lower = reference.clone().with(axis, 290.0).unwrap();
    let upper = reference.clone().with(axis, 310.0).unwrap();
    let keys = [key(0), key(1), key(2)];
    let strains = if isotropic {
        FreeExpansionStrainSources::Isotropic {
            strain: &keys[0],
            selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        }
    } else {
        FreeExpansionStrainSources::FixedPrincipal {
            frame: PRINCIPAL,
            strains: [&keys[0], &keys[1], &keys[2]],
            selections: [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 3],
        }
    };
    PreparedFreeExpansion::try_new(
        card,
        FreeExpansionGeometryQuery {
            reference: &reference,
            strains,
            temperature_axis: axis,
            lower: &lower,
            upper: &upper,
        },
    )
    .unwrap()
}

fn frame(q: [[f64; 3]; 3]) -> ExpansionThermalFrame {
    ExpansionThermalFrame::FixedPrincipal {
        principal_frame: PRINCIPAL,
        mesh_frame: MESH,
        principal_to_mesh: q,
    }
}

fn bind_expansion(card: &MaterialCard, q: [[f64; 3]; 3], temperature: f64) -> FiniteSolidMaterial {
    let thermal = ExpansionThermalReference::try_new(
        card,
        path(card, false, "T", point(300.0)),
        temperature,
        frame(q),
    )
    .unwrap();
    FiniteSolidMaterial::resolve_thermal(
        card,
        0,
        ThermalReference::Expansion(thermal),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn det(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

// Independent Ogden-alpha=2 energy through invariants, using F's projections
// onto principal directions rather than the production pullback tensors.
fn energy(f: &[f64; 9], s: [f64; 3], q: [[f64; 3]; 3]) -> f64 {
    let jth = s.iter().product::<f64>();
    let je = det(f) / jth;
    let mut norm2 = 0.0;
    for i in 0..3 {
        for a in 0..3 {
            let component = (0..3).map(|j| f[3 * i + j] * q[j][a]).sum::<f64>() / s[a];
            norm2 += component * component;
        }
    }
    jth * (7.5 * (norm2 / je.powf(2.0 / 3.0) - 3.0) + 20.0 * (je - 1.0).powi(2))
}

#[test]
fn g1_oriented_expansion_energy_piola_tangent_and_native_assembly_are_consistent() {
    with_cx(|_, cx| {
        let source = card(STRETCHES, "T", false);
        let material = bind_expansion(&source, Q, 310.0);
        let f = [1.09, 0.04, -0.02, 0.03, 0.97, 0.05, 0.01, -0.06, 1.08];
        let p = material.first_piola(&f).unwrap();
        let a = material.consistent_tangent(&f).unwrap();
        near(
            material.energy_density(&f),
            energy(&f, STRETCHES, Q),
            2.0e-13,
        );
        let h = 2.0e-6;
        for column in 0..9 {
            let mut plus = f;
            plus[column] += h;
            let mut minus = f;
            minus[column] -= h;
            near(
                (energy(&plus, STRETCHES, Q) - energy(&minus, STRETCHES, Q)) / (2.0 * h),
                p[column],
                3.0e-8,
            );
            let pp = material.first_piola(&plus).unwrap();
            let pm = material.first_piola(&minus).unwrap();
            for row in 0..9 {
                near((pp[row] - pm[row]) / (2.0 * h), a[row][column], 3.0e-7);
                near(a[row][column], a[column][row], 3.0e-12);
            }
        }
        let u = NODES.map(|x| {
            std::array::from_fn::<_, 3, _>(|i| {
                (0..3).map(|j| f[3 * i + j] * x[j]).sum::<f64>() - x[i]
            })
        });
        let flat: Vec<_> = u.iter().flatten().copied().collect();
        let problem = problem(&material, &[], &[]);
        near(
            problem.stored_energy_j(&flat, cx).unwrap(),
            VOLUME * energy(&f, STRETCHES, Q),
            1.0e-13,
        );
        let force = problem.internal_force(&flat, cx).unwrap();
        let stiffness = problem.tangent_matrix(&flat, cx).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                near(force[3 * (j + 1) + i], VOLUME * p[3 * i + j], 2.0e-13);
                for k in 0..3 {
                    for l in 0..3 {
                        near(
                            stiffness.get(3 * (j + 1) + i, 3 * (l + 1) + k),
                            VOLUME * a[3 * i + j][3 * k + l],
                            2.0e-12,
                        );
                    }
                }
            }
        }
    });
}

#[test]
fn g1_directional_free_expansion_rotation_cooling_and_restraint_keep_original_mass() {
    with_cx(|_, cx| {
        let source = card(STRETCHES, "T", false);
        for temperature in [290.0, 300.0, 310.0] {
            let material = bind_expansion(&source, Q, temperature);
            let ThermalReference::Expansion(reference) = material.thermal_reference().unwrap()
            else {
                panic!("engineering length source required")
            };
            let fth = *reference.thermal_gradient();
            let s = reference.geometry().principal_stretches();
            for axis in 0..3 {
                near(
                    s[axis],
                    1.0 + (temperature - 300.0) / 10.0 * (STRETCHES[axis] - 1.0),
                    3.0e-16,
                );
            }
            let rotation = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
            let f: [f64; 9] = std::array::from_fn(|index| {
                (0..3)
                    .map(|k| rotation[index / 3][k] * fth[3 * k + index % 3])
                    .sum()
            });
            for stress in material.first_piola(&f).unwrap() {
                near(stress, 0.0, 1.0e-12);
            }
            let u = NODES.map(|x| {
                std::array::from_fn(|i| (0..3).map(|j| f[3 * i + j] * x[j]).sum::<f64>() - x[i])
            });
            let body = FiniteSolidBody::new(
                problem(&material, &[], &[]),
                &u,
                &[[0.0; 3]; 4],
                0.0,
                Default::default(),
                cx,
            )
            .unwrap();
            let (next, receipt) = body
                .step_candidate(0.01, &[], [0.0; 3], controls(1.0e-8), cx)
                .unwrap();
            for (actual, expected) in next
                .displacement_m()
                .iter()
                .flatten()
                .zip(u.iter().flatten())
            {
                near(*actual, *expected, 2.0e-10);
            }
            near(next.reference_mass_kg(), 1.0, 2.0e-15);
            near(receipt.strain_after_j, 0.0, 2.0e-12);
            near(reference.geometry().relative_volume(), det(&f), 1.0e-15);
            assert_eq!(material.density().query(), point(300.0));
            assert_eq!(material.law().context(), &point(temperature));
            let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
            if temperature != 300.0 {
                assert!(material.energy_density(&identity) > 0.0);
                assert!(
                    material
                        .first_piola(&identity)
                        .unwrap()
                        .iter()
                        .any(|p| p.abs() > 0.1)
                );
            }
        }
    });
}

#[test]
fn g3_equal_principal_lengths_reduce_to_original_isotropic_mean_reference() {
    let source = card([1.1; 3], "T", false);
    let directional = bind_expansion(&source, Q, 310.0);
    let isotropic = ExpansionThermalReference::try_new(
        &source,
        path(&source, true, "T", point(300.0)),
        310.0,
        ExpansionThermalFrame::Isotropic,
    )
    .unwrap();
    let isotropic = FiniteSolidMaterial::resolve_thermal(
        &source,
        0,
        ThermalReference::Expansion(isotropic),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let mean_source = super::thermal::thermal_card(0.01, 1.0, false);
    let mean = super::thermal::prepared(&mean_source);
    let f = [1.2, 0.1, 0.0, 0.0, 0.9, -0.04, 0.01, 0.0, 1.05];
    for material in [&directional, &isotropic] {
        near(
            material.energy_density(&f),
            mean.energy_density(&f),
            5.0e-13,
        );
        for (actual, expected) in material
            .first_piola(&f)
            .unwrap()
            .iter()
            .zip(mean.first_piola(&f).unwrap())
        {
            near(*actual, expected, 3.0e-12);
        }
        for (actual, expected) in material
            .consistent_tangent(&f)
            .unwrap()
            .iter()
            .flatten()
            .zip(mean.consistent_tangent(&f).unwrap().iter().flatten())
        {
            near(*actual, *expected, 3.0e-11);
        }
    }
}

#[test]
fn g4_finite_expansion_refuses_wrong_sources_frames_dynamic_axes_and_ambiguous_temperature() {
    let source = card(STRETCHES, "T", false);
    let original = path(&source, false, "T", point(300.0));
    let prepare =
        |frame, t| ExpansionThermalReference::try_new(&source, original.clone(), t, frame);
    assert!(prepare(ExpansionThermalFrame::Isotropic, 310.0).is_err());
    for q in [
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
        [[1.0, 0.1, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[f64::NAN; 3]; 3],
    ] {
        assert!(prepare(frame(q), 310.0).is_err());
    }
    for (principal_frame, mesh_frame) in [
        (MESH, MESH),
        (PRINCIPAL, PRINCIPAL),
        (PRINCIPAL, ContentHash([0; 32])),
    ] {
        assert!(
            prepare(
                ExpansionThermalFrame::FixedPrincipal {
                    principal_frame,
                    mesh_frame,
                    principal_to_mesh: Q,
                },
                310.0
            )
            .is_err()
        );
    }
    for t in [280.0, 320.0, f64::NAN] {
        assert!(prepare(frame(Q), t).is_err());
    }
    let other = card([1.2; 3], "T", false);
    assert!(ExpansionThermalReference::try_new(&other, original.clone(), 310.0, frame(Q)).is_err());
    let reference = prepare(frame(Q), 310.0).unwrap();
    assert!(
        FiniteSolidMaterial::resolve_thermal(
            &other,
            0,
            ThermalReference::Expansion(reference),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let stressed = card(STRETCHES, "T", true);
    let query = point(300.0)
        .with_quantity("stress", QuantitySpec::dimensional(Pressure::DIMS), 10.0)
        .unwrap();
    assert!(
        ExpansionThermalReference::try_new(
            &stressed,
            path(&stressed, false, "T", query),
            310.0,
            frame(Q)
        )
        .is_err()
    );
    let custom = card(STRETCHES, "theta", false);
    let query = QueryPoint::new()
        .with_quantity("theta", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap();
    let clean = path(&custom, false, "theta", query.clone());
    let reference = ExpansionThermalReference::try_new(&custom, clean, 310.0, frame(Q)).unwrap();
    assert_eq!(reference.reference_context(), &point(300.0));
    assert_eq!(reference.current_context(), &point(310.0));
    let ambiguous = query
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 305.0)
        .unwrap();
    assert!(
        ExpansionThermalReference::try_new(
            &custom,
            path(&custom, false, "theta", ambiguous),
            310.0,
            frame(Q)
        )
        .is_err()
    );
}
