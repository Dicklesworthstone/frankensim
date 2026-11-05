//! Kinetic history -> original material sources -> native continuum mechanics.
//! Manufactured constitutive data test the integration, not a calibrated polymer.
//! Frozen stress snapshots do not supply thermomechanical parameter-work feedback.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::IsotropicElastic;
use fs_material::phase::crystallization::CrystallizationCaloricPair;
use fs_material::phase::crystallization::mechanics::{
    CrystallizationLock, CrystallizationStrainLimits,
};
use fs_material::phase::crystallization::properties::{
    CrystallizationMaterialPoint, CrystallizationPropertySource,
};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase};
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::state_point::{
    IsotropicElasticStatePoint, MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_isotropic_elastic_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::{
    Formulation, IncrementSettings, Mesh2, PlaneAssemblyBudget, PlaneKind, SmallStrainLawField,
    SmallStrainPlaneProblem,
};
use fs_thermal::crystallization::CrystallizationSolve;

fn pressure_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}

fn card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, axis, knots) in [
        (
            "young_modulus",
            Pressure::DIMS,
            "Xm",
            vec![(0.0, 1e9), (0.75, 2.5e9)],
        ),
        // A separately supported fit exists beyond the first source's domain.
        // Continuing a pinned history must not silently switch to this claim.
        (
            "young_modulus",
            Pressure::DIMS,
            "Xm",
            vec![(0.76, 2.6e9), (0.8, 2.7e9)],
        ),
        (
            "density",
            Density::DIMS,
            "Xv",
            vec![(0.0, 1000.0), (1.0, 1030.0)],
        ),
        (
            "poisson_ratio",
            Dims::NONE,
            "Xm",
            vec![(0.0, 0.3), (1.0, 0.3)],
        ),
    ] {
        let lower = knots[0].0;
        let upper = knots[knots.len() - 1].0;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: axis.into(),
                    abscissa_dims: Dims::NONE,
                    knots,
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        280.0,
                        380.0,
                    )
                    .with_quantity("p", pressure_quantity(), 100_000.0, 100_000.0)
                    .with_quantity("Xm", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0)
                    .with_quantity("Xv", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0)
                    .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), lower, upper),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured crystallinity-dependent mechanics".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "load-bearing amorphous parent and crystal".into(),
            process: "kinetic source integration fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn branch(card: &MaterialCard, datum: f64, rho: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        datum,
        &[280.0, 380.0].map(|temperature_k| HeatCapacityKnot {
            temperature_k,
            specific_heat_capacity_j_kg_k: 1000.0,
            bulk_density_kg_m3: rho,
        }),
        1e-6,
    )
    .unwrap()
}

fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 330.0)
        .unwrap()
        .with_quantity("p", pressure_quantity(), 100_000.0)
        .unwrap()
        .with_quantity("Xm", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap()
        .with_quantity("Xv", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap()
}

fn requirements() -> [ScalarPropertyRequirement; 3] {
    [
        (
            "density",
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        ),
        (
            "young_modulus",
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        ),
        (
            "poisson_ratio",
            Dims::NONE,
            ScalarAdmissibility::OpenInterval {
                lower: -1.0,
                upper: 0.5,
            },
        ),
    ]
    .map(|(name, dims, range)| ScalarPropertyRequirement::try_new(name, dims, range).unwrap())
}

fn point<'a>(card: &'a MaterialCard, context: &'a QueryPoint) -> CrystallizationMaterialPoint<'a> {
    CrystallizationMaterialPoint {
        material: card,
        context,
        temperature_axis: "T",
        crystalline_mass_fraction_axis: Some("Xm"),
        crystalline_volume_fraction_axis: Some("Xv"),
    }
}

fn limits() -> CrystallizationStrainLimits {
    CrystallizationStrainLimits {
        maximum_absolute_total_mandel_strain: 0.02,
        maximum_absolute_elastic_mandel_strain: 0.02,
        maximum_absolute_free_isotropic_strain: 0.01,
    }
}

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 73,
                kernel_id: 73,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_same_temperature_history_drives_source_stiffness_shrinkage_and_native_stress() {
    let card = card();
    let parent = branch(&card, 100_000.0, 1000.0);
    let crystal = branch(&card, 0.0, 1030.0);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 380.0]).unwrap();
    let rate = |_: f64| Ok(1.0);
    let kinetic = pair
        .law(&rate, ContentHash([0x73; 32]), 1.0, 0.8, 400.0)
        .unwrap();
    let initial = kinetic.state_at_temperature(330.0, 0.0).unwrap();
    let context = context();
    let source = CrystallizationPropertySource::new(
        pair,
        kinetic,
        initial,
        point(&card, &context),
        &requirements(),
        MaterialPropertySelection::SingleClaimOnly,
        CrystallizationSolve::default(),
    )
    .unwrap();
    // An amorphous solid can already carry load. This is an explicitly chosen
    // reference, not an inferred crystalline percolation threshold.
    let lock = CrystallizationLock::new(&pair, 330.0, 0.0, [0.0; 6]).unwrap();
    let initial_source = source.evaluate(initial).unwrap();
    let original_claim = initial_source
        .material()
        .property("young_modulus")
        .unwrap()
        .answer()
        .receipt
        .selected;
    let mut tensions = Vec::new();
    for elapsed_s in [0.5, 1.5] {
        let growth = kinetic.advance_isothermal(initial, elapsed_s).unwrap();
        let snapshot = source.evaluate(growth.state).unwrap();
        assert_eq!(snapshot.phase(), growth.state);
        assert_eq!(snapshot.phase().temperature_k(), 330.0);
        let x = growth.state.crystalline_mass_fraction();
        let expected_x = 0.8 * (1.0 - (-elapsed_s).exp());
        close(x, expected_x, 2e-15);
        close(growth.specific_heat_j_kg, -100_000.0 * x, 1e-9);
        let volume = snapshot.volume();
        let expected_volume = (1.0 - x) / 1000.0 + x / 1030.0;
        close(volume.specific_volume_m3_kg, expected_volume, 1e-18);
        assert!(volume.crystalline_volume_fraction < x);
        let material =
            IsotropicElasticStatePoint::try_from_resolved(snapshot.material().clone()).unwrap();
        // rho(phi) is affine in VOLUME fraction, while E(X) uses MASS fraction.
        close(material.density_kg_m3(), expected_volume.recip(), 5e-12);
        let expected_modulus = 1e9 + 2e9 * expected_x;
        close(material.young_modulus_pa(), expected_modulus, 1e-5);
        assert_eq!(
            material
                .resolved()
                .property("young_modulus")
                .unwrap()
                .answer()
                .receipt
                .selected,
            original_claim
        );
        let elastic =
            IsotropicElastic::new(material.young_modulus_pa(), material.poisson_ratio(), 0.02)
                .unwrap();
        let law = lock.freeze(&pair, 330.0, x, &elastic, limits()).unwrap();
        let expected_strain = (1000.0 * expected_volume).cbrt() - 1.0;
        close(law.free_isotropic_strain(), expected_strain, 3e-16);
        let expected_tension = -expected_modulus * expected_strain / (1.0 - 2.0 * 0.3);
        tensions.push(expected_tension);
        with_cx(|cx| {
            for mesh in [
                Mesh2::triangles(2.0, 1.0, 2, 2),
                Mesh2::quads(2.0, 1.0, 2, 2),
            ] {
                let right = mesh.nodes.iter().position(|&p| p == [2.0, 0.0]).unwrap();
                for constrained in [false, true] {
                    let prescribed: Vec<_> = if constrained {
                        (0..2 * mesh.nodes.len()).map(|dof| (dof, 0.0)).collect()
                    } else {
                        vec![(0, 0.0), (1, 0.0), (2 * right + 1, 0.0)]
                    };
                    let problem = SmallStrainPlaneProblem {
                        mesh: &mesh,
                        laws: SmallStrainLawField::Uniform(&law),
                        plane: if constrained {
                            PlaneKind::Strain
                        } else {
                            PlaneKind::Stress
                        },
                        formulation: Formulation::Standard,
                        thickness_m: 0.2,
                        reduction: ReductionConfig {
                            max_iterations: 8,
                            stress_tolerance_pa: 1e-6,
                            strain_limit: 0.02,
                        },
                        rotation_limit_rad: 0.02,
                        prescribed_m: &prescribed,
                        nodal_forces_n: &[],
                        load_path: &[1.0],
                        budget: PlaneAssemblyBudget::default(),
                        settings: IncrementSettings::default(),
                    };
                    let solved = problem.solve(cx).unwrap();
                    for (position, displacement) in mesh.nodes.iter().zip(&solved.displacement_m) {
                        for axis in 0..2 {
                            close(
                                displacement[axis],
                                if constrained {
                                    0.0
                                } else {
                                    expected_strain * position[axis]
                                },
                                2e-10,
                            );
                        }
                    }
                    for stress in &solved.quadrature_stress_pa {
                        for (component, value) in stress.iter().enumerate() {
                            close(
                                *value,
                                if constrained && component < 3 {
                                    expected_tension
                                } else {
                                    0.0
                                },
                                0.05,
                            );
                        }
                    }
                }
            }
        });
    }
    assert!(tensions[1] > tensions[0]);
    assert_eq!(
        source.evaluate(initial).unwrap().material(),
        initial_source.material()
    );
}

#[test]
fn g0_mechanical_source_history_refuses_switching_after_fraction_support_ends() {
    let card = card();
    let parent = branch(&card, 100_000.0, 1000.0);
    let crystal = branch(&card, 0.0, 1030.0);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 380.0]).unwrap();
    let rate = |_: f64| Ok(1.0);
    let kinetic = pair
        .law(&rate, ContentHash([0x73; 32]), 1.0, 0.8, 400.0)
        .unwrap();
    let initial = kinetic.state_at_temperature(330.0, 0.0).unwrap();
    let context = context();
    let source = CrystallizationPropertySource::new(
        pair,
        kinetic,
        initial,
        point(&card, &context),
        &requirements(),
        MaterialPropertySelection::SingleClaimOnly,
        CrystallizationSolve::default(),
    )
    .unwrap();
    let accepted = source.evaluate(initial).unwrap();
    let later = kinetic.advance_isothermal(initial, 4.0).unwrap().state;
    let x = later.crystalline_mass_fraction();
    assert!(x > 0.76 && x < 0.8);
    let volume = pair.volume_state(later.temperature_k(), x).unwrap();
    let query = context
        .clone()
        .with("Xm", x)
        .unwrap()
        .with("Xv", volume.crystalline_volume_fraction)
        .unwrap();
    let replacement = resolve_isotropic_elastic_state_point(
        &card,
        &query,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_ne!(
        replacement
            .resolved()
            .property("young_modulus")
            .unwrap()
            .answer()
            .receipt
            .selected,
        accepted
            .material()
            .property("young_modulus")
            .unwrap()
            .answer()
            .receipt
            .selected,
    );
    assert!(source.evaluate(later).is_err());
    assert_eq!(
        source.evaluate(initial).unwrap().material(),
        accepted.material()
    );
    assert_eq!(source.evaluate(initial).unwrap().phase(), initial);
}
