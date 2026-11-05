//! G0/G3: material substitutions propagate through the actual LP and sizing.
//! Synthetic curves are independently chosen numerical fixtures, not measured
//! alloy data or an experimental claim about a building.

use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Pressure, Temperature};
use fs_truss::material::{MaterialLayout, TrussMaterial};
use fs_truss::{
    GroundLimits, GroundStructure, LayoutCase, LayoutCertificateLimits, LayoutCertificateStatus,
    LayoutLimits, PdhgSettings,
};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 1,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn card(density: [f64; 2], young: [f64; 2], strength: [f64; 2]) -> MaterialCard {
    card_on_axis(density, young, strength, "T", Temperature::DIMS)
}

fn card_on_axis(
    density: [f64; 2],
    young: [f64; 2],
    strength: [f64; 2],
    axis: &str,
    axis_dims: fs_qty::Dims,
) -> MaterialCard {
    property_card(
        &[
            (DENSITY_PROPERTY, Density::DIMS, density),
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, young),
            (YIELD_STRESS_PROPERTY, Pressure::DIMS, strength),
        ],
        axis,
        axis_dims,
    )
}

fn property_card(
    properties: &[(&str, fs_qty::Dims, [f64; 2])],
    axis: &str,
    axis_dims: fs_qty::Dims,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for &(key, dims, values) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Curve {
                    abscissa: axis.into(),
                    abscissa_dims: axis_dims,
                    knots: vec![(300.0, values[0]), (600.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with(axis, 300.0, 600.0),
                uncertainty: UncertaintyModel::RelativeHalfWidth {
                    fraction: 0.01,
                    confidence: 0.95,
                },
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic truss property curve".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic truss material".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn static_sizing_refuses_unowned_frequency_and_rate_laws() {
    for axis in ["omega", "strain_rate"] {
        let source = card_on_axis(
            [8000.0; 2],
            [200e9; 2],
            [300e6; 2],
            axis,
            fs_qty::Frequency::DIMS,
        );
        let point = QueryPoint::new().with(axis, 300.0).unwrap();
        assert!(
            TrussMaterial::resolve(&source, &point, MaterialPropertySelection::SingleClaimOnly)
                .is_err()
        );
    }
}

#[test]
fn unrepresentable_physical_volume_is_not_published() {
    with_cx(|cx, _| {
        let source = card([8000.0; 2], [200e9; 2], [1e-320; 2]);
        let states = [material(&source, 300.0)];
        let (ground, case) = two_members(cx, 1.0);
        let layout = MaterialLayout::try_new(
            &ground,
            &case,
            &states,
            &[0, 0],
            LayoutLimits::default(),
            cx,
        )
        .unwrap();
        assert!(matches!(
            layout.solve(PdhgSettings::default()),
            Err(fs_truss::material::MaterialLayoutSolveError::UnrepresentableVolume)
        ));
    });
}

fn material(card: &MaterialCard, temperature: f64) -> TrussMaterial {
    TrussMaterial::resolve(
        card,
        &QueryPoint::new().with("T", temperature).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn two_members(cx: &Cx<'_>, force: f64) -> (GroundStructure, LayoutCase) {
    let ground = GroundStructure::try_from_parts(
        &[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]],
        &[(0, 1), (1, 2)],
        &[1.0, 1.0],
        GroundLimits::default(),
        cx,
    )
    .unwrap();
    let case = LayoutCase::try_new(
        vec![[true, true], [false, true], [false, true]],
        vec![[0.0, 0.0], [0.0, 0.0], [force, 0.0]],
        3,
    )
    .unwrap();
    (ground, case)
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "{actual} != {expected}"
    );
}

const TENSILE: &str = "test_tensile_limit";
const COMPRESSIVE: &str = "test_compressive_limit";

fn axial_card(tensile: [f64; 2], compressive: [f64; 2]) -> MaterialCard {
    property_card(
        &[
            (DENSITY_PROPERTY, Density::DIMS, [2.0, 8.0]),
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, [1e12; 2]),
            (TENSILE, Pressure::DIMS, tensile),
            (COMPRESSIVE, Pressure::DIMS, compressive),
        ],
        "T",
        Temperature::DIMS,
    )
}

fn axial_material(card: &MaterialCard, temperature: f64) -> TrussMaterial {
    TrussMaterial::resolve_axial_limits(
        card,
        &QueryPoint::new().with("T", temperature).unwrap(),
        TENSILE,
        COMPRESSIVE,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn explicit_equal_yield_limits_retain_the_symmetric_material() {
    let source = card([8000.0; 2], [200e9; 2], [300e6; 2]);
    let symmetric = material(&source, 300.0);
    let named = TrussMaterial::resolve_axial_limits(
        &source,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        YIELD_STRESS_PROPERTY,
        YIELD_STRESS_PROPERTY,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(named, symmetric);
    assert_eq!(named.resolved().properties().len(), 3);
    assert_eq!(named.tensile_limit_pa(), named.compressive_limit_pa());
    assert_eq!(named.tension_property(), YIELD_STRESS_PROPERTY);
    assert_eq!(named.compression_property(), YIELD_STRESS_PROPERTY);
}

#[test]
fn asymmetric_limits_reach_objective_refitted_sizing_mass_and_certificate() {
    with_cx(|cx, _| {
        let source = axial_card([100.0, 50.0], [400.0, 200.0]);
        let states = [axial_material(&source, 300.0), axial_material(&source, 600.0)];
        for force in [1.0, -1.0] {
            let (ground, case) = two_members(cx, force);
            let layout = MaterialLayout::try_new(
                &ground, &case, &states, &[0, 1], LayoutLimits::default(), cx,
            )
            .unwrap();
            assert_eq!(layout.lp().c(), &[1.0, 2.0, 0.25, 0.5]);
            assert_eq!(layout.strength_scale_pa(), 100.0);
            let settings = PdhgSettings::default();
            let solution = layout.solve(settings).unwrap();
            let (limits, catalog_areas, volume) = if force > 0.0 {
                ([100.0, 50.0], [0.011, 0.021], 0.03)
            } else {
                ([400.0, 200.0], [0.003, 0.006], 0.0075)
            };
            relative(solution.volume_m3, volume, 1e-4);
            // Deliberately supply the opposite initial force sign: the common
            // equilibrium refit must choose the actual load case's limit.
            let wrong_sign = if force > 0.0 {
                [0.0, 0.0, 1.0, 1.0]
            } else {
                [1.0, 1.0, 0.0, 0.0]
            };
            let sized = layout
                .size_and_snap(&wrong_sign, &[0.003, 0.006, 0.011, 0.021], 1e-3, cx)
                .unwrap();
            assert!(sized.audit.all_pass);
            for (i, member) in sized.audit.members.iter().enumerate() {
                relative(member.force, force, 1e-12);
                relative(member.area_yield, 1.0 / limits[i], 1e-12);
                assert_eq!(member.area_catalog, catalog_areas[i]);
                assert_eq!(member.area_buckling > 0.0, force < 0.0);
            }
            relative(sized.mass_kg.unwrap(), 2.0 * catalog_areas[0] + 8.0 * catalog_areas[1], 1e-12);
            let LayoutCertificateStatus::Certified(certificate) = layout.lp()
                .certify_optimum(&solution.forces, &solution.dual, settings,
                    LayoutCertificateLimits::default(), cx).unwrap()
            else { panic!("the well-conditioned two-member oracle must certify") };
            let oracle = volume * layout.strength_scale_pa();
            let bounds = certificate.bounds();
            assert!(bounds.lower() <= oracle && oracle <= bounds.upper());
            assert!(bounds.lower() > 0.45 * oracle);
            assert!(bounds.upper() < oracle * (1.0 + 1e-10));
        }
        assert_eq!(states[0].tension_property(), TENSILE);
        assert_eq!(states[1].compression_property(), COMPRESSIVE);
        assert_eq!(states[0].resolved().property(TENSILE).unwrap().value_si(), 100.0);
        assert_eq!(states[1].resolved().property(COMPRESSIVE).unwrap().value_si(), 200.0);
    });
}

#[test]
fn reversal_changes_the_optimal_load_path_with_the_same_material_assignments() {
    with_cx(|cx, _| {
        let ground = GroundStructure::try_from_parts(
            &[[0.0, -1.0], [0.0, 0.0], [0.0, 1.0], [1.0, 0.0]],
            &[(0, 3), (1, 3), (2, 3)],
            &[2.0_f64.sqrt(), 1.0, 2.0_f64.sqrt()],
            GroundLimits::default(), cx,
        ).unwrap();
        let source = axial_card([1000.0, 100.0], [100.0, 1000.0]);
        let states = [axial_material(&source, 300.0), axial_material(&source, 600.0)];
        for force in [1.0, -1.0] {
            let case = LayoutCase::try_new(
                vec![[true, true], [true, true], [true, true], [false, false]],
                vec![[0.0; 2], [0.0; 2], [0.0; 2], [force, 0.0]], 4,
            ).unwrap();
            let layout = MaterialLayout::try_new(
                &ground, &case, &states, &[1, 0, 1], LayoutLimits::default(), cx,
            ).unwrap();
            let solution = layout.solve(PdhgSettings::default()).unwrap();
            if force > 0.0 {
                assert!(solution.forces[1] > 0.99);
                relative(solution.volume_m3, 1.0 / 1000.0, 1e-4);
            } else {
                assert!(solution.forces[4] < 1e-4);
                assert!(solution.forces[3] > 0.7 && solution.forces[5] > 0.7);
                relative(solution.volume_m3, 2.0 / 1000.0, 1e-4);
            }
        }
    });
}

#[test]
fn asymmetric_limits_persist_through_self_weight_rebuilds() {
    with_cx(|cx, _| {
        let source = axial_card([100.0; 2], [1000.0; 2]);
        let states = [axial_material(&source, 300.0)];
        for force in [1.0, -1.0] {
            let (ground, case) = two_members(cx, force);
            let layout = MaterialLayout::try_new(
                &ground, &case, &states, &[0, 0], LayoutLimits::default(), cx,
            ).unwrap();
            let acceleration = [10.0 * force, 0.0];
            let design = layout.solve_with_body_loads(
                acceleration, &[0.002, 0.006, 0.014, 0.02, 0.04], 1e-3,
                PdhgSettings { max_iters: 50_000, gap_tol: 1e-9, check_every: 50 },
                8, cx,
            ).unwrap();
            let (areas, limit) = if force > 0.0 {
                (vec![0.02, 0.014], 100.0)
            } else {
                (vec![0.002; 2], 1000.0)
            };
            assert_eq!(design.areas_m2, areas);
            assert_eq!(design.body_loads, layout.body_loads(&areas, acceleration, cx).unwrap());
            for member in &design.design.audit.members {
                assert_eq!(member.force > 0.0, force > 0.0);
                relative(member.area_yield, member.force.abs() / limit, 1e-12);
            }
            assert!(design.design.audit.all_pass);
        }
    });
}

#[test]
fn named_axial_limits_require_both_static_positive_pressure_sources() {
    let resolve = |source: &MaterialCard, point: &QueryPoint| {
        TrussMaterial::resolve_axial_limits(
            source, point, TENSILE, COMPRESSIVE, MaterialPropertySelection::SingleClaimOnly,
        )
    };
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    let incomplete = property_card(&[
        (DENSITY_PROPERTY, Density::DIMS, [2.0; 2]),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, [1e12; 2]),
        (TENSILE, Pressure::DIMS, [100.0; 2]),
    ], "T", Temperature::DIMS);
    assert!(resolve(&incomplete, &point).is_err());
    let mistyped = property_card(&[
        (DENSITY_PROPERTY, Density::DIMS, [2.0; 2]),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, [1e12; 2]),
        (TENSILE, Pressure::DIMS, [100.0; 2]),
        (COMPRESSIVE, Temperature::DIMS, [1000.0; 2]),
    ], "T", Temperature::DIMS);
    assert!(resolve(&mistyped, &point).is_err());
    for (tensile, compressive) in [(0.0, 1000.0), (100.0, 0.0), (100.0, -1.0)] {
        assert!(resolve(&axial_card([tensile; 2], [compressive; 2]), &point).is_err());
    }
    let dynamic = property_card(&[
        (DENSITY_PROPERTY, Density::DIMS, [2.0; 2]),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, [1e12; 2]),
        (TENSILE, Pressure::DIMS, [100.0; 2]),
        (COMPRESSIVE, Pressure::DIMS, [1000.0; 2]),
    ], "strain_rate", fs_qty::Frequency::DIMS);
    assert!(resolve(&dynamic, &QueryPoint::new().with("strain_rate", 300.0).unwrap()).is_err());
}

#[test]
fn temperature_and_member_materials_change_yield_cost_sizing_and_mass() {
    with_cx(|cx, _| {
        let source = card([8000.0, 7900.0], [200e9, 50e9], [300e6, 100e6]);
        let states = [material(&source, 300.0), material(&source, 600.0)];
        assert_ne!(
            states[0].resolved().identity(),
            states[1].resolved().identity()
        );
        let (ground, case) = two_members(cx, 1.0);
        let layout = MaterialLayout::try_new(
            &ground,
            &case,
            &states,
            &[0, 1],
            LayoutLimits::default(),
            cx,
        )
        .unwrap();
        assert_eq!(layout.lp().c(), &[1.0, 3.0, 1.0, 3.0]);
        let result = layout.solve(PdhgSettings::default()).unwrap();
        relative(result.volume_m3, 1.0 / 300e6 + 1.0 / 100e6, 1e-4);
        let sized = layout
            .size_and_snap(&result.forces, &[4e-9, 1.1e-8, 2e-8], 1e-3, cx)
            .unwrap();
        assert!(sized.audit.all_pass);
        assert!(sized.audit.rows.iter().all(|row| !row.contains("aisc")));
        assert_eq!(sized.audit.members[0].area_catalog, 4e-9);
        assert_eq!(sized.audit.members[1].area_catalog, 1.1e-8);
        relative(
            sized.mass_kg.unwrap(),
            4e-9 * 8000.0 + 1.1e-8 * 7900.0,
            1e-12,
        );
        relative(sized.volume_m3.unwrap(), 4e-9 + 1.1e-8, 1e-12);
        assert_eq!(
            layout.material(1).unwrap().resolved().identity(),
            states[1].resolved().identity()
        );
    });
}

#[test]
fn compression_buckling_uses_each_members_modulus_not_global_steel() {
    with_cx(|cx, _| {
        let source = card([8000.0, 8000.0], [200e9, 50e9], [300e6, 300e6]);
        let states = [material(&source, 300.0), material(&source, 600.0)];
        let (ground, case) = two_members(cx, -1.0);
        let layout = MaterialLayout::try_new(
            &ground,
            &case,
            &states,
            &[0, 1],
            LayoutLimits::default(),
            cx,
        )
        .unwrap();
        let sized = layout
            .size_and_snap(&[0.0, 0.0, 1.0, 1.0], &[3e-6, 6e-6, 1e-5], 1e-3, cx)
            .unwrap();
        assert!(sized.audit.all_pass);
        relative(
            sized.audit.members[1].area_buckling,
            2.0 * sized.audit.members[0].area_buckling,
            1e-12,
        );
        assert_eq!(sized.audit.members[0].area_catalog, 3e-6);
        assert_eq!(sized.audit.members[1].area_catalog, 6e-6);
        relative(sized.mass_kg.unwrap(), 0.072, 1e-12);
    });
}

#[test]
fn swapping_member_assignments_changes_actual_optimum_load_path() {
    with_cx(|cx, _| {
        // Direct tie and the two inclined ties are alternative load paths.
        // Lower yield costs must select the newly assigned stronger path.
        let ground = GroundStructure::try_from_parts(
            &[[0.0, -1.0], [0.0, 0.0], [0.0, 1.0], [1.0, 0.0]],
            &[(0, 3), (1, 3), (2, 3)],
            &[2.0_f64.sqrt(), 1.0, 2.0_f64.sqrt()],
            GroundLimits::default(),
            cx,
        )
        .unwrap();
        let case = LayoutCase::try_new(
            vec![[true, true], [true, true], [true, true], [false, false]],
            vec![[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [1.0, 0.0]],
            4,
        )
        .unwrap();
        let source = card([8000.0; 2], [200e9; 2], [100e6, 1e9]);
        let states = [material(&source, 300.0), material(&source, 600.0)];
        let direct = MaterialLayout::try_new(
            &ground,
            &case,
            &states,
            &[0, 1, 0],
            LayoutLimits::default(),
            cx,
        )
        .unwrap()
        .solve(PdhgSettings::default())
        .unwrap();
        let inclined = MaterialLayout::try_new(
            &ground,
            &case,
            &states,
            &[1, 0, 1],
            LayoutLimits::default(),
            cx,
        )
        .unwrap()
        .solve(PdhgSettings::default())
        .unwrap();
        assert!(direct.forces[1] > 0.99);
        assert!(inclined.forces[1] < 1e-4);
        relative(direct.volume_m3, 1e-9, 1e-4);
        relative(inclined.volume_m3, 2e-9, 1e-4);
    });
}

#[test]
fn incomplete_sources_invalid_assignments_and_cancelled_publication_refuse() {
    let source = card([8000.0; 2], [200e9; 2], [300e6; 2]);
    assert!(
        TrussMaterial::resolve(
            &source,
            &QueryPoint::new().with("T", 700.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    with_cx(|cx, gate| {
        let states = [material(&source, 300.0)];
        let (ground, case) = two_members(cx, -1.0);
        assert!(
            MaterialLayout::try_new(&ground, &case, &states, &[0], LayoutLimits::default(), cx)
                .is_err()
        );
        assert!(
            MaterialLayout::try_new(
                &ground,
                &case,
                &states,
                &[0, 1],
                LayoutLimits::default(),
                cx
            )
            .is_err()
        );
        let layout = MaterialLayout::try_new(
            &ground,
            &case,
            &states,
            &[0, 0],
            LayoutLimits::default(),
            cx,
        )
        .unwrap();
        assert!(layout.size_and_snap(&[], &[1e-5], 0.0, cx).is_err());
        assert!(
            layout
                .size_and_snap(&[0.0, 0.0, 1.0, 1.0], &[2e-5, 1e-5], 0.0, cx)
                .is_err()
        );
        let missing = layout
            .size_and_snap(&[0.0, 0.0, 1.0, 1.0], &[1e-12], 0.0, cx)
            .unwrap();
        assert!(!missing.audit.all_pass);
        assert_eq!(missing.mass_kg, None);
        assert_eq!(missing.member_mass_kg, vec![None, None]);
        gate.request();
        assert!(
            layout
                .size_and_snap(&[0.0, 0.0, 1.0, 1.0], &[1e-5], 0.0, cx)
                .is_err()
        );
    });
}
