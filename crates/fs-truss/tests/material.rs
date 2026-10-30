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
use fs_truss::{GroundLimits, GroundStructure, LayoutCase, LayoutLimits, PdhgSettings};

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
    let mut claims = ClaimSet::new();
    for (key, dims, values) in [
        (DENSITY_PROPERTY, Density::DIMS, density),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, young),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, strength),
    ] {
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
