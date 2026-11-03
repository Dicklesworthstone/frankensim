//! Reversible adhesion data reach the existing sparse equilibrium wall operator.
use super::*;
use fs_material::liquid::surface::adhesion::{
    EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY, EquilibriumAngleOrigin, resolve_young_dupre_wetting,
};

fn adhesion_card(
    liquid: &ResolvedLiquidSurfaceState,
    work: Option<PropertyClaim>,
) -> InterfaceSystemCard {
    let context = interface_card(liquid, None, Dims::NONE, false, None);
    let mut claims = ClaimSet::new();
    if let Some(work) = work {
        claims.insert_claim(work).unwrap();
    }
    InterfaceSystemCard::assemble(
        context.surface_a().clone(),
        context.surface_b().clone(),
        context.context().clone(),
        claims,
        vec![],
    )
    .unwrap()
}
fn work(value: f64) -> PropertyClaim {
    claim(
        EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY,
        value,
        SURFACE_TENSION_DIMS,
    )
}

#[test]
fn original_adhesion_and_tension_drive_native_capillary_response() {
    let liquid = liquid(0, 300.0, Some(0.05));
    let cards = [
        adhesion_card(&liquid, Some(work(0.025))),
        adhesion_card(&liquid, Some(work(0.075))),
    ];
    let boundaries: Vec<_> = cards
        .iter()
        .map(|card| {
            resolve_young_dupre_wetting(card, &liquid, MaterialPropertySelection::SingleClaimOnly)
                .unwrap()
        })
        .collect();
    for boundary in &boundaries {
        assert_eq!(boundary.angle_origin(), EquilibriumAngleOrigin::YoungDupre);
        assert!(
            boundary
                .resolved()
                .property(EQUILIBRIUM_CONTACT_ANGLE_PROPERTY)
                .is_none()
        );
        assert_eq!(
            boundary.surface_tension_source(),
            liquid.material().property(SURFACE_TENSION_PROPERTY)
        );
    }
    assert!((boundaries[0].contact_angle_radians() - 2.0 * PI / 3.0).abs() < 2e-15);
    assert!((boundaries[1].contact_angle_radians() - PI / 3.0).abs() < 2e-15);
    let mut flows: Vec<_> = boundaries
        .iter()
        .map(|boundary| {
            plan(&liquid)
                .with_planar_wetting(PlanarWall3::ZMin, boundary)
                .unwrap()
                .free_surface3(domain(), drop_shape)
                .unwrap()
        })
        .collect();
    let initial = flows[0].ledger_mass_kg().unwrap();
    for _ in 0..12 {
        for flow in &mut flows {
            flow.step().unwrap();
        }
    }
    assert_ne!(
        flows[0].lattice().grid().state_bits(),
        flows[1].lattice().grid().state_bits()
    );
    for (flow, boundary) in flows.iter().zip(&boundaries) {
        assert_eq!(
            flow.lattice().planar_contact_angle(),
            Some((PlanarWall3::ZMin, boundary.contact_angle_radians()))
        );
        assert!((flow.ledger_mass_kg().unwrap() - initial).abs() / initial < 1e-10);
    }
}

#[test]
fn liquid_tension_substitution_and_missing_partial_wetting_data_refuse() {
    let original = liquid(0, 300.0, Some(0.05));
    let card = adhesion_card(&original, Some(work(0.06)));
    let boundary =
        resolve_young_dupre_wetting(&card, &original, MaterialPropertySelection::SingleClaimOnly)
            .unwrap();
    for other in [
        liquid(0, 300.0, Some(0.06)),
        liquid(1, 300.0, Some(0.05)),
        liquid(0, 301.0, Some(0.05)),
        liquid(0, 300.0, None),
    ] {
        assert!(
            plan(&other)
                .with_planar_wetting(PlanarWall3::ZMin, &boundary)
                .is_err()
        );
    }
    for bad in [
        None,
        Some(work(0.0)),
        Some(work(0.1)),
        Some(work(0.11)),
        Some(claim(
            EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY,
            0.06,
            Density::DIMS,
        )),
    ] {
        let card = adhesion_card(&original, bad);
        assert!(
            resolve_young_dupre_wetting(
                &card,
                &original,
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    let neglected = liquid(0, 300.0, None);
    assert!(
        resolve_young_dupre_wetting(
            &adhesion_card(&neglected, Some(work(0.06))),
            &neglected,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}

#[test]
fn temperature_refresh_requires_both_original_sources_at_the_same_point() {
    let original = liquid(0, 300.0, Some(0.05));
    let mut source = work(0.05);
    source.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: fs_qty::Temperature::DIMS,
        knots: vec![(290.0, 0.02), (310.0, 0.08)],
        dims: SURFACE_TENSION_DIMS,
    };
    source.interpolation = InterpolationPolicy::LinearInside;
    let card = adhesion_card(&original, Some(source));
    let old =
        resolve_young_dupre_wetting(&card, &original, MaterialPropertySelection::SingleClaimOnly)
            .unwrap();
    let warm = liquid(0, 310.0, Some(0.05));
    let refreshed = old
        .at_query_point_with_liquid(&card, &point(310.0), &warm)
        .unwrap();
    assert!(refreshed.contact_angle_radians() < old.contact_angle_radians());
    assert!((old.contact_angle_radians() - FRAC_PI_2).abs() < 2e-15);
    assert!((refreshed.young_dupre_law().unwrap().work_of_adhesion_j_m2() - 0.08).abs() < 1e-16);
    assert!(old.at_query_point(&card, &point(310.0)).is_err());
    assert!(
        old.at_query_point_with_liquid(&card, &point(310.0), &original)
            .is_err()
    );
    let substituted = liquid(0, 310.0, Some(0.06));
    assert!(
        old.at_query_point_with_liquid(&card, &point(310.0), &substituted)
            .is_err()
    );
    assert_eq!(old.angle_origin(), EquilibriumAngleOrigin::YoungDupre);
}

#[test]
fn directly_sourced_angles_keep_original_refresh_behavior() {
    let liquid = liquid(0, 300.0, Some(0.05));
    let card = interface_card(&liquid, Some(FRAC_PI_2), Dims::NONE, false, None);
    let boundary = resolve_liquid_wetting_state(
        &card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(
        boundary.angle_origin(),
        EquilibriumAngleOrigin::SourcedAngle
    );
    assert!(boundary.young_dupre_law().is_none());
    assert!(boundary.surface_tension_source().is_none());
    let refreshed = boundary.at_query_point(&card, &point(310.0)).unwrap();
    assert_eq!(refreshed.contact_angle_radians(), FRAC_PI_2);
}

fn domain2() -> fs_lbm::material::free_surface::MaterialFreeSurfaceDomain2 {
    fs_lbm::material::free_surface::MaterialFreeSurfaceDomain2 {
        cells: [18, 14],
        gravity_m_s2: [0.0; 2],
        contact: fs_lbm::freesurface::ContactModel::Neutral,
        max_cells: 18 * 14,
    }
}

fn drop2(x: usize, y: usize) -> bool {
    let x = x as f64 - 8.5;
    let y = y as f64 - 0.5;
    x * x + y * y < 25.0
}

#[test]
fn g1_g3_adhesion_drives_dense_capillary_flow_with_mass_per_depth() {
    use fs_lbm::freesurface::PlanarWall2;
    let liquid = liquid(0, 300.0, Some(0.05));
    let cards = [adhesion_card(&liquid, Some(work(0.025))),
        adhesion_card(&liquid, Some(work(0.075)))];
    let boundaries: Vec<_> = cards.iter().map(|card| resolve_young_dupre_wetting(
        card, &liquid, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()).collect();
    let make = |boundary| plan(&liquid).with_planar_wetting(PlanarWall3::YMin, boundary)
        .unwrap().free_surface2(domain2(), drop2).unwrap();
    let mut flows: Vec<_> = boundaries.iter().map(make).collect();
    let mut replay = make(&boundaries[0]);
    let initial = flows[0].ledger_mass_kg_per_m().unwrap();
    let expected = flows[0].lattice().ledger_mass() * liquid.density_kg_m3() * 1e-6;
    assert!((initial - expected).abs() < 1e-14);
    for _ in 0..12 {
        for flow in &mut flows { flow.step(); }
        replay.step();
    }
    assert_ne!(flows[0].lattice().grid.f, flows[1].lattice().grid.f);
    assert_eq!(flows[0].lattice().grid.f, replay.lattice().grid.f);
    assert_eq!(flows[0].lattice().mass, replay.lattice().mass);
    for (flow, boundary) in flows.iter().zip(&boundaries) {
        assert_eq!(flow.lattice().planar_contact_angle(),
            Some((PlanarWall2::YMin, boundary.contact_angle_radians())));
        assert!((flow.ledger_mass_kg_per_m().unwrap() - initial).abs() / initial < 1e-10);
        let sample = flow.sample_at_si(8, 1).unwrap().unwrap();
        assert!(sample.density_kg_m3 > 0.0);
        assert!(sample.velocity_m_s.iter().all(|v| v.is_finite()));
        assert_eq!(flow.sample_at_si(0, 0).unwrap(), None);
        assert_eq!(flow.sample_at_si(usize::MAX, 1).unwrap(), None);
    }
}

#[test]
fn g3_dense_equilibrium_ninety_degrees_is_exactly_neutral() {
    let liquid = liquid(0, 300.0, Some(0.05));
    let boundary = wetting(&liquid, FRAC_PI_2);
    let mut neutral = plan(&liquid).free_surface2(domain2(), drop2).unwrap();
    let mut sourced = plan(&liquid).with_planar_wetting(PlanarWall3::YMin, &boundary)
        .unwrap().free_surface2(domain2(), drop2).unwrap();
    for _ in 0..12 { neutral.step(); sourced.step(); }
    assert_eq!(neutral.lattice().grid.f, sourced.lattice().grid.f);
    assert_eq!(neutral.lattice().mass, sourced.lattice().mass);
    assert_eq!(neutral.ledger_mass_kg_per_m().unwrap(), sourced.ledger_mass_kg_per_m().unwrap());
}

#[test]
fn dense_source_and_domain_refusals_precede_geometry_evaluation() {
    use fs_lbm::freesurface::ContactModel;
    let liquid = liquid(0, 300.0, Some(0.05));
    let boundary = wetting(&liquid, PI / 3.0);
    let fail = |_: usize, _: usize| -> bool { panic!("refused input evaluated geometry") };
    for wall in [PlanarWall3::ZMin, PlanarWall3::ZMax] {
        assert!(plan(&liquid).with_planar_wetting(wall, &boundary).unwrap()
            .free_surface2(domain2(), fail).is_err());
    }
    let scaling = plan(&liquid).with_planar_wetting(PlanarWall3::YMin, &boundary).unwrap();
    for domain in [
        fs_lbm::material::free_surface::MaterialFreeSurfaceDomain2 { contact: ContactModel::Wetting, ..domain2() },
        fs_lbm::material::free_surface::MaterialFreeSurfaceDomain2 { max_cells: 10, ..domain2() },
        fs_lbm::material::free_surface::MaterialFreeSurfaceDomain2 { cells: [3, 10], ..domain2() },
        fs_lbm::material::free_surface::MaterialFreeSurfaceDomain2 { cells: [usize::MAX, 10], ..domain2() },
        fs_lbm::material::free_surface::MaterialFreeSurfaceDomain2 { gravity_m_s2: [f64::NAN, 0.0], ..domain2() },
    ] {
        assert!(scaling.free_surface2(domain, fail).is_err());
    }
}
