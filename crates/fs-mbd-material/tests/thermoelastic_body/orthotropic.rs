//! G1/G3: one data-selected stiffness law reaches real heat, FEM and vibration.
use super::*;
use fs_matdb::SelectionPolicy;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_mbd_material::thermal::heating::{
    ThermalRegion,
    constrained::{ConstrainedHeatedBody, ConstrainedRegion},
};
use fs_thermal::enthalpy::*;

const IDENTITY: [[f64; 3]; 3] = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TURN: [[f64; 3]; 3] = [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]];
fn model(rotation: [[f64; 3]; 3]) -> ThermalStiffnessModel {
    ThermalStiffnessModel::OrientedOrthotropic {
        principal_to_mesh: rotation,
        orientation_identity: ContentHash([47; 32]),
    }
}
fn source(base: MaterialCard, e: [f64; 3], nu: [f64; 3], g: [f64; 3], upper: f64) -> MaterialCard {
    let mut claims = base.claims().clone();
    for (names, dims, values, scale) in [
        (ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES, Pressure::DIMS, e, 0.5),
        (ORTHOTROPIC_POISSON_RATIO_PROPERTIES, Dims::NONE, nu, 1.),
        (ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES, Pressure::DIMS, g, 0.5),
    ] {
        for (name, value) in names.into_iter().zip(values) {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(300., value), (upper, value * scale)],
                        dims,
                    },
                    validity: ValidityDomain::unconstrained()
                        .with("T", 300., upper)
                        .with("pressure", 1e5, 1e5),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::LinearInside,
                    observations: vec![],
                    provenance: Provenance {
                        source: format!("synthetic orthotropic {name}; not measurements"),
                        license: "CC0-1.0".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

#[test]
fn isotropic_limit_matches_shape_reactions_mass_and_prestressed_pencil() {
    let cards = [source(
        card("a", 6., [1000., 500.], 1e-4),
        [1000.; 3],
        [0.25; 3],
        [400.; 3],
        400.,
    )];
    let original = reference(&cards);
    let mesh =
        TetrahedralMass::bind(&original, &NODES[..4], &TETS[..1], &KEYS[..1], 4, || false).unwrap();
    with_cx(|cx| {
        let mut cells = inputs(&cards, &[400.]);
        let forces = [(0, 2.), (3, 1.)];
        let a = mesh
            .solve_thermal_with_nodal_forces(&cells, &SUPPORT, &forces, config(), cx)
            .unwrap();
        cells[0].stiffness_model = model(IDENTITY);
        let b = mesh
            .solve_thermal_with_nodal_forces(&cells, &SUPPORT, &forces, config(), cx)
            .unwrap();
        for (x, y) in a
            .body()
            .nodes_m()
            .iter()
            .flatten()
            .zip(b.body().nodes_m().iter().flatten())
        {
            near(*x, *y);
        }
        for (x, y) in a
            .support_reactions_n()
            .iter()
            .flatten()
            .zip(b.support_reactions_n().iter().flatten())
        {
            near(*x, *y);
        }
        near(a.body().inertia().mass(), b.body().inertia().mass());
        near(a.stress().elastic_energy_j(), b.stress().elastic_energy_j());
        let (pa, pb) = (
            a.initial_stress_pencil(cx).unwrap(),
            b.initial_stress_pencil(cx).unwrap(),
        );
        assert_eq!(pa.pencil().mass.to_dense(), pb.pencil().mass.to_dense());
        for (x, y) in pa
            .pencil()
            .stiffness
            .to_dense()
            .iter()
            .zip(pb.pencil().stiffness.to_dense())
        {
            near(*x, y);
        }
        assert_eq!(b.sources()[0].stiffness().properties().len(), 9);
        assert!(
            b.sources()[0]
                .stiffness()
                .property(DENSITY_PROPERTY)
                .is_none()
        );
        assert_eq!(b.sources()[0].stiffness_model(), model(IDENTITY));
        assert_eq!(
            b,
            mesh.solve_thermal_with_nodal_forces(&cells, &SUPPORT, &forces, config(), cx)
                .unwrap()
        );
    });
}

#[test]
fn axes_and_temperature_change_the_same_loaded_shape_inertia_and_vibration() {
    let cards = [source(
        card("a", 6., [1000., 500.], 1e-4),
        [1000., 2000., 3000.],
        [0.; 3],
        [500., 600., 700.],
        400.,
    )];
    let original = reference(&cards);
    let mesh =
        TetrahedralMass::bind(&original, &NODES[..4], &TETS[..1], &KEYS[..1], 4, || false).unwrap();
    let mut snapshots = Vec::new();
    with_cx(|cx| {
        for (rotation, e) in [(IDENTITY, 1000.), (TURN, 2000.)] {
            for t in [300., 400.] {
                let mut cells = inputs(&cards, &[t]);
                cells[0].stiffness_model = model(rotation);
                let body = mesh
                    .solve_thermal_with_nodal_forces(&cells, &SUPPORT, &[(3, 1.)], config(), cx)
                    .unwrap();
                let current_e = e * if t == 300. { 1. } else { 0.5 };
                near(
                    body.displacement().displacement_m()[1][0],
                    (t - 300.) * 1e-4 + 6. / current_e,
                );
                near(body.stress().elements()[0].stress_mandel_pa()[0], 6.);
                assert_eq!(body.body().inertia().mass(), original.inertia().mass());
                let pencil = body.initial_stress_pencil(cx).unwrap();
                assert_eq!(pencil.stress(), body.stress());
                snapshots.push((body.body().inertia(), pencil));
            }
        }
    });
    assert_ne!(snapshots[0].0, snapshots[1].0);
    assert_ne!(snapshots[1].0, snapshots[3].0);
    assert_ne!(
        snapshots[1].1.pencil().stiffness.to_dense(),
        snapshots[3].1.pencil().stiffness.to_dense()
    );
    for (_, pencil) in &snapshots {
        assert_eq!(
            pencil.pencil().mass.to_dense(),
            snapshots[0].1.pencil().mass.to_dense()
        );
    }
}

fn heat_source(upper: f64) -> MaterialCard {
    let base = card("a", 6., [1000., 500.], 1e-4);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == DENSITY_PROPERTY {
            claim.validity = ValidityDomain::unconstrained()
                .with("T", 300., 400.)
                .with("pressure", 1e5, 1e5);
        }
        claims.insert_claim(claim).unwrap();
    }
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS),
            value: PropertyValue::Scalar {
                value: 10.,
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 300., 400.)
                .with("pressure", 1e5, 1e5),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic caloric fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    source(
        MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap(),
        [1000., 2000., 3000.],
        [0.; 3],
        [500., 600., 700.],
        upper,
    )
}
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    let p = |t| {
        QueryPoint::new()
            .with("T", t)
            .unwrap()
            .with("pressure", 1e5)
            .unwrap()
    };
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &p(300.),
            upper: &p(400.),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.,
        1e-8,
    )
    .unwrap()
}
fn heated<'a>(
    original: &'a ResolvedBody,
    card: &'a MaterialCard,
    curve: &'a ResolvedHeatCapacityCurve,
    cx: &Cx<'_>,
) -> ConstrainedHeatedBody<'a> {
    let mesh =
        TetrahedralMass::bind(original, &NODES[..4], &TETS[..1], &KEYS[..1], 4, || false).unwrap();
    let boundary = |source_w| Boundary {
        ambient_k: 300.,
        htc_w_m2_k: 0.,
        reservoir: None,
        source_w,
    };
    let storage = ThermalRegion {
        key: "a",
        card,
        caloric: curve,
        cell_id: CellId(0),
        material_id: MaterialId(0),
        segment_id: 0,
        body_id: 0,
        geometry: CellGeometry {
            start_m: [0., 0., 0.],
            end_m: [0.5, 0.5, 0.5],
            volume_m3: original.regions()[0].geometry().volume_m3(),
            surface_area_m2: 2.4,
        },
        initial_specific_enthalpy_j_kg: 0.,
        boundary: boundary(100.),
        expansion_selection: MaterialPropertySelection::SingleClaimOnly,
        fixed_point_m: Vec3::default(),
    };
    ConstrainedHeatedBody::new_with_nodal_forces(
        mesh,
        vec![ConstrainedRegion {
            storage,
            stiffness_selection: MaterialPropertySelection::SingleClaimOnly,
            stiffness_model: model(TURN),
        }],
        SUPPORT.to_vec(),
        vec![(3, 1.)],
        config(),
        vec![TimedEvent {
            time_s: 1.,
            event: Event::Boundary(CellId(0), boundary(-100.)),
        }],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 2,
            max_substeps: 100,
        },
        step(),
        cx,
    )
    .unwrap()
}
fn step() -> StepConfig {
    StepConfig {
        max_dt_s: 0.25,
        energy_tolerance_j: 1e-8,
    }
}

#[test]
fn accepted_thermal_substeps_and_restarts_keep_the_selected_material_axes() {
    let cards = [heat_source(400.)];
    let original = reference(&cards);
    let curve = caloric(&cards[0]);
    with_cx(|cx| {
        let body = heated(&original, &cards[0], &curve, cx);
        let (hot, report) = body.advance_to(1., step(), 4, cx).unwrap();
        assert_eq!(report.shape_solves, 4);
        near(report.thermal.external_heat_j, 100.);
        let t = hot.thermal().temperature_k(CellId(0)).unwrap();
        near(t, 310.);
        let mut cells = inputs(&cards, &[t]);
        cells[0].stiffness_model = model(TURN);
        let direct = body
            .mesh()
            .solve_thermal_with_nodal_forces(&cells, &SUPPORT, &[(3, 1.)], config(), cx)
            .unwrap();
        assert_eq!(&direct, hot.configuration());
        let (half, _) = body.advance_to(0.5, step(), 2, cx).unwrap();
        let (resumed, _) = half.advance_to(1., step(), 2, cx).unwrap();
        assert_eq!(hot.configuration(), resumed.configuration());
        let a = hot.configuration().initial_stress_pencil(cx).unwrap();
        let b = resumed.configuration().initial_stress_pencil(cx).unwrap();
        assert_eq!(a.identity(), b.identity());
        assert_eq!(
            a.pencil().stiffness.to_dense(),
            b.pencil().stiffness.to_dense()
        );
        let (cool, _) = hot.advance_to(2., step(), 4, cx).unwrap();
        assert_eq!(cool.configuration(), body.configuration());
        near(body.thermal().time_s(), 0.);
    });
}

#[test]
fn unsupported_orthotropy_refuses_without_isotropic_fallback_or_partial_heat() {
    let cards = [card("a", 6., [1000., 500.], 1e-4)];
    let original = reference(&cards);
    let mesh =
        TetrahedralMass::bind(&original, &NODES[..4], &TETS[..1], &KEYS[..1], 4, || false).unwrap();
    let mut cells = inputs(&cards, &[300.]);
    cells[0].stiffness_model = model(IDENTITY);
    assert!(with_cx(|cx| mesh.solve_thermal(&cells, &SUPPORT, config(), cx)).is_err());
    let cards = [heat_source(305.)];
    let original = reference(&cards);
    let curve = caloric(&cards[0]);
    with_cx(|cx| {
        let body = heated(&original, &cards[0], &curve, cx);
        let before = body.configuration().clone();
        assert!(body.advance_to(2., step(), 8, cx).is_err());
        near(body.thermal().time_s(), 0.);
        assert_eq!(body.configuration(), &before);
        let mut cells = inputs(&cards, &[300.]);
        cells[0].stiffness_model = model([[1., 0., 0.], [0., 2., 0.], [0., 0., 1.]]);
        assert!(
            body.mesh()
                .solve_thermal(&cells, &SUPPORT, config(), cx)
                .is_err()
        );
        cells[0].stiffness_model = ThermalStiffnessModel::OrientedOrthotropic {
            principal_to_mesh: IDENTITY,
            orientation_identity: ContentHash([0; 32]),
        };
        assert!(
            body.mesh()
                .solve_thermal(&cells, &SUPPORT, config(), cx)
                .is_err()
        );
    });
}
