//! Original source/latent charts on physical layered control volumes.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve},
    state_point::THERMAL_CONDUCTIVITY_DIMS,
    thermal::{
        conduction::phase::PreparedPhaseConductivity,
        source::ThermalMaterialPoint,
        transport::stack::{
            LayerStackGeometry,
            phase_transient::{
                BoundPhaseLayeredThermalBody, PhaseLayerSource, PreparedPhaseLayeredThermalBody,
            },
            transient::LayeredThermalIds,
        },
    },
};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_thermal::enthalpy::{
    DepositionBudget, DepositionError, DepositionSimulation, Event, StepConfig,
    conduction::implicit::{
        ImplicitConductionControls, ImplicitConductionWorkspace,
        adaptive::AdaptiveImplicitConductionControls,
    },
};
const PI: f64 = core::f64::consts::PI;

fn material(solid_k: f64, liquid_k: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("k", THERMAL_CONDUCTIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "f".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, solid_k), (1.0, liquid_k)],
                dims: THERMAL_CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    190.0,
                    900.0,
                )
                .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic phase transport oracle, not empirical data".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured conductivity".into(),
            phase: "explicit h/T/f chart".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn caloric(card: &MaterialCard, tm: f64, latent: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-100.0, tm - 100.0, 0.0),
            (0.0, tm, 0.0),
            (latent, tm, 1.0),
            (latent + 400.0, tm + 400.0, 1.0),
        ]
        .into_iter()
        .map(|(h, t, f)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: 1.0,
        })
        .collect(),
    )
    .unwrap()
}
fn prepare<'a>(
    card: &MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
) -> PreparedPhaseConductivity<'a> {
    let context = QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            curve.knots()[0].temperature_k,
        )
        .unwrap()
        .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap();
    PreparedPhaseConductivity::new(
        curve,
        ThermalMaterialPoint {
            material: card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: Some("f"),
        },
        &PropertyKey::new("k", THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        [
            curve.knots()[0].specific_enthalpy_j_kg,
            curve.knots().last().unwrap().specific_enthalpy_j_kg,
        ],
        || Ok(()),
    )
    .unwrap()
}
fn shapes() -> [LayerStackGeometry; 3] {
    [
        LayerStackGeometry::Planar { area_m2: 2.0 },
        LayerStackGeometry::Cylindrical {
            inner_radius_m: 1.0,
            length_m: 0.7,
        },
        LayerStackGeometry::Spherical {
            inner_radius_m: 1.0,
        },
    ]
}
fn start(shape: LayerStackGeometry) -> f64 {
    match shape {
        LayerStackGeometry::Planar { .. } => 0.0,
        LayerStackGeometry::Cylindrical { inner_radius_m, .. }
        | LayerStackGeometry::Spherical { inner_radius_m } => inner_radius_m,
    }
}
fn area(shape: LayerStackGeometry, r: f64) -> f64 {
    match shape {
        LayerStackGeometry::Planar { area_m2 } => area_m2,
        LayerStackGeometry::Cylindrical { length_m, .. } => 2.0 * PI * length_m * r,
        LayerStackGeometry::Spherical { .. } => 4.0 * PI * r * r,
    }
}
fn volume(shape: LayerStackGeometry, a: f64, b: f64) -> f64 {
    match shape {
        LayerStackGeometry::Planar { area_m2 } => area_m2 * (b - a),
        LayerStackGeometry::Cylindrical { length_m, .. } => PI * length_m * (b * b - a * a),
        LayerStackGeometry::Spherical { .. } => 4.0 * PI * (b.powi(3) - a.powi(3)) / 3.0,
    }
}
// Independent dimensional resistance integrals, with no production metrics.
fn metric(shape: LayerStackGeometry, a: f64, b: f64) -> f64 {
    match shape {
        LayerStackGeometry::Planar { area_m2 } => (b - a) / area_m2,
        LayerStackGeometry::Cylindrical { length_m, .. } => (b / a).ln() / (2.0 * PI * length_m),
        LayerStackGeometry::Spherical { .. } => (1.0 / a - 1.0 / b) / (4.0 * PI),
    }
}
fn centroid(shape: LayerStackGeometry, a: f64, b: f64) -> f64 {
    match shape {
        LayerStackGeometry::Planar { .. } => (a + b) / 2.0,
        LayerStackGeometry::Cylindrical { .. } => {
            2.0 * (b.powi(3) - a.powi(3)) / (3.0 * (b * b - a * a))
        }
        LayerStackGeometry::Spherical { .. } => {
            3.0 * (b.powi(4) - a.powi(4)) / (4.0 * (b.powi(3) - a.powi(3)))
        }
    }
}
fn body<'a>(
    sources: &[&'a PreparedPhaseConductivity<'a>],
    widths: &[f64],
    counts: &[usize],
    contacts: &[f64],
    geometry: LayerStackGeometry,
) -> PreparedPhaseLayeredThermalBody<'a> {
    let mut x = start(geometry);
    let layers: Vec<_> = sources
        .iter()
        .enumerate()
        .map(|(i, &conductivity)| {
            let mass_kg = volume(geometry, x, x + widths[i]);
            x += widths[i];
            PhaseLayerSource {
                id: 100 + i as u64,
                thickness_m: widths[i],
                mass_kg,
                cell_count: counts[i],
                conductivity,
            }
        })
        .collect();
    PreparedPhaseLayeredThermalBody::try_new(&layers, contacts, geometry, || Ok(())).unwrap()
}
fn ids() -> LayeredThermalIds {
    LayeredThermalIds {
        first_cell: 11,
        first_material: 21,
        first_link: 31,
        body_id: 41,
    }
}
fn initial<'a>(
    bound: &BoundPhaseLayeredThermalBody<'_, 'a>,
    h: &[f64],
    temperatures: [Option<f64>; 2],
) -> DepositionSimulation<'a> {
    let n = bound.prepared_body().cell_count();
    let events = bound.birth_events(0.0, h, temperatures, || Ok(())).unwrap();
    DepositionSimulation::new(
        events,
        &bound.materials(),
        DepositionBudget {
            max_cells: n,
            max_links: n - 1,
            max_events: 2 * n - 1,
            max_substeps: 100_000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn controls(dt: f64) -> ImplicitConductionControls {
    ImplicitConductionControls {
        maximum_dt_s: dt,
        minimum_dt_s: 1e-12,
        maximum_steps: 100_000,
        maximum_newton_iterations: 48,
        maximum_backtracks: 32,
        maximum_halvings: 20,
        maximum_residual_evaluations: 2_000_000,
        cell_tolerance_j: 1e-10,
        total_tolerance_j: 1e-8,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}, allowance={tolerance}"
    );
}
fn root(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) <= 0.0 && residual(hi) >= 0.0);
    for _ in 0..90 {
        let mid = f64::midpoint(lo, hi);
        if residual(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    f64::midpoint(lo, hi)
}

#[test]
fn exact_latent_births_and_physical_inventory_survive_the_native_adaptive_owner() {
    let card = material(1.0, 4.0);
    let curve = caloric(&card, 400.0, 100.0);
    let law = prepare(&card, &curve);
    let original_identity = law.identity();
    for shape in shapes() {
        let body = body(&[&law], &[0.8], &[4], &[], shape);
        let a = start(shape);
        let v = volume(shape, a, a + 0.8);
        close(body.total_volume_m3(), v, 1e-12);
        close(body.total_mass_kg(), v, 1e-12);
        for (cell, x) in body
            .cell_geometries()
            .zip(body.cell_temperature_positions_m())
        {
            close(cell.volume_m3, v / 4.0, 1e-12);
            close(
                volume(shape, cell.start_m[0], cell.end_m[0]),
                v / 4.0,
                1e-12,
            );
            close(x, centroid(shape, cell.start_m[0], cell.end_m[0]), 1e-12);
        }
        let bound = body.bind(ids(), [None; 2]).unwrap();
        let before = initial(&bound, &[10.0, 30.0, 50.0, 70.0], [None; 2]);
        assert!(before.uses_material_law(bound.materials()[0].0, &curve));
        let chain = bound.implicit_chain(&before, || Ok(())).unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(4).unwrap();
        let policy = AdaptiveImplicitConductionControls {
            step: controls(0.5),
            initial_dt_s: 0.5,
            maximum_trials: 16,
            enthalpy_absolute_tolerance_j_kg: 1e-6,
            temperature_absolute_tolerance_k: 1e-6,
            phase_fraction_absolute_tolerance: Some(1e-8),
            relative_tolerance: 1e-6,
        };
        let (next, report) = chain
            .advance_adaptive_to(&before, 1.0, policy, &mut workspace, || Ok(()))
            .unwrap();
        for (a, b) in before.cells().values().zip(next.cells().values()) {
            assert_eq!(a.specific_h_j_kg.to_bits(), b.specific_h_j_kg.to_bits());
            assert_eq!(next.temperature_k(b.id).unwrap(), 400.0);
        }
        close(report.solve.energy.external_heat_j, 0.0, 1e-12);
        close(report.solve.energy.residual_j, 0.0, 1e-12);
        assert_eq!(law.identity(), original_identity);
        assert!(std::ptr::eq(
            body.layers()[0].conductivity.caloric_curve(),
            &curve
        ));
    }
}

#[test]
fn two_original_phase_sources_and_contact_match_an_independent_latent_energy_balance() {
    let card_a = material(1.0, 4.0);
    let card_b = material(2.0, 2.0);
    let curve_a = caloric(&card_a, 400.0, 100.0);
    let curve_b = caloric(&card_b, 400.0, 100.0);
    let law_a = prepare(&card_a, &curve_a);
    let law_b = prepare(&card_b, &curve_b);
    for shape in shapes() {
        let body = body(&[&law_a, &law_b], &[0.2, 0.3], &[1, 1], &[0.05], shape);
        let bound = body.bind(ids(), [None; 2]).unwrap();
        let before = initial(&bound, &[80.0, -20.0], [None; 2]);
        let a = start(shape);
        let face = a + 0.2;
        let end = face + 0.3;
        let ma = volume(shape, a, face);
        let mb = volume(shape, face, end);
        let xa = centroid(shape, a, face);
        let xb = centroid(shape, face, end);
        let ra = metric(shape, xa, face);
        let rb = metric(shape, face, xb) / 2.0;
        let contact = 0.05 / area(shape, face);
        let dt = 0.01;
        let expected = root(0.0, (80.0 * ma).min(20.0 * mb), |q| {
            let k_a = 1.0 + 3.0 * (80.0 - q / ma) / 100.0;
            q - dt * (20.0 - q / mb) / (ra / k_a + contact + rb)
        });
        let frozen_g = 1.0 / (ra / 3.4 + contact + rb);
        let frozen = dt * frozen_g * 20.0 / (1.0 + dt * frozen_g / mb);
        assert!((expected - frozen).abs() > 1e-4);
        let chain = bound.implicit_chain(&before, || Ok(())).unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        let (next, report) = chain
            .advance_to(&before, dt, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        let states: Vec<_> = next.cells().values().collect();
        close(states[0].specific_h_j_kg, 80.0 - expected / ma, 2e-8);
        close(states[1].specific_h_j_kg, -20.0 + expected / mb, 2e-8);
        close(next.temperature_k(states[0].id).unwrap(), 400.0, 1e-12);
        close(
            next.temperature_k(states[1].id).unwrap(),
            380.0 + expected / mb,
            2e-8,
        );
        close(
            ma * (states[0].specific_h_j_kg - 80.0) + mb * (states[1].specific_h_j_kg + 20.0),
            0.0,
            1e-8,
        );
        close(report.energy.external_heat_j, 0.0, 1e-12);
        close(report.energy.residual_j, 0.0, 1e-8);
    }
}

#[test]
fn real_face_reservoir_and_phase_dependent_half_cell_transport_close_one_latent_store() {
    let card = material(1.0, 4.0);
    let curve = caloric(&card, 400.0, 100.0);
    let law = prepare(&card, &curve);
    for shape in shapes() {
        let body = body(&[&law], &[0.5], &[1], &[], shape);
        let bound = body.bind(ids(), [None, Some(0.1)]).unwrap();
        let before = initial(&bound, &[20.0], [None, Some(450.0)]);
        let a = start(shape);
        let b = a + 0.5;
        let mass = volume(shape, a, b);
        let bulk = metric(shape, centroid(shape, a, b), b);
        let film = 0.1 / area(shape, b);
        let dt = 0.01;
        let expected = root(0.0, 80.0 * mass, |q| {
            let k = 1.0 + 3.0 * (20.0 + q / mass) / 100.0;
            q - dt * 50.0 / (bulk / k + film)
        });
        let chain = bound.implicit_chain(&before, || Ok(())).unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        let (next, report) = chain
            .advance_to(&before, dt, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        let cell = next.cells().values().next().unwrap();
        close(cell.specific_h_j_kg, 20.0 + expected / mass, 2e-8);
        close(next.temperature_k(cell.id).unwrap(), 400.0, 1e-12);
        close(report.energy.external_heat_j, expected, 1e-8);
        close(report.energy.residual_j, 0.0, 1e-8);
        close(
            cell.boundary.reservoir.unwrap().area_m2,
            area(shape, b),
            1e-12,
        );
    }
}

#[test]
fn source_bound_freezing_front_and_heat_refine_to_the_stefan_similarity_solution() {
    // rho=Cp=k=L=1, lambda=1/2. The only applied drive is a cold physical
    // face; the exact front s(t)=sqrt(t) is used solely as an observation.
    let card = material(1.0, 1.0);
    let curve = caloric(&card, 300.0, 1.0);
    let law = prepare(&card, &curve);
    let lambda = 0.5_f64;
    let delta = lambda * PI.sqrt() * (lambda * lambda).exp() * fs_math::det::erf(lambda);
    let time = 0.16_f64;
    let exact_front = time.sqrt();
    let exact_heat = -(0.25_f64).exp() * time.sqrt();
    let shape = LayerStackGeometry::Planar { area_m2: 1.0 };
    let mut errors = Vec::new();
    for n in [16, 32, 64] {
        let dx = 1.0 / n as f64;
        let body = body(&[&law], &[1.0], &[n], &[], shape);
        let bound = body.bind(ids(), [Some(0.0), None]).unwrap();
        let before = initial(&bound, &vec![1.0; n], [Some(300.0 - delta), None]);
        let chain = bound.implicit_chain(&before, || Ok(())).unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(n).unwrap();
        let (next, report) = chain
            .advance_to(
                &before,
                time,
                controls(0.25 * dx * dx),
                &mut workspace,
                || Ok(()),
            )
            .unwrap();
        let mut frozen_length = 0.0;
        let mut temperature_error = 0.0;
        for (i, cell) in next.cells().values().enumerate() {
            let state = law
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .unwrap();
            frozen_length += dx * (1.0 - state.liquid_mass_fraction());
            let x = (i as f64 + 0.5) * dx;
            let exact = if x < exact_front {
                300.0 - delta
                    + delta * fs_math::det::erf(x / (2.0 * time.sqrt())) / fs_math::det::erf(lambda)
            } else {
                300.0
            };
            temperature_error += dx * (state.temperature_k() - exact).powi(2);
            assert!(
                state.temperature_k() >= 300.0 - delta - 1e-9
                    && state.temperature_k() <= 300.0 + 1e-9
            );
        }
        close(report.energy.residual_j, 0.0, 1e-8);
        errors.push([
            (frozen_length - exact_front).abs(),
            (report.energy.external_heat_j - exact_heat).abs(),
            temperature_error.sqrt(),
        ]);
    }
    for metric in 0..3 {
        assert!(errors[2][metric] < errors[0][metric], "{errors:?}");
    }
    assert!(
        errors[2][0] < 0.008 && errors[2][1] < 0.01 && errors[2][2] < 0.005,
        "{errors:?}"
    );
}

#[test]
fn source_domain_binding_and_cancellation_refuse_without_changing_phase_history() {
    let card = material(1.0, 4.0);
    let curve = caloric(&card, 400.0, 100.0);
    let law = prepare(&card, &curve);
    let shape = LayerStackGeometry::Planar { area_m2: 1.0 };
    let body = body(&[&law], &[1.0], &[2], &[], shape);
    let bound = body.bind(ids(), [None; 2]).unwrap();
    for h in [vec![0.0], vec![f64::NAN, 20.0], vec![-101.0, 20.0]] {
        assert!(bound.birth_events(0.0, &h, [None; 2], || Ok(())).is_err());
    }
    assert!(matches!(
        bound.birth_events(0.0, &[10.0, 20.0], [None; 2], || Err(
            DepositionError::Cancelled
        )),
        Err(DepositionError::Cancelled)
    ));
    let before = initial(&bound, &[10.0, 20.0], [None; 2]);
    let original = before.cells().clone();
    let chain = bound.implicit_chain(&before, || Ok(())).unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    assert!(
        chain
            .advance_to(&before, 0.1, controls(0.1), &mut workspace, || Err(
                DepositionError::Cancelled
            ))
            .is_err()
    );
    assert_eq!(before.cells(), &original);
    let one = self::body(&[&law], &[1.0], &[1], &[], shape);
    let both_faces = one.bind(ids(), [Some(0.0), Some(0.0)]).unwrap();
    assert!(
        both_faces
            .birth_events(0.0, &[20.0], [Some(400.0), None], || Ok(()))
            .is_err()
    );
    assert!(body.bind(ids(), [Some(-0.1), None]).is_err());
    let mut events = bound
        .birth_events(0.0, &[10.0, 20.0], [None; 2], || Ok(()))
        .unwrap();
    if let Event::Birth(cell) = &mut events[0].event {
        cell.mass_kg *= 2.0;
    }
    let foreign = DepositionSimulation::new(
        events,
        &bound.materials(),
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 10,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    assert!(bound.conduction_observer(&foreign, || Ok(())).is_err());
    assert_eq!(body.layers()[0].conductivity.identity(), law.identity());
}

#[path = "phase_layered_transport/boundaries.rs"]
mod boundaries;

#[path = "phase_layered_transport/solid.rs"]
mod solid;
