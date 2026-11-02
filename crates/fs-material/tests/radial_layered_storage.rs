//! Source-backed shell inventory and transport on the original thermal owner.
//! Authored numerical materials; these are implementation, not measured-data tests.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::{
    conduction::LayerConductionConfig,
    transport::{
        ThermalTransportError,
        stack::{
            LayerStackGeometry, LayerStackSourceQuery, PreparedLayerStackConduction,
            transient::{
                BoundLayeredThermalBody, LayeredThermalIds, PreparedLayeredThermalBody,
                TransientLayerSourceQuery,
                implicit::{
                    LayeredThermalAdaptiveControls, LayeredThermalSolveControls,
                    LayeredThermalStep, LayeredThermalWorkspace,
                },
            },
        },
        storage::{HeatStorageAmount, HeatStorageComponentQuery},
    },
};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::{
    DepositionBudget, DepositionError, DepositionSimulation, Event, StepConfig,
};

const K_DIMS: Dims = Dims([1, 1, -3, -1, 0, 0]);
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const CV_DIMS: Dims = Dims([-1, 1, -2, -1, 0, 0]);
const RHO: f64 = 1000.0;
const CP: f64 = 200.0;
const PI: f64 = core::f64::consts::PI;

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            t,
        )
        .unwrap()
}
fn material(name: &str, conductivity: [(f64, f64); 3]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (
            "k",
            K_DIMS,
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: conductivity.to_vec(),
                dims: K_DIMS,
            },
        ),
        (
            "cp",
            CP_DIMS,
            PropertyValue::Scalar {
                value: CP,
                dims: CP_DIMS,
            },
        ),
        (
            "cv",
            CV_DIMS,
            PropertyValue::Scalar {
                value: RHO * CP,
                dims: CV_DIMS,
            },
        ),
    ] {
        let interpolation = if matches!(&value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    200.0,
                    800.0,
                ),
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored radial finite-volume oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "fixed reference shell".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn constant(name: &str, k: f64) -> MaterialCard {
    material(name, [(200.0, k), (450.0, k), (800.0, k)])
}
fn shape(cylinder: bool, inner: f64) -> LayerStackGeometry {
    if cylinder {
        LayerStackGeometry::Cylindrical {
            inner_radius_m: inner,
            length_m: 0.7,
        }
    } else {
        LayerStackGeometry::Spherical {
            inner_radius_m: inner,
        }
    }
}
fn inner(geometry: LayerStackGeometry) -> f64 {
    match geometry {
        LayerStackGeometry::Planar { .. } => 0.0,
        LayerStackGeometry::Cylindrical { inner_radius_m, .. }
        | LayerStackGeometry::Spherical { inner_radius_m } => inner_radius_m,
    }
}
fn area(geometry: LayerStackGeometry, r: f64) -> f64 {
    match geometry {
        LayerStackGeometry::Planar { area_m2 } => area_m2,
        LayerStackGeometry::Cylindrical { length_m, .. } => 2.0 * PI * r * length_m,
        LayerStackGeometry::Spherical { .. } => 4.0 * PI * r * r,
    }
}
fn volume(geometry: LayerStackGeometry, a: f64, b: f64) -> f64 {
    match geometry {
        LayerStackGeometry::Planar { area_m2 } => area_m2 * (b - a),
        LayerStackGeometry::Cylindrical { length_m, .. } => PI * length_m * (b * b - a * a),
        LayerStackGeometry::Spherical { .. } => 4.0 * PI / 3.0 * (b.powi(3) - a.powi(3)),
    }
}
fn metric(geometry: LayerStackGeometry, a: f64, b: f64) -> f64 {
    match geometry {
        LayerStackGeometry::Planar { area_m2 } => (b - a) / area_m2,
        LayerStackGeometry::Cylindrical { length_m, .. } => (b / a).ln() / (2.0 * PI * length_m),
        LayerStackGeometry::Spherical { .. } => (1.0 / a - 1.0 / b) / (4.0 * PI),
    }
}
fn prepare(
    cards: &[&MaterialCard],
    widths: &[f64],
    counts: &[usize],
    contacts: &[f64],
    geometry: LayerStackGeometry,
    volumetric: bool,
) -> Result<PreparedLayeredThermalBody, ThermalTransportError> {
    let lower = point(200.0);
    let upper = point(800.0);
    let k = PropertyKey::new("k", K_DIMS);
    let capacity = PropertyKey::new(
        if volumetric { "cv" } else { "cp" },
        if volumetric { CV_DIMS } else { CP_DIMS },
    );
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let mut start = inner(geometry);
    let queries: Vec<_> = cards
        .iter()
        .enumerate()
        .map(|(i, &card)| {
            let end = start + widths[i];
            let v = volume(geometry, start, end);
            start = end;
            TransientLayerSourceQuery {
                conduction: LayerStackSourceQuery {
                    id: 100 + i as u64,
                    card,
                    conductivity: &k,
                    thickness_m: widths[i],
                    temperature_axis: "T",
                    lower: &lower,
                    upper: &upper,
                    selection,
                },
                capacity: HeatStorageComponentQuery {
                    id: 200 + i as u64,
                    card,
                    capacity: &capacity,
                    amount: if volumetric {
                        HeatStorageAmount::ReferenceVolumeM3(v)
                    } else {
                        HeatStorageAmount::MassKg(RHO * v)
                    },
                    temperature_axis: "T",
                    lower: &lower,
                    upper: &upper,
                    selection,
                },
                mass_kg: RHO * v,
                cell_count: counts[i],
            }
        })
        .collect();
    PreparedLayeredThermalBody::try_new_geometric(&queries, contacts, geometry, 300.0)
}
fn ids() -> LayeredThermalIds {
    LayeredThermalIds {
        first_cell: 10,
        first_material: 30,
        first_link: 50,
        body_id: 70,
    }
}
fn conduction() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 96,
        absolute_heat_tolerance_w: 1e-10,
        relative_heat_tolerance: 1e-12,
    }
}
fn initial<'a>(
    bound: &'a BoundLayeredThermalBody<'_>,
    temperatures: &[f64],
    sources: &[f64],
    htc: f64,
    ambient: f64,
) -> DepositionSimulation<'a> {
    let n = bound.prepared_body().cell_count();
    let mut events = bound
        .birth_events(0.0, temperatures, conduction(), || Ok(()))
        .unwrap();
    let mut i = 0;
    for event in &mut events {
        if let Event::Birth(cell) = &mut event.event {
            cell.boundary.source_w = sources[i];
            cell.boundary.htc_w_m2_k = htc;
            cell.boundary.ambient_k = ambient;
            i += 1;
        }
    }
    let materials: Vec<_> = bound.materials().into_iter().collect();
    DepositionSimulation::new(
        events,
        &materials,
        DepositionBudget {
            max_cells: n,
            max_links: n - 1,
            max_events: 2 * n - 1,
            max_substeps: 64,
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
    .0
}
fn close(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "{actual:.16e} != {expected:.16e}, allowance={absolute}"
    );
}

#[test]
fn g1_radial_stack_uses_actual_contact_area_and_original_conductivity_integral() {
    let a = constant("inner", 20.0);
    let b = constant("outer", 0.6);
    let nonlinear = material("affine", [(200.0, 4.0), (450.0, 6.5), (800.0, 10.0)]);
    for cylinder in [false, true] {
        let geometry = shape(cylinder, 0.01);
        let body = prepare(&[&a, &b], &[0.01, 0.02], &[1, 1], &[0.03], geometry, false).unwrap();
        let stack = body.conduction();
        let resistance = metric(geometry, 0.01, 0.02) / 20.0
            + 0.03 / area(geometry, 0.02)
            + metric(geometry, 0.02, 0.04) / 0.6;
        for temperatures in [[600.0, 300.0], [300.0, 600.0], [450.0; 2]] {
            let response = stack
                .evaluate(temperatures, conduction(), || Ok(()))
                .unwrap();
            let q = (temperatures[0] - temperatures[1]) / resistance;
            close(response.power_a_to_b_w, q, 1e-8);
            close(response.endpoint_derivative_w_k[0], 1.0 / resistance, 1e-10);
            close(
                response.endpoint_derivative_w_k[1],
                -1.0 / resistance,
                1e-10,
            );
            let faces = response.layer_temperatures_k();
            close(
                faces[0][1] - faces[1][0],
                q * 0.03 / area(geometry, 0.02),
                1e-8,
            );
        }
        assert_eq!(stack.layers()[0].thickness_m(), 0.01);
        assert_eq!(stack.contact_resistances_m2_k_w(), &[0.03]);
        let nonlinear_body = prepare(&[&nonlinear], &[0.03], &[1], &[], geometry, false).unwrap();
        let response = nonlinear_body
            .conduction()
            .evaluate([600.0, 300.0], conduction(), || Ok(()))
            .unwrap();
        let integral = 2.0 * 300.0 + 0.005 * (600.0_f64.powi(2) - 300.0_f64.powi(2));
        close(
            response.power_a_to_b_w,
            integral / metric(geometry, 0.01, 0.04),
            1e-8,
        );
    }
}

#[test]
fn g1_equal_volume_shells_keep_real_inventory_centroids_and_exposed_faces() {
    let card = constant("inventory", 8.0);
    for cylinder in [false, true] {
        let geometry = shape(cylinder, 0.01);
        let expected_volume = volume(geometry, 0.01, 0.04);
        for volumetric in [false, true] {
            let body = prepare(&[&card], &[0.03], &[7], &[], geometry, volumetric).unwrap();
            close(
                body.total_volume_m3(),
                expected_volume,
                2e-15 * expected_volume,
            );
            close(
                body.total_mass_kg(),
                RHO * expected_volume,
                2e-15 * RHO * expected_volume,
            );
            assert_eq!(body.storages().len(), 1);
            let laws = body.thermal_laws(false).unwrap();
            let bound = body.bind(&laws, ids()).unwrap();
            let state = initial(&bound, &[450.0; 7], &[0.0; 7], 0.0, 300.0);
            let d = if cylinder { 2 } else { 3 };
            for (i, (cell, center)) in body
                .cell_geometries()
                .zip(body.cell_temperature_positions_m())
                .enumerate()
            {
                let [a, b] = [cell.start_m[0], cell.end_m[0]];
                close(
                    volume(geometry, a, b),
                    expected_volume / 7.0,
                    2e-14 * expected_volume,
                );
                close(
                    cell.volume_m3,
                    expected_volume / 7.0,
                    2e-15 * expected_volume,
                );
                let centroid = d as f64 / (d + 1) as f64 * (b.powi(d + 1) - a.powi(d + 1))
                    / (b.powi(d) - a.powi(d));
                close(center, centroid, 2e-15);
                close(
                    cell.surface_area_m2,
                    area(geometry, a) + area(geometry, b),
                    2e-15,
                );
                let exposed = if i == 0 {
                    area(geometry, 0.01)
                } else if i == 6 {
                    area(geometry, 0.04)
                } else {
                    0.0
                };
                close(
                    state.exposed_area_m2(bound.cell_ids()[i]).unwrap(),
                    exposed,
                    2e-15,
                );
            }
            close(state.energy_j(), RHO * expected_volume * CP * 150.0, 2e-10);
        }
    }
}

#[test]
fn g1_radial_nonlinear_flux_tangent_and_implicit_energy_use_one_owner() {
    let card = material("nonlinear", [(200.0, 4.0), (450.0, 9.0), (800.0, 10.0)]);
    for cylinder in [false, true] {
        let body = prepare(
            &[&card, &card],
            &[0.01, 0.02],
            &[2, 3],
            &[0.002],
            shape(cylinder, 0.01),
            false,
        )
        .unwrap();
        let laws = body.thermal_laws(false).unwrap();
        let bound = body.bind(&laws, ids()).unwrap();
        let old = initial(
            &bound,
            &[600.0, 500.0, 430.0, 370.0, 300.0],
            &[0.0; 5],
            0.0,
            300.0,
        );
        let old_cells = old.cells().clone();
        let step = LayeredThermalStep::new(&bound, &old, 0.1, conduction(), || Ok(())).unwrap();
        let energy = [0.02, -0.01, 0.01, 0.02, -0.01];
        let direction = [0.2, -0.3, 0.5, -0.1, 0.4];
        let trial = step.evaluate(&energy, || Ok(())).unwrap();
        let mut action = [0.0; 5];
        trial
            .payload()
            .apply_into(&direction, &mut action, || Ok(()))
            .unwrap();
        let epsilon = 1e-4;
        let plus = step
            .evaluate(
                &std::array::from_fn::<_, 5, _>(|i| energy[i] + epsilon * direction[i]),
                || Ok(()),
            )
            .unwrap();
        let minus = step
            .evaluate(
                &std::array::from_fn::<_, 5, _>(|i| energy[i] - epsilon * direction[i]),
                || Ok(()),
            )
            .unwrap();
        for i in 0..5 {
            close(
                action[i],
                (plus.residual_j()[i] - minus.residual_j()[i]) / (2.0 * epsilon),
                2e-6 * action[i].abs().max(1.0),
            );
        }
        let mut workspace = LayeredThermalWorkspace::new(5).unwrap();
        let (next, _, report) = step
            .solve(
                LayeredThermalSolveControls::default(),
                &mut workspace,
                || Ok(()),
            )
            .unwrap();
        close(next.energy_j(), old.energy_j(), 2e-8);
        close(report.energy.external_heat_j, 0.0, 0.0);
        assert!(report.energy.residual_j.abs() <= 1e-8);
        assert!(next.temperature_k(bound.cell_ids()[0]).unwrap() < 600.0);
        assert!(next.temperature_k(bound.cell_ids()[4]).unwrap() > 300.0);
        assert!(matches!(
            step.evaluate(&energy, || Err(DepositionError::Cancelled)),
            Err(DepositionError::Cancelled)
        ));
        assert_eq!(old.cells(), &old_cells);
    }
}

#[test]
fn g1_radial_boundary_and_volumetric_source_energy_use_physical_measures() {
    let card = constant("boundary", 8.0);
    for cylinder in [false, true] {
        let geometry = shape(cylinder, 0.01);
        let body = prepare(&[&card], &[0.03], &[1], &[], geometry, true).unwrap();
        let laws = body.thermal_laws(false).unwrap();
        let bound = body.bind(&laws, ids()).unwrap();
        let volume = body.total_volume_m3();
        let source = 10000.0 * volume;
        let old = initial(&bound, &[400.0], &[source], 30.0, 320.0);
        let dt = 0.5;
        let g = 30.0 * (area(geometry, 0.01) + area(geometry, 0.04));
        let c = RHO * CP * volume;
        let expected = (c * 400.0 + dt * (g * 320.0 + source)) / (c + dt * g);
        let step = LayeredThermalStep::new(&bound, &old, dt, conduction(), || Ok(())).unwrap();
        let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
        let (next, _, report) = step
            .solve(
                LayeredThermalSolveControls::default(),
                &mut workspace,
                || Ok(()),
            )
            .unwrap();
        close(
            next.temperature_k(bound.cell_ids()[0]).unwrap(),
            expected,
            2e-10,
        );
        close(
            report.energy.external_heat_j,
            dt * (source - g * (expected - 320.0)),
            2e-9,
        );
        close(
            report.energy.stored_energy_change_j,
            c * (expected - 400.0),
            2e-9,
        );
    }
}

#[test]
fn g1_existing_adaptive_owner_conserves_distributed_radial_source_heat() {
    let card = constant("adaptive", 8.0);
    for cylinder in [false, true] {
        let body = prepare(&[&card], &[0.03], &[4], &[], shape(cylinder, 0.01), true).unwrap();
        let laws = body.thermal_laws(true).unwrap();
        let bound = body.bind(&laws, ids()).unwrap();
        let sources: Vec<_> = body
            .cell_geometries()
            .map(|c| 10000.0 * c.volume_m3)
            .collect();
        let old = initial(&bound, &[350.0; 4], &sources, 0.0, 300.0);
        let mut controls = LayeredThermalAdaptiveControls::default();
        controls.initial_dt_s = 0.5;
        controls.step.maximum_dt_s = 0.5;
        controls.temperature_absolute_tolerance_k = 1e-5;
        let mut workspace = LayeredThermalWorkspace::new(4).unwrap();
        let (next, report) = bound
            .advance_adaptive_to(&old, 1.0, conduction(), controls, &mut workspace, || Ok(()))
            .unwrap();
        for id in bound.cell_ids() {
            close(
                next.temperature_k(id).unwrap(),
                350.0 + 10000.0 / (RHO * CP),
                2e-9,
            );
        }
        close(
            report.solve.energy.external_heat_j,
            10000.0 * body.total_volume_m3(),
            2e-9,
        );
        assert!(report.solve.energy.residual_j.abs() <= 1e-8);
        assert!(report.accepted_comparisons > 0);
    }
}

#[test]
fn g1_radial_diffusion_refines_toward_manufactured_parabolic_temperature() {
    // T=300+a*r²+2*d*alpha*a*t exactly satisfies radial diffusion. Constant
    // boundary powers -k*A*dT/dr supply the entire known external heat. Time
    // dependence is affine, so backward Euler introduces no temporal error in
    // this reference; the computed discrepancy tests the spatial reduction.
    let card = constant("manufactured", 8.0);
    for cylinder in [false, true] {
        let geometry = shape(cylinder, 0.02);
        let d = if cylinder { 2.0 } else { 3.0 };
        let mut errors = Vec::new();
        for n in [8, 32, 128] {
            let body = prepare(&[&card], &[0.08], &[n], &[], geometry, false).unwrap();
            let centers: Vec<_> = body.cell_temperature_positions_m().collect();
            let temperatures: Vec<_> = centers.iter().map(|r| 300.0 + 1000.0 * r * r).collect();
            let mut sources = vec![0.0; n];
            sources[0] = -2.0 * 8.0 * 1000.0 * 0.02 * area(geometry, 0.02);
            sources[n - 1] = 2.0 * 8.0 * 1000.0 * 0.1 * area(geometry, 0.1);
            let laws = body.thermal_laws(false).unwrap();
            let bound = body.bind(&laws, ids()).unwrap();
            let old = initial(&bound, &temperatures, &sources, 0.0, 300.0);
            let step = LayeredThermalStep::new(&bound, &old, 1.0, conduction(), || Ok(())).unwrap();
            let mut workspace = LayeredThermalWorkspace::new(n).unwrap();
            let (next, _, report) = step
                .solve(
                    LayeredThermalSolveControls::default(),
                    &mut workspace,
                    || Ok(()),
                )
                .unwrap();
            let increase = 2.0 * d * 8.0 / (RHO * CP) * 1000.0;
            let error = bound
                .cell_ids()
                .iter()
                .zip(&temperatures)
                .map(|(&id, &t)| (next.temperature_k(id).unwrap() - t - increase).abs())
                .fold(0.0_f64, f64::max);
            errors.push(error);
            close(
                report.energy.external_heat_j,
                sources.iter().sum::<f64>(),
                2e-8,
            );
            assert!(report.energy.residual_j.abs() <= 1e-8);
        }
        assert!(errors[1] < 0.7 * errors[0], "{errors:?}");
        assert!(errors[2] < 0.7 * errors[1], "{errors:?}");
        assert!(errors[2] < 0.007, "{errors:?}");
    }
}

#[test]
fn g0_planar_constructor_and_source_identity_remain_compatible() {
    let card = constant("planar", 8.0);
    let lower = point(200.0);
    let upper = point(800.0);
    let key = PropertyKey::new("k", K_DIMS);
    let queries = [LayerStackSourceQuery {
        id: 1,
        card: &card,
        conductivity: &key,
        thickness_m: 0.03,
        temperature_axis: "T",
        lower: &lower,
        upper: &upper,
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }];
    let old = PreparedLayerStackConduction::try_new(&queries, &[], 0.7).unwrap();
    let explicit = PreparedLayerStackConduction::try_new_geometric(
        &queries,
        &[],
        LayerStackGeometry::Planar { area_m2: 0.7 },
    )
    .unwrap();
    assert_eq!(old, explicit);
    close(
        old.evaluate([600.0, 300.0], conduction(), || Ok(()))
            .unwrap()
            .power_a_to_b_w,
        0.7 * 8.0 * 300.0 / 0.03,
        1e-8,
    );
    assert!(
        PreparedLayerStackConduction::try_new_geometric(&queries, &[], shape(false, 0.0)).is_err()
    );
    assert!(prepare(&[&card], &[0.03], &[513], &[], shape(true, 0.01), false).is_err());
    assert!(matches!(
        explicit.evaluate([900.0, 300.0], conduction(), || Ok(())),
        Err(DepositionError::MaterialDomain)
    ));
    assert!(
        explicit
            .evaluate([600.0, 300.0], conduction(), || Err(
                DepositionError::Cancelled
            ))
            .is_err()
    );
}

#[test]
fn g3_radial_capacity_refuses_planar_proxy_volume_and_unresolved_shell_cells() {
    let card = constant("admission", 8.0);
    let lower = point(200.0);
    let upper = point(800.0);
    let conductivity = PropertyKey::new("k", K_DIMS);
    let capacity = PropertyKey::new("cv", CV_DIMS);
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let geometry = shape(false, 0.01);
    let physical_volume = volume(geometry, 0.01, 0.04);
    let planar_proxy = area(geometry, 0.01) * 0.03;
    assert!(physical_volume > 5.0 * planar_proxy);
    let mut query = TransientLayerSourceQuery {
        conduction: LayerStackSourceQuery {
            id: 1,
            card: &card,
            conductivity: &conductivity,
            thickness_m: 0.03,
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
            selection,
        },
        capacity: HeatStorageComponentQuery {
            id: 2,
            card: &card,
            capacity: &capacity,
            amount: HeatStorageAmount::ReferenceVolumeM3(planar_proxy),
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
            selection,
        },
        mass_kg: RHO * physical_volume,
        cell_count: 3,
    };
    assert!(PreparedLayeredThermalBody::try_new_geometric(&[query], &[], geometry, 300.0).is_err());
    query.capacity.amount = HeatStorageAmount::ReferenceVolumeM3(physical_volume);
    assert!(PreparedLayeredThermalBody::try_new_geometric(&[query], &[], geometry, 300.0).is_ok());
    for cylinder in [false, true] {
        assert!(
            prepare(
                &[&card],
                &[f64::EPSILON],
                &[2],
                &[],
                shape(cylinder, 1.0),
                false
            )
            .is_err()
        );
    }
}
