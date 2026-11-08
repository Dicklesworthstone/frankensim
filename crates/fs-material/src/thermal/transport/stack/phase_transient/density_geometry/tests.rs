//! Manufactured liquid sources exercise the real thermal timeline and radial integrals.
use super::*;
use crate::{
    phase::{
        EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacitySourceQuery,
        ResolvedHeatCapacityCurve, SolidLiquidPhase,
    },
    state_point::THERMAL_CONDUCTIVITY_DIMS,
    thermal::source::ThermalMaterialPoint,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_thermal::enthalpy::{DepositionBudget, StepConfig};

const PI: f64 = core::f64::consts::PI;
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
fn point(t: f64, f: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
        .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), 100_000.0)
        .unwrap()
        .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), f)
        .unwrap()
}
fn material(k: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, knots) in [
        ("cp", CP, vec![(300.0, 10.0), (400.0, 10.0)]),
        ("rho", RHO, vec![(300.0, 1000.0), (400.0, 500.0)]),
        ("k", THERMAL_CONDUCTIVITY_DIMS, vec![(280.0, k), (420.0, k)]),
    ] {
        let lo = knots[0].0;
        let hi = knots[1].0;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots,
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), lo, hi)
                    .with_quantity(
                        "p",
                        QuantitySpec::dimensional(Pressure::DIMS),
                        90_000.0,
                        110_000.0,
                    )
                    .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured liquid; large density slope is a numerical oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic spherical liquid".into(),
            phase: "declared caloric chart".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new_exact(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("cp", CP),
            density: &PropertyKey::new("rho", RHO),
            temperature_axis: "T",
            lower: &point(300.0, 1.0),
            upper: &point(400.0, 1.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        0.0,
    )
    .unwrap()
}
fn prepare<'a>(
    card: &MaterialCard,
    chart: &'a EquilibriumEnthalpyPhaseCurve,
) -> PreparedPhaseConductivity<'a> {
    PreparedPhaseConductivity::new(
        chart,
        ThermalMaterialPoint {
            material: card,
            context: &point(
                chart.knots()[0].temperature_k,
                chart.knots()[0].liquid_mass_fraction,
            ),
            temperature_axis: "T",
            liquid_fraction_axis: Some("f"),
        },
        &PropertyKey::new("k", THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        [
            chart.knots()[0].specific_enthalpy_j_kg,
            chart.knots().last().unwrap().specific_enthalpy_j_kg,
        ],
        || Ok(()),
    )
    .unwrap()
}
fn ids() -> LayeredThermalIds {
    LayeredThermalIds {
        first_cell: 11,
        first_material: 21,
        first_link: 31,
        body_id: 41,
    }
}
fn volume(r: f64) -> f64 {
    4.0 * PI * r.powi(3) / 3.0
}
fn centroid(a: f64, b: f64) -> f64 {
    3.0 * (b.powi(4) - a.powi(4)) / (4.0 * (b.powi(3) - a.powi(3)))
}
fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-9,
    }
}
fn sphere<'a>(
    law: &'a PreparedPhaseConductivity<'a>,
    mass: f64,
    count: usize,
) -> PreparedPhaseLayeredThermalBody<'a> {
    PreparedPhaseLayeredThermalBody::try_new_solid_spherical(
        &[PhaseLayerSource {
            id: 1,
            thickness_m: 0.1,
            mass_kg: mass,
            cell_count: count,
            conductivity: law,
        }],
        &[],
        || Ok(()),
    )
    .unwrap()
}
fn initial<'a>(
    bound: &BoundPhaseLayeredThermalBody<'_, 'a>,
    h: &[f64],
    bath: Option<f64>,
    mut edit: impl FnMut(&mut DepositedCell),
) -> DepositionSimulation<'a> {
    let mut events = bound.birth_events(0.0, h, [None, bath], || Ok(())).unwrap();
    for timed in &mut events {
        if let Event::Birth(cell) = &mut timed.event {
            edit(cell);
        }
    }
    DepositionSimulation::new(
        events,
        &bound.materials(),
        DepositionBudget {
            max_cells: h.len(),
            max_links: h.len() - 1,
            max_events: 2 * h.len() - 1,
            max_substeps: 10_000,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}, tolerance={tolerance}"
    );
}

#[test]
fn uniform_source_heating_conserves_inventory_and_expands_actual_radial_transport() {
    let card = material(1.0);
    let law = caloric(&card);
    let k = prepare(&card, law.curve());
    let mass = 1000.0 * volume(0.1);
    let body = sphere(&k, mass, 3);
    let bound = body.bind(ids(), [None; 2]).unwrap();
    let before = initial(&bound, &[0.0; 3], None, |cell| {
        cell.boundary.source_w = 500.0 * cell.mass_kg
    });
    let observer = bound
        .spherical_liquid_density_observer(&before, || Ok(()))
        .unwrap();
    let mut cold = observer.clone();
    let g_cold = cold.conductances(&before, &mut || Ok(())).unwrap();
    let (after, mut warm, audit) = before
        .advance_observed(&observer, 1.0, config(0.25), || Ok(()))
        .unwrap();
    let shape = warm.current_geometry().clone();
    let scale = (1000.0_f64 / 750.0).cbrt();
    close(shape.total_mass_kg, mass, 1e-14);
    close(shape.total_volume_m3, mass / 750.0, 1e-16);
    close(shape.outer_radius_m, 0.1 * scale, 1e-15);
    close(shape.outer_area_m2, 4.0 * PI * (0.1 * scale).powi(2), 1e-15);
    close(audit.external_heat_j, 500.0 * mass, 1e-9);
    close(audit.stored_energy_change_j, audit.external_heat_j, 1e-9);
    close(audit.residual_j, 0.0, 1e-9);
    assert_eq!(audit.substeps, 4);
    let g_warm = warm.conductances(&after, &mut || Ok(())).unwrap();
    for id in g_cold.keys() {
        close(g_warm[id] / g_cold[id], scale, 2e-13);
    }
    for (index, cell) in after.cells().values().enumerate() {
        let old = &before.cells()[&cell.id];
        assert_eq!(cell.mass_kg, old.mass_kg);
        assert_eq!(cell.geometry, old.geometry);
        assert!(after.uses_material_law(cell.material_id, law.curve()));
        close(after.temperature_k(cell.id).unwrap(), 350.0, 1e-11);
        close(shape.cells[index].density_kg_m3, 750.0, 1e-10);
        let enclosed = (index + 1) as f64 * mass / (3.0 * 750.0);
        let radius = (3.0 * enclosed / (4.0 * PI)).cbrt();
        close(shape.cells[index].outer_radius_m, radius, 1e-15);
    }
    assert!(!law.curve().source_receipts().is_empty());
    for receipt in law.curve().source_receipts() {
        card.claims().verify_receipt(receipt).unwrap();
    }
}

#[test]
fn unequal_liquid_inventories_update_the_next_substep_resistance_and_share_one_heat_transfer() {
    let (card_a, card_b) = (material(1.0), material(2.0));
    let (law_a, law_b) = (caloric(&card_a), caloric(&card_b));
    let (ka, kb) = (
        prepare(&card_a, law_a.curve()),
        prepare(&card_b, law_b.curve()),
    );
    let ma = 750.0 * volume(0.04);
    let mb = 1000.0 * (volume(0.1) - volume(0.04));
    let contact = 0.002;
    let body = PreparedPhaseLayeredThermalBody::try_new_solid_spherical(
        &[
            PhaseLayerSource {
                id: 1,
                thickness_m: 0.04,
                mass_kg: ma,
                cell_count: 1,
                conductivity: &ka,
            },
            PhaseLayerSource {
                id: 2,
                thickness_m: 0.06,
                mass_kg: mb,
                cell_count: 1,
                conductivity: &kb,
            },
        ],
        &[contact],
        || Ok(()),
    )
    .unwrap();
    let bound = body.bind(ids(), [None; 2]).unwrap();
    let before = initial(&bound, &[500.0, 0.0], None, |_| {});
    let observer = bound
        .spherical_liquid_density_observer(&before, || Ok(()))
        .unwrap();
    let mut fixed = observer.clone();
    let original_g = fixed.conductances(&before, &mut || Ok(())).unwrap()[&31];
    let dt = 0.25;
    let (middle, mid_observer, _) = before
        .advance_observed(&observer, dt, config(dt), || Ok(()))
        .unwrap();
    let (ta, tb) = (
        middle.temperature_k(CellId(11)).unwrap(),
        middle.temperature_k(CellId(12)).unwrap(),
    );
    let va = ma / (1000.0 - 5.0 * (ta - 300.0));
    let vb = mb / (1000.0 - 5.0 * (tb - 300.0));
    let face = (3.0 * va / (4.0 * PI)).cbrt();
    let outer = (3.0 * (va + vb) / (4.0 * PI)).cbrt();
    let (xa, xb) = (centroid(0.0, face), centroid(face, outer));
    // Independent integrals of dr/(4*pi*k*r^2), plus the original areal film.
    let resistance = (1.0 / xa - 1.0 / face) / (4.0 * PI)
        + contact / (4.0 * PI * face * face)
        + (1.0 / face - 1.0 / xb) / (8.0 * PI);
    let heat = dt * (ta - tb) / resistance;
    assert!((heat - dt * original_g * (ta - tb)).abs() > 1e-5);
    let (after, geometry, audit) = middle
        .advance_observed(&mid_observer, 2.0 * dt, config(dt), || Ok(()))
        .unwrap();
    close(
        after.cells()[&CellId(11)].specific_h_j_kg,
        middle.cells()[&CellId(11)].specific_h_j_kg - heat / ma,
        1e-10,
    );
    close(
        after.cells()[&CellId(12)].specific_h_j_kg,
        middle.cells()[&CellId(12)].specific_h_j_kg + heat / mb,
        1e-10,
    );
    close(audit.external_heat_j, 0.0, 1e-13);
    close(audit.stored_energy_change_j, 0.0, 1e-10);
    close(audit.residual_j, 0.0, 1e-10);
    close(geometry.current_geometry().total_mass_kg, ma + mb, 1e-14);
    assert_eq!(audit.substeps, 1);
}

#[test]
fn outer_ambient_and_series_bath_use_current_surface_and_bulk_distance() {
    let card = material(1.0);
    let law = caloric(&card);
    let k = prepare(&card, law.curve());
    let mass = 1000.0 * volume(0.1);
    let body = sphere(&k, mass, 1);
    for film in [None, Some(0.01)] {
        let bound = body.bind(ids(), [None, film]).unwrap();
        let before = initial(&bound, &[0.0], film.map(|_| 360.0), |cell| {
            if film.is_none() {
                cell.boundary.htc_w_m2_k = 100.0;
                cell.boundary.ambient_k = 360.0;
            }
        });
        let observer = bound
            .spherical_liquid_density_observer(&before, || Ok(()))
            .unwrap();
        let dt = 0.5;
        let (middle, mid_observer, _) = before
            .advance_observed(&observer, dt, config(dt), || Ok(()))
            .unwrap();
        let temperature = middle.temperature_k(CellId(11)).unwrap();
        let density = 1000.0 - 5.0 * (temperature - 300.0);
        let radius = (3.0 * mass / (4.0 * PI * density)).cbrt();
        let area = 4.0 * PI * radius * radius;
        let conductance = film.map_or(100.0 * area, |r| {
            1.0 / ((1.0 / (0.75 * radius) - 1.0 / radius) / (4.0 * PI) + r / area)
        });
        let heat = dt * conductance * (360.0 - temperature);
        let (after, _, audit) = middle
            .advance_observed(&mid_observer, 2.0 * dt, config(dt), || Ok(()))
            .unwrap();
        close(audit.external_heat_j, heat, 1e-10);
        close(
            after.cells()[&CellId(11)].specific_h_j_kg,
            middle.cells()[&CellId(11)].specific_h_j_kg + heat / mass,
            1e-10,
        );
        assert_eq!(
            after.cells()[&CellId(11)].boundary,
            before.cells()[&CellId(11)].boundary
        );
        assert_eq!(audit.substeps, 1);
    }
}

#[test]
fn unsupported_constraints_inventory_and_exposed_interior_refuse_at_admission() {
    let card = material(1.0);
    let law = caloric(&card);
    let k = prepare(&card, law.curve());
    let layer = PhaseLayerSource {
        id: 1,
        thickness_m: 0.1,
        mass_kg: 1000.0 * volume(0.1),
        cell_count: 1,
        conductivity: &k,
    };
    for shape in [
        LayerStackGeometry::Spherical {
            inner_radius_m: 0.1,
        },
        LayerStackGeometry::Cylindrical {
            inner_radius_m: 0.1,
            length_m: 1.0,
        },
        LayerStackGeometry::Planar { area_m2: 1.0 },
    ] {
        let body =
            PreparedPhaseLayeredThermalBody::try_new(&[layer], &[], shape, || Ok(())).unwrap();
        let bound = body.bind(ids(), [None; 2]).unwrap();
        let thermal = initial(&bound, &[0.0], None, |_| {});
        assert!(
            bound
                .spherical_liquid_density_observer(&thermal, || Ok(()))
                .is_err()
        );
    }
    let wrong_mass = sphere(&k, layer.mass_kg * 1.01, 1);
    let wrong = wrong_mass.bind(ids(), [None; 2]).unwrap();
    let thermal = initial(&wrong, &[0.0], None, |_| {});
    assert!(
        wrong
            .spherical_liquid_density_observer(&thermal, || Ok(()))
            .is_err()
    );
    let body = sphere(&k, layer.mass_kg, 2);
    let bound = body.bind(ids(), [None; 2]).unwrap();
    let thermal = initial(&bound, &[0.0; 2], None, |cell| {
        if cell.id == CellId(11) {
            cell.boundary.htc_w_m2_k = 1.0;
        }
    });
    assert!(
        bound
            .spherical_liquid_density_observer(&thermal, || Ok(()))
            .is_err()
    );
    let body = sphere(&k, layer.mass_kg, 1);
    let bound = body.bind(ids(), [None, Some(0.01)]).unwrap();
    let thermal = initial(&bound, &[0.0], Some(350.0), |cell| {
        cell.boundary.htc_w_m2_k = 1.0
    });
    assert!(
        bound
            .spherical_liquid_density_observer(&thermal, || Ok(()))
            .is_err()
    );
}

#[test]
fn leaving_liquid_phase_or_cancelling_rolls_back_thermal_and_geometric_candidates() {
    let card = material(1.0);
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-100.0, 290.0, 0.0, 1000.0),
            (0.0, 300.0, 0.0, 1000.0),
            (100.0, 300.0, 1.0, 1000.0),
            (1100.0, 400.0, 1.0, 500.0),
        ]
        .into_iter()
        .map(|(h, t, f, rho)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: rho,
        })
        .collect(),
    )
    .unwrap();
    let k = prepare(&card, &curve);
    let body = sphere(&k, 1000.0 * volume(0.1), 1);
    let bound = body.bind(ids(), [None; 2]).unwrap();
    for h in [-10.0, 50.0] {
        let thermal = initial(&bound, &[h], None, |_| {});
        assert!(
            bound
                .spherical_liquid_density_observer(&thermal, || Ok(()))
                .is_err()
        );
    }
    let density = curve
        .values_at_specific_enthalpy(110.0)
        .unwrap()
        .bulk_density_kg_m3;
    let body = sphere(&k, density * volume(0.1), 1);
    let bound = body.bind(ids(), [None; 2]).unwrap();
    let before = initial(&bound, &[110.0], None, |cell| {
        cell.boundary.source_w = -20.0 * cell.mass_kg
    });
    let observer = bound
        .spherical_liquid_density_observer(&before, || Ok(()))
        .unwrap();
    let original_cells = before.cells().clone();
    let original_geometry = observer.current_geometry().clone();
    assert!(matches!(
        before.advance_observed(&observer, 1.0, config(0.1), || Ok(())),
        Err(DepositionError::InvalidInput(_))
    ));
    assert_eq!(before.time_s(), 0.0);
    assert_eq!(before.cells(), &original_cells);
    assert_eq!(observer.current_geometry(), &original_geometry);
    let mut calls = 0;
    assert!(matches!(
        before.advance_observed(&observer, 0.1, config(0.001), || {
            calls += 1;
            if calls == 100 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(before.cells(), &original_cells);
    assert_eq!(observer.current_geometry(), &original_geometry);
}
