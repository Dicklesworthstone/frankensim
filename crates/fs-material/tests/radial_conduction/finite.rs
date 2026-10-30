use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_material::thermal::conduction::network::{
    GeometricLayerLinkBinding, LayerConductionObserver, PreparedLayerConduction,
};
use fs_material::thermal::conduction::source::{LayerSourceError, ResolvedConductionLayer};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Density, Temperature};
use fs_thermal::enthalpy::{backward_euler::BackwardEulerStep, *};

fn card(phase: SolidLiquidPhase, rho: f64, cp: f64, k: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, rho),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            cp,
        ),
        (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, k),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value), (1000.0, value)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 200.0, 1000.0),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic radial material fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: format!("{phase:?}"),
            process: "manufactured test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn chart(card: &MaterialCard, phase: SolidLiquidPhase) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &QueryPoint::new().with("T", 200.0).unwrap(),
            upper: &QueryPoint::new().with("T", 1000.0).unwrap(),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        phase,
        0.0,
        1e-8,
    )
    .unwrap()
}
fn layer<'a>(
    card: &'a MaterialCard,
    chart: &'a ResolvedHeatCapacityCurve,
    context: &'a QueryPoint,
    phase: SolidLiquidPhase,
) -> ResolvedConductionLayer<'a> {
    ResolvedConductionLayer::new(
        chart.curve(),
        ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        &PropertyKey::new(THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        phase,
        [200.0, 1000.0],
        || Ok(()),
    )
    .unwrap()
}
fn geometry() -> LayerGeometry {
    LayerGeometry::spherical([0.015, 0.025], 0.02).unwrap()
}
fn thermal<'a>(charts: [&'a ResolvedHeatCapacityCurve; 2]) -> DepositionSimulation<'a> {
    let radii = [0.01_f64, 0.02, 0.03];
    let pi = core::f64::consts::PI;
    let mut events = Vec::new();
    for (i, t) in [600.0, 300.0].into_iter().enumerate() {
        let state = charts[i].state_at_temperature(t).unwrap();
        let volume = 4.0 * pi / 3.0 * (radii[i + 1].powi(3) - radii[i].powi(3));
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                body_id: 1,
                segment_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: state.bulk_density_kg_m3() * volume,
                specific_h_j_kg: state.specific_enthalpy_j_kg(),
                geometry: CellGeometry {
                    start_m: [radii[i], 0.0, 0.0],
                    end_m: [radii[i + 1], 0.0, 0.0],
                    volume_m3: volume,
                    surface_area_m2: 4.0 * pi * (radii[i].powi(2) + radii[i + 1].powi(2)),
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: None,
                },
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: CellId(0),
            b: CellId(1),
            area_m2: geometry().contact_area_m2(),
            conductance_w_k: 0.0,
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), charts[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), charts[1].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(0.0, step(1.0), || Ok(()))
    .unwrap()
    .0
}
fn binding<'a, 'm>(
    layers: [&'a ResolvedConductionLayer<'m>; 2],
) -> GeometricLayerLinkBinding<'a, 'm> {
    GeometricLayerLinkBinding {
        link: 7,
        layers,
        geometry: geometry(),
        contact_area_resistance_m2_k_w: 0.0,
    }
}
fn step(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-7,
    }
}
fn capacities(state: &DepositionSimulation<'_>, cp: [f64; 2]) -> [f64; 2] {
    core::array::from_fn(|i| state.cells()[&CellId(i as u64)].mass_kg * cp[i])
}
fn conductance(k: [f64; 2]) -> f64 {
    4.0 * core::f64::consts::PI
        / ((1.0 / 0.015 - 1.0 / 0.02) / k[0] + (1.0 / 0.02 - 1.0 / 0.025) / k[1])
}

#[test]
fn g1_finite_solid_and_liquid_exchange_once_with_native_energy_entropy_and_refinement() {
    let a = card(SolidLiquidPhase::Solid, 7800.0, 400.0, 20.0);
    let b = card(SolidLiquidPhase::Liquid, 1000.0, 4200.0, 0.6);
    let ca = chart(&a, SolidLiquidPhase::Solid);
    let cb = chart(&b, SolidLiquidPhase::Liquid);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let la = layer(&a, &ca, &context, SolidLiquidPhase::Solid);
    let lb = layer(&b, &cb, &context, SolidLiquidPhase::Liquid);
    let old = thermal([&ca, &cb]);
    let model =
        PreparedLayerConduction::new_geometric(&old, vec![binding([&la, &lb])], 1, || Ok(()))
            .unwrap();
    let observer = LayerConductionObserver::new(&model, &old, config(), || Ok(())).unwrap();
    let c = capacities(&old, [400.0, 4200.0]);
    let g = conductance([20.0, 0.6]);
    let mean = (c[0] * 600.0 + c[1] * 300.0) / (c[0] + c[1]);
    let target = 25.0;
    let difference = 300.0 * (-g * (1.0 / c[0] + 1.0 / c[1]) * target).exp();
    let exact = [
        mean + c[1] / (c[0] + c[1]) * difference,
        mean - c[0] / (c[0] + c[1]) * difference,
    ];
    let mut errors = Vec::new();
    for dt in [0.25, 0.125] {
        let (next, observed, energy, entropy) = old
            .advance_observed_with_entropy(&observer, target, step(dt), || Ok(()))
            .unwrap();
        let t = [
            next.temperature_k(CellId(0)).unwrap(),
            next.temperature_k(CellId(1)).unwrap(),
        ];
        errors.push((t[0] - exact[0]).abs().max((t[1] - exact[1]).abs()));
        assert!(t[0] < 600.0 && t[1] > 300.0);
        close(next.energy_j(), old.energy_j(), 1e-7);
        assert_eq!(energy.external_heat_j, 0.0);
        assert!(energy.residual_j.abs() < 1e-7);
        assert!(entropy.conduction_generation_j_k > 0.0);
        assert!(!entropy.has_unclassified_sources);
        assert_eq!(entropy.reservoir_entropy_input_j_k, 0.0);
        close(
            entropy.thermal_storage_change_j_k,
            c[0] * (t[0] / 600.0).ln() + c[1] * (t[1] / 300.0).ln(),
            1e-10,
        );
        assert_eq!(observed.last_samples()[0].geometry, geometry());
        assert_eq!(next.links()[&7], old.links()[&7]);
        assert_eq!(
            next.cells()[&CellId(0)].mass_kg,
            old.cells()[&CellId(0)].mass_kg
        );
    }
    assert!(errors[1] < 0.51 * errors[0], "time refinement: {errors:?}");
    assert_eq!(old.time_s(), 0.0);
    // The original implicit residual accepts the same source/geometry model.
    let dt = 2.0;
    let q = dt * g * 300.0 / (1.0 + dt * g * (1.0 / c[0] + 1.0 / c[1]));
    let implicit = BackwardEulerStep::new(&old, dt, || Ok(())).unwrap();
    let trial = implicit
        .evaluate(
            &[-q, q],
            |s| model.resolve(s, config(), || Ok(())),
            || Ok(()),
        )
        .unwrap();
    let (next, used, audit) = trial.accept(1e-7, 1e-7, || Ok(())).unwrap();
    close(used[0].applied_power_a_to_b_w * dt, q, 1e-7);
    close(
        next.temperature_k(CellId(1)).unwrap(),
        300.0 + q / c[1],
        1e-10,
    );
    assert_eq!(audit.external_heat_j, 0.0);
}

#[test]
fn g3_material_only_capacity_conductivity_and_density_changes_reach_both_finite_states() {
    let mut answers = Vec::new();
    for (rho, cp, k) in [
        (1000.0, 4200.0, 0.6),
        (1000.0, 8400.0, 0.6),
        (1000.0, 4200.0, 1.2),
        (2000.0, 4200.0, 0.6),
    ] {
        let a = card(SolidLiquidPhase::Solid, 7800.0, 400.0, 20.0);
        let b = card(SolidLiquidPhase::Liquid, rho, cp, k);
        let ca = chart(&a, SolidLiquidPhase::Solid);
        let cb = chart(&b, SolidLiquidPhase::Liquid);
        let context = QueryPoint::new().with("T", 300.0).unwrap();
        let la = layer(&a, &ca, &context, SolidLiquidPhase::Solid);
        let lb = layer(&b, &cb, &context, SolidLiquidPhase::Liquid);
        let old = thermal([&ca, &cb]);
        let model =
            PreparedLayerConduction::new_geometric(&old, vec![binding([&la, &lb])], 1, || Ok(()))
                .unwrap();
        let observer = LayerConductionObserver::new(&model, &old, config(), || Ok(())).unwrap();
        let (next, _, audit) = old
            .advance_observed(&observer, 25.0, step(0.125), || Ok(()))
            .unwrap();
        assert_eq!(audit.external_heat_j, 0.0);
        assert!(audit.residual_j.abs() < 1e-7);
        answers.push([
            next.temperature_k(CellId(0)).unwrap(),
            next.temperature_k(CellId(1)).unwrap(),
        ]);
    }
    assert!(answers[1][1] < answers[0][1]);
    assert!(answers[2][0] < answers[0][0] && answers[2][1] > answers[0][1]);
    close(answers[3][0], answers[1][0], 1e-10);
    close(answers[3][1], answers[1][1], 1e-10);
}

#[test]
fn g4_late_cancel_foreign_source_and_changed_area_publish_no_thermal_history() {
    let a = card(SolidLiquidPhase::Solid, 7800.0, 400.0, 20.0);
    let b = card(SolidLiquidPhase::Liquid, 1000.0, 4200.0, 0.6);
    let ca = chart(&a, SolidLiquidPhase::Solid);
    let cb = chart(&b, SolidLiquidPhase::Liquid);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let la = layer(&a, &ca, &context, SolidLiquidPhase::Solid);
    let lb = layer(&b, &cb, &context, SolidLiquidPhase::Liquid);
    let old = thermal([&ca, &cb]);
    let model =
        PreparedLayerConduction::new_geometric(&old, vec![binding([&la, &lb])], 1, || Ok(()))
            .unwrap();
    let observer = LayerConductionObserver::new(&model, &old, config(), || Ok(())).unwrap();
    let original = observer.last_samples().to_vec();
    let mut polls = 0;
    let expected = old
        .advance_observed(&observer, 0.25, step(0.125), || {
            polls += 1;
            Ok(())
        })
        .unwrap()
        .0;
    let mut calls = 0;
    assert!(matches!(
        old.advance_observed(&observer, 0.25, step(0.125), || {
            calls += 1;
            if calls == polls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(LayerSourceError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(observer.last_samples(), original);
    let replay = old
        .advance_observed(&observer, 0.25, step(0.125), || Ok(()))
        .unwrap()
        .0;
    assert_eq!(replay.cells(), expected.cells());
    assert_eq!(old.time_s(), 0.0);
    let mut bad = binding([&la, &lb]);
    bad.geometry = LayerGeometry::spherical([0.015, 0.025], 0.021).unwrap();
    assert!(PreparedLayerConduction::new_geometric(&old, vec![bad], 1, || Ok(())).is_err());
    let foreign = ca.clone();
    let foreign_layer = layer(&a, &foreign, &context, SolidLiquidPhase::Solid);
    assert!(
        PreparedLayerConduction::new_geometric(
            &old,
            vec![binding([&foreign_layer, &lb])],
            1,
            || Ok(())
        )
        .is_err()
    );
    let old_patch = old.contact_patch(7).unwrap();
    let mut changed = old_patch;
    changed.covered_area_m2 = [0.002; 2];
    let changed = old
        .replace_contact_patches(
            &[ContactPatchUpdate {
                id: 7,
                expected: Some(old_patch),
                replacement: Some(changed),
            }],
            || Ok(()),
        )
        .unwrap();
    assert!(
        changed
            .advance_observed(&observer, 0.25, step(0.125), || Ok(()))
            .is_err()
    );
}
