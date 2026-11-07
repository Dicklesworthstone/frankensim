//! Live finite-cell heat, original H(T,z), actual slip and three normal accounts.

#[path = "depth/removal.rs"]
mod removal;

#[path = "depth/resistance.rs"]
mod resistance;

use super::*;
use fs_contact::normal_patch::FiniteGapPlasticLinearSolver;
use fs_contact::normal_patch::plastic_source::depth_wear::SourceDepthHardnessWear;
use fs_contact::normal_patch::plastic_source::wear::sliding::hardness::{
    HardnessSlidingControl, SourceDepthSlidingWear,
};
use fs_couple::friction_heat::plastic_wear::depth::DepthSlidingHeatContact;
use fs_material::hardness::depth::{EXPOSED_MATERIAL_DEPTH_AXIS as DEPTH, ResolvedDepthHardness};

const MAX_DEPTH: f64 = 1e-7;
const GRADIENT: [f64; 2] = [-5e12, 2e12];
const WEAR_K: [f64; 2] = [1e-6, 2e-6];

struct DepthFixture {
    cards: [MaterialCard; 2],
    caloric: [ResolvedHeatCapacityCurve; 2],
    source: SourceDepthSlidingWear,
}

fn depth_material(side: usize, maximum_temperature_k: f64) -> (MaterialCard, PropertyKey) {
    let (base, original_key) = material(&format!("depth-heat-body-{side}"), side);
    let mut claims = base.claims().clone();
    let key = PropertyKey::with_quantity("depth-hardness", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(original_key.hardness_test().unwrap().clone())
        .unwrap();
    let h = |t: f64, z: f64| {
        let surface = if side == 0 {
            1e6 * (1.0 - 0.01 * (t - T0))
        } else {
            2e6
        };
        surface + GRADIENT[side] * z
    };
    let length = Dims([1, 0, 0, 0, 0, 0]);
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Surface2D {
                axes: [DEPTH.into(), "T".into()],
                axis_dims: [length, Temperature::DIMS],
                samples: [
                    [0.0, 280.0],
                    [MAX_DEPTH, 280.0],
                    [MAX_DEPTH, maximum_temperature_k],
                    [0.0, maximum_temperature_k],
                ]
                .into_iter()
                .map(|[z, t]| ([z, t], h(t, z)))
                .collect(),
                triangles: vec![[0, 1, 2], [0, 2, 3]],
                dims: INDENTATION_HARDNESS_QUANTITY.dims(),
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 280.0, maximum_temperature_k)
                .with_quantity(DEPTH, QuantitySpec::dimensional(length), 0.0, MAX_DEPTH),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![key.hardness_test().unwrap().observation()],
            provenance: provenance(),
        })
        .unwrap();
    (
        MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap(),
        key,
    )
}

fn depth_fixture(maximum_temperature_k: f64, cx: &Cx<'_>) -> DepthFixture {
    depth_fixture_with_extra_claim(maximum_temperature_k, cx, None)
}

fn depth_fixture_with_extra_claim(
    maximum_temperature_k: f64,
    cx: &Cx<'_>,
    extra: Option<PropertyClaim>,
) -> DepthFixture {
    let [(a, ha), (b, hb)] = [0, 1].map(|side| depth_material(side, maximum_temperature_k));
    let cards = [a, b];
    let keys = [ha, hb];
    let caloric = cards.each_ref().map(|card| {
        ResolvedHeatCapacityCurve::try_new_exact(
            card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new(
                    "specific_heat_capacity",
                    SPECIFIC_HEAT_CAPACITY_DIMS,
                ),
                density: &PropertyKey::new("density", Density::DIMS),
                temperature_axis: "T",
                lower: &point(280.0),
                upper: &point(320.0),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            0.0,
        )
        .unwrap()
    });
    let base = source_common::bind(&cards[0], &cards[1], None, 0.0);
    let original = base.normal_state().elastic().interface().card();
    let mut claims = original.claims().clone();
    let wear_keys = ["depth_k_a", "depth_k_b"].map(|name| PropertyKey::new(name, Dims::NONE));
    for side in 0..2 {
        claims
            .insert_claim(claim(
                wear_keys[side].clone(),
                [280.0, 320.0].map(|t| WEAR_K[side] * (1.0 + 0.02 * (t - T0))),
            ))
            .unwrap();
    }
    for (key, values) in [
        (PropertyKey::new("mu_static", Dims::NONE), [0.9; 2]),
        (
            PropertyKey::new("mu_kinetic", Dims::NONE),
            [280.0, 320.0].map(|t| 0.4 + 0.01 * (t - T0)),
        ),
    ] {
        claims.insert_claim(claim(key, values)).unwrap();
    }
    if let Some(extra) = extra {
        claims.insert_claim(extra).unwrap();
    }
    let interface = InterfaceSystemCard::assemble(
        original.surface_a().clone(),
        original.surface_b().clone(),
        original.context().clone(),
        claims,
        original.models().to_vec(),
    )
    .unwrap();
    let query = point(T0);
    let select = MaterialPropertySelection::SingleClaimOnly;
    let solids = cards
        .each_ref()
        .map(|card| resolve_isotropic_solid_state_point(card, &query, select.clone()).unwrap());
    let dry = bind_dry_interface_system_card(&interface, InputAuthority::SyntheticFixture).unwrap();
    let elastic = bind_isotropic_elastic_interface(&dry, &solids[0], &solids[1]).unwrap();
    let adhesion = resolve_interface_state_point(
        &interface,
        &query,
        &[ScalarPropertyRequirement::try_new(
            ADHESION_ENERGY_PROPERTY,
            ADHESION_ENERGY_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        select.clone(),
    )
    .unwrap();
    let normal = bind_normal_interface_state(&elastic, &adhesion).unwrap();
    let profiles = [0, 1].map(|side| {
        ResolvedDepthHardness::new(
            &cards[side],
            &query,
            keys[side].clone(),
            select.clone(),
            MAX_DEPTH,
        )
        .unwrap()
    });
    let coefficients = resolve_interface_state_point(
        &interface,
        &query,
        &wear_keys.each_ref().map(|key| {
            ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative).unwrap()
        }),
        select.clone(),
    )
    .unwrap();
    let mut gap = vec![1e-3; 25];
    gap[12] = 0.0;
    let wear = SourceDepthHardnessWear::new(
        normal,
        coefficients,
        wear_keys,
        profiles,
        FiniteGapGrid {
            cells_x: 5,
            cells_y: 5,
            cell_width_m: SIDE,
            cell_depth_m: SIDE,
        },
        gap,
        PlasticContactControls {
            maximum_active_set_iterations: 128,
            complementarity_tolerance_m: 1e-12,
            boundary_clearance_cells: 1,
        },
        PlasticContactEnvelope {
            maximum_approach_m: 1e-4,
            maximum_plastic_depth_m: 1e-4,
            half_space_depth_m: [LENGTH; 2],
            maximum_patch_to_depth: [0.1; 2],
        },
        [MAX_DEPTH; 2],
        FiniteGapPlasticLinearSolver::Dense,
        cx,
    )
    .unwrap();
    let requirement = |name| {
        ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative)
            .unwrap()
    };
    let friction = SourceThermalCoulomb::resolve(
        dry,
        query,
        requirement("mu_static"),
        requirement("mu_kinetic"),
        select,
        ThermalFrictionAxes {
            temperature: "T".into(),
            nominal_pressure: None,
        },
    )
    .unwrap();
    DepthFixture {
        cards,
        caloric,
        source: SourceDepthSlidingWear::new(wear, friction).unwrap(),
    }
}

fn depth_loaded(data: &DepthFixture, cx: &Cx<'_>) -> SourceDepthSlidingWear {
    data.source
        .advance_normal(HardnessSlidingControl::Approach(5e-5), cx)
        .unwrap()
        .candidate
}

fn bind_depth<'a>(
    data: &'a DepthFixture,
    thermal: &DepositionSimulation<'_>,
    route: &FrictionHeatRoute,
    source: SourceDepthSlidingWear,
    limits: PlasticSlidingHeatLimits,
    cx: &Cx<'_>,
) -> DepthSlidingHeatContact<'a> {
    DepthSlidingHeatContact::new(
        thermal,
        source,
        [&data.cards[0], &data.cards[1]],
        [&data.caloric[0], &data.caloric[1]],
        route,
        limits,
        cx,
    )
    .unwrap()
}

#[test]
fn g1_live_depth_hardness_heat_matches_independent_quadratic_and_three_work_accounts() {
    with_cx(|cx, _| {
        let data = depth_fixture(320.0, cx);
        // A preconditioned, cooled specimen starts with actual removed depth.
        let preworn = depth_loaded(&data, cx)
            .advance(
                &motion(0.005),
                0.1,
                controls().friction,
                HardnessSlidingControl::Approach(5e-5),
                cx,
            )
            .unwrap()
            .candidate;
        let saved_source = preworn.clone();
        let thermal = network(&data.caloric, T0);
        let saved_cells = thermal.cells().clone();
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let contact = bind_depth(&data, &thermal, &route, preworn, limits(), cx);
        let distance = 0.01;
        let step = contact
            .step_candidate(
                &thermal,
                0.1,
                &motion(distance),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                &[0.0; 2],
                T0,
                transport,
                cx,
            )
            .unwrap();
        let report = step.report();
        let capacity = RHO * AREA * LENGTH * CP;
        let half_g = CONDUCTIVITY * AREA / (0.5 * LENGTH);
        let b = 0.1 * half_g / (capacity + 0.1 * half_g);
        let z0 = saved_source.wear_state().recession_m()[0][12];
        let cold_h = 1e6 + GRADIENT[0] * z0;
        let measure = AREA * distance;
        // 2 C b theta = A s (H_cold - 1e4 b theta) (0.4 + 0.01 theta).
        let qa = measure * 1e4 * b * 0.01;
        let qb = 2.0 * capacity * b - measure * (0.01 * cold_h - 0.4 * 1e4 * b);
        let qc = measure * cold_h * 0.4;
        let theta = 2.0 * qc / (qb + (qb * qb + 4.0 * qa * qc).sqrt());
        let body_k = T0 + b * theta;
        let pressure = cold_h - 1e4 * b * theta;
        let friction_work = measure * pressure * (0.4 + 0.01 * theta);
        close(report.heat.interface_temperature_k, T0 + theta, 0.0, 2e-6);
        close(report.frictional_work_j, friction_work, 2e-6, 1e-13);
        close(
            report
                .mechanics
                .temperature
                .candidate
                .wear_state()
                .contact_state()
                .pressure_pa()[12],
            pressure,
            2e-6,
            1e-3,
        );
        let mut stored = 0.0;
        for (side, id) in CELLS.into_iter().enumerate() {
            close(step.thermal().temperature_k(id).unwrap(), body_k, 0.0, 2e-6);
            close(report.source[side], body_k, 0.0, 2e-6);
            let old = &thermal.cells()[&id];
            let new = &step.thermal().cells()[&id];
            stored += old.mass_kg * (new.specific_h_j_kg - old.specific_h_j_kg);
            let z = saved_source.wear_state().recession_m()[side][12];
            let h = if side == 0 {
                pressure
            } else {
                2e6 + GRADIENT[side] * z
            };
            let k = WEAR_K[side] * (1.0 + 0.02 * theta);
            let q = k * pressure * distance;
            let increment = 2.0 * q / (h + (h * h + 2.0 * GRADIENT[side] * q).sqrt());
            close(
                report.mechanics.sliding.recession_increment_m[side][12],
                increment,
                3e-6,
                1e-15,
            );
            close(
                report.mechanics.sliding.wear[side].increment_m3,
                AREA * increment,
                3e-6,
                1e-23,
            );
            close(
                report.removed_volume_m3[side],
                saved_source.wear_state().volumes_m3()[side] + AREA * increment,
                3e-6,
                1e-23,
            );
        }
        close(stored, friction_work, 2e-6, 3e-11);
        close(
            report.heat.body_heat_j.iter().sum(),
            friction_work,
            2e-6,
            3e-11,
        );
        assert_eq!(report.mechanics.wear_normals().len(), 2);
        let normals = std::iter::once(report.mechanics.temperature_normal())
            .chain(report.mechanics.wear_normals().iter().map(|r| &r.normal));
        let (mut normal_work, mut plastic_work, mut parameter_work) = (0.0, 0.0, 0.0);
        for normal in normals {
            normal_work += normal.normal_work_j;
            plastic_work += normal.response.plastic_work_j;
            parameter_work += normal.parameter_work_j;
        }
        close(report.normal_work_j, normal_work, 1e-12, 1e-20);
        close(report.normal_plastic_work_j, plastic_work, 1e-12, 1e-20);
        close(
            report.temperature_parameter_work_j + report.wear_parameter_work_j,
            parameter_work,
            1e-12,
            1e-20,
        );
        assert!(
            report.mechanics.wear_normals()[1]
                .normal
                .response
                .plastic_work_j
                > 0.0
        );
        assert_eq!(report.normal_plastic_heat_j, [0.0; 2]);
        close(report.nonthermal_plastic_work_j, plastic_work, 1e-12, 1e-20);
        assert!(report.heat.total_energy_residual_j.abs() <= controls().heat.total_j);
        assert_eq!(contact.source(), &saved_source);
        assert_eq!(thermal.cells(), &saved_cells);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(step.thermal().time_s(), 0.1);
    });
}

#[test]
fn g1_depth_heat_force_spring_and_reversal_preserve_actual_history_and_thermal_feedback() {
    with_cx(|cx, _| {
        let data = depth_fixture(320.0, cx);
        let force = 0.004;
        let source = data
            .source
            .advance_normal(
                HardnessSlidingControl::Load {
                    force_n: force,
                    controls: normal_controls(),
                },
                cx,
            )
            .unwrap()
            .candidate;
        let approach = source.wear_state().contact_state().approach_m();
        for boundary in [
            PlasticSlidingHeatBoundary::MaintainedLoad {
                force_n: force,
                controls: normal_controls(),
            },
            PlasticSlidingHeatBoundary::FixedDrive {
                drive: FiniteGapPlasticDrive {
                    stiffness_n_per_m: 1000.0,
                    position_m: approach + force / 1000.0,
                },
                controls: normal_controls(),
            },
        ] {
            let thermal = network(&data.caloric, T0);
            let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
            let contact = bind_depth(&data, &thermal, &route, source.clone(), limits(), cx);
            let distance = 0.01;
            let forward = contact
                .step_candidate(
                    &thermal,
                    0.1,
                    &motion(distance),
                    boundary,
                    controls(),
                    &[0.0; 2],
                    T0,
                    transport,
                    cx,
                )
                .unwrap();
            let capacity = RHO * AREA * LENGTH * CP;
            let g = CONDUCTIVITY * AREA / (0.5 * LENGTH);
            let b = 0.1 * g / (capacity + 0.1 * g);
            let theta = force * distance * 0.4 / (2.0 * capacity * b - force * distance * 0.01);
            close(
                forward.report().heat.interface_temperature_k,
                T0 + theta,
                0.0,
                2e-6,
            );
            close(
                forward.report().frictional_work_j,
                force * distance * (0.4 + 0.01 * theta),
                2e-6,
                1e-13,
            );
            let reverse = forward
                .contact()
                .step_candidate(
                    forward.thermal(),
                    0.2,
                    &motion(-distance),
                    boundary,
                    controls(),
                    &[0.0; 2],
                    forward.report().heat.interface_temperature_k,
                    transport,
                    cx,
                )
                .unwrap();
            assert_eq!(
                reverse.contact().source().tangential_states()[12].slip_m(),
                [0.0; 2]
            );
            assert!(reverse.report().frictional_work_j > 0.0);
            assert!(
                reverse.thermal().temperature_k(CELLS[0]).unwrap()
                    > forward.thermal().temperature_k(CELLS[0]).unwrap()
            );
            for side in 0..2 {
                assert!(
                    reverse.report().removed_volume_m3[side]
                        > forward.report().removed_volume_m3[side]
                );
            }
            if matches!(boundary, PlasticSlidingHeatBoundary::MaintainedLoad { .. }) {
                close(
                    reverse
                        .contact()
                        .source()
                        .wear_state()
                        .contact_state()
                        .pressure_pa()[12]
                        * AREA,
                    force,
                    1e-10,
                    1e-13,
                );
            }
            assert!(
                reverse.report().contact_work_residual_j <= controls().contact_work_tolerance_j
            );
            assert_eq!(reverse.thermal().time_s(), 0.2);
            assert_eq!(thermal.time_s(), 0.0);
        }
    });
}

#[test]
fn g4_depth_heat_source_inventory_budget_and_cancellation_refuse_without_partial_state() {
    with_cx(|cx, gate| {
        for maximum_temperature in [T0 + 0.01, 320.0] {
            let data = depth_fixture(maximum_temperature, cx);
            let thermal = network(&data.caloric, T0);
            let saved_cells = thermal.cells().clone();
            let source = depth_loaded(&data, cx);
            let saved_source = source.clone();
            let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
            let mut bound = limits();
            if maximum_temperature == 320.0 {
                bound.maximum_removed_volume_fraction = [1e-8; 2];
            }
            let contact = bind_depth(&data, &thermal, &route, source, bound, cx);
            assert!(
                contact
                    .step_candidate(
                        &thermal,
                        0.1,
                        &motion(0.01),
                        PlasticSlidingHeatBoundary::HeldApproach,
                        controls(),
                        &[0.0; 2],
                        T0,
                        transport,
                        cx
                    )
                    .is_err()
            );
            assert_eq!(contact.source(), &saved_source);
            assert_eq!(thermal.cells(), &saved_cells);
            assert_eq!(thermal.time_s(), 0.0);
        }
        let data = depth_fixture(320.0, cx);
        let thermal = network(&data.caloric, T0);
        let source = depth_loaded(&data, cx);
        let saved_source = source.clone();
        let saved_cells = thermal.cells().clone();
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let contact = bind_depth(&data, &thermal, &route, source, limits(), cx);
        let mut limited = controls();
        limited.maximum_trials = 1;
        assert!(
            contact
                .step_candidate(
                    &thermal,
                    0.1,
                    &motion(0.01),
                    PlasticSlidingHeatBoundary::HeldApproach,
                    limited,
                    &[0.0; 2],
                    T0,
                    transport,
                    cx
                )
                .is_err()
        );
        gate.request();
        assert!(
            contact
                .step_candidate(
                    &thermal,
                    0.1,
                    &motion(0.01),
                    PlasticSlidingHeatBoundary::HeldApproach,
                    controls(),
                    &[0.0; 2],
                    T0,
                    transport,
                    cx
                )
                .is_err()
        );
        assert_eq!(contact.source(), &saved_source);
        assert_eq!(thermal.cells(), &saved_cells);
    });
}
