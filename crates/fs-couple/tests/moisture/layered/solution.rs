//! Authored finite-solution/coating laws checked against independent balances.
use super::*;
use crate::solution as aqueous;
use fs_couple::moisture::{
    radial::{
        MoistureRadial, MoistureRadialGeometry, RadialMoistureSymmetry,
        layered::{MoistureLayeredRadial, MoistureRadialLayer},
    },
    solution::{
        FiniteMoistureSolution, LayeredMoistureParticipant, LayeredRadialMoistureParticipant,
        SolutionMoistureError, SolutionMoistureStep,
    },
};
use fs_material::osmosis::{
    AqueousSolution, OsmosisError,
    nonideal::{NonidealSolutionError, VirialSolutionModel},
    solubility::kinetics::KineticSolubilitySolution,
};
use fs_qty::{Amount, AmountConcentration, Temperature, Volume, VolumetricFlowRate};
use fs_thermochem::{
    SpeciesId,
    solution::solubility::{
        IdealSolubilityV1,
        kinetics::{PrecipitationPolicyV1, SolubilityKineticStateV1, SolubilityKineticsV1},
    },
};
use std::sync::Arc;

#[path = "solution/sdirk.rs"]
mod sdirk;

const RHO: [f64; 2] = [300.0, 900.0];
const SLOPE: [f64; 2] = [0.2, 0.4];
const DIFFUSION: [f64; 2] = [0.001, 0.0005];
const JOIN_H: f64 = 0.03;
const OUTER_H: f64 = 0.04;
const WATER_DENSITY: f64 = 0.018 / aqueous::WATER_VOLUME;

fn slab_body(initial: [f64; 2], cx: &Cx<'_>) -> MoistureLayeredSlab {
    MoistureLayeredSlab::new(
        (0..2)
            .map(|i| {
                layer(
                    material(RHO[i], vec![(0.0, 0.0), (1.0, SLOPE[i])], [DIFFUSION[i]; 2]),
                    [0.03, 0.02][i],
                    vec![SLOPE[i] * initial[i]],
                    cx,
                )
            })
            .collect(),
        2,
        cx,
    )
    .unwrap()
    .with_transfer(
        [
            MoistureActivityTransfer::IdealContact,
            MoistureActivityTransfer::ConductanceKgM2S(OUTER_H),
        ],
        vec![MoistureActivityTransfer::ConductanceKgM2S(JOIN_H)],
        cx,
    )
    .unwrap()
}

fn radial_body(
    symmetry: RadialMoistureSymmetry,
    initial: [f64; 2],
    cx: &Cx<'_>,
) -> MoistureLayeredRadial {
    MoistureLayeredRadial::new(
        (0..2)
            .map(|i| {
                let card = material(RHO[i], vec![(0.0, 0.0), (1.0, SLOPE[i])], [DIFFUSION[i]; 2]);
                let sorption = isotherm(&card);
                let field = MoistureRadial::new(
                    card,
                    point(0.1),
                    "moisture_diffusivity",
                    MaterialPropertySelection::SingleClaimOnly,
                    MoistureRadialGeometry {
                        symmetry,
                        faces_m: if i == 0 {
                            vec![0.0, 0.03]
                        } else {
                            vec![0.03, 0.05]
                        },
                        maximum_moisture: 0.6,
                        maximum_cells: 2,
                    },
                    vec![SLOPE[i] * initial[i]],
                    cx,
                )
                .unwrap();
                MoistureRadialLayer::new(field, sorption, cx).unwrap()
            })
            .collect(),
        2,
        cx,
    )
    .unwrap()
    .with_transfer(
        [
            MoistureActivityTransfer::IdealContact,
            MoistureActivityTransfer::ConductanceKgM2S(OUTER_H),
        ],
        vec![MoistureActivityTransfer::ConductanceKgM2S(JOIN_H)],
        cx,
    )
    .unwrap()
}

struct Physics {
    capacity: [f64; 2],
    half: [f64; 2],
    join: f64,
    outer: f64,
    outer_bulk: f64,
}

fn coefficients(symmetry: Option<RadialMoistureSymmetry>) -> Physics {
    let volume = |a: f64, b: f64| match symmetry {
        None => 0.001 * (b - a),
        Some(RadialMoistureSymmetry::Sphere) => {
            4.0 * core::f64::consts::PI * (b.powi(3) - a.powi(3)) / 3.0
        }
        Some(RadialMoistureSymmetry::Cylinder { length_m }) => {
            core::f64::consts::PI * length_m * (b * b - a * a)
        }
    };
    let area = |r: f64| match symmetry {
        None => 0.001,
        Some(RadialMoistureSymmetry::Sphere) => 4.0 * core::f64::consts::PI * r * r,
        Some(RadialMoistureSymmetry::Cylinder { length_m }) => {
            2.0 * core::f64::consts::PI * length_m * r
        }
    };
    let resistance = |a: f64, b: f64| match symmetry {
        None => (b - a) / 0.001,
        Some(RadialMoistureSymmetry::Sphere) => (1.0 / a - 1.0 / b) / (4.0 * core::f64::consts::PI),
        Some(RadialMoistureSymmetry::Cylinder { length_m }) => {
            (b / a).ln() / (2.0 * core::f64::consts::PI * length_m)
        }
    };
    let capacity = [
        RHO[0] * SLOPE[0] * volume(0.0, 0.03),
        RHO[1] * SLOPE[1] * volume(0.03, 0.05),
    ];
    let half = [
        RHO[0] * DIFFUSION[0] * SLOPE[0] / resistance(0.015, 0.03),
        RHO[1] * DIFFUSION[1] * SLOPE[1] / resistance(0.03, 0.04),
    ];
    let join = 1.0 / (1.0 / half[0] + 1.0 / (JOIN_H * area(0.03)) + 1.0 / half[1]);
    let outer_bulk = RHO[1] * DIFFUSION[1] * SLOPE[1] / resistance(0.04, 0.05);
    let outer = 1.0 / (1.0 / outer_bulk + 1.0 / (OUTER_H * area(0.05)));
    Physics {
        capacity,
        half,
        join,
        outer,
        outer_bulk,
    }
}

// Independent two-cell BE equations with actual water capacities and the
// three material/interface resistances. The target bath activity is prescribed
// analytically, then the ORIGINAL bath volume is fixed by closed water balance.
fn endpoint(p: &Physics, initial: [f64; 2], dt: f64, activity: f64) -> ([f64; 2], f64) {
    let a = p.capacity[0] + dt * p.join;
    let b = dt * p.join;
    let d = p.capacity[1] + dt * (p.join + p.outer);
    let rhs = [
        p.capacity[0] * initial[0],
        p.capacity[1] * initial[1] + dt * p.outer * activity,
    ];
    let determinant = a * d - b * b;
    let next = [
        (d * rhs[0] + b * rhs[1]) / determinant,
        (b * rhs[0] + a * rhs[1]) / determinant,
    ];
    let water = (0..2).map(|i| p.capacity[i] * (next[i] - initial[i])).sum();
    (next, water)
}

fn check_endpoint(
    step: &SolutionMoistureStep<MoistureLayeredReport>,
    p: &Physics,
    expected: [f64; 2],
    activity: f64,
    water: f64,
    dt: f64,
) {
    near(step.water_activity[1], activity, 2e-11);
    near(step.water_to_solid_kg, water, 2e-13);
    near(step.water_balance_residual_kg, 0.0, 1e-13);
    near(step.transport.boundary_water_kg[1], water, 2e-13);
    assert_eq!(step.transport.boundary_water_kg[0], 0.0);
    let flux = p.join * (expected[0] - expected[1]);
    let interface = &step.transport.interfaces[0];
    near(interface.water_left_to_right_kg, dt * flux, 2e-13);
    near(
        interface.relative_activity,
        expected[0] - flux / p.half[0],
        2e-11,
    );
    near(
        interface.right_relative_activity,
        expected[1] + flux / p.half[1],
        2e-11,
    );
    for i in 0..2 {
        let face_activity = [
            interface.relative_activity,
            interface.right_relative_activity,
        ][i];
        near(interface.moisture[i], SLOPE[i] * face_activity, 2e-12);
    }
    near(
        step.transport.surface_activities[1].unwrap(),
        expected[1] + water / (dt * p.outer_bulk),
        2e-11,
    );
    assert_eq!(step.transport.sorption_sources.len(), 2);
    assert_eq!(step.transport.diffusivity_sources.len(), 2);
}

#[test]
fn layered_solution_slab_keeps_finite_films_partition_and_analytic_water_feedback() {
    with_cx(|cx| {
        let final_volume = 1e-5;
        let solute = 0.02;
        let activity = (-aqueous::WATER_VOLUME * solute / final_volume).exp();
        let p = coefficients(None);
        for initial in [[0.1, 0.2], [0.98, 0.99]] {
            let body = slab_body(initial, cx);
            let original_activity = body.activity().to_vec();
            let (expected, water) = endpoint(&p, initial, 2.0, activity);
            let initial_volume = final_volume + water / WATER_DENSITY;
            let bath = aqueous::reservoir(initial_volume, vec![solute]);
            let saved = bath.solution().clone();
            let result = bath
                .advance_layered(
                    &body,
                    2.0,
                    [false, true],
                    config(),
                    nonlinear(),
                    aqueous::controls(initial_volume),
                    cx,
                )
                .unwrap();
            for i in 0..2 {
                near(result.moisture().activity()[i], expected[i], 2e-11);
                near(
                    result.moisture().layers()[i].slab().moisture()[0],
                    SLOPE[i] * expected[i],
                    2e-12,
                );
                assert_eq!(
                    result.moisture().layers()[i].sorption().identity(),
                    body.layers()[i].sorption().identity()
                );
            }
            check_endpoint(result.step(), &p, expected, activity, water, 2.0);
            near(
                result.reservoir().solution().solvent_volume_m3(),
                final_volume,
                2e-16,
            );
            near(
                result.reservoir().solvent_mass_kg() + result.moisture().water_mass_kg(),
                bath.solvent_mass_kg() + body.water_mass_kg(),
                2e-13,
            );
            assert_eq!(
                result.reservoir().solution().solute_moles(),
                saved.solute_moles()
            );
            assert_eq!(result.moisture().time_s(), 2.0);
            assert_eq!(bath.solution(), &saved);
            assert_eq!(body.activity(), original_activity);
            assert_eq!(body.time_s(), 0.0);
        }
    });
}

#[test]
fn coated_solution_uses_each_curved_area_and_actual_unequal_dry_inventory() {
    with_cx(|cx| {
        let final_volume = 1e-5;
        let solute = 0.02;
        let activity = (-aqueous::WATER_VOLUME * solute / final_volume).exp();
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
        ] {
            let p = coefficients(Some(symmetry));
            for initial in [[0.1, 0.2], [0.98, 0.99]] {
                let body = radial_body(symmetry, initial, cx);
                let original_activity = body.activity().to_vec();
                let (expected, water) = endpoint(&p, initial, 2.0, activity);
                let initial_volume = final_volume + water / WATER_DENSITY;
                let bath = aqueous::reservoir(initial_volume, vec![solute]);
                let saved = bath.solution().clone();
                let result = bath
                    .advance_layered_radial(
                        &body,
                        2.0,
                        [false, true],
                        config(),
                        nonlinear(),
                        aqueous::controls(initial_volume),
                        cx,
                    )
                    .unwrap();
                for i in 0..2 {
                    near(result.moisture().activity()[i], expected[i], 2e-11);
                    near(
                        result.moisture().layers()[i].radial().moisture()[0],
                        SLOPE[i] * expected[i],
                        2e-12,
                    );
                    near(
                        result.moisture().layers()[i].radial().dry_mass_kg(),
                        p.capacity[i] / SLOPE[i],
                        2e-14,
                    );
                    assert_eq!(
                        result.moisture().layers()[i].sorption().identity(),
                        body.layers()[i].sorption().identity()
                    );
                }
                check_endpoint(result.step(), &p, expected, activity, water, 2.0);
                near(
                    result.reservoir().solution().solvent_volume_m3(),
                    final_volume,
                    2e-16,
                );
                near(
                    result.reservoir().solvent_mass_kg() + result.moisture().water_mass_kg(),
                    bath.solvent_mass_kg() + body.water_mass_kg(),
                    2e-13,
                );
                assert_eq!(
                    result.reservoir().solution().solute_moles(),
                    saved.solute_moles()
                );
                assert_eq!(body.activity(), original_activity);
                assert_eq!(body.time_s(), 0.0);
                assert_eq!(bath.solution(), &saved);
            }
        }
    });
}

fn slab_participant(id: u64, moisture: &MoistureLayeredSlab) -> LayeredMoistureParticipant<'_> {
    LayeredMoistureParticipant {
        id,
        moisture,
        exposed: [false, true],
        transport: config(),
        nonlinear: nonlinear(),
    }
}
fn radial_participant(
    id: u64,
    moisture: &MoistureLayeredRadial,
) -> LayeredRadialMoistureParticipant<'_> {
    LayeredRadialMoistureParticipant {
        id,
        moisture,
        exposed: [false, true],
        transport: config(),
        nonlinear: nonlinear(),
    }
}

#[test]
fn shared_layered_solution_preserves_signed_exchange_and_canonical_physical_bodies() {
    with_cx(|cx| {
        let dt = 2.0;
        let volume = 1e-5;
        let solute = 0.02;
        let activity = (-aqueous::WATER_VOLUME * solute / volume).exp();
        let initial = [[0.1, 0.2], [0.98, 0.99]];
        let p = coefficients(None);
        let expected = initial.map(|a| endpoint(&p, a, dt, activity));
        assert!(expected[0].1 > 0.0 && expected[1].1 < 0.0);
        let water = expected[0].1 + expected[1].1;
        let start = volume + water / WATER_DENSITY;
        let bath = aqueous::reservoir(start, vec![solute]);
        let bodies = initial.map(|a| slab_body(a, cx));
        let run = |reverse| {
            let participants = if reverse {
                [
                    slab_participant(23, &bodies[1]),
                    slab_participant(7, &bodies[0]),
                ]
            } else {
                [
                    slab_participant(7, &bodies[0]),
                    slab_participant(23, &bodies[1]),
                ]
            };
            bath.advance_shared_layered(&participants, dt, 4, aqueous::controls(start), cx)
                .unwrap()
        };
        let result = run(false);
        let reversed = run(true);
        assert_eq!(
            result.reservoir().solution(),
            reversed.reservoir().solution()
        );
        assert_eq!(
            result.solids().iter().map(|s| s.id).collect::<Vec<_>>(),
            vec![7, 23]
        );
        for (i, (body, reverse)) in result.solids().iter().zip(reversed.solids()).enumerate() {
            assert_eq!(body.moisture.activity(), reverse.moisture.activity());
            for j in 0..2 {
                near(body.moisture.activity()[j], expected[i].0[j], 2e-11);
            }
            near(body.transport.boundary_water_kg[1], expected[i].1, 2e-13);
        }
        near(result.step().water_to_solids_kg, water, 2e-13);
        near(result.step().direct_water_balance_residual_kg, 0.0, 1e-13);
        near(
            result.reservoir().solution().solvent_volume_m3(),
            volume,
            2e-16,
        );
        assert_eq!(result.reservoir().solution().solute_moles(), &[solute]);

        let symmetries = [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
        ];
        let p = symmetries.map(|s| coefficients(Some(s)));
        let expected: [_; 2] = std::array::from_fn(|i| endpoint(&p[i], initial[i], dt, activity));
        assert!(expected[0].1 > 0.0 && expected[1].1 < 0.0);
        let water = expected[0].1 + expected[1].1;
        let start = volume + water / WATER_DENSITY;
        let bath = aqueous::reservoir(start, vec![solute]);
        let bodies: [_; 2] = std::array::from_fn(|i| radial_body(symmetries[i], initial[i], cx));
        let run = |reverse| {
            let participants = if reverse {
                [
                    radial_participant(23, &bodies[1]),
                    radial_participant(7, &bodies[0]),
                ]
            } else {
                [
                    radial_participant(7, &bodies[0]),
                    radial_participant(23, &bodies[1]),
                ]
            };
            bath.advance_shared_layered_radial(&participants, dt, 4, aqueous::controls(start), cx)
                .unwrap()
        };
        let result = run(false);
        let reversed = run(true);
        assert_eq!(
            result.reservoir().solution(),
            reversed.reservoir().solution()
        );
        assert_eq!(
            result.solids().iter().map(|s| s.id).collect::<Vec<_>>(),
            vec![7, 23]
        );
        for (i, (body, reverse)) in result.solids().iter().zip(reversed.solids()).enumerate() {
            assert_eq!(body.moisture.activity(), reverse.moisture.activity());
            for j in 0..2 {
                near(body.moisture.activity()[j], expected[i].0[j], 2e-11);
            }
            near(body.transport.boundary_water_kg[1], expected[i].1, 2e-13);
        }
        near(result.step().water_to_solids_kg, water, 2e-13);
        near(result.step().direct_water_balance_residual_kg, 0.0, 1e-13);
        near(
            result.reservoir().solution().solvent_volume_m3(),
            volume,
            2e-16,
        );
        near(
            result.step().solid_water_kg[0],
            bodies
                .iter()
                .map(MoistureLayeredRadial::water_mass_kg)
                .sum(),
            2e-13,
        );
        assert_eq!(result.reservoir().solution().solute_moles(), &[solute]);
        for (body, a) in bodies.iter().zip(initial) {
            for i in 0..2 {
                near(body.activity()[i], a[i], 2e-16);
            }
            assert_eq!(body.time_s(), 0.0);
        }
    });
}

fn kinetic_bath() -> FiniteMoistureSolution {
    let law = IdealSolubilityV1::new(
        SpeciesId::new("water").unwrap(),
        SpeciesId::new("authored-layered-solute").unwrap(),
        Temperature::new(293.15),
        AmountConcentration::new(1000.0),
        0.0,
        [Temperature::new(293.15), Temperature::new(293.15)],
    )
    .unwrap();
    FiniteMoistureSolution::with_kinetic_solubility(
        KineticSolubilitySolution::new(
            Temperature::new(293.15),
            aqueous::WATER_VOLUME,
            Volume::new(5e-5),
            vec![(
                SolubilityKineticsV1::new(
                    law,
                    VolumetricFlowRate::new(2.5e-5),
                    PrecipitationPolicyV1::ExistingSolidOnly,
                )
                .unwrap(),
                SolubilityKineticStateV1::from_phases(Amount::new(0.0), Amount::new(0.2)).unwrap(),
            )],
        )
        .unwrap(),
        0.018,
    )
    .unwrap()
}

fn check_phase(bath: &FiniteMoistureSolution, duration: f64) {
    let phase = bath.kinetic_solubility().unwrap().last_steps()[0].state();
    near(
        phase.dissolved().value(),
        0.05 * (1.0 - (-0.5 * duration).exp()),
        3e-15,
    );
    near(
        phase.dissolved().value() + phase.solid().value(),
        0.2,
        2e-16,
    );
    assert_eq!(bath.solution().solvent_volume_m3(), 5e-5);
}

#[test]
fn effective_seals_keep_internal_diffusion_and_actual_kinetic_phase_continuation() {
    with_cx(|cx| {
        let bath = kinetic_bath();
        let initial = slab_body([0.8, 0.1], cx)
            .with_transfer(
                [MoistureActivityTransfer::Sealed; 2],
                vec![MoistureActivityTransfer::ConductanceKgM2S(JOIN_H)],
                cx,
            )
            .unwrap();
        let first = bath
            .advance_layered(
                &initial,
                2.0,
                [true; 2],
                config(),
                nonlinear(),
                aqueous::controls(5e-5),
                cx,
            )
            .unwrap();
        check_phase(first.reservoir(), 2.0);
        assert_eq!(first.step().iterations, 0);
        assert_eq!(first.step().transport.boundary_water_kg, [0.0; 2]);
        assert!(first.moisture().activity()[0] < 0.8);
        assert!(first.moisture().activity()[1] > 0.1);
        near(
            first.moisture().water_mass_kg(),
            initial.water_mass_kg(),
            1e-13,
        );
        let second = first
            .reservoir()
            .advance_layered(
                first.moisture(),
                2.0,
                [true; 2],
                config(),
                nonlinear(),
                aqueous::controls(5e-5),
                cx,
            )
            .unwrap();
        check_phase(second.reservoir(), 4.0);
        assert_eq!(second.moisture().time_s(), 4.0);

        let initial = radial_body(RadialMoistureSymmetry::Sphere, [0.8, 0.1], cx)
            .with_transfer(
                [MoistureActivityTransfer::Sealed; 2],
                vec![MoistureActivityTransfer::ConductanceKgM2S(JOIN_H)],
                cx,
            )
            .unwrap();
        let first = bath
            .advance_layered_radial(
                &initial,
                2.0,
                [false, true],
                config(),
                nonlinear(),
                aqueous::controls(5e-5),
                cx,
            )
            .unwrap();
        check_phase(first.reservoir(), 2.0);
        assert_eq!(first.step().iterations, 0);
        assert_eq!(first.step().transport.boundary_water_kg, [0.0; 2]);
        assert!(first.moisture().activity()[0] < 0.8);
        assert!(first.moisture().activity()[1] > 0.1);
        near(
            first.moisture().water_mass_kg(),
            initial.water_mass_kg(),
            1e-13,
        );
        let second = first
            .reservoir()
            .advance_layered_radial(
                first.moisture(),
                2.0,
                [false, true],
                config(),
                nonlinear(),
                aqueous::controls(5e-5),
                cx,
            )
            .unwrap();
        check_phase(second.reservoir(), 4.0);
        assert_eq!(second.moisture().time_s(), 4.0);
        assert_eq!(
            bath.kinetic_solubility().unwrap().last_steps()[0]
                .state()
                .dissolved()
                .value(),
            0.0
        );
        assert_eq!(
            bath.kinetic_solubility().unwrap().last_steps()[0]
                .state()
                .solid()
                .value(),
            0.2
        );
        near(initial.activity()[0], 0.8, 2e-16);
        near(initial.activity()[1], 0.1, 2e-16);
    });
}

fn pressure_stack(cx: &Cx<'_>) -> MoistureLayeredSlab {
    let card = material(500.0, vec![(0.0, 0.0), (1.0, 0.2)], [0.001; 2]);
    let pressure = |point: QueryPoint| {
        point
            .with_quantity("pressure", qty(QuantityKind::Pressure), 100_000.0)
            .unwrap()
    };
    let slab = MoistureSlab::new_kirchhoff(
        card.clone(),
        pressure(point(0.1)),
        "moisture_diffusivity",
        MaterialPropertySelection::SingleClaimOnly,
        MoistureSlabGeometry {
            length_m: 0.05,
            cross_section_m2: 0.001,
            maximum_moisture: 0.6,
            maximum_cells: 2,
        },
        vec![0.16, 0.02],
        nonlinear(),
        cx,
    )
    .unwrap();
    let sorption = isotherm_at(
        &card,
        pressure(
            QueryPoint::new()
                .with_quantity(
                    "temperature",
                    qty(QuantityKind::AbsoluteTemperature),
                    293.15,
                )
                .unwrap()
                .with_quantity(
                    "relative_humidity",
                    qty(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
                    0.5,
                )
                .unwrap(),
        ),
    );
    MoistureLayeredSlab::new(vec![MoistureLayer::new(slab, sorption, cx).unwrap()], 2, cx).unwrap()
}

#[test]
fn layered_solution_source_center_and_shared_budget_refusals_are_atomic() {
    with_cx(|cx| {
        let body = radial_body(RadialMoistureSymmetry::Sphere, [0.1, 0.2], cx);
        let volume = 1e-4;
        let model = Arc::new(
            VirialSolutionModel::new(
                293.15,
                aqueous::WATER_VOLUME,
                1.0,
                vec![1100.0],
                vec![0.002],
            )
            .unwrap(),
        );
        let bath = FiniteMoistureSolution::new(
            AqueousSolution::nonideal(model, volume, vec![0.1]).unwrap(),
            0.018,
        )
        .unwrap();
        let saved = bath.solution().clone();
        assert!(matches!(
            bath.advance_layered_radial(
                &body,
                2.0,
                [false, true],
                config(),
                nonlinear(),
                aqueous::controls(volume),
                cx
            ),
            Err(SolutionMoistureError::Solution(OsmosisError::Nonideal(
                NonidealSolutionError::OutsideDomain { .. }
            )))
        ));
        let mut admitted = aqueous::controls(volume);
        admitted.solvent_volume_bracket_m3 = [0.95 * volume, 1.05 * volume];
        let good = bath
            .advance_layered_radial(
                &body,
                2.0,
                [false, true],
                config(),
                nonlinear(),
                admitted,
                cx,
            )
            .unwrap();
        assert!(good.step().water_to_solid_kg > 0.0);
        assert!(matches!(
            bath.advance_layered_radial(
                &body,
                2.0,
                [true, false],
                config(),
                nonlinear(),
                admitted,
                cx
            ),
            Err(SolutionMoistureError::Invalid(_))
        ));
        assert!(matches!(
            bath.advance_shared_layered_radial(
                &[radial_participant(7, &body), radial_participant(7, &body)],
                2.0,
                4,
                admitted,
                cx
            ),
            Err(SolutionMoistureError::Invalid(_))
        ));
        assert!(matches!(
            bath.advance_shared_layered_radial(
                &[radial_participant(7, &body), radial_participant(9, &body)],
                2.0,
                3,
                admitted,
                cx
            ),
            Err(SolutionMoistureError::Invalid(_))
        ));
        let pressure = pressure_stack(cx);
        assert!(matches!(
            bath.advance_layered(
                &pressure,
                2.0,
                [false, true],
                config(),
                nonlinear(),
                admitted,
                cx
            ),
            Err(SolutionMoistureError::Moisture(_))
        ));
        let closed = pressure
            .with_transfer([MoistureActivityTransfer::Sealed; 2], vec![], cx)
            .unwrap();
        let sealed = bath
            .advance_layered(&closed, 2.0, [true; 2], config(), nonlinear(), admitted, cx)
            .unwrap();
        assert_eq!(sealed.step().transport.boundary_water_kg, [0.0; 2]);
        near(
            sealed.moisture().water_mass_kg(),
            closed.water_mass_kg(),
            1e-13,
        );
        assert!(sealed.moisture().activity()[0] < closed.activity()[0]);
        assert_eq!(bath.solution(), &saved);
        near(body.activity()[0], 0.1, 2e-16);
        near(body.activity()[1], 0.2, 2e-16);
        assert_eq!(body.time_s(), 0.0);
    });
}
