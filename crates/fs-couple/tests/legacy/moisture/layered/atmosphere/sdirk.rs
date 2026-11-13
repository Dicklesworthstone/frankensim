//! Independent finite-gas graph and exact-time references for both native stages.
#[path = "sdirk/adaptive.rs"]
mod adaptive;

use super::*;
use fs_couple::moisture::atmosphere::FiniteMoistureAtmosphere;
use fs_couple::moisture::radial::{
    MoistureRadial, MoistureRadialGeometry, RadialMoistureSymmetry,
    layered::{
        MoistureLayeredRadial, MoistureRadialLayer, atmosphere::LayeredRadialMoistureParticipant,
    },
};

const GAMMA: f64 = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
const SLAB_LINKS: [(usize, usize, f64); 3] = [(0, 1, 0.012), (0, 2, 0.016), (1, 2, 0.048)];

/// Independent conservative physical storage graph; production uses neither
/// this dense matrix nor this elimination.
fn rate_matrix<const N: usize>(capacity: [f64; N], links: &[(usize, usize, f64)]) -> [[f64; N]; N] {
    let mut rate = [[0.0; N]; N];
    for &(left, right, conductance) in links {
        rate[left][left] -= conductance / capacity[left];
        rate[left][right] += conductance / capacity[left];
        rate[right][right] -= conductance / capacity[right];
        rate[right][left] += conductance / capacity[right];
    }
    rate
}
fn product<const N: usize>(matrix: &[[f64; N]; N], state: [f64; N]) -> [f64; N] {
    std::array::from_fn(|i| (0..N).map(|j| matrix[i][j] * state[j]).sum())
}
fn solve<const N: usize>(mut matrix: [[f64; N]; N], mut rhs: [f64; N]) -> [f64; N] {
    for row in 0..N {
        let pivot = (row..N)
            .max_by(|&a, &b| matrix[a][row].abs().total_cmp(&matrix[b][row].abs()))
            .unwrap();
        matrix.swap(row, pivot);
        rhs.swap(row, pivot);
        assert!(matrix[row][row].abs() > 0.0);
        for next in row + 1..N {
            let ratio = matrix[next][row] / matrix[row][row];
            for column in row..N {
                matrix[next][column] -= ratio * matrix[row][column];
            }
            rhs[next] -= ratio * rhs[row];
        }
    }
    let mut result = [0.0; N];
    for row in (0..N).rev() {
        result[row] = (rhs[row]
            - (row + 1..N)
                .map(|column| matrix[row][column] * result[column])
                .sum::<f64>())
            / matrix[row][row];
    }
    result
}
fn rk_stages<const N: usize>(
    rate: &[[f64; N]; N],
    initial: [f64; N],
    dt: f64,
) -> ([f64; N], [f64; N]) {
    let implicit: [[f64; N]; N] = std::array::from_fn(|i| {
        std::array::from_fn(|j| (if i == j { 1.0 } else { 0.0 }) - GAMMA * dt * rate[i][j])
    });
    let first = solve(implicit, initial);
    let first_rate = product(rate, first);
    let rhs = std::array::from_fn(|i| initial[i] + (1.0 - GAMMA) * dt * first_rate[i]);
    (first, solve(implicit, rhs))
}
fn weighted(dt: f64, conductance: f64, first_difference: f64, final_difference: f64) -> f64 {
    dt * conductance * ((1.0 - GAMMA) * first_difference + GAMMA * final_difference)
}

/// Closed form of a conservative three-store graph: one equilibrium mode and
/// two exponential relaxation modes. This is independent of either time stepper.
fn exact_three(capacity: [f64; 3], rate: &[[f64; 3]; 3], initial: [f64; 3], time: f64) -> [f64; 3] {
    let mean = (0..3).map(|i| capacity[i] * initial[i]).sum::<f64>() / capacity.iter().sum::<f64>();
    let trace = (0..3).map(|i| rate[i][i]).sum::<f64>();
    let trace_square = (0..3)
        .map(|i| (0..3).map(|j| rate[i][j] * rate[j][i]).sum::<f64>())
        .sum::<f64>();
    let product = 0.5 * (trace * trace - trace_square);
    let gap = (trace * trace - 4.0 * product).sqrt();
    let slow = (-trace - gap) / 2.0;
    let fast = (-trace + gap) / 2.0;
    let derivative = self::product(rate, initial);
    std::array::from_fn(|i| {
        let difference = initial[i] - mean;
        let a = (derivative[i] + fast * difference) / gap;
        mean + a * (-slow * time).exp() + (difference - a) * (-fast * time).exp()
    })
}

fn gas_energy(atmosphere: &FiniteMoistureAtmosphere) -> f64 {
    atmosphere.gas().mass_kg() * atmosphere.gas().state().specific_internal_energy_j_kg()
}
fn slab_participant(id: u64, moisture: &MoistureLayeredSlab) -> LayeredMoistureParticipant<'_> {
    LayeredMoistureParticipant {
        id,
        moisture,
        exposed: [true; 2],
        transport: config(),
        nonlinear: nonlinear(),
    }
}

#[test]
fn finite_slab_stages_match_coupled_rk_and_whole_interval_water_and_gas_energy() {
    with_cx(|cx| {
        for (rh, initial) in [(0.8, [0.2, 0.1]), (0.1, [0.8, 0.5])] {
            let mut atmosphere = air::atmosphere(rh, 1.0);
            let mut body = stack(initial, cx);
            let gas_capacity = atmosphere.vapor_mass_kg().unwrap() / rh;
            let rate = rate_matrix([CAPACITY[0], CAPACITY[1], gas_capacity], &SLAB_LINKS);
            let carrier = air::carrier_mass(&atmosphere);
            for _ in 0..3 {
                let dt = 0.25;
                let old = [
                    body.activity()[0],
                    body.activity()[1],
                    atmosphere.relative_humidity().unwrap(),
                ];
                let (first, last) = rk_stages(&rate, old, dt);
                let candidate = atmosphere
                    .advance_layered_sdirk2(
                        &body,
                        dt,
                        [true; 2],
                        config(),
                        nonlinear(),
                        air::controls(),
                        cx,
                    )
                    .unwrap();
                let report = candidate.step();
                let stage = &report.transport.first_stage;
                let endpoint = &report.transport.endpoint;
                near(
                    candidate.atmosphere().relative_humidity().unwrap(),
                    last[2],
                    3e-10,
                );
                assert_eq!(endpoint.interval_s, [body.time_s(), body.time_s() + dt]);
                near(stage.interval_s[1], body.time_s() + GAMMA * dt, 2e-15);
                for i in 0..2 {
                    let conductance = [0.016, 0.048][i];
                    near(candidate.moisture().activity()[i], last[i], 3e-10);
                    near(stage.surface_activities[i].unwrap(), first[2], 3e-10);
                    near(
                        stage.boundary_water_kg[i],
                        GAMMA * dt * conductance * (first[2] - first[i]),
                        2e-12,
                    );
                    near(
                        endpoint.boundary_water_kg[i],
                        weighted(dt, conductance, first[2] - first[i], last[2] - last[i]),
                        2e-12,
                    );
                    let source = report.sorption[i].as_ref().unwrap();
                    near(source.query().axes()["relative_humidity"], last[2], 3e-10);
                    near(
                        source.property("sorption").unwrap().value_si(),
                        SORPTION[i] * last[2],
                        3e-10,
                    );
                    assert_eq!(
                        source.card_identity(),
                        body.layers()[i].slab().card().content_hash()
                    );
                }
                let internal = weighted(dt, 0.012, first[0] - first[1], last[0] - last[1]);
                near(
                    endpoint.interfaces[0].water_left_to_right_kg,
                    internal,
                    2e-12,
                );
                near(
                    stage.interfaces[0].water_left_to_right_kg,
                    GAMMA * dt * 0.012 * (first[0] - first[1]),
                    2e-12,
                );
                near(
                    CAPACITY[0] * (last[0] - old[0]),
                    endpoint.boundary_water_kg[0] - internal,
                    2e-12,
                );
                near(
                    CAPACITY[1] * (last[1] - old[1]),
                    endpoint.boundary_water_kg[1] + internal,
                    2e-12,
                );
                near(
                    report.water_to_solid_kg,
                    endpoint.boundary_water_kg.iter().sum(),
                    2e-12,
                );
                near(
                    candidate.atmosphere().vapor_mass_kg().unwrap()
                        + candidate.moisture().water_mass_kg(),
                    atmosphere.vapor_mass_kg().unwrap() + body.water_mass_kg(),
                    2e-13,
                );
                near(
                    candidate.atmosphere().gas().mass_kg() + candidate.moisture().wet_mass_kg(),
                    atmosphere.gas().mass_kg() + body.wet_mass_kg(),
                    2e-13,
                );
                near(air::carrier_mass(candidate.atmosphere()), carrier, 2e-13);
                near(
                    report.gas_internal_energy_change_j,
                    gas_energy(candidate.atmosphere()) - gas_energy(&atmosphere),
                    2e-10,
                );
                let (next_air, next_body, _) = candidate.into_parts();
                atmosphere = next_air;
                body = next_body;
            }
        }
    });
}

#[test]
fn finite_slab_sdirk_refines_at_second_order_against_the_continuous_closed_system() {
    with_cx(|cx| {
        let initial = [0.2, 0.1, 0.8];
        let gas = air::atmosphere(initial[2], 1.0);
        let capacity = [
            CAPACITY[0],
            CAPACITY[1],
            gas.vapor_mass_kg().unwrap() / initial[2],
        ];
        let rate = rate_matrix(capacity, &SLAB_LINKS);
        let exact = exact_three(capacity, &rate, initial, 1.0);
        let errors = [4, 8, 16].map(|count| {
            let mut atmosphere = gas.clone();
            let mut body = stack([initial[0], initial[1]], cx);
            for _ in 0..count {
                let next = atmosphere
                    .advance_layered_sdirk2(
                        &body,
                        1.0 / f64::from(count),
                        [true; 2],
                        config(),
                        nonlinear(),
                        air::controls(),
                        cx,
                    )
                    .unwrap();
                (atmosphere, body, _) = next.into_parts();
            }
            let state = [
                body.activity()[0],
                body.activity()[1],
                atmosphere.relative_humidity().unwrap(),
            ];
            (0..3)
                .map(|i| (state[i] - exact[i]).abs())
                .fold(0.0, f64::max)
        });
        for pair in errors.windows(2) {
            assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "{errors:?}");
        }
        assert!(errors[2] < 4e-5, "{errors:?}");
    });
}

const RHO: [f64; 2] = [300.0, 900.0];
const SLOPE: [f64; 2] = [0.2, 0.4];
const DIFFUSION: [f64; 2] = [0.001, 0.0005];
fn coated(
    symmetry: RadialMoistureSymmetry,
    activity: [f64; 2],
    cx: &Cx<'_>,
) -> MoistureLayeredRadial {
    let layers = (0..2)
        .map(|i| {
            let card = material(RHO[i], vec![(0.0, 0.0), (1.0, SLOPE[i])], [DIFFUSION[i]; 2]);
            let sorption = isotherm(&card);
            let radial = MoistureRadial::new(
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
                vec![SLOPE[i] * activity[i]],
                cx,
            )
            .unwrap();
            MoistureRadialLayer::new(radial, sorption, cx).unwrap()
        })
        .collect();
    MoistureLayeredRadial::new(layers, 2, cx)
        .unwrap()
        .with_transfer(
            [
                MoistureActivityTransfer::IdealContact,
                MoistureActivityTransfer::ConductanceKgM2S(0.04),
            ],
            vec![MoistureActivityTransfer::ConductanceKgM2S(0.03)],
            cx,
        )
        .unwrap()
}
struct Curved {
    capacity: [f64; 2],
    join: f64,
    outer: f64,
    outer_bulk: f64,
}
fn curved(symmetry: RadialMoistureSymmetry) -> Curved {
    let pi = core::f64::consts::PI;
    let volume = |a: f64, b: f64| match symmetry {
        RadialMoistureSymmetry::Sphere => 4.0 * pi * (b.powi(3) - a.powi(3)) / 3.0,
        RadialMoistureSymmetry::Cylinder { length_m } => pi * length_m * (b * b - a * a),
    };
    let path = |a: f64, b: f64| match symmetry {
        RadialMoistureSymmetry::Sphere => (1.0 / a - 1.0 / b) / (4.0 * pi),
        RadialMoistureSymmetry::Cylinder { length_m } => (b / a).ln() / (2.0 * pi * length_m),
    };
    let area = |r: f64| match symmetry {
        RadialMoistureSymmetry::Sphere => 4.0 * pi * r * r,
        RadialMoistureSymmetry::Cylinder { length_m } => 2.0 * pi * length_m * r,
    };
    let half = [
        RHO[0] * DIFFUSION[0] * SLOPE[0] / path(0.015, 0.03),
        RHO[1] * DIFFUSION[1] * SLOPE[1] / path(0.03, 0.04),
    ];
    let outer_bulk = RHO[1] * DIFFUSION[1] * SLOPE[1] / path(0.04, 0.05);
    Curved {
        capacity: [
            RHO[0] * SLOPE[0] * volume(0.0, 0.03),
            RHO[1] * SLOPE[1] * volume(0.03, 0.05),
        ],
        join: 1.0 / (1.0 / half[0] + 1.0 / (area(0.03) * 0.03) + 1.0 / half[1]),
        outer: 1.0 / (1.0 / outer_bulk + 1.0 / (area(0.05) * 0.04)),
        outer_bulk,
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
fn finite_coated_cylinder_and_sphere_retain_curved_resistance_and_weighted_storage() {
    with_cx(|cx| {
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
        ] {
            let p = curved(symmetry);
            for (rh, initial) in [(0.8, [0.1, 0.2]), (0.15, [0.8, 0.6])] {
                let mut atmosphere = air::atmosphere(rh, 0.5);
                let mut body = coated(symmetry, initial, cx);
                let rate = rate_matrix(
                    [
                        p.capacity[0],
                        p.capacity[1],
                        atmosphere.vapor_mass_kg().unwrap() / rh,
                    ],
                    &[(0, 1, p.join), (1, 2, p.outer)],
                );
                for _ in 0..2 {
                    let dt = 2.0;
                    let old = [
                        body.activity()[0],
                        body.activity()[1],
                        atmosphere.relative_humidity().unwrap(),
                    ];
                    let (first, last) = rk_stages(&rate, old, dt);
                    let candidate = atmosphere
                        .advance_layered_radial_sdirk2(
                            &body,
                            dt,
                            [false, true],
                            config(),
                            nonlinear(),
                            air::controls(),
                            cx,
                        )
                        .unwrap();
                    let report = candidate.step();
                    let endpoint = &report.transport.endpoint;
                    near(
                        candidate.atmosphere().relative_humidity().unwrap(),
                        last[2],
                        4e-11,
                    );
                    for i in 0..2 {
                        near(candidate.moisture().activity()[i], last[i], 4e-12);
                    }
                    let exchange = weighted(dt, p.outer, first[2] - first[1], last[2] - last[1]);
                    let internal = weighted(dt, p.join, first[0] - first[1], last[0] - last[1]);
                    near(endpoint.boundary_water_kg[1], exchange, 2e-12);
                    assert_eq!(endpoint.boundary_water_kg[0], 0.0);
                    near(
                        endpoint.interfaces[0].water_left_to_right_kg,
                        internal,
                        2e-12,
                    );
                    near(p.capacity[0] * (last[0] - old[0]), -internal, 2e-12);
                    near(
                        p.capacity[1] * (last[1] - old[1]),
                        internal + exchange,
                        2e-12,
                    );
                    near(
                        report.transport.first_stage.boundary_water_kg[1],
                        GAMMA * dt * p.outer * (first[2] - first[1]),
                        2e-12,
                    );
                    let surface = last[1] + p.outer * (last[2] - last[1]) / p.outer_bulk;
                    near(endpoint.surface_activities[1].unwrap(), surface, 4e-12);
                    near(
                        report.transport.first_stage.surface_activities[1].unwrap(),
                        first[1] + p.outer * (first[2] - first[1]) / p.outer_bulk,
                        4e-12,
                    );
                    let source = report.sorption[1].as_ref().unwrap();
                    near(source.query().axes()["relative_humidity"], surface, 4e-12);
                    assert_eq!(
                        source
                            .property("sorption")
                            .unwrap()
                            .answer()
                            .receipt
                            .selected,
                        body.layers()[1].sorption().selected_claim()
                    );
                    assert!(report.sorption[0].is_none());
                    near(report.water_to_solid_kg, exchange, 2e-12);
                    near(
                        candidate.atmosphere().vapor_mass_kg().unwrap()
                            + candidate.moisture().water_mass_kg(),
                        atmosphere.vapor_mass_kg().unwrap() + body.water_mass_kg(),
                        2e-13,
                    );
                    near(
                        candidate.atmosphere().gas().mass_kg() + candidate.moisture().wet_mass_kg(),
                        atmosphere.gas().mass_kg() + body.wet_mass_kg(),
                        2e-13,
                    );
                    near(
                        air::carrier_mass(candidate.atmosphere()),
                        air::carrier_mass(&atmosphere),
                        2e-13,
                    );
                    near(
                        report.gas_internal_energy_change_j,
                        gas_energy(candidate.atmosphere()) - gas_energy(&atmosphere),
                        2e-10,
                    );
                    let (next_air, next_body, _) = candidate.into_parts();
                    atmosphere = next_air;
                    body = next_body;
                }
            }
        }
    });
}

#[test]
fn shared_slab_stages_match_five_stores_and_keep_canonical_opposing_transfers() {
    with_cx(|cx| {
        let atmosphere = air::atmosphere(0.5, 1.0);
        let dry = stack([0.1, 0.2], cx);
        let wet = stack([0.8, 0.7], cx);
        let rate = rate_matrix(
            [
                CAPACITY[0],
                CAPACITY[1],
                CAPACITY[0],
                CAPACITY[1],
                atmosphere.vapor_mass_kg().unwrap() / 0.5,
            ],
            &[
                (0, 1, 0.012),
                (0, 4, 0.016),
                (1, 4, 0.048),
                (2, 3, 0.012),
                (2, 4, 0.016),
                (3, 4, 0.048),
            ],
        );
        let (_, expected) = rk_stages(&rate, [0.1, 0.2, 0.8, 0.7, 0.5], 0.25);
        let candidate = atmosphere
            .advance_shared_layered_sdirk2(
                &[slab_participant(7, &dry), slab_participant(2, &wet)],
                0.25,
                4,
                air::controls(),
                cx,
            )
            .unwrap();
        let reversed = atmosphere
            .advance_shared_layered_sdirk2(
                &[slab_participant(2, &wet), slab_participant(7, &dry)],
                0.25,
                4,
                air::controls(),
                cx,
            )
            .unwrap();
        assert_eq!(candidate.atmosphere(), reversed.atmosphere());
        assert_eq!(
            candidate.solids().iter().map(|s| s.id).collect::<Vec<_>>(),
            [2, 7]
        );
        near(
            candidate.atmosphere().relative_humidity().unwrap(),
            expected[4],
            3e-10,
        );
        for (solid, reverse) in candidate.solids().iter().zip(reversed.solids()) {
            assert_eq!(solid.moisture.activity(), reverse.moisture.activity());
            let (offset, original) = if solid.id == 7 { (0, &dry) } else { (2, &wet) };
            for i in 0..2 {
                near(solid.moisture.activity()[i], expected[offset + i], 3e-10);
            }
            let change = solid.moisture.water_mass_kg() - original.water_mass_kg();
            assert_eq!(change.is_sign_positive(), solid.id == 7);
            near(
                change,
                solid.transport.endpoint.boundary_water_kg.iter().sum(),
                2e-12,
            );
        }
        near(
            candidate.atmosphere().vapor_mass_kg().unwrap()
                + candidate
                    .solids()
                    .iter()
                    .map(|s| s.moisture.water_mass_kg())
                    .sum::<f64>(),
            atmosphere.vapor_mass_kg().unwrap() + dry.water_mass_kg() + wet.water_mass_kg(),
            2e-13,
        );
        near(candidate.step().total_mass_residual_kg, 0.0, 1e-13);
        near(
            candidate.step().gas_internal_energy_change_j,
            gas_energy(candidate.atmosphere()) - gas_energy(&atmosphere),
            2e-10,
        );
        assert_eq!(dry.time_s(), 0.0);
        assert_eq!(wet.time_s(), 0.0);
    });
}

#[test]
fn shared_curved_stages_keep_actual_shell_capacities_and_one_original_gas() {
    with_cx(|cx| {
        let symmetries = [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
        ];
        let p = symmetries.map(curved);
        let dry = coated(symmetries[0], [0.1, 0.2], cx);
        let wet = coated(symmetries[1], [0.8, 0.7], cx);
        let atmosphere = air::atmosphere(0.4, 1.0);
        let rate = rate_matrix(
            [
                p[0].capacity[0],
                p[0].capacity[1],
                p[1].capacity[0],
                p[1].capacity[1],
                atmosphere.vapor_mass_kg().unwrap() / 0.4,
            ],
            &[
                (0, 1, p[0].join),
                (1, 4, p[0].outer),
                (2, 3, p[1].join),
                (3, 4, p[1].outer),
            ],
        );
        let (_, expected) = rk_stages(&rate, [0.1, 0.2, 0.8, 0.7, 0.4], 2.0);
        let candidate = atmosphere
            .advance_shared_layered_radial_sdirk2(
                &[radial_participant(7, &dry), radial_participant(2, &wet)],
                2.0,
                4,
                air::controls(),
                cx,
            )
            .unwrap();
        let reversed = atmosphere
            .advance_shared_layered_radial_sdirk2(
                &[radial_participant(2, &wet), radial_participant(7, &dry)],
                2.0,
                4,
                air::controls(),
                cx,
            )
            .unwrap();
        assert_eq!(candidate.atmosphere(), reversed.atmosphere());
        assert_eq!(
            candidate.solids().iter().map(|s| s.id).collect::<Vec<_>>(),
            [2, 7]
        );
        near(
            candidate.atmosphere().relative_humidity().unwrap(),
            expected[4],
            4e-11,
        );
        for (solid, reverse) in candidate.solids().iter().zip(reversed.solids()) {
            assert_eq!(solid.moisture.activity(), reverse.moisture.activity());
            let (offset, original) = if solid.id == 7 { (0, &dry) } else { (2, &wet) };
            for i in 0..2 {
                near(solid.moisture.activity()[i], expected[offset + i], 4e-12);
            }
            let change = solid.moisture.water_mass_kg() - original.water_mass_kg();
            assert_eq!(change.is_sign_positive(), solid.id == 7);
            near(change, solid.transport.endpoint.boundary_water_kg[1], 2e-12);
        }
        near(
            candidate.atmosphere().vapor_mass_kg().unwrap()
                + candidate
                    .solids()
                    .iter()
                    .map(|s| s.moisture.water_mass_kg())
                    .sum::<f64>(),
            atmosphere.vapor_mass_kg().unwrap() + dry.water_mass_kg() + wet.water_mass_kg(),
            2e-13,
        );
        near(
            air::carrier_mass(candidate.atmosphere()),
            air::carrier_mass(&atmosphere),
            2e-13,
        );
        near(candidate.step().total_mass_residual_kg, 0.0, 1e-13);
        near(
            candidate.step().gas_internal_energy_change_j,
            gas_energy(candidate.atmosphere()) - gas_energy(&atmosphere),
            2e-10,
        );
    });
}

#[test]
fn finite_air_late_stage_and_participant_refusals_preserve_every_original_owner() {
    with_cx(|cx| {
        let atmosphere = air::atmosphere(0.5, 0.5);
        let saved = atmosphere.clone();
        let body = coated(RadialMoistureSymmetry::Sphere, [0.8, 0.1], cx);
        let original_activity = body.activity().to_vec();
        let p = curved(RadialMoistureSymmetry::Sphere);
        let dt = 10.0 / (p.join * (1.0 / p.capacity[0] + 1.0 / p.capacity[1]));
        atmosphere
            .advance_layered_radial(
                &body,
                GAMMA * dt,
                [false; 2],
                config(),
                nonlinear(),
                air::controls(),
                cx,
            )
            .unwrap();
        assert!(matches!(
            atmosphere.advance_layered_radial_sdirk2(
                &body,
                dt,
                [false; 2],
                config(),
                nonlinear(),
                air::controls(),
                cx,
            ),
            Err(AtmosphereMoistureError::Moisture(
                MoistureError::ActivityStageRange { .. }
            ))
        ));
        let retry = atmosphere
            .advance_layered_radial_sdirk2(
                &body,
                dt / 50.0,
                [false; 2],
                config(),
                nonlinear(),
                air::controls(),
                cx,
            )
            .unwrap();
        assert_eq!(retry.atmosphere(), &atmosphere);
        near(
            retry.moisture().water_mass_kg(),
            body.water_mass_kg(),
            2e-13,
        );
        assert!(matches!(
            atmosphere.advance_layered_radial_sdirk2(
                &body,
                2.0,
                [true, false],
                config(),
                nonlinear(),
                air::controls(),
                cx,
            ),
            Err(AtmosphereMoistureError::Invalid(_))
        ));
        assert!(matches!(
            atmosphere.advance_shared_layered_radial_sdirk2(
                &[radial_participant(1, &body), radial_participant(1, &body)],
                2.0,
                4,
                air::controls(),
                cx,
            ),
            Err(AtmosphereMoistureError::Invalid(_))
        ));
        assert!(matches!(
            atmosphere.advance_shared_layered_radial_sdirk2(
                &[radial_participant(1, &body)],
                2.0,
                1,
                air::controls(),
                cx,
            ),
            Err(AtmosphereMoistureError::Invalid(_))
        ));
        let slab = stack([0.2, 0.1], cx);
        let mut limited = air::controls();
        limited.maximum_iterations = 1;
        assert!(matches!(
            atmosphere.advance_layered_sdirk2(
                &slab,
                0.25,
                [true; 2],
                config(),
                nonlinear(),
                limited,
                cx,
            ),
            Err(AtmosphereMoistureError::NotConverged { .. })
        ));
        assert_eq!(atmosphere, saved);
        assert_eq!(body.activity(), original_activity);
        assert_eq!(body.time_s(), 0.0);
        assert_eq!(slab.time_s(), 0.0);
    });
}
