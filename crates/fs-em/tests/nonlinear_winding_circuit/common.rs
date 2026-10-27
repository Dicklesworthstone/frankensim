#![allow(dead_code)]
//! Manufactured two-cell magnetic modes and actual material/circuit fixtures.

#[path = "../nonlinear_eddy/common.rs"]
pub mod material;

use fs_em::{
    planar::{eddy::nonlinear::PlanarNonlinearEddy, nonlinear::PlanarNonlinearMagnetostatics},
    winding::{
        nonlinear::{
            NonlinearWindingCircuit, NonlinearWindingConfig, NonlinearWindingControls,
            NonlinearWindingInitial, NonlinearWindings,
        },
        ohmic::WindingConductor,
    },
};
use fs_matdb::MaterialCard;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::NewtonKrylovConfig;

pub fn source(saturating: bool, sigma: f64) -> MaterialCard {
    material::source(
        if saturating {
            &[(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)]
        } else {
            &[(0.0, 0.0), (4.0, 4.0)]
        },
        sigma,
    )
}

pub fn problem(card: &MaterialCard, depth: f64) -> PlanarNonlinearMagnetostatics {
    material::problem(card, mesh(depth))
}
pub fn mesh(depth: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "two-square-stranded-magnetic-section",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [2.0, 1.0],
            [0.5, 0.5],
            [1.5, 0.5],
        ],
        vec![
            [0, 1, 6],
            [1, 4, 6],
            [4, 3, 6],
            [3, 0, 6],
            [1, 2, 7],
            [2, 5, 7],
            [5, 4, 7],
            [4, 1, 7],
        ],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap()
}
pub fn boundary() -> Vec<(usize, f64)> {
    (0..6).map(|v| (v, 0.0)).collect()
}
pub fn profiles(count: usize) -> Vec<Vec<f64>> {
    match count {
        1 => vec![vec![6.0, 6.0, 6.0, 6.0, -6.0, -6.0, -6.0, -6.0]],
        2 => vec![
            vec![12.0, 12.0, 0.0, 0.0, -12.0, -12.0, 0.0, 0.0],
            vec![0.0, 0.0, 12.0, 12.0, 0.0, 0.0, -12.0, -12.0],
        ],
        _ => panic!("fixture supports one or two disjoint physical windings"),
    }
}
pub fn config() -> NonlinearWindingConfig {
    NonlinearWindingConfig {
        max_windings: 4,
        max_profile_values: 1000,
        current_balance_tolerance_turns: 1e-10,
    }
}
pub fn controls() -> NonlinearWindingControls {
    NonlinearWindingControls {
        field: material::controls(),
        maximum_field_trials: 1000,
        maximum_newton_iterations: 30,
        difference_step_a: 1e-5,
        terminal_voltage_tolerance_v: 1e-8,
        power_tolerance_w: 1e-8,
        energy_tolerance_j: 2e-8,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}
pub fn conductors(card: &MaterialCard, temperature: f64) -> Vec<WindingConductor> {
    let state = material::conductor(card, temperature, ElectricalPropertyBasis::Resistivity);
    (0..8)
        .map(|triangle| WindingConductor {
            triangle,
            conductor: state.clone(),
        })
        .collect()
}
pub fn initial_field<'m>(
    problem: &'m PlanarNonlinearMagnetostatics,
    profiles: &[Vec<f64>],
    currents: &[f64],
) -> PlanarNonlinearEddy<'m> {
    let mut j = vec![0.0; 8];
    for (profile, &i) in profiles.iter().zip(currents) {
        for (density, &p) in j.iter_mut().zip(profile) {
            *density = p.mul_add(i, *density);
        }
    }
    let initial = problem
        .solve(
            &j,
            &boundary(),
            &vec![0.0; 8],
            material::controls().nonlinear,
            &mut || false,
        )
        .unwrap();
    PlanarNonlinearEddy::new_with_impressed_current(
        &initial,
        vec![],
        &(0..8).collect::<Vec<_>>(),
        &boundary(),
        &j,
        0.0,
        material::controls(),
        &mut || false,
    )
    .unwrap()
}

pub fn single<'m>(
    problem: &'m PlanarNonlinearMagnetostatics,
    card: &MaterialCard,
    series_ohm: f64,
    initial_current_a: f64,
    initial_source_v: f64,
) -> NonlinearWindingCircuit<'m> {
    let profiles = profiles(1);
    NonlinearWindingCircuit::new(
        initial_field(problem, &profiles, &[initial_current_a]),
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: 123.0 }),
                (2, 0, Branch::Resistor { ohms: series_ohm }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: profiles,
            end_inductance_h: vec![0.0],
        },
        conductors(card, 300.0),
        NonlinearWindingInitial {
            currents_a: vec![initial_current_a],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![initial_source_v],
        },
        config(),
        &mut || false,
    )
    .unwrap()
}

pub fn h(b: f64, saturating: bool) -> f64 {
    b.signum()
        * if saturating && b.abs() > 1.0 {
            1.0 + 4.0 * (b.abs() - 1.0)
        } else {
            b.abs()
        }
}
pub fn b_from_current(i: f64, saturating: bool) -> f64 {
    i.signum()
        * if saturating && i.abs() > 1.0 {
            1.0 + (i.abs() - 1.0) / 4.0
        } else {
            i.abs()
        }
}
pub fn density_energy(b: f64, saturating: bool) -> f64 {
    let b = b.abs();
    if saturating && b > 1.0 {
        0.5 + (b - 1.0) + 2.0 * (b - 1.0) * (b - 1.0)
    } else {
        0.5 * b * b
    }
}
pub fn linkage(i: f64, depth: f64, saturating: bool) -> f64 {
    2.0 * depth * b_from_current(i, saturating)
}
pub fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-7 * (1.0 + expected.abs()),
        "{actual} != {expected}"
    );
}
