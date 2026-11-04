//! Independent binary mixing solutions through the public material parcel API.
//! All caloric/transport cards are synthetic; no measured-gas claim is made.

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
#[path = "../examples/material_gas_spring/source.rs"]
mod source;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::{
    GasSpeciesExchangeSolve,
    diffusion::{
        BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY, BinaryDiffusionMaterial,
        GasDiffusionGeometry, GasDiffusionSolve,
    },
};
use fs_material::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve, parcel::ClosedGasParcel,
};
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    Composition, CompositionBasis, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

const T: f64 = 400.0;
const P: f64 = 100_000.0;
const D: f64 = 2.0e-5;

fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn controls() -> GasDiffusionSolve {
    let gas = GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    };
    GasDiffusionSolve {
        maximum_temperature_difference_k: 1e-5,
        maximum_pressure_difference_pa: 0.01,
        gradient_tolerance_per_m: 1e-9,
        molar_flux_tolerance_mol_m2_s: 1e-12,
        entropy_tolerance_j_k: 1e-10,
        exchange: GasSpeciesExchangeSolve {
            gas: [gas; 2],
            maximum_species_fraction: 0.9,
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-8,
        },
    }
}

fn geometry() -> GasDiffusionGeometry {
    GasDiffusionGeometry {
        area_m2: 0.01,
        distance_m: 0.01,
    }
}

fn gas(x: f64, volume: f64) -> ClosedGasParcel {
    let model = source::load(
        include_str!("../examples/material_gas_spring/mixture-a.tsv"),
        "binary exponential diffusion regression",
    )
    .unwrap()
    .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![x, 1.0 - x]).unwrap())
    .unwrap();
    ClosedGasParcel::new(model, source::point(T, P).unwrap(), volume).unwrap()
}

fn first(parcel: &ClosedGasParcel) -> f64 {
    parcel
        .model()
        .thermodynamics()
        .components()
        .iter()
        .zip(
            parcel
                .model()
                .thermodynamics()
                .mole_composition()
                .fractions(),
        )
        .find(|(s, _)| s.species().as_str() == "A")
        .map_or(0.0, |(_, x)| *x)
}

fn card(d: f64, range: [f64; 2]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(range[0], d), (range[1], d)],
                dims: BINARY_DIFFUSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    typed(QuantityKind::AbsoluteTemperature),
                    range[0],
                    range[1],
                )
                .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e7),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured binary diffusion coefficient".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "binary test pair".into(),
            phase: "gas".into(),
            process: "analytic regression".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn pair<'a>(
    names: [&str; 2],
    material: &'a MaterialCard,
    point: &'a QueryPoint,
) -> BinaryDiffusionMaterial<'a> {
    BinaryDiffusionMaterial {
        species: names.map(|s| SpeciesId::new(s).unwrap()),
        material,
        context: point,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e} (allowance {tolerance:e})"
    );
}

#[test]
fn unequal_volume_binary_pair_has_analytic_endpoint_beyond_explicit_stability() {
    let (left, right) = (gas(0.75, 0.001), gas(0.25, 0.004));
    let material = card(D, [200.0, 1800.0]);
    let point = source::point(T, P).unwrap();
    let pairs = [pair(["A", "B"], &material, &point)];
    let g = geometry();
    let rate = D * g.area_m2 / g.distance_m * (1.0 / 0.001 + 1.0 / 0.004);
    let equilibrium = (0.001 * 0.75 + 0.004 * 0.25) / 0.005;
    for duration in [1e-4, 25.0, 1e8] {
        let candidate = left
            .diffuse_binary_exponential_with(&right, g, duration, &pairs, controls(), || false)
            .unwrap();
        let decay = (-rate * duration).exp();
        near(
            first(candidate.exchange().left()),
            equilibrium + 0.4 * decay,
            2e-13,
        );
        near(
            first(candidate.exchange().right()),
            equilibrium - 0.1 * decay,
            2e-13,
        );
        near(
            candidate.step().relaxation_fraction,
            -(-rate * duration).exp_m1(),
            2e-14,
        );
        assert_eq!(candidate.step().duration_s, duration);
        assert_eq!(candidate.initial_flux().pairs[0].diffusivity_m2_s, D);
        assert!(candidate.entropy_change_j_k() >= -controls().entropy_tolerance_j_k);
        assert!(candidate.exchange().step().maximum_species_mass_residual_kg < 1e-14);
        assert!(candidate.exchange().step().energy_residual_j.abs() < 1e-8);
        for p in [candidate.exchange().left(), candidate.exchange().right()] {
            near(p.state().state().temperature, T, 1e-7);
            near(p.state().state().pressure, P, 1e-5);
        }
        assert!(
            candidate.step().first_species_transfer_mol
                < candidate.initial_flux().species[0].molar_flux_mol_m2_s * g.area_m2 * duration
        );
    }
    assert!(
        left.diffuse_with(&right, g, 1e8, &pairs, controls(), || false)
            .is_err()
    );
}

#[test]
fn isothermal_semigroup_and_complete_material_observations_are_preserved() {
    let (left, right) = (gas(0.75, 0.001), gas(0.25, 0.001));
    let material = card(D, [200.0, 1800.0]);
    let point = source::point(T, P).unwrap();
    let pairs = [pair(["A", "B"], &material, &point)];
    let whole = left
        .diffuse_binary_exponential_with(&right, geometry(), 40.0, &pairs, controls(), || false)
        .unwrap();
    let half = left
        .diffuse_binary_exponential_with(&right, geometry(), 20.0, &pairs, controls(), || false)
        .unwrap();
    let second = half
        .exchange()
        .left()
        .diffuse_binary_exponential_with(
            half.exchange().right(),
            geometry(),
            20.0,
            &pairs,
            controls(),
            || false,
        )
        .unwrap();
    for (a, b) in [
        (whole.exchange().left(), second.exchange().left()),
        (whole.exchange().right(), second.exchange().right()),
    ] {
        near(first(a), first(b), 2e-13);
        near(a.mass_kg(), b.mass_kg(), 1e-14);
        near(
            a.state().state().sound_speed,
            b.state().state().sound_speed,
            1e-7,
        );
    }
    assert_ne!(
        whole.exchange().left().state().state().sound_speed,
        left.state().state().sound_speed
    );
    assert_ne!(
        whole.exchange().left().state().state().thermal_conductivity,
        left.state().state().thermal_conductivity
    );
    // Unequal species molar masses mean molar counterflow is NOT zero mass flow.
    assert!(whole.exchange().step().net_mass_left_to_right_kg.abs() > 1e-5);
}

#[test]
fn pure_species_mix_to_finite_entropy_without_a_seed_or_clipping() {
    let (left, right) = (gas(1.0, 0.001), gas(0.0, 0.001));
    let material = card(D, [200.0, 1800.0]);
    let point = source::point(T, P).unwrap();
    let pairs = [pair(["B", "A"], &material, &point)];
    let candidate = left
        .diffuse_binary_exponential_with(&right, geometry(), 1e8, &pairs, controls(), || false)
        .unwrap();
    near(first(candidate.exchange().left()), 0.5, 2e-14);
    near(first(candidate.exchange().right()), 0.5, 2e-14);
    near(
        candidate.entropy_change_j_k(),
        2.0 * P * 0.001 / T * 2.0_f64.ln(),
        1e-10,
    );
    near(
        candidate.step().first_species_transfer_mol,
        0.5 * P * 0.001 / (R * T),
        1e-15,
    );
    assert_eq!(candidate.step().relaxation_fraction, 1.0);
    assert_eq!(
        candidate
            .exchange()
            .left()
            .model()
            .thermodynamics()
            .components()
            .len(),
        2
    );
    assert_eq!(
        candidate
            .exchange()
            .right()
            .model()
            .thermodynamics()
            .components()
            .len(),
        2
    );
}

#[test]
fn reverse_orientation_diffusivity_time_scaling_and_weak_exchange_match() {
    let (left, right) = (gas(0.75, 0.001), gas(0.25, 0.004));
    let point = source::point(T, P).unwrap();
    let material = card(D, [200.0, 1800.0]);
    let faster = card(2.0 * D, [200.0, 1800.0]);
    let run = |a: &ClosedGasParcel, b: &ClosedGasParcel, d: &MaterialCard, dt| {
        a.diffuse_binary_exponential_with(
            b,
            geometry(),
            dt,
            &[pair(["A", "B"], d, &point)],
            controls(),
            || false,
        )
        .unwrap()
    };
    let ordinary = run(&left, &right, &material, 30.0);
    let scaled = run(&left, &right, &faster, 15.0);
    let reversed = run(&right, &left, &material, 30.0);
    near(
        first(ordinary.exchange().left()),
        first(scaled.exchange().left()),
        1e-14,
    );
    near(
        first(ordinary.exchange().left()),
        first(reversed.exchange().right()),
        1e-14,
    );
    near(
        ordinary.step().first_species_transfer_mol,
        -reversed.step().first_species_transfer_mol,
        1e-15,
    );
    let weak = run(&left, &right, &material, 1e-12);
    let linear = P / (R * T) * D * geometry().area_m2 / geometry().distance_m * 0.5 * 1e-12;
    assert!(weak.step().first_species_transfer_mol > 0.0);
    near(weak.step().first_species_transfer_mol / linear, 1.0, 5e-13);
}

#[test]
fn closed_and_equilibrated_faces_are_exact_identities() {
    let (left, right) = (gas(0.75, 0.001), gas(0.25, 0.004));
    let closed = left
        .diffuse_binary_exponential_with(
            &right,
            GasDiffusionGeometry {
                area_m2: 0.0,
                ..geometry()
            },
            1e8,
            &[],
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(closed.exchange().left(), &left);
    assert_eq!(closed.exchange().right(), &right);
    assert_eq!(closed.step().first_species_transfer_mol, 0.0);
    let material = card(D, [200.0, 1800.0]);
    let point = source::point(T, P).unwrap();
    let equilibrium = left
        .diffuse_binary_exponential_with(
            &left,
            geometry(),
            1e8,
            &[pair(["A", "B"], &material, &point)],
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(equilibrium.exchange().left(), &left);
    assert_eq!(equilibrium.exchange().right(), &left);
    assert_eq!(equilibrium.entropy_change_j_k(), 0.0);
}

#[test]
fn original_withdrawal_source_and_cancellation_gates_remain_atomic() {
    let (left, right) = (gas(0.75, 0.001), gas(0.25, 0.001));
    let saved = (left.clone(), right.clone());
    let material = card(D, [200.0, 1800.0]);
    let point = source::point(T, P).unwrap();
    let pairs = [pair(["A", "B"], &material, &point)];
    let mut limited = controls();
    limited.exchange.maximum_species_fraction = 0.1;
    assert!(
        left.diffuse_binary_exponential_with(&right, geometry(), 1e8, &pairs, limited, || false)
            .is_err()
    );
    for duration in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            left.diffuse_binary_exponential_with(
                &right,
                geometry(),
                duration,
                &pairs,
                controls(),
                || false
            )
            .is_err()
        );
    }
    let hot = left.at_temperature(T + 0.02).unwrap().0;
    let narrow = card(D, [399.99, 400.011]);
    let mut permissive = controls();
    permissive.maximum_temperature_difference_k = 1.0;
    permissive.maximum_pressure_difference_pa = 10.0;
    assert!(
        left.diffuse_binary_exponential_with(
            &hot,
            geometry(),
            1.0,
            &[pair(["A", "B"], &narrow, &point)],
            permissive,
            || false
        )
        .is_err()
    );
    let mut calls = 0;
    left.diffuse_binary_exponential_with(&right, geometry(), 25.0, &pairs, controls(), || {
        calls += 1;
        false
    })
    .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut polls = 0;
        assert_eq!(
            left.diffuse_binary_exponential_with(
                &right,
                geometry(),
                25.0,
                &pairs,
                controls(),
                || {
                    polls += 1;
                    polls == stop
                }
            ),
            Err(GasEnergyError::Cancelled)
        );
        assert_eq!((&left, &right), (&saved.0, &saved.1));
    }
}

#[test]
fn multicomponent_faces_keep_the_general_maxwell_stefan_owner() {
    let mixture = source::load(
        "A 0.4 0.028 3.5 0.001 1.7e-5 0.025\nB 0.3 0.004 2.5 0 2e-5 0.15\nC 0.3 0.044 4 0 1.5e-5 0.02",
        "ternary refusal fixture",
    ).unwrap();
    let point = source::point(T, P).unwrap();
    let left = ClosedGasParcel::new(mixture, point.clone(), 0.001).unwrap();
    let material = card(D, [200.0, 1800.0]);
    let pairs = [["A", "B"], ["A", "C"], ["B", "C"]].map(|names| pair(names, &material, &point));
    assert!(
        left.diffusion_flux_to(&left, geometry(), &pairs, controls(), || false)
            .is_ok()
    );
    assert!(
        left.diffuse_binary_exponential_with(&left, geometry(), 1.0, &pairs, controls(), || false)
            .is_err()
    );
}
