//! G0/G1/G3/G4: molecular recoil, finite inventories and inertial-frame accounts.
//! Genuine synthetic NASA/material sources; no empirical gas accuracy claim.
use super::*;
use crate::rocket_blowdown::GasRocketState;
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::{
    effusion::{
        EffusiveInterface, EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport,
        momentum::EffusiveMomentumInterface,
    },
    network::{SpeciesNetworkSolve, transport::GasTransportSolve},
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_mbd::Vec3;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
#[path = "../../../../../fs-material/examples/material_gas_spring/source.rs"]
mod source;

const PURE: &str = "A 1 0.028 2.5 0 1.7e-5 0.025";
const MIXTURE: &str = "A 0.5 0.028 2.5 0 1.7e-5 0.025\nB 0.5 0.004 2.5 0 2e-5 0.15";
const T0: f64 = 400.0;
const P0: f64 = 100_000.0;
const VOLUME: f64 = 0.001;
const AREA: f64 = 1e-7;
const DT: f64 = 0.1;
const FRAME: u64 = 7;

fn initial(text: &str, velocity: Vec3) -> GasRocketState {
    let model = source::load(text, "effusive recoil analytic fixture").unwrap();
    let gas = ClosedGasParcel::new(model, source::point(T0, P0).unwrap(), VOLUME).unwrap();
    GasRocketState::new(gas, 0.01, Vec3::ZERO, velocity, FRAME).unwrap()
}

fn interface(model: &ThermochemicalGasModel) -> EffusiveMomentumInterface {
    let rows = model
        .thermodynamics()
        .components()
        .iter()
        .map(|species| {
            EffusiveSpeciesTransmission::try_new(
                species.species().clone(),
                1.0,
                "synthetic reciprocal transmission independent of energy/internal state",
            )
            .unwrap()
        })
        .collect();
    EffusiveMomentumInterface::try_new(
        EffusiveInterface::try_new(
            model.clone(),
            rows,
            EffusiveSupport {
                temperature_k: (200.0, 1800.0),
                pressure_pa: (1000.0, 1e7),
            },
            EffusiveRegime::try_new(
                1e-9,
                0.0,
                1e-7,
                "synthetic collisionless thin-array qualification, not measured gas validation",
            )
            .unwrap(),
        )
        .unwrap(),
        "synthetic uncollimated angle-independent nondeflecting Maxwellian aperture",
    )
    .unwrap()
}

fn outlet(interface: &EffusiveMomentumInterface, normal: Vec3) -> EffusionOutlet<'_> {
    EffusionOutlet {
        id: 4,
        interface,
        area_m2: AREA,
        outward_normal_world: normal,
    }
}

fn endpoint() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}

fn controls() -> EffusionRocketSolve {
    EffusionRocketSolve {
        transport: GasTransportSolve {
            exchange: SpeciesNetworkSolve {
                maximum_species_fraction: 0.1,
                mass_tolerance_kg: 1e-13,
                energy_tolerance_j: 1e-8,
            },
            entropy_tolerance_j_k: 1e-10,
        },
        momentum_tolerance_n_s: 1e-12,
    }
}

fn near(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "{actual:.17e} != {expected:.17e} (absolute allowance {absolute:.3e})"
    );
}

fn vector_near(actual: Vec3, expected: Vec3, absolute: f64) {
    for (a, b) in [
        (actual.x, expected.x),
        (actual.y, expected.y),
        (actual.z, expected.z),
    ] {
        near(a, b, absolute);
    }
}

fn mass(state: &GasRocketState) -> f64 {
    state.vehicle().dry_mass_kg + state.reservoir().mass_kg()
}

fn kinetic(state: &GasRocketState) -> f64 {
    0.5 * mass(state) * state.vehicle().velocity_m_s.norm_squared()
}

fn internal(state: &GasRocketState) -> f64 {
    state.reservoir().mass_kg() * state.reservoir().state().specific_internal_energy_j_kg()
}

fn escaped_moles(molar_mass: f64, fraction: f64) -> f64 {
    DT * AREA * fraction * P0 / (2.0 * core::f64::consts::PI * molar_mass * R * T0).sqrt()
}

fn check_accounts(before: &GasRocketState, candidate: &EffusionRocketCandidate) {
    let after = candidate.endpoint();
    let step = candidate.step();
    let body_momentum_change = after
        .vehicle()
        .velocity_m_s
        .scale(mass(after))
        .sub(before.vehicle().velocity_m_s.scale(mass(before)));
    vector_near(
        body_momentum_change
            .add(step.exhaust_momentum_world_kg_m_s)
            .sub(step.exhaust.external_impulse_n_s),
        Vec3::ZERO,
        2e-12,
    );
    near(
        internal(after) - internal(before) + kinetic(after) - kinetic(before)
            + step.exhaust_energy_world_j
            - step.exhaust.external_work_j,
        0.0,
        2e-8,
    );
    vector_near(step.momentum_residual_kg_m_s, Vec3::ZERO, 2e-12);
    near(
        step.energy_residual_j,
        0.0,
        controls().transport.exchange.energy_tolerance_j,
    );
    near(
        step.remaining_mass_residual_kg,
        0.0,
        controls().transport.exchange.mass_tolerance_kg,
    );
    near(
        step.consumed_mass_residual_kg,
        0.0,
        controls().transport.exchange.mass_tolerance_kg,
    );
    near(
        after.vehicle().propellant_mass_kg,
        after.reservoir().mass_kg(),
        1e-13,
    );
    assert_eq!(after.vehicle().dry_mass_kg, before.vehicle().dry_mass_kg);
    assert_eq!(after.reservoir().volume_m3(), VOLUME);
}

#[test]
fn monatomic_effusion_accelerates_the_real_vehicle_with_half_pressure_recoil_once() {
    let before = initial(PURE, Vec3::ZERO);
    let saved = before.clone();
    let law = interface(before.reservoir().model());
    let normal = Vec3::new(0.6, 0.8, 0.0);
    let candidate = before
        .advance_effusion(
            FRAME,
            &[outlet(&law, normal)],
            DT,
            endpoint(),
            controls(),
            || false,
        )
        .unwrap();
    let escaped = escaped_moles(0.028, 1.0);
    let dm = escaped * 0.028;
    let impulse = AREA * P0 * DT / 2.0;
    let speed = (core::f64::consts::PI * R * T0 / (2.0 * 0.028)).sqrt();
    near(candidate.step().expelled_mass_kg, dm, 1e-18);
    vector_near(
        candidate.step().relative_momentum_kg_m_s,
        normal.scale(impulse),
        2e-16,
    );
    vector_near(
        candidate.step().exhaust.mechanics.thrust_impulse_n_s,
        normal.scale(-impulse),
        2e-16,
    );
    let dv = speed * (mass(&before) / mass(candidate.endpoint())).ln();
    vector_near(
        candidate.endpoint().vehicle().velocity_m_s,
        normal.scale(-dv),
        1e-11,
    );
    assert!(candidate.endpoint().vehicle().velocity_m_s.dot(normal) < 0.0);
    // The incident half-space moment is pA/2. Adding a separate pA term would
    // triple this recoil and double-count momentum already in the molecules.
    near(
        candidate.step().relative_momentum_kg_m_s.dot(normal) / (AREA * P0 * DT),
        0.5,
        2e-13,
    );
    let n0 = P0 * VOLUME / (R * T0);
    let expected_t = T0 * (n0 - 4.0 * escaped / 3.0) / (n0 - escaped);
    near(
        candidate.endpoint().reservoir().state().state().temperature,
        expected_t,
        2e-8,
    );
    near(
        candidate.step().relative_energy_j,
        2.0 * escaped * R * T0,
        2e-10,
    );
    check_accounts(&before, &candidate);
    assert_eq!(before, saved);
}

#[test]
fn mixture_recoil_uses_actual_species_losses_and_preserves_fractionation() {
    let before = initial(MIXTURE, Vec3::ZERO);
    let law = interface(before.reservoir().model());
    let candidate = before
        .advance_effusion(
            FRAME,
            &[outlet(&law, Vec3::new(1.0, 0.0, 0.0))],
            DT,
            endpoint(),
            controls(),
            || false,
        )
        .unwrap();
    let amounts = [escaped_moles(0.028, 0.5), escaped_moles(0.004, 0.5)];
    near(amounts[1] / amounts[0], 7.0_f64.sqrt(), 1e-14);
    let boundaries = &candidate.step().gas.exchange.boundaries;
    assert_eq!(boundaries.len(), 2);
    for (row, (id, molar_mass, amount)) in boundaries
        .iter()
        .zip([("A", 0.028, amounts[0]), ("B", 0.004, amounts[1])])
    {
        assert_eq!(row.transfer.species.as_str(), id);
        near(row.transfer.outgoing_mass_kg, amount * molar_mass, 1e-18);
        assert_eq!(row.transfer.incoming_mass_kg, 0.0);
    }
    let dn = amounts[0] + amounts[1];
    let n0 = P0 * VOLUME / (R * T0);
    let remaining = candidate
        .endpoint()
        .reservoir()
        .model()
        .thermodynamics()
        .mole_composition();
    near(
        remaining.fractions()[1],
        (0.5 * n0 - amounts[1]) / (n0 - dn),
        2e-13,
    );
    assert!(
        remaining.fractions()[1] < 0.5,
        "lighter molecules leave more quickly"
    );
    let dm = amounts[0] * 0.028 + amounts[1] * 0.004;
    near(candidate.step().expelled_mass_kg, dm, 1e-18);
    vector_near(
        candidate.step().relative_momentum_kg_m_s,
        Vec3::new(AREA * P0 * DT / 2.0, 0.0, 0.0),
        2e-16,
    );
    near(
        candidate.endpoint().reservoir().state().state().temperature,
        T0 * (n0 - 4.0 * dn / 3.0) / (n0 - dn),
        2e-8,
    );
    near(candidate.step().relative_energy_j, 2.0 * dn * R * T0, 2e-10);
    check_accounts(&before, &candidate);
}

#[test]
fn opposed_apertures_cancel_thrust_while_both_beams_carry_mass_and_energy() {
    let before = initial(MIXTURE, Vec3::ZERO);
    let law = interface(before.reservoir().model());
    let first = outlet(&law, Vec3::new(1.0, 0.0, 0.0));
    let second = EffusionOutlet {
        id: 9,
        outward_normal_world: Vec3::new(-1.0, 0.0, 0.0),
        ..first
    };
    let candidate = before
        .advance_effusion(FRAME, &[first, second], DT, endpoint(), controls(), || {
            false
        })
        .unwrap();
    let reverse = before
        .advance_effusion(FRAME, &[second, first], DT, endpoint(), controls(), || {
            false
        })
        .unwrap();
    assert_eq!(candidate, reverse);
    assert_eq!(candidate.step().relative_momentum_kg_m_s, Vec3::ZERO);
    assert_eq!(
        candidate.step().exhaust.mechanics.thrust_impulse_n_s,
        Vec3::ZERO
    );
    assert_eq!(candidate.endpoint().vehicle().velocity_m_s, Vec3::ZERO);
    assert_eq!(candidate.endpoint().vehicle().position_m, Vec3::ZERO);
    let dm = 2.0 * (0.028 * escaped_moles(0.028, 0.5) + 0.004 * escaped_moles(0.004, 0.5));
    near(candidate.step().expelled_mass_kg, dm, 2e-18);
    assert!(candidate.step().exhaust_energy_world_j > 0.0);
    assert!(candidate.endpoint().reservoir().state().state().temperature < T0);
    assert_eq!(candidate.step().gas.exchange.boundaries.len(), 4);
    check_accounts(&before, &candidate);
}

#[test]
fn exported_world_momentum_and_energy_transform_with_the_same_galilean_boost() {
    let velocity = Vec3::new(3.0, -2.0, 1.0);
    let boost = Vec3::new(30.0, -40.0, 10.0);
    let before = initial(MIXTURE, velocity);
    let boosted = initial(MIXTURE, velocity.add(boost));
    let law = interface(before.reservoir().model());
    let paths = [outlet(&law, Vec3::new(0.6, 0.8, 0.0))];
    let a = before
        .advance_effusion(FRAME, &paths, DT, endpoint(), controls(), || false)
        .unwrap();
    let b = boosted
        .advance_effusion(FRAME, &paths, DT, endpoint(), controls(), || false)
        .unwrap();
    assert_eq!(a.endpoint().reservoir(), b.endpoint().reservoir());
    assert_eq!(
        a.step().relative_momentum_kg_m_s,
        b.step().relative_momentum_kg_m_s
    );
    vector_near(
        b.endpoint()
            .vehicle()
            .velocity_m_s
            .sub(a.endpoint().vehicle().velocity_m_s),
        boost,
        1e-13,
    );
    vector_near(
        b.endpoint()
            .vehicle()
            .position_m
            .sub(a.endpoint().vehicle().position_m),
        boost.scale(DT),
        1e-13,
    );
    let exported_a = a.step();
    let exported_b = b.step();
    let dm = a.step().expelled_mass_kg;
    vector_near(
        exported_b.exhaust_momentum_world_kg_m_s,
        exported_a
            .exhaust_momentum_world_kg_m_s
            .add(boost.scale(dm)),
        2e-15,
    );
    near(
        exported_b.exhaust_energy_world_j,
        exported_a.exhaust_energy_world_j
            + boost.dot(exported_a.exhaust_momentum_world_kg_m_s)
            + 0.5 * dm * boost.norm_squared(),
        2e-12,
    );
    check_accounts(&before, &a);
    check_accounts(&boosted, &b);
}

#[test]
fn closed_coast_source_refusal_and_late_cancellation_publish_neither_owner() {
    let before = initial(PURE, Vec3::new(3.0, 2.0, -1.0));
    let original = before.clone();
    let law = interface(before.reservoir().model());
    let path = outlet(&law, Vec3::new(1.0, 0.0, 0.0));
    let closed = EffusionOutlet {
        area_m2: 0.0,
        ..path
    };
    for paths in [&[][..], &[closed][..]] {
        let coast = before
            .advance_effusion(FRAME, paths, DT, endpoint(), controls(), || false)
            .unwrap();
        assert_eq!(coast.endpoint().reservoir(), before.reservoir());
        assert_eq!(
            coast.endpoint().vehicle().velocity_m_s,
            before.vehicle().velocity_m_s
        );
        vector_near(
            coast.endpoint().vehicle().position_m,
            before.vehicle().velocity_m_s.scale(DT),
            1e-15,
        );
        assert_eq!(coast.step().expelled_mass_kg, 0.0);
        assert_eq!(coast.step().exhaust_energy_world_j, 0.0);
        assert_eq!(coast.step().exhaust_momentum_world_kg_m_s, Vec3::ZERO);
        check_accounts(&before, &coast);
    }
    assert!(
        before
            .advance_effusion(FRAME + 1, &[path], DT, endpoint(), controls(), || false)
            .is_err()
    );
    assert!(
        before
            .advance_effusion(FRAME, &[path, path], DT, endpoint(), controls(), || false)
            .is_err()
    );
    assert!(
        before
            .advance_effusion(
                FRAME,
                &[EffusionOutlet {
                    outward_normal_world: Vec3::ZERO,
                    ..path
                }],
                DT,
                endpoint(),
                controls(),
                || false
            )
            .is_err()
    );
    assert!(
        before
            .advance_effusion(FRAME, &[path], 0.0, endpoint(), controls(), || false)
            .is_err()
    );
    let changed = initial("A 1 0.028 3.5 0 1.7e-5 0.025", Vec3::ZERO);
    let foreign = interface(changed.reservoir().model());
    assert!(
        before
            .advance_effusion(
                FRAME,
                &[EffusionOutlet {
                    interface: &foreign,
                    ..path
                }],
                DT,
                endpoint(),
                controls(),
                || false
            )
            .is_err()
    );
    let narrow = GasEnergySolve {
        minimum_temperature_k: T0 - 1e-6,
        ..endpoint()
    };
    assert!(
        before
            .advance_effusion(FRAME, &[path], DT, narrow, controls(), || false)
            .is_err()
    );
    assert!(
        before
            .advance_effusion(
                FRAME,
                &[EffusionOutlet {
                    area_m2: 1.0,
                    ..path
                }],
                DT,
                endpoint(),
                controls(),
                || false
            )
            .is_err()
    );
    // This boost makes a real recoil smaller than the representable velocity
    // spacing. A conservative candidate must refuse the lost momentum.
    let huge_boost = initial(PURE, Vec3::new(1e16, 0.0, 0.0));
    let saved_boost = huge_boost.clone();
    assert!(matches!(
        huge_boost.advance_effusion(FRAME, &[path], DT, endpoint(), controls(), || false),
        Err(EffusionRocketError::MomentumBalance { .. })
    ));
    assert_eq!(huge_boost, saved_boost);
    let mut polls = 0;
    before
        .advance_effusion(FRAME, &[path], DT, endpoint(), controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls - 1, polls] {
        let mut count = 0;
        assert!(
            before
                .advance_effusion(FRAME, &[path], DT, endpoint(), controls(), || {
                    count += 1;
                    count == stop
                })
                .is_err()
        );
        assert_eq!(before, original);
    }
}

fn loaded_environment() -> EffusionRocketEnvironment {
    EffusionRocketEnvironment {
        inertial_frame: FRAME,
        gravity_world_m_s2: Vec3::new(0.0, -9.81, 0.0),
        external_force_world_n: Vec3::new(0.002, 0.003, -0.001),
    }
}

#[test]
fn explicit_zero_loads_retain_the_original_free_entry_exactly() {
    let before = initial(MIXTURE, Vec3::new(2.0, -1.0, 0.5));
    let law = interface(before.reservoir().model());
    let outlets = [outlet(&law, Vec3::new(0.6, 0.8, 0.0))];
    let original = before
        .advance_effusion(FRAME, &outlets, DT, endpoint(), controls(), || false)
        .unwrap();
    let explicit = before
        .advance_effusion_with_environment(
            EffusionRocketEnvironment {
                inertial_frame: FRAME,
                gravity_world_m_s2: Vec3::ZERO,
                external_force_world_n: Vec3::ZERO,
            },
            &outlets,
            DT,
            endpoint(),
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(original, explicit);
    assert_eq!(explicit.step().exhaust.external_impulse_n_s, Vec3::ZERO);
    assert_eq!(explicit.step().exhaust.external_work_j, 0.0);
    assert_eq!(explicit.step().gravity_mass_history_uncertainty_kg, 0.0);
}

#[test]
fn loaded_escaping_gas_uses_changing_mass_and_the_same_external_work_under_a_boost() {
    let before = initial(MIXTURE, Vec3::new(2.0, -1.0, 0.5));
    let law = interface(before.reservoir().model());
    let outlets = [outlet(&law, Vec3::new(0.6, 0.8, 0.0))];
    let environment = loaded_environment();
    let candidate = before
        .advance_effusion_with_environment(
            environment,
            &outlets,
            DT,
            endpoint(),
            controls(),
            || false,
        )
        .unwrap();
    let free = before
        .advance_effusion(FRAME, &outlets, DT, endpoint(), controls(), || false)
        .unwrap();
    assert_eq!(
        candidate.endpoint().reservoir(),
        free.endpoint().reservoir()
    );
    let step = candidate.step();
    let dm = step.expelled_mass_kg;
    let m0 = mass(&before);
    let fraction = dm / m0;
    let logarithm = -(-fraction).ln_1p();
    let position_factor = (fraction - (1.0 - fraction) * logarithm) / (fraction * fraction);
    let thrustward_speed = step.relative_momentum_kg_m_s.scale(-1.0 / dm);
    let expected_velocity = before
        .vehicle()
        .velocity_m_s
        .add(thrustward_speed.scale(logarithm))
        .add(
            environment
                .external_force_world_n
                .scale(DT * logarithm / dm),
        )
        .add(environment.gravity_world_m_s2.scale(DT));
    let expected_position = before
        .vehicle()
        .position_m
        .add(before.vehicle().velocity_m_s.scale(DT))
        .add(thrustward_speed.scale(DT * fraction * position_factor))
        .add(
            environment
                .external_force_world_n
                .scale(DT * DT * position_factor / m0),
        )
        .add(environment.gravity_world_m_s2.scale(0.5 * DT * DT));
    vector_near(
        candidate.endpoint().vehicle().velocity_m_s,
        expected_velocity,
        2e-12,
    );
    vector_near(
        candidate.endpoint().vehicle().position_m,
        expected_position,
        2e-10,
    );
    // Gravity follows the linearly decreasing held-rate mass, unlike the fixed force.
    vector_near(
        step.exhaust.external_impulse_n_s,
        environment
            .gravity_world_m_s2
            .scale((m0 - 0.5 * dm) * DT)
            .add(environment.external_force_world_n.scale(DT)),
        2e-14,
    );
    assert!(step.exhaust.external_work_j.abs() > 1e-5);
    assert!(step.gravity_impulse_uncertainty_n_s.y >= 0.0);
    assert!(
        step.energy_residual_j.abs() + step.gravity_work_uncertainty_j
            <= controls().transport.exchange.energy_tolerance_j
    );
    check_accounts(&before, &candidate);

    let boost = Vec3::new(30.0, -40.0, 10.0);
    let shifted = initial(MIXTURE, before.vehicle().velocity_m_s.add(boost));
    let shifted_candidate = shifted
        .advance_effusion_with_environment(
            environment,
            &outlets,
            DT,
            endpoint(),
            controls(),
            || false,
        )
        .unwrap();
    let shifted_step = shifted_candidate.step();
    vector_near(
        shifted_candidate
            .endpoint()
            .vehicle()
            .velocity_m_s
            .sub(candidate.endpoint().vehicle().velocity_m_s),
        boost,
        1e-13,
    );
    vector_near(
        shifted_step.exhaust.external_impulse_n_s,
        step.exhaust.external_impulse_n_s,
        2e-14,
    );
    near(
        shifted_step.exhaust.external_work_j,
        step.exhaust.external_work_j + boost.dot(step.exhaust.external_impulse_n_s),
        2e-12,
    );
    vector_near(
        shifted_step.exhaust_momentum_world_kg_m_s,
        step.exhaust_momentum_world_kg_m_s.add(boost.scale(dm)),
        2e-15,
    );
    near(
        shifted_step.exhaust_energy_world_j,
        step.exhaust_energy_world_j
            + boost.dot(step.exhaust_momentum_world_kg_m_s)
            + 0.5 * dm * boost.norm_squared(),
        2e-12,
    );
    check_accounts(&shifted, &shifted_candidate);
}

#[test]
fn closed_loaded_vessel_follows_constant_mass_ballistics_without_a_gas_query() {
    let before = initial(PURE, Vec3::new(3.0, 2.0, -1.0));
    let environment = loaded_environment();
    let foreign = initial("A 1 0.028 3.5 0 1.7e-5 0.025", Vec3::ZERO);
    let law = interface(foreign.reservoir().model());
    let closed = EffusionOutlet {
        area_m2: 0.0,
        ..outlet(&law, Vec3::new(1.0, 0.0, 0.0))
    };
    for outlets in [&[][..], &[closed][..]] {
        let candidate = before
            .advance_effusion_with_environment(
                environment,
                outlets,
                DT,
                endpoint(),
                controls(),
                || false,
            )
            .unwrap();
        let force = environment
            .external_force_world_n
            .add(environment.gravity_world_m_s2.scale(mass(&before)));
        let acceleration = force.scale(1.0 / mass(&before));
        let displacement = before
            .vehicle()
            .velocity_m_s
            .scale(DT)
            .add(acceleration.scale(0.5 * DT * DT));
        vector_near(
            candidate.endpoint().vehicle().velocity_m_s,
            before.vehicle().velocity_m_s.add(acceleration.scale(DT)),
            2e-14,
        );
        vector_near(
            candidate.endpoint().vehicle().position_m,
            before.vehicle().position_m.add(displacement),
            2e-14,
        );
        vector_near(
            candidate.step().exhaust.external_impulse_n_s,
            force.scale(DT),
            2e-14,
        );
        near(
            candidate.step().exhaust.external_work_j,
            force.dot(displacement),
            2e-14,
        );
        assert_eq!(candidate.endpoint().reservoir(), before.reservoir());
        assert_eq!(candidate.step().expelled_mass_kg, 0.0);
        assert_eq!(candidate.step().exhaust_momentum_world_kg_m_s, Vec3::ZERO);
        assert_eq!(candidate.step().exhaust_energy_world_j, 0.0);
        assert_eq!(candidate.step().gravity_mass_history_uncertainty_kg, 0.0);
        check_accounts(&before, &candidate);
    }
}

#[test]
fn initial_inventory_rounding_spends_physical_budget_before_canonicalization() {
    let mut before = initial(PURE, Vec3::new(1e6, 0.0, 0.0));
    // Manufacture the admitted small inventory mismatch a previous native owner
    // may leave. It is small in kg but observable as momentum/work under a boost.
    before.vehicle.propellant_mass_kg += 1e-15;
    let saved = before.clone();
    let quiet = EffusionRocketEnvironment {
        inertial_frame: FRAME,
        gravity_world_m_s2: Vec3::ZERO,
        external_force_world_n: Vec3::ZERO,
    };
    assert!(matches!(
        before.advance_effusion_with_environment(quiet, &[], DT, endpoint(), controls(), || false),
        Err(EffusionRocketError::RoundingBudget { .. })
    ));
    let mut admitted = controls();
    admitted.momentum_tolerance_n_s = 1e-8;
    admitted.transport.exchange.energy_tolerance_j = 1e-2;
    let candidate = before
        .advance_effusion_with_environment(quiet, &[], DT, endpoint(), admitted, || false)
        .unwrap();
    assert!(candidate.step().initial_mass_residual_kg > 0.0);
    assert!(candidate.step().initial_inventory_momentum_correction_n_s.x < 0.0);
    assert!(candidate.step().initial_inventory_kinetic_correction_j < 0.0);
    assert_eq!(
        candidate.endpoint().vehicle().propellant_mass_kg,
        candidate.endpoint().reservoir().mass_kg()
    );
    assert_eq!(before, saved);
}

#[test]
fn loaded_refusals_and_late_cancellation_leave_both_original_owners_unchanged() {
    let before = initial(MIXTURE, Vec3::new(2.0, -1.0, 0.5));
    let saved = before.clone();
    let law = interface(before.reservoir().model());
    let outlets = [outlet(&law, Vec3::new(1.0, 0.0, 0.0))];
    let environment = loaded_environment();
    for invalid in [
        EffusionRocketEnvironment {
            inertial_frame: FRAME + 1,
            ..environment
        },
        EffusionRocketEnvironment {
            gravity_world_m_s2: Vec3::new(f64::NAN, 0.0, 0.0),
            ..environment
        },
        EffusionRocketEnvironment {
            external_force_world_n: Vec3::new(0.0, f64::INFINITY, 0.0),
            ..environment
        },
        EffusionRocketEnvironment {
            external_force_world_n: Vec3::new(1e308, 0.0, 0.0),
            ..environment
        },
    ] {
        assert!(
            before
                .advance_effusion_with_environment(
                    invalid,
                    &outlets,
                    DT,
                    endpoint(),
                    controls(),
                    || false
                )
                .is_err()
        );
        assert_eq!(before, saved);
    }
    let narrow = GasEnergySolve {
        minimum_temperature_k: T0 - 1e-6,
        ..endpoint()
    };
    assert!(
        before
            .advance_effusion_with_environment(
                environment,
                &outlets,
                DT,
                narrow,
                controls(),
                || false
            )
            .is_err()
    );
    let mut polls = 0;
    before
        .advance_effusion_with_environment(
            environment,
            &outlets,
            DT,
            endpoint(),
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, polls - 1, polls] {
        let mut count = 0;
        assert!(
            before
                .advance_effusion_with_environment(
                    environment,
                    &outlets,
                    DT,
                    endpoint(),
                    controls(),
                    || {
                        count += 1;
                        count == stop
                    }
                )
                .is_err()
        );
        assert_eq!(before, saved);
    }
}

#[test]
fn gravity_reserve_keeps_representable_products_across_extreme_scales() {
    // epsilon*dt underflows, but multiplication with finite gravity and RMS
    // speed yields real, representable impulse and work reserves.
    let (impulse, work) =
        gravity_rounding_bounds(1e-250, Vec3::new(1e200, 0.0, 0.0), 1e-100, 1e200).unwrap();
    near(impulse.x / 1e-150, 1.0, 2e-15);
    assert_eq!(impulse.y, 0.0);
    assert_eq!(impulse.z, 0.0);
    near(work / 1e-50, 1.0, 2e-15);
    // A mathematically positive reserve must refuse when its FINAL value
    // cannot be represented; it must never silently become an exact zero.
    assert!(gravity_rounding_bounds(1e-300, Vec3::new(1e-100, 0.0, 0.0), 1e-100, 1.0,).is_err());
    assert!(gravity_rounding_bounds(1e-100, Vec3::new(1e-100, 0.0, 0.0), 1.0, 1e-300,).is_err());

    let before = initial(PURE, Vec3::ZERO);
    let saved = before.clone();
    let law = interface(before.reservoir().model());
    let candidate = before
        .advance_effusion_with_environment(
            loaded_environment(),
            &[outlet(&law, Vec3::new(1.0, 0.0, 0.0))],
            DT,
            endpoint(),
            controls(),
            || false,
        )
        .unwrap();
    assert!(candidate.step().gravity_work_uncertainty_j.is_finite());
    assert!(candidate.step().gravity_work_uncertainty_j >= 0.0);
    assert_eq!(before, saved);
}

mod trajectory;

mod exchange;

mod thermal;

mod channel;
