//! G1/G3 finite source-driven tank evolution and G4 atomic physical limits.
use super::*;
use super::super::nozzle::{
    EquilibriumNozzleControl, equilibrium_nozzle_flow,
    blowdown::{
        EquilibriumNozzleBlowdownControl, EquilibriumNozzleBlowdownError,
        equilibrium_nozzle_blowdown,
    },
};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::equilibrium::blowdown::{
    EquilibriumBlowdownError, EquilibriumBlowdownSolve, EquilibriumGasReservoir,
};
use fs_qty::Area;
use std::sync::Arc;

const AREA: f64 = 1e-5;
const VOLUME: f64 = 0.01;

fn blowdown_controls() -> EquilibriumNozzleBlowdownControl {
    let mut nozzle = EquilibriumNozzleControl::new(controls().isentropic, Pressure::new(100.0));
    nozzle.relative_sonic_tolerance = 1e-10;
    nozzle.relative_pressure_tolerance = 1e-10;
    EquilibriumNozzleBlowdownControl {
        reservoir: EquilibriumBlowdownSolve {
            minimum_pressure_pa: 100.0,
            relative_density_tolerance: 1e-12,
            mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 1e-6,
            maximum_density_iterations: 80,
        },
        nozzle,
        maximum_outgoing_fraction: 0.25,
        maximum_mass_iterations: 64,
    }
}
fn feed(chemistry: &IdealGasEquilibriumV1) -> Vec<(SpeciesId, Amount)> {
    vec![(
        id(if chemistry.models().len() == 1 {
            "A"
        } else {
            "A2"
        }),
        Amount::new(1.0),
    )]
}
fn tank(chemistry: &IdealGasEquilibriumV1) -> EquilibriumGasReservoir {
    EquilibriumGasReservoir::new(
        Arc::new(chemistry.clone()),
        feed(chemistry),
        VOLUME,
        Temperature::new(T0),
        Pressure::new(P0),
        controls().isentropic,
        || false,
    )
    .unwrap()
}
fn initial_rate(chemistry: &IdealGasEquilibriumV1, back: f64) -> f64 {
    equilibrium_nozzle_flow(
        chemistry,
        Temperature::new(T0),
        Pressure::new(P0),
        &feed(chemistry),
        Pressure::new(back),
        Area::new(AREA),
        blowdown_controls().nozzle,
        || false,
    )
    .unwrap()
    .mass_flow
    .value()
}

#[test]
fn g1_g3_frozen_limit_has_analytic_tank_cooling_conservative_energy_and_second_order_time() {
    let chemistry = single();
    let initial = tank(&chemistry);
    let gamma = 4.0_f64 / 3.0;
    let r = R / 0.01;
    let rate = AREA * P0 / (r * T0).sqrt()
        * gamma.sqrt()
        * (2.0 / (gamma + 1.0)).powf((gamma + 1.0) / (2.0 * (gamma - 1.0)));
    close(initial_rate(&chemistry, 0.0), rate, 1e-12);
    let duration = 0.5;
    let exponent = 0.5 * (gamma - 1.0);
    let exact_mass = initial.mass_kg()
        * (1.0 + exponent * rate * duration / initial.mass_kg()).powf(-1.0 / exponent);
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let mut current = initial.clone();
        let mut exhausted = 0.0;
        for _ in 0..steps {
            let before = current.mass_kg();
            let candidate = equilibrium_nozzle_blowdown(
                &current,
                Pressure::new(0.0),
                Area::new(AREA),
                duration / steps as f64,
                blowdown_controls(),
                || false,
            )
            .unwrap();
            let step = candidate.step();
            let midpoint_mass = before - 0.5 * step.gas.expelled_mass_kg;
            let expected_rate =
                rate * (midpoint_mass / initial.mass_kg()).powf(0.5 * (gamma + 1.0));
            close(step.midpoint_nozzle.mass_flow.value(), expected_rate, 1e-11);
            assert!(
                step.mass_quadrature_residual_kg.abs()
                    <= blowdown_controls().reservoir.mass_tolerance_kg
            );
            assert!(step.initial_choked && step.midpoint_nozzle.choked && step.endpoint_choked);
            assert_eq!(step.isentrope_preparations, 1);
            exhausted += step.gas.exhausted_enthalpy_j;
            current = candidate.into_parts().0;
        }
        let ratio = current.mass_kg() / initial.mass_kg();
        let temperature = T0 * ratio.powf(gamma - 1.0);
        close(
            current.isentropic_sample().temperature.value(),
            temperature,
            2e-6,
        );
        close(
            current.isentropic_sample().pressure.value(),
            P0 * ratio.powf(gamma),
            2e-4,
        );
        let expected_energy = (initial.mass_kg() * T0 - current.mass_kg() * temperature) * 3.0 * r;
        close(exhausted, expected_energy, 1e-5);
        assert_eq!(current.volume_m3(), VOLUME);
        errors.push((current.mass_kg() - exact_mass).abs());
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 3.5 && pair[0] / pair[1] < 4.5,
            "{errors:?}"
        );
    }
}

#[test]
fn g1_reacting_tank_crosses_from_choked_to_subsonic_with_one_element_and_energy_account() {
    let chemistry = dissociation();
    let original = tank(&chemistry);
    let mut current = original.clone();
    let duration = 0.015 * original.mass_kg() / initial_rate(&chemistry, 0.4 * P0);
    let mut exported_mass = 0.0;
    let mut exported_energy = 0.0;
    let mut saw_choked = false;
    let mut saw_subsonic = false;
    for _ in 0..80 {
        if current.isentropic_sample().pressure.value() <= 0.6 * P0 {
            break;
        }
        let candidate = equilibrium_nozzle_blowdown(
            &current,
            Pressure::new(0.4 * P0),
            Area::new(AREA),
            duration,
            blowdown_controls(),
            || false,
        )
        .unwrap();
        let step = candidate.step();
        saw_choked |= step.initial_choked;
        saw_subsonic |= !step.endpoint_choked;
        exported_mass += step.gas.expelled_mass_kg;
        exported_energy += step.gas.exhausted_enthalpy_j;
        check_products(&step.midpoint_nozzle.inlet, &chemistry);
        check_products(&step.midpoint_nozzle.outlet, &chemistry);
        assert!(
            step.midpoint_nozzle.stagnation_enthalpy_residual_j_kg.abs()
                <= blowdown_controls().nozzle.absolute_energy_tolerance_j_kg
        );
        assert!(candidate.endpoint().isentropic_sample().pressure.value() >= 0.4 * P0);
        current = candidate.into_parts().0;
    }
    assert!(saw_choked && saw_subsonic);
    assert!(current.isentropic_sample().pressure.value() <= 0.6 * P0);
    close(current.mass_kg() + exported_mass, original.mass_kg(), 1e-14);
    close(
        exported_energy,
        original.internal_energy_j() - current.internal_energy_j(),
        2e-5,
    );
    // The original A2 feed has two moles of elements per .002 kg. Products
    // react, so an individual molecular inventory difference is not its flux.
    let remaining_elements = current
        .thermodynamics()
        .equilibrium()
        .species()
        .iter()
        .map(|s| s.amount.value() * if s.species == id("A2") { 2.0 } else { 1.0 })
        .sum::<f64>()
        * current.physical_amount_scale();
    close(
        remaining_elements + exported_mass * 1000.0,
        original.mass_kg() * 1000.0,
        1e-9,
    );
    let a = |tank: &EquilibriumGasReservoir| {
        tank.thermodynamics()
            .equilibrium()
            .species()
            .iter()
            .find(|s| s.species == id("A"))
            .unwrap()
            .amount
            .value()
    };
    assert!((a(&current) - a(&original)).abs() > 1e-3);
    assert!(
        current.isentropic_sample().entropy_change_j_k.abs()
            <= controls().isentropic.entropy_tolerance_j_k
    );
}

#[test]
fn g4_outgoing_inventory_and_back_pressure_crossing_refuse_without_clipping() {
    let chemistry = single();
    let reservoir = tank(&chemistry);
    let before = reservoir.clone();
    let back = 0.99 * P0;
    let mass_to_equilibrium = reservoir.mass_kg() * (1.0 - 0.99_f64.powf(0.75));
    let crossing_time = 2.0 * mass_to_equilibrium / initial_rate(&chemistry, back);
    assert!(matches!(
        equilibrium_nozzle_blowdown(
            &reservoir,
            Pressure::new(back),
            Area::new(AREA),
            crossing_time,
            blowdown_controls(),
            || false
        ),
        Err(EquilibriumNozzleBlowdownError::BackPressureCrossing { .. })
    ));
    assert!(matches!(
        equilibrium_nozzle_blowdown(
            &reservoir,
            Pressure::new(0.0),
            Area::new(AREA),
            10.0 * reservoir.mass_kg() / initial_rate(&chemistry, 0.0),
            blowdown_controls(),
            || false
        ),
        Err(EquilibriumNozzleBlowdownError::InventoryLimit { .. })
    ));
    let mut changed = blowdown_controls();
    changed.nozzle.isentropic.equilibrium.maximum_iterations -= 1;
    assert!(matches!(
        equilibrium_nozzle_blowdown(
            &reservoir,
            Pressure::new(0.0),
            Area::new(AREA),
            0.01,
            changed,
            || false
        ),
        Err(EquilibriumNozzleBlowdownError::Invalid(_))
    ));
    assert_eq!(reservoir, before);
}

#[test]
fn g0_zero_drive_or_area_preserves_a_previously_depleted_tank_exactly() {
    let chemistry = dissociation();
    let initial = tank(&chemistry);
    let reservoir = initial
        .exhaust_mass_adiabatic(
            0.1 * initial.mass_kg(),
            blowdown_controls().reservoir,
            || false,
        )
        .unwrap()
        .into_parts()
        .0;
    assert!(reservoir.isentropic_sample().specific_enthalpy_change_j_kg < 0.0);
    for (back, area) in [
        (reservoir.isentropic_sample().pressure.value(), AREA),
        (0.0, 0.0),
    ] {
        let candidate = equilibrium_nozzle_blowdown(
            &reservoir,
            Pressure::new(back),
            Area::new(area),
            1.0,
            blowdown_controls(),
            || false,
        )
        .unwrap();
        assert_eq!(candidate.endpoint(), &reservoir);
        assert_eq!(candidate.step().gas.expelled_mass_kg, 0.0);
        assert_eq!(candidate.step().gas.exhausted_enthalpy_j, 0.0);
        assert_eq!(candidate.step().midpoint_nozzle.mass_flow.value(), 0.0);
        assert_eq!(
            candidate
                .step()
                .midpoint_nozzle
                .specific_enthalpy_change_j_kg,
            0.0
        );
        assert_eq!(
            candidate.step().midpoint_nozzle.inlet,
            candidate.step().midpoint_nozzle.outlet
        );
        assert_eq!(candidate.step().mass_iterations, 0);
        assert_eq!(candidate.step().density_pressure_queries, 0);
        assert_eq!(candidate.step().nozzle_pressure_queries, 0);
    }
}

#[test]
fn g4_final_finite_step_cancellation_preserves_tank_sources_and_deterministic_retry() {
    let chemistry = dissociation();
    let reservoir = tank(&chemistry);
    let before = reservoir.clone();
    let call = |cancelled: &mut dyn FnMut() -> bool| {
        equilibrium_nozzle_blowdown(
            &reservoir,
            Pressure::new(0.4 * P0),
            Area::new(AREA),
            0.003,
            blowdown_controls(),
            cancelled,
        )
    };
    let mut total = 0;
    let expected = call(&mut || {
        total += 1;
        false
    })
    .unwrap();
    let mut count = 0;
    assert!(matches!(
        call(&mut || {
            count += 1;
            count == total
        }),
        Err(EquilibriumNozzleBlowdownError::Cancelled)
    ));
    assert_eq!(count, total);
    assert_eq!(call(&mut || false).unwrap(), expected);
    assert_eq!(reservoir, before);
}

#[test]
fn g4_replacing_the_prepared_path_cannot_rebase_birth_entropy_or_loosen_policy() {
    let chemistry = single();
    let reservoir = tank(&chemistry);
    let before = reservoir.clone();
    let mut changed = controls().isentropic;
    changed.entropy_tolerance_j_k *= 2.0;
    for (temperature, policy) in [(0.9 * T0, controls().isentropic), (T0, changed)] {
        let mut prepared = reservoir.prepare_depletion(
            blowdown_controls().reservoir, || false,
        ).unwrap();
        // A mutable native path can be safely replaced in Rust. It must not
        // become the physical reservoir's new entropy reference or work policy.
        *prepared.isentrope() = EquilibriumIsentropeV1::new(
            &chemistry, Temperature::new(temperature), Pressure::new(P0),
            &feed(&chemistry), policy, || false,
        ).unwrap();
        assert!(matches!(
            prepared.sample_remaining_mass(0.9 * reservoir.mass_kg(), || false),
            Err(EquilibriumBlowdownError::Invalid(_))
        ));
        assert!(matches!(
            prepared.finish_remaining_mass(0.9 * reservoir.mass_kg(), || false),
            Err(EquilibriumBlowdownError::Invalid(_))
        ));
    }
    assert_eq!(reservoir, before);
}
