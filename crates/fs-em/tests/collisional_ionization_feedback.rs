//! Atomic populations and both thermal collision channels evolve together.
//! Initial FEEC thermal powers are explicitly held over the reacting interval;
//! the refreshed voltage solve is subsequent feedback, not a continuously
//! fixed-voltage transient. All source formation energy is counted once.

mod common;
#[path = "../../fs-material/tests/common/ionized_gas.rs"]
mod fixture;

use common::{config, ends, layered, limits};
use fs_em::conduction::{ConductionControls, PlanarConduction};
use fs_gas::chemistry::energy::ElectronHeavyCollisionSources;
use fs_gas::chemistry::ionization::{
    AtomicIonizationSolve, advance_atomic_ionization_with_collisions,
};
use fs_material::ionized_gas::{CarrierMomentumRelaxation, IonizedGasConductivity};
use fs_qty::{Power, Temperature, Time};
use fs_thermochem::coulomb::CoulombModelV1;
use fs_thermochem::electron_collision::{
    ElectronCrossSectionV1, MaxwellianRateControlsV1, MaxwellianRateWeightingV1, nist,
};
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateV1;
use fs_thermochem::kinetics::ElectronImpactRateV1;
use fs_thermochem::kinetics::atomic_ionization::{AtomicIonizationV1, AtomicStatisticalWeightsV1};
use fs_thermochem::{ReactionId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

const MOMENTUM_CROSS_SECTION_M2: f64 = 1e-22;

fn rate_controls() -> MaxwellianRateControlsV1 {
    MaxwellianRateControlsV1 {
        relative_tolerance: 1e-9,
        maximum_subdivisions: 4096,
    }
}

fn background(
    state: &FrozenIonizedGasStateV1,
    neutral: &ElectronCrossSectionV1,
) -> IonizedGasConductivity {
    let evaluated = neutral
        .evaluate(
            state.electron_temperature(),
            MaxwellianRateWeightingV1::WeakDriftMomentum,
            rate_controls(),
            || false,
        )
        .unwrap();
    let electron = CarrierMomentumRelaxation::from_electron_neutral_cross_sections(
        state,
        vec![evaluated],
        0.1,
        || false,
    )
    .unwrap();
    let ion = fixture::source_two_temperature(
        "H+",
        state.heavy_temperature().value(),
        state.electron_temperature().value(),
        state.pressure().value(),
        1e6,
    );
    IonizedGasConductivity::from_state(
        state.clone(),
        vec![electron, ion],
        fixture::limits(),
        || false,
    )
    .unwrap()
}

fn channel(state: &FrozenIonizedGasStateV1) -> AtomicIonizationV1 {
    AtomicIonizationV1::new(
        state.thermodynamics().clone(),
        ReactionId::new("coupled-H-impact-and-three-body-reverse").unwrap(),
        fixture::id("H+"),
        ElectronImpactRateV1::new(
            fixture::id("e-"),
            nist::hydrogen(fixture::id("H")),
            rate_controls(),
            0.1,
        )
        .unwrap(),
        AtomicStatisticalWeightsV1 { neutral: 2, ion: 1 },
        1e-12,
    )
    .unwrap()
}

fn field_controls() -> ConductionControls {
    ConductionControls {
        linear: config(),
        current_tolerance_a: 1e-6,
        power_tolerance_w: 1e-6,
    }
}

#[test]
fn g1_reacting_coulomb_feedback_conserves_both_pool_budgets() {
    // Both remain inside the explicit classical/weak-coupling source limits.
    // Recombination is slower here; its longer interval also exercises stiff
    // electron/heavy equilibration alongside changing charged populations.
    for (name, fraction, th, te, volume, duration, voltage, sign) in [
        ("ionization", 1e-4, 6000.0, 16000.0, 0.2, 1e-8, 2.0, 1.0),
        ("recombination", 0.002, 6000.0, 6000.0, 1.0, 1e-4, 1.0, -1.0),
    ] {
        let initial = fixture::frozen_state(fraction, th, te, volume);
        let channel = channel(&initial);
        let neutral = fixture::momentum_cross_section("H", MOMENTUM_CROSS_SECTION_M2);
        let original_coulomb = fixture::coulomb_model(&initial);
        let collisions = ElectronHeavyCollisionSources {
            momentum_sources: vec![neutral.clone()],
            coulomb: Some(original_coulomb.clone()),
            maximum_charged_mole_fraction: 0.2,
            maximum_target_thermal_speed_ratio: 0.1,
            rate_controls: rate_controls(),
        };
        let before = background(&initial, &neutral)
            .with_coulomb_drag(original_coulomb.clone(), || false)
            .unwrap();
        assert!(before.admits_electric_field([voltage, 0.0]), "{name}");
        let sigma_before = before.conductivity_s_m();
        let original = PlanarConduction::new(
            layered(0.5, 1.0, volume),
            vec![before; 4],
            limits(),
            &mut || false,
        )
        .unwrap();
        let field = original
            .solve(&ends(voltage), field_controls(), &mut || false)
            .unwrap();
        let resistance_before = original
            .resistance(&[2, 5], &[0, 3], field_controls(), &mut || false)
            .unwrap()
            .resistance_ohm();
        fixture::relative(
            field.current_density_a_m2()[0][0],
            -voltage * sigma_before,
            2e-12,
        );
        fixture::relative(
            field.power_w(),
            voltage * voltage * volume * sigma_before,
            2e-12,
        );
        let mut powers = [0.0, 0.0];
        let mut field_volume = 0.0;
        for (cell, electric_field) in field.electric_field_v_m().iter().enumerate() {
            let heat = original.conductors()[cell]
                .joule_heating_w_m3(*electric_field)
                .unwrap();
            let cell_volume = original.mesh().face_measure(cell).unwrap();
            powers[0] += cell_volume * heat.heavy_power_density_w_m3();
            powers[1] += cell_volume * heat.electron_power_density_w_m3();
            field_volume += cell_volume;
            fixture::relative(
                cell_volume * heat.total_power_density_w_m3(),
                field.cell_power_w()[cell],
                2e-12,
            );
        }
        fixture::relative(field_volume, volume, 1e-14);
        fixture::relative(powers[0] + powers[1], field.power_w(), 2e-12);
        let mut controls =
            AtomicIonizationSolve::new(Temperature::new(5500.0), Temperature::new(18000.0));
        controls.relative_tolerance = 1e-6;
        controls.energy_tolerance_j = 1e-8;
        controls.initial_step_s = 1e-10;
        controls.maximum_step_s = duration / 10.0;
        let result = advance_atomic_ionization_with_collisions(
            &initial,
            &channel,
            &collisions,
            Power::new(powers[1]),
            Power::new(powers[0]),
            Time::new(duration),
            controls,
            || false,
        )
        .unwrap_or_else(|error| panic!("{name}: {error}"));
        let next = result.state();
        let n0 = initial.amounts();
        let n1 = next.amounts();
        let extent = result.reaction_extent().value();
        assert!(sign * extent > 1e-8, "{name}: extent={extent:e}");
        fixture::relative(n1[1].value() - n0[1].value(), extent, 1e-9);
        fixture::relative(n1[2].value() - n0[2].value(), extent, 1e-9);
        assert!((n1[0].value() - n0[0].value() + extent).abs() < 1e-15);
        fixture::relative(n1[0].value() + n1[1].value(), 1.0, 1e-14);
        assert_eq!(n1[1], n1[2]);
        fixture::relative(next.mass().value(), initial.mass().value(), 1e-14);
        assert_eq!(next.volume(), initial.volume());

        let balance = result.balance();
        assert!(balance.electron_to_heavy_neutral_energy_j > 0.0, "{name}");
        assert!(balance.electron_to_heavy_coulomb_energy_j > 0.0, "{name}");
        fixture::relative(
            balance.electron_to_heavy_energy_j,
            balance.electron_to_heavy_neutral_energy_j + balance.electron_to_heavy_coulomb_energy_j,
            1e-12,
        );
        // Independent monatomic calorics retain both electron number and Te.
        // Formation energy, held electrical work, and each internal exchange
        // enter once. The NIST event threshold is not a second energy sink.
        let electron_energy = 1.5
            * R
            * (n1[2].value() * next.electron_temperature().value()
                - n0[2].value() * initial.electron_temperature().value());
        let heavy_energy =
            1.5 * R * (next.heavy_temperature().value() - initial.heavy_temperature().value());
        let chemical_storage = channel.ionization_energy_j_mol() * extent;
        assert!(
            (electron_energy + chemical_storage + balance.electron_to_heavy_energy_j
                - powers[1] * duration)
                .abs()
                < 2e-8,
            "{name}"
        );
        assert!(
            (heavy_energy - balance.electron_to_heavy_energy_j - powers[0] * duration).abs() < 2e-8,
            "{name}"
        );
        assert!(
            (electron_energy + heavy_energy + chemical_storage - field.power_w() * duration).abs()
                < 3e-8,
            "{name}"
        );
        fixture::relative(balance.chemical_storage_change_j, chemical_storage, 1e-12);
        assert!(
            balance.energy_residual_j.abs() <= controls.energy_tolerance_j,
            "{name}"
        );
        assert!(
            balance.minimum_step_entropy_excess_j_k >= -controls.entropy_tolerance_j_k,
            "{name}"
        );

        let retained = result.collision_sources().unwrap();
        assert_eq!(retained.momentum_sources, vec![neutral.clone()]);
        assert_eq!(retained.coulomb.as_ref(), Some(&original_coulomb));
        assert!(!original_coulomb.is_bound_to_populations(next));
        let after_background = background(next, &neutral);
        assert!(
            after_background
                .clone()
                .with_coulomb_drag(original_coulomb.clone(), || false)
                .is_err()
        );
        let refreshed_coulomb = CoulombModelV1::new(
            next,
            original_coulomb.screening(),
            original_coulomb.limits(),
            || false,
        )
        .unwrap();
        let after = after_background
            .with_coulomb_drag(refreshed_coulomb, || false)
            .unwrap();
        let electron = after.carriers().iter().find(|c| c.is_electron()).unwrap();
        let ion = after.carriers().iter().find(|c| !c.is_electron()).unwrap();
        fixture::relative(
            electron.number_density_per_m3(),
            n1[2].value() * fixture::NA / volume,
            1e-14,
        );
        let neutral_density = n1[0].value() * fixture::NA / volume;
        let expected_rate = (4.0 / 3.0)
            * neutral_density
            * MOMENTUM_CROSS_SECTION_M2
            * (8.0 * fixture::KB * next.electron_temperature().value()
                / (core::f64::consts::PI * fixture::ME))
                .sqrt();
        fixture::relative(electron.relaxation().rate_per_s(), expected_rate, 2e-9);
        let ae = electron.number_density_per_m3()
            * electron.mass_kg()
            * electron.relaxation().rate_per_s();
        let ai = ion.number_density_per_m3() * ion.mass_kg() * ion.relaxation().rate_per_s();
        let force = electron.number_density_per_m3() * fixture::E;
        let k = after.coulomb_pairs()[0].friction_coefficient_kg_m3_s();
        let independent_sigma = force * force * (ae + ai) / (ae * ai + k * (ae + ai));
        let sigma_after = after.conductivity_s_m();
        fixture::relative(sigma_after, independent_sigma, 2e-12);
        assert!((sigma_after / sigma_before - 1.0).abs() > 1e-8, "{name}");

        // Holding the old populations at the final two temperatures gives a
        // distinct response: the refresh actually consumes the reacting amounts.
        let old_populations = FrozenIonizedGasStateV1::new(
            initial.thermodynamics(),
            initial.amounts(),
            initial.volume(),
            next.heavy_temperature(),
            next.electron_temperature(),
            1e-12,
            || false,
        )
        .unwrap();
        let old_coulomb = CoulombModelV1::new(
            &old_populations,
            original_coulomb.screening(),
            original_coulomb.limits(),
            || false,
        )
        .unwrap();
        let unchanged_composition = background(&old_populations, &neutral)
            .with_coulomb_drag(old_coulomb, || false)
            .unwrap();
        assert!(
            (sigma_after / unchanged_composition.conductivity_s_m() - 1.0).abs() > 1e-8,
            "{name}"
        );

        let updated = PlanarConduction::new(
            layered(0.5, 1.0, volume),
            vec![after; 4],
            limits(),
            &mut || false,
        )
        .unwrap();
        let new_field = updated
            .solve(&ends(voltage), field_controls(), &mut || false)
            .unwrap();
        fixture::relative(
            new_field.current_density_a_m2()[0][0],
            -voltage * independent_sigma,
            2e-12,
        );
        fixture::relative(
            new_field.power_w(),
            voltage * voltage * volume * independent_sigma,
            2e-12,
        );
        let resistance = updated
            .resistance(&[2, 5], &[0, 3], field_controls(), &mut || false)
            .unwrap()
            .resistance_ohm();
        fixture::relative(resistance, 1.0 / (independent_sigma * volume), 2e-12);
        fixture::relative(
            resistance / resistance_before,
            sigma_before / sigma_after,
            2e-12,
        );
    }
}
