//! Rigid open-inventory energy/entropy through original phase calorics.
//! References below come from an independent PR cubic/fugacity solve plus the
//! species U,V balance, not this intake wrapper or a constant-Cp dense surrogate.
#[path = "common/nonideal_flash.rs"]
mod common;

use common::{METHANE, model, source};
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
    flash::{
        ThermochemicalFlashError as NativeError,
        equilibrium::{PhaseEnthalpyControls, volume::PhaseEnergyVolumeControls},
        parcel::{
            IsochoricMixtureEquilibriumParcel,
            transfer::{
                MixtureEquilibriumInflow, MixtureEquilibriumInflowControl,
                MixtureEquilibriumInflowError as Error, MixtureEquilibriumInflowStep,
            },
        },
    },
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::{
    Composition, CompositionBasis, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}; difference {}, allowance {tolerance}",
        actual - expected
    );
}
fn composition(
    base: &ThermochemicalPengRobinsonModel,
    methane: f64,
) -> ThermochemicalPengRobinsonModel {
    let source = base
        .source()
        .with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![1.0 - methane, methane]).unwrap(),
        )
        .unwrap();
    let bindings = [
        ("methane", 190.564, 4.5992e6, 0.01142),
        ("ethane", 305.322, 4.8722e6, 0.0995),
    ]
    .into_iter()
    .filter_map(|(name, tc, pc, omega)| {
        let species = SpeciesId::new(name).unwrap();
        source
            .thermodynamics()
            .components()
            .iter()
            .any(|s| s.species() == &species)
            .then_some(PengRobinsonSpeciesBinding {
                species,
                critical_temperature: Temperature::new(tc),
                critical_pressure: Pressure::new(pc),
                acentric_factor: omega,
            })
    })
    .collect();
    ThermochemicalPengRobinsonModel::new(source, bindings, vec![]).unwrap()
}
fn controls(p: [f64; 2], t: [f64; 2]) -> MixtureEquilibriumInflowControl {
    MixtureEquilibriumInflowControl::new(PhaseEnergyVolumeControls::new(
        p.map(Pressure::new),
        PhaseEnthalpyControls::new(t.map(Temperature::new)),
    ))
}
fn dense_controls() -> MixtureEquilibriumInflowControl {
    controls([3.1e6, 4.3e6], [220.0, 255.0])
}
fn tank(
    base: ThermochemicalPengRobinsonModel,
    t: f64,
    p: f64,
) -> IsochoricMixtureEquilibriumParcel {
    IsochoricMixtureEquilibriumParcel::from_pressure_temperature(
        base,
        1.0,
        Temperature::new(t),
        Pressure::new(p),
        dense_controls().equilibrium.phase.phase_options(),
        || false,
    )
    .unwrap()
}
fn portion(
    model: &ThermochemicalPengRobinsonModel,
    t: f64,
    p: f64,
    mass: f64,
) -> MixtureEquilibriumInflow<'_> {
    MixtureEquilibriumInflow {
        model,
        temperature: Temperature::new(t),
        pressure: Pressure::new(p),
        mass_kg: mass,
    }
}
fn audit(
    tank: &IsochoricMixtureEquilibriumParcel,
    step: &MixtureEquilibriumInflowStep,
    c: MixtureEquilibriumInflowControl,
) {
    assert!(
        step.energy_residual_j.abs()
            <= tank.mass_kg() * c.equilibrium.internal_energy_tolerance_j_kg
    );
    assert!(
        step.volume_residual_m3.abs()
            <= tank.mass_kg() * c.equilibrium.specific_volume_tolerance_m3_kg
    );
    assert!(step.entropy_production_j_k >= -c.entropy_tolerance_j_k);
    assert!(
        step.mass_residual_kg.abs()
            <= c.absolute_mass_tolerance_kg + c.relative_mass_tolerance * tank.mass_kg()
    );
    for s in &step.species {
        assert!(
            s.residual_kg.abs()
                <= c.absolute_mass_tolerance_kg
                    + c.relative_mass_tolerance
                        * s.final_mass_kg.max(s.initial_mass_kg + s.incoming_mass_kg)
        );
    }
    assert!(step.pressure_evaluations <= c.equilibrium.maximum_pressure_evaluations);
    assert!(step.pt_evaluations <= c.equilibrium.maximum_pt_evaluations);
    assert!(
        step.eos_evaluation_upper_bound <= c.equilibrium.maximum_eos_evaluation_bound().unwrap()
    );
    assert_eq!(tank.state().constitutive_law(), tank.model().mechanics());
}

#[test]
fn distinct_inlet_enthalpy_fills_rigid_two_phase_inventory_with_independent_uv_oracle() {
    let base = source();
    let feed = composition(&base, 0.8);
    let mut vessel = tank(base, 230.0, 3e6);
    let before = vessel.clone();
    let c = dense_controls();
    let result = vessel
        .receive_stream(portion(&feed, 280.0, 8e6, 0.1), None, c, || false)
        .unwrap();
    close(vessel.mass_kg(), 1.1, 0.0);
    assert_eq!(vessel.volume_m3().to_bits(), before.volume_m3().to_bits());
    close(vessel.volume_m3(), 0.011_114_793_676_556_982, 2e-12);
    close(
        vessel
            .model()
            .source()
            .thermodynamics()
            .mole_composition()
            .fractions()[METHANE],
        0.532_698_010_825_802_5,
        3e-15,
    );
    close(vessel.state().temperature_k(), 235.373_785_330_774_07, 3e-6);
    close(vessel.state().pressure_pa(), 3_606_532.867_079_617, 0.05);
    close(vessel.state().density_kg_m3(), 98.967_199_213_071_21, 4e-7);
    close(
        vessel.state().vapor_mass_fraction().unwrap(),
        0.460_805_431_027_524_1,
        2e-8,
    );
    close(
        vessel
            .state()
            .phase()
            .two_phase()
            .unwrap()
            .vapor_fraction_mol(),
        0.521_712_091_624_367_9,
        2e-8,
    );
    close(result.entropy_production_j_k, 29.948_027_310_295_856, 2e-5);
    close(result.incoming_enthalpy_j, 32_760.374_503_398_75, 2e-5);
    close(
        result.incoming_pressure_flow_work_j,
        8_825.272_854_210_845,
        2e-5,
    );
    close(
        result.stored_internal_energy_change_j,
        vessel.internal_energy_j() - before.internal_energy_j(),
        2e-8,
    );
    close(
        result.stored_internal_energy_change_j - result.incoming_enthalpy_j,
        result.energy_residual_j,
        2e-8,
    );
    close(
        result.stored_entropy_change_j_k,
        vessel.mass_kg() * vessel.state().specific_entropy_j_kg_k()
            - before.mass_kg() * before.state().specific_entropy_j_kg_k(),
        2e-8,
    );
    close(
        result.stored_entropy_change_j_k - result.incoming_entropy_j_k,
        result.entropy_production_j_k,
        2e-8,
    );
    // Omitting pressure-flow work would lose more than8.8kJ in this small fill.
    let inlet = result.inlet_state.as_ref().unwrap();
    close(
        result.incoming_enthalpy_j - 0.1 * inlet.specific_internal_energy_j_kg(),
        result.incoming_pressure_flow_work_j,
        2e-8,
    );
    assert!(result.ideal_mixing_entropy_j_k > 0.0);
    assert_eq!(result.species.len(), 2);
    assert!(result.pt_evaluations >= result.pressure_evaluations + 2);
    audit(&vessel, &result, c);
}

#[test]
fn incoming_flow_enthalpy_can_convert_a_split_tank_to_one_equilibrium_phase() {
    let base = source();
    let feed = composition(&base, 0.8);
    let mut vessel = tank(base, 230.0, 3e6);
    assert!(vessel.state().phase().two_phase().is_some());
    let volume = vessel.volume_m3();
    let c = controls([6e6, 9e6], [230.0, 330.0]);
    let result = vessel
        .receive_stream(portion(&feed, 350.0, 12e6, 0.5), None, c, || false)
        .unwrap();
    close(vessel.mass_kg(), 1.5, 0.0);
    assert_eq!(vessel.volume_m3().to_bits(), volume.to_bits());
    close(
        vessel
            .model()
            .source()
            .thermodynamics()
            .mole_composition()
            .fractions()[METHANE],
        0.613_853_115_444_340_6,
        3e-15,
    );
    close(vessel.state().temperature_k(), 268.893_916_189_584, 4e-6);
    close(vessel.state().pressure_pa(), 7_376_609.639_042_512_5, 0.1);
    close(
        vessel.state().specific_internal_energy_j_kg(),
        169_991.941_643_232_42,
        2e-4,
    );
    close(result.entropy_production_j_k, 218.757_309_203_603_65, 3e-5);
    close(
        result.incoming_pressure_flow_work_j,
        65_316.190_721_704_03,
        3e-5,
    );
    assert!(vessel.state().phase().two_phase().is_none());
    assert_eq!(vessel.state().vapor_mass_fraction(), None);
    audit(&vessel, &result, c);
    // The accepted source is reusable by the existing fixed-inventory U,V owner.
    let temperature = vessel.state().temperature_k();
    let heated = vessel.add_heat(1000.0, c.equilibrium, || false).unwrap();
    close(heated.stored_energy_change_j, 1000.0, 1.5e-5);
    assert!(vessel.state().temperature_k() > temperature);
    assert_eq!(vessel.volume_m3().to_bits(), volume.to_bits());
    close(vessel.mass_kg(), 1.5, 0.0);
}

#[test]
fn dilute_fill_uses_tank_cv_but_incoming_cp_and_conserves_distinct_species() {
    let base = source();
    let mut c = controls([1.5, 2.5], [290.0, 430.0]);
    c.equilibrium.specific_volume_tolerance_m3_kg = 1e-6;
    c.equilibrium.pressure_tolerance_pa = 1e-9;
    for (methane, pr_temperature, pr_pressure, entropy) in [
        (
            0.5,
            377.777_754_933_523_2,
            1.888_888_782_557_249,
            362.458_785_475_907_46,
        ),
        (
            0.8,
            388.100_820_488_529_96,
            2.084_917_110_337_802_4,
            454.944_189_716_135_95,
        ),
    ] {
        let feed = composition(&base, methane);
        let mut vessel = tank(base.clone(), 300.0, 1.0);
        let volume = vessel.volume_m3();
        let result = vessel
            .receive_stream(portion(&feed, 400.0, 10.0, 0.5), None, c, || false)
            .unwrap();
        let mass0 = common::MIXTURE_MASS;
        let mass_in = methane * 0.016_042_46 + (1.0 - methane) * 0.030_069_04;
        let cp_in = (methane * 3.5 + (1.0 - methane) * 4.5) * R / mass_in;
        let cv0 = 3.0 * R / mass0;
        let cv_in = cp_in - R / mass_in;
        let expected_t = (cv0 * 300.0 + 0.5 * cp_in * 400.0) / (cv0 + 0.5 * cv_in);
        let n0 = 1.0 / mass0;
        let nin = 0.5 / mass_in;
        let expected_x = (0.5 * n0 + methane * nin) / (n0 + nin);
        let expected_p = (n0 + nin) * R * expected_t / volume;
        close(vessel.state().temperature_k(), expected_t, 4e-5);
        close(vessel.state().pressure_pa(), expected_p, 4e-7);
        close(vessel.state().temperature_k(), pr_temperature, 4e-6);
        close(vessel.state().pressure_pa(), pr_pressure, 3e-8);
        close(result.entropy_production_j_k, entropy, 3e-5);
        close(
            vessel
                .model()
                .source()
                .thermodynamics()
                .mole_composition()
                .fractions()[METHANE],
            expected_x,
            3e-15,
        );
        assert!(result.incoming_pressure_flow_work_j > 0.0);
        assert_eq!(vessel.volume_m3().to_bits(), volume.to_bits());
        audit(&vessel, &result, c);
    }
}

#[test]
fn caloric_reference_temperature_and_huge_shared_species_datums_leave_physics_unchanged() {
    let plain = source();
    let shifted = model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]);
    let run = |base: &ThermochemicalPengRobinsonModel, reference| {
        let feed = composition(base, 0.8);
        let mut vessel = tank(base.clone(), 230.0, 3e6);
        let mut c = dense_controls();
        c.reference_temperature = reference;
        let report = vessel
            .receive_stream(portion(&feed, 280.0, 8e6, 0.1), None, c, || false)
            .unwrap();
        (vessel, report)
    };
    let (a, ra) = run(&plain, None);
    let (b, rb) = run(&shifted, None);
    assert_eq!(
        a.state().temperature_k().to_bits(),
        b.state().temperature_k().to_bits()
    );
    assert_eq!(
        a.state().pressure_pa().to_bits(),
        b.state().pressure_pa().to_bits()
    );
    assert_eq!(a.state().phase(), b.state().phase());
    assert_eq!(
        ra.energy_residual_j.to_bits(),
        rb.energy_residual_j.to_bits()
    );
    assert_eq!(
        ra.entropy_production_j_k.to_bits(),
        rb.entropy_production_j_k.to_bits()
    );
    assert_eq!(
        ra.ideal_mixing_entropy_j_k.to_bits(),
        rb.ideal_mixing_entropy_j_k.to_bits()
    );
    assert_eq!(ra.pt_evaluations, rb.pt_evaluations);
    assert!(rb.incoming_enthalpy_j.abs() > 1e19);
    assert!(rb.incoming_entropy_j_k.abs() > 1e15);
    let (c, rc) = run(&plain, Some(Temperature::new(240.0)));
    close(c.state().temperature_k(), a.state().temperature_k(), 4e-6);
    close(c.state().pressure_pa(), a.state().pressure_pa(), 0.05);
    close(rc.entropy_production_j_k, ra.entropy_production_j_k, 2e-5);
    close(
        rc.stored_internal_energy_change_j,
        ra.stored_internal_energy_change_j,
        3e-5,
    );
}

#[test]
fn subdividing_a_held_reservoir_portion_preserves_the_conserved_endpoint() {
    let base = source();
    let feed = composition(&base, 0.8);
    let mut whole = tank(base, 230.0, 3e6);
    let mut halves = whole.clone();
    let c = dense_controls();
    let all = whole
        .receive_stream(portion(&feed, 280.0, 8e6, 0.1), None, c, || false)
        .unwrap();
    let first = halves
        .receive_stream(portion(&feed, 280.0, 8e6, 0.05), None, c, || false)
        .unwrap();
    let second = halves
        .receive_stream(portion(&feed, 280.0, 8e6, 0.05), None, c, || false)
        .unwrap();
    close(whole.mass_kg(), halves.mass_kg(), 3e-16);
    assert_eq!(whole.volume_m3().to_bits(), halves.volume_m3().to_bits());
    close(
        whole.state().temperature_k(),
        halves.state().temperature_k(),
        4e-6,
    );
    close(
        whole.state().pressure_pa(),
        halves.state().pressure_pa(),
        0.06,
    );
    close(
        all.entropy_production_j_k,
        first.entropy_production_j_k + second.entropy_production_j_k,
        3e-5,
    );
    close(
        all.stored_internal_energy_change_j,
        first.stored_internal_energy_change_j + second.stored_internal_energy_change_j,
        4e-5,
    );
    for (x, y) in whole
        .model()
        .source()
        .thermodynamics()
        .mole_composition()
        .fractions()
        .iter()
        .zip(
            halves
                .model()
                .source()
                .thermodynamics()
                .mole_composition()
                .fractions(),
        )
    {
        close(*x, *y, 3e-15);
    }
}

#[test]
fn zero_portion_preserves_every_stored_bit_and_spends_no_source_work() {
    let mut vessel = tank(source(), 230.0, 3e6);
    let before = vessel.clone();
    let incompatible = model([1.0, 0.0], [0.0; 2], [100.0, 450.0]);
    let mut c = dense_controls();
    c.equilibrium.maximum_pt_evaluations = 1;
    c.equilibrium.maximum_pressure_evaluations = 1;
    c.equilibrium.phase.enthalpy.maximum_flash_evaluations = 1;
    // No incoming source selection, catalog join or pressure-path admission is
    // needed for zero mass, even though600K is outside that unused source.
    let result = vessel
        .receive_stream(portion(&incompatible, 600.0, 1.0, 0.0), None, c, || false)
        .unwrap();
    assert_eq!(vessel, before);
    assert_eq!(result.pt_evaluations, 0);
    assert_eq!(result.pressure_evaluations, 0);
    assert_eq!(result.eos_evaluation_upper_bound, 0);
    assert_eq!(result.energy_residual_j, 0.0);
    assert_eq!(result.entropy_production_j_k, 0.0);
    assert!(result.inlet_state.is_none());
    assert!(result.species.is_empty());
}

#[test]
fn positive_entropy_does_not_waive_reservoir_pressure_admission_or_atomicity() {
    let base = source();
    let feed = composition(&base, 0.8);
    let mut vessel = tank(base, 230.0, 3e6);
    let before = vessel.clone();
    // Independent unconstrained U,V root:350.682757K,20.0401727MPa,
    // positive Sgen226.9524456J/K. It cannot be supplied passively by12MPa.
    let c = controls([16e6, 24e6], [270.0, 420.0]);
    let refusal = vessel
        .receive_stream(portion(&feed, 350.0, 12e6, 1.0), None, c, || false)
        .unwrap_err();
    match refusal {
        Error::SupplyPressure {
            supply_pressure_pa,
            tank_pressure_pa,
        } => {
            close(supply_pressure_pa, 12e6, 0.0);
            close(tank_pressure_pa, 20_040_172.663_135_927, 0.3);
        }
        other => panic!("expected final supply-pressure refusal: {other:?}"),
    }
    assert_eq!(vessel, before);
    assert!(matches!(
        vessel.receive_stream(
            portion(&feed, 280.0, 2e6, 0.1),
            None,
            dense_controls(),
            || false
        ),
        Err(Error::SupplyPressure { .. })
    ));
    assert_eq!(vessel, before);
}

#[test]
fn original_source_catalogs_are_required_and_conflicts_are_atomic() {
    let base = source();
    let mut vessel = tank(base.clone(), 230.0, 3e6);
    let before = vessel.clone();
    let conflicting = composition(&model([1.0, 0.0], [0.0; 2], [100.0, 450.0]), 0.8);
    assert!(matches!(
        vessel.receive_stream(
            portion(&conflicting, 280.0, 8e6, 0.1),
            None,
            dense_controls(),
            || false
        ),
        Err(Error::Catalog(_))
    ));
    assert_eq!(vessel, before);
    // Disjoint original species require an explicit source for the cross pair;
    // the full binary catalog admits the same finite pure-species intake.
    let ethane = composition(&base, 0.0);
    let methane = composition(&base, 1.0);
    let mut pure = tank(ethane, 300.0, 1.0);
    let pure_before = pure.clone();
    let mut c = controls([1.5, 3.5], [290.0, 440.0]);
    c.equilibrium.specific_volume_tolerance_m3_kg = 1e-6;
    c.equilibrium.pressure_tolerance_pa = 1e-9;
    assert!(matches!(
        pure.receive_stream(portion(&methane, 400.0, 10.0, 0.5), None, c, || false),
        Err(Error::Catalog(_))
    ));
    assert_eq!(pure, pure_before);
    let accepted = pure
        .receive_stream(portion(&methane, 400.0, 10.0, 0.5), Some(&base), c, || {
            false
        })
        .unwrap();
    assert!(accepted.ideal_mixing_entropy_j_k > 0.0);
    assert_eq!(accepted.species.len(), 2);
    let expected_x = (0.5 / 0.016_042_46) / (1.0 / 0.030_069_04 + 0.5 / 0.016_042_46);
    close(
        pure.model()
            .source()
            .thermodynamics()
            .mole_composition()
            .fractions()[METHANE],
        expected_x,
        3e-15,
    );
    audit(&pure, &accepted, c);
}

#[test]
fn shared_native_budget_and_final_cancellation_leave_complete_inventory_unchanged() {
    let base = source();
    let feed = composition(&base, 0.8);
    let initial = tank(base, 230.0, 3e6);
    let c = dense_controls();
    let mut accepted = initial.clone();
    let mut polls = 0usize;
    let report = accepted
        .receive_stream(portion(&feed, 280.0, 8e6, 0.1), None, c, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut exact = c;
    exact.equilibrium.maximum_pt_evaluations = report.pt_evaluations;
    exact.equilibrium.maximum_pressure_evaluations = report.pressure_evaluations;
    let mut repeated = initial.clone();
    let exact_report = repeated
        .receive_stream(portion(&feed, 280.0, 8e6, 0.1), None, exact, || false)
        .unwrap();
    assert_eq!(repeated, accepted);
    assert_eq!(exact_report.pt_evaluations, report.pt_evaluations);
    let mut one_less = exact;
    one_less.equilibrium.maximum_pt_evaluations -= 1;
    let mut refused = initial.clone();
    let error = refused
        .receive_stream(portion(&feed, 280.0, 8e6, 0.1), None, one_less, || false)
        .unwrap_err();
    match error {
        Error::Equilibrium {
            completed_pt_evaluations,
            source: NativeError::WorkLimit { flash_evaluations },
            ..
        } => assert_eq!(
            completed_pt_evaluations + flash_evaluations,
            one_less.equilibrium.maximum_pt_evaluations
        ),
        other => panic!("expected shared PT work refusal: {other:?}"),
    }
    assert_eq!(refused, initial);
    let mut seen = 0usize;
    let error = refused
        .receive_stream(portion(&feed, 280.0, 8e6, 0.1), None, c, || {
            seen += 1;
            seen == polls
        })
        .unwrap_err();
    assert_eq!(error, Error::Cancelled);
    assert_eq!(seen, polls);
    assert_eq!(refused, initial);
    assert!(matches!(
        refused.receive_stream(portion(&feed, 280.0, 8e6, 1e-300), None, c, || false),
        Err(Error::Invalid(_))
    ));
    assert_eq!(refused, initial);
    let mut no_queries = c;
    no_queries.equilibrium.maximum_pt_evaluations = 1;
    assert_eq!(
        refused
            .receive_stream(portion(&feed, 280.0, 8e6, 0.1), None, no_queries, || false)
            .unwrap_err(),
        Error::WorkLimit { pt_evaluations: 1 }
    );
    assert_eq!(refused, initial);
}
