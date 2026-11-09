//! Source-backed reacting rigid storage through the existing thermal timeline.
use super::conduction::GasIsochoricEquilibriumConductionLaw;
use super::*;
use crate::thermal::conduction::ConductivityCurve;
use crate::thermochemical_gas::energy::equilibrium::tests::{fixture, isomers};
use crate::thermochemical_gas::tests::{
    point, species, species_with_energy_reference, transport, transport_in_domain,
};
use fs_thermal::enthalpy::conduction::implicit::ImplicitConductionLaw;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositionBudget, DepositionSimulation, Event, MaterialId,
    StepConfig, TimedEvent,
};
use fs_thermochem::equilibrium::composition::{
    EquilibriumSolveV1,
    volume::{EquilibriumInternalEnergySolveV1, EquilibriumVolumeSolveV1},
};
use fs_thermochem::{
    Nasa9RegionV1, Nasa9StandardStateModelV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

fn controls() -> DrivenEquilibriumGasSolve {
    let mut chemistry = EquilibriumSolveV1::default();
    chemistry.relative_balance_tolerance = 1e-13;
    chemistry.dimensionless_potential_tolerance = 1e-12;
    DrivenEquilibriumGasSolve {
        equilibrium: EquilibriumInternalEnergySolveV1 {
            volume: EquilibriumVolumeSolveV1 {
                equilibrium: chemistry,
                relative_volume_tolerance: 1e-12,
                maximum_pressure_iterations: 96,
            },
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1500.0,
            energy_tolerance_j: 1e-7,
            maximum_temperature_iterations: 128,
        },
        mass_tolerance_kg: 1e-12,
        energy_tolerance_j: 1e-5,
        entropy_tolerance_j_k: 1e-8,
    }
}
fn feed() -> [(SpeciesId, Amount); 1] {
    [(SpeciesId::new("A").unwrap(), Amount::new(2.0))]
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn fraction(t: f64) -> f64 {
    1.0 / (1.0 + (-24000.0 / (R * t)).exp())
}
fn u(t: f64) -> f64 {
    (3.0 * R * (t - 200.0) - 24000.0 * (fraction(t) - fraction(200.0))) / 0.01
}
fn cv(t: f64) -> f64 {
    let y = fraction(t);
    (3.0 * R + 24000.0_f64.powi(2) * y * (1.0 - y) / (R * t * t)) / 0.01
}
fn entropy(t: f64) -> f64 {
    let y = fraction(t);
    R / 0.01 * (3.0 * t.ln() - y * y.ln() - (1.0 - y) * (1.0 - y).ln())
}
fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: 400.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    }
}
fn cell(mass: f64, volume: f64, energy: f64, source: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(1),
        segment_id: 1,
        body_id: 1,
        material_id: MaterialId(1),
        mass_kg: mass,
        specific_h_j_kg: energy,
        geometry: CellGeometry {
            start_m: [0.0, 0.0, 0.0],
            end_m: [1.0, 0.0, 0.0],
            volume_m3: volume,
            surface_area_m2: 1.0,
        },
        boundary: boundary(source),
    }
}

#[test]
fn g1_rigid_reacting_storage_retains_cv_pressure_entropy_and_original_materials() {
    let (bank, chemistry) = isomers();
    let original = bank.clone();
    let context = point(900.0, 1e5).with("sample_id", 7.0).unwrap();
    let law = GasIsochoricEquilibriumThermalStorage::new(
        &bank,
        &chemistry,
        &context,
        &feed(),
        Volume::new(0.1),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(law.storage_basis(), ThermalStorageBasis::InternalEnergy);
    assert_eq!(law.feed_mass_kg(), 0.02);
    assert_eq!(law.density_kg_m3(), 0.02 / 0.1);
    let source_bound = 0.01 / (3.0 * R);
    assert!(law.maximum_temperature_slope() >= source_bound);
    assert!(law.maximum_temperature_slope() <= source_bound * (1.0 + 1e-10));
    for t in [200.0, 400.0, 900.0, 1300.0, 1500.0] {
        let (actual_u, actual_cv) = law
            .internal_energy_and_heat_capacity_at_temperature(Temperature::new(t), || false)
            .unwrap();
        close(actual_u, u(t), 2e-5);
        close(actual_cv, cv(t), 2e-7);
        close(law.temperature_k(actual_u).unwrap(), t, 2e-8);
        assert_eq!(law.bulk_density_kg_m3(actual_u).unwrap(), Some(0.02 / 0.1));
        let state = law
            .caloric_state_at_specific_internal_energy(actual_u, || false)
            .unwrap();
        close(state.pressure.value(), 2.0 * R * t / 0.1, 2e-6);
        assert!(1.0 / actual_cv <= law.maximum_temperature_slope());
    }
    let result = law
        .resolve_specific_internal_energy(u(900.0), || false)
        .unwrap();
    let gas = result.chamber().parcel();
    assert_eq!(gas.mass_kg(), 0.02);
    assert_eq!(gas.volume_m3(), 0.1);
    assert_eq!(
        gas.model().thermodynamics().components(),
        chemistry.models()
    );
    assert_eq!(
        gas.model().transport_materials(),
        bank.transport_materials()
    );
    assert_eq!(gas.point().axes()["sample_id"], 7.0);
    close(
        result.internal_energy_change_j(),
        0.02 * u(900.0),
        controls().energy_tolerance_j,
    );
    assert!(result.energy_residual_j().abs() <= controls().energy_tolerance_j);
    let frozen_cv = gas
        .state()
        .thermodynamics()
        .properties()
        .mass_specific_cv()
        .value();
    assert!(result.native_sample().response.heat_capacity_v_j_k > 0.02 * frozen_cv);
    close(
        law.thermal_entropy_change_j_kg_k(u(400.0), u(1300.0), || false)
            .unwrap(),
        entropy(1300.0) - entropy(400.0),
        2e-7,
    );
    close(
        law.thermal_entropy_change_j_kg_k(u(1300.0), u(400.0), || false)
            .unwrap(),
        entropy(400.0) - entropy(1300.0),
        2e-7,
    );
    assert!(law.advected_composition().is_none());
    assert_eq!(bank, original);
}

#[test]
fn g1_existing_timeline_heats_and_cools_reacting_rigid_cells_without_boundary_work() {
    let (bank, chemistry) = isomers();
    let law = GasIsochoricEquilibriumThermalStorage::new(
        &bank,
        &chemistry,
        &point(900.0, 1e5),
        &feed(),
        Volume::new(0.1),
        controls(),
        || false,
    )
    .unwrap();
    // Scale both mass and volume; the same intensive reacting law applies.
    let mass = 0.03;
    let volume = mass / law.density_kg_m3();
    let heat = mass * (u(1200.0) - u(400.0));
    let birth = cell(mass, volume, u(400.0), heat);
    law.validate_cell(&birth).unwrap();
    let original = DepositionSimulation::new(
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(birth.clone()),
            },
            TimedEvent {
                time_s: 1.0,
                event: Event::Boundary(CellId(1), boundary(-heat)),
            },
        ],
        &[(MaterialId(1), &law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 2,
            max_substeps: 128,
        },
    )
    .unwrap();
    let config = StepConfig {
        max_dt_s: 0.125,
        energy_tolerance_j: 1e-6,
    };
    let (initial, _) = original.advance_to(0.0, config, || Ok(())).unwrap();
    assert_eq!(
        initial.material_storage_basis(MaterialId(1)).unwrap(),
        ThermalStorageBasis::InternalEnergy
    );
    let (hot, audit) = initial.advance_to(1.0, config, || Ok(())).unwrap();
    close(hot.temperature_k(CellId(1)).unwrap(), 1200.0, 2e-7);
    close(audit.stored_energy_change_j, heat, 1e-7);
    close(audit.external_heat_j, heat, 1e-7);
    assert!(audit.residual_j.abs() <= config.energy_tolerance_j);
    assert_eq!(hot.mass_kg(), mass);
    assert_eq!(hot.cells()[&CellId(1)].geometry, birth.geometry);
    let material = law
        .resolve_specific_internal_energy(hot.cells()[&CellId(1)].specific_h_j_kg, || false)
        .unwrap();
    close(
        material.chamber().parcel().state().state().pressure,
        2.0 * R * 1200.0 / 0.1,
        2e-5,
    );
    let (cooled, audit) = hot.advance_to(2.0, config, || Ok(())).unwrap();
    close(cooled.temperature_k(CellId(1)).unwrap(), 400.0, 2e-7);
    close(audit.external_heat_j, -heat, 1e-7);
    assert!(audit.residual_j.abs() <= config.energy_tolerance_j);
    assert_eq!(cooled.cells()[&CellId(1)].geometry, birth.geometry);
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn g1_rigid_conduction_uses_dissociation_cv_without_isobaric_expansion() {
    let mut models = Vec::new();
    for (name, mass, cp_r) in [("A", 0.001, 2.5), ("A2", 0.002, 3.5)] {
        let source = species(name, mass, cp_r, 0.0);
        let mut card = source.regions()[0].card().clone();
        card.parameters.get_mut("a8").unwrap().value = cp_r * (1.0 - 1000.0_f64.ln());
        models.push(
            Nasa9StandardStateModelV1::new(
                source.species().clone(),
                source.molar_mass(),
                source.convention().clone(),
                vec![Nasa9RegionV1::from_card(card).unwrap()],
            )
            .unwrap(),
        );
    }
    let (bank, chemistry) = fixture(
        models,
        vec![transport("A", 1.7e-5, 0.025), transport("A2", 2.0e-5, 0.15)],
        &["H"],
        vec![vec![1, 2]],
    );
    let law = GasIsochoricEquilibriumThermalStorage::new(
        &bank,
        &chemistry,
        &point(900.0, 1e5),
        &[(SpeciesId::new("A2").unwrap(), Amount::new(1.0))],
        Volume::new(0.1),
        controls(),
        || false,
    )
    .unwrap();
    let conductor = GasIsochoricEquilibriumConductionLaw::try_new(
        &law,
        ConductivityCurve::new(vec![[200.0, 0.1], [800.0, 0.22], [1500.0, 0.5]]).unwrap(),
    )
    .unwrap();
    let (energy, _) = law
        .internal_energy_and_heat_capacity_at_temperature(Temperature::new(1000.0), || false)
        .unwrap();
    let state = conductor.state_at_specific_internal_energy(energy).unwrap();
    let q = R * 1000.0 / (0.1 * 1e5);
    let alpha = 2.0 / (1.0 + (1.0 + 16.0 * q).sqrt());
    let da = 0.5 * alpha * (1.0 - alpha) / (1000.0 * (2.0 - alpha));
    let capacity = R * (2.5 + 0.5 * alpha + 500.0 * da);
    close(state.temperature_slope_kg_k_j, 0.002 / capacity, 2e-13);
    close(state.conductivity_w_m_k, 0.3, 1e-11);
    close(
        state.conductivity_slope_kg_w_m_k_j,
        0.0004 * 0.002 / capacity,
        2e-14,
    );
    assert_eq!(
        conductor
            .specific_volume_derivative_at_specific_enthalpy(energy)
            .unwrap(),
        Some(0.0)
    );
    assert_eq!(
        conductor.caloric_law().storage_basis(),
        ThermalStorageBasis::InternalEnergy
    );
    assert_eq!(
        conductor.enthalpy_bounds_j_kg(),
        law.specific_internal_energy_domain_j_kg()
    );
    assert!(
        GasIsochoricEquilibriumConductionLaw::try_new(
            &law,
            ConductivityCurve::new(vec![[300.0, 0.1], [1500.0, 0.5]]).unwrap()
        )
        .is_err()
    );
}

#[test]
fn g4_rigid_density_and_evolving_pressure_preserve_original_admission() {
    let (bank, chemistry) = isomers();
    let law = GasIsochoricEquilibriumThermalStorage::new(
        &bank,
        &chemistry,
        &point(900.0, 1e5),
        &feed(),
        Volume::new(0.1),
        controls(),
        || false,
    )
    .unwrap();
    let valid = cell(0.03, 0.15, u(400.0), 0.0);
    law.validate_cell(&valid).unwrap();
    for (mass, volume) in [
        (0.03, 0.16),
        (0.04, 0.15),
        (0.0, 0.15),
        (0.03, f64::INFINITY),
    ] {
        assert!(
            law.validate_cell(&cell(mass, volume, u(400.0), 0.0))
                .is_err()
        );
    }
    let (limited, _) = fixture(
        chemistry.models().to_vec(),
        vec![
            transport_in_domain("A", 1.7e-5, 0.025, [200.0, 2000.0], [1000.0, 1e5]),
            transport_in_domain("B", 2.0e-5, 0.15, [200.0, 2000.0], [1000.0, 1e5]),
        ],
        &["H"],
        vec![vec![1, 1]],
    );
    // Cold pressure is supported, but the sealed gas's hot pressure is not.
    // The context pressure must not masquerade as a fixed-pressure constraint.
    assert!(
        GasIsochoricEquilibriumThermalStorage::new(
            &limited,
            &chemistry,
            &point(900.0, 1e5),
            &feed(),
            Volume::new(0.1),
            controls(),
            || false,
        )
        .is_err()
    );
}

#[test]
fn g4_failed_and_cancelled_queries_do_not_mutate_the_rigid_material_law() {
    let (bank, chemistry) = isomers();
    let original = bank.clone();
    let make = || {
        GasIsochoricEquilibriumThermalStorage::new(
            &bank,
            &chemistry,
            &point(900.0, 1e5),
            &feed(),
            Volume::new(0.1),
            controls(),
            || false,
        )
        .unwrap()
    };
    let law = make();
    let reference = law.reference_state().clone();
    let expected = law
        .caloric_state_at_specific_internal_energy(u(900.0), || false)
        .unwrap();
    for energy in [
        -1.0,
        f64::NAN,
        f64::INFINITY,
        law.specific_internal_energy_domain_j_kg()[1].next_up(),
    ] {
        assert!(law.temperature_k(energy).is_err());
        assert!(law.bulk_density_kg_m3(energy).is_err());
        assert_eq!(
            law.caloric_state_at_specific_internal_energy(u(900.0), || false)
                .unwrap(),
            expected
        );
    }
    let mut count = 0;
    law.resolve_specific_internal_energy(u(900.0), || {
        count += 1;
        false
    })
    .unwrap();
    let cancelled = make();
    let mut polls = 0;
    assert!(matches!(
        cancelled.resolve_specific_internal_energy(u(900.0), || {
            polls += 1;
            polls == count
        }),
        Err(Error::Cancelled)
    ));
    assert_eq!(polls, count);
    assert_eq!(cancelled.reference_state(), &reference);
    assert_eq!(
        cancelled
            .caloric_state_at_specific_internal_energy(u(900.0), || false)
            .unwrap(),
        expected
    );
    assert_eq!(bank, original);
}

#[test]
fn g3_shared_freezing_retains_amount_correction_before_large_formation_energy() {
    let source = species_with_energy_reference("A", 0.01, 4.0, 0.0, 1e13);
    let (bank, chemistry) = fixture(
        vec![source.clone()],
        vec![transport("A", 1.7e-5, 0.025)],
        &["H"],
        vec![vec![1]],
    );
    let amount = 1.1;
    let feed = [(SpeciesId::new("A").unwrap(), Amount::new(amount))];
    let mass = amount * source.molar_mass().value();
    let mut c = controls();
    c.energy_tolerance_j = 1e-3; // Admit the separately reported absolute-U rounding.
    let native = chemistry
        .equilibrate_tv(
            Temperature::new(900.0),
            Volume::new(0.1),
            &feed,
            c.equilibrium.volume,
            || false,
        )
        .unwrap();
    let native_amount = native.equilibrium().species()[0].amount.value();
    let (_, defect, _) = freeze(
        &bank,
        &chemistry,
        &point(900.0, 1e5),
        &feed,
        native,
        mass,
        c,
        &mut || false,
    )
    .unwrap();
    // Independent fused-product remainder: the rounded mass/M equals amount,
    // but their exact represented mass difference is nonzero. Rounding that
    // difference before multiplying by formation energy would erase it.
    let molar_mass = source.molar_mass().value();
    assert_eq!(mass / molar_mass, amount);
    let amount_correction = (-native_amount).mul_add(molar_mass, mass) / molar_mass;
    let expected = amount_correction
        * source
            .evaluate_properties(Temperature::new(900.0))
            .unwrap()
            .u()
            .value();
    let rounded_defect = (mass / molar_mass - native_amount)
        * source
            .evaluate_properties(Temperature::new(900.0))
            .unwrap()
            .u()
            .value();
    assert!((expected - rounded_defect).abs() > 1e-7);
    close(defect, expected, 1e-16);
}
