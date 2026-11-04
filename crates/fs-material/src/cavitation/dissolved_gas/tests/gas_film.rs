//! Independent two-resistance flux/partition oracles; synthetic coefficients
//! exercise the native gas/caloric path without claiming empirical gas data.
use super::*;
use crate::cavitation::dissolved_gas::mixture::advance_mixture;
use crate::thermochemical_gas::tests::{species, transport};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

const DG: f64 = 2e-10;

fn mole_fraction_type() -> QuantitySpec {
    typed(QuantityKind::Composition(
        fs_qty::semantic::CompositionBasis::MoleFraction,
    ))
}

fn gas_source(name: &str, diffusivity: f64, pressures: [f64; 2]) -> GasFilmSource {
    let mut property = claim(
        GAS_FILM_DIFFUSIVITY_PROPERTY,
        DIFFUSIVITY_DIMS,
        diffusivity,
        pressures,
    );
    property.validity = ValidityDomain::unconstrained()
        .with_quantity(
            "temperature",
            typed(QuantityKind::AbsoluteTemperature),
            200.0,
            800.0,
        )
        .with_quantity(
            "pressure",
            typed(QuantityKind::Pressure),
            pressures[0],
            pressures[1],
        )
        .with_quantity(GAS_FILM_MOLE_FRACTION_AXIS, mole_fraction_type(), 0.0, 0.1);
    let mut claims = ClaimSet::new();
    claims.insert_claim(property).unwrap();
    GasFilmSource {
        species: id(name),
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: format!("synthetic {name} dilute in declared carrier"),
                phase: "gas".into(),
                process: "authored constant gas-phase diffusivity".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        point: point(T, 100_000.0)
            .with_quantity(GAS_FILM_MOLE_FRACTION_AXIS, mole_fraction_type(), 0.01)
            .unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
        maximum_mole_fraction: 0.05,
    }
}

fn equal_resistance_thickness(h: f64, d: f64) -> f64 {
    DG * film().thickness_m / (h * R * T * d)
}

#[test]
fn two_films_retain_finite_equilibrium_and_match_both_interface_fluxes() {
    let catalog = model(vec![0.01, 0.99], 0.0);
    let thickness = equal_resistance_thickness(H, D);
    let base = law(&catalog);
    let law = base
        .clone()
        .with_gas_film(gas_source("A", DG, [1000.0, 1e7]), thickness)
        .unwrap();
    let authored_film = DissolutionFilm {
        gas_side_transfer: Some(MolarPermeance::try_new(DG / (R * T * thickness)).unwrap()),
        ..film()
    };
    let coefficient = 1.0 / (film().thickness_m / D + H * R * T * thickness / DG);
    let rate = film().area_m2 * coefficient * (1.0 / VL + H * R * T / VG);
    let original_rate = D * film().area_m2 / film().thickness_m * (1.0 / VL + H * R * T / VG);
    near(rate, 0.5 * original_rate, 1e-16);
    for initial_dissolved in [0.0, 0.00008] {
        let mut gas = gas(&catalog);
        let inert = gas.species_mass_kg(&id("B")).unwrap();
        let total = gas.species_mass_kg(&id("A")).unwrap() / MW + initial_dissolved;
        let equilibrium = total * (H * R * T * VL) / (VG + H * R * T * VL);
        let mut liquid = DissolvedGasInventory::new(VL, initial_dissolved, 0.0).unwrap();
        let mut authored_gas = gas.clone();
        let mut authored_liquid = liquid;
        // One physical resistance must not be described twice or selected by
        // precedence. The refusal preserves both original finite inventories.
        assert_eq!(
            law.advance(
                &mut authored_gas,
                &mut authored_liquid,
                authored_film,
                200.0,
                controls(),
                || false,
            ),
            Err(DissolvedGasError::Invalid(
                "choose either a gas-film source or a gas-side coefficient, not both"
            ))
        );
        assert_eq!(authored_gas, gas);
        assert_eq!(authored_liquid, liquid);
        let expected = equilibrium + (initial_dissolved - equilibrium) * (-rate * 200.0).exp();
        let step = law
            .advance(&mut gas, &mut liquid, film(), 200.0, controls(), || false)
            .unwrap();
        let authored = base
            .advance(
                &mut authored_gas,
                &mut authored_liquid,
                authored_film,
                200.0,
                controls(),
                || false,
            )
            .unwrap();
        near(
            authored_liquid.dissolved_moles(),
            liquid.dissolved_moles(),
            1e-16,
        );
        near(
            authored_gas.state().state().pressure,
            gas.state().state().pressure,
            1e-8,
        );
        near(
            authored.relaxation_rate_per_s,
            step.relaxation_rate_per_s,
            1e-18,
        );
        near(
            authored.film_relaxation_ratio,
            step.film_relaxation_ratio,
            1e-18,
        );
        near(
            authored.thermostat_heat_into_system_j,
            step.thermostat_heat_into_system_j,
            1e-12,
        );
        assert_eq!(authored.gas_film_relaxation_ratio, None);
        for (source, coefficient) in step.sources.iter().zip(&authored.sources) {
            near(
                source.molar_flux_mol_m2_s,
                coefficient.molar_flux_mol_m2_s,
                1e-18,
            );
            near(
                source.interface_concentration_mol_m3,
                coefficient.interface_concentration_mol_m3,
                1e-17,
            );
            near(
                source.interface_partial_pressure_pa,
                coefficient.interface_partial_pressure_pa,
                1e-10,
            );
            assert!(coefficient.gas_film_material.is_none());
        }
        near(liquid.dissolved_moles(), expected, 1e-16);
        near(step.relaxation_rate_per_s, rate, 1e-16);
        near(
            step.film_relaxation_ratio,
            rate * film().thickness_m.powi(2) / D,
            1e-16,
        );
        near(
            step.gas_film_relaxation_ratio.unwrap(),
            rate * thickness.powi(2) / DG,
            1e-15,
        );
        for sample in &step.sources {
            let liquid_flux = D / film().thickness_m
                * (sample.interface_concentration_mol_m3 - sample.concentration_mol_m3);
            let gas_flux = DG / thickness
                * (sample.partial_pressure_pa - sample.interface_partial_pressure_pa)
                / (R * T);
            near(sample.molar_flux_mol_m2_s, liquid_flux, 1e-18);
            near(sample.molar_flux_mol_m2_s, gas_flux, 1e-18);
            near(
                sample.interface_concentration_mol_m3,
                H * sample.interface_partial_pressure_pa,
                1e-17,
            );
            near(
                sample.interface_concentration_mol_m3,
                0.5 * (sample.equilibrium_concentration_mol_m3 + sample.concentration_mol_m3),
                1e-17,
            );
            assert!(sample.gas_film_material.is_some());
            assert_eq!(sample.molar_flux_mol_m2_s > 0.0, initial_dissolved == 0.0);
        }
        near(
            gas.species_mass_kg(&id("A")).unwrap() / MW + liquid.dissolved_moles(),
            total,
            1e-15,
        );
        near(gas.species_mass_kg(&id("B")).unwrap(), inert, 1e-16);
        near(
            step.thermostat_heat_into_system_j,
            (expected - initial_dissolved) * (R * T - SOLUTION_HEAT),
            1e-12,
        );
        near(step.energy_residual_j, 0.0, 1e-8);
        assert_eq!(gas.state().state().temperature, T);
        law.advance(&mut gas, &mut liquid, film(), 1e6, controls(), || false)
            .unwrap();
        near(liquid.dissolved_moles(), equilibrium, 1e-15);
    }
    // A finite gas resistance lowers the actual liquid-side trace. The
    // bulk-gas Henry proxy need not be a physically present concentration.
    let mut dilute = controls();
    dilute.maximum_concentration_mol_m3 = 0.012;
    let mut gas = gas(&catalog);
    let mut liquid = DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap();
    let step = law
        .advance(&mut gas, &mut liquid, film(), 200.0, dilute, || false)
        .unwrap();
    assert!(step.sources[0].equilibrium_concentration_mol_m3 > dilute.maximum_concentration_mol_m3);
    assert!(step.sources[0].interface_concentration_mol_m3 < dilute.maximum_concentration_mol_m3);
}

#[test]
fn two_film_subdivision_seals_and_liquid_only_limit_use_the_original_owner() {
    let catalog = model(vec![0.01, 0.99], 0.0);
    let base = law(&catalog);
    let two = base
        .clone()
        .with_gas_film(
            gas_source("A", DG, [1000.0, 1e7]),
            equal_resistance_thickness(H, D),
        )
        .unwrap();
    let run = |law: &HenryDissolution, times: &[f64]| {
        let mut gas = gas(&catalog);
        let mut liquid = DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap();
        let mut steps = Vec::new();
        for &end in times {
            steps.push(
                law.advance(&mut gas, &mut liquid, film(), end, controls(), || false)
                    .unwrap(),
            );
        }
        (gas, liquid, steps)
    };
    let full = run(&two, &[400.0]);
    let split = run(&two, &[200.0, 400.0]);
    near(full.1.dissolved_moles(), split.1.dissolved_moles(), 1e-16);
    near(
        full.0.state().state().pressure,
        split.0.state().state().pressure,
        1e-8,
    );
    assert_eq!(split, run(&two, &[200.0, 400.0]));

    let thin = base
        .clone()
        .with_gas_film(gas_source("A", 2e-5, [1000.0, 1e7]), 1e-5)
        .unwrap();
    let original = run(&base, &[200.0]);
    let nearly_original = run(&thin, &[200.0]);
    near(
        nearly_original.2[0].relaxation_rate_per_s / original.2[0].relaxation_rate_per_s,
        1.0,
        6e-6,
    );
    assert_eq!(original.2[0].gas_film_relaxation_ratio, None);

    let sealed = base
        .with_gas_film(gas_source("A", 0.0, [1000.0, 1e7]), 1e-5)
        .unwrap();
    let no_transfer = run(&sealed, &[200.0]);
    assert_eq!(no_transfer.1.dissolved_moles(), 0.0);
    assert_eq!(no_transfer.0, gas(&catalog));
    assert_eq!(
        no_transfer.2[0].sources[0].interface_concentration_mol_m3,
        0.0
    );
    assert_eq!(
        no_transfer.2[0].sources[0].interface_partial_pressure_pa,
        0.0
    );
    assert_eq!(no_transfer.2[0].sources[0].molar_flux_mol_m2_s, 0.0);

    let sealed_liquid = HenryDissolution::new(
        &catalog,
        id("A"),
        source(0.0, [1000.0, 1e7]),
        calorics(&catalog),
    )
    .unwrap();
    assert!(
        sealed_liquid
            .clone()
            .with_gas_film(gas_source("A", 0.0, [1000.0, 1e7]), 1e-5)
            .is_err()
    );
    let sealed_liquid = sealed_liquid
        .with_gas_film(gas_source("A", DG, [1000.0, 1e7]), 1e-5)
        .unwrap();
    let no_transfer = run(&sealed_liquid, &[200.0]);
    assert_eq!(no_transfer.1.dissolved_moles(), 0.0);
    near(
        no_transfer.2[0].sources[0].interface_partial_pressure_pa,
        1000.0,
        1e-10,
    );
}

#[test]
fn gas_film_source_limits_timescales_and_final_cancellation_are_atomic() {
    let catalog = model(vec![0.01, 0.99], 0.0);
    let base = law(&catalog);
    let thickness = equal_resistance_thickness(H, D);
    // Both coefficients can be finite while their positive ratio falls below
    // f64's range. A conducting film must not disappear from the interface.
    for (liquid_d, gas_d) in [(1e200, 1e-200), (1e-200, 1e200)] {
        let extreme = HenryDissolution::new(
            &catalog,
            id("A"),
            source(liquid_d, [1000.0, 1e7]),
            calorics(&catalog),
        )
        .unwrap()
        .with_gas_film(gas_source("A", gas_d, [1000.0, 1e7]), 1.0)
        .unwrap();
        let mut gas = gas(&catalog);
        let before_gas = gas.clone();
        let mut liquid = DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap();
        let before_liquid = liquid;
        assert_eq!(
            extreme.advance(
                &mut gas,
                &mut liquid,
                DissolutionFilm {
                    area_m2: 1e-7,
                    thickness_m: 1.0,
                    gas_side_transfer: None,
                },
                200.0,
                controls(),
                || false,
            ),
            Err(DissolvedGasError::Invalid(
                "positive two-film interface weight underflows"
            ))
        );
        assert_eq!(gas, before_gas);
        assert_eq!(liquid, before_liquid);
        assert!(matches!(
            extreme
                .gas_film
                .as_ref()
                .unwrap()
                .sample(100_000.0, f64::from_bits(1)),
            Err(DissolvedGasError::Invalid(
                "positive gas-film mole fraction underflows"
            ))
        ));
    }
    assert!(
        base.clone()
            .with_gas_film(gas_source("B", DG, [1000.0, 1e7]), thickness)
            .is_err()
    );
    let mut wrong_basis = gas_source("A", DG, [1000.0, 1e7]);
    wrong_basis.point = point(T, 100_000.0)
        .with_quantity(
            GAS_FILM_MOLE_FRACTION_AXIS,
            typed(QuantityKind::Composition(
                fs_qty::semantic::CompositionBasis::MassFraction,
            )),
            0.01,
        )
        .unwrap();
    assert!(base.clone().with_gas_film(wrong_basis, thickness).is_err());
    let initial_gas = gas(&catalog);
    let initial_liquid = DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap();
    let mut gas = initial_gas.clone();
    let mut liquid = initial_liquid;
    let narrow = base
        .clone()
        .with_gas_film(gas_source("A", DG, [99_995.0, 100_010.0]), thickness)
        .unwrap();
    assert!(matches!(
        narrow.advance(&mut gas, &mut liquid, film(), 200.0, controls(), || false),
        Err(DissolvedGasError::Material(_))
    ));
    assert_eq!(gas, initial_gas);
    assert_eq!(liquid, initial_liquid);
    let thick = base
        .clone()
        .with_gas_film(gas_source("A", DG, [1000.0, 1e7]), 0.1)
        .unwrap();
    assert!(matches!(
        thick.advance(&mut gas, &mut liquid, film(), 200.0, controls(), || false),
        Err(DissolvedGasError::GasFilmRelaxation { .. })
    ));
    let law = base
        .with_gas_film(gas_source("A", DG, [1000.0, 1e7]), thickness)
        .unwrap();
    let mut concentrated = DissolvedGasInventory::new(VL, 0.001, 0.0).unwrap();
    let before_concentrated = concentrated;
    assert!(matches!(
        law.advance(
            &mut gas,
            &mut concentrated,
            film(),
            200.0,
            controls(),
            || false
        ),
        Err(DissolvedGasError::GasFilmConcentration { .. })
    ));
    assert_eq!(concentrated, before_concentrated);
    assert_eq!(gas, initial_gas);
    let mut polls = 0;
    let expected = law
        .advance(&mut gas, &mut liquid, film(), 200.0, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut gas = initial_gas.clone();
    let mut liquid = initial_liquid;
    let mut cancelled_polls = 0;
    assert_eq!(
        law.advance(&mut gas, &mut liquid, film(), 200.0, controls(), || {
            cancelled_polls += 1;
            cancelled_polls == polls
        }),
        Err(DissolvedGasError::Cancelled)
    );
    assert_eq!(gas, initial_gas);
    assert_eq!(liquid, initial_liquid);
    assert_eq!(
        law.advance(&mut gas, &mut liquid, film(), 200.0, controls(), || false)
            .unwrap(),
        expected
    );
}

#[test]
fn simultaneous_two_film_exchange_preserves_carrier_and_bounds_combined_dilution() {
    let catalog = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species("A", MW, 3.5, 0.0),
                species("B", 0.004, 2.5, 0.0),
                species("C", 0.04, 3.5, 0.0),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.01, 0.015, 0.975]).unwrap(),
        )
        .unwrap(),
        vec![
            transport("A", 1.7e-5, 0.025),
            transport("B", 2e-5, 0.15),
            transport("C", 1.8e-5, 0.04),
        ],
        None,
    )
    .unwrap();
    let specs = [("A", H, D, 6000.0), ("B", 0.5 * H, 3.0 * D, 2500.0)];
    let laws: Vec<_> = specs
        .iter()
        .map(|&(name, h, d, heat)| {
            super::mixture::species_law(&catalog, name, h, d, heat, [1000.0, 1e7])
                .with_gas_film(
                    gas_source(name, DG, [1000.0, 1e7]),
                    equal_resistance_thickness(h, d),
                )
                .unwrap()
        })
        .collect();
    let initial_gas = gas(&catalog);
    let initial_liquid = [
        DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap(),
        DissolvedGasInventory::new(VL, 0.00012, 0.0).unwrap(),
    ];
    let mut gas = initial_gas.clone();
    let mut liquids = initial_liquid;
    let report = advance_mixture(
        &laws,
        &mut gas,
        &mut liquids,
        &[film(); 2],
        200.0,
        controls(),
        || false,
    )
    .unwrap();
    assert!(report.exchanges[0].dissolved_moles_change > 0.0);
    assert!(report.exchanges[1].dissolved_moles_change < 0.0);
    for (index, &(name, h, d, _)) in specs.iter().enumerate() {
        let mw = [MW, 0.004][index];
        let total = initial_gas.species_mass_kg(&id(name)).unwrap() / mw
            + initial_liquid[index].dissolved_moles();
        let equilibrium = total * (h * R * T * VL) / (VG + h * R * T * VL);
        let coefficient =
            1.0 / (film().thickness_m / d + h * R * T * equal_resistance_thickness(h, d) / DG);
        let rate = coefficient * film().area_m2 * (1.0 / VL + h * R * T / VG);
        let expected = equilibrium
            + (initial_liquid[index].dissolved_moles() - equilibrium) * (-rate * 200.0).exp();
        near(liquids[index].dissolved_moles(), expected, 1e-16);
        near(
            gas.species_mass_kg(&id(name)).unwrap() / mw + liquids[index].dissolved_moles(),
            total,
            1e-15,
        );
    }
    near(
        gas.species_mass_kg(&id("C")).unwrap(),
        initial_gas.species_mass_kg(&id("C")).unwrap(),
        1e-16,
    );
    near(report.energy_residual_j, 0.0, 1e-8);

    let mut limited = laws.clone();
    let mut source = limited[0].gas_film_source().unwrap().clone();
    source.maximum_mole_fraction = 0.04;
    limited[0] = limited[0]
        .clone()
        .with_gas_film(source, equal_resistance_thickness(H, D))
        .unwrap();
    let mut gas = initial_gas.clone();
    let mut liquids = initial_liquid;
    assert!(matches!(
        advance_mixture(
            &limited,
            &mut gas,
            &mut liquids,
            &[film(); 2],
            200.0,
            controls(),
            || false
        ),
        Err(DissolvedGasError::GasFilmConcentration { .. })
    ));
    assert_eq!(gas, initial_gas);
    assert_eq!(liquids, initial_liquid);

    let mut gas_reversed = initial_gas;
    let mut liquids_reversed = [initial_liquid[1], initial_liquid[0]];
    let reversed = advance_mixture(
        &[laws[1].clone(), laws[0].clone()],
        &mut gas_reversed,
        &mut liquids_reversed,
        &[film(); 2],
        200.0,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(reversed, report);
}
