//! A source-bound atomic ionization substep changes the original FEEC current.
//! NIST BEB supplies the event rate; the momentum cross section is explicitly
//! manufactured. The two powers calculated from the initial field are HELD
//! over the chemistry interval. The final fixed-voltage field is a subsequent
//! feedback evaluation, not the time integral of an evolving voltage source.
//! Elastic exchange is a separate owner; this test isolates chemical storage.

mod common;
#[path = "../../fs-material/tests/common/ionized_gas.rs"]
mod fixture;

use common::{config, ends, layered, limits};
use fs_em::conduction::{ConductionControls, PlanarConduction};
use fs_gas::chemistry::ionization::{AtomicIonizationSolve, advance_atomic_ionization};
use fs_material::ionized_gas::IonizedGasConductivity;
use fs_qty::{Power, Temperature, Time};
use fs_thermochem::electron_collision::{MaxwellianRateControlsV1, nist};
use fs_thermochem::equilibrium::composition::state::FrozenIonizedGasStateV1;
use fs_thermochem::kinetics::ElectronImpactRateV1;
use fs_thermochem::kinetics::atomic_ionization::{AtomicIonizationV1, AtomicStatisticalWeightsV1};
use fs_thermochem::{ReactionId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn material(state: &FrozenIonizedGasStateV1) -> IonizedGasConductivity {
    let electron = fixture::electron_relaxation(state, 1e-22);
    let ion = fixture::source_two_temperature(
        "H+",
        state.heavy_temperature().value(),
        state.electron_temperature().value(),
        state.pressure().value(),
        1e9,
    );
    IonizedGasConductivity::from_state(
        state.clone(),
        vec![electron, ion],
        fixture::limits(),
        || false,
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
fn g1_spectral_ionization_preserves_band_energy_through_current_and_heating() {
    use fs_gas::chemistry::photoionization::spectrum::{
        SpectralPhotoionizationSolve, advance_spectral_photoionization,
    };
    use fs_math::c64::C64;
    use fs_qty::{Area, Energy};
    use fs_thermochem::kinetics::atomic_ionization::photoionization::{
        AtomicPhotoionizationV1, PhotonNumberFluxV1,
    };

    // Different photon energies see different original H I opacities. The
    // deliberately broad beam keeps both transmitted bands measurable.
    let initial = fixture::frozen_state(0.0, 6000.0, 16000.0, 0.2);
    let channel = AtomicIonizationV1::new(
        initial.thermodynamics().clone(),
        ReactionId::new("spectral-hydrogen-electrical-feedback").unwrap(),
        fixture::id("H+"),
        ElectronImpactRateV1::new(
            fixture::id("e-"),
            nist::hydrogen(fixture::id("H")),
            MaxwellianRateControlsV1 {
                relative_tolerance: 1e-10,
                maximum_subdivisions: 4096,
            },
            0.1,
        )
        .unwrap(),
        AtomicStatisticalWeightsV1 { neutral: 2, ion: 1 },
        1e-12,
    )
    .unwrap();
    let source = AtomicPhotoionizationV1::verner_hydrogen(channel.clone()).unwrap();
    let area = 1000.0;
    let duration = 1e-7;
    let groups = [
        PhotonNumberFluxV1 {
            photon_energy: Energy::new(14.0 * fixture::E),
            number_flux_per_m2_s: 3e-3 * fixture::NA / (area * duration),
        },
        PhotonNumberFluxV1 {
            photon_energy: Energy::new(17.0 * fixture::E),
            number_flux_per_m2_s: 4e-3 * fixture::NA / (area * duration),
        },
    ];
    let mut controls =
        SpectralPhotoionizationSolve::new(Temperature::new(2000.0), Temperature::new(20000.0));
    controls.relative_tolerance = 1e-8;
    controls.initial_step_s = duration;
    controls.maximum_step_s = duration;
    let photo = advance_spectral_photoionization(
        &initial,
        &source,
        &groups,
        Area::new(area),
        Time::new(duration),
        controls,
        || false,
    )
    .unwrap();
    let state = photo.state();
    assert_eq!(state.thermodynamics(), initial.thermodynamics());
    assert_eq!(state.heavy_temperature(), initial.heavy_temperature());
    assert_eq!(state.amounts()[1], state.amounts()[2]);
    fixture::relative(
        state.amounts()[2].value(),
        photo.reaction_extent().value(),
        1e-13,
    );
    fixture::relative(
        state.amounts()[0].value() + state.amounts()[1].value(),
        1.0,
        1e-13,
    );

    let band = photo.group_balances();
    assert!(band[0].transmitted_photons > 0.0 && band[1].transmitted_photons > 0.0);
    assert!(
        band[1].transmitted_photons / band[1].incident_photons
            > band[0].transmitted_photons / band[0].incident_photons
    );
    let absorbed = band
        .iter()
        .zip(groups)
        .map(|(b, g)| {
            fixture::relative(
                b.absorbed_photons + b.transmitted_photons,
                b.incident_photons,
                1e-12,
            );
            b.absorbed_photons * g.photon_energy.value()
        })
        .sum::<f64>();
    // Independently reconstruct the original chemical and electron sensible
    // stores. No average incident photon energy or second ionization debit.
    let electrons = state.amounts()[2].value();
    let chemical = electrons * channel.ionization_energy_j_mol();
    let sensible = 1.5 * R * electrons * state.electron_temperature().value();
    assert!((chemical + sensible - absorbed).abs() < controls.energy_tolerance_j);
    assert!(
        (state.internal_energy().value() - initial.internal_energy().value() - absorbed).abs()
            < 2e-7
    );

    let conductor = material(state);
    assert_eq!(conductor.state(), state);
    assert!(conductor.equilibrium().is_none());
    let electron = conductor
        .carriers()
        .iter()
        .find(|c| c.is_electron())
        .unwrap();
    fixture::relative(
        electron.number_density_per_m3(),
        electrons * fixture::NA / 0.2,
        1e-13,
    );
    let sigma = conductor.conductivity_s_m();
    assert!(sigma > 0.0);
    let harmonic = conductor.prepare_harmonic(|| false).unwrap();
    let response = harmonic
        .response(1e12, [C64::from_re(0.001), C64::ZERO, C64::ZERO], || false)
        .unwrap();
    fixture::relative(
        response.mean_electrical_power_w_m3,
        response.mean_heating_w_m3.total_power_density_w_m3(),
        1e-12,
    );
    assert!(response.mean_drift_kinetic_energy_j_m3 > 0.0);

    let operator = PlanarConduction::new(
        layered(0.5, 1.0, 0.2),
        vec![conductor; 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let field = operator
        .solve(&ends(2.0), field_controls(), &mut || false)
        .unwrap();
    fixture::relative(field.current_density_a_m2()[0][0], -2.0 * sigma, 1e-12);
    fixture::relative(field.power_w(), 0.8 * sigma, 1e-12);
    let mut supplied_heat = 0.0;
    for (cell, electric) in field.electric_field_v_m().iter().enumerate() {
        let heat = operator.conductors()[cell]
            .joule_heating_w_m3(*electric)
            .unwrap();
        supplied_heat += operator.mesh().face_measure(cell).unwrap()
            * (heat.electron_power_density_w_m3() + heat.heavy_power_density_w_m3());
    }
    fixture::relative(supplied_heat, field.power_w(), 1e-12);
}

#[test]
fn g1_two_temperature_atomic_ionization_changes_current_and_closes_electrical_work() {
    let initial = fixture::frozen_state(1e-4, 6000.0, 16000.0, 0.2);
    let forward = ElectronImpactRateV1::new(
        fixture::id("e-"),
        nist::hydrogen(fixture::id("H")),
        MaxwellianRateControlsV1 {
            relative_tolerance: 1e-10,
            maximum_subdivisions: 4096,
        },
        0.1,
    )
    .unwrap();
    let channel = AtomicIonizationV1::new(
        initial.thermodynamics().clone(),
        ReactionId::new("source-bound-hydrogen-impact").unwrap(),
        fixture::id("H+"),
        forward,
        AtomicStatisticalWeightsV1 { neutral: 2, ion: 1 },
        1e-12,
    )
    .unwrap();
    let before = material(&initial);
    let sigma_before = before.conductivity_s_m();
    assert!(before.admits_electric_field([2.0, 0.0]));
    let original = PlanarConduction::new(
        layered(0.5, 1.0, 0.2),
        vec![before; 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let field = original
        .solve(&ends(2.0), field_controls(), &mut || false)
        .unwrap();
    let resistance_before = original
        .resistance(&[2, 5], &[0, 3], field_controls(), &mut || false)
        .unwrap()
        .resistance_ohm();
    let mut powers = [0.0, 0.0];
    let mut volume = 0.0;
    for (cell, e) in field.electric_field_v_m().iter().enumerate() {
        let work = original.conductors()[cell].joule_heating_w_m3(*e).unwrap();
        let cell_volume = original.mesh().face_measure(cell).unwrap();
        powers[0] += cell_volume * work.heavy_power_density_w_m3();
        powers[1] += cell_volume * work.electron_power_density_w_m3();
        volume += cell_volume;
        fixture::relative(
            cell_volume * work.total_power_density_w_m3(),
            field.cell_power_w()[cell],
            1e-13,
        );
    }
    fixture::relative(volume, initial.volume().value(), 1e-14);
    fixture::relative(powers[0] + powers[1], field.power_w(), 1e-13);
    let mut controls =
        AtomicIonizationSolve::new(Temperature::new(2000.0), Temperature::new(20000.0));
    controls.relative_tolerance = 1e-7;
    controls.energy_tolerance_j = 1e-8;
    controls.initial_step_s = 1e-10;
    controls.maximum_step_s = 1e-9;
    let duration = 1e-8;
    let result = advance_atomic_ionization(
        &initial,
        &channel,
        Power::new(powers[1]),
        Power::new(powers[0]),
        Time::new(duration),
        controls,
        || false,
    )
    .unwrap();
    let next = result.state();
    let n0 = initial.amounts();
    let n1 = next.amounts();
    let extent = result.reaction_extent().value();
    assert!(
        extent > 1e-8,
        "the source must produce observable electron-impact growth"
    );
    fixture::relative(n1[1].value() - n0[1].value(), extent, 1e-10);
    fixture::relative(n1[2].value() - n0[2].value(), extent, 1e-10);
    fixture::relative(n1[0].value() + n1[1].value(), 1.0, 1e-12);
    fixture::relative(n1[1].value(), n1[2].value(), 1e-12);
    fixture::relative(next.mass().value(), initial.mass().value(), 1e-12);
    assert_eq!(next.volume(), initial.volume());
    assert!(next.electron_temperature().value() < initial.electron_temperature().value());
    assert_ne!(next.electron_temperature(), next.heavy_temperature());

    // Independent monatomic energy account: electron number changes as well as
    // Te. The source ionization gap enters chemical storage once, without an
    // additional subtraction of the NIST cross-section threshold.
    let electron_sensible = 1.5
        * R
        * (n1[2].value() * next.electron_temperature().value()
            - n0[2].value() * initial.electron_temperature().value());
    let heavy_sensible =
        1.5 * R * (next.heavy_temperature().value() - initial.heavy_temperature().value());
    let chemical_storage = extent * channel.ionization_energy_j_mol();
    let supplied = (powers[0] + powers[1]) * duration;
    assert!((electron_sensible + chemical_storage - powers[1] * duration).abs() < 2e-8);
    assert!((heavy_sensible - powers[0] * duration).abs() < 2e-8);
    assert!((electron_sensible + heavy_sensible + chemical_storage - supplied).abs() < 3e-8);
    let balance = result.balance();
    assert!((balance.electron_deposited_energy_j - powers[1] * duration).abs() < 1e-15);
    assert!((balance.heavy_deposited_energy_j - powers[0] * duration).abs() < 1e-15);
    fixture::relative(balance.chemical_storage_change_j, chemical_storage, 1e-12);
    assert!(balance.energy_residual_j.abs() <= controls.energy_tolerance_j);

    let after = material(next);
    assert!(after.equilibrium().is_none());
    assert_eq!(after.state().amounts(), n1);
    assert!(after.conductivity_s_m() > sigma_before);
    assert!(after.admits_electric_field([2.0, 0.0]));
    let sigma_after = after.conductivity_s_m();
    let updated = PlanarConduction::new(
        layered(0.5, 1.0, 0.2),
        vec![after; 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let new_field = updated
        .solve(&ends(2.0), field_controls(), &mut || false)
        .unwrap();
    fixture::relative(
        new_field.current_density_a_m2()[0][0],
        -2.0 * sigma_after,
        1e-12,
    );
    fixture::relative(new_field.power_w(), 0.8 * sigma_after, 1e-12);
    assert!(new_field.power_w() > field.power_w());
    let resistance_after = updated
        .resistance(&[2, 5], &[0, 3], field_controls(), &mut || false)
        .unwrap()
        .resistance_ohm();
    assert!(resistance_after < resistance_before);
}

#[test]
fn g1_photons_seed_a_neutral_gas_and_feed_existing_current_and_impact_owners() {
    use fs_gas::chemistry::photoionization::{
        PhotoionizationControls, PhotoionizationExposure, advance_monochromatic_photoionization,
    };
    use fs_qty::{Area, Energy};
    use fs_thermochem::kinetics::atomic_ionization::photoionization::{
        AtomicPhotoionizationV1, PhotonNumberFluxV1,
    };

    let initial = fixture::frozen_state(0.0, 6000.0, 16000.0, 0.2);
    let channel = AtomicIonizationV1::new(
        initial.thermodynamics().clone(),
        ReactionId::new("photo-seeded-hydrogen-impact").unwrap(),
        fixture::id("H+"),
        ElectronImpactRateV1::new(
            fixture::id("e-"), nist::hydrogen(fixture::id("H")),
            MaxwellianRateControlsV1 { relative_tolerance: 1e-10, maximum_subdivisions: 4096 },
            0.1,
        ).unwrap(),
        AtomicStatisticalWeightsV1 { neutral: 2, ion: 1 },
        1e-12,
    ).unwrap();
    let photon_source = AtomicPhotoionizationV1::verner_hydrogen(channel.clone()).unwrap();
    let photo = advance_monochromatic_photoionization(
        &initial, &photon_source,
        PhotoionizationExposure {
            beam: PhotonNumberFluxV1 {
                photon_energy: Energy::new(16.0 * fixture::E),
                number_flux_per_m2_s: 1e-3 * fixture::NA / 1e-7,
            },
            illuminated_area: Area::new(1.0), duration: Time::new(1e-7),
        },
        PhotoionizationControls { relative_charge_tolerance: 1e-12, energy_tolerance_j: 1e-7 },
        || false,
    ).unwrap();
    let seeded = photo.state();
    assert!(seeded.amounts()[2].value() > 0.0);
    fixture::relative(seeded.amounts()[2].value(), 1e-3, 1e-12);
    assert_eq!(seeded.amounts()[1], seeded.amounts()[2]);
    assert!(photo.balance().absorbed_photons <= photo.balance().incident_photons * (1.0 + 1e-14));
    assert_eq!(seeded.thermodynamics(), initial.thermodynamics());

    // Original transport sources are evaluated at the actual new composition
    // and temperatures; no arbitrary conductivity is attached to illumination.
    let conductor = material(seeded);
    assert_eq!(conductor.state().amounts(), seeded.amounts());
    let sigma = conductor.conductivity_s_m();
    assert!(sigma > 0.0);
    {
        use fs_math::c64::C64;
        let omega = 1e12;
        let transport = conductor.prepare_harmonic(|| false).unwrap();
        let response = transport.response(
            omega, [C64::from_re(0.001), C64::ZERO, C64::ZERO], || false,
        ).unwrap();
        let expected = conductor.carriers().iter().fold(C64::ZERO, |sum, carrier| {
            let q = carrier.charge_number() as f64 * fixture::E;
            sum + C64::from_re(carrier.number_density_per_m3() * q * q / carrier.mass_kg())
                / C64::new(carrier.relaxation().rate_per_s(), omega)
        });
        assert!((response.conductivity_s_m - expected).abs() <= 3e-12 * expected.abs());
        fixture::relative(response.mean_electrical_power_w_m3,
            response.mean_heating_w_m3.total_power_density_w_m3(), 3e-12);
    }
    let operator = PlanarConduction::new(
        layered(0.5, 1.0, 0.2), vec![conductor; 4], limits(), &mut || false,
    ).unwrap();
    let field = operator.solve(&ends(2.0), field_controls(), &mut || false).unwrap();
    fixture::relative(field.current_density_a_m2()[0][0], -2.0 * sigma, 1e-12);
    fixture::relative(field.power_w(), 0.8 * sigma, 1e-12);

    // A subsequent collision-only interval consumes the photo-created
    // electrons through the existing channel, with no artificial seed or LTE reset.
    let mut controls = AtomicIonizationSolve::new(
        Temperature::new(2000.0), Temperature::new(20000.0),
    );
    controls.relative_tolerance = 1e-7;
    controls.energy_tolerance_j = 1e-8;
    controls.initial_step_s = 1e-12;
    controls.maximum_step_s = 1e-11;
    let reacted = advance_atomic_ionization(
        seeded, &channel, Power::new(0.0), Power::new(0.0), Time::new(1e-11),
        controls, || false,
    ).unwrap();
    assert!(reacted.reaction_extent().value() > 0.0);
    assert_eq!(reacted.state().amounts()[1], reacted.state().amounts()[2]);
    let total_gain = reacted.state().internal_energy().value() - initial.internal_energy().value();
    assert!((total_gain - photo.balance().absorbed_energy_j).abs() < 2e-7);

    // Radiation subsequently removes carriers and chemical energy from the
    // same original gas; the existing material owner sees their actual loss.
    use fs_gas::chemistry::radiative_recombination::{
        RadiativeRecombinationSolve, advance_radiative_recombination,
    };
    use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::AtomicRadiativeRecombinationV1;
    let reciprocal = AtomicRadiativeRecombinationV1::new(
        photon_source.clone(),
        Area::new(6.4e-22),
        fs_matdb::Provenance {
            source: "test-declared global hydrogen cross-section envelope".into(),
            license: "analytic test envelope assumption".into(),
            artifact: None,
        },
    )
    .unwrap();
    let conductivity_before_capture = material(reacted.state()).conductivity_s_m();
    let capture = advance_radiative_recombination(
        reacted.state(),
        &reciprocal,
        Time::new(1e-5),
        RadiativeRecombinationSolve::new(Temperature::new(2000.0), Temperature::new(20000.0)),
        || false,
    )
    .unwrap();
    assert!(capture.recombination_extent().value() > 0.0);
    assert!(capture.state().amounts()[2].value() < reacted.state().amounts()[2].value());
    assert_eq!(capture.state().amounts()[1], capture.state().amounts()[2]);
    assert!(material(capture.state()).conductivity_s_m() < conductivity_before_capture);
    let retained = capture.state().internal_energy().value() - initial.internal_energy().value();
    assert!(
        (retained + capture.balance().emitted_radiation_energy_j
            - photo.balance().absorbed_energy_j)
            .abs()
            < 3e-7
    );
}
