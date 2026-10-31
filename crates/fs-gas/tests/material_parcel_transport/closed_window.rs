//! G0/G1/G3: original endpoint accounts, including sub-ULP absolute energies.

use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::diffusion::network::DiffusionNetworkError;
use transport::audit_closed_window;

#[test]
fn closed_window_uses_actual_endpoints_and_cancels_formation_references() {
    let run = |formation| {
        let (bank, mut op) = fixture(2.0 * MASS, formation, 0.0);
        let original = BTreeMap::from([
            (1, parcel(&bank, 0.7, 390.0, 1e5)),
            (2, parcel(&bank, 0.3, 410.0, 1e5)),
        ]);
        let mut current = original.clone();
        for _ in 0..8 {
            current = solve(&mut op, &current, &[link()], 0.1)
                .unwrap()
                .into_parts()
                .0;
        }
        let audit =
            audit_closed_window(&original, &current, controls().physical, || false).unwrap();
        // Both manufactured species have Cv=3R per mole. Independently recover
        // their actual sensible storage, irrespective of unequal molar masses.
        let sensible = |gases: &BTreeMap<u64, ClosedGasParcel>| {
            gases
                .values()
                .map(|p| {
                    (mass(p, "A") / MASS + mass(p, "B") / (2.0 * MASS))
                        * 3.0
                        * R
                        * p.state().state().temperature
                })
                .sum::<f64>()
        };
        close(
            audit.energy_change_j,
            sensible(&current) - sensible(&original),
            2e-10,
        );
        assert!(audit.entropy_change_j_k > 0.0);
        assert!(
            audit.maximum_species_mass_residual_kg
                <= controls().physical.exchange.mass_tolerance_kg
        );
        audit
    };
    let base = run(0.0);
    let shifted = run(1e14);
    close(base.energy_change_j, shifted.energy_change_j, 1e-10);
    close(base.entropy_change_j_k, shifted.entropy_change_j_k, 1e-12);
}

#[test]
fn closed_window_detects_heat_hidden_below_absolute_formation_energy_resolution() {
    let (bank, _) = fixture(MASS, 1e24, 0.0);
    let t0 = 400.0_f64;
    let t1 = t0 + 1e-8;
    let original = BTreeMap::from([(1, parcel(&bank, 0.5, t0, 1e5))]);
    let heated = BTreeMap::from([(1, parcel(&bank, 0.5, t1, 1e5 * (t1 / t0)))]);
    let old_absolute =
        original[&1].mass_kg() * original[&1].state().specific_internal_energy_j_kg();
    let new_absolute = heated[&1].mass_kg() * heated[&1].state().specific_internal_energy_j_kg();
    let expected = heated[&1].mass_kg() / MASS * 3.0 * R * (t1 - t0);
    assert!(expected > 0.0);
    // The naive answer is either zero or a formation-weighted rounding artifact.
    assert!(((new_absolute - old_absolute) - expected).abs() >= expected * 0.99);
    let mut allowance = controls().physical;
    allowance.exchange.energy_tolerance_j = 1.0;
    let audit = audit_closed_window(&original, &heated, allowance, || false).unwrap();
    close(audit.energy_change_j, expected, expected * 1e-10);
    allowance.exchange.energy_tolerance_j = expected * 0.5;
    assert!(matches!(
        audit_closed_window(&original, &heated, allowance, || false),
        Err(GasTransportError::Invalid(_))
    ));
}

#[test]
fn closed_window_handles_absent_species_and_types_actual_entropy_decrease() {
    let (bank, _) = fixture(MASS, 0.0, 0.0);
    let separated = BTreeMap::from([
        (1, parcel(&bank, 1.0, 400.0, 1e5)),
        (2, parcel(&bank, 0.0, 400.0, 1e5)),
    ]);
    let mixed = BTreeMap::from([
        (1, parcel(&bank, 0.5, 400.0, 1e5)),
        (2, parcel(&bank, 0.5, 400.0, 1e5)),
    ]);
    let audit = audit_closed_window(&separated, &mixed, controls().physical, || false).unwrap();
    let expected = 2.0 * separated[&1].mass_kg() / MASS * R * 2.0_f64.ln();
    close(audit.entropy_change_j_k, expected, 1e-12);
    close(audit.energy_change_j, 0.0, 1e-15);
    assert!(matches!(
        audit_closed_window(&mixed, &separated, controls().physical, || false),
        Err(GasTransportError::Diffusion(DiffusionNetworkError::EntropyDecrease {
            change_j_k,
            tolerance_j_k,
        })) if change_j_k < -tolerance_j_k
    ));
}

#[test]
fn closed_window_rejects_species_creation_source_changes_and_late_cancellation() {
    let (bank, _) = fixture(MASS, 0.0, 0.0);
    let original = BTreeMap::from([
        (1, parcel(&bank, 0.5, 400.0, 1e5)),
        (2, parcel(&bank, 0.5, 400.0, 1e5)),
    ]);
    let mut altered = original.clone();
    altered.insert(1, parcel(&bank, 0.6, 400.0, 1e5));
    assert!(matches!(
        audit_closed_window(&original, &altered, controls().physical, || false),
        Err(GasTransportError::Diffusion(
            DiffusionNetworkError::Invalid(_)
        ))
    ));
    let (different, _) = fixture(MASS, 1e5, 0.0);
    altered.insert(1, parcel(&different, 0.5, 400.0, 1e5));
    assert!(matches!(
        audit_closed_window(&original, &altered, controls().physical, || false),
        Err(GasTransportError::Diffusion(
            DiffusionNetworkError::Invalid(_)
        ))
    ));
    // Even identical source cards cannot silently switch the caloric evaluator.
    let reconciled = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            bank.thermodynamics()
                .components()
                .iter()
                .map(|source| {
                    source
                        .with_caloric_continuity(
                            Temperature::new(400.0),
                            fs_thermochem::MolarEnergyQuantityV1::new(0.0),
                            fs_thermochem::MolarThermalQuantityV1::new(0.0),
                        )
                        .unwrap()
                })
                .collect(),
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
        )
        .unwrap(),
        bank.transport_materials().to_vec(),
        None,
    )
    .unwrap();
    altered.insert(1, parcel(&reconciled, 0.5, 400.0, 1e5));
    assert!(matches!(
        audit_closed_window(&original, &altered, controls().physical, || false),
        Err(GasTransportError::Diffusion(
            DiffusionNetworkError::Invalid(_)
        ))
    ));
    let mut polls = 0;
    audit_closed_window(&original, &original, controls().physical, || {
        polls += 1;
        false
    })
    .unwrap();
    let mut cancelled_polls = 0;
    assert_eq!(
        audit_closed_window(&original, &original, controls().physical, || {
            cancelled_polls += 1;
            cancelled_polls == polls
        }),
        Err(GasTransportError::Cancelled)
    );
    assert_eq!(original[&1].state().state().temperature, 400.0);
}

#[test]
fn porous_counterdiffusion_uses_the_same_reference_invariant_window_account() {
    use fs_material::porous::permeability::{
        DarcyGasBarrier, IntrinsicPermeability, PorousGasLayer, PorousGasSupport,
    };
    use fs_material::thermochemical_gas::energy::parcel::transfer::species::diffusion::{
        BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY, BinaryDiffusionMaterial,
        dusty_gas::{DustyGasFaceSolve, DustyGasMedium, DustyGasPores},
    };
    use transport::pore_domain::{PoreCellGeometry, PoreDomain, PoreDomainSolve, PoreFace};
    let medium = |length| {
        DustyGasMedium::try_new(
            DarcyGasBarrier::try_new(
                vec![
                    PorousGasLayer::try_new(
                        IntrinsicPermeability::isotropic(0.0).unwrap(),
                        length,
                        0.0,
                        "manufactured closed diffusion window",
                    )
                    .unwrap(),
                ],
                [1.0, 0.0, 0.0],
                PorousGasSupport {
                    temperature_k: (200.0, 1000.0),
                    pressure_pa: (100.0, 1e8),
                    maximum_superficial_speed_m_s: 100.0,
                },
            )
            .unwrap(),
            DustyGasPores::try_new(1e-6, 0.4, 2.0, "manufactured geometry").unwrap(),
        )
        .unwrap()
    };
    let coefficient = card(
        "pore molecular diffusion",
        [scalar(
            BINARY_DIFFUSIVITY_PROPERTY,
            BINARY_DIFFUSIVITY_DIMS,
            DIFFUSIVITY,
        )],
    );
    let context = point(400.0, 1e5);
    let pairs = [BinaryDiffusionMaterial {
        species: [id("A"), id("B")],
        material: &coefficient,
        context: &context,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }];
    let throat = medium(0.001);
    let face = PoreFace {
        medium: &throat,
        area_m2: 0.001,
        pairs: &pairs,
        gas_context: &context,
        solve: DustyGasFaceSolve {
            maximum_temperature_difference_k: 1e-6,
            maximum_relative_pressure_difference: 1e-6,
            gradient_tolerance_per_m: 1e-7,
        },
    };
    let run = |formation| {
        let (bank, _) = fixture(MASS, formation, 0.0);
        let initial = PoreDomain::try_new(
            0.1,
            vec![PoreCellGeometry::try_new(medium(0.025)).unwrap()],
            parcel(&bank, 0.8, 400.0, 1e5),
            vec![parcel(&bank, 0.2, 400.0, 1e5)],
            parcel(&bank, 0.2, 400.0, 1e5),
        )
        .unwrap();
        initial
            .step_candidate(
                1.0,
                &[face, face],
                &[endpoint(); 3],
                PoreDomainSolve {
                    maximum_step_s: 0.1,
                    minimum_step_s: 1e-6,
                    maximum_substeps: 100,
                    maximum_attempts: 200,
                    maximum_throat_volume_fraction: 0.01,
                    transport: controls().physical,
                },
                || false,
            )
            .unwrap()
    };
    let base = run(0.0);
    let shifted = run(1e14);
    assert!(
        base.step()
            .faces
            .iter()
            .any(|f| f.mass_left_to_right_kg > 0.0)
    );
    close(
        base.step().total_energy_residual_j,
        shifted.step().total_energy_residual_j,
        1e-10,
    );
    close(
        base.step().entropy_change_j_k,
        shifted.step().entropy_change_j_k,
        1e-12,
    );
    close(
        mass(base.domain().pore_gas(0).unwrap(), "A"),
        mass(shifted.domain().pore_gas(0).unwrap(), "A"),
        1e-18,
    );
}
