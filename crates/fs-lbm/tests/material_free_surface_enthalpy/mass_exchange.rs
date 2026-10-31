//! G0/G1/G3/G4: phase transfer changes actual native volume and storage.

use super::*;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceMassExchange, FreeSurfaceMassExchangeControls, MaterialFreeSurfaceEnthalpy3,
};

#[allow(dead_code)]
#[path = "../nsf3_thermochemical/source.rs"]
mod gas_source;

const INTERFACE: [usize; 3] = [3, 3, 3];

fn controls() -> FreeSurfaceMassExchangeControls {
    FreeSurfaceMassExchangeControls {
        mass_tolerance_kg: 1e-18,
        energy_tolerance_j: 1e-12,
        momentum_tolerance_kg_m_s: 1e-18,
    }
}

fn with_liquid<R>(
    datum: f64,
    f: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    with_liquid_at(datum, 400.0, f)
}

fn with_liquid_at<R>(
    datum: f64,
    temperature_k: f64,
    f: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let card = card(0.02, 1000.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point(temperature_k),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        &card,
        &point(temperature_k),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Liquid,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: datum,
                temperature_k: 300.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: datum + 600_000.0,
                temperature_k: 600.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
        ],
    )
    .unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut cells = initial(false);
        cells
            .iter_mut()
            .find(|c| c.coordinates == INTERFACE)
            .unwrap()
            .fill_fraction = 0.4;
        let flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), &cells, cx)
            .unwrap();
        let mut liquid = flow
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                |_| datum + (temperature_k - 300.0) * 2000.0,
                &gate,
            )
            .unwrap();
        f(&mut liquid, &gate)
    })
}

#[test]
fn g1_shared_saturation_owner_transfers_its_film_once_into_native_liquid() {
    use fs_material::phase::condensation::{WaterCondensationModel, WetGasParcel};
    use fs_material::thermochemical_gas::energy::IsochoricGasParcel;

    let temperature = 310.0;
    let source = gas_source::moist_model(0.08, 2e6);
    let water = source.water_species().unwrap().clone();
    let model = WaterCondensationModel::new(source.clone(), water, temperature).unwrap();
    let gas =
        IsochoricGasParcel::new(source, gas_source::point(temperature, 101325.0), 1e-6).unwrap();
    let inlet = WetGasParcel::new(model, gas, 0.0).unwrap();
    let phase = inlet.at_temperature(temperature, || false).unwrap();
    assert!(phase.condensed_mass_kg > 0.0);
    assert!(phase.total_internal_energy_change_j < 0.0);
    assert!(phase.saturation_residual_pa.abs() < 1e-8);
    let sat = phase.parcel.liquid_saturation();

    with_liquid_at(0.0, temperature, |liquid, gate| {
        let old_mass = liquid.flow().ledger_mass_kg().unwrap();
        let old_energy = liquid.total_enthalpy_j().unwrap();
        let native_h = liquid.specific_enthalpy_at(INTERFACE).unwrap();
        // Both charts describe this explicitly imposed transfer temperature.
        // Map the original IAPWS enthalpy to the native synthetic caloric datum;
        // NASA's chemical formation reference is never a liquid temperature.
        let native_offset = native_h - sat.liquid_specific_enthalpy().value();
        let report = liquid
            .exchange_interface_mass(
                &[FreeSurfaceMassExchange {
                    coordinates: INTERFACE,
                    mass_kg: phase.condensed_mass_kg,
                    specific_enthalpy_j_kg: sat.liquid_specific_enthalpy().value() + native_offset,
                    // This is an isothermal separator. Its original phase-owner heat
                    // duty already removes latent energy, so none is added again here.
                    heat_j: 0.0,
                }],
                controls(),
                gate,
            )
            .unwrap();
        let new_mass = liquid.flow().ledger_mass_kg().unwrap();
        near(
            new_mass + phase.parcel.gas().mass_kg(),
            old_mass + inlet.gas().mass_kg(),
            2e-13,
        );
        near(
            liquid.specific_enthalpy_at(INTERFACE).unwrap(),
            native_h,
            2e-13,
        );
        near(
            liquid.total_enthalpy_j().unwrap(),
            old_energy + report.mass_change_kg * native_h,
            2e-13,
        );

        // Export ONLY the proposed gas. The phase proposal's film has moved
        // into native liquid; retaining it as a second store would double mass.
        // Reconcile U (fixed gas volume) versus native H and both source datums.
        let alignment = inlet.model().liquid_energy_offset_j_kg() - native_offset;
        let pressure_work =
            report.mass_change_kg * sat.pressure().value() / sat.liquid_density().value();
        let native_aligned_u =
            report.enthalpy_change_j + report.mass_change_kg * alignment - pressure_work;
        assert!(
            (phase.gas_internal_energy_change_j + native_aligned_u
                - phase.total_internal_energy_change_j)
                .abs()
                < 1e-11
        );
        assert_eq!(report.heat_j, 0.0);
        assert!(phase.latent_enthalpy_j > 0.0);
        let continued = liquid.step(gate).unwrap();
        near(
            continued.mass_after_kg + phase.parcel.gas().mass_kg(),
            old_mass + inlet.gas().mass_kg(),
            2e-12,
        );
        near(
            continued.enthalpy_after_j,
            old_energy + report.carried_enthalpy_j,
            2e-11,
        );
    });
}

#[test]
fn g1_condensate_adds_volume_enthalpy_and_carried_momentum_before_native_flow() {
    with_liquid(0.0, |liquid, gate| {
        let mass = liquid.flow().ledger_mass_kg().unwrap();
        let energy = liquid.total_enthalpy_j().unwrap();
        let before = liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
        let populations = liquid.flow().lattice().grid().state_bits();
        let report = liquid
            .exchange_interface_mass(
                &[FreeSurfaceMassExchange {
                    coordinates: INTERFACE,
                    mass_kg: 1e-7,
                    specific_enthalpy_j_kg: 100_000.0,
                    heat_j: 0.02,
                }],
                controls(),
                gate,
            )
            .unwrap();
        let after = liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
        near(after.fill_fraction, 0.5, 2e-14);
        near(report.mass_change_kg, 1e-7, 2e-14);
        near(liquid.flow().ledger_mass_kg().unwrap(), mass + 1e-7, 2e-14);
        near(
            liquid.specific_enthalpy_at(INTERFACE).unwrap(),
            220_000.0,
            2e-14,
        );
        near(report.enthalpy_change_j, 0.03, 2e-14);
        near(report.carried_enthalpy_j, 0.01, 2e-14);
        near(liquid.total_enthalpy_j().unwrap(), energy + 0.03, 2e-14);
        assert_eq!(before.velocity_m_s, after.velocity_m_s);
        assert_eq!(populations, liquid.flow().lattice().grid().state_bits());
        for axis in 0..3 {
            near(
                report.momentum_change_kg_m_s[axis],
                1e-7 * before.velocity_m_s[axis],
                2e-14,
            );
        }
        assert_eq!(liquid.physical_time_s(), 0.0);
        assert_eq!(report.physical_time_s, 0.0);
        let endpoint = liquid.material_at(INTERFACE).unwrap();
        near(endpoint.dynamic_viscosity_pa_s(), 0.056, 2e-14);
        let advanced = liquid.step(gate).unwrap();
        near(advanced.mass_before_kg, mass + 1e-7, 2e-14);
        near(advanced.mass_after_kg, mass + 1e-7, 2e-12);
        near(advanced.enthalpy_after_j, energy + 0.03, 2e-11);
        assert_eq!(liquid.physical_time_s(), 0.001);
    });
}

#[test]
fn g1_evaporation_removes_actual_liquid_mass_and_applies_latent_cooling_once() {
    with_liquid(0.0, |liquid, gate| {
        let before = liquid.total_enthalpy_j().unwrap();
        let report = liquid
            .exchange_interface_mass(
                &[FreeSurfaceMassExchange {
                    coordinates: INTERFACE,
                    mass_kg: -1e-8,
                    specific_enthalpy_j_kg: 200_000.0,
                    heat_j: -0.02,
                }],
                controls(),
                gate,
            )
            .unwrap();
        near(
            liquid
                .flow()
                .sample_at_si(3, 3, 3)
                .unwrap()
                .unwrap()
                .fill_fraction,
            0.39,
            2e-14,
        );
        near(report.mass_change_kg, -1e-8, 2e-14);
        near(report.carried_enthalpy_j, -0.002, 2e-14);
        near(report.enthalpy_change_j, -0.022, 2e-14);
        near(liquid.total_enthalpy_j().unwrap(), before - 0.022, 2e-14);
        near(
            liquid.specific_enthalpy_at(INTERFACE).unwrap(),
            200_000.0 - 0.02 / 3.9e-7,
            2e-14,
        );
    });
}

#[test]
fn g3_common_caloric_reference_does_not_change_transferred_temperature_or_fill() {
    let run = |datum| {
        with_liquid(datum, |liquid, gate| {
            let report = liquid
                .exchange_interface_mass(
                    &[FreeSurfaceMassExchange {
                        coordinates: INTERFACE,
                        mass_kg: 1e-7,
                        specific_enthalpy_j_kg: datum + 100_000.0,
                        heat_j: 0.02,
                    }],
                    controls(),
                    gate,
                )
                .unwrap();
            (
                liquid.specific_enthalpy_at(INTERFACE).unwrap() - datum,
                liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap(),
                report,
            )
        })
    };
    let base = run(0.0);
    let shifted = run(1e9);
    near(base.0, shifted.0, 1e-12);
    assert_eq!(base.1, shifted.1);
    near(
        shifted.2.carried_enthalpy_j - base.2.carried_enthalpy_j,
        shifted.2.mass_change_kg * 1e9,
        2e-14,
    );
    assert!(shifted.2.absolute_energy_defect_j <= controls().energy_tolerance_j);
}

#[test]
fn g4_phase_transfer_refusals_keep_native_fill_heat_sources_and_clock() {
    with_liquid(0.0, |liquid, gate| {
        let sample = liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
        let mass = liquid.flow().ledger_mass_kg().unwrap();
        let energy = liquid.total_enthalpy_j().unwrap();
        let material = liquid.material_at(INTERFACE).unwrap().clone();
        let bits = liquid.flow().lattice().grid().state_bits();
        let transfer = FreeSurfaceMassExchange {
            coordinates: INTERFACE,
            mass_kg: 1e-8,
            specific_enthalpy_j_kg: 200_000.0,
            heat_j: 0.0,
        };
        for bad in [
            FreeSurfaceMassExchange {
                mass_kg: 1e-6,
                ..transfer
            },
            FreeSurfaceMassExchange {
                mass_kg: -4e-7,
                ..transfer
            },
            FreeSurfaceMassExchange {
                mass_kg: 1e-30,
                ..transfer
            },
            FreeSurfaceMassExchange {
                heat_j: 0.3,
                ..transfer
            },
            FreeSurfaceMassExchange {
                coordinates: [2, 3, 3],
                ..transfer
            },
        ] {
            assert!(
                liquid
                    .exchange_interface_mass(&[bad], controls(), gate)
                    .is_err()
            );
        }
        // The first cell's private candidate cannot leak when the second fails.
        assert!(
            liquid
                .exchange_interface_mass(
                    &[
                        transfer,
                        FreeSurfaceMassExchange {
                            coordinates: [3, 3, 4],
                            mass_kg: 1e-6,
                            ..transfer
                        }
                    ],
                    controls(),
                    gate
                )
                .is_err()
        );
        assert!(
            liquid
                .exchange_interface_mass(&[transfer, transfer], controls(), gate)
                .is_err()
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            liquid
                .exchange_interface_mass(&[transfer], controls(), &cancelled)
                .is_err()
        );
        assert_eq!(
            liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap(),
            sample
        );
        assert_eq!(liquid.flow().ledger_mass_kg().unwrap(), mass);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), energy);
        assert_eq!(liquid.material_at(INTERFACE).unwrap(), &material);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}
