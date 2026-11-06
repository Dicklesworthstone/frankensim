//! G0/G3/G4: one continuous material law on the native superfluid heat timeline.
//! Authored Cp(T)=8T is an analytic fixture, not a helium calibration.

use fs_blake3::ContentHash;
use fs_couple::superfluid::{FrictionHeatControls, advance_closed_friction_heat};
use fs_duct::superfluid::{
    LinearizationLimits, MutualFriction, TwoFluidChannel, VortexTransportBoundaries,
    VortexTransportBoundary, VortexTransportControls,
};
use fs_material::superfluid::caloric::{
    IsobaricSuperfluidLaw, SuperfluidCaloricKnot, SuperfluidCaloricReference,
};
use fs_material::superfluid::vinen::{VinenParameters, VinenStepControls, VinenTangle};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, ThermalStorageBasis, TimedEvent,
};

const N: usize = 4;
const MASS: f64 = 0.01;
const FRICTION: MutualFriction = MutualFriction::EvolvingVinenTangle;

fn law(upper_k: f64) -> IsobaricSuperfluidLaw {
    IsobaricSuperfluidLaw::try_new(
        ContentHash([0x72; 32]),
        SuperfluidCaloricReference {
            pressure_pa: 1000.0,
            density_at_min_temperature_kg_m3: 4.0,
            specific_enthalpy_at_min_temperature_j_kg: 0.0,
            specific_entropy_at_min_temperature_j_kg_k: 4.0,
        },
        [1.0, 1.7, upper_k]
            .into_iter()
            .map(|temperature_k| SuperfluidCaloricKnot {
                temperature_k,
                specific_heat_cp_j_kg_k: 8.0 * temperature_k,
                superfluid_fraction: 0.95 - 0.2 * temperature_k,
                adiabatic_sound_speed_m_s: 10.0,
                thermal_expansion_per_k: 0.0,
                normal_viscosity_pa_s: 0.03,
            })
            .collect(),
    )
    .unwrap()
}

fn controls() -> FrictionHeatControls {
    FrictionHeatControls {
        thermal: StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-12,
        },
        cell_heat_tolerance_j: 1e-14,
        heat_energy_tolerance_j: 1e-13,
        combined_energy_tolerance_j: 1e-12,
    }
}

fn ids() -> [CellId; N] {
    [CellId(10), CellId(11), CellId(12), CellId(13)]
}

fn channel(law: &IsobaricSuperfluidLaw, seed: f64, pulse: bool) -> TwoFluidChannel {
    let material = *law.state_at_temperature_k(2.0).unwrap().two_fluid();
    let mut channel = TwoFluidChannel::try_new(
        material,
        1.0,
        0.01,
        N,
        LinearizationLimits {
            max_pressure_perturbation_pa: 1.0,
            max_temperature_perturbation_k: 0.02,
            max_component_speed_m_s: 0.1,
        },
    )
    .unwrap();
    if pulse {
        let relative = [0.0, 0.025, 0.04, 0.025, 0.0];
        let normal = relative.map(|w| material.superfluid_fraction() * w);
        let superfluid = relative.map(|w| -material.normal_fraction() * w);
        channel
            .initialize(&[0.0; N], &[0.0; N], &normal, &superfluid, || false)
            .unwrap();
    }
    let model = VinenTangle::try_new(
        material,
        VinenParameters {
            circulation_quantum_m2_s: 1e-5,
            dissipative_b: 1.0,
            production_coefficient: 0.02,
            decay_coefficient_m2_s: 1e-6,
            line_energy_j_m: 1e-8,
        },
        VinenStepControls {
            maximum_fractional_step: 0.04,
            maximum_substeps: 512,
        },
    )
    .unwrap()
    .with_counterflow_onset_m_s(0.005)
    .unwrap();
    channel
        .configure_vinen_tangle(model, &[seed; N - 1], || false)
        .unwrap();
    channel
}

fn thermal(law: &IsobaricSuperfluidLaw, source_w: f64) -> DepositionSimulation<'_> {
    let initial_h = law
        .state_at_temperature_k(2.0)
        .unwrap()
        .specific_enthalpy_j_kg();
    let events = ids()
        .into_iter()
        .enumerate()
        .map(|(i, id)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: i as u64,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: MASS,
                specific_h_j_kg: initial_h,
                geometry: CellGeometry {
                    start_m: [i as f64 / N as f64, 0.0, 0.0],
                    end_m: [(i + 1) as f64 / N as f64, 0.0, 0.0],
                    volume_m3: 0.01 / N as f64,
                    surface_area_m2: 0.01,
                },
                boundary: Boundary {
                    ambient_k: 2.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w,
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: N,
            max_links: 0,
            max_events: N,
            max_substeps: 128,
        },
    )
    .unwrap()
    .advance_to(0.0, controls().thermal, || Ok(()))
    .unwrap()
    .0
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tol={tolerance:e}"
    );
}

fn snapshot(channel: &TwoFluidChannel, thermal: &DepositionSimulation<'_>) -> Vec<u64> {
    let mut values = vec![
        channel.time_s(),
        channel.represented_energy_j().unwrap(),
        channel.cumulative_mutual_friction_heat_j(),
        thermal.time_s(),
        thermal.energy_j(),
    ];
    for i in 0..=N {
        let f = channel.face(i).unwrap();
        values.extend([f.normal_velocity_m_s, f.superfluid_velocity_m_s]);
        if let Some(l) = channel.vortex_line_density_per_m2(i) {
            values.push(l);
        }
        if i < N {
            let c = channel.cell(i).unwrap();
            values.extend([c.pressure_perturbation_pa, c.temperature_perturbation_k]);
            values.push(thermal.cells()[&ids()[i]].specific_h_j_kg);
        }
    }
    values.extend(channel.last_dissipated_heat_j());
    values.into_iter().map(f64::to_bits).collect()
}

#[test]
fn native_heating_and_cooling_use_the_same_integrated_caloric_law() {
    let law = law(3.0);
    assert_eq!(
        ThermalLaw::storage_basis(&law),
        ThermalStorageBasis::Enthalpy
    );
    let slope_bound = ThermalLaw::maximum_temperature_slope(&law);
    assert!(slope_bound >= 0.125 && slope_bound < 0.125 + 1e-14);
    assert!(ThermalLaw::caloric_identity(&law).is_some());
    for source_w in [-0.1, 0.1] {
        let mut channel = channel(&law, 4e5, true);
        let mut thermal = thermal(&law, source_w);
        let initial = channel.represented_energy_j().unwrap() + thermal.energy_j();
        let initial_h = thermal.cells()[&ids()[0]].specific_h_j_kg;
        let mut external = 0.0;
        let mut friction_heat = 0.0;
        for _ in 0..12 {
            let report = advance_closed_friction_heat(
                &mut channel,
                &mut thermal,
                &ids(),
                0.005,
                FRICTION,
                controls(),
                || false,
            )
            .unwrap();
            external += report.thermal_march.external_heat_j;
            friction_heat += report.channel.dissipated_heat_j();
            close(report.transfer_residual_j, 0.0, 1e-13);
            close(report.combined_storage_residual_j, 0.0, 1e-12);
        }
        assert!(friction_heat > 0.0);
        close(
            channel.represented_energy_j().unwrap() + thermal.energy_j() - initial,
            external,
            1e-12,
        );
        for id in ids() {
            let h = thermal.cells()[&id].specific_h_j_kg;
            let temperature = thermal.temperature_k(id).unwrap();
            // Independent exact Cp=8T integral: h=4(T^2-1), delta s=8 delta T.
            close(temperature, (1.0 + h / 4.0).sqrt(), 2e-14);
            let entropy = ThermalLaw::specific_entropy_change_j_kg_k(&law, initial_h, h)
                .unwrap()
                .unwrap();
            close(entropy, 8.0 * (temperature - 2.0), 2e-13);
            assert_eq!(temperature > 2.0, source_w > 0.0);
        }
        // The finite thermal owner changes T; the wave linearization remains
        // the explicitly selected initial reference, with no hidden remap.
        assert_eq!(channel.material().properties().temperature_k, 2.0);
    }
}

#[test]
fn finite_caloric_domain_refusal_rolls_back_tangle_waves_heat_and_clock() {
    let law = law(2.000000001);
    let mut channel = channel(&law, 4e5, true);
    let mut thermal = thermal(&law, 0.0);
    let before = snapshot(&channel, &thermal);
    assert!(
        advance_closed_friction_heat(
            &mut channel,
            &mut thermal,
            &ids(),
            0.005,
            FRICTION,
            controls(),
            || false,
        )
        .is_err()
    );
    assert_eq!(snapshot(&channel, &thermal), before);
    let report = advance_closed_friction_heat(
        &mut channel,
        &mut thermal,
        &ids(),
        1e-8,
        FRICTION,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(channel.time_s(), 1e-8);
    assert_eq!(thermal.time_s(), 1e-8);
    assert!(report.channel.dissipated_heat_j() > 0.0);
    close(report.combined_storage_residual_j, 0.0, 1e-12);
}

#[test]
fn imported_vortex_energy_is_split_between_line_storage_and_native_caloric_heat_once() {
    let law = law(3.0);
    let mut channel = channel(&law, 0.0, false);
    let mut thermal = thermal(&law, 0.0);
    channel
        .configure_vortex_line_transport(
            &[0.1; N],
            VortexTransportBoundaries {
                left: VortexTransportBoundary::Reservoir {
                    incoming_line_density_per_m2: 1e8,
                },
                right: VortexTransportBoundary::Reservoir {
                    incoming_line_density_per_m2: 0.0,
                },
            },
            VortexTransportControls {
                maximum_outgoing_courant: 0.4,
                maximum_substeps: 64,
            },
            || false,
        )
        .unwrap();
    let initial = thermal.energy_j();
    let mut imported = 0.0;
    let mut dissipated = 0.0;
    for _ in 0..8 {
        let report = advance_closed_friction_heat(
            &mut channel,
            &mut thermal,
            &ids(),
            0.005,
            FRICTION,
            controls(),
            || false,
        )
        .unwrap();
        imported += report
            .channel
            .vortex_transport_boundary_energy_j
            .iter()
            .sum::<f64>();
        dissipated += report.channel.dissipated_heat_j();
        assert!(report.channel.vortex_transport_boundary_line_m[0] > 0.0);
        assert!(report.channel.vortex_transport_boundary_line_m[1] <= 0.0);
        close(report.combined_storage_residual_j, 0.0, 1e-12);
    }
    assert!(imported > 0.0);
    assert!(dissipated > 0.0 && dissipated < imported);
    assert!(channel.vortex_energy_j().unwrap() > 0.0);
    close(channel.wave_energy_j().unwrap(), 0.0, 1e-30);
    close(thermal.energy_j() - initial, dissipated, 1e-13);
    close(
        channel.vortex_energy_j().unwrap() + dissipated,
        imported,
        1e-12,
    );
    assert!(thermal.temperature_k(ids()[0]).unwrap() > 2.0);
}
