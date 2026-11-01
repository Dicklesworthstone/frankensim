//! Separate visible illumination and infrared exchange share native enthalpy.

#[path = "joint/refinement.rs"]
mod refinement;

use super::*;
use fs_blake3::ContentHash;
use fs_render::volumes::heating::exchange::feedback::{
    RadiativeFeedbackBudget, TemperatureRadiativeExchange,
};
use fs_render::volumes::heating::exchange::{RadiativeBath, RadiativeExchangeTolerance};
use fs_render::volumes::heating::feedback::TemperatureAbsorptionTable;
use fs_render::volumes::heating::illumination::{
    IlluminationBudget, IncidentIlluminationSource, TemperatureIlluminationHeating,
};
use fs_render::volumes::heating::joint::JointOpticalHeating;
use fs_render::volumes::heating::scattering::ScatteringHeatingTolerance;
use fs_render::volumes::heating::scattering_feedback::IncidentScatteringBeam;
use fs_render::volumes::heating::temperature_field::{
    TemperatureOpticalField, TemperatureOpticalFieldBudget, TemperatureScatteringTable,
    ThermalScatteringCell,
};
use fs_render::volumes::transport::scattering::deposition::ScatteringDepositionConfig;
use fs_render::volumes::transport::scattering::exchange::RadiativeExchangeConfig;
use fs_thermal::enthalpy::{DepositionObserver, ThermalLink};
use std::sync::Arc;

struct JointCaloric;
impl ThermalLaw for JointCaloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = h / 100.0;
        if t.is_finite() && (100.0..=1000.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.01
    }
}
static JOINT_CALORIC: JointCaloric = JointCaloric;

fn joint_cell(id: u64, temperature_k: f64) -> DepositedCell {
    let mut result = cell(id, 1.0, 1.0, 0.0);
    result.specific_h_j_kg = 100.0 * temperature_k;
    result.geometry.start_m = [(id - 11) as f64, 0.5, 0.5];
    result.geometry.end_m = [(id - 10) as f64, 0.5, 0.5];
    result
}

#[allow(clippy::too_many_arguments)]
fn joint_field<'g, const N: usize>(
    cx: &Cx<'_>,
    density: &'g [f64],
    thermal: &DepositionSimulation<'_>,
    wavelengths: [f64; N],
    support: (f64, f64),
    absorption: [f64; 2],
    scattering: f64,
    tag: u8,
) -> Arc<TemperatureOpticalField<'g, N>> {
    let a: Vec<_> = [100.0, 1000.0]
        .into_iter()
        .zip(absorption)
        .map(|(t, k)| {
            (
                t,
                MassAbsorptionSpectrum::source_table(
                    &[(support.0, k), (support.1, k)],
                    ContentHash([tag; 32]),
                )
                .unwrap(),
            )
        })
        .collect();
    let s = MassScatteringSpectrum::source_table(
        &[(support.0, scattering), (support.1, scattering)],
        ContentHash([tag + 1; 32]),
    )
    .unwrap();
    let bindings: Vec<_> = (0..density.len())
        .map(|i| ThermalScatteringCell {
            grid_cell: i,
            thermal_cell: CellId(11 + i as u64),
            absorption_table: 0,
            scattering_table: 0,
            phase: 0,
        })
        .collect();
    Arc::new(
        TemperatureOpticalField::try_new(
            grid(cx, density, MassAbsorptionSpectrum::constant(0.0).unwrap()),
            wavelengths,
            &[TemperatureAbsorptionTable::try_new(ContentHash([tag + 2; 32]), &a, 2).unwrap()],
            &[TemperatureScatteringTable::try_new(
                ContentHash([tag + 3; 32]),
                &[(100.0, s), (1000.0, s)],
                2,
            )
            .unwrap()],
            &[PhaseFunction::Isotropic],
            &bindings,
            thermal,
            TemperatureOpticalFieldBudget {
                max_grid_cells: 8,
                max_tables: 2,
                max_total_temperature_knots: 4,
            },
            || Ok(()),
        )
        .unwrap(),
    )
}

fn joint_fixture<'g>(
    cx: &Cx<'_>,
    density: &'g [f64],
    thermal: &DepositionSimulation<'_>,
    incident_power_w: f64,
    scattering: f64,
    temperature_dependent: bool,
) -> JointOpticalHeating<'g, 2, 1> {
    let illumination_field = joint_field(
        cx,
        density,
        thermal,
        [450.0, 650.0],
        (400.0, 700.0),
        if temperature_dependent {
            [0.2, 0.8]
        } else {
            [0.4; 2]
        },
        scattering,
        181,
    );
    let radiation_field = joint_field(
        cx,
        density,
        thermal,
        [10_000.0],
        (8000.0, 12_000.0),
        if temperature_dependent {
            [0.1, 0.5]
        } else {
            [0.25; 2]
        },
        scattering,
        191,
    );
    let sources = [
        IncidentIlluminationSource {
            source_id: 17,
            beam: IncidentScatteringBeam {
                origin: Point3::new(-1.0, 0.5, 0.5),
                unit_direction: Vec3::new(1.0, 0.0, 0.0),
                wavelengths_nm: [450.0, 650.0],
                incident_power_w: [incident_power_w * 0.45, incident_power_w * 0.15],
            },
        },
        IncidentIlluminationSource {
            source_id: 23,
            beam: IncidentScatteringBeam {
                origin: Point3::new(density.len() as f64 + 1.0, 0.5, 0.5),
                unit_direction: Vec3::new(-1.0, 0.0, 0.0),
                wavelengths_nm: [450.0, 650.0],
                incident_power_w: [incident_power_w * 0.1, incident_power_w * 0.3],
            },
        },
    ];
    let illumination = TemperatureIlluminationHeating::from_field(
        illumination_field,
        &sources,
        thermal,
        IlluminationBudget {
            max_sources: 2,
            max_source_evaluations: 5000,
            max_traced_segments: 10_000_000,
            max_retained_deposition_cells: 16,
        },
        ScatteringDepositionConfig {
            histories_per_lane: 16,
            max_collisions_per_history: 96,
            max_total_segments: 100_000,
            max_deposition_cells: 8,
            seed: 74,
            source_id: 31,
        },
        ScatteringHeatingTolerance {
            power_roundoff_w: 1e-9,
            max_unresolved_w: 1e-7,
            max_total_standard_error_w: 1000.0,
            max_cell_standard_error_w: 1000.0,
        },
        || Ok(()),
    )
    .unwrap();
    let radiation = TemperatureRadiativeExchange::try_new(
        radiation_field,
        [20.0],
        RadiativeBath::Blackbody {
            temperature_k: 300.0,
        },
        thermal,
        RadiativeFeedbackBudget {
            max_source_evaluations: 5000,
            max_traced_segments: 10_000_000,
            max_thermal_pairs: 16,
        },
        RadiativeExchangeConfig {
            histories_per_emitter_lane: 16,
            max_collisions_per_history: 96,
            max_total_segments: 100_000,
            max_material_cells: 8,
            max_emitter_cells: 8,
            max_pairs: 16,
            seed: 75,
            source_id: 32,
        },
        RadiativeExchangeTolerance {
            coefficient_roundoff_m2_sr: 1e-9,
            power_roundoff_w: 1e-9,
            max_cell_standard_error_w: 1000.0,
            max_total_standard_error_w: 1000.0,
            max_unresolved_exchange_w: 1e-6,
            max_unresolved_standard_error_w: 1e-6,
        },
        || Ok(()),
    )
    .unwrap();
    JointOpticalHeating::try_new(illumination, radiation, thermal, || Ok(())).unwrap()
}

#[test]
fn different_visible_and_infrared_packets_drive_one_conservative_thermal_trajectory() {
    with_cx(|cx| {
        let original = thermal(
            vec![joint_cell(11, 500.0), joint_cell(12, 350.0)],
            &JOINT_CALORIC,
            vec![TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 1,
                    a: CellId(11),
                    b: CellId(12),
                    area_m2: 0.2,
                    conductance_w_k: 0.3,
                }),
            }],
        );
        let joint = joint_fixture(cx, &[1.0; 2], &original, 80.0, 0.0, true);
        assert_eq!(
            joint.illumination().optical_field().wavelengths_nm(),
            &[450.0, 650.0]
        );
        assert_eq!(joint.radiation().field().wavelengths_nm(), &[10_000.0]);
        let (final_state, accepted, audit) = original
            .advance_observed(&joint, 0.25, config(0.01), || Ok(()))
            .unwrap();
        assert_eq!(accepted.time_s(), final_state.time_s());
        assert_eq!(accepted.radiation().time_s(), final_state.time_s());
        assert!(accepted.illumination().energy().thermal_input_j > 0.0);
        assert!(accepted.radiation().energy().thermal_external_j < 0.0);
        assert!(accepted.radiation().operator_builds() > 1);
        let optical_heat = accepted.illumination().energy().thermal_input_j
            + accepted.radiation().energy().thermal_external_j;
        close(audit.external_heat_j, optical_heat);
        close(audit.stored_energy_change_j, optical_heat);
        close(accepted.illumination().energy().duration_s, 0.25);
        let (dark, _, _) = original
            .advance_observed(joint.radiation(), 0.25, config(0.01), || Ok(()))
            .unwrap();
        assert!(final_state.energy_j() > dark.energy_j());
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(joint.time_s(), 0.0);
    });
}

#[test]
fn joint_candidate_rolls_back_both_sources_when_a_later_owner_refuses() {
    with_cx(|cx| {
        let original = thermal(vec![joint_cell(11, 450.0)], &JOINT_CALORIC, vec![]);
        let pair = joint_fixture(cx, &[1.0], &original, 30.0, 0.4, true);
        let (illumination, radiation) = pair.into_parts();
        let pair = JointOpticalHeating::try_new(
            illumination,
            radiation.with_work_limits(1, 1_000_000),
            &original,
            || Ok(()),
        )
        .unwrap();
        assert!(
            original
                .advance_observed(&pair, 1.0, config(0.1), || Ok(()))
                .is_err()
        );
        assert_eq!(pair.illumination().energy().thermal_input_j, 0.0);
        assert_eq!(pair.radiation().energy().thermal_external_j, 0.0);
        assert_eq!(pair.illumination().source_evaluations(), 0);
        assert_eq!(pair.radiation().source_evaluations(), 0);
        assert_eq!(pair.time_s(), 0.0);
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn admission_refuses_independently_advanced_optical_clocks() {
    with_cx(|cx| {
        let original = thermal(vec![joint_cell(11, 450.0)], &JOINT_CALORIC, vec![]);
        let pair = joint_fixture(cx, &[1.0], &original, 30.0, 0.0, false);
        let (later, illumination, _) = original
            .advance_observed(pair.illumination(), 0.125, config(0.1), || Ok(()))
            .unwrap();
        assert!(
            JointOpticalHeating::try_new(illumination, pair.radiation().clone(), &later, || Ok(()))
                .is_err()
        );
    });
}
