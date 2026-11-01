//! G0/G1/G3/G4: native latent inventory drives the original optical owners.

use super::*;
use fs_blake3::ContentHash;
use fs_render::volumes::heating::exchange::feedback::{
    RadiativeFeedbackBudget, TemperatureRadiativeExchange,
};
use fs_render::volumes::heating::exchange::{RadiativeBath, RadiativeExchangeTolerance};
use fs_render::volumes::heating::feedback::OpticalFeedbackError;
use fs_render::volumes::heating::illumination::{
    IlluminationBudget, IncidentIlluminationSource, TemperatureIlluminationHeating,
};
use fs_render::volumes::heating::joint::JointOpticalHeating;
use fs_render::volumes::heating::scattering::ScatteringHeatingTolerance;
use fs_render::volumes::heating::scattering_feedback::IncidentScatteringBeam;
use fs_render::volumes::heating::temperature_field::{
    EnthalpyAbsorptionTable, EnthalpyOpticalFieldBudget, EnthalpyScatteringCell,
    EnthalpyScatteringTable, TemperatureOpticalField,
};
use fs_render::volumes::transport::scattering::deposition::ScatteringDepositionConfig;
use fs_render::volumes::transport::scattering::exchange::RadiativeExchangeConfig;
use fs_thermal::enthalpy::ThermalStorageBasis;
use std::sync::Arc;

/// A sharp fusion plateau with an explicit caloric identity and arbitrary datum.
/// Optical tables may resolve h inside [datum, datum+100] although T=500 K.
struct Plateau {
    datum: f64,
    identity: u8,
    internal_energy: bool,
}
impl ThermalLaw for Plateau {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        (self.identity != 0).then_some([self.identity; 32])
    }
    fn storage_basis(&self) -> ThermalStorageBasis {
        if self.internal_energy {
            ThermalStorageBasis::InternalEnergy
        } else {
            ThermalStorageBasis::Enthalpy
        }
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let x = h - self.datum;
        if !x.is_finite() || !(-100.0..=200.0).contains(&x) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if x < 0.0 {
            500.0 + x
        } else if x <= 100.0 {
            500.0
        } else {
            500.0 + x - 100.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
static PLATEAU: Plateau = Plateau {
    datum: 0.0,
    identity: 51,
    internal_energy: false,
};
static SHIFTED: Plateau = Plateau {
    datum: -1000.0,
    identity: 52,
    internal_energy: false,
};
static INTERNAL: Plateau = Plateau {
    datum: 0.0,
    identity: 51,
    internal_energy: true,
};
static UNIDENTIFIED: Plateau = Plateau {
    datum: 0.0,
    identity: 0,
    internal_energy: false,
};

fn specimen(law: &'static Plateau, x: f64, applied_w: f64) -> DepositionSimulation<'static> {
    let mut body = cell(11, 1.0, 1.0, applied_w);
    body.specific_h_j_kg = law.datum + x;
    thermal(vec![body], law, vec![])
}

#[derive(Clone, Copy)]
struct OpticalData {
    datum: f64,
    caloric_identity: u8,
    material_id: MaterialId,
    wavelength_nm: f64,
    span_j_kg: (f64, f64),
    absorption: [f64; 2],
    scattering: [f64; 2],
}
fn optical_data() -> OpticalData {
    OpticalData {
        datum: 0.0,
        caloric_identity: 51,
        material_id: MaterialId(7),
        wavelength_nm: 550.0,
        span_j_kg: (0.0, 100.0),
        absorption: [0.2, 1.2],
        scattering: [0.0; 2],
    }
}

fn enthalpy_field(
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'_>,
    data: OpticalData,
) -> Result<Arc<TemperatureOpticalField<'static, 1>>, OpticalFeedbackError> {
    let hs = [data.datum + data.span_j_kg.0, data.datum + data.span_j_kg.1];
    let support = (0.8 * data.wavelength_nm, 1.2 * data.wavelength_nm);
    let absorption: Vec<_> = hs
        .into_iter()
        .zip(data.absorption)
        .map(|(h, k)| {
            (
                h,
                MassAbsorptionSpectrum::source_table(
                    &[(support.0, k), (support.1, k)],
                    ContentHash([61; 32]),
                )
                .unwrap(),
            )
        })
        .collect();
    let scattering: Vec<_> = hs
        .into_iter()
        .zip(data.scattering)
        .map(|(h, k)| {
            (
                h,
                MassScatteringSpectrum::source_table(
                    &[(support.0, k), (support.1, k)],
                    ContentHash([62; 32]),
                )
                .unwrap(),
            )
        })
        .collect();
    let caloric = ContentHash([data.caloric_identity; 32]);
    TemperatureOpticalField::try_new_enthalpy(
        grid(cx, &[1.0], MassAbsorptionSpectrum::constant(0.0).unwrap()),
        [data.wavelength_nm],
        &[EnthalpyAbsorptionTable::try_new(
            ContentHash([63; 32]),
            caloric,
            (450.0, 550.0),
            &absorption,
            2,
        )?],
        &[EnthalpyScatteringTable::try_new(
            ContentHash([64; 32]),
            caloric,
            (450.0, 550.0),
            &scattering,
            2,
        )?],
        &[PhaseFunction::Isotropic],
        &[EnthalpyScatteringCell {
            grid_cell: 0,
            thermal_cell: CellId(11),
            material_id: data.material_id,
            absorption_table: 0,
            scattering_table: 0,
            phase: 0,
        }],
        thermal,
        EnthalpyOpticalFieldBudget {
            max_grid_cells: 1,
            max_tables: 2,
            max_total_enthalpy_knots: 4,
        },
        || Ok(()),
    )
    .map(Arc::new)
}

fn illuminate(
    field: Arc<TemperatureOpticalField<'static, 1>>,
    thermal: &DepositionSimulation<'_>,
    power_w: f64,
) -> TemperatureIlluminationHeating<'static, 1> {
    let wavelengths_nm = *field.wavelengths_nm();
    TemperatureIlluminationHeating::from_field(
        field,
        &[IncidentIlluminationSource {
            source_id: 11,
            beam: IncidentScatteringBeam {
                origin: Point3::new(-1.0, 0.5, 0.5),
                unit_direction: Vec3::new(1.0, 0.0, 0.0),
                wavelengths_nm,
                incident_power_w: [power_w],
            },
        }],
        thermal,
        IlluminationBudget {
            max_sources: 1,
            max_source_evaluations: 1000,
            max_traced_segments: 1_000_000,
            max_retained_deposition_cells: 1,
        },
        ScatteringDepositionConfig {
            histories_per_lane: 64,
            max_collisions_per_history: 96,
            max_total_segments: 100_000,
            max_deposition_cells: 1,
            seed: 314,
            source_id: 159,
        },
        ScatteringHeatingTolerance {
            power_roundoff_w: 1e-10,
            max_unresolved_w: 1e-7,
            max_total_standard_error_w: 100.0,
            max_cell_standard_error_w: 100.0,
        },
        || Ok(()),
    )
    .unwrap()
}

// One metre, rho=1 kg/m³ and m=1 kg. With κ=a+b*h the Beer ODE is
// h'=P*(1-exp(-a-b*h)); exp(b*h) solves a constant-coefficient linear ODE.
fn exact_plateau_h(initial_h: f64, power_w: f64, duration_s: f64) -> f64 {
    let a = 0.2;
    let b = 0.01;
    initial_h + (-(-(a + b * initial_h)).exp_m1() * (b * power_w * duration_s).exp_m1()).ln_1p() / b
}

#[test]
fn isothermal_melting_changes_absorption_and_converges_to_independent_beer_ode() {
    with_cx(|cx| {
        let original = specimen(&PLATEAU, 10.0, 0.0);
        let field = enthalpy_field(cx, &original, optical_data()).unwrap();
        let input = illuminate(Arc::clone(&field), &original, 20.0);
        let exact = exact_plateau_h(10.0, 20.0, 4.0);
        let mut errors = Vec::new();
        for dt in [0.25, 0.125, 0.0625] {
            let (heated, accepted, audit) = original
                .advance_observed(&input, 4.0, config(dt), || Ok(()))
                .unwrap();
            let h = heated.cells()[&CellId(11)].specific_h_j_kg;
            errors.push(exact - h);
            assert!(h > 10.0 && h < exact && h < 100.0);
            assert_eq!(heated.temperature_k(CellId(11)).unwrap(), 500.0);
            close(h - 10.0, audit.external_heat_j);
            close(
                audit.stored_energy_change_j,
                accepted.energy().thermal_input_j,
            );
            close(accepted.energy().absorbed_j[0], h - 10.0);
            close(accepted.energy().incident_j[0], 80.0);
            close(
                accepted.energy().absorbed_j[0] + accepted.energy().ballistic_escape_j[0],
                80.0,
            );
            assert_eq!(accepted.transport_builds(), audit.substeps);
            assert_eq!(accepted.traced_segments(), audit.substeps as u64);
            assert_eq!(accepted.deposition(11).unwrap().sampled_histories, [0]);
            let packet = field.initial_state(&heated, || Ok(())).unwrap();
            close(packet.absorption_m2_kg()[0][0], 0.2 + 0.01 * h);
            assert!(
                accepted.deposition(11).unwrap().absorbed.mean_w[0] > 20.0 * -(-0.3_f64).exp_m1()
            );
        }
        for pair in errors.windows(2) {
            assert!(
                pair[1] > 0.45 * pair[0] && pair[1] < 0.55 * pair[0],
                "{errors:?}"
            );
        }
        assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 10.0);
        assert_eq!(input.energy().thermal_input_j, 0.0);
    });
}

#[test]
fn visible_absorption_and_infrared_exchange_follow_the_same_isothermal_phase_state() {
    with_cx(|cx| {
        let original = specimen(&PLATEAU, 10.0, 0.0);
        let visible = enthalpy_field(cx, &original, optical_data()).unwrap();
        let infrared = enthalpy_field(
            cx,
            &original,
            OpticalData {
                wavelength_nm: 10_000.0,
                absorption: [0.1, 0.6],
                ..optical_data()
            },
        )
        .unwrap();
        let illumination = illuminate(visible, &original, 20.0);
        let radiation = TemperatureRadiativeExchange::try_new(
            Arc::clone(&infrared),
            [20.0],
            RadiativeBath::Blackbody {
                temperature_k: 500.0,
            },
            &original,
            RadiativeFeedbackBudget {
                max_source_evaluations: 1000,
                max_traced_segments: 1_000_000,
                max_thermal_pairs: 1,
            },
            RadiativeExchangeConfig {
                histories_per_emitter_lane: 16,
                max_collisions_per_history: 96,
                max_total_segments: 100_000,
                max_material_cells: 1,
                max_emitter_cells: 1,
                max_pairs: 1,
                seed: 271,
                source_id: 828,
            },
            RadiativeExchangeTolerance {
                coefficient_roundoff_m2_sr: 1e-9,
                power_roundoff_w: 1e-9,
                max_cell_standard_error_w: 1000.0,
                max_total_standard_error_w: 1000.0,
                max_unresolved_exchange_w: 1e-7,
                max_unresolved_standard_error_w: 1e-7,
            },
            || Ok(()),
        )
        .unwrap();
        let input =
            JointOpticalHeating::try_new(illumination, radiation, &original, || Ok(())).unwrap();
        let (heated, accepted, audit) = original
            .advance_observed(&input, 1.0, config(0.125), || Ok(()))
            .unwrap();
        let h = heated.cells()[&CellId(11)].specific_h_j_kg;
        assert!(h > 10.0);
        assert_eq!(heated.temperature_k(CellId(11)).unwrap(), 500.0);
        assert_eq!(accepted.radiation().energy().thermal_external_j, 0.0);
        assert_eq!(accepted.radiation().energy().net_bath_to_material_j, [0.0]);
        assert_eq!(
            accepted.radiation().energy().internal_transfer_magnitude_j,
            [0.0]
        );
        assert!(accepted.radiation().operator_builds() > 1);
        assert_eq!(accepted.illumination().transport_builds(), audit.substeps);
        close(
            audit.external_heat_j,
            accepted.illumination().energy().thermal_input_j,
        );
        close(audit.stored_energy_change_j, h - 10.0);
        let packet = infrared.initial_state(&heated, || Ok(())).unwrap();
        close(packet.absorption_m2_kg()[0][0], 0.1 + 0.005 * h);
    });
}

#[test]
fn negative_enthalpy_datum_preserves_physical_absorption_and_latent_heat() {
    with_cx(|cx| {
        let original = specimen(&PLATEAU, 10.0, 0.0);
        let shifted = specimen(&SHIFTED, 10.0, 0.0);
        let input = illuminate(
            enthalpy_field(cx, &original, optical_data()).unwrap(),
            &original,
            20.0,
        );
        let shifted_input = illuminate(
            enthalpy_field(
                cx,
                &shifted,
                OpticalData {
                    datum: -1000.0,
                    caloric_identity: 52,
                    ..optical_data()
                },
            )
            .unwrap(),
            &shifted,
            20.0,
        );
        assert_ne!(input.model_identity(), shifted_input.model_identity());
        let (a, heat_a, audit_a) = original
            .advance_observed(&input, 2.0, config(0.125), || Ok(()))
            .unwrap();
        let (b, heat_b, audit_b) = shifted
            .advance_observed(&shifted_input, 2.0, config(0.125), || Ok(()))
            .unwrap();
        close(
            a.cells()[&CellId(11)].specific_h_j_kg,
            b.cells()[&CellId(11)].specific_h_j_kg + 1000.0,
        );
        close(heat_a.energy().absorbed_j[0], heat_b.energy().absorbed_j[0]);
        close(audit_a.external_heat_j, audit_b.external_heat_j);
        assert_eq!(
            a.temperature_k(CellId(11)).unwrap(),
            b.temperature_k(CellId(11)).unwrap()
        );
    });
}

#[test]
fn incompatible_caloric_identity_material_and_internal_energy_are_refused() {
    with_cx(|cx| {
        let original = specimen(&PLATEAU, 10.0, 0.0);
        for data in [
            OpticalData {
                caloric_identity: 52,
                ..optical_data()
            },
            OpticalData {
                material_id: MaterialId(8),
                ..optical_data()
            },
        ] {
            assert!(matches!(
                enthalpy_field(cx, &original, data),
                Err(OpticalFeedbackError::Invalid(_))
            ));
        }
        for law in [&INTERNAL, &UNIDENTIFIED] {
            let incompatible = specimen(law, 10.0, 0.0);
            assert!(matches!(
                enthalpy_field(cx, &incompatible, optical_data()),
                Err(OpticalFeedbackError::Invalid(_))
            ));
        }
        assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 10.0);
    });
}

#[test]
fn transparent_material_still_checks_endpoint_enthalpy_and_rolls_back_the_native_step() {
    with_cx(|cx| {
        let original = specimen(&PLATEAU, 1.0, 10.0);
        let field = enthalpy_field(
            cx,
            &original,
            OpticalData {
                span_j_kg: (0.0, 5.0),
                absorption: [0.0; 2],
                ..optical_data()
            },
        )
        .unwrap();
        let input = illuminate(field, &original, 0.0);
        assert!(matches!(
            original.advance_observed(&input, 1.0, config(0.25), || Ok(())),
            Err(OpticalFeedbackError::EnthalpyDomain {
                thermal_cell: CellId(11),
                ..
            })
        ));
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 1.0);
        assert_eq!(input.time_s(), 0.0);
        assert_eq!(input.energy().thermal_input_j, 0.0);
        let (admitted, accepted, audit) = original
            .advance_observed(&input, 0.25, config(0.25), || Ok(()))
            .unwrap();
        assert_eq!(admitted.temperature_k(CellId(11)).unwrap(), 500.0);
        close(admitted.cells()[&CellId(11)].specific_h_j_kg, 3.5);
        close(audit.external_heat_j, 2.5);
        assert_eq!(accepted.energy().thermal_input_j, 0.0);
    });
}

#[test]
fn constant_enthalpy_coefficients_reuse_transport_and_continue_transactionally() {
    with_cx(|cx| {
        let original = specimen(&PLATEAU, 10.0, 0.0);
        let field = enthalpy_field(
            cx,
            &original,
            OpticalData {
                absorption: [0.4; 2],
                ..optical_data()
            },
        )
        .unwrap();
        let input = illuminate(field, &original, 20.0);
        let (mid, accepted, _) = original
            .advance_observed(&input, 1.0, config(0.25), || Ok(()))
            .unwrap();
        assert_eq!(accepted.transport_builds(), 1);
        let cached = accepted.with_work_limits(1000, 0);
        let (end, continued, _) = mid
            .advance_observed(&cached, 2.0, config(0.25), || Ok(()))
            .unwrap();
        assert_eq!(continued.transport_builds(), 0);
        assert_eq!(continued.traced_segments(), 0);
        let (whole, whole_input, _) = original
            .advance_observed(&input, 2.0, config(0.25), || Ok(()))
            .unwrap();
        assert_eq!(end.cells(), whole.cells());
        assert_eq!(continued.energy(), whole_input.energy());
        let mut checkpoints = 0;
        assert!(matches!(
            mid.advance_observed(&cached, 2.0, config(0.25), || {
                checkpoints += 1;
                if checkpoints == 25 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(OpticalFeedbackError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(cached.time_s(), 1.0);
        assert_eq!(cached.energy(), accepted.energy());
        assert_eq!(mid.time_s(), 1.0);
    });
}

#[test]
fn isothermal_enthalpy_change_can_activate_scattering_from_an_exact_absorption_lane() {
    with_cx(|cx| {
        let original = specimen(&PLATEAU, 0.0, 10.0);
        let field = enthalpy_field(
            cx,
            &original,
            OpticalData {
                absorption: [0.4; 2],
                scattering: [0.0, 100.0],
                ..optical_data()
            },
        )
        .unwrap();
        assert!(field.has_scattering());
        let input = illuminate(Arc::clone(&field), &original, 5.0);
        let (mid, first, _) = original
            .advance_observed(&input, 0.125, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(first.deposition(11).unwrap().sampled_histories, [0]);
        let packet = field.initial_state(&mid, || Ok(())).unwrap();
        assert!(packet.scattering_m2_kg()[0][0] > 1.0);
        let (end, continued, audit) = mid
            .advance_observed(&first, 0.25, config(0.125), || Ok(()))
            .unwrap();
        let report = continued.deposition(11).unwrap();
        assert_eq!(report.sampled_histories, [64]);
        assert!(report.diffuse_escape.mean_w[0] > 0.0);
        assert_eq!(end.temperature_k(CellId(11)).unwrap(), 500.0);
        close(
            report.absorbed.mean_w[0]
                + report.ballistic_escape.mean_w[0]
                + report.diffuse_escape.mean_w[0]
                + report.unresolved.mean_w[0],
            5.0,
        );
        close(
            audit.external_heat_j,
            10.0 * 0.125 + report.absorbed.mean_w[0] * 0.125,
        );
        close(audit.stored_energy_change_j, audit.external_heat_j);
    });
}
