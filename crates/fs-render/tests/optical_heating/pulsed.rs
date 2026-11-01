//! G0/G1/G3/G4: finite spectral fluence on the original native thermal clock.

use super::*;
use fs_blake3::ContentHash;
use fs_render::volumes::heating::exchange::feedback::{
    RadiativeFeedbackBudget, TemperatureRadiativeExchange,
};
use fs_render::volumes::heating::exchange::{RadiativeBath, RadiativeExchangeTolerance};
use fs_render::volumes::heating::feedback::{OpticalFeedbackError, TemperatureAbsorptionTable};
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_render::volumes::heating::illumination::{
    IlluminationBudget, ScheduledIlluminationSource, TemperatureIlluminationHeating,
};
use fs_render::volumes::heating::joint::JointOpticalHeating;
use fs_render::volumes::heating::scattering::ScatteringHeatingTolerance;
use fs_render::volumes::heating::temperature_field::{
    EnthalpyAbsorptionTable, EnthalpyOpticalFieldBudget, EnthalpyScatteringCell,
    EnthalpyScatteringTable, TemperatureOpticalField, TemperatureOpticalFieldBudget,
    TemperatureScatteringTable, ThermalScatteringCell,
};
use fs_render::volumes::transport::scattering::deposition::ScatteringDepositionConfig;
use fs_render::volumes::transport::scattering::exchange::RadiativeExchangeConfig;
use fs_thermal::enthalpy::{DepositionObserver, ThermalLink};
use fs_thermal::observers::CoupledObservers;
use std::collections::BTreeMap;
use std::sync::Arc;

fn limits() -> IlluminationBudget {
    IlluminationBudget {
        max_sources: 4,
        max_source_evaluations: 1000,
        max_traced_segments: 1_000_000,
        max_retained_deposition_cells: 8,
    }
}

fn photons() -> ScatteringDepositionConfig {
    ScatteringDepositionConfig {
        histories_per_lane: 64,
        max_collisions_per_history: 64,
        max_total_segments: 100_000,
        max_deposition_cells: 8,
        seed: 826,
        source_id: 1,
    }
}

fn tolerance() -> ScatteringHeatingTolerance {
    ScatteringHeatingTolerance {
        power_roundoff_w: 1e-10,
        max_unresolved_w: 100.0,
        max_total_standard_error_w: 100.0,
        max_cell_standard_error_w: 100.0,
    }
}

fn source<const N: usize>(segments: &[SpectralPowerSegment<N>]) -> ScheduledIlluminationSource<N> {
    ScheduledIlluminationSource {
        source_id: 11,
        origin: Point3::new(-1.0, 0.5, 0.5),
        unit_direction: Vec3::new(1.0, 0.0, 0.0),
        power: SpectralPowerSchedule::try_new(segments, 16).unwrap(),
    }
}

fn segment<const N: usize>(
    start_s: f64,
    end_s: f64,
    start_power_w: [f64; N],
    end_power_w: [f64; N],
) -> SpectralPowerSegment<N> {
    SpectralPowerSegment {
        start_s,
        end_s,
        start_power_w,
        end_power_w,
    }
}

fn scheduled<'g, const N: usize>(
    field: Arc<TemperatureOpticalField<'g, N>>,
    source: ScheduledIlluminationSource<N>,
    thermal: &DepositionSimulation<'_>,
) -> TemperatureIlluminationHeating<'g, N> {
    TemperatureIlluminationHeating::from_scheduled_field(
        field,
        &[source],
        thermal,
        limits(),
        photons(),
        tolerance(),
        16,
        || Ok(()),
    )
    .unwrap()
}

fn field<'g, const N: usize>(
    cx: &Cx<'_>,
    density: &'g [f64],
    thermal: &DepositionSimulation<'_>,
    wavelengths: [f64; N],
    absorption: [f64; N],
    scattering: f64,
    temperatures: (f64, f64),
) -> Arc<TemperatureOpticalField<'g, N>> {
    let samples: Vec<_> = if N == 1 {
        vec![
            (wavelengths[0] * 0.8, absorption[0]),
            (wavelengths[0] * 1.2, absorption[0]),
        ]
    } else {
        wavelengths.into_iter().zip(absorption).collect()
    };
    let a = MassAbsorptionSpectrum::source_table(&samples, ContentHash([231; 32])).unwrap();
    let s = MassScatteringSpectrum::source_table(
        &[
            (samples[0].0, scattering),
            (samples.last().unwrap().0, scattering),
        ],
        ContentHash([232; 32]),
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
            &[TemperatureAbsorptionTable::try_new(
                ContentHash([233; 32]),
                &[(temperatures.0, a), (temperatures.1, a)],
                2,
            )
            .unwrap()],
            &[TemperatureScatteringTable::try_new(
                ContentHash([234; 32]),
                &[(temperatures.0, s), (temperatures.1, s)],
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

#[test]
fn short_boxes_triangles_and_spectral_changes_deliver_their_exact_fluence() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let input = scheduled(
            field(
                cx,
                &[1.0],
                &original,
                [450.0, 650.0],
                [0.4, 1.2],
                0.0,
                (290.0, 600.0),
            ),
            source(&[
                segment(0.125, 0.25, [8.0, 0.0], [8.0, 0.0]),
                segment(0.5, 0.75, [0.0; 2], [8.0, 12.0]),
                segment(0.75, 1.0, [8.0, 12.0], [0.0; 2]),
                segment(1.0, 1.125, [4.0, 6.0], [4.0, 6.0]),
            ]),
            &original,
        );
        // The requested step is longer than the entire source program.
        let (hot, accepted, audit) = original
            .advance_observed(&input, 1.5, config(2.0), || Ok(()))
            .unwrap();
        let fluence = [3.5, 3.75];
        let mut absorbed = 0.0;
        for (lane, kappa) in [0.4_f64, 1.2].into_iter().enumerate() {
            let q = fluence[lane] * -(-kappa).exp_m1();
            close(accepted.energy().incident_j[lane], fluence[lane]);
            close(accepted.energy().absorbed_j[lane], q);
            close(
                accepted.energy().ballistic_escape_j[lane],
                fluence[lane] * (-kappa).exp(),
            );
            absorbed += q;
        }
        close(hot.cells()[&CellId(11)].specific_h_j_kg, absorbed);
        close(audit.external_heat_j, absorbed);
        close(accepted.energy().thermal_input_j, absorbed);
        assert_eq!(accepted.energy().diffuse_escape_j, [0.0; 2]);
        assert_eq!(accepted.energy().unresolved_j, [0.0; 2]);
        assert_eq!(accepted.energy().duration_s, 1.5);
        assert_eq!(audit.events, 0);
        assert_eq!(audit.substeps, 7);
        assert_eq!(accepted.source_evaluations(), audit.substeps);
        assert_eq!(accepted.transport_builds(), 1);
    });
}

#[test]
fn dark_intervals_draw_no_photons_and_power_only_changes_reuse_scattering() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let input = scheduled(
            field(cx, &[1.0], &original, [550.0], [0.4], 1.5, (290.0, 600.0)),
            source(&[
                segment(1.0, 2.0, [0.0], [8.0]),
                segment(3.0, 4.0, [8.0], [0.0]),
            ]),
            &original,
        );
        let (dark, no_photons, _) = original
            .advance_observed(&input, 1.0, config(0.125), || Ok(()))
            .unwrap();
        // Work-limit setters only tighten. Prove zero-work dark admission
        // separately so the continuing observer retains its future allowance.
        let (_, zero_budget, _) = original
            .advance_observed(&input.with_work_limits(1000, 0), 1.0, config(0.125), || {
                Ok(())
            })
            .unwrap();
        assert_eq!(zero_budget.energy(), no_photons.energy());
        assert_eq!(dark.cells(), original.cells());
        assert_eq!(no_photons.traced_segments(), 0);
        assert_eq!(no_photons.transport_builds(), 0);
        assert!(no_photons.deposition(11).is_none());
        let (mid, sampled, _) = dark
            .advance_observed(&no_photons, 1.5, config(0.125), || Ok(()))
            .unwrap();
        let reference = sampled.deposition(11).unwrap().clone();
        assert!(reference.sampled_histories[0] > 0);
        let (gap, cached, _) = mid
            .advance_observed(
                &sampled.with_work_limits(1000, 0),
                2.5,
                config(0.125),
                || Ok(()),
            )
            .unwrap();
        assert_eq!(cached.transport_builds(), 0);
        assert_eq!(cached.traced_segments(), 0);
        assert!(
            cached
                .sources_w()
                .unwrap()
                .values()
                .all(|watts| *watts == 0.0)
        );
        let (end, accepted, _) = gap
            .advance_observed(&cached, 4.5, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(accepted.deposition(11), Some(&reference));
        assert_eq!(accepted.transport_builds(), 0);
        assert_eq!(accepted.traced_segments(), 0);
        // Two 4 J triangles equal one second at the 8 W reference power.
        close(accepted.energy().incident_j[0], 8.0);
        close(
            accepted.energy().absorbed_j[0],
            reference.absorbed.mean_w[0],
        );
        close(
            accepted.energy().diffuse_escape_j[0],
            reference.diffuse_escape.mean_w[0],
        );
        close(
            accepted.energy().unresolved_j[0],
            reference.unresolved.mean_w[0],
        );
        let (whole, uninterrupted, _) = original
            .advance_observed(&input, 4.5, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(end.cells(), whole.cells());
        assert_eq!(accepted.energy(), uninterrupted.energy());
    });
}

struct PulsePlateau;
impl ThermalLaw for PulsePlateau {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([235; 32])
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && (0.0..=100.0).contains(&h) {
            Ok(500.0)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.0
    }
}
static PULSE_PLATEAU: PulsePlateau = PulsePlateau;

#[test]
fn pulsed_latent_absorption_converges_to_an_independent_fluence_solution() {
    with_cx(|cx| {
        let mut body = cell(11, 1.0, 1.0, 0.0);
        body.specific_h_j_kg = 10.0;
        let original = thermal(vec![body], &PULSE_PLATEAU, vec![]);
        let caloric = ContentHash([235; 32]);
        let zero = MassScatteringSpectrum::constant(0.0).unwrap();
        let field = Arc::new(
            TemperatureOpticalField::try_new_enthalpy(
                grid(cx, &[1.0], MassAbsorptionSpectrum::constant(0.0).unwrap()),
                [550.0],
                &[EnthalpyAbsorptionTable::try_new(
                    ContentHash([236; 32]),
                    caloric,
                    (500.0, 500.0),
                    &[
                        (0.0, MassAbsorptionSpectrum::constant(0.2).unwrap()),
                        (100.0, MassAbsorptionSpectrum::constant(1.2).unwrap()),
                    ],
                    2,
                )
                .unwrap()],
                &[EnthalpyScatteringTable::try_new(
                    ContentHash([237; 32]),
                    caloric,
                    (500.0, 500.0),
                    &[(0.0, zero), (100.0, zero)],
                    2,
                )
                .unwrap()],
                &[PhaseFunction::Isotropic],
                &[EnthalpyScatteringCell {
                    grid_cell: 0,
                    thermal_cell: CellId(11),
                    material_id: MaterialId(7),
                    absorption_table: 0,
                    scattering_table: 0,
                    phase: 0,
                }],
                &original,
                EnthalpyOpticalFieldBudget {
                    max_grid_cells: 1,
                    max_tables: 2,
                    max_total_enthalpy_knots: 4,
                },
                || Ok(()),
            )
            .unwrap(),
        );
        let input = scheduled(
            field,
            source(&[
                segment(0.25, 0.75, [0.0], [40.0]),
                segment(0.75, 1.25, [40.0], [40.0]),
                segment(1.25, 1.75, [40.0], [0.0]),
            ]),
            &original,
        );
        // κ(h)=0.2+0.01*h, rho=L=m=1. Integrating the Beer ODE in
        // incident fluence F gives this closed form for arbitrary P(t).
        let fluence = 40.0;
        let exact = 10.0 + (-(-0.3_f64).exp_m1() * (0.01_f64 * fluence).exp_m1()).ln_1p() / 0.01;
        let mut errors = Vec::new();
        for dt in [0.125, 0.03125] {
            let (hot, accepted, audit) = original
                .advance_observed(&input, 2.0, config(dt), || Ok(()))
                .unwrap();
            let h = hot.cells()[&CellId(11)].specific_h_j_kg;
            assert_eq!(hot.temperature_k(CellId(11)).unwrap(), 500.0);
            close(accepted.energy().incident_j[0], fluence);
            close(h - 10.0, accepted.energy().thermal_input_j);
            close(audit.external_heat_j, h - 10.0);
            assert!(h > 10.0 && h < exact);
            errors.push(exact - h);
        }
        assert!(errors[1] < 0.3 * errors[0], "{errors:?}");
    });
}

#[derive(Clone, Default)]
struct Endpoints(Vec<f64>);
impl DepositionObserver for Endpoints {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        self.0.push(thermal.time_s());
        Ok(())
    }
}

#[test]
fn joint_radiation_and_conduction_reach_pulse_edges_and_cool_after_switch_off() {
    with_cx(|cx| {
        let mut hot = cell(11, 1.0, 1.0, 0.0);
        hot.specific_h_j_kg = 200.0;
        let mut cold = cell(12, 1.0, 1.0, 0.0);
        cold.specific_h_j_kg = 80.0;
        let original = thermal(
            vec![hot, cold],
            &SENSIBLE,
            vec![TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 1,
                    a: CellId(11),
                    b: CellId(12),
                    area_m2: 1.0,
                    conductance_w_k: 0.3,
                }),
            }],
        );
        let illumination = scheduled(
            field(
                cx,
                &[1.0; 2],
                &original,
                [550.0],
                [0.4],
                0.0,
                (290.0, 600.0),
            ),
            source(&[segment(0.25, 0.5, [8.0], [8.0])]),
            &original,
        );
        let radiation = TemperatureRadiativeExchange::try_new(
            field(
                cx,
                &[1.0; 2],
                &original,
                [10_000.0],
                [0.25],
                0.0,
                (290.0, 600.0),
            ),
            [20.0],
            RadiativeBath::Blackbody {
                temperature_k: 300.0,
            },
            &original,
            RadiativeFeedbackBudget {
                max_source_evaluations: 1000,
                max_traced_segments: 1_000_000,
                max_thermal_pairs: 1,
            },
            RadiativeExchangeConfig {
                histories_per_emitter_lane: 32,
                max_collisions_per_history: 64,
                max_total_segments: 100_000,
                max_material_cells: 2,
                max_emitter_cells: 2,
                max_pairs: 1,
                seed: 827,
                source_id: 2,
            },
            RadiativeExchangeTolerance {
                coefficient_roundoff_m2_sr: 1e-10,
                power_roundoff_w: 1e-10,
                max_cell_standard_error_w: 100.0,
                max_total_standard_error_w: 100.0,
                max_unresolved_exchange_w: 1e-6,
                max_unresolved_standard_error_w: 1e-6,
            },
            || Ok(()),
        )
        .unwrap();
        let joint =
            JointOpticalHeating::try_new(illumination, radiation, &original, || Ok(())).unwrap();
        let (off, accepted, audit) = original
            .advance_observed(
                &CoupledObservers::new(joint, Endpoints::default()),
                0.75,
                config(2.0),
                || Ok(()),
            )
            .unwrap();
        assert!(accepted.second.0.contains(&0.25));
        assert!(accepted.second.0.contains(&0.5));
        assert!(audit.substeps < 20);
        let joint = accepted.first;
        close(joint.illumination().energy().incident_j[0], 2.0);
        close(
            audit.external_heat_j,
            joint.illumination().energy().thermal_input_j
                + joint.radiation().energy().thermal_external_j,
        );
        close(audit.stored_energy_change_j, audit.external_heat_j);
        assert_eq!(joint.time_s(), 0.75);
        assert_eq!(joint.radiation().time_s(), 0.75);
        let light_at_switch_off = joint.illumination().energy().thermal_input_j;
        let (cooled, dark, audit) = off
            .advance_observed(&joint, 1.5, config(2.0), || Ok(()))
            .unwrap();
        assert_eq!(
            dark.illumination().energy().thermal_input_j,
            light_at_switch_off
        );
        assert_eq!(dark.illumination().traced_segments(), 0);
        assert_eq!(dark.illumination().transport_builds(), 0);
        assert!(cooled.energy_j() < off.energy_j());
        assert!(cooled.temperature_k(CellId(11)).unwrap() < off.temperature_k(CellId(11)).unwrap());
        assert!(audit.external_heat_j < 0.0);
        close(
            audit.external_heat_j,
            dark.radiation().energy().thermal_external_j
                - joint.radiation().energy().thermal_external_j,
        );
    });
}

#[test]
fn pulse_storage_work_and_cancellation_refuse_without_partial_history() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let field = field(cx, &[1.0], &original, [550.0], [0.4], 0.0, (290.0, 600.0));
        let beam = source(&[segment(0.25, 0.75, [8.0], [8.0])]);
        assert!(matches!(
            TemperatureIlluminationHeating::from_scheduled_field(
                Arc::clone(&field),
                &[beam.clone()],
                &original,
                limits(),
                photons(),
                tolerance(),
                0,
                || Ok(()),
            ),
            Err(OpticalFeedbackError::Budget)
        ));
        let input = scheduled(field, beam, &original);
        assert!(
            original
                .advance_observed(
                    &input.with_work_limits(1, 1_000_000),
                    1.0,
                    config(0.125),
                    || Ok(())
                )
                .is_err()
        );
        let mut checks = 0;
        assert!(
            original
                .advance_observed(&input, 1.0, config(0.125), || {
                    checks += 1;
                    if checks == 40 {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                })
                .is_err()
        );
        assert_eq!(checks, 40);
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 0.0);
        assert_eq!(input.time_s(), 0.0);
        assert_eq!(input.energy().incident_j, [0.0]);
        assert_eq!(input.source_evaluations(), 0);
        assert!(input.deposition(11).is_none());
    });
}

#[test]
fn a_dark_transparent_interval_still_checks_the_accepted_material_domain() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 2.0)], &SENSIBLE, vec![]);
        let input = scheduled(
            field(cx, &[1.0], &original, [550.0], [0.0], 0.0, (300.0, 300.25)),
            source(&[segment(2.0, 3.0, [8.0], [8.0])]),
            &original,
        );
        assert!(
            original
                .advance_observed(&input, 1.0, config(1.0), || Ok(()))
                .is_err()
        );
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(input.time_s(), 0.0);
        assert_eq!(input.energy().thermal_input_j, 0.0);
        let (_, accepted, _) = original
            .advance_observed(&input, 0.25, config(1.0), || Ok(()))
            .unwrap();
        assert_eq!(accepted.traced_segments(), 0);
        assert_eq!(accepted.energy().incident_j, [0.0]);
    });
}
