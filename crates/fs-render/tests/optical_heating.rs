//! Incident optical absorption reaches the original conservative enthalpy cells.

#[path = "optical_heating/scattering.rs"]
mod scattering;

#[path = "optical_heating/feedback.rs"]
mod feedback;

#[path = "optical_heating/scattering_feedback.rs"]
mod scattering_feedback;

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::volumes::VolumeGrid;
use fs_render::volumes::heating::{OpticalHeating, OpticalHeatingError};
use fs_render::volumes::transport::deposition::BeamPowerDeposition;
use fs_render::volumes::transport::scattering::{
    MassScatteringSpectrum, PhaseFunction, ScatteringGrid,
};
use fs_render::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
};

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        operation(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 173,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn grid<'a>(
    cx: &Cx<'_>,
    density: &'a [f64],
    absorption: MassAbsorptionSpectrum,
) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([density.len(), 1, 1], density, [0.0; 3], [1.0; 3]),
        absorption,
        None,
        AbsorbingGridBudget {
            max_grid_cells: 8,
            max_steps_per_segment: 16,
        },
    )
    .unwrap()
}

fn report<const N: usize>(
    cx: &Cx<'_>,
    density: &[f64],
    absorption: MassAbsorptionSpectrum,
    wavelengths: [f64; N],
    powers: [f64; N],
) -> BeamPowerDeposition<N> {
    grid(cx, density, absorption)
        .deposit_beam_power(
            cx,
            Point3::new(-1.0, 0.5, 0.5),
            Vec3::new(1.0, 0.0, 0.0),
            density.len() as f64 + 2.0,
            &wavelengths,
            &powers,
            8,
        )
        .unwrap()
}

struct CaloricLaw {
    latent: bool,
}
impl ThermalLaw for CaloricLaw {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(0.0..=1000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if self.latent {
            if h <= 10.0 {
                300.0 + h
            } else if h <= 110.0 {
                310.0
            } else {
                310.0 + h - 110.0
            }
        } else {
            300.0 + h / 4.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        if self.latent { 1.0 } else { 0.25 }
    }
}
static SENSIBLE: CaloricLaw = CaloricLaw { latent: false };
static LATENT: CaloricLaw = CaloricLaw { latent: true };

fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    }
}
fn cell(id: u64, mass_kg: f64, volume_m3: f64, source_w: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: 1,
        material_id: MaterialId(7),
        mass_kg,
        specific_h_j_kg: 0.0,
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [1.0, 0.0, 0.0],
            volume_m3,
            surface_area_m2: 6.0,
        },
        boundary: boundary(source_w),
    }
}
fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1.0e-10,
    }
}
fn thermal(
    cells: Vec<DepositedCell>,
    law: &'static dyn ThermalLaw,
    events: Vec<TimedEvent>,
) -> DepositionSimulation<'static> {
    let mut schedule: Vec<_> = cells
        .into_iter()
        .map(|c| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(c),
        })
        .collect();
    schedule.extend(events);
    let empty = DepositionSimulation::new(
        schedule,
        &[(MaterialId(7), law)],
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 16,
            max_substeps: 1000,
        },
    )
    .unwrap();
    empty.advance_to(0.0, config(0.1), || Ok(())).unwrap().0
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2.0e-10,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn sourced_uv_and_infrared_power_heat_native_cells_without_visible_extrapolation() {
    with_cx(|cx| {
        let wavelengths = [250.0, 1064.0, 10_000.0];
        let powers = [7.0, 3.0, 11.0];
        let coefficients: [f64; 3] = [0.5, 0.25, 1.0];
        let spectrum = MassAbsorptionSpectrum::source_table(
            &[(250.0, 0.5), (1064.0, 0.25), (10_000.0, 1.0)],
            fs_blake3::ContentHash([29; 32]),
        )
        .unwrap();
        let deposition = report(cx, &[1.0], spectrum, wavelengths, powers);
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let heating = OpticalHeating::try_new(
            deposition.clone(),
            &[(0, CellId(11))],
            &original,
            1e-12,
            || Ok(()),
        )
        .unwrap();
        let (heated, accepted, audit) = original
            .advance_observed(&heating, 2.0, config(0.1), || Ok(()))
            .unwrap();
        let absorbed_w: f64 = (0..3)
            .map(|lane| powers[lane] * -(-coefficients[lane]).exp_m1())
            .sum();
        close(
            heated.temperature_k(CellId(11)).unwrap(),
            300.0 + absorbed_w / 2.0,
        );
        close(audit.external_heat_j, 2.0 * absorbed_w);
        close(accepted.energy().thermal_input_j, audit.external_heat_j);
        assert_eq!(accepted.power().wavelength_bounds_nm, (250.0, 10_000.0));
        for bounds in [(380.0, 780.0), (0.0, 10_000.0), (10_000.0, 250.0)] {
            let mut invalid = deposition.clone();
            invalid.wavelength_bounds_nm = bounds;
            assert!(
                OpticalHeating::try_new(invalid, &[(0, CellId(11))], &original, 1e-12, || Ok(()),)
                    .is_err()
            );
        }
    });
}

#[test]
fn spectral_beam_heats_original_cells_by_their_exact_absorbed_watts() {
    with_cx(|cx| {
        let incident = [7.0, 3.0];
        let coefficients = [0.5, 1.5];
        let spectrum =
            MassAbsorptionSpectrum::table(&[(380.0, coefficients[0]), (780.0, coefficients[1])])
                .unwrap();
        let deposition = report(cx, &[1.0, 2.0], spectrum, [380.0, 780.0], incident);
        let original = thermal(
            vec![cell(11, 1.0, 1.0, 0.0), cell(22, 2.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let heating = OpticalHeating::try_new(
            deposition,
            &[(0, CellId(11)), (1, CellId(22))],
            &original,
            1.0e-12,
            || Ok(()),
        )
        .unwrap();
        let (heated, accepted, audit) = original
            .advance_observed(&heating, 2.0, config(0.17), || Ok(()))
            .unwrap();
        let mut first_w = 0.0;
        let mut second_w = 0.0;
        for lane in 0..2 {
            let attenuation = (-coefficients[lane]).exp();
            first_w += incident[lane] * (1.0 - attenuation);
            second_w += incident[lane] * attenuation * (1.0 - (-2.0 * coefficients[lane]).exp());
            close(accepted.energy().incident_j[lane], 2.0 * incident[lane]);
            close(
                accepted.energy().transmitted_j[lane],
                2.0 * incident[lane] * (-3.0 * coefficients[lane]).exp(),
            );
            close(accepted.energy().scattered_j[lane], 0.0);
        }
        close(
            heated.temperature_k(CellId(11)).unwrap(),
            300.0 + 2.0 * first_w / 4.0,
        );
        close(
            heated.temperature_k(CellId(22)).unwrap(),
            300.0 + 2.0 * second_w / 8.0,
        );
        close(audit.external_heat_j, 2.0 * (first_w + second_w));
        close(audit.stored_energy_change_j, audit.external_heat_j);
        close(accepted.energy().thermal_input_j, audit.external_heat_j);
        close(accepted.energy().duration_s, 2.0);
        assert_eq!(heated.cells()[&CellId(11)].boundary.source_w, 0.0);
        assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 0.0);
        assert_eq!(heating.energy().thermal_input_j, 0.0);
    });
}

#[test]
fn native_latent_storage_absorbs_heat_without_inventing_a_temperature_rise() {
    with_cx(|cx| {
        let deposition = report(
            cx,
            &[1.0],
            MassAbsorptionSpectrum::constant(std::f64::consts::LN_2).unwrap(),
            [550.0],
            [100.0],
        );
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &LATENT, vec![]);
        let heating = OpticalHeating::try_new(
            deposition,
            &[(0, CellId(11))],
            &original,
            1.0e-12,
            || Ok(()),
        )
        .unwrap();
        let (heated, accepted, audit) = original
            .advance_observed(&heating, 1.0, config(0.07), || Ok(()))
            .unwrap();
        close(heated.cells()[&CellId(11)].specific_h_j_kg, 50.0);
        close(heated.temperature_k(CellId(11)).unwrap(), 310.0);
        close(audit.stored_energy_change_j, 50.0);
        close(accepted.energy().absorbed_j[0], 50.0);
        close(accepted.energy().transmitted_j[0], 50.0);
    });
}

#[test]
fn mapping_aggregation_events_and_split_advances_do_not_count_sources_twice() {
    with_cx(|cx| {
        let deposition = report(
            cx,
            &[1.0, 2.0],
            MassAbsorptionSpectrum::constant(0.7).unwrap(),
            [550.0],
            [10.0],
        );
        let absorbed = deposition.absorbed_power_w[0];
        let original = thermal(
            vec![cell(11, 3.0, 2.0, 2.0)],
            &SENSIBLE,
            vec![TimedEvent {
                time_s: 0.375,
                event: Event::Boundary(CellId(11), boundary(5.0)),
            }],
        );
        let heating = OpticalHeating::try_new(
            deposition,
            &[(0, CellId(11)), (1, CellId(11))],
            &original,
            1.0e-12,
            || Ok(()),
        )
        .unwrap();
        assert_eq!(heating.sources_w().len(), 1);
        let (whole, whole_optics, whole_audit) = original
            .advance_observed(&heating, 1.0, config(0.2), || Ok(()))
            .unwrap();
        let (middle, middle_optics, first) = original
            .advance_observed(&heating, 0.375, config(0.2), || Ok(()))
            .unwrap();
        let (split, split_optics, second) = middle
            .advance_observed(&middle_optics, 1.0, config(0.2), || Ok(()))
            .unwrap();
        let authored_heat = 2.0 * 0.375 + 5.0 * (1.0 - 0.375);
        close(whole_audit.external_heat_j, authored_heat + absorbed);
        close(
            first.external_heat_j + second.external_heat_j,
            whole_audit.external_heat_j,
        );
        close(
            split.temperature_k(CellId(11)).unwrap(),
            whole.temperature_k(CellId(11)).unwrap(),
        );
        close(whole_optics.energy().thermal_input_j, absorbed);
        close(split_optics.energy().thermal_input_j, absorbed);
        assert_eq!(whole.cells()[&CellId(11)].boundary.source_w, 5.0);
        assert_eq!(original.cells()[&CellId(11)].boundary.source_w, 2.0);
        assert_eq!(split_optics.time_s(), 1.0);
    });
}

#[test]
fn pure_scattering_does_not_heat_the_thermal_receiver() {
    with_cx(|cx| {
        let density = [1.0, 2.0];
        let medium = ScatteringGrid::try_new(
            cx,
            grid(cx, &density, MassAbsorptionSpectrum::constant(0.0).unwrap()),
            MassScatteringSpectrum::constant(2.0).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let deposition = medium
            .deposit_beam_power(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                4.0,
                &[550.0],
                &[12.0],
                8,
            )
            .unwrap();
        let original = thermal(vec![cell(11, 3.0, 2.0, 0.0)], &SENSIBLE, vec![]);
        let heating =
            OpticalHeating::try_new(deposition, &[], &original, 1.0e-12, || Ok(())).unwrap();
        assert!(heating.sources_w().is_empty());
        let (same, accepted, audit) = original
            .advance_observed(&heating, 2.0, config(0.3), || Ok(()))
            .unwrap();
        close(same.temperature_k(CellId(11)).unwrap(), 300.0);
        assert_eq!(audit.external_heat_j, 0.0);
        assert_eq!(accepted.energy().thermal_input_j, 0.0);
        assert_eq!(accepted.energy().absorbed_j, [0.0]);
        close(
            accepted.energy().scattered_j[0] + accepted.energy().transmitted_j[0],
            24.0,
        );
    });
}

#[test]
fn invalid_reports_mappings_and_stale_clocks_refuse_without_physical_changes() {
    with_cx(|cx| {
        let deposition = report(
            cx,
            &[1.0, 2.0],
            MassAbsorptionSpectrum::constant(1.0).unwrap(),
            [550.0],
            [10.0],
        );
        let original = thermal(vec![cell(11, 3.0, 2.0, 0.0)], &SENSIBLE, vec![]);
        assert!(matches!(
            OpticalHeating::try_new(
                deposition.clone(),
                &[(0, CellId(11))],
                &original,
                1.0e-12,
                || Ok(()),
            ),
            Err(OpticalHeatingError::UnmappedCell(1))
        ));
        assert!(matches!(
            OpticalHeating::try_new(
                deposition.clone(),
                &[(0, CellId(11)), (0, CellId(11))],
                &original,
                1.0e-12,
                || Ok(()),
            ),
            Err(OpticalHeatingError::DuplicateCell(0))
        ));
        assert!(matches!(
            OpticalHeating::try_new(
                deposition.clone(),
                &[(0, CellId(99)), (1, CellId(11))],
                &original,
                1.0e-12,
                || Ok(()),
            ),
            Err(OpticalHeatingError::InactiveCell(CellId(99)))
        ));
        let bindings = [(0, CellId(11)), (1, CellId(11))];
        let mut nonconservative = deposition.clone();
        nonconservative.transmitted_power_w[0] += 1.0;
        assert!(
            OpticalHeating::try_new(nonconservative, &bindings, &original, 1.0e-12, || Ok(()))
                .is_err()
        );
        let mut inconsistent = deposition.clone();
        inconsistent.cells[0].absorbed_power_w[0] *= 0.5;
        assert!(
            OpticalHeating::try_new(inconsistent, &bindings, &original, 1.0e-12, || Ok(()))
                .is_err()
        );
        let mut active = deposition.clone();
        active.cells[0].absorbed_power_w[0] = -1.0;
        assert!(OpticalHeating::try_new(active, &bindings, &original, 1.0e-12, || Ok(())).is_err());
        let heating =
            OpticalHeating::try_new(deposition, &bindings, &original, 1.0e-12, || Ok(())).unwrap();
        let advanced_without_light = original.advance_to(0.1, config(0.1), || Ok(())).unwrap().0;
        assert!(matches!(
            advanced_without_light.advance_observed(&heating, 0.2, config(0.1), || Ok(())),
            Err(OpticalHeatingError::Invalid(_))
        ));
        assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 0.0);
        assert_eq!(heating.time_s(), 0.0);
        assert_eq!(heating.energy().thermal_input_j, 0.0);
    });
}

#[test]
fn cancellation_and_caloric_refusal_roll_back_heat_and_optical_accounts_together() {
    with_cx(|cx| {
        let deposition = report(
            cx,
            &[1.0],
            MassAbsorptionSpectrum::constant(10.0).unwrap(),
            [550.0],
            [20.0],
        );
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let heating = OpticalHeating::try_new(
            deposition,
            &[(0, CellId(11))],
            &original,
            1.0e-12,
            || Ok(()),
        )
        .unwrap();
        let mut checks = 0;
        let cancelled = original.advance_observed(&heating, 1.0, config(0.01), || {
            checks += 1;
            if checks == 100 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            cancelled,
            Err(OpticalHeatingError::Thermal(DepositionError::Cancelled))
        ));
        assert!(matches!(
            original.advance_observed(&heating, 100.0, config(0.5), || Ok(())),
            Err(OpticalHeatingError::Thermal(
                DepositionError::MaterialDomain
            ))
        ));
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 0.0);
        assert_eq!(heating.time_s(), 0.0);
        assert_eq!(heating.energy().duration_s, 0.0);
        let (heated, accepted, audit) = original
            .advance_observed(&heating, 1.0, config(0.1), || Ok(()))
            .unwrap();
        close(
            heated.cells()[&CellId(11)].specific_h_j_kg,
            audit.external_heat_j,
        );
        close(accepted.energy().thermal_input_j, audit.external_heat_j);
    });
}
