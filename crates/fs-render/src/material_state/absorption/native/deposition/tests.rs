//! Accepted source laws become conservative incident heat in original cells.

use super::super::tests::{LENGTH, budget, close, direction, geometry, origin, palette, with_cx};
use super::super::{MaterialAbsorptionEmission, MaterialAbsorptionField};
use super::*;
use crate::material_state::absorption::MaterialAbsorptionBinding;
use crate::material_state::absorption::tests::{complex_fixture, density_key, point};

#[test]
fn g1_original_material_laws_deposit_serial_power_without_counting_lte_emission() {
    with_cx(|cx, _| {
        let (card, cold, cold_optics) = complex_fixture([0.01, 0.03]);
        let hot = cold.at_query_point(&card, &point(400.0)).unwrap();
        let hot_optics = cold_optics.at_query_point(&card, &point(400.0)).unwrap();
        let materials = palette(
            cx,
            vec![
                MaterialAbsorptionBinding::try_spectral(&cold, &cold_optics, &density_key())
                    .unwrap(),
                MaterialAbsorptionBinding::try_spectral(&hot, &hot_optics, &density_key()).unwrap(),
            ],
            &[0, 1],
        );
        let density = [1000.0, 800.0];
        let geometry = geometry(cx, &density, 16);
        let wavelengths = [8000.0, 10000.0, 12000.0];
        let powers = [2.0, 3.0, 4.0];
        let field = MaterialAbsorptionField::try_new(
            cx,
            &geometry,
            materials.clone(),
            wavelengths,
            MaterialAbsorptionEmission::AcceptedLte,
            budget(),
        )
        .unwrap();
        let no_emission = MaterialAbsorptionField::try_new(
            cx,
            &geometry,
            materials,
            wavelengths,
            MaterialAbsorptionEmission::None,
            budget(),
        )
        .unwrap();
        assert!(
            field
                .transport(cx, origin(), direction(), f64::INFINITY)
                .unwrap()
                .emission_w_m2_sr_nm
                .iter()
                .all(|v| *v > 0.0)
        );
        let report = field
            .deposit_beam_power(cx, origin(), direction(), f64::INFINITY, &powers, 2)
            .unwrap();
        let cold_report = no_emission
            .deposit_beam_power(cx, origin(), direction(), f64::INFINITY, &powers, 2)
            .unwrap();
        assert_eq!(report.cells, cold_report.cells);
        assert_eq!(report.absorbed_power_w, cold_report.absorbed_power_w);
        assert_eq!(report.medium_identity, field.content_identity());
        assert_eq!(report.wavelength_bounds_nm, (8000.0, 12000.0));
        assert_eq!(report.wavelengths_nm, wavelengths);
        assert_eq!(
            report
                .cells
                .iter()
                .map(|c| c.cell_index)
                .collect::<Vec<_>>(),
            vec![0, 1]
        );
        let reverse = field
            .deposit_beam_power(
                cx,
                Point3::new(0.0, 0.0, 3.0 * LENGTH),
                Vec3::new(0.0, 0.0, -1.0),
                f64::INFINITY,
                &powers,
                2,
            )
            .unwrap();
        let finite = field
            .deposit_beam_power(cx, origin(), direction(), 2.0 * LENGTH, &powers, 1)
            .unwrap();
        assert_eq!(
            reverse
                .cells
                .iter()
                .map(|c| c.cell_index)
                .collect::<Vec<_>>(),
            vec![1, 0]
        );
        assert_eq!(finite.cells.len(), 1);
        for lane in 0..3 {
            let tau = 4.0 * core::f64::consts::PI * [0.01, 0.02, 0.03][lane] * LENGTH
                / (wavelengths[lane] * 1e-9);
            let first = powers[lane] * -(-tau).exp_m1();
            let second = powers[lane] * (-tau).exp() * -(-2.0 * tau).exp_m1();
            close(report.cells[0].absorbed_power_w[lane], first);
            close(report.cells[1].absorbed_power_w[lane], second);
            close(
                report.transmitted_power_w[lane],
                powers[lane] * (-3.0 * tau).exp(),
            );
            close(report.absorbed_power_w[lane], first + second);
            close(
                report.transmitted_power_w[lane] + report.absorbed_power_w[lane],
                powers[lane],
            );
            close(
                reverse.cells[0].absorbed_power_w[lane],
                powers[lane] * -(-2.0 * tau).exp_m1(),
            );
            close(
                reverse.cells[1].absorbed_power_w[lane],
                powers[lane] * (-2.0 * tau).exp() * -(-tau).exp_m1(),
            );
            close(
                reverse.transmitted_power_w[lane],
                report.transmitted_power_w[lane],
            );
            close(finite.absorbed_power_w[lane], first);
            close(
                finite.transmitted_power_w[lane],
                powers[lane] * (-tau).exp(),
            );
            assert_eq!(report.scattered_power_w[lane], 0.0);
            assert_eq!(report.cells[0].scattered_power_w[lane], 0.0);
            assert_eq!(report.cells[1].scattered_power_w[lane], 0.0);
        }
        let dark = field
            .deposit_beam_power(cx, origin(), direction(), f64::INFINITY, &[0.0; 3], 0)
            .unwrap();
        assert_eq!(dark.absorbed_power_w, [0.0; 3]);
        assert!(dark.cells.is_empty());
    });
}

#[test]
fn g1_extreme_source_columns_keep_representable_absorbed_and_transmitted_watts() {
    with_cx(|cx, _| {
        for k in [1e-300, 800.0 / (4.0 * core::f64::consts::PI)] {
            let (_, accepted, optical) = complex_fixture([k; 2]);
            let material =
                MaterialAbsorptionBinding::try_spectral(&accepted, &optical, &density_key())
                    .unwrap();
            let density = [1000.0];
            let geometry = geometry(cx, &density, 16);
            let field = MaterialAbsorptionField::try_new(
                cx,
                &geometry,
                palette(cx, vec![material], &[0]),
                [10000.0],
                MaterialAbsorptionEmission::None,
                budget(),
            )
            .unwrap();
            let report = field
                .deposit_beam_power(cx, origin(), direction(), f64::INFINITY, &[1e300], 1)
                .unwrap();
            let tau = 4.0 * core::f64::consts::PI * k;
            close(report.absorbed_power_w[0], 1e300 * -(-tau).exp_m1());
            close(report.transmitted_power_w[0], (1e300_f64.ln() - tau).exp());
            assert!(report.absorbed_power_w[0] > 0.0);
            if k < 1.0 {
                assert_eq!(
                    field
                        .transmittance(cx, origin(), direction(), f64::INFINITY)
                        .unwrap(),
                    [1.0]
                );
            } else {
                assert_eq!(
                    field
                        .transmittance(cx, origin(), direction(), f64::INFINITY)
                        .unwrap(),
                    [0.0]
                );
                assert!(report.transmitted_power_w[0] > 0.0);
            }
        }
    });
}

#[test]
fn g0_vacuum_limits_and_invalid_inputs_keep_native_admission() {
    with_cx(|cx, gate| {
        let (_, clear_state, clear_optical) = complex_fixture([0.0; 2]);
        let (_, absorbing_state, absorbing_optical) = complex_fixture([0.01, 0.03]);
        let clear =
            MaterialAbsorptionBinding::try_spectral(&clear_state, &clear_optical, &density_key())
                .unwrap();
        let absorber = MaterialAbsorptionBinding::try_spectral(
            &absorbing_state,
            &absorbing_optical,
            &density_key(),
        )
        .unwrap();
        let density = [1000.0; 2];
        let one_step = geometry(cx, &density, 1);
        let clear_field = MaterialAbsorptionField::try_new(
            cx,
            &one_step,
            palette(cx, vec![clear, absorber], &[0, 0]),
            [10000.0],
            MaterialAbsorptionEmission::None,
            budget(),
        )
        .unwrap();
        let report = clear_field
            .deposit_beam_power(cx, origin(), direction(), f64::INFINITY, &[3.0], 0)
            .unwrap();
        assert_eq!(report.transmitted_power_w, [3.0]);
        assert!(report.cells.is_empty());
        for distance in [f64::NAN, -1.0] {
            assert_eq!(
                clear_field.deposit_beam_power(cx, origin(), direction(), distance, &[3.0], 0),
                Err(VolumeTransportError::Ray)
            );
        }
        assert_eq!(
            clear_field.deposit_beam_power(cx, origin(), Vec3::new(0.0, 0.0, 0.0), 0.0, &[3.0], 0),
            Err(VolumeTransportError::Ray)
        );
        assert_eq!(
            clear_field.deposit_beam_power(cx, origin(), direction(), 0.0, &[f64::NAN], 0),
            Err(VolumeTransportError::BeamPower { lane: 0 })
        );
        let absorbing = MaterialAbsorptionField::try_new(
            cx,
            &one_step,
            palette(cx, vec![absorber], &[0, 0]),
            [10000.0],
            MaterialAbsorptionEmission::None,
            budget(),
        )
        .unwrap();
        assert_eq!(
            absorbing.deposit_beam_power(cx, origin(), direction(), f64::INFINITY, &[3.0], 2),
            Err(VolumeTransportError::Budget)
        );
        assert_eq!(
            absorbing.deposit_beam_power(cx, origin(), direction(), 2.0 * LENGTH, &[3.0], 0),
            Err(VolumeTransportError::Budget)
        );
        gate.request();
        assert_eq!(
            clear_field.deposit_beam_power(cx, origin(), direction(), f64::INFINITY, &[0.0], 0),
            Err(VolumeTransportError::Cancelled)
        );
    });
}

#[cfg(feature = "optical-heating")]
mod heating {
    use super::*;
    use crate::volumes::heating::OpticalHeating;
    use fs_thermal::enthalpy::{
        Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
        DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
    };

    struct ConstantCapacity;
    impl ThermalLaw for ConstantCapacity {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            if h.is_finite() && (0.0..=1e8).contains(&h) {
                Ok(300.0 + h / 1000.0)
            } else {
                Err(DepositionError::MaterialDomain)
            }
        }
        fn maximum_temperature_slope(&self) -> f64 {
            0.001
        }
    }
    static CAPACITY: ConstantCapacity = ConstantCapacity;

    #[test]
    fn g3_accepted_spectral_absorption_heats_original_cells_on_native_enthalpy_clock() {
        with_cx(|cx, _| {
            let (_, accepted, optical) = complex_fixture([0.01, 0.03]);
            let material =
                MaterialAbsorptionBinding::try_spectral(&accepted, &optical, &density_key())
                    .unwrap();
            let density = [1000.0; 2];
            let geometry = geometry(cx, &density, 16);
            let field = MaterialAbsorptionField::try_new(
                cx,
                &geometry,
                palette(cx, vec![material], &[0, 0]),
                [8000.0, 10000.0, 12000.0],
                MaterialAbsorptionEmission::None,
                budget(),
            )
            .unwrap();
            let report = field
                .deposit_beam_power(
                    cx,
                    origin(),
                    direction(),
                    f64::INFINITY,
                    &[0.002, 0.003, 0.005],
                    2,
                )
                .unwrap();
            let heat: [f64; 2] =
                std::array::from_fn(|cell| report.cells[cell].absorbed_power_w.iter().sum());
            let mass = 1000.0 * (4.0 * LENGTH); // Same optical density and physical cell volume.
            let events = [11_u64, 22]
                .into_iter()
                .enumerate()
                .map(|(index, id)| TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(DepositedCell {
                        id: CellId(id),
                        segment_id: id,
                        body_id: 1,
                        material_id: MaterialId(1),
                        mass_kg: mass,
                        specific_h_j_kg: 0.0,
                        geometry: CellGeometry {
                            start_m: [0.0, 0.0, index as f64 * LENGTH],
                            end_m: [0.0, 0.0, (index + 1) as f64 * LENGTH],
                            volume_m3: 4.0 * LENGTH,
                            surface_area_m2: 8.0 + 8.0 * LENGTH,
                        },
                        boundary: Boundary {
                            ambient_k: 300.0,
                            htc_w_m2_k: 0.0,
                            reservoir: None,
                            source_w: if index == 0 { 0.001 } else { 0.0 },
                        },
                    }),
                })
                .collect();
            let initial = DepositionSimulation::new(
                events,
                &[(MaterialId(1), &CAPACITY as &dyn ThermalLaw)],
                DepositionBudget {
                    max_cells: 2,
                    max_links: 0,
                    max_events: 2,
                    max_substeps: 100,
                },
            )
            .unwrap();
            let config = StepConfig {
                max_dt_s: 0.1,
                energy_tolerance_j: 1e-12,
            };
            let original = initial.advance_to(0.0, config, || Ok(())).unwrap().0;
            let observer = OpticalHeating::try_new(
                report,
                &[(0, CellId(11)), (1, CellId(22))],
                &original,
                1e-13,
                || Ok(()),
            )
            .unwrap();
            assert_eq!(observer.power().medium_identity, field.content_identity());
            assert_eq!(observer.power().wavelength_bounds_nm, (8000.0, 12000.0));
            assert_eq!(observer.power().scattered_w, [0.0; 3]);
            let (after, observed, audit) = original
                .advance_observed(&observer, 2.0, config, || Ok(()))
                .unwrap();
            for (index, id) in [CellId(11), CellId(22)].into_iter().enumerate() {
                let boundary = if index == 0 { 0.001 } else { 0.0 };
                let expected = 300.0 + 2.0 * (heat[index] + boundary) / (1000.0 * mass);
                assert!((after.temperature_k(id).unwrap() - expected).abs() < 2e-10);
                assert_eq!(original.temperature_k(id).unwrap(), 300.0);
            }
            close(audit.external_heat_j, 2.0 * (heat[0] + heat[1] + 0.001));
            close(observed.energy().thermal_input_j, 2.0 * (heat[0] + heat[1]));
            assert!(heat[0] > heat[1]);
        });
    }
}
