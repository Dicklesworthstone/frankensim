//! Shares and uncertainty of unchanged native histories on original Mie sources.

use super::super::*;
use crate::material_state::optical_mixture::deposition::ConstituentDepositionBudget;
use crate::volumes::transport::scattering::deposition::ScatteringDepositionConfig;

fn scattering_config(histories: u32) -> ScatteringDepositionConfig {
    ScatteringDepositionConfig {
        histories_per_lane: histories,
        max_collisions_per_history: 32,
        max_total_segments: u64::from(histories) * 32 * 3,
        max_deposition_cells: 2,
        seed: 83,
        source_id: 19,
    }
}

fn partition_budget() -> ConstituentDepositionBudget {
    ConstituentDepositionBudget {
        max_rows: 4,
        max_component_evaluations: 12,
    }
}

#[test]
fn g1_constituent_scattering_retains_native_histories_and_scales_cell_means_and_errors() {
    with_cx(|cx, _| {
        let a = particles(cx, 0, 300.0);
        let b = particles(cx, 1, 300.0);
        let components = [
            MaterialOpticalComponent::Particles(&a),
            MaterialOpticalComponent::Particles(&b),
        ];
        let composition = composition(&[0.25, 0.75]);
        let recipe =
            MaterialOpticalMixture::try_new(cx, &components, &composition, mixture_budget())
                .unwrap();
        let density = [0.003, 0.001];
        let field = MaterialOpticalMixtureGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [0.0; 3], [2.0, 2.0, 0.5]),
            palette(cx, vec![recipe], &[0, 0]),
            None,
            None,
            0.01,
            grid_budget(),
        )
        .unwrap();
        let wavelengths = [430.0, 650.0, 760.0];
        let powers = [1.0, 2.0, 0.0];
        let origin = Point3::new(1.0, 1.0, -0.25);
        let direction = Vec3::new(0.0, 0.0, 1.0);
        let config = scattering_config(256);
        let native = field
            .medium()
            .deposit_multiple_scattering_power(cx, origin, direction, &wavelengths, &powers, config)
            .unwrap();
        let partition = field
            .deposit_multiple_scattering_power_by_component(
                cx,
                origin,
                direction,
                &wavelengths,
                &powers,
                config,
                partition_budget(),
            )
            .unwrap();
        let replay = field
            .deposit_multiple_scattering_power_by_component(
                cx,
                origin,
                direction,
                &wavelengths,
                &powers,
                config,
                partition_budget(),
            )
            .unwrap();
        assert_eq!(partition, replay);
        assert_eq!(partition.scattering(), &native);
        assert_eq!(partition.mixture_identity(), field.content_identity());
        assert_eq!(native.sampled_histories, [256, 256, 0]);
        assert!(
            native.diffuse_escape.mean_w[..2]
                .iter()
                .any(|power| *power > 0.0)
        );
        assert!(native.traced_segments > 512);
        assert!(!partition.cells().is_empty());
        assert_eq!(
            partition.component_evaluations(),
            native.cells.len() * 2 * 3
        );
        let mut total = [0.0; 3];
        for row in partition.cells() {
            let native_cell = native
                .cells
                .iter()
                .find(|cell| cell.cell_index == row.optical_cell_index())
                .unwrap();
            let (particle, weight) = if row.component_index() == 0 {
                (&a, 0.25)
            } else {
                (&b, 0.75)
            };
            assert_eq!(row.source_identity(), particle.binding_identity());
            for lane in 0..3 {
                let denominator = 0.25 * a.spectrum().absorption.at(wavelengths[lane]).unwrap()
                    + 0.75 * b.spectrum().absorption.at(wavelengths[lane]).unwrap();
                let share = weight
                    * particle
                        .spectrum()
                        .absorption
                        .at(wavelengths[lane])
                        .unwrap()
                    / denominator;
                close(
                    row.absorbed_power_w()[lane],
                    share * native_cell.absorbed.mean_w[lane],
                );
                if native.sampled_histories[lane] == 0 {
                    assert_eq!(row.standard_error_w()[lane], Some(0.0));
                } else if native_cell.absorbed.mean_w[lane] == 0.0 {
                    assert_eq!(row.standard_error_w()[lane], None);
                } else {
                    close(
                        row.standard_error_w()[lane].unwrap(),
                        share * native_cell.absorbed.standard_error_w[lane],
                    );
                }
                total[lane] += row.absorbed_power_w()[lane];
            }
        }
        for lane in 0..3 {
            close(
                total[lane] + partition.partition_residual_w()[lane],
                native.absorbed.mean_w[lane],
            );
        }
    });
}

#[test]
fn g0_unobserved_absorbing_lanes_do_not_claim_a_zero_standard_error() {
    with_cx(|cx, _| {
        let a = particles(cx, 0, 300.0);
        let b = particles(cx, 1, 300.0);
        let components = [
            MaterialOpticalComponent::Particles(&a),
            MaterialOpticalComponent::Particles(&b),
        ];
        let composition = composition(&[0.4, 0.6]);
        let recipe =
            MaterialOpticalMixture::try_new(cx, &components, &composition, mixture_budget())
                .unwrap();
        let wavelengths = [430.0, 760.0];
        let largest_extinction = wavelengths
            .iter()
            .map(|&wavelength| {
                recipe.absorption_spectrum().at(wavelength).unwrap()
                    + recipe
                        .scattering_optics()
                        .scattering
                        .at(wavelength)
                        .unwrap()
            })
            .fold(0.0_f64, f64::max);
        let density = [0.4 / largest_extinction];
        let field = MaterialOpticalMixtureGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0; 3]),
            palette(cx, vec![recipe], &[0]),
            None,
            None,
            0.01,
            grid_budget(),
        )
        .unwrap();
        // Select a deterministic small-sample fixture with an unobserved lane.
        // This tests representation of missing evidence, not a confidence bound
        // or a precision-driven stopping rule for a production measurement.
        let mut found = false;
        for seed in 0..64 {
            let config = ScatteringDepositionConfig {
                seed,
                ..scattering_config(2)
            };
            let report = field
                .deposit_multiple_scattering_power_by_component(
                    cx,
                    Point3::new(0.5, 0.5, -0.25),
                    Vec3::new(0.0, 0.0, 1.0),
                    &wavelengths,
                    &[1.0, 1.0],
                    config,
                    partition_budget(),
                )
                .unwrap();
            for cell in report.scattering().cells.iter() {
                for lane in 0..2 {
                    if cell.absorbed.mean_w[lane] == 0.0 {
                        assert_eq!(report.scattering().sampled_histories[lane], 2);
                        for row in report
                            .cells()
                            .iter()
                            .filter(|row| row.optical_cell_index() == cell.cell_index)
                        {
                            assert_eq!(row.absorbed_power_w()[lane], 0.0);
                            assert_eq!(row.standard_error_w()[lane], None);
                            found = true;
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        assert!(
            found,
            "the bounded deterministic fixture set must exercise an unobserved absorbing lane"
        );
    });
}

#[test]
fn g1_pure_molecular_scattering_has_exactly_zero_constituent_heat() {
    with_cx(|cx, _| {
        let state = gas(41).evaluate(300.0, 100_000.0).unwrap();
        let density = [state.gas().density];
        let carrier = MaterialGasScatteringBinding::try_new(state.clone()).unwrap();
        let components = [MaterialOpticalComponent::Gas(&carrier)];
        let composition = composition(&[1.0]);
        let recipe =
            MaterialOpticalMixture::try_new(cx, &components, &composition, mixture_budget())
                .unwrap();
        let states = [state];
        let field = MaterialOpticalMixtureGrid::try_new(
            cx,
            grid(&density),
            palette(cx, vec![recipe], &[0]),
            Some(&states),
            None,
            0.01,
            grid_budget(),
        )
        .unwrap();
        let wavelengths = [420.0, 650.0];
        let powers = [1.0, 2.0];
        let config = scattering_config(32);
        let origin = Point3::new(0.5, 0.5, -1.0);
        let direction = Vec3::new(0.0, 0.0, 1.0);
        let native = field
            .medium()
            .deposit_multiple_scattering_power(cx, origin, direction, &wavelengths, &powers, config)
            .unwrap();
        let report = field
            .deposit_multiple_scattering_power_by_component(
                cx,
                origin,
                direction,
                &wavelengths,
                &powers,
                config,
                ConstituentDepositionBudget {
                    max_rows: 0,
                    max_component_evaluations: 0,
                },
            )
            .unwrap();
        assert_eq!(report.scattering(), &native);
        assert_eq!(native.absorbed.mean_w, [0.0; 2]);
        assert_eq!(native.absorbed.standard_error_w, [0.0; 2]);
        assert!(report.cells().is_empty());
        assert_eq!(report.component_evaluations(), 0);
        assert_eq!(report.partition_residual_w(), &[0.0; 2]);
    });
}

#[test]
fn g0_constituent_scattering_respects_native_and_partition_budgets_and_cancellation() {
    with_cx(|cx, gate| {
        let a = particles(cx, 0, 300.0);
        let b = particles(cx, 1, 300.0);
        let components = [
            MaterialOpticalComponent::Particles(&a),
            MaterialOpticalComponent::Particles(&b),
        ];
        let composition = composition(&[0.5, 0.5]);
        let recipe =
            MaterialOpticalMixture::try_new(cx, &components, &composition, mixture_budget())
                .unwrap();
        let density = [0.1];
        let field = MaterialOpticalMixtureGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [0.0; 3], [2.0; 3]),
            palette(cx, vec![recipe], &[0]),
            None,
            None,
            0.01,
            grid_budget(),
        )
        .unwrap();
        let origin = Point3::new(1.0, 1.0, -0.25);
        let direction = Vec3::new(0.0, 0.0, 1.0);
        let run = |config, budget| {
            field.deposit_multiple_scattering_power_by_component(
                cx,
                origin,
                direction,
                &[550.0],
                &[1.0],
                config,
                budget,
            )
        };
        for budget in [
            ConstituentDepositionBudget {
                max_rows: 1,
                max_component_evaluations: 2,
            },
            ConstituentDepositionBudget {
                max_rows: 2,
                max_component_evaluations: 1,
            },
        ] {
            assert!(matches!(
                run(scattering_config(8), budget),
                Err(MaterialOpticalMixtureError::Transport(
                    VolumeTransportError::Budget
                ))
            ));
        }
        assert!(matches!(
            run(
                ScatteringDepositionConfig {
                    max_total_segments: 0,
                    ..scattering_config(8)
                },
                partition_budget()
            ),
            Err(MaterialOpticalMixtureError::Transport(
                VolumeTransportError::Budget
            ))
        ));
        gate.request();
        assert!(matches!(
            run(scattering_config(8), partition_budget()),
            Err(MaterialOpticalMixtureError::Transport(
                VolumeTransportError::Cancelled
            ))
        ));
    });
}

#[cfg(feature = "optical-heating")]
mod thermal_bridge {
    use super::*;
    use crate::material_state::optical_mixture::deposition::heating::ConstituentHeatingError;
    use crate::volumes::heating::OpticalHeatingError;
    use crate::volumes::heating::scattering::ScatteringHeatingTolerance;
    use fs_thermal::enthalpy::{
        Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
        DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
    };

    struct ConstantCapacity;
    impl ThermalLaw for ConstantCapacity {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            if h.is_finite() && (0.0..=1e9).contains(&h) {
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

    fn thermal() -> DepositionSimulation<'static> {
        let events = [CellId(11), CellId(22)]
            .into_iter()
            .map(|id| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: id.0,
                    body_id: id.0,
                    material_id: MaterialId(1),
                    mass_kg: 1.0,
                    specific_h_j_kg: 0.0,
                    geometry: CellGeometry {
                        start_m: [0.0; 3],
                        end_m: [0.0, 0.0, 1.0],
                        volume_m3: 1e-3,
                        surface_area_m2: 1.0,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            })
            .collect();
        DepositionSimulation::new(
            events,
            &[(MaterialId(1), &CAPACITY as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: 2,
                max_links: 0,
                max_events: 2,
                max_substeps: 128,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 0.1,
                energy_tolerance_j: 1e-12,
            },
            || Ok(()),
        )
        .unwrap()
        .0
    }

    fn tolerance<const N: usize>(
        report: &crate::material_state::optical_mixture::deposition::scattering::ConstituentScatteringDeposition<N>,
    ) -> ScatteringHeatingTolerance {
        let native = report.scattering();
        let unresolved = native.unresolved.mean_w.iter().sum::<f64>();
        let total_error = native
            .absorbed
            .standard_error_w
            .iter()
            .chain(native.unresolved.standard_error_w.iter())
            .copied()
            .sum::<f64>();
        let cell_error = report
            .cells()
            .iter()
            .flat_map(|row| row.standard_error_w().iter().copied().flatten())
            .sum::<f64>();
        ScatteringHeatingTolerance {
            power_roundoff_w: 1e-10,
            max_unresolved_w: unresolved + 1e-12,
            max_total_standard_error_w: total_error + 1e-12,
            max_cell_standard_error_w: cell_error + 1e-12,
        }
    }

    #[test]
    fn g1_scattered_constituent_heat_advances_distinct_native_thermal_cells() {
        with_cx(|cx, _| {
            let a = particles(cx, 0, 300.0);
            let b = particles(cx, 1, 300.0);
            let components = [
                MaterialOpticalComponent::Particles(&a),
                MaterialOpticalComponent::Particles(&b),
            ];
            let composition = composition(&[0.3, 0.7]);
            let recipe =
                MaterialOpticalMixture::try_new(cx, &components, &composition, mixture_budget())
                    .unwrap();
            let density = [0.004, 0.002];
            let field = MaterialOpticalMixtureGrid::try_new(
                cx,
                VolumeGrid::new([1, 1, 2], &density, [0.0; 3], [2.0, 2.0, 0.5]),
                palette(cx, vec![recipe], &[0, 0]),
                None,
                None,
                0.01,
                grid_budget(),
            )
            .unwrap();
            let report = field
                .deposit_multiple_scattering_power_by_component(
                    cx,
                    Point3::new(1.0, 1.0, -0.25),
                    Vec3::new(0.0, 0.0, 1.0),
                    &[430.0, 650.0],
                    &[0.01, 0.02],
                    scattering_config(1024),
                    ConstituentDepositionBudget {
                        max_rows: 4,
                        max_component_evaluations: 8,
                    },
                )
                .unwrap();
            assert!(report.cells().iter().all(|row| {
                row.standard_error_w().iter().all(Option::is_some)
            }));
            let tol = tolerance(&report);
            let bindings: Vec<_> = report
                .cells()
                .iter()
                .map(|row| {
                    (
                        row.optical_cell_index(),
                        row.component_index(),
                        if row.component_index() == 0 {
                            CellId(11)
                        } else {
                            CellId(22)
                        },
                    )
                })
                .collect();
            let initial = thermal();
            let mut heating = report
                .into_heating(&bindings, &initial, tol, || Ok(()))
                .unwrap();
            let source_a = heating.sources_w()[&CellId(11)];
            let source_b = heating.sources_w()[&CellId(22)];
            assert!(source_a > 0.0 && source_b > 0.0);
            let native_absorbed = heating.deposition().absorbed.mean_w.iter().sum::<f64>();
            close(source_a + source_b, native_absorbed);
            let (after, _, audit) = initial
                .advance_observed(
                    &mut heating,
                    1.0,
                    StepConfig {
                        max_dt_s: 0.1,
                        energy_tolerance_j: 1e-12,
                    },
                    || Ok(()),
                )
                .unwrap();
            close(after.temperature_k(CellId(11)).unwrap(), 300.0 + source_a / 1000.0);
            close(after.temperature_k(CellId(22)).unwrap(), 300.0 + source_b / 1000.0);
            close(audit.external_heat_j, source_a + source_b);
            close(heating.energy().thermal_input_j, source_a + source_b);
            assert!(heating.energy().diffuse_escape_j.iter().any(|v| *v > 0.0));
        });
    }

    #[test]
    fn g0_unobserved_or_incomplete_constituent_heat_refuses_before_native_advance() {
        with_cx(|cx, _| {
            let a = particles(cx, 0, 300.0);
            let b = particles(cx, 1, 300.0);
            let components = [
                MaterialOpticalComponent::Particles(&a),
                MaterialOpticalComponent::Particles(&b),
            ];
            let composition = composition(&[0.4, 0.6]);
            let recipe =
                MaterialOpticalMixture::try_new(cx, &components, &composition, mixture_budget())
                    .unwrap();
            let wavelengths = [430.0, 760.0];
            let largest_extinction = wavelengths
                .iter()
                .map(|&wavelength| {
                    recipe.absorption_spectrum().at(wavelength).unwrap()
                        + recipe
                            .scattering_optics()
                            .scattering
                            .at(wavelength)
                            .unwrap()
                })
                .fold(0.0_f64, f64::max);
            let density = [0.4 / largest_extinction];
            let field = MaterialOpticalMixtureGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0; 3]),
                palette(cx, vec![recipe], &[0]),
                None,
                None,
                0.01,
                grid_budget(),
            )
            .unwrap();
            let initial = thermal();
            let mut found = false;
            for seed in 0..64 {
                let report = field
                    .deposit_multiple_scattering_power_by_component(
                        cx,
                        Point3::new(0.5, 0.5, -0.25),
                        Vec3::new(0.0, 0.0, 1.0),
                        &wavelengths,
                        &[1.0, 1.0],
                        ScatteringDepositionConfig {
                            seed,
                            ..scattering_config(2)
                        },
                        partition_budget(),
                    )
                    .unwrap();
                if report
                    .cells()
                    .iter()
                    .any(|row| row.standard_error_w().iter().any(Option::is_none))
                {
                    let tol = tolerance(&report);
                    let bindings: Vec<_> = report
                        .cells()
                        .iter()
                        .map(|row| {
                            (
                                row.optical_cell_index(),
                                row.component_index(),
                                CellId(11),
                            )
                        })
                        .collect();
                    assert!(matches!(
                        report.into_heating(&bindings, &initial, tol, || Ok(())),
                        Err(ConstituentHeatingError::Heating(
                            OpticalHeatingError::Invalid(_)
                        ))
                    ));
                    found = true;
                    break;
                }
            }
            assert!(found);
            assert_eq!(initial.temperature_k(CellId(11)).unwrap(), 300.0);
            assert_eq!(initial.temperature_k(CellId(22)).unwrap(), 300.0);
        });
    }
}
