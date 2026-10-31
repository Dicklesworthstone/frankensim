//! Original material fixtures, exact spectral shares and native thermal storage.

use super::*;
use crate::material_state::optical_mixture::deposition::{
    ConstituentBeamDeposition, ConstituentDepositionBudget,
};

fn with_two_populations(cx: &Cx<'_>, run: impl FnOnce(&MaterialOpticalMixtureGrid<'_>)) {
    let a = particles(cx, 0, 300.0);
    let b = particles(cx, 1, 300.0);
    let components = [
        MaterialOpticalComponent::Particles(&a),
        MaterialOpticalComponent::Particles(&b),
    ];
    let first = composition(&[0.25, 0.75]);
    let second = composition(&[0.7, 0.3]);
    let recipes = vec![
        MaterialOpticalMixture::try_new(cx, &components, &first, mixture_budget()).unwrap(),
        MaterialOpticalMixture::try_new(cx, &components, &second, mixture_budget()).unwrap(),
    ];
    let density = [1e-4, 2e-4];
    let medium = MaterialOpticalMixtureGrid::try_new(
        cx,
        VolumeGrid::new([1, 1, 2], &density, [0.0; 3], [1.0; 3]),
        palette(cx, recipes, &[1, 0]),
        None,
        None,
        0.01,
        grid_budget(),
    )
    .unwrap();
    run(&medium);
}

fn partition(
    cx: &Cx<'_>,
    medium: &MaterialOpticalMixtureGrid<'_>,
    power: &[f64; 3],
    budget: ConstituentDepositionBudget,
) -> Result<ConstituentBeamDeposition<3>, MaterialOpticalMixtureError> {
    medium.deposit_beam_power_by_component(
        cx,
        Point3::new(0.5, 0.5, -1.0),
        Vec3::new(0.0, 0.0, 1.0),
        f64::INFINITY,
        &[380.0, 550.0, 780.0],
        power,
        2,
        budget,
    )
}

fn partition_budget() -> ConstituentDepositionBudget {
    ConstituentDepositionBudget {
        max_rows: 4,
        max_component_evaluations: 12,
    }
}

#[test]
fn original_constituent_spectra_partition_one_unchanged_native_beam() {
    with_cx(|cx, _| {
        with_two_populations(cx, |medium| {
            let report = partition(cx, medium, &[2.0, 3.0, 5.0], partition_budget()).unwrap();
            let reference = medium
                .medium()
                .deposit_beam_power(
                    cx,
                    Point3::new(0.5, 0.5, -1.0),
                    Vec3::new(0.0, 0.0, 1.0),
                    f64::INFINITY,
                    &[380.0, 550.0, 780.0],
                    &[2.0, 3.0, 5.0],
                    2,
                )
                .unwrap();
            assert_eq!(report.beam(), &reference);
            assert_eq!(report.mixture_identity(), medium.content_identity());
            assert_eq!(report.cells().len(), 4);
            assert_eq!(report.component_evaluations(), 12);
            let mut reduced = [0.0; 3];
            for row in report.cells() {
                let cell = row.optical_cell_index();
                let mixture = medium.mixtures().at_cell(cell).unwrap();
                let original = &reference.cells[cell];
                assert_eq!(original.cell_index, cell);
                let MaterialOpticalComponent::Particles(owner) =
                    mixture.components()[row.component_index()]
                else {
                    panic!("fixture contains original particle owners");
                };
                assert_eq!(row.source_identity(), owner.binding_identity());
                for (lane, wavelength) in [380.0, 550.0, 780.0].into_iter().enumerate() {
                    let mut denominator = 0.0;
                    for (component, weight) in mixture
                        .components()
                        .iter()
                        .zip(mixture.composition().fractions())
                    {
                        let MaterialOpticalComponent::Particles(particles) = component else {
                            unreachable!()
                        };
                        denominator +=
                            weight * particles.spectrum().absorption.at(wavelength).unwrap();
                    }
                    let expected = original.absorbed_power_w[lane]
                        * mixture.composition().fractions()[row.component_index()]
                        * owner.spectrum().absorption.at(wavelength).unwrap()
                        / denominator;
                    close(row.absorbed_power_w()[lane], expected);
                    reduced[lane] += row.absorbed_power_w()[lane];
                }
            }
            for lane in 0..3 {
                close(reduced[lane], reference.absorbed_power_w[lane]);
                assert!(report.partition_residual_w()[lane].abs() < 1e-12);
            }
            let first_component = &report.cells()[0];
            let blue_share =
                first_component.absorbed_power_w()[0] / reference.cells[0].absorbed_power_w[0];
            let red_share =
                first_component.absorbed_power_w()[2] / reference.cells[0].absorbed_power_w[2];
            assert!(
                (blue_share - red_share).abs() > 1e-6,
                "absorption shares must retain spectral dependence"
            );
            assert!(
                (blue_share - 0.7).abs() > 1e-3,
                "heat must not be split by mass alone"
            );
        })
    });
}

#[test]
fn pure_molecular_scattering_and_zero_incident_power_create_no_constituent_heat() {
    with_cx(|cx, _| {
        let state = gas(3).evaluate(300.0, 100_000.0).unwrap();
        let gas = MaterialGasScatteringBinding::try_new(state).unwrap();
        let components = [MaterialOpticalComponent::Gas(&gas)];
        let fractions = composition(&[1.0]);
        let recipe =
            MaterialOpticalMixture::try_new(cx, &components, &fractions, mixture_budget()).unwrap();
        let density = [state.gas().density];
        let states = [state];
        let medium = MaterialOpticalMixtureGrid::try_new(
            cx,
            grid(&density),
            palette(cx, vec![recipe], &[0]),
            Some(&states),
            None,
            0.01,
            grid_budget(),
        )
        .unwrap();
        let report = partition(
            cx,
            &medium,
            &[2.0, 3.0, 5.0],
            ConstituentDepositionBudget {
                max_rows: 0,
                max_component_evaluations: 0,
            },
        )
        .unwrap();
        assert!(
            report
                .beam()
                .scattered_power_w
                .iter()
                .all(|power| *power > 0.0)
        );
        assert_eq!(report.beam().absorbed_power_w, [0.0; 3]);
        assert!(report.cells().is_empty());
        assert_eq!(report.component_evaluations(), 0);

        with_two_populations(cx, |medium| {
            let zero = partition(
                cx,
                medium,
                &[0.0; 3],
                ConstituentDepositionBudget {
                    max_rows: 0,
                    max_component_evaluations: 0,
                },
            )
            .unwrap();
            assert!(zero.cells().is_empty());
            assert!(zero.beam().cells.is_empty());
            assert_eq!(zero.partition_residual_w(), &[0.0; 3]);
        });
    });
}

#[test]
fn constituent_output_work_and_cancellation_refuse_before_publishing() {
    with_cx(|cx, gate| {
        with_two_populations(cx, |medium| {
            for budget in [
                ConstituentDepositionBudget {
                    max_rows: 3,
                    max_component_evaluations: 12,
                },
                ConstituentDepositionBudget {
                    max_rows: 4,
                    max_component_evaluations: 11,
                },
            ] {
                assert!(matches!(
                    partition(cx, medium, &[2.0, 3.0, 5.0], budget),
                    Err(MaterialOpticalMixtureError::Transport(
                        VolumeTransportError::Budget
                    ))
                ));
            }
            gate.request();
            assert!(matches!(
                partition(cx, medium, &[2.0, 3.0, 5.0], partition_budget()),
                Err(MaterialOpticalMixtureError::Transport(
                    VolumeTransportError::Cancelled
                ))
            ));
        })
    });
}

#[cfg(feature = "optical-heating")]
mod heating {
    use super::*;
    use crate::material_state::optical_mixture::deposition::heating::ConstituentHeatingError;
    use crate::volumes::heating::OpticalHeatingError;
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
    fn separate_constituents_heat_separate_native_states_and_accumulate_original_cells() {
        with_cx(|cx, _| {
            with_two_populations(cx, |medium| {
                let report =
                    partition(cx, medium, &[0.002, 0.003, 0.005], partition_budget()).unwrap();
                let mut masses = [0.0; 2];
                for (cell, density) in medium.density_kg_m3().iter().enumerate() {
                    for (component, weight) in medium
                        .mixtures()
                        .at_cell(cell)
                        .unwrap()
                        .composition()
                        .fractions()
                        .iter()
                        .enumerate()
                    {
                        masses[component] += density * weight; // Original one-m³ optical cells.
                    }
                }
                let events = [11_u64, 22]
                    .into_iter()
                    .enumerate()
                    .map(|(component, id)| TimedEvent {
                        time_s: 0.0,
                        event: Event::Birth(DepositedCell {
                            id: CellId(id),
                            segment_id: id,
                            body_id: id,
                            material_id: MaterialId(1),
                            mass_kg: masses[component],
                            specific_h_j_kg: 0.0,
                            geometry: CellGeometry {
                                start_m: [0.0; 3],
                                end_m: [0.0, 0.0, 2.0],
                                volume_m3: masses[component] / (1000.0 * (component + 1) as f64),
                                surface_area_m2: 1.0,
                            },
                            boundary: Boundary {
                                ambient_k: 300.0,
                                htc_w_m2_k: 0.0,
                                reservoir: None,
                                source_w: if component == 0 { 0.001 } else { 0.0 },
                            },
                        }),
                    })
                    .collect();
                let config = StepConfig {
                    max_dt_s: 0.1,
                    energy_tolerance_j: 1e-12,
                };
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
                let original = initial.advance_to(0.0, config, || Ok(())).unwrap().0;
                let bindings = [
                    (0, 0, CellId(11)),
                    (0, 1, CellId(22)),
                    (1, 0, CellId(11)),
                    (1, 1, CellId(22)),
                ];
                let heating = report
                    .prepare_heating(&bindings, &original, 1e-13, || Ok(()))
                    .unwrap();
                let heat = [
                    heating.sources_w()[&CellId(11)],
                    heating.sources_w()[&CellId(22)],
                ];
                assert!(heat.iter().all(|power| *power > 0.0));
                let (after, accepted, audit) = original
                    .advance_observed(&heating, 2.0, config, || Ok(()))
                    .unwrap();
                for (component, id) in [CellId(11), CellId(22)].into_iter().enumerate() {
                    let boundary = if component == 0 { 0.001 } else { 0.0 };
                    let expected =
                        300.0 + 2.0 * (heat[component] + boundary) / (1000.0 * masses[component]);
                    assert!((after.temperature_k(id).unwrap() - expected).abs() < 2e-10);
                    assert_eq!(original.temperature_k(id).unwrap(), 300.0);
                }
                assert_ne!(
                    after.temperature_k(CellId(11)).unwrap(),
                    after.temperature_k(CellId(22)).unwrap()
                );
                close(audit.external_heat_j, 2.0 * (heat[0] + heat[1] + 0.001));
                close(accepted.energy().thermal_input_j, 2.0 * (heat[0] + heat[1]));
                close(
                    accepted.power().thermal_input_w,
                    report.beam().absorbed_power_w.iter().sum(),
                );
                assert!(
                    accepted
                        .power()
                        .scattered_w
                        .iter()
                        .any(|power| *power > 0.0)
                );

                assert!(matches!(
                    report.prepare_heating(&bindings[..3], &original, 1e-13, || Ok(())),
                    Err(ConstituentHeatingError::MissingBinding {
                        optical_cell: 1,
                        component: 1
                    })
                ));
                let duplicate = [bindings[0], bindings[1], bindings[2], bindings[2]];
                assert!(matches!(
                    report.prepare_heating(&duplicate, &original, 1e-13, || Ok(())),
                    Err(ConstituentHeatingError::DuplicateBinding {
                        optical_cell: 1,
                        component: 0
                    })
                ));
                assert!(matches!(
                    report.prepare_heating(&bindings, &original, 1e-13, || Err(
                        DepositionError::Cancelled
                    )),
                    Err(ConstituentHeatingError::Heating(
                        OpticalHeatingError::Thermal(DepositionError::Cancelled)
                    ))
                ));
            })
        });
    }
}
