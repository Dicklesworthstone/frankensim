//! An appended elastic fiber must not bypass the original material domains.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::LongitudinalStringResolution;
use fs_couple::fiber_ribbon::{FiberRibbonGeometry, FiberRibbonResolution, TensionedFiberRibbon};
use fs_couple::render::schedule::force::coupled::contact::friction::{
    multiple::MultiFrictionConfig, nonlinear::NonlinearModalConfig,
};

fn elastic_fibers(gate: &CancelGate) -> TensionedFiberRibbon {
    // Manufactured source, deliberately without hereditary fiber memory.
    // Both fiber losses are explicitly zero in the existing material fixture.
    let section = ConcentricElasticSection::try_new(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: 40e-6,
            material: material(0.0, 0.004, 200_000.0),
        }],
    )
    .unwrap();
    TensionedFiberRibbon::try_new(
        FiberRibbonGeometry {
            length_m: 0.6,
            width_m: 0.01,
            fiber_count: 180,
            tension_n: 40.0,
        },
        FiberRibbonResolution {
            cells: 1,
            axial_modes: 2,
            transverse_modes: 1,
        },
        section.clone().with_axial_viscosity(ETA).unwrap(),
        section.with_kelvin_voigt_bending(ETA).unwrap(),
        48_000,
        ModalAcousticTimeBudget::audible_reference(),
        gate,
    )
    .unwrap()
}

#[test]
fn elastic_fibers_revalidate_original_bending_layered_and_torsional_memory() {
    with_cx(|cx, gate| {
        let fibers = elastic_fibers(gate);
        assert!(fibers.material_relaxation().is_none());
        let joint = || MultiFrictionConfig {
            max_ports: 1,
            max_sweeps: 128,
            max_setup_terms: 1_048_576,
        };
        // Exercise both source-validation paths, and each existing string
        // history owner. The restricted band is derived from the actual
        // original and enlarged operators, not an assumed modal frequency.
        for longitudinal in [false, true] {
            for kind in 0..3 {
                let make = |band| {
                    let original = match kind {
                        0 => bank_with_sources(
                            cx,
                            gate,
                            vec![original_source(40e9, 0.004, band, true)],
                        )
                        .unwrap(),
                        1 => bank(cx, gate, 40e9, 0.004, band, 2, 1)
                            .with_material_string_bending_relaxation(selection(), None, 16, gate)
                            .unwrap(),
                        _ => {
                            bank_with_sources(cx, gate, vec![shear_source(40e9, 0.004, band, true)])
                                .unwrap()
                        }
                    };
                    if longitudinal {
                        original
                            .with_longitudinal_waves(
                                vec![LongitudinalStringResolution {
                                    axial_modes: vec![1],
                                    maximum_slope: 0.2,
                                    maximum_strain: 0.01,
                                }],
                                NonlinearModalConfig {
                                    max_iterations: 64,
                                    maximum_generalized_force: 1e5,
                                    force_absolute_tolerance: 1e-11,
                                    force_relative_tolerance: 1e-10,
                                },
                                gate,
                            )
                            .unwrap()
                    } else {
                        original
                    }
                };
                let original = make(200_000.0);
                let original_band = original.system().linear_frequency_upper_bound_rad_s();
                let memory = original.system().material_memory().to_vec();
                assert!(!memory.is_empty());
                let admitted = original
                    .with_tensioned_fiber_ribbon(fibers.clone(), 0.4, joint(), gate)
                    .unwrap();
                let enlarged_band = admitted.system().linear_frequency_upper_bound_rad_s();
                assert!(enlarged_band > original_band);
                assert_eq!(admitted.system().material_memory(), memory);
                assert_eq!(admitted.system().samples_rendered(), 0);

                let restricted_band = 0.5 * (original_band + enlarged_band);
                let restricted = make(restricted_band);
                assert_eq!(
                    restricted.system().linear_frequency_upper_bound_rad_s(),
                    original_band,
                );
                match restricted.with_tensioned_fiber_ribbon(fibers.clone(), 0.4, joint(), gate) {
                    Ok(_) => panic!(
                        "elastic fibers escaped original material support: kind={kind}, longitudinal={longitudinal}",
                    ),
                    Err(error) => assert!(error.to_string().contains("support"), "{error}"),
                }
            }
        }
    });
}
