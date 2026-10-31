//! Independent source sums and native transport for unequal constituent temperatures.

use super::{ThermalEmissionComponent, ThermalEmissionModel};
use crate::thermal_surface::planck_radiance_w_m2_sr_nm;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::scattering::{
    CellPhaseFunction, CellScatteringGrid, MassScatteringSpectrum, PhaseFunction, ScatteringGrid,
};
use crate::volumes::transport::{
    AbsorbingGrid, AbsorbingGridBudget, CellPalette, CellPaletteBudget, MassAbsorptionSpectrum,
    OpticalMixtureBudget, VolumeTransportError,
};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use std::sync::Arc;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 43,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}

fn budget() -> AbsorbingGridBudget {
    AbsorbingGridBudget {
        max_grid_cells: 8,
        max_steps_per_segment: 16,
    }
}

fn model(cx: &Cx<'_>, components: &[ThermalEmissionComponent]) -> ThermalEmissionModel {
    ThermalEmissionModel::try_new(cx, components, OpticalMixtureBudget { max_components: 8 })
        .unwrap()
}

fn term(weight: f64, absorption: f64, temperature_k: f64) -> ThermalEmissionComponent {
    ThermalEmissionComponent {
        mass_weight: weight,
        absorption: MassAbsorptionSpectrum::constant(absorption).unwrap(),
        temperature_k: Some(temperature_k),
    }
}

fn palette(
    cx: &Cx<'_>,
    models: Vec<ThermalEmissionModel>,
    indices: &[u32],
) -> CellPalette<ThermalEmissionModel> {
    CellPalette::try_new(
        cx,
        models,
        Arc::from(indices),
        CellPaletteBudget {
            max_models: 8,
            max_cells: 8,
        },
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-100),
        "{actual:e} != {expected:e}"
    );
}

// Independent SI Planck equation; these moderate temperatures/wavelengths do
// not require the production owner's logarithmic stabilization.
fn planck(temperature: f64, wavelength_nm: f64) -> f64 {
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0;
    let k = 1.380_649e-23;
    let wavelength = wavelength_nm * 1e-9;
    2.0 * h * c * c / wavelength.powi(5) / ((h * c / (wavelength * k * temperature)).exp_m1())
        * 1e-9
}

#[test]
fn g1_two_temperature_slab_retains_each_original_spectral_source() {
    with_cx(|cx, _| {
        let source = model(
            cx,
            &[
                ThermalEmissionComponent {
                    mass_weight: 0.25,
                    absorption: MassAbsorptionSpectrum::table(&[(380.0, 0.4), (780.0, 1.2)])
                        .unwrap(),
                    temperature_k: Some(900.0),
                },
                ThermalEmissionComponent {
                    mass_weight: 0.75,
                    absorption: MassAbsorptionSpectrum::table(&[(380.0, 2.0), (780.0, 0.4)])
                        .unwrap(),
                    temperature_k: Some(1800.0),
                },
            ],
        );
        let density = [2.0];
        let sources = palette(cx, vec![source.clone()], &[0]);
        let medium = AbsorbingGrid::try_new_thermal_sources(
            cx,
            VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0, 1.0, 0.3]),
            sources,
            budget(),
        )
        .unwrap();
        assert!(medium.has_thermal_emission());
        assert!(medium.has_emission_proposal());
        assert!(core::ptr::eq(medium.grid.data.as_ptr(), density.as_ptr()));
        assert!(Arc::ptr_eq(
            medium.thermal_sources().unwrap().cell_indices(),
            medium.absorption_palette().unwrap().cell_indices()
        ));
        let wavelengths = [380.0, 523.0, 780.0];
        let result = medium
            .transport(
                cx,
                Point3::new(0.5, 0.5, -1.0),
                Vec3::new(0.0, 0.0, 1.0),
                f64::INFINITY,
                &wavelengths,
            )
            .unwrap();
        for (lane, &wavelength) in wavelengths.iter().enumerate() {
            let first = 0.25 * (0.4 + 0.002 * (wavelength - 380.0));
            let second = 0.75 * (2.0 - 0.004 * (wavelength - 380.0));
            let absorption = first + second;
            let emitted = first * planck(900.0, wavelength) + second * planck(1800.0, wavelength);
            let transmission = (-0.6 * absorption).exp();
            close(
                source.absorption_spectrum().at(wavelength).unwrap(),
                absorption,
            );
            close(
                source.source_radiance_at(wavelength).unwrap(),
                emitted / absorption,
            );
            close(result.transmittance[lane], transmission);
            close(
                result.emission_w_m2_sr_nm[lane],
                (1.0 - transmission) * emitted / absorption,
            );
        }
        // A source change must not be hidden by unchanged total opacity.
        let warmer = model(cx, &[term(0.25, 0.4, 1000.0), term(0.75, 2.0, 1800.0)]);
        let colder = model(cx, &[term(0.25, 0.4, 900.0), term(0.75, 2.0, 1800.0)]);
        assert_eq!(warmer.absorption_spectrum(), colder.absorption_spectrum());
        assert_ne!(warmer.content_identity(), colder.content_identity());
    });
}

#[test]
fn g1_layer_order_attenuates_each_constituent_source_with_the_actual_foreground() {
    with_cx(|cx, _| {
        let density = [1.0, 2.0];
        let medium = AbsorbingGrid::try_new_thermal_sources(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [0.0; 3], [1.0, 1.0, 0.4]),
            palette(
                cx,
                vec![
                    model(cx, &[term(0.4, 0.5, 1000.0), term(0.6, 2.0, 1600.0)]),
                    model(cx, &[term(0.6, 1.5, 2100.0), term(0.4, 0.25, 800.0)]),
                ],
                &[0, 1],
            ),
            budget(),
        )
        .unwrap();
        for (z, direction, order) in [(-1.0, 1.0, [0, 1]), (2.0, -1.0, [1, 0])] {
            let wavelengths = [550.0, 700.0];
            let result = medium
                .transport(
                    cx,
                    Point3::new(0.5, 0.5, z),
                    Vec3::new(0.0, 0.0, direction),
                    f64::INFINITY,
                    &wavelengths,
                )
                .unwrap();
            for (lane, &wavelength) in wavelengths.iter().enumerate() {
                let source = [
                    (0.2 * planck(1000.0, wavelength) + 1.2 * planck(1600.0, wavelength)) / 1.4,
                    0.9 * planck(2100.0, wavelength) + 0.1 * planck(800.0, wavelength),
                ];
                let tau = [0.56_f64, 0.8];
                let (near, far) = (order[0], order[1]);
                let expected = -(-tau[near]).exp_m1() * source[near]
                    + (-tau[near]).exp() * -(-tau[far]).exp_m1() * source[far];
                close(result.emission_w_m2_sr_nm[lane], expected);
                close(result.transmittance[lane], (-1.36_f64).exp());
            }
        }
    });
}

#[test]
fn g3_equal_constituent_temperatures_reduce_to_the_original_lte_source() {
    with_cx(|cx, _| {
        let source = model(
            cx,
            &[
                ThermalEmissionComponent {
                    mass_weight: 0.25,
                    absorption: MassAbsorptionSpectrum::table(&[(380.0, 0.3), (780.0, 1.1)])
                        .unwrap(),
                    temperature_k: Some(1800.0),
                },
                term(0.75, 0.7, 1800.0),
                ThermalEmissionComponent {
                    mass_weight: 0.2,
                    absorption: MassAbsorptionSpectrum::constant(0.0).unwrap(),
                    temperature_k: None,
                },
            ],
        );
        let density = [0.1, 1.0, 2.0];
        let temperature = [1800.0; 3];
        let field = || VolumeGrid::new([1, 1, 3], &density, [0.0; 3], [1.0, 1.0, 0.2]);
        let old = AbsorbingGrid::try_new(
            cx,
            field(),
            source.absorption_spectrum(),
            Some(&temperature),
            budget(),
        )
        .unwrap();
        let new = AbsorbingGrid::try_new_thermal_sources(
            cx,
            field(),
            palette(cx, vec![source.clone()], &[0, 0, 0]),
            budget(),
        )
        .unwrap();
        for wavelength in [380.0, 497.0, 780.0] {
            assert_eq!(
                source.source_radiance_at(wavelength).unwrap(),
                planck_radiance_w_m2_sr_nm(1800.0, wavelength).unwrap()
            );
            let evaluate = |grid: &AbsorbingGrid<'_>| {
                grid.transport(
                    cx,
                    Point3::new(0.5, 0.5, -1.0),
                    Vec3::new(0.0, 0.0, 1.0),
                    f64::INFINITY,
                    &[wavelength],
                )
                .unwrap()
            };
            let (a, b) = (evaluate(&old), evaluate(&new));
            assert_eq!(a.transmittance, b.transmittance);
            close(a.emission_w_m2_sr_nm[0], b.emission_w_m2_sr_nm[0]);
        }
    });
}

#[test]
fn g1_scattering_free_flight_retains_the_complete_unequal_temperature_source() {
    with_cx(|cx, _| {
        let density = [2.0];
        let absorption = AbsorbingGrid::try_new_thermal_sources(
            cx,
            VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0, 1.0, 0.5]),
            palette(
                cx,
                vec![model(cx, &[term(0.3, 0.4, 800.0), term(0.7, 1.2, 1700.0)])],
                &[0],
            ),
            budget(),
        )
        .unwrap();
        let medium = ScatteringGrid::try_new(
            cx,
            absorption,
            MassScatteringSpectrum::constant(2.0).unwrap(),
            PhaseFunction::henyey_greenstein(0.6).unwrap(),
        )
        .unwrap();
        let wavelengths = [420.0, 560.0, 760.0];
        let origin = Point3::new(0.5, 0.5, -0.2);
        let direction = Vec3::new(0.0, 0.0, 1.0);
        let transported = medium
            .transport(cx, origin, direction, 1.0, &wavelengths)
            .unwrap();
        let sampled = medium
            .sample_segment(cx, origin, direction, 1.0, &wavelengths, 0.1, 0.99)
            .unwrap();
        assert!(sampled.event.is_some());
        assert_eq!(sampled.emission_w_m2_sr_nm, transported.emission_w_m2_sr_nm);
        for (lane, &wavelength) in wavelengths.iter().enumerate() {
            let emitted = 0.12 * planck(800.0, wavelength) + 0.84 * planck(1700.0, wavelength);
            close(transported.transmittance[lane], (-2.96_f64).exp());
            close(
                transported.emission_w_m2_sr_nm[lane],
                -(-2.96_f64).exp_m1() * emitted / 2.96,
            );
            close(sampled.weight[lane], 2.0 / 2.96);
        }
    });
}

#[test]
fn g0_prepared_packets_keep_matched_constituent_absorption_and_refuse_substitution() {
    with_cx(|cx, _| {
        let density = [1.0, 2.0];
        let geometry = AbsorbingGrid::try_new_thermal_sources(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [0.0; 3], [1.0, 1.0, 0.5]),
            palette(
                cx,
                vec![
                    model(cx, &[term(0.5, 0.2, 800.0), term(0.5, 1.0, 1900.0)]),
                    model(cx, &[term(1.0, 0.4, 1100.0)]),
                ],
                &[0, 1],
            ),
            budget(),
        )
        .unwrap();
        let wavelengths = [450.0, 700.0];
        let absorption: [[f64; 2]; 2] =
            std::array::from_fn(|cell| geometry.coefficients_at_cell(cell, &wavelengths).unwrap());
        let scattering = [[0.8; 2], [0.3; 2]];
        let phase = PhaseFunction::Isotropic;
        let field = CellScatteringGrid::try_new(
            cx,
            &geometry,
            wavelengths,
            &absorption,
            &scattering,
            CellPhaseFunction::Uniform(&phase),
            (380.0, 780.0),
            ContentHash([7; 32]),
        )
        .unwrap();
        let transported = field
            .transport(
                cx,
                Point3::new(0.5, 0.5, -1.0),
                Vec3::new(0.0, 0.0, 1.0),
                f64::INFINITY,
            )
            .unwrap();
        for (lane, &wavelength) in wavelengths.iter().enumerate() {
            let first = (0.1 * planck(800.0, wavelength) + 0.5 * planck(1900.0, wavelength)) / 1.4;
            let second = 0.4 * planck(1100.0, wavelength) / 0.7;
            close(
                transported.emission_w_m2_sr_nm[lane],
                -(-0.7_f64).exp_m1() * (first + (-0.7_f64).exp() * second),
            );
            close(transported.transmittance[lane], (-1.4_f64).exp());
        }
        let mut unrelated = absorption;
        unrelated[1][0] += 0.1;
        assert!(matches!(
            CellScatteringGrid::try_new(
                cx,
                &geometry,
                wavelengths,
                &unrelated,
                &scattering,
                CellPhaseFunction::Uniform(&phase),
                (380.0, 780.0),
                ContentHash([7; 32])
            ),
            Err(VolumeTransportError::ThermalSource)
        ));
        assert!(matches!(
            CellScatteringGrid::try_new(
                cx,
                &geometry,
                wavelengths,
                &absorption,
                &scattering,
                CellPhaseFunction::Uniform(&phase),
                (300.0, 900.0),
                ContentHash([7; 32])
            ),
            Err(VolumeTransportError::ThermalSource)
        ));
    });
}

#[test]
fn g0_source_admission_requires_each_absorbers_temperature_and_binds_original_terms() {
    with_cx(|cx, gate| {
        let options = OpticalMixtureBudget { max_components: 8 };
        let absorbing = MassAbsorptionSpectrum::constant(1.0).unwrap();
        let transparent = MassAbsorptionSpectrum::constant(0.0).unwrap();
        let unknown = ThermalEmissionComponent {
            mass_weight: 1.0,
            absorption: absorbing,
            temperature_k: None,
        };
        assert_eq!(
            ThermalEmissionModel::try_new(cx, &[unknown], options),
            Err(VolumeTransportError::Temperature)
        );
        for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let invalid = ThermalEmissionComponent {
                temperature_k: Some(temperature),
                ..unknown
            };
            assert_eq!(
                ThermalEmissionModel::try_new(cx, &[invalid], options),
                Err(VolumeTransportError::Temperature)
            );
        }
        for inactive in [
            ThermalEmissionComponent {
                mass_weight: 0.0,
                ..unknown
            },
            ThermalEmissionComponent {
                absorption: transparent,
                ..unknown
            },
        ] {
            let admitted = model(cx, &[inactive]);
            assert_eq!(admitted.maximum_temperature_k(), None);
            assert_eq!(admitted.maximum_absorption_m2_kg(), 0.0);
            assert_eq!(admitted.source_radiance_at(550.0), Ok(0.0));
        }
        let original = [
            term(1.0, 1.0, 1500.0),
            ThermalEmissionComponent {
                mass_weight: 0.0,
                absorption: transparent,
                temperature_k: Some(9000.0),
            },
        ];
        let base = model(cx, &original);
        assert_eq!(base.components(), &original);
        assert_eq!(base.maximum_temperature_k(), Some(1500.0));
        let mut changed = original;
        changed[1].temperature_k = Some(5000.0);
        let changed_temperature = model(cx, &changed);
        changed[1].mass_weight = 0.25;
        let changed_weight = model(cx, &changed);
        assert_ne!(
            base.content_identity(),
            changed_temperature.content_identity()
        );
        assert_ne!(
            changed_temperature.content_identity(),
            changed_weight.content_identity()
        );
        assert_eq!(
            base.source_radiance_at(550.0),
            changed_temperature.source_radiance_at(550.0)
        );
        assert_eq!(
            base.source_radiance_at(550.0),
            changed_weight.source_radiance_at(550.0)
        );
        assert_eq!(
            base.source_radiance_at(800.0),
            Err(VolumeTransportError::Wavelength)
        );
        assert_eq!(
            ThermalEmissionModel::try_new(cx, &[], options),
            Err(VolumeTransportError::Mixture)
        );
        assert_eq!(
            ThermalEmissionModel::try_new(
                cx,
                &original,
                OpticalMixtureBudget { max_components: 1 }
            ),
            Err(VolumeTransportError::Budget)
        );
        let foreign = ThermalEmissionComponent {
            mass_weight: 1.0,
            absorption: MassAbsorptionSpectrum::source_table(
                &[(250.0, 1.0), (780.0, 1.0)],
                ContentHash([3; 32]),
            )
            .unwrap(),
            temperature_k: Some(1000.0),
        };
        assert_eq!(
            ThermalEmissionModel::try_new(cx, &[original[0], foreign], options),
            Err(VolumeTransportError::Wavelength)
        );
        gate.request();
        assert_eq!(
            ThermalEmissionModel::try_new(cx, &original, options),
            Err(VolumeTransportError::Cancelled)
        );
    });
}
