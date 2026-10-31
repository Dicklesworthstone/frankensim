use super::super::tests::{budget, close, with_cx};
use super::*;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::CellPaletteBudget;
use std::sync::Arc;

fn palette<T>(cx: &Cx<'_>, values: Vec<T>, map: &[u32]) -> CellPalette<T> {
    CellPalette::try_new(
        cx,
        values,
        Arc::from(map),
        CellPaletteBudget {
            max_models: 4,
            max_cells: 4,
        },
    )
    .unwrap()
}

#[test]
fn g1_local_spectra_and_phase_match_independent_layer_collision_integrals() {
    with_cx(|cx, _| {
        let densities = [2.0, 1.0];
        let a = [
            MassAbsorptionSpectrum::table(&[(380.0, 0.2), (780.0, 0.6)]).unwrap(),
            MassAbsorptionSpectrum::table(&[(380.0, 0.8), (780.0, 0.1)]).unwrap(),
        ];
        let s = [
            MassScatteringSpectrum::table(&[(380.0, 0.1), (780.0, 1.4)]).unwrap(),
            MassScatteringSpectrum::table(&[(380.0, 1.1), (780.0, 0.3)]).unwrap(),
        ];
        let medium = ScatteringGrid::try_new_indexed(
            cx,
            AbsorbingGrid::try_new_indexed(
                cx,
                VolumeGrid::new([1, 1, 2], &densities, [0.0; 3], [1.0, 1.0, 0.4]),
                // Absorption and scattering deliberately have different maps.
                palette(cx, a.to_vec(), &[1, 0]),
                None,
                budget(8),
            )
            .unwrap(),
            palette(
                cx,
                vec![
                    ScatteringOptics {
                        scattering: s[0],
                        phase: PhaseFunction::HenyeyGreenstein { g: -0.6 },
                    },
                    ScatteringOptics {
                        scattering: s[1],
                        phase: PhaseFunction::HenyeyGreenstein { g: 0.7 },
                    },
                ],
                &[0, 1],
            ),
        )
        .unwrap();
        assert!(medium.phase_ref().is_none());
        let origin = Point3::new(0.5, 0.5, -1.0);
        let direction = Vec3::new(0.0, 0.0, 1.0);
        let wavelengths = [380.0, 580.0, 780.0];
        let transmitted = medium
            .transmittance(cx, origin, direction, 3.0, &wavelengths)
            .unwrap();
        let mut expected_scatter = [0.0; 3];
        let mut expected_angular = [0.0; 3];
        for (lane, &wavelength) in wavelengths.iter().enumerate() {
            let a = [a[1].at(wavelength).unwrap(), a[0].at(wavelength).unwrap()];
            let s = [s[0].at(wavelength).unwrap(), s[1].at(wavelength).unwrap()];
            let extinction = [a[0] + s[0], a[1] + s[1]];
            let tau = [0.8 * extinction[0], 0.4 * extinction[1]];
            close(transmitted[lane], (-(tau[0] + tau[1])).exp(), 3e-14);
            let probabilities = [
                s[0] / extinction[0] * -(-tau[0]).exp_m1(),
                (-tau[0]).exp() * s[1] / extinction[1] * -(-tau[1]).exp_m1(),
            ];
            expected_scatter[lane] = probabilities.iter().sum();
            for (cell, g) in [-0.6_f64, 0.7].into_iter().enumerate() {
                let phase = (1.0 - g * g)
                    / (4.0 * core::f64::consts::PI * (1.0 + g * g - 0.6 * g).powf(1.5));
                expected_angular[lane] += probabilities[cell] * phase;
            }
        }
        let (mut escaped, mut scattered, mut angular) = ([0.0; 3], [0.0; 3], [0.0; 3]);
        const SAMPLES: usize = 16384;
        for hero in 0..3 {
            for sample in 0..SAMPLES {
                let flight = medium
                    .sample_segment(
                        cx,
                        origin,
                        direction,
                        3.0,
                        &wavelengths,
                        (hero as f64 + 0.5) / 3.0,
                        (sample as f64 + 0.5) / SAMPLES as f64,
                    )
                    .unwrap();
                if let Some(event) = flight.event {
                    assert_eq!(event.cell_index, usize::from(event.position.z >= 0.4));
                    let phase = medium.phase_at_cell(event.cell_index).unwrap();
                    for lane in 0..3 {
                        scattered[lane] += flight.weight[lane];
                        angular[lane] +=
                            flight.weight[lane] * phase.eval_at(wavelengths[lane], 0.3).unwrap();
                    }
                } else {
                    for lane in 0..3 {
                        escaped[lane] += flight.weight[lane];
                    }
                }
            }
        }
        for lane in 0..3 {
            let count = (3 * SAMPLES) as f64;
            close(escaped[lane] / count, transmitted[lane], 4e-4);
            close(scattered[lane] / count, expected_scatter[lane], 4e-4);
            close(angular[lane] / count, expected_angular[lane], 6e-4);
        }
    });
}

#[test]
fn g0_indexed_uniform_reduction_preserves_free_flight_and_complete_lte_source() {
    with_cx(|cx, gate| {
        let density = [0.5, 2.0];
        let temperatures = [1600.0, 2100.0];
        let grid = || VolumeGrid::new([1, 1, 2], &density, [0.0; 3], [1.0; 3]);
        let absorption = MassAbsorptionSpectrum::constant(0.3).unwrap();
        let scattering = MassScatteringSpectrum::constant(0.7).unwrap();
        let optics = ScatteringOptics {
            scattering,
            phase: PhaseFunction::Rayleigh,
        };
        let uniform = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(cx, grid(), absorption, Some(&temperatures), budget(8)).unwrap(),
            scattering,
            PhaseFunction::Rayleigh,
        )
        .unwrap();
        let indexed = ScatteringGrid::try_new_indexed(
            cx,
            AbsorbingGrid::try_new_indexed(
                cx,
                grid(),
                palette(cx, vec![absorption; 2], &[1, 0]),
                Some(&temperatures),
                budget(8),
            )
            .unwrap(),
            palette(cx, vec![optics.clone(), optics], &[0, 1]),
        )
        .unwrap();
        let o = Point3::new(0.5, 0.5, -1.0);
        let d = Vec3::new(0.0, 0.0, 1.0);
        let wavelengths = [380.0, 550.0, 780.0];
        assert_eq!(
            uniform.transport(cx, o, d, 4.0, &wavelengths),
            indexed.transport(cx, o, d, 4.0, &wavelengths)
        );
        for u in [0.0, 0.1, 0.4, 0.8, 1.0 - f64::EPSILON] {
            assert_eq!(
                uniform.sample_segment(cx, o, d, 4.0, &wavelengths, 0.5, u),
                indexed.sample_segment(cx, o, d, 4.0, &wavelengths, 0.5, u)
            );
        }
        assert_eq!(indexed.phase_at_cell(2), Err(VolumeTransportError::Palette));
        assert_ne!(uniform.content_identity(), indexed.content_identity());
        gate.request();
        assert_eq!(
            indexed.transmittance(cx, o, d, 4.0, &wavelengths),
            Err(VolumeTransportError::Cancelled)
        );
    });
}
