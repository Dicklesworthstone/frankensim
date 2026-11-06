//! Independent geometric densities, support and bounded refusal checks.

use super::*;
use crate::volumes::transport::{
    AbsorbingGridBudget, CellPaletteBudget, MassAbsorptionSpectrum, OpticalMixtureBudget,
    ThermalEmissionComponent,
};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use std::sync::Arc;

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        operation(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 731,
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

fn grid<'a>(
    cx: &Cx<'_>,
    field: VolumeGrid<'a>,
    temperature: Option<&'a [f64]>,
    absorption: f64,
    steps: u32,
) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(
        cx,
        field,
        MassAbsorptionSpectrum::constant(absorption).unwrap(),
        temperature,
        AbsorbingGridBudget {
            max_grid_cells: 64,
            max_steps_per_segment: steps,
        },
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual / expected - 1.0).abs() <= relative,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn thermal_model(cx: &Cx<'_>, terms: &[(f64, f64, Option<f64>)]) -> ThermalEmissionModel {
    let components: Vec<_> = terms
        .iter()
        .map(
            |&(mass_weight, absorption, temperature_k)| ThermalEmissionComponent {
                mass_weight,
                absorption: MassAbsorptionSpectrum::constant(absorption).unwrap(),
                temperature_k,
            },
        )
        .collect();
    ThermalEmissionModel::try_new(cx, &components, OpticalMixtureBudget { max_components: 4 })
        .unwrap()
}

fn source_palette(
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

fn indexed_grid<'a>(
    cx: &Cx<'_>,
    field: VolumeGrid<'a>,
    temperatures: &'a [f64],
    models: Vec<MassAbsorptionSpectrum>,
    indices: &[u32],
) -> AbsorbingGrid<'a> {
    let absorption = CellPalette::try_new(
        cx,
        models,
        Arc::from(indices),
        CellPaletteBudget {
            max_models: 8,
            max_cells: 8,
        },
    )
    .unwrap();
    AbsorbingGrid::try_new_indexed(
        cx,
        field,
        absorption,
        Some(temperatures),
        AbsorbingGridBudget {
            max_grid_cells: 8,
            max_steps_per_segment: 16,
        },
    )
    .unwrap()
}

#[test]
fn g1_indexed_proposal_weights_absorption_and_never_samples_hot_transparent_cells() {
    with_cx(|cx, _| {
        let density = [1e300, 1.0, 1e300, 1.0, 1e300, 0.0];
        let temperature = [f64::MAX, 2000.0, f64::MAX, 2000.0, f64::MAX, f64::MAX];
        let medium = indexed_grid(
            cx,
            VolumeGrid::new([1, 1, 6], &density, [-0.5, -0.5, 0.0], [1.0; 3]),
            &temperature,
            vec![
                MassAbsorptionSpectrum::constant(0.0).unwrap(),
                MassAbsorptionSpectrum::table(&[(380.0, 0.0), (600.0, 1.0), (780.0, 0.0)]).unwrap(),
                MassAbsorptionSpectrum::table(&[(380.0, 0.0), (600.0, 2.0), (780.0, 0.0)]).unwrap(),
                MassAbsorptionSpectrum::constant(1e300).unwrap(),
            ],
            &[0, 1, 0, 2, 0, 3],
        );
        // Both active spectra vanish at the source-band Wien maximum (780 nm).
        // Their interior knots still emit, with absorption maxima in ratio 1:2.
        // Huge transparent temperatures must neither overflow admission nor
        // take any of the support mixture away from these actual absorbers.
        let probability = [0.0, 11.0 / 32.0, 0.0, 21.0 / 32.0, 0.0, 0.0];
        for (cell, &expected) in probability.iter().enumerate() {
            if expected == 0.0 {
                assert_eq!(medium.emission.probability(cell), 0.0);
            } else {
                near(medium.emission.probability(cell), expected, 2e-13);
            }
        }
        let mut count = [0usize; 6];
        for index in 0..4096 {
            let direction = medium
                .sample_emission_direction(
                    cx,
                    Point3::new(-2.0, 0.0, 0.0),
                    [(f64::from(index) + 0.5) / 4096.0, 0.5, 0.5, 0.5],
                )
                .unwrap()
                .unwrap();
            let cell = (2.0 * direction.z / direction.x).floor() as usize;
            count[cell] += 1;
        }
        assert_eq!(count, [0, 1408, 0, 2688, 0, 0]);
        let origin = Point3::new(0.0, 0.0, -1.0);
        let direction = Vec3::new(0.0, 0.0, 1.0);
        // The two unit-volume source cells span radii [2,3] and [4,5].
        // Thus q = (11/32)*(27-8)/3 + (21/32)*(125-64)/3.
        near(
            medium
                .emission_direction_pdf(cx, origin, direction)
                .unwrap(),
            745.0 / 48.0,
            2e-13,
        );
        assert_eq!(
            medium
                .emission_direction_pdf(cx, Point3::new(0.0, -1.0, 0.5), Vec3::new(0.0, 1.0, 0.0),)
                .unwrap(),
            0.0,
        );
        // The ordinary thermal solution remains Beer--Lambert/Kirchhoff: the
        // same-temperature emitting cells have total depth 1+2 at 600 nm.
        let transport = medium
            .transport(cx, origin, direction, f64::INFINITY, &[600.0, 780.0])
            .unwrap();
        let lambda_m = 600e-9_f64;
        let planck = 2.0 * 6.626_070_15e-34 * 299_792_458.0_f64.powi(2) * 1e-9
            / (lambda_m.powi(5)
                * (6.626_070_15e-34 * 299_792_458.0 / (lambda_m * 1.380_649e-23 * 2000.0))
                    .exp_m1());
        near(transport.transmittance[0], (-3.0_f64).exp(), 3e-14);
        near(
            transport.emission_w_m2_sr_nm[0],
            planck * -(-3.0_f64).exp_m1(),
            2e-13,
        );
        assert_eq!(transport.transmittance[1], 1.0);
        assert_eq!(transport.emission_w_m2_sr_nm[1], 0.0);
    });
}

#[test]
fn g1_indexed_direction_density_has_exact_octant_shape_and_unit_normalization() {
    with_cx(|cx, _| {
        let density = [1.0; 8];
        let temperature = [2000.0; 8];
        let medium = indexed_grid(
            cx,
            VolumeGrid::new([2; 3], &density, [-1.0; 3], [1.0; 3]),
            &temperature,
            vec![
                MassAbsorptionSpectrum::constant(0.0).unwrap(),
                MassAbsorptionSpectrum::constant(1.0).unwrap(),
                MassAbsorptionSpectrum::constant(3.0).unwrap(),
            ],
            &[0, 1, 0, 0, 0, 0, 2, 0],
        );
        let probability = [0.0, 17.0 / 64.0, 0.0, 0.0, 0.0, 0.0, 47.0 / 64.0, 0.0];
        assert_eq!(medium.emission.cdf.last(), Some(&1.0));
        let mut integral = 0.0;
        const SAMPLES: usize = 16_384;
        for index in 0..SAMPLES {
            let z = 1.0 - 2.0 * (index as f64 + 0.5) / SAMPLES as f64;
            let phi = index as f64 * (std::f64::consts::PI * (3.0 - 5.0_f64.sqrt()));
            let radius = (1.0 - z * z).sqrt();
            let direction = Vec3::new(radius * phi.cos(), radius * phi.sin(), z);
            let cell = usize::from(direction.x > 0.0)
                + 2 * usize::from(direction.y > 0.0)
                + 4 * usize::from(direction.z > 0.0);
            // Each unit cube owns one octant from the receiver at its corner.
            // Its exact radial Jacobian is r_exit^3/3, with no occlusion term.
            let expected = probability[cell]
                / (3.0
                    * direction
                        .x
                        .abs()
                        .max(direction.y.abs())
                        .max(direction.z.abs())
                        .powi(3));
            let actual = medium
                .emission_direction_pdf(cx, Point3::new(0.0, 0.0, 0.0), direction)
                .unwrap();
            if expected == 0.0 {
                assert_eq!(actual, 0.0);
            } else {
                near(actual, expected, 2e-13);
            }
            integral += actual;
        }
        near(
            integral * (4.0 * std::f64::consts::PI / SAMPLES as f64),
            1.0,
            2e-3,
        );
    });
}

#[test]
fn g3_indexed_absorption_ratios_retain_cold_infrared_support_without_extrapolation() {
    with_cx(|cx, _| {
        let density = [1e-300, 1e300, 1e300];
        let temperature = [20.0, 10_000.0, f64::MAX];
        let spectrum = |coefficient| {
            MassAbsorptionSpectrum::source_table(
                &[(2000.0, coefficient), (150_000.0, coefficient)],
                fs_blake3::ContentHash([19; 32]),
            )
            .unwrap()
        };
        let medium = indexed_grid(
            cx,
            VolumeGrid::new([1, 1, 3], &density, [0.0; 3], [1.0; 3]),
            &temperature,
            vec![spectrum(1e-300), spectrum(1e300), spectrum(0.0)],
            &[0, 1, 2],
        );
        // Density*absorption would underflow/overflow before normalization.
        // Log weights plus a support mixture preserve the cold IR contributor.
        assert_eq!(medium.emission.probability(0), 1.0 / 32.0);
        assert_eq!(medium.emission.probability(1), 31.0 / 32.0);
        assert_eq!(medium.emission.probability(2), 0.0);
        let origin = Point3::new(0.5, -1.0, 0.5);
        let direction = medium
            .sample_emission_direction(cx, origin, [0.01, 0.5, 0.5, 0.5])
            .unwrap()
            .unwrap();
        assert_eq!(direction, Vec3::new(0.0, 1.0, 0.0));
        near(
            medium
                .emission_direction_pdf(cx, origin, direction)
                .unwrap(),
            7.0 / 96.0,
            3e-14,
        );
        assert_eq!(
            medium.transmittance(cx, origin, direction, 3.0, &[550.0]),
            Err(VolumeTransportError::Wavelength),
        );
        // Cold-only IR sources also retain a proposal, unlike their visible
        // tails. The nonzero opacity ratio still controls its source sampling.
        let cold_density = [1.0; 2];
        let cold_temperature = [20.0; 2];
        let cold = indexed_grid(
            cx,
            VolumeGrid::new([1, 1, 2], &cold_density, [0.0; 3], [1.0; 3]),
            &cold_temperature,
            vec![spectrum(1.0), spectrum(2.0)],
            &[0, 1],
        );
        near(cold.emission.probability(0), 11.0 / 32.0, 2e-13);
        near(cold.emission.probability(1), 21.0 / 32.0, 2e-13);
        assert!(
            cold.transport(cx, origin, direction, 3.0, &[150_000.0])
                .unwrap()
                .emission_w_m2_sr_nm[0]
                > 0.0,
        );
    });
}

#[test]
fn g1_uniform_absorption_preserves_temperature_field_proposal_bits() {
    with_cx(|cx, _| {
        let density = [1.0, 2.0, 0.01];
        let temperature = [1200.0, 2000.0, 10_000.0];
        let field = || VolumeGrid::new([1, 1, 3], &density, [0.0; 3], [1.0; 3]);
        let control = grid(cx, field(), Some(&temperature), 1.0, 8);
        for absorption in [1e-300, 0.4, 1e300] {
            let medium = grid(cx, field(), Some(&temperature), absorption, 8);
            assert_eq!(medium.emission.cdf, control.emission.cdf);
        }
    });
}

#[test]
fn g1_cube_density_matches_radial_jacobian_and_extreme_metric_rescaling() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperature = [2000.0];
        let diagonal = 1.0 / 3.0_f64.sqrt();
        for scale in [1.0e-120, 1.0, 1.0e120] {
            let medium = grid(
                cx,
                VolumeGrid::new([1; 3], &density, [-scale; 3], [2.0 * scale; 3]),
                Some(&temperature),
                0.4,
                8,
            );
            // Uniform source cube: q(w)=r_exit³/(3V), so axis=1/24 and
            // body-diagonal=sqrt(3)/8. These are independent of length units.
            near(
                medium
                    .emission_direction_pdf(
                        cx,
                        Point3::new(0.0, 0.0, 0.0),
                        Vec3::new(1.0, 0.0, 0.0),
                    )
                    .unwrap(),
                1.0 / 24.0,
                3.0e-12,
            );
            near(
                medium
                    .emission_direction_pdf(
                        cx,
                        Point3::new(0.0, 0.0, 0.0),
                        Vec3::new(diagonal, diagonal, diagonal),
                    )
                    .unwrap(),
                3.0_f64.sqrt() / 8.0,
                3.0e-12,
            );
            // Outside the cube, (4³-2³)/(3*8)=7/3; the entry distance
            // cannot be dropped from the marginal radial probability.
            near(
                medium
                    .emission_direction_pdf(
                        cx,
                        Point3::new(0.0, 0.0, 3.0 * scale),
                        Vec3::new(0.0, 0.0, -1.0),
                    )
                    .unwrap(),
                7.0 / 3.0,
                3.0e-12,
            );
            let sampled = medium
                .sample_emission_direction(cx, Point3::new(0.0, 0.0, 0.0), [0.2, 0.75, 0.75, 0.75])
                .unwrap()
                .unwrap();
            near(sampled.x, diagonal, 2.0e-14);
            near(sampled.dot(sampled), 1.0, 2.0e-14);
        }
    });
}

#[test]
fn g1_direction_density_normalizes_and_survives_uniform_grid_refinement() {
    with_cx(|cx, _| {
        let density = [1.0; 8];
        let temperature = [2000.0; 8];
        let coarse = grid(
            cx,
            VolumeGrid::new([1; 3], &density[..1], [-1.0; 3], [2.0; 3]),
            Some(&temperature[..1]),
            0.4,
            8,
        );
        let fine = grid(
            cx,
            VolumeGrid::new([2; 3], &density, [-1.0; 3], [1.0; 3]),
            Some(&temperature),
            0.4,
            8,
        );
        let origin = Point3::new(0.0, 0.0, 0.0);
        let mut integral = 0.0;
        const SAMPLES: usize = 16_384;
        for index in 0..SAMPLES {
            // Independent equal-area Fibonacci sphere quadrature.
            let z = 1.0 - 2.0 * (index as f64 + 0.5) / SAMPLES as f64;
            let phi = index as f64 * (std::f64::consts::PI * (3.0 - 5.0_f64.sqrt()));
            let radius = (1.0 - z * z).sqrt();
            let direction = Vec3::new(radius * phi.cos(), radius * phi.sin(), z);
            let analytic = 1.0
                / (24.0
                    * direction
                        .x
                        .abs()
                        .max(direction.y.abs())
                        .max(direction.z.abs())
                        .powi(3));
            let value = coarse
                .emission_direction_pdf(cx, origin, direction)
                .unwrap();
            near(value, analytic, 3.0e-14);
            near(
                fine.emission_direction_pdf(cx, origin, direction).unwrap(),
                value,
                3.0e-14,
            );
            integral += value;
        }
        near(
            integral * (4.0 * std::f64::consts::PI / SAMPLES as f64),
            1.0,
            1.0e-3,
        );
    });
}

#[test]
fn g1_sample_frequencies_and_full_ray_pdf_share_the_same_cell_mixture() {
    with_cx(|cx, _| {
        let density = [1.0, 2.0, 0.0];
        let temperature = [2000.0; 3];
        let medium = grid(
            cx,
            VolumeGrid::new([1, 1, 3], &density, [-0.5, -0.5, 0.0], [1.0; 3]),
            Some(&temperature),
            0.4,
            8,
        );
        // 15/16 density-proportional plus 1/16 uniform gives 11/32 and 21/32.
        let mut count = [0usize; 2];
        for index in 0..4096 {
            let direction = medium
                .sample_emission_direction(
                    cx,
                    Point3::new(-2.0, 0.0, 0.0),
                    [(index as f64 + 0.5) / 4096.0, 0.5, 0.5, 0.5],
                )
                .unwrap()
                .unwrap();
            count[usize::from(direction.z / direction.x > 0.5)] += 1;
        }
        assert_eq!(count, [1408, 2688]);
        // Both source cells contribute, even though the first lies in front
        // of the second. Geometric visibility is deliberately not a PDF input.
        near(
            medium
                .emission_direction_pdf(cx, Point3::new(0.0, 0.0, -1.0), Vec3::new(0.0, 0.0, 1.0))
                .unwrap(),
            119.0 / 24.0,
            3.0e-14,
        );
    });
}

#[test]
fn g3_extreme_source_ratios_keep_cold_spectral_tail_support() {
    with_cx(|cx, _| {
        let density = [1.0e-300, 1.0e300];
        let temperature = [50.0, 10_000.0];
        let medium = grid(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [0.0; 3], [1.0; 3]),
            Some(&temperature),
            1.0e300,
            8,
        );
        // The importance ratio underflows, but the uniform component retains
        // probability 1/32 for the cold, positive-density cell.
        let origin = Point3::new(0.5, -1.0, 0.5);
        let direction = medium
            .sample_emission_direction(cx, origin, [0.01, 0.5, 0.5, 0.5])
            .unwrap()
            .unwrap();
        assert_eq!(direction, Vec3::new(0.0, 1.0, 0.0));
        near(
            medium
                .emission_direction_pdf(cx, origin, direction)
                .unwrap(),
            7.0 / 96.0,
            3.0e-14,
        );
    });
}

#[test]
fn g4_empty_proposals_invalid_samples_traversal_budgets_and_cancel_are_explicit() {
    with_cx(|cx, gate| {
        let density = [1.0; 2];
        let hot = [2000.0; 2];
        let cold = [1.0e-6; 2];
        let origin = Point3::new(0.5, 0.5, -1.0);
        let direction = Vec3::new(0.0, 0.0, 1.0);
        for (absorption, temperature) in
            [(0.0, Some(&hot[..])), (1.0, None), (1.0, Some(&cold[..]))]
        {
            let medium = grid(
                cx,
                VolumeGrid::new([1, 1, 2], &density, [0.0; 3], [1.0; 3]),
                temperature,
                absorption,
                8,
            );
            assert!(!medium.has_emission_proposal());
            assert_eq!(
                medium
                    .sample_emission_direction(cx, origin, [0.5; 4])
                    .unwrap(),
                None
            );
            assert_eq!(
                medium
                    .emission_direction_pdf(cx, origin, direction)
                    .unwrap(),
                0.0
            );
        }
        let medium = grid(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [0.0; 3], [1.0; 3]),
            Some(&hot),
            1.0,
            1,
        );
        assert_eq!(
            medium.sample_emission_direction(cx, origin, [1.0, 0.5, 0.5, 0.5]),
            Err(VolumeTransportError::Sample)
        );
        assert_eq!(
            medium.emission_direction_pdf(cx, origin, direction.scale(2.0)),
            Err(VolumeTransportError::Ray)
        );
        assert_eq!(
            medium.emission_direction_pdf(cx, origin, direction),
            Err(VolumeTransportError::Budget)
        );
        gate.request();
        assert_eq!(
            medium.sample_emission_direction(cx, origin, [0.5; 4]),
            Err(VolumeTransportError::Cancelled)
        );
        assert_eq!(
            medium.emission_direction_pdf(cx, origin, direction),
            Err(VolumeTransportError::Cancelled)
        );
    });
}

#[test]
fn g1_constituent_proposal_uses_absorption_and_preserves_the_complete_direction_pdf() {
    with_cx(|cx, _| {
        let density = [1.0, 1.0, 0.0, 1.0];
        let sources = source_palette(
            cx,
            vec![
                thermal_model(cx, &[(0.5, 1.0, Some(1000.0)), (0.5, 3.0, Some(2000.0))]),
                thermal_model(cx, &[(1.0, 4.0, Some(2000.0))]),
                thermal_model(cx, &[(1.0, 0.0, None)]),
                thermal_model(cx, &[(1.0, 100.0, Some(10_000.0))]),
            ],
            &[0, 1, 3, 2],
        );
        let field = || VolumeGrid::new([1, 1, 4], &density, [-0.5, -0.5, 0.0], [1.0; 3]);
        let distribution = EmissionDistribution::try_new_sources(cx, &field(), &sources).unwrap();
        // The two emitting cells have the same maximum temperature and total
        // absorption 2 and 4. Thus the exact proposal is 15/16*(1/3,2/3)+1/32.
        near(distribution.probability(0), 11.0 / 32.0, 2e-13);
        near(distribution.probability(1), 21.0 / 32.0, 2e-13);
        assert_eq!(distribution.probability(2), 0.0);
        assert_eq!(distribution.probability(3), 0.0);
        let medium = AbsorbingGrid::try_new_thermal_sources(
            cx,
            field(),
            sources,
            AbsorbingGridBudget {
                max_grid_cells: 8,
                max_steps_per_segment: 8,
            },
        )
        .unwrap();
        assert!(medium.has_emission_proposal());
        near(
            medium
                .emission_direction_pdf(cx, Point3::new(0.0, 0.0, -1.0), Vec3::new(0.0, 0.0, 1.0))
                .unwrap(),
            119.0 / 24.0,
            2e-13,
        );
        let mut count = [0_usize; 2];
        for index in 0..4096 {
            let direction = medium
                .sample_emission_direction(
                    cx,
                    Point3::new(-2.0, 0.0, 0.0),
                    [(f64::from(index) + 0.5) / 4096.0, 0.5, 0.5, 0.5],
                )
                .unwrap()
                .unwrap();
            count[usize::from(direction.z / direction.x > 0.5)] += 1;
        }
        assert_eq!(count, [1408, 2688]);
    });
}

#[test]
fn g3_constituent_proposal_keeps_cold_support_and_excludes_unoccupied_sources() {
    with_cx(|cx, gate| {
        let models = vec![
            thermal_model(cx, &[(1.0, 1e-300, Some(50.0))]),
            thermal_model(cx, &[(1.0, 1e300, Some(10_000.0))]),
        ];
        let sources = source_palette(cx, models, &[0, 1]);
        let density = [1e-300, 1e300];
        let field = VolumeGrid::new([1, 1, 2], &density, [0.0; 3], [1.0; 3]);
        let distribution = EmissionDistribution::try_new_sources(cx, &field, &sources).unwrap();
        assert_eq!(distribution.probability(0), 1.0 / 32.0);
        assert_eq!(distribution.probability(1), 31.0 / 32.0);
        let zero = [0.0; 2];
        let empty = EmissionDistribution::try_new_sources(
            cx,
            &VolumeGrid::new([1, 1, 2], &zero, [0.0; 3], [1.0; 3]),
            &sources,
        )
        .unwrap();
        assert!(empty.cdf.is_empty());
        let cold = source_palette(
            cx,
            vec![thermal_model(cx, &[(1.0, 1.0, Some(1e-6))])],
            &[0, 0],
        );
        assert!(
            EmissionDistribution::try_new_sources(cx, &field, &cold)
                .unwrap()
                .cdf
                .is_empty()
        );
        gate.request();
        assert!(matches!(
            EmissionDistribution::try_new_sources(cx, &field, &sources),
            Err(VolumeTransportError::Cancelled)
        ));
    });
}
