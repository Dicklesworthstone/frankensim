//! Original-law attenuation and accepted-temperature layered-source references.

use super::super::tests::{complex_fixture, density_key, point};
use super::*;
use crate::camera::{Aperture, CameraProjection};
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGridBudget, CellPaletteBudget, MassAbsorptionSpectrum};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use std::sync::Arc;

const LENGTH: f64 = 1e-5;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 17,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    });
}

fn geometry<'a>(cx: &Cx<'_>, density: &'a [f64], steps: u32) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new(
            [1, 1, density.len()],
            density,
            [-1.0, -1.0, 0.0],
            [2.0, 2.0, LENGTH],
        ),
        // An unrelated visible-only optical law must not become the material
        // field's absorption, thermal source or infrared-support authority.
        MassAbsorptionSpectrum::constant(9.0).unwrap(),
        None,
        AbsorbingGridBudget {
            max_grid_cells: density.len(),
            max_steps_per_segment: steps,
        },
    )
    .unwrap()
}

fn palette<'a>(
    cx: &Cx<'_>,
    values: Vec<MaterialAbsorptionBinding<'a>>,
    indices: &[u32],
) -> CellPalette<MaterialAbsorptionBinding<'a>> {
    let budget = CellPaletteBudget {
        max_models: values.len(),
        max_cells: indices.len(),
    };
    CellPalette::try_new(cx, values, Arc::from(indices), budget).unwrap()
}

fn budget() -> MaterialAbsorptionPreparationBudget {
    MaterialAbsorptionPreparationBudget {
        max_models: 8,
        max_spectral_evaluations: 128,
    }
}

fn settings() -> VolumeCameraSettings {
    VolumeCameraSettings {
        resolution: [1, 1],
        samples_per_pixel: 1,
        seed: 71,
        max_pixels: 1,
        max_samples: 1,
        max_cell_visits: 16,
    }
}

fn origin() -> Point3 {
    Point3::new(0.0, 0.0, -LENGTH)
}
fn direction() -> Vec3 {
    Vec3::new(0.0, 0.0, 1.0)
}
fn camera() -> PhysicalCamera {
    PhysicalCamera::try_look_at(
        origin(),
        Point3::new(0.0, 0.0, LENGTH),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_half_tangent(0.1).unwrap(),
        2.0 * LENGTH,
        Aperture::try_circular(0.0).unwrap(),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8e-12 * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

// Independent SI Planck equation at moderate test wavelengths/temperatures.
fn planck(temperature: f64, wavelength_nm: f64) -> f64 {
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0;
    let k = 1.380_649e-23;
    let wavelength = wavelength_nm * 1e-9;
    2.0 * h * c * c / wavelength.powi(5) / (h * c / (wavelength * k * temperature)).exp_m1() * 1e-9
}

#[test]
fn g1_original_index_law_drives_native_packets_and_physical_camera() {
    with_cx(|cx, _| {
        let (_, accepted, optical) = complex_fixture([0.01, 0.03]);
        let material =
            MaterialAbsorptionBinding::try_spectral(&accepted, &optical, &density_key()).unwrap();
        let density = [1000.0];
        let geometry = geometry(cx, &density, 16);
        let materials = palette(cx, vec![material], &[0]);
        let original_map = Arc::clone(materials.cell_indices());
        let wavelengths = [8000.0, 10000.0, 12000.0];
        let field = MaterialAbsorptionField::try_new(
            cx,
            &geometry,
            materials,
            wavelengths,
            MaterialAbsorptionEmission::None,
            budget(),
        )
        .unwrap();
        assert!(Arc::ptr_eq(field.materials().cell_indices(), &original_map));
        assert_eq!(field.source_wavelength_bounds_nm(), [8000.0, 12000.0]);
        let native = field
            .transport(cx, origin(), direction(), f64::INFINITY)
            .unwrap();
        let image = field.render_camera(cx, &camera(), settings()).unwrap();
        assert_eq!(image.wavelengths_nm(), &wavelengths);
        assert_eq!(image.field_identity(), field.content_identity());
        assert_eq!(image.pixels()[0].emitted_radiance_w_m2_sr_nm(), &[0.0; 3]);
        let illuminated = field
            .render_camera_with_background(cx, &camera(), settings(), |_, _, _, nodes| {
                assert_eq!(nodes, &wavelengths);
                Ok(VolumeCameraBackground {
                    radiance: [2.0, 3.0, 4.0],
                    distance_m: f64::INFINITY,
                })
            })
            .unwrap();
        for (lane, (wavelength, k)) in wavelengths.into_iter().zip([0.01, 0.02, 0.03]).enumerate() {
            let tau = 4.0 * core::f64::consts::PI * k / (wavelength * 1e-9) * LENGTH;
            close(native.transmittance[lane], (-tau).exp());
            assert_eq!(
                image.pixels()[0].transmittance()[lane],
                native.transmittance[lane]
            );
            close(image.pixels()[0].effective_optical_depth()[lane], tau);
            close(
                illuminated.radiance_w_m2_sr_nm()[0][lane],
                [2.0, 3.0, 4.0][lane] * (-tau).exp(),
            );
        }
        let converted_knot_midpoint = 0.5
            * (material.mass_absorption_at(8000.0).unwrap()
                + material.mass_absorption_at(12000.0).unwrap());
        assert!(
            (material.mass_absorption_at(10000.0).unwrap() / converted_knot_midpoint - 1.0).abs()
                > 0.06
        );
    });
}

#[test]
fn g1_each_materials_accepted_lte_source_obeys_layer_order_and_finite_endpoint() {
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
        // A single exact node is admitted without inventing a second knot.
        let field = MaterialAbsorptionField::try_new(
            cx,
            &geometry,
            materials,
            [10000.0],
            MaterialAbsorptionEmission::AcceptedLte,
            budget(),
        )
        .unwrap();
        let a = 4.0 * core::f64::consts::PI * 0.02;
        let b = 4.0 * core::f64::consts::PI * 0.04;
        let expected = planck(300.0, 10000.0) * -(-a).exp_m1()
            + (-a).exp() * planck(400.0, 10000.0) * -(-b).exp_m1();
        let native = field
            .transport(cx, origin(), direction(), f64::INFINITY)
            .unwrap();
        close(native.emission_w_m2_sr_nm[0], expected);
        close(native.transmittance[0], (-a - b).exp());
        assert!((expected / (planck(350.0, 10000.0) * -(-a - b).exp_m1()) - 1.0).abs() > 0.05);
        let image = field.render_camera(cx, &camera(), settings()).unwrap();
        assert_eq!(
            image.pixels()[0].emitted_radiance_w_m2_sr_nm()[0],
            native.emission_w_m2_sr_nm[0]
        );
        let finite = field
            .render_camera_with_background(cx, &camera(), settings(), |_, _, _, _| {
                Ok(VolumeCameraBackground {
                    radiance: [7.0],
                    distance_m: 2.0 * LENGTH,
                })
            })
            .unwrap();
        close(
            finite.radiance_w_m2_sr_nm()[0][0],
            planck(300.0, 10000.0) * -(-a).exp_m1() + 7.0 * (-a).exp(),
        );
        let reversed = field
            .transport(
                cx,
                Point3::new(0.0, 0.0, 3.0 * LENGTH),
                Vec3::new(0.0, 0.0, -1.0),
                f64::INFINITY,
            )
            .unwrap();
        close(
            reversed.emission_w_m2_sr_nm[0],
            planck(400.0, 10000.0) * -(-b).exp_m1()
                + (-b).exp() * planck(300.0, 10000.0) * -(-a).exp_m1(),
        );
    });
}

#[test]
fn g0_cold_absorption_needs_no_temperature_and_geometry_cannot_supply_one() {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyValue, QueryPoint,
    };
    use fs_material::state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        SpectralComplexIndexStatePoint, resolve_material_state_point,
    };
    with_cx(|cx, _| {
        let (source, accepted, optical) = complex_fixture([0.01, 0.03]);
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        for (_, original) in source.claims().claims_ordered() {
            let mut claim = original.clone();
            claim.value = PropertyValue::Scalar {
                value: accepted.property_by_key(&claim.key).unwrap().value_si(),
                dims: claim.key.quantity().dims(),
            };
            claim.validity = ValidityDomain::unconstrained();
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
            requirements.push(
                ScalarPropertyRequirement::try_with_key(&claim.key, ScalarAdmissibility::Finite)
                    .unwrap(),
            );
            claims.insert_claim(claim).unwrap();
        }
        let source = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic temperature-free absorber".into(),
                phase: "solid".into(),
                process: "constant analytic test".into(),
                revision: 0,
            },
            claims,
            Vec::new(),
        )
        .unwrap();
        let accepted = resolve_material_state_point(
            &source,
            &QueryPoint::new(),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let optical =
            SpectralComplexIndexStatePoint::try_from_resolved(accepted.clone(), optical.bindings())
                .unwrap();
        let material =
            MaterialAbsorptionBinding::try_spectral(&accepted, &optical, &density_key()).unwrap();
        assert_eq!(material.accepted_temperature_k(), None);
        let density = [1000.0];
        let temperature = [6000.0];
        let geometry = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, LENGTH]),
            MassAbsorptionSpectrum::constant(1.0).unwrap(),
            Some(&temperature),
            AbsorbingGridBudget {
                max_grid_cells: 1,
                max_steps_per_segment: 16,
            },
        )
        .unwrap();
        let materials = palette(cx, vec![material], &[0]);
        let cold = MaterialAbsorptionField::try_new(
            cx,
            &geometry,
            materials.clone(),
            [10000.0],
            MaterialAbsorptionEmission::None,
            budget(),
        )
        .unwrap();
        assert_eq!(
            cold.transport(cx, origin(), direction(), f64::INFINITY)
                .unwrap()
                .emission_w_m2_sr_nm,
            [0.0]
        );
        assert!(matches!(
            MaterialAbsorptionField::try_new(
                cx,
                &geometry,
                materials,
                [10000.0],
                MaterialAbsorptionEmission::AcceptedLte,
                budget()
            ),
            Err(MaterialAbsorptionFieldError::Temperature { model: 0 })
        ));
    });
}

#[test]
fn g0_source_map_preparation_and_traversal_limits_precede_publication() {
    with_cx(|cx, gate| {
        let (_, accepted, optical) = complex_fixture([0.01, 0.03]);
        let material =
            MaterialAbsorptionBinding::try_spectral(&accepted, &optical, &density_key()).unwrap();
        let density = [1000.0; 2];
        let geometry = geometry(cx, &density, 1);
        let materials = palette(cx, vec![material], &[0, 0]);
        let insufficient = MaterialAbsorptionPreparationBudget {
            max_models: 1,
            max_spectral_evaluations: 0,
        };
        assert!(matches!(
            MaterialAbsorptionField::try_new(
                cx,
                &geometry,
                materials.clone(),
                [10000.0],
                MaterialAbsorptionEmission::None,
                insufficient
            ),
            Err(MaterialAbsorptionFieldError::Transport(
                VolumeTransportError::Budget
            ))
        ));
        assert!(matches!(
            MaterialAbsorptionField::try_new(
                cx,
                &geometry,
                palette(cx, vec![material], &[0]),
                [10000.0],
                MaterialAbsorptionEmission::None,
                budget()
            ),
            Err(MaterialAbsorptionFieldError::Transport(
                VolumeTransportError::Palette
            ))
        ));
        assert!(matches!(
            MaterialAbsorptionField::try_new(
                cx,
                &geometry,
                materials.clone(),
                [7000.0],
                MaterialAbsorptionEmission::None,
                budget()
            ),
            Err(MaterialAbsorptionFieldError::Absorption(
                MaterialAbsorptionError::WavelengthOutsideSource { .. }
            ))
        ));
        assert!(matches!(
            MaterialAbsorptionField::try_new(
                cx,
                &geometry,
                materials.clone(),
                [10000.0, 10000.0],
                MaterialAbsorptionEmission::None,
                budget()
            ),
            Err(MaterialAbsorptionFieldError::Absorption(
                MaterialAbsorptionError::InvalidWavelengthSamples
            ))
        ));
        let field = MaterialAbsorptionField::try_new(
            cx,
            &geometry,
            materials.clone(),
            [10000.0],
            MaterialAbsorptionEmission::None,
            budget(),
        )
        .unwrap();
        assert_eq!(
            field.transmittance(cx, origin(), direction(), f64::INFINITY),
            Err(VolumeTransportError::Budget)
        );
        let mut zero_cells = settings();
        zero_cells.max_cell_visits = 0;
        assert!(matches!(
            field.render_camera(cx, &camera(), zero_cells),
            Err(VolumeCameraError::Transport(VolumeTransportError::Budget))
        ));
        gate.request();
        assert_eq!(
            field.transport(cx, origin(), direction(), 0.0),
            Err(VolumeTransportError::Cancelled)
        );
        assert!(matches!(
            MaterialAbsorptionField::try_new(
                cx,
                &geometry,
                materials,
                [10000.0],
                MaterialAbsorptionEmission::None,
                budget()
            ),
            Err(MaterialAbsorptionFieldError::Transport(
                VolumeTransportError::Cancelled
            ))
        ));
    });
}

#[test]
fn g1_unreferenced_absorber_cannot_create_haze_or_thermal_emission() {
    with_cx(|cx, _| {
        let (_, clear_state, clear_optical) = complex_fixture([0.0, 0.0]);
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
        let geometry = geometry(cx, &density, 1);
        let field = MaterialAbsorptionField::try_new(
            cx,
            &geometry,
            palette(cx, vec![clear, absorber], &[0, 0]),
            [10000.0],
            MaterialAbsorptionEmission::AcceptedLte,
            budget(),
        )
        .unwrap();
        assert!(field.is_vacuum());
        let mut zero_cells = settings();
        zero_cells.max_cell_visits = 0;
        let observed = field
            .render_camera_with_background(cx, &camera(), zero_cells, |_, _, _, _| {
                Ok(VolumeCameraBackground {
                    radiance: [3.0],
                    distance_m: f64::INFINITY,
                })
            })
            .unwrap();
        assert_eq!(observed.radiance_w_m2_sr_nm(), &[[3.0]]);
        let sample = field
            .transport(cx, origin(), direction(), f64::INFINITY)
            .unwrap();
        assert_eq!(sample.transmittance, [1.0]);
        assert_eq!(sample.emission_w_m2_sr_nm, [0.0]);
    });
}
