use super::*;
use crate::gradient_index::GradientIndexField;
use crate::material_state::gas_scattering::GasRefractiveFieldError;
use crate::volumes::transport::scattering::PhaseFunction;
use fs_geom::Vec3;
use fs_material::gas::GasSpec;
use fs_material::gas_scattering::composition::{
    MolecularCauchyDispersion, MolecularRefractionModel, PartialDensityBasis,
};
use fs_material::gas_scattering::{AVOGADRO_PER_MOL, MolecularRayleigh, MolecularRayleighMixture};
use fs_qty::semantic::{Composition, CompositionBasis};
use std::sync::atomic::{AtomicU8, Ordering};

fn model() -> MolecularRefractionModel {
    let mut specs = [GasSpec::dry_air_ussa1976(); 2];
    specs[1].molar_mass *= 2.0;
    let optics = [(0.00015, 0.01), (0.0006, 0.2)].map(|(r, delta)| {
        MolecularRayleigh::try_new(r, 2.5e25, 550.0, delta, ContentHash([11; 32])).unwrap()
    });
    let original = MolecularRayleighMixture::try_new(
        &specs,
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
        &optics,
        None,
    )
    .unwrap();
    let reference = MolecularRefractionModel::from_mixture(
        &original,
        original.evaluate(300.0, 100_000.0).unwrap().gas(),
    )
    .unwrap();
    let shapes: Vec<_> = reference
        .components()
        .iter()
        .enumerate()
        .map(|(i, row)| {
            MolecularCauchyDispersion::try_new(
                row.optical(),
                [1.0e-6, 2.0e-5][i],
                [1.0e-7, 3.0e-7][i],
                [400.0, 700.0],
                ContentHash([i as u8 + 12; 32]),
            )
            .unwrap()
        })
        .collect();
    reference.with_dispersion(&shapes).unwrap()
}

// The source fixture has constant total mass density and changing composition.
// Metadata mutation below is only a refusal fixture, never an admitted snapshot.
struct CompositionField {
    identity: ContentHash,
    basis: PartialDensityBasis,
    scale: [f64; 2],
    total_mass_density: f64,
    written: usize,
    metadata_mode: AtomicU8,
}

impl CompositionField {
    fn new(model: &MolecularRefractionModel, basis: PartialDensityBasis) -> Self {
        Self {
            identity: model.content_identity(),
            basis,
            scale: core::array::from_fn(|i| match basis {
                PartialDensityBasis::MassKgPerM3 => 1.0,
                PartialDensityBasis::MoleculesPerM3 => {
                    AVOGADRO_PER_MOL / model.components()[i].molar_mass_kg_mol()
                }
            }),
            total_mass_density: 1.0,
            written: 2,
            metadata_mode: AtomicU8::new(0),
        }
    }
}

impl GasCompositionGradientField for CompositionField {
    fn model_identity(&self) -> ContentHash {
        if self.metadata_mode.load(Ordering::Relaxed) == 1 {
            ContentHash([99; 32])
        } else {
            self.identity
        }
    }
    fn component_count(&self) -> usize {
        if self.metadata_mode.load(Ordering::Relaxed) == 2 {
            1
        } else {
            2
        }
    }
    fn density_basis(&self) -> PartialDensityBasis {
        if self.metadata_mode.load(Ordering::Relaxed) == 3 {
            PartialDensityBasis::MoleculesPerM3
        } else {
            self.basis
        }
    }
    fn sample_partial_densities(
        &self,
        p: Point3,
        output: &mut [GasPartialDensitySample],
    ) -> Result<(), GradientIndexError> {
        if !(0.0..=1.0).contains(&p.x) {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        for (i, row) in output.iter_mut().take(self.written).enumerate() {
            let factor = self.total_mass_density * self.scale[i];
            *row = GasPartialDensitySample {
                density: factor * (if i == 0 { 1.0 - p.x } else { p.x }),
                gradient_per_m: Vec3::new(if i == 0 { -factor } else { factor }, 0.0, 0.0),
            };
        }
        Ok(())
    }
}

fn bind<'a>(
    model: &'a MolecularRefractionModel,
    provider: &'a CompositionField,
) -> MaterialGasCompositionRefractiveField<'a, CompositionField> {
    MaterialGasCompositionRefractiveField::try_new(model, provider, ContentHash([17; 32])).unwrap()
}

fn close(got: f64, expected: f64) {
    assert!(
        (got - expected).abs() <= 3.0e-14 * expected.abs() + 1.0e-30,
        "{got} != {expected}"
    );
}

#[test]
fn g1_dispersive_rates_and_local_phases_use_the_same_original_composition() {
    let model = model();
    let provider = CompositionField::new(&model, PartialDensityBasis::MassKgPerM3);
    let refractive = bind(&model, &provider);
    let field = refractive.bind_scattering_field(&[1.0, 1.0]).unwrap();
    let native: &dyn CurvedRayleighField = &field;
    for wavelength in [410.0, 550.0, 690.0] {
        let spectral = model.at_wavelength(wavelength, provider.basis).unwrap();
        let majorant = native.scattering_majorant_per_m(wavelength).unwrap();
        assert_eq!(
            majorant,
            spectral.scattering(&[1.0, 1.0]).unwrap().scattering_per_m
        );
        for x in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let point = Point3::new(x, 0.0, 0.0);
            let expected = spectral.scattering(&[1.0 - x, x]).unwrap();
            assert_eq!(
                field.sample_scattering_optics(point, wavelength).unwrap(),
                expected
            );
            let actual = native.sample_scattering(point, wavelength).unwrap();
            assert_eq!(actual.scattering_per_m, expected.scattering_per_m);
            assert_eq!(actual.depolarization_ratio, expected.depolarization_ratio);
            assert!(actual.scattering_per_m > 0.0 && actual.scattering_per_m <= majorant);
            assert_eq!(expected.mass_density_kg_m3, 1.0);
        }
        close(
            native
                .sample_scattering(Point3::new(0.0, 0.0, 0.0), wavelength)
                .unwrap()
                .depolarization_ratio,
            0.01,
        );
        close(
            native
                .sample_scattering(Point3::new(1.0, 0.0, 0.0), wavelength)
                .unwrap()
                .depolarization_ratio,
            0.2,
        );
    }
    let point = Point3::new(0.25, 0.0, 0.0);
    let blue = native.sample_scattering(point, 410.0).unwrap();
    let red = native.sample_scattering(point, 690.0).unwrap();
    assert!(blue.scattering_per_m / red.scattering_per_m > (690.0_f64 / 410.0).powi(4));
    assert!((blue.depolarization_ratio - red.depolarization_ratio).abs() > 1.0e-5);
    let index_a = refractive
        .sample(Point3::new(0.0, 0.0, 0.0), 550.0)
        .unwrap();
    let index_b = refractive
        .sample(Point3::new(1.0, 0.0, 0.0), 550.0)
        .unwrap();
    assert!(index_b.phase_index() > index_a.phase_index());
}

#[test]
fn g1_mass_and_number_bases_share_physics_but_retain_source_and_bound_identities() {
    let model = model();
    let mass_provider = CompositionField::new(&model, PartialDensityBasis::MassKgPerM3);
    let number_provider = CompositionField::new(&model, PartialDensityBasis::MoleculesPerM3);
    let mass = bind(&model, &mass_provider);
    let number = bind(&model, &number_provider);
    let mass_field = mass.bind_scattering_field(&mass_provider.scale).unwrap();
    let number_field = number
        .bind_scattering_field(&number_provider.scale)
        .unwrap();
    let broader = mass.bind_scattering_field(&[1.1, 1.1]).unwrap();
    assert!(core::ptr::eq(mass_field.refractive_field(), &mass));
    assert_eq!(mass_field.source_field_identity(), mass.content_identity());
    assert_eq!(broader.source_field_identity(), mass.content_identity());
    assert_ne!(mass_field.content_identity(), broader.content_identity());
    assert_ne!(
        mass_field.source_field_identity(),
        number_field.source_field_identity()
    );
    assert_eq!(mass_field.component_upper_bounds(), &[1.0, 1.0]);
    for wavelength in [410.0, 550.0, 690.0] {
        let p = Point3::new(0.37, 0.0, 0.0);
        let a = mass_field.sample_scattering_optics(p, wavelength).unwrap();
        let b = number_field
            .sample_scattering_optics(p, wavelength)
            .unwrap();
        close(a.scattering_per_m, b.scattering_per_m);
        close(a.depolarization_ratio, b.depolarization_ratio);
        close(a.mass_density_kg_m3, b.mass_density_kg_m3);
        assert_eq!(a, broader.sample_scattering_optics(p, wavelength).unwrap());
    }
}

#[test]
fn g1_effective_phase_preserves_the_complete_molecular_mueller_mixture() {
    let model = model();
    let provider = CompositionField::new(&model, PartialDensityBasis::MassKgPerM3);
    let refractive = bind(&model, &provider);
    let field = refractive.bind_scattering_field(&[1.0, 1.0]).unwrap();
    let original_phases = model
        .components()
        .iter()
        .map(|component| {
            PhaseFunction::rayleigh_depolarized(component.optical().depolarization_ratio()).unwrap()
        })
        .collect::<Vec<_>>();
    for wavelength in [410.0, 550.0, 690.0] {
        let spectral = model.at_wavelength(wavelength, provider.basis).unwrap();
        for x in [0.2, 0.65] {
            let local = field
                .sample_scattering_optics(Point3::new(x, 0.0, 0.0), wavelength)
                .unwrap();
            let beta = [
                spectral
                    .scattering(&[1.0 - x, 0.0])
                    .unwrap()
                    .scattering_per_m,
                spectral.scattering(&[0.0, x]).unwrap().scattering_per_m,
            ];
            let effective =
                PhaseFunction::rayleigh_depolarized(local.depolarization_ratio).unwrap();
            for cosine in [-1.0, -0.37, 0.0, 0.43, 1.0] {
                let actual = effective.rayleigh_mueller(cosine).unwrap().0;
                let a = original_phases[0].rayleigh_mueller(cosine).unwrap().0;
                let b = original_phases[1].rayleigh_mueller(cosine).unwrap().0;
                for row in 0..4 {
                    for column in 0..4 {
                        close(
                            actual[row][column],
                            (beta[0] * a[row][column] + beta[1] * b[row][column])
                                / local.scattering_per_m,
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn g1_native_curved_scene_consumes_prepared_optics_and_rejects_a_foreign_snapshot() {
    use crate::charts::TriMesh;
    use crate::gradient_index::event::GradientEventOptions;
    use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
    use crate::polarization::Stokes;
    use crate::polarization::dielectric::scene::curved::participating::{
        CurvedRayleighBinding, CurvedRayleighScene, CurvedRayleighTraceSettings,
    };
    use crate::polarization::dielectric::scene::curved::{
        CurvedPolarizedScene, CurvedSceneError, CurvedTraceSettings, SceneMedium,
    };
    use crate::polarization::dielectric::scene::{
        Interaction, PathTermination, SceneError, SurfaceBinding,
    };
    use crate::polarization::dielectric::transport::PolarizationFrame;
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 73,
                kernel_id: 19,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let model = model();
        let provider = CompositionField::new(&model, PartialDensityBasis::MassKgPerM3);
        let refractive = bind(&model, &provider);
        let scattering = refractive.bind_scattering_field(&[1.0, 1.0]).unwrap();
        let prepared = refractive.at_wavelength(550.0).unwrap();
        let geometry = InstanceScene::try_new(vec![
            GeometryInstance::try_new(
                71,
                ContentHash([71; 32]),
                SharedGeometry::mesh(TriMesh::new(
                    vec![[-10.0, -10.0, 1.0], [10.0, -10.0, 1.0], [0.0, 10.0, 1.0]],
                    vec![[0, 1, 2]],
                )),
                RigidTransform::identity(),
            )
            .unwrap(),
        ])
        .unwrap();
        let media = [SceneMedium::GradientIndex {
            field: &prepared,
            identity: prepared.content_identity(),
        }];
        let bindings = vec![SurfaceBinding {
            object_id: 71,
            interaction: Interaction::Opaque,
        }];
        let scene =
            CurvedPolarizedScene::try_new(&cx, &geometry, &media, bindings.clone(), 100_000)
                .unwrap();
        let fields = [Some(CurvedRayleighBinding::Field(&scattering))];
        let participating = CurvedRayleighScene::try_new_with_fields(&cx, &scene, &fields).unwrap();
        assert!(core::ptr::eq(participating.surface_scene(), &scene));
        let beam = scene
            .launch_beam(
                Point3::new(0.25, 0.0, 0.0),
                PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 1.0, 0.0))
                    .unwrap(),
                550.0,
                0,
                Stokes::unpolarized(1.0),
            )
            .unwrap();
        let rate = scattering.scattering_majorant_per_m(550.0).unwrap();
        let cdf = -(-rate * 0.25).exp_m1();
        let path = participating
            .trace_with_tracking_samples(
                &cx,
                beam,
                0,
                CurvedRayleighTraceSettings {
                    surfaces: CurvedTraceSettings {
                        max_surface_events: 4,
                        max_distance_m: 2.0,
                        max_total_steps: 16_384,
                        max_query_work: 1_000_000,
                        event: GradientEventOptions::default(),
                    },
                    max_volume_events: 1,
                },
                |_| panic!("candidate precedes the opaque surface"),
                |_| [cdf, 0.4, 0.7, 0.0],
            )
            .unwrap();
        assert_eq!(path.volume_events(), 1);
        assert_eq!(path.free_flight_draws(), 1);
        assert_eq!(path.path().termination(), PathTermination::VolumeLimit);
        assert!((path.path().traveled_m() - 0.25).abs() < 1.0e-7);
        assert_eq!(
            path.path().medium_snapshot_identity(),
            Some(refractive.content_identity())
        );
        let wrong_media = [SceneMedium::GradientIndex {
            field: &refractive,
            identity: scattering.content_identity(),
        }];
        let wrong_scene =
            CurvedPolarizedScene::try_new(&cx, &geometry, &wrong_media, bindings, 100_000).unwrap();
        assert!(matches!(
            CurvedRayleighScene::try_new_with_fields(&cx, &wrong_scene, &fields),
            Err(CurvedSceneError::Scene(SceneError::InvalidBindings))
        ));
    });
}

#[test]
fn g0_invalid_or_exceeded_component_envelopes_and_unwritten_rows_refuse() {
    let model = model();
    let provider = CompositionField::new(&model, PartialDensityBasis::MassKgPerM3);
    let refractive = bind(&model, &provider);
    assert!(matches!(
        refractive.bind_scattering_field(&[1.0]),
        Err(GasCompositionScatteringError::Components)
    ));
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            refractive.bind_scattering_field(&[1.0, bad]),
            Err(GasCompositionScatteringError::InvalidBound { component: 1 })
        ));
    }
    assert!(matches!(
        refractive.bind_scattering_field(&[1000.0, 1000.0]),
        Err(GasCompositionScatteringError::Gas(
            GasScatteringError::RefractivityDomain
        ))
    ));
    let too_small = refractive.bind_scattering_field(&[0.5, 1.0]).unwrap();
    assert!(matches!(
        too_small.sample_scattering_optics(Point3::new(0.0, 0.0, 0.0), 550.0),
        Err(GasCompositionScatteringError::BoundExceeded { component: 0 })
    ));
    assert_eq!(
        too_small.sample_scattering(Point3::new(0.0, 0.0, 0.0), 550.0),
        Err(VolumeTransportError::Density { index: 0 })
    );
    let mut partial = CompositionField::new(&model, PartialDensityBasis::MassKgPerM3);
    partial.written = 1;
    let incomplete = bind(&model, &partial);
    let field = incomplete.bind_scattering_field(&[1.0, 1.0]).unwrap();
    assert!(matches!(
        field.sample_scattering_optics(Point3::new(0.5, 0.0, 0.0), 550.0),
        Err(GasCompositionScatteringError::InvalidSample { component: 1 })
    ));
}

#[test]
fn g0_original_snapshot_source_order_and_units_cannot_change_silently() {
    let model = model();
    let provider = CompositionField::new(&model, PartialDensityBasis::MassKgPerM3);
    assert!(matches!(
        MaterialGasCompositionRefractiveField::try_new(&model, &provider, ContentHash([0; 32])),
        Err(GasRefractiveFieldError::MissingSnapshotIdentity)
    ));
    let refractive = bind(&model, &provider);
    let field = refractive.bind_scattering_field(&[1.0, 1.0]).unwrap();
    for mode in [1, 2, 3] {
        provider.metadata_mode.store(mode, Ordering::Relaxed);
        let error = if mode == 2 {
            GasCompositionScatteringError::Components
        } else {
            GasCompositionScatteringError::SourceMismatch
        };
        assert_eq!(
            field.sample_scattering_optics(Point3::new(0.5, 0.0, 0.0), 550.0),
            Err(error.clone())
        );
        assert_eq!(field.scattering_majorant_optics(550.0), Err(error));
        assert!(refractive.bind_scattering_field(&[1.0, 1.0]).is_err());
    }
}

#[test]
fn g0_vacuum_is_exact_and_original_domain_refusals_are_not_vacuum() {
    let model = model();
    let mut provider = CompositionField::new(&model, PartialDensityBasis::MassKgPerM3);
    provider.total_mass_density = 0.0;
    let refractive = bind(&model, &provider);
    let field = refractive.bind_scattering_field(&[0.0, 0.0]).unwrap();
    assert_eq!(field.scattering_majorant_per_m(550.0).unwrap(), 0.0);
    assert_eq!(
        field
            .sample_scattering(Point3::new(0.5, 0.0, 0.0), 550.0)
            .unwrap(),
        CurvedRayleighSample {
            scattering_per_m: 0.0,
            depolarization_ratio: 0.0
        }
    );
    assert!(matches!(
        field.sample_scattering_optics(Point3::new(2.0, 0.0, 0.0), 550.0),
        Err(GasCompositionScatteringError::Field(
            GradientIndexError::OutsideFieldDomain
        ))
    ));
    for wavelength in [399.0, 701.0] {
        assert!(field.scattering_majorant_optics(wavelength).is_err());
        assert!(
            field
                .sample_scattering_optics(Point3::new(0.5, 0.0, 0.0), wavelength)
                .is_err()
        );
    }
}
