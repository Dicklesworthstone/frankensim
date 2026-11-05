use super::*;
use crate::charts::Ray;
use crate::gradient_index::bundle::{BundleOptions, ProjectionOptions, propagate_bundle};
use crate::gradient_index::{GradientIndexField, GradientIndexOptions, propagate};
use crate::material_state::thermo_optic::transient::P1ConductionTemperatureSnapshot;
use crate::material_state::thermo_optic::transient::mesh::water::pressure_temperature::{
    NodalPressureValues, P1PressureSnapshot, P1PressureTemperatureSnapshot,
};
use crate::material_state::thermo_optic::{TemperatureGradientSample, TemperatureJet};
use crate::material_state::water_refraction::pressure_temperature::PressureTemperatureJet;
use fs_conduction::ConductionMesh;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::gas::{GasSpec, R_USSA_1976};
use fs_material::gas_scattering::{AVOGADRO_PER_MOL, MolecularRayleigh};
use fs_qty::semantic::{Composition, CompositionBasis};
use fs_rep_mesh::TetComplex;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance {tolerance}"
    );
}
fn point(p: [f64; 3]) -> Point3 {
    Point3::new(p[0], p[1], p[2])
}
fn optical(id: u8) -> MolecularRayleigh {
    MolecularRayleigh::try_new(0.0003, 2.5e25, 550.0, 0.03, ContentHash([id; 32])).unwrap()
}
fn model(water: bool, id: u8) -> MolecularRayleighMixture {
    let spec = if water {
        GasSpec::water_vapor_nist()
    } else {
        GasSpec::dry_air_ussa1976()
    };
    MolecularRayleighMixture::try_new(
        &[spec],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        &[optical(id)],
        water.then_some(0),
    )
    .unwrap()
}
fn bind<'a, F: PressureTemperatureGradientField + ?Sized>(
    source: &'a F,
    model: &'a MolecularRayleighMixture,
    id: u8,
) -> PressureTemperatureGasField<'a, F> {
    PressureTemperatureGasField::try_new(
        source,
        model,
        Temperature::new(300.0),
        Pressure::new(100_000.0),
        ContentHash([id; 32]),
    )
    .unwrap()
}
fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 117,
                kernel_id: 0,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}
fn uniform_jet(t: f64, p: f64) -> PressureTemperatureJet {
    PressureTemperatureJet {
        temperature: TemperatureJet {
            sample: TemperatureGradientSample {
                temperature_k: t,
                gradient_k_per_m: Vec3::new(0.0, 0.0, 0.0),
            },
            diagonal_k_per_m2: [0.0; 3],
            off_diagonal_k_per_m2: [0.0; 3],
        },
        pressure_pa: p,
        pressure_gradient_pa_per_m: Vec3::new(0.0, 0.0, 0.0),
        pressure_diagonal_pa_per_m2: [0.0; 3],
        pressure_off_diagonal_pa_per_m2: [0.0; 3],
    }
}
struct Reconstruction<F>(F);
impl<F: Fn(Point3) -> Result<PressureTemperatureJet, GradientIndexError>>
    PressureTemperatureGradientField for Reconstruction<F>
{
    fn pressure_temperature_gradient(
        &self,
        p: Point3,
    ) -> Result<PressureTemperatureGradientSample, GradientIndexError> {
        (self.0)(p).map(|jet| jet.sample())
    }
}
impl<F: Fn(Point3) -> Result<PressureTemperatureJet, GradientIndexError>>
    PressureTemperatureSecondOrderField for Reconstruction<F>
{
    fn pressure_temperature_jet(
        &self,
        p: Point3,
    ) -> Result<PressureTemperatureJet, GradientIndexError> {
        (self.0)(p)
    }
}

#[test]
fn g1_original_joint_p1_source_drives_density_and_full_molecular_hessian() {
    with_cx(|cx| {
        for reverse in [false, true] {
            let mut tet = [0, 1, 2, 3];
            if reverse {
                tet.swap(0, 1);
            }
            let mesh = ConductionMesh::new(
                TetComplex::from_tets(4, vec![tet]),
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
            )
            .unwrap();
            let gt = [5.0, -3.0, 2.0];
            let gp = [20_000.0, -13_000.0, 7_000.0];
            let temperatures = [300.0, 305.0, 297.0, 302.0];
            let pressures = [100_000.0, 120_000.0, 87_000.0, 107_000.0];
            let thermal =
                P1ConductionTemperatureSnapshot::try_new(cx, &mesh, &temperatures).unwrap();
            let pressure = P1PressureSnapshot::try_new(
                cx,
                &mesh,
                NodalPressureValues::Pascals(&pressures),
                ContentHash([117; 32]),
            )
            .unwrap();
            let snapshot =
                P1PressureTemperatureSnapshot::try_new(thermal, ContentHash([118; 32]), pressure)
                    .unwrap();
            let source = snapshot.element(0).unwrap();
            let original = model(false, 119);
            let field = PressureTemperatureGasField::try_new(
                &source,
                &original,
                Temperature::new(300.0),
                Pressure::new(100_000.0),
                snapshot.content_identity(),
            )
            .unwrap();
            let x = [0.18, 0.23, 0.16];
            let native = source.sample_jet(point(x)).unwrap();
            assert_eq!(
                field.sample_state(point(x)).unwrap(),
                *original
                    .evaluate(native.temperature.sample.temperature_k, native.pressure_pa)
                    .unwrap()
                    .gas()
            );
            assert!(core::ptr::eq(field.state_source(), &source));
            assert!(core::ptr::eq(field.model(), &original));
            assert_eq!(field.snapshot_identity(), snapshot.content_identity());
            let density = field.sample_jet(point(x)).unwrap();
            assert_eq!(density.sample, field.sample_gradient(point(x)).unwrap());
            let optics = field.bind_refraction().unwrap();
            let index = optics.sample_optical_jet(point(x), 550.0).unwrap();
            let dg = density.sample.gradient_kg_m4;
            let ng = index.optical().phase_gradient_per_m();
            // Independent scalar original EOS + original molecular source;
            // finite differences belong only to this falsifier, not sampling.
            let scalar = |y: [f64; 3]| {
                let t = 300.0 + gt[0] * y[0] + gt[1] * y[1] + gt[2] * y[2];
                let p = 100_000.0 + gp[0] * y[0] + gp[1] * y[1] + gp[2] * y[2];
                let gas = original.evaluate(t, p).unwrap();
                [
                    gas.gas().density,
                    gas.refractive_index(550.0).unwrap().phase_index,
                ]
            };
            let center = scalar(x);
            let h = 0.005;
            let gradients = [[dg.x, dg.y, dg.z], [ng.x, ng.y, ng.z]];
            let diagonals = [density.hessian_diagonal_kg_m5, index.hessian().diagonal()];
            let mixed = [
                density.hessian_off_diagonal_kg_m5,
                index.hessian().off_diagonal(),
            ];
            for i in 0..3 {
                let (mut a, mut b) = (x, x);
                a[i] += h;
                b[i] -= h;
                let (a, b) = (scalar(a), scalar(b));
                for (k, tolerance) in [2e-9, 5e-11].into_iter().enumerate() {
                    near(gradients[k][i], (a[k] - b[k]) / (2.0 * h), tolerance);
                    near(
                        diagonals[k][i],
                        (a[k] - 2.0 * center[k] + b[k]) / (h * h),
                        tolerance,
                    );
                }
            }
            for (i, (a, b)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
                let (mut pp, mut pm, mut mp, mut mm) = (x, x, x, x);
                pp[a] += h;
                pp[b] += h;
                pm[a] += h;
                pm[b] -= h;
                mp[a] -= h;
                mp[b] += h;
                mm[a] -= h;
                mm[b] -= h;
                let (pp, pm, mp, mm) = (scalar(pp), scalar(pm), scalar(mp), scalar(mm));
                for (k, tolerance) in [2e-9, 5e-11].into_iter().enumerate() {
                    near(
                        mixed[k][i],
                        (pp[k] - pm[k] - mp[k] + mm[k]) / (4.0 * h * h),
                        tolerance,
                    );
                }
            }
            assert!(
                density
                    .hessian_off_diagonal_kg_m5
                    .into_iter()
                    .all(|h| h.abs() > 1e-4)
            );
            assert_eq!(
                optics.sample(point([1.0; 3]), 550.0),
                Err(GradientIndexError::OutsideFieldDomain)
            );
        }
    });
}

#[test]
fn g1_pressure_lens_bends_and_focuses_native_rays_and_angular_bundles() {
    let source = Reconstruction(|p: Point3| {
        if p.x.abs() > 1.0 {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        // Prescribed frozen pressure lens, not an asserted mechanical equilibrium.
        let mut jet = uniform_jet(300.0, 100_000.0 * (1.0 - 0.1 * p.x * p.x));
        jet.pressure_gradient_pa_per_m.x = -20_000.0 * p.x;
        jet.pressure_diagonal_pa_per_m2[0] = -20_000.0;
        Ok(jet)
    });
    let original = model(false, 119);
    let field = bind(&source, &original, 120);
    let optics = field.bind_refraction().unwrap();
    let options = GradientIndexOptions {
        max_step_m: 0.1,
        relative_tolerance: 1e-12,
        position_tolerance_m: 1e-12,
        direction_tolerance: 1e-12,
        optical_path_tolerance_m: 1e-12,
        ..GradientIndexOptions::default()
    };
    // Independent source polarizability and ideal-gas number density. On the
    // symmetry axis the Jacobi equations are X''=-omega² X and Y''=0.
    let n2_minus_one = 0.0003 * (2.0 + 0.0003);
    let volume = (n2_minus_one / (3.0 + n2_minus_one)) / 2.5e25;
    let l = volume * (100_000.0 * AVOGADRO_PER_MOL / (R_USSA_1976 * 300.0));
    let n = ((1.0 + 2.0 * l) / (1.0 - l)).sqrt();
    let dn_dl = 3.0 / (2.0 * n * (1.0 - l).powi(2));
    let omega = (0.2 * dn_dl * l / n).sqrt();
    let length = 20.0;
    let expected_area = length * (omega * length).sin() / omega;
    for wavelength in [380.0, 550.0, 780.0] {
        let bundle = propagate_bundle(
            &optics,
            Ray {
                origin: point([0.0; 3]),
                dir: Vec3::new(0.0, 0.0, 1.0),
            },
            Vec3::new(1.0, 0.0, 0.0),
            wavelength,
            length,
            BundleOptions {
                integration: options,
                position_derivative_tolerance_m: 1e-11,
                direction_derivative_tolerance: 1e-12,
            },
        )
        .unwrap();
        let area = bundle
            .project_at_endpoint(Vec3::new(0.0, 0.0, 1.0), ProjectionOptions::default())
            .unwrap()
            .area_per_solid_angle_m2();
        near(area, expected_area, 2e-8);
        assert!(area < length * length - 0.1);
        near(
            bundle.propagation().optical_path_length_m,
            n * length,
            2e-11,
        );
        near(
            bundle.propagation().group_optical_path_length_m.unwrap(),
            n * length,
            2e-11,
        );
        assert_eq!(bundle.propagation().power_transmittance(), 1.0);
    }
    let bent = propagate(
        &optics,
        Ray {
            origin: point([0.4, 0.0, 0.0]),
            dir: Vec3::new(0.0, 0.0, 1.0),
        },
        Vec3::new(1.0, 0.0, 0.0),
        550.0,
        2.0,
        options,
    )
    .unwrap();
    assert!(
        bent.ray.origin.x < 0.399999,
        "bend toward the higher-pressure symmetry axis"
    );
    near(
        bent.final_sample.phase_index() * bent.ray.dir.z,
        bent.initial_sample.phase_index(),
        2e-12,
    );
}

#[test]
fn g3_correlated_pressure_temperature_jets_cancel_without_extra_source_queries() {
    let calls = core::cell::Cell::new(0usize);
    let source = Reconstruction(|p: Point3| {
        calls.set(calls.get() + 1);
        let t = 320.0 + 5.0 * p.x - 4.0 * p.y + 3.0 * p.z + 2.0 * p.x * p.x - 1.5 * p.y * p.y
            + 0.7 * p.z * p.z
            + 1.1 * p.x * p.y
            - 0.9 * p.x * p.z
            + 0.6 * p.y * p.z;
        let mut jet = uniform_jet(t, 400.0 * t);
        let g = Vec3::new(
            5.0 + 4.0 * p.x + 1.1 * p.y - 0.9 * p.z,
            -4.0 - 3.0 * p.y + 1.1 * p.x + 0.6 * p.z,
            3.0 + 1.4 * p.z - 0.9 * p.x + 0.6 * p.y,
        );
        jet.temperature.sample.gradient_k_per_m = g;
        jet.temperature.diagonal_k_per_m2 = [4.0, -3.0, 1.4];
        jet.temperature.off_diagonal_k_per_m2 = [1.1, -0.9, 0.6];
        jet.pressure_gradient_pa_per_m = g.scale(400.0);
        jet.pressure_diagonal_pa_per_m2 = jet.temperature.diagonal_k_per_m2.map(|h| 400.0 * h);
        jet.pressure_off_diagonal_pa_per_m2 =
            jet.temperature.off_diagonal_k_per_m2.map(|h| 400.0 * h);
        Ok(jet)
    });
    let original = model(false, 119);
    let field = bind(&source, &original, 121);
    assert_eq!(calls.get(), 0);
    for p in [point([0.0; 3]), point([0.2, -0.1, 0.3])] {
        let before = calls.get();
        let jet = field.sample_jet(p).unwrap();
        assert_eq!(calls.get(), before + 1);
        let before = calls.get();
        assert_eq!(jet.sample, field.sample_gradient(p).unwrap());
        assert_eq!(calls.get(), before + 1);
        near(
            jet.sample.density_kg_m3,
            400.0 * GasSpec::dry_air_ussa1976().molar_mass / R_USSA_1976,
            1e-15,
        );
        let g = jet.sample.gradient_kg_m4;
        for derivative in [g.x, g.y, g.z]
            .into_iter()
            .chain(jet.hessian_diagonal_kg_m5)
            .chain(jet.hessian_off_diagonal_kg_m5)
        {
            near(derivative, 0.0, 5e-17);
        }
    }
    assert_ne!(
        field.content_identity(),
        bind(&source, &original, 122).content_identity()
    );
    let other_optical_source = model(false, 122);
    assert_ne!(
        field.content_identity(),
        bind(&source, &other_optical_source, 121).content_identity()
    );
}

#[test]
fn g0_original_gas_domains_labels_source_refusals_and_derivative_range_survive() {
    let p = point([0.0; 3]);
    for water in [false, true] {
        let original = model(water, 119);
        for (t, pressure) in [(50.0, 10_000.0), (300.0, 123_456.0), (2000.0, 10_000_000.0)] {
            let source = Reconstruction(|_| Ok(uniform_jet(t, pressure)));
            let field = bind(&source, &original, 123);
            assert_eq!(
                field.sample_state(p).unwrap(),
                *original.evaluate(t, pressure).unwrap().gas()
            );
            assert_eq!(field.temperature_k(p).unwrap(), t);
            assert_eq!(
                field.sample_jet(p).unwrap().hessian_diagonal_kg_m5,
                [0.0; 3]
            );
        }
        for (t, pressure) in [
            (49.9, 1e5),
            (2000.1, 1e5),
            (f64::NAN, 1e5),
            (f64::INFINITY, 1e5),
            (300.0, 0.0),
            (300.0, -1.0),
            (300.0, 10_000_001.0),
            (300.0, f64::NAN),
            (300.0, f64::INFINITY),
        ] {
            let source = Reconstruction(|_| Ok(uniform_jet(t, pressure)));
            let field = bind(&source, &original, 123);
            assert!(matches!(
                field.sample_state(p),
                Err(PressureTemperatureGasError::Thermodynamics(_))
            ));
            assert_eq!(
                field.sample_density(p),
                Err(GradientIndexError::OutsideFieldDomain)
            );
            assert_eq!(
                field.sample_density_second_order(p),
                Err(GradientIndexError::OutsideFieldDomain)
            );
        }
    }
    let original = model(false, 119);
    let dense = Reconstruction(|_| Ok(uniform_jet(50.0, 10_000_000.0)));
    let dense = bind(&dense, &original, 123);
    assert!(dense.sample_state(p).is_ok());
    assert_eq!(
        dense.bind_refraction().unwrap().sample(p, 550.0),
        Err(GradientIndexError::OutsideFieldDomain)
    );
    for i in 0..4 {
        let mut jet = uniform_jet(300.0, 100_000.0);
        match i {
            0 => jet.temperature.sample.gradient_k_per_m.x = f64::NAN,
            1 => jet.pressure_gradient_pa_per_m.z = f64::INFINITY,
            2 => jet.temperature.off_diagonal_k_per_m2[1] = f64::NAN,
            _ => jet.pressure_diagonal_pa_per_m2[2] = f64::INFINITY,
        }
        let source = Reconstruction(|_| Ok(jet));
        assert_eq!(
            bind(&source, &original, 123).sample_density_second_order(p),
            Err(GradientIndexError::InvalidSample)
        );
    }
    let overflow = Reconstruction(|_| {
        let mut jet = uniform_jet(50.0, 10_000_000.0);
        jet.temperature.sample.gradient_k_per_m.x = f64::MAX;
        Ok(jet)
    });
    assert_eq!(
        bind(&overflow, &original, 123).sample_density(p),
        Err(GradientIndexError::NumericalRange)
    );
    let cancelled = Reconstruction(|_| Err(GradientIndexError::Cancelled));
    let field = bind(&cancelled, &original, 123);
    assert_eq!(field.sample_density(p), Err(GradientIndexError::Cancelled));
    assert_eq!(
        field.sample_density_second_order(p),
        Err(GradientIndexError::Cancelled)
    );
    assert_eq!(
        field.sample_density(point([f64::NAN, 0.0, 0.0])),
        Err(GradientIndexError::OutsideFieldDomain)
    );
    assert!(matches!(
        PressureTemperatureGasField::try_new(
            &cancelled,
            &original,
            Temperature::new(300.0),
            Pressure::new(100_000.0),
            ContentHash([0; 32])
        ),
        Err(PressureTemperatureGasError::MissingSnapshotIdentity)
    ));
    assert!(
        PressureTemperatureGasField::try_new(
            &cancelled,
            &original,
            Temperature::new(300.0),
            Pressure::new(0.0),
            ContentHash([123; 32])
        )
        .is_err()
    );
    let multiple = MolecularRayleighMixture::try_new(
        &[GasSpec::dry_air_ussa1976(), GasSpec::water_vapor_nist()],
        Composition::new(CompositionBasis::MoleFraction, vec![0.99, 0.01]).unwrap(),
        &[optical(119); 2],
        Some(1),
    )
    .unwrap();
    assert!(matches!(
        PressureTemperatureGasField::try_new(
            &cancelled,
            &multiple,
            Temperature::new(300.0),
            Pressure::new(100_000.0),
            ContentHash([123; 32])
        ),
        Err(PressureTemperatureGasError::UnsupportedComposition)
    ));
}
