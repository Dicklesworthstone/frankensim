use super::*;
use crate::charts::Ray;
use crate::gradient_index::bundle::{BundleOptions, ProjectionOptions, propagate_bundle};
use crate::gradient_index::{GradientIndexField, GradientIndexOptions, propagate};
use crate::material_state::thermo_optic::{
    TemperatureJet, conduction::PlanarConductionTemperatureField,
};
use fs_conduction::{
    ConductivityTable,
    layered::{GeneratedLayeredWallSolution, LayeredWall, PlanarThermalLayer},
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::Vec3;
use fs_material::gas::{GasSpec, R_USSA_1976};
use fs_material::gas_scattering::{AVOGADRO_PER_MOL, MolecularRayleigh};
use fs_qty::semantic::{Composition, CompositionBasis};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}
fn model(water: bool, source: u8) -> MolecularRayleighMixture {
    let spec = if water {
        GasSpec::water_vapor_nist()
    } else {
        GasSpec::dry_air_ussa1976()
    };
    MolecularRayleighMixture::try_new(
        &[spec],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        &[
            MolecularRayleigh::try_new(0.0003, 2.5e25, 550.0, 0.03, ContentHash([source; 32]))
                .unwrap(),
        ],
        water.then_some(0),
    )
    .unwrap()
}
fn solved_temperature() -> GeneratedLayeredWallSolution {
    // Original steady heat equation: k=4, q=240, T(0)=T(2)=300,
    // hence T(x)=300+60x-30x². This is a prescribed thermal profile of a
    // quiescent gas; no momentum or convective equilibrium is asserted.
    let wall = LayeredWall::new(
        vec![
            PlanarThermalLayer::new("heated gas", 2.0, ConductivityTable::declared(4.0).unwrap())
                .unwrap(),
        ],
        vec![],
        [0.0; 2],
    )
    .unwrap();
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        wall.solve_with_generation(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 67,
                    kernel_id: 0,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            [300.0; 2],
            1.0,
            &[240.0],
        )
    })
    .unwrap()
}

#[test]
fn original_thermal_solution_drives_exact_rotated_density_gradient_and_hessian() {
    let solution = solved_temperature();
    let axis = Vec3::new(0.6, 0.8, 0.0);
    let temperature =
        PlanarConductionTemperatureField::try_new(&solution, 0, Point3::new(0.0, 0.0, 0.0), axis)
            .unwrap();
    let model = model(false, 67);
    let field = IsobaricGasField::try_new(
        &temperature,
        &model,
        Pressure::new(100_000.0),
        Temperature::new(300.0),
        ContentHash([68; 32]),
    )
    .unwrap();
    let p = Point3::new(0.3, 0.4, 0.0);
    let native = *model.evaluate(322.5, 100_000.0).unwrap().gas();
    assert_eq!(field.sample_state(p).unwrap(), native);
    assert_eq!(field.sample_state(p).unwrap().temperature, 322.5);
    let jet = field.sample_jet(p).unwrap();
    let slope = -native.density * 30.0 / 322.5;
    let curvature = native.density * (60.0 / 322.5 + 2.0 * (30.0 / 322.5_f64).powi(2));
    near(jet.sample.gradient_kg_m4.x, slope * 0.6, 3e-16);
    near(jet.sample.gradient_kg_m4.y, slope * 0.8, 3e-16);
    for (actual, expected) in jet
        .hessian_diagonal_kg_m5
        .into_iter()
        .zip([0.36, 0.64, 0.0])
    {
        near(actual, curvature * expected, 3e-16);
    }
    for (actual, expected) in jet
        .hessian_off_diagonal_kg_m5
        .into_iter()
        .zip([0.48, 0.0, 0.0])
    {
        near(actual, curvature * expected, 3e-16);
    }
    assert_eq!(jet.sample, field.sample_gradient(p).unwrap());
    assert!(core::ptr::eq(field.model(), &model));
    assert!(core::ptr::eq(field.temperature_source(), &temperature));
}

#[test]
fn original_conduction_to_molecular_refraction_bends_and_defocuses_native_rays() {
    let solution = solved_temperature();
    let temperature = PlanarConductionTemperatureField::try_new(
        &solution,
        0,
        Point3::new(-1.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
    )
    .unwrap();
    let model = model(false, 67);
    let field = IsobaricGasField::try_new(
        &temperature,
        &model,
        Pressure::new(100_000.0),
        Temperature::new(300.0),
        ContentHash([69; 32]),
    )
    .unwrap();
    let optics = field.bind_refraction().unwrap();
    let integration = GradientIndexOptions {
        max_step_m: 0.1,
        relative_tolerance: 1e-12,
        position_tolerance_m: 1e-12,
        direction_tolerance: 1e-12,
        optical_path_tolerance_m: 1e-12,
        ..GradientIndexOptions::default()
    };
    // Independent Lorentz--Lorenz source calculation at the hot symmetry axis.
    let refractivity = 0.0003;
    let n2_minus_one = refractivity * (2.0 + refractivity);
    let molecular_volume = (n2_minus_one / (3.0 + n2_minus_one)) / 2.5e25;
    let l = molecular_volume * (100_000.0 * AVOGADRO_PER_MOL / (R_USSA_1976 * 330.0));
    let n = ((1.0 + 2.0 * l) / (1.0 - l)).sqrt();
    let dn_dl = 3.0 / (2.0 * n * (1.0 - l).powi(2));
    let omega = (dn_dl * l * (60.0 / 330.0) / n).sqrt();
    let length = 20.0;
    let expected_area = length * (omega * length).sinh() / omega;
    for wavelength in [380.0, 550.0, 780.0] {
        let bundle = propagate_bundle(
            &optics,
            Ray {
                origin: Point3::new(0.0, 0.0, 0.0),
                dir: Vec3::new(0.0, 0.0, 1.0),
            },
            Vec3::new(1.0, 0.0, 0.0),
            wavelength,
            length,
            BundleOptions {
                integration,
                position_derivative_tolerance_m: 1e-11,
                direction_derivative_tolerance: 1e-12,
            },
        )
        .unwrap();
        // T=330-30x² gives X''=omega² X and Y''=0: a diverging thermal lens.
        let area = bundle
            .project_at_endpoint(Vec3::new(0.0, 0.0, 1.0), ProjectionOptions::default())
            .unwrap()
            .area_per_solid_angle_m2();
        near(area, expected_area, 2e-8);
        assert!(area > length * length + 0.1);
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
            origin: Point3::new(0.4, 0.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        },
        Vec3::new(1.0, 0.0, 0.0),
        550.0,
        2.0,
        integration,
    )
    .unwrap();
    assert!(
        bent.ray.origin.x > 0.400001,
        "ray must bend toward colder, denser gas"
    );
    near(
        bent.final_sample.phase_index() * bent.ray.dir.z,
        bent.initial_sample.phase_index(),
        2e-12,
    );
    assert!(!optics.is_spatially_uniform());
    assert_eq!(
        optics.sample(Point3::new(1.01, 0.0, 0.0), 550.0),
        Err(GradientIndexError::OutsideFieldDomain)
    );
}

struct FixedTemperature {
    t: f64,
    error: Option<GradientIndexError>,
}
impl TemperatureField for FixedTemperature {
    fn temperature_k(&self, _: Point3) -> Result<f64, GradientIndexError> {
        self.error.map_or(Ok(self.t), Err)
    }
}
impl TemperatureGradientField for FixedTemperature {
    fn temperature_gradient(
        &self,
        p: Point3,
    ) -> Result<TemperatureGradientSample, GradientIndexError> {
        Ok(TemperatureGradientSample {
            temperature_k: self.temperature_k(p)?,
            gradient_k_per_m: Vec3::new(0.0, 0.0, 0.0),
        })
    }
}
impl TemperatureSecondOrderField for FixedTemperature {
    fn temperature_jet(&self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        Ok(TemperatureJet {
            sample: self.temperature_gradient(p)?,
            diagonal_k_per_m2: [0.0; 3],
            off_diagonal_k_per_m2: [0.0; 3],
        })
    }
}

#[test]
fn exact_original_source_windows_labels_and_thermal_refusals_survive() {
    let p = Point3::new(0.0, 0.0, 0.0);
    for water in [false, true] {
        let model = model(water, 67);
        for t in [50.0, 300.0, 400.0, 2000.0] {
            let temperature = FixedTemperature { t, error: None };
            let field = IsobaricGasField::try_new(
                &temperature,
                &model,
                Pressure::new(10_000.0),
                Temperature::new(300.0),
                ContentHash([70; 32]),
            )
            .unwrap();
            assert_eq!(
                field.sample_state(p).unwrap(),
                *model.evaluate(t, 10_000.0).unwrap().gas()
            );
            assert_eq!(field.temperature_k(p).unwrap(), t);
            assert_eq!(
                field.sample_jet(p).unwrap().hessian_diagonal_kg_m5,
                [0.0; 3]
            );
        }
        for t in [49.9, 2000.1, f64::NAN] {
            let temperature = FixedTemperature { t, error: None };
            let field = IsobaricGasField::try_new(
                &temperature,
                &model,
                Pressure::new(10_000.0),
                Temperature::new(300.0),
                ContentHash([70; 32]),
            )
            .unwrap();
            assert_eq!(
                field.bind_refraction().unwrap().sample(p, 550.0),
                Err(GradientIndexError::OutsideFieldDomain)
            );
        }
    }
    let model = model(false, 67);
    let temperature = FixedTemperature {
        t: 300.0,
        error: Some(GradientIndexError::Cancelled),
    };
    let field = IsobaricGasField::try_new(
        &temperature,
        &model,
        Pressure::new(100_000.0),
        Temperature::new(300.0),
        ContentHash([70; 32]),
    )
    .unwrap();
    assert_eq!(
        field.bind_refraction().unwrap().sample(p, 550.0),
        Err(GradientIndexError::Cancelled)
    );
    assert_eq!(
        field.sample_density_second_order(p),
        Err(GradientIndexError::Cancelled)
    );
    for pressure in [0.0, -1.0, 10_000_001.0, f64::INFINITY, f64::NAN] {
        assert!(
            IsobaricGasField::try_new(
                &temperature,
                &model,
                Pressure::new(pressure),
                Temperature::new(300.0),
                ContentHash([70; 32])
            )
            .is_err()
        );
    }
    assert!(matches!(
        IsobaricGasField::try_new(
            &temperature,
            &model,
            Pressure::new(100_000.0),
            Temperature::new(300.0),
            ContentHash([0; 32])
        ),
        Err(IsobaricGasError::MissingSnapshotIdentity)
    ));
}
