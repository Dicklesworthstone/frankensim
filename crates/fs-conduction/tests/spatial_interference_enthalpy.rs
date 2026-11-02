//! Original-owner regressions for spatial spectral heating and exact source joules.
mod support;

use fs_blake3::ContentHash;
use fs_conduction::transient::enthalpy::interference::{
    NodalPowerSpectrum, PreparedSpatialInterference, SpatialInterferenceBudget,
    SpatialInterferenceControls, SpatialInterferenceError, SpatialInterferenceRadiation,
};
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyError, EnthalpyStepConfig,
};
use fs_conduction::{
    AmbientRadiationPatch, ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    ConductivityTable, EMISSIVITY_DIMS, PointwiseAmbientRadiation, SURFACE_EMISSIVITY_PROPERTY,
    ScalarField, SurfaceEmissivity, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder,
};
use fs_evidence::ValidityDomain;
use fs_fft::harmonic_window::{
    AngularHarmonicSpectrumTerm, HarmonicSpectrumTermLimits, PreparedHarmonicProductSpectrum,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_math::dd::Dd;
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

const NAMES: [&str; 4] = ["vertex-a", "vertex-b", "vertex-c", "vertex-d"];

fn mesh() -> ConductionMesh {
    ConductionMesh::new(
        TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    )
    .unwrap()
}

fn adiabatic(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap()
}

fn fixed(mesh: &ConductionMesh, temperature: f64) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region("bath", |_| true, ThermalBc::dirichlet(temperature).unwrap())
        .unwrap()
        .finish()
        .unwrap()
}

fn curve(latent: bool) -> EquilibriumEnthalpyPhaseCurve {
    let points = if latent {
        vec![
            (0.0, 300.0, 0.0, 4.0),
            (100.0, 350.0, 0.0, 3.0),
            (300.0, 350.0, 1.0, 2.0),
            (600.0, 500.0, 1.0, 1.0),
        ]
    } else {
        vec![(0.0, 280.0, 0.0, 2.0), (1100.0, 500.0, 0.0, 2.0)]
    };
    let knots = points
        .into_iter()
        .map(|(h, t, f, rho)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: rho,
        })
        .collect();
    if latent {
        EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x71; 32]), knots).unwrap()
    } else {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x72; 32]),
            SolidLiquidPhase::Solid,
            knots,
        )
        .unwrap()
    }
}

fn storage_budget() -> EnthalpyBudget {
    EnthalpyBudget {
        max_vertices: 4,
        max_elements: 1,
    }
}

fn source_budget() -> SpatialInterferenceBudget {
    SpatialInterferenceBudget {
        maximum_bindings: 4,
        maximum_spectral_terms: 32,
        maximum_source_workspace_bytes: 16_384,
    }
}

fn config() -> EnthalpyStepConfig {
    EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-12,
            linear_restart: 4,
            max_linear_cycles: 8,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-3,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 32,
        energy_tolerance_j: 1e-8,
    }
}

fn controls() -> SpatialInterferenceControls {
    SpatialInterferenceControls {
        maximum_thermal_step_s: 2.0,
        maximum_source_ripple_j: 100.0,
    }
}

fn problem<'a>(
    mesh: &'a ConductionMesh,
    boundary: &'a ThermalBoundary,
    material: &'a ConductivityModel,
    source: &'a ScalarField,
) -> ConductionProblem<'a> {
    ConductionProblem {
        mesh,
        boundary,
        material,
        source,
        element_materials: None,
    }
}

fn spectrum(dc: f64, omega: f64, real: f64, imaginary: f64) -> PreparedHarmonicProductSpectrum {
    PreparedHarmonicProductSpectrum::from_angular_terms(
        [
            AngularHarmonicSpectrumTerm {
                angular_frequency_rad_s: Dd::ZERO,
                peak_amplitude: [Dd::from_f64(dc), Dd::ZERO],
            },
            AngularHarmonicSpectrumTerm {
                angular_frequency_rad_s: Dd::from_f64(omega),
                peak_amplitude: [Dd::from_f64(real), Dd::from_f64(imaginary)],
            },
        ],
        HarmonicSpectrumTermLimits {
            maximum_terms: 2,
            maximum_workspace_bytes: 4096,
        },
    )
    .unwrap()
}

fn bindings(
    spectra: &[PreparedHarmonicProductSpectrum; 4],
    scale: f64,
) -> Vec<NodalPowerSpectrum<'_>> {
    spectra
        .iter()
        .enumerate()
        .map(|(vertex, spectrum)| NodalPowerSpectrum {
            name: NAMES[vertex],
            vertex,
            spectrum,
            input_unit: "declared spectral power density",
            watts_per_m3_per_input_unit: scale,
        })
        .collect()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn prescribed_nodes_use_consistent_p1_load_and_exact_reservoir_reaction() {
    let mesh = mesh();
    let boundary = fixed(&mesh, 300.0);
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let zero = ScalarField::Uniform(0.0);
    let curve = curve(false);
    let spectra = [1.0, 2.0, 3.0, 4.0].map(|dc| spectrum(dc, 0.0, 0.0, 0.0));
    with_cx(|cx| {
        let engine =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 2.0, storage_budget()).unwrap();
        let prepared = PreparedSpatialInterference::new(
            cx,
            &engine,
            &bindings(&spectra, 1.0),
            "laboratory clock",
            source_budget(),
        )
        .unwrap();
        let got = prepared
            .advance_interval(
                cx,
                problem(&mesh, &boundary, &material, &zero),
                None,
                &[100.0; 4],
                Dd::ZERO,
                1.0,
                config(),
                controls(),
                None,
            )
            .unwrap();
        // V=1/6. The exact P1 dual row is V*(sum(E)+E_i)/20,
        // while the independently integrated physical total is V*sum(E)/4.
        for (i, q) in got.source_load_j.iter().enumerate() {
            close(q.to_f64(), (11 + i) as f64 / 120.0, 2e-17);
        }
        close(got.energy.source_j.to_f64(), 5.0 / 12.0, 6e-17);
        close(got.energy.dirichlet_in_j.to_f64(), -5.0 / 12.0, 1e-12);
        assert_eq!(got.energy.stored_energy_change_j.to_f64(), 0.0);
        assert!(got.nodal_residual_j.iter().all(|r| r.to_f64() == 0.0));
        assert_eq!(got.specific_enthalpy_j_kg, vec![100.0; 4]);
        assert!(got.energy.residual_j.to_f64().abs() < 1e-12);
        assert_eq!(got.newton.iterations, 0);
    });
}

#[test]
fn spatially_nonuniform_dc_reuses_original_nonlinear_conductivity_endpoint() {
    let mesh = mesh();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(280.0, 1.0), (500.0, 11.0)]).unwrap(),
    );
    let curve = curve(false);
    let q = [20.0, -4.0, 12.0, 0.0];
    let static_source = ScalarField::Nodal(q.to_vec());
    let zero = ScalarField::Uniform(0.0);
    let old = [150.0, 250.0, 100.0, 300.0];
    let spectra = q.map(|dc| spectrum(dc, 0.0, 0.0, 0.0));
    with_cx(|cx| {
        let engine =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 2.0, storage_budget()).unwrap();
        let native = engine
            .advance(
                cx,
                problem(&mesh, &boundary, &material, &static_source),
                None,
                &old,
                0.03,
                config(),
            )
            .unwrap();
        let prepared = PreparedSpatialInterference::new(
            cx,
            &engine,
            &bindings(&spectra, 1.0),
            "source",
            source_budget(),
        )
        .unwrap();
        let got = prepared
            .advance_interval(
                cx,
                problem(&mesh, &boundary, &material, &zero),
                None,
                &old,
                Dd::ZERO,
                0.03,
                config(),
                controls(),
                None,
            )
            .unwrap();
        for (a, b) in got
            .specific_enthalpy_j_kg
            .iter()
            .zip(&native.specific_enthalpy_j_kg)
        {
            close(*a, *b, 1e-8);
        }
        close(
            got.energy.source_j.to_f64(),
            0.03 * (20.0 - 4.0 + 12.0) / 24.0,
            1e-17,
        );
        close(
            got.energy.stored_energy_change_j.to_f64(),
            got.energy.source_j.to_f64(),
            1e-8,
        );
        assert!(
            got.temperature
                .windows(2)
                .any(|pair| (pair[0] - pair[1]).abs() > 0.1)
        );
        assert_eq!(old, [150.0, 250.0, 100.0, 300.0]);
    });
}

#[test]
fn local_heating_changes_phase_at_constant_temperature_without_changing_reference_mass() {
    let mesh = mesh();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let curve = curve(true);
    let zero = ScalarField::Uniform(0.0);
    let spectra = [100.0, 0.0, 0.0, 0.0].map(|dc| spectrum(dc, 0.0, 0.0, 0.0));
    with_cx(|cx| {
        let engine =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 2.0, storage_budget()).unwrap();
        let before = engine.reference_nodal_masses_kg().to_vec();
        let prepared = PreparedSpatialInterference::new(
            cx,
            &engine,
            &bindings(&spectra, 1.0),
            "source",
            source_budget(),
        )
        .unwrap();
        let got = prepared
            .advance_interval(
                cx,
                problem(&mesh, &boundary, &material, &zero),
                None,
                &[150.0; 4],
                Dd::ZERO,
                1.0,
                config(),
                controls(),
                None,
            )
            .unwrap();
        for (i, expected) in [170.0, 160.0, 160.0, 160.0].into_iter().enumerate() {
            close(got.specific_enthalpy_j_kg[i], expected, 1e-9);
            close(got.temperature[i], 350.0, 1e-12);
            close(
                got.liquid_mass_fraction[i],
                (expected - 100.0) / 200.0,
                1e-11,
            );
        }
        close(got.energy.source_j.to_f64(), 100.0 / 24.0, 1e-15);
        close(
            got.energy.stored_energy_change_j.to_f64(),
            100.0 / 24.0,
            1e-9,
        );
        assert_eq!(engine.reference_nodal_masses_kg(), before);
    });
}

#[test]
fn opposite_local_carriers_cannot_hide_behind_zero_total_energy_or_full_cycles() {
    let mesh = mesh();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(0.01).unwrap();
    let curve = curve(false);
    let zero = ScalarField::Uniform(0.0);
    let omega = 4.0 * std::f64::consts::PI;
    let spectra = [
        spectrum(0.0, omega, 1.0, 0.0),
        spectrum(0.0, omega, -1.0, 0.0),
        spectrum(0.0, 0.0, 0.0, 0.0),
        spectrum(0.0, 0.0, 0.0, 0.0),
    ];
    with_cx(|cx| {
        let engine =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 2.0, storage_budget()).unwrap();
        let prepared = PreparedSpatialInterference::new(
            cx,
            &engine,
            &bindings(&spectra, 1.0),
            "source",
            source_budget(),
        )
        .unwrap();
        let controls = SpatialInterferenceControls {
            maximum_thermal_step_s: 1.0,
            maximum_source_ripple_j: 2e-5,
        };
        assert!(matches!(
            prepared.advance_interval(
                cx,
                problem(&mesh, &boundary, &material, &zero),
                None,
                &[100.0; 4],
                Dd::ZERO,
                1.0,
                config(),
                controls,
                None,
            ),
            Err(SpatialInterferenceError::SourceRipple { .. })
        ));
        let fine = prepared
            .advance_interval(
                cx,
                problem(&mesh, &boundary, &material, &zero),
                None,
                &[100.0; 4],
                Dd::ZERO,
                0.01,
                config(),
                controls,
                None,
            )
            .unwrap();
        assert_eq!(fine.energy.source_j.to_f64(), 0.0);
        assert!(fine.specific_enthalpy_j_kg[0] > fine.specific_enthalpy_j_kg[1]);
        assert!(fine.source_ripple_bound_j > 0.0 && fine.source_ripple_bound_j < 2e-5);
    });
}

#[test]
fn high_frequency_input_is_admitted_only_when_its_integrated_ripple_is_physically_small() {
    let mesh = mesh();
    let boundary = fixed(&mesh, 300.0);
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let curve = curve(false);
    let zero = ScalarField::Uniform(0.0);
    let spectra = std::array::from_fn(|_| spectrum(0.0, 1e9, 1e-3, 0.0));
    with_cx(|cx| {
        let engine =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 2.0, storage_budget()).unwrap();
        let prepared = PreparedSpatialInterference::new(
            cx,
            &engine,
            &bindings(&spectra, 1.0),
            "source",
            source_budget(),
        )
        .unwrap();
        let got = prepared
            .advance_interval(
                cx,
                problem(&mesh, &boundary, &material, &zero),
                None,
                &[100.0; 4],
                Dd::ZERO,
                1.0,
                config(),
                SpatialInterferenceControls {
                    maximum_thermal_step_s: 1.0,
                    maximum_source_ripple_j: 1e-10,
                },
                None,
            )
            .unwrap();
        assert!(got.source_ripple_bound_j > 0.0 && got.source_ripple_bound_j < 1e-10);
        assert_eq!(got.specific_enthalpy_j_kg, vec![100.0; 4]);
    });
}

#[test]
fn dd_clock_chains_late_intervals_without_losing_source_phase() {
    let mesh = mesh();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let curve = curve(false);
    let zero = ScalarField::Uniform(0.0);
    let spectra = std::array::from_fn(|_| spectrum(20.0, 1.0, 5.0, 3.0));
    let epoch = Dd::from_f64(2_f64.powi(80));
    with_cx(|cx| {
        let engine =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 2.0, storage_budget()).unwrap();
        let prepared = PreparedSpatialInterference::new(
            cx,
            &engine,
            &bindings(&spectra, 1.0),
            "absolute source epoch",
            source_budget(),
        )
        .unwrap();
        let p = problem(&mesh, &boundary, &material, &zero);
        let first = prepared
            .advance_interval(
                cx,
                p,
                None,
                &[100.0; 4],
                epoch,
                0.25,
                config(),
                controls(),
                None,
            )
            .unwrap();
        let second = prepared
            .advance_interval(
                cx,
                p,
                None,
                &first.specific_enthalpy_j_kg,
                first.end_clock_s,
                0.25,
                config(),
                controls(),
                None,
            )
            .unwrap();
        assert_eq!(first.end_clock_s.hi, epoch.hi);
        close(first.end_clock_s.lo, 0.25, 0.0);
        close(second.end_clock_s.lo, 0.5, 0.0);
        let whole = spectra[0]
            .causal_integral_from_clock_dd(epoch, 0.5, 0.0)
            .unwrap()
            .integral_dd()
            / Dd::from_f64(6.0);
        let parts = first.energy.source_j + second.energy.source_j;
        assert!((parts - whole).to_f64().abs() < 1e-13);
        let expected_second = spectra[0]
            .causal_integral_from_clock_dd(Dd::from_pair(epoch.hi, 0.25), 0.25, 0.0)
            .unwrap()
            .integral_dd()
            / Dd::from_f64(6.0);
        assert!((second.energy.source_j - expected_second).to_f64().abs() < 1e-13);
    });
}

#[test]
fn physical_unit_conversion_precedes_tiny_duration_and_survives_the_spatial_reservoir_ledger() {
    let mesh = mesh();
    let boundary = fixed(&mesh, 300.0);
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let curve = curve(false);
    let zero = ScalarField::Uniform(0.0);
    let spectra = std::array::from_fn(|_| spectrum(1e-200, 0.0, 0.0, 0.0));
    with_cx(|cx| {
        let engine =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 2.0, storage_budget()).unwrap();
        let prepared = PreparedSpatialInterference::new(
            cx,
            &engine,
            &bindings(&spectra, 1e200),
            "source",
            source_budget(),
        )
        .unwrap();
        let mut config = config();
        config.energy_tolerance_j = 1e-210;
        let got = prepared
            .advance_interval(
                cx,
                problem(&mesh, &boundary, &material, &zero),
                None,
                &[100.0; 4],
                Dd::ZERO,
                1e-200,
                config,
                controls(),
                None,
            )
            .unwrap();
        let expected = 1e-200 / 6.0;
        assert!(got.energy.source_j.to_f64() > 0.0);
        close(got.energy.source_j.to_f64() / expected, 1.0, 2e-14);
        close(got.energy.dirichlet_in_j.to_f64() / expected, -1.0, 1e-12);
        assert_eq!(got.energy.stored_energy_change_j.to_f64(), 0.0);
        assert!(got.energy.residual_j.to_f64().abs() < 1e-210);
    });
}

#[test]
fn original_mesh_and_single_source_ownership_are_mandatory() {
    let mesh = mesh();
    let other_mesh = self::mesh();
    let boundary = adiabatic(&mesh);
    let other_boundary = adiabatic(&other_mesh);
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let curve = curve(false);
    let zero = ScalarField::Uniform(0.0);
    let duplicate = ScalarField::Uniform(1.0);
    let spectra = std::array::from_fn(|_| spectrum(1.0, 0.0, 0.0, 0.0));
    let old = [100.0; 4];
    with_cx(|cx| {
        let engine =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 2.0, storage_budget()).unwrap();
        let prepared = PreparedSpatialInterference::new(
            cx,
            &engine,
            &bindings(&spectra, 1.0),
            "source",
            source_budget(),
        )
        .unwrap();
        for p in [
            problem(&other_mesh, &other_boundary, &material, &zero),
            problem(&mesh, &boundary, &material, &duplicate),
        ] {
            assert!(matches!(
                prepared.advance_interval(
                    cx,
                    p,
                    None,
                    &old,
                    Dd::ZERO,
                    0.1,
                    config(),
                    controls(),
                    None,
                ),
                Err(SpatialInterferenceError::Enthalpy(
                    EnthalpyError::InvalidInput(_)
                ))
            ));
        }
        assert_eq!(old, [100.0; 4]);
    });
}

#[test]
fn mapping_units_budgets_clocks_and_step_domains_refuse_without_an_endpoint() {
    let mesh = mesh();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let curve = curve(false);
    let zero = ScalarField::Uniform(0.0);
    let spectra = std::array::from_fn(|_| spectrum(1.0, 0.0, 0.0, 0.0));
    with_cx(|cx| {
        let engine =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 2.0, storage_budget()).unwrap();
        let good = bindings(&spectra, 1.0);
        for which in 0..5 {
            let mut bad = good.clone();
            match which {
                0 => bad[1].vertex = 0,
                1 => bad[1].name = bad[0].name,
                2 => bad[1].watts_per_m3_per_input_unit = 0.0,
                3 => bad[1].input_unit = "",
                _ => bad[1].vertex = 4,
            }
            assert!(
                PreparedSpatialInterference::new(cx, &engine, &bad, "source", source_budget())
                    .is_err()
            );
        }
        for budget in [
            SpatialInterferenceBudget {
                maximum_bindings: 3,
                ..source_budget()
            },
            SpatialInterferenceBudget {
                maximum_spectral_terms: 3,
                ..source_budget()
            },
            SpatialInterferenceBudget {
                maximum_source_workspace_bytes: 1,
                ..source_budget()
            },
        ] {
            assert!(
                PreparedSpatialInterference::new(cx, &engine, &good, "source", budget).is_err()
            );
        }
        let prepared =
            PreparedSpatialInterference::new(cx, &engine, &good, "source", source_budget())
                .unwrap();
        for (clock, dt) in [
            (Dd::ZERO, 0.0),
            (Dd::ZERO, 3.0),
            (Dd::ZERO, f64::NAN),
            (Dd::from_pair(f64::NAN, 0.0), 0.1),
        ] {
            assert!(
                prepared
                    .advance_interval(
                        cx,
                        problem(&mesh, &boundary, &material, &zero),
                        None,
                        &[100.0; 4],
                        clock,
                        dt,
                        config(),
                        controls(),
                        None,
                    )
                    .is_err()
            );
        }
    });
}

#[test]
fn empty_mapping_repeats_deterministically_and_native_domain_failure_is_atomic() {
    let mesh = mesh();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let curve = curve(false);
    let zero = ScalarField::Uniform(0.0);
    let excessive = std::array::from_fn(|_| spectrum(1e10, 0.0, 0.0, 0.0));
    let old = [150.0, 200.0, 100.0, 250.0];
    with_cx(|cx| {
        let engine =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 2.0, storage_budget()).unwrap();
        let p = problem(&mesh, &boundary, &material, &zero);
        let prepared =
            PreparedSpatialInterference::new(cx, &engine, &[], "source", source_budget()).unwrap();
        let one = prepared
            .advance_interval(cx, p, None, &old, Dd::ZERO, 0.1, config(), controls(), None)
            .unwrap();
        let two = prepared
            .advance_interval(cx, p, None, &old, Dd::ZERO, 0.1, config(), controls(), None)
            .unwrap();
        assert_eq!(
            one.specific_enthalpy_j_kg
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>(),
            two.specific_enthalpy_j_kg
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>()
        );
        assert_eq!(one.energy.source_j.to_f64(), 0.0);
        let hot = PreparedSpatialInterference::new(
            cx,
            &engine,
            &bindings(&excessive, 1.0),
            "source",
            source_budget(),
        )
        .unwrap();
        assert!(
            hot.advance_interval(cx, p, None, &old, Dd::ZERO, 1.0, config(), controls(), None)
                .is_err()
        );
        assert_eq!(old, [150.0, 200.0, 100.0, 250.0]);
    });
}

#[test]
fn cancellation_preserves_the_original_history_and_prepared_sources() {
    let mesh = mesh();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let curve = curve(false);
    let zero = ScalarField::Uniform(0.0);
    let spectra = std::array::from_fn(|_| spectrum(10.0, 1.0, 2.0, 0.0));
    let old = [100.0; 4];
    with_cx(|cx| {
        let engine =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 2.0, storage_budget()).unwrap();
        let prepared = PreparedSpatialInterference::new(
            cx,
            &engine,
            &bindings(&spectra, 1.0),
            "source",
            source_budget(),
        )
        .unwrap();
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                prepared.advance_interval(
                    cancelled,
                    problem(&mesh, &boundary, &material, &zero),
                    None,
                    &old,
                    Dd::ZERO,
                    0.1,
                    config(),
                    controls(),
                    None,
                ),
                Err(SpatialInterferenceError::Enthalpy(
                    EnthalpyError::Conduction(ConductionError::Cancelled { .. })
                ))
            ));
        });
        assert_eq!(old, [100.0; 4]);
        assert_eq!(prepared.bindings().len(), 4);
    });
}

fn radiation_patch() -> AmbientRadiationPatch {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: 0.4,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 240.0, 600.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "synthetic spatial interference regression".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic gray surface".into(),
            phase: "solid".into(),
            process: "numerical reference".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    let emissivity =
        SurfaceEmissivity::from_card("exterior", &card, 350.0, SelectionPolicy::SingleClaimOnly)
            .unwrap();
    AmbientRadiationPatch::new("exterior", emissivity, 280.0).unwrap()
}

#[test]
fn original_pointwise_radiation_and_robin_remain_separate_in_the_same_newton_step() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .remainder("exterior", ThermalBc::robin(1.0, 290.0).unwrap())
        .unwrap()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(280.0, 1.0), (500.0, 5.0)]).unwrap(),
    );
    let curve = curve(false);
    let q = [2.0, 4.0, 0.0, 0.0];
    let static_source = ScalarField::Nodal(q.to_vec());
    let zero = ScalarField::Uniform(0.0);
    let spectra = q.map(|dc| spectrum(dc, 0.0, 0.0, 0.0));
    let patches = [radiation_patch()];
    let old = [150.0, 200.0, 250.0, 300.0];
    with_cx(|cx| {
        let engine =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 2.0, storage_budget()).unwrap();
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let native = engine
            .advance_with_pointwise_ambient_radiation(
                cx,
                problem(&mesh, &boundary, &material, &static_source),
                None,
                &old,
                0.002,
                config(),
                &radiation,
            )
            .unwrap();
        let prepared = PreparedSpatialInterference::new(
            cx,
            &engine,
            &bindings(&spectra, 1.0),
            "source",
            source_budget(),
        )
        .unwrap();
        let got = prepared
            .advance_interval(
                cx,
                problem(&mesh, &boundary, &material, &zero),
                None,
                &old,
                Dd::ZERO,
                0.002,
                config(),
                controls(),
                Some(SpatialInterferenceRadiation::Ambient(&radiation)),
            )
            .unwrap();
        for (a, b) in got.temperature.iter().zip(&native.temperature) {
            close(*a, *b, 1e-8);
        }
        assert!(got.energy.robin_out_j.to_f64() > 0.0);
        assert!(got.energy.radiation_out_j.to_f64() > got.energy.robin_out_j.to_f64());
        close(
            got.energy.radiation_out_j.to_f64(),
            0.002 * got.pointwise_radiation.as_ref().unwrap().outward_heat_w,
            1e-12,
        );
        close(got.energy.source_j.to_f64(), 0.002 * 6.0 / 24.0, 1e-18);
        assert!(got.energy.residual_j.to_f64().abs() < config().energy_tolerance_j);
        let other_boundary = ThermalBoundaryBuilder::new(&mesh)
            .remainder("exterior", ThermalBc::robin(1.0, 290.0).unwrap())
            .unwrap()
            .finish()
            .unwrap();
        let foreign = PointwiseAmbientRadiation::new(cx, &mesh, &other_boundary, &patches).unwrap();
        assert!(
            prepared
                .advance_interval(
                    cx,
                    problem(&mesh, &boundary, &material, &zero),
                    None,
                    &old,
                    Dd::ZERO,
                    0.002,
                    config(),
                    controls(),
                    Some(SpatialInterferenceRadiation::Ambient(&foreign)),
                )
                .is_err()
        );
    });
}
