//! Original source-temperature refresh, coherent banks and native cell energy.
use super::*;
use fs_couple::layered_wave_heating::{
    interference::{
        ThermalInterferenceDrive, ThermalInterferenceFactor, ThermalInterferenceLimits,
        ThermalInterferenceLine,
    },
    spatial_interference::feedback::{
        QuasiSteadySpatialInterferenceHeating, SpatialInterferenceFeedbackError,
        SpatialInterferenceFeedbackLimits,
    },
};
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumFactorStatistics, AngularSpectrumIlluminationLimits,
};
use fs_material::field_storage::dielectric::lorentz::{
    LorentzPoleProperties, ResolvedLorentzMaterial,
};
use fs_qty::semantic::FrequencyConvention;

fn temperature_card(upper: f64, lorentz: bool) -> MaterialCard {
    let base = card("temperature-dependent causal film", 2.7, 0.6);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        if ![MU_NAME, "delta", "tau"].contains(&claim.key.name()) {
            claims.insert_claim(claim.clone()).unwrap();
        }
    }
    for (name, dims, values) in [
        (MU_NAME, MU, [0.9 * M, 1.1 * M]),
        ("delta", EPS, [0.4 * E, 1.6 * E]),
        ("tau", Time::DIMS, [0.8 / OMEGA, 1.6 / OMEGA]),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: fs_qty::Temperature::DIMS,
                    knots: vec![(280.0, values[0]), (360.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        280.0,
                        upper,
                    )
                    .with_quantity(
                        "pressure",
                        quantity(QuantityKind::Pressure),
                        90_000.0,
                        200_000.0,
                    ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored spatial feedback oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    if lorentz {
        for (name, q, values) in [
            (
                "omega0",
                quantity(QuantityKind::Frequency(FrequencyConvention::Angular)),
                [1.35 * OMEGA, 1.05 * OMEGA],
            ),
            (
                "gamma",
                QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
                [0.12 * OMEGA, 0.25 * OMEGA],
            ),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::with_quantity(name, q),
                    value: PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: fs_qty::Temperature::DIMS,
                        knots: vec![(280.0, values[0]), (360.0, values[1])],
                        dims: q.dims(),
                    },
                    validity: ValidityDomain::unconstrained()
                        .with_quantity(
                            "T",
                            quantity(QuantityKind::AbsoluteTemperature),
                            280.0,
                            upper,
                        )
                        .with_quantity(
                            "pressure",
                            quantity(QuantityKind::Pressure),
                            90_000.0,
                            200_000.0,
                        ),
                    interpolation: InterpolationPolicy::LinearInside,
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: "authored spatial Lorentz feedback oracle".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}
fn initial_at<'t>(
    bound: &'t BoundLayeredThermalBody<'_>,
    temperatures: &[f64],
) -> DepositionSimulation<'t> {
    let mut events = bound
        .birth_events(0.0, temperatures, conduction(), || Ok(()))
        .unwrap();
    for event in &mut events {
        if let Event::Birth(cell) = &mut event.event {
            if cell.id == bound.cell_ids()[0] {
                cell.boundary.source_w = 5.0;
            }
        }
    }
    let materials: Vec<_> = bound.materials().into_iter().collect();
    DepositionSimulation::new(
        events,
        &materials,
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn feedback_limits() -> SpatialInterferenceFeedbackLimits {
    let spatial = limits();
    SpatialInterferenceFeedbackLimits {
        maximum_cells: 2,
        fields: ThermalInterferenceLimits {
            maximum_lines: 8,
            maximum_preparation_workspace_bytes: 128 << 20,
            angular: AngularSpectrumLimits {
                maximum_grid_points: 16,
                maximum_workspace_bytes: 8 << 20,
                maximum_fft_work: 1 << 20,
            },
            illumination: AngularSpectrumIlluminationLimits {
                maximum_layer_mode_pairs: 128,
                maximum_workspace_bytes: 8 << 20,
                maximum_fft_work: 1 << 20,
            },
            exposure: AngularSpectrumExposureLimits {
                maximum_lines: 8,
                maximum_point_line_pairs: 1024,
                maximum_point_frequency_pairs: 1024,
                maximum_workspace_bytes: 8 << 20,
                maximum_fft_work: 1 << 20,
            },
            ensemble: AngularSpectrumExposureEnsembleLimits {
                maximum_factors: 4,
                maximum_point_line_pairs: 1024,
                maximum_point_frequency_pairs: 1024,
                maximum_workspace_bytes: 16 << 20,
                maximum_fft_work: 1 << 20,
            },
            layer: spatial.layer,
            heat: spatial.cell_heat,
        },
    }
}

fn exercise(endpoint_domain_refusal: bool, carrier: [f64; 2], lorentz: bool) {
    let outside = state(&card("feedback-exterior", 1.0, 0.0), false, 100_000.0);
    let source_card = temperature_card(
        if endpoint_domain_refusal { 300.000001 } else { 360.0 },
        lorentz,
    );
    let source = debye(&source_card, 100_000.0);
    let resonant = lorentz.then(|| {
        let material = source
            .material()
            .with_additional_requirements(
                &source_card,
                &[
                    ScalarPropertyRequirement::try_with_quantity(
                        "omega0",
                        quantity(QuantityKind::Frequency(FrequencyConvention::Angular)),
                        ScalarAdmissibility::StrictlyPositive,
                    )
                    .unwrap(),
                    ScalarPropertyRequirement::try_new(
                        "gamma",
                        Dims([0, 0, -1, 0, 0, 0]),
                        ScalarAdmissibility::NonNegative,
                    )
                    .unwrap(),
                ],
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
        ResolvedLorentzMaterial::try_from_material(
            material,
            EPS_NAME,
            &[LorentzPoleProperties {
                permittivity_increment: "delta".into(),
                angular_frequency: "omega0".into(),
                damping_rate: "gamma".into(),
            }],
        )
        .unwrap()
    });
    let response = resonant
        .as_ref()
        .map(DielectricResponse::Lorentz)
        .unwrap_or(DielectricResponse::Debye(&source));
    let width = 0.6 / k();
    let native = body(&[&source_card], &[width], &[2], 0.2);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let original = LayeredWaveHeating::try_new_with_incident_response(
        &bound,
        DielectricResponse::Lossless(&outside),
        DielectricResponse::Lossless(&outside),
        OMEGA,
        &[response],
        || false,
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 4,
        height: 1,
        spacing_m: [1.0, 1.0],
    };
    let left = [[C64::new(20.0, 10.0), C64::from_re(100.0)]; 4];
    let right = [[C64::from_re(15.0), C64::new(-20.0, 30.0)]; 4];
    let oblique: Vec<_> = (0..4)
        .map(|i| {
            let angle = -(i as f64) * core::f64::consts::FRAC_PI_2;
            [C64::ZERO, C64::new(60.0 * angle.cos(), 60.0 * angle.sin())]
        })
        .collect();
    // Authored order differs from sorted frequency order: the correct original
    // source stack must remain the common ensemble anchor after rebuilding.
    let lines = [
        ThermalInterferenceLine {
            angular_frequency_rad_s: OMEGA + 50.0,
            left_tangential_e_v_per_m: &left,
            right_tangential_e_v_per_m: &right,
        },
        ThermalInterferenceLine {
            angular_frequency_rad_s: OMEGA,
            left_tangential_e_v_per_m: &oblique,
            right_tangential_e_v_per_m: &right,
        },
    ];
    let factors = [
        ThermalInterferenceFactor {
            lines: &[],
            statistics: AngularSpectrumFactorStatistics::Real,
        },
        ThermalInterferenceFactor {
            lines: &lines,
            statistics: AngularSpectrumFactorStatistics::Real,
        },
        ThermalInterferenceFactor {
            lines: &lines,
            statistics: AngularSpectrumFactorStatistics::ProperComplex,
        },
    ];
    let drive = ThermalInterferenceDrive {
        source_clock: "shared experiment",
        mean: &lines,
        factors: &factors,
    };
    let depth = depth(width / 2.0);
    let quadrature = [SpatialInterferenceCellQuadrature {
        depth_samples: &depth,
    }; 2];
    let footprint = AngularSpectrumAperture::try_new([0.15, -0.2], [0.5, 0.4]).unwrap();
    let cards = [&source_card];
    let owner = QuasiSteadySpatialInterferenceHeating::try_new_with_carrier(
        &original,
        &cards,
        grid,
        carrier,
        drive,
        &quadrature,
        footprint,
        feedback_limits(),
        || false,
    )
    .unwrap();
    assert_eq!(owner.carrier_transverse_wavenumber_rad_m(), carrier);
    assert_eq!(owner.frequencies_rad_s(), [OMEGA, OMEGA + 50.0]);
    let temperatures = if endpoint_domain_refusal {
        [300.0; 2]
    } else {
        [295.0, 315.0]
    };
    let initial = owner
        .initial_state(
            initial_at(&bound, &temperatures),
            Dd::from_pair(1024.0, 0.017),
            || false,
        )
        .unwrap();
    let before = initial.thermal().cells().clone();
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    let dt = 0.02;
    let temporal = SpatialInterferenceHeatingControls {
        maximum_thermal_step_s: dt,
        maximum_source_ripple_j: 100.0,
    };
    let mut polls = 0;
    let result = owner.advance_to(
        &initial,
        dt,
        2.0,
        temporal,
        conduction(),
        LayeredThermalSolveControls::default(),
        &mut workspace,
        || {
            polls += 1;
            false
        },
    );
    if endpoint_domain_refusal {
        assert!(matches!(
            result,
            Err(SpatialInterferenceFeedbackError::Sources(_))
        ));
        assert_eq!(initial.thermal().cells(), &before);
        assert_eq!(initial.thermal().time_s(), 0.0);
        return;
    }
    let result = result.unwrap();
    let legacy = QuasiSteadySpatialInterferenceHeating::try_new(
        &original,
        &cards,
        grid,
        drive,
        &quadrature,
        footprint,
        feedback_limits(),
        || false,
    )
    .unwrap();
    assert_eq!(legacy.carrier_transverse_wavenumber_rad_m(), [0.0; 2]);
    let zero_carrier = legacy
        .advance_to(
            &initial,
            dt,
            2.0,
            temporal,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    if carrier == [0.0; 2] {
        assert_eq!(result.step.cell_energy_j, zero_carrier.step.cell_energy_j);
        assert_eq!(
            result.step.state.thermal().cells(),
            zero_carrier.step.state.thermal().cells()
        );
    } else {
        assert!(
            result
                .step
                .cell_energy_j
                .iter()
                .zip(&zero_carrier.step.cell_energy_j)
                .any(|((_, tilted), (_, normal))| (tilted - normal).abs() > 1e-4 * tilted.abs())
        );
    }
    // Same geometry, drive and phase epoch, but different accepted native
    // temperatures: the applied optical joules themselves must change. This
    // falsifies a refresh that updates metadata while leaving fields frozen.
    let uniform = owner
        .initial_state(
            initial_at(&bound, &[300.0; 2]),
            initial.source_time_s(),
            || false,
        )
        .unwrap();
    let uniform_step = owner
        .advance_to(
            &uniform,
            dt,
            2.0,
            temporal,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    assert!(
        result
            .step
            .cell_energy_j
            .iter()
            .zip(&uniform_step.step.cell_energy_j)
            .any(|((_, warm), (_, held))| (warm - held).abs() > 1e-5 * warm.abs())
    );
    for (index, &temperature) in result.applied_temperatures_k.iter().enumerate() {
        let reference = response
            .material()
            .at_query_point(&source_card, &at(temperature, 100_000.0))
            .unwrap();
        let applied = result.applied_response(index).unwrap();
        assert_eq!(
            applied.material().identity(),
            reference.identity()
        );
        assert_eq!(applied.material().state_coordinate("T"), Some(temperature));
        assert_eq!(
            applied.material().state_coordinate("pressure"),
            Some(100_000.0)
        );
        near(
            fs_material::field_storage::linear_magnetic_permeability(applied.material()).unwrap(),
            fs_material::field_storage::linear_magnetic_permeability(&reference).unwrap(),
            1e-20,
        );
        if lorentz {
            // Original authored coefficients, independently eliminated at each
            // accepted local temperature; neither a layer mean nor T0 is used.
            let fraction = (temperature - 280.0) / 80.0;
            let omega0 = (1.35 - 0.3 * fraction) * OMEGA;
            let gamma = (0.12 + 0.13 * fraction) * OMEGA;
            let delta = (0.4 + 1.2 * fraction) * E;
            for omega in [OMEGA, OMEGA + 50.0] {
                let expected = C64::from_re(2.1 * E)
                    + C64::from_re(delta * omega0.powi(2))
                        / C64::new(omega0.powi(2) - omega.powi(2), gamma * omega);
                let epsilon = applied.permittivity_f_m(omega).unwrap();
                near(epsilon.re, expected.re, 1e-24);
                near(epsilon.im, expected.im, 1e-24);
            }
        }
        let endpoint = result
            .step
            .state
            .thermal()
            .temperature_k(original.cell_ids()[index])
            .unwrap();
        assert_eq!(
            result
                .endpoint_response(index)
                .unwrap()
                .material()
                .state_coordinate("T"),
            Some(endpoint)
        );
    }
    assert_ne!(
        result.applied_response(0).unwrap().material().identity(),
        result.applied_response(1).unwrap().material().identity()
    );
    // Independent two-cell constant-capacity backward-Euler balance, driven by
    // actual integrated optical joules plus the unchanged authored source.
    let capacity = RHO * 0.2 * width / 2.0 * CP;
    let conductance_dt = 0.2 * K / (width / 2.0) * dt;
    let q = [
        result.step.cell_energy_j[0].1 + 5.0 * dt,
        result.step.cell_energy_j[1].1,
    ];
    let t = &result.applied_temperatures_k;
    let mean = 0.5 * (t[0] + t[1]) + 0.5 * (q[0] + q[1]) / capacity;
    let contrast =
        (0.5 * capacity * (t[0] - t[1]) + 0.5 * (q[0] - q[1])) / (capacity + 2.0 * conductance_dt);
    for (index, expected) in [mean + contrast, mean - contrast].into_iter().enumerate() {
        near(result.endpoint_temperatures_k[index], expected, 2e-8);
    }
    near(
        result.step.heat_report.energy.external_heat_j,
        result.step.absorbed_energy_j.to_f64() + 5.0 * dt,
        1e-8,
    );
    near(
        result.step.heat_report.energy.stored_energy_change_j,
        result.step.heat_report.energy.external_heat_j,
        2e-8,
    );
    near(
        (result.step.state.source_time_s() - initial.source_time_s()).to_f64(),
        dt,
        1e-17,
    );
    let next = owner
        .advance_to(
            &result.step.state,
            2.0 * dt,
            2.0,
            temporal,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    for index in 0..2 {
        assert_eq!(
            next.applied_response(index).unwrap().material().identity(),
            result
                .endpoint_response(index)
                .unwrap()
                .material()
                .identity()
        );
    }
    near(
        (next.step.state.source_time_s() - initial.source_time_s()).to_f64(),
        2.0 * dt,
        1e-17,
    );
    assert_ne!(next.step.cell_energy_j, result.step.cell_energy_j);
    assert!(matches!(
        owner.advance_to(
            &initial,
            dt,
            0.0,
            temporal,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || false
        ),
        Err(SpatialInterferenceFeedbackError::TemperatureExcursion { .. })
    ));
    let mut count = 0;
    assert!(
        owner
            .advance_to(
                &initial,
                dt,
                2.0,
                temporal,
                conduction(),
                LayeredThermalSolveControls::default(),
                &mut workspace,
                || {
                    count += 1;
                    count == polls
                }
            )
            .is_err()
    );
    assert_eq!(initial.thermal().cells(), &before);
    assert_eq!(initial.thermal().time_s(), 0.0);
}

#[test]
fn nonuniform_native_cells_refresh_causal_sources_without_resetting_mixed_bank_clocks() {
    exercise(false, [0.0; 2], false);
}

#[test]
fn off_axis_envelopes_heat_finite_cells_and_keep_their_carrier_through_temperature_rebuilds() {
    exercise(false, [0.65, -0.4], false);
}

#[test]
fn endpoint_optical_domain_refusal_leaves_the_original_cell_history_unchanged() {
    exercise(true, [0.0; 2], false);
}

#[test]
fn local_lorentz_resonances_refract_and_heat_coherent_fields_on_the_native_thermal_clock() {
    exercise(false, [0.65, -0.4], true);
}

#[test]
fn late_lorentz_source_refusal_preserves_both_native_cells_and_the_original_phase_clock() {
    exercise(true, [0.0; 2], true);
}
