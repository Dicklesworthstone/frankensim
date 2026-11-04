#![cfg(feature = "layered-wave-heating")]
//! Actual-owner G1/G3/G4 feedback tests with authored constitutive data.
//! These are numerical model oracles, not measured-material qualification.

use fs_conduction::lumped::{
    BiotGate, LumpedEnthalpyBody, LumpedEnthalpyMarchConfig,
    interference::enthalpy::LumpedEnthalpySpectrumControls,
};
use fs_couple::layered_wave_heating::interference::{
    QuasiSteadyFilmHeating, ThermalInterferenceDrive, ThermalInterferenceFactor,
    ThermalInterferenceLimits, ThermalInterferenceLine,
};
use fs_em::plane_wave::{
    DielectricResponse, Polarization,
    angular_spectrum::{
        AngularSpectrumExposureEnsembleLimits, AngularSpectrumExposureLayerLimits,
        AngularSpectrumExposureLimits, AngularSpectrumFactorStatistics, AngularSpectrumGrid,
        AngularSpectrumHeatAreaQuadrature, AngularSpectrumHeatDepthSample,
        AngularSpectrumHeatVolumeQuadrature, AngularSpectrumIlluminationLimits,
        AngularSpectrumIntegratedHeatLimits, AngularSpectrumLimits,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS as EPS, ELECTRIC_PERMITTIVITY_PROPERTY as EPS_NAME,
        MAGNETIC_PERMEABILITY_DIMS as MU, MAGNETIC_PERMEABILITY_PROPERTY as MU_NAME,
        dielectric::{DebyePoleProperties, ResolvedDebyeCapacitor},
    },
    phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase},
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_math::{c64::C64, dd::Dd};
use fs_qty::{
    Dims, QuantitySpec, Time,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const WIDTH: f64 = 0.25;
const RHO: f64 = 2.0;
const CP: f64 = 4.0;
const GRID: AngularSpectrumGrid = AngularSpectrumGrid {
    width: 1,
    height: 1,
    spacing_m: [1.0, 1.0],
};
const LEFT: [[[C64; 2]; 1]; 6] = [
    [[C64::new(5.0, 0.5), C64::new(1.0, -0.2)]],
    [[C64::new(2.0, -0.3), C64::new(-0.4, 0.7)]],
    [[C64::new(1.1, 0.4), C64::new(0.0, 0.3)]],
    [[C64::new(-0.2, 0.6), C64::new(0.5, 0.1)]],
    [[C64::new(0.7, -0.2), C64::new(0.2, 0.3)]],
    [[C64::new(0.4, 0.1), C64::new(-0.3, 0.4)]],
];
const RIGHT: [[[C64; 2]; 1]; 6] = [
    [[C64::new(1.0, -0.4), C64::new(0.3, 0.2)]],
    [[C64::new(-0.6, 0.2), C64::new(0.1, 0.4)]],
    [[C64::new(0.2, 0.1), C64::new(-0.1, 0.2)]],
    [[C64::new(0.1, -0.3), C64::new(0.2, 0.0)]],
    [[C64::new(-0.2, 0.3), C64::new(0.1, -0.2)]],
    [[C64::new(0.3, -0.1), C64::new(0.2, 0.1)]],
];

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn at(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            quantity(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
        .unwrap()
}
fn specs(film: bool) -> Vec<(&'static str, Dims, [f64; 2])> {
    if !film {
        return vec![(EPS_NAME, EPS, [1.0; 2]), (MU_NAME, MU, [1.0; 2])];
    }
    vec![
        (EPS_NAME, EPS, [1.8, 2.6]),
        (MU_NAME, MU, [0.9, 1.3]),
        ("delta", EPS, [0.4, 1.2]),
        ("tau", Time::DIMS, [0.1, 0.3]),
        (
            "electrical_conductivity",
            ElectricalPropertyBasis::Conductivity.dimensions(),
            [0.2, 1.0],
        ),
    ]
}
fn card(name: &str, upper: f64, film: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in specs(film) {
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
                        110_000.0,
                    ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored interference/material-feedback oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "numerical constitutive oracle".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn state(card: &MaterialCard, film: bool) -> ResolvedMaterialStatePoint {
    let requirements: Vec<_> = specs(film)
        .into_iter()
        .map(|(name, dims, _)| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap()
        })
        .collect();
    resolve_material_state_point(
        card,
        &at(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn sources(card: &MaterialCard) -> (ResolvedDebyeCapacitor, ResolvedConductorState) {
    let material = state(card, true);
    let conductor = ResolvedConductorState::try_from_material_with_basis(
        material.clone(),
        ElectricalPropertyBasis::Conductivity,
    )
    .unwrap();
    let dielectric = ResolvedDebyeCapacitor::try_from_material(
        material,
        EPS_NAME,
        &[DebyePoleProperties {
            permittivity_increment: "delta".into(),
            relaxation_time: "tau".into(),
        }],
        1.0,
        WIDTH,
    )
    .unwrap();
    (dielectric, conductor)
}
fn chart(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        [280.0, 300.0, 360.0]
            .map(|temperature| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: CP * (temperature - 300.0),
                temperature_k: temperature,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: RHO,
            })
            .to_vec(),
    )
    .unwrap()
}
fn body(chart: &EquilibriumEnthalpyPhaseCurve) -> LumpedEnthalpyBody<'_> {
    LumpedEnthalpyBody::try_new(
        "original complete film",
        RHO * WIDTH,
        2.0,
        0.0,
        0.0,
        WIDTH / 2.0,
        100.0,
        chart,
    )
    .unwrap()
}
fn limits() -> ThermalInterferenceLimits {
    ThermalInterferenceLimits {
        maximum_lines: 16,
        maximum_preparation_workspace_bytes: 128 << 20,
        angular: AngularSpectrumLimits {
            maximum_grid_points: 16,
            maximum_workspace_bytes: 8 << 20,
            maximum_fft_work: 1 << 20,
        },
        illumination: AngularSpectrumIlluminationLimits {
            maximum_layer_mode_pairs: 64,
            maximum_workspace_bytes: 8 << 20,
            maximum_fft_work: 1 << 20,
        },
        exposure: AngularSpectrumExposureLimits {
            maximum_lines: 16,
            maximum_point_line_pairs: 64,
            maximum_point_frequency_pairs: 256,
            maximum_workspace_bytes: 8 << 20,
            maximum_fft_work: 1 << 20,
        },
        ensemble: AngularSpectrumExposureEnsembleLimits {
            maximum_factors: 8,
            maximum_point_line_pairs: 64,
            maximum_point_frequency_pairs: 256,
            maximum_workspace_bytes: 8 << 20,
            maximum_fft_work: 1 << 20,
        },
        layer: AngularSpectrumExposureLayerLimits {
            maximum_frequency_groups: 16,
            maximum_frequency_layer_mode_pairs: 256,
            maximum_point_line_pairs: 64,
            maximum_point_reservoir_frequency_pairs: 2048,
            maximum_workspace_bytes: 8 << 20,
            maximum_fft_work: 1 << 20,
        },
        heat: AngularSpectrumIntegratedHeatLimits {
            maximum_depth_samples: 32,
            maximum_point_line_pairs: 4096,
            maximum_fft_work: 1 << 20,
            maximum_spatial_point_frequency_pairs: 4096,
            maximum_reservoir_frequency_pairs: 256,
            maximum_generated_terms: 128,
            maximum_workspace_bytes: 8 << 20,
        },
    }
}
fn depth() -> Vec<AngularSpectrumHeatDepthSample> {
    (0..=8)
        .map(|i| AngularSpectrumHeatDepthSample {
            depth_m: WIDTH * i as f64 / 8.0,
            weight_m: WIDTH / 24.0
                * if i == 0 || i == 8 {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                },
        })
        .collect()
}
fn volume(nodes: &[AngularSpectrumHeatDepthSample]) -> AngularSpectrumHeatVolumeQuadrature<'_> {
    AngularSpectrumHeatVolumeQuadrature {
        depth_span_m: [0.0, WIDTH],
        depth_samples: nodes,
        transverse: AngularSpectrumHeatAreaQuadrature::PeriodicCell,
    }
}
fn lines() -> [ThermalInterferenceLine<'static>; 6] {
    std::array::from_fn(|i| ThermalInterferenceLine {
        angular_frequency_rad_s: if i % 2 == 0 { 2.0 } else { 2.35 },
        left_tangential_e_v_per_m: &LEFT[i],
        right_tangential_e_v_per_m: &RIGHT[i],
    })
}
fn factors<'a>(lines: &'a [ThermalInterferenceLine<'a>]) -> [ThermalInterferenceFactor<'a>; 2] {
    [
        ThermalInterferenceFactor {
            lines: &lines[2..4],
            statistics: AngularSpectrumFactorStatistics::Real,
        },
        ThermalInterferenceFactor {
            lines: &lines[4..6],
            statistics: AngularSpectrumFactorStatistics::ProperComplex,
        },
    ]
}
fn config(initial_h: f64, duration: f64) -> LumpedEnthalpyMarchConfig {
    LumpedEnthalpyMarchConfig {
        initial_specific_enthalpy_j_kg: initial_h,
        ambient_temperature_k: 300.0,
        radiation_temperature_k: 300.0,
        internal_power_w: 0.0,
        duration_s: duration,
        maximum_step_s: 0.01,
        maximum_steps: 1000,
        enthalpy_tolerance_j_kg: 1e-12,
    }
}
fn controls() -> LumpedEnthalpySpectrumControls {
    LumpedEnthalpySpectrumControls {
        start_s: -999.0,
        minimum_step_s: 1e-6,
        source_ripple_tolerance_j: 1e-3,
        maximum_step_energy_residual_j: 1e-10,
        maximum_attempts: 2000,
        maximum_spectral_work: 1 << 20,
        maximum_output_samples: 1001,
        maximum_workspace_bytes: 1 << 20,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

// Independent physical oracle: direct native two-port fields, no angular FFT,
// no spectral-product contraction and no PoleTransfer helper. The original
// Debye ODE supplies Pdot = i*omega*Delta*E/(1+i*omega*tau). Enumerating eight
// coefficient states exactly realizes the mixed real/proper second moments.
fn reference_power(
    medium: &ResolvedMaterialStatePoint,
    dielectric: &ResolvedDebyeCapacitor,
    conductor: &ResolvedConductorState,
    lines: &[ThermalInterferenceLine<'_>],
    nodes: &[AngularSpectrumHeatDepthSample],
    time: f64,
) -> f64 {
    let pole = dielectric.poles()[0];
    let sigma = conductor.resistivity_ohm_m().recip();
    let response = DielectricResponse::DebyeWithOhmic {
        dielectric,
        conductor,
    };
    let stacks: Vec<_> = lines
        .iter()
        .map(|line| {
            PlanarDielectricLayers::try_new(
                medium,
                DielectricResponse::Lossless(medium),
                line.angular_frequency_rad_s,
                &[DielectricLayer {
                    response,
                    thickness_m: WIDTH,
                }],
            )
            .unwrap()
        })
        .collect();
    let mut total = 0.0;
    for node in nodes {
        let fields: Vec<_> = stacks
            .iter()
            .zip(lines)
            .map(|(stack, line)| {
                let mut electric = [C64::ZERO; 3];
                for (axis, polarization) in [(0, Polarization::Tm), (1, Polarization::Te)] {
                    let illumination = stack.prepare_illumination(1.0, polarization).unwrap();
                    let point = illumination.prepare_point(0, node.depth_m).unwrap();
                    let field = point
                        .evaluate([
                            line.left_tangential_e_v_per_m[0][axis],
                            line.right_tangential_e_v_per_m[0][axis],
                        ])
                        .unwrap();
                    for axis in 0..3 {
                        electric[axis] = electric[axis] + field.electric_v_m[axis];
                    }
                }
                electric
            })
            .collect();
        for real in [-1.0, 1.0] {
            for proper in [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)] {
                let mut electric = [0.0; 3];
                let mut current = [0.0; 3];
                for (index, (line, field)) in lines.iter().zip(&fields).enumerate() {
                    let coefficient = if index < 2 {
                        C64::ONE
                    } else if index < 4 {
                        C64::from_re(real)
                    } else {
                        proper
                    };
                    let omega = line.angular_frequency_rad_s;
                    let (s, c) = (omega * time).sin_cos();
                    let temporal = coefficient * C64::new(c, s);
                    let filter = C64::new(0.0, omega * pole.permittivity_increment_f_m)
                        / C64::new(1.0, omega * pole.relaxation_time_s);
                    for axis in 0..3 {
                        electric[axis] += (field[axis] * temporal).re;
                        current[axis] += (filter * field[axis] * temporal).re;
                    }
                }
                total += node.weight_m / 8.0
                    * (0..3)
                        .map(|axis| {
                            sigma * electric[axis] * electric[axis]
                                + pole.relaxation_time_s / pole.permittivity_increment_f_m
                                    * current[axis]
                                    * current[axis]
                        })
                        .sum::<f64>();
            }
        }
    }
    total
}

#[test]
fn g1_refreshed_two_face_mixed_heat_matches_direct_native_fields_and_debye_ode() {
    let exterior_card = card("external medium", 360.0, false);
    let medium = state(&exterior_card, false);
    let card = card("temperature dependent film", 360.0, true);
    let (dielectric, conductor) = sources(&card);
    let chart = chart(&card);
    let body = body(&chart);
    let stack = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&medium),
        2.0,
        &[DielectricLayer {
            response: DielectricResponse::DebyeWithOhmic {
                dielectric: &dielectric,
                conductor: &conductor,
            },
            thickness_m: WIDTH,
        }],
    )
    .unwrap();
    let lines = lines();
    let factors = factors(&lines);
    let drive = ThermalInterferenceDrive {
        source_clock: "original oscillator",
        mean: &lines[..2],
        factors: &factors,
    };
    let nodes = depth();
    let owner = QuasiSteadyFilmHeating::try_new(
        &stack,
        0,
        &card,
        "T",
        &body,
        GRID,
        drive,
        volume(&nodes),
        limits(),
        || false,
    )
    .unwrap();
    let original_identity = dielectric.material().identity();
    let load_300 = owner
        .heat_load(owner.initial_state(Dd::ZERO, 0.0).unwrap(), || false)
        .unwrap();
    let load_320 = owner
        .heat_load(owner.initial_state(Dd::ZERO, 20.0 * CP).unwrap(), || false)
        .unwrap();
    let warm = dielectric.at_query_point(&card, &at(320.0)).unwrap();
    let warm_conductor = conductor.with_material(warm.material().clone()).unwrap();
    assert_eq!(
        load_320.response().material().identity(),
        warm.material().identity()
    );
    assert_eq!(
        load_320.response().material().state_coordinate("pressure"),
        Some(100_000.0)
    );
    close(
        load_320
            .response()
            .material()
            .property(MU_NAME)
            .unwrap()
            .value_si(),
        1.1,
        1e-14,
    );
    for time in [0.0, 0.17, 0.71, 1.4] {
        close(
            load_300.spectrum().evaluate(time).unwrap(),
            reference_power(&medium, &dielectric, &conductor, &lines, &nodes, time),
            2e-11,
        );
        close(
            load_320.spectrum().evaluate(time).unwrap(),
            reference_power(&medium, &warm, &warm_conductor, &lines, &nodes, time),
            2e-11,
        );
    }
    assert!(
        (load_300.spectrum().evaluate(0.17).unwrap() - load_320.spectrum().evaluate(0.17).unwrap())
            .abs()
            > 1e-3
    );
    assert_eq!(dielectric.material().identity(), original_identity);
    close(load_300.represented_volume_m3_dd().to_f64(), WIDTH, 1e-14);
}

#[test]
fn g3_native_heat_feedback_rebuilds_next_interval_and_preserves_late_dd_phase_clock() {
    let exterior_card = card("external medium", 360.0, false);
    let medium = state(&exterior_card, false);
    let card = card("temperature dependent film", 360.0, true);
    let (dielectric, conductor) = sources(&card);
    let chart = chart(&card);
    let body = body(&chart);
    let stack = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&medium),
        2.0,
        &[DielectricLayer {
            response: DielectricResponse::DebyeWithOhmic {
                dielectric: &dielectric,
                conductor: &conductor,
            },
            thickness_m: WIDTH,
        }],
    )
    .unwrap();
    let lines = lines();
    let factors = factors(&lines);
    let drive = ThermalInterferenceDrive {
        source_clock: "original oscillator",
        mean: &lines[..2],
        factors: &factors,
    };
    let nodes = depth();
    let continuous_volume = AngularSpectrumHeatVolumeQuadrature {
        transverse: AngularSpectrumHeatAreaQuadrature::ContinuousPeriodicCell,
        ..volume(&nodes)
    };
    let mut continuous_limits = limits();
    continuous_limits.heat.maximum_fft_work = 0;
    let owner = QuasiSteadyFilmHeating::try_new(
        &stack,
        0,
        &card,
        "T",
        &body,
        GRID,
        drive,
        continuous_volume,
        continuous_limits,
        || false,
    )
    .unwrap();
    let clock = Dd::from_pair(9_007_199_254_740_992.0, 0.25);
    let initial = owner.initial_state(clock, 0.0).unwrap();
    let duration = 0.125;
    let first = owner
        .advance(
            initial,
            BiotGate::corpus_default(),
            config(0.0, duration),
            controls(),
            20.0,
            || false,
        )
        .unwrap();
    let source_j = first
        .applied_load()
        .spectrum()
        .causal_integral_from_clock_dd(clock, duration, 0.0)
        .unwrap()
        .integral_dd()
        .to_f64();
    close(first.energy().source_input_j_dd().to_f64(), source_j, 2e-12);
    close(
        first.endpoint().phase_state().temperature_k(),
        300.0 + source_j / (body.mass_kg() * CP),
        2e-12,
    );
    close(first.energy().residual_j_dd().to_f64(), 0.0, 1e-9);
    assert_eq!(
        (first.endpoint().source_time_s_dd() - clock).to_f64(),
        duration
    );
    assert_eq!(
        first.endpoint_response().material().state_coordinate("T"),
        Some(first.endpoint().phase_state().temperature_k())
    );
    let next = owner
        .advance(
            first.endpoint(),
            BiotGate::corpus_default(),
            config(
                first.endpoint().phase_state().specific_enthalpy_j_kg(),
                duration,
            ),
            controls(),
            20.0,
            || false,
        )
        .unwrap();
    assert_eq!(
        next.applied_load().response().material().identity(),
        first.endpoint_response().material().identity()
    );
    assert_ne!(
        next.applied_load().response().material().identity(),
        first.applied_load().response().material().identity()
    );
    assert_eq!(
        (next.endpoint().source_time_s_dd() - clock).to_f64(),
        2.0 * duration
    );
    assert_eq!(initial.source_time_s_dd(), clock);
    assert_eq!(initial.phase_state().temperature_k(), 300.0);
}

#[test]
fn g4_domain_excursion_volume_budgets_and_cancellation_preserve_original_state() {
    let exterior_card = card("external medium", 360.0, false);
    let medium = state(&exterior_card, false);
    let card = card("bounded original film", 300.001, true);
    let (dielectric, conductor) = sources(&card);
    let chart = chart(&card);
    let body = body(&chart);
    let stack = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&medium),
        2.0,
        &[DielectricLayer {
            response: DielectricResponse::DebyeWithOhmic {
                dielectric: &dielectric,
                conductor: &conductor,
            },
            thickness_m: WIDTH,
        }],
    )
    .unwrap();
    let lines = lines();
    let factors = factors(&lines);
    let drive = ThermalInterferenceDrive {
        source_clock: "original oscillator",
        mean: &lines[..2],
        factors: &factors,
    };
    let nodes = depth();
    let owner = QuasiSteadyFilmHeating::try_new(
        &stack,
        0,
        &card,
        "T",
        &body,
        GRID,
        drive,
        volume(&nodes),
        limits(),
        || false,
    )
    .unwrap();
    let initial = owner.initial_state(Dd::ZERO, 0.0).unwrap();
    let before = owner.heat_load(initial, || false).unwrap();
    assert!(
        owner
            .advance(
                initial,
                BiotGate::corpus_default(),
                config(0.0, 0.1),
                controls(),
                20.0,
                || false
            )
            .is_err()
    );
    assert!(
        owner
            .advance(
                initial,
                BiotGate::corpus_default(),
                config(0.0, 0.0001),
                controls(),
                0.0,
                || false
            )
            .is_err()
    );
    assert!(owner.heat_load(initial, || true).is_err());
    let retry = owner.heat_load(initial, || false).unwrap();
    assert_eq!(
        before.spectrum().evaluate(0.3).unwrap(),
        retry.spectrum().evaluate(0.3).unwrap()
    );
    let small = owner
        .advance(
            initial,
            BiotGate::corpus_default(),
            config(0.0, 1e-6),
            controls(),
            0.001,
            || false,
        )
        .unwrap();
    assert!(small.endpoint().phase_state().temperature_k() > 300.0);
    assert_eq!(initial.phase_state().temperature_k(), 300.0);
    assert_eq!(dielectric.material().state_coordinate("T"), Some(300.0));

    let areas = [1.0];
    let masked = AngularSpectrumHeatVolumeQuadrature {
        transverse: AngularSpectrumHeatAreaQuadrature::CellAreasM2(&areas),
        ..volume(&nodes)
    };
    assert!(
        QuasiSteadyFilmHeating::try_new(
            &stack,
            0,
            &card,
            "T",
            &body,
            GRID,
            drive,
            masked,
            limits(),
            || false,
        )
        .is_err()
    );
    let wrong_mass = LumpedEnthalpyBody::try_new(
        "different mass",
        2.0 * body.mass_kg(),
        2.0,
        0.0,
        0.0,
        WIDTH / 2.0,
        100.0,
        &chart,
    )
    .unwrap();
    assert!(
        QuasiSteadyFilmHeating::try_new(
            &stack,
            0,
            &card,
            "T",
            &wrong_mass,
            GRID,
            drive,
            volume(&nodes),
            limits(),
            || false,
        )
        .is_err()
    );
    let mut bounded = limits();
    bounded.maximum_lines = 5;
    assert!(
        QuasiSteadyFilmHeating::try_new(
            &stack,
            0,
            &card,
            "T",
            &body,
            GRID,
            drive,
            volume(&nodes),
            bounded,
            || false,
        )
        .is_err()
    );
    bounded = limits();
    bounded.maximum_preparation_workspace_bytes = 0;
    let bounded_owner = QuasiSteadyFilmHeating::try_new(
        &stack,
        0,
        &card,
        "T",
        &body,
        GRID,
        drive,
        volume(&nodes),
        bounded,
        || false,
    )
    .unwrap();
    assert!(bounded_owner.heat_load(initial, || false).is_err());
}
