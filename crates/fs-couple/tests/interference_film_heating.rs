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
        dielectric::{
            DebyePoleProperties, ResolvedDebyeCapacitor,
            lorentz::{LorentzPoleProperties, ResolvedLorentzMaterial},
        },
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
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
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

/// Every f64 a refresh must reproduce: the held state coordinates and each
/// resolved property value.
fn physical_bits(state: &ResolvedMaterialStatePoint) -> Vec<u64> {
    let mut values: Vec<f64> = state
        .query_point()
        .iter()
        .map(|(_, value)| *value)
        .collect();
    values.extend(state.properties().iter().map(|p| p.value_si()));
    values.into_iter().map(f64::to_bits).collect()
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
    carrier: [f64; 2],
    time: f64,
) -> f64 {
    let pole = dielectric.poles()[0];
    let sigma = conductor.resistivity_ohm_m().recip();
    let response = DielectricResponse::DebyeWithOhmic {
        dielectric,
        conductor,
    };
    reference_reservoir_power(
        medium,
        response,
        sigma,
        lines,
        nodes,
        carrier,
        time,
        |omega| {
            C64::new(0.0, omega * pole.permittivity_increment_f_m)
                / C64::new(1.0, omega * pole.relaxation_time_s)
        },
        pole.relaxation_time_s / pole.permittivity_increment_f_m,
    )
}

// Shared independent field/ensemble integral. The caller supplies Pdot/E and
// the physical positive heat coefficient from its original oscillator ODE.
#[allow(clippy::too_many_arguments)]
fn reference_reservoir_power(
    medium: &ResolvedMaterialStatePoint,
    response: DielectricResponse<'_>,
    sigma: f64,
    lines: &[ThermalInterferenceLine<'_>],
    nodes: &[AngularSpectrumHeatDepthSample],
    carrier: [f64; 2],
    time: f64,
    transfer: impl Fn(f64) -> C64,
    dissipation: f64,
) -> f64 {
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
    let kt = carrier[0].hypot(carrier[1]);
    let p = if kt == 0.0 {
        [1.0, 0.0]
    } else {
        carrier.map(|k| k / kt)
    };
    // Integral of exp(-2 i k_carrier.x) over the authored unit-square cell.
    // For a common carrier, integral Re(z exp(-i k.x))² equals
    // (|z|² + Re(z² integral exp(-2 i k.x)))/2, including optical sum tones.
    let spatial = carrier.into_iter().fold(C64::ONE, |value, k| {
        if k == 0.0 {
            value
        } else {
            value * C64::new(k.cos(), -k.sin()).scale(k.sin() / k)
        }
    });
    let mean_square = |value: C64| {
        if kt == 0.0 {
            value.re * value.re
        } else {
            0.5 * (value.norm_sq() + (value * value * spatial).re)
        }
    };
    let mut total = 0.0;
    for node in nodes {
        let fields: Vec<_> = stacks
            .iter()
            .zip(lines)
            .map(|(stack, line)| {
                let mut electric = [C64::ZERO; 3];
                let k = line.angular_frequency_rad_s
                    * (stack.incident_permittivity_f_m() * stack.incident_permeability_h_m())
                        .sqrt();
                let cosine = (1.0 - (kt / k).powi(2)).sqrt();
                for polarization in [Polarization::Tm, Polarization::Te] {
                    let component = |e: [C64; 2]| match polarization {
                        Polarization::Tm => e[0].scale(p[0]) + e[1].scale(p[1]),
                        Polarization::Te => e[1].scale(p[0]) - e[0].scale(p[1]),
                    };
                    let illumination = stack.prepare_illumination(cosine, polarization).unwrap();
                    let point = illumination.prepare_point(0, node.depth_m).unwrap();
                    let field = point
                        .evaluate([
                            component(line.left_tangential_e_v_per_m[0]),
                            component(line.right_tangential_e_v_per_m[0]),
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
                let mut electric = [C64::ZERO; 3];
                let mut current = [C64::ZERO; 3];
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
                    let filter = transfer(omega);
                    for axis in 0..3 {
                        electric[axis] = electric[axis] + field[axis] * temporal;
                        current[axis] = current[axis] + filter * field[axis] * temporal;
                    }
                }
                total += node.weight_m / 8.0
                    * (0..3)
                        .map(|axis| {
                            sigma * mean_square(electric[axis])
                                + dissipation * mean_square(current[axis])
                        })
                        .sum::<f64>();
            }
        }
    }
    total
}

fn lorentz_card() -> MaterialCard {
    let base = card("temperature-dependent original Lorentz film", 360.0, true);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        claims.insert_claim(claim.clone()).unwrap();
    }
    for (name, q, values) in [
        (
            "omega0",
            quantity(QuantityKind::Frequency(FrequencyConvention::Angular)),
            [2.8, 2.2],
        ),
        (
            "gamma",
            QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
            [0.16, 0.48],
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
                        360.0,
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
                    source: "authored temperature-dependent inertial polarization".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

#[test]
fn g1_hot_lorentz_refraction_and_mixed_heat_follow_the_original_reservoir_ode() {
    let medium = state(&card("Lorentz feedback exterior", 360.0, false), false);
    let card = lorentz_card();
    let material = state(&card, true)
        .with_additional_requirements(
            &card,
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
    let dielectric = ResolvedLorentzMaterial::try_from_material(
        material.clone(),
        EPS_NAME,
        &[LorentzPoleProperties {
            permittivity_increment: "delta".into(),
            angular_frequency: "omega0".into(),
            damping_rate: "gamma".into(),
        }],
    )
    .unwrap();
    let conductor = ResolvedConductorState::try_from_material_with_basis(
        material.clone(),
        ElectricalPropertyBasis::Conductivity,
    )
    .unwrap();
    let warm = ResolvedLorentzMaterial::try_from_material(
        material.at_query_point(&card, &at(320.0)).unwrap(),
        dielectric.high_frequency_property(),
        dielectric.pole_properties(),
    )
    .unwrap();
    let warm_conductor = conductor.with_material(warm.material().clone()).unwrap();
    let chart = chart(&card);
    let body = body(&chart);
    let lines = lines();
    let factors = factors(&lines);
    let nodes = depth();
    let carrier = [0.7, -0.45];
    for leaky in [false, true] {
        let response = if leaky {
            DielectricResponse::LorentzWithOhmic {
                dielectric: &dielectric,
                conductor: &conductor,
            }
        } else {
            DielectricResponse::Lorentz(&dielectric)
        };
        let warm_response = if leaky {
            DielectricResponse::LorentzWithOhmic {
                dielectric: &warm,
                conductor: &warm_conductor,
            }
        } else {
            DielectricResponse::Lorentz(&warm)
        };
        let stack = PlanarDielectricLayers::try_new(
            &medium,
            DielectricResponse::Lossless(&medium),
            2.0,
            &[DielectricLayer {
                response,
                thickness_m: WIDTH,
            }],
        )
        .unwrap();
        let owner = QuasiSteadyFilmHeating::try_new_with_carrier(
            &stack,
            0,
            &card,
            "T",
            &body,
            GRID,
            carrier,
            ThermalInterferenceDrive {
                source_clock: "original Lorentz oscillator",
                mean: &lines[..2],
                factors: &factors,
            },
            AngularSpectrumHeatVolumeQuadrature {
                transverse: AngularSpectrumHeatAreaQuadrature::ContinuousPeriodicCell,
                ..volume(&nodes)
            },
            limits(),
            || false,
        )
        .unwrap();
        let cold = owner
            .heat_load(owner.initial_state(Dd::ZERO, 0.0).unwrap(), || false)
            .unwrap();
        let hot = owner
            .heat_load(owner.initial_state(Dd::ZERO, 20.0 * CP).unwrap(), || false)
            .unwrap();
        assert_eq!(hot.response().material(), warm.material());
        // A heat load refreshes the film through at_query_point, which pins
        // each selected claim; `material` was resolved with SingleClaimOnly.
        // "A refresh records pinned-claim selection; an independent
        // resolution records SingleClaimOnly. Their physical tuples, not
        // policy receipts, must agree exactly" (fs-material
        // src/liquid/density/tests.rs:298-300). So the cold load equals the
        // original refreshed through the same path at its own point, whole,
        // and every physical value equals the original's bit for bit.
        assert_eq!(
            cold.response().material(),
            &material.at_query_point(&card, &material.query()).unwrap()
        );
        assert_eq!(
            physical_bits(cold.response().material()),
            physical_bits(&material)
        );
        if let DielectricResponse::LorentzWithOhmic { dielectric, conductor } = hot.response() {
            assert_eq!(dielectric.material(), conductor.material());
        }
        // Independent elimination of Pddot + gamma Pdot + omega0² P
        // = Delta_epsilon omega0² E. Damping heat is gamma Pdot²/(Delta omega0²).
        // Inertial/restoring storage and net electromagnetic work are not heat.
        let omega0: f64 = 2.5;
        let gamma: f64 = 0.32;
        let delta: f64 = 0.8;
        let sigma = if leaky { 0.6 } else { 0.0 };
        let susceptibility = |omega: f64| {
            C64::from_re(delta * omega0.powi(2))
                / C64::new(omega0.powi(2) - omega.powi(2), gamma * omega)
        };
        for omega in [2.0, 2.35] {
            let expected = C64::new(2.2, -sigma / omega) + susceptibility(omega);
            let actual = hot.response().permittivity_f_m(omega).unwrap();
            close(actual.re, expected.re, 2e-14);
            close(actual.im, expected.im, 2e-14);
        }
        for time in [0.0, 0.17, 0.71, 1.4] {
            close(
                hot.spectrum().evaluate(time).unwrap(),
                reference_reservoir_power(
                    &medium,
                    warm_response,
                    sigma,
                    &lines,
                    &nodes,
                    carrier,
                    time,
                    |omega| C64::new(0.0, omega) * susceptibility(omega),
                    gamma / (delta * omega0.powi(2)),
                ),
                3e-11,
            );
        }
        assert!(
            (hot.spectrum().evaluate(0.17).unwrap() - cold.spectrum().evaluate(0.17).unwrap())
                .abs()
                > 1e-3
        );
        assert_eq!(dielectric.material(), &material);
    }
}

#[test]
fn g1_refreshed_two_face_mixed_heat_matches_direct_native_fields_and_debye_ode() {
    check_refreshed_carrier_heat([0.0; 2]);
}

#[test]
fn off_axis_two_frequency_heat_retains_physical_carrier_and_continuous_cell_interference() {
    check_refreshed_carrier_heat([1.1, 0.6]);
}

fn check_refreshed_carrier_heat(carrier: [f64; 2]) {
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
    let quadrature = AngularSpectrumHeatVolumeQuadrature {
        transverse: AngularSpectrumHeatAreaQuadrature::ContinuousPeriodicCell,
        ..volume(&nodes)
    };
    let owner = QuasiSteadyFilmHeating::try_new_with_carrier(
        &stack,
        0,
        &card,
        "T",
        &body,
        GRID,
        carrier,
        drive,
        quadrature,
        limits(),
        || false,
    )
    .unwrap();
    assert_eq!(owner.carrier_transverse_wavenumber_rad_m(), carrier);
    let legacy = QuasiSteadyFilmHeating::try_new(
        &stack,
        0,
        &card,
        "T",
        &body,
        GRID,
        drive,
        quadrature,
        limits(),
        || false,
    )
    .unwrap();
    assert_eq!(legacy.carrier_transverse_wavenumber_rad_m(), [0.0; 2]);
    let normal = legacy
        .heat_load(legacy.initial_state(Dd::ZERO, 0.0).unwrap(), || false)
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
            reference_power(
                &medium,
                &dielectric,
                &conductor,
                &lines,
                &nodes,
                carrier,
                time,
            ),
            2e-11,
        );
        close(
            load_320.spectrum().evaluate(time).unwrap(),
            reference_power(
                &medium,
                &warm,
                &warm_conductor,
                &lines,
                &nodes,
                carrier,
                time,
            ),
            2e-11,
        );
        if carrier == [0.0; 2] {
            assert_eq!(
                load_300.spectrum().evaluate(time).unwrap(),
                normal.spectrum().evaluate(time).unwrap()
            );
        }
    }
    if carrier != [0.0; 2] {
        assert!(
            (load_300.spectrum().evaluate(0.17).unwrap()
                - normal.spectrum().evaluate(0.17).unwrap())
            .abs()
                > 1e-2
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
fn g3_off_axis_heat_feedback_rebuilds_next_interval_and_preserves_late_dd_phase_clock() {
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
    let carrier = [0.7, -0.45];
    let owner = QuasiSteadyFilmHeating::try_new_with_carrier(
        &stack,
        0,
        &card,
        "T",
        &body,
        GRID,
        carrier,
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

    for carrier in [[f64::NAN, 0.0], [0.0, f64::INFINITY]] {
        assert!(
            QuasiSteadyFilmHeating::try_new_with_carrier(
                &stack,
                0,
                &card,
                "T",
                &body,
                GRID,
                carrier,
                drive,
                volume(&nodes),
                limits(),
                || false,
            )
            .is_err()
        );
    }
    // A finite carrier still has to propagate at EVERY source frequency.
    // It is grazing at omega=2 while the higher-frequency line propagates;
    // it must not be silently dropped or reinterpreted as a fixed angle.
    let cutoff = QuasiSteadyFilmHeating::try_new_with_carrier(
        &stack,
        0,
        &card,
        "T",
        &body,
        GRID,
        [2.0, 0.0],
        drive,
        volume(&nodes),
        limits(),
        || false,
    )
    .unwrap();
    assert!(
        cutoff
            .heat_load(cutoff.initial_state(Dd::ZERO, 0.0).unwrap(), || false)
            .is_err()
    );

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
