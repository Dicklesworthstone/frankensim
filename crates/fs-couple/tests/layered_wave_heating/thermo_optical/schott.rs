//! G1/G3/G4: reference-only catalog optics on accepted spatial heat states.
//! Authored coefficients test implementation; they do not qualify a glass.

use super::*;
use fs_blake3::ContentHash;
use fs_couple::layered_wave_heating::{
    BroadbandThermoOpticalLimits, BroadbandThermoOpticalWaveHeating, SchottThermoOpticalSource,
    ThermoOpticalSource, ThermoOpticalSpectralLine,
};
use fs_material::{
    state_point::{
        SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES, SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
        SPECTRAL_DIELECTRIC_IOR_PROPERTIES, SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES,
        SpectralDielectricAbsorptionStatePoint, SpectralOpticalStatePoint,
        SpectralTabulatedDielectricStatePoint, resolve_spectral_dielectric_absorption_state_point,
        resolve_spectral_tabulated_dielectric_state_point,
    },
    thermo_optic::{SchottThermoOpticLaw, ThermoOpticCoefficients},
};

const GAP: f64 = 26e-9;
const DEPTH: f64 = 600e-9;
const ALPHA: f64 = 2e5;

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
}

fn reference_card(alpha: f64) -> MaterialCard {
    let base = card("T0-only glass optical spectrum", 1.0, 0.0);
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        claim.validity = ValidityDomain::unconstrained().with_quantity(
            "T",
            quantity(QuantityKind::AbsoluteTemperature),
            280.0,
            360.0,
        );
        claims.insert_claim(claim).unwrap();
    }
    for i in 0..9 {
        for (name, q, value) in [
            (
                SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[i],
                QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                (500.0 + 25.0 * i as f64) * 1e-9,
            ),
            (
                SPECTRAL_DIELECTRIC_IOR_PROPERTIES[i],
                QuantitySpec::dimensional(Dims::NONE),
                1.5,
            ),
            (
                SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i],
                SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
                alpha,
            ),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::with_quantity(name, q),
                    value: PropertyValue::Scalar {
                        value,
                        dims: q.dims(),
                    },
                    validity: ValidityDomain::unconstrained().with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        300.0,
                        300.0,
                    ),
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: "authored reference-only thermal glass".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

fn actual(card: &MaterialCard, t: f64) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(
        card,
        &point(t),
        &[
            ScalarPropertyRequirement::try_new(MU_NAME, MU, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn optical(card: &MaterialCard) -> SpectralOpticalStatePoint {
    let optical = resolve_spectral_tabulated_dielectric_state_point(
        card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let complete = optical
        .resolved()
        .with_additional_requirements(
            card,
            &[ScalarPropertyRequirement::try_new(
                MU_NAME,
                MU,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
    SpectralOpticalStatePoint::TabulatedDielectric(
        SpectralTabulatedDielectricStatePoint::try_from_material(complete).unwrap(),
    )
}

fn law(upper: f64) -> SchottThermoOpticLaw {
    SchottThermoOpticLaw::try_new(
        ThermoOpticCoefficients {
            d: [1e-4, 0.0, 0.0],
            e: [2e-5, 0.0],
            resonance_um: 0.1,
        },
        300.0,
        [280.0, upper],
        [500.0, 700.0],
        ContentHash([81; 32]),
    )
    .unwrap()
}

fn schott_source(card: &MaterialCard, upper: f64) -> SchottThermoOpticalSource {
    SchottThermoOpticalSource::try_with_reference_absorption(
        optical(card),
        law(upper),
        actual(card, 300.0),
        "T",
        GAP,
    )
    .unwrap()
}

fn thermal_body(card: &MaterialCard, count: usize) -> PreparedLayeredThermalBody {
    let lower = point(280.0);
    let upper = point(360.0);
    let conductivity = PropertyKey::new(K_NAME, THERMAL_CONDUCTIVITY_DIMS);
    let capacity = PropertyKey::new(CP_NAME, SPECIFIC_HEAT_CAPACITY_DIMS);
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let mass = RHO * AREA * DEPTH;
    PreparedLayeredThermalBody::try_new(
        &[TransientLayerSourceQuery {
            conduction: LayerStackSourceQuery {
                id: 100,
                card,
                conductivity: &conductivity,
                thickness_m: DEPTH,
                temperature_axis: "T",
                lower: &lower,
                upper: &upper,
                selection,
            },
            capacity: HeatStorageComponentQuery {
                id: 200,
                card,
                capacity: &capacity,
                amount: HeatStorageAmount::MassKg(mass),
                temperature_axis: "T",
                lower: &lower,
                upper: &upper,
                selection,
            },
            mass_kg: mass,
            cell_count: count,
        }],
        &[],
        AREA,
        300.0,
    )
    .unwrap()
}

fn phase(wavelength_m: f64, t: f64) -> f64 {
    1.5 + (1.5 * 1.5 - 1.0) / (2.0 * 1.5)
        * (t - 300.0)
        * (1e-4 + 2e-5 / ((wavelength_m * 1e6).powi(2) - 0.1_f64.powi(2)))
}

#[test]
fn reference_only_optics_survive_hot_cells_without_discarding_actual_source_admission() {
    let card = reference_card(ALPHA);
    let source = schott_source(&card, 350.0);
    assert!(
        source
            .reference_source()
            .at_query_point(&card, &point(330.0))
            .is_err()
    );
    assert_eq!(source.resolved().properties().len(), 1);
    assert_eq!(source.reference_source().resolved().properties().len(), 28);
    let hot = source.at_query_point(&card, &point(330.0)).unwrap();
    assert!(std::ptr::eq(
        source.reference_source(),
        hot.reference_source()
    ));
    assert_eq!(
        hot.reference_source().resolved().state_coordinate("T"),
        Some(300.0)
    );
    assert_eq!(hot.resolved().state_coordinate("T"), Some(330.0));
    let index = hot
        .response()
        .unwrap()
        .complex_index(OPTICAL_OMEGA)
        .unwrap();
    near(index.re, phase(550e-9, 330.0), 3e-15);
    near(
        index.im,
        -ALPHA * 550e-9 / (2.0 * core::f64::consts::TAU),
        2e-17,
    );
    assert_eq!(hot.at_query_point(&card, &point(300.0)).unwrap(), source);
    assert!(
        SchottThermoOpticalSource::try_lossless(
            optical(&card),
            law(350.0),
            actual(&card, 300.0),
            "T",
            GAP,
        )
        .is_err()
    );
    assert!(
        source
            .at_query_point(&card, &point(351.0))
            .unwrap()
            .wave_response(OPTICAL_OMEGA)
            .is_err()
    );
    assert!(source.at_query_point(&card, &point(361.0)).is_err());
    assert!(
        source
            .wave_response(core::f64::consts::TAU * LIGHT_SPEED / 710e-9)
            .is_err()
    );
    let wrong = reference_card(ALPHA * 2.0);
    assert!(source.at_query_point(&wrong, &point(310.0)).is_err());
    let transparent = reference_card(0.0);
    let lossless = SchottThermoOpticalSource::try_lossless(
        optical(&transparent),
        law(350.0),
        actual(&transparent, 330.0),
        "T",
        GAP,
    )
    .unwrap();
    assert!(!lossless.retains_reference_absorption());
    let index = lossless
        .response()
        .unwrap()
        .complex_index(OPTICAL_OMEGA)
        .unwrap();
    near(index.re, phase(550e-9, 330.0), 3e-15);
    assert_eq!(index.im, 0.0);
}

#[test]
fn nonuniform_accepted_temperatures_drive_native_coherent_heat_and_next_step() {
    let card = reference_card(ALPHA);
    let source = schott_source(&card, 350.0);
    let outside = exterior();
    let native = thermal_body(&card, 2);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let initial = initial_at(&bound, &[290.0, 330.0]);
    let before = initial.cells().clone();
    let temperatures: Vec<_> = bound
        .cell_ids()
        .iter()
        .map(|&id| initial.temperature_k(id).unwrap())
        .collect();
    let bindings = [ThermoOpticalLayer {
        source: &source,
        card: &card,
    }];
    let owner = ThermoOpticalWaveHeating::try_new(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        OPTICAL_OMEGA,
        0.8,
        Polarization::Tm,
        &bindings,
        || false,
    )
    .unwrap();
    let refreshed: Vec<_> = temperatures
        .iter()
        .map(|&t| source.at_query_point(&card, &point(t)).unwrap())
        .collect();
    let responses: Vec<_> = refreshed
        .iter()
        .map(|s| s.wave_response(OPTICAL_OMEGA).unwrap())
        .collect();
    let original = LayeredWaveHeating::try_new_cellwise(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        OPTICAL_OMEGA,
        &responses,
        || false,
    )
    .unwrap();
    let drive = LayeredWaveDrive::Coherent([C64::new(2e6, 0.0), C64::new(3e5, 2e5)]);
    let mut powers = [0.0; 2];
    original
        .prepare(0.8, Polarization::Tm, || false)
        .unwrap()
        .cell_power_into(drive, &mut powers, || false)
        .unwrap();
    let frozen = LayeredWaveHeating::try_new(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        OPTICAL_OMEGA,
        &[source.wave_response(OPTICAL_OMEGA).unwrap()],
        || false,
    )
    .unwrap();
    let mut frozen_powers = [0.0; 2];
    frozen
        .prepare(0.8, Polarization::Tm, || false)
        .unwrap()
        .cell_power_into(drive, &mut frozen_powers, || false)
        .unwrap();
    assert!(
        powers
            .iter()
            .zip(frozen_powers)
            .any(|(&p, f)| (p - f).abs() > 1e-5 * p)
    );
    let dt = 1e-10;
    let controls = LayeredThermalSolveControls {
        cell_tolerance_j: 1e-11,
        total_tolerance_j: 1e-11,
        ..LayeredThermalSolveControls::default()
    };
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    let mut polls = 0;
    let result = owner
        .advance_to(
            &initial,
            dt,
            drive,
            conduction(),
            controls,
            &mut workspace,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(result.applied_optical_states, refreshed);
    // Independent two-cell constant-capacity backward Euler solution. Optical
    // work is the original field owner's local positive volume dissipation.
    let capacity = RHO * AREA * (DEPTH / 2.0) * CP;
    let g_dt = AREA * K / (DEPTH / 2.0) * dt;
    let mean =
        0.5 * (temperatures[0] + temperatures[1]) + 0.5 * dt * (powers[0] + powers[1]) / capacity;
    let contrast = (0.5 * capacity * (temperatures[0] - temperatures[1])
        + 0.5 * dt * (powers[0] - powers[1]))
        / (capacity + 2.0 * g_dt);
    for (i, expected) in [mean + contrast, mean - contrast].into_iter().enumerate() {
        let t = result
            .step
            .thermal
            .temperature_k(bound.cell_ids()[i])
            .unwrap();
        near(t, expected, 3e-9);
        near(result.step.cell_energy_j[i].1, powers[i] * dt, 1e-12);
        assert_eq!(
            result.endpoint_optical_states[i],
            source.at_query_point(&card, &point(t)).unwrap()
        );
    }
    near(
        result.step.heat_report.energy.stored_energy_change_j,
        result.step.optical.absorbed_j,
        3e-11,
    );
    near(result.step.optical.residual_j, 0.0, 3e-11);
    let next = owner
        .advance_to(
            &result.step.thermal,
            2.0 * dt,
            drive,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(next.applied_optical_states, result.endpoint_optical_states);
    assert_ne!(next.applied_optical_states, result.applied_optical_states);
    let mut cancelled_polls = 0;
    assert!(matches!(
        owner.advance_to(
            &initial,
            dt,
            drive,
            conduction(),
            controls,
            &mut workspace,
            || {
                cancelled_polls += 1;
                cancelled_polls == polls
            }
        ),
        Err(LayeredWaveHeatingError::Cancelled)
    ));
    assert_eq!(initial.cells(), &before);
    assert_eq!(initial.time_s(), 0.0);
}

// Independent equal-exterior complex Airy series, using standard trig/exp.
fn airy_absorption(wavelength: f64, t: f64) -> f64 {
    airy_absorption_with_alpha(wavelength, t, ALPHA)
}

fn airy_absorption_with_alpha(wavelength: f64, t: f64, alpha: f64) -> f64 {
    let n = C64::new(
        phase(wavelength, t),
        -alpha * wavelength / (2.0 * core::f64::consts::TAU),
    );
    let r = (C64::ONE - n) / (C64::ONE + n);
    let angle = -core::f64::consts::TAU * n.re * DEPTH / wavelength;
    let p = C64::new(angle.cos(), angle.sin()).scale((-alpha * DEPTH / 2.0).exp());
    let denominator = C64::ONE - r * r * p * p;
    let reflected = r * (C64::ONE - p * p) / denominator;
    let transmitted = (C64::ONE - r * r) * p / denominator;
    1.0 - reflected.norm_sq() - transmitted.norm_sq()
}

#[test]
fn broadband_thermal_glass_matches_independent_airy_heat_and_refuses_endpoint_domain_exit() {
    let card = reference_card(ALPHA);
    let source = schott_source(&card, 350.0);
    let outside = exterior();
    let native = thermal_body(&card, 1);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let initial = initial_at(&bound, &[320.0]);
    let t0 = initial.temperature_k(bound.cell_ids()[0]).unwrap();
    let wavelengths = [550e-9, 650e-9];
    let lines = wavelengths.map(|w| ThermoOpticalSpectralLine {
        angular_frequency_rad_s: core::f64::consts::TAU * LIGHT_SPEED / w,
        incidence_cosine: 1.0,
        polarization: Polarization::Te,
    });
    let limits = BroadbandThermoOpticalLimits {
        maximum_lines: 2,
        maximum_cell_line_pairs: 2,
        maximum_cell_factor_pairs: 2,
    };
    let bindings = [ThermoOpticalLayer {
        source: &source,
        card: &card,
    }];
    let owner = BroadbandThermoOpticalWaveHeating::try_new(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        &lines,
        &bindings,
        limits,
        || false,
    )
    .unwrap();
    let amplitudes = [1.0e6, 1.4e6];
    let drives = amplitudes.map(|a| LayeredWaveDrive::Coherent([C64::new(a, 0.0), C64::ZERO]));
    let dt = 1e-10;
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    let controls = LayeredThermalSolveControls::default();
    let result = owner
        .advance_to(
            &initial,
            dt,
            &drives,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    let mut expected = 0.0;
    for (i, &w) in wavelengths.iter().enumerate() {
        let flux = amplitudes[i].powi(2) / (2.0 * M * LIGHT_SPEED);
        let absorption = flux * airy_absorption(w, t0);
        near(
            result.line_power[i].absorbed_w_m2,
            absorption,
            2e-11 * absorption,
        );
        expected += AREA * absorption * dt;
    }
    near(result.step.optical.absorbed_j, expected, 2e-11 * expected);
    near(
        result
            .step
            .thermal
            .temperature_k(bound.cell_ids()[0])
            .unwrap(),
        t0 + expected / (RHO * AREA * DEPTH * CP),
        3e-9,
    );
    // Thermal support admits this endpoint; the independent SCHOTT bound
    // refuses it before the joint optical/thermal result can be published.
    let narrow = schott_source(&card, t0 + expected / (2.0 * RHO * AREA * DEPTH * CP));
    let narrow_bindings = [ThermoOpticalLayer {
        source: &narrow,
        card: &card,
    }];
    let narrow_owner = BroadbandThermoOpticalWaveHeating::try_new(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        &lines,
        &narrow_bindings,
        limits,
        || false,
    )
    .unwrap();
    let before = initial.cells().clone();
    assert!(matches!(
        narrow_owner.advance_to(
            &initial,
            dt,
            &drives,
            conduction(),
            controls,
            &mut workspace,
            || false
        ),
        Err(LayeredWaveHeatingError::Electromagnetic(_))
    ));
    assert_eq!(initial.cells(), &before);
    assert_eq!(initial.time_s(), 0.0);
}

fn current_alpha_card(upper: f64) -> MaterialCard {
    let base = reference_card(ALPHA);
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        let name = claim.key.name();
        if SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES.contains(&name) {
            claim.validity = ValidityDomain::unconstrained().with_quantity(
                "T",
                quantity(QuantityKind::AbsoluteTemperature),
                280.0,
                upper,
            );
        } else if SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES.contains(&name) {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: vec![(280.0, 0.5 * ALPHA), (360.0, 2.5 * ALPHA)],
                dims: claim.key.dims(),
            };
            claim.validity = ValidityDomain::unconstrained().with_quantity(
                "T",
                quantity(QuantityKind::AbsoluteTemperature),
                280.0,
                upper,
            );
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

fn current_alpha_source(card: &MaterialCard) -> SchottThermoOpticalSource {
    let source = resolve_spectral_dielectric_absorption_state_point(
        card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let complete = source
        .resolved()
        .with_additional_requirements(
            card,
            &[ScalarPropertyRequirement::try_new(
                MU_NAME,
                MU,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
    SchottThermoOpticalSource::try_with_current_absorption(
        optical(card),
        law(350.0),
        SpectralDielectricAbsorptionStatePoint::try_from_material(complete).unwrap(),
        "T",
        GAP,
        GAP,
    )
    .unwrap()
}

fn alpha_at(t: f64) -> f64 {
    ALPHA * (0.5 + (t - 280.0) / 40.0)
}

#[test]
fn current_alpha_refresh_keeps_reference_phase_and_original_actual_requirements_separate() {
    let card = current_alpha_card(360.0);
    let source = current_alpha_source(&card);
    let hot = source.at_query_point(&card, &point(330.0)).unwrap();
    assert!(
        source
            .reference_source()
            .at_query_point(&card, &point(330.0))
            .is_err()
    );
    assert!(std::ptr::eq(source.reference_source(), hot.reference_source()));
    assert!(!hot.retains_reference_absorption());
    assert_eq!(hot.resolved().properties().len(), 19);
    assert_eq!(hot.current_absorption().unwrap().resolved(), hot.resolved());
    assert_eq!(
        hot.reference_source().resolved().state_coordinate("T"),
        Some(300.0)
    );
    assert_eq!(hot.resolved().state_coordinate("T"), Some(330.0));
    let index = hot
        .response()
        .unwrap()
        .complex_index(OPTICAL_OMEGA)
        .unwrap();
    near(index.re, phase(550e-9, 330.0), 3e-15);
    near(
        index.im,
        -alpha_at(330.0) * 550e-9 / (4.0 * core::f64::consts::PI),
        3e-17,
    );
    assert_eq!(hot.at_query_point(&card, &point(300.0)).unwrap(), source);
    assert!(source.at_query_point(&card, &point(361.0)).is_err());
    let wrong = current_alpha_card(340.0);
    assert!(source.at_query_point(&wrong, &point(330.0)).is_err());
}

#[test]
fn current_alpha_drives_native_mono_and_broadband_heat_with_atomic_source_bounds() {
    for upper in [360.0, 320.2] {
        let card = current_alpha_card(upper);
        let source = current_alpha_source(&card);
        let outside = exterior();
        let native = thermal_body(&card, 1);
        let laws = native.thermal_laws(false).unwrap();
        let bound = native.bind(&laws, ids()).unwrap();
        let initial = initial_at(&bound, &[320.0]);
        let before = initial.cells().clone();
        let t0 = initial.temperature_k(bound.cell_ids()[0]).unwrap();
        let wavelengths = [550e-9, 650e-9];
        let lines = wavelengths.map(|w| ThermoOpticalSpectralLine {
            angular_frequency_rad_s: core::f64::consts::TAU * LIGHT_SPEED / w,
            incidence_cosine: 1.0,
            polarization: Polarization::Te,
        });
        let limits = BroadbandThermoOpticalLimits {
            maximum_lines: 2,
            maximum_cell_line_pairs: 2,
            maximum_cell_factor_pairs: 2,
        };
        let bindings = [ThermoOpticalLayer {
            source: &source,
            card: &card,
        }];
        let owner = BroadbandThermoOpticalWaveHeating::try_new(
            &bound,
            &outside,
            DielectricResponse::Lossless(&outside),
            &lines,
            &bindings,
            limits,
            || false,
        )
        .unwrap();
        let amplitudes: [f64; 2] = [1e6, 1.4e6];
        let drives = amplitudes
            .map(|a| LayeredWaveDrive::Coherent([C64::new(a, 0.0), C64::ZERO]));
        let dt = 1e-10;
        let mut powers = [0.0; 2];
        let mut frozen_power = 0.0;
        for (i, &w) in wavelengths.iter().enumerate() {
            let flux = amplitudes[i].powi(2) / (2.0 * M * LIGHT_SPEED);
            powers[i] = AREA * flux * airy_absorption_with_alpha(w, t0, alpha_at(t0));
            frozen_power += AREA * flux * airy_absorption(w, t0);
        }
        let expected_energy = dt * powers.iter().sum::<f64>();
        let expected_t = t0 + expected_energy / (RHO * AREA * DEPTH * CP);
        assert!(expected_energy > 1.4 * dt * frozen_power);
        let controls = LayeredThermalSolveControls::default();
        let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
        if upper < expected_t {
            assert!(
                source
                    .at_query_point(&card, &point(t0))
                    .unwrap()
                    .wave_response(lines[0].angular_frequency_rad_s)
                    .is_ok()
            );
            assert!(matches!(
                owner.advance_to(
                    &initial,
                    dt,
                    &drives,
                    conduction(),
                    controls,
                    &mut workspace,
                    || false,
                ),
                Err(LayeredWaveHeatingError::Material(_))
            ));
        } else {
            let result = owner
                .advance_to(
                    &initial,
                    dt,
                    &drives,
                    conduction(),
                    controls,
                    &mut workspace,
                    || false,
                )
                .unwrap();
            near(
                result.step.optical.absorbed_j,
                expected_energy,
                2e-11 * expected_energy,
            );
            near(
                result
                    .step
                    .thermal
                    .temperature_k(bound.cell_ids()[0])
                    .unwrap(),
                expected_t,
                3e-9,
            );
            near(
                result.step.heat_report.energy.stored_energy_change_j,
                expected_energy,
                3e-11,
            );
            let mono = ThermoOpticalWaveHeating::try_new(
                &bound,
                &outside,
                DielectricResponse::Lossless(&outside),
                lines[0].angular_frequency_rad_s,
                1.0,
                Polarization::Te,
                &bindings,
                || false,
            )
            .unwrap();
            let one_line = mono
                .advance_to(
                    &initial,
                    dt,
                    drives[0],
                    conduction(),
                    controls,
                    &mut workspace,
                    || false,
                )
                .unwrap();
            near(
                one_line.step.optical.absorbed_j,
                dt * powers[0],
                2e-11 * dt * powers[0],
            );
            assert_eq!(one_line.applied_optical_states, result.applied_optical_states);
            let next = owner
                .advance_to(
                    &result.step.thermal,
                    2.0 * dt,
                    &drives,
                    conduction(),
                    controls,
                    &mut workspace,
                    || false,
                )
                .unwrap();
            assert_eq!(next.applied_optical_states, result.endpoint_optical_states);
            assert_ne!(next.applied_optical_states, result.applied_optical_states);
            let accepted_alpha = next.applied_optical_states[0]
                .current_absorption()
                .unwrap()
                .absorption_per_m_at(550e-9)
                .unwrap();
            near(accepted_alpha, alpha_at(expected_t), 1e-5);
        }
        assert_eq!(initial.cells(), &before);
        assert_eq!(initial.time_s(), 0.0);
    }
}
