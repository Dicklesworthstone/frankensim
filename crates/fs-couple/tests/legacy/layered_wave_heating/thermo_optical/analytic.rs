//! Original analytic phase and authored current alpha(T) in native cell heat.
//! The authored absorption/density/permeability curves are numerical controls.
use super::*;
use fs_couple::layered_wave_heating::{
    AnalyticThermoOpticalSource, BroadbandThermoOpticalLimits, BroadbandThermoOpticalWaveHeating,
    ThermoOpticalSource, ThermoOpticalSpectralLine,
};
use fs_material::{
    silica_refraction::Corning7980Silica,
    state_point::{
        DENSITY_PROPERTY, SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
        SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY, SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES,
        SpectralDielectricAbsorptionStatePoint, resolve_spectral_dielectric_absorption_state_point,
    },
    water_refraction::IapwsWaterRefraction,
};
use fs_qty::{Density, Temperature};

const DEPTH: f64 = 600e-9;
const GAP: f64 = 26e-9;

fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), t).unwrap()
}
fn alpha(t: f64, wavelength_m: f64) -> f64 {
    (100_000.0 + 5000.0 * (t - 280.0)) * (0.8 + 0.002 * (wavelength_m * 1e9 - 500.0))
}
fn rho(t: f64) -> f64 {
    999.0 - 0.3 * (t - 280.0)
}
fn mu(t: f64) -> f64 {
    M * (1.0 + 0.001 * (t - 280.0))
}
fn source_card(alpha_upper: f64) -> MaterialCard {
    let base = card("authored current loss with original analytic water/silica phase", 1.0, 0.0);
    let validity = ValidityDomain::unconstrained().with_quantity(
        "T", quantity(QuantityKind::AbsoluteTemperature), 280.0, 360.0,
    );
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        claim.validity = validity.clone();
        if claim.key.name() == MU_NAME {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, mu(280.0)), (360.0, mu(360.0))], dims: MU,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    let mut extra = Vec::new();
    extra.push((DENSITY_PROPERTY, QuantitySpec::dimensional(Density::DIMS), [rho(280.0), rho(360.0)], false));
    for i in 0..9 {
        let wavelength = (500.0 + 25.0 * i as f64) * 1e-9;
        extra.push((SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[i],
            QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])), [wavelength; 2], false));
        extra.push((SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i],
            SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
            [alpha(280.0, wavelength), alpha(360.0, wavelength)], true));
    }
    for (name, q, values, bounded_alpha) in extra {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(name, q),
            value: PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, values[0]), (360.0, values[1])], dims: q.dims(),
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T", quantity(QuantityKind::AbsoluteTemperature), 280.0,
                if bounded_alpha { alpha_upper } else { 360.0 },
            ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "authored original-alpha/analytic-phase thermal control".into(),
                license: "test-only".into(), artifact: None,
            },
        }).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}
fn source(card: &MaterialCard, water: bool) -> AnalyticThermoOpticalSource {
    let selected = resolve_spectral_dielectric_absorption_state_point(
        card, &point(290.0), MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let complete = selected.resolved().with_additional_requirements(
        card,
        &[
            ScalarPropertyRequirement::try_new(MU_NAME, MU, ScalarAdmissibility::StrictlyPositive).unwrap(),
            ScalarPropertyRequirement::try_new(DENSITY_PROPERTY, Density::DIMS, ScalarAdmissibility::StrictlyPositive).unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let selected = SpectralDielectricAbsorptionStatePoint::try_from_material(complete).unwrap();
    if water {
        AnalyticThermoOpticalSource::try_iapws_water(selected, "T", GAP).unwrap()
    } else {
        AnalyticThermoOpticalSource::try_corning_7980_silica(selected, "T", GAP).unwrap()
    }
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
                id: 100, card, conductivity: &conductivity, thickness_m: DEPTH,
                temperature_axis: "T", lower: &lower, upper: &upper, selection,
            },
            capacity: HeatStorageComponentQuery {
                id: 200, card, capacity: &capacity, amount: HeatStorageAmount::MassKg(mass),
                temperature_axis: "T", lower: &lower, upper: &upper, selection,
            },
            mass_kg: mass, cell_count: count,
        }], &[], AREA, 300.0,
    ).unwrap()
}

#[test]
fn nonuniform_native_cells_refresh_original_alpha_density_mu_and_feed_the_next_step() {
    let card = source_card(360.0);
    let outside = exterior();
    for water in [true, false] {
        let source = source(&card, water);
        let native = thermal_body(&card, 2);
        let laws = native.thermal_laws(false).unwrap();
        let bound = native.bind(&laws, ids()).unwrap();
        let initial = initial_at(&bound, &[285.0, 295.0]);
        let before = initial.cells().clone();
        let temperatures = bound.cell_ids().iter().map(|&id| initial.temperature_k(id).unwrap()).collect::<Vec<_>>();
        let bindings = [ThermoOpticalLayer { source: &source, card: &card }];
        let owner = ThermoOpticalWaveHeating::try_new(
            &bound, &outside, DielectricResponse::Lossless(&outside),
            OPTICAL_OMEGA, 0.8, Polarization::Tm, &bindings, || false,
        ).unwrap();
        let refreshed: Vec<_> = temperatures.iter()
            .map(|&t| source.at_query_point(&card, &point(t)).unwrap()).collect();
        for (state, &temperature) in refreshed.iter().zip(&temperatures) {
            assert_eq!(state.resolved().card_identity(), source.resolved().card_identity());
            assert_eq!(state.resolved().properties().len(), 20);
            near(state.source().absorption_per_m_at(550e-9).unwrap(), alpha(temperature, 550e-9), 1e-8);
            near(state.resolved().property(DENSITY_PROPERTY).unwrap().value_si(), rho(temperature), 3e-12);
            near(state.resolved().property(MU_NAME).unwrap().value_si(), mu(temperature), 1e-20);
        }
        let responses: Vec<_> = refreshed.iter().map(|s| s.wave_response(OPTICAL_OMEGA).unwrap()).collect();
        let original = LayeredWaveHeating::try_new_cellwise(
            &bound, &outside, DielectricResponse::Lossless(&outside),
            OPTICAL_OMEGA, &responses, || false,
        ).unwrap();
        let drive = LayeredWaveDrive::Coherent([C64::new(1e6, 0.0), C64::new(2e5, 1e5)]);
        let mut powers = [0.0; 2];
        original.prepare(0.8, Polarization::Tm, || false).unwrap()
            .cell_power_into(drive, &mut powers, || false).unwrap();
        let frozen = LayeredWaveHeating::try_new(
            &bound, &outside, DielectricResponse::Lossless(&outside), OPTICAL_OMEGA,
            &[source.wave_response(OPTICAL_OMEGA).unwrap()], || false,
        ).unwrap();
        let mut frozen_powers = [0.0; 2];
        frozen.prepare(0.8, Polarization::Tm, || false).unwrap()
            .cell_power_into(drive, &mut frozen_powers, || false).unwrap();
        assert!(powers.iter().zip(frozen_powers).any(|(&p, f)| (p - f).abs() > 0.01 * p));
        let dt = 1e-10;
        let controls = LayeredThermalSolveControls {
            cell_tolerance_j: 1e-11, total_tolerance_j: 1e-11,
            ..LayeredThermalSolveControls::default()
        };
        let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
        let mut polls = 0;
        let result = owner.advance_to(
            &initial, dt, drive, conduction(), controls, &mut workspace,
            || { polls += 1; false },
        ).unwrap();
        assert_eq!(result.applied_optical_states, refreshed);
        let capacity = RHO * AREA * (DEPTH / 2.0) * CP;
        let g_dt = AREA * K / (DEPTH / 2.0) * dt;
        let mean = 0.5 * (temperatures[0] + temperatures[1])
            + 0.5 * dt * (powers[0] + powers[1]) / capacity;
        let contrast = (0.5 * capacity * (temperatures[0] - temperatures[1])
            + 0.5 * dt * (powers[0] - powers[1])) / (capacity + 2.0 * g_dt);
        for (i, expected) in [mean + contrast, mean - contrast].into_iter().enumerate() {
            let t = result.step.thermal.temperature_k(bound.cell_ids()[i]).unwrap();
            near(t, expected, 3e-9);
            near(result.step.cell_energy_j[i].1, powers[i] * dt, 1e-12);
            assert_eq!(result.endpoint_optical_states[i], source.at_query_point(&card, &point(t)).unwrap());
        }
        near(result.step.heat_report.energy.stored_energy_change_j, result.step.optical.absorbed_j, 3e-11);
        let next = owner.advance_to(
            &result.step.thermal, 2.0 * dt, drive, conduction(), controls, &mut workspace, || false,
        ).unwrap();
        assert_eq!(next.applied_optical_states, result.endpoint_optical_states);
        assert_ne!(next.applied_optical_states, result.applied_optical_states);
        let mut cancelled_polls = 0;
        assert!(matches!(owner.advance_to(
            &initial, dt, drive, conduction(), controls, &mut workspace,
            || { cancelled_polls += 1; cancelled_polls == polls },
        ), Err(LayeredWaveHeatingError::Cancelled)));
        assert_eq!(initial.cells(), &before);
        assert_eq!(initial.time_s(), 0.0);
    }
}

fn airy_absorption(water: bool, wavelength_m: f64, t: f64) -> f64 {
    let phase = if water {
        IapwsWaterRefraction::prepare(wavelength_m * 1e9).unwrap().sample(t, rho(t)).unwrap().phase_index
    } else {
        Corning7980Silica::at_temperature(t).unwrap().evaluate(wavelength_m * 1e9).unwrap().phase_index
    };
    let a = alpha(t, wavelength_m);
    let n = C64::new(phase, -a * wavelength_m / (4.0 * core::f64::consts::PI));
    let admittance_ratio = n.scale(M / mu(t));
    let r = (C64::ONE - admittance_ratio) / (C64::ONE + admittance_ratio);
    let angle = -core::f64::consts::TAU * n.re * DEPTH / wavelength_m;
    let p = C64::new(angle.cos(), angle.sin()).scale((-a * DEPTH / 2.0).exp());
    let denominator = C64::ONE - r * r * p * p;
    let reflected = r * (C64::ONE - p * p) / denominator;
    let transmitted = (C64::ONE - r * r) * p / denominator;
    1.0 - reflected.norm_sq() - transmitted.norm_sq()
}

#[test]
fn broadband_heat_matches_original_phase_airy_and_source_domain_exits_are_transactional() {
    let card = source_card(360.0);
    let outside = exterior();
    let wavelengths = [550e-9, 650e-9];
    let lines = wavelengths.map(|w| ThermoOpticalSpectralLine {
        angular_frequency_rad_s: core::f64::consts::TAU * LIGHT_SPEED / w,
        incidence_cosine: 1.0, polarization: Polarization::Te,
    });
    let limits = BroadbandThermoOpticalLimits {
        maximum_lines: 2, maximum_cell_line_pairs: 2, maximum_cell_factor_pairs: 2,
    };
    let amplitudes = [1e6_f64, 1.4e6];
    let drives = amplitudes.map(|a| LayeredWaveDrive::Coherent([C64::new(a, 0.0), C64::ZERO]));
    let dt = 1e-10;
    for water in [true, false] {
        let source = source(&card, water);
        let native = thermal_body(&card, 1);
        let laws = native.thermal_laws(false).unwrap();
        let bound = native.bind(&laws, ids()).unwrap();
        let initial = initial_at(&bound, &[295.0]);
        let t0 = initial.temperature_k(bound.cell_ids()[0]).unwrap();
        let bindings = [ThermoOpticalLayer { source: &source, card: &card }];
        let owner = BroadbandThermoOpticalWaveHeating::try_new(
            &bound, &outside, DielectricResponse::Lossless(&outside),
            &lines, &bindings, limits, || false,
        ).unwrap();
        let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
        let controls = LayeredThermalSolveControls::default();
        let result = owner.advance_to(
            &initial, dt, &drives, conduction(), controls, &mut workspace, || false,
        ).unwrap();
        let mut energy = 0.0;
        for (i, &w) in wavelengths.iter().enumerate() {
            let power = amplitudes[i].powi(2) / (2.0 * M * LIGHT_SPEED) * airy_absorption(water, w, t0);
            near(result.line_power[i].absorbed_w_m2, power, 2e-11 * power);
            energy += AREA * power * dt;
        }
        near(result.step.optical.absorbed_j, energy, 2e-11 * energy);
        let endpoint = result.step.thermal.temperature_k(bound.cell_ids()[0]).unwrap();
        near(endpoint, t0 + energy / (RHO * AREA * DEPTH * CP), 3e-9);
        near(result.endpoint_optical_states[0].source().absorption_per_m_at(550e-9).unwrap(),
            alpha(endpoint, 550e-9), 1e-8);
        assert!(endpoint > t0 + 0.1 && endpoint < 300.0);
        if !water {
            // The thermal and alpha sources allow >300 K; the original CHARMS
            // phase law must still refuse BEFORE publishing the heated state.
            let near_boundary = initial_at(&bound, &[299.95]);
            let before = near_boundary.cells().clone();
            assert!(matches!(owner.advance_to(
                &near_boundary, dt, &drives, conduction(), controls, &mut workspace, || false,
            ), Err(LayeredWaveHeatingError::Electromagnetic(_))));
            assert_eq!(near_boundary.cells(), &before);
        }
    }
    // An independently selected alpha domain is equally binding for water.
    let bounded_card = source_card(295.05);
    let source = source(&bounded_card, true);
    let native = thermal_body(&bounded_card, 1);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let initial = initial_at(&bound, &[295.0]);
    let before = initial.cells().clone();
    let bindings = [ThermoOpticalLayer { source: &source, card: &bounded_card }];
    let owner = BroadbandThermoOpticalWaveHeating::try_new(
        &bound, &outside, DielectricResponse::Lossless(&outside), &lines, &bindings, limits, || false,
    ).unwrap();
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    assert!(owner.advance_to(
        &initial, dt, &drives, conduction(), LayeredThermalSolveControls::default(), &mut workspace, || false,
    ).is_err());
    assert_eq!(initial.cells(), &before);
    assert_eq!(initial.time_s(), 0.0);
}
