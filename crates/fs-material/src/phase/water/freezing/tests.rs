use super::*;
use fs_thermochem::water::ice::{NORMAL_PRESSURE_PA, melting_pressure_derivative};

fn controls(minimum: f64, maximum: f64) -> FreezingWaterControls {
    FreezingWaterControls {
        minimum_temperature_k: minimum,
        maximum_temperature_k: maximum,
        enthalpy_tolerance_j_kg: 1e-6,
        maximum_iterations: 96,
        maximum_entropy_seam_j_kg_k: 0.1,
        maximum_gibbs_mismatch_j_kg: 30.0,
        melting: IceMeltingSolve {
            pressure_tolerance_pa: 1e-5,
            temperature_tolerance_k: 1e-10,
            maximum_iterations: 96,
        },
    }
}
fn chart() -> FreezingWater {
    FreezingWater::new(
        Pressure::new(NORMAL_PRESSURE_PA),
        controls(243.15, 313.15),
        || false,
    )
    .unwrap()
}
fn at_t(chart: &FreezingWater, t: f64) -> FreezingWaterState {
    chart
        .state_at_temperature(Temperature::new(t), || false)
        .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}, tolerance {tolerance:.4e}"
    );
}

#[test]
fn original_sources_meet_on_pressure_dependent_fusion_with_explicit_matching_defects() {
    let c = chart();
    let pair = c.melting();
    // Released normal-pressure melting value rounded to nine decimals.
    near(pair.temperature().value(), 273.152519, 5e-7);
    assert!((333_000.0..334_000.0).contains(&pair.latent_enthalpy_j_kg()));
    assert!((916.0..918.0).contains(&pair.ice().density().value()));
    assert!((999.0..1001.0).contains(&pair.liquid().density().value()));
    let raw_ice = IceIhState::at_pressure_temperature(pair.pressure(), pair.temperature()).unwrap();
    let raw_liquid = If97State::at_pressure_temperature(
        pair.pressure(),
        pair.temperature(),
        If97Region::Region1,
    )
    .unwrap();
    assert_eq!(*pair.ice(), raw_ice);
    assert_eq!(*pair.liquid(), raw_liquid);
    near(
        pair.entropy_defect_j_kg_k(),
        pair.latent_enthalpy_j_kg() / pair.temperature().value() - pair.latent_entropy_j_kg_k(),
        0.0,
    );
    near(
        pair.gibbs_mismatch_j_kg(),
        raw_liquid.specific_gibbs_energy().value() - raw_ice.specific_gibbs_energy().value(),
        0.0,
    );
    assert!(pair.entropy_defect_j_kg_k() != 0.0);
    assert!(pair.gibbs_mismatch_j_kg() != 0.0);
    let p = pair.pressure().value();
    assert!(
        pair.melting_solution().pressure_residual_pa().abs()
            <= c.controls.melting.pressure_tolerance_pa
    );
    let low =
        FreezingWater::new(Pressure::new(10_000.0), controls(243.15, 313.15), || false).unwrap();
    assert!(low.melting().temperature().value() > pair.temperature().value());
    let dv = raw_liquid.specific_volume().value() - raw_ice.specific_volume().value();
    let clapeyron = pair.latent_enthalpy_j_kg() / (pair.temperature().value() * dv);
    let correlation_slope = melting_pressure_derivative(pair.temperature())
        .unwrap()
        .value();
    assert!(clapeyron < 0.0);
    near(clapeyron / correlation_slope, 1.0, 5e-4);
    for t in [243.15, 253.15, 270.0, 280.0, 293.15, 313.15] {
        let s = at_t(&c, t);
        let solved = c
            .state_at_specific_enthalpy(s.specific_enthalpy_j_kg(), || false)
            .unwrap();
        near(solved.state().temperature_k(), t, 1e-8);
        near(s.pressure_closure_defect_j_kg(), 0.0, 2e-9);
        assert_eq!(s.pressure_pa(), p);
    }
}

#[test]
fn freezing_thawing_retains_mass_physical_phase_volumes_and_independent_energy() {
    let c = chart();
    let (warm, cold) = (at_t(&c, 293.15), at_t(&c, 253.15));
    let mass = 0.2;
    let heat = mass * (cold.specific_enthalpy_j_kg() - warm.specific_enthalpy_j_kg());
    let frozen = c.advance_heat(warm, mass, heat, 2e-6, || false).unwrap();
    near(frozen.state.temperature_k(), 253.15, 1e-8);
    assert_eq!(frozen.initial_phase_masses_kg, [0.0, mass]);
    assert_eq!(frozen.phase_masses_kg, [mass, 0.0]);
    assert!(heat < 0.0 && frozen.sensible_enthalpy_change_j < 0.0);
    near(
        frozen.fusion_enthalpy_change_j,
        -mass * c.melting.latent_enthalpy_j_kg(),
        1e-10,
    );
    near(
        frozen.sensible_enthalpy_change_j + frozen.fusion_enthalpy_change_j,
        heat,
        2e-6,
    );
    near(
        frozen.internal_energy_change_j + frozen.boundary_work_j,
        heat,
        2e-6,
    );
    assert!(frozen.volume_change_m3 > 0.0);
    assert!(frozen.boundary_work_j > 0.0);
    assert!((1.08..1.10).contains(&(frozen.volumes_m3[1] / frozen.volumes_m3[0])));
    near(frozen.enthalpy_residual_j, 0.0, 2e-6);
    near(frozen.energy_residual_j, 0.0, 2e-6);
    near(frozen.represented_enthalpy_residual_j, 0.0, 2e-6);
    near(frozen.represented_energy_residual_j, 0.0, 2e-6);
    near(
        frozen.thermal_entropy_change_j_k,
        frozen.source_entropy_change_j_k + frozen.entropy_seam_correction_change_j_k,
        1e-12,
    );
    let thawed = c
        .advance_heat(frozen.state, mass, -heat, 2e-6, || false)
        .unwrap();
    near(thawed.state.temperature_k(), warm.temperature_k(), 1e-8);
    near(
        thawed.volume_change_m3 + frozen.volume_change_m3,
        0.0,
        2e-15,
    );
    near(
        thawed.thermal_entropy_change_j_k + frozen.thermal_entropy_change_j_k,
        0.0,
        1e-8,
    );

    let mixed = c.state_at_liquid_fraction(0.25, || false).unwrap();
    assert_eq!(mixed.phase(), FreezingWaterPhase::IceLiquid);
    assert_eq!(mixed.liquid_mass_fraction(), 0.25);
    assert_eq!(mixed.ice_mass_fraction(), 0.75);
    let v = 0.75 * c.melting.ice.specific_volume().value()
        + 0.25 * c.melting.liquid.specific_volume().value();
    near(mixed.specific_volume_m3_kg(), v, 1e-18);
    near(
        mixed.liquid_volume_fraction() + mixed.ice_volume_fraction(),
        1.0,
        3e-16,
    );
    assert!(mixed.liquid_volume_fraction() < mixed.liquid_mass_fraction());
}

fn simpson(a: f64, b: f64, read: impl Fn(f64) -> f64) -> f64 {
    let n = 512;
    let dx = (b - a) / n as f64;
    let mut sum = read(a) + read(b);
    for i in 1..n {
        sum += (if i % 2 == 0 { 2.0 } else { 4.0 }) * read(a + i as f64 * dx);
    }
    sum * dx / 3.0
}

#[test]
fn independent_cp_and_cp_over_t_integrals_close_sensible_and_fusion_accounts() {
    let c = chart();
    let (a, b) = (at_t(&c, 253.15), at_t(&c, 303.15));
    let tm = c.melting.temperature().value();
    let ice_cp = |t| {
        IceIhState::at_pressure_temperature(c.melting.pressure(), Temperature::new(t))
            .unwrap()
            .cp()
            .value()
    };
    let liquid_cp = |t| {
        If97State::at_pressure_temperature(
            c.melting.pressure(),
            Temperature::new(t),
            If97Region::Region1,
        )
        .unwrap()
        .cp()
        .value()
    };
    let sensible_h =
        simpson(a.temperature_k(), tm, ice_cp) + simpson(tm, b.temperature_k(), liquid_cp);
    let sensible_s = simpson(a.temperature_k(), tm, |t| ice_cp(t) / t)
        + simpson(tm, b.temperature_k(), |t| liquid_cp(t) / t);
    let heat = b.specific_enthalpy_j_kg() - a.specific_enthalpy_j_kg();
    let step = c.advance_heat(a, 1.0, heat, 2e-6, || false).unwrap();
    near(step.sensible_enthalpy_change_j, sensible_h, 2e-5);
    near(
        step.fusion_enthalpy_change_j,
        c.melting.latent_enthalpy_j_kg(),
        1e-10,
    );
    near(
        step.thermal_entropy_change_j_k,
        sensible_s + c.melting.latent_enthalpy_j_kg() / tm,
        2e-7,
    );
    near(
        c.specific_entropy_change_j_kg_k(
            a.specific_enthalpy_j_kg(),
            b.specific_enthalpy_j_kg(),
            || false,
        )
        .unwrap(),
        step.thermal_entropy_change_j_k,
        2e-8,
    );
    near(
        c.specific_entropy_change_j_kg_k(
            b.specific_enthalpy_j_kg(),
            a.specific_enthalpy_j_kg(),
            || false,
        )
        .unwrap(),
        -step.thermal_entropy_change_j_k,
        2e-8,
    );
}

#[test]
fn subdivisions_latent_only_windows_and_small_heat_preserve_the_original_path() {
    let c = chart();
    let (a, b) = (at_t(&c, 253.15), at_t(&c, 293.15));
    let mass = 0.125;
    let heat = mass * (b.specific_enthalpy_j_kg() - a.specific_enthalpy_j_kg());
    let direct = c.advance_heat(a, mass, heat, 2e-6, || false).unwrap();
    let mut state = a;
    let mut total = [0.0; 4];
    for target in [
        c.state_at_liquid_fraction(0.0, || false).unwrap(),
        c.state_at_liquid_fraction(0.3, || false).unwrap(),
        c.state_at_liquid_fraction(1.0, || false).unwrap(),
        b,
    ] {
        let q = mass * (target.specific_enthalpy_j_kg() - state.specific_enthalpy_j_kg());
        let step = c.advance_heat(state, mass, q, 2e-6, || false).unwrap();
        for (sum, x) in total.iter_mut().zip([
            step.enthalpy_change_j,
            step.internal_energy_change_j,
            step.boundary_work_j,
            step.thermal_entropy_change_j_k,
        ]) {
            *sum += x;
        }
        state = step.state;
    }
    near(state.temperature_k(), direct.state.temperature_k(), 1e-9);
    for (sum, x) in total.into_iter().zip([
        direct.enthalpy_change_j,
        direct.internal_energy_change_j,
        direct.boundary_work_j,
        direct.thermal_entropy_change_j_k,
    ]) {
        near(sum, x, 2e-7);
    }
    for state in [a, b, c.state_at_liquid_fraction(0.4, || false).unwrap()] {
        let zero = c.advance_heat(state, mass, 0.0, 0.0, || false).unwrap();
        assert_eq!(zero.state, state);
        assert_eq!(zero.inversion_iterations, 0);
        assert_eq!(zero.enthalpy_change_j, 0.0);
        let small = c.advance_heat(state, mass, 0.001, 2e-7, || false).unwrap();
        assert!(small.enthalpy_change_j > 0.0);
        near(small.enthalpy_change_j, 0.001, 2e-7);
    }
    let tm = c.melting.temperature().value();
    let plateau = FreezingWater::new(c.melting.pressure(), controls(tm, tm), || false).unwrap();
    let solid = plateau.state_at_liquid_fraction(0.0, || false).unwrap();
    for x in [1e-10, 0.1, 0.5, 1.0 - 1e-10, 1.0] {
        let target = plateau.state_at_liquid_fraction(x, || false).unwrap();
        let q = target.specific_enthalpy_j_kg() - solid.specific_enthalpy_j_kg();
        let step = plateau.advance_heat(solid, 1.0, q, 2e-7, || false).unwrap();
        assert_eq!(step.state.temperature_k(), tm);
        assert_eq!(step.sensible_enthalpy_change_j, 0.0);
        near(step.state.liquid_mass_fraction(), x, 3e-16);
        near(
            step.thermal_entropy_change_j_k,
            step.enthalpy_change_j / tm,
            2e-13,
        );
    }
}

#[test]
fn source_overlap_limits_matching_and_work_budgets_refuse_without_extrapolation() {
    let c = chart();
    let pm = melting_pressure(Temperature::new(273.15)).unwrap();
    let edge = FreezingWater::new(pm, controls(253.15, 293.15), || false).unwrap();
    assert_eq!(edge.melting.temperature().value(), 273.15);
    assert!(matches!(
        FreezingWater::new(
            Pressure::new(pm.value().next_up()),
            controls(253.15, 293.15),
            || false
        ),
        Err(FreezingWaterError::OutsideDomain)
    ));
    assert!(matches!(
        FreezingWater::new(
            Pressure::new(TRIPLE_POINT_PRESSURE_PA.next_down()),
            controls(253.15, 293.15),
            || false
        ),
        Err(FreezingWaterError::OutsideDomain)
    ));
    assert!(FreezingWater::new(c.melting.pressure(), controls(49.0, 293.15), || false).is_err());
    assert!(FreezingWater::new(c.melting.pressure(), controls(280.0, 293.15), || false).is_err());
    assert!(matches!(
        FreezingWater::new(c.melting.pressure(), controls(253.15, 400.0), || false),
        Err(FreezingWaterError::LiquidSource(_))
    ));
    let mut limits = controls(253.15, 293.15);
    limits.maximum_entropy_seam_j_kg_k = 0.0;
    assert!(matches!(
        FreezingWater::new(c.melting.pressure(), limits, || false),
        Err(FreezingWaterError::EntropySeam { .. })
    ));
    limits.maximum_entropy_seam_j_kg_k = 0.1;
    limits.maximum_gibbs_mismatch_j_kg = 0.0;
    assert!(matches!(
        FreezingWater::new(c.melting.pressure(), limits, || false),
        Err(FreezingWaterError::GibbsMismatch { .. })
    ));
    limits.maximum_gibbs_mismatch_j_kg = 30.0;
    limits.melting.maximum_iterations = 1;
    assert!(matches!(
        FreezingWater::new(c.melting.pressure(), limits, || false),
        Err(FreezingWaterError::IceSource(
            IceIhError::NotConverged { .. }
        ))
    ));
    assert_eq!(
        c.state_at_temperature(c.melting.temperature(), || false),
        Err(FreezingWaterError::AmbiguousMelting)
    );
    for x in [-0.01, 1.01, f64::NAN] {
        assert!(c.state_at_liquid_fraction(x, || false).is_err());
    }
    let [lo, hi] = c.enthalpy_bounds_j_kg();
    for h in [lo.next_down(), hi.next_up(), f64::NAN] {
        assert_eq!(
            c.state_at_specific_enthalpy(h, || false),
            Err(FreezingWaterError::OutsideDomain)
        );
    }
    let mut short = controls(243.15, 313.15);
    short.maximum_iterations = 1;
    short.enthalpy_tolerance_j_kg = 0.0;
    let short = FreezingWater::new(c.melting.pressure(), short, || false).unwrap();
    assert!(matches!(
        short.state_at_specific_enthalpy(at_t(&c, 256.77).specific_enthalpy_j_kg(), || false),
        Err(FreezingWaterError::NotConverged { .. })
    ));
}

#[test]
fn energy_representability_chart_identity_and_early_late_cancellation_are_atomic() {
    let c = chart();
    let a = at_t(&c, 253.15);
    let b = at_t(&c, 293.15);
    let q = b.specific_enthalpy_j_kg() - a.specific_enthalpy_j_kg();
    assert_eq!(
        c.advance_heat(a, 1.0, q, 2e-6, || true),
        Err(FreezingWaterError::Cancelled)
    );
    assert_eq!(
        FreezingWater::new(c.melting.pressure(), controls(243.15, 313.15), || true),
        Err(FreezingWaterError::Cancelled)
    );
    let mut polls = 0;
    c.advance_heat(a, 1.0, q, 2e-6, || {
        polls += 1;
        false
    })
    .unwrap();
    let mut actual = 0;
    assert_eq!(
        c.advance_heat(a, 1.0, q, 2e-6, || {
            actual += 1;
            actual == polls
        }),
        Err(FreezingWaterError::Cancelled)
    );
    let other =
        FreezingWater::new(Pressure::new(100_000.0), controls(243.15, 313.15), || false).unwrap();
    assert_eq!(
        other.advance_heat(a, 1.0, 0.0, 0.0, || false),
        Err(FreezingWaterError::ChartMismatch)
    );
    assert!(matches!(
        c.advance_heat(a, 1.0, f64::from_bits(1), 1.0, || false),
        Err(FreezingWaterError::Unrepresentable(_))
    ));
    let mixed = c.state_at_liquid_fraction(0.5, || false).unwrap();
    assert!(matches!(
        c.advance_heat(mixed, f64::from_bits(1), 0.0, 0.0, || false),
        Err(FreezingWaterError::Unrepresentable(_))
    ));
    assert!(c.advance_heat(a, 1.0, q, 0.0, || false).is_err());
    assert_eq!(a, at_t(&c, 253.15));
    assert_eq!(b, at_t(&c, 293.15));
}
