//! Conservative finite material storage with quasi-steady original-phase ports.
#![cfg(all(feature = "material", feature = "experimental-euler"))]

#[allow(dead_code)]
#[path = "../../fs-material/tests/common/nonideal_flash.rs"]
mod common;

use common::{MIXTURE_MASS, source};
use fs_gas::chemistry::ThermalConductance;
use fs_gas::material::phase_heat_exchanger::{
    MixtureHeatExchangerError as StreamError, MixtureHeatExchangerInlet as Inlet,
    wall::{
        FiniteWallStreamControl as Control, FiniteWallStreamError as E,
        MixtureWallPort as Port, WallStreamMechanicalConstraint as Mechanical,
        advance_wall_with_mixture_streams,
    },
    wall_stream::{MixtureWallStreamControl, exchange_mixture_stream_with_wall},
};
use fs_gas::multicomponent::wall_exchange::{
    MaterialWallCaloric, SharedWallStorage, WallDrag, WallThermal,
};
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::flash::{
    ThermochemicalFlashError, equilibrium::PhaseEnthalpyControls,
};
use fs_qty::{MassFlowRate, Pressure, Temperature, Time};
use fs_thermal::enthalpy::{
    Boundary as ThermalBoundary, CellGeometry, CellId, DepositedCell, DepositionError,
    MaterialId, ThermalLaw, ThermalStorageBasis,
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "{actual} != {expected}; difference {}", actual - expected);
}
fn calorics() -> PhaseEnthalpyControls {
    PhaseEnthalpyControls::new([Temperature::new(295.0), Temperature::new(405.0)])
}
fn port(model: &ThermochemicalPengRobinsonModel, temperature: f64, conductance: f64) -> Port<'_> {
    Port {
        inlet: Inlet { model, temperature: Temperature::new(temperature),
            pressure: Pressure::new(1.0), mass_flow: MassFlowRate::new(0.5) },
        conductance: ThermalConductance::new(conductance),
        controls: MixtureWallStreamControl::new([-50_000.0, 50_000.0], calorics()),
    }
}
fn control() -> Control { Control::new([-1e6, 1e6]) }
fn effective_conductance() -> f64 {
    500.0 / (1.0 + 500.0 / (0.5 * 4.0 * R / MIXTURE_MASS))
}

/// Independent analytic original material: constant heat capacity or an exact
/// latent interval, with an optional pressure-constrained density and entropy.
#[derive(Clone)]
struct AnalyticWall {
    basis: ThermalStorageBasis,
    latent: bool,
    entropy: bool,
    density: bool,
}
impl AnalyticWall {
    fn constant(basis: ThermalStorageBasis) -> Self {
        Self { basis, latent: false, entropy: true, density: true }
    }
    fn range(&self) -> [f64; 2] {
        if self.latent { [-40_000.0, 200_000.0] } else { [-750_000.0, 750_000.0] }
    }
    fn entropy_anchor(&self, h: f64) -> Result<f64, DepositionError> {
        let t = self.temperature_k(h)?;
        Ok(if !self.latent { 5000.0 * (t / 300.0).ln() }
        else if h < 0.0 { 1000.0 * (t / 340.0).ln() }
        else if h <= 100_000.0 { h / 340.0 }
        else { 100_000.0 / 340.0 + 1000.0 * (t / 340.0).ln() })
    }
}
impl ThermalLaw for AnalyticWall {
    fn storage_basis(&self) -> ThermalStorageBasis { self.basis }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let [lo, hi] = self.range();
        if !(lo..=hi).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(if !self.latent { 300.0 + h / 5000.0 }
        else if h < 0.0 { 340.0 + h / 1000.0 }
        else if h <= 100_000.0 { 340.0 }
        else { 340.0 + (h - 100_000.0) / 1000.0 })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        if self.latent { 0.001 } else { 1.0 / 5000.0 }
    }
    fn bulk_density_kg_m3(&self, h: f64) -> Result<Option<f64>, DepositionError> {
        if !self.density { return Ok(None); }
        Ok(Some(1000.0 / (1.0 + 0.001 * (self.temperature_k(h)? - 300.0))))
    }
    fn specific_entropy_change_j_kg_k(&self, a: f64, b: f64) -> Result<Option<f64>, DepositionError> {
        if !self.entropy { return Ok(None); }
        if !self.latent {
            let ta = self.temperature_k(a)?;
            self.temperature_k(b)?;
            return Ok(Some(5000.0 * ((b - a) / (5000.0 * ta)).ln_1p()));
        }
        Ok(Some(self.entropy_anchor(b)? - self.entropy_anchor(a)?))
    }
}
fn wall(law: AnalyticWall, initial_h: f64) -> SharedWallStorage<'static> {
    let range = law.range();
    let original = MaterialWallCaloric::from_thermal_law(law,
        DepositedCell {
            id: CellId(109), segment_id: 3, body_id: 4, material_id: MaterialId(17),
            mass_kg: 2.0, specific_h_j_kg: initial_h,
            geometry: CellGeometry { start_m: [0.0; 3], end_m: [1.0, 0.0, 0.0],
                volume_m3: 0.002, surface_area_m2: 1.0 },
            boundary: ThermalBoundary { ambient_k: 300.0, htc_w_m2_k: 0.0,
                reservoir: None, source_w: 0.0 },
        }, range).unwrap();
    SharedWallStorage {
        thermal: original.initial_contact(0.0, WallDrag::None).unwrap().thermal,
        heat_storage: Some(original),
    }
}
fn constant_wall(t: f64, basis: ThermalStorageBasis) -> SharedWallStorage<'static> {
    wall(AnalyticWall::constant(basis), 5000.0 * (t - 300.0))
}
fn wall_h(wall: &SharedWallStorage<'_>) -> f64 {
    match wall.thermal {
        WallThermal::Caloric { mass_kg, specific_energy_j_kg, conductance_w_per_k } => {
            assert_eq!(mass_kg, 2.0);
            assert_eq!(conductance_w_per_k, 0.0);
            specific_energy_j_kg
        }
        _ => panic!("original energy-coordinate wall was lost"),
    }
}

#[test]
fn g1_original_wall_and_stream_match_the_analytic_backward_euler_balance() {
    let original = source();
    let k = effective_conductance();
    let dt = 20.0;
    let cw = 10_000.0;
    let expected_wall = (cw * 400.0 + dt * k * 300.0) / (cw + dt * k);
    let expected_heat = k * (expected_wall - 300.0);
    let expected_out = 300.0 + expected_heat / (0.5 * 4.0 * R / MIXTURE_MASS);
    for (basis, mechanics) in [
        (ThermalStorageBasis::InternalEnergy, Mechanical::IsochoricInternalEnergy),
        (ThermalStorageBasis::Enthalpy, Mechanical::NegligibleExpansionEnthalpy),
    ] {
        let mut storage = constant_wall(400.0, basis);
        let result = advance_wall_with_mixture_streams(&mut storage, &[port(&original, 300.0, 500.0)],
            Time::new(dt), mechanics, control(), || false).unwrap();
        // The reference is the exact ideal Cp reduction; finite 1-Pa PR
        // departures remain in the production gas rather than being removed.
        close(result.wall_temperature_k[1], expected_wall, 1e-4);
        close(result.streams[0].enthalpy_flow_change_w, expected_heat, 0.01);
        close(result.streams[0].outlet_equilibrium.state().temperature_k(), expected_out, 1e-4);
        close(result.wall_storage_change_j, cw * (expected_wall - 400.0), 0.1);
        close(result.wall_storage_change_j, 2.0 * (wall_h(&storage) - 500_000.0), 1e-9);
        close(result.wall_entropy_change_j_k, cw * (expected_wall / 400.0).ln(), 1e-4);
        let stream_s = dt * 0.5 * 4.0 * R / MIXTURE_MASS * (expected_out / 300.0).ln();
        close(result.stream_entropy_change_j_k[0], stream_s, 1e-4);
        assert_eq!(result.wall_storage_basis, basis);
        assert_eq!(result.wall_pressure_volume_work_j, 0.0);
        assert_eq!(result.wall_internal_energy_change_j, result.wall_storage_change_j);
        assert!(result.energy_residual_j.abs() < 1e-3);
        assert!(result.maximum_contact_residual_j < 1e-3);
        assert!(result.entropy_production_j_k > 0.0);
        assert!(result.wall_temporal_entropy_defect_j_k > 49.0);
        close(result.entropy_production_j_k,
            result.wall_temporal_entropy_defect_j_k + dt * result.streams[0].entropy_production_w_k,
            1e-5);
    }
}


#[test]
fn g1_finite_u_wall_couples_to_the_independent_two_phase_flow_reference() {
    let original = source();
    let target_wall = 228.945_170_130_525_04;
    let target_heat = 8401.552_172_252_039;
    // Choose the time interval from the independent PR mixed-stream endpoint
    // and the exact original wall storage. Both energy equations then share
    // this reference, including the source vapor fraction and entropy.
    let dt = 10_000.0 * (230.0 - target_wall) / target_heat;
    let mut storage = constant_wall(230.0, ThermalStorageBasis::InternalEnergy);
    let cryogenic = Port {
        inlet: Inlet { model: &original, temperature: Temperature::new(210.0),
            pressure: Pressure::new(3e6), mass_flow: MassFlowRate::new(0.5) },
        conductance: ThermalConductance::new(500.0),
        controls: MixtureWallStreamControl::new([-10_000.0, 10_000.0],
            PhaseEnthalpyControls::new([Temperature::new(200.0), Temperature::new(240.0)])),
    };
    let result = advance_wall_with_mixture_streams(
        &mut storage, &[cryogenic], Time::new(dt),
        Mechanical::IsochoricInternalEnergy, control(), || false).unwrap();
    close(result.wall_temperature_k[1], target_wall, 3e-7);
    close(result.streams[0].enthalpy_flow_change_w, target_heat, 1e-3);
    let outlet = result.streams[0].outlet_equilibrium.state();
    close(outlet.temperature_k(), 212.142_065_786_020_96, 3e-7);
    close(outlet.vapor_mass_fraction().unwrap(), 0.120_681_519_207_042_33, 3e-8);
    close(result.entropy_production_j_k,
        10_000.0 * (target_wall / 230.0).ln() + dt * 39.806_095_914_368_5, 3e-6);
    assert!(result.wall_temporal_entropy_defect_j_k > 0.1);
    assert_eq!(result.wall_storage_basis, ThermalStorageBasis::InternalEnergy);
    assert!(result.energy_residual_j.abs() < 5e-5);
}

#[test]
fn g1_time_refinement_recovers_continuous_wall_cooling_and_reduces_temporal_entropy() {
    let original = source();
    let k = effective_conductance();
    let exact = 300.0 + 100.0 * (-k * 20.0 / 10_000.0).exp();
    let mut errors = Vec::new();
    let mut defects = Vec::new();
    for n in [1_u32, 2, 4] {
        let mut storage = constant_wall(400.0, ThermalStorageBasis::InternalEnergy);
        let dt = 20.0 / f64::from(n);
        let mut total_heat = 0.0;
        let mut defect = 0.0;
        for _ in 0..n {
            let result = advance_wall_with_mixture_streams(
                &mut storage, &[port(&original, 300.0, 500.0)], Time::new(dt),
                Mechanical::IsochoricInternalEnergy, control(), || false).unwrap();
            total_heat += result.stream_heat_j[0];
            defect += result.wall_temporal_entropy_defect_j_k;
        }
        let actual = storage.heat_storage.as_ref().unwrap().temperature_k(storage.thermal).unwrap();
        let discrete = 300.0 + 100.0 / (1.0 + k * dt / 10_000.0).powi(n as i32);
        close(actual, discrete, 1e-4);
        close(2.0 * (wall_h(&storage) - 500_000.0) + total_heat, 0.0, 0.01);
        errors.push((actual - exact).abs());
        defects.push(defect);
    }
    assert!(errors[2] < errors[1] && errors[1] < errors[0]);
    assert!(errors[2] < 0.35 * errors[0]);
    assert!(defects[2] < defects[1] && defects[1] < defects[0]);
}

#[test]
fn g1_latent_wall_inventory_advances_at_constant_temperature_and_crosses_the_whole_plateau() {
    let original = source();
    let law = AnalyticWall { latent: true, ..AnalyticWall::constant(ThermalStorageBasis::Enthalpy) };
    let mut plateau = wall(law.clone(), 50_000.0);
    let result = advance_wall_with_mixture_streams(
        &mut plateau, &[port(&original, 300.0, 500.0)], Time::new(1.0),
        Mechanical::NegligibleExpansionEnthalpy, control(), || false).unwrap();
    assert_eq!(result.wall_temperature_k, [340.0, 340.0]);
    assert!((0.0..50_000.0).contains(&wall_h(&plateau)));
    close(result.stream_heat_j[0], effective_conductance() * 40.0, 0.01);
    close(result.wall_entropy_change_j_k, result.wall_storage_change_j / 340.0, 1e-10);
    close(result.wall_temporal_entropy_defect_j_k, 0.0, 1e-10);
    let mut crossing = wall(law, 110_000.0);
    let result = advance_wall_with_mixture_streams(
        &mut crossing, &[port(&original, 300.0, 500.0)], Time::new(20.0),
        Mechanical::NegligibleExpansionEnthalpy, control(), || false).unwrap();
    let expected = (900_000.0 + 20.0 * effective_conductance() * 300.0)
        / (2000.0 + 20.0 * effective_conductance());
    assert!(wall_h(&crossing) < 0.0);
    close(result.wall_temperature_k[1], expected, 1e-4);
    close(result.wall_storage_change_j, 2000.0 * (expected - 340.0) - 220_000.0, 0.1);
    assert!(result.entropy_production_j_k > 0.0);
}

#[test]
fn g1_isobaric_h_wall_reports_original_density_work_and_internal_energy_separately() {
    let original = source();
    let mut storage = constant_wall(400.0, ThermalStorageBasis::Enthalpy);
    let result = advance_wall_with_mixture_streams(
        &mut storage, &[port(&original, 300.0, 500.0)], Time::new(20.0),
        Mechanical::IsobaricEnthalpy { pressure: Pressure::new(1e8) },
        control(), || false).unwrap();
    let delta_t = result.wall_temperature_k[1] - 400.0;
    let expected_work = 200.0 * delta_t; // p m alpha / rho_ref * Delta T.
    close(result.wall_pressure_volume_work_j, expected_work, 1e-6);
    assert!(result.wall_pressure_volume_work_j < -7000.0);
    close(result.wall_internal_energy_change_j,
        result.wall_storage_change_j - expected_work, 1e-6);
    close(result.wall_internal_energy_change_j + result.wall_pressure_volume_work_j
        + result.stream_heat_j[0], 0.0, 1e-3);
    let volume = result.wall_volume_m3.unwrap();
    close(volume[0], 0.0022, 1e-15);
    close(volume[1], 0.002 * (1.0 + 0.001 * (result.wall_temperature_k[1] - 300.0)), 1e-15);
    // Original reference geometry is not silently expanded by a thermal adapter.
    assert_eq!(storage.heat_storage.as_ref().unwrap().initial_cell().unwrap().geometry.volume_m3, 0.002);
}

#[test]
fn g0_zero_duration_returns_instantaneous_outlets_without_changing_wall_inventory() {
    let original = source();
    let mut storage = constant_wall(400.0, ThermalStorageBasis::InternalEnergy);
    let before = storage.clone();
    let result = advance_wall_with_mixture_streams(
        &mut storage, &[port(&original, 300.0, 500.0)], Time::new(0.0),
        Mechanical::IsochoricInternalEnergy, control(), || false).unwrap();
    assert_eq!(storage, before);
    assert_eq!(result.wall_storage_change_j, 0.0);
    assert_eq!(result.stream_heat_j, vec![0.0]);
    assert_eq!(result.stream_entropy_change_j_k, vec![0.0]);
    assert_eq!(result.entropy_production_j_k, 0.0);
    assert_eq!(result.wall_evaluations, 1);
    assert!(result.streams[0].enthalpy_flow_change_w > 29_000.0);
    assert!(result.streams[0].outlet_equilibrium.state().temperature_k() > 340.0);
}

#[test]
fn g1_nonzero_stream_throughflow_can_leave_the_shared_wall_exactly_at_rest() {
    let original = source();
    let cold = port(&original, 300.0, 500.0);
    let heat = exchange_mixture_stream_with_wall(cold.inlet, Temperature::new(350.0),
        cold.conductance, cold.controls, || false).unwrap().enthalpy_flow_change_w;
    // Construct a physically balanced hot feed using the same original PH
    // increment, then derive its supplied conductance from the contact law.
    // This tests the no-storage-change semantics independently of outer-root
    // iteration; it does not assert an idealized Cp symmetry for the PR gas.
    let hot_path = original.prepare_phase_enthalpy_path(
        Temperature::new(400.0), Pressure::new(1.0), Pressure::new(1.0),
        calorics(), || false).unwrap();
    let hot_out = hot_path.solve_increment(-heat / 0.5, 128, || false).unwrap();
    let hot = port(&original, 400.0, heat / (hot_out.state().temperature_k() - 350.0));
    let mut storage = constant_wall(350.0, ThermalStorageBasis::InternalEnergy);
    let before = storage.clone();
    let mut controls = control();
    controls.absolute_energy_tolerance_j = 1e-3;
    controls.relative_energy_tolerance = 0.0;
    let result = advance_wall_with_mixture_streams(
        &mut storage, &[cold, hot], Time::new(1.0),
        Mechanical::IsochoricInternalEnergy, controls, || false).unwrap();
    assert_eq!(storage, before);
    assert_eq!(result.wall_storage_change_j, 0.0);
    assert_eq!(result.wall_evaluations, 1);
    assert!(result.streams[0].enthalpy_flow_change_w > 14_000.0);
    assert!(result.streams[1].enthalpy_flow_change_w < -14_000.0);
    close(result.stream_heat_j.iter().sum(), 0.0, 1e-3);
    assert!(result.entropy_production_j_k > 0.0);
    assert_eq!(result.wall_entropy_change_j_k, 0.0);
}

#[test]
fn g0_unsupported_storage_entropy_and_volume_refuse_without_publication() {
    let original = source();
    for (law, mode) in [
        (AnalyticWall { entropy: false, ..AnalyticWall::constant(ThermalStorageBasis::InternalEnergy) },
            Mechanical::IsochoricInternalEnergy),
        (AnalyticWall { density: false, ..AnalyticWall::constant(ThermalStorageBasis::Enthalpy) },
            Mechanical::IsobaricEnthalpy { pressure: Pressure::new(1e5) }),
        (AnalyticWall::constant(ThermalStorageBasis::Enthalpy), Mechanical::IsochoricInternalEnergy),
        (AnalyticWall::constant(ThermalStorageBasis::DeclaredSensibleStorage),
            Mechanical::NegligibleExpansionEnthalpy),
    ] {
        let mut storage = wall(law, 500_000.0);
        let before = storage.clone();
        assert!(matches!(advance_wall_with_mixture_streams(
            &mut storage, &[port(&original, 300.0, 500.0)], Time::new(1.0),
            mode, control(), || false), Err(E::Unsupported(_))));
        assert_eq!(storage, before);
    }
}

#[test]
fn g0_shared_work_source_limits_and_final_cancellation_leave_wall_atomic() {
    let original = source();
    let baseline = constant_wall(400.0, ThermalStorageBasis::InternalEnergy);
    let solve = |storage: &mut SharedWallStorage<'_>, ports: &[Port<'_>],
        controls: Control, cancelled: &mut dyn FnMut() -> bool| {
        advance_wall_with_mixture_streams(storage, ports, Time::new(1.0),
            Mechanical::IsochoricInternalEnergy, controls, cancelled)
    };
    let ports = [port(&original, 300.0, 500.0)];
    let mut storage = baseline.clone();
    let mut polls = 0;
    let accepted = solve(&mut storage, &ports, control(), &mut || { polls += 1; false }).unwrap();
    let mut replay = baseline.clone();
    let mut count = 0;
    assert_eq!(solve(&mut replay, &ports, control(), &mut || {
        count += 1; count == polls
    }).unwrap_err(), E::Cancelled);
    assert_eq!(replay, baseline);
    let mut controls = control();
    controls.maximum_wall_evaluations = 1;
    assert!(matches!(solve(&mut replay, &ports, controls, &mut || false), Err(E::WorkLimit { wall_evaluations: 1, .. })));
    assert_eq!(replay, baseline);
    controls = control();
    controls.maximum_pt_evaluations = accepted.pt_evaluations;
    controls.maximum_duty_evaluations = accepted.duty_evaluations;
    let exact = solve(&mut replay, &ports, controls, &mut || false).unwrap();
    assert_eq!(exact.wall_storage_change_j, accepted.wall_storage_change_j);
    assert_eq!(exact.pt_evaluations, accepted.pt_evaluations);
    replay = baseline.clone();
    controls.maximum_pt_evaluations -= 1;
    match solve(&mut replay, &ports, controls, &mut || false).unwrap_err() {
        E::Stream { completed_pt_evaluations: outer, source: StreamError::Equilibrium {
            completed_pt_evaluations: inner,
            source: ThermochemicalFlashError::WorkLimit { flash_evaluations }, .. }, .. } =>
            assert_eq!(outer + inner + flash_evaluations, controls.maximum_pt_evaluations),
        E::Stream { completed_pt_evaluations: outer, source: StreamError::WorkLimit { pt_evaluations, .. }, .. } =>
            assert_eq!(outer + pt_evaluations, controls.maximum_pt_evaluations),
        E::WorkLimit { pt_evaluations, .. } => assert_eq!(pt_evaluations, controls.maximum_pt_evaluations),
        other => panic!("unexpected shared work refusal: {other:?}"),
    }
    assert_eq!(replay, baseline);
    let mut narrow = port(&original, 300.0, 500.0);
    narrow.controls.equilibrium.enthalpy.temperature_bracket =
        [Temperature::new(300.0), Temperature::new(301.0)];
    assert!(matches!(solve(&mut replay, &[narrow], control(), &mut || false), Err(E::Stream { .. })));
    assert_eq!(replay, baseline);
    assert!(matches!(solve(&mut replay, &ports, Control::new([-1.0, 1.0]), &mut || false),
        Err(E::NotBracketed { .. })));
    assert_eq!(replay, baseline);
}

#[path = "phase_equilibrium_finite_wall_streams/spatial.rs"]
mod spatial;
