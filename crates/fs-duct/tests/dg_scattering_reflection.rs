//! Direct wave excitation of original DG fields; manufactured equation oracles.

use fs_duct::dg::{
    DgBoundary, DgDuct, DgFluidMaterial, DgMesh, DgScatteringResponse, PortId,
};
use fs_duct::{Duct, DuctError, Segment};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::gas::{GasSpec, GasState};
use fs_material::visco::GeneralizedMaxwell;
use fs_math::{c64::C64, det};
use fs_phs::WallPin;

const INLET: PortId = PortId(7);
const OUTLET: PortId = PortId(11);
const RADIUS: f64 = 0.0073;
const LENGTH: f64 = 0.5;
const WORK: usize = 100_000_000;

fn geometry() -> Duct {
    Duct { segments: vec![Segment::Cylinder { radius: RADIUS, length: LENGTH }] }
}

fn duct(
    medium: FluidAcousticState,
    elements: usize,
    order: usize,
) -> DgDuct<DgFluidMaterial<FluidAcousticState>> {
    DgDuct::from_material(
        &geometry(), &medium, DgMesh { elements, order },
        DgBoundary::Port(INLET), DgBoundary::Port(OUTLET),
    ).unwrap()
}

fn impedance(medium: FluidAcousticState) -> f64 {
    medium.characteristic_impedance_pa_s_m() / (core::f64::consts::PI * RADIUS * RADIUS)
}

fn cis(angle: f64) -> C64 {
    C64::new(det::cos(angle), det::sin(angle))
}

fn power_identity(response: &DgScatteringResponse) {
    let scale: f64 = response.waves.iter().map(|wave|
        wave.incident_power_w + wave.reflected_power_w).sum();
    for wave in &response.waves {
        let port = response.physical.ports.iter().find(|port| port.id == wave.id).unwrap();
        assert_eq!(wave.net_input_power_w.to_bits(), port.mean_input_power_w.to_bits());
        assert!(wave.power_residual_w.abs() < 2e-10 * scale.max(1e-30));
        let p = wave.incident_pressure_pa + wave.reflected_pressure_pa;
        let q = (wave.incident_pressure_pa - wave.reflected_pressure_pa)
            .scale(1.0 / wave.reference_impedance_pa_s_m3);
        assert!((p - port.pressure_pa).abs() < 1e-11 * p.abs().max(1.0));
        assert!((q - port.inflow_m3_s).abs() < 1e-11 * q.abs().max(1e-20));
    }
    assert!(response.physical.power_residual_w.abs() < 2e-9 * scale.max(1e-30));
    assert!(response.wave_power_residual_w.abs() < 2e-9 * scale.max(1e-30));
}

#[test]
fn two_port_reflection_and_transmission_converge_to_the_analytic_pipe_in_gas_and_liquid() {
    for medium in [
        FluidAcousticState::try_new(1.2, 343.0).unwrap(),
        FluidAcousticState::try_new(1000.0, 1480.0).unwrap(),
    ] {
        let z = impedance(medium);
        let phase = 2.4;
        let omega = phase * medium.sound_speed_m_s() / LENGTH;
        let mut errors = Vec::new();
        for elements in [1, 2, 4] {
            let dg = duct(medium, elements, 3);
            let factor = dg.scattering_factor(
                omega, &[(INLET, z), (OUTLET, z)], WORK, 1e-10, &mut || false,
            ).unwrap();
            let response = factor.response(&[(INLET, C64::from_re(1.0))], &mut || false).unwrap();
            assert_eq!(response.waves[1].incident_pressure_pa, C64::ZERO);
            let error = response.waves[0].reflected_pressure_pa.abs().max(
                (response.waves[1].reflected_pressure_pa - cis(phase)).abs(),
            );
            errors.push(error);
            power_identity(&response);
        }
        assert!(errors[1] < 0.6 * errors[0], "{errors:?}");
        assert!(errors[2] < 0.6 * errors[1], "{errors:?}");
        assert!(errors[2] < 1e-4, "{errors:?}");
    }
}

#[test]
fn direct_scattering_stays_finite_at_the_actual_discrete_duct_impedance_pole() {
    let medium = FluidAcousticState::try_new(1.2, 343.0).unwrap();
    let dg = duct(medium, 1, 1);
    let z = impedance(medium);
    // The original one-element linear closed-flow DG image has exactly one
    // nonzero oscillator pair. Its frequency follows from ||L||_F^2=2*omega^2,
    // without assuming that the finite mesh's pole equals continuum pi*c/L.
    let omega = det::sqrt(0.5 * dg.continuous_operator().generator_per_s()
        .iter().map(|value| value * value).sum::<f64>());
    assert!(omega > 0.0);
    let bare = dg.harmonic_response(
        omega, &[(INLET, C64::from_re(1.0))], WORK, 1e-10, &mut || false,
    );
    if let Ok(bare) = bare {
        // At a floating-point approximation to the pole LU may return a huge
        // impedance instead of an exact-zero-pivot refusal.
        assert!(bare.ports[0].pressure_pa.abs() > 1e8 * z);
    }
    let factor = dg.scattering_factor(
        omega, &[(INLET, z)], WORK, 1e-10, &mut || false,
    ).unwrap();
    let response = factor.response(&[(INLET, C64::from_re(1.0))], &mut || false).unwrap();
    assert!((response.waves[0].reflected_pressure_pa - C64::from_re(1.0)).abs() < 1e-9);
    assert_eq!(response.physical.energy_coordinates.len(), dg.total_dofs());
    assert!(response.physical.energy_coordinates.iter().all(|z| z.re.is_finite() && z.im.is_finite()));
    assert_eq!(response.physical.ports[1].inflow_m3_s, C64::ZERO);
    assert_eq!(response.physical.wall_dissipation_w, 0.0);
    assert_eq!(response.physical.molecular_dissipation_w, 0.0);
    assert_eq!(response.physical.matched_absorption_w, 0.0);
    // Flow vanishes at this pole; the independent wave identity is checked
    // against the supplied amplitude, not a relative tolerance on zero flow.
    let scale = response.waves[0].incident_power_w;
    assert!(response.physical.power_residual_w.abs() < 1e-9 * scale);
    assert!(response.waves[0].power_residual_w.abs() < 1e-9 * scale);
}

#[test]
fn one_factor_retains_coherent_drives_reference_order_and_signed_port_work() {
    let medium = FluidAcousticState::try_new(1.2, 343.0).unwrap();
    let dg = duct(medium, 3, 3);
    let z = impedance(medium);
    let references = [(OUTLET, 0.7 * z), (INLET, 1.4 * z)];
    let factor = dg.scattering_factor(
        1.3 * medium.sound_speed_m_s() / LENGTH, &references, WORK, 1e-10, &mut || false,
    ).unwrap();
    assert_eq!(factor.reference_impedances(), &references);
    assert_eq!(factor.operator_identity(), dg.operator_identity());
    let a = C64::new(1.0, 0.2);
    let b = C64::new(0.05, 0.1);
    let left = factor.response(&[(INLET, a)], &mut || false).unwrap();
    let right = factor.response(&[(OUTLET, b)], &mut || false).unwrap();
    let coherent = factor.response(&[(INLET, a), (OUTLET, b)], &mut || false).unwrap();
    let state_scale = coherent.physical.energy_coordinates.iter()
        .map(|x| x.abs()).fold(0.0_f64, f64::max);
    for ((x, y), z) in left.physical.energy_coordinates.iter()
        .zip(&right.physical.energy_coordinates).zip(&coherent.physical.energy_coordinates)
    {
        assert!((*x + *y - *z).abs() < 1e-10 * state_scale);
    }
    for k in 0..2 {
        assert!((left.waves[k].reflected_pressure_pa + right.waves[k].reflected_pressure_pa
            - coherent.waves[k].reflected_pressure_pa).abs() < 1e-10);
    }
    assert_eq!(coherent.waves[0].id, OUTLET);
    assert_eq!(coherent.waves[1].id, INLET);
    assert_eq!(coherent.physical.ports[0].id, INLET);
    assert_eq!(coherent.physical.ports[1].id, OUTLET);
    assert!(coherent.waves[0].net_input_power_w < 0.0);
    assert!(coherent.waves[1].net_input_power_w > 0.0);
    power_identity(&coherent);
    let repeated = factor.response(&[(OUTLET, b), (INLET, a)], &mut || false).unwrap();
    assert_eq!(repeated, coherent);
}

#[test]
fn unselected_original_ports_keep_zero_flow_instead_of_becoming_absorbers() {
    let medium = FluidAcousticState::try_new(1.2, 343.0).unwrap();
    let dg = duct(medium, 4, 4);
    let phase = 0.9;
    let factor = dg.scattering_factor(
        phase * medium.sound_speed_m_s() / LENGTH,
        &[(INLET, impedance(medium))], WORK, 1e-10, &mut || false,
    ).unwrap();
    let response = factor.response(&[(INLET, C64::from_re(1.0))], &mut || false).unwrap();
    assert_eq!(response.physical.ports[1].inflow_m3_s, C64::ZERO);
    assert!((response.waves[0].reflected_pressure_pa - cis(2.0 * phase)).abs() < 2e-5);
    assert!(factor.response(&[(OUTLET, C64::from_re(1.0))], &mut || false).is_err());
    power_identity(&response);
}

#[test]
fn original_viscothermal_and_maxwell_states_supply_the_absorbed_power() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let dg = DgDuct::new(
        &geometry(), &gas, DgMesh { elements: 1, order: 2 },
        DgBoundary::Port(INLET), DgBoundary::Port(OUTLET),
    ).unwrap()
        .with_viscothermal(16, [50.0, 2000.0], 0.01).unwrap()
        .with_compliant_wall_relaxation(
            WallPin { surface_density: 0.2, stiffness_per_area: 8e7, resistance: 30.0 },
            &GeneralizedMaxwell::new(2e9, vec![(8e8, 0.002)]).unwrap(),
        ).unwrap();
    let z = gas.density * gas.sound_speed / (core::f64::consts::PI * RADIUS * RADIUS);
    let factor = dg.scattering_factor(
        3000.0, &[(INLET, z), (OUTLET, z)], WORK, 1e-10, &mut || false,
    ).unwrap();
    let response = factor.response(&[(INLET, C64::new(1.0, 0.3))], &mut || false).unwrap();
    assert_eq!(response.physical.energy_coordinates.len(), dg.total_dofs());
    assert!(dg.total_dofs() > 2 * dg.field_dofs());
    assert!(response.physical.mean_wall_energy_j > 0.0);
    assert!(response.physical.wall_dissipation_w > 0.0);
    assert!(response.physical.input_power_w > 0.0);
    assert_eq!(response.physical.matched_absorption_w, 0.0);
    let loss = response.physical.wall_dissipation_w + response.physical.molecular_dissipation_w;
    assert!((loss - response.physical.input_power_w).abs() < 2e-8 * loss);
    power_identity(&response);
}

#[test]
fn invalid_sources_loaded_overflow_and_cancellation_preserve_a_reusable_factor() {
    let medium = FluidAcousticState::try_new(1.2, 343.0).unwrap();
    let dg = duct(medium, 1, 2);
    let z = impedance(medium);
    for references in [
        vec![], vec![(PortId(99), z)], vec![(INLET, z), (INLET, z)],
        vec![(INLET, 0.0)], vec![(INLET, -1.0)], vec![(INLET, f64::NAN)],
        vec![(INLET, f64::INFINITY)], vec![(INLET, f64::MIN_POSITIVE)],
    ] {
        assert!(dg.scattering_factor(1000.0, &references, WORK, 1e-10, &mut || false).is_err());
    }
    assert!(dg.scattering_factor(1000.0, &[(INLET, z)], 0, 1e-10, &mut || false).is_err());
    assert!(dg.scattering_factor(0.0, &[(INLET, z)], WORK, 1e-10, &mut || false).is_err());
    assert!(dg.scattering_factor(1000.0, &[(INLET, z)], WORK, 0.0, &mut || false).is_err());
    assert!(matches!(
        dg.scattering_factor(1000.0, &[(INLET, z)], WORK, 1e-10, &mut || true),
        Err(DuctError::Cancelled { .. })
    ));
    let factor = dg.scattering_factor(
        1000.0, &[(INLET, z), (OUTLET, z)], WORK, 1e-10, &mut || false,
    ).unwrap();
    let input = [(INLET, C64::from_re(1.0))];
    let before = factor.response(&input, &mut || false).unwrap();
    for invalid in [
        vec![(INLET, C64::new(f64::NAN, 0.0))],
        vec![(PortId(99), C64::ZERO)],
        vec![(INLET, C64::ZERO), (INLET, C64::ZERO)],
    ] {
        assert!(factor.response(&invalid, &mut || false).is_err());
    }
    let mut calls = 0;
    assert!(matches!(
        factor.response(&input, &mut || { calls += 1; calls == 4 }),
        Err(DuctError::Cancelled { .. })
    ));
    assert_eq!(factor.response(&input, &mut || false).unwrap(), before);
    let silent = factor.response(&[], &mut || false).unwrap();
    assert!(silent.physical.energy_coordinates.iter().all(|&x| x == C64::ZERO));
}
