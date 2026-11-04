//! Manufactured resistive-boundary equation oracles, not measured material claims.

use fs_duct::dg::{DgBoundary, DgDuct, DgFluidMaterial, DgMesh, DgResistance, PortId};
use fs_duct::{Duct, Segment};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, det};

const RADIUS: f64 = 0.01;
const LENGTH: f64 = 0.4;
const WORK: usize = 20_000_000;

fn duct(
    medium: FluidAcousticState,
    inlet: DgBoundary,
    outlet: DgBoundary,
) -> DgDuct<DgFluidMaterial<FluidAcousticState>> {
    DgDuct::from_material(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: LENGTH,
            }],
        },
        &medium,
        DgMesh {
            elements: 4,
            order: 4,
        },
        inlet,
        outlet,
    )
    .unwrap()
}

fn impedance(medium: FluidAcousticState) -> f64 {
    medium.characteristic_impedance_pa_s_m() / (core::f64::consts::PI * RADIUS * RADIUS)
}

#[test]
fn refinement_converges_to_physical_reflection_and_removes_boundary_residual_loss() {
    let medium = FluidAcousticState::try_new(1.2, 343.0).unwrap();
    let zc = impedance(medium);
    let mut errors = Vec::new();
    let mut penalties = Vec::new();
    for mesh in [
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgMesh {
            elements: 4,
            order: 4,
        },
    ] {
        let dg = DgDuct::from_material(
            &Duct {
                segments: vec![Segment::Cylinder {
                    radius: RADIUS,
                    length: LENGTH,
                }],
            },
            &medium,
            mesh,
            DgBoundary::Port(PortId(0)),
            DgBoundary::resistive(4.0 * zc).unwrap(),
        )
        .unwrap();
        let response = dg
            .harmonic_response(
                2.0 * medium.sound_speed_m_s() / LENGTH,
                &[(PortId(0), C64::from_re(1e-7))],
                WORK,
                1e-10,
                &mut || false,
            )
            .unwrap();
        let zin = response.ports[0].pressure_pa.scale(1e7);
        let reflection = (zin - C64::from_re(zc)) / (zin + C64::from_re(zc));
        errors.push((reflection - C64::new(0.6 * det::cos(4.0), 0.6 * det::sin(4.0))).abs());
        let p = *response.pressure_pa.last().unwrap();
        let q = *response.volume_flow_m3_s.last().unwrap();
        let load_power = 2.0 * zc * (p.scale(1.0 / zc) + q).scale(0.2).norm_sq();
        let penalty = (p.scale(0.2) - q.scale(0.8 * zc)).norm_sq() / (2.0 * zc);
        let total = response.matched_absorption_w;
        assert!((total - load_power - penalty).abs() < 1e-12 * total);
        penalties.push(penalty / total);
    }
    assert!(errors[1] < errors[0] * 0.01);
    assert!(penalties[1] < penalties[0] * 0.01);
}

#[test]
fn finite_load_reflection_has_the_correct_amplitude_phase_and_inlet_orientation() {
    for medium in [
        FluidAcousticState::try_new(1.2, 343.0).unwrap(),
        FluidAcousticState::try_new(1000.0, 1480.0).unwrap(),
    ] {
        let zc = impedance(medium);
        for ratio in [0.2, 1.0, 4.0] {
            let boundary = DgBoundary::resistive(ratio * zc).unwrap();
            for (inlet, outlet) in [
                (DgBoundary::Port(PortId(0)), boundary),
                (boundary, DgBoundary::Port(PortId(0))),
            ] {
                let dg = duct(medium, inlet, outlet);
                for phase in [0.8, 2.0] {
                    let response = dg
                        .harmonic_response(
                            phase * medium.sound_speed_m_s() / LENGTH,
                            &[(PortId(0), C64::from_re(1e-7))],
                            WORK,
                            1e-10,
                            &mut || false,
                        )
                        .unwrap();
                    let zin = response.ports[0].pressure_pa.scale(1e7);
                    let reflection = (zin - C64::from_re(zc)) / (zin + C64::from_re(zc));
                    let r = (ratio - 1.0) / (ratio + 1.0);
                    let expected = C64::new(r * det::cos(2.0 * phase), r * det::sin(2.0 * phase));
                    assert!(
                        (reflection - expected).abs() < 2e-5,
                        "ratio={ratio}, phase={phase}: {reflection:?} != {expected:?}"
                    );
                    assert!(response.input_power_w > 0.0);
                    assert!(response.matched_absorption_w > 0.0);
                    assert!(response.power_residual_w.abs() < 2e-8 * response.input_power_w);
                }
            }
        }
    }
}

#[test]
fn resistance_equal_to_the_original_medium_reproduces_the_matched_generator() {
    let medium = FluidAcousticState::try_new(1000.0, 1480.0).unwrap();
    let matched = duct(medium, DgBoundary::Closed, DgBoundary::Matched);
    let resistive = duct(
        medium,
        DgBoundary::Closed,
        DgBoundary::resistive(impedance(medium)).unwrap(),
    );
    assert_eq!(
        matched.continuous_operator().generator_per_s(),
        resistive.continuous_operator().generator_per_s(),
    );
    assert_ne!(matched.operator_identity(), resistive.operator_identity());
    let other = duct(
        medium,
        DgBoundary::Closed,
        DgBoundary::resistive(2.0 * impedance(medium)).unwrap(),
    );
    assert_ne!(other.operator_identity(), resistive.operator_identity());
}

#[test]
fn boundary_absorption_frequency_tangent_uses_the_actual_resistance() {
    let medium = FluidAcousticState::try_new(1.2, 343.0).unwrap();
    let dg = duct(
        medium,
        DgBoundary::Port(PortId(0)),
        DgBoundary::resistive(0.2 * impedance(medium)).unwrap(),
    );
    let omega = 1.3 * medium.sound_speed_m_s() / LENGTH;
    let inputs = [(PortId(0), C64::from_re(1e-7))];
    let tangent = dg
        .harmonic_response_with_frequency_derivative(omega, &inputs, None, WORK, 1e-10, &mut || {
            false
        })
        .unwrap();
    let delta = omega * 1e-4;
    let before = dg
        .harmonic_response(omega - delta, &inputs, WORK, 1e-10, &mut || false)
        .unwrap();
    let after = dg
        .harmonic_response(omega + delta, &inputs, WORK, 1e-10, &mut || false)
        .unwrap();
    let measured = (after.matched_absorption_w - before.matched_absorption_w) / (2.0 * delta);
    let predicted = tangent.derivative.matched_absorption_w;
    assert!(predicted.abs() > 1e-15);
    assert!((measured - predicted).abs() < 2e-6 * predicted.abs());
    assert!(tangent.derivative.power_residual_w.abs() < 2e-8 * predicted.abs());
}

#[test]
fn transient_resistive_end_absorbs_without_producing_energy() {
    let medium = FluidAcousticState::try_new(1.2, 343.0).unwrap();
    let dg = duct(
        medium,
        DgBoundary::Closed,
        DgBoundary::resistive(4.0 * impedance(medium)).unwrap(),
    );
    let propagator = dg
        .propagator(LENGTH / medium.sound_speed_m_s() / 16.0, 1e-12)
        .unwrap();
    let mut state = dg.project_initial_state(
        |x| det::exp(-0.5 * ((x - LENGTH * 0.4) / (LENGTH * 0.1)).powi(2)),
        |_| 0.0,
    );
    let initial = state.stored_energy_j();
    let mut previous = initial;
    for _ in 0..64 {
        let output = propagator.step(&mut state, &[]).unwrap();
        assert!(!output.ledger.balance_checked);
        let energy = state.stored_energy_j();
        assert!(energy <= previous + initial * 1e-11);
        previous = energy;
    }
    assert!(previous < initial * 0.5);
}

#[test]
fn resistance_admission_rejects_active_and_nonfinite_loads() {
    for invalid in [0.0, -0.0, -1.0, f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
        assert!(DgResistance::new(invalid).is_err());
        assert!(DgBoundary::resistive(invalid).is_err());
    }
    for resistance in [f64::MIN_POSITIVE, 1.0, f64::MAX] {
        assert_eq!(
            DgResistance::new(resistance).unwrap().resistance_pa_s_m3(),
            resistance
        );
    }
}
