//! `NetworkLoad::Cavity`: a closed cavity's lumped compliance `C = V / (rho c^2)`, `Z = i / (omega C)`
//! (exp(-i omega t)). The dimensions are test geometry: a 16.1 cm3 cavity behind a 0.45 cm neck of 0.076 cm2, the
//! size of a maxillary sinus behind its ostium, not a specimen claim.

use super::*;
use crate::network::fluid::drive_sources;
use fs_material::fluid_acoustics::FluidAcousticState;

const VOLUME_M3: f64 = 16.1e-6;
const NECK_RADIUS_M: f64 = 1.555_4e-3;
const NECK_LENGTH_M: f64 = 4.5e-3;

fn compliance_load(gas: &GasState, omega: f64) -> C64 {
    // The oracle, from the volume and the gas alone: C = V / (rho c^2), Z = 1 / (-i omega C).
    let c = VOLUME_M3 / (gas.density * gas.sound_speed * gas.sound_speed);
    C64::new(0.0, 1.0 / (omega * c))
}

fn behind_neck(load: NetworkLoad) -> DuctNetwork {
    DuctNetwork {
        nodes: vec![port("input"), node("cavity", NodeBoundary::Termination(load))],
        edges: vec![edge("neck", 0, 1, NECK_RADIUS_M, NECK_LENGTH_M)],
    }
}

#[test]
fn cavity_is_its_volume_compliance_behind_a_neck() {
    let gas = gas();
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        // Below, near and above the neck-and-cavity resonance (559.5 Hz for these sizes in this dry air, lumped).
        for (rung, f) in [53.0, 211.0, 537.0, 913.0, 2_917.0].into_iter().enumerate() {
            let start = Instant::now();
            let w = core::f64::consts::TAU * f;
            let net = behind_neck(NetworkLoad::Cavity { volume_m3: VOLUME_M3 });
            let response = network_port_impedance(&net, &gas, w, loss).unwrap();
            let oracle = cylinder_input(&gas, NECK_RADIUS_M, NECK_LENGTH_M, w, loss, compliance_load(&gas, w));
            let error = relative(response.z[0], oracle);
            row("network-cavity-behind-neck", &net, &response, error, 1.0e-12, rung, start);
            assert!(error <= 1.0e-12, "cavity behind a neck, {loss:?} at {f} Hz: {error:.3e}");
            // The same cavity as a labelled port's termination takes the same law.
            let line = DuctNetwork {
                nodes: vec![port("input"), port("cavity")],
                edges: vec![edge("neck", 0, 1, NECK_RADIUS_M, NECK_LENGTH_M)],
            };
            let driven = |load: NetworkLoad| {
                network_drive(&line, &gas, w, loss, "input", &[PortTermination { label: "cavity".into(), load }])
                    .unwrap()
            };
            let as_port = driven(NetworkLoad::Cavity { volume_m3: VOLUME_M3 });
            let supplied = driven(NetworkLoad::Impedance(compliance_load(&gas, w)));
            let port_error = relative(as_port.impedance, supplied.impedance);
            assert!(port_error <= 1.0e-14, "cavity as a port load, {loss:?} at {f} Hz: {port_error:.3e}");
        }
    }
}

#[test]
fn cavity_is_an_equal_volume_closed_tube_at_low_ka() {
    // A closed tube of the same volume is a compliance to first order: its input impedance is the cavity's times
    // (kL)/tan(kL). Both sit behind the same neck, so the two inputs differ by that tube's own O((kL)^2), and a sign
    // error in the cavity's law would put them apart by a factor near -1.
    let gas = gas();
    let (tube_radius, loss) = (0.012, LossModel::Lossless);
    let tube_length = VOLUME_M3 / (core::f64::consts::PI * tube_radius * tube_radius);
    for f in [5.0, 20.0, 80.0] {
        let w = core::f64::consts::TAU * f;
        let kl = w / gas.sound_speed * tube_length;
        let cavity = network_port_impedance(&behind_neck(NetworkLoad::Cavity { volume_m3: VOLUME_M3 }), &gas, w, loss)
            .unwrap()
            .z[0];
        let tube = DuctNetwork {
            nodes: vec![
                port("input"),
                node("ostium", NodeBoundary::Junction),
                node("end", NodeBoundary::Termination(NetworkLoad::Closed)),
            ],
            edges: vec![
                edge("neck", 0, 1, NECK_RADIUS_M, NECK_LENGTH_M),
                edge("tube", 1, 2, tube_radius, tube_length),
            ],
        };
        let closed = network_port_impedance(&tube, &gas, w, loss).unwrap().z[0];
        let error = relative(cavity, closed);
        eprintln!(
            "{{\"suite\":\"fs-duct\",\"case\":\"network-cavity-closed-tube\",\"hz\":{f},\"kl\":{kl:.4e},\"rel_err\":{error:.4e},\"bound\":{:.4e}}}",
            kl * kl
        );
        assert!(cavity.im > 0.0 && closed.im > 0.0, "a compliance at {f} Hz: {cavity:?} and {closed:?}");
        assert!(error <= kl * kl, "cavity against an equal-volume closed tube at {f} Hz: {error:.3e} > (kL)^2");
    }
}

#[test]
fn cavity_refuses_what_it_cannot_hold() {
    let gas = gas();
    let w = core::f64::consts::TAU * 211.0;
    for volume_m3 in [0.0, -VOLUME_M3, f64::NAN, f64::INFINITY] {
        assert!(
            network_port_impedance(&behind_neck(NetworkLoad::Cavity { volume_m3 }), &gas, w, LossModel::Lossless)
                .is_err(),
            "volume {volume_m3} admitted"
        );
    }
}

#[test]
fn cavity_on_a_shared_fluid_chain_is_its_bulk_modulus_compliance() {
    // The generic-fluid route takes the chain's own adiabatic bulk modulus. Explicit analytic coefficients, not a
    // measured-liquid claim.
    let state = FluidAcousticState::try_new(1000.0, 1500.0).unwrap();
    let line = DuctNetwork {
        nodes: vec![port("input"), port("cavity")],
        edges: vec![NetworkEdge {
            id: "neck".into(),
            from: 0,
            to: 1,
            duct: duct(NECK_RADIUS_M, NECK_LENGTH_M),
            corrections: JunctionCorrections::Disabled,
        }],
    };
    let zc = state.characteristic_impedance_pa_s_m() / (core::f64::consts::PI * NECK_RADIUS_M * NECK_RADIUS_M);
    for f in [53.0, 911.0, 4_001.0] {
        let w = core::f64::consts::TAU * f;
        let source = [PortSource {
            label: "input".into(),
            excitation: PortExcitation::Pressure {
                open_circuit_pressure_pa: C64::ONE,
                impedance_pa_s_m3: C64::from_re(zc),
            },
        }];
        let driven = |load: NetworkLoad| {
            drive_sources(&line, &state, w, &source, None, &[PortTermination { label: "cavity".into(), load }], None)
                .unwrap()
        };
        let cavity = driven(NetworkLoad::Cavity { volume_m3: VOLUME_M3 });
        let supplied = driven(NetworkLoad::Impedance(C64::new(
            0.0,
            state.adiabatic_bulk_modulus_pa() / (w * VOLUME_M3),
        )));
        for (a, b) in cavity.fields.port_pressures.iter().chain(&cavity.fields.port_inflows).zip(
            supplied.fields.port_pressures.iter().chain(&supplied.fields.port_inflows),
        ) {
            let error = relative(*a, *b);
            assert!(error <= 1.0e-14, "fluid-chain cavity at {f} Hz: {error:.3e}");
        }
    }
}
