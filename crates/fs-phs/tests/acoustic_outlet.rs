//! G1/G3: physical outlet flow, coordinate orientation and complete storage.
use fs_math::c64::C64;
use fs_phs::acoustic_outlet::acoustic_chain_with_outlet_flow;
use fs_phs::harmonic::{HarmonicBudget, HarmonicModel};
use fs_phs::{
    AcousticSection, ChainOutlet, MouthFlange, ViscothermalPin, acoustic_chain_outlet_wall,
};

fn section(cells: usize, cone: bool) -> AcousticSection {
    AcousticSection {
        length: 0.3,
        radius: 0.01,
        outlet_radius: if cone { 0.02 } else { 0.01 },
        cells,
    }
}

fn observed_sections(cells: usize, tapered: bool) -> Vec<AcousticSection> {
    if tapered {
        // Two distinct slopes exercise the original frustum kernel, whose
        // coordinates retain physical compression volume throughout the bore.
        vec![
            AcousticSection {
                length: 0.1,
                radius: 0.01,
                outlet_radius: 0.015,
                cells: 1,
            },
            AcousticSection {
                length: 0.2,
                radius: 0.015,
                outlet_radius: 0.02,
                cells: cells - 1,
            },
        ]
    } else {
        vec![section(cells, false)]
    }
}

#[test]
fn outward_flow_closes_independent_compression_volume_balance_for_odd_and_even_chains() {
    for cells in [2, 3, 4, 5] {
        for outlet in [
            ChainOutlet::PressureRelease,
            ChainOutlet::Radiating(MouthFlange::Flanged),
        ] {
            let (system, flow) = acoustic_chain_with_outlet_flow(
                &[section(cells, false)],
                1.2,
                340.0,
                outlet,
                1,
                &[],
                None,
                None,
            )
            .unwrap();
            let mut initial = vec![0.0; system.state_dim()];
            for i in 0..cells {
                initial[2 * i] = (i as f64 + 1.0) * 1e-8;
                initial[2 * i + 1] = (2 * i + 1) as f64 * 1e-4;
            }
            let input = 2e-4;
            let dt = 1e-5;
            let step = fs_phs::step(&system, &initial, &[input], dt).unwrap();
            let midpoint: Vec<_> = initial
                .iter()
                .zip(&step.x)
                .map(|(a, b)| f64::midpoint(*a, *b))
                .collect();
            let volume_change: f64 = (0..cells)
                .map(|i| {
                    let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
                    sign * (step.x[2 * i] - initial[2 * i])
                })
                .sum();
            let boundary_change = dt * (input - flow.volume_flow_m3_s(&midpoint).unwrap());
            assert!(
                (volume_change - boundary_change).abs() < 1e-11 * boundary_change.abs().max(1e-12),
                "cells={cells}, {outlet:?}: stored={volume_change:e}, boundary={boundary_change:e}"
            );
        }
    }
}

#[test]
fn cylinder_and_frustum_preserve_physical_low_frequency_through_flow() {
    let pin = ViscothermalPin {
        dynamic_viscosity: 1.8e-5,
        gamma: 1.4,
        prandtl: 0.71,
        foster_branches: 2,
    };
    for cells in [3, 4] {
        for tapered in [false, true] {
            for loss in [None, Some(&pin)] {
                let (system, flow) = acoustic_chain_with_outlet_flow(
                    &observed_sections(cells, tapered),
                    1.2,
                    340.0,
                    ChainOutlet::Radiating(MouthFlange::Flanged),
                    1,
                    &[],
                    loss,
                    None,
                )
                .unwrap();
                let model = HarmonicModel::from_ode(
                    &system,
                    HarmonicBudget {
                        max_dimension: 256,
                        max_dense_bytes: 32_000_000,
                        max_dense_work: 100_000_000,
                        relative_residual_tolerance: 1e-10,
                    },
                    || false,
                )
                .unwrap();
                let response = model
                    .factor(0.01, || false)
                    .unwrap()
                    .response(&[C64::ONE], || false)
                    .unwrap();
                let outward = flow
                    .coefficients()
                    .iter()
                    .zip(&response.state)
                    .fold(C64::ZERO, |sum, (a, x)| sum + x.scale(*a));
                // Independently, every nonleaking duct transports its inlet
                // flow to the outlet as omega -> 0, regardless of taper or
                // the arbitrary alternating orientation of its native states.
                assert!(
                    (outward - C64::ONE).abs() < 2e-5,
                    "cells={cells}, tapered={tapered}, loss={}: {outward:?}",
                    loss.is_some()
                );
                let real: Vec<_> = response.state.iter().map(|x| x.re).collect();
                assert!((flow.volume_flow_m3_s(&real).unwrap() - outward.re).abs() < 1e-12);
            }
        }
    }
}

#[test]
fn observed_terminal_flow_reproduces_independent_physical_radiation_work() {
    let density = 1.2;
    let speed = 340.0;
    let dt = 1e-5;
    for cells in [3, 4] {
        for tapered in [false, true] {
            let sections = observed_sections(cells, tapered);
            let (system, flow) = acoustic_chain_with_outlet_flow(
                &sections,
                density,
                speed,
                ChainOutlet::Radiating(MouthFlange::Flanged),
                1,
                &[],
                None,
                None,
            )
            .unwrap();
            let radius = sections.last().unwrap().outlet_radius;
            let length: f64 = sections.iter().map(|section| section.length).sum();
            // Independently, the compact flanged load at the original pin
            // frequency dissipates R * U_out^2. The accepted step's original
            // material resistance account must equal this physical work.
            let omega = core::f64::consts::PI * speed / (2.0 * length);
            let area = core::f64::consts::PI * radius * radius;
            let resistance = density * speed / area * (omega * radius / speed).powi(2) / 2.0;
            let mut initial = vec![0.0; system.state_dim()];
            for i in 0..cells {
                initial[2 * i] = (i as f64 + 1.0) * 1e-8;
                initial[2 * i + 1] = (2 * i + 1) as f64 * 1e-4;
            }
            let step = fs_phs::step(&system, &initial, &[2e-4], dt).unwrap();
            let midpoint: Vec<_> = initial
                .iter()
                .zip(&step.x)
                .map(|(a, b)| f64::midpoint(*a, *b))
                .collect();
            let outward = flow.volume_flow_m3_s(&midpoint).unwrap();
            let expected = dt * resistance * outward * outward;
            assert!(expected > 0.0);
            assert!(
                (step.dissipated - expected).abs() < 1e-11 * expected,
                "cells={cells}, tapered={tapered}: native={}, physical={expected}",
                step.dissipated
            );
        }
    }
}

#[test]
fn spherical_source_remains_available_but_its_approximate_terminal_is_not_a_physical_observer() {
    for cells in [3, 4] {
        for outlet in [
            ChainOutlet::PressureRelease,
            ChainOutlet::Radiating(MouthFlange::Flanged),
        ] {
            let sections = [section(cells, true)];
            let original =
                acoustic_chain_outlet_wall(&sections, 1.2, 340.0, outlet, 1, &[], None, None)
                    .unwrap();
            assert_eq!(original.state_dim(), 2 * cells + 1);
            let step =
                fs_phs::step(&original, &vec![0.0; original.state_dim()], &[2e-4], 1e-5).unwrap();
            assert!(step.x.iter().any(|x| *x != 0.0));
            assert!(
                acoustic_chain_with_outlet_flow(&sections, 1.2, 340.0, outlet, 1, &[], None, None,)
                    .is_err()
            );
        }
    }
}

#[test]
fn rigid_terminal_never_radiates_its_disconnected_flux_and_invalid_histories_refuse() {
    let (system, flow) = acoustic_chain_with_outlet_flow(
        &[section(3, true)],
        1.2,
        340.0,
        ChainOutlet::Rigid,
        1,
        &[],
        None,
        None,
    )
    .unwrap();
    let mut state = vec![1.0; system.state_dim()];
    assert_eq!(flow.volume_flow_m3_s(&state).unwrap(), 0.0);
    assert!(flow.volume_flow_m3_s(&state[1..]).is_err());
    state[0] = f64::NAN;
    assert!(flow.volume_flow_m3_s(&state).is_err());
}
