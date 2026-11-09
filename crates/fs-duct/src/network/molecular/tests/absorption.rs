//! Independent constitutive work and coherent finite-volume absorption.

use super::*;
use crate::molecular::local::absorption::MolecularCylinder;
use crate::network::field::{NetworkFieldBudget, NetworkFieldProbe, NetworkFieldSection};

fn tube(state: &RelaxationState, radius: f64, length: f64, wall: bool) -> MolecularNetwork {
    model(
        DuctNetwork {
            nodes: vec![port("left"), port("right")],
            edges: vec![edge("tube", 0, 1, radius, length)],
        },
        &[state.clone()],
        wall.then_some(MolecularWallResolution {
            poles_per_field: 8,
            band_rad_s: [0.1, 7000.0],
            maximum_relative_error: 0.02,
        }),
    )
}

fn native(network: &MolecularNetwork, omega: f64) -> MolecularCylinder<'_> {
    network
        .edge("tube")
        .unwrap()
        .prepare_spatial_controlled(omega, &mut || false)
        .unwrap()
        .cylinder(0)
        .unwrap()
}

fn points(offsets: &[f64]) -> Vec<NetworkFieldProbe> {
    offsets
        .iter()
        .map(|&offset_m| NetworkFieldProbe {
            edge_id: "tube".into(),
            section: NetworkFieldSection::Cylinder {
                segment: 0,
                offset_m,
            },
        })
        .collect()
}

fn request<'a>(
    omega: f64,
    sources: &'a [PortSource],
    probes: &'a [NetworkFieldProbe],
) -> MolecularNetworkFieldRequest<'a> {
    MolecularNetworkFieldRequest {
        omega,
        source_clock: "original-molecular-loss-clock",
        sources,
        source_bank: None,
        loads: &[],
        coupled_load: None,
        probes,
        budget: NetworkFieldBudget::default(),
    }
}

fn matched(z: C64, a: C64, b: C64) -> [PortSource; 2] {
    ["left", "right"]
        .into_iter()
        .zip([a, b])
        .map(|(label, amplitude)| PortSource {
            label: label.into(),
            excitation: PortExcitation::Pressure {
                open_circuit_pressure_pa: amplitude.scale(2.0),
                impedance_pa_s_m3: z,
            },
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
}

fn phase(k: C64, x: f64) -> C64 {
    C64::new((k.re * x).cos(), (k.re * x).sin()).scale((-k.im * x).exp())
}

fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-290),
        "{a:e} != {b:e}"
    );
}

// Derive work from actual thermal-mode and reservoir-temperature amplitudes.
// This reference does not use Cp.im to assign molecular loss, the new loss
// helper, wave-integral weights, or boundary-work subtraction.
fn temperature_work(point: MolecularCylinder<'_>) -> [f64; 2] {
    let gas = point.state.gas();
    let area = core::f64::consts::PI * point.radius_m * point.radius_m;
    let mut result = [0.0; 2];
    let modes: Vec<_> = point.wall.map_or_else(
        || vec![(1.0, 0.0)],
        |wall| {
            wall.shunt_terms()
                .iter()
                .map(|term| {
                    (
                        term.gain_per_s / (term.rate_per_s * (gas.gamma - 1.0)),
                        gas.specific_heat_cp * term.rate_per_s,
                    )
                })
                .collect()
        },
    );
    for (weight, diffusion) in modes {
        let temperature = (point.cp + C64::new(0.0, diffusion / point.omega))
            .recip()
            .scale(weight.sqrt() / gas.density);
        result[0] += area * gas.density / gas.temperature * diffusion * temperature.norm_sq();
        for reservoir in point.state.processes() {
            let x = point.omega * reservoir.relaxation_time_s;
            // T-T_j with T_j=T/(1-i omega tau), evaluated without cancellation.
            let slip = temperature * C64::new(0.0, -x) / C64::new(1.0, -x);
            result[1] += area * gas.density / gas.temperature
                * (reservoir.heat_capacity_j_kg_k / reservoir.relaxation_time_s)
                * slip.norm_sq();
        }
    }
    result
}

#[test]
fn original_positive_fourier_and_reservoir_forms_match_native_shunt_and_temperature_work() {
    for (count, rate_scale) in [(0, 1.0), (2, 1e-7), (2, 1.0), (2, 1e7)] {
        let state = law(gas(293.15), count, rate_scale);
        for wall in [false, true] {
            let network = tube(&state, 0.001, 0.09, wall);
            for omega in [1.0, 1200.0, 6500.0] {
                let point = native(&network, omega);
                let weights = point.loss_weights().unwrap();
                let independent = temperature_work(point);
                close(weights[1], independent[0], 4e-11);
                close(weights[2], independent[1], 4e-11);
                let wave = point.wave;
                close(
                    weights[1] + weights[2],
                    (wave.wavenumber / wave.characteristic_impedance).im,
                    4e-11,
                );
                if wall {
                    let model = point.wall.unwrap();
                    let area = core::f64::consts::PI * point.radius_m * point.radius_m;
                    let resistance = model.series_drag_per_s()
                        + model
                            .series_terms()
                            .iter()
                            .map(|term| {
                                term.gain_per_s * omega * omega
                                    / (omega * omega + term.rate_per_s * term.rate_per_s)
                            })
                            .sum::<f64>();
                    close(
                        weights[0],
                        resistance * point.state.gas().density / area,
                        2e-13,
                    );
                } else {
                    assert_eq!(weights[0], 0.0);
                    assert_eq!(weights[1], 0.0);
                }
                if count == 0 {
                    assert_eq!(weights[2], 0.0);
                }
                assert!(weights.into_iter().all(|v| v >= 0.0));
            }
        }
    }
}

#[test]
fn molecular_cells_close_physical_work_partition_and_keep_cubic_pressure_node_loss() {
    let state = law(gas(293.15), 2, 1.0);
    let (radius, length, omega) = (0.001, 0.09, 5107.0);
    let network = tube(&state, radius, length, true);
    let original = native(&network, omega);
    let z =
        original.wave.specific_impedance / C64::from_re(core::f64::consts::PI * radius * radius);
    let a = C64::new(1.3, -0.4);
    let sources = matched(z, a, -a); // exact pressure node at the centre.
    let centre = length * 0.5;
    let h = 0.0002;
    let probes = points(&[
        0.0,
        0.02,
        centre,
        0.07,
        length,
        centre - h,
        centre + h,
        centre - h * 0.5,
        centre + h * 0.5,
    ]);
    let owner = network
        .prepare_fields(&request(omega, &sources, &probes), || false)
        .unwrap();
    let pairs = [
        [0, 4],
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 4],
        [5, 6],
        [7, 8],
        [2, 2],
    ];
    let absorption = owner
        .prepare_cylinder_absorption(&pairs, pairs.len(), || false)
        .unwrap();
    assert!(core::ptr::eq(absorption.fields(), &owner));
    let excitation = owner.fields().port_response().default_excitations();
    let losses = absorption.sample(excitation, || false).unwrap();
    let fields = owner.fields().sample(excitation, || false).unwrap();
    for ((pair, cell), loss) in pairs.iter().zip(absorption.cells()).zip(&losses) {
        let [left, right] = *pair;
        close(
            cell.volume_m3,
            core::f64::consts::PI * radius * radius * (cell.end_m - cell.start_m),
            2e-15,
        );
        close(
            loss.absorbed_power_w,
            fields[left].mean_power_w - fields[right].mean_power_w,
            3e-9,
        );
        close(
            loss.absorbed_power_w,
            loss.viscous_power_w + loss.thermal_power_w + loss.molecular_power_w,
            3e-15,
        );
    }
    let components: [fn(&MolecularNetworkCylinderAbsorption) -> f64; 3] = [
        |v: &MolecularNetworkCylinderAbsorption| v.viscous_power_w,
        |v: &MolecularNetworkCylinderAbsorption| v.thermal_power_w,
        |v: &MolecularNetworkCylinderAbsorption| v.molecular_power_w,
    ];
    for component in components {
        close(
            component(&losses[0]),
            losses[1..5].iter().map(component).sum(),
            5e-12,
        );
    }
    assert!(fields[2].pressure_pa.abs() < 1e-12);
    assert!(losses[6].thermal_power_w > 0.0 && losses[6].molecular_power_w > 0.0);
    close(
        losses[5].thermal_power_w / losses[6].thermal_power_w,
        8.0,
        3e-4,
    );
    close(
        losses[5].molecular_power_w / losses[6].molecular_power_w,
        8.0,
        3e-4,
    );
    close(
        losses[5].viscous_power_w / losses[6].viscous_power_w,
        2.0,
        3e-4,
    );
    assert_eq!(losses[7], MolecularNetworkCylinderAbsorption::default());
    // The generic equilibrium route remains unavailable; only the retained
    // molecular law can authorize these three-channel absorption forms.
    assert!(
        owner
            .fields()
            .prepare_cylinder_absorption(&[[0, 4]], 1, || false)
            .is_err()
    );
}

#[test]
fn molecular_absorption_preserves_covariance_phase_and_independent_interior_quadrature() {
    let state = law(gas(293.15), 2, 1.0);
    let (radius, length, omega) = (0.001, 0.09, 5107.0);
    let network = tube(&state, radius, length, true);
    let point = native(&network, omega);
    let area = core::f64::consts::PI * radius * radius;
    let z = point.wave.specific_impedance / C64::from_re(area);
    let sources = matched(z, C64::ZERO, C64::ZERO);
    let probes = points(&[0.0, length]);
    let owner = network
        .prepare_fields(&request(omega, &sources, &probes), || false)
        .unwrap();
    let absorption = owner
        .prepare_cylinder_absorption(&[[0, 1]], 1, || false)
        .unwrap();
    let (a, b) = (C64::new(1.2, -0.3), C64::new(0.4, 0.6));
    let coherent = [a.scale(2.0), b.scale(2.0)];
    let left = [coherent[0], C64::ZERO];
    let right = [C64::ZERO, coherent[1]];
    let combined = absorption.sample(&coherent, || false).unwrap()[0];
    let factor = absorption
        .observe_coherence(None, &[&coherent], || false)
        .unwrap()[0];
    let independent = absorption
        .observe_coherence(None, &[&left, &right], || false)
        .unwrap()[0];
    close(combined.absorbed_power_w, factor.absorbed_power_w, 2e-15);
    assert!(
        (combined.molecular_power_w - independent.molecular_power_w).abs()
            > 1e-3 * combined.molecular_power_w
    );
    let independent_weights = temperature_work(point);
    let viscous = point.loss_weights().unwrap()[0];
    let mut integrals = [0.0; 3];
    let panels = 1024;
    for i in 0..=panels {
        let x = length * (i as f64 / panels as f64);
        let f = a * phase(point.wave.wavenumber, x);
        let b = b * phase(point.wave.wavenumber, length - x);
        let p = f + b;
        let u = (f - b) / z;
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for (sum, value) in integrals.iter_mut().zip([
            viscous * u.norm_sq(),
            independent_weights[0] * p.norm_sq(),
            independent_weights[1] * p.norm_sq(),
        ]) {
            *sum += weight * value;
        }
    }
    for (actual, integral) in [
        combined.viscous_power_w,
        combined.thermal_power_w,
        combined.molecular_power_w,
    ]
    .into_iter()
    .zip(integrals)
    {
        close(actual, integral * length / (6.0 * panels as f64), 3e-11);
    }
    let mixed = absorption
        .observe_coherence(Some(&left), &[&right], || false)
        .unwrap()[0];
    close(mixed.absorbed_power_w, independent.absorbed_power_w, 2e-15);
}

#[test]
fn bare_molecular_opaque_absorption_and_empty_reservoir_limits_remain_source_bound() {
    for reservoirs in [0, 2] {
        let state = law(gas(293.15), reservoirs, 1.0);
        let (radius, omega) = (0.008, 17_000.0);
        let k = state.relaxing_wavenumber(omega).unwrap();
        let length = if reservoirs == 0 { 0.27 } else { 600.0 / k.im };
        let network = tube(&state, radius, length, false);
        let area = core::f64::consts::PI * radius * radius;
        let z = C64::from_re(omega * state.gas().density / area) / k;
        let amplitude = C64::new(1.3, -0.4);
        let sources = matched(z, amplitude, C64::ZERO);
        let probes = points(&[0.0, length * 0.5, length]);
        let owner = network
            .prepare_fields(&request(omega, &sources, &probes), || false)
            .unwrap();
        let absorption = owner
            .prepare_cylinder_absorption(&[[0, 2], [0, 1], [1, 2]], 3, || false)
            .unwrap();
        let result = absorption
            .sample(owner.fields().port_response().default_excitations(), || {
                false
            })
            .unwrap();
        assert_eq!(owner.network().edge("tube").unwrap().state(), &state);
        for value in &result {
            assert_eq!(value.viscous_power_w, 0.0);
            assert_eq!(value.thermal_power_w, 0.0);
        }
        if reservoirs == 0 {
            assert_eq!(result[0], MolecularNetworkCylinderAbsorption::default());
        } else {
            close(
                result[0].molecular_power_w,
                0.5 * z.recip().re * amplitude.norm_sq(),
                5e-12,
            );
            close(
                result[0].molecular_power_w,
                result[1].molecular_power_w + result[2].molecular_power_w,
                5e-12,
            );
            assert!(
                result[2].molecular_power_w > 0.0,
                "opaque decay must retain weak positive local loss"
            );
        }
    }
}

#[test]
fn molecular_absorption_refuses_non_cylinders_wrong_cells_bad_sources_and_cancellation() {
    let state = law(gas(293.15), 2, 1.0);
    let network = tube(&state, 0.001, 0.09, true);
    let sources = matched(C64::from_re(1e8), C64::ONE, C64::ZERO);
    let probes = points(&[0.0, 0.09]);
    let owner = network
        .prepare_fields(&request(5107.0, &sources, &probes), || false)
        .unwrap();
    assert!(
        owner
            .prepare_cylinder_absorption(&[[0, 1]], 0, || false)
            .is_err()
    );
    assert!(
        owner
            .prepare_cylinder_absorption(&[[1, 0]], 1, || false)
            .is_err()
    );
    assert!(
        owner
            .prepare_cylinder_absorption(&[[0, 2]], 1, || false)
            .is_err()
    );
    assert!(matches!(
        owner.prepare_cylinder_absorption(&[[0, 1]], 1, || true),
        Err(DuctError::Cancelled { .. })
    ));
    let prepared = owner
        .prepare_cylinder_absorption(&[[0, 1]], 1, || false)
        .unwrap();
    assert!(prepared.sample(&[C64::ONE], || false).is_err());
    assert!(
        prepared
            .sample(&[C64::new(f64::NAN, 0.0), C64::ZERO], || false)
            .is_err()
    );
    let values = owner.fields().port_response().default_excitations();
    let first = prepared.sample(values, || false).unwrap();
    let mut polls = 0;
    assert!(matches!(
        prepared.sample(values, || {
            polls += 1;
            polls == 2
        }),
        Err(DuctError::Cancelled { .. })
    ));
    assert_eq!(first, prepared.sample(values, || false).unwrap());
    let mut graph = network.network().clone();
    graph.edges[0].duct.segments = vec![Segment::Cone {
        inlet_radius: 0.001,
        outlet_radius: 0.0015,
        length: 0.09,
    }];
    let cone = model(graph, &[state], None);
    let probes = [
        NetworkFieldProbe {
            edge_id: "tube".into(),
            section: NetworkFieldSection::Inlet { segment: 0 },
        },
        NetworkFieldProbe {
            edge_id: "tube".into(),
            section: NetworkFieldSection::Outlet { segment: 0 },
        },
    ];
    let owner = cone
        .prepare_fields(&request(1700.0, &sources, &probes), || false)
        .unwrap();
    assert!(
        owner
            .prepare_cylinder_absorption(&[[0, 1]], 1, || false)
            .is_err()
    );
}
