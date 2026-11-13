//! Original wall loading, physical work and independent line/geometry equations.
//! Manufactured local Maxwell walls are implementation fixtures, not measured
//! material qualification or a claim of nonlocal shell dynamics.

use super::*;
use crate::WallLaw;
use crate::molecular::local::absorption::MolecularCylinder;
use crate::network::field::{NetworkFieldBudget, NetworkFieldProbe, NetworkFieldSection};
use std::sync::atomic::{AtomicUsize, Ordering};

struct Surface {
    stiffness_scale: f64,
    calls: AtomicUsize,
}

impl Surface {
    fn new(stiffness_scale: f64) -> Self {
        Self {
            stiffness_scale,
            calls: AtomicUsize::new(0),
        }
    }

    // Thin local hoop stiffness, one retained Maxwell arm, mass and damper.
    fn coefficients(&self, radius: f64) -> (f64, f64, f64, f64, f64) {
        let geometry = self.stiffness_scale * (0.002 / radius).powi(2);
        (0.08, 8e5 * geometry, 4e5 * geometry, 4e-4, 80.0)
    }

    fn impedance(&self, radius: f64, omega: f64) -> C64 {
        let (mass, equilibrium, arm, tau, resistance) = self.coefficients(radius);
        let s = C64::new(0.0, -omega);
        C64::from_re(resistance) + s.scale(mass)
            + C64::from_re(equilibrium) / s
            + C64::from_re(arm * tau) / (C64::ONE + s.scale(tau))
    }
}

impl WallLaw for Surface {
    fn specific_impedance(&self, radius: f64, omega: f64) -> Result<C64, DuctError> {
        self.calls.fetch_add(1, Ordering::Relaxed);
        if !(100.0..=6000.0).contains(&omega) || !(0.0005..=0.01).contains(&radius) {
            return Err(bad("synthetic Maxwell wall source domain"));
        }
        Ok(self.impedance(radius, omega))
    }
}

fn cp(state: &RelaxationState, omega: f64) -> C64 {
    let s = C64::new(0.0, -omega);
    state.processes().iter().fold(
        C64::from_re(state.frozen_cv_j_kg_k() + state.gas().specific_gas_constant),
        |sum, arm| {
            sum + C64::from_re(arm.heat_capacity_j_kg_k)
                / (C64::ONE + s.scale(arm.relaxation_time_s))
        },
    )
}

// Independent telegraph coefficients under exp(-i omega t). The original
// radial spectral measure is data; no production wave/composition/loss map is
// used here. The wall is an additional reciprocal *mechanical* displacement.
fn gas_pair(point: MolecularCylinder<'_>) -> (C64, C64) {
    let gas = point.state.gas();
    let area = core::f64::consts::PI * point.radius_m * point.radius_m;
    let s = C64::new(0.0, -point.omega);
    let cp = cp(point.state, point.omega);
    match point.wall {
        Some(wall) => {
            let series = wall.response_at(point.omega).unwrap().0.conj().scale(gas.density / area);
            let temperature = wall.shunt_terms().iter().fold(C64::ZERO, |sum, term| {
                let weight = term.gain_per_s / (term.rate_per_s * (gas.gamma - 1.0));
                sum + C64::from_re(weight)
                    / (cp + C64::new(0.0, gas.specific_heat_cp * term.rate_per_s / point.omega))
            });
            let compressibility = (C64::ONE - temperature.scale(gas.specific_gas_constant))
                .scale(1.0 / (gas.density * gas.specific_gas_constant * gas.temperature));
            (series, (s * compressibility).scale(area))
        }
        None => (
            s.scale(gas.density / area),
            (s * (cp - C64::from_re(gas.specific_gas_constant)) / cp)
                .scale(area / gas.pressure),
        ),
    }
}

fn bare_pair(state: &RelaxationState, radius: f64, omega: f64) -> (C64, C64) {
    let gas = state.gas();
    let area = core::f64::consts::PI * radius * radius;
    let s = C64::new(0.0, -omega);
    let cp = cp(state, omega);
    (
        s.scale(gas.density / area),
        (s * (cp - C64::from_re(gas.specific_gas_constant)) / cp)
            .scale(area / gas.pressure),
    )
}

fn telegraph(z: C64, y: C64, length: f64) -> [C64; 4] {
    let k = (-z * y).sqrt();
    let (sin, cos) = sin_cos(k.scale(length));
    [cos, z * sin / k, y * sin / k, cos]
}

fn probe(id: &str, x: f64) -> NetworkFieldProbe {
    NetworkFieldProbe {
        edge_id: id.into(),
        section: NetworkFieldSection::Cylinder { segment: 0, offset_m: x },
    }
}

fn request<'a>(
    omega: f64,
    sources: &'a [PortSource],
    loads: &'a [PortTermination],
    probes: &'a [NetworkFieldProbe],
) -> MolecularNetworkFieldRequest<'a> {
    MolecularNetworkFieldRequest {
        omega,
        source_clock: "molecular-mechanical-wall-clock",
        sources,
        source_bank: None,
        loads,
        coupled_load: None,
        probes,
        budget: NetworkFieldBudget::default(),
    }
}

#[track_caller]
fn relative(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-30),
        "{a:?} != {b:?}"
    );
}

#[track_caller]
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-290),
        "{a:e} != {b:e}"
    );
}

#[test]
fn source_owned_wall_closes_native_molecular_fields_and_separate_physical_work() {
    let state = law(gas(293.15), 2, 1.0);
    let (radius, length) = (0.002, 0.08);
    let net = DuctNetwork {
        nodes: vec![port("drive"), port("load")],
        edges: vec![edge("tube", 0, 1, radius, length)],
    };
    for diffusion in [false, true] {
        let surface = Arc::new(Surface::new(1.0));
        let wall: MolecularCompliantWall = surface.clone();
        let original = MolecularNetwork::new_with_compliant_walls(
            net.clone(),
            records(&net, &[state.clone()]),
            bindings(&net, &[state.clone()], diffusion.then_some(MolecularWallResolution {
                poles_per_field: 8,
                band_rad_s: [100.0, 6000.0],
                maximum_relative_error: 0.02,
            })),
            vec![("tube".into(), wall.clone())],
            64,
        ).unwrap();
        assert!(std::ptr::eq(
            original.edge("tube").unwrap().compliant_wall().unwrap(),
            wall.as_ref(),
        ));
        assert_eq!(surface.calls.load(Ordering::Relaxed), 0);
        let reference = state.gas().characteristic_impedance
            / (core::f64::consts::PI * radius * radius);
        let source_z = C64::new(0.25, -0.1).scale(reference);
        let load_z = C64::new(0.8, 0.2).scale(reference);
        let amplitude = C64::new(1.3, -0.4);
        let sources = [PortSource {
            label: "drive".into(),
            excitation: PortExcitation::Pressure {
                open_circuit_pressure_pa: amplitude,
                impedance_pa_s_m3: source_z,
            },
        }];
        let loads = [PortTermination {
            label: "load".into(),
            load: NetworkLoad::Impedance(load_z),
        }];
        let probes = [probe("tube", 0.0), probe("tube", 0.031), probe("tube", length)];
        for omega in [400.0, 1500.0, 4000.0] {
            let calls = surface.calls.load(Ordering::Relaxed);
            let fields = original.prepare_fields(
                &request(omega, &sources, &loads, &probes), || false,
            ).unwrap();
            assert_eq!(surface.calls.load(Ordering::Relaxed), calls + 1);
            let snapshot = fields.cylinder("tube", 0).unwrap();
            let (z, gas_y) = gas_pair(snapshot);
            let wall_z = surface.impedance(radius, omega);
            let wall_y = wall_z.recip().scale(core::f64::consts::TAU * radius);
            let y = gas_y + wall_y;
            let [a, b, c, d] = telegraph(z, y, length);
            let outlet_flow = amplitude / (a * load_z + b + source_z * (c * load_z + d));
            let at = |x: f64| {
                let [a, b, c, d] = telegraph(z, y, length - x);
                [(a * load_z + b) * outlet_flow, (c * load_z + d) * outlet_flow]
            };
            let actual = fields.fields().sample(&[amplitude], || false).unwrap();
            for (index, x) in [0.0, 0.031, length].into_iter().enumerate() {
                let [p, u] = at(x);
                relative(actual[index].pressure_pa, p, 3e-9);
                relative(actual[index].volume_flow_m3_s, u, 3e-9);
            }
            let losses = fields.prepare_cylinder_absorption(
                &[[0, 2], [0, 1], [1, 2]], 3, || false,
            ).unwrap();
            let observed = losses.sample(&[amplitude], || false).unwrap();
            let (_, _, arm, tau, resistance) = surface.coefficients(radius);
            // Dissipation in the actual dashpot and Maxwell arm, evaluated
            // from the wall velocity rather than a pressure-work residual.
            let mechanical_resistance = resistance + arm * tau / (1.0 + (omega * tau).powi(2));
            let mut integrated = [0.0; 2];
            // Simpson converges as panels^-4 onto the native integral; at 4000 rad/s 512 panels
            // are 3.1e-8 short of it and 2048 1.2e-10, so 4096 keeps the reference well inside 3e-8.
            let panels = 4096;
            for j in 0..=panels {
                let [p, u] = at(length * j as f64 / panels as f64);
                let weight = if j == 0 || j == panels { 1.0 } else if j % 2 == 0 { 2.0 } else { 4.0 };
                integrated[0] += weight * 0.5 * (z.re * u.norm_sq() + gas_y.re * p.norm_sq());
                integrated[1] += weight * 0.5 * core::f64::consts::TAU * radius
                    * mechanical_resistance * (p / wall_z).norm_sq();
            }
            for value in &mut integrated {
                *value *= length / (3 * panels) as f64;
            }
            close(observed[0].absorbed_power_w, integrated[0], 3e-8);
            close(observed[0].wall_work_w, integrated[1], 3e-8);
            close(
                observed[0].total_removed_power_w(),
                actual[0].mean_power_w - actual[2].mean_power_w,
                3e-9,
            );
            close(observed[0].wall_work_w, observed[1].wall_work_w + observed[2].wall_work_w, 3e-10);
            close(observed[0].absorbed_power_w, observed[1].absorbed_power_w + observed[2].absorbed_power_w, 3e-10);
            assert!(observed[0].wall_work_w > 0.0 && observed[0].molecular_power_w > 0.0);
            assert_eq!(surface.calls.load(Ordering::Relaxed), calls + 1);
        }
    }
}

#[test]
fn unlike_source_owned_walls_preserve_branch_reciprocity_and_coherent_work() {
    let states = [law(gas(293.15), 2, 1.0), law(gas(343.15), 2, 0.2)];
    let radii = [0.002, 0.003];
    let lengths = [0.08, 0.12];
    let net = DuctNetwork {
        nodes: vec![port("left"), port("right"), port("receiver")],
        edges: vec![
            edge("left-arm", 0, 2, radii[0], lengths[0]),
            edge("right-arm", 1, 2, radii[1], lengths[1]),
        ],
    };
    let surfaces = [Arc::new(Surface::new(1.0)), Arc::new(Surface::new(1.8))];
    let walls: Vec<(String, MolecularCompliantWall)> = vec![
        ("right-arm".into(), surfaces[1].clone()),
        ("left-arm".into(), surfaces[0].clone()),
    ];
    let original = MolecularNetwork::new_with_compliant_walls(
        net.clone(), records(&net, &states), bindings(&net, &states, None), walls, 64,
    ).unwrap();
    let omega = 1100.0;
    let q = [C64::new(1.0, 0.2).scale(1e-8), C64::new(-0.3, 0.8).scale(0.7e-8)];
    let sources = ["left", "right"].into_iter().zip(q).map(|(label, value)| PortSource {
        label: label.into(), excitation: PortExcitation::VolumeFlow(value),
    }).collect::<Vec<_>>();
    let load_z = C64::new(3e6, -7e5);
    let loads = [PortTermination {
        label: "receiver".into(), load: NetworkLoad::Impedance(load_z),
    }];
    let probes = [
        probe("left-arm", 0.0), probe("left-arm", lengths[0]),
        probe("right-arm", 0.0), probe("right-arm", lengths[1]),
    ];
    let fields = original.prepare_fields(
        &request(omega, &sources, &loads, &probes), || false,
    ).unwrap();
    let actual = fields.fields().sample(&q, || false).unwrap();
    let matrices: Vec<_> = (0..2).map(|i| {
        let (z, gas_y) = bare_pair(&states[i], radii[i], omega);
        let wall_y = surfaces[i].impedance(radii[i], omega).recip()
            .scale(core::f64::consts::TAU * radii[i]);
        telegraph(z, gas_y + wall_y, lengths[i])
    }).collect();
    // The common receiver pressure follows continuity and the finite load;
    // each specified source is the actual inward edge flow.
    let receiver = (q[0] / matrices[0][3] + q[1] / matrices[1][3])
        / (load_z.recip() + matrices[0][2] / matrices[0][3]
            + matrices[1][2] / matrices[1][3]);
    for i in 0..2 {
        let [a, b, c, d] = matrices[i];
        let outlet = (q[i] - c * receiver) / d;
        relative(actual[2 * i].pressure_pa, a * receiver + b * outlet, 3e-9);
        relative(actual[2 * i].volume_flow_m3_s, q[i], 3e-9);
        relative(actual[2 * i + 1].pressure_pa, receiver, 3e-9);
        relative(actual[2 * i + 1].volume_flow_m3_s, outlet, 3e-9);
    }
    let pressure = fields.fields().pressure_transfer();
    relative(pressure[1], pressure[4], 3e-10); // Z_left,right = Z_right,left
    let cells = fields.prepare_cylinder_absorption(&[[0, 1], [2, 3]], 2, || false).unwrap();
    let coherent = cells.sample(&q, || false).unwrap();
    let factor = cells.observe_coherence(None, &[&q], || false).unwrap();
    let left = [q[0], C64::ZERO];
    let right = [C64::ZERO, q[1]];
    let independent = cells.observe_coherence(None, &[&left, &right], || false).unwrap();
    let a = cells.sample(&left, || false).unwrap();
    let b = cells.sample(&right, || false).unwrap();
    let mut interference = 0.0;
    let mut scale = 0.0;
    for i in 0..2 {
        close(coherent[i].wall_work_w, factor[i].wall_work_w, 3e-13);
        close(independent[i].wall_work_w, a[i].wall_work_w + b[i].wall_work_w, 3e-13);
        close(
            coherent[i].total_removed_power_w(),
            actual[2 * i].mean_power_w - actual[2 * i + 1].mean_power_w,
            3e-9,
        );
        interference += (coherent[i].wall_work_w - independent[i].wall_work_w).abs();
        scale += coherent[i].wall_work_w + independent[i].wall_work_w;
        assert_eq!(surfaces[i].calls.load(Ordering::Relaxed), 1);
    }
    assert!(interference > 1e-4 * scale);
}

// Independent RK4 integration of the physical variable-area p/U equations.
// No original duct chain, spherical-wave split or Magnus integrator is used.
#[allow(clippy::too_many_arguments)]
fn variable_reference(
    state: &RelaxationState,
    surface: &Surface,
    profile: bool,
    radii: [f64; 2],
    length: f64,
    omega: f64,
    load: C64,
    at: f64,
) -> [C64; 2] {
    let gas = state.gas();
    let reference = gas.characteristic_impedance
        / (core::f64::consts::PI * radii[0] * radii[0]);
    let capacity = cp(state, omega);
    let beta = (capacity - C64::from_re(gas.specific_gas_constant))
        / capacity / C64::from_re(gas.pressure);
    let s = C64::new(0.0, -omega);
    let growth = (radii[1] / radii[0]).ln() / length;
    let rhs = |x: f64, q: [C64; 2]| {
        let (radius, slope) = if profile {
            let radius = radii[0] * (growth * x).exp();
            (radius, growth * radius)
        } else {
            let slope = (radii[1] - radii[0]) / length;
            (radii[0] + slope * x, slope)
        };
        let area = core::f64::consts::PI * radius * radius;
        let z = s.scale(gas.density / area);
        let wall = surface.impedance(radius, omega).recip()
            .scale(core::f64::consts::TAU * radius * (1.0 + slope * slope).sqrt());
        let y = (s * beta).scale(area) + wall;
        [(-z).scale(1.0 / reference) * q[1], (-y).scale(reference) * q[0]]
    };
    let mut q = [load.scale(1.0 / reference), C64::ONE];
    let panels = 2048;
    let h = (at - length) / panels as f64;
    for i in 0..panels {
        let x = length + h * i as f64;
        let a = rhs(x, q);
        let b = rhs(x + h * 0.5, [q[0] + a[0].scale(h * 0.5), q[1] + a[1].scale(h * 0.5)]);
        let c = rhs(x + h * 0.5, [q[0] + b[0].scale(h * 0.5), q[1] + b[1].scale(h * 0.5)]);
        let d = rhs(x + h, [q[0] + c[0].scale(h), q[1] + c[1].scale(h)]);
        for j in 0..2 {
            q[j] = q[j] + (a[j] + b[j].scale(2.0) + c[j].scale(2.0) + d[j]).scale(h / 6.0);
        }
    }
    [q[0].scale(reference), q[1]]
}

#[test]
fn original_composed_waves_and_slants_reach_cone_and_profile_fields() {
    let state = law(gas(313.15), 2, 0.7);
    let (radii, length, omega) = ([0.002, 0.008], 0.04, 1100.0);
    let amplitude = C64::new(1.0, -0.2).scale(1e-8);
    let load = C64::new(2e6, -4e5);
    let sources = [PortSource {
        label: "drive".into(), excitation: PortExcitation::VolumeFlow(amplitude),
    }];
    let loads = [PortTermination {
        label: "load".into(), load: NetworkLoad::Impedance(load),
    }];
    for profile in [false, true] {
        let surface = Arc::new(Surface::new(1.0));
        let middle = if profile { 7.0 / 16.0 } else { 0.431 };
        let positions = [0.0, middle * length, length];
        let reference_in = variable_reference(
            &state, &surface, profile, radii, length, omega, load, 0.0,
        );
        let scale = amplitude / reference_in[1];
        let expected: Vec<_> = positions.iter().map(|&x| {
            let q = variable_reference(&state, &surface, profile, radii, length, omega, load, x);
            [q[0] * scale, q[1] * scale]
        }).collect();
        let mut errors = Vec::new();
        for n in [16, 32] {
            let segment = if profile {
                Segment::Profile {
                    // exponential() takes the area flare, r(x) = r0 exp(area_flare x / 2), so the
                    // reference's radius growth ln(r1/r0)/L is half of it and the bore ends at r1.
                    profile: Arc::new(crate::RadiusProfile::exponential(
                        radii[0], 2.0 * (radii[1] / radii[0]).ln() / length, length,
                    ).unwrap()),
                    steps: crate::profile::RiccatiSteps(n),
                }
            } else {
                Segment::Cone { inlet_radius: radii[0], outlet_radius: radii[1], length }
            };
            let net = DuctNetwork {
                nodes: vec![port("drive"), port("load")],
                edges: vec![NetworkEdge {
                    id: "taper".into(), from: 0, to: 1,
                    duct: Duct { segments: vec![segment] },
                    corrections: JunctionCorrections::Disabled,
                }],
            };
            let mut edges = bindings(&net, &[state.clone()], None);
            edges[0].cones = ConeSlices::Fixed(n);
            let original = MolecularNetwork::new_with_compliant_walls(
                net.clone(), records(&net, &[state.clone()]), edges,
                vec![("taper".into(), surface.clone())], 128,
            ).unwrap();
            let probes: Vec<_> = (0..3).map(|i| NetworkFieldProbe {
                edge_id: "taper".into(),
                section: if profile {
                    NetworkFieldSection::Profile {
                        segment: 0, boundary: [0, 7 * n / 16, n][i],
                    }
                } else {
                    NetworkFieldSection::Cone { segment: 0, offset_m: positions[i] }
                },
            }).collect();
            let fields = original.prepare_fields(
                &request(omega, &sources, &loads, &probes), || false,
            ).unwrap();
            let actual = fields.fields().sample(&[amplitude], || false).unwrap();
            let mut error = 0.0_f64;
            for i in 0..3 {
                error = error.max((actual[i].pressure_pa - expected[i][0]).abs()
                    / expected[i][0].abs());
                error = error.max((actual[i].volume_flow_m3_s - expected[i][1]).abs()
                    / expected[i][1].abs());
            }
            errors.push(error);
            assert!(fields.prepare_cylinder_absorption(&[[0, 2]], 1, || false).is_err());
        }
        assert!(errors[1] < 0.7 * errors[0] && errors[1] < 2e-3, "{errors:?}");
    }
}

#[test]
fn mechanical_bindings_keep_rigid_bits_original_pins_and_atomic_refusals() {
    struct Rigid;
    impl WallLaw for Rigid {
        fn specific_impedance(&self, _: f64, _: f64) -> Result<C64, DuctError> {
            Err(bad("an exact rigid wall has no finite impedance"))
        }
        fn specific_admittance(&self, _: f64, _: f64) -> Result<C64, DuctError> {
            Ok(C64::ZERO)
        }
    }
    let state = law(gas(293.15), 2, 1.0);
    let net = DuctNetwork {
        nodes: vec![port("drive"), port("load")],
        edges: vec![edge("tube", 0, 1, 0.002, 0.08)],
    };
    let make = |media: NetworkMedia, walls: Vec<(String, MolecularCompliantWall)>| {
        MolecularNetwork::new_with_compliant_walls(
            net.clone(), media, bindings(&net, &[state.clone()], None), walls, 64,
        )
    };
    let plain = model(net.clone(), &[state.clone()], None);
    let rigid = make(records(&net, &[state.clone()]), vec![("tube".into(), Arc::new(Rigid))]).unwrap();
    let amplitude = C64::new(1.0, 0.2).scale(1e-8);
    let sources = [PortSource {
        label: "drive".into(), excitation: PortExcitation::VolumeFlow(amplitude),
    }];
    let loads = [PortTermination {
        label: "load".into(), load: NetworkLoad::Impedance(C64::new(2e6, -4e5)),
    }];
    let probes = [probe("tube", 0.0), probe("tube", 0.04), probe("tube", 0.08)];
    let mut request = request(1100.0, &sources, &loads, &probes);
    let old = plain.prepare_fields(&request, || false).unwrap();
    let zero = rigid.prepare_fields(&request, || false).unwrap();
    assert_eq!(old.fields().pressure_transfer(), zero.fields().pressure_transfer());
    assert_eq!(old.fields().flow_transfer(), zero.fields().flow_transfer());
    let old_cells = old.prepare_cylinder_absorption(&[[0, 2]], 1, || false).unwrap();
    let zero_cells = zero.prepare_cylinder_absorption(&[[0, 2]], 1, || false).unwrap();
    let old_loss = old_cells.sample(&[amplitude], || false).unwrap();
    let zero_loss = zero_cells.sample(&[amplitude], || false).unwrap();
    assert_eq!(old_loss, zero_loss);
    assert_eq!(zero_loss[0].wall_work_w, 0.0);

    let pin = fs_phs::WallPin {
        surface_density: 0.08, stiffness_per_area: 8e5, resistance: 80.0,
    };
    let media = records(&net, &[state.clone()])
        .with_edge_walls(&net, vec![("tube".into(), pin)]).unwrap();
    let pinned = MolecularNetwork::new(
        net.clone(), media.clone(), bindings(&net, &[state.clone()], None), 64,
    ).unwrap();
    assert_eq!(pinned.media().edge_wall("tube"), Some(&pin));
    let explicit = make(
        records(&net, &[state.clone()]), vec![("tube".into(), Arc::new(pin))],
    ).unwrap();
    let native = pinned.prepare_fields(&request, || false).unwrap();
    let declared = explicit.prepare_fields(&request, || false).unwrap();
    assert_eq!(native.fields().pressure_transfer(), declared.fields().pressure_transfer());
    assert_eq!(native.fields().flow_transfer(), declared.fields().flow_transfer());
    assert_ne!(native.fields().pressure_transfer(), old.fields().pressure_transfer());

    let source = Arc::new(Surface::new(1.0));
    let wall: MolecularCompliantWall = source.clone();
    assert!(make(media, vec![("tube".into(), wall.clone())]).is_err());
    assert!(make(
        records(&net, &[state.clone()]), vec![("missing".into(), wall.clone())],
    ).is_err());
    assert!(make(
        records(&net, &[state.clone()]),
        vec![("tube".into(), wall.clone()), ("tube".into(), wall.clone())],
    ).is_err());
    assert_eq!(source.calls.load(Ordering::Relaxed), 0);
    assert!(matches!(
        MolecularNetwork::new_with_compliant_walls_controlled(
            net.clone(), records(&net, &[state.clone()]), bindings(&net, &[state.clone()], None),
            vec![("tube".into(), wall.clone())], 64, &mut || true,
        ),
        Err(MolecularDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    let owned = make(
        records(&net, &[state.clone()]), vec![("tube".into(), wall)],
    ).unwrap();
    let accepted = owned.prepare_fields(&request, || false).unwrap();
    assert!(matches!(
        owned.prepare_fields(&request, || true),
        Err(MolecularDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    request.omega = 7000.0;
    assert!(matches!(
        owned.prepare_fields(&request, || false),
        Err(MolecularDuctError::Duct(DuctError::BadParameter {
            what: "synthetic Maxwell wall source domain"
        }))
    ));
    request.omega = 1100.0;
    let replay = owned.prepare_fields(&request, || false).unwrap();
    assert_eq!(accepted.fields().pressure_transfer(), replay.fields().pressure_transfer());
    assert_eq!(accepted.fields().flow_transfer(), replay.fields().flow_transfer());
}
