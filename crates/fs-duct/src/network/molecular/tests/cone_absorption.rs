//! Original spherical fields, independent physical volume work and coherence.
//! Manufactured gas reservoirs and a local wall are implementation fixtures,
//! not an empirical gas/material qualification.

use super::*;
use crate::molecular::local::absorption::MolecularCone;
use crate::network::field::{NetworkFieldBudget, NetworkFieldProbe, NetworkFieldSection};
use std::sync::atomic::{AtomicUsize, Ordering};

struct Surface(AtomicUsize);
impl crate::WallLaw for Surface {
    fn specific_impedance(&self, radius: f64, omega: f64) -> Result<C64, DuctError> {
        self.0.fetch_add(1, Ordering::Relaxed);
        // Passive, locally reacting mass/damper/hoop spring in exp(-i omega t).
        Ok(C64::new(
            80.0,
            3e5 * (0.002 / radius).powi(2) / omega - 0.02 * omega,
        ))
    }
}

fn cone(
    state: &RelaxationState,
    radii: [f64; 2],
    length: f64,
    slices: usize,
    diffusion: bool,
    mechanical: bool,
) -> (MolecularNetwork, Arc<Surface>) {
    let network = DuctNetwork {
        nodes: vec![port("left"), port("right")],
        edges: vec![NetworkEdge {
            id: "cone".into(),
            from: 0,
            to: 1,
            duct: Duct {
                segments: vec![Segment::Cone {
                    inlet_radius: radii[0],
                    outlet_radius: radii[1],
                    length,
                }],
            },
            corrections: JunctionCorrections::Disabled,
        }],
    };
    let media = records(&network, &[state.clone()]);
    let mut edges = bindings(
        &network,
        &[state.clone()],
        // 12 poles: at 8, a 3.2 mm cone end's slice (2.95 mm) is refused at 2 % on this band.
        diffusion.then_some(MolecularWallResolution {
            poles_per_field: 12,
            band_rad_s: [0.1, 7000.0],
            maximum_relative_error: 0.02,
        }),
    );
    edges[0].cones = ConeSlices::Fixed(slices);
    let wall = Arc::new(Surface(AtomicUsize::new(0)));
    let walls = if mechanical {
        vec![("cone".into(), wall.clone() as MolecularCompliantWall)]
    } else {
        Vec::new()
    };
    (
        MolecularNetwork::new_with_compliant_walls(network, media, edges, walls, 4096).unwrap(),
        wall,
    )
}

fn probes(offsets: &[f64]) -> Vec<NetworkFieldProbe> {
    offsets
        .iter()
        .map(|&offset_m| NetworkFieldProbe {
            edge_id: "cone".into(),
            section: NetworkFieldSection::Cone {
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
        source_clock: "original-spherical-work-clock",
        sources,
        source_bank: None,
        loads: &[],
        coupled_load: None,
        probes,
        budget: NetworkFieldBudget::default(),
    }
}
fn controls() -> ConeAbsorptionControls {
    ConeAbsorptionControls {
        relative_tolerance: 1e-11,
        maximum_panels_per_cell: 8192,
        ..ConeAbsorptionControls::default()
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-280),
        "actual={actual:e}, expected={expected:e}"
    );
}
fn phase(k: C64, x: f64) -> C64 {
    C64::new((k.re * x).cos(), (k.re * x).sin()).scale((-k.im * x).exp())
}

// Independent spherical solution and its derivative. A and B are local
// pressure amplitudes at opposite endpoints; no production TwoPort or
// loss/quadrature helper is used in this physical reference.
fn spherical(source: MolecularCone<'_>, x: f64, a: C64, b: C64) -> [C64; 2] {
    let slope = (source.outlet_radius_m - source.inlet_radius_m) / source.length_m;
    let radius = source.inlet_radius_m + slope * x;
    let k = source.original.wave.wavenumber;
    let ik = C64::new(0.0, 1.0) * k;
    let f = a.scale(source.inlet_radius_m / radius) * phase(k, x);
    let b = b.scale(source.outlet_radius_m / radius) * phase(k, source.length_m - x);
    let p = f + b;
    let derivative = ik * (f - b) - p.scale(slope / radius);
    let u = derivative.scale(core::f64::consts::PI * radius * radius)
        / (ik * source.original.wave.specific_impedance);
    [p, u]
}

fn driven(source: MolecularCone<'_>, a: C64, b: C64) -> [PortSource; 2] {
    let z = source.original.wave.characteristic_impedance.re;
    let [p0, u0] = spherical(source, 0.0, a, b);
    let [p1, u1] = spherical(source, source.length_m, a, b);
    [("left", p0 + u0.scale(z)), ("right", p1 - u1.scale(z))].map(|(label, amplitude)| PortSource {
        label: label.into(),
        excitation: PortExcitation::Pressure {
            open_circuit_pressure_pa: amplitude,
            impedance_pa_s_m3: C64::from_re(z),
        },
    })
}

// Direct radial spectral work, before either native wave square root. This
// uses the retained original measure, not a signed flux difference.
fn coefficients(source: MolecularCone<'_>) -> [f64; 4] {
    let original = source.original;
    let gas = original.state.gas();
    let omega = original.omega;
    let area = core::f64::consts::PI * original.radius_m.powi(2);
    let cp = original.state.processes().iter().fold(
        C64::from_re(original.state.frozen_cv_j_kg_k() + gas.specific_gas_constant),
        |sum, reservoir| {
            sum + C64::from_re(reservoir.heat_capacity_j_kg_k)
                / C64::new(1.0, -omega * reservoir.relaxation_time_s)
        },
    );
    let mut weights = [0.0; 4];
    if let Some(wall) = original.wall {
        weights[0] = gas.density / area
            * (wall.series_drag_per_s()
                + wall
                    .series_terms()
                    .iter()
                    .map(|term| {
                        term.gain_per_s * omega * omega
                            / (term.rate_per_s * term.rate_per_s + omega * omega)
                    })
                    .sum::<f64>());
        for term in wall.shunt_terms() {
            let weight = term.gain_per_s / (term.rate_per_s * (gas.gamma - 1.0));
            let diffusion = gas.specific_heat_cp * term.rate_per_s;
            let denominator = cp + C64::new(0.0, diffusion / omega);
            let scale = area / (gas.density * gas.temperature) * weight / denominator.norm_sq();
            weights[1] += scale * diffusion;
            weights[2] += scale * omega * cp.im;
        }
    } else {
        weights[2] = area * gas.specific_gas_constant / gas.pressure * omega * cp.im / cp.norm_sq();
    }
    weights[3] = original.mechanical_admittance_per_metre.re;
    weights
}

fn powers(value: MolecularNetworkCylinderAbsorption) -> [f64; 4] {
    [
        value.viscous_power_w,
        value.thermal_power_w,
        value.molecular_power_w,
        value.wall_work_w,
    ]
}

#[test]
fn spherical_induction_matches_independent_volume_work_in_both_directions() {
    let state = law(gas(293.15), 2, 1.0);
    let (length, omega) = (0.06, 400.0);
    let (a, b) = (C64::new(1.3, -0.4), C64::new(-0.2, 0.6));
    for radii in [[0.0008, 0.0032], [0.0032, 0.0008]] {
        let (network, wall) = cone(&state, radii, length, 1, true, true);
        let native = network
            .edge("cone")
            .unwrap()
            .prepare_spatial_controlled(omega, &mut || false)
            .unwrap()
            .cone(0)
            .unwrap();
        let weights = coefficients(native);
        let sources = driven(native, a, b);
        let points = probes(&[0.0, 0.37 * length, length]);
        wall.0.store(0, Ordering::Relaxed);
        let owner = network
            .prepare_fields(&request(omega, &sources, &points), || false)
            .unwrap();
        assert_eq!(wall.0.load(Ordering::Relaxed), 1);
        let bank = owner.fields().port_response().default_excitations();
        let cells = owner
            .prepare_cone_absorption(&[[0, 2], [0, 1], [1, 2]], controls(), || false)
            .unwrap();
        let result = cells.sample(bank, || false).unwrap();
        let values = owner.fields().sample(bank, || false).unwrap();
        for (&x, actual) in [0.0, 0.37 * length, length].iter().zip(&values) {
            let [p, u] = spherical(native, x, a, b);
            near(actual.pressure_pa, p, p.abs());
            near(actual.volume_flow_m3_s, u, u.abs());
        }
        let mut integral = [0.0; 4];
        let mut without_induction = 0.0;
        let steps = 4096;
        for j in 0..=steps {
            let x = length * j as f64 / steps as f64;
            let radius = radii[0] + (radii[1] - radii[0]) * (x / length);
            let area_ratio = (radius / native.original.radius_m).powi(2);
            let [p, u] = spherical(native, x, a, b);
            let weight = if j == 0 || j == steps {
                1.0
            } else if j % 2 == 0 {
                2.0
            } else {
                4.0
            };
            integral[0] += weight * weights[0] * u.norm_sq() / area_ratio;
            for channel in 1..4 {
                integral[channel] += weight * weights[channel] * p.norm_sq() * area_ratio;
            }
            let f = a.scale(radii[0] / native.original.radius_m)
                * phase(native.original.wave.wavenumber, x);
            let b = b.scale(radii[1] / native.original.radius_m)
                * phase(native.original.wave.wavenumber, length - x);
            without_induction += weight
                * weights[0]
                * ((f - b) / native.original.wave.characteristic_impedance).norm_sq();
        }
        for channel in 0..4 {
            integral[channel] *= 0.5 * length / (3.0 * steps as f64);
            close(powers(result[0])[channel], integral[channel], 3e-8);
            close(
                powers(result[1])[channel] + powers(result[2])[channel],
                powers(result[0])[channel],
                3e-9,
            );
        }
        without_induction *= 0.5 * length / (3.0 * steps as f64);
        assert!(
            (without_induction - result[0].viscous_power_w).abs() > 0.1 * result[0].viscous_power_w
        );
        close(
            result[0].total_removed_power_w(),
            values[0].mean_power_w - values[2].mean_power_w,
            2e-8,
        );
        close(
            cells.cells()[0].volume_m3,
            core::f64::consts::PI
                * length
                * (radii[0].powi(2) + radii[0] * radii[1] + radii[1].powi(2))
                / 3.0,
            1e-14,
        );
        assert!(
            cells
                .quadrature()
                .iter()
                .all(|q| q.panels > 0 && q.relative_change <= controls().relative_tolerance)
        );
        assert_eq!(wall.0.load(Ordering::Relaxed), 1);
    }
}

#[test]
fn native_slice_partition_and_coherence_retain_the_original_material_bank() {
    let state = law(gas(293.15), 2, 1.0);
    let length = 0.12;
    let (network, wall) = cone(&state, [0.0012, 0.0032], length, 4, true, true);
    let sources = [
        PortSource {
            label: "left".into(),
            excitation: PortExcitation::Pressure {
                open_circuit_pressure_pa: C64::new(1.0, 0.2),
                impedance_pa_s_m3: C64::from_re(8e7),
            },
        },
        PortSource {
            label: "right".into(),
            excitation: PortExcitation::Pressure {
                open_circuit_pressure_pa: C64::new(-0.3, 0.5),
                impedance_pa_s_m3: C64::from_re(2e7),
            },
        },
    ];
    let offsets: Vec<_> = (0..=8).map(|j| length * (j as f64 / 8.0)).collect();
    let points = probes(&offsets);
    let owner = network
        .prepare_fields(&request(2600.0, &sources, &points), || false)
        .unwrap();
    assert!(core::ptr::eq(owner.network(), &network));
    let original_transfer = owner.fields().pressure_transfer().to_vec();
    assert_eq!(wall.0.load(Ordering::Relaxed), 4);
    let whole: Vec<_> = (0..4).map(|j| [2 * j, 2 * j + 2]).collect();
    let split: Vec<_> = (0..8).map(|j| [j, j + 1]).collect();
    let whole = owner
        .prepare_cone_absorption(&whole, controls(), || false)
        .unwrap();
    let split = owner
        .prepare_cone_absorption(&split, controls(), || false)
        .unwrap();
    assert!(
        owner
            .prepare_cone_absorption(&[[0, 8]], controls(), || false)
            .is_err()
    );
    let factors = [
        [C64::new(0.7, -0.2), C64::new(0.3, 0.6)],
        [C64::new(-0.1, 0.4), C64::new(0.9, -0.3)],
    ];
    let total = whole
        .observe_coherence(None, &[&factors[0], &factors[1]], || false)
        .unwrap();
    let mut summed = [[0.0; 4]; 4];
    for factor in factors {
        let values = owner.fields().sample(&factor, || false).unwrap();
        let loss = whole.sample(&factor, || false).unwrap();
        let parts = split.sample(&factor, || false).unwrap();
        for (j, cell) in loss.iter().enumerate() {
            for channel in 0..4 {
                summed[j][channel] += powers(*cell)[channel];
                close(
                    powers(parts[2 * j])[channel] + powers(parts[2 * j + 1])[channel],
                    powers(*cell)[channel],
                    4e-9,
                );
            }
            close(
                cell.total_removed_power_w(),
                values[2 * j].mean_power_w - values[2 * j + 2].mean_power_w,
                3e-8,
            );
        }
    }
    for (j, cell) in total.iter().enumerate() {
        for channel in 0..4 {
            close(powers(*cell)[channel], summed[j][channel], 1e-13);
        }
    }
    close(
        whole.cells().iter().map(|cell| cell.volume_m3).sum(),
        core::f64::consts::PI
            * length
            * (0.0012_f64.powi(2) + 0.0012 * 0.0032 + 0.0032_f64.powi(2))
            / 3.0,
        1e-14,
    );
    assert_eq!(owner.fields().pressure_transfer(), original_transfer);
    assert_eq!(wall.0.load(Ordering::Relaxed), 4);
}

#[test]
fn zero_opaque_and_refused_cells_preserve_physical_fields_and_work_bounds() {
    let state = law(gas(293.15), 0, 1.0);
    let (network, _) = cone(&state, [0.001, 0.003], 0.09, 1, false, false);
    let source = network
        .edge("cone")
        .unwrap()
        .prepare_spatial_controlled(3100.0, &mut || false)
        .unwrap()
        .cone(0)
        .unwrap();
    let sources = driven(source, C64::ONE, C64::new(-0.2, 0.6));
    let points = probes(&[0.0, 0.09]);
    let owner = network
        .prepare_fields(&request(3100.0, &sources, &points), || false)
        .unwrap();
    let no_quadrature = ConeAbsorptionControls {
        maximum_quadrature_evaluations: 0,
        ..controls()
    };
    let zero = owner
        .prepare_cone_absorption(&[[0, 1], [1, 1]], no_quadrature, || false)
        .unwrap();
    assert_eq!(
        zero.sample(owner.fields().port_response().default_excitations(), || {
            false
        })
        .unwrap(),
        vec![MolecularNetworkCylinderAbsorption::default(); 2]
    );
    assert!(
        zero.quadrature()
            .iter()
            .all(|q| *q == ConeAbsorptionQuadrature::default())
    );
    assert_eq!(zero.cells()[1].volume_m3, 0.0);

    let state = law(gas(293.15), 2, 1.0);
    let (short, _) = cone(&state, [0.001, 0.003], 1.0, 1, true, false);
    let native = short
        .edge("cone")
        .unwrap()
        .prepare_spatial_controlled(6500.0, &mut || false)
        .unwrap()
        .cone(0)
        .unwrap();
    let length = 100.0 / native.original.wave.wavenumber.im;
    let (network, _) = cone(&state, [0.001, 0.003], length, 1, true, false);
    let native = network
        .edge("cone")
        .unwrap()
        .prepare_spatial_controlled(6500.0, &mut || false)
        .unwrap()
        .cone(0)
        .unwrap();
    let sources = driven(native, C64::new(1.0, 0.2), C64::new(-0.1, 0.3));
    let points = probes(&[0.0, length]);
    let owner = network
        .prepare_fields(&request(6500.0, &sources, &points), || false)
        .unwrap();
    let bank = owner.fields().port_response().default_excitations();
    let fields = owner.fields().sample(bank, || false).unwrap();
    for budget in [
        ConeAbsorptionControls {
            maximum_panels_per_cell: 2,
            ..controls()
        },
        ConeAbsorptionControls {
            maximum_quadrature_evaluations: 1,
            ..controls()
        },
    ] {
        assert!(
            owner
                .prepare_cone_absorption(&[[0, 1]], budget, || false)
                .is_err()
        );
    }
    assert!(
        owner
            .prepare_cone_absorption(&[[1, 0]], controls(), || false)
            .is_err()
    );
    assert!(
        owner
            .prepare_cone_absorption(&[[0, 3]], controls(), || false)
            .is_err()
    );
    assert!(
        owner
            .prepare_cone_absorption(&[[0, 1]], controls(), || true)
            .is_err()
    );
    let mut calls = 0;
    assert!(
        owner
            .prepare_cone_absorption(&[[0, 1]], controls(), || {
                calls += 1;
                calls == 12
            })
            .is_err()
    );
    let accepted = owner
        .prepare_cone_absorption(&[[0, 1], [0, 0]], controls(), || false)
        .unwrap();
    assert!(accepted.sample(&[C64::ONE], || false).is_err());
    let loss = accepted.sample(bank, || false).unwrap();
    assert!(loss[0].viscous_power_w > 0.0 && loss[0].molecular_power_w > 0.0);
    assert_eq!(loss[1], MolecularNetworkCylinderAbsorption::default());
    close(
        loss[0].total_removed_power_w(),
        fields[0].mean_power_w - fields[1].mean_power_w,
        3e-8,
    );
    assert_eq!(owner.fields().sample(bank, || false).unwrap(), fields);
}
