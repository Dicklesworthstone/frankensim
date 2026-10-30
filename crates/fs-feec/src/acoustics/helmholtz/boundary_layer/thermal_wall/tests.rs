use super::super::HelmholtzFe;
use super::super::tests::{
    AIR, GEOMETRY, TRANSPORT, box_mesh, complete_basis, end_faces, quadratic, side_faces,
};
use super::*;

const SOLID: WallThermalTransport = WallThermalTransport {
    conductivity_w_m_k: 0.01,
    volumetric_heat_capacity_j_m3_k: 1000.0,
};

fn wall(depth: ThermalWallDepth) -> ThermalWall {
    ThermalWall::Conductive(ConductiveThermalWall {
        transport: SOLID,
        depth,
        minimum_tangential_scale_m: 0.1,
    })
}
fn norm(z: [f64; 2]) -> f64 {
    z[0].hypot(z[1])
}
fn sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}
fn mul(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}
fn div(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    let d = b[0] * b[0] + b[1] * b[1];
    [
        (a[0] * b[0] + a[1] * b[1]) / d,
        (a[1] * b[0] - a[0] * b[1]) / d,
    ]
}

// Independent P1 weak heat equation through the actual solid depth. The Robin
// boundary is the gas half-space flux; the remote boundary is enforced directly.
// No tanh, slab admittance or production thermal transfer factor appears here.
fn resolved_slab(n: usize, d: f64, f: f64, backing: ThermalBacking) -> ([f64; 2], f64) {
    let omega = 2.0 * core::f64::consts::PI * f;
    let h = d / n as f64;
    let k = SOLID.conductivity_w_m_k;
    let capacity = SOLID.volumetric_heat_capacity_j_m3_k;
    let delta = (2.0 * TRANSPORT.thermal_conductivity_w_m_k
        / (AIR.density_kg_m3 * TRANSPORT.cp_j_kg_k * omega))
        .sqrt();
    let ag = TRANSPORT.thermal_conductivity_w_m_k / delta;
    let count = if backing == ThermalBacking::Isothermal {
        n
    } else {
        n + 1
    };
    let off = [-k / h, omega * capacity * h / 6.0];
    let mut diagonal = vec![[2.0 * k / h, 2.0 * omega * capacity * h / 3.0]; count];
    diagonal[0] = [k / h + ag, omega * capacity * h / 3.0 + ag];
    if backing == ThermalBacking::Adiabatic {
        diagonal[n] = [k / h, omega * capacity * h / 3.0];
    }
    let mut rhs = vec![[0.0, 0.0]; count];
    rhs[0] = [ag, ag];
    for i in 1..count {
        let factor = div(off, diagonal[i - 1]);
        diagonal[i] = sub(diagonal[i], mul(factor, off));
        rhs[i] = sub(rhs[i], mul(factor, rhs[i - 1]));
    }
    let mut temperature = vec![[0.0, 0.0]; n + 1];
    for i in (0..count).rev() {
        let right = if i + 1 < count {
            sub(rhs[i], mul(off, temperature[i + 1]))
        } else {
            rhs[i]
        };
        temperature[i] = div(right, diagonal[i]);
    }
    let solid_loss = temperature
        .windows(2)
        .map(|pair| {
            let jump = sub(pair[1], pair[0]);
            k / h * (jump[0] * jump[0] + jump[1] * jump[1]) / ag
        })
        .sum();
    (temperature[0], solid_loss)
}

#[test]
fn conductive_slab_matches_independent_normal_heat_fe_refinement() {
    let f = 300.0;
    let delta = (SOLID.conductivity_w_m_k
        / (core::f64::consts::PI * f * SOLID.volumetric_heat_capacity_j_m3_k))
        .sqrt();
    let d = 1.3 * delta;
    for backing in [ThermalBacking::Isothermal, ThermalBacking::Adiabatic] {
        let expected = wall(ThermalWallDepth::Slab {
            thickness_m: d,
            backing,
        })
        .response(TRANSPORT, AIR, GEOMETRY, f)
        .unwrap();
        let mut errors = Vec::new();
        for n in [32, 64, 128] {
            let (temperature, loss) = resolved_slab(n, d, f, backing);
            let error = (norm(sub(temperature, expected.wall_temperature))
                / norm(expected.wall_temperature))
            .max((loss / expected.solid_loss_factor - 1.0).abs());
            println!("conductive thermal slab backing={backing:?} cells={n} error={error:e}");
            errors.push(error);
        }
        assert!(errors[2] < 2e-5, "{errors:?}");
        for pair in errors.windows(2) {
            let order = (pair[0] / pair[1]).log2();
            assert!((1.8..=2.2).contains(&order), "{order}: {errors:?}");
        }
    }
}

#[test]
fn conductive_wall_limits_entropy_partition_and_thin_slab_loss_are_passive() {
    let f = 300.0;
    let eg =
        (TRANSPORT.thermal_conductivity_w_m_k * AIR.density_kg_m3 * TRANSPORT.cp_j_kg_k).sqrt();
    let ew = (SOLID.conductivity_w_m_k * SOLID.volumetric_heat_capacity_j_m3_k).sqrt();
    let fraction = ew / (eg + ew);
    let half = wall(ThermalWallDepth::HalfSpace)
        .response(TRANSPORT, AIR, GEOMETRY, f)
        .unwrap();
    assert!((half.gas_temperature_deficit[0] - fraction).abs() < 3e-16);
    assert!(half.gas_temperature_deficit[1].abs() < 1e-16);
    assert!((half.gas_loss_factor - fraction * fraction).abs() < 3e-16);
    assert!((half.solid_loss_factor - fraction * (1.0 - fraction)).abs() < 3e-16);
    let delta = (SOLID.conductivity_w_m_k
        / (core::f64::consts::PI * f * SOLID.volumetric_heat_capacity_j_m3_k))
        .sqrt();
    for backing in [ThermalBacking::Isothermal, ThermalBacking::Adiabatic] {
        let thick = wall(ThermalWallDepth::Slab {
            thickness_m: 30.0 * delta,
            backing,
        })
        .response(TRANSPORT, AIR, GEOMETRY, f)
        .unwrap();
        assert_eq!(thick, half);
        for exponent in -14..=3 {
            let x = 10.0_f64.powi(exponent);
            let r = wall(ThermalWallDepth::Slab {
                thickness_m: x * delta,
                backing,
            })
            .response(TRANSPORT, AIR, GEOMETRY, f)
            .unwrap();
            assert!(r.gas_loss_factor >= 0.0 && r.solid_loss_factor >= 0.0);
            let total = r.gas_loss_factor + r.solid_loss_factor;
            assert!(
                (total - (r.gas_temperature_deficit[0] - r.gas_temperature_deficit[1])).abs()
                    < 2e-15
            );
            assert!(
                norm(sub(
                    [1.0, 0.0],
                    [
                        r.wall_temperature[0] + r.gas_temperature_deficit[0],
                        r.wall_temperature[1] + r.gas_temperature_deficit[1]
                    ]
                )) < 5e-16
            );
            if x < 1e-10 && backing == ThermalBacking::Adiabatic {
                // The leading reactive term cancels in ordinary complex subtraction.
                // The independently positive gas+solid entropy form must survive.
                assert!(total > 0.0);
                assert!(norm(sub(r.wall_temperature, [1.0, 0.0])) < 1e-9);
            }
            if x < 1e-10 && backing == ThermalBacking::Isothermal {
                assert!(norm(sub(r.gas_temperature_deficit, [1.0, 0.0])) < 1e-9);
            }
        }
    }
}

#[test]
fn original_full_fe_response_closes_conductive_thermal_work_and_retains_default_bits() {
    let (complex, positions) = box_mesh(1, [0.5, 0.04, 0.04]);
    let fe = HelmholtzFe::new(&complex, &positions, 2, AIR).unwrap();
    let source = end_faces(&fe, 0.0);
    let outlet = end_faces(&fe, 0.5);
    let faces = side_faces(&fe, 0.5);
    let modes = complete_basis(&fe);
    let ports = [source.as_slice(), outlet.as_slice()];
    let build = || {
        ReducedHelmholtz::new_with_boundary_layers(
            &fe, &modes, &ports, false, &faces, TRANSPORT, GEOMETRY,
        )
        .unwrap()
    };
    let baseline = build();
    let default = build().with_thermal_wall(ThermalWall::Isothermal).unwrap();
    let c = fe.boundary_mass(&outlet).unwrap();
    for f in [200.0, 600.0, 1200.0] {
        let loads = [[0.0, 0.0], [0.7, 0.0]];
        let old = baseline.response(f, 0, 0.01, &loads).unwrap();
        let same = default.response(f, 0, 0.01, &loads).unwrap();
        assert_eq!(
            old.amplitudes
                .iter()
                .flatten()
                .map(|x| x.to_bits())
                .collect::<Vec<_>>(),
            same.amplitudes
                .iter()
                .flatten()
                .map(|x| x.to_bits())
                .collect::<Vec<_>>()
        );
        for thermal in [
            ThermalWall::Adiabatic,
            wall(ThermalWallDepth::HalfSpace),
            wall(ThermalWallDepth::Slab {
                thickness_m: 1e-4,
                backing: ThermalBacking::Isothermal,
            }),
            wall(ThermalWallDepth::Slab {
                thickness_m: 1e-4,
                backing: ThermalBacking::Adiabatic,
            }),
        ] {
            let model = build().with_thermal_wall(thermal).unwrap();
            let result = model.response(f, 0, 0.01, &loads).unwrap();
            let loss = model.boundary_layer_power(f, &result.amplitudes).unwrap();
            let input = 0.5 * result.mean_pressure[0][0] * result.inflow;
            let re: Vec<_> = result.amplitudes.iter().map(|p| p[0]).collect();
            let im: Vec<_> = result.amplitudes.iter().map(|p| p[1]).collect();
            let outgoing = 0.35 / (AIR.density_kg_m3 * AIR.sound_speed_m_s)
                * (quadratic(&c, &re) + quadratic(&c, &im));
            assert!(
                (input - outgoing - loss.viscous_w - loss.thermal_w).abs() < 5e-10 * input.abs()
            );
            assert!(loss.viscous_w > 0.0 && loss.thermal_w >= 0.0);
            assert_ne!(result.amplitudes, old.amplitudes);
            if thermal == ThermalWall::Adiabatic {
                assert_eq!(loss.thermal_w, 0.0);
            }
            assert_eq!(model.thermal_wall(), Some(thermal));
        }
    }
}

#[test]
fn conductive_wall_refuses_inconsistent_transport_and_unresolved_lateral_profiles() {
    let mut config = ConductiveThermalWall {
        transport: SOLID,
        depth: ThermalWallDepth::HalfSpace,
        minimum_tangential_scale_m: 0.1,
    };
    config.transport.conductivity_w_m_k = 0.0;
    assert!(
        ThermalWall::Conductive(config)
            .response(TRANSPORT, AIR, GEOMETRY, 300.0)
            .is_err()
    );
    config.transport = SOLID;
    config.minimum_tangential_scale_m = 1e-6;
    assert!(
        ThermalWall::Conductive(config)
            .response(TRANSPORT, AIR, GEOMETRY, 300.0)
            .is_err()
    );
    config.minimum_tangential_scale_m = 0.1;
    config.depth = ThermalWallDepth::Slab {
        thickness_m: f64::NAN,
        backing: ThermalBacking::Adiabatic,
    };
    assert!(
        ThermalWall::Conductive(config)
            .response(TRANSPORT, AIR, GEOMETRY, 300.0)
            .is_err()
    );
    let (complex, positions) = box_mesh(1, [0.5, 0.04, 0.04]);
    let fe = HelmholtzFe::new(&complex, &positions, 1, AIR).unwrap();
    let source = end_faces(&fe, 0.0);
    assert!(
        ReducedHelmholtz::new(&fe, &complete_basis(&fe), &[&source], false)
            .unwrap()
            .with_thermal_wall(ThermalWall::Adiabatic)
            .is_err()
    );
}
