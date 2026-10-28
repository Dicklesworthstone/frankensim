//! G0/G1: the original spatial wall measure drives the coupled gas law.
use super::*;
use crate::profile::{RadiusProfile, RiccatiSteps};
use fs_la::eigen_complex::lu_complex;
use std::sync::Arc;

const WORK: usize = 1_000_000_000;

fn resolution(points: usize) -> DgWallResolution {
    DgWallResolution {
        poles_per_field: 8,
        quadrature_points: points,
        band_rad_s: [50.0, 2000.0].map(|f| core::f64::consts::TAU * f),
        maximum_relative_error: 0.05,
    }
}

fn curve() -> Duct {
    Duct {
        segments: vec![Segment::Profile {
            profile: Arc::new(RadiusProfile::exponential(0.0007, 40.0, 0.05).unwrap()),
            steps: RiccatiSteps(32),
        }],
    }
}

fn radius(shape: &Duct, x: f64) -> f64 {
    match &shape.segments[0] {
        Segment::Cylinder { radius, .. } => *radius,
        Segment::Cone {
            inlet_radius,
            outlet_radius,
            length,
        } => inlet_radius + (outlet_radius - inlet_radius) * x / length,
        Segment::Profile { profile, .. } => profile.radius(x).unwrap(),
        _ => panic!("this independent oracle uses one uninterrupted bore"),
    }
}

fn base(shape: &Duct, elements: usize, order: usize) -> DgDuct {
    DgDuct::new(
        shape,
        &gas(),
        DgMesh { elements, order },
        DgBoundary::Port(PortId(7)),
        DgBoundary::Closed,
    )
    .unwrap()
}

fn relative(a: &[C64], b: &[C64]) -> f64 {
    let error = a
        .iter()
        .zip(b)
        .map(|(a, b)| (*a - *b).abs())
        .fold(0.0_f64, f64::max);
    error / b.iter().map(|z| z.abs()).fold(0.0_f64, f64::max)
}

// Eliminate ALL actual auxiliary coordinates by an independent dense solve.
// In particular, do not use a diagonal-pole shortcut: molecular/thermal
// auxiliaries are coupled and no longer form independent scalar arms.
fn pressure_schur(duct: &DgDuct, omega: f64) -> Vec<C64> {
    let m = duct.field_dofs;
    let n = duct.total_dofs();
    let start = 2 * m;
    let extra = n - start;
    let s = C64::new(0.0, -omega);
    let mut memory = vec![C64::ZERO; extra * extra];
    for i in 0..extra {
        for j in 0..extra {
            memory[i * extra + j] = C64::from_re(-duct.l_tilde[(start + i) * n + start + j]);
            if i == j {
                memory[i * extra + j] = memory[i * extra + j] + s;
            }
        }
    }
    let factor = lu_complex(&memory, extra).unwrap();
    let mut result = vec![C64::ZERO; m * m];
    for j in 0..m {
        let mut rhs: Vec<_> = (start..n)
            .map(|k| C64::from_re(duct.l_tilde[k * n + j]))
            .collect();
        factor.solve(&mut rhs);
        for i in 0..m {
            let correction = (start..n).fold(C64::ZERO, |sum, k| {
                sum + rhs[k - start].scale(duct.l_tilde[i * n + k])
            });
            result[i * m + j] = C64::from_re(-duct.l_tilde[i * n + j]) - correction;
            if i == j {
                result[i * m + j] = result[i * m + j] + s;
            }
        }
    }
    result
}

// Independent weak heat-equation evaluation from geometry, original mass,
// the shared Bessel owner and material cp(s); no augmented coefficients or
// retained channel projections are consumed here. At finite axial quadrature
// the effective map is I + sum_q (B_inf*kappa_q(s)-1) v_q v_q^T. Replacing I
// by sum_q v_q v_q^T would silently change the original pressure mass.
fn expected_pressure(
    original: &DgDuct,
    shape: &Duct,
    material: &RelaxationState,
    omega: f64,
    points: usize,
) -> Vec<C64> {
    let m = original.field_dofs;
    let local = original.mesh.order + 1;
    let gas = material.gas();
    let bulk = gas.density * material.frozen_sound_speed_m_s().powi(2);
    let equilibrium_bulk = gas.density * gas.sound_speed.powi(2);
    let s = C64::new(0.0, -omega);
    let cp = material.processes().iter().fold(
        C64::from_re(material.frozen_cv_j_kg_k() + gas.specific_gas_constant),
        |sum, p| {
            sum + C64::from_re(p.heat_capacity_j_kg_k) / (C64::ONE + s.scale(p.relaxation_time_s))
        },
    );
    let mut result = vec![C64::ZERO; m * m];
    for i in 0..m {
        result[i * m + i] = s;
    }
    let r = resolution(points);
    let (nodes, weights) = gauss_legendre_rule(points);
    for (e, element) in original.elements.iter().enumerate() {
        let jac = 0.5 * element.length;
        for (&xi, &weight) in nodes.iter().zip(&weights) {
            let x = element.x_left + jac * (xi + 1.0);
            let radius = radius(shape, x);
            let model = fs_phs::ZkDiffusive::new(
                radius,
                gas.density,
                gas.dynamic_viscosity,
                gas.gamma,
                gas.prandtl,
                r.poles_per_field,
                r.band_rad_s,
                r.maximum_relative_error,
            )
            .unwrap();
            let heat = model.shunt_terms().iter().fold(C64::ZERO, |sum, t| {
                let w = t.gain_per_s / (t.rate_per_s * (gas.gamma - 1.0));
                sum + s.scale(w) / (s * cp + C64::from_re(gas.specific_heat_cp * t.rate_per_s))
            });
            let kappa = (C64::ONE - heat.scale(gas.specific_gas_constant))
                .scale(1.0 / (gas.density * gas.specific_gas_constant * gas.temperature));
            let coefficient = s * (kappa.scale(bulk) - C64::ONE);
            let phi = original.lgl.eval_basis(xi);
            let root = det::sqrt(
                jac * weight * core::f64::consts::PI * radius * radius / equilibrium_bulk,
            );
            let projection: Vec<f64> = (0..local)
                .map(|i| {
                    root * (0..local)
                        .map(|j| element.mp_inv_sqrt[i * local + j] * phi[j])
                        .sum::<f64>()
                })
                .collect();
            for i in 0..local {
                for j in 0..local {
                    let index = (e * local + i) * m + e * local + j;
                    result[index] =
                        result[index] + coefficient.scale(projection[i] * projection[j]);
                }
            }
        }
    }
    result
}

#[test]
fn uniform_molecular_spatial_and_exact_mass_routes_agree_in_the_coupled_response() {
    let shape = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.001,
            length: 0.05,
        }],
    };
    let original = base(&shape, 1, 2);
    let material = law(gas(), 2);
    let r = resolution(3);
    let wall = fs_phs::WallPin {
        surface_density: 0.4,
        stiffness_per_area: 8e8,
        resistance: 100.0,
    };
    let uniform = original
        .clone()
        .with_viscothermal(8, [50.0, 2000.0], 0.05)
        .unwrap()
        .with_molecular_relaxation(material.clone())
        .unwrap()
        .with_compliant_wall(wall)
        .unwrap();
    let spatial = original
        .with_spatial_viscothermal(r)
        .unwrap()
        .with_molecular_relaxation(material)
        .unwrap()
        .with_compliant_wall(wall)
        .unwrap();
    for omega in [400.0, 1200.0, 6000.0] {
        let solve = |duct: &DgDuct| {
            duct.harmonic_response(
                omega,
                &[(PortId(7), C64::from_re(1e-9))],
                WORK,
                1e-10,
                &mut || false,
            )
            .unwrap()
        };
        let a = solve(&uniform);
        let b = solve(&spatial);
        assert!(relative(&a.pressure_pa, &b.pressure_pa) < 2e-8);
        assert!(relative(&a.volume_flow_m3_s, &b.volume_flow_m3_s) < 2e-8);
        for (a, b) in [
            (a.mean_stored_energy_j, b.mean_stored_energy_j),
            (a.wall_dissipation_w, b.wall_dissipation_w),
            (a.molecular_dissipation_w, b.molecular_dissipation_w),
        ] {
            assert!(a > 0.0 && b > 0.0 && (a / b - 1.0).abs() < 2e-8, "{a} {b}");
        }
    }
}

#[test]
fn tapered_molecular_operator_matches_independent_local_thermodynamics() {
    let material = law(gas(), 1);
    for shape in [
        curve(),
        Duct {
            segments: vec![Segment::Cone {
                inlet_radius: 0.0007,
                outlet_radius: 0.002,
                length: 0.05,
            }],
        },
    ] {
        let original = base(&shape, 1, 1);
        for points in [2, 5] {
            let resolved = original
                .clone()
                .with_spatial_viscothermal(resolution(points))
                .unwrap()
                .with_molecular_relaxation(material.clone())
                .unwrap();
            for omega in [700.0, 4000.0] {
                let actual = pressure_schur(&resolved, omega);
                let expected = expected_pressure(&original, &shape, &material, omega, points);
                assert!(
                    relative(&actual, &expected) < 3e-9,
                    "q={points} omega={omega}"
                );
            }
            // Pressure is normalized by frozen physical storage, including
            // the curved area integral rather than a quadrature-lumped mass.
            let state = resolved.project_initial_state(|x| 0.2 + x, |_| 0.0);
            let (x, w) = gauss_legendre_rule(64);
            let bulk = gas().density * material.frozen_sound_speed_m_s().powi(2);
            let energy: f64 = x
                .iter()
                .zip(w)
                .map(|(x, w)| {
                    let x = 0.025 * (x + 1.0);
                    let r = radius(&shape, x);
                    0.5 * 0.025 * w * core::f64::consts::PI * r * r * (0.2 + x).powi(2) / bulk
                })
                .sum();
            assert!((state.stored_energy_j() / energy - 1.0).abs() < 2e-10);
        }
    }
}

#[test]
fn curved_molecular_quadrature_refines_and_material_rates_change_pressure_admittance() {
    let shape = curve();
    let original = base(&shape, 1, 1);
    let material = law(gas(), 1);
    let reference = expected_pressure(&original, &shape, &material, 700.0, 64);
    let mut previous = f64::INFINITY;
    for points in [2, 4, 8] {
        let base = original
            .clone()
            .with_spatial_viscothermal(resolution(points))
            .unwrap();
        let resolved = base
            .clone()
            .with_molecular_relaxation(material.clone())
            .unwrap();
        assert_eq!(resolved.total_dofs() - base.total_dofs(), points * 8);
        let actual = pressure_schur(&resolved, 700.0);
        let error = relative(&actual, &reference);
        assert!(error < previous, "q={points}: {error} >= {previous}");
        previous = error;
        if points == 4 {
            let changed = base
                .with_molecular_relaxation(law_rate(gas(), 1, 0.1))
                .unwrap();
            assert!(relative(&pressure_schur(&changed, 700.0), &actual) > 1e-7);
        }
    }
    assert!(previous < 1e-6, "curved quadrature error {previous}");
}

#[test]
fn spatial_molecular_loss_ports_and_exact_time_evolution_remain_consistent() {
    let shape = Duct {
        segments: vec![Segment::Cone {
            inlet_radius: 0.0007,
            outlet_radius: 0.002,
            length: 0.05,
        }],
    };
    let base = base(&shape, 2, 1)
        .with_spatial_viscothermal(resolution(3))
        .unwrap();
    let junction = [DgSidePort {
        id: PortId(11),
        downstream_element: 1,
    }];
    let material = law(gas(), 1);
    let a = base
        .clone()
        .with_side_ports(&junction)
        .unwrap()
        .with_molecular_relaxation(material.clone())
        .unwrap();
    let b = base
        .clone()
        .with_molecular_relaxation(material)
        .unwrap()
        .with_side_ports(&junction)
        .unwrap();
    for ((ia, pa), (ib, pb)) in a.ports.iter().zip(&b.ports) {
        assert_eq!(ia, ib);
        let scale = pa.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        for (a, b) in pa.iter().zip(pb) {
            assert!((a - b).abs() < 2e-14 * scale);
        }
    }
    let mut arbitrary = a.zero_state();
    for (i, x) in arbitrary.tilde.iter_mut().enumerate() {
        *x = 1e-4 * det::sin(0.37 * i as f64);
    }
    let n = a.total_dofs();
    let derivative: f64 = (0..n)
        .map(|i| {
            arbitrary.tilde[i]
                * (0..n)
                    .map(|j| a.l_tilde[i * n + j] * arbitrary.tilde[j])
                    .sum::<f64>()
        })
        .sum();
    let wall = a.wall_dissipated_power_w(&arbitrary).unwrap();
    let molecular = a.molecular_dissipated_power_w(&arbitrary).unwrap();
    assert!(wall > 0.0 && molecular > 0.0);
    assert!((derivative + wall + molecular).abs() < 2e-11 * (wall + molecular));
    let mut coarse = a.project_initial_state(|x| det::cos(20.0 * x), |_| 0.0);
    let mut fine = coarse.clone();
    let initial = coarse.stored_energy_j();
    let large = a.propagator(2e-6, 1e-11).unwrap();
    let small = a.propagator(1e-6, 1e-11).unwrap();
    for _ in 0..4 {
        let before = coarse.stored_energy_j();
        large.step(&mut coarse, &[]).unwrap();
        small.step(&mut fine, &[]).unwrap();
        small.step(&mut fine, &[]).unwrap();
        assert!(coarse.stored_energy_j() <= before * (1.0 + 1e-10));
    }
    let difference: f64 = coarse
        .tilde
        .iter()
        .zip(&fine.tilde)
        .map(|(a, b)| (a - b) * (a - b))
        .sum();
    assert!(difference < 1e-16 * initial);
    let classical = base
        .clone()
        .with_molecular_relaxation(law(gas(), 0))
        .unwrap();
    assert_eq!(classical.l_tilde, base.l_tilde);
    assert_eq!(classical.operator_identity(), base.operator_identity());
    let field = crate::profile::gas_field::tests::field([300.0; 2], [0.8; 2], 101325.0);
    let mut live = a;
    let before = coarse.clone();
    assert!(
        live.rebind_spatial_viscothermal_gas_field(&mut coarse, &field, 0.05, 1.0, &mut || false)
            .is_err()
    );
    assert_eq!(coarse, before);
    let varying = crate::profile::gas_field::tests::field([300.0, 330.0], [0.8, 0.2], 101325.0);
    let shape = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.001,
            length: varying.length_m(),
        }],
    };
    let varying_duct = DgDuct::new_with_gas_field(
        &shape,
        &varying,
        DgMesh {
            elements: 1,
            order: 1,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap()
    .with_spatial_viscothermal(resolution(3))
    .unwrap();
    let material = law(*varying.state_at(0.0).unwrap().state(), 1);
    assert!(varying_duct.with_molecular_relaxation(material).is_err());
}
