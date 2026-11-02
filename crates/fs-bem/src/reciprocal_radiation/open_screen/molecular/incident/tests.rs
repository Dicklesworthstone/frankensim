//! G1/G3: original pressure-jump forcing, work transpose and total field.
//! Manufactured gas data verify the implementation, not measured shell sound.
use super::super::tests as fixture;
use super::*;
use crate::helmholtz::bc::{
    Admittance, ForceDipole, IncidentSource, Monopole, PlaneWave, SurfaceBc,
};

fn limits(surface: &OpenPanels) -> MolecularScreenIlluminationLimits {
    MolecularScreenIlluminationLimits {
        boundary: fixture::limits(surface),
        maximum_incident_cases: 3,
        maximum_incident_jump_entries: 3 * surface.vertices().len(),
        maximum_source_point_products: 100_000,
        minimum_source_clearance_m: 0.01,
        maximum_live_bytes: 8_000_000,
    }
}

fn gauss(n: usize) -> Vec<(f64, f64)> {
    (0..n)
        .map(|i| {
            let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
            let mut derivative = 0.0;
            for _ in 0..32 {
                let (mut a, mut b) = (1.0, x);
                for degree in 2..=n {
                    let next =
                        ((2 * degree - 1) as f64 * x * b - (degree - 1) as f64 * a) / degree as f64;
                    a = b;
                    b = next;
                }
                derivative = n as f64 * (x * b - a) / (x * x - 1.0);
                let update = b / derivative;
                x -= update;
                if update.abs() < 2e-16 {
                    break;
                }
            }
            (
                0.5 * (x + 1.0),
                1.0 / ((1.0 - x * x) * derivative * derivative),
            )
        })
        .collect()
}

// A different triangle parameterization and independently generated Gauss
// points reconstruct the incident weak RHS and the solved P1 field.
fn integrate(
    surface: &OpenPanels,
    rule: &[(f64, f64)],
    mut sample: impl FnMut(usize, [f64; 3], Point, Point, f64),
) {
    for (facet, ids) in surface.triangles().iter().enumerate() {
        let (area, normal) = fixture::area_normal(surface, facet);
        let t = ids.map(|v| surface.vertices()[v]);
        for &(a, wa) in rule {
            for &(b, wb) in rule {
                let shape = [1.0 - a, a * (1.0 - b), a * b];
                let point =
                    core::array::from_fn(|axis| (0..3).map(|i| shape[i] * t[i][axis]).sum());
                sample(facet, shape, point, normal, 2.0 * area * a * wa * wb);
            }
        }
    }
}

fn plane_pressure(plane: PlaneWave, k: C64, point: Point) -> C64 {
    let coordinate = (0..3).map(|a| plane.direction[a] * point[a]).sum::<f64>();
    let phase = k.re * coordinate;
    plane.amplitude_pa * C64::new(phase.cos(), phase.sin()).scale((-k.im * coordinate).exp())
}

fn double_layer(surface: &OpenPanels, k: C64, jump: &[C64], point: Point, order: usize) -> C64 {
    let mut result = C64::ZERO;
    integrate(surface, &gauss(order), |facet, shape, y, normal, weight| {
        let r_vec: Point = core::array::from_fn(|a| point[a] - y[a]);
        let r = r_vec.iter().map(|v| v * v).sum::<f64>().sqrt();
        let n_dot_r = (0..3).map(|a| normal[a] * r_vec[a]).sum::<f64>() / r;
        let phase = k.re * r;
        let green = C64::new(phase.cos(), phase.sin())
            .scale((-k.im * r).exp() / (4.0 * std::f64::consts::PI * r));
        let derivative = green * C64::new(1.0 / r + k.im, -k.re).scale(n_dot_r);
        let p = (0..3).fold(C64::ZERO, |sum, a| {
            sum + jump[surface.triangles()[facet][a]].scale(shape[a])
        });
        result = result + (derivative * p).scale(weight);
    });
    result
}

#[test]
fn empty_reservoir_incident_drive_recovers_original_rigid_screen_and_force_pairing() {
    let curved = fixture::surface();
    let surface = OpenPanels::new(
        curved
            .vertices()
            .iter()
            .map(|p| [p[0], p[1], 0.0])
            .collect(),
        curved.triangles().to_vec(),
    )
    .unwrap();
    let gas = fixture::exterior(false);
    let omega = 650.0;
    let k = gas.source().relaxing_wavenumber(omega).unwrap();
    let medium = Medium {
        density: gas.source().gas().density,
        sound_speed: omega / k.re,
    };
    let plane = PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::new(0.7, -0.2),
    };
    let field = IncidentField::PlaneWave(plane);
    let owner = PreparedScreen::new(&surface, 6.0, fixture::options(1), || false).unwrap();
    let rows = fixture::velocity_rows(&surface);
    let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
    let actual = owner
        .solve_molecular_illumination(&gas, omega, &refs, &[field], limits(&surface), || false)
        .unwrap();
    let drive = vec![
        -plane
            .amplitude_pa
            .scale(1.0 / (medium.density * medium.sound_speed));
        surface.vertices().len()
    ];
    let original = owner
        .solve(
            k.re,
            medium,
            SurfaceBc {
                drive_velocity: &drive,
                admittance: Admittance::Rigid,
            },
            || false,
        )
        .unwrap();
    fixture::close(
        &actual.blocked_pressure_jumps_pa()[0],
        original.pressure_jump_pa(),
        5e-13,
    );
    let original_motion = owner
        .solve_molecular_velocity_ports(&gas, omega, &refs, fixture::limits(&surface), || false)
        .unwrap();
    fixture::close(
        actual.motion().impedance(),
        original_motion.impedance(),
        5e-13,
    );
    for (mode, row) in rows.iter().enumerate() {
        fixture::close(
            &[actual.blocked_forces()[0][mode]],
            &[-fixture::effort(&surface, row, original.pressure_jump_pa())],
            5e-13,
        );
    }
    let rates = [C64::new(0.001, -0.0002), C64::new(-0.0003, 0.0006)];
    let motion_jump = actual.motion().pressure_jump(&rates, || false).unwrap();
    let expected: Vec<_> = original
        .pressure_jump_pa()
        .iter()
        .zip(motion_jump)
        .map(|(&a, b)| a + b)
        .collect();
    fixture::close(
        &actual.total_pressure_jump(0, &rates, || false).unwrap(),
        &expected,
        5e-13,
    );
    let reaction: Vec<_> = rows
        .iter()
        .map(|row| -fixture::effort(&surface, row, &expected))
        .collect();
    fixture::close(
        &actual.total_fluid_reaction(0, &rates, || false).unwrap(),
        &reaction,
        5e-13,
    );
    assert!(actual.maximum_backward_residual() < 1e-10);
}

#[test]
fn relaxing_screen_keeps_physical_weak_drive_and_coherent_total_receivers() {
    let surface = fixture::surface();
    let gas = fixture::exterior(true);
    let omega = 650.0;
    let k = gas.source().relaxing_wavenumber(omega).unwrap();
    let rho_omega = gas.source().gas().density * omega;
    let plane = PlaneWave {
        direction: [0.3, -0.4, 0.75_f64.sqrt()],
        amplitude_pa: C64::new(0.7, -0.2),
    };
    let finite = [
        IncidentSource::Monopole(Monopole {
            position_m: [0.2, -0.1, 1.2],
            volume_velocity_m3_s: C64::new(0.004, 0.002),
        }),
        IncidentSource::ForceDipole(ForceDipole {
            position_m: [-0.2, 0.3, -1.2],
            force_on_fluid_n: [C64::new(0.03, -0.02), C64::ZERO, C64::new(0.08, 0.01)],
        }),
    ];
    let all = [IncidentSource::PlaneWave(plane), finite[0], finite[1]];
    let incidents = [
        IncidentField::PlaneWave(plane),
        IncidentField::Sources(&finite),
        IncidentField::Sources(&all),
    ];
    let owner = PreparedScreen::new(&surface, 6.0, fixture::options(1), || false).unwrap();
    let rows = fixture::velocity_rows(&surface);
    let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
    let actual = owner
        .solve_molecular_illumination(&gas, omega, &refs, &incidents, limits(&surface), || false)
        .unwrap();
    let mut weak = vec![C64::ZERO; owner.unknowns()];
    let dofs = owner.pressure_dof_map();
    integrate(
        &surface,
        &gauss(12),
        |facet, shape, point, normal, weight| {
            let normal_projection = (0..3).map(|a| plane.direction[a] * normal[a]).sum::<f64>();
            let drive = -(k * plane_pressure(plane, k, point)).scale(normal_projection / rho_omega);
            for a in 0..3 {
                if let Some(u) = dofs[surface.triangles()[facet][a]] {
                    weak[u] = weak[u] + drive.scale(weight * shape[a]);
                }
            }
        },
    );
    let h = owner
        .complex_matrix(k, fixture::operator_limits(), || false)
        .unwrap();
    let mut active = vec![C64::ZERO; owner.unknowns()];
    for (vertex, &dof) in dofs.iter().enumerate() {
        if let Some(u) = dof {
            active[u] = actual.blocked_pressure_jumps_pa()[0][vertex];
        }
    }
    let lhs: Vec<_> = h
        .chunks(owner.unknowns())
        .map(|row| {
            row.iter()
                .zip(&active)
                .fold(C64::ZERO, |sum, (&a, &b)| sum + a * b)
        })
        .collect();
    let rhs: Vec<_> = weak.iter().map(|&v| C64::new(0.0, rho_omega) * v).collect();
    fixture::close(&lhs, &rhs, 2e-9);
    let false_momentum =
        C64::new(0.0, 1.0) * k.scale(gas.source().gas().density * gas.source().gas().sound_speed);
    let false_rhs: Vec<_> = weak.iter().map(|&v| false_momentum * v).collect();
    assert!(
        lhs.iter()
            .zip(false_rhs)
            .any(|(&a, b)| (a - b).abs() > 1e-4 * a.abs())
    );
    let summed: Vec<_> = actual.blocked_pressure_jumps_pa()[0]
        .iter()
        .zip(&actual.blocked_pressure_jumps_pa()[1])
        .map(|(&a, &b)| a + b)
        .collect();
    fixture::close(&actual.blocked_pressure_jumps_pa()[2], &summed, 2e-12);
    let rates = [
        [C64::new(0.001, -0.0002), C64::new(-0.0003, 0.0006)],
        [C64::new(-0.0002, 0.0004), C64::new(0.0005, -0.0003)],
    ];
    let total_rates = [rates[0][0] + rates[1][0], rates[0][1] + rates[1][1]];
    let points = [[0.1, -0.05, 0.65], [-0.08, 0.1, -0.65]];
    let receiver = actual
        .prepare_receivers(
            &points,
            0.01,
            crate::near_field::Options::default(),
            1_000_000,
            || false,
        )
        .unwrap();
    let work = receiver.kernel_evaluations();
    let a = receiver.evaluate(0, &rates[0], || false).unwrap();
    let b = receiver.evaluate(1, &rates[1], || false).unwrap();
    let combined = receiver.evaluate(2, &total_rates, || false).unwrap();
    let jump = actual.total_pressure_jump(0, &rates[0], || false).unwrap();
    let total_pressure =
        |point| plane_pressure(plane, k, point) + double_layer(&surface, k, &jump, point, 20);
    for (i, &point) in points.iter().enumerate() {
        fixture::close(
            &[combined.scalar.pressure[i]],
            &[a.scalar.pressure[i] + b.scalar.pressure[i]],
            5e-12,
        );
        fixture::close(
            &combined.particle_velocity_m_s[i],
            &core::array::from_fn::<_, 3, _>(|axis| {
                a.particle_velocity_m_s[i][axis] + b.particle_velocity_m_s[i][axis]
            }),
            5e-12,
        );
        let expected = total_pressure(point);
        fixture::close(&[a.scalar.pressure[i]], &[expected], 5e-7);
        let lower_order =
            plane_pressure(plane, k, point) + double_layer(&surface, k, &jump, point, 12);
        fixture::close(&[lower_order], &[expected], 2e-10);
        let step = 4e-5;
        let gradient: [C64; 3] = core::array::from_fn(|axis| {
            let value = |offset: f64| {
                let mut x = point;
                x[axis] += offset * step;
                total_pressure(x)
            };
            (-value(2.0) + value(1.0).scale(8.0) - value(-1.0).scale(8.0) + value(-2.0))
                .scale(1.0 / (12.0 * step))
        });
        fixture::close(
            &a.particle_velocity_m_s[i],
            &gradient.map(|v| v * C64::new(0.0, -1.0 / rho_omega)),
            8e-7,
        );
        assert!(a.scalar.quadrature_error_estimate_pa[i] >= 0.0);
    }
    assert_eq!(receiver.kernel_evaluations(), work);
    assert_eq!(receiver.points_m(), points);
}

#[test]
fn source_clearance_whole_batch_limits_cancellation_and_replay_remain_owned() {
    let surface = fixture::surface();
    let gas = fixture::exterior(true);
    let owner = PreparedScreen::new(&surface, 6.0, fixture::options(1), || false).unwrap();
    let rows = fixture::velocity_rows(&surface);
    let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
    let plane = IncidentField::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::ONE,
    });
    let incidents = [plane, plane];
    let controls = limits(&surface);
    for kind in 0..5 {
        let mut short = controls;
        match kind {
            0 => short.maximum_incident_cases = 1,
            1 => short.maximum_incident_jump_entries = 2 * surface.vertices().len() - 1,
            2 => short.maximum_source_point_products = 1,
            3 => short.maximum_live_bytes = 0,
            _ => short.boundary.ports.maximum_matrix_entries = 3,
        }
        assert!(
            owner
                .solve_molecular_illumination(&gas, 650.0, &refs, &incidents, short, || false)
                .is_err()
        );
    }
    // Near a free edge, far from either panel centroid: face/edge geometry must
    // own this refusal, not a source-to-centroid distance approximation.
    let t = surface.triangles()[0].map(|v| surface.vertices()[v]);
    let (_, normal) = fixture::area_normal(&surface, 0);
    let near_edge: Point = core::array::from_fn(|a| 0.5 * (t[0][a] + t[1][a]) + 0.005 * normal[a]);
    for position in [t[0], near_edge] {
        let emitters = [Monopole {
            position_m: position,
            volume_velocity_m3_s: C64::ONE,
        }];
        assert!(
            owner
                .solve_molecular_illumination(
                    &gas,
                    650.0,
                    &refs,
                    &[IncidentField::Monopoles(&emitters)],
                    controls,
                    || false
                )
                .is_err()
        );
    }
    assert!(matches!(
        owner.solve_molecular_illumination(&gas, 650.0, &refs, &incidents, controls, || true),
        Err(MolecularError::Boundary(HelmholtzError::Cancelled))
    ));
    let mut polls = 0;
    let original = owner
        .solve_molecular_illumination(&gas, 650.0, &refs, &incidents, controls, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut cancelled = 0;
    assert!(matches!(
        owner.solve_molecular_illumination(&gas, 650.0, &refs, &incidents, controls, || {
            cancelled += 1;
            cancelled >= polls - 1
        }),
        Err(MolecularError::Boundary(HelmholtzError::Cancelled))
    ));
    let replay = owner
        .solve_molecular_illumination(&gas, 650.0, &refs, &incidents, controls, || false)
        .unwrap();
    assert_eq!(original.motion().impedance(), replay.motion().impedance());
    assert_eq!(
        original.blocked_pressure_jumps_pa(),
        replay.blocked_pressure_jumps_pa()
    );
    assert_eq!(original.blocked_forces(), replay.blocked_forces());
    let rates = [C64::ZERO; 2];
    assert!(original.total_pressure_jump(2, &rates, || false).is_err());
    assert!(
        original
            .total_fluid_reaction(0, &[C64::ONE], || false)
            .is_err()
    );
    assert!(matches!(
        original.total_pressure_jump(0, &rates, || true),
        Err(HelmholtzError::Cancelled)
    ));
    let points = [[0.0, 0.0, 0.65]];
    assert!(
        original
            .prepare_receivers(
                &points,
                0.01,
                crate::near_field::Options::default(),
                1,
                || false
            )
            .is_err()
    );
    let receiver = original
        .prepare_receivers(
            &points,
            0.01,
            crate::near_field::Options::default(),
            1_000_000,
            || false,
        )
        .unwrap();
    assert!(receiver.evaluate(2, &rates, || false).is_err());
    assert!(matches!(
        receiver.evaluate(0, &rates, || true),
        Err(HelmholtzError::Cancelled)
    ));
    let a = receiver.evaluate(0, &rates, || false).unwrap();
    let b = receiver.evaluate(0, &rates, || false).unwrap();
    assert_eq!(a.scalar.pressure, b.scalar.pressure);
    assert_eq!(a.particle_velocity_m_s, b.particle_velocity_m_s);
}
