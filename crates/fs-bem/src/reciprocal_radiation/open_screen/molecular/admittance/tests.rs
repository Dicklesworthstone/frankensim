//! G0/G1/G3: physical sheet work in the original relaxing-gas screen.
//! Cubic barycentric moments are independent of the production Duffy rule.
use super::super::tests as fixture;
use super::*;
use crate::helmholtz::bc::{IncidentField, PlaneWave};

fn limits(surface: &OpenPanels, ports: usize) -> MolecularScreenPortLimits {
    let mut limits = fixture::limits(surface);
    limits.ports.maximum_ports = ports;
    limits.ports.maximum_matrix_entries = ports * ports;
    limits.ports.maximum_jump_entries = (ports + 1) * surface.vertices().len();
    limits
}

fn illumination_limits(surface: &OpenPanels, ports: usize) -> MolecularScreenIlluminationLimits {
    MolecularScreenIlluminationLimits {
        boundary: limits(surface, ports),
        maximum_incident_cases: 2,
        maximum_incident_jump_entries: 2 * surface.vertices().len(),
        maximum_source_point_products: 100_000,
        minimum_source_clearance_m: 0.01,
        maximum_live_bytes: 8_000_000,
    }
}

fn sheet(surface: &OpenPanels) -> Vec<C64> {
    surface
        .vertices()
        .iter()
        .map(|&[x, y, _]| C64::new(0.002 * (1.0 + 0.3 * x), 0.0007 * (1.0 - 0.2 * y)))
        .collect()
}

// Integral phi_a*phi_b*phi_c dS = A/10, A/30, A/60 according to
// the multiplicities of three barycentric coordinates. Free-edge Y is kept.
fn sheet_matrix(surface: &OpenPanels, dofs: &[Option<usize>], y: &[C64]) -> Vec<C64> {
    let n = surface.unknowns();
    let mut matrix = vec![C64::ZERO; n * n];
    for (facet, ids) in surface.triangles().iter().enumerate() {
        let (area, _) = fixture::area_normal(surface, facet);
        for a in 0..3 {
            if let Some(i) = dofs[ids[a]] {
                for b in 0..3 {
                    if let Some(j) = dofs[ids[b]] {
                        for c in 0..3 {
                            let weight = if a == b && b == c {
                                area / 10.0
                            } else if a == b || a == c || b == c {
                                area / 30.0
                            } else {
                                area / 60.0
                            };
                            matrix[i * n + j] = matrix[i * n + j] + y[ids[c]].scale(weight);
                        }
                    }
                }
            }
        }
    }
    matrix
}

fn active(dofs: &[Option<usize>], jump: &[C64]) -> Vec<C64> {
    let mut p = vec![C64::ZERO; dofs.iter().flatten().count()];
    for (v, &dof) in dofs.iter().enumerate() {
        if let Some(i) = dof {
            p[i] = jump[v];
        }
    }
    p
}

fn image(matrix: &[C64], p: &[C64]) -> Vec<C64> {
    matrix
        .chunks(p.len())
        .map(|row| {
            row.iter()
                .zip(p)
                .fold(C64::ZERO, |sum, (&a, &b)| sum + a * b)
        })
        .collect()
}

fn quadratic(matrix: &[C64], p: &[C64]) -> C64 {
    p.iter()
        .zip(image(matrix, p))
        .fold(C64::ZERO, |sum, (&a, b)| sum + a.conj() * b)
}

fn sheet_heat(matrix: &[C64], p: &[C64]) -> f64 {
    let conductance: Vec<_> = matrix.iter().map(|y| C64::from_re(y.re)).collect();
    0.5 * quadratic(&conductance, p).re
}

fn close_real(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-300),
        "actual {a}, expected {b}"
    );
}

#[test]
fn no_reservoir_matches_original_robin_and_zero_sheet_keeps_rigid_bits() {
    let surface = fixture::surface();
    let material = fixture::exterior(false);
    let omega = 650.0;
    let k = material.source().relaxing_wavenumber(omega).unwrap();
    let medium = Medium {
        density: material.source().gas().density,
        sound_speed: omega / k.re,
    };
    let owner = PreparedScreen::new(&surface, 2.0, fixture::options(1), || false).unwrap();
    let rows = fixture::velocity_rows(&surface);
    let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
    let mut y = sheet(&surface);
    let original_y = y.clone();
    let beta: Vec<_> = y
        .iter()
        .map(|y| y.scale(medium.density * medium.sound_speed))
        .collect();
    let controls = limits(&surface, 2);
    let original = owner
        .solve_velocity_ports_with_local_admittance(
            k.re,
            medium,
            &refs,
            &beta,
            controls.ports,
            || false,
        )
        .unwrap();
    let actual = owner
        .solve_molecular_velocity_ports_with_local_admittance(
            &material,
            omega,
            &refs,
            &y,
            controls,
            || false,
        )
        .unwrap();
    y.fill(C64::ZERO);
    assert_eq!(
        actual.local_admittance_m_per_pa_s(),
        Some(original_y.as_slice())
    );
    assert!(std::ptr::eq(actual.source().source(), material.source()));
    fixture::close(actual.impedance(), original.impedance(), 5e-13);
    for (a, b) in actual.pressure_jumps_pa().iter().zip(original.solutions()) {
        fixture::close(a, b.pressure_jump_pa(), 5e-13);
    }
    let rates = [C64::new(0.7, -0.2), C64::new(-0.3, 0.6)];
    let power = actual.power(&rates, || false).unwrap();
    close_real(
        power.mechanical_power_w,
        original.input_power_w(&rates).unwrap(),
        2e-12,
    );
    close_real(
        power.sheet_absorbed_power_w,
        original.absorbed_power_w(&rates, || false).unwrap(),
        2e-12,
    );
    close_real(
        power.fluid_power_w,
        original.radiated_power_w(&rates, || false).unwrap(),
        2e-11,
    );

    let rigid = owner
        .solve_molecular_velocity_ports(&material, omega, &refs, controls, || false)
        .unwrap();
    let zero = owner
        .solve_molecular_velocity_ports_with_local_admittance(
            &material,
            omega,
            &refs,
            &y,
            controls,
            || false,
        )
        .unwrap();
    assert_eq!(zero.impedance(), rigid.impedance());
    assert_eq!(zero.pressure_jumps_pa(), rigid.pressure_jumps_pa());
    assert_eq!(
        zero.power(&rates, || false).unwrap(),
        rigid.power(&rates, || false).unwrap()
    );
}

#[test]
fn relaxing_gas_and_spatial_sheet_close_the_original_weak_work_independently() {
    let surface = fixture::surface();
    let material = fixture::exterior(true);
    let omega = 650.0;
    let owner = PreparedScreen::new(&surface, 2.0, fixture::options(1), || false).unwrap();
    let rows = fixture::velocity_rows(&surface);
    let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
    let y = sheet(&surface);
    let actual = owner
        .solve_molecular_velocity_ports_with_local_admittance(
            &material,
            omega,
            &refs,
            &y,
            limits(&surface, 2),
            || false,
        )
        .unwrap();
    let k = actual.wavenumber_per_m();
    assert!(k.im > 0.0);
    assert_eq!(k, material.source().relaxing_wavenumber(omega).unwrap());
    let rho_omega = material.source().gas().density * omega;
    let h = owner
        .complex_matrix(k, fixture::operator_limits(), || false)
        .unwrap();
    let mass = sheet_matrix(&surface, owner.pressure_dof_map(), &y);
    let physical = C64::new(0.0, rho_omega);
    let robin: Vec<_> = h
        .iter()
        .zip(&mass)
        .map(|(&a, &b)| a + physical * b)
        .collect();
    let weak = owner
        .velocity_port_pairing(&refs, 1024, u64::MAX, || false)
        .unwrap();
    for (j, jump) in actual.pressure_jumps_pa().iter().enumerate() {
        let p = active(owner.pressure_dof_map(), jump);
        let rhs: Vec<_> = weak[j].iter().map(|&b| physical.scale(b)).collect();
        fixture::close(&image(&robin, &p), &rhs, 8e-13);
        // A complex wave impedance must not replace inviscid physical momentum.
        let wrong_scale = C64::new(0.0, 1.0)
            * k.scale(material.source().gas().density * material.source().gas().sound_speed);
        assert!((wrong_scale - physical).abs() > 1e-4 * physical.abs());
    }
    fixture::close(&[actual.impedance()[1]], &[actual.impedance()[2]], 2e-12);
    let rates = [C64::new(0.7, -0.2), C64::new(-0.3, 0.6)];
    let jump = actual.pressure_jump(&rates, || false).unwrap();
    let p = active(owner.pressure_dof_map(), &jump);
    let power = actual.power(&rates, || false).unwrap();
    let fluid = 0.5 * (quadratic(&h, &p) * C64::new(0.0, -1.0 / rho_omega)).re;
    let heat = sheet_heat(&mass, &p);
    let mechanical = 0.5
        * rows
            .iter()
            .zip(rates)
            .map(|(row, rate)| (rate.conj() * fixture::effort(&surface, row, &jump)).re)
            .sum::<f64>();
    assert!(fluid > 0.0 && heat > 0.0);
    close_real(power.sheet_absorbed_power_w, heat, 3e-13);
    close_real(power.fluid_power_w, fluid, 4e-12);
    close_real(power.mechanical_power_w, mechanical, 4e-13);
    close_real(mechanical, fluid + heat, 4e-12);
    let incoherent: f64 = actual
        .pressure_jumps_pa()
        .iter()
        .zip(rates)
        .map(|(jump, rate)| {
            rate.norm_sq() * sheet_heat(&mass, &active(owner.pressure_dof_map(), jump))
        })
        .sum();
    assert!(
        (incoherent - heat).abs() > 0.01 * heat,
        "coherent cross-port sheet work must survive"
    );
}

#[test]
fn illuminated_sheet_keeps_incident_motion_interference_and_transparent_tracking() {
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
    let material = fixture::exterior(true);
    let omega = 650.0;
    let owner = PreparedScreen::new(&surface, 2.0, fixture::options(1), || false).unwrap();
    let uniform = vec![[1.0; 3]; surface.triangles().len()];
    let y = sheet(&surface);
    let plane = PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::new(0.7, -0.2),
    };
    let incident = IncidentField::PlaneWave(plane);
    let actual = owner
        .solve_molecular_illumination_with_local_admittance(
            &material,
            omega,
            &[&uniform],
            &y,
            &[incident],
            illumination_limits(&surface, 1),
            || false,
        )
        .unwrap();
    let k = actual.motion().wavenumber_per_m();
    let incident_velocity =
        (k * plane.amplitude_pa).scale(1.0 / (material.source().gas().density * omega));
    let blocked = actual.power(0, &[C64::ZERO], || false).unwrap();
    assert_eq!(blocked.mechanical_power_w, 0.0);
    assert!(blocked.sheet_absorbed_power_w > 0.0);
    assert_eq!(blocked.fluid_power_w, -blocked.sheet_absorbed_power_w);
    let rates = [incident_velocity.scale(0.25)];
    let coupled = actual.power(0, &rates, || false).unwrap();
    assert!(
        coupled.mechanical_power_w < 0.0,
        "incident sound supplies structural work"
    );
    close_real(
        coupled.sheet_absorbed_power_w,
        0.75_f64.powi(2) * blocked.sheet_absorbed_power_w,
        4e-12,
    );
    let jump = actual.total_pressure_jump(0, &rates, || false).unwrap();
    let mass = sheet_matrix(&surface, owner.pressure_dof_map(), &y);
    close_real(
        coupled.sheet_absorbed_power_w,
        sheet_heat(&mass, &active(owner.pressure_dof_map(), &jump)),
        3e-13,
    );
    let force = actual.total_fluid_reaction(0, &rates, || false).unwrap()[0];
    close_real(
        coupled.mechanical_power_w,
        -0.5 * (rates[0].conj() * force).re,
        4e-13,
    );

    // Moving at the local incident velocity creates no pressure jump or leakage:
    // total sound is the unchanged incident p/u, even with molecular and sheet loss.
    let tracking = [incident_velocity];
    let jump = actual.total_pressure_jump(0, &tracking, || false).unwrap();
    let scale = actual.blocked_pressure_jumps_pa()[0]
        .iter()
        .map(|p| p.abs())
        .fold(0.0, f64::max);
    assert!(jump.iter().all(|p| p.abs() < 2e-12 * scale));
    let power = actual.power(0, &tracking, || false).unwrap();
    assert!(power.sheet_absorbed_power_w < 1e-23 * blocked.sheet_absorbed_power_w);
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
    let evaluations = receiver.kernel_evaluations();
    let total = receiver.evaluate(0, &tracking, || false).unwrap();
    let incoming = material
        .prepare_incident(omega, incident, || false)
        .unwrap();
    for (i, point) in points.into_iter().enumerate() {
        let expected = incoming.pressure_velocity_at(point, || false).unwrap();
        fixture::close(&[total.scalar.pressure[i]], &[expected.pressure_pa], 2e-12);
        fixture::close(
            &total.particle_velocity_m_s[i],
            &expected.velocity_m_s,
            2e-12,
        );
    }
    assert_eq!(receiver.kernel_evaluations(), evaluations);
}

#[test]
fn high_q_sheet_heat_keeps_phase_and_passive_resource_admission() {
    let surface = fixture::surface();
    let material = fixture::exterior(true);
    let omega = 650.0;
    let owner = PreparedScreen::new(&surface, 2.0, fixture::options(1), || false).unwrap();
    let rows = fixture::velocity_rows(&surface);
    let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
    let rates = [C64::new(0.7, -0.2), C64::new(-0.3, 0.6)];
    let phase = C64::new(0.71_f64.cos(), 0.71_f64.sin());
    for conductance in [0.0, 1e-20] {
        let y = vec![C64::new(conductance, 0.01); surface.vertices().len()];
        let actual = owner
            .solve_molecular_velocity_ports_with_local_admittance(
                &material,
                omega,
                &refs,
                &y,
                limits(&surface, 2),
                || false,
            )
            .unwrap();
        let jump = actual.pressure_jump(&rates, || false).unwrap();
        let expected = sheet_heat(
            &sheet_matrix(&surface, owner.pressure_dof_map(), &y),
            &active(owner.pressure_dof_map(), &jump),
        );
        let heat = actual.sheet_absorbed_power_w(&rates, || false).unwrap();
        let rotated = rates.map(|rate| rate * phase);
        let rotated_heat = actual.sheet_absorbed_power_w(&rotated, || false).unwrap();
        if conductance == 0.0 {
            assert_eq!(heat, 0.0);
            assert_eq!(rotated_heat, 0.0);
        } else {
            assert!(heat > 0.0);
            close_real(heat, expected, 4e-13);
            close_real(rotated_heat, heat, 4e-13);
        }
    }
    let y = sheet(&surface);
    for invalid in [
        vec![C64::new(-1e-20, 0.0); y.len()],
        vec![C64::new(f64::NAN, 0.0); y.len()],
        vec![C64::new(0.0, f64::INFINITY); y.len()],
        y[..y.len() - 1].to_vec(),
    ] {
        assert!(
            owner
                .solve_molecular_velocity_ports_with_local_admittance(
                    &material,
                    omega,
                    &refs,
                    &invalid,
                    limits(&surface, 2),
                    || false,
                )
                .is_err()
        );
    }
    for kind in 0..3 {
        let mut short = limits(&surface, 2);
        match kind {
            0 => short.ports.maximum_jump_entries -= 1,
            1 => {
                short.ports.maximum_work = owner
                    .molecular_velocity_port_work_with_local_admittance(2)
                    .unwrap()
                    - 1
            }
            _ => {
                short.maximum_workspace_bytes =
                    owner.molecular_port_boundary_payload(2, true).unwrap() as usize - 1
            }
        }
        assert!(
            owner
                .solve_molecular_velocity_ports_with_local_admittance(
                    &material,
                    omega,
                    &refs,
                    &y,
                    short,
                    || false,
                )
                .is_err()
        );
    }
    assert!(
        owner
            .solve_molecular_velocity_ports_with_local_admittance(
                &material,
                omega,
                &refs,
                &y,
                limits(&surface, 2),
                || true,
            )
            .is_err()
    );
}
