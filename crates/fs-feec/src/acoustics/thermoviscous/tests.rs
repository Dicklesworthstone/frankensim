use super::*;

mod conductive;
mod fluid;

fn air() -> ThermoviscousMedium {
    let temperature = 293.15;
    let pressure = 101325.0;
    ThermoviscousMedium {
        density_kg_m3: pressure / (287.05 * temperature),
        pressure_pa: pressure,
        temperature_k: temperature,
        cp_j_kg_k: 1005.0,
        dynamic_viscosity_pa_s: 1.81e-5,
        thermal_conductivity_w_m_k: 0.0257,
        bulk_viscosity_pa_s: 0.0,
    }
}

fn box_mesh(cells: [usize; 3], size: [f64; 3]) -> (TetComplex, Vec<[f64; 3]>) {
    let index =
        |i: usize, j: usize, k: usize| ((i * (cells[1] + 1) + j) * (cells[2] + 1) + k) as u32;
    let mut positions = Vec::new();
    for i in 0..=cells[0] {
        for j in 0..=cells[1] {
            for k in 0..=cells[2] {
                positions.push([
                    size[0] * i as f64 / cells[0] as f64,
                    size[1] * j as f64 / cells[1] as f64,
                    size[2] * k as f64 / cells[2] as f64,
                ]);
            }
        }
    }
    let mut tets = Vec::new();
    for i in 0..cells[0] {
        for j in 0..cells[1] {
            for k in 0..cells[2] {
                for permutation in [
                    [0, 1, 2],
                    [0, 2, 1],
                    [1, 0, 2],
                    [1, 2, 0],
                    [2, 0, 1],
                    [2, 1, 0],
                ] {
                    let mut corner = [i, j, k];
                    let mut tet = [0; 4];
                    tet[0] = index(corner[0], corner[1], corner[2]);
                    for (step, axis) in permutation.into_iter().enumerate() {
                        corner[axis] += 1;
                        tet[step + 1] = index(corner[0], corner[1], corner[2]);
                    }
                    tets.push(tet);
                }
            }
        }
    }
    (TetComplex::from_tets(positions.len(), tets), positions)
}

fn boundaries(
    mesh: &TetComplex,
    pos: &[[f64; 3]],
    size: [f64; 3],
    isothermal: bool,
) -> ThermoviscousBoundary {
    let mut left = Vec::new();
    let mut right = Vec::new();
    let mut walls = Vec::new();
    for &face in &mesh.faces {
        let on = |axis: usize, value: f64| {
            face.iter()
                .all(|&v| (pos[v as usize][axis] - value).abs() < 1e-13 * size[axis])
        };
        if on(0, 0.0) {
            left.push(face);
        } else if on(0, size[0]) {
            right.push(face);
        } else if on(1, 0.0) || on(1, size[1]) || on(2, 0.0) || on(2, size[2]) {
            walls.push(face);
        }
    }
    ThermoviscousBoundary {
        isothermal_faces: if isothermal {
            walls.clone()
        } else {
            Vec::new()
        },
        no_slip_faces: walls,
        ports: vec![BoundaryMode::uniform(&left), BoundaryMode::uniform(&right)],
    }
}

fn limits() -> ThermoviscousLimits {
    ThermoviscousLimits {
        max_dofs: MAX_LNS_DOFS,
        max_dense_bytes: 128 * 1024 * 1024,
        max_elements: MAX_LNS_ELEMENTS,
        max_element_aspect_ratio: 1000.0,
    }
}

fn assert_work(field: &ThermoviscousField) {
    let p = field.power;
    assert!(
        p.viscous_w >= 0.0
            && p.bulk_w >= 0.0
            && p.conductive_w >= 0.0
            && p.solid_conductive_w >= 0.0
            && p.caloric_w >= 0.0,
        "{p:?}"
    );
    let scale = p.input_w.abs().max(p.dissipated_w());
    assert!(
        p.balance_residual_w().abs() < 2e-8 * scale.max(1e-30),
        "{p:?}"
    );
}

#[test]
fn full_lns_work_reciprocity_complex_caloric_loss_and_power_scattering() {
    let size = [0.004, 0.0008, 0.001];
    let (mesh, pos) = box_mesh([2, 2, 2], size);
    let boundary = boundaries(&mesh, &pos, size, true);
    let mut gas = air();
    gas.bulk_viscosity_pa_s = 1.2e-5;
    let model = ThermoviscousFe::new(&mesh, &pos, gas, &boundary, limits(), &mut || false).unwrap();
    let cp = C64::new(970.0, -45.0);
    let prepared = model
        .prepare_with_heat_capacity(1700.0, cp, &mut || false)
        .unwrap();
    let drive = [C64::new(1.0, 0.2), C64::new(-0.3, 0.1)];
    let field = prepared.response(&drive, &mut || false).unwrap();
    assert_work(&field);
    assert!(
        field.power.viscous_w > 0.0
            && field.power.bulk_w > 0.0
            && field.power.conductive_w > 0.0
            && field.power.caloric_w > 0.0
    );
    let repeat = prepared.response(&drive, &mut || false).unwrap();
    assert_eq!(field.pressure_pa, repeat.pressure_pa);
    assert_eq!(field.velocity_m_s, repeat.velocity_m_s);
    let y = model
        .admittance_with_heat_capacity(1700.0, cp, &mut || false)
        .unwrap();
    assert!((y[1] - y[2]).abs() < 1e-10 * y[1].abs());
    for i in 0..2 {
        let q = y[2 * i] * drive[0] + y[2 * i + 1] * drive[1];
        assert!((q - field.inward_flow_m3_s[i]).abs() < 1e-11 * q.abs());
    }
    let reference = [5e8, 3e8];
    let s = model
        .scattering_with_heat_capacity(1700.0, cp, &reference, &mut || false)
        .unwrap();
    let incident = [C64::new(0.7, 0.2), C64::new(-0.3, 0.6)];
    let reflected = [
        s[0] * incident[0] + s[1] * incident[1],
        s[2] * incident[0] + s[3] * incident[1],
    ];
    let pressure = std::array::from_fn::<_, 2, _>(|i| {
        (incident[i] + reflected[i]).scale(det::sqrt(reference[i]))
    });
    let actual = prepared.response(&pressure, &mut || false).unwrap();
    assert_work(&actual);
    let wave_power = 0.5
        * (incident.iter().map(|v| v.norm_sq()).sum::<f64>()
            - reflected.iter().map(|v| v.norm_sq()).sum::<f64>());
    assert!(wave_power > 0.0);
    assert!((wave_power - actual.input_power_w()).abs() < 1e-9 * wave_power);
    for i in 0..2 {
        let q = (incident[i] - reflected[i]).scale(1.0 / det::sqrt(reference[i]));
        assert!((q - actual.inward_flow_m3_s[i]).abs() < 1e-10 * q.abs());
    }
}

#[test]
fn full_lns_low_frequency_isothermal_adiabatic_and_relaxing_compliance() {
    // Independent ideal-gas closed-volume limits. The open traction surfaces
    // supply the small compression flow; tangential no-slip remains resolved.
    let size = [0.003, 0.001, 0.0012];
    let volume = size.iter().product::<f64>();
    let (mesh, pos) = box_mesh([2, 2, 2], size);
    let gas = air();
    for (isothermal, cp) in [
        (true, C64::new(gas.cp_j_kg_k, 0.0)),
        (false, C64::new(gas.cp_j_kg_k, 0.0)),
        (false, C64::new(910.0, -80.0)),
    ] {
        // The isothermal limit needs a period long compared with transverse
        // diffusion. An insulated volume has no such slow thermal limit; use
        // 10 Hz to avoid subtracting nearly balanced static tractions to infer
        // a 1e-17 m^3/s compression flow. Both remain far below the first mode.
        let frequency = if isothermal { 1e-3 } else { 10.0 };
        let omega = std::f64::consts::TAU * frequency;
        let boundary = boundaries(&mesh, &pos, size, isothermal);
        let model =
            ThermoviscousFe::new(&mesh, &pos, gas, &boundary, limits(), &mut || false).unwrap();
        let solution = model
            .prepare_with_heat_capacity(frequency, cp, &mut || false)
            .unwrap()
            .response(&[C64::ONE, C64::ONE], &mut || false)
            .unwrap();
        let compliance = if isothermal {
            C64::new(volume / gas.pressure_pa, 0.0)
        } else {
            (C64::ONE - C64::new(gas.gas_constant_j_kg_k(), 0.0) / cp)
                .scale(volume / gas.pressure_pa)
        };
        let expected = C64::new(0.0, omega) * compliance;
        let q = solution.inward_flow_m3_s[0] + solution.inward_flow_m3_s[1];
        let error = (q - expected).abs() / expected.abs();
        eprintln!("LNS compliance isothermal={isothermal} cp={cp:?} relative_error={error:e}");
        assert!(error < 2e-5);
        if cp.im < 0.0 {
            assert_work(&solution);
            assert!(solution.power.caloric_w > 0.0);
        }
    }
}

#[test]
fn full_lns_rectangular_gap_mesh_ladder_against_poiseuille() {
    // Long finite rectangular duct, quasi-steady limit. The analytic target is
    // fully developed Poiseuille flow, so the remaining finite-end defect is
    // explicitly present; this is not the 1e-4 infinite-slit ZK qualification.
    let size = [0.05, 0.0005, 0.001];
    let gas = air();
    let [length, height, width] = size;
    let sum = (1..=501)
        .step_by(2)
        .map(|m| {
            let m = m as f64;
            (m * std::f64::consts::PI * width / (2.0 * height)).tanh() / m.powi(5)
        })
        .sum::<f64>();
    let conductance = width * height.powi(3) / (12.0 * gas.dynamic_viscosity_pa_s * length)
        * (1.0 - 192.0 * height / (std::f64::consts::PI.powi(5) * width) * sum);
    let mut errors = Vec::new();
    for across in [2, 3, 4] {
        let (mesh, pos) = box_mesh([1, across, across], size);
        let boundary = boundaries(&mesh, &pos, size, true);
        let model =
            ThermoviscousFe::new(&mesh, &pos, gas, &boundary, limits(), &mut || false).unwrap();
        let field = model
            .prepare(1e-4, &mut || false)
            .unwrap()
            .response(&[C64::new(0.5, 0.0), C64::new(-0.5, 0.0)], &mut || false)
            .unwrap();
        assert_work(&field);
        let error = (field.inward_flow_m3_s[0].re - conductance).abs() / conductance;
        eprintln!(
            "LNS rectangular gap across={across} dofs={} aspect={} conductance_error={error:e}",
            model.ndof(),
            model.max_element_aspect_ratio()
        );
        errors.push(error);
    }
    assert!(errors[1] < errors[0] && errors[2] < errors[1], "{errors:?}");
    assert!(errors[2] < 0.015, "{errors:?}");
}

#[test]
fn full_lns_refuses_invalid_partition_material_budgets_and_cancellation() {
    let size = [0.002, 0.001, 0.001];
    let (mesh, pos) = box_mesh([1, 2, 2], size);
    let boundary = boundaries(&mesh, &pos, size, true);
    let gas = air();
    assert!(matches!(
        ThermoviscousFe::new(&mesh, &pos, gas, &boundary, limits(), &mut || true),
        Err(ThermoviscousError::Cancelled)
    ));
    let mut bad = boundary.clone();
    bad.no_slip_faces.pop();
    assert!(ThermoviscousFe::new(&mesh, &pos, gas, &bad, limits(), &mut || false).is_err());
    let mut bad_gas = gas;
    bad_gas.cp_j_kg_k = gas.gas_constant_j_kg_k();
    assert!(
        ThermoviscousFe::new(&mesh, &pos, bad_gas, &boundary, limits(), &mut || false).is_err()
    );
    let model = ThermoviscousFe::new(&mesh, &pos, gas, &boundary, limits(), &mut || false).unwrap();
    let mut cap = limits();
    cap.max_dofs = model.ndof() - 1;
    assert!(matches!(
        ThermoviscousFe::new(&mesh, &pos, gas, &boundary, cap, &mut || false),
        Err(ThermoviscousError::Budget {
            what: "mixed unknowns",
            ..
        })
    ));
    cap = limits();
    cap.max_dense_bytes = 48 * model.ndof() * model.ndof() - 1;
    assert!(matches!(
        ThermoviscousFe::new(&mesh, &pos, gas, &boundary, cap, &mut || false),
        Err(ThermoviscousError::Budget {
            what: "dense bytes",
            ..
        })
    ));
    cap = limits();
    cap.max_element_aspect_ratio = model.max_element_aspect_ratio() * (1.0 - 1e-9);
    assert!(matches!(
        ThermoviscousFe::new(&mesh, &pos, gas, &boundary, cap, &mut || false),
        Err(ThermoviscousError::Aspect { .. })
    ));
    assert!(
        model
            .prepare_with_heat_capacity(100.0, C64::new(1000.0, 1.0), &mut || false)
            .is_err()
    );
    assert!(model.prepare(0.0, &mut || false).is_err());
    assert!(model.prepare(100.0, &mut || true).is_err());
    assert!(
        model
            .scattering_with_heat_capacity(100.0, C64::new(1000.0, 0.0), &[0.0, 1.0], &mut || false)
            .is_err()
    );
    let prepared = model.prepare(100.0, &mut || false).unwrap();
    assert!(prepared.response(&[C64::ONE], &mut || false).is_err());
    assert!(
        prepared
            .response(&[C64::ONE, C64::ZERO], &mut || true)
            .is_err()
    );
}
