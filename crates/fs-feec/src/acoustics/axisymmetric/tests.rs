//! G0/G1/G3: independent swept moments, cylinder/annulus fields and boundary work.
use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

const LIMITS: WeightedAssemblyLimits = WeightedAssemblyLimits {
    max_cells: 2048,
    max_dofs: 2048,
    max_triplets: 18432,
};
const RHO: f64 = 1.3;
const SPEED: f64 = 2.0;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 91,
                kernel_id: 21,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 256,
        maximum_payload_bytes: 32 * 1024 * 1024,
        maximum_relative_residual: 5e-12,
    }
}
fn fluid() -> AxisymmetricFluid {
    AxisymmetricFluid {
        density_kg_m3: RHO,
        compressibility_pa_inv: 1.0 / (RHO * SPEED * SPEED),
    }
}
fn rectangle(
    inner: f64,
    outer: f64,
    length: f64,
    nr: usize,
    nz: usize,
    reversed: bool,
) -> TriComplex2 {
    let mut vertices = Vec::new();
    for iz in 0..=nz {
        for ir in 0..=nr {
            vertices.push([
                inner + (outer - inner) * ir as f64 / nr as f64,
                length * iz as f64 / nz as f64,
            ]);
        }
    }
    let mut faces = Vec::new();
    for iz in 0..nz {
        for ir in 0..nr {
            let a = (iz * (nr + 1) + ir) as u32;
            let b = a + 1;
            let c = a + nr as u32 + 1;
            let d = c + 1;
            for mut face in [[a, b, d], [a, d, c]] {
                if reversed {
                    face.swap(1, 2);
                }
                faces.push(face);
            }
        }
    }
    TriComplex2::from_indexed_triangles(
        "axisymmetric-pressure-fixture",
        vertices,
        faces,
        Metric2::axisymmetric(core::f64::consts::TAU).unwrap(),
    )
    .unwrap()
}
fn edges(model: &AxisymmetricPressureFe<'_>, axis: usize, coordinate: f64) -> Vec<[u32; 2]> {
    model
        .boundary_edges()
        .iter()
        .filter(|edge| {
            edge.vertices
                .iter()
                .all(|&v| (model.mesh().vertices()[v as usize][axis] - coordinate).abs() < 1e-12)
        })
        .map(|edge| edge.vertices)
        .collect()
}
fn constant_drive(edges: &[[u32; 2]], velocity: C64) -> Vec<AxisymmetricVelocityTrace> {
    edges.iter().map(|&edge| (edge, [velocity; 2])).collect()
}
fn form(matrix: &Csr, left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .enumerate()
        .map(|(row, &x)| {
            let (columns, values) = matrix.row(row);
            x * columns
                .iter()
                .zip(values)
                .map(|(&col, &value)| value * right[col])
                .sum::<f64>()
        })
        .sum()
}
fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}
fn complex_close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:?}, expected={expected:?}, tolerance={tolerance:e}"
    );
}
fn relative_error(actual: &[C64], expected: &[C64]) -> f64 {
    let numerator = actual
        .iter()
        .zip(expected)
        .map(|(&a, &b)| (a - b).abs().powi(2))
        .sum::<f64>();
    let denominator = expected.iter().map(|v| v.abs().powi(2)).sum::<f64>();
    (numerator / denominator).sqrt()
}

#[test]
fn exact_cylinder_annulus_moments_keep_material_jumps_and_both_spatial_directions() {
    with_cx(&CancelGate::new(), |cx| {
        for inner in [0.0, 0.35] {
            for reversed in [false, true] {
                let outer = 1.2;
                let length = 0.8;
                let middle = 0.5 * (inner + outer);
                let mesh = rectangle(inner, outer, length, 2, 3, reversed);
                let choices = [
                    AxisymmetricFluid {
                        density_kg_m3: 1.1,
                        compressibility_pa_inv: 0.4,
                    },
                    AxisymmetricFluid {
                        density_kg_m3: 2.3,
                        compressibility_pa_inv: 0.7,
                    },
                ];
                let materials = mesh
                    .faces()
                    .iter()
                    .map(|face| {
                        let radius = face
                            .iter()
                            .map(|&v| mesh.vertices()[v as usize][0])
                            .sum::<f64>()
                            / 3.0;
                        if radius > middle {
                            choices[1]
                        } else {
                            choices[0]
                        }
                    })
                    .collect::<Vec<_>>();
                let model = AxisymmetricPressureFe::new(&mesh, &materials, LIMITS, cx).unwrap();
                let one = vec![1.0; model.ndof()];
                let radial = mesh.vertices().iter().map(|p| p[0]).collect::<Vec<_>>();
                let axial = mesh.vertices().iter().map(|p| p[1]).collect::<Vec<_>>();
                let mut mass = 0.0;
                let mut radial_mass = 0.0;
                let mut axial_mass = 0.0;
                let mut mixed_mass = 0.0;
                let mut gradient_energy = 0.0;
                for ((a, b), material) in
                    [(inner, middle), (middle, outer)].into_iter().zip(choices)
                {
                    let area = core::f64::consts::PI * (b * b - a * a);
                    let beta = material.compressibility_pa_inv;
                    mass += beta * area * length;
                    radial_mass +=
                        beta * core::f64::consts::PI * length * (b.powi(4) - a.powi(4)) / 2.0;
                    axial_mass += beta * area * length.powi(3) / 3.0;
                    mixed_mass +=
                        beta * core::f64::consts::PI * length.powi(2) * (b.powi(3) - a.powi(3))
                            / 3.0;
                    gradient_energy += area * length / material.density_kg_m3;
                }
                close(
                    model.volume_m3(),
                    core::f64::consts::PI * (outer * outer - inner * inner) * length,
                    model.volume_m3(),
                );
                close(form(model.mass(), &one, &one), mass, mass);
                close(
                    form(model.mass(), &radial, &radial),
                    radial_mass,
                    radial_mass,
                );
                close(form(model.mass(), &axial, &axial), axial_mass, axial_mass);
                close(form(model.mass(), &radial, &axial), mixed_mass, mixed_mass);
                close(
                    form(model.stiffness(), &radial, &radial),
                    gradient_energy,
                    gradient_energy,
                );
                close(
                    form(model.stiffness(), &axial, &axial),
                    gradient_energy,
                    gradient_energy,
                );
                close(
                    form(model.stiffness(), &radial, &axial),
                    0.0,
                    gradient_energy,
                );
                close(form(model.stiffness(), &one, &one), 0.0, gradient_energy);
                for row in 0..model.ndof() {
                    assert!(
                        model.mass().get(row, row) > 0.0,
                        "axis vertices must retain compression mass"
                    );
                    for col in 0..model.ndof() {
                        assert_eq!(
                            model.mass().get(row, col).to_bits(),
                            model.mass().get(col, row).to_bits()
                        );
                        assert_eq!(
                            model.stiffness().get(row, col).to_bits(),
                            model.stiffness().get(col, row).to_bits()
                        );
                    }
                }
                assert_eq!(
                    model.boundary_edges().len(),
                    if inner == 0.0 { 7 } else { 10 }
                );
                for edge in model.boundary_edges() {
                    let a = mesh.vertices()[edge.vertices[0] as usize];
                    let b = mesh.vertices()[edge.vertices[1] as usize];
                    let normal = if a[1].abs() < 1e-12 && b[1].abs() < 1e-12 {
                        [0.0, -1.0]
                    } else if (a[1] - length).abs() < 1e-12 && (b[1] - length).abs() < 1e-12 {
                        [0.0, 1.0]
                    } else if (a[0] - outer).abs() < 1e-12 && (b[0] - outer).abs() < 1e-12 {
                        [1.0, 0.0]
                    } else {
                        [-1.0, 0.0]
                    };
                    for axis in 0..2 {
                        close(edge.outward_normal_r_z[axis], normal[axis], 1.0);
                    }
                }
            }
        }
    });
}

#[test]
fn original_endcap_moments_distinguish_center_and_rim_without_extra_area_factors() {
    with_cx(&CancelGate::new(), |cx| {
        let radius = 0.7;
        let mesh = TriComplex2::from_indexed_triangles(
            "axisymmetric-endcap",
            vec![[0.0, 0.0], [radius, 0.0], [0.0, 0.8]],
            vec![[0, 1, 2]],
            Metric2::axisymmetric(core::f64::consts::TAU).unwrap(),
        )
        .unwrap();
        let model = AxisymmetricPressureFe::new(&mesh, &[fluid()], LIMITS, cx).unwrap();
        let edge = model.edge([0, 1]).unwrap();
        let area = core::f64::consts::PI * radius * radius;
        close(edge.surface_area_m2, area, area);
        close(edge.mass[0][0], area / 6.0, area);
        close(edge.mass[0][1], area / 6.0, area);
        close(edge.mass[1][1], area / 2.0, area);
        let velocities = [C64::new(0.3, -0.1), C64::new(-0.2, 0.7)];
        let expected = [
            (velocities[0] + velocities[1]).scale(area / 6.0),
            velocities[0].scale(area / 6.0) + velocities[1].scale(area / 2.0),
            C64::ZERO,
        ];
        for trace in [
            ([0, 1], velocities),
            ([1, 0], [velocities[1], velocities[0]]),
        ] {
            let load = model.velocity_load(&[trace], &BTreeSet::new(), cx).unwrap();
            for (&actual, &expected) in load.iter().zip(&expected) {
                complex_close(actual, expected, 3e-13 * area);
            }
        }
        let admittance = 0.13;
        let prepared = model
            .prepare(
                0.3,
                &[AxisymmetricAdmittance {
                    edges: vec![[1, 0]],
                    admittance_m_pa_s: admittance,
                }],
                &[],
                controls(),
                cx,
            )
            .unwrap();
        close(
            prepared.conductance.get(0, 0),
            admittance * area / 6.0,
            area,
        );
        close(
            prepared.conductance.get(0, 1),
            admittance * area / 6.0,
            area,
        );
        close(
            prepared.conductance.get(1, 1),
            admittance * area / 2.0,
            area,
        );
        let first = prepared.solve(&[([0, 1], velocities)], cx).unwrap();
        let reversed = prepared
            .solve(&[([1, 0], [velocities[1], velocities[0]])], cx)
            .unwrap();
        for (&a, &b) in first.pressure_pa().iter().zip(reversed.pressure_pa()) {
            complex_close(a, b, 2e-12 * b.abs().max(1.0));
        }
        close(
            first.input_power_w,
            first.boundary_dissipation_w,
            first.boundary_dissipation_w,
        );
        assert!(first.boundary_dissipation_w > 0.0);
        // The two velocity families remain linked by p=i*omega*psi, under the
        // same +i convention as the other full pressure-acoustic owners.
        let iw = C64::new(0.0, core::f64::consts::TAU * prepared.frequency_hz());
        for (&p, &psi) in first
            .pressure_pa()
            .iter()
            .zip(first.pressure_potential_pa_s())
        {
            complex_close(p, iw * psi, 2e-12 * p.abs().max(1.0));
        }
    });
}

#[test]
fn cylinder_and_annulus_duct_response_and_wall_work_converge_together() {
    with_cx(&CancelGate::new(), |cx| {
        for inner in [0.0, 0.1] {
            let radius = 0.35;
            let length = 0.9;
            let k = 1.1 / length;
            let omega = SPEED * k;
            let y = 0.35 / (RHO * SPEED);
            let velocity = C64::new(0.03, -0.017);
            let outlet_pressure =
                velocity / C64::new(y * (k * length).cos(), (k * length).sin() / (RHO * SPEED));
            let area = core::f64::consts::PI * (radius * radius - inner * inner);
            let expected_heat = 0.5 * y * area * outlet_pressure.abs().powi(2);
            let mut errors = Vec::new();
            let mut heat_errors = Vec::new();
            for (nr, nz) in [(1, 4), (2, 8), (4, 16)] {
                let mesh = rectangle(inner, radius, length, nr, nz, false);
                let model = AxisymmetricPressureFe::new(
                    &mesh,
                    &vec![fluid(); mesh.faces().len()],
                    LIMITS,
                    cx,
                )
                .unwrap();
                let inlet = edges(&model, 1, 0.0);
                let outlet = edges(&model, 1, length);
                let prepared = model
                    .prepare(
                        omega / core::f64::consts::TAU,
                        &[AxisymmetricAdmittance {
                            edges: outlet,
                            admittance_m_pa_s: y,
                        }],
                        &[],
                        controls(),
                        cx,
                    )
                    .unwrap();
                assert_eq!(prepared.free_dofs().len(), model.ndof());
                let field = prepared
                    .solve(&constant_drive(&inlet, velocity), cx)
                    .unwrap();
                // Independent 1-D continuum solution embedded in the complete
                // meridian space: p=P_L(cos(k(L-z))+i rho c Y sin(k(L-z))).
                let expected = mesh
                    .vertices()
                    .iter()
                    .map(|x| {
                        let phase = k * (length - x[1]);
                        outlet_pressure * C64::new(phase.cos(), RHO * SPEED * y * phase.sin())
                    })
                    .collect::<Vec<_>>();
                errors.push(relative_error(field.pressure_pa(), &expected));
                heat_errors.push((field.boundary_dissipation_w / expected_heat - 1.0).abs());
                assert!(field.boundary_dissipation_w > 0.0);
                close(
                    field.input_power_w,
                    field.boundary_dissipation_w,
                    expected_heat,
                );
                assert!(field.relative_residual < controls().maximum_relative_residual);
                let point = field.probe(0, [1.0 / 3.0; 3]).unwrap();
                assert!(point.particle_velocity_r_z_m_s[1].abs() > 0.0);
                assert!(point.active_intensity_r_z_w_m2[1] > 0.0);
            }
            assert!(
                errors[1] < 0.35 * errors[0] && errors[2] < 0.35 * errors[1],
                "{errors:?}"
            );
            assert!(
                heat_errors[1] < 0.35 * heat_errors[0] && heat_errors[2] < 0.35 * heat_errors[1],
                "{heat_errors:?}"
            );
            assert!(errors[2] < 8e-4 && heat_errors[2] < 8e-4);
        }
    });
}

// Independent convergent series for J0 and J1 on the test interval |x|<=0.9.
// These are continuum field oracles, not a numerical eigensolve of the FE pencil.
fn j0(x: f64) -> f64 {
    let mut term = 1.0;
    let mut sum = term;
    for n in 1..25 {
        term *= -x * x / (4.0 * (n * n) as f64);
        sum += term;
    }
    sum
}
fn j1(x: f64) -> f64 {
    let mut term = 0.5 * x;
    let mut sum = term;
    for n in 1..25 {
        term *= -x * x / (4.0 * (n * (n + 1)) as f64);
        sum += term;
    }
    sum
}

#[test]
fn radial_bessel_field_and_natural_axis_regularity_converge_without_axis_constraints() {
    with_cx(&CancelGate::new(), |cx| {
        let radius = 0.6;
        let length = 0.7;
        let k = 0.9 / radius;
        let omega = SPEED * k;
        let amplitude = C64::new(0.7, -0.2);
        // p=A J0(kr), v_r=-i A J1(kr)/(rho c), so positive INWARD motion
        // at r=R has the opposite sign. Rigid end caps admit this radial field.
        let velocity = C64::new(0.0, j1(k * radius) / (RHO * SPEED)) * amplitude;
        let mut errors = Vec::new();
        let mut axis_velocity = Vec::new();
        for n in [2, 4, 8] {
            let mesh = rectangle(0.0, radius, length, n, n, false);
            let model =
                AxisymmetricPressureFe::new(&mesh, &vec![fluid(); mesh.faces().len()], LIMITS, cx)
                    .unwrap();
            let outer = edges(&model, 0, radius);
            assert_eq!(outer.len(), n);
            assert!(
                edges(&model, 0, 0.0).is_empty(),
                "the axis must not be treated as a swept wall"
            );
            let prepared = model
                .prepare(omega / core::f64::consts::TAU, &[], &[], controls(), cx)
                .unwrap();
            assert_eq!(prepared.free_dofs().len(), model.ndof());
            let field = prepared
                .solve(&constant_drive(&outer, velocity), cx)
                .unwrap();
            let expected = mesh
                .vertices()
                .iter()
                .map(|x| amplitude.scale(j0(k * x[0])))
                .collect::<Vec<_>>();
            errors.push(relative_error(field.pressure_pa(), &expected));
            let axis = field.probe(0, [1.0, 0.0, 0.0]).unwrap();
            assert_eq!(axis.position_r_z_m[0], 0.0);
            assert!(axis.pressure_pa.abs() > 0.5);
            axis_velocity.push(axis.particle_velocity_r_z_m_s[0].abs());
            assert!(
                (field.pressure_pa()[0] - field.pressure_pa()[n]).abs() > 0.1,
                "a one-dimensional or uniform-pressure surrogate erases this radial field"
            );
            assert!(field.kinetic_energy_j > 0.0 && field.compressional_energy_j > 0.0);
            close(field.boundary_dissipation_w, 0.0, 1.0);
            close(field.input_power_w, 0.0, 1.0);
        }
        assert!(
            errors[1] < 0.35 * errors[0] && errors[2] < 0.35 * errors[1],
            "{errors:?}"
        );
        assert!(errors[2] < 1.5e-3);
        assert!(
            axis_velocity[1] < 0.55 * axis_velocity[0]
                && axis_velocity[2] < 0.55 * axis_velocity[1],
            "axis-adjacent P1 derivatives {axis_velocity:?}"
        );
    });
}

#[test]
fn constant_compression_is_retained_and_pressure_release_restricts_only_its_trace() {
    with_cx(&CancelGate::new(), |cx| {
        let radius = 0.35;
        let length = 0.9;
        let velocity = C64::new(0.03, -0.017);
        let mesh = rectangle(0.0, radius, length, 2, 8, false);
        let model =
            AxisymmetricPressureFe::new(&mesh, &vec![fluid(); mesh.faces().len()], LIMITS, cx)
                .unwrap();
        let inlet = edges(&model, 1, 0.0);
        let omega = 0.2;
        let closed = model
            .prepare(omega / core::f64::consts::TAU, &[], &[], controls(), cx)
            .unwrap();
        let field = closed.solve(&constant_drive(&inlet, velocity), cx).unwrap();
        // Testing continuity with the retained constant function gives exactly
        // i*omega*beta*V*mean(p)=Q_in, at any spatial resolution.
        let mut integral = C64::ZERO;
        let mut measure = 0.0;
        for (row, &p) in field.pressure_pa().iter().enumerate() {
            let weight = model.mass().row(row).1.iter().sum::<f64>();
            integral = integral + p.scale(weight);
            measure += weight;
        }
        let expected = C64::new(
            0.0,
            -1.0 / (omega * fluid().compressibility_pa_inv * length),
        ) * velocity;
        complex_close(
            integral.scale(1.0 / measure),
            expected,
            3e-8 * expected.abs(),
        );
        assert_eq!(closed.free_dofs().len(), model.ndof());
        assert!(field.kinetic_energy_j > 0.0 && field.compressional_energy_j > 0.0);

        let k = 0.6 / length;
        let omega = SPEED * k;
        let amplitude = C64::new(0.0, RHO * SPEED / (k * length).cos()) * velocity;
        let mut errors = Vec::new();
        for (nr, nz) in [(1, 4), (2, 8), (4, 16)] {
            let mesh = rectangle(0.0, radius, length, nr, nz, false);
            let model =
                AxisymmetricPressureFe::new(&mesh, &vec![fluid(); mesh.faces().len()], LIMITS, cx)
                    .unwrap();
            let inlet = edges(&model, 1, 0.0);
            let release = edges(&model, 1, length);
            let prepared = model
                .prepare(
                    omega / core::f64::consts::TAU,
                    &[],
                    &release,
                    controls(),
                    cx,
                )
                .unwrap();
            let field = prepared
                .solve(&constant_drive(&inlet, velocity), cx)
                .unwrap();
            assert_eq!(prepared.free_dofs().len(), model.ndof() - nr - 1);
            let expected = mesh
                .vertices()
                .iter()
                .map(|x| amplitude.scale((k * (length - x[1])).sin()))
                .collect::<Vec<_>>();
            errors.push(relative_error(field.pressure_pa(), &expected));
            for edge in &release {
                for &vertex in edge {
                    assert_eq!(field.pressure_pa()[vertex as usize], C64::ZERO);
                    assert_eq!(field.pressure_potential_pa_s()[vertex as usize], C64::ZERO);
                }
            }
            assert!(
                prepared.free_dofs().contains(&0),
                "release must not erase the rest of the rotation axis"
            );
            close(field.input_power_w, 0.0, 1.0);
            assert!(
                prepared
                    .solve(&constant_drive(&release, velocity), cx)
                    .is_err()
            );
            assert!(
                model
                    .prepare(
                        omega / core::f64::consts::TAU,
                        &[AxisymmetricAdmittance {
                            edges: release.clone(),
                            admittance_m_pa_s: 0.0,
                        }],
                        &release,
                        controls(),
                        cx
                    )
                    .is_err()
            );
        }
        assert!(
            errors[1] < 0.4 * errors[0] && errors[2] < 0.4 * errors[1],
            "{errors:?}"
        );
        assert!(errors[2] < 2e-4);
    });
}

#[test]
fn invalid_physics_boundary_data_budgets_and_cancellation_refuse_before_publication() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let mesh = rectangle(0.0, 0.6, 0.7, 1, 1, false);
        let materials = vec![fluid(); mesh.faces().len()];
        let model = AxisymmetricPressureFe::new(&mesh, &materials, LIMITS, cx).unwrap();
        let planar = TriComplex2::from_indexed_triangles(
            "axisymmetric-wrong-metric",
            mesh.vertices().to_vec(),
            mesh.faces().to_vec(),
            Metric2::planar(1.0).unwrap(),
        )
        .unwrap();
        assert!(AxisymmetricPressureFe::new(&planar, &materials, LIMITS, cx).is_err());
        assert!(AxisymmetricPressureFe::new(&mesh, &[], LIMITS, cx).is_err());
        let mut invalid_materials = materials.clone();
        invalid_materials[0].compressibility_pa_inv = 0.0;
        assert!(AxisymmetricPressureFe::new(&mesh, &invalid_materials, LIMITS, cx).is_err());
        assert!(
            AxisymmetricPressureFe::new(
                &mesh,
                &materials,
                WeightedAssemblyLimits {
                    max_dofs: 3,
                    ..LIMITS
                },
                cx
            )
            .is_err()
        );
        for edge in [[0, 2], [0, 3]] {
            assert!(
                model.prepare(0.3, &[], &[edge], controls(), cx).is_err(),
                "axis and internal edges cannot prescribe physical surface pressure"
            );
        }
        assert!(
            model
                .prepare(
                    0.3,
                    &[AxisymmetricAdmittance {
                        edges: vec![[0, 1]],
                        admittance_m_pa_s: -0.1,
                    }],
                    &[],
                    controls(),
                    cx
                )
                .is_err()
        );
        assert!(
            model
                .prepare(
                    0.3,
                    &[AxisymmetricAdmittance {
                        edges: vec![[0, 1], [1, 0]],
                        admittance_m_pa_s: 0.1,
                    }],
                    &[],
                    controls(),
                    cx
                )
                .is_err()
        );
        assert!(
            model
                .prepare(0.3, &[], &[[0, 1], [1, 3], [2, 3]], controls(), cx)
                .is_err()
        );
        assert!(matches!(
            model.prepare(
                0.3,
                &[],
                &[],
                HarmonicResponseControls {
                    maximum_coordinates: 3,
                    ..controls()
                },
                cx
            ),
            Err(AxisymmetricError::Harmonic(
                HarmonicResponseError::Budget { .. }
            ))
        ));
        assert!(matches!(
            model.prepare(
                0.3,
                &[],
                &[],
                HarmonicResponseControls {
                    maximum_payload_bytes: 1,
                    ..controls()
                },
                cx
            ),
            Err(AxisymmetricError::Harmonic(
                HarmonicResponseError::Budget { .. }
            ))
        ));
        let prepared = model.prepare(0.3, &[], &[], controls(), cx).unwrap();
        assert!(prepared.solve(&[([0, 2], [C64::ZERO; 2])], cx).is_err());
        assert!(
            prepared
                .solve(&[([0, 1], [C64::new(f64::NAN, 0.0); 2])], cx)
                .is_err()
        );
        assert!(
            prepared
                .solve(&[([0, 1], [C64::ONE; 2]), ([1, 0], [C64::ONE; 2]),], cx)
                .is_err()
        );
        let rest = prepared.solve(&[], cx).unwrap();
        assert!(rest.pressure_pa().iter().all(|p| *p == C64::ZERO));
        gate.request();
        assert!(matches!(
            prepared.solve(&[], cx),
            Err(AxisymmetricError::Cancelled)
        ));
        assert!(matches!(
            AxisymmetricPressureFe::new(&mesh, &materials, LIMITS, cx),
            Err(AxisymmetricError::Cancelled)
        ));
    });
}
