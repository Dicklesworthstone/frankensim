use super::*;
use crate::acoustics::thermoviscous::moving_wall::{MovingWall, modal_impedance};
use fs_modal::time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};

fn structure(omega: f64, damping: f64) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new_with_free_coordinates(
        48_000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: damping,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap()
}

fn liquid() -> (ThermoviscousMedium, ThermoviscousThermodynamics) {
    // Manufactured zero-expansion fluid separates mechanical storage from heat.
    (
        ThermoviscousMedium {
            density_kg_m3: 1000.0,
            pressure_pa: 101325.0,
            temperature_k: 300.0,
            cp_j_kg_k: 4000.0,
            dynamic_viscosity_pa_s: 0.001,
            thermal_conductivity_w_m_k: 0.6,
            bulk_viscosity_pa_s: 0.002,
        },
        ThermoviscousThermodynamics {
            isothermal_compressibility_pa_inv: 1.0 / 2.0e9,
            isobaric_expansion_per_k: 0.0,
        },
    )
}

fn closed<'c>(
    mesh: &'c TetComplex,
    pos: &[[f64; 3]],
    medium: ThermoviscousMedium,
    thermo: ThermoviscousThermodynamics,
    source: ModalAcousticTimeModel,
    shape: impl Fn([f64; 3]) -> [f64; 3],
) -> ThermoviscousFe<'c> {
    let geometry = HelmholtzFe::new(
        mesh,
        pos,
        2,
        Medium {
            density_kg_m3: medium.density_kg_m3,
            sound_speed_m_s: thermo.sound_speed_m_s(medium).unwrap(),
        },
    )
    .unwrap();
    let faces = geometry.boundary_faces();
    let space = SimplexSpace::new(mesh, 2);
    let mut phi = vec![[0.0; 3]; space.ndof];
    for face in &faces {
        for &v in face {
            phi[v as usize] = shape(pos[v as usize]);
        }
    }
    ThermoviscousFe::new_with_moving_wall(
        mesh,
        pos,
        medium,
        thermo,
        &ThermoviscousBoundary {
            no_slip_faces: faces.clone(),
            isothermal_faces: vec![],
            ports: vec![],
        },
        MovingWall {
            faces,
            velocity_shapes: vec![phi],
            structure: source,
        },
        limits(),
        &mut || false,
    )
    .unwrap()
}

fn assert_boundary_and_reaction(
    model: &ThermoviscousFe<'_>,
    field: &ThermoviscousField,
    force: C64,
    frequency: f64,
) {
    let wall = model.moving_wall().unwrap();
    let u = field.moving_wall.modal_velocity_m_sqrt_kg_per_s[0];
    let mut boundary_dofs = BTreeSet::new();
    for face in &wall.faces {
        for &v in face {
            boundary_dofs.insert(v as usize);
        }
        let face = canonical(*face);
        for a in 0..3 {
            for b in a + 1..3 {
                let edge = model
                    .space
                    .complex
                    .edges
                    .binary_search(&[face[a], face[b]])
                    .unwrap();
                boundary_dofs.insert(model.space.edge_off + edge);
            }
        }
    }
    for i in boundary_dofs {
        for c in 0..3 {
            let expected = u.scale(wall.velocity_shapes[0][i][c]);
            assert!((field.velocity_m_s[i][c] - expected).abs() < 1e-13 * u.abs().max(1e-30));
        }
    }
    let fluid = field.moving_wall.fluid_force_n_per_sqrt_kg[0];
    let structural = modal_impedance(
        &wall.structure.modes()[0],
        std::f64::consts::TAU * frequency,
    ) * u;
    assert!(
        (force + fluid - structural).abs() < 2e-9 * force.abs().max(structural.abs()).max(1e-30)
    );
    let fluid_losses = field.power.dissipated_w() - field.power.structural_damping_w;
    let port_work: f64 = field
        .pressure_drive_pa
        .iter()
        .zip(&field.inward_flow_m3_s)
        .map(|(&p, &q)| 0.5 * (p.conj() * q).re)
        .sum();
    let scale = field.power.dissipated_w().abs().max(1e-30);
    assert!((port_work + field.moving_wall.power_into_fluid_w - fluid_losses).abs() < 2e-8 * scale);
    assert_work(field);
}

#[test]
fn g1_normal_wall_load_has_actual_fluid_compliance_and_material_swap() {
    let size = [0.002, 0.001, 0.001];
    let (mesh, pos) = box_mesh([2, 2, 2], size);
    let (medium, thermo) = liquid();
    let amplitude = 100.0;
    let frequency = 10.0;
    let omega = std::f64::consts::TAU * frequency;
    let mut observed = Vec::new();
    for bulk_scale in [1.0, 4.0] {
        let eos = ThermoviscousThermodynamics {
            isothermal_compressibility_pa_inv: thermo.isothermal_compressibility_pa_inv
                / bulk_scale,
            ..thermo
        };
        let model = closed(&mesh, &pos, medium, eos, structure(300.0, 0.02), |p| {
            [amplitude * p[0] / size[0], 0.0, 0.0]
        });
        let field = model
            .prepare(frequency, &mut || false)
            .unwrap()
            .response_with_wall_forces(&[], &[C64::ONE], &mut || false)
            .unwrap();
        assert_boundary_and_reaction(&model, &field, C64::ONE, frequency);
        let u = field.moving_wall.modal_velocity_m_sqrt_kg_per_s[0];
        let fluid_z = -field.moving_wall.fluid_force_n_per_sqrt_kg[0] / u;
        // Uniform adiabatic compression: Kf=(dV/dq)^2/(beta_S*V).
        let expected = amplitude * amplitude * size[1] * size[2]
            / (size[0] * eos.isothermal_compressibility_pa_inv);
        assert!((-omega * fluid_z.im / expected - 1.0).abs() < 2e-5);
        observed.push(u.abs());
    }
    assert!(observed[1] < 0.251 * observed[0]);
}

#[test]
fn g1_accelerating_closed_volume_adds_the_actual_fluid_mass_to_structure() {
    let size = [0.002, 0.001, 0.001];
    let (mesh, pos) = box_mesh([2, 2, 2], size);
    let (medium, thermo) = liquid();
    let mass = medium.density_kg_m3 * size.iter().product::<f64>();
    let amplitude = 1.0 / det::sqrt(mass);
    let frequency = 30.0;
    let omega = std::f64::consts::TAU * frequency;
    let model = closed(&mesh, &pos, medium, thermo, structure(20.0, 0.05), |_| {
        [amplitude, 0.0, 0.0]
    });
    let field = model
        .prepare(frequency, &mut || false)
        .unwrap()
        .response_with_wall_forces(&[], &[C64::ONE], &mut || false)
        .unwrap();
    assert_boundary_and_reaction(&model, &field, C64::ONE, frequency);
    let u = field.moving_wall.modal_velocity_m_sqrt_kg_per_s[0];
    let z = -field.moving_wall.fluid_force_n_per_sqrt_kg[0] / u;
    // Wavelength >> box length: the liquid translates with the box. Its mass
    // in this mass-normalized trace is rho*V*amplitude^2 = 1.
    assert!((z.im / omega - 1.0).abs() < 2e-5);
    let unloaded =
        C64::ONE / modal_impedance(&model.moving_wall().unwrap().structure.modes()[0], omega);
    assert!(u.abs() < 0.51 * unloaded.abs());
}

#[test]
fn g1_tangential_stokes_shear_uses_full_no_slip_and_physical_viscosity() {
    let size = [0.002, 0.0002, 0.001];
    let (mesh, pos) = box_mesh([2, 2, 2], size);
    let (mut medium, thermo) = liquid();
    medium.dynamic_viscosity_pa_s = 0.1;
    medium.bulk_viscosity_pa_s = 0.2;
    let amplitude = 100.0;
    let frequency = 0.01;
    let model = closed(&mesh, &pos, medium, thermo, structure(100.0, 0.02), |p| {
        [amplitude * p[1] / size[1], 0.0, 0.0]
    });
    let field = model
        .prepare(frequency, &mut || false)
        .unwrap()
        .response_with_wall_forces(&[], &[C64::ONE], &mut || false)
        .unwrap();
    assert_boundary_and_reaction(&model, &field, C64::ONE, frequency);
    let u = field.moving_wall.modal_velocity_m_sqrt_kg_per_s[0];
    // Affine Couette/Stokes limit, dv_x/dy=amplitude*u/H. This is shear
    // dissipation even though the upper-wall velocity has zero normal component.
    let expected =
        0.5 * medium.dynamic_viscosity_pa_s * size.iter().product::<f64>() * amplitude * amplitude
            / (size[1] * size[1])
            * u.norm_sq();
    assert!((field.power.viscous_w / expected - 1.0).abs() < 2e-4);
    assert!(field.power.bulk_w < 1e-6 * expected);
    assert!(field.power.structural_damping_w > 0.0);
}

#[test]
fn g1_traction_and_solid_drives_are_reciprocal_and_shape_sign_is_physical() {
    let size = [0.004, 0.0008, 0.001];
    let (mesh, pos) = box_mesh([2, 2, 2], size);
    let boundary = boundaries(&mesh, &pos, size, false);
    let medium = air();
    let moving: Vec<_> = boundary
        .no_slip_faces
        .iter()
        .copied()
        .filter(|f| f.iter().all(|&v| pos[v as usize][1] == size[1]))
        .collect();
    let space = SimplexSpace::new(&mesh, 2);
    let build = |sign: f64| {
        let mut phi = vec![[0.0; 3]; space.ndof];
        for (i, p) in pos.iter().enumerate() {
            if p[1] == size[1] && p[0] > 0.0 && p[0] < size[0] && p[2] > 0.0 && p[2] < size[2] {
                phi[i][1] = sign * 100.0;
            }
        }
        ThermoviscousFe::new_with_moving_wall(
            &mesh,
            &pos,
            medium,
            ThermoviscousThermodynamics::ideal_gas(medium).unwrap(),
            &boundary,
            MovingWall {
                faces: moving.clone(),
                velocity_shapes: vec![phi],
                structure: structure(3000.0, 0.02),
            },
            limits(),
            &mut || false,
        )
        .unwrap()
    };
    let frequency = 700.0;
    let model = build(1.0);
    let prepared = model.prepare(frequency, &mut || false).unwrap();
    let pressure = prepared
        .response(&[C64::ONE, C64::ZERO], &mut || false)
        .unwrap();
    let force = prepared
        .response_with_wall_forces(&[C64::ZERO; 2], &[C64::ONE], &mut || false)
        .unwrap();
    assert_boundary_and_reaction(&model, &pressure, C64::ZERO, frequency);
    assert_boundary_and_reaction(&model, &force, C64::ONE, frequency);
    let coupling = pressure.moving_wall.modal_velocity_m_sqrt_kg_per_s[0];
    assert!(coupling.abs() > 0.0);
    assert!((coupling - force.inward_flow_m3_s[0]).abs() < 2e-9 * coupling.abs());
    let reversed = build(-1.0);
    let opposite = reversed
        .prepare(frequency, &mut || false)
        .unwrap()
        .response(&[C64::ONE, C64::ZERO], &mut || false)
        .unwrap();
    assert!(
        (coupling + opposite.moving_wall.modal_velocity_m_sqrt_kg_per_s[0]).abs()
            < 2e-10 * coupling.abs()
    );
    for (a, b) in pressure
        .inward_flow_m3_s
        .iter()
        .zip(&opposite.inward_flow_m3_s)
    {
        assert!((*a - *b).abs() < 2e-10 * a.abs().max(1e-30));
    }
    assert!(
        prepared
            .response_with_wall_forces(&[C64::ZERO; 2], &[], &mut || false)
            .is_err()
    );
    assert!(matches!(
        prepared.response_with_wall_forces(&[C64::ZERO; 2], &[C64::ONE], &mut || true),
        Err(ThermoviscousError::Cancelled)
    ));
}

#[test]
fn g0_moving_trace_refuses_interior_and_stationary_wall_motion() {
    let size = [0.004, 0.0008, 0.001];
    let (mesh, pos) = box_mesh([1, 1, 1], size);
    let boundary = boundaries(&mesh, &pos, size, false);
    let medium = air();
    let faces: Vec<_> = boundary
        .no_slip_faces
        .iter()
        .copied()
        .filter(|f| f.iter().all(|&v| pos[v as usize][1] == size[1]))
        .collect();
    let space = SimplexSpace::new(&mesh, 2);
    for point in [0usize, faces[0][0] as usize] {
        let mut phi = vec![[0.0; 3]; space.ndof];
        phi[point][0] = 1.0;
        let result = ThermoviscousFe::new_with_moving_wall(
            &mesh,
            &pos,
            medium,
            ThermoviscousThermodynamics::ideal_gas(medium).unwrap(),
            &boundary,
            MovingWall {
                faces: faces.clone(),
                velocity_shapes: vec![phi],
                structure: structure(100.0, 0.02),
            },
            limits(),
            &mut || false,
        );
        assert!(matches!(
            result,
            Err(ThermoviscousError::InvalidInput(
                "wall shape must vanish outside its moving trace and on stationary intersections"
            ))
        ));
    }
}
