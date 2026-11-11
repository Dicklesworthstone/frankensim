//! Native element heating and availability, with independent polynomial moments.
use super::*;
use crate::acoustics::helmholtz::QuadraticGeometryControls;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn tetrahedron() -> (TetComplex, Vec<[f64; 3]>, ThermoviscousBoundary) {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let positions = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let boundary = ThermoviscousBoundary {
        no_slip_faces: vec![],
        isothermal_faces: vec![],
        ports: vec![BoundaryMode::uniform(&mesh.faces)],
    };
    (mesh, positions, boundary)
}

fn relative(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-12 * expected.abs().max(1e-30),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn affine_physical_gradients_separate_all_four_native_channels_and_mass() {
    let (mesh, positions, boundary) = tetrahedron();
    let mut medium = air();
    medium.bulk_viscosity_pa_s = 2.7e-5;
    let model = ThermoviscousFe::new(&mesh, &positions, medium, &boundary, limits(), &mut || {
        false
    })
    .unwrap();
    let frequency = 137.0;
    let cp = C64::new(970.0, -43.0);
    let prepared = model
        .prepare_with_heat_capacity(frequency, cp, &mut || false)
        .unwrap();
    let a = C64::new(0.7, -0.4);
    let b = C64::new(-0.2, 0.6);
    let mut velocity = vec![[C64::ZERO; 3]; model.space.ndof];
    let mut temperature = vec![C64::ZERO; model.space.ndof];
    for (i, &[x, y, z]) in positions.iter().enumerate() {
        velocity[i] = [
            a.scale(2.0 * x - y),
            a.scale(3.0 * y + z),
            a.scale(-x + 4.0 * z),
        ];
        temperature[i] = b.scale(1.0 + 2.0 * x - 3.0 * y + z);
    }
    let elements = prepared
        .element_dissipation(&velocity, &temperature, &mut || false)
        .unwrap();
    assert_eq!(elements.len(), 1);
    let power = elements[0];
    let volume = 1.0 / 6.0;
    relative(power.volume_m3, volume);
    relative(power.mass_kg, medium.density_kg_m3 * volume);
    // sym(grad v) has deviatoric diagonal [-1,0,1] and three
    // off-diagonal entries of magnitude 1/2; div(v)=9.
    relative(
        power.shear_heating_w,
        3.5 * medium.dynamic_viscosity_pa_s * volume * a.norm_sq(),
    );
    relative(
        power.bulk_heating_w,
        40.5 * medium.bulk_viscosity_pa_s * volume * a.norm_sq(),
    );
    relative(
        power.thermal_availability_loss_w,
        7.0 * medium.thermal_conductivity_w_m_k / medium.temperature_k * volume * b.norm_sq(),
    );
    relative(
        power.conductive_entropy_w_k,
        power.thermal_availability_loss_w / medium.temperature_k,
    );
    // Exact tetrahedron moments give mean (1+2x-3y+z)^2 = 17/10.
    relative(
        power.caloric_availability_loss_w,
        -0.5 * std::f64::consts::TAU * frequency * medium.density_kg_m3 * cp.im
            / medium.temperature_k
            * volume
            * 1.7
            * b.norm_sq(),
    );
    relative(
        power.viscous_heat_w(),
        power.shear_heating_w + power.bulk_heating_w,
    );
    assert!(power.disturbance_dissipation_w() > power.viscous_heat_w());

    // Global phase is immaterial; doubling peak amplitude multiplies all loss
    // channels by four while original geometry and source mass stay unchanged.
    let phase = C64::new(0.0, 2.0);
    let scaled_velocity = velocity
        .iter()
        .map(|v| v.map(|z| phase * z))
        .collect::<Vec<_>>();
    let scaled_temperature = temperature.iter().map(|&z| phase * z).collect::<Vec<_>>();
    let scaled = prepared
        .element_dissipation(&scaled_velocity, &scaled_temperature, &mut || false)
        .unwrap()[0];
    assert_eq!(scaled.volume_m3, power.volume_m3);
    assert_eq!(scaled.mass_kg, power.mass_kg);
    relative(scaled.viscous_heat_w(), 4.0 * power.viscous_heat_w());
    relative(
        scaled.thermal_availability_loss_w,
        4.0 * power.thermal_availability_loss_w,
    );
    relative(
        scaled.caloric_availability_loss_w,
        4.0 * power.caloric_availability_loss_w,
    );
}

#[test]
fn quadratic_edge_coefficients_recover_independent_tetrahedron_moments() {
    let (mesh, positions, boundary) = tetrahedron();
    let mut medium = air();
    medium.bulk_viscosity_pa_s = 2.7e-5;
    let model = ThermoviscousFe::new(&mesh, &positions, medium, &boundary, limits(), &mut || {
        false
    })
    .unwrap();
    let cp = C64::new(970.0, -43.0);
    let frequency = 137.0;
    let prepared = model
        .prepare_with_heat_capacity(frequency, cp, &mut || false)
        .unwrap();
    // Native edge bubbles are lambda_a*lambda_b, so the coefficient of x^2
    // is -(x_a-x_b)^2, not the midpoint's nodal value.
    let mut temperature = positions
        .iter()
        .map(|p| C64::new(p[0] * p[0], 0.0))
        .collect::<Vec<_>>();
    for &[a, b] in &mesh.edges {
        let difference = positions[a as usize][0] - positions[b as usize][0];
        temperature.push(C64::new(-difference * difference, 0.0));
    }
    let velocity = temperature
        .iter()
        .map(|&t| [t, C64::ZERO, C64::ZERO])
        .collect::<Vec<_>>();
    let power = prepared
        .element_dissipation(&velocity, &temperature, &mut || false)
        .unwrap()[0];
    let volume = 1.0 / 6.0;
    // On the unit tetrahedron, integral x^2=V/10 and integral x^4=V/35.
    relative(
        power.shear_heating_w,
        4.0 / 15.0 * medium.dynamic_viscosity_pa_s * volume,
    );
    relative(
        power.bulk_heating_w,
        medium.bulk_viscosity_pa_s * volume / 5.0,
    );
    relative(
        power.thermal_availability_loss_w,
        medium.thermal_conductivity_w_m_k * volume / (5.0 * medium.temperature_k),
    );
    relative(
        power.caloric_availability_loss_w,
        -0.5 * std::f64::consts::TAU * frequency * medium.density_kg_m3 * cp.im
            / medium.temperature_k
            * volume
            / 35.0,
    );
}

#[test]
fn rigid_spin_and_thermal_availability_do_not_create_viscous_heat() {
    let (mesh, positions, boundary) = tetrahedron();
    let mut medium = air();
    medium.bulk_viscosity_pa_s = 2.7e-5;
    let model = ThermoviscousFe::new(&mesh, &positions, medium, &boundary, limits(), &mut || {
        false
    })
    .unwrap();
    let prepared = model.prepare(137.0, &mut || false).unwrap();
    let mut velocity = vec![[C64::ZERO; 3]; model.space.ndof];
    let mut temperature = vec![C64::ZERO; model.space.ndof];
    for (i, &[x, y, z]) in positions.iter().enumerate() {
        velocity[i] = [C64::new(-y, 0.0), C64::new(x, 0.0), C64::ZERO];
        temperature[i] = C64::new(z, 0.0);
    }
    let power = prepared
        .element_dissipation(&velocity, &temperature, &mut || false)
        .unwrap()[0];
    assert_eq!(power.viscous_heat_w(), 0.0);
    assert_eq!(power.caloric_availability_loss_w, 0.0);
    relative(
        power.thermal_availability_loss_w,
        medium.thermal_conductivity_w_m_k / (12.0 * medium.temperature_k),
    );
    assert!(power.disturbance_dissipation_w() > 0.0);

    for (i, &position) in positions.iter().enumerate() {
        velocity[i] = position.map(|v| C64::new(v, 0.0));
    }
    let dilation = prepared
        .element_dissipation(&velocity, &temperature, &mut || false)
        .unwrap()[0];
    assert_eq!(dilation.shear_heating_w, 0.0);
    relative(dilation.bulk_heating_w, 0.75 * medium.bulk_viscosity_pa_s);
}

#[test]
fn original_curved_geometry_sets_actual_volume_mass_and_physical_gradients() {
    let (mesh, original, _) = tetrahedron();
    let a = 0.6;
    let map = |[x, y, z]: [f64; 3]| [x, y, z * (1.0 + a * x)];
    let positions = original.iter().copied().map(map).collect::<Vec<_>>();
    let midpoints = mesh
        .edges
        .iter()
        .map(|&[i, j]| {
            map(std::array::from_fn(|axis| {
                f64::midpoint(original[i as usize][axis], original[j as usize][axis])
            }))
        })
        .collect::<Vec<_>>();
    let medium = air();
    let thermo = ThermoviscousThermodynamics::ideal_gas(medium).unwrap();
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 19,
                kernel_id: 95,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let geometry = HelmholtzFe::new_quadratic(
            &mesh,
            &positions,
            &midpoints,
            2,
            7,
            Medium {
                density_kg_m3: medium.density_kg_m3,
                sound_speed_m_s: thermo.sound_speed_m_s(medium).unwrap(),
            },
            QuadraticGeometryControls::default(),
            &cx,
        )
        .unwrap();
        // A port needs a planar face (`boundary_geometry` refuses a nonaffine one, pinned by `curved.rs`): the two faces
        // this map bends are no-slip walls, the two that stay planar carry the port. The element observations below
        // take their fields as arguments and do not depend on the boundary data.
        let (planar, bent): (Vec<_>, Vec<_>) = geometry
            .boundary_faces()
            .into_iter()
            .partition(|&face| geometry.boundary_geometry(face).is_ok());
        assert_eq!((planar.len(), bent.len()), (2, 2));
        let boundary = ThermoviscousBoundary {
            no_slip_faces: bent,
            isothermal_faces: vec![],
            ports: vec![BoundaryMode::uniform(&planar)],
        };
        let model =
            ThermoviscousFe::from_geometry(geometry, medium, &boundary, limits(), &mut || false)
                .unwrap();
        let prepared = model.prepare(137.0, &mut || false).unwrap();
        let mut velocity = vec![[C64::ZERO; 3]; model.space.ndof];
        let mut temperature = vec![C64::ZERO; model.space.ndof];
        for (i, &[_, y, z]) in positions.iter().enumerate() {
            velocity[i][0] = C64::new(y, 0.0);
            temperature[i] = C64::new(z, 0.0);
        }
        for (edge, (&[i, j], midpoint)) in mesh.edges.iter().zip(&midpoints).enumerate() {
            temperature[mesh.vertex_count + edge] = C64::new(
                4.0 * (midpoint[2] - 0.5 * (positions[i as usize][2] + positions[j as usize][2])),
                0.0,
            );
        }
        let power = prepared
            .element_dissipation(&velocity, &temperature, &mut || false)
            .unwrap()[0];
        let volume = 1.0 / 6.0 + a / 24.0;
        relative(power.volume_m3, volume);
        relative(power.mass_kg, medium.density_kg_m3 * volume);
        relative(
            power.shear_heating_w,
            0.5 * medium.dynamic_viscosity_pa_s * volume,
        );
        relative(
            power.thermal_availability_loss_w,
            0.5 * medium.thermal_conductivity_w_m_k / medium.temperature_k * volume,
        );
    });
}

#[test]
fn solved_fields_match_element_sums_and_invalid_observations_are_atomic() {
    let size = [0.004, 0.001, 0.0012];
    let (mesh, positions) = box_mesh([1, 2, 2], size);
    let boundary = boundaries(&mesh, &positions, size, true);
    let mut medium = air();
    medium.bulk_viscosity_pa_s = 2.7e-5;
    let model = ThermoviscousFe::new(&mesh, &positions, medium, &boundary, limits(), &mut || {
        false
    })
    .unwrap();
    let cp = C64::new(970.0, -43.0);
    let prepared = model
        .prepare_with_heat_capacity(1700.0, cp, &mut || false)
        .unwrap();
    let loaded = model
        .prepare_scattering(1700.0, cp, &[5e8, 3e8], &mut || false)
        .unwrap();
    let drive = [C64::new(1.0, 0.2), C64::new(-0.3, 0.1)];
    let direct = prepared.response(&drive, &mut || false).unwrap();
    let incident = drive.map(|value| value.scale(1e-5));
    let waves = loaded.response(&incident, &mut || false).unwrap();
    for (field, elements) in [
        (
            &direct,
            prepared
                .element_dissipation(&direct.velocity_m_s, &direct.temperature_k, &mut || false)
                .unwrap(),
        ),
        (
            &waves,
            loaded
                .element_dissipation(&waves.velocity_m_s, &waves.temperature_k, &mut || false)
                .unwrap(),
        ),
    ] {
        assert_work(field);
        assert_eq!(elements.len(), mesh.tets.len());
        relative(
            elements.iter().map(|e| e.volume_m3).sum(),
            size.iter().product(),
        );
        relative(
            elements.iter().map(|e| e.mass_kg).sum(),
            medium.density_kg_m3 * size.iter().product::<f64>(),
        );
        relative(
            elements.iter().map(|e| e.shear_heating_w).sum(),
            field.power.viscous_w,
        );
        relative(
            elements.iter().map(|e| e.bulk_heating_w).sum(),
            field.power.bulk_w,
        );
        relative(
            elements.iter().map(|e| e.thermal_availability_loss_w).sum(),
            field.power.conductive_w,
        );
        relative(
            elements.iter().map(|e| e.caloric_availability_loss_w).sum(),
            field.power.caloric_w,
        );
        for (i, power) in elements.iter().enumerate() {
            relative(power.volume_m3, model.geometry.vol_signed[i].abs());
            assert!(power.viscous_heat_w() >= 0.0 && power.conductive_entropy_w_k >= 0.0);
        }
    }
    let expected = prepared
        .element_dissipation(&direct.velocity_m_s, &direct.temperature_k, &mut || false)
        .unwrap();
    assert!(matches!(
        prepared.element_dissipation(
            &direct.velocity_m_s[1..],
            &direct.temperature_k,
            &mut || false
        ),
        Err(ThermoviscousError::InvalidInput(_))
    ));
    assert!(matches!(
        prepared.element_dissipation(
            &direct.velocity_m_s,
            &direct.temperature_k[1..],
            &mut || false
        ),
        Err(ThermoviscousError::InvalidInput(_))
    ));
    let mut invalid_velocity = direct.velocity_m_s.clone();
    invalid_velocity[0][1].im = f64::NAN;
    assert!(matches!(
        prepared.element_dissipation(&invalid_velocity, &direct.temperature_k, &mut || false),
        Err(ThermoviscousError::NonFinite)
    ));
    let mut invalid_temperature = direct.temperature_k.clone();
    invalid_temperature[0].re = f64::INFINITY;
    assert!(matches!(
        loaded.element_dissipation(&direct.velocity_m_s, &invalid_temperature, &mut || false),
        Err(ThermoviscousError::NonFinite)
    ));
    let mut polls = 0;
    assert!(matches!(
        prepared.element_dissipation(&direct.velocity_m_s, &direct.temperature_k, &mut || {
            polls += 1;
            polls == 8
        }),
        Err(ThermoviscousError::Cancelled)
    ));
    assert_eq!(
        prepared
            .element_dissipation(&direct.velocity_m_s, &direct.temperature_k, &mut || false)
            .unwrap(),
        expected
    );
}
