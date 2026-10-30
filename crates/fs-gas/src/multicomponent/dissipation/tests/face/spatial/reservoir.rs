use super::*;

#[test]
fn physical_reservoir_uses_original_bath_enthalpy_and_boundary_mechanical_work() {
    let mut op = prepared3(1e6);
    let left = input3(&[0.8, 0.2], [0.3, -0.4, 0.8], 1.5);
    let right = input3(&[0.2, 0.8], [1.3, 1.6, -0.2], 2.5);
    let gradients = [[6.0, -3.0, 0.0], [-2.0, 0.0, 1.5]];
    let mut interior = [0.0; 6];
    op.face_from_primitives(
        left,
        right,
        0.5,
        gradients,
        MaxwellStefanSolveV1::default(),
        &mut interior,
        || false,
    )
    .unwrap();
    for boundary_on_left in [true, false] {
        let bath = if boundary_on_left { left } else { right };
        let mut flux = [0.0; 6];
        op.reservoir_face_from_primitives(
            left,
            right,
            0.5,
            gradients,
            boundary_on_left,
            MaxwellStefanSolveV1::default(),
            &mut flux,
            || false,
        )
        .unwrap();
        assert_eq!(
            flux[..5].iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            interior[..5]
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>()
        );
        // Equal molecular masses and rho=1 make original bath T equal p.
        let enthalpy_difference = (2.5 - 4.0) * (bath.pressure - 0.001);
        close(
            flux[5],
            flux[0] * enthalpy_difference - 0.09 * 2.0
                + bath.normal_velocity * flux[2]
                + bath.transverse_velocity[0] * flux[3]
                + bath.transverse_velocity[1] * flux[4],
        );
    }
    // A stationary isothermal reservoir bears shear but receives zero stress
    // work, even when the adjacent cell has both transverse velocities.
    let moving = input3(&[1.0, 0.0], [0.0, 1.0, -2.0], 2.0);
    let stationary = input3(&[1.0, 0.0], [0.0; 3], 2.0);
    let mut spatial = [0.0; 6];
    op.reservoir_face_from_primitives(
        moving,
        stationary,
        0.25,
        [[0.0; 3]; 2],
        false,
        MaxwellStefanSolveV1::default(),
        &mut spatial,
        || false,
    )
    .unwrap();
    close(spatial[3], 0.07 / 0.25);
    close(spatial[4], -2.0 * 0.07 / 0.25);
    assert_eq!(spatial[5], 0.0);
    let mut planar = prepared(1e6);
    let mut plane = [0.0; 5];
    planar
        .reservoir_face_from_primitives(
            input(&[1.0, 0.0], [0.0, 1.0], 2.0),
            input(&[1.0, 0.0], [0.0; 2], 2.0),
            0.25,
            [0.0; 2],
            false,
            MaxwellStefanSolveV1::default(),
            &mut plane,
            || false,
        )
        .unwrap();
    close(plane[3], spatial[3]);
    assert_eq!(plane[4], 0.0);
}
