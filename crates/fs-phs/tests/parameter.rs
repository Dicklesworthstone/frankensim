//! Original coefficient directions retain reciprocal maps and signed loss slopes.

use fs_phs::{PhsStructureParameterDirection as Direction, PortHamiltonian, QuadraticStorage};

fn source(p: f64) -> PortHamiltonian {
    PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0 + 0.2 * p, -1.0 - 0.2 * p, 0.0],
        vec![0.4 - 0.1 * p, 0.1 + 0.03 * p, 0.1 + 0.03 * p, 0.5 + 0.2 * p],
        vec![1.0 + 0.3 * p, -0.2 + 0.1 * p, 0.4 - 0.2 * p, 0.7 + 0.5 * p],
        Box::new(
            QuadraticStorage::new(
                vec![2.0 + 0.4 * p, 0.3 - 0.1 * p, 0.3 - 0.1 * p, 1.0 + 0.2 * p],
                2,
            )
            .unwrap(),
        ),
    )
    .unwrap()
}

#[test]
fn full_original_flow_and_reciprocal_output_match_independent_source_changes() {
    let system = source(0.0);
    let prepared = system
        .prepare_parameter_direction(Direction {
            interconnection: Some(&[0.0, 0.2, -0.2, 0.0]),
            dissipation: Some(&[-0.1, 0.03, 0.03, 0.2]),
            ports: Some(&[0.3, 0.1, -0.2, 0.5]),
        })
        .unwrap();
    let state = [0.2, -0.4];
    let effort = system.effort(&state);
    let mut flow = [0.0; 2];
    let mut output = [0.0; 2];
    prepared
        .vector_field_from_effort_into(
            &effort,
            &[0.12, -0.10],
            &[1.2, -0.7],
            &[0.3, 0.2],
            &mut flow,
        )
        .unwrap();
    prepared
        .output_from_effort_into(&effort, &[0.12, -0.10], &mut output)
        .unwrap();
    let observe = |p: f64| {
        let system = source(p);
        let mut field = [0.0; 2];
        system
            .vector_field_from_effort_into(
                &system.effort(&state),
                &[1.2 + 0.3 * p, -0.7 + 0.2 * p],
                &mut field,
            )
            .unwrap();
        (field, system.output(&state))
    };
    let step = 1e-4;
    let (fp, yp) = observe(step);
    let (fm, ym) = observe(-step);
    for i in 0..2 {
        assert!((flow[i] - (fp[i] - fm[i]) / (2.0 * step)).abs() < 2e-11);
        assert!((output[i] - (yp[i] - ym[i]) / (2.0 * step)).abs() < 2e-11);
    }
}

#[test]
fn coherent_direct_sources_and_changing_output_map_keep_the_weak_remainder() {
    let large = 2.0_f64.powi(54);
    let system = PortHamiltonian::new(
        1,
        3,
        vec![0.0],
        vec![1.0],
        vec![large, 1.0, large],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
    )
    .unwrap();
    let dr = [-large]; // A signed local slope, not an admitted new resistor.
    let prepared = system
        .prepare_parameter_direction(Direction {
            dissipation: Some(&dr),
            ..Default::default()
        })
        .unwrap();
    let mut flow = [0.0];
    prepared
        .vector_field_from_effort_into(&[1.0], &[large], &[0.0; 3], &[1.0, 1.0, -1.0], &mut flow)
        .unwrap();
    assert_eq!(flow, [1.0]);

    let identity = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let system = PortHamiltonian::new(
        3,
        1,
        vec![0.0; 9],
        identity.clone(),
        vec![1.0; 3],
        Box::new(QuadraticStorage::new(identity, 3).unwrap()),
    )
    .unwrap();
    let prepared = system
        .prepare_parameter_direction(Direction {
            ports: Some(&[-1.0, 0.0, -1.0]),
            ..Default::default()
        })
        .unwrap();
    let mut output = [0.0];
    prepared
        .output_from_effort_into(&[large, 1.0, large], &[large, 1.0, large], &mut output)
        .unwrap();
    assert_eq!(output, [1.0]);
}

#[test]
fn malformed_directions_refuse_without_imposing_positivity_on_signed_slopes() {
    let system = source(0.0);
    for direction in [
        Direction {
            interconnection: Some(&[0.0, 0.2, 0.2, 0.0]),
            ..Default::default()
        },
        Direction {
            dissipation: Some(&[0.0, 0.2, 0.0, 0.0]),
            ..Default::default()
        },
        Direction {
            ports: Some(&[0.0]),
            ..Default::default()
        },
        Direction {
            dissipation: Some(&[0.0, 0.0, 0.0, f64::NAN]),
            ..Default::default()
        },
    ] {
        assert!(system.prepare_parameter_direction(direction).is_err());
    }
    let prepared = system
        .prepare_parameter_direction(Direction {
            dissipation: Some(&[-1.0, 0.0, 0.0, 0.0]),
            ..Default::default()
        })
        .unwrap();
    let mut out = [7.0; 2];
    assert!(
        prepared
            .vector_field_from_effort_into(&[0.0], &[0.0; 2], &[0.0; 2], &[0.0; 2], &mut out)
            .is_err()
    );
    assert_eq!(out, [7.0; 2]);
    assert!(
        prepared
            .output_from_effort_into(&[0.0; 2], &[f64::NAN, 0.0], &mut out)
            .is_err()
    );
    assert_eq!(out, [7.0; 2]);
    prepared
        .output_from_effort_into(&[0.0; 2], &[0.0; 2], &mut out)
        .unwrap();
    assert_eq!(out, [0.0; 2]);
}
