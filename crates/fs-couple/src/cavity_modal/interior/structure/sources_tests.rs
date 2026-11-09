//! Exact uniform compression, phasor convention and original source admission.

use super::*;
use fs_feec::acoustics::helmholtz::heterogeneous::{HarmonicLimits, MAX_VOLUME_SOURCES};
use std::cell::Cell;

#[test]
fn original_volume_injection_drives_complete_pressure_potential_with_exact_compression() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 2, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let calls = Cell::new(0);
    let volume = interior
        .volume_velocity_source(
            4,
            &mut |_, _| {
                calls.set(calls.get() + 1);
                Ok::<_, HeterogeneousError>(6.0)
            },
            &mut || false,
        )
        .unwrap();
    assert!((volume.integrated_shape() - 1.0).abs() < 2e-14);
    let sampled = calls.get();
    let mass = dense(&[2.0], 1);
    let stiffness = dense(&[7.0], 1);
    let damping = dense(&[0.3], 1);
    with_cx(|cx| {
        // A detached structural coordinate: the cavity is exactly rigid, and
        // the complete fluid constant-compression coordinate must remain live.
        let model = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: Some(&damping),
            },
            &[Vec::new()],
            &[],
            controls(),
            cx,
        )
        .unwrap();
        let frequency = 0.7;
        let omega = core::f64::consts::TAU * frequency;
        let prepared = model.prepare(frequency, cx).unwrap();
        let q = C64::new(0.02, -0.007);
        let response = prepared
            .solve_with_sources(&[C64::ZERO], &[], &[], &[(&volume, q)], cx)
            .unwrap();
        // -i*w*beta*V*p=Q, with V=1/6. This is not the +i pressure convention.
        let exact = C64::new(-q.im, q.re)
            .scale(6.0 * medium.density_kg_m3 * medium.sound_speed_m_s.powi(2) / omega);
        for (i, &p) in response.pressure_pa.iter().enumerate() {
            let expected = if i < 4 { exact } else { C64::ZERO };
            assert!((p - expected).abs() < 2e-10 * exact.abs());
            assert!(
                (C64::new(
                    omega * response.pressure_potential_pa_s[i].im,
                    -omega * response.pressure_potential_pa_s[i].re,
                ) - p)
                    .abs()
                    < 2e-12 * exact.abs()
            );
        }
        let expected_work = 0.5 * (exact.conj() * q).re;
        assert!(
            (response.average_acoustic_input_power_w - expected_work).abs()
                < 2e-11 * exact.abs() * q.abs()
        );
        assert_eq!(response.structural_displacement, vec![C64::ZERO]);
        let opposite = interior
            .prepare(frequency, &[], HarmonicLimits::default(), &mut || false)
            .unwrap()
            .solve_with_boundary_sources(&[], &[], &[(&volume, q.conj())], &mut || false)
            .unwrap();
        for (&minus, &plus) in response.pressure_pa.iter().zip(&opposite.pressure_pa) {
            assert!((minus - plus.conj()).abs() < 2e-10 * exact.abs());
        }

        let other = uniform(&geometry, medium);
        let foreign = other
            .volume_velocity_source(4, &mut |_, _| Ok::<_, HeterogeneousError>(6.0), &mut || {
                false
            })
            .unwrap();
        assert!(
            prepared
                .solve_with_sources(&[C64::ZERO], &[], &[], &[(&foreign, C64::ZERO)], cx)
                .is_err()
        );
        let boundary = interior
            .quadratic_velocity_source(&[([0, 1, 2], [0.0; 6])], &mut || false)
            .unwrap();
        assert!(
            interior
                .project_source_load(
                    &[[0, 1, 2]],
                    &[],
                    &[(&boundary, C64::ZERO)],
                    &[],
                    &mut || false
                )
                .is_err()
        );
        assert!(
            interior
                .project_source_load(
                    &[],
                    &[],
                    &[],
                    &[(&volume, C64::new(f64::NAN, 0.0))],
                    &mut || false
                )
                .is_err()
        );
        assert!(
            interior
                .project_source_load(
                    &[],
                    &[],
                    &[],
                    &vec![(&volume, C64::ZERO); MAX_VOLUME_SOURCES + 1],
                    &mut || false
                )
                .is_err()
        );
        assert!(
            interior
                .project_source_load(&[], &[], &[], &[(&volume, q)], &mut || true)
                .is_err()
        );
        assert_eq!(
            calls.get(),
            sampled,
            "prepared shapes are never resampled by the coupled factor"
        );
    });
}
