use super::*;

fn film() -> ResistiveFilm {
    let gap = GapPort {
        reference_m: 0.0002,
        closure: vec![1.0, -1.0],
    };
    ResistiveFilm::new(
        vec![FilmCell {
            area_m2: 0.0001,
            gap: gap.clone(),
        }],
        vec![FilmChannel {
            from: 0,
            to: None,
            width_m: 0.01,
            length_m: 0.02,
            gap,
        }],
        2,
        0.01,
        FilmLimits {
            minimum_cell_gap_m: 0.00001,
            maximum_gap_m: 0.0005,
            maximum_pressure_pa: 10000.0,
        },
    )
    .unwrap()
}

#[test]
fn changed_viscosity_drives_the_same_pressure_force_work_and_tangent() {
    let old = film();
    let next = old.with_viscosity(0.004).unwrap();
    let q = [0.00001, 0.0];
    let v = [0.001, 0.0];
    let mut old_f = [0.0; 2];
    let mut new_f = [0.0; 2];
    let mut old_p = [0.0];
    let mut new_p = [0.0];
    let a = old.evaluate_into(&q, &v, &mut old_f, &mut old_p).unwrap();
    let b = next.evaluate_into(&q, &v, &mut new_f, &mut new_p).unwrap();
    // Independent Poiseuille balance with the actual, already-closed aperture.
    let g = 0.01 * 0.00019_f64.powi(3) / (12.0 * 0.004 * 0.02);
    let p = 0.0001 * 0.001 / g;
    assert!((new_p[0] - p).abs() < 1e-10);
    assert!((new_f[0] - 0.0001 * p).abs() < 1e-14);
    assert_eq!(new_f[1], -new_f[0]);
    assert!((b.dissipated_power_w - 0.4 * a.dissipated_power_w).abs() < 1e-17);
    assert_eq!(a.minimum_gap_m, b.minimum_gap_m);
    let mut da = [0.0; 2];
    let mut db = [0.0; 2];
    old.tangent_into(&q, &v, &[0.000001, 0.0], &[0.0001, 0.0], &mut da)
        .unwrap();
    next.tangent_into(&q, &v, &[0.000001, 0.0], &[0.0001, 0.0], &mut db)
        .unwrap();
    for i in 0..2 {
        assert!((db[i] - 0.4 * da[i]).abs() < 1e-14);
    }
}

#[test]
fn coefficient_recompilation_is_nonaccumulating_and_keeps_the_original_limits() {
    let original = film();
    let mut refreshed = original.clone();
    for _ in 0..64 {
        refreshed = refreshed
            .with_viscosity(0.013)
            .unwrap()
            .with_viscosity(0.01)
            .unwrap();
    }
    assert_eq!(refreshed.viscosity_pa_s(), original.viscosity_pa_s());
    assert_eq!(refreshed.factors, original.factors);
    let mut f = [42.0; 2];
    let mut p = [42.0];
    assert!(
        refreshed
            .evaluate_into(&[0.0002, 0.0], &[0.0; 2], &mut f, &mut p)
            .is_err()
    );
    assert_eq!(f, [42.0; 2]);
    assert_eq!(p, [42.0]);
    assert!(
        refreshed
            .evaluate_into(&[0.0; 2], &[1.0, 0.0], &mut f, &mut p)
            .is_err()
    );
    assert_eq!(f, [42.0; 2]);
    assert_eq!(p, [42.0]);
    for mu in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(original.with_viscosity(mu).is_err());
        assert_eq!(original.viscosity_pa_s(), 0.01);
    }
}
