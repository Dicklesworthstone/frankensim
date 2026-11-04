use super::*;
use asupersync::types::Budget;
use fs_exec::{CancelGate, ExecMode, StreamKey};

fn with_cx(operation: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 7,
                kernel_id: 0x4D41_534B,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx, &gate);
    });
}

fn clear() -> GoboTransmission {
    GoboTransmission::try_new(2.0, 2.0, 1, 1, vec![1.0]).unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5.0e-12 * expected.abs().max(1.0),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_shutter_cells_match_independent_rectangle_and_rotated_triangle_areas() {
    with_cx(|cx, _| {
        for (angle, offset, area) in [
            (0.0, 0.3, 2.6),
            (
                core::f64::consts::FRAC_PI_4,
                -core::f64::consts::FRAC_1_SQRT_2,
                0.5,
            ),
        ] {
            let mask =
                GoboMask::try_new(None, vec![FramingShutter::try_new(angle, offset).unwrap()])
                    .unwrap();
            for subdivisions in [[1, 1], [7, 5]] {
                let result = clear().masked(&mask, subdivisions, cx).unwrap();
                close(result.mean_transmission() * 4.0, area);
                assert!(
                    result
                        .data
                        .values
                        .iter()
                        .all(|value| (0.0..=1.0).contains(value))
                );
            }
        }
    });
}

#[test]
fn g1_polygonal_iris_area_and_shutter_intersection_have_analytic_references() {
    with_cx(|cx, _| {
        for blades in [3, 4, 8, MAX_IRIS_BLADES] {
            let radius = 0.4;
            let mask = GoboMask::try_new(
                Some(PolygonalIris::try_new(radius, blades, 0.37).unwrap()),
                Vec::new(),
            )
            .unwrap();
            let area =
                blades as f64 * radius * radius * (core::f64::consts::PI / blades as f64).tan();
            for subdivisions in [[1, 1], [7, 5]] {
                close(
                    clear()
                        .masked(&mask, subdivisions, cx)
                        .unwrap()
                        .mean_transmission()
                        * 4.0,
                    area,
                );
            }
        }
        // The four-sided iris is [-0.5,0.5]^2. Two blades leave the rectangle
        // [-0.5,0] x [-0.5,0.25], with independently known area 0.375 m².
        let mask = GoboMask::try_new(
            Some(PolygonalIris::try_new(0.5, 4, 0.0).unwrap()),
            vec![
                FramingShutter::try_new(0.0, 0.0).unwrap(),
                FramingShutter::try_new(core::f64::consts::FRAC_PI_2, 0.25).unwrap(),
            ],
        )
        .unwrap();
        close(
            clear()
                .masked(&mask, [9, 7], cx)
                .unwrap()
                .mean_transmission()
                * 4.0,
            0.375,
        );
        let closed = GoboMask::try_new(
            Some(PolygonalIris::try_new(0.0, 8, 0.0).unwrap()),
            Vec::new(),
        )
        .unwrap();
        assert_eq!(
            clear()
                .masked(&closed, [3, 4], cx)
                .unwrap()
                .mean_transmission(),
            0.0
        );
    });
}

#[test]
fn g3_masking_retains_authored_texels_calibration_and_identity_state() {
    with_cx(|cx, _| {
        let original = GoboTransmission::try_new(2.0, 2.0, 2, 1, vec![0.2, 0.8]).unwrap();
        let empty = GoboMask::try_new(None, Vec::new()).unwrap();
        let identity = original.masked(&empty, [1, 1], cx).unwrap();
        assert!(Arc::ptr_eq(&original.data, &identity.data));
        let mask =
            GoboMask::try_new(None, vec![FramingShutter::try_new(0.0, 0.5).unwrap()]).unwrap();
        // Entire left texel plus half the right: (0.2*2 + 0.8*1)/4 = 0.3.
        let first = original.masked(&mask, [3, 2], cx).unwrap();
        let replay = original.masked(&mask, [3, 2], cx).unwrap();
        close(first.mean_transmission(), 0.3);
        assert_eq!(first, replay);
        assert_eq!((first.width_m(), first.height_m()), (2.0, 2.0));
        assert_eq!((first.columns(), first.rows()), (6, 2));
        assert_eq!(original.data.values, [0.2, 0.8]);
        close(original.mean_transmission(), 0.5);
        close(
            original
                .masked(&empty, [5, 3], cx)
                .unwrap()
                .mean_transmission(),
            0.5,
        );
    });
}

#[test]
fn g3_subdivision_preserves_power_and_reduces_binary_edge_profile_error() {
    with_cx(|cx, _| {
        let boundary = 0.1;
        let mask =
            GoboMask::try_new(None, vec![FramingShutter::try_new(0.0, boundary).unwrap()]).unwrap();
        for (columns, expected_error) in [(1, 1.98), (8, 0.24), (32, 0.06)] {
            let map = clear().masked(&mask, [columns, 1], cx).unwrap();
            close(map.mean_transmission(), 0.55);
            // Integrate |represented transmission - ideal binary shutter|
            // against its exact geometric step, without sampling that step.
            let width = 2.0 / columns as f64;
            let mut error = 0.0;
            for (column, &value) in map.data.values.iter().enumerate() {
                let left = -1.0 + column as f64 * width;
                let open_width = (boundary - left).clamp(0.0, width);
                error += 2.0 * ((1.0 - value) * open_width + value * (width - open_width));
            }
            close(error, expected_error);
            if columns > 1 {
                assert_eq!(map.transmission_at(-0.5, 0.0), Ok(1.0));
                assert_eq!(map.transmission_at(0.75, 0.0), Ok(0.0));
            }
        }
    });
}

#[test]
fn g0_g4_mask_admission_underflow_and_cancellation_are_explicit() {
    for radius in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(PolygonalIris::try_new(radius, 8, 0.0).is_err());
    }
    for blades in [0, 2, MAX_IRIS_BLADES + 1] {
        assert!(PolygonalIris::try_new(1.0, blades, 0.0).is_err());
    }
    assert!(PolygonalIris::try_new(1.0, 8, f64::INFINITY).is_err());
    assert!(FramingShutter::try_new(f64::NAN, 0.0).is_err());
    assert!(FramingShutter::try_new(0.0, f64::INFINITY).is_err());
    let shutter = FramingShutter::try_new(0.0, 0.0).unwrap();
    assert!(GoboMask::try_new(None, vec![shutter; MAX_FRAMING_SHUTTERS + 1]).is_err());
    with_cx(|cx, gate| {
        let original = clear();
        let mask = GoboMask::try_new(None, vec![shutter]).unwrap();
        for subdivisions in [[0, 1], [usize::MAX, 2], [MAX_GOBO_TEXELS + 1, 1]] {
            assert!(matches!(
                original.masked(&mask, subdivisions, cx),
                Err(GoboMaskError::Source(_))
            ));
        }
        let tiny = GoboMask::try_new(
            Some(PolygonalIris::try_new(1.0e-300, 8, 0.0).unwrap()),
            Vec::new(),
        )
        .unwrap();
        assert!(matches!(
            original.masked(&tiny, [1, 1], cx),
            Err(GoboMaskError::Source(SourceError::NumericalRange))
        ));
        let narrow =
            GoboTransmission::try_new(f64::from_bits(2), 1.0e308, 1, 1, vec![1.0]).unwrap();
        assert!(matches!(
            narrow.masked(&mask, [4, 1], cx),
            Err(GoboMaskError::Source(SourceError::NumericalRange))
        ));
        let faint = GoboTransmission::try_new(2.0, 2.0, 1, 1, vec![f64::from_bits(1)]).unwrap();
        assert!(matches!(
            faint.masked(&mask, [1, 1], cx),
            Err(GoboMaskError::Source(SourceError::NumericalRange))
        ));
        gate.request();
        assert!(matches!(
            original.masked(&mask, [8, 8], cx),
            Err(GoboMaskError::Cancelled)
        ));
        assert_eq!(original.mean_transmission(), 1.0);
        assert_eq!(original.data.values, [1.0]);
    });
}
