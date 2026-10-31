//! Multiple rectangles retain every original reciprocal Rayleigh coordinate.
use super::*;

#[test]
fn every_rectangular_pair_and_terminal_pair_has_its_original_mutual_image() {
    let (first, terminals, options) = fixture();
    let second = DrivenMouthSurface { center_m: [0.012, 0.004], ..first };
    let mouths = [first, second];
    let gas = gas(293.15, 101_325.0);
    let gate = CancelGate::new_clock_free();
    let frequencies = [100.0, 1000.0, 4000.0];
    let surface = MixedSurface::new_many(&mouths, &terminals, &gas,
        [1.0, 4000.0], options, &gate).unwrap();
    let rows = surface.rows(&gas, &frequencies, &gate).unwrap();
    let centers = [first.center_m, second.center_m, terminals[0].center_m];
    for (&frequency, row) in frequencies.iter().zip(&rows) {
        assert_eq!(row.len(), 9);
        for i in 0..3 {
            for j in i + 1..3 {
                assert_eq!(row[3 * i + j], row[3 * j + i]);
                let distance = (centers[i][0] - centers[j][0])
                    .hypot(centers[i][1] - centers[j][1]);
                let phase = core::f64::consts::TAU * frequency * distance / gas.sound_speed;
                let expected = C64::new(det::sin(phase), det::cos(phase))
                    .scale(gas.density * frequency / distance);
                let relative = det::sqrt((row[3 * i + j] - expected).norm_sq()
                    / expected.norm_sq());
                assert!(relative < 0.002, "pair {i},{j} at {frequency}: {relative}");
            }
        }
        let q = [1.0, -0.4, -0.6];
        let power: f64 = (0..3).flat_map(|i| (0..3).map(move |j|
            q[i] * row[3 * i + j].re * q[j])).sum();
        let scale: f64 = row.iter().map(|z| z.re.abs()).sum();
        assert!(power >= -1e-12 * scale);
    }
    let reordered = MixedSurface::new_many(&[second, first], &terminals, &gas,
        [1.0, 4000.0], options, &gate).unwrap()
        .rows(&gas, &frequencies, &gate).unwrap();
    let permutation = [1, 0, 2];
    for (a, b) in rows.iter().zip(&reordered) {
        for i in 0..3 { for j in 0..3 {
            let expected = a[3 * permutation[i] + permutation[j]];
            let actual = b[3 * i + j];
            assert!(det::sqrt((actual - expected).norm_sq())
                <= 1e-12 * det::sqrt(expected.norm_sq()).max(1e-30));
        }}
    }
}

#[test]
fn empty_overlap_quadrature_and_singleton_contracts_remain_explicit() {
    let (mouth, terminals, options) = fixture();
    let gas = gas(293.15, 101_325.0);
    let gate = CancelGate::new_clock_free();
    let build = |mouths: &[DrivenMouthSurface]| MixedSurface::new_many(
        mouths, &terminals, &gas, [1.0, 4000.0], options, &gate);
    assert!(build(&[]).is_err());
    assert!(build(&[mouth, mouth]).is_err());
    assert!(build(&[mouth, DrivenMouthSurface { center_m: [0.0; 2], ..mouth }]).is_err());
    assert!(build(&[mouth, DrivenMouthSurface {
        center_m: [0.012, 0.004], terminal_angular_points: 64, ..mouth
    }]).is_err());
    let frequencies = [100.0, 1000.0, 4000.0];
    assert_eq!(build(&[mouth]).unwrap().rows(&gas, &frequencies, &gate).unwrap(),
        MixedSurface::new(mouth, &terminals, &gas, [1.0, 4000.0], options, &gate)
            .unwrap().rows(&gas, &frequencies, &gate).unwrap());
}
