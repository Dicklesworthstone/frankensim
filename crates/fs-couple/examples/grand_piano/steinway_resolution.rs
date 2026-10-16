//! Bounded rib-normal subdivision of the source drawing's existing sweep.
//! Coordinates stay in drawing units here. The caller converts the requested
//! SI normal spacing; all original feature stations are retained exactly.
//! This changes mechanical resolution, not acoustic quadrature or a mode cap.

pub const MAX_NODES: usize = 20_000;
pub const MAX_TRIANGLES: usize = 40_000;

/// Convert distance normal to parallel ribs to drawing-station increments.
/// For t = y + slope*x, the normal distance is dt*scale/sqrt(1+slope^2).
pub fn drawing_gap(normal_step_m: f64, metres_per_pixel: f64, slope: f64)
    -> Result<f64, String> {
    if !normal_step_m.is_finite() || normal_step_m <= 0.0
        || !metres_per_pixel.is_finite() || metres_per_pixel <= 0.0 || !slope.is_finite() {
        return Err("rib-normal step and drawing scale must be finite and positive".into());
    }
    let unit = metres_per_pixel / slope.hypot(1.0);
    let gap = normal_step_m / unit;
    if unit <= 0.0 || !unit.is_finite() || !gap.is_finite() || gap <= 0.0 {
        return Err("rib-normal step is not representable in the source drawing".into());
    }
    Ok(gap)
}

/// Split every supplied interval, preserving feature stations and stable order.
/// Admission precedes allocation, with at most divisions+4 vertices per row:
/// divisions+1 regular positions and the three existing structural curves.
/// No feature or requested interval is dropped to obtain a smaller mesh.
pub fn refine(stations: &[f64], divisions: usize, maximum_gap: f64)
    -> Result<Vec<f64>, String> {
    if !(4..=64).contains(&divisions) || stations.len() < 2
        || stations.iter().any(|t| !t.is_finite())
        || stations.windows(2).any(|p| p[0] >= p[1])
        || !maximum_gap.is_finite() || maximum_gap <= 0.0 {
        return Err("resolved mesh requires divisions 4..64, increasing finite stations and a positive finite gap".into());
    }
    let maximum_rows = (MAX_NODES / (divisions + 4))
        .min(MAX_TRIANGLES / (2 * (divisions + 4)) + 1);
    if stations.len() > maximum_rows {
        return Err("source stations exceed resolved mesh node/triangle budget".into());
    }
    let mut rows = stations.len();
    // Count before allocating either the refined list or any mesh vertices.
    for pair in stations.windows(2) {
        let pieces = ((pair[1] - pair[0]) / maximum_gap).ceil().max(1.0);
        if !pieces.is_finite() || pieces > (maximum_rows - rows + 1) as f64 {
            return Err("requested rib-normal resolution exceeds the 20000-node/40000-triangle mesh budget".into());
        }
        rows += pieces as usize - 1;
    }
    let mut result = Vec::with_capacity(rows);
    for pair in stations.windows(2) {
        result.push(pair[0]);
        let pieces = ((pair[1] - pair[0]) / maximum_gap).ceil().max(1.0) as usize;
        for index in 1..pieces {
            let value = pair[0] + (pair[1] - pair[0]) * index as f64 / pieces as f64;
            if value <= *result.last().expect("left feature retained") || value >= pair[1] {
                return Err("requested rib-normal resolution has indistinguishable stations".into());
            }
            result.push(value);
        }
    }
    result.push(*stations.last().expect("at least two stations"));
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn si_step_is_normal_to_ribs_not_vertical_or_along_them() {
        for slope in [0.0_f64, 0.8, -1.2, 4.0] {
            let scale = 0.002173_f64;
            let gap = drawing_gap(0.02, scale, slope).unwrap();
            // Independent projection on the unit normal (slope,1).
            let normal = [slope / (1.0+slope*slope).sqrt(), 1.0/(1.0+slope*slope).sqrt()];
            let separation = [0.0, gap*scale];
            let distance = normal[0]*separation[0]+normal[1]*separation[1];
            assert!((distance - 0.02).abs() < 2e-17);
        }
    }

    #[test]
    fn every_original_station_is_exact_and_every_interval_is_resolved() {
        let stations = [1.3, 9.7, 41.2, 55.0, 93.1];
        for gap in [100.0, 11.0, 3.0, 0.9] {
            let actual = refine(&stations, 8, gap).unwrap();
            assert!(stations.iter().all(|v| actual.iter().any(|x| x.to_bits()==v.to_bits())));
            assert!(actual.windows(2).all(|p| p[0] < p[1] && p[1]-p[0] <= gap*(1.0+1e-13)));
            assert_eq!(actual, refine(&stations, 8, gap).unwrap());
            assert_eq!(actual.first(), stations.first());
            assert_eq!(actual.last(), stations.last());
        }
        assert_eq!(refine(&stations, 8, 100.0).unwrap(), stations);
    }

    #[test]
    fn positive_gap_cannot_bypass_work_or_coordinate_admission() {
        for gap in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MIN_POSITIVE] {
            assert!(refine(&[1.0,2.0], 8, gap).is_err());
        }
        for divisions in [0,3,65,usize::MAX] {
            assert!(refine(&[1.0,2.0], divisions, 0.1).is_err());
        }
        assert!(refine(&[0.0,1.0,1.0],8,0.1).is_err());
        assert!(refine(&[f64::NEG_INFINITY,1.0],8,0.1).is_err());
        assert!(refine(&[-f64::MAX,f64::MAX],8,f64::MAX).is_err());
        assert!(drawing_gap(f64::MAX, 0.001, 0.8).is_err());
        assert!(drawing_gap(0.01, 0.0, 0.8).is_err());
    }

    #[test]
    fn exact_capacity_passes_and_one_more_row_refuses() {
        let limit = MAX_NODES / 8;
        let actual = refine(&[0.0, (limit-1) as f64],4,1.0).unwrap();
        assert_eq!(actual.len(),limit);
        assert!(refine(&[0.0,limit as f64],4,1.0).unwrap_err().contains("budget"));
        let crowded: Vec<_> = (0..=limit).map(|i|i as f64).collect();
        assert!(refine(&crowded,4,10.0).is_err());
    }
}
