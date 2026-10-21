//! Closed lathed exterior bodies carrying named aperture patches, for the
//! exterior observation of solved aperture volume flows.
//!
//! The body is the existing [`lathe_profile`] closed soup with axial stations
//! refined around every wall aperture. An aperture is the union of panels whose
//! centroid falls inside its footprint: a circle on the developed wall, or a
//! centred disc on the +x end cap. Its edge is therefore resolved only to the
//! panel scale, and that staircase is part of the panel-ladder error. Each
//! aperture's trace is a uniform outward normal velocity that integrates to
//! exactly one m^3/s over the SAME panel areas the solver uses. Every other
//! panel is rigid. No chimney, riser, key cup, pad or player is represented.
use crate::{BemError, panel3d::SpherePanels, radiation_bake::lathe_profile};
use fs_math::det;

/// A circular opening in the lathed wall. Azimuth is measured about +x from
/// the +y axis toward +z.
#[derive(Debug, Clone, PartialEq)]
pub struct WallAperture {
    pub name: String,
    pub axial_m: f64,
    pub azimuth_rad: f64,
    pub radius_m: f64,
}

/// Complete body declaration. The profile gives the OUTER radius at increasing
/// axial stations; both ends are capped, the -x end rigid.
#[derive(Debug, Clone, PartialEq)]
pub struct LathedApertureSpec {
    pub profile: Vec<(f64, f64)>,
    pub around: usize,
    pub axial_step_m: f64,
    pub aperture_step_m: f64,
    pub wall: Vec<WallAperture>,
    /// A centred disc on the +x end cap: name and radius.
    pub end: Option<(String, f64)>,
}

/// The closed surface and one unit-flow trace per aperture, in declaration
/// order: wall apertures, then the end aperture.
#[derive(Debug, Clone)]
pub struct LathedApertureBody {
    pub surface: SpherePanels,
    pub names: Vec<String>,
    pub traces: Vec<Vec<f64>>,
    /// Summed area of each aperture's panels [m^2].
    pub patch_areas_m2: Vec<f64>,
}

fn bad(name: &'static str, value: f64, requirement: &'static str) -> BemError {
    BemError::InvalidScalar {
        name,
        value,
        requirement,
    }
}

fn radius_at(profile: &[(f64, f64)], x: f64) -> f64 {
    let i = profile
        .windows(2)
        .position(|w| x <= w[1].0)
        .unwrap_or(profile.len() - 2);
    let (a, b) = (profile[i], profile[i + 1]);
    a.1 + (b.1 - a.1) * (x - a.0) / (b.0 - a.0)
}

/// Mesh the body and bind every aperture.
///
/// # Errors
/// A degenerate profile or step, an aperture outside the wall, an empty or
/// overlapping aperture, or a surface the panel constructor refuses.
pub fn lathed_aperture_body(spec: &LathedApertureSpec) -> Result<LathedApertureBody, BemError> {
    let profile = &spec.profile;
    if profile.len() < 2
        || spec.around < 8
        || !profile
            .windows(2)
            .all(|w| w[1].0 > w[0].0 && w[0].1 > 0.0 && w[1].1 > 0.0)
    {
        return Err(bad(
            "profile",
            profile.len() as f64,
            ">= 2 stations, increasing x, positive radii, >= 8 around",
        ));
    }
    for step in [spec.axial_step_m, spec.aperture_step_m] {
        if !(step.is_finite() && step > 0.0) {
            return Err(bad("axial step", step, "positive and finite"));
        }
    }
    let (x0, x1) = (profile[0].0, profile[profile.len() - 1].0);
    for a in &spec.wall {
        if a.name.is_empty()
            || !(a.radius_m.is_finite() && a.radius_m > 0.0 && a.azimuth_rad.is_finite())
            || !(a.axial_m - a.radius_m > x0 && a.axial_m + a.radius_m < x1)
        {
            return Err(bad(
                "wall aperture",
                a.axial_m,
                "named, finite, positive radius, inside the wall",
            ));
        }
    }
    // Axial stations: the fine step wherever a station interval touches an
    // aperture footprint widened by one fine step, the coarse step elsewhere.
    let zones: Vec<(f64, f64)> = spec
        .wall
        .iter()
        .map(|a| {
            (
                a.axial_m - a.radius_m - spec.aperture_step_m,
                a.axial_m + a.radius_m + spec.aperture_step_m,
            )
        })
        .collect();
    let mut xs = vec![x0];
    loop {
        let x = xs[xs.len() - 1];
        let fine = zones
            .iter()
            .any(|&(a, b)| x + spec.axial_step_m > a && x < b);
        let step = if fine {
            spec.aperture_step_m
        } else {
            spec.axial_step_m
        };
        let next = x + step;
        if next >= x1 - 0.25 * step {
            xs.push(x1);
            break;
        }
        xs.push(next);
    }
    // Equal-area faceting: each ring polygon encloses the circle's area.
    let angle = core::f64::consts::TAU / spec.around as f64;
    let facet = det::sqrt(angle / det::sin(angle));
    let stations: Vec<(f64, f64)> = xs
        .iter()
        .map(|&x| (x, facet * radius_at(profile, x)))
        .collect();
    let lathed = lathe_profile(&stations, spec.around)
        .map_err(|_| bad("profile", f64::NAN, "the lathe refused the stations"))?;
    let wall_panels = 2 * spec.around * (stations.len() - 1);
    let mouth_start = lathed.mouth_start;
    let surface = SpherePanels::from_triangles(lathed.triangles)?;

    let mut names = Vec::new();
    let mut masks: Vec<Vec<bool>> = Vec::new();
    for a in &spec.wall {
        let mask = surface
            .centroids()
            .iter()
            .enumerate()
            .map(|(i, c)| {
                if i >= wall_panels {
                    return false;
                }
                let mut d = det::atan2(c[2], c[1]) - a.azimuth_rad;
                d -= core::f64::consts::TAU * (d / core::f64::consts::TAU).round();
                let arc = radius_at(profile, c[0]) * d;
                let dx = c[0] - a.axial_m;
                dx * dx + arc * arc <= a.radius_m * a.radius_m
            })
            .collect();
        names.push(a.name.clone());
        masks.push(mask);
    }
    if let Some((name, radius)) = &spec.end {
        if name.is_empty() || !(radius.is_finite() && *radius > 0.0) {
            return Err(bad("end aperture", *radius, "named, positive radius"));
        }
        let mask = surface
            .centroids()
            .iter()
            .enumerate()
            .map(|(i, c)| i >= mouth_start && det::hypot(c[1], c[2]) <= *radius)
            .collect();
        names.push(name.clone());
        masks.push(mask);
    }
    let mut owner = vec![false; surface.areas().len()];
    let mut traces = Vec::with_capacity(masks.len());
    let mut patch_areas_m2 = Vec::with_capacity(masks.len());
    for mask in &masks {
        let area: f64 = mask
            .iter()
            .zip(surface.areas())
            .filter(|(m, _)| **m)
            .map(|(_, a)| a)
            .sum();
        if !(area > 0.0) {
            return Err(bad(
                "aperture",
                area,
                "selects at least one panel; refine the steps",
            ));
        }
        for (taken, &m) in owner.iter_mut().zip(mask) {
            if m && *taken {
                return Err(bad("aperture", area, "apertures must not share panels"));
            }
            *taken |= m;
        }
        traces.push(
            mask.iter()
                .map(|&m| if m { 1.0 / area } else { 0.0 })
                .collect(),
        );
        patch_areas_m2.push(area);
    }
    Ok(LathedApertureBody {
        surface,
        names,
        traces,
        patch_areas_m2,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec(around: usize, fine: f64) -> LathedApertureSpec {
        LathedApertureSpec {
            profile: vec![(0.0, 0.01), (0.3, 0.01)],
            around,
            axial_step_m: 0.006,
            aperture_step_m: fine,
            wall: vec![
                WallAperture {
                    name: "hole-a".into(),
                    axial_m: 0.1,
                    azimuth_rad: core::f64::consts::FRAC_PI_2,
                    radius_m: 0.004,
                },
                WallAperture {
                    name: "hole-b".into(),
                    axial_m: 0.2,
                    azimuth_rad: core::f64::consts::FRAC_PI_2,
                    radius_m: 0.004,
                },
            ],
            end: Some(("end".into(), 0.0095)),
        }
    }

    /// Unit flow, disjoint apertures, the right count, and patch areas that
    /// approach the declared footprints as the mesh refines.
    #[test]
    fn apertures_carry_unit_flow_and_converge_to_their_footprints() {
        let mut errors = Vec::new();
        for (around, fine) in [(16, 0.002), (32, 0.001)] {
            let body = lathed_aperture_body(&spec(around, fine)).unwrap();
            assert_eq!(body.names, ["hole-a", "hole-b", "end"]);
            for trace in &body.traces {
                let flow: f64 = trace
                    .iter()
                    .zip(body.surface.areas())
                    .map(|(v, a)| v * a)
                    .sum();
                assert!((flow - 1.0).abs() < 1e-12, "unit flow {flow}");
            }
            let hole = core::f64::consts::PI * 0.004 * 0.004;
            errors.push((body.patch_areas_m2[0] - hole).abs() / hole);
        }
        // Mirrored offline: 10.1 % then 1.7 % for these two rungs.
        assert!(errors[1] < errors[0] && errors[1] < 0.05, "{errors:?}");
    }

    #[test]
    fn overlapping_or_outside_apertures_refuse() {
        let mut s = spec(16, 0.002);
        s.wall[1].axial_m = 0.103;
        assert!(lathed_aperture_body(&s).is_err());
        let mut s = spec(16, 0.002);
        s.wall[0].axial_m = 0.002;
        assert!(lathed_aperture_body(&s).is_err());
    }
}
