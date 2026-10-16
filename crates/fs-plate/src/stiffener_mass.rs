//! Inertia of the existing two-node Hermite stiffener field.
//! The coordinate order is `(w1, dw/ds1, w2, dw/ds2)`, with both slopes
//! measured along the same directed segment. The translational option exactly
//! integrates `rho*A*w_dot(s)^2`. The eccentric option also retains the supplied
//! bending moment of area and motion of the beam's offset centroid. Saint-Venant
//! torsion J is not a polar mass moment and is never used as one here.

/// Mass law for the existing Euler–Bernoulli beams.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum StiffenerMass {
    /// Half the segment mass at each endpoint displacement; the original law.
    #[default]
    Lumped,
    /// Exact mass integral of the cubic Hermite displacement used by bending.
    ConsistentHermite,
    /// Hermite translation, bending rotary inertia and eccentric-centroid
    /// motion from the supplied area, bending inertia and centroid offset.
    /// Torsional section rotary inertia is not inferred from Saint-Venant J.
    ConsistentEccentric,
}

/// Exact `integral rho*A*N(s)^T*N(s) ds`. Zero line density admits a
/// deliberately massless brace, as in the existing stiffness-isolation tests.
pub(crate) fn hermite_mass(length: f64, line_density: f64) -> Option<[[f64; 4]; 4]> {
    if !length.is_finite() || length <= 0.0 || !line_density.is_finite() || line_density < 0.0 {
        return None;
    }
    let l = length;
    let l2 = l * l;
    let scale = line_density * l / 420.0;
    let coefficients = [
        [156.0, 22.0 * l, 54.0, -13.0 * l],
        [22.0 * l, 4.0 * l2, 13.0 * l, -3.0 * l2],
        [54.0, 13.0 * l, 156.0, -22.0 * l],
        [-13.0 * l, -3.0 * l2, -22.0 * l, 4.0 * l2],
    ];
    let mass = coefficients.map(|row| row.map(|value| scale * value));
    mass.iter()
        .flatten()
        .all(|value| value.is_finite())
        .then_some(mass)
}

/// Additional inertia beyond `hermite_mass`. The along-beam block acts on
/// `(w1,slope1,w2,slope2)`; the across-beam block acts on the two nodal cross
/// slopes interpolated linearly, as in the existing torsional stiffness.
pub(crate) struct EccentricInertia {
    pub along: [[f64; 4]; 4],
    pub across: [[f64; 2]; 2],
}

/// Exact integrals `rho*(I+A*e^2)*N'(s)^T*N'(s)` and
/// `rho*A*e^2*L(s)^T*L(s)`. With the plate's physical rotations and a vertical
/// centroid arm e, the latter two A*e^2 contributions are centroid translation,
/// not a second copy of transverse mass. I is the supplied centroidal bending
/// moment of area; no missing torsional polar inertia is invented.
pub(crate) fn eccentric_inertia(
    length: f64, density: f64, area: f64, inertia: f64, eccentricity: f64,
) -> Option<EccentricInertia> {
    if !length.is_finite() || length <= 0.0 || !density.is_finite() || density < 0.0
        || !area.is_finite() || area <= 0.0 || !inertia.is_finite() || inertia < 0.0
        || !eccentricity.is_finite() {
        return None;
    }
    let offset_inertia = area * eccentricity * eccentricity;
    let line_rotary = density * (inertia + offset_inertia);
    let line_offset = density * offset_inertia;
    if !line_rotary.is_finite() || !line_offset.is_finite() { return None; }
    let l = length;
    let l2 = l * l;
    let derivative_scale = line_rotary / (30.0 * l);
    let derivative = [
        [36.0, 3.0*l, -36.0, 3.0*l],
        [3.0*l, 4.0*l2, -3.0*l, -l2],
        [-36.0, -3.0*l, 36.0, -3.0*l],
        [3.0*l, -l2, -3.0*l, 4.0*l2],
    ];
    let along = derivative.map(|row| row.map(|value| derivative_scale * value));
    let cross_scale = line_offset * l / 6.0;
    let across = [[2.0*cross_scale, cross_scale], [cross_scale, 2.0*cross_scale]];
    along.iter().flatten().chain(across.iter().flatten()).all(|v| v.is_finite())
        .then_some(EccentricInertia { along, across })
}

#[cfg(test)]
mod tests {
    use super::{eccentric_inertia, hermite_mass};

    fn bilinear(m: &[[f64; 4]; 4], a: [f64; 4], b: [f64; 4]) -> f64 {
        (0..4)
            .flat_map(|i| (0..4).map(move |j| a[i] * m[i][j] * b[j]))
            .sum()
    }

    // Independent polynomial endpoints, slopes and continuum product integral.
    // These do not evaluate the element shape functions or repeat its matrix.
    fn endpoints(p: [f64; 4], length: f64) -> [f64; 4] {
        let value = p[0] + length * (p[1] + length * (p[2] + length * p[3]));
        let slope = p[1] + length * (2.0 * p[2] + 3.0 * length * p[3]);
        [p[0], p[1], value, slope]
    }

    fn integral(a: [f64; 4], b: [f64; 4], length: f64) -> f64 {
        (0..4)
            .flat_map(|i| {
                (0..4).map(move |j| {
                    a[i] * b[j] * length.powi((i + j + 1) as i32) / (i + j + 1) as f64
                })
            })
            .sum()
    }

    fn derivative(p: [f64; 4]) -> [f64; 4] { [p[1], 2.0*p[2], 3.0*p[3], 0.0] }
    fn cross_bilinear(m: &[[f64; 2]; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
        (0..2).flat_map(|i| (0..2).map(move |j| a[i]*m[i][j]*b[j])).sum()
    }

    #[test]
    fn eccentric_inertia_integrates_cubic_bending_and_linear_cross_slope_products() {
        let fields = [[1.0,0.0,0.0,0.0], [0.0,1.0,0.0,0.0],
            [0.0,0.0,1.0,0.0], [0.0,0.0,0.0,1.0], [0.7,-0.4,0.9,0.2]];
        let (density, area, inertia, offset) = (530.0, 0.0003, 2.3e-8, -0.025);
        for length in [0.03, 0.7, 2.0] {
            let m = eccentric_inertia(length, density, area, inertia, offset).unwrap();
            let line_rotary = density * (inertia + area*offset*offset);
            for a in fields { for b in fields {
                let actual = bilinear(&m.along, endpoints(a,length), endpoints(b,length));
                let expected = line_rotary*integral(derivative(a),derivative(b),length);
                assert!((actual-expected).abs() < 4e-13*line_rotary.max(expected.abs()));
            } }
            for a in [[1.0,0.0,0.0,0.0], [0.0,1.0,0.0,0.0], [-0.4,0.7,0.0,0.0]] {
                for b in [[1.0,0.0,0.0,0.0], [0.0,1.0,0.0,0.0], [0.2,-0.9,0.0,0.0]] {
                    let actual = cross_bilinear(&m.across,[a[0],a[0]+a[1]*length],
                        [b[0],b[0]+b[1]*length]);
                    let expected = density*area*offset*offset*integral(a,b,length);
                    assert!((actual-expected).abs() < 3e-14*line_rotary.max(expected.abs()));
                }
            }
            assert_eq!(bilinear(&m.along,[1.0,0.0,1.0,0.0],[1.0,0.0,1.0,0.0]),0.0);
        }
    }

    #[test]
    fn eccentric_rigid_tilt_has_the_supplied_physical_moment_and_translation_mass() {
        let (length, density, area, inertia, offset) = (0.73, 480.0, 0.0004, 3e-8, 0.031);
        let m = eccentric_inertia(length,density,area,inertia,offset).unwrap();
        let translation = hermite_mass(length,density*area).unwrap();
        let (shift, slope, cross) = (0.4, -0.7, 0.2);
        let q = [shift,slope,shift+slope*length,slope];
        let actual = bilinear(&translation,q,q) + bilinear(&m.along,q,q)
            + cross_bilinear(&m.across,[cross;2],[cross;2]);
        let expected = density*area*integral([shift,slope,0.0,0.0],[shift,slope,0.0,0.0],length)
            + density*length*(inertia*slope*slope + area*offset*offset*(slope*slope+cross*cross));
        assert!((actual-expected).abs() < 2e-14*expected);
        let rigid = [1.0,0.0,1.0,0.0];
        assert!((bilinear(&translation,rigid,rigid)+bilinear(&m.along,rigid,rigid)
            -density*area*length).abs() < 2e-14*density*area*length);
        let centered = eccentric_inertia(length,density,area,inertia,0.0).unwrap();
        assert_eq!(centered.across,[[0.0;2];2]);
        assert!((bilinear(&centered.along,q,q)-density*inertia*length*slope*slope).abs()
            < 2e-14*density*inertia*length);
        let absent = eccentric_inertia(length,density,area,0.0,0.0).unwrap();
        assert!(absent.along.iter().flatten().all(|&v|v==0.0));
        assert_eq!(absent.across,[[0.0;2];2]);
    }

    #[test]
    fn eccentric_inertia_is_direction_independent_and_exact_under_subdivision() {
        let (length,density,area,inertia,offset) = (1.7,510.0,0.0003,2e-8,0.026);
        let p = [0.4,-0.2,0.7,-0.1];let cross = [0.2,-0.4,0.0,0.0];
        let full = eccentric_inertia(length,density,area,inertia,offset).unwrap();
        let q = endpoints(p,length);let qc = [cross[0],cross[0]+length*cross[1]];
        let energy = bilinear(&full.along,q,q)+cross_bilinear(&full.across,qc,qc);
        let reverse = [q[2],-q[3],q[0],-q[1]];let rc = [-qc[1],-qc[0]];
        assert!((bilinear(&full.along,reverse,reverse)+cross_bilinear(&full.across,rc,rc)-energy).abs()
            < 2e-14*energy);
        let opposite = eccentric_inertia(length,density,area,inertia,-offset).unwrap();
        assert_eq!(full.along,opposite.along);assert_eq!(full.across,opposite.across);
        for i in 0..4 {for j in 0..4 {assert_eq!(full.along[i][j],full.along[j][i]);}}
        for segments in [1,2,7,16] {
            let h = length/f64::from(segments);
            let local = eccentric_inertia(h,density,area,inertia,offset).unwrap();
            let mut actual = 0.0;
            for segment in 0..segments {
                let x = f64::from(segment)*h;
                let shifted = [p[0]+x*(p[1]+x*(p[2]+x*p[3])),p[1]+x*(2.0*p[2]+3.0*x*p[3]),
                    p[2]+3.0*x*p[3],p[3]];
                let q = endpoints(shifted,h);
                let qc = [cross[0]+cross[1]*x,cross[0]+cross[1]*(x+h)];
                actual += bilinear(&local.along,q,q)+cross_bilinear(&local.across,qc,qc);
            }
            assert!((actual-energy).abs() < 3e-13*energy);
        }
    }

    #[test]
    fn eccentric_inertia_refuses_nonphysical_or_overflowing_inputs() {
        for [length,density,area,inertia,offset] in [
            [0.0,500.0,0.001,1e-8,0.02], [1.0,-1.0,0.001,1e-8,0.02],
            [1.0,500.0,0.0,1e-8,0.02], [1.0,500.0,0.001,-1e-8,0.02],
            [f64::NAN,500.0,0.001,1e-8,0.02], [1.0,f64::INFINITY,0.001,1e-8,0.02],
            [1.0,500.0,f64::NAN,1e-8,0.02], [1.0,500.0,0.001,f64::INFINITY,0.02],
            [1.0,500.0,0.001,1e-8,f64::NAN], [1e200,500.0,0.001,1e-8,0.02],
            [1.0,500.0,0.001,1e-8,1e200],
        ] { assert!(eccentric_inertia(length,density,area,inertia,offset).is_none()); }
        let zero = eccentric_inertia(0.7,0.0,0.001,1e-8,0.02).unwrap();
        assert!(zero.along.iter().flatten().all(|&v|v==0.0));
        assert_eq!(zero.across,[[0.0;2];2]);
    }

    #[test]
    fn consistent_hermite_mass_integrates_every_cubic_velocity_product() {
        let fields = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.7, -1.1, 0.4, 0.9],
        ];
        for length in [0.03, 0.7, 2.0] {
            let density = 3.7;
            let m = hermite_mass(length, density).unwrap();
            for a in fields {
                for b in fields {
                    let actual = bilinear(&m, endpoints(a, length), endpoints(b, length));
                    let expected = density * integral(a, b, length);
                    assert!((actual - expected).abs() < 3e-13 * expected.abs().max(1.0));
                }
                assert!(bilinear(&m, endpoints(a, length), endpoints(a, length)) > 0.0);
            }
            assert!(
                (bilinear(&m, [1.0, 0.0, 1.0, 0.0], [1.0, 0.0, 1.0, 0.0]) - density * length).abs()
                    < 2e-14
            );
        }
    }

    #[test]
    fn consistent_hermite_mass_preserves_energy_when_beam_direction_is_reversed() {
        let q = [0.3, 0.8, -1.2, 0.4];
        let reversed = [q[2], -q[3], q[0], -q[1]];
        let m = hermite_mass(0.73, 4.2).unwrap();
        assert!((bilinear(&m, q, q) - bilinear(&m, reversed, reversed)).abs() < 1e-14);
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(m[i][j], m[j][i]);
            }
        }
    }

    #[test]
    fn consistent_hermite_mass_preserves_polynomial_energy_under_subdivision() {
        let p = [0.4, -0.2, 0.7, -0.1];
        let length = 1.7;
        let density = 2.3;
        let expected = density * integral(p, p, length);
        for segments in [1, 2, 7, 16] {
            let h = length / f64::from(segments);
            let m = hermite_mass(h, density).unwrap();
            let mut actual = 0.0;
            for segment in 0..segments {
                let x = f64::from(segment) * h;
                let shifted = [
                    p[0] + x * (p[1] + x * (p[2] + x * p[3])),
                    p[1] + x * (2.0 * p[2] + 3.0 * x * p[3]),
                    p[2] + 3.0 * x * p[3],
                    p[3],
                ];
                let q = endpoints(shifted, h);
                actual += bilinear(&m, q, q);
            }
            assert!((actual - expected).abs() < 2e-13 * expected);
        }
    }

    #[test]
    fn consistent_hermite_mass_refuses_nonphysical_and_overflowing_inputs() {
        for (length, density) in [
            (0.0, 1.0),
            (-1.0, 1.0),
            (1.0, -1.0),
            (f64::NAN, 1.0),
            (1.0, f64::INFINITY),
            (1e200, 1e200),
        ] {
            assert!(hermite_mass(length, density).is_none());
        }
        assert_eq!(hermite_mass(0.5, 0.0).unwrap(), [[0.0; 4]; 4]);
    }
}
