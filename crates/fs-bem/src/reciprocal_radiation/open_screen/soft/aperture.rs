//! Finite apertures in an infinite rigid baffle by exact acoustic Babinet.
//!
//! For the original soft-screen solution u_D=u_i+S*sigma with trace
//! S*sigma=-u_i on the aperture, the complementary rigid-baffle field is
//! -S*sigma on the transmitted side and u_i+u_i(reflected point)+S*sigma
//! on the incident side. Differentiating the same identity supplies the
//! physical Euler velocity; reflection acts on its vector, not only pressure.
//! See Gibbs et al. (2020), Remark 1, equation (4),
//! <https://doi.org/10.1007/s42985-020-00013-3>, and the 3D screen/aperture
//! equivalence in Chandler-Wilde et al. (2021), §1,
//! <https://doi.org/10.1007/s00211-021-01182-y>.
//!
//! This is the full-wave complementary boundary-value problem, not a
//! Kirchhoff prescription of unobstructed velocity in each hole. Disconnected
//! apertures retain the original single-layer solve's mutual interactions.

use super::*;
use crate::helmholtz::bc::{IncidentSample, IncidentSource, PlaneWave};
use crate::near_field;
use fs_ivl::{Sign, orient2d, orient3d};
use fs_math::dd::Dd;

const SCOPE: &str = "Finite coplanar apertures in an infinite stationary zero-thickness \
rigid acoustic baffle, with all prescribed sources in one declared incident half-space. \
Exact Babinet conversion of the original convergent pressure-release screen solve; \
no Kirchhoff opening approximation, finite wall thickness, finite baffle, material \
impedance, moving medium, source feedback or compact multipole-array qualification.";

/// Source-bound aperture response borrowing the original soft-screen density.
///
/// The surface mesh is the OPEN AREA, not the rigid complement. Its exact
/// plane separates the incident side from the transmitted side. Neither the
/// original LU, density nor Green rows are copied or re-solved by this adapter.
pub struct RigidApertureScattering<'b, 'p, 's, 'i> {
    source: &'b SoftScreenScattering<'p, 's, 'i>,
    plane: [Point; 3],
    transmission_normal: Point,
    transmitted_orientation: Sign,
}

/// Original adaptive rows plus the explicitly reflected source-side field.
/// Evaluation performs no Green quadrature or analytic source sampling.
pub struct PreparedRigidApertureReceivers<'b, 'p, 's, 'i> {
    rows: PreparedSoftScreenReceivers<'b, 'p, 's, 'i>,
    // None on the transmitted side; Some retains the already vector-reflected
    // incident field, including finite-source reactive terms and source phase.
    reflected: Vec<Option<IncidentSample>>,
}

impl<'p, 's, 'i> SoftScreenScattering<'p, 's, 'i> {
    /// Interpret this solved soft-screen support as holes in a rigid baffle.
    ///
    /// `transmission_normal` is unit and points FROM the source half-space TO
    /// the transmitted half-space. It must agree with one of the exact mesh
    /// plane's normals; the stored geometric normal is used for reflection.
    /// Every vertex must be EXACTLY coplanar under fs-ivl's orientation
    /// predicate. A rounded arbitrary rotation may therefore refuse; author
    /// geometry exactly on a plane rather than silently flattening a mesh.
    /// Finite sources must lie strictly on the incident side, and all plane
    /// waves must have a resolved positive component along the given normal.
    pub fn as_rigid_aperture(
        &self,
        transmission_normal: Point,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RigidApertureScattering<'_, 'p, 's, 'i>, HelmholtzError> {
        poll(&mut cancelled)?;
        let first = &self.prepared.facets[0];
        let plane = first.triangle;
        admit_plane(plane)?;
        let normal = first.normal;
        let alignment = dot(normal, transmission_normal);
        if transmission_normal.iter().any(|v| !v.is_finite())
            || (norm(transmission_normal) - 1.0).abs() > 1e-10
            || norm(cross(normal, transmission_normal)) > 64.0 * f64::EPSILON
        {
            return Err(bad(
                "rigid aperture needs a unit normal parallel to its exact plane",
            ));
        }
        for (index, &point) in self.prepared.surface.vertices.iter().enumerate() {
            if index % 128 == 0 {
                poll(&mut cancelled)?;
            }
            if plane_side(plane, point)? != Sign::Zero {
                return Err(bad("rigid-baffle apertures must be exactly coplanar"));
            }
        }
        let forward = alignment > 0.0;
        let out = RigidApertureScattering {
            source: self,
            plane,
            transmission_normal: normal.map(|v| if forward { v } else { -v }),
            transmitted_orientation: if forward {
                Sign::Negative
            } else {
                Sign::Positive
            },
        };
        let admit_plane = |plane: PlaneWave| {
            if precise_dot(plane.direction, out.transmission_normal).to_f64() <= 64.0 * f64::EPSILON
            {
                Err(bad(
                    "aperture plane waves must travel toward the declared baffle side",
                ))
            } else {
                Ok(())
            }
        };
        let admit_source = |position| {
            if out.side(position)? {
                Err(bad(
                    "all aperture sources must lie in the declared incident half-space",
                ))
            } else {
                Ok(())
            }
        };
        match self.incident {
            IncidentField::PlaneWave(plane) => admit_plane(plane)?,
            IncidentField::Monopoles(sources) => {
                for source in sources {
                    poll(&mut cancelled)?;
                    admit_source(source.position_m)?;
                }
            }
            IncidentField::ForceDipoles(sources) => {
                for source in sources {
                    poll(&mut cancelled)?;
                    admit_source(source.position_m)?;
                }
            }
            IncidentField::Sources(sources) => {
                for source in sources {
                    poll(&mut cancelled)?;
                    match source {
                        IncidentSource::PlaneWave(plane) => admit_plane(*plane)?,
                        IncidentSource::Monopole(s) => admit_source(s.position_m)?,
                        IncidentSource::ForceDipole(s) => admit_source(s.position_m)?,
                        IncidentSource::StressQuadrupole(s) => admit_source(s.position_m)?,
                    }
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(out)
    }
}

impl<'b, 'p, 's, 'i> RigidApertureScattering<'b, 'p, 's, 'i> {
    pub const fn scope(&self) -> &'static str {
        SCOPE
    }

    pub const fn transmission_normal(&self) -> Point {
        self.transmission_normal
    }

    /// Total outgoing transmitted power [W]. A planar single layer has equal
    /// outgoing powers in its two mirror hemispheres, so Babinet transmission
    /// is half the original scattered power. The same mesh/angular-rung
    /// convergence obligations apply; this does not certify a source load.
    pub fn transmitted_power_w(&self) -> f64 {
        0.5 * self.source.scattered_power_w
    }

    /// r*exp(-ikr)*p_transmitted(r*d) [Pa m] in the transmission hemisphere.
    /// The original global phase and complete mutually coupled density survive.
    pub fn transmitted_far_field(
        &self,
        direction: Point,
        cancelled: impl FnMut() -> bool,
    ) -> Result<C64, HelmholtzError> {
        if precise_dot(direction, self.transmission_normal).to_f64() <= 64.0 * f64::EPSILON {
            return Err(bad(
                "aperture far field needs a resolved transmission-side direction",
            ));
        }
        Ok(-self.source.scattered_far_field(direction, cancelled)?)
    }

    /// Bind receivers strictly off the baffle plane on either side. Exact
    /// on-plane points refuse instead of selecting an unspecified trace.
    /// Original aperture clearance and adaptive pressure/gradient gates apply.
    /// The byte ceiling includes retained original rows plus reflected source
    /// samples; temporary/output allocations retain the original row contract.
    /// Source work charges every original point and every reflected-side point
    /// before sampling, against the bound of the original scattering solution.
    pub fn prepare_receivers(
        &self,
        points: &[Point],
        separation: f64,
        options: near_field::Options,
        maximum_row_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedRigidApertureReceivers<'b, 'p, 's, 'i>, HelmholtzError> {
        poll(&mut cancelled)?;
        if points.is_empty() || points.len() > near_field::MAX_RECEIVERS {
            return Err(bad(
                "rigid-aperture receiver cardinality is outside its bound",
            ));
        }
        let extra = points
            .len()
            .checked_mul(core::mem::size_of::<Option<IncidentSample>>())
            .and_then(|n| n.checked_add(core::mem::size_of::<Vec<Option<IncidentSample>>>()));
        let core_bytes = extra
            .and_then(|extra| maximum_row_bytes.checked_sub(extra))
            .ok_or_else(|| bad("rigid-aperture reflected receiver storage exceeds its budget"))?;
        let mut incident_points = 0;
        for &point in points {
            poll(&mut cancelled)?;
            incident_points += usize::from(!self.side(point)?);
        }
        if self.source.incident.source_count() as u128 * (points.len() + incident_points) as u128
            > self.source.maximum_source_point_products as u128
        {
            return Err(bad(
                "rigid-aperture direct and reflected source work exceeds its budget",
            ));
        }
        let rho_omega = self.source.k * self.source.medium.sound_speed * self.source.medium.density;
        let mut reflected = Vec::with_capacity(points.len());
        for &point in points {
            poll(&mut cancelled)?;
            if self.side(point)? {
                reflected.push(None);
            } else {
                let reflected_point = self.reflect_point(point)?;
                let mut sample = self.source.incident.sample_wave(
                    self.source.k,
                    rho_omega,
                    Some(self.source.medium.density * self.source.medium.sound_speed),
                    reflected_point,
                    &mut cancelled,
                )?;
                let mut projection = ComplexSum::default();
                for (u, n) in sample
                    .velocity_m_s
                    .into_iter()
                    .zip(self.transmission_normal)
                {
                    projection.add_product(u, C64::from_re(n));
                }
                let projection = projection.finish();
                sample.velocity_m_s = core::array::from_fn(|axis| {
                    let mut sum = ComplexSum::default();
                    sum.add(sample.velocity_m_s[axis]);
                    sum.add_product(
                        projection,
                        C64::from_re(-2.0 * self.transmission_normal[axis]),
                    );
                    sum.finish()
                });
                reflected.push(Some(sample));
            }
        }
        let rows = self.source.prepare_receivers(
            points,
            separation,
            options,
            core_bytes,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(PreparedRigidApertureReceivers { rows, reflected })
    }

    fn side(&self, point: Point) -> Result<bool, HelmholtzError> {
        if point.iter().any(|v| !v.is_finite() || v.abs() > 1e9) {
            return Err(bad(
                "rigid aperture needs finite resolved point coordinates",
            ));
        }
        let sign = plane_side(self.plane, point)?;
        if sign == Sign::Zero {
            return Err(bad(
                "rigid-aperture source or receiver must lie strictly off its plane",
            ));
        }
        Ok(sign == self.transmitted_orientation)
    }

    fn reflect_point(&self, point: Point) -> Result<Point, HelmholtzError> {
        let mut distance = Dd::ZERO;
        for axis in 0..3 {
            distance = distance
                + Dd::from_f64(self.transmission_normal[axis])
                    * (Dd::from_f64(point[axis]) - Dd::from_f64(self.plane[0][axis]));
        }
        let reflected = core::array::from_fn(|axis| {
            (Dd::from_f64(point[axis])
                - Dd::from_f64(2.0 * self.transmission_normal[axis]) * distance)
                .to_f64()
        });
        if !self.side(reflected)? {
            return Err(bad(
                "rigid-aperture reflected point is unresolved across its plane",
            ));
        }
        Ok(reflected)
    }
}

impl PreparedRigidApertureReceivers<'_, '_, '_, '_> {
    pub fn points_m(&self) -> &[Point] {
        self.rows.points_m()
    }

    pub fn kernel_evaluations(&self) -> usize {
        self.rows.kernel_evaluations()
    }

    /// Total physical pressure and Euler velocity on the retained side of each
    /// receiver. Signed intensity and microphones consume these same p/u;
    /// quadrature discrepancies remain those of the original scattered rows.
    pub fn evaluate(
        &self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<near_field::VectorObservation, HelmholtzError> {
        let mut field = self.rows.evaluate_scattered(&mut cancelled)?;
        for index in 0..self.reflected.len() {
            poll(&mut cancelled)?;
            if let Some(reflected) = self.reflected[index] {
                let direct = self.rows.incident_samples()[index];
                let mut pressure = ComplexSum::default();
                pressure.add(field.scalar.pressure[index]);
                pressure.add(direct.pressure_pa);
                pressure.add(reflected.pressure_pa);
                field.scalar.pressure[index] = pressure.finish();
                for axis in 0..3 {
                    let mut velocity = ComplexSum::default();
                    velocity.add(field.particle_velocity_m_s[index][axis]);
                    velocity.add(direct.velocity_m_s[axis]);
                    velocity.add(reflected.velocity_m_s[axis]);
                    field.particle_velocity_m_s[index][axis] = velocity.finish();
                }
            } else {
                field.scalar.pressure[index] = -field.scalar.pressure[index];
                field.particle_velocity_m_s[index] = field.particle_velocity_m_s[index].map(|u| -u);
            }
            if !finite(field.scalar.pressure[index])
                || field.particle_velocity_m_s[index]
                    .iter()
                    .any(|&v| !finite(v))
            {
                return Err(bad("rigid-aperture total pressure or velocity overflows"));
            }
        }
        poll(&mut cancelled)?;
        Ok(field)
    }
}

fn precise_dot(a: Point, b: Point) -> Dd {
    a.into_iter().zip(b).fold(Dd::ZERO, |sum, (a, b)| {
        sum + Dd::from_f64(a) * Dd::from_f64(b)
    })
}

// Preserve the original predicate's no-intermediate-underflow premise.
// A common power-of-two division is exact. Admitted scaled coordinates lie
// in [-2,2], and each nonzero coordinate is at least 2^-128: every exact
// coordinate difference is then an integer multiple of 2^-180 and every
// degree-three determinant product of 2^-540, safely inside normal f64.
// Extreme coordinate dynamic range refuses; no tolerance declares a warp flat.
fn scaled_points(plane: [Point; 3], point: Point) -> Result<[Point; 4], HelmholtzError> {
    let points = [plane[0], plane[1], plane[2], point];
    let maximum = points.iter().flatten().fold(0.0_f64, |m, v| m.max(v.abs()));
    if !maximum.is_finite() || maximum == 0.0 {
        return Err(bad(
            "aperture plane predicate requires finite resolved coordinates",
        ));
    }
    let exponent = maximum.to_bits() & 0x7ff0_0000_0000_0000;
    let power = if exponent == 0 {
        f64::MIN_POSITIVE
    } else {
        f64::from_bits(exponent)
    };
    let mut scaled = [[0.0; 3]; 4];
    for i in 0..4 {
        for j in 0..3 {
            scaled[i][j] = points[i][j] / power;
            if !points[i][j].is_finite()
                || (points[i][j] != 0.0 && scaled[i][j].abs() < 2.938_735_877_055_719e-39)
            {
                return Err(bad(
                    "aperture plane predicate exceeds its exact exponent range",
                ));
            }
        }
    }
    Ok(scaled)
}

// Float cross products can produce a tiny nonzero vector even when three
// dyadic points are exactly collinear. Such a triple would make every 3D
// coplanarity test vacuously zero, so establish exact rank two first.
fn admit_plane(plane: [Point; 3]) -> Result<(), HelmholtzError> {
    let p = scaled_points(plane, plane[0])?;
    for [i, j] in [[0, 1], [1, 2], [2, 0]] {
        if orient2d([p[0][i], p[0][j]], [p[1][i], p[1][j]], [p[2][i], p[2][j]]) != Sign::Zero {
            return Ok(());
        }
    }
    Err(bad(
        "rigid aperture needs an exactly noncollinear defining triangle",
    ))
}

fn plane_side(plane: [Point; 3], point: Point) -> Result<Sign, HelmholtzError> {
    let p = scaled_points(plane, point)?;
    Ok(orient3d(p[0], p[1], p[2], p[3]))
}

#[cfg(test)]
mod tests;
