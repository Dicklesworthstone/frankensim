//! Geometry and finite-control-volume momentum balance; no hydraulic integrator.

use fs_qty::Area;

#[derive(Clone, Debug, PartialEq)]
pub enum PipeLoadError {
    Invalid(&'static str),
    Cancelled,
}
impl core::fmt::Display for PipeLoadError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pipe load: {self:?}")
    }
}
impl std::error::Error for PipeLoadError {}
type Result<T> = core::result::Result<T, PipeLoadError>;

/// A fixed slender pipe's actual centreline, uniform flow area and cached
/// geometric integrals. End sections are normal to the first/last path segment.
/// Bulk velocity is axial and uniform over each section, with no swirl or spin.
/// Bends retain their full path for angular momentum and weight moments.
#[derive(Clone, Debug)]
pub struct PipeLoadGeometry {
    pub(super) vertices_m: Vec<[f64; 3]>,
    pub(super) area_m2: f64,
    pub(super) length_m: f64,
    sections: Vec<(f64, [f64; 3])>,
    tangents: [[f64; 3]; 2],
    chord: [f64; 3],
    line_moment: [f64; 3],
}

/// Interval-average fluid load on the enclosing solid, including its exterior
/// uniform pressure. Force and moment use the geometry's inertial frame. The
/// moment is about the retained `about_m`, not an implicit centre of mass.
/// Divide neither impulse by mass nor torque by inertia before passing them to
/// the mechanical owner: that owner retains the actual tensor and kinematics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PipeLoad {
    pub about_m: [f64; 3],
    pub duration_s: f64,
    pub force_n: [f64; 3],
    pub torque_nm: [f64; 3],
    pub impulse_ns: [f64; 3],
    pub angular_impulse_nms: [f64; 3],
}

pub(super) struct FlowInterval {
    pub duration_s: f64,
    pub volume_flows_m3_s: [f64; 2],
    pub density_kg_m3: f64,
    pub pressures_pa: [f64; 2],
    pub exterior_pressure_pa: f64,
    // Internal assembly cuts carry equal/opposite interface tractions; omit
    // them analytically from its total, retaining all fluid storage and weight.
    pub exposed_ports: [bool; 2],
    pub gravity_m_s2: [f64; 3],
}

impl PipeLoadGeometry {
    /// Compile the physical path once, with an explicit positive section cap.
    /// The original native branch independently admits its length and area.
    pub fn try_new(
        vertices_m: Vec<[f64; 3]>,
        flow_area: Area,
        maximum_sections: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        check(&mut cancelled)?;
        let area = flow_area.value();
        if vertices_m.len() < 2
            || vertices_m.len() - 1 > maximum_sections
            || !positive(area)
            || vertices_m.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(PipeLoadError::Invalid(
                "bounded finite pipe path and positive area required",
            ));
        }
        let first = vertices_m[0];
        let chord = vector(sub(*vertices_m.last().unwrap(), first))?;
        let mut sections = Vec::with_capacity(vertices_m.len() - 1);
        let (mut length, mut correction) = (0.0, 0.0);
        let mut moment = [0.0; 3];
        let mut tangents = [[0.0; 3]; 2];
        for (i, pair) in vertices_m.windows(2).enumerate() {
            check(&mut cancelled)?;
            let delta = vector(sub(pair[1], pair[0]))?;
            let scale = delta.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
            if !positive(scale) {
                return Err(PipeLoadError::Invalid(
                    "positive pipe section length required",
                ));
            }
            let n = delta.map(|x| x / scale);
            let ds = number(scale * fs_math::det::sqrt(dot(n, n)))?;
            if !positive(ds) {
                return Err(PipeLoadError::Invalid(
                    "representable pipe section required",
                ));
            }
            let tangent = vector(delta.map(|x| x / ds))?;
            if i == 0 {
                tangents[0] = tangent;
            }
            tangents[1] = tangent;
            let start = vector(sub(pair[0], first))?;
            moment = vector(add(moment, cross(start, delta)))?;
            sections.push((ds, vector(add(start, mul(delta, 0.5)))?));
            let (sum, error) = fs_math::eft::two_sum(length, ds);
            length = sum;
            correction += error;
        }
        let length = number(length + correction)?;
        if !positive(area * length) {
            return Err(PipeLoadError::Invalid(
                "representable positive pipe volume required",
            ));
        }
        check(&mut cancelled)?;
        Ok(Self {
            vertices_m,
            area_m2: area,
            length_m: length,
            sections,
            tangents,
            chord,
            line_moment: moment,
        })
    }

    // Normal points OUT of the point junction and INTO this adjacent pipe.
    pub(super) fn junction_port(&self, at_start: bool) -> ([f64; 3], [f64; 3]) {
        if at_start {
            (self.vertices_m[0], self.tangents[0])
        } else {
            (
                *self.vertices_m.last().unwrap(),
                mul(self.tangents[1], -1.0),
            )
        }
    }

    // The hydraulic trial supplies its actual interval, datum and refreshed
    // body-force density. No source query, flow root or second step occurs here.
    pub(super) fn evaluate(
        &self,
        flow: FlowInterval,
        about_m: [f64; 3],
        mut gravity_density: impl FnMut(usize) -> f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PipeLoad> {
        check(&mut cancelled)?;
        let FlowInterval {
            duration_s: dt,
            volume_flows_m3_s: [q0, q1],
            density_kg_m3: rho,
            pressures_pa: [p0, p1],
            exterior_pressure_pa: exterior,
            exposed_ports: exposed,
            gravity_m_s2: g,
        } = flow;
        if !positive(dt)
            || !positive(rho)
            || [q0, q1, p0, p1, exterior].iter().any(|x| !x.is_finite())
            || g.iter().chain(&about_m).any(|x| !x.is_finite())
        {
            return Err(PipeLoadError::Invalid(
                "finite original flow/pressure and positive time/density required",
            ));
        }
        // The flux is EVEN in Q: reversing both velocity and mass flux preserves
        // the outward momentum flux. Squaring average Q loses reversal loading.
        let q2 = affine_mean_square([q0, q1])?;
        let flux = number(dt * rho * (q2 / self.area_m2))?;
        let stored = number(rho * (q1 - q0))?;
        let pressure_in = if exposed[0] {
            number(dt * number((p0 - exterior) * self.area_m2)?)?
        } else {
            0.0
        };
        let pressure_out = if exposed[1] {
            number(dt * number((p1 - exterior) * self.area_m2)?)?
        } else {
            0.0
        };
        let flux_direction = std::array::from_fn(|axis| {
            (if exposed[1] {
                self.tangents[1][axis]
            } else {
                0.0
            }) - (if exposed[0] {
                self.tangents[0][axis]
            } else {
                0.0
            })
        });
        let outlet_moment = cross(self.chord, self.tangents[1]);
        let mut impulse = sub(
            sub(
                mul(self.tangents[0], pressure_in),
                mul(self.tangents[1], pressure_out),
            ),
            add(mul(flux_direction, flux), mul(self.chord, stored)),
        );
        let mut angular = sub(
            mul(
                outlet_moment,
                -pressure_out - if exposed[1] { flux } else { 0.0 },
            ),
            mul(self.line_moment, stored),
        );
        // Local source densities supply weight, while the native reference rho
        // supplies stored/convected momentum. This preserves Boussinesq scope.
        for (i, &(length, arm)) in self.sections.iter().enumerate() {
            check(&mut cancelled)?;
            let density = gravity_density(i);
            if !positive(density) {
                return Err(PipeLoadError::Invalid(
                    "positive original gravity density required",
                ));
            }
            let weight = vector(mul(g, number(dt * density * (self.area_m2 * length))?))?;
            impulse = vector(add(impulse, weight))?;
            angular = vector(add(angular, cross(arm, weight)))?;
        }
        // Translate a wrench, not its components' individual reference frames.
        angular = vector(add(
            angular,
            cross(vector(sub(self.vertices_m[0], about_m))?, impulse),
        ))?;
        let result = PipeLoad {
            about_m,
            duration_s: dt,
            force_n: vector(impulse.map(|x| x / dt))?,
            torque_nm: vector(angular.map(|x| x / dt))?,
            impulse_ns: vector(impulse)?,
            angular_impulse_nms: angular,
        };
        check(&mut cancelled)?;
        Ok(result)
    }
}

pub(super) fn affine_mean_square([a, b]: [f64; 2]) -> Result<f64> {
    let mean = 0.5 * a + 0.5 * b;
    let half_change = 0.5 * b - 0.5 * a;
    number(mean * mean + half_change * (half_change / 3.0))
}

pub(super) fn affine_mean_cube([a, b]: [f64; 2]) -> Result<f64> {
    let mean = 0.5 * a + 0.5 * b;
    let half_change = 0.5 * b - 0.5 * a;
    number(mean * number(mean * mean + half_change * half_change)?)
}

pub(super) fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
pub(super) fn number(x: f64) -> Result<f64> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(PipeLoadError::Invalid(
            "unrepresentable pipe momentum arithmetic",
        ))
    }
}
pub(super) fn vector(v: [f64; 3]) -> Result<[f64; 3]> {
    for x in v {
        number(x)?;
    }
    Ok(v)
}
pub(super) fn check(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(PipeLoadError::Cancelled)
    } else {
        Ok(())
    }
}
pub(super) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
pub(super) fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] + b[i])
}
pub(super) fn mul(a: [f64; 3], s: f64) -> [f64; 3] {
    a.map(|x| x * s)
}
pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
pub(super) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[cfg(test)]
mod tests;
