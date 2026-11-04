//! Uniform spatial averages of the original harmonic wall observation.
//!
//! A half edge `a` gives the homogeneous gain `sinc(k.a)` and the wall gain
//! `sinh(z)/z`, where `z = -q*a_z + i*(kx*a_x + ky*a_y)`. Products of these
//! exact one-dimensional integrals describe any translated parallelepiped,
//! including collapsed face/line/point footprints. Analytic differentiation
//! commutes with this fixed observation operator, preserving divergence.
//!
//! The wall gain is stored multiplied by `exp(-q*r)`, where
//! `r = sum(abs(a_z))`. At sampling it receives `exp(-q*(h-r))`. Both factors
//! are bounded for a footprint wholly above the wall. Thus a tall box touching
//! the wall retains its wall-layer average even when its center correction
//! underflows. No variance is restored and no nonlinear flux is averaged.

use super::{MannPlaneWall, ModePoint};
use crate::mann_field::filter::MannBoxFilter;
use crate::mann_field::{
    FlowSample, MannSpectralModel, Refusal, SpectralMode, det, dot, finite_covariance,
    finite_vector, refuse,
};

/// Source-bound gains for one fixed footprint and slip plane. Sampling and
/// covariance read the existing realization without allocations or new draws.
#[derive(Clone, Debug, PartialEq)]
pub struct MannWallBoxFilter {
    homogeneous: MannBoxFilter,
    wall: MannPlaneWall,
    vertical_radius_m: f64,
    // Scaled wall gain (real, imaginary), and bulk minus scaled real gain.
    // The last value retains tiny near-wall differences during preparation.
    correction: Vec<[f64; 3]>,
}

fn error(code: &'static str, message: &str) -> Refusal {
    refuse(
        code,
        message.into(),
        "use a matching Mann model, finite footprint wholly above its plane and explicit gain-storage budget",
    )
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(error(
            "mann-wall-filter-cancelled",
            "filter preparation cancelled",
        ))
    } else {
        Ok(())
    }
}

fn radius(half_edges: [[f64; 3]; 3]) -> Result<f64, Refusal> {
    let r = half_edges.iter().map(|e| e[2].abs()).sum::<f64>();
    if !half_edges.iter().flatten().all(|x| x.is_finite()) || !r.is_finite() {
        return Err(error("mann-wall-filter-invalid", "footprint is not finite"));
    }
    Ok(r)
}

fn sinc_minus_one(x: f64) -> f64 {
    if x.abs() < 1e-3 {
        let x2 = x * x;
        x2 * (-1.0 / 6.0 + x2 * (1.0 / 120.0 - x2 / 5040.0))
    } else {
        det::sin(x) / x - 1.0
    }
}

// Return (real, imaginary, real minus one) of exp(-|Re z|) sinh(z)/z.
// Retain both the gain and increment: deep-box gains can be tiny, while
// tiny-footprint increments must not be lost by subtracting two rounded ones.
fn scaled_edge_gain(real: f64, imaginary: f64) -> [f64; 3] {
    let alpha = real.abs();
    if alpha == 0.0 {
        let delta = sinc_minus_one(imaginary);
        let gain = if imaginary.abs() < 1e-3 {
            1.0 + delta
        } else {
            det::sin(imaginary) / imaginary
        };
        return [gain, 0.0, delta];
    }
    // sinh(z)/z is even. Choose Re z <= 0, with no large positive exponent.
    let beta = if real > 0.0 { -imaginary } else { imaginary };
    let scale = alpha.max(beta.abs());
    if scale < 1e-3 {
        let z2 = [alpha * alpha - beta * beta, -2.0 * alpha * beta];
        let mul = |a: [f64; 2], b: [f64; 2]| [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
        let z4 = mul(z2, z2);
        let z6 = mul(z4, z2);
        let increment =
            core::array::from_fn::<_, 2, _>(|i| z2[i] / 6.0 + z4[i] / 120.0 + z6[i] / 5040.0);
        let decay = det::exp(-alpha);
        let delta = det::expm1(-alpha) + decay * increment[0];
        return [1.0 + delta, decay * increment[1], delta];
    }
    let decay = det::exp(-2.0 * alpha);
    let numerator = [
        0.5 * det::expm1(-2.0 * alpha) * det::cos(beta),
        0.5 * (1.0 + decay) * det::sin(beta),
    ];
    let re = -alpha / scale;
    let im = beta / scale;
    let denominator = re * re + im * im;
    let gain = ((numerator[0] * re + numerator[1] * im) / denominator) / scale;
    [
        gain,
        ((numerator[1] * re - numerator[0] * im) / denominator) / scale,
        gain - 1.0,
    ]
}

impl MannWallBoxFilter {
    /// Numerical gain bytes: zero for a point, 8 per mode for a tangential
    /// footprint, 32 per mode with vertical extent. Fixed metadata is excluded.
    pub fn required_workspace_bytes(
        model: &MannSpectralModel,
        half_edges_m: [[f64; 3]; 3],
    ) -> Result<usize, Refusal> {
        let r = radius(half_edges_m)?;
        let values = if half_edges_m.iter().flatten().all(|&x| x == 0.0) {
            0
        } else if r == 0.0 {
            1
        } else {
            4
        };
        Ok(model.mode_count() * values * core::mem::size_of::<f64>())
    }

    /// Prepare analytic gains in O(modes). The entire translated footprint
    /// must fit above the wall at each later observation; no clipping occurs.
    pub fn prepare(
        model: &MannSpectralModel,
        wall: MannPlaneWall,
        half_edges_m: [[f64; 3]; 3],
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Refusal> {
        poll(&mut cancelled)?;
        let vertical_radius_m = radius(half_edges_m)?;
        if Self::required_workspace_bytes(model, half_edges_m)? > maximum_workspace_bytes {
            return Err(error(
                "mann-wall-filter-budget",
                "filter exceeds its gain-storage budget",
            ));
        }
        let homogeneous =
            MannBoxFilter::prepare(model, half_edges_m, maximum_workspace_bytes, &mut cancelled)?;
        let mut correction = Vec::new();
        if vertical_radius_m != 0.0 {
            correction
                .try_reserve_exact(model.mode_count())
                .map_err(|_| {
                    error(
                        "mann-wall-filter-allocation",
                        "unable to reserve wall gains",
                    )
                })?;
            for (index, mode) in model.modes.iter().enumerate() {
                poll(&mut cancelled)?;
                let q = det::hypot(mode.k[0], mode.k[1]);
                let mut bulk_delta = 0.0;
                let mut wall_gain = [1.0, 0.0];
                let mut wall_delta = 0.0;
                for edge in half_edges_m {
                    let bulk_angle = dot(mode.k, edge);
                    let real = -q * edge[2];
                    let imaginary = mode.k[0] * edge[0] + mode.k[1] * edge[1];
                    if !bulk_angle.is_finite() || !real.is_finite() || !imaginary.is_finite() {
                        return Err(error(
                            "mann-wall-filter-invalid",
                            "filter phase is not representable",
                        ));
                    }
                    let bulk = sinc_minus_one(bulk_angle);
                    bulk_delta += bulk + bulk_delta * bulk;
                    let w = scaled_edge_gain(real, imaginary);
                    wall_delta += w[2] + wall_delta * w[2] - wall_gain[1] * w[1];
                    wall_gain = [
                        wall_gain[0] * w[0] - wall_gain[1] * w[1],
                        wall_gain[1] * w[0] + wall_gain[0] * w[1],
                    ];
                }
                let bulk = homogeneous
                    .gains(model)?
                    .expect("vertical extent is not a point")[index];
                let gap = if bulk.abs().max(wall_gain[0].abs()) < 0.5 {
                    bulk - wall_gain[0]
                } else {
                    bulk_delta - wall_delta
                };
                let gain = [wall_gain[0], wall_gain[1], gap];
                if !gain.iter().all(|x| x.is_finite()) {
                    return Err(error(
                        "mann-wall-filter-invalid",
                        "wall gain is not representable",
                    ));
                }
                correction.push(gain);
            }
        }
        let filter = Self {
            homogeneous,
            wall,
            vertical_radius_m,
            correction,
        };
        if filter.workspace_bytes() > maximum_workspace_bytes {
            return Err(error(
                "mann-wall-filter-budget",
                "allocated gain capacities exceed the budget",
            ));
        }
        poll(&mut cancelled)?;
        Ok(filter)
    }

    pub const fn wall(&self) -> MannPlaneWall {
        self.wall
    }
    pub fn half_edges_m(&self) -> [[f64; 3]; 3] {
        self.homogeneous.half_edges_m()
    }
    pub const fn vertical_radius_m(&self) -> f64 {
        self.vertical_radius_m
    }
    pub fn workspace_bytes(&self) -> usize {
        self.homogeneous.workspace_bytes()
            + self.correction.capacity() * core::mem::size_of::<[f64; 3]>()
    }

    fn admit(&self, model: &MannSpectralModel, point: [f64; 3]) -> Result<f64, Refusal> {
        self.homogeneous.gains(model)?;
        let height = self.wall.distance(point)?;
        if height < self.vertical_radius_m {
            return Err(error(
                "mann-wall-filter-domain",
                "footprint crosses below the slip plane",
            ));
        }
        Ok(height)
    }

    fn geometry(
        &self,
        model: &MannSpectralModel,
        index: usize,
        time: f64,
        point: [f64; 3],
        wall: MannPlaneWall,
        height: f64,
    ) -> Result<ModePoint, Refusal> {
        let mode = &model.modes[index];
        let bulk = self.homogeneous.gains(model)?.map_or(1.0, |g| g[index]);
        if self.vertical_radius_m == 0.0 {
            let mut point = ModePoint::new(mode, model.advection_mps, time, point, wall, height)?;
            point.phase = point.phase.map(|x| x * bulk);
            point.wall_phase = point.wall_phase.map(|x| x * bulk);
            point.difference = point.difference.map(|x| x * bulk);
            return Ok(point);
        }
        averaged_geometry(
            mode,
            model.advection_mps,
            time,
            point,
            wall,
            height,
            self.vertical_radius_m,
            bulk,
            self.correction[index],
        )
    }
}

fn averaged_geometry(
    mode: &SpectralMode,
    advection: [f64; 3],
    time: f64,
    point: [f64; 3],
    wall: MannPlaneWall,
    height: f64,
    radius: f64,
    bulk: f64,
    gain: [f64; 3],
) -> Result<ModePoint, Refusal> {
    let k = mode.k;
    let position = core::array::from_fn(|i| point[i] - advection[i] * time);
    let angle = dot(k, position);
    let q = det::hypot(k[0], k[1]);
    if !angle.is_finite() || !q.is_finite() {
        return Err(error(
            "mann-wall-filter-invalid",
            "observation phase is not representable",
        ));
    }
    let attenuation = q * (height - radius);
    let decay = det::exp(-attenuation);
    let mut phase = [bulk * det::cos(angle), bulk * det::sin(angle)];
    if q == 0.0 || decay == 0.0 {
        if q == 0.0 && mode.polarization.iter().any(|p| p[2] != 0.0) {
            return Err(error(
                "mann-wall-filter-invalid",
                "vertical mode has a normal amplitude",
            ));
        }
        return Ok(ModePoint {
            k,
            q,
            horizontal: if q == 0.0 {
                [0.0; 2]
            } else {
                [k[0] / q, k[1] / q]
            },
            phase,
            wall_phase: [0.0; 2],
            difference: phase,
        });
    }
    let wall_angle = dot(
        k,
        [
            position[0],
            position[1],
            wall.height_m - advection[2] * time,
        ],
    );
    if !wall_angle.is_finite() {
        return Err(error(
            "mann-wall-filter-invalid",
            "wall phase is not representable",
        ));
    }
    let wall_trig = [det::cos(wall_angle), det::sin(wall_angle)];
    let rotate = |v: [f64; 2]| {
        [
            wall_trig[0] * v[0] - wall_trig[1] * v[1],
            wall_trig[1] * v[0] + wall_trig[0] * v[1],
        ]
    };
    let wall_phase = rotate([decay * gain[0], decay * gain[1]]);
    let delta = k[2] * height;
    let difference = if delta.abs() <= 0.5 && attenuation <= 0.5 {
        let sine = det::sin(delta);
        let half_sine = det::sin(0.5 * delta);
        let cos_minus_one = -2.0 * half_sine * half_sine;
        phase = rotate([bulk * (1.0 + cos_minus_one), bulk * sine]);
        rotate([
            bulk * cos_minus_one + gain[2] - gain[0] * det::expm1(-attenuation),
            bulk * sine - decay * gain[1],
        ])
    } else {
        core::array::from_fn(|i| phase[i] - wall_phase[i])
    };
    Ok(ModePoint {
        k,
        q,
        horizontal: [k[0] / q, k[1] / q],
        phase,
        wall_phase,
        difference,
    })
}

pub(in crate::mann_field) fn sample_modes(
    model: &MannSpectralModel,
    time: f64,
    amplitudes: [&[f64]; 4],
    point: [f64; 3],
    filter: &MannWallBoxFilter,
) -> Result<FlowSample, Refusal> {
    let height = filter.admit(model, point)?;
    // Preserve the established point observation's arithmetic exactly.
    if filter.homogeneous.gains(model)?.is_none() {
        return super::sample_modes(model, time, amplitudes, point, filter.wall);
    }
    let mut sample = FlowSample {
        u: [0.0; 3],
        grad: [[0.0; 3]; 3],
    };
    for (index, mode) in model.modes.iter().enumerate() {
        let geometry = filter.geometry(model, index, time, point, filter.wall, height)?;
        let a = core::array::from_fn(|i| {
            mode.polarization[0][i] * amplitudes[0][index]
                + mode.polarization[1][i] * amplitudes[2][index]
        });
        let b = core::array::from_fn(|i| {
            mode.polarization[0][i] * amplitudes[1][index]
                + mode.polarization[1][i] * amplitudes[3][index]
        });
        geometry.add_sample(a, b, &mut sample);
    }
    if !sample
        .u
        .iter()
        .chain(sample.grad.iter().flatten())
        .all(|x| x.is_finite())
    {
        return Err(error(
            "mann-wall-filter-invalid",
            "averaged field is not representable",
        ));
    }
    Ok(sample)
}

pub(in crate::mann_field) fn covariance_modes(
    model: &MannSpectralModel,
    times: Option<&[f64]>,
    first: [f64; 3],
    second: [f64; 3],
    lag: f64,
    filter: &MannWallBoxFilter,
) -> Result<[[f64; 3]; 3], Refusal> {
    let first_height = filter.admit(model, first)?;
    let second_height = filter.admit(model, second)?;
    if !lag.is_finite() {
        return Err(error(
            "mann-wall-filter-invalid",
            "covariance lag must be finite",
        ));
    }
    if filter.homogeneous.gains(model)?.is_none() {
        return super::covariance_modes(model, times, first, second, lag, filter.wall);
    }
    let displacement = core::array::from_fn(|i| second[i] - first[i]);
    finite_vector(displacement)?;
    let wall = MannPlaneWall {
        height_m: -first_height,
    };
    let mut covariance = [[0.0; 3]; 3];
    for (index, mode) in model.modes.iter().enumerate() {
        let a = filter.geometry(model, index, 0.0, [0.0; 3], wall, first_height)?;
        let b = filter.geometry(model, index, lag, displacement, wall, second_height)?;
        let decay = times.map_or(1.0, |t| det::exp(-lag.abs() / t[index]));
        for &polarization in &mode.polarization {
            let ac = a.values(polarization, [0.0; 3]);
            let as_ = a.values([0.0; 3], polarization);
            let bc = b.values(polarization, [0.0; 3]);
            let bs = b.values([0.0; 3], polarization);
            for i in 0..3 {
                for j in 0..3 {
                    covariance[i][j] += decay * (ac[i] * bc[j] + as_[i] * bs[j]);
                }
            }
        }
    }
    finite_covariance(covariance)
}

#[cfg(test)]
mod tests;
