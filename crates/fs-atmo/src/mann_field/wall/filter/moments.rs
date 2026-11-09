//! Instantaneous spatial stress of the harmonic wall-blocked Mann realization.
//!
//! Each mode is the sum of a Fourier function and its decaying wall correction.
//! Their products integrate analytically over the three independent uniform
//! footprint coordinates. All bulk/bulk, bulk/wall and wall/wall cross terms
//! are retained; this is a realization's spatial covariance, not an ensemble
//! stress or a new turbulence closure. Small phases use central uniform moments.
//! Wall functions are normalized at the footprint's lowest point, so neither
//! preparation nor sampling forms an exponentially growing wall amplitude.

use super::{MannWallBoxFilter, error, poll, scaled_edge_gain};
use crate::mann_field::filter::moments::MannFilteredMoments;
use crate::mann_field::wall::{MannPlaneWall, ModePoint};
use crate::mann_field::{MannSpectralModel, det, dot};
use crate::{MAX_MODES, Refusal};

type Complex = [f64; 2];
type Covariance = [[f64; 2]; 2];

fn multiply(a: Complex, b: Complex) -> Complex {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

fn gain(z: Complex) -> Complex {
    let value = scaled_edge_gain(z[0], z[1]);
    [value[0], value[1]]
}

// Covariance of the real/imaginary parts of exp(z*s-|Re z|), s uniform
// on [-1,1]. Both functions are bounded on the interval. On small phases,
// integrate their centered Taylor series without subtracting rounded means.
fn edge_covariance(a: Complex, b: Complex) -> Result<Covariance, Refusal> {
    if !a.into_iter().chain(b).all(f64::is_finite) {
        return Err(error(
            "mann-wall-moments-invalid",
            "modal exponent is not finite",
        ));
    }
    if a == [0.0; 2] || b == [0.0; 2] {
        return Ok([[0.0; 2]; 2]);
    }
    if a.into_iter().chain(b).all(|x| x.abs() <= 0.5) {
        // |z| <= sqrt(1/2); the first omitted term is below 8e-18, and
        // its contribution to a central second moment is smaller still.
        const TERMS: usize = 16;
        let mut powers_a = [[0.0; 2]; TERMS];
        let mut powers_b = [[0.0; 2]; TERMS];
        let mut ta = [1.0, 0.0];
        let mut tb = [1.0, 0.0];
        for n in 0..TERMS {
            ta = multiply(ta, a).map(|x| x / (n + 1) as f64);
            tb = multiply(tb, b).map(|x| x / (n + 1) as f64);
            powers_a[n] = ta;
            powers_b[n] = tb;
        }
        let mut covariance = [[0.0; 2]; 2];
        for (n, pa) in powers_a.iter().enumerate() {
            // Odd/even cross moments vanish by symmetry. For both even,
            // Cov(s^n,s^m)=nm/[(n+m+1)(n+1)(m+1)]; for both odd it is
            // 1/(n+m+1). Here the polynomial degrees start at one.
            for m in (n % 2..TERMS).step_by(2) {
                let dn = (n + 1) as f64;
                let dm = (m + 1) as f64;
                let weight = if n % 2 == 0 {
                    1.0 / (dn + dm + 1.0)
                } else {
                    dn * dm / ((dn + dm + 1.0) * (dn + 1.0) * (dm + 1.0))
                };
                for i in 0..2 {
                    for j in 0..2 {
                        covariance[i][j] += weight * pa[i] * powers_b[m][j];
                    }
                }
            }
        }
        let scale = det::exp(-a[0].abs() - b[0].abs());
        return Ok(covariance.map(|row| row.map(|x| scale * x)));
    }
    let real = a[0] + b[0];
    let plus = a[1] + b[1];
    let minus = a[1] - b[1];
    let normalization = a[0].abs() + b[0].abs();
    if ![real, plus, minus, normalization]
        .into_iter()
        .all(f64::is_finite)
    {
        return Err(error(
            "mann-wall-moments-invalid",
            "modal product is not representable",
        ));
    }
    // In the actual wall basis both real parts have the same sign for
    // each edge. Keeping this scale also makes the identity general.
    let scale = det::exp(real.abs() - normalization);
    let plus = gain([real, plus]).map(|x| x * scale);
    let minus = gain([real, minus]).map(|x| x * scale);
    let ma = gain(a);
    let mb = gain(b);
    Ok([
        [
            0.5 * (minus[0] + plus[0]) - ma[0] * mb[0],
            0.5 * (plus[1] - minus[1]) - ma[0] * mb[1],
        ],
        [
            0.5 * (plus[1] + minus[1]) - ma[1] * mb[0],
            0.5 * (minus[0] - plus[0]) - ma[1] * mb[1],
        ],
    ])
}

// Terms of real/imaginary complex multiplication: (component, component, sign).
const PRODUCT: [[(usize, usize, f64); 2]; 2] =
    [[(0, 0, 1.0), (1, 1, -1.0)], [(0, 1, 1.0), (1, 0, 1.0)]];

pub(in crate::mann_field) fn pair_covariance(
    a: [Complex; 3],
    b: [Complex; 3],
) -> Result<Covariance, Refusal> {
    let mut mean_a = [1.0, 0.0];
    let mut mean_b = [1.0, 0.0];
    let mut covariance = [[0.0; 2]; 2];
    for axis in 0..3 {
        let edge = edge_covariance(a[axis], b[axis])?;
        let ea = gain(a[axis]);
        let eb = gain(b[axis]);
        let mut next = [[0.0; 2]; 2];
        // For independent axes, E[XY]=Cov(X,Y)+E[X]E[Y]. Expanding
        // the central product before evaluating avoids subtracting unit
        // means and preserves quartic cosine variance on tiny footprints.
        for i in 0..2 {
            for j in 0..2 {
                for (ip, ie, si) in PRODUCT[i] {
                    for (jp, je, sj) in PRODUCT[j] {
                        next[i][j] += si
                            * sj
                            * (covariance[ip][jp] * ea[ie] * eb[je]
                                + mean_a[ip] * mean_b[jp] * edge[ie][je]
                                + covariance[ip][jp] * edge[ie][je]);
                    }
                }
            }
        }
        covariance = next;
        mean_a = multiply(mean_a, ea);
        mean_b = multiply(mean_b, eb);
    }
    if !covariance.iter().flatten().all(|x| x.is_finite()) {
        return Err(error(
            "mann-wall-moments-invalid",
            "central modal product is not finite",
        ));
    }
    Ok(covariance)
}

/// Prepared quadratic observation of one wall-blocked spectrum and footprint.
/// Sampling includes every cross-mode product, draws no noise and allocates
/// nothing. A tangential footprint combines each mode's wall and bulk parts
/// before accumulation; this preserves exact normal impermeability at the wall.
#[derive(Clone, Debug, PartialEq)]
pub struct MannWallBoxMoments {
    filter: MannWallBoxFilter,
    parts: usize,
    coefficients: Vec<Covariance>,
}

impl MannWallBoxMoments {
    /// Prepare analytic products with explicit numerical storage and pair-work
    /// limits. Nonpoint tangential footprints use n functions; a footprint
    /// with vertical extent uses 2n. Each unordered function pair stores four
    /// f64 covariances. Preparation additionally uses at most 24 KiB of stack
    /// exponent scratch, independent of the byte budget for retained storage.
    pub fn prepare(
        model: &MannSpectralModel,
        wall: MannPlaneWall,
        half_edges_m: [[f64; 3]; 3],
        maximum_workspace_bytes: usize,
        maximum_function_pairs: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Refusal> {
        let filter = MannWallBoxFilter::prepare(
            model,
            wall,
            half_edges_m,
            maximum_workspace_bytes,
            &mut cancelled,
        )?;
        let point = half_edges_m.iter().flatten().all(|&x| x == 0.0);
        let parts = if point {
            0
        } else if filter.vertical_radius_m == 0.0 {
            1
        } else {
            2
        };
        let count = parts * model.mode_count();
        let pairs = count * (count + 1) / 2;
        let bytes = filter.workspace_bytes() + pairs * core::mem::size_of::<Covariance>();
        if pairs > maximum_function_pairs || bytes > maximum_workspace_bytes {
            return Err(error(
                "mann-wall-moments-budget",
                "quadratic wall observation exceeds its pair or byte budget",
            ));
        }
        let mut coefficients = Vec::new();
        coefficients.try_reserve_exact(pairs).map_err(|_| {
            error(
                "mann-wall-moments-allocation",
                "unable to reserve quadratic wall coefficients",
            )
        })?;
        if filter.workspace_bytes() + coefficients.capacity() * core::mem::size_of::<Covariance>()
            > maximum_workspace_bytes
        {
            return Err(error(
                "mann-wall-moments-budget",
                "allocated quadratic wall capacities exceed the byte budget",
            ));
        }
        let mut exponents = [[[0.0; 2]; 3]; 2 * MAX_MODES];
        if parts != 0 {
            for (m, mode) in model.modes.iter().enumerate() {
                poll(&mut cancelled)?;
                let q = det::hypot(mode.k[0], mode.k[1]);
                for (axis, edge) in half_edges_m.iter().enumerate() {
                    exponents[parts * m][axis] = [0.0, dot(mode.k, *edge)];
                    if parts == 2 {
                        exponents[parts * m + 1][axis] =
                            [-q * edge[2], mode.k[0] * edge[0] + mode.k[1] * edge[1]];
                    }
                }
            }
        }
        for m in 0..count {
            for n in m..count {
                poll(&mut cancelled)?;
                coefficients.push(pair_covariance(exponents[m], exponents[n])?);
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            filter,
            parts,
            coefficients,
        })
    }

    /// Original linear observation; its mean/gradient arithmetic is unchanged.
    pub const fn mean_filter(&self) -> &MannWallBoxFilter {
        &self.filter
    }

    /// Allocated numerical capacities, including the original linear filter.
    pub fn workspace_bytes(&self) -> usize {
        self.filter.workspace_bytes()
            + self.coefficients.capacity() * core::mem::size_of::<Covariance>()
    }

    /// Unordered bulk/wall function pairs; zero for an exact point.
    pub fn function_pairs(&self) -> usize {
        self.coefficients.len()
    }
}

pub(in crate::mann_field) fn sample(
    model: &MannSpectralModel,
    time: f64,
    amplitudes: [&[f64]; 4],
    point: [f64; 3],
    observation: &MannWallBoxMoments,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MannFilteredMoments, Refusal> {
    poll(&mut cancelled)?;
    let height = observation.filter.admit(model, point)?;
    let mean = super::sample_modes(model, time, amplitudes, point, &observation.filter)?;
    let mut stress = [[0.0; 3]; 3];
    if observation.parts != 0 {
        // At most 24 KiB: two vector quadratures per bulk/wall function.
        let mut vectors = [[[0.0; 3]; 2]; 2 * MAX_MODES];
        let position = core::array::from_fn(|i| point[i] - model.advection_mps[i] * time);
        for (m, mode) in model.modes.iter().enumerate() {
            poll(&mut cancelled)?;
            let a = core::array::from_fn(|i| {
                mode.polarization[0][i] * amplitudes[0][m]
                    + mode.polarization[1][i] * amplitudes[2][m]
            });
            let b = core::array::from_fn(|i| {
                mode.polarization[0][i] * amplitudes[1][m]
                    + mode.polarization[1][i] * amplitudes[3][m]
            });
            if observation.parts == 1 {
                let geometry = ModePoint::new(
                    mode,
                    model.advection_mps,
                    time,
                    point,
                    observation.filter.wall,
                    height,
                )?;
                vectors[m] = [geometry.values(a, b), geometry.values(b, a.map(|x| -x))];
                continue;
            }
            let phase = dot(mode.k, position);
            let (s, c) = det::sin_cos(phase);
            vectors[2 * m] = [
                core::array::from_fn(|i| a[i] * c + b[i] * s),
                core::array::from_fn(|i| -a[i] * s + b[i] * c),
            ];
            let q = det::hypot(mode.k[0], mode.k[1]);
            let decay = det::exp(-q * (height - observation.filter.vertical_radius_m));
            if q == 0.0 || decay == 0.0 {
                continue;
            }
            let wall_phase = dot(
                mode.k,
                [
                    position[0],
                    position[1],
                    observation.filter.wall.height_m() - model.advection_mps[2] * time,
                ],
            );
            if !wall_phase.is_finite() {
                return Err(error(
                    "mann-wall-moments-invalid",
                    "wall product phase is not representable",
                ));
            }
            let (s, c) = det::sin_cos(wall_phase);
            let ac = a[2] * c + b[2] * s;
            let bc = -a[2] * s + b[2] * c;
            vectors[2 * m + 1] = [
                [mode.k[0] / q * bc, mode.k[1] / q * bc, -ac].map(|x| x * decay),
                [-mode.k[0] / q * ac, -mode.k[1] / q * ac, -bc].map(|x| x * decay),
            ];
        }
        let count = observation.parts * model.mode_count();
        let mut coefficients = observation.coefficients.iter();
        for m in 0..count {
            for n in m..count {
                poll(&mut cancelled)?;
                let pair = coefficients
                    .next()
                    .expect("prepared wall function-pair count");
                for r in 0..2 {
                    for s in 0..2 {
                        if pair[r][s] == 0.0 {
                            continue;
                        }
                        for i in 0..3 {
                            for j in i..3 {
                                let mut product = vectors[m][r][i] * vectors[n][s][j];
                                if m != n {
                                    product += vectors[n][s][i] * vectors[m][r][j];
                                }
                                stress[i][j] += pair[r][s] * product;
                            }
                        }
                    }
                }
            }
        }
        for i in 0..3 {
            for j in 0..i {
                stress[i][j] = stress[j][i];
            }
        }
    }
    if !stress.iter().flatten().all(|x| x.is_finite()) {
        return Err(error(
            "mann-wall-moments-invalid",
            "wall subfilter stress is not finite",
        ));
    }
    poll(&mut cancelled)?;
    Ok(MannFilteredMoments {
        mean,
        subfilter_stress_m2_s2: stress,
    })
}

#[cfg(test)]
mod tests;
