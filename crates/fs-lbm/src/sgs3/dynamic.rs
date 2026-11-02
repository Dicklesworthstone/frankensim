//! Germano--Lilly dynamic viscosity on a uniform, periodic Cartesian field.
//!
//! The test filter is the equal-volume average of the 27 adjacent cells.
//! Interpreting the resolved values as cell averages gives a grid box of width
//! `h` and a test box of width `3 h`; the modeled squared-width ratio is NINE,
//! not the often-assumed four. Density-weighted (Favre) test velocities retain
//! variable-density covariance. Only the deviatoric stress is modeled.
//!
//! With `A = rho |S| S_dev`, `L = hat(rho u u) - hat(rho u)hat(rho u)/hat(rho)`
//! and `M = 2 [hat(A) - 9 hat(rho) |S_hat| S_hat_dev]`, the least-squares
//! coefficient is `C_s^2 = <L_dev:M>/<M:M>`. Lengths in this identity are in
//! cells; the returned viscosity restores physical `h`. Averaging is either
//! over the entire periodic domain or over the same local 27-cell test box.
//! Negative coefficients are explicitly set to zero (no backscatter). Positive
//! coefficients above the caller's bound REFUSE rather than being truncated.
//!
//! References: Lilly (1992), doi:10.1063/1.858280; Moin, Squires, Cabot & Lee
//! (1991), doi:10.1063/1.858164. This implements the deviatoric dynamic viscosity,
//! not their complete compressible SGS energy/heat model. It supplies no wall
//! treatment, shock model, scale-dependent correction or DNS validation.

use super::SgsUse;
use fs_math::det;

/// Numerator and denominator are averaged separately, never their ratio.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DynamicAveraging {
    /// One coefficient for an explicitly homogeneous periodic domain.
    Domain,
    /// A spatially varying coefficient, stabilized over the same test box.
    LocalBox,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DynamicSgsConfig {
    pub averaging: DynamicAveraging,
    /// Upper admission bound on C_s squared, not on C_s.
    pub maximum_coefficient_squared: f64,
    /// Used by the consuming solver; never a viscosity clamp.
    pub maximum_eddy_to_molecular: f64,
    pub usage: SgsUse,
}

impl DynamicSgsConfig {
    pub fn validate(self) -> Result<(), DynamicSgsError> {
        if self.usage == SgsUse::AcceptanceRender {
            return Err(DynamicSgsError::Invalid("acceptance use is not supported"));
        }
        if !self.maximum_coefficient_squared.is_finite()
            || self.maximum_coefficient_squared < 0.0
            || !self.maximum_eddy_to_molecular.is_finite()
            || self.maximum_eddy_to_molecular < 0.0
        {
            return Err(DynamicSgsError::Invalid(
                "finite nonnegative bounds required",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum DynamicSgsError {
    Invalid(&'static str),
    WorkspaceBytes {
        required: usize,
        maximum: usize,
    },
    Allocation,
    Cell {
        cell: usize,
        reason: &'static str,
    },
    Coefficient {
        cell: usize,
        observed: f64,
        maximum: f64,
    },
    Cancelled,
}

impl core::fmt::Display for DynamicSgsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LBM-DYNAMIC-SGS: {self:?}")
    }
}
impl std::error::Error for DynamicSgsError {}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DynamicSgsReport {
    pub maximum_coefficient_squared: f64,
    pub minimum_raw_coefficient_squared: f64,
    /// Cells where a negative least-squares coefficient was set to zero.
    pub backscatter_limited_cells: usize,
    /// Cells with exactly zero least-squares denominator.
    pub undetermined_cells: usize,
    pub maximum_eddy_viscosity: f64,
}

/// Reusable nonlocal filter storage. Preparation allocates nothing. Scratch
/// and returned viscosities are consumable only after successful preparation.
/// Cell order is x-fastest, followed by y and z; every dimension is periodic.
pub struct PeriodicDynamicSgs {
    shape: [usize; 3],
    spacing: f64,
    velocity: Vec<[f64; 3]>,
    filtered: Vec<[f64; 4]>,
    basis: Vec<[f64; 6]>,
    strain_rate: Vec<f64>,
    contraction: Vec<[f64; 2]>,
    viscosity: Vec<f64>,
}

fn allocated<T: Clone>(n: usize, value: T) -> Result<Vec<T>, DynamicSgsError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| DynamicSgsError::Allocation)?;
    values.resize(n, value);
    Ok(values)
}

fn poll(cancelled: &mut impl FnMut() -> bool, cell: usize) -> Result<(), DynamicSgsError> {
    if cell % 64 == 0 && cancelled() {
        Err(DynamicSgsError::Cancelled)
    } else {
        Ok(())
    }
}

fn index(shape: [usize; 3], p: [usize; 3]) -> usize {
    (p[2] * shape[1] + p[1]) * shape[0] + p[0]
}

fn coordinates(shape: [usize; 3], i: usize) -> [usize; 3] {
    [
        i % shape[0],
        (i / shape[0]) % shape[1],
        i / (shape[0] * shape[1]),
    ]
}

fn neighbors(n: usize, i: usize) -> [usize; 3] {
    [
        if i == 0 { n - 1 } else { i - 1 },
        i,
        if i + 1 == n { 0 } else { i + 1 },
    ]
}

fn each_box(shape: [usize; 3], cell: usize, mut visit: impl FnMut(usize)) {
    let p = coordinates(shape, cell);
    for z in neighbors(shape[2], p[2]) {
        for y in neighbors(shape[1], p[1]) {
            for x in neighbors(shape[0], p[0]) {
                visit(index(shape, [x, y, z]));
            }
        }
    }
}

// Packed symmetric tensor: xx, yy, zz, xy, xz, yz. Off-diagonals count twice.
fn contract(a: [f64; 6], b: [f64; 6]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + 2.0 * (a[3] * b[3] + a[4] * b[4] + a[5] * b[5])
}

fn deviator(mut a: [f64; 6]) -> [f64; 6] {
    let trace = (a[0] + a[1] + a[2]) / 3.0;
    for v in &mut a[..3] {
        *v -= trace;
    }
    a
}

fn strain<const N: usize>(shape: [usize; 3], cell: usize, u: &[[f64; N]]) -> (f64, [f64; 6]) {
    let at = coordinates(shape, cell);
    let mut g = [[0.0; 3]; 3];
    for axis in 0..3 {
        let sides = neighbors(shape[axis], at[axis]);
        let mut low = at;
        let mut high = at;
        low[axis] = sides[0];
        high[axis] = sides[2];
        for component in 0..3 {
            g[component][axis] =
                0.5 * (u[index(shape, high)][component] - u[index(shape, low)][component]);
        }
    }
    let s = [
        g[0][0],
        g[1][1],
        g[2][2],
        0.5 * (g[0][1] + g[1][0]),
        0.5 * (g[0][2] + g[2][0]),
        0.5 * (g[1][2] + g[2][1]),
    ];
    (det::sqrt(2.0 * contract(s, s)), deviator(s))
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let y = value - self.correction;
        let t = self.value + y;
        self.correction = (t - self.value) - y;
        self.value = t;
    }
}

impl PeriodicDynamicSgs {
    /// Numerical payload, excluding the fixed-size owner and allocator metadata.
    pub fn required_bytes(cells: usize) -> Result<usize, DynamicSgsError> {
        cells
            .checked_mul(17 * core::mem::size_of::<f64>())
            .ok_or(DynamicSgsError::Invalid("workspace size overflow"))
    }

    pub fn new(
        shape: [usize; 3],
        spacing: f64,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        if shape.iter().any(|&n| n < 3) || !spacing.is_finite() || spacing <= 0.0 {
            return Err(DynamicSgsError::Invalid(
                "at least three cells per periodic axis and positive isotropic spacing required",
            ));
        }
        let n = shape
            .into_iter()
            .try_fold(1usize, |n, d| n.checked_mul(d))
            .ok_or(DynamicSgsError::Invalid("cell count overflow"))?;
        let required = Self::required_bytes(n)?;
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            });
        }
        let result = Self {
            shape,
            spacing,
            velocity: allocated(n, [0.0; 3])?,
            filtered: allocated(n, [0.0; 4])?,
            basis: allocated(n, [0.0; 6])?,
            strain_rate: allocated(n, 0.0)?,
            contraction: allocated(n, [0.0; 2])?,
            viscosity: allocated(n, 0.0)?,
        };
        let required = result.workspace_bytes();
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            });
        }
        Ok(result)
    }

    /// Actual numerical vector capacities, excluding the fixed-size owner and
    /// allocator metadata. Consumers include this in their complete byte gate.
    pub fn workspace_bytes(&self) -> usize {
        (3 * self.velocity.capacity()
            + 4 * self.filtered.capacity()
            + 6 * self.basis.capacity()
            + self.strain_rate.capacity()
            + 2 * self.contraction.capacity()
            + self.viscosity.capacity())
            * core::mem::size_of::<f64>()
    }

    pub const fn shape(&self) -> [usize; 3] {
        self.shape
    }
    pub fn viscosities(&self) -> &[f64] {
        &self.viscosity
    }

    /// Evaluate on the same-time physical density and velocity fields.
    /// Normalizing by density and velocity scales keeps the fourth-degree
    /// least-squares contractions representable. Subtracting a reference
    /// velocity before covariance avoids cancellation under uniform boosts.
    pub fn prepare(
        &mut self,
        config: DynamicSgsConfig,
        density: &[f64],
        velocity: &[[f64; 3]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DynamicSgsReport, DynamicSgsError> {
        poll(&mut cancelled, 0)?;
        config.validate()?;
        let n = self.viscosity.len();
        if density.len() != n || velocity.len() != n {
            return Err(DynamicSgsError::Invalid(
                "field and workspace lengths differ",
            ));
        }
        let reference = velocity[0];
        let mut density_scale: f64 = 0.0;
        let mut velocity_scale: f64 = 0.0;
        for cell in 0..n {
            poll(&mut cancelled, cell)?;
            if !density[cell].is_finite() || density[cell] <= 0.0 {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "positive finite density required",
                });
            }
            density_scale = density_scale.max(density[cell]);
            for (component, base) in reference.iter().enumerate() {
                let value = velocity[cell][component] - base;
                if !value.is_finite() {
                    return Err(DynamicSgsError::Cell {
                        cell,
                        reason: "finite representable velocity differences required",
                    });
                }
                self.velocity[cell][component] = value;
                velocity_scale = velocity_scale.max(value.abs());
            }
        }
        if velocity_scale == 0.0 {
            for (tile, values) in self.viscosity.chunks_mut(64).enumerate() {
                poll(&mut cancelled, tile * 64)?;
                values.fill(0.0);
            }
            poll(&mut cancelled, 0)?;
            return Ok(DynamicSgsReport {
                undetermined_cells: n,
                ..DynamicSgsReport::default()
            });
        }
        for (cell, v) in self.velocity.iter_mut().enumerate() {
            poll(&mut cancelled, cell)?;
            for component in v {
                *component /= velocity_scale;
            }
        }
        for cell in 0..n {
            poll(&mut cancelled, cell)?;
            let (rate, s) = strain(self.shape, cell, &self.velocity);
            self.strain_rate[cell] = rate;
            self.basis[cell] = s.map(|v| (density[cell] / density_scale) * rate * v);
            let mut filtered = [0.0; 4];
            each_box(self.shape, cell, |j| {
                let rho = density[j] / density_scale / 27.0;
                filtered[3] += rho;
                for (k, v) in filtered[..3].iter_mut().enumerate() {
                    *v += rho * self.velocity[j][k];
                }
            });
            if !(filtered[3] > 0.0) {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "test density underflow",
                });
            }
            for k in 0..3 {
                filtered[k] /= filtered[3];
            }
            self.filtered[cell] = filtered;
        }
        let mut numerator = Sum::default();
        let mut denominator = Sum::default();
        for cell in 0..n {
            poll(&mut cancelled, cell)?;
            let mut a = [0.0; 6];
            let mut l = [0.0; 6];
            let mean = self.filtered[cell];
            each_box(self.shape, cell, |j| {
                let rho = density[j] / density_scale / 27.0;
                let d = core::array::from_fn::<_, 3, _>(|k| self.velocity[j][k] - mean[k]);
                let product = [
                    d[0] * d[0],
                    d[1] * d[1],
                    d[2] * d[2],
                    d[0] * d[1],
                    d[0] * d[2],
                    d[1] * d[2],
                ];
                for k in 0..6 {
                    l[k] += rho * product[k];
                    a[k] += self.basis[j][k] / 27.0;
                }
            });
            let (rate, s) = strain(self.shape, cell, &self.filtered);
            let m = core::array::from_fn(|k| 2.0 * (a[k] - 9.0 * mean[3] * rate * s[k]));
            let pair = [contract(deviator(l), m), contract(m, m)];
            self.contraction[cell] = pair;
            numerator.add(pair[0]);
            denominator.add(pair[1]);
        }
        let domain = [numerator.value, denominator.value];
        let mut report = DynamicSgsReport {
            minimum_raw_coefficient_squared: f64::INFINITY,
            ..DynamicSgsReport::default()
        };
        for cell in 0..n {
            poll(&mut cancelled, cell)?;
            let pair = match config.averaging {
                DynamicAveraging::Domain => domain,
                DynamicAveraging::LocalBox => {
                    let mut pair = [0.0; 2];
                    each_box(self.shape, cell, |j| {
                        for k in 0..2 {
                            pair[k] += self.contraction[j][k] / 27.0;
                        }
                    });
                    pair
                }
            };
            let raw = if pair[1] == 0.0 {
                report.undetermined_cells += 1;
                0.0
            } else {
                pair[0] / pair[1]
            };
            if !raw.is_finite() {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "unrepresentable least-squares coefficient",
                });
            }
            report.minimum_raw_coefficient_squared =
                report.minimum_raw_coefficient_squared.min(raw);
            if raw < 0.0 {
                report.backscatter_limited_cells += 1;
            }
            let coefficient = raw.max(0.0);
            if coefficient > config.maximum_coefficient_squared {
                return Err(DynamicSgsError::Coefficient {
                    cell,
                    observed: coefficient,
                    maximum: config.maximum_coefficient_squared,
                });
            }
            let nu = coefficient * self.strain_rate[cell] * velocity_scale * self.spacing;
            if !nu.is_finite() {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "unrepresentable modeled viscosity",
                });
            }
            self.viscosity[cell] = nu;
            report.maximum_coefficient_squared =
                report.maximum_coefficient_squared.max(coefficient);
            report.maximum_eddy_viscosity = report.maximum_eddy_viscosity.max(nu);
        }
        poll(&mut cancelled, 0)?;
        Ok(report)
    }
}

#[cfg(test)]
#[path = "dynamic/tests.rs"]
mod tests;
