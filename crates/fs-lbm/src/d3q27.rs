//! Full tensor-product D3Q27 central-moment and cumulant collisions.
//!
//! Geier, Schoenherr, Pasquali & Krafczyk (2015), CAMWA70:507--547,
//! doi:10.1016/j.camwa.2015.05.001, section4, equations43--96. All 27
//! independent moments are represented, including C111, C211 permutations,
//! C221 permutations and C222. The multivariate moment/cumulant recurrence
//! below is the log-generating-function definition, not a fourth-order
//! projection. Relaxation uses the paper's ten rotational groups and its
//! diagonal cubic-alias correction (58--63). No empirical rate defaults.
//!
//! Fixed tensor sweep and recurrence order; safe, allocation-free per cell.
//! Isothermal weakly-compressible model only. This does not supply total
//! energy, qualified high-Re jets, SGS, wall interpolation or thermal noise.

use core::fmt;

/// Version of tables, transforms, relaxation grouping and arithmetic order.
/// V2 supplies the third/fifth central force moments; unforced bits are unchanged.
pub const D3Q27_BIT_SEMANTICS_VERSION: u32 = 2;
/// Population count.
pub const Q27: usize = 27;
/// Lattice sound speed squared.
pub const CS2_27: f64 = 1.0 / 3.0;
/// Rest, faces, edges (same direction order as D3Q19), then corners.
pub const E27: [[i32; 3]; Q27] = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 0],
    [-1, -1, 0],
    [1, -1, 0],
    [-1, 1, 0],
    [1, 0, 1],
    [-1, 0, -1],
    [1, 0, -1],
    [-1, 0, 1],
    [0, 1, 1],
    [0, -1, -1],
    [0, 1, -1],
    [0, -1, 1],
    [1, 1, 1],
    [-1, -1, -1],
    [1, 1, -1],
    [-1, -1, 1],
    [1, -1, 1],
    [-1, 1, -1],
    [-1, 1, 1],
    [1, -1, -1],
];
/// Exact integer weights, denominator216.
pub const W216: [i64; Q27] = [
    64, 16, 16, 16, 16, 16, 16, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1,
];
/// Normalized weights, derived from the integer table.
pub const W27: [f64; Q27] = weights();
const fn weights() -> [f64; Q27] {
    let mut out = [0.0; Q27];
    let mut q = 0;
    while q < Q27 {
        out[q] = W216[q] as f64 / 216.0;
        q += 1;
    }
    out
}
/// Opposite directions.
pub const OPP27: [usize; Q27] = [
    0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26,
    25,
];

/// Explicit BGK equilibrium choice; used to isolate the lattice in controls.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Equilibrium27 {
    Quadratic,
    Factorized,
}

/// Geier omega1..omega10 in publication order: deviatoric stress, trace,
/// third-order sums/differences, C111, fourth-order diagonal differences/sum,
/// mixed fourth order, fifth order, sixth order. Every entry is required.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rates27 {
    pub omega: [f64; 10],
}

impl Rates27 {
    /// All rates must be finite and strictly in (0,2).
    pub fn validate(self) -> Result<(), CollisionError27> {
        for (group, rate) in self.omega.into_iter().enumerate() {
            if !(rate.is_finite() && rate > 0.0 && rate < 2.0) {
                return Err(CollisionError27::Rate {
                    group: group + 1,
                    rate,
                });
            }
        }
        Ok(())
    }
    /// Shear viscosity in lattice units; no bulk/higher rate is inferred.
    pub fn viscosity(self) -> Result<f64, CollisionError27> {
        self.validate()?;
        Ok(CS2_27 * (1.0 / self.omega[0] - 0.5))
    }
}

/// Separate, explicitly selected collision models on the complete lattice.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CollisionModel27 {
    Bgk {
        tau: f64,
        equilibrium: Equilibrium27,
    },
    CentralMoment {
        rates: Rates27,
    },
    Cumulant {
        rates: Rates27,
    },
}

/// Cell admission includes an explicit low-Mach envelope.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CollisionConfig27 {
    pub model: CollisionModel27,
    pub maximum_mach: f64,
}

impl CollisionConfig27 {
    pub fn validate(self) -> Result<(), CollisionError27> {
        if !(self.maximum_mach.is_finite()
            && self.maximum_mach > 0.0
            && self.maximum_mach <= crate::MACH_LIMIT)
        {
            return Err(CollisionError27::MachLimit {
                limit: self.maximum_mach,
            });
        }
        match self.model {
            CollisionModel27::Bgk { tau, .. } if tau.is_finite() && tau > 0.5 => Ok(()),
            CollisionModel27::Bgk { tau, .. } => Err(CollisionError27::Tau { tau }),
            CollisionModel27::CentralMoment { rates } | CollisionModel27::Cumulant { rates } => {
                rates.validate()
            }
        }
    }
}

/// Structured refusal; populations are never clamped or silently repaired.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CollisionError27 {
    Tau { tau: f64 },
    Rate { group: usize, rate: f64 },
    MachLimit { limit: f64 },
    Mach { observed: f64, limit: f64 },
    NonFinite { field: &'static str, index: usize },
    Density { rho: f64 },
}
impl fmt::Display for CollisionError27 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FS-LBM-D3Q27: {self:?}")
    }
}
impl std::error::Error for CollisionError27 {}

/// D3Q27 equilibrium with exact tensor-product second moments. The
/// factorized form retains mixed third through sixth equilibrium moments;
/// the quadratic form is an explicit same-polynomial D3Q19 control.
#[must_use]
pub fn equilibrium27(rho: f64, u: [f64; 3], form: Equilibrium27) -> [f64; Q27] {
    core::array::from_fn(|q| match form {
        Equilibrium27::Quadratic => {
            let eu = (0..3).map(|a| f64::from(E27[q][a]) * u[a]).sum::<f64>();
            let u2 = u.iter().map(|v| v * v).sum::<f64>();
            W27[q] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2)
        }
        Equilibrium27::Factorized => {
            let mut p = rho;
            for a in 0..3 {
                let second = CS2_27 + u[a] * u[a];
                p *= match E27[q][a] {
                    -1 => 0.5 * (second - u[a]),
                    0 => 1.0 - second,
                    _ => 0.5 * (second + u[a]),
                };
            }
            p
        }
    })
}

/// Density and physical lattice velocity, including the half-force shift.
pub fn macroscopic27(
    pop: &[f64; Q27],
    force: [f64; 3],
    maximum_mach: f64,
) -> Result<(f64, [f64; 3]), CollisionError27> {
    let mut rho = 0.0;
    let mut j = [0.0; 3];
    for (q, &v) in pop.iter().enumerate() {
        if !v.is_finite() {
            return Err(CollisionError27::NonFinite {
                field: "population",
                index: q,
            });
        }
        rho += v;
        for a in 0..3 {
            j[a] += f64::from(E27[q][a]) * v;
        }
    }
    if !(rho.is_finite() && rho > 0.0) {
        return Err(CollisionError27::Density { rho });
    }
    for (a, f) in force.into_iter().enumerate() {
        if !f.is_finite() {
            return Err(CollisionError27::NonFinite {
                field: "force",
                index: a,
            });
        }
    }
    let u = core::array::from_fn(|a| (j[a] + 0.5 * force[a]) / rho);
    let mach = (u.iter().map(|x| x * x).sum::<f64>() / CS2_27).sqrt();
    if !(maximum_mach.is_finite() && maximum_mach > 0.0 && maximum_mach <= crate::MACH_LIMIT) {
        return Err(CollisionError27::MachLimit {
            limit: maximum_mach,
        });
    }
    if !mach.is_finite() || mach > maximum_mach {
        return Err(CollisionError27::Mach {
            observed: mach,
            limit: maximum_mach,
        });
    }
    Ok((rho, u))
}

const fn at(a: usize, b: usize, c: usize) -> usize {
    a + 3 * b + 9 * c
}
const fn exponents(i: usize) -> [usize; 3] {
    [i % 3, (i / 3) % 3, i / 9]
}

/// Separable three-point transform. Coordinate and exponent indices both run
/// 0..2, with speed index0=-1,1=0,2=1. No 27x27 solve or fitted inverse.
fn transform(mut values: [f64; Q27], u: [f64; 3], inverse: bool) -> [f64; Q27] {
    let axes: [usize; 3] = if inverse { [2, 1, 0] } else { [0, 1, 2] };
    for axis in axes {
        let stride = [1, 3, 9][axis];
        for base in 0..Q27 {
            if exponents(base)[axis] != 0 {
                continue;
            }
            let [a, b, c] = [
                values[base],
                values[base + stride],
                values[base + 2 * stride],
            ];
            let result = if inverse {
                let raw1 = b + u[axis] * a;
                let raw2 = c + 2.0 * u[axis] * b + u[axis] * u[axis] * a;
                [0.5 * (raw2 - raw1), a - raw2, 0.5 * (raw2 + raw1)]
            } else {
                let zero = (a + b) + c;
                let first = c - a;
                [
                    zero,
                    first - u[axis] * zero,
                    (a + c) - 2.0 * u[axis] * first + u[axis] * u[axis] * zero,
                ]
            };
            for k in 0..3 {
                values[base + k * stride] = result[k];
            }
        }
    }
    values
}

fn choose(n: usize, k: usize) -> f64 {
    if n == 2 && k == 1 { 2.0 } else { 1.0 }
}

/// Normalized multivariate moment/cumulant recurrence from dM = M d(log M).
/// Componentwise truncation at (2,2,2) is closed under this recurrence.
/// mu000=1, k000=0; density stays outside the logarithm. First-order entries
/// are retained, including the nonzero half-force frame shift.
fn cumulants(mu: [f64; Q27]) -> [f64; Q27] {
    let mut k = [0.0; Q27];
    for i in 1..Q27 {
        let alpha = exponents(i);
        let d = (0..3).find(|&a| alpha[a] > 0).unwrap_or(0);
        let mut n = alpha;
        n[d] -= 1;
        let mut value = mu[i];
        for z in 0..=n[2] {
            for y in 0..=n[1] {
                for x in 0..=n[0] {
                    let beta = [x, y, z];
                    if beta == n {
                        continue;
                    }
                    let mut b = beta;
                    b[d] += 1;
                    let coefficient = choose(n[0], x) * choose(n[1], y) * choose(n[2], z);
                    value -= coefficient
                        * k[at(b[0], b[1], b[2])]
                        * mu[at(n[0] - x, n[1] - y, n[2] - z)];
                }
            }
        }
        k[i] = value;
    }
    k
}

fn moments(k: [f64; Q27]) -> [f64; Q27] {
    let mut mu = [0.0; Q27];
    mu[0] = 1.0;
    for i in 1..Q27 {
        let alpha = exponents(i);
        let d = (0..3).find(|&a| alpha[a] > 0).unwrap_or(0);
        let mut n = alpha;
        n[d] -= 1;
        let mut value = 0.0;
        for z in 0..=n[2] {
            for y in 0..=n[1] {
                for x in 0..=n[0] {
                    let mut b = [x, y, z];
                    b[d] += 1;
                    let coefficient = choose(n[0], x) * choose(n[1], y) * choose(n[2], z);
                    value += coefficient
                        * k[at(b[0], b[1], b[2])]
                        * mu[at(n[0] - x, n[1] - y, n[2] - z)];
                }
            }
        }
        mu[i] = value;
    }
    mu
}

fn relax_pair(v: &mut [f64; Q27], a: usize, b: usize, sum_rate: f64, diff_rate: f64) {
    let sum = (1.0 - sum_rate) * (v[a] + v[b]);
    let diff = (1.0 - diff_rate) * (v[a] - v[b]);
    v[a] = 0.5 * (sum + diff);
    v[b] = 0.5 * (sum - diff);
}

fn relax(mut v: [f64; Q27], eq: [f64; Q27], u: [f64; 3], rates: Rates27) -> [f64; Q27] {
    let w = rates.omega;
    // Normalized forms of Geier58--63, using pre-collision values.
    let [xx, yy, zz] = [at(2, 0, 0), at(0, 2, 0), at(0, 0, 2)];
    let trace = v[xx] + v[yy] + v[zz];
    let dx = -0.5 * w[0] * (2.0 * v[xx] - v[yy] - v[zz]) - 0.5 * w[1] * (trace - 1.0);
    let dy = dx + 1.5 * w[0] * (v[xx] - v[yy]);
    let dz = dx + 1.5 * w[0] * (v[xx] - v[zz]);
    let g = [u[0] * u[0] * dx, u[1] * u[1] * dy, u[2] * u[2] * dz];
    let dxy = (1.0 - w[0]) * (v[xx] - v[yy]) - 3.0 * (1.0 - 0.5 * w[0]) * (g[0] - g[1]);
    let dxz = (1.0 - w[0]) * (v[xx] - v[zz]) - 3.0 * (1.0 - 0.5 * w[0]) * (g[0] - g[2]);
    let tr = w[1] + (1.0 - w[1]) * trace - 3.0 * (1.0 - 0.5 * w[1]) * (g[0] + g[1] + g[2]);
    v[xx] = (tr + dxy + dxz) / 3.0;
    v[yy] = v[xx] - dxy;
    v[zz] = v[xx] - dxz;
    for i in [at(1, 1, 0), at(1, 0, 1), at(0, 1, 1)] {
        v[i] *= 1.0 - w[0];
    }
    for (a, b) in [
        (at(1, 2, 0), at(1, 0, 2)),
        (at(2, 1, 0), at(0, 1, 2)),
        (at(2, 0, 1), at(0, 2, 1)),
    ] {
        relax_pair(&mut v, a, b, w[2], w[3]);
    }
    v[at(1, 1, 1)] *= 1.0 - w[4];
    // Both independent diagonal fourth-order differences share omega6.
    let ids = [at(2, 2, 0), at(2, 0, 2), at(0, 2, 2)];
    let mean = ids.iter().map(|&i| v[i] - eq[i]).sum::<f64>() / 3.0;
    for i in ids {
        v[i] = eq[i] + (1.0 - w[5]) * (v[i] - eq[i] - mean) + (1.0 - w[6]) * mean;
    }
    for i in [at(2, 1, 1), at(1, 2, 1), at(1, 1, 2)] {
        v[i] *= 1.0 - w[7];
    }
    for i in [at(2, 2, 1), at(2, 1, 2), at(1, 2, 2)] {
        v[i] *= 1.0 - w[8];
    }
    v[at(2, 2, 2)] = eq[at(2, 2, 2)] + (1.0 - w[9]) * (v[at(2, 2, 2)] - eq[at(2, 2, 2)]);
    v
}

/// Collision with explicit force density in lattice units. Cumulant forcing
/// changes the first cumulants from -F/(2rho) to +F/(2rho) in the midpoint
/// frame (Geier85--87). Central-moment forcing also supplies the third/fifth
/// Gaussian source moments with their own grouped trapezoidal prefactors
/// (Fei, Luo & Li2018, PRE97:053309, equations5,18).
/// No positivity limiter or automatic rate selection.
/// Finite-force Navier--Stokes accuracy requires its own downstream ladder.
pub fn collide_cell27(
    pop: [f64; Q27],
    config: CollisionConfig27,
    force: [f64; 3],
) -> Result<[f64; Q27], CollisionError27> {
    config.validate()?;
    let (rho, u) = macroscopic27(&pop, force, config.maximum_mach)?;
    let post = match config.model {
        CollisionModel27::Bgk { tau, equilibrium } => {
            let eq = equilibrium27(rho, u, equilibrium);
            core::array::from_fn(|q| {
                let e = E27[q].map(f64::from);
                let eu = (0..3).map(|a| e[a] * u[a]).sum::<f64>();
                let source = (0..3)
                    .map(|a| ((e[a] - u[a]) / CS2_27 + eu * e[a] / (CS2_27 * CS2_27)) * force[a])
                    .sum::<f64>();
                (1.0 / tau).mul_add(eq[q] - pop[q], pop[q]) + (1.0 - 0.5 / tau) * W27[q] * source
            })
        }
        CollisionModel27::CentralMoment { rates } | CollisionModel27::Cumulant { rates } => {
            let mut tensor = [0.0; Q27];
            for q in 0..Q27 {
                let e = E27[q].map(|c| (c + 1) as usize);
                tensor[at(e[0], e[1], e[2])] = pop[q] / rho;
            }
            let mu = transform(tensor, u, false);
            let mut gaussian = [0.0; Q27];
            for (i, v) in gaussian.iter_mut().enumerate() {
                let p = exponents(i);
                if p.iter().all(|&a| a % 2 == 0) {
                    *v = p
                        .into_iter()
                        .filter(|&a| a == 2)
                        .fold(1.0, |v, _| v * CS2_27);
                }
            }
            let mut updated = match config.model {
                CollisionModel27::Cumulant { .. } => {
                    let mut equilibrium = [0.0; Q27];
                    for index in [at(2, 0, 0), at(0, 2, 0), at(0, 0, 2)] {
                        equilibrium[index] = CS2_27;
                    }
                    let mut k = relax(cumulants(mu), equilibrium, u, rates);
                    for (a, index) in [at(1, 0, 0), at(0, 1, 0), at(0, 0, 1)]
                        .into_iter()
                        .enumerate()
                    {
                        k[index] = 0.5 * force[a] / rho;
                    }
                    moments(k)
                }
                _ => {
                    let mut central = relax(mu, gaussian, u, rates);
                    if force.iter().any(|&f| f != 0.0) {
                        // The Gaussian source has C120=C102=cs² Fx/rho,
                        // its axis permutations, and C122=cs⁴ Fx/rho etc.
                        // Third-order pair differences have zero source, so
                        // omega3 (the sum rate), not omega4, sets the prefactor.
                        // Retain the untouched path for exactly zero force.
                        let third = (1.0 - 0.5 * rates.omega[2]) * CS2_27;
                        let fifth = (1.0 - 0.5 * rates.omega[8]) * CS2_27 * CS2_27;
                        for (axis, (a, b, c)) in [
                            (at(1, 2, 0), at(1, 0, 2), at(1, 2, 2)),
                            (at(2, 1, 0), at(0, 1, 2), at(2, 1, 2)),
                            (at(2, 0, 1), at(0, 2, 1), at(2, 2, 1)),
                        ]
                        .into_iter()
                        .enumerate()
                        {
                            let acceleration = force[axis] / rho;
                            central[a] += third * acceleration;
                            central[b] += third * acceleration;
                            central[c] += fifth * acceleration;
                        }
                    }
                    central
                }
            };
            updated[0] = 1.0;
            for (a, index) in [at(1, 0, 0), at(0, 1, 0), at(0, 0, 1)]
                .into_iter()
                .enumerate()
            {
                updated[index] = 0.5 * force[a] / rho;
            }
            let tensor = transform(updated, u, true);
            core::array::from_fn(|q| {
                let e = E27[q].map(|c| (c + 1) as usize);
                rho * tensor[at(e[0], e[1], e[2])]
            })
        }
    };
    macroscopic27(&post, force, config.maximum_mach)?;
    Ok(post)
}

/// Marker for the shared lattice-generic grid.
#[derive(Clone, Copy, Debug)]
pub struct D3q27;
impl crate::lattice3::Lattice3<Q27> for D3q27 {
    type Collision = CollisionConfig27;
    type Error = CollisionError27;
    const E: [[i32; 3]; Q27] = E27;
    const W: [f64; Q27] = W27;
    const OPP: [usize; Q27] = OPP27;
    const CS2: f64 = CS2_27;
    fn validate(c: Self::Collision) -> Result<(), Self::Error> {
        c.validate()
    }
    fn collide(f: [f64; Q27], c: Self::Collision, g: [f64; 3]) -> Result<[f64; Q27], Self::Error> {
        collide_cell27(f, c, g)
    }
}

/// D3Q27 consumer alias. Construction always requires a complete config.
pub type Grid27 = crate::lattice3::Grid3<D3q27, Q27>;
