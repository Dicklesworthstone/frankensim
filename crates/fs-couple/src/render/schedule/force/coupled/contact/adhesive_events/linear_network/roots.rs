//! Bounded dense-polynomial event isolation, with every derivative extremum.
use super::{Error, ORDER, finite, poll};
use fs_exec::CancelGate;

pub(super) struct RootBudget {
    maximum: usize,
    iterations: usize,
    used: usize,
}
impl RootBudget {
    pub fn new(maximum: usize, iterations: usize) -> Self {
        Self {
            maximum,
            iterations,
            used: 0,
        }
    }
    pub fn reset(&mut self) {
        self.used = 0;
    }
    pub const fn used(&self) -> usize {
        self.used
    }
    fn eval(
        &mut self,
        c: &[f64; ORDER + 1],
        degree: usize,
        x: f64,
        gate: &CancelGate,
    ) -> Result<f64, Error> {
        if self.used >= self.maximum {
            return Err(Error::Budget("dense event evaluations"));
        }
        self.used += 1;
        if self.used % 32 == 1 {
            poll(gate)?;
        }
        finite(c[..=degree].iter().rev().fold(0.0, |v, a| v * x + a))
    }
}

#[derive(Clone, Copy)]
pub(super) struct EventRoot {
    pub fraction: f64,
    pub direction: i8,
}

struct Roots {
    values: [f64; ORDER],
    len: usize,
}
impl Roots {
    fn new() -> Self {
        Self {
            values: [0.0; ORDER],
            len: 0,
        }
    }
    fn add(&mut self, x: f64) -> Result<(), Error> {
        if self.len > 0 && x <= self.values[self.len - 1] {
            return Ok(());
        }
        if self.len == ORDER {
            return Err(Error::EventUnresolved);
        }
        self.values[self.len] = x;
        self.len += 1;
        Ok(())
    }
}

fn derivative(c: &[f64; ORDER + 1], degree: usize) -> [f64; ORDER + 1] {
    let mut d = [0.0; ORDER + 1];
    for i in 1..=degree {
        d[i - 1] = i as f64 * c[i];
    }
    d
}
fn degree(c: &[f64; ORDER + 1]) -> usize {
    let mut n = ORDER;
    while n > 0 && c[n] == 0.0 {
        n -= 1;
    }
    n
}
fn opposite(a: f64, b: f64) -> bool {
    a != 0.0 && b != 0.0 && a.is_sign_negative() != b.is_sign_negative()
}

fn bracket(
    c: &[f64; ORDER + 1],
    degree: usize,
    mut lo: f64,
    mut hi: f64,
    mut flo: f64,
    tolerance: f64,
    budget: &mut RootBudget,
    gate: &CancelGate,
) -> Result<f64, Error> {
    for _ in 0..budget.iterations {
        let mid = f64::midpoint(lo, hi);
        let fm = budget.eval(c, degree, mid, gate)?;
        if fm == 0.0 {
            return Ok(mid);
        }
        if hi - lo <= tolerance || mid == lo || mid == hi {
            return Ok(mid);
        }
        if opposite(flo, fm) {
            hi = mid;
        } else {
            lo = mid;
            flo = fm;
        }
    }
    Err(Error::Budget("event root iterations"))
}

/// Derivative recursion isolates all sign-changing polynomial roots. Even
/// derivative roots do not change monotonicity and need no extra partition.
fn roots(
    c: &[f64; ORDER + 1],
    mut n: usize,
    tolerance: f64,
    budget: &mut RootBudget,
    gate: &CancelGate,
) -> Result<Roots, Error> {
    while n > 0 && c[n] == 0.0 {
        n -= 1;
    }
    let mut result = Roots::new();
    if n == 0 {
        return Ok(result);
    }
    let critical = if n == 1 {
        Roots::new()
    } else {
        roots(&derivative(c, n), n - 1, tolerance, budget, gate)?
    };
    let mut left = 0.0;
    let mut fl = budget.eval(c, n, left, gate)?;
    if fl == 0.0 {
        result.add(0.0)?;
    }
    for i in 0..=critical.len {
        let right = if i == critical.len {
            1.0
        } else {
            critical.values[i]
        };
        if right <= left {
            continue;
        }
        let fr = budget.eval(c, n, right, gate)?;
        if opposite(fl, fr) {
            result.add(bracket(c, n, left, right, fl, tolerance, budget, gate)?)?;
        } else if fr == 0.0 {
            result.add(right)?;
        }
        left = right;
        fl = fr;
    }
    Ok(result)
}

pub(super) fn first_zero(
    c: &[f64; ORDER + 1],
    direction: i8,
    tolerance: f64,
    budget: &mut RootBudget,
    gate: &CancelGate,
) -> Result<Option<EventRoot>, Error> {
    let n = degree(c);
    if c[0] == 0.0 {
        let first = if n == 0 {
            0.0
        } else {
            c[1..=n].iter().find(|x| **x != 0.0).copied().unwrap_or(0.0)
        };
        if first == 0.0 || first * (direction as f64) < 0.0 {
            return Ok(Some(EventRoot {
                fraction: 0.0,
                direction,
            }));
        }
    }
    let all = roots(c, n, tolerance, budget, gate)?;
    let found = all.values[..all.len].iter().copied().find(|r| *r > 0.0);
    // A near-zero interior extremum without an isolated crossing can be an
    // unresolved grazing stop. Refuse instead of silently skipping it.
    if n > 1 {
        let extrema = roots(&derivative(c, n), n - 1, tolerance, budget, gate)?;
        let roundoff = 64.0 * f64::EPSILON * c.iter().map(|v| v.abs()).sum::<f64>();
        for &z in &extrema.values[..extrema.len] {
            if z <= 0.0 || z >= found.unwrap_or(1.0) {
                continue;
            }
            let value = budget.eval(c, n, z, gate)?;
            if value != 0.0 && value.abs() <= roundoff {
                return Err(Error::EventUnresolved);
            }
        }
    }
    Ok(found.map(|fraction| EventRoot {
        fraction,
        direction,
    }))
}

fn first_positive(
    c: &[f64; ORDER + 1],
    tolerance: f64,
    budget: &mut RootBudget,
    gate: &CancelGate,
) -> Result<Option<f64>, Error> {
    let n = degree(c);
    let extrema = if n == 0 {
        Roots::new()
    } else {
        roots(&derivative(c, n), n - 1, tolerance, budget, gate)?
    };
    let mut left = 0.0;
    let mut fl = budget.eval(c, n, left, gate)?;
    if fl > 0.0 {
        return Ok(Some(0.0));
    }
    for i in 0..=extrema.len {
        let right = if i == extrema.len {
            1.0
        } else {
            extrema.values[i]
        };
        if right <= left {
            continue;
        }
        let fr = budget.eval(c, n, right, gate)?;
        if fr > 0.0 {
            if fl == 0.0 {
                return Ok(Some(left));
            }
            return Ok(Some(bracket(
                c, n, left, right, fl, tolerance, budget, gate,
            )?));
        }
        left = right;
        fl = fr;
    }
    Ok(None)
}

pub(super) fn first_exit(
    c: &[f64; ORDER + 1],
    capacity: f64,
    tolerance: f64,
    budget: &mut RootBudget,
    gate: &CancelGate,
) -> Result<Option<EventRoot>, Error> {
    let mut upper = *c;
    upper[0] -= capacity;
    let mut lower = c.map(|v| -v);
    lower[0] -= capacity;
    let up = first_positive(&upper, tolerance, budget, gate)?;
    let down = first_positive(&lower, tolerance, budget, gate)?;
    Ok(match (up, down) {
        (None, None) => None,
        (Some(z), None) => Some(EventRoot {
            fraction: z,
            direction: 1,
        }),
        (None, Some(z)) => Some(EventRoot {
            fraction: z,
            direction: -1,
        }),
        (Some(a), Some(b)) => {
            if a <= b {
                Some(EventRoot {
                    fraction: a,
                    direction: 1,
                })
            } else {
                Some(EventRoot {
                    fraction: b,
                    direction: -1,
                })
            }
        }
    })
}
