//! Fixed analytic derivatives for two symmetric source directions and one
//! receiving direction. This is a 40-scalar local algebra, not a numerical
//! source stencil or an allocation-backed automatic differentiation graph.
use super::*;

pub(super) const PAIRS: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

#[derive(Clone, Copy, Debug)]
pub(super) struct Derivatives {
    pub value: C64,
    pub source: [C64; 3],
    pub source2: [C64; 6],
    pub receiver: [C64; 3],
    pub mixed: [[C64; 3]; 3],
    pub mixed2: [[C64; 6]; 3],
}
impl Derivatives {
    pub fn constant(value: C64) -> Self {
        Self {
            value,
            source: [C64::ZERO; 3],
            source2: [C64::ZERO; 6],
            receiver: [C64::ZERO; 3],
            mixed: [[C64::ZERO; 3]; 3],
            mixed2: [[C64::ZERO; 6]; 3],
        }
    }
    pub fn real(value: f64) -> Self {
        Self::constant(C64::from_re(value))
    }
    pub fn coordinate(value: f64, source: [f64; 3], receiver: [f64; 3]) -> Self {
        let mut result = Self::real(value);
        result.source = source.map(C64::from_re);
        result.receiver = receiver.map(C64::from_re);
        result
    }
    pub fn map(self, f: impl Fn(C64) -> C64) -> Self {
        Self {
            value: f(self.value),
            source: self.source.map(&f),
            source2: self.source2.map(&f),
            receiver: self.receiver.map(&f),
            mixed: self.mixed.map(|row| row.map(&f)),
            mixed2: self.mixed2.map(|row| row.map(&f)),
        }
    }
    pub fn scale(self, value: f64) -> Self {
        self.map(|v| v.scale(value))
    }
    pub fn times(self, value: C64) -> Self {
        self.map(|v| v * value)
    }
    pub fn chain(self, value: C64, first: C64, second: C64, third: C64) -> Self {
        let mut result = Self::constant(value);
        for j in 0..3 {
            result.source[j] = first * self.source[j];
            result.receiver[j] = first * self.receiver[j];
            for i in 0..3 {
                result.mixed[i][j] =
                    first * self.mixed[i][j] + second * self.receiver[i] * self.source[j];
            }
        }
        for (a, (j, k)) in PAIRS.into_iter().enumerate() {
            result.source2[a] = first * self.source2[a] + second * self.source[j] * self.source[k];
            for i in 0..3 {
                result.mixed2[i][a] = first * self.mixed2[i][a]
                    + second
                        * (self.receiver[i] * self.source2[a]
                            + self.mixed[i][j] * self.source[k]
                            + self.mixed[i][k] * self.source[j])
                    + third * self.receiver[i] * self.source[j] * self.source[k];
            }
        }
        result
    }
    pub fn recip(self) -> Self {
        let r = self.value.recip();
        let r2 = r * r;
        self.chain(r, -r2, (r2 * r).scale(2.0), (r2 * r2).scale(-6.0))
    }
    pub fn sqrt(self) -> Self {
        let r = self.value.sqrt();
        let inverse = r.recip();
        let inverse2 = inverse * inverse;
        self.chain(
            r,
            inverse.scale(0.5),
            (inverse * inverse2).scale(-0.25),
            (inverse * inverse2 * inverse2).scale(0.375),
        )
    }
    pub fn argument_derivatives(self) -> Self {
        let r = self.value.recip();
        self.chain(C64::ZERO, r, -r * r, (r * r * r).scale(2.0))
            .map(|z| C64::from_re(z.im))
    }
    pub fn all_finite(self) -> bool {
        finite(self.value)
            && self.source.into_iter().all(finite)
            && self.source2.into_iter().all(finite)
            && self.receiver.into_iter().all(finite)
            && self.mixed.into_iter().flatten().all(finite)
            && self.mixed2.into_iter().flatten().all(finite)
    }
}
impl core::ops::Add for Derivatives {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        let mut r = Self::constant(self.value + b.value);
        for j in 0..3 {
            r.source[j] = self.source[j] + b.source[j];
            r.receiver[j] = self.receiver[j] + b.receiver[j];
            for i in 0..3 {
                r.mixed[i][j] = self.mixed[i][j] + b.mixed[i][j];
            }
        }
        for j in 0..6 {
            r.source2[j] = self.source2[j] + b.source2[j];
            for i in 0..3 {
                r.mixed2[i][j] = self.mixed2[i][j] + b.mixed2[i][j];
            }
        }
        r
    }
}
impl core::ops::Sub for Derivatives {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + b.scale(-1.0)
    }
}
impl core::ops::Mul for Derivatives {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        let a = self;
        let mut r = Self::constant(a.value * b.value);
        for j in 0..3 {
            r.source[j] = a.source[j] * b.value + a.value * b.source[j];
            r.receiver[j] = a.receiver[j] * b.value + a.value * b.receiver[j];
            for i in 0..3 {
                r.mixed[i][j] = a.mixed[i][j] * b.value
                    + a.receiver[i] * b.source[j]
                    + a.source[j] * b.receiver[i]
                    + a.value * b.mixed[i][j];
            }
        }
        for (v, (j, k)) in PAIRS.into_iter().enumerate() {
            r.source2[v] = a.source2[v] * b.value
                + a.source[j] * b.source[k]
                + a.source[k] * b.source[j]
                + a.value * b.source2[v];
            for i in 0..3 {
                r.mixed2[i][v] = a.mixed2[i][v] * b.value
                    + a.source2[v] * b.receiver[i]
                    + a.mixed[i][j] * b.source[k]
                    + a.source[j] * b.mixed[i][k]
                    + a.mixed[i][k] * b.source[j]
                    + a.source[k] * b.mixed[i][j]
                    + a.receiver[i] * b.source2[v]
                    + a.value * b.mixed2[i][v];
            }
        }
        r
    }
}
