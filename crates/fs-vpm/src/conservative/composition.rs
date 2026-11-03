//! Signed composition of the existing reversible pair/source flows.

/// Time order of conservative inviscid evolution at a fixed step size.
///
/// Both choices use the same exact pair maps and retain their linear/angular
/// vortex-impulse invariants. Fourth order uses a negative internal interval;
/// it is available only for reversible inviscid and held-affine evolution,
/// never for PSE diffusion or remapping. It is not an adaptive error bound.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum CompositionOrder {
    /// The original symmetric forward/reverse pair sweep.
    #[default]
    Second,
    /// Three symmetric sweeps with Yoshida's fourth-order coefficients.
    Fourth,
}

impl CompositionOrder {
    pub(super) const fn sweeps(self) -> usize {
        match self {
            Self::Second => 1,
            Self::Fourth => 3,
        }
    }

    pub(super) const fn source_passes(self) -> u64 {
        match self {
            Self::Second => 2,
            Self::Fourth => 4,
        }
    }
}

// Yoshida, Phys. Lett. A 150 (1990), equations (4.1), (4.5)-(4.6),
// DOI 10.1016/0375-9601(90)90092-3: w1=1/(2-cbrt(2)), w0=1-2*w1.
// Thus 2*w1+w0=1 and 2*w1^3+w0^3=0. Fixed correctly rounded w1 avoids
// a platform-dependent cube root in the method's numerical identity.
const FORWARD: f64 = 1.351_207_191_959_657_5;
const REVERSE: f64 = 1.0 - 2.0 * FORWARD;

pub(super) struct Composition {
    weights: [f64; 3],
    count: usize,
}

impl Composition {
    pub(super) fn new(order: CompositionOrder) -> Self {
        let weights = match order {
            CompositionOrder::Second => [1.0, 0.0, 0.0],
            CompositionOrder::Fourth => [FORWARD, REVERSE, FORWARD],
        };
        Self {
            weights,
            count: order.sweeps(),
        }
    }

    // Do not premultiply dt: fractional subnormal time can otherwise round
    // away before circulation or source velocity makes the motion resolvable.
    pub(super) fn weights(&self) -> &[f64] {
        &self.weights[..self.count]
    }

    // Adjacent maps of one held affine field commute. Their half intervals
    // combine, giving four source passes for three pair sweeps rather than six.
    pub(super) fn bridge_weight(&self) -> Option<f64> {
        (self.count == 3).then(|| self.weights[0] + self.weights[1])
    }
}
