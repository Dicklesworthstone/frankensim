//! Conservative scalar transport shared by the D2Q5 and D3Q7 slabs.
//!
//! Linear equilibrium alone recovers an unwanted
//! div[(tau - 1/2) d(T u)/dt] term. Add the zero-mass first-moment source
//! of Chopard, Falcone & Latt (2009), section 2.4, equation (14):
//! w_i (1 - 1/(2 tau)) e_i . [T u - (T u)_previous] / c_T^2.
//! <https://doi.org/10.1140/epjst/e2009-01035-5>
//!
//! All velocities belong to the beginning of the thermal interval. A fresh
//! history uses zero correction for one startup step; no fictitious previous
//! temperature or acceleration is extrapolated. This does not establish the
//! temporal order of the complete partitioned buoyancy/source solver.

use core::ops::Range;

/// One vector per cell; retained by transactional snapshots, never by a
/// step-local relaxation field. Wall entries are unused.
#[derive(Clone, Debug, Default)]
pub(crate) struct FluxHistory<const D: usize> {
    flux: Vec<[f64; D]>,
    ready: bool,
}

impl<const D: usize> FluxHistory<D> {
    pub(crate) fn reset(&mut self) {
        self.ready = false;
    }

    fn prepare(&mut self, cells: usize) {
        if self.flux.len() != cells {
            self.flux.resize(cells, [0.0; D]);
            self.ready = false;
        }
    }
}

/// Rest plus one opposite pair per coordinate axis. The two supported
/// lattices have w_plus = w_minus = c_T^2 / 2.
pub(crate) struct AxisLattice<const D: usize, const Q: usize> {
    pub(crate) weights: [f64; Q],
    pub(crate) cs2: f64,
    positive: [usize; D],
    negative: [usize; D],
}

pub(crate) const D2Q5: AxisLattice<2, 5> = AxisLattice {
    weights: [1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0],
    cs2: 1.0 / 3.0,
    positive: [1, 2],
    negative: [3, 4],
};

pub(crate) const D3Q7: AxisLattice<3, 7> = AxisLattice {
    weights: [0.25, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125],
    cs2: 0.25,
    positive: [1, 3, 5],
    negative: [2, 4, 6],
};

impl<const D: usize, const Q: usize> AxisLattice<D, Q> {
    pub(crate) fn equilibrium(&self, temperature: f64, velocity: [f64; D]) -> [f64; Q] {
        let mut populations = [0.0; Q];
        populations[0] = self.weights[0] * temperature;
        for axis in 0..D {
            let positive = self.positive[axis];
            let negative = self.negative[axis];
            let mean = self.weights[positive] * temperature;
            let advective = 0.5 * temperature * velocity[axis];
            populations[positive] = mean + advective;
            populations[negative] = mean - advective;
        }
        populations
    }

    /// Write tightly packed post-collision populations into the original flow
    /// scratch. Opposite correction sources cancel without a thermal source
    /// or a face allocation. Rates are the admitted, step-local rates.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn collide(
        &self,
        populations: &[[f64; Q]],
        post: &mut [f64],
        active: Range<usize>,
        background_tau: f64,
        local_tau: Option<&[f64]>,
        history: &mut FluxHistory<D>,
        velocity: impl Fn(usize) -> [f64; D],
    ) {
        history.prepare(populations.len());
        for cell in active {
            let before = &populations[cell];
            let temperature: f64 = before.iter().sum();
            let u = velocity(cell);
            let flux: [f64; D] = core::array::from_fn(|axis| temperature * u[axis]);
            let tau = local_tau.map_or(background_tau, |field| field[cell]);
            let half_prefactor = 0.5 * ((tau - 0.5) / tau);
            let out = &mut post[Q * cell..Q * (cell + 1)];
            out[0] = before[0] + (self.weights[0] * temperature - before[0]) / tau;
            for axis in 0..D {
                let positive = self.positive[axis];
                let negative = self.negative[axis];
                let mean = self.weights[positive] * temperature;
                let advective = 0.5 * flux[axis];
                let correction = if history.ready {
                    half_prefactor * (flux[axis] - history.flux[cell][axis])
                } else {
                    0.0
                };
                out[positive] =
                    before[positive] + (mean + advective - before[positive]) / tau + correction;
                out[negative] =
                    before[negative] + (mean - advective - before[negative]) / tau - correction;
            }
            history.flux[cell] = flux;
        }
        history.ready = true;
    }

    /// Tensor counterpart of the scalar collision. The zero-mass pair
    /// correction replaces only the first-moment relaxation block. It is
    /// Eq. (5), with K2=k0, and Eq. (13) of Zhao et al. (2020):
    /// <https://arxiv.org/abs/1902.09813>.
    /// Existing scalar calls keep their original arithmetic path.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn collide_tensor(
        &self,
        populations: &[[f64; Q]],
        post: &mut [f64],
        active: Range<usize>,
        background_tau: f64,
        diffusion: super::tensor::DiffusivityField<'_, D>,
        history: &mut FluxHistory<D>,
        velocity: impl Fn(usize) -> [f64; D],
    ) {
        history.prepare(populations.len());
        let even_rate = 1.0 / background_tau;
        for cell in active {
            let before = &populations[cell];
            let temperature: f64 = before.iter().sum();
            let u = velocity(cell);
            let flux: [f64; D] = core::array::from_fn(|axis| temperature * u[axis]);
            let nonequilibrium: [f64; D] = core::array::from_fn(|axis| {
                before[self.positive[axis]] - before[self.negative[axis]] - flux[axis]
            });
            let derivative: [f64; D] = core::array::from_fn(|axis| {
                if history.ready { flux[axis] - history.flux[cell][axis] } else { 0.0 }
            });
            let relaxation = &diffusion.at(cell).relaxation;
            let out = &mut post[Q * cell..Q * (cell + 1)];
            out[0] = before[0] + (self.weights[0] * temperature - before[0]) / background_tau;
            for axis in 0..D {
                let mut correction = even_rate * nonequilibrium[axis] + derivative[axis];
                for other in 0..D {
                    correction -= relaxation[axis][other]
                        * (nonequilibrium[other] + 0.5 * derivative[other]);
                }
                correction *= 0.5;
                let positive = self.positive[axis];
                let negative = self.negative[axis];
                let mean = self.weights[positive] * temperature;
                let advective = 0.5 * flux[axis];
                out[positive] = before[positive]
                    + (mean + advective - before[positive]) / background_tau + correction;
                out[negative] = before[negative]
                    + (mean - advective - before[negative]) / background_tau - correction;
            }
            history.flux[cell] = flux;
        }
        history.ready = true;
    }

    /// Pull streaming, periodic in the in-plane coordinates. If walls are
    /// supplied, the first/last planes are fixed halfway-temperature walls;
    /// otherwise all axes wrap. The same pair operations serve both cases.
    pub(crate) fn stream(
        &self,
        populations: &mut [[f64; Q]],
        post: &[f64],
        shape: [usize; D],
        walls: Option<[f64; 2]>,
    ) {
        let mut stride = [1; D];
        for axis in 1..D {
            stride[axis] = stride[axis - 1] * shape[axis - 1];
        }
        let cells = populations.len();
        let plane = stride[D - 1];
        let active = if walls.is_some() {
            plane..cells - plane
        } else {
            0..cells
        };
        for cell in active {
            let mut incoming = [0.0; Q];
            incoming[0] = post[Q * cell];
            for axis in 0..D {
                let positive = self.positive[axis];
                let negative = self.negative[axis];
                if let Some([bottom, top]) = walls.filter(|_| axis == D - 1) {
                    incoming[positive] = if cell < 2 * plane {
                        2.0f64.mul_add(self.weights[positive] * bottom, -post[Q * cell + negative])
                    } else {
                        post[Q * (cell - plane) + positive]
                    };
                    incoming[negative] = if cell >= cells - 2 * plane {
                        2.0f64.mul_add(self.weights[negative] * top, -post[Q * cell + positive])
                    } else {
                        post[Q * (cell + plane) + negative]
                    };
                } else {
                    let coordinate = (cell / stride[axis]) % shape[axis];
                    let from_negative = if coordinate == 0 {
                        cell + (shape[axis] - 1) * stride[axis]
                    } else {
                        cell - stride[axis]
                    };
                    let from_positive = if coordinate + 1 == shape[axis] {
                        cell - (shape[axis] - 1) * stride[axis]
                    } else {
                        cell + stride[axis]
                    };
                    incoming[positive] = post[Q * from_negative + positive];
                    incoming[negative] = post[Q * from_positive + negative];
                }
            }
            populations[cell] = incoming;
        }
    }

    /// One prescribed source half-kick at fixed physical velocity.
    ///
    /// Adding equilibrium(delta T, u), and delta T * u to the retained
    /// history, preserves both g - g_eq and T*u - history. Isotropic kicks
    /// would change the nonconserved flux; unshifted history would mistake
    /// the source increment for transport. This is a source-substep rebind,
    /// not a claim of second-order coupled time integration.
    pub(crate) fn source_half_kick(
        &self,
        populations: &mut [[f64; Q]],
        heat: &[f64],
        active: Range<usize>,
        history: &mut FluxHistory<D>,
        velocity: impl Fn(usize) -> [f64; D],
    ) {
        history.prepare(populations.len());
        for cell in active {
            if heat[cell] == 0.0 {
                continue;
            }
            let increment = 0.5 * heat[cell];
            let u = velocity(cell);
            let kick = self.equilibrium(increment, u);
            for (population, delta) in populations[cell].iter_mut().zip(kick) {
                *population += delta;
            }
            if history.ready {
                for axis in 0..D {
                    history.flux[cell][axis] += increment * u[axis];
                }
            }
        }
    }
}

#[cfg(test)]
#[path = "advection/tests.rs"]
mod tests;

#[cfg(test)]
#[path = "advection/tensor_tests.rs"]
mod tensor_tests;
