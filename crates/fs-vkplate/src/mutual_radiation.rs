//! Reciprocal mutual radiation resistance between distinct mechanical bodies.
//!
//! This is an explicitly declared narrowband approximation: one complete
//! symmetric Re Z sample, on the producer's original generalized coordinates.
//! It retains every signed cross term. It is not a causal fit of a frequency
//! sweep, a retarded load, an added-mass model, or a mutual-scattering solver.
//!
//! For a length coordinate, Z is kg/s. Given surface velocity shapes Phi and
//! area quadrature W, the work-conjugate nodal projection is
//! P = (Phi^T W Phi)^-1 Phi^T W. The caller supplies those rows, not Phi itself.
//! Then qdot = P v and f = -P^T R qdot, so nodal work is -qdot^T R qdot.
//! A BEM ModalLoadingRow supplies omega/(2 pi) and its unchanged real impedance;
//! the declared density and sound speed must be the producer's medium.

use fs_math::dd::Dd;
use std::collections::{BTreeMap, BTreeSet};

/// Deliberate dense-port and mapped-operator limits; no port is truncated.
pub const MAX_RESISTANCE_PORTS: usize = 128;
pub const MAX_FOOTPRINT_TERMS: usize = 65_536;

/// One nodal projection term: body, node, dimensionless generalized-coordinate weight.
pub type RadiationFootprint = Vec<(usize, usize, f64)>;

/// Complete, immutable narrowband resistance with a positive-semidefinite factor.
#[derive(Clone, Debug)]
pub struct ResistanceSample {
    ports: usize,
    reference_hz: f64,
    density: f64,
    sound_speed: f64,
    resistance: Vec<f64>,
    rows: Vec<Vec<Dd>>,
    pivots: Vec<Dd>,
}

fn finite(x: Dd) -> bool {
    x.hi.is_finite() && x.lo.is_finite()
}

fn invalid(what: &str) -> String {
    format!("mutual radiation: {what}")
}

impl ResistanceSample {
    /// Admit the full row-major Re Z in its original work-conjugate basis.
    ///
    /// Exact symmetry is required. Pivoted double-double LDL admits zero
    /// eigenvalues without adding a loss floor. A negative or unresolved
    /// Schur complement refuses; no eigenvalue clipping, averaging of
    /// nonreciprocal entries, or diagonal-only fallback changes the input.
    /// The factor retains about 106 significand bits; arithmetic outside that
    /// finite range refuses. A frequency sample alone does not certify a
    /// broadband causal impedance.
    pub fn new(
        ports: usize,
        reference_hz: f64,
        density: f64,
        sound_speed: f64,
        resistance: Vec<f64>,
    ) -> Result<Self, String> {
        if !(1..=MAX_RESISTANCE_PORTS).contains(&ports) || resistance.len() != ports * ports {
            return Err(invalid(
                "requires 1..128 ports and the complete square matrix",
            ));
        }
        if [reference_hz, density, sound_speed]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(invalid(
                "frequency, density and sound speed must be finite and positive",
            ));
        }
        for i in 0..ports {
            for j in 0..ports {
                if !resistance[i * ports + j].is_finite()
                    || resistance[i * ports + j] != resistance[j * ports + i]
                {
                    return Err(invalid("resistance must be finite and exactly symmetric"));
                }
            }
        }
        let scale = resistance.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
        let mut rows = Vec::new();
        let mut pivots = Vec::new();
        if scale > 0.0 {
            let scale_dd = Dd::from_f64(scale);
            let mut a: Vec<Dd> = resistance
                .iter()
                .map(|x| Dd::from_f64(*x) / scale_dd)
                .collect();
            if a.iter()
                .zip(&resistance)
                .any(|(a, r)| !finite(*a) || (*r != 0.0 && a.hi == 0.0))
            {
                return Err(invalid(
                    "matrix dynamic range exceeds finite factor arithmetic",
                ));
            }
            let mut active: Vec<usize> = (0..ports).collect();
            while !active.is_empty() {
                let mut at = 0;
                for k in 1..active.len() {
                    if a[active[at] * ports + active[at]].lt(a[active[k] * ports + active[k]]) {
                        at = k;
                    }
                }
                let p = active[at];
                let d = a[p * ports + p];
                if d.lt(Dd::ZERO) {
                    return Err(invalid(
                        "negative resistance pivot; sampled load is not passive",
                    ));
                }
                if d == Dd::ZERO {
                    if active
                        .iter()
                        .any(|i| active.iter().any(|j| a[i * ports + j] != Dd::ZERO))
                    {
                        return Err(invalid("zero resistance pivot has a nonzero Schur block"));
                    }
                    break;
                }
                let mut row = vec![Dd::ZERO; ports];
                for &j in &active {
                    row[j] = a[p * ports + j] / d;
                    if !finite(row[j]) {
                        return Err(invalid("unresolved semidefinite resistance factor"));
                    }
                }
                row[p] = Dd::ONE;
                active.swap_remove(at);
                for (ii, &i) in active.iter().enumerate() {
                    for &j in &active[ii..] {
                        let next = a[i * ports + j] - row[i] * a[p * ports + j];
                        if !finite(next) {
                            return Err(invalid("resistance factor exceeds finite arithmetic"));
                        }
                        a[i * ports + j] = next;
                        a[j * ports + i] = next;
                    }
                }
                let pivot = d * scale_dd;
                if !finite(pivot) || pivot.hi <= 0.0 {
                    return Err(invalid("resistance pivot exceeds finite arithmetic"));
                }
                rows.push(row);
                pivots.push(pivot);
            }
        }
        Ok(Self {
            ports,
            reference_hz,
            density,
            sound_speed,
            resistance,
            rows,
            pivots,
        })
    }

    pub const fn ports(&self) -> usize {
        self.ports
    }
    pub const fn reference_hz(&self) -> f64 {
        self.reference_hz
    }
    pub const fn density(&self) -> f64 {
        self.density
    }
    pub const fn sound_speed(&self) -> f64 {
        self.sound_speed
    }
    /// Original complete matrix, including negative cross terms.
    pub fn resistance(&self) -> &[f64] {
        &self.resistance
    }
    /// Number of admitted dissipative channels; zero is an exact dark sample.
    pub fn rank(&self) -> usize {
        self.pivots.len()
    }

    pub(crate) fn bind(
        self,
        footprints: &[RadiationFootprint],
        mass: impl Fn(usize, usize) -> Option<f64>,
    ) -> Result<BoundResistance, String> {
        if footprints.len() != self.ports
            || footprints.iter().any(Vec::is_empty)
            || footprints
                .iter()
                .try_fold(0usize, |n, p| n.checked_add(p.len()))
                .is_none_or(|n| n > MAX_FOOTPRINT_TERMS)
        {
            return Err(invalid(
                "requires a nonempty bounded footprint for every port",
            ));
        }
        let mut nodes = BTreeSet::new();
        for footprint in footprints {
            for &(body, node, weight) in footprint {
                if !weight.is_finite()
                    || mass(body, node).is_none_or(|m| !m.is_finite() || m <= 0.0)
                {
                    return Err(invalid("invalid node, projection weight or nodal mass"));
                }
                nodes.insert((body, node));
            }
        }
        let nodes: Vec<_> = nodes.into_iter().collect();
        let index: BTreeMap<_, _> = nodes
            .iter()
            .copied()
            .enumerate()
            .map(|(i, n)| (n, i))
            .collect();
        // Retain sparse P; expanding every factor over all bodies would make
        // callback work scale as rank * all_nodes for disjoint source ports.
        let projections: Vec<Vec<(usize, Dd)>> = footprints
            .iter()
            .map(|p| {
                p.iter()
                    .map(|&(b, n, w)| (index[&(b, n)], Dd::from_f64(w)))
                    .collect()
            })
            .collect();
        // diag(P^T R P) is enough for trace(M^-1 C), even after a mass edit.
        // One temporary mapped row keeps preparation storage linear in nodes.
        let mut diagonal = vec![Dd::ZERO; nodes.len()];
        let mut mapped = vec![Dd::ZERO; nodes.len()];
        for (factor, pivot) in self.rows.iter().zip(&self.pivots) {
            mapped.fill(Dd::ZERO);
            for (port, projection) in projections.iter().enumerate() {
                for &(node, weight) in projection {
                    mapped[node] = mapped[node] + factor[port] * weight;
                }
            }
            for (d, weight) in diagonal.iter_mut().zip(&mapped) {
                *d = *d + *pivot * *weight * *weight;
                if !finite(*d) {
                    return Err(invalid("mapped resistance exceeds finite arithmetic"));
                }
            }
        }
        let load = BoundResistance {
            channels: vec![Dd::ZERO; self.rank()],
            nodal_forces: vec![Dd::ZERO; nodes.len()],
            source: self,
            nodes,
            projections,
            diagonal,
            dissipated_j: Dd::ZERO,
        };
        if !load
            .rate_bound(|b, n| mass(b, n).expect("validated nodal mass"))
            .is_finite()
        {
            return Err(invalid(
                "mass-normalized resistance exceeds finite arithmetic",
            ));
        }
        Ok(load)
    }
}

/// A prepared reciprocal load. All hot-loop storage is bounded or allocated by bind().
#[derive(Clone)]
pub(crate) struct BoundResistance {
    pub(crate) source: ResistanceSample,
    nodes: Vec<(usize, usize)>,
    projections: Vec<Vec<(usize, Dd)>>,
    diagonal: Vec<Dd>,
    channels: Vec<Dd>,
    nodal_forces: Vec<Dd>,
    dissipated_j: Dd,
}

impl BoundResistance {
    pub(crate) fn bodies(&self) -> impl Iterator<Item = usize> + '_ {
        self.nodes.iter().map(|n| n.0)
    }

    /// trace(M^-1 P^T R P), a conservative upper bound for lambda_max.
    /// The physical diagonal is cached, so updated masses cost O(nodes).
    pub(crate) fn rate_bound(&self, mass: impl Fn(usize, usize) -> f64) -> f64 {
        let mut trace = Dd::ZERO;
        for (&(b, n), &diagonal) in self.nodes.iter().zip(&self.diagonal) {
            let m = mass(b, n);
            if !m.is_finite() || m <= 0.0 {
                return f64::INFINITY;
            }
            trace = trace + diagonal / Dd::from_f64(m);
        }
        trace.to_f64()
    }

    fn port_velocity(&self, velocity: impl Fn(usize, usize) -> f64) -> [Dd; MAX_RESISTANCE_PORTS] {
        let mut ports = [Dd::ZERO; MAX_RESISTANCE_PORTS];
        for (flow, projection) in ports.iter_mut().zip(&self.projections) {
            for &(node, weight) in projection {
                let (b, n) = self.nodes[node];
                *flow = *flow + weight * Dd::from_f64(velocity(b, n));
            }
        }
        ports
    }

    pub(crate) fn sample_backward_velocity(&mut self, velocity: impl Fn(usize, usize) -> f64) {
        let ports = self.port_velocity(velocity);
        for (channel, row) in self.channels.iter_mut().zip(&self.source.rows) {
            let mut sum = Dd::ZERO;
            for (&weight, &flow) in row.iter().zip(&ports) {
                sum = sum + weight * flow;
            }
            *channel = sum;
        }
    }

    /// Sparse reciprocal scatter, O(footprint_terms + rank * ports).
    /// DD is retained across channel, port and shared-node interference.
    pub(crate) fn scatter_forces(&mut self, mut add_force: impl FnMut(usize, usize, f64)) {
        let mut port_forces = [Dd::ZERO; MAX_RESISTANCE_PORTS];
        for ((row, pivot), channel) in self
            .source
            .rows
            .iter()
            .zip(&self.source.pivots)
            .zip(&self.channels)
        {
            let force = -*pivot * *channel;
            for (sum, &weight) in port_forces.iter_mut().zip(row) {
                *sum = *sum + weight * force;
            }
        }
        self.nodal_forces.fill(Dd::ZERO);
        for (projection, &force) in self.projections.iter().zip(&port_forces) {
            for &(node, weight) in projection {
                self.nodal_forces[node] = self.nodal_forces[node] + weight * force;
            }
        }
        for (&(b, n), &force) in self.nodes.iter().zip(&self.nodal_forces) {
            add_force(b, n, force.to_f64());
        }
    }

    fn power(&self, velocity: impl Fn(usize, usize) -> f64) -> Dd {
        let ports = self.port_velocity(velocity);
        let mut power = Dd::ZERO;
        for (row, pivot) in self.source.rows.iter().zip(&self.source.pivots) {
            let mut flow = Dd::ZERO;
            for (&weight, &port) in row.iter().zip(&ports) {
                flow = flow + weight * port;
            }
            power = power + *pivot * flow * flow;
        }
        power
    }

    /// Nonnegative centered-velocity dissipation of the declared resistance.
    pub(crate) fn record_centered_loss(&mut self, dt: f64, velocity: impl Fn(usize, usize) -> f64) {
        self.dissipated_j = self.dissipated_j + Dd::from_f64(dt) * self.power(velocity);
    }

    pub(crate) fn dissipated_j(&self) -> f64 {
        self.dissipated_j.to_f64()
    }

    /// Backward-velocity force has a finite-step storage correction:
    /// E_modified = E_mechanical - dt/4 v_half^T C v_half.
    /// Its decrease equals dt v_center^T C v_center in a free linear system.
    /// This is a discretization correction, not energy stored in physical air.
    pub(crate) fn energy_correction(&self, dt: f64, velocity: impl Fn(usize, usize) -> f64) -> f64 {
        (-Dd::from_f64(0.25 * dt) * self.power(velocity)).to_f64()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(r: Vec<f64>, n: usize) -> ResistanceSample {
        ResistanceSample::new(n, 200.0, 1.2, 343.0, r).unwrap()
    }

    fn separate(n: usize) -> Vec<RadiationFootprint> {
        (0..n).map(|i| vec![(i, 0, 1.0)]).collect()
    }

    #[test]
    fn full_signed_psd_load_has_exact_bright_and_dark_collective_modes() {
        let source = sample(vec![2.0, -2.0, -2.0, 2.0], 2);
        assert_eq!(source.rank(), 1);
        let mut load = source.bind(&separate(2), |_, _| Some(1.0)).unwrap();
        let mut forces = [0.0; 2];
        load.sample_backward_velocity(|_, _| 3.0);
        load.scatter_forces(|b, _, f| forces[b] = f);
        assert_eq!(forces, [0.0; 2]);
        load.sample_backward_velocity(|b, _| if b == 0 { 3.0 } else { -3.0 });
        load.scatter_forces(|b, _, f| forces[b] = f);
        assert_eq!(forces, [-12.0, 12.0]);
        assert_eq!(
            load.power(|b, _| if b == 0 { 3.0 } else { -3.0 }).to_f64(),
            72.0
        );
        assert_eq!(load.rate_bound(|_, _| 1.0), 4.0);
    }

    #[test]
    fn weak_channel_survives_coherent_strong_weak_opposite_sum() {
        let source = sample(vec![1.0; 9], 3);
        assert_eq!(source.rank(), 1);
        let mut load = source.bind(&separate(3), |_, _| Some(1.0)).unwrap();
        let weak = 2.0_f64.powi(-60);
        let velocities = [1.0, weak, -1.0];
        load.sample_backward_velocity(|b, _| velocities[b]);
        let mut forces = [0.0; 3];
        load.scatter_forces(|b, _, f| forces[b] = f);
        assert_eq!(forces, [-weak; 3]);
        assert_eq!(load.power(|b, _| velocities[b]).to_f64(), weak * weak);
    }

    #[test]
    fn projection_and_force_are_work_conjugate_and_admission_counts_cross_terms() {
        let source = sample(vec![4.0, -1.0, -1.0, 2.0], 2);
        let p = vec![
            vec![(0, 0, 0.5), (1, 0, 0.5)],
            vec![(0, 0, 1.0), (1, 0, -1.0)],
        ];
        let mut load = source
            .bind(&p, |b, _| Some(if b == 0 { 2.0 } else { 3.0 }))
            .unwrap();
        let v = [0.3, -0.7];
        let mut f = [0.0; 2];
        load.sample_backward_velocity(|b, _| v[b]);
        load.scatter_forces(|b, _, x| f[b] = x);
        let q = [0.5 * (v[0] + v[1]), v[0] - v[1]];
        let rq = [4.0 * q[0] - q[1], -q[0] + 2.0 * q[1]];
        assert!((f[0] + 0.5 * rq[0] + rq[1]).abs() < 1e-15);
        assert!((f[1] + 0.5 * rq[0] - rq[1]).abs() < 1e-15);
        assert!((v[0] * f[0] + v[1] * f[1] + load.power(|b, _| v[b]).to_f64()).abs() < 1e-15);
        // P^T R P = [[2,-1],[-1,4]], trace(M^-1 C)=2/2+4/3.
        assert!((load.rate_bound(|b, _| if b == 0 { 2.0 } else { 3.0 }) - 7.0 / 3.0).abs() < 1e-15);
    }

    #[test]
    fn explicit_mutual_force_balances_modified_energy_and_centered_dissipation() {
        let mut load = sample(vec![3.0, -1.0, -1.0, 2.0], 2)
            .bind(&separate(2), |b, _| Some([2.0, 0.7][b]))
            .unwrap();
        let dt = 0.1;
        let mass = [2.0, 0.7];
        assert!(0.5 * dt * load.rate_bound(|b, _| mass[b]) < 1.0);
        let mut v = [0.8, -0.3];
        let energy = |v: [f64; 2]| 0.5 * (mass[0] * v[0] * v[0] + mass[1] * v[1] * v[1]);
        let initial = energy(v) + load.energy_correction(dt, |b, _| v[b]);
        for _ in 0..120 {
            let mut f = [0.0; 2];
            load.sample_backward_velocity(|b, _| v[b]);
            load.scatter_forces(|b, _, x| f[b] = x);
            let next = [v[0] + dt * f[0] / mass[0], v[1] + dt * f[1] / mass[1]];
            load.record_centered_loss(dt, |b, _| 0.5 * (v[b] + next[b]));
            let balance =
                energy(next) + load.energy_correction(dt, |b, _| next[b]) + load.dissipated_j();
            assert!(
                (balance - initial).abs() < 2e-15,
                "{balance} versus {initial}"
            );
            v = next;
        }
        assert!(load.dissipated_j() > 0.0);
    }

    #[test]
    fn active_nonreciprocal_unbounded_and_invalid_footprints_refuse() {
        for r in [
            vec![1.0, 2.0, 2.0, 1.0],
            vec![0.0, 1.0, 1.0, 0.0],
            vec![1.0, 0.1, 0.2, 1.0],
            vec![f64::NAN; 4],
        ] {
            assert!(ResistanceSample::new(2, 200.0, 1.2, 343.0, r).is_err());
        }
        assert!(ResistanceSample::new(129, 200.0, 1.2, 343.0, vec![]).is_err());
        assert!(
            sample(vec![1.0], 1)
                .bind(&[vec![(0, 0, 1.0)]], |_, _| None)
                .is_err()
        );
        assert!(
            sample(vec![1.0], 1)
                .bind(&[vec![(0, 0, f64::INFINITY)]], |_, _| Some(1.0))
                .is_err()
        );
        assert!(
            sample(vec![1.0], 1)
                .bind(&[vec![(0, 0, 1.0)]], |_, _| Some(0.0))
                .is_err()
        );
        let zero = sample(vec![0.0; 4], 2);
        assert_eq!(zero.rank(), 0);
        assert_eq!(
            zero.bind(&separate(2), |_, _| Some(1.0))
                .unwrap()
                .rate_bound(|_, _| 1.0),
            0.0
        );
    }
}
