//! Experimental Galerkin restriction of the actual fine-grid body.
//!
//! This retains the original fine-grid nonlinear energy/Airy solve and
//! retarded microphone taps. It is a reference reduction, not a fast modal
//! surrogate or a bandwidth certificate. Basis preparation belongs to the
//! existing modal solvers; no eigenvectors or missing frequencies are invented.

use super::{Band, Body, Contact, Law, contact_force, contact_phi, dot};
use std::collections::BTreeMap;

#[derive(Clone)]
pub(crate) struct Reduction {
    // Column-major physical displacement shapes, Phi' M Phi = I.
    pub(super) basis: Vec<Vec<f64>>,
    base: Option<Band>,
    k: f64,
    pg: Vec<f64>,
    pu: Vec<f64>,
    zg: Vec<f64>,
    zu: Vec<f64>,
    increment: Vec<f64>,
    // Inverse of C^-1 + U' B^-1 U, U=[Phi'g,Phi'u].
    correction: [f64; 3],
}

impl Body {
    /// Restrict subsequent motion to an explicit mass-orthonormal basis.
    ///
    /// Install before loading or striking the body: nonzero initial motion is
    /// refused, never silently projected away. Every shape uses this body's
    /// exact node order and mass. All retained coordinates, including supplied
    /// rigid/null shapes, participate in force and contact projection.
    /// Direct edits to public nodal state after installation must themselves
    /// lie in the retained space; ordinary force/strike stepping preserves it.
    ///
    /// This bounded experimental path retains the full nonlinear force solve,
    /// fine-grid clock/CFL, mounts, forces and microphone propagation delays.
    /// It does NOT establish that the supplied basis resolves an audible band,
    /// the omitted contact compliance, or that rendering becomes faster.
    /// Normal contact forces are solved jointly; friction retains the source's
    /// post-normal backward-slip split, without a joint cone-passivity claim.
    /// A contact-budget or convergence failure stops the experimental step.
    /// `Instrument::try_image` refuses external-stepper export of this path.
    pub fn set_galerkin_basis(&mut self, basis: Vec<Vec<f64>>) -> Result<(), String> {
        if self.reduction.is_some()
            || self
                .w
                .iter()
                .chain(&self.wp)
                .chain(&self.wn)
                .any(|&x| x != 0.0)
            || self.psi != 0.0
            || self.psi_next != 0.0
        {
            return Err(
                "Galerkin basis must be installed once, before initial motion or preload".into(),
            );
        }
        let m = basis.len();
        if m == 0
            || m > self.n
            || m > 256
            || self.n.checked_mul(m).is_none_or(|n| n > 4_194_304)
            || self.mass.len() != self.n
            || self.mass.iter().any(|x| !x.is_finite() || *x <= 0.0)
            || basis
                .iter()
                .any(|v| v.len() != self.n || v.iter().any(|x| !x.is_finite()))
        {
            return Err("Galerkin basis requires 1..=min(nodes,256) finite shapes and positive masses; maximum 4194304 entries".into());
        }
        for i in 0..m {
            for j in 0..=i {
                let gram: f64 = (0..self.n)
                    .map(|k| basis[i][k] * self.mass[k] * basis[j][k])
                    .sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                if !gram.is_finite() || (gram - expected).abs() > 1e-8 {
                    return Err(format!(
                        "Galerkin basis is not mass-orthonormal at ({i},{j}): {gram}"
                    ));
                }
            }
        }
        self.reduction = Some(Reduction {
            basis,
            base: None,
            k: 0.0,
            pg: vec![0.0; m],
            pu: vec![0.0; m],
            zg: vec![0.0; m],
            zu: vec![0.0; m],
            increment: vec![0.0; self.n],
            correction: [0.0; 3],
        });
        Ok(())
    }
}

impl Reduction {
    /// Galerkin solve of A wn = rhs, A=D+g g'/4+mu u u'.
    /// The existing fine solver supplies wn=A^-1 rhs, so reconstruct its
    /// actual rhs before restricting. Solve for the increment from current
    /// motion (which is already in span Phi), then add it to that carrier.
    /// This is the same Galerkin equation while avoiding repeated rounding of
    /// the entire displacement before microphone second differences.
    pub(crate) fn predict(
        &mut self,
        k: f64,
        cf: &[[f64; 4]],
        g: &[f64],
        h2: f64,
        mu: f64,
        current: &[f64],
        wn: &mut [f64],
    ) {
        let m = self.basis.len();
        if self.base.is_none() || self.k != k {
            let rows: Vec<BTreeMap<usize, f64>> = self
                .basis
                .iter()
                .map(|a| {
                    self.basis
                        .iter()
                        .enumerate()
                        .map(|(j, b)| (j, (0..wn.len()).map(|i| a[i] * b[i] / cf[i][0]).sum()))
                        .collect()
                })
                .collect();
            // Positive D and independent admitted shapes imply SPD. No
            // regularization or dropped mode is permitted on factor failure.
            self.base = Some(
                Band::factor(&rows)
                    .expect("admitted Galerkin step matrix must be positive definite"),
            );
            self.k = k;
        }
        for i in 0..wn.len() {
            self.increment[i] = wn[i] - current[i];
        }
        let gw = dot(g, &self.increment);
        let uw = h2 * self.increment.iter().sum::<f64>();
        let mut rhs = vec![0.0; m];
        for j in 0..m {
            let phi = &self.basis[j];
            self.pg[j] = dot(phi, g);
            self.pu[j] = h2 * phi.iter().sum::<f64>();
            rhs[j] = (0..wn.len())
                .map(|i| phi[i] * self.increment[i] / cf[i][0])
                .sum::<f64>()
                + 0.25 * self.pg[j] * gw
                + mu * self.pu[j] * uw;
        }
        let base = self.base.as_ref().unwrap();
        self.zg.copy_from_slice(&self.pg);
        self.zu.copy_from_slice(&self.pu);
        base.solve(&mut self.zg);
        base.solve(&mut self.zu);
        let a = 4.0 + dot(&self.pg, &self.zg);
        if mu > 0.0 {
            let b = dot(&self.pg, &self.zu);
            let c = 1.0 / mu + dot(&self.pu, &self.zu);
            let det = a * c - b * b;
            self.correction = [c / det, -b / det, a / det];
        } else {
            self.correction = [1.0 / a, 0.0, 0.0];
        }
        self.solve(&mut rhs);
        wn.copy_from_slice(current);
        self.lift(&rhs, wn);
    }

    fn solve(&self, x: &mut [f64]) {
        self.base
            .as_ref()
            .expect("body prediction precedes contact")
            .solve(x);
        let g = dot(&self.pg, x);
        let u = dot(&self.pu, x);
        let [a, b, c] = self.correction;
        for j in 0..x.len() {
            x[j] -= self.zg[j] * (a * g + b * u) + self.zu[j] * (b * g + c * u);
        }
    }

    fn project(&self, body: usize, port: &[(usize, usize, f64)]) -> Vec<f64> {
        self.basis
            .iter()
            .map(|phi| {
                port.iter()
                    .filter(|t| t.0 == body)
                    .map(|t| phi[t.1] * t.2)
                    .sum()
            })
            .collect()
    }

    fn lift(&self, q: &[f64], out: &mut [f64]) {
        for (phi, &x) in self.basis.iter().zip(q) {
            for (y, &p) in out.iter_mut().zip(phi) {
                *y += p * x;
            }
        }
    }
}

/// Full reciprocal mobility, including cross-node terms on each reduced body.
pub(crate) fn mobility(
    bodies: &[Body],
    x: &[(usize, usize, f64)],
    y: &[(usize, usize, f64)],
    k: f64,
) -> f64 {
    let mut value = 0.0;
    for (b, body) in bodies.iter().enumerate() {
        if !x.iter().any(|t| t.0 == b) || !y.iter().any(|t| t.0 == b) {
            continue;
        }
        if let Some(r) = &body.reduction {
            let px = r.project(b, x);
            let mut py = r.project(b, y);
            r.solve(&mut py);
            value += dot(&px, &py);
        } else {
            for &(_, i, wx) in x.iter().filter(|t| t.0 == b) {
                for &(_, j, wy) in y.iter().filter(|t| t.0 == b && t.1 == i) {
                    value += wx * wy * body.compliance(j, k);
                }
            }
        }
    }
    value
}

pub(crate) fn apply(bodies: &mut [Body], port: &[(usize, usize, f64)], force: f64, k: f64) {
    for (b, body) in bodies.iter_mut().enumerate() {
        if !port.iter().any(|t| t.0 == b) {
            continue;
        }
        if let Some(r) = &body.reduction {
            let mut p = r.project(b, port);
            r.solve(&mut p);
            for v in &mut p {
                *v *= force;
            }
            r.lift(&p, &mut body.wn);
        } else {
            for &(_, i, w) in port.iter().filter(|t| t.0 == b) {
                body.wn[i] += w * force * body.compliance(i, k);
            }
        }
    }
}

/// Simultaneously enforce normal contact laws in the projected mobility.
/// Later contacts can move earlier ports even when their grid nodes differ.
/// No body correction is applied until the final force residual converges.
pub(crate) fn contacts(
    bodies: &mut [Body],
    contacts: &mut [Contact],
    k: f64,
) -> Result<(), String> {
    let active: Vec<_> = contacts
        .iter()
        .enumerate()
        .filter(|(_, c)| c.active)
        .map(|(j, _)| j)
        .collect();
    let n = active.len();
    // A grid2 hi-hat has about 300 rim pairs. Bound the dense matrix to
    // 512² doubles (2 MiB), rather than imposing the unrelated mode budget.
    if n > 512 {
        return Err("experimental Galerkin contact island exceeds 512 active contacts (2 MiB matrix budget)".into());
    }
    let ports: Vec<Vec<_>> = active
        .iter()
        .map(|&j| {
            contacts[j]
                .a
                .iter()
                .copied()
                .chain(contacts[j].b.iter().map(|&(b, i, w)| (b, i, -w)))
                .collect()
        })
        .collect();
    let previous: Vec<f64> = active
        .iter()
        .zip(&ports)
        .map(|(&j, p)| {
            p.iter().map(|&(b, i, w)| w * bodies[b].wp[i]).sum::<f64>() - contacts[j].gap
        })
        .collect();
    let mut d: Vec<f64> = active
        .iter()
        .zip(&ports)
        .map(|(&j, p)| {
            p.iter().map(|&(b, i, w)| w * bodies[b].wn[i]).sum::<f64>() - contacts[j].gap
        })
        .collect();
    let mut h = vec![vec![0.0; n]; n];
    for (b, body) in bodies.iter().enumerate() {
        if let Some(r) = &body.reduction {
            // O(contacts*modes² + contacts²*modes), rather than solving
            // the same projected port again for every contact pair.
            let p: Vec<_> = ports.iter().map(|port| r.project(b, port)).collect();
            let mut response = p.clone();
            for z in &mut response {
                r.solve(z);
            }
            for i in 0..n {
                for j in 0..=i {
                    h[i][j] += dot(&p[i], &response[j]);
                }
            }
        } else {
            for i in 0..n {
                for j in 0..=i {
                    for &(_, node, wi) in ports[i].iter().filter(|t| t.0 == b) {
                        for &(_, _, wj) in ports[j].iter().filter(|t| t.0 == b && t.1 == node) {
                            h[i][j] += wi * wj * body.compliance(node, k);
                        }
                    }
                }
            }
        }
    }
    for i in 0..n {
        for j in 0..=i {
            if !h[i][j].is_finite() {
                return Err("nonfinite Galerkin contact mobility".into());
            }
            h[j][i] = h[i][j];
        }
        if h[i][i] < 0.0 {
            return Err("negative Galerkin contact compliance".into());
        }
    }
    let mut force = vec![0.0; n];
    for _ in 0..512 {
        for i in 0..n {
            let c = &contacts[active[i]];
            let f = contact_force(
                d[i] + h[i][i] * force[i],
                h[i][i],
                force[i],
                previous[i],
                c.law,
                c.fixed_point_exit,
                k,
            );
            let change = f - force[i];
            force[i] = f;
            for j in 0..n {
                d[j] -= h[j][i] * change;
            }
        }
        let converged = (0..n).all(|i| {
            let Law::Power { k: kk, alpha, mu } = contacts[active[i]].law;
            let expected = contact_phi(d[i], previous[i], kk, alpha, mu, k).0;
            force[i].is_finite()
                && expected.is_finite()
                && (force[i] - expected).abs() <= 1e-12 + 1e-10 * force[i].abs().max(expected.abs())
        });
        if converged {
            for i in 0..n {
                if force[i] != 0.0 {
                    apply(bodies, &ports[i], -force[i], k);
                }
                let c = &mut contacts[active[i]];
                c.force = force[i];
                c.peak_force = c.peak_force.max(force[i]);
                c.d_now = d[i];
            }
            return Ok(());
        }
    }
    Err("experimental Galerkin normal contact solve did not converge in 512 sweeps".into())
}

/// The reduced slip Gram matrix can be rank one: solve that physical
/// one-direction Coulomb disk instead of silently disabling friction.
pub(crate) fn friction(dx: f64, dy: f64, xx: f64, yy: f64, xy: f64, cap: f64) -> (f64, f64, bool) {
    if cap == 0.0 {
        return (0.0, 0.0, false);
    }
    if dx == 0.0 && dy == 0.0 {
        return (0.0, 0.0, true);
    }
    let trace = xx + yy;
    let det = xx * yy - xy * xy;
    if trace > 0.0 && det <= 64.0 * f64::EPSILON * trace * trace {
        let (fx, fy) = (-dx / trace, -dy / trace);
        let residual = (dx + xx * fx + xy * fy).hypot(dy + xy * fx + yy * fy);
        if residual <= 1e-12 * dx.hypot(dy) {
            let norm = fx.hypot(fy);
            if norm <= cap {
                return (fx, fy, true);
            }
            return (fx * cap / norm, fy * cap / norm, false);
        }
    }
    // Positive-definite or a sliding solution with a nonzero null-direction
    // slip: lambda>0 makes the same existing disk solver nonsingular.
    super::friction_force(dx, dy, xx, yy, xy, cap)
}
