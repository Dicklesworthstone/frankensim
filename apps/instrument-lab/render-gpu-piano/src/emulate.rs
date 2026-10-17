//! The phases of `plan.rs` run on one CPU core, in the order the GPU kernels run them. This is the
//! oracle the kernels are checked against (`--check`) and the fallback when no GPU is usable.
//!
//! The dampers run in station space. A chain sees its string only through the station speeds u = S' v
//! and adds S c to it, and a unit response is w = S uc; so through both damper calls of a step and the
//! next one a pad carries just u and the accumulated station changes Delta (u moves by G d with each
//! change d), and its string velocities take S Delta once, before the string loss reads them.
use crate::plan::{Mat, Program};

/// The per-mode string-loss update of one group at bridge speed `b`; returns its kinetic-energy decrease.
fn string_loss(p: &Program, range: (usize, usize), b: f64, v: &mut [f64]) -> f64 {
    let mut loss = 0.0;
    for (offset, velocity) in v[range.0..range.1].iter_mut().enumerate() {
        let k = range.0 + offset;
        let old = *velocity;
        let mid = p.inverse[k] * old + p.cross[k] * b;
        let rel = mid - p.beta[k] * b;
        *velocity = old - p.duration * (p.mass[k] * mid + p.strain[k] * rel);
        loss += 0.5 * (old - *velocity) * (old + *velocity);
    }
    loss
}

fn dot(m: &Mat, i: usize, x: &[f64]) -> f64 {
    m.dot_row(i, x)
}

pub struct Emulator {
    q: Vec<f64>,
    v: Vec<f64>,
    dvb: Vec<f64>,
    fq: Vec<f64>,
    fv: Vec<f64>,
    vbt: Vec<f64>,
    qbo: Vec<f64>,
    src: Vec<f64>,
    dd: Vec<f64>,
    rea: Vec<f64>,
    b0: Vec<f64>,
    rea2: Vec<f64>,
    src2: Vec<f64>,
    dd2: Vec<f64>,
    a0: Vec<f64>,
    a1: Vec<f64>,
    l: Vec<f64>,
    /// per pad: station speeds, accumulated station changes, the present chain's changes
    u: Vec<Vec<f64>>,
    delta: Vec<Vec<f64>>,
    c: Vec<Vec<f64>>,
}

impl Emulator {
    pub fn new(p: &Program) -> Self {
        let (n, r, ng, np) = (p.n, p.r, p.ng, p.np());
        let st = |pad: &crate::plan::Pad| vec![0.0; pad.stations()];
        Self {
            q: vec![0.0; n + r],
            v: vec![0.0; n + r],
            dvb: vec![0.0; r],
            fq: vec![0.0; n + r],
            fv: vec![0.0; n + r],
            vbt: vec![0.0; r],
            qbo: vec![0.0; r],
            src: vec![0.0; ng],
            dd: vec![0.0; ng],
            rea: vec![0.0; ng],
            b0: vec![0.0; ng],
            rea2: vec![0.0; ng],
            src2: vec![0.0; ng],
            dd2: vec![0.0; ng],
            a0: vec![0.0; np],
            a1: vec![0.0; np],
            l: vec![0.0; np],
            u: p.pads.iter().map(st).collect(),
            delta: p.pads.iter().map(st).collect(),
            c: p.pads.iter().map(st).collect(),
        }
    }

    /// Every engaged pad's station speeds from its string velocities; no change accumulated yet.
    fn stations(&mut self, p: &Program) {
        for (i, pad) in p.pads.iter().enumerate() {
            if !p.dampers.active[i] {
                continue;
            }
            let m = pad.len();
            let v = &self.v[pad.modes.0..pad.modes.1];
            for s in 0..pad.stations() {
                self.u[i][s] = pad.shape[s * m..(s + 1) * m]
                    .iter()
                    .zip(v)
                    .map(|(g, x)| g * x)
                    .sum();
                self.delta[i][s] = 0.0;
            }
        }
    }

    /// The accumulated station changes into the string velocities.
    fn apply(&mut self, p: &Program) {
        for (i, pad) in p.pads.iter().enumerate() {
            if !p.dampers.active[i] {
                continue;
            }
            let m = pad.len();
            let v = &mut self.v[pad.modes.0..pad.modes.1];
            for (k, x) in v.iter_mut().enumerate() {
                *x += (0..pad.stations())
                    .map(|s| pad.shape[s * m + k] * self.delta[i][s])
                    .sum::<f64>();
            }
        }
    }

    /// The chains of one direction (plan::ChainOp) from the pads' station speeds; returns their loss.
    fn chains(&mut self, p: &Program, dir: usize) -> f64 {
        let mut loss = 0.0;
        for (i, pad) in p.pads.iter().enumerate() {
            let a = if dir == 0 {
                &mut self.a0[i]
            } else {
                &mut self.a1[i]
            };
            if !p.dampers.active[i] {
                *a = 0.0;
                continue;
            }
            let (ns, op, u) = (pad.stations(), &pad.op[dir], &self.u[i]);
            let dotu = |row: &[f64]| row.iter().zip(u).map(|(a, b)| a * b).sum::<f64>();
            for s in 0..ns {
                self.c[i][s] = dotu(&op.t[s * ns..(s + 1) * ns]);
            }
            *a = dotu(&op.mu);
            self.l[i] = dotu(&op.rho);
            loss += (0..ns)
                .map(|s| u[s] * dotu(&op.lambda[s * ns..(s + 1) * ns]))
                .sum::<f64>();
        }
        loss
    }

    /// The superposition of direction `dir` (its chains already ran): every pad's bridge speed, its
    /// station changes c + beta uc (accumulated, and moving its station speeds) and loss; the board.
    fn superpose(&mut self, p: &Program, dir: usize) -> f64 {
        let (n, r) = (p.n, p.r);
        let d = &p.dampers;
        let mut loss = 0.0;
        let a = if dir == 0 { &self.a0 } else { &self.a1 };
        let board_in: Vec<f64> = if dir == 0 {
            self.v[n..].to_vec()
        } else {
            self.vbt.clone()
        };
        for (i, pad) in p.pads.iter().enumerate() {
            if !d.active[i] {
                continue;
            }
            let beta = dot(&d.p[dir], i, &board_in) + dot(&d.q[dir], i, a);
            let ns = pad.stations();
            let change: Vec<f64> = (0..ns)
                .map(|s| self.c[i][s] + beta * pad.uc[dir][s])
                .collect();
            for s in 0..ns {
                self.delta[i][s] += change[s];
                self.u[i][s] += pad.gram[s * ns..(s + 1) * ns]
                    .iter()
                    .zip(&change)
                    .map(|(g, x)| g * x)
                    .sum::<f64>();
            }
            loss += beta * self.l[i] + beta * beta * pad.unit_loss[dir];
        }
        for j in 0..r {
            let out = dot(&d.v[dir], j, &board_in) + dot(&d.a[dir], j, a);
            if dir == 0 {
                self.vbt[j] = out;
            } else {
                self.v[n + j] = out;
            }
        }
        loss
    }

    /// Run `steps` contact-free steps from (q, v). Writes each step's board velocity (r values) into
    /// `trace` and each step's contact positions into `positions`; returns (modal loss, damper loss).
    pub fn run(
        &mut self,
        p: &Program,
        q: &mut [f64],
        v: &mut [f64],
        steps: usize,
        trace: &mut [f64],
        positions: &mut [f64],
    ) -> (f64, f64) {
        let (n, r, ng, nc) = (p.n, p.r, p.ng, p.contacts.len());
        self.q.copy_from_slice(q);
        self.v.copy_from_slice(v);
        let damp = p.dampers.any;
        let (mut modal, mut damper) = (0.0, 0.0);
        if damp {
            self.stations(p);
            damper += self.chains(p, 0);
        }
        for step in 0..steps {
            if damp {
                damper += self.superpose(p, 0);
                damper += self.chains(p, 1);
                damper += self.superpose(p, 1);
                self.apply(p);
            }
            // SL1: the string-loss sources of the predicted velocity
            for g in 0..ng {
                let (s, e) = p.group_modes[g];
                self.src[g] = (s..e).map(|k| p.cross[k] * self.v[k]).sum();
                self.dd[g] = (s..e)
                    .map(|k| p.strain[k] * p.beta[k] * p.inverse[k] * self.v[k])
                    .sum();
            }
            // Y1: string loss + board flow (predict's damped velocity), free motion, reactions
            for j in 0..r {
                let vb = &self.v[n..];
                self.dvb[j] =
                    dot(&p.k1, j, vb) + dot(&p.k2, j, &self.src) + dot(&p.k3, j, &self.dd);
            }
            for g in 0..ng {
                let b = dot(&p.r1, g, &self.v[n..]) + dot(&p.r2, g, &self.src);
                modal += string_loss(p, p.group_modes[g], b, &mut self.v);
                let (s, e) = p.group_modes[g];
                let b0 = dot(&p.gm, g, &self.q[n..]);
                let mut reaction = 0.0;
                for k in s..e {
                    self.fq[k] = p.qq[k] * self.q[k] + p.qv[k] * self.v[k];
                    self.fv[k] = p.vq[k] * self.q[k] + p.vv[k] * self.v[k];
                    reaction +=
                        0.5 * p.a[k] * (self.q[k] + self.fq[k] + 0.5 * p.bq[k] * p.a[k] * b0);
                }
                self.b0[g] = b0;
                self.rea[g] = reaction;
            }
            for j in 0..r {
                let k = n + j;
                modal += 0.5 * (self.v[k] * self.v[k] - self.dvb[j] * self.dvb[j]);
                self.v[k] = self.dvb[j];
                self.fq[k] = p.qq[k] * self.q[k] + p.qv[k] * self.v[k];
                self.fv[k] = p.vq[k] * self.q[k] + p.vv[k] * self.v[k];
            }
            // P2: board end positions, string end state with the finish ledger, contact positions
            for g in 0..ng {
                let b1 = dot(&p.r3, g, &self.fq[n..]) + dot(&p.r4, g, &self.rea);
                let bbar = 0.5 * (self.b0[g] + b1);
                let (s, e) = p.group_modes[g];
                let (mut rea2, mut src2, mut dd2) = (0.0, 0.0, 0.0);
                for k in s..e {
                    let f = p.a[k] * bbar;
                    let nq = self.fq[k] + p.bq[k] * f;
                    let nv = self.fv[k] + p.bv[k] * f;
                    rea2 += p.a[k] * 0.5 * (self.q[k] + nq);
                    modal += f * (nq - self.q[k])
                        - 0.5
                            * (nv * nv - self.v[k] * self.v[k]
                                + p.omega2[k] * (nq * nq - self.q[k] * self.q[k]));
                    self.q[k] = nq;
                    self.v[k] = nv;
                    src2 += p.cross[k] * nv;
                    dd2 += p.strain[k] * p.beta[k] * p.inverse[k] * nv;
                }
                self.rea2[g] = rea2;
                self.src2[g] = src2;
                self.dd2[g] = dd2;
                for &c in &p.group_contacts[g] {
                    let contact = &p.contacts[c];
                    positions[step * nc + c] = (contact.modes.0..contact.modes.1)
                        .map(|k| p.hammer[k] * self.q[k])
                        .sum::<f64>()
                        + contact.lift * b1;
                }
            }
            let board_end: Vec<f64> = (0..r)
                .map(|j| dot(&p.k4, j, &self.fq[n..]) + dot(&p.k5, j, &self.rea))
                .collect();
            self.qbo[..r].copy_from_slice(&self.q[n..n + r]);
            self.q[n..n + r].copy_from_slice(&board_end);
            // P3: board finish with its ledger, flow and string loss; the groups' string loss; commit
            for j in 0..r {
                let k = n + j;
                let f = dot(&p.gt, j, &self.rea2);
                let nv = self.fv[k] + p.bv[k] * f;
                let (qn, qo) = (self.q[k], self.qbo[j]);
                modal += f * (qn - qo)
                    - 0.5
                        * (nv * nv - self.dvb[j] * self.dvb[j] + p.omega2[k] * (qn * qn - qo * qo));
                let fin = dot(&p.k6, j, &self.fv[n..])
                    + dot(&p.k7, j, &self.rea2)
                    + dot(&p.e2, j, &self.src2)
                    + dot(&p.e3, j, &self.dd2);
                modal += 0.5 * (nv * nv - fin * fin);
                self.vbt[j] = fin;
            }
            for g in 0..ng {
                let b = dot(&p.r5, g, &self.fv[n..])
                    + dot(&p.r6, g, &self.rea2)
                    + dot(&p.r2, g, &self.src2);
                modal += string_loss(p, p.group_modes[g], b, &mut self.v);
            }
            let fin = self.vbt.clone();
            self.v[n..].copy_from_slice(&fin);
            if damp {
                self.stations(p);
                damper += self.chains(p, 0);
                damper += self.superpose(p, 0);
                damper += self.chains(p, 1);
                damper += self.superpose(p, 1);
                // the next step's first damper call continues from these station speeds
                if step + 1 < steps {
                    damper += self.chains(p, 0);
                } else {
                    self.apply(p);
                }
            }
            trace[step * r..(step + 1) * r].copy_from_slice(&self.v[n..]);
        }
        q.copy_from_slice(&self.q);
        v.copy_from_slice(&self.v);
        (modal, damper)
    }
}
