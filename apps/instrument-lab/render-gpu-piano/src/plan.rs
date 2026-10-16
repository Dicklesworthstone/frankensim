//! The contact-free mechanics step of a grand_piano instrument, reorganised for parallel hardware.
//!
//! One step of `Instrument::step_linear_reference` (damp, predict, finish with zero contact force,
//! commit, damp) is linear in the bank state. Its work splits into items that own their data:
//!
//! * a **group** item per bridge group (`Bank::groups`, strings sharing one bridge row): its string
//!   modes, its strings' damper pads and their contact sites;
//! * a **board** item per loaded-board coordinate.
//!
//! Items meet only through small vectors (per group: bridge sources and reactions; per pad: the board
//! amount a damper chain moved). Every coupling between them is a product with a matrix composed here
//! from the bank's own operators (`BankImage`, `DamperImage`), so each phase is one dot product per item
//! and a step needs 7 grid-wide barriers with dampers engaged, 3 without. In exact arithmetic the step
//! equals the CPU's; the rounding differs. `emulate.rs` runs these phases on the CPU, and the GPU
//! kernels run the same phases.
//!
//! Notation (all from the CPU code): `Gm` the group bridge rows (ng x r), `Sl^-1` the string-loss board
//! Schur inverse, `S^-1` the bank's Schur inverse, `F` the board half-step damping flow, `c_g` the string
//! loss's board reaction per unit bridge speed. A string-loss pass maps the board velocity to
//! `E1 v_b + E2 src + E3 d` where per group `src_g = sum cross_k v_k`, `d_g = sum strain_k beta_k inverse_k v_k`,
//! and gives each group the bridge speed `b_g = R1 v_b + R2 src`.
use crate::grand_piano::linear::{BankImage, dampers::DamperImage};

/// A row-major matrix.
#[derive(Clone, Debug, Default)]
pub struct Mat {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<f64>,
}
impl Mat {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }
    pub fn identity(n: usize) -> Self {
        let mut m = Self::zeros(n, n);
        for i in 0..n {
            m.data[i * n + i] = 1.0;
        }
        m
    }
    fn from(rows: usize, cols: usize, data: Vec<f64>) -> Result<Self, String> {
        if data.len() != rows * cols {
            return Err("matrix image has the wrong size".into());
        }
        Ok(Self { rows, cols, data })
    }
    pub fn row(&self, i: usize) -> &[f64] {
        &self.data[i * self.cols..(i + 1) * self.cols]
    }
    fn mul(&self, b: &Mat) -> Mat {
        assert_eq!(self.cols, b.rows);
        let mut out = Mat::zeros(self.rows, b.cols);
        for i in 0..self.rows {
            let o = &mut out.data[i * b.cols..(i + 1) * b.cols];
            for (k, &a) in self.row(i).iter().enumerate() {
                if a == 0.0 {
                    continue;
                }
                for (x, &y) in o.iter_mut().zip(b.row(k)) {
                    *x += a * y;
                }
            }
        }
        out
    }
    fn transpose(&self) -> Mat {
        let mut out = Mat::zeros(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                out.data[j * self.rows + i] = self.data[i * self.cols + j];
            }
        }
        out
    }
    /// diag(d) * self
    fn scale_rows(&self, d: &[f64]) -> Mat {
        let mut out = self.clone();
        for (i, &s) in d.iter().enumerate() {
            for x in &mut out.data[i * self.cols..(i + 1) * self.cols] {
                *x *= s;
            }
        }
        out
    }
    fn add(&self, b: &Mat, s: f64) -> Mat {
        let mut out = self.clone();
        for (x, y) in out.data.iter_mut().zip(&b.data) {
            *x += s * y;
        }
        out
    }
    fn scaled(&self, s: f64) -> Mat {
        Mat {
            rows: self.rows,
            cols: self.cols,
            data: self.data.iter().map(|x| x * s).collect(),
        }
    }
    pub fn dot_row(&self, i: usize, x: &[f64]) -> f64 {
        self.row(i).iter().zip(x).map(|(a, b)| a * b).sum()
    }
}

/// A damper pad: its string's modes, its quadrature stations and its unit responses per direction.
#[derive(Clone, Debug)]
pub struct Pad {
    pub group: usize,
    pub course: usize,
    pub modes: (usize, usize),
    /// station-major: stations x modes
    pub shape: Vec<f64>,
    pub lift: Vec<f64>,
    pub norm: Vec<f64>,
    pub factor: Vec<f64>,
    /// per direction: w (modes), alpha, loss, speed and change (stations, in pass order)
    pub w: [Vec<f64>; 2],
    pub alpha: [f64; 2],
    pub unit_loss: [f64; 2],
    pub speed: [Vec<f64>; 2],
    pub change: [Vec<f64>; 2],
    /// station Gram matrix s_i . s_j + |b|^2 lift_i lift_j (stations x stations, static)
    pub h: Vec<f64>,
    /// the shapes' own Gram matrix s_i . s_j (static): a change sum_i d_i s_i of the string velocities
    /// moves the station speeds u = S' v by G d
    pub gram: Vec<f64>,
    /// per direction, each station's change in the unit response (station-indexed): w = S uc
    pub uc: [Vec<f64>; 2],
    /// per direction, the chain as matrices (`ChainOp`), for the present pedal
    pub op: [ChainOp; 2],
}

/// One damper chain from a zero bridge speed, exactly, as products: with u_i = s_i . v (station i's
/// speed from the string velocities it starts from), the chain's station changes are c = T u, it leaves
/// v + sum_i c_i s_i, moves the board by mu . u, and has loss u' Lambda u and cross term rho . u with
/// the unit response. Station-indexed; T is triangular in the direction's pass order.
#[derive(Clone, Debug, Default)]
pub struct ChainOp {
    pub t: Vec<f64>,
    pub lambda: Vec<f64>,
    pub rho: Vec<f64>,
    pub mu: Vec<f64>,
    /// G T and G uc: what the chain's changes and the unit response do to the station speeds
    pub gt: Vec<f64>,
    pub guc: Vec<f64>,
}

impl Pad {
    /// The chain of direction `dir` as a `ChainOp` (`factor`, `speed`, `change` must be current).
    ///
    /// Station i's speed is u_i + sum over the stations j before it of H_ij c_j (the string modes and
    /// the own bridge speed the earlier changes left), and c_i = f_i speed_i. So (I - F H~) c = F u with
    /// H~ the part of H before each station, T = (I - F H~)^-1 F, speed = (I + H~ T) u, and the chain's
    /// loss sum(-c_i speed_i - c_i^2 norm_i / 2) and cross term are quadratic and linear forms in u.
    pub fn chain_op(&self, dir: usize) -> ChainOp {
        let ns = self.stations();
        let order: Vec<usize> = if dir == 0 {
            (0..ns).collect()
        } else {
            (0..ns).rev().collect()
        };
        let mut before = vec![false; ns * ns];
        for (t, &i) in order.iter().enumerate() {
            for &j in &order[..t] {
                before[i * ns + j] = true;
            }
        }
        let hm: Vec<f64> = (0..ns * ns)
            .map(|x| if before[x] { self.h[x] } else { 0.0 })
            .collect();
        let mut t = vec![0.0; ns * ns];
        for (k, &i) in order.iter().enumerate() {
            let mut row = vec![0.0; ns];
            row[i] = 1.0;
            for &j in &order[..k] {
                let hij = hm[i * ns + j];
                if hij == 0.0 {
                    continue;
                }
                for x in 0..ns {
                    row[x] += hij * t[j * ns + x];
                }
            }
            for x in 0..ns {
                t[i * ns + x] = self.factor[i] * row[x];
            }
        }
        // speed map (I + H~ T)
        let mut sp = vec![0.0; ns * ns];
        for i in 0..ns {
            sp[i * ns + i] = 1.0;
            for j in 0..ns {
                let hij = hm[i * ns + j];
                if hij == 0.0 {
                    continue;
                }
                for x in 0..ns {
                    sp[i * ns + x] += hij * t[j * ns + x];
                }
            }
        }
        let (mut us, mut uc) = (vec![0.0; ns], vec![0.0; ns]);
        for (k, &i) in order.iter().enumerate() {
            us[i] = self.speed[dir][k];
            uc[i] = self.change[dir][k];
        }
        let mut lambda = vec![0.0; ns * ns];
        let (mut rho, mut mu) = (vec![0.0; ns], vec![0.0; ns]);
        for x in 0..ns {
            for y in 0..ns {
                let mut acc = 0.0;
                for i in 0..ns {
                    acc += -t[i * ns + x] * sp[i * ns + y]
                        - 0.5 * t[i * ns + x] * self.norm[i] * t[i * ns + y];
                }
                lambda[x * ns + y] = acc;
            }
            for i in 0..ns {
                rho[x] += -t[i * ns + x] * us[i]
                    - sp[i * ns + x] * uc[i]
                    - t[i * ns + x] * uc[i] * self.norm[i];
                mu[x] += t[i * ns + x] * self.lift[i];
            }
        }
        for x in 0..ns {
            for y in x + 1..ns {
                let a = 0.5 * (lambda[x * ns + y] + lambda[y * ns + x]);
                lambda[x * ns + y] = a;
                lambda[y * ns + x] = a;
            }
        }
        ChainOp {
            t,
            lambda,
            rho,
            mu,
            gt: Vec::new(),
            guc: Vec::new(),
        }
    }
}
impl Pad {
    pub fn stations(&self) -> usize {
        self.lift.len()
    }
    pub fn len(&self) -> usize {
        self.modes.1 - self.modes.0
    }
}

/// The damper operators for one pedal value and one set of lifted courses.
#[derive(Clone, Debug, Default)]
pub struct Dampers {
    pub active: Vec<bool>,
    pub any: bool,
    /// per direction: beta = P v_b + Q a (np x r, np x np); new board = V v_b + A a (r x r, r x np)
    pub p: [Mat; 2],
    pub q: [Mat; 2],
    pub v: [Mat; 2],
    pub a: [Mat; 2],
}

#[derive(Clone, Debug)]
pub struct Contact {
    pub modes: (usize, usize),
    pub lift: f64,
}

#[derive(Clone, Debug)]
pub struct Program {
    pub n: usize,
    pub r: usize,
    pub ng: usize,
    /// per coordinate (string modes, then board)
    pub qq: Vec<f64>,
    pub qv: Vec<f64>,
    pub vq: Vec<f64>,
    pub vv: Vec<f64>,
    pub bq: Vec<f64>,
    pub bv: Vec<f64>,
    pub omega2: Vec<f64>,
    /// per string mode
    pub a: Vec<f64>,
    pub beta: Vec<f64>,
    pub hammer: Vec<f64>,
    pub inverse: Vec<f64>,
    pub cross: Vec<f64>,
    pub mass: Vec<f64>,
    pub strain: Vec<f64>,
    pub duration: f64,
    /// per group: string-mode range and bridge row
    pub group_modes: Vec<(usize, usize)>,
    pub gm: Mat,
    pub r1: Mat,
    pub r2: Mat,
    pub r3: Mat,
    pub r4: Mat,
    pub r5: Mat,
    pub r6: Mat,
    pub k1: Mat,
    pub k2: Mat,
    pub k3: Mat,
    pub k4: Mat,
    pub k5: Mat,
    pub k6: Mat,
    pub k7: Mat,
    pub e2: Mat,
    pub e3: Mat,
    pub gt: Mat,
    pub contacts: Vec<Contact>,
    pub pads: Vec<Pad>,
    /// per group: its pads (contiguous) and its contacts
    pub group_pads: Vec<(usize, usize)>,
    pub group_contacts: Vec<Vec<usize>>,
    pub gram: Mat,
    pub dampers: Dampers,
}

impl Program {
    pub fn new(bank: &BankImage, dampers: Option<&DamperImage>) -> Result<Self, String> {
        let (n, r) = (bank.n, bank.r);
        let ng = bank.groups.len();
        if bank.transition.len() != n + r
            || bank.modes.len() != n
            || bank.omega2.len() != n + r
            || bank.q.len() != n + r
            || bank.v.len() != n + r
            || bank.schur_inverse.len() != r * r
        {
            return Err("bank image dimensions disagree".into());
        }
        let col = |i: usize| bank.transition.iter().map(|t| t[i]).collect::<Vec<f64>>();
        let (qq, qv, vq, vv, bq, bv) = (col(0), col(1), col(2), col(3), col(4), col(5));
        // groups: contiguous strings with contiguous modes, one bridge row
        let mut group_modes = Vec::with_capacity(ng);
        let mut gm = Mat::zeros(ng, r);
        let mut string_group = vec![usize::MAX; bank.strings.len()];
        let mut next = 0;
        for (g, range) in bank.groups.iter().enumerate() {
            let first = &bank.strings[range.start];
            if first.bridge.len() != r {
                return Err("bridge row length".into());
            }
            gm.data[g * r..(g + 1) * r].copy_from_slice(&first.bridge);
            let start = first.modes.start;
            let mut end = start;
            for si in range.clone() {
                let s = &bank.strings[si];
                if s.modes.start != end {
                    return Err("a bridge group's string modes are not contiguous".into());
                }
                if s.bridge != first.bridge {
                    return Err("a bridge group's strings have different bridge rows".into());
                }
                end = s.modes.end;
                string_group[si] = g;
            }
            if start != next {
                return Err("bridge groups do not tile the string modes in order".into());
            }
            next = end;
            group_modes.push((start, end));
        }
        if next != n {
            return Err("bridge groups do not cover every string mode".into());
        }
        let a: Vec<f64> = bank.modes.iter().map(|m| m.a).collect();
        let beta: Vec<f64> = bank.modes.iter().map(|m| m.beta).collect();
        let hammer: Vec<f64> = bank.modes.iter().map(|m| m.hammer_shape).collect();
        // string loss, or the identity pass when none is configured
        let (inverse, cross, mass, strain, drag, duration, sl_inv) = match &bank.string_loss {
            Some(l) => (
                l.inverse.clone(),
                l.cross_response.clone(),
                l.mass_rate.clone(),
                l.strain_rate.clone(),
                l.endpoint_drag.clone(),
                l.duration_s,
                Mat::from(r, r, l.schur_inverse.clone())?,
            ),
            None => (
                vec![1.0; n],
                vec![0.0; n],
                vec![0.0; n],
                vec![0.0; n],
                vec![0.0; bank.strings.len()],
                0.0,
                Mat::identity(r),
            ),
        };
        let mut c = vec![0.0; ng];
        for (g, range) in bank.groups.iter().enumerate() {
            for si in range.clone() {
                c[g] += drag[si];
                for k in bank.strings[si].modes.clone() {
                    c[g] -= strain[k] * beta[k] * (cross[k] - beta[k]);
                }
            }
        }
        let flow = match &bank.board_half_damping {
            Some(f) => Mat::from(r, r, f.clone())?,
            None => Mat::identity(r),
        };
        let s_inv = Mat::from(r, r, bank.schur_inverse.clone())?;
        let gt = gm.transpose();
        let r1 = gm.mul(&sl_inv);
        let r2 = r1.mul(&gt);
        let gtc = gt.mul(&Mat::identity(ng).scale_rows(&c));
        let e1 = Mat::identity(r).add(&gtc.mul(&r1), -duration);
        let e2 = gtc.mul(&r2).scaled(-duration);
        let e3 = gt.scaled(duration);
        let k1 = flow.mul(&e1);
        let k2 = flow.mul(&e2);
        let k3 = flow.mul(&e3);
        let inv_bq: Vec<f64> = bq[n..].iter().map(|b| 1.0 / b).collect();
        let k4 = s_inv.transpose().scale_rows(&inv_bq).transpose();
        let k5 = s_inv.mul(&gt);
        let r3 = gm.mul(&k4);
        let r4 = gm.mul(&k5);
        let bv_gt = gt.scale_rows(&bv[n..]);
        let k6 = e1.mul(&flow);
        let k7 = k6.mul(&bv_gt);
        let r5 = r1.mul(&flow);
        let r6 = r5.mul(&bv_gt);
        let mut contacts = Vec::with_capacity(bank.contact_strings.len());
        let mut group_contacts = vec![Vec::new(); ng];
        for (ci, &si) in bank.contact_strings.iter().enumerate() {
            let s = &bank.strings[si];
            let g = string_group[si];
            contacts.push(Contact {
                modes: (s.modes.start, s.modes.end),
                lift: s.hammer_lift,
            });
            group_contacts[g].push(ci);
        }
        let mut program = Program {
            n,
            r,
            ng,
            qq,
            qv,
            vq,
            vv,
            bq,
            bv,
            omega2: bank.omega2.clone(),
            a,
            beta,
            hammer,
            inverse,
            cross,
            mass,
            strain,
            duration,
            group_modes,
            gm,
            r1,
            r2,
            r3,
            r4,
            r5,
            r6,
            k1,
            k2,
            k3,
            k4,
            k5,
            k6,
            k7,
            e2,
            e3,
            gt,
            contacts,
            pads: Vec::new(),
            group_pads: vec![(0, 0); ng],
            group_contacts,
            gram: Mat::zeros(0, 0),
            dampers: Dampers::default(),
        };
        if let Some(d) = dampers {
            program.set_pads(bank, d, &string_group)?;
        }
        Ok(program)
    }

    fn set_pads(
        &mut self,
        bank: &BankImage,
        d: &DamperImage,
        string_group: &[usize],
    ) -> Result<(), String> {
        if d.board_start != self.n {
            return Err("damper image addresses another bank".into());
        }
        let np = d.pads.len();
        self.pads.clear();
        let mut last_group = 0;
        for (p, pad) in d.pads.iter().enumerate() {
            let si = bank.modes[pad.modes.start].string;
            let g = string_group[si];
            if g < last_group {
                return Err("damper pads are not in bridge-group order".into());
            }
            last_group = g;
            if bank.strings[si].modes != pad.modes || pad.bridge != bank.strings[si].bridge {
                return Err("a damper pad does not match its string".into());
            }
            if self.group_pads[g].1 == 0 {
                self.group_pads[g] = (p, p + 1);
            } else {
                self.group_pads[g].1 = p + 1;
            }
            let m = pad.modes.len();
            let mut shape = Vec::with_capacity(pad.shape.len() * m);
            for s in &pad.shape {
                if s.len() != m {
                    return Err("damper station shape length".into());
                }
                shape.extend_from_slice(s);
            }
            let ns = pad.lift.len();
            let mut gram = vec![0.0; ns * ns];
            for i in 0..ns {
                for j in 0..ns {
                    gram[i * ns + j] = shape[i * m..(i + 1) * m]
                        .iter()
                        .zip(&shape[j * m..(j + 1) * m])
                        .map(|(a, b)| a * b)
                        .sum::<f64>();
                }
            }
            let h = (0..ns * ns)
                .map(|x| gram[x] + pad.bridge_norm2 * pad.lift[x / ns] * pad.lift[x % ns])
                .collect();
            self.pads.push(Pad {
                h,
                gram,
                uc: [vec![0.0; ns], vec![0.0; ns]],
                op: [ChainOp::default(), ChainOp::default()],
                group: g,
                course: pad.course,
                modes: (pad.modes.start, pad.modes.end),
                shape,
                lift: pad.lift.clone(),
                norm: pad.norm.clone(),
                factor: pad.factor.clone(),
                w: [pad.unit[0].w.clone(), pad.unit[1].w.clone()],
                alpha: [pad.unit[0].alpha, pad.unit[1].alpha],
                unit_loss: [pad.unit[0].loss, pad.unit[1].loss],
                speed: [pad.unit[0].speed.clone(), pad.unit[1].speed.clone()],
                change: [pad.unit[0].change.clone(), pad.unit[1].change.clone()],
            });
        }
        // groups without pads keep an empty range at their position
        let mut cursor = 0;
        for g in 0..self.ng {
            if self.group_pads[g].1 == 0 {
                self.group_pads[g] = (cursor, cursor);
            } else {
                cursor = self.group_pads[g].1;
            }
        }
        self.gram = Mat::from(np, np, d.gram.clone())?;
        Ok(())
    }

    /// Refresh the pads' factors and unit responses (they depend on the pedal) and compose the damper
    /// operators for this set of lifted courses.
    pub fn set_dampers(
        &mut self,
        d: Option<&DamperImage>,
        lifted: &[bool],
        sustain: f64,
    ) -> Result<(), String> {
        let (n, r, np) = (self.n, self.r, self.pads.len());
        let Some(d) = d else {
            self.dampers = Dampers {
                active: vec![],
                any: false,
                ..Dampers::default()
            };
            return Ok(());
        };
        if d.pads.len() != np {
            return Err("damper image changed shape".into());
        }
        for (pad, image) in self.pads.iter_mut().zip(&d.pads) {
            pad.factor.clone_from(&image.factor);
            for dir in 0..2 {
                pad.w[dir].clone_from(&image.unit[dir].w);
                pad.alpha[dir] = image.unit[dir].alpha;
                pad.unit_loss[dir] = image.unit[dir].loss;
                pad.speed[dir].clone_from(&image.unit[dir].speed);
                pad.change[dir].clone_from(&image.unit[dir].change);
            }
            pad.op = [pad.chain_op(0), pad.chain_op(1)];
            let ns = pad.stations();
            pad.uc = [
                pad.change[0].clone(),
                (0..ns).map(|i| pad.change[1][ns - 1 - i]).collect(),
            ];
            for dir in 0..2 {
                let (gram, op) = (&pad.gram, &mut pad.op[dir]);
                op.gt = (0..ns * ns)
                    .map(|x| {
                        (0..ns)
                            .map(|y| gram[(x / ns) * ns + y] * op.t[y * ns + x % ns])
                            .sum()
                    })
                    .collect();
                op.guc = (0..ns)
                    .map(|i| (0..ns).map(|y| gram[i * ns + y] * pad.uc[dir][y]).sum())
                    .collect();
            }
        }
        let active: Vec<bool> = self
            .pads
            .iter()
            .map(|p| sustain < 1.0 && !lifted.get(p.course).copied().unwrap_or(false))
            .collect();
        let any = active.iter().any(|&x| x);
        let mut bridge = Mat::zeros(np, r);
        for (p, pad) in self.pads.iter().enumerate() {
            if active[p] {
                bridge.data[p * r..(p + 1) * r].copy_from_slice(self.gm.row(pad.group));
            }
        }
        let _ = n;
        let bt = bridge.transpose();
        let mut out = Dampers {
            active: active.clone(),
            any,
            ..Dampers::default()
        };
        for dir in 0..2 {
            let order: Vec<usize> = if dir == 0 {
                (0..np).filter(|&p| active[p]).collect()
            } else {
                (0..np).rev().filter(|&p| active[p]).collect()
            };
            let alpha: Vec<f64> = self.pads.iter().map(|p| p.alpha[dir]).collect();
            // forward substitution in pass order: X_p = B_p + sum_{q before p} gram_pq alpha_q X_q, and the
            // same for Y with right-hand side T_p (gram_pq for q before p)
            let mut x = Mat::zeros(np, r);
            let mut y = Mat::zeros(np, np);
            for (k, &p) in order.iter().enumerate() {
                let mut xp = bridge.row(p).to_vec();
                let mut yp = vec![0.0; np];
                for &q in &order[..k] {
                    let t = self.gram.data[p * np + q];
                    yp[q] += t;
                    let s = t * alpha[q];
                    if s == 0.0 {
                        continue;
                    }
                    for (o, v) in xp.iter_mut().zip(x.row(q)) {
                        *o += s * v;
                    }
                    for (o, v) in yp.iter_mut().zip(y.row(q)) {
                        *o += s * v;
                    }
                }
                x.data[p * r..(p + 1) * r].copy_from_slice(&xp);
                y.data[p * np..(p + 1) * np].copy_from_slice(&yp);
            }
            let bta = bt.mul(&Mat::identity(np).scale_rows(&alpha));
            out.v[dir] = Mat::identity(r).add(&bta.mul(&x), 1.0);
            out.a[dir] = bt.add(&bta.mul(&y), 1.0);
            out.p[dir] = x;
            out.q[dir] = y;
        }
        self.dampers = out;
        Ok(())
    }

    pub fn np(&self) -> usize {
        self.pads.len()
    }
}
