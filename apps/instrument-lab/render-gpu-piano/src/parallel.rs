//! The phase emulator (emulate.rs) on several CPU cores: each thread owns a contiguous run of bridge
//! groups (with their damper pads) and of board rows, the phases are those of the GPU kernel, and the
//! threads meet at a spin barrier between them. Every item's arithmetic is the single-core emulator's,
//! so the state and the traces are bit-identical to it; only the loss sums are added in another order.
//! On Apple silicon the performance cores are used (FS_PIANO_CPU_THREADS overrides).
use crate::accel::Backend;
use crate::plan::{Mat, Program};
use std::sync::atomic::{AtomicUsize, Ordering};

/// A shared f64 array. Within a phase, threads touch disjoint ranges mutably and read only ranges no
/// thread writes in that phase (the phase plan in `step` guarantees it), so the references never alias.
#[derive(Clone, Copy)]
struct Shared(*mut f64, usize);
unsafe impl Send for Shared {}
unsafe impl Sync for Shared {}
impl Shared {
    fn of(v: &mut Vec<f64>) -> Self {
        Self(v.as_mut_ptr(), v.len())
    }
    /// SAFETY: no thread writes [a, b) during this phase.
    unsafe fn r(&self, a: usize, b: usize) -> &[f64] {
        assert!(a <= b && b <= self.1);
        unsafe { std::slice::from_raw_parts(self.0.add(a), b - a) }
    }
    /// SAFETY: only this thread touches [a, b) during this phase.
    #[allow(clippy::mut_from_ref)]
    unsafe fn w(&self, a: usize, b: usize) -> &mut [f64] {
        assert!(a <= b && b <= self.1);
        unsafe { std::slice::from_raw_parts_mut(self.0.add(a), b - a) }
    }
}

struct SpinBarrier {
    count: AtomicUsize,
    generation: AtomicUsize,
    n: usize,
}
impl SpinBarrier {
    fn wait(&self) {
        let g = self.generation.load(Ordering::Acquire);
        if self.count.fetch_add(1, Ordering::AcqRel) + 1 == self.n {
            self.count.store(0, Ordering::Relaxed);
            self.generation.fetch_add(1, Ordering::Release);
        } else {
            while self.generation.load(Ordering::Acquire) == g {
                std::hint::spin_loop();
            }
        }
    }
}

fn dot(m: &Mat, i: usize, x: &[f64]) -> f64 {
    m.dot_row(i, x)
}

struct State {
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
    u: Vec<f64>,
    delta: Vec<f64>,
    c: Vec<f64>,
}

#[derive(Clone, Copy)]
struct Ptrs {
    q: Shared,
    v: Shared,
    dvb: Shared,
    fq: Shared,
    fv: Shared,
    vbt: Shared,
    qbo: Shared,
    src: Shared,
    dd: Shared,
    rea: Shared,
    b0: Shared,
    rea2: Shared,
    src2: Shared,
    dd2: Shared,
    a0: Shared,
    a1: Shared,
    l: Shared,
    u: Shared,
    delta: Shared,
    c: Shared,
    trace: Shared,
    positions: Shared,
}

/// One thread's share: groups [g0, g1) with pads [p0, p1), board rows [j0, j1).
#[derive(Clone, Copy)]
struct Part {
    g0: usize,
    g1: usize,
    p0: usize,
    p1: usize,
    j0: usize,
    j1: usize,
}

pub struct Parallel {
    s: State,
    /// first station of each pad (and the total)
    st: Vec<usize>,
    parts: Vec<Part>,
}

/// Performance cores (macOS) or available parallelism, unless FS_PIANO_CPU_THREADS says otherwise.
pub fn threads() -> usize {
    if let Some(t) = std::env::var("FS_PIANO_CPU_THREADS")
        .ok()
        .and_then(|s| s.parse().ok())
    {
        return t;
    }
    #[cfg(target_os = "macos")]
    if let Ok(out) = std::process::Command::new("/usr/sbin/sysctl")
        .args(["-n", "hw.perflevel0.physicalcpu"])
        .output()
    {
        if let Ok(n) = String::from_utf8_lossy(&out.stdout).trim().parse::<usize>() {
            return n.max(1);
        }
    }
    std::thread::available_parallelism().map_or(1, |n| n.get())
}

impl Parallel {
    pub fn new(p: &Program, threads: usize) -> Self {
        let (n, r, ng, np) = (p.n, p.r, p.ng, p.np());
        let mut st = vec![0usize];
        for pad in &p.pads {
            st.push(st.last().unwrap() + pad.stations());
        }
        let ns = *st.last().unwrap();
        let s = State {
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
            u: vec![0.0; ns],
            delta: vec![0.0; ns],
            c: vec![0.0; ns],
        };
        // groups by weight (modes, damper station-modes, the dense rows), rows evenly
        let t = threads.clamp(1, ng.max(1));
        let weight = |g: usize| {
            let (a, b) = p.group_modes[g];
            let (pb, pe) = p.group_pads[g];
            (b - a) * 4
                + p.pads[pb..pe]
                    .iter()
                    .map(|x| x.stations() * x.len())
                    .sum::<usize>()
                + 3 * r
                + 4 * ng
        };
        let total: usize = (0..ng).map(weight).sum();
        let mut parts = Vec::with_capacity(t);
        let (mut g0, mut acc) = (0usize, 0usize);
        for k in 0..t {
            let target = total * (k + 1) / t;
            let mut g1 = g0;
            while g1 < ng && (acc + weight(g1) <= target || g1 == g0) && ng - g1 > t - k - 1 {
                acc += weight(g1);
                g1 += 1;
            }
            if k + 1 == t {
                while g1 < ng {
                    acc += weight(g1);
                    g1 += 1;
                }
            }
            let (p0, p1) = if g1 > g0 {
                (p.group_pads[g0].0, p.group_pads[g1 - 1].1)
            } else {
                (0, 0)
            };
            parts.push(Part {
                g0,
                g1,
                p0,
                p1,
                j0: r * k / t,
                j1: r * (k + 1) / t,
            });
            g0 = g1;
        }
        Self { s, st, parts }
    }

    pub fn threads(&self) -> usize {
        self.parts.len()
    }
}

/// Station speeds of the part's engaged pads from their string velocities; no change accumulated.
unsafe fn stations(p: &Program, x: &Ptrs, st: &[usize], part: Part) {
    for i in part.p0..part.p1 {
        if !p.dampers.active[i] {
            continue;
        }
        let pad = &p.pads[i];
        let m = pad.len();
        let v = unsafe { x.v.r(pad.modes.0, pad.modes.1) };
        let (u, d) = unsafe { (x.u.w(st[i], st[i + 1]), x.delta.w(st[i], st[i + 1])) };
        for s in 0..pad.stations() {
            u[s] = pad.shape[s * m..(s + 1) * m]
                .iter()
                .zip(v)
                .map(|(g, y)| g * y)
                .sum();
            d[s] = 0.0;
        }
    }
}

unsafe fn apply(p: &Program, x: &Ptrs, st: &[usize], part: Part) {
    for i in part.p0..part.p1 {
        if !p.dampers.active[i] {
            continue;
        }
        let pad = &p.pads[i];
        let m = pad.len();
        let d = unsafe { x.delta.r(st[i], st[i + 1]) };
        let v = unsafe { x.v.w(pad.modes.0, pad.modes.1) };
        for (k, y) in v.iter_mut().enumerate() {
            *y += (0..pad.stations())
                .map(|s| pad.shape[s * m + k] * d[s])
                .sum::<f64>();
        }
    }
}

unsafe fn chains(p: &Program, x: &Ptrs, st: &[usize], part: Part, dir: usize) -> f64 {
    let mut loss = 0.0;
    let a = if dir == 0 { x.a0 } else { x.a1 };
    for i in part.p0..part.p1 {
        let ai = unsafe { a.w(i, i + 1) };
        if !p.dampers.active[i] {
            ai[0] = 0.0;
            continue;
        }
        let pad = &p.pads[i];
        let (ns, op) = (pad.stations(), &pad.op[dir]);
        let u = unsafe { x.u.r(st[i], st[i + 1]) };
        let c = unsafe { x.c.w(st[i], st[i + 1]) };
        let dotu = |row: &[f64]| row.iter().zip(u).map(|(a, b)| a * b).sum::<f64>();
        for s in 0..ns {
            c[s] = dotu(&op.t[s * ns..(s + 1) * ns]);
        }
        ai[0] = dotu(&op.mu);
        unsafe { x.l.w(i, i + 1)[0] = dotu(&op.rho) };
        loss += (0..ns)
            .map(|s| u[s] * dotu(&op.lambda[s * ns..(s + 1) * ns]))
            .sum::<f64>();
    }
    loss
}

/// Direction `dir`'s superposition for the part's pads (from `board_in`, every pad's a), and the
/// part's board rows into `board_out`.
unsafe fn superpose(
    p: &Program,
    x: &Ptrs,
    st: &[usize],
    part: Part,
    dir: usize,
    board_in: (Shared, usize),
    board_out: (Shared, usize),
) -> f64 {
    let d = &p.dampers;
    let (r, np) = (p.r, p.np());
    let mut loss = 0.0;
    let all = if dir == 0 { x.a0 } else { x.a1 };
    let a = unsafe { all.r(0, np) };
    let bin = unsafe { board_in.0.r(board_in.1, board_in.1 + r) };
    for i in part.p0..part.p1 {
        if !d.active[i] {
            continue;
        }
        let pad = &p.pads[i];
        let beta = dot(&d.p[dir], i, bin) + dot(&d.q[dir], i, a);
        let ns = pad.stations();
        let c = unsafe { x.c.r(st[i], st[i + 1]) };
        let change: Vec<f64> = (0..ns).map(|s| c[s] + beta * pad.uc[dir][s]).collect();
        let (u, dl) = unsafe { (x.u.w(st[i], st[i + 1]), x.delta.w(st[i], st[i + 1])) };
        for s in 0..ns {
            dl[s] += change[s];
            u[s] += pad.gram[s * ns..(s + 1) * ns]
                .iter()
                .zip(&change)
                .map(|(g, y)| g * y)
                .sum::<f64>();
        }
        loss += beta * unsafe { x.l.r(i, i + 1)[0] } + beta * beta * pad.unit_loss[dir];
    }
    for j in part.j0..part.j1 {
        let o = board_out.1 + j;
        unsafe { board_out.0.w(o, o + 1)[0] = dot(&d.v[dir], j, bin) + dot(&d.a[dir], j, a) };
    }
    loss
}

unsafe fn sources(p: &Program, x: &Ptrs, part: Part) {
    for g in part.g0..part.g1 {
        let (s, e) = p.group_modes[g];
        let v = unsafe { x.v.r(s, e) };
        unsafe {
            x.src.w(g, g + 1)[0] = (s..e).map(|k| p.cross[k] * v[k - s]).sum();
            x.dd.w(g, g + 1)[0] = (s..e)
                .map(|k| p.strain[k] * p.beta[k] * p.inverse[k] * v[k - s])
                .sum();
        }
    }
}

fn string_loss(p: &Program, range: (usize, usize), b: f64, v: &mut [f64]) -> f64 {
    let mut loss = 0.0;
    for k in range.0..range.1 {
        let i = k - range.0;
        let old = v[i];
        let mid = p.inverse[k] * old + p.cross[k] * b;
        let rel = mid - p.beta[k] * b;
        v[i] = old - p.duration * (p.mass[k] * mid + p.strain[k] * rel);
        loss += 0.5 * (old - v[i]) * (old + v[i]);
    }
    loss
}

impl Parallel {
    fn ptrs(&mut self, trace: &mut [f64], positions: &mut [f64]) -> Ptrs {
        let s = &mut self.s;
        Ptrs {
            q: Shared::of(&mut s.q),
            v: Shared::of(&mut s.v),
            dvb: Shared::of(&mut s.dvb),
            fq: Shared::of(&mut s.fq),
            fv: Shared::of(&mut s.fv),
            vbt: Shared::of(&mut s.vbt),
            qbo: Shared::of(&mut s.qbo),
            src: Shared::of(&mut s.src),
            dd: Shared::of(&mut s.dd),
            rea: Shared::of(&mut s.rea),
            b0: Shared::of(&mut s.b0),
            rea2: Shared::of(&mut s.rea2),
            src2: Shared::of(&mut s.src2),
            dd2: Shared::of(&mut s.dd2),
            a0: Shared::of(&mut s.a0),
            a1: Shared::of(&mut s.a1),
            l: Shared::of(&mut s.l),
            u: Shared::of(&mut s.u),
            delta: Shared::of(&mut s.delta),
            c: Shared::of(&mut s.c),
            trace: Shared(trace.as_mut_ptr(), trace.len()),
            positions: Shared(positions.as_mut_ptr(), positions.len()),
        }
    }

    pub fn run(
        &mut self,
        p: &Program,
        q: &mut [f64],
        v: &mut [f64],
        steps: usize,
        trace: &mut [f64],
        positions: &mut [f64],
    ) -> (f64, f64) {
        self.s.q.copy_from_slice(q);
        self.s.v.copy_from_slice(v);
        let x = self.ptrs(trace, positions);
        let barrier = SpinBarrier {
            count: AtomicUsize::new(0),
            generation: AtomicUsize::new(0),
            n: self.parts.len(),
        };
        let st = &self.st;
        let totals: Vec<(f64, f64)> = std::thread::scope(|scope| {
            let handles: Vec<_> = self.parts[1..]
                .iter()
                .map(|&part| {
                    let barrier = &barrier;
                    scope.spawn(move || unsafe { worker(p, &x, st, part, steps, barrier) })
                })
                .collect();
            let mut all = vec![unsafe { worker(p, &x, st, self.parts[0], steps, &barrier) }];
            all.extend(
                handles
                    .into_iter()
                    .map(|h| h.join().expect("emulator thread")),
            );
            all
        });
        q.copy_from_slice(&self.s.q);
        v.copy_from_slice(&self.s.v);
        totals
            .iter()
            .fold((0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1))
    }
}

/// One thread's phases for `steps` steps (the order of emulate.rs, split at the GPU kernel's barriers).
unsafe fn worker(
    p: &Program,
    x: &Ptrs,
    st: &[usize],
    part: Part,
    steps: usize,
    barrier: &SpinBarrier,
) -> (f64, f64) {
    let (n, r, ng, nc) = (p.n, p.r, p.ng, p.contacts.len());
    let damp = p.dampers.any;
    let (mut modal, mut damper) = (0.0, 0.0);
    unsafe {
        if damp {
            stations(p, x, st, part);
            damper += chains(p, x, st, part, 0);
        } else {
            sources(p, x, part);
        }
        barrier.wait();
        for step in 0..steps {
            if damp {
                damper += superpose(p, x, st, part, 0, (x.v, n), (x.vbt, 0));
                damper += chains(p, x, st, part, 1);
                barrier.wait();
                damper += superpose(p, x, st, part, 1, (x.vbt, 0), (x.v, n));
                apply(p, x, st, part);
                sources(p, x, part);
                barrier.wait();
            }
            // Y1: rows keep the damped board velocity aside (groups still read v_b)
            for j in part.j0..part.j1 {
                let k = n + j;
                let dvb = dot(&p.k1, j, x.v.r(n, n + r))
                    + dot(&p.k2, j, x.src.r(0, ng))
                    + dot(&p.k3, j, x.dd.r(0, ng));
                let vb = x.v.r(k, k + 1)[0];
                modal += 0.5 * (vb * vb - dvb * dvb);
                x.dvb.w(j, j + 1)[0] = dvb;
                let qk = x.q.r(k, k + 1)[0];
                x.fq.w(k, k + 1)[0] = p.qq[k] * qk + p.qv[k] * dvb;
                x.fv.w(k, k + 1)[0] = p.vq[k] * qk + p.vv[k] * dvb;
            }
            for g in part.g0..part.g1 {
                let (s, e) = p.group_modes[g];
                let b = dot(&p.r1, g, x.v.r(n, n + r)) + dot(&p.r2, g, x.src.r(0, ng));
                let vg = x.v.w(s, e);
                modal += string_loss(p, (s, e), b, vg);
                let b0 = dot(&p.gm, g, x.q.r(n, n + r));
                let qg = x.q.r(s, e);
                let (fq, fv) = (x.fq.w(s, e), x.fv.w(s, e));
                let mut reaction = 0.0;
                for k in s..e {
                    let i = k - s;
                    fq[i] = p.qq[k] * qg[i] + p.qv[k] * vg[i];
                    fv[i] = p.vq[k] * qg[i] + p.vv[k] * vg[i];
                    reaction += 0.5 * p.a[k] * (qg[i] + fq[i] + 0.5 * p.bq[k] * p.a[k] * b0);
                }
                x.b0.w(g, g + 1)[0] = b0;
                x.rea.w(g, g + 1)[0] = reaction;
            }
            barrier.wait();
            // P2
            for g in part.g0..part.g1 {
                let (s, e) = p.group_modes[g];
                let b1 = dot(&p.r3, g, x.fq.r(n, n + r)) + dot(&p.r4, g, x.rea.r(0, ng));
                let bbar = 0.5 * (x.b0.r(g, g + 1)[0] + b1);
                let (qg, vg) = (x.q.w(s, e), x.v.w(s, e));
                let (fq, fv) = (x.fq.r(s, e), x.fv.r(s, e));
                let (mut rea2, mut src2, mut dd2) = (0.0, 0.0, 0.0);
                for k in s..e {
                    let i = k - s;
                    let f = p.a[k] * bbar;
                    let nq = fq[i] + p.bq[k] * f;
                    let nv = fv[i] + p.bv[k] * f;
                    rea2 += p.a[k] * 0.5 * (qg[i] + nq);
                    modal += f * (nq - qg[i])
                        - 0.5 * (nv * nv - vg[i] * vg[i] + p.omega2[k] * (nq * nq - qg[i] * qg[i]));
                    qg[i] = nq;
                    vg[i] = nv;
                    src2 += p.cross[k] * nv;
                    dd2 += p.strain[k] * p.beta[k] * p.inverse[k] * nv;
                }
                x.rea2.w(g, g + 1)[0] = rea2;
                x.src2.w(g, g + 1)[0] = src2;
                x.dd2.w(g, g + 1)[0] = dd2;
                for &c in &p.group_contacts[g] {
                    let contact = &p.contacts[c];
                    let pos = (contact.modes.0..contact.modes.1)
                        .map(|k| p.hammer[k] * qg[k - s])
                        .sum::<f64>()
                        + contact.lift * b1;
                    x.positions.w(step * nc + c, step * nc + c + 1)[0] = pos;
                }
            }
            for j in part.j0..part.j1 {
                let k = n + j;
                x.v.w(k, k + 1)[0] = x.dvb.r(j, j + 1)[0];
                let end = dot(&p.k4, j, x.fq.r(n, n + r)) + dot(&p.k5, j, x.rea.r(0, ng));
                x.qbo.w(j, j + 1)[0] = x.q.r(k, k + 1)[0];
                x.q.w(k, k + 1)[0] = end;
            }
            barrier.wait();
            // P3: rows write the board velocity at once (groups read fv_b, rea2, src2, never v_b here)
            for j in part.j0..part.j1 {
                let k = n + j;
                let f = dot(&p.gt, j, x.rea2.r(0, ng));
                let nv = x.fv.r(k, k + 1)[0] + p.bv[k] * f;
                let (qn, qo, dvb) = (
                    x.q.r(k, k + 1)[0],
                    x.qbo.r(j, j + 1)[0],
                    x.dvb.r(j, j + 1)[0],
                );
                modal +=
                    f * (qn - qo) - 0.5 * (nv * nv - dvb * dvb + p.omega2[k] * (qn * qn - qo * qo));
                let fin = dot(&p.k6, j, x.fv.r(n, n + r))
                    + dot(&p.k7, j, x.rea2.r(0, ng))
                    + dot(&p.e2, j, x.src2.r(0, ng))
                    + dot(&p.e3, j, x.dd2.r(0, ng));
                modal += 0.5 * (nv * nv - fin * fin);
                x.v.w(k, k + 1)[0] = fin;
                if !damp {
                    x.trace.w(step * r + j, step * r + j + 1)[0] = fin;
                }
            }
            for g in part.g0..part.g1 {
                let (s, e) = p.group_modes[g];
                let b = dot(&p.r5, g, x.fv.r(n, n + r))
                    + dot(&p.r6, g, x.rea2.r(0, ng))
                    + dot(&p.r2, g, x.src2.r(0, ng));
                modal += string_loss(p, (s, e), b, x.v.w(s, e));
            }
            if damp {
                stations(p, x, st, part);
                damper += chains(p, x, st, part, 0);
            } else {
                sources(p, x, part);
            }
            barrier.wait();
            if damp {
                damper += superpose(p, x, st, part, 0, (x.v, n), (x.vbt, 0));
                damper += chains(p, x, st, part, 1);
                barrier.wait();
                damper += superpose(p, x, st, part, 1, (x.vbt, 0), (x.v, n));
                for j in part.j0..part.j1 {
                    x.trace.w(step * r + j, step * r + j + 1)[0] = x.v.r(n + j, n + j + 1)[0];
                }
                // the next step's first damper call continues from these station speeds
                if step + 1 < steps {
                    damper += chains(p, x, st, part, 0);
                } else {
                    apply(p, x, st, part);
                }
                barrier.wait();
            }
        }
    }
    (modal, damper)
}

pub struct Cores(pub Parallel);
impl Backend for Cores {
    fn name(&self) -> String {
        format!("CPU ({} threads, phase emulator)", self.0.threads())
    }
    fn load_dampers(&mut self, _: &Program) -> Result<(), String> {
        Ok(())
    }
    fn run(
        &mut self,
        p: &Program,
        q: &mut [f64],
        v: &mut [f64],
        steps: usize,
        trace: &mut [f64],
        positions: &mut [f64],
    ) -> Result<(f64, f64), String> {
        Ok(self.0.run(p, q, v, steps, trace, positions))
    }
}
