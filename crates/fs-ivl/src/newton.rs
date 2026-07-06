//! Interval Newton and Krawczyk root certification (plan §6.4), plus
//! Lipschitz extraction — the primitives that turn "the solver found a
//! root" into "a root EXISTS, is UNIQUE in this box, and here is the box"
//! (what the word certified means for roots).
//!
//! Semantics: `Certified` is issued ONLY when the contraction lands
//! strictly inside the box (the Newton/Krawczyk existence-uniqueness
//! theorem); everything else is `Possible` (may contain roots, could not
//! certify) — a double root can never be falsely certified (tested).

use crate::Interval;

/// A root-search result box with its certification status.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RootBox {
    /// Exactly one root exists in this box (Newton/Krawczyk contraction
    /// strictly interior — the classical existence + uniqueness theorem).
    Certified(Interval),
    /// The box may contain roots; certification did not succeed at the
    /// subdivision limit (multiple/tangent roots land here — honestly).
    Possible(Interval),
}

impl RootBox {
    /// The underlying box.
    #[must_use]
    pub fn interval(&self) -> Interval {
        match *self {
            RootBox::Certified(iv) | RootBox::Possible(iv) => iv,
        }
    }

    /// Is this a certified (exists + unique) box?
    #[must_use]
    pub fn is_certified(&self) -> bool {
        matches!(self, RootBox::Certified(_))
    }
}

/// One interval-Newton step: N(X) = m − f(m)/F′(X), intersected with X.
/// Returns `None` when N(X) ∩ X is empty (NO root in X — a certificate of
/// absence) and the contraction plus a strict-interior flag otherwise.
fn newton_step<F, D>(f: &F, fp: &D, x: Interval) -> Option<(Interval, bool)>
where
    F: Fn(Interval) -> Interval,
    D: Fn(Interval) -> Interval,
{
    let m = x.midpoint();
    let fm = f(Interval::point(m));
    let d = fp(x);
    if d.contains_zero() {
        // Division yields the whole line: no contraction information.
        return Some((x, false));
    }
    let n = Interval::point(m) - fm / d;
    let contracted = n.intersect(x)?;
    let strict_interior = n.lo() > x.lo() && n.hi() < x.hi();
    Some((contracted, strict_interior))
}

/// Find all roots of `f` in `domain` by recursive bisection with interval
/// Newton contraction. `f` and `fp` are interval extensions of the
/// function and its derivative (fs-ivl arithmetic keeps them rigorous).
/// `min_width` bounds subdivision; boxes still ambiguous at that width
/// come back [`RootBox::Possible`].
#[must_use]
pub fn newton_roots<F, D>(f: &F, fp: &D, domain: Interval, min_width: f64) -> Vec<RootBox>
where
    F: Fn(Interval) -> Interval,
    D: Fn(Interval) -> Interval,
{
    let mut out = Vec::new();
    let mut stack = vec![domain];
    while let Some(x) = stack.pop() {
        // Exclusion test first: 0 ∉ F(X) means no root here.
        if !f(x).contains_zero() {
            continue;
        }
        // Newton contraction loop.
        let mut cur = x;
        let mut certified = false;
        let mut absent = false;
        for _ in 0..64 {
            match newton_step(f, fp, cur) {
                None => {
                    // Empty intersection: certificate of ABSENCE.
                    absent = true;
                    break;
                }
                Some((next, strict)) => {
                    if strict {
                        certified = true;
                    }
                    let stalled = next.width() >= cur.width() * 0.9;
                    cur = next;
                    if stalled {
                        break;
                    }
                }
            }
        }
        if absent {
            continue;
        }
        if certified {
            // Polish: iterate to a tight certified box.
            for _ in 0..64 {
                match newton_step(f, fp, cur) {
                    Some((next, _)) if next.width() < cur.width() => cur = next,
                    _ => break,
                }
            }
            out.push(RootBox::Certified(cur));
        } else if cur.width() <= min_width {
            out.push(RootBox::Possible(cur));
        } else {
            let m = cur.midpoint();
            stack.push(Interval::new(cur.lo(), m));
            stack.push(Interval::new(m, cur.hi()));
        }
    }
    // Deterministic presentation: sort by lower bound; merge adjacent
    // Possible boxes that share an endpoint (bisection artifacts).
    out.sort_by(|a, b| a.interval().lo().partial_cmp(&b.interval().lo()).unwrap());
    let mut merged: Vec<RootBox> = Vec::new();
    for r in out {
        if let (Some(RootBox::Possible(prev)), RootBox::Possible(cur)) = (merged.last().copied(), r)
            && prev.hi() >= cur.lo()
        {
            let joined = prev.hull(cur);
            *merged.last_mut().unwrap() = RootBox::Possible(joined);
            continue;
        }
        merged.push(r);
    }
    merged
}

/// One Krawczyk step: K(X) = m − y·f(m) + (1 − y·F′(X))·(X − m) with
/// y = 1/f′(m). Same certification semantics as interval Newton.
#[must_use]
pub fn krawczyk_step<F, D>(f: &F, fp: &D, x: Interval) -> Option<(Interval, bool)>
where
    F: Fn(Interval) -> Interval,
    D: Fn(Interval) -> Interval,
{
    let m = x.midpoint();
    let fm = f(Interval::point(m));
    let dm = fp(Interval::point(m));
    if dm.contains_zero() {
        return Some((x, false));
    }
    let y = Interval::point(1.0) / dm;
    let k =
        Interval::point(m) - y * fm + (Interval::point(1.0) - y * fp(x)) * (x - Interval::point(m));
    let contracted = k.intersect(x)?;
    let strict = k.lo() > x.lo() && k.hi() < x.hi();
    Some((contracted, strict))
}

/// A CERTIFIED Lipschitz constant for f over `domain`: the magnitude of
/// the derivative enclosure, rounded up — the primitive fs-geom's
/// certified-Lipschitz chart contract consumes. Returns +∞ when the
/// derivative enclosure is unbounded (honest, never understated).
#[must_use]
pub fn lipschitz_bound<D>(fp: &D, domain: Interval) -> f64
where
    D: Fn(Interval) -> Interval,
{
    let d = fp(domain);
    let mag = d.lo().abs().max(d.hi().abs());
    if mag.is_finite() {
        fs_math::next_up(mag)
    } else {
        f64::INFINITY
    }
}
