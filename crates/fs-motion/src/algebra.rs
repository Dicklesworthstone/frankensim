//! Taylor-model multivectors over PGA structure constants extracted
//! from fs-ga at runtime.
//!
//! The geometric product on `fs_ga::Pga` is bilinear with structure
//! constants in {−1, 0, +1}. Rather than transcribing the Cayley table
//! (a sign-error trap), [`gp_table`] extracts it once by multiplying
//! basis blades through fs-ga itself; the conformance suite then
//! falsifies the extracted table against `Motor::transform_point`.
//! With the table in hand, a motor whose components are univariate
//! Taylor models in time supports rigorous products, reversal, and
//! sandwich actions — every rounding and truncation lands in the
//! models' remainders, so enclosures stay sound.

use crate::MotionError;
use fs_ga::Pga;
use fs_ivl::{Interval, TaylorModel1};
use std::sync::OnceLock;

/// Number of PGA basis blades.
pub const BLADES: usize = 16;

/// Blade index of the scalar part.
pub const SCALAR: usize = 0;

/// Blade index of the pseudoscalar `e0123`.
pub const PSEUDO: usize = 15;

/// Homogeneous point blade indices in fs-ga's packed order, and the
/// component mapping used by `Motor::transform_point`: a Euclidean
/// point `(x, y, z)` embeds as `[7] = −z, [11] = y, [13] = −x,
/// [14] = 1` (weight).
pub const POINT_BLADE_NZ: usize = 7;
/// See [`POINT_BLADE_NZ`].
pub const POINT_BLADE_Y: usize = 11;
/// See [`POINT_BLADE_NZ`].
pub const POINT_BLADE_NX: usize = 13;
/// See [`POINT_BLADE_NZ`].
pub const POINT_BLADE_W: usize = 14;

/// One nonzero structure-constant term: `out[k] += sign · a[i] · b[j]`.
#[derive(Debug, Clone, Copy)]
pub struct GpTerm {
    /// Left blade index.
    pub i: u8,
    /// Right blade index.
    pub j: u8,
    /// Output blade index.
    pub k: u8,
    /// ±1.
    pub sign: i8,
}

/// The extracted product table plus per-blade reverse signs.
#[derive(Debug)]
pub struct GaTables {
    /// Sparse geometric-product terms (at most one output per pair in
    /// a degenerate Clifford algebra; null products are absent).
    pub gp: Vec<GpTerm>,
    /// `reverse(e_k) = rev_sign[k] · e_k`.
    pub rev_sign: [i8; BLADES],
}

fn basis(i: usize) -> Pga {
    Pga::blade(i, 1.0)
}

fn extract_tables() -> GaTables {
    let mut gp = Vec::new();
    for i in 0..BLADES {
        for j in 0..BLADES {
            let prod = basis(i).gp(&basis(j));
            for (k, &coeff) in prod.0.iter().enumerate() {
                if coeff == 1.0 || coeff == -1.0 {
                    gp.push(GpTerm {
                        i: i as u8,
                        j: j as u8,
                        k: k as u8,
                        sign: if coeff == 1.0 { 1 } else { -1 },
                    });
                } else {
                    debug_assert!(
                        coeff == 0.0,
                        "fs-ga basis product produced non-unit coefficient {coeff}"
                    );
                }
            }
        }
    }
    let mut rev_sign = [1i8; BLADES];
    for (k, slot) in rev_sign.iter_mut().enumerate() {
        let rev = basis(k).reverse();
        let coeff = rev.0[k];
        debug_assert!(coeff == 1.0 || coeff == -1.0);
        *slot = if coeff < 0.0 { -1 } else { 1 };
    }
    GaTables { gp, rev_sign }
}

/// The extracted structure constants (built once, deterministically).
pub fn ga_tables() -> &'static GaTables {
    static TABLES: OnceLock<GaTables> = OnceLock::new();
    TABLES.get_or_init(extract_tables)
}

/// A PGA multivector whose sixteen components are univariate Taylor
/// models sharing one time domain and order. `nonzero` tracks
/// structural zeros so products skip dead terms.
#[derive(Debug, Clone)]
pub struct TmMv {
    comp: Vec<TaylorModel1>,
    nonzero: [bool; BLADES],
    domain: Interval,
    order: usize,
}

impl TmMv {
    /// The zero multivector on `domain`.
    pub fn zero(domain: Interval, order: usize) -> Result<TmMv, MotionError> {
        let zero = TaylorModel1::constant(0.0, domain, order)?;
        Ok(TmMv {
            comp: vec![zero; BLADES],
            nonzero: [false; BLADES],
            domain,
            order,
        })
    }

    /// A constant multivector (every component a constant model).
    pub fn constant(value: &Pga, domain: Interval, order: usize) -> Result<TmMv, MotionError> {
        let mut out = TmMv::zero(domain, order)?;
        for (k, &v) in value.0.iter().enumerate() {
            if v != 0.0 {
                out.set(k, TaylorModel1::constant(v, domain, order)?)?;
            }
        }
        Ok(out)
    }

    /// The shared time domain.
    #[must_use]
    pub fn domain(&self) -> Interval {
        self.domain
    }

    /// The shared Taylor order.
    #[must_use]
    pub fn order(&self) -> usize {
        self.order
    }

    /// Install a component model (must share domain and order).
    pub fn set(&mut self, k: usize, tm: TaylorModel1) -> Result<(), MotionError> {
        if tm.order() != self.order || tm.domain() != self.domain {
            return Err(MotionError::MixedModelShape {
                blade: k,
                expected_order: self.order,
                got_order: tm.order(),
            });
        }
        let bound = tm.bound();
        self.nonzero[k] = !(bound.lo() == 0.0 && bound.hi() == 0.0);
        self.comp[k] = tm;
        Ok(())
    }

    /// Borrow a component model.
    #[must_use]
    pub fn component(&self, k: usize) -> &TaylorModel1 {
        &self.comp[k]
    }

    /// Whether a component is structurally nonzero.
    #[must_use]
    pub fn is_nonzero(&self, k: usize) -> bool {
        self.nonzero[k]
    }

    /// Rigorous geometric product through the extracted table.
    pub fn gp(&self, rhs: &TmMv) -> Result<TmMv, MotionError> {
        if rhs.domain != self.domain || rhs.order != self.order {
            return Err(MotionError::MixedModelShape {
                blade: 0,
                expected_order: self.order,
                got_order: rhs.order,
            });
        }
        let mut acc: Vec<Option<TaylorModel1>> = vec![None; BLADES];
        for term in &ga_tables().gp {
            let (i, j, k) = (term.i as usize, term.j as usize, term.k as usize);
            if !self.nonzero[i] || !rhs.nonzero[j] {
                continue;
            }
            let mut prod = self.comp[i].try_mul(&rhs.comp[j])?;
            if term.sign < 0 {
                prod = prod.scale(-1.0)?;
            }
            acc[k] = Some(match acc[k].take() {
                None => prod,
                Some(prev) => prev.try_add(&prod)?,
            });
        }
        let mut out = TmMv::zero(self.domain, self.order)?;
        for (k, slot) in acc.into_iter().enumerate() {
            if let Some(tm) = slot {
                out.set(k, tm)?;
            }
        }
        Ok(out)
    }

    /// Reversal (per-blade sign flips from the extracted table).
    pub fn reverse(&self) -> Result<TmMv, MotionError> {
        let mut out = self.clone();
        for k in 0..BLADES {
            if out.nonzero[k] && ga_tables().rev_sign[k] < 0 {
                let flipped = out.comp[k].scale(-1.0)?;
                out.set(k, flipped)?;
            }
        }
        Ok(out)
    }

    /// Negate every component (double-cover sign choice).
    pub fn negate(&self) -> Result<TmMv, MotionError> {
        let mut out = self.clone();
        for k in 0..BLADES {
            if out.nonzero[k] {
                let flipped = out.comp[k].scale(-1.0)?;
                out.set(k, flipped)?;
            }
        }
        Ok(out)
    }

    /// Interval enclosures of all components over a subinterval of the
    /// domain.
    pub fn eval_all(&self, t: Interval) -> Result<[Interval; BLADES], MotionError> {
        if !self.domain.encloses(t) {
            return Err(MotionError::OutOfDomain {
                lo: t.lo(),
                hi: t.hi(),
                domain_lo: self.domain.lo(),
                domain_hi: self.domain.hi(),
            });
        }
        let mut out = [Interval::point(0.0); BLADES];
        for (k, slot) in out.iter_mut().enumerate() {
            if self.nonzero[k] {
                *slot = self.comp[k].eval_interval(t);
            }
        }
        Ok(out)
    }

    /// Upper bound of the versor defect `‖M M̃ − 1‖∞` over the whole
    /// domain, computed by rigorous model arithmetic. For a true unit
    /// versor path this is rounding-plus-truncation noise; a broken
    /// construction (non-unit generator, sign mixing) shows up here
    /// instead of being silently absorbed.
    pub fn versor_defect(&self) -> Result<f64, MotionError> {
        let mm = self.gp(&self.reverse()?)?;
        let mut worst = 0.0f64;
        for k in 0..BLADES {
            let b = if mm.is_nonzero(k) {
                mm.component(k).bound()
            } else {
                Interval::point(0.0)
            };
            let dev = if k == SCALAR {
                (b - Interval::point(1.0)).abs_bound()
            } else {
                b.abs_bound()
            };
            if dev.hi() > worst {
                worst = dev.hi();
            }
        }
        Ok(worst)
    }
}

/// Embed a Euclidean point into the homogeneous point blades (the
/// component mapping `Motor::transform_point` uses).
#[must_use]
pub fn point_to_mv(x: f64, y: f64, z: f64) -> Pga {
    let mut p = Pga::zero();
    p.0[POINT_BLADE_NZ] = -z;
    p.0[POINT_BLADE_Y] = y;
    p.0[POINT_BLADE_NX] = -x;
    p.0[POINT_BLADE_W] = 1.0;
    p
}

/// Extract interval point coordinates from sandwich output enclosures,
/// dividing out the homogeneous weight. Refuses when the weight
/// enclosure contains zero. Uniform versor scaling cancels here, so
/// the result is exact for the constructed component path even when
/// its norm drifts from one.
pub fn homogeneous_point(
    enclosures: &[Interval; BLADES],
) -> Result<[Interval; 3], MotionError> {
    let w = enclosures[POINT_BLADE_W];
    if w.lo() <= 0.0 && w.hi() >= 0.0 {
        return Err(MotionError::DegenerateWeight {
            lo: w.lo(),
            hi: w.hi(),
        });
    }
    let x = (Interval::point(0.0) - enclosures[POINT_BLADE_NX]) / w;
    let y = enclosures[POINT_BLADE_Y] / w;
    let z = (Interval::point(0.0) - enclosures[POINT_BLADE_NZ]) / w;
    Ok([x, y, z])
}
