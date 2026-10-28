use super::{Conserved, Error, FrozenGas, Primitive, poll, positive};

/// One shared oriented face flux, per unit cross-sectional area.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FaceFlux {
    /// Mass flux [kg/(m² s)].
    pub mass: f64,
    /// Momentum flux [Pa].
    pub momentum: f64,
    /// Referenced total-energy flux [W/m²]. Add `e_ref * mass` for physical energy.
    pub energy: f64,
    /// Largest magnitude of the actual numerical face wave speeds [m/s].
    pub signal_speed: f64,
    /// Whether this face used the more dissipative HLL fallback.
    pub used_hll_fallback: bool,
}
impl FaceFlux {
    pub(super) fn array(self) -> [f64; 3] {
        [self.mass, self.momentum, self.energy]
    }
}

impl FrozenGas {
    /// HLLC flux between two conservative states of this same frozen mixture.
    ///
    /// Left/right Davis wave bounds and the contact Rankine–Hugoniot relations
    /// use each state's actual frozen sound speed. A negative-pressure,
    /// out-of-domain or degenerate contact state falls back to HLL. If necessary,
    /// HLL signal speeds expand by bounded factors of two until its intermediate
    /// state is admitted. These checks do not assert a general entropy theorem.
    pub fn hllc_flux(
        &self,
        left: Conserved,
        right: Conserved,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let pl = self.decode(left, None, &mut cancelled)?;
        let pr = self.decode(right, None, &mut cancelled)?;
        self.interface(left, pl, right, pr, &mut cancelled)
    }

    pub(super) fn interface(
        &self,
        left: Conserved,
        pl: Primitive,
        right: Conserved,
        pr: Primitive,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(cancelled)?;
        let fl = physical_flux(left, pl);
        let fr = physical_flux(right, pr);
        let sl = (pl.velocity - pl.sound_speed).min(pr.velocity - pr.sound_speed);
        let sr = (pl.velocity + pl.sound_speed).max(pr.velocity + pr.sound_speed);
        let speed = sl.abs().max(sr.abs());
        if !positive(speed) || !fl.iter().chain(&fr).all(|x| x.is_finite()) {
            return Err(Error::Invalid("unrepresentable Euler flux or signal speed"));
        }
        if sl >= 0.0 {
            return flux(fl, speed, false);
        }
        if sr <= 0.0 {
            return flux(fr, speed, false);
        }
        let al = left.density * (sl - pl.velocity);
        let ar = right.density * (sr - pr.velocity);
        let contact = (pr.pressure - pl.pressure + al * pl.velocity - ar * pr.velocity) / (al - ar);
        let pstar = pl.pressure + al * (contact - pl.velocity);
        if contact.is_finite() && sl < contact && contact < sr && positive(pstar) {
            let star = |q: Conserved, p: Primitive, s: f64| Conserved {
                density: q.density * (s - p.velocity) / (s - contact),
                momentum: q.density * (s - p.velocity) / (s - contact) * contact,
                energy: ((s - p.velocity) * q.energy - p.pressure * p.velocity + pstar * contact)
                    / (s - contact),
            };
            let ls = star(left, pl, sl);
            let rs = star(right, pr, sr);
            if self.admissible_state(ls) && self.admissible_state(rs) {
                let (q, qs, f, s) = if contact >= 0.0 {
                    (left, ls, fl, sl)
                } else {
                    (right, rs, fr, sr)
                };
                let (q, qs) = (q.array(), qs.array());
                return flux(
                    core::array::from_fn(|k| f[k] + s * (qs[k] - q[k])),
                    speed,
                    false,
                );
            }
        }
        let (ql, qr) = (left.array(), right.array());
        for expansion in 0..=8 {
            poll(cancelled)?;
            let (l, r) = if expansion == 0 {
                (sl, sr)
            } else {
                let a = speed * (1_u32 << (expansion - 1)) as f64;
                (-a, a)
            };
            if !l.is_finite() || !r.is_finite() {
                return Err(Error::Invalid("HLL wave-speed overflow"));
            }
            let middle = Conserved::from_array(core::array::from_fn(|k| {
                (r * qr[k] - l * ql[k] + fl[k] - fr[k]) / (r - l)
            }));
            if self.admissible_state(middle) {
                return flux(
                    core::array::from_fn(|k| {
                        (r * fl[k] - l * fr[k] + l * r * (qr[k] - ql[k])) / (r - l)
                    }),
                    l.abs().max(r.abs()),
                    true,
                );
            }
        }
        Err(Error::Invalid(
            "HLL intermediate state outside admitted source domain",
        ))
    }
}

fn physical_flux(q: Conserved, p: Primitive) -> [f64; 3] {
    [
        q.momentum,
        q.momentum * p.velocity + p.pressure,
        (q.energy + p.pressure) * p.velocity,
    ]
}
fn flux(f: [f64; 3], signal_speed: f64, used_hll_fallback: bool) -> Result<FaceFlux, Error> {
    if !f.iter().all(|x| x.is_finite()) {
        return Err(Error::Invalid("non-finite numerical flux"));
    }
    Ok(FaceFlux {
        mass: f[0],
        momentum: f[1],
        energy: f[2],
        signal_speed,
        used_hll_fallback,
    })
}
