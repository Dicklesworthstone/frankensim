use super::{
    Cell, Error, MAX_COMPONENTS, MAX_SPATIAL_COMPONENTS, MixtureGas, NormalCell, Primitive,
    SpatialCell, poll, positive, sum,
};

/// Wave-speed receipt accompanying an allocation-free species/energy face flux.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FaceFlux {
    /// Wave-speed bound [m/s], relative to the face for moving-wall APIs and
    /// in the laboratory frame for stationary-face APIs.
    pub signal_speed: f64,
    /// Whether an inadmissible HLLC contact fell back to an admitted HLL state.
    pub used_hll_fallback: bool,
}
pub(super) type Face = FaceFlux;

impl MixtureGas {
    /// Compute a shared oriented flux into caller-provided storage.
    ///
    /// `out` has exactly `species_count()+2` entries: every species mass flux
    /// [kg/(m² s)], momentum flux [Pa], and species-referenced energy flux [W/m²].
    /// The physical energy flux adds `sum(e_ref_i * species_flux_i)`. This
    /// reusable closure and its successful numerical path allocate nothing.
    /// `out` is usable only on success; cancellation may leave scratch values.
    pub fn hllc_flux(
        &self,
        left: Cell<'_>,
        right: Cell<'_>,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let n = self.species_count();
        let m = n + 2;
        if out.len() != m {
            return Err(Error::Invalid("face flux output dimension"));
        }
        let pl = self.primitive(left, &mut cancelled)?;
        let pr = self.primitive(right, &mut cancelled)?;
        let mut ql = [0.0; MAX_COMPONENTS];
        let mut qr = [0.0; MAX_COMPONENTS];
        ql[..n].copy_from_slice(left.partial_densities);
        qr[..n].copy_from_slice(right.partial_densities);
        ql[n] = left.momentum;
        qr[n] = right.momentum;
        ql[n + 1] = left.energy;
        qr[n + 1] = right.energy;
        self.interface(&ql[..m], pl, &qr[..m], pr, out, &mut cancelled)
    }

    /// Evaluate all species, momentum and energy fluxes with one shared set of
    /// HLLC wave speeds. The species energy reference is transported by the
    /// corresponding species flux, not by an independently rounded mass flux.
    pub(super) fn interface(
        &self,
        left: &[f64],
        pl: Primitive,
        right: &[f64],
        pr: Primitive,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Face, Error> {
        self.interface_components::<0>(left, pl, right, pr, out, cancelled)
    }

    /// Laboratory ALE flux through a face translating at the supplied axial
    /// speed. Every species retains its own original energy reference.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn moving_interface(
        &self,
        left: &[f64],
        pl: Primitive,
        right: &[f64],
        pr: Primitive,
        velocity: f64,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Face, Error> {
        if velocity == 0.0 {
            return self.interface(left, pl, right, pr, out, cancelled);
        }
        let n = self.species_count();
        let m = n + 2;
        if !velocity.is_finite() || left.len() != m || right.len() != m || out.len() != m {
            return Err(Error::Invalid("moving mixture face speed or dimensions"));
        }
        let mut ql = [0.0; MAX_COMPONENTS];
        let mut qr = [0.0; MAX_COMPONENTS];
        poll(cancelled)?;
        let pl = self.state_at(
            &left[..n],
            pl.velocity - velocity,
            pl.pressure,
            &mut ql[..m],
        )?;
        poll(cancelled)?;
        let pr = self.state_at(
            &right[..n],
            pr.velocity - velocity,
            pr.pressure,
            &mut qr[..m],
        )?;
        let face = self.interface(&ql[..m], pl, &qr[..m], pr, out, cancelled)?;
        // Total mass has no independently evolved field. Use the same ordered
        // species-flux sum in both Galilean momentum and energy transformations.
        let mass_flux = sum(out[..n].iter().copied());
        let relative_momentum = out[n];
        out[n] = relative_momentum + velocity * mass_flux;
        out[n + 1] += velocity * relative_momentum + (0.5 * velocity * velocity) * mass_flux;
        if out.iter().any(|value| !value.is_finite()) {
            return Err(Error::Invalid("moving mixture laboratory flux overflow"));
        }
        poll(cancelled)?;
        Ok(face)
    }

    // Axial, planar and volumetric consumers use these same wave speeds, star
    // states and bounded HLL fallback. Every transverse momentum is transported
    // by the SAME waves, with its kinetic energy in every source-domain gate.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn interface_components<const TANGENTS: usize>(
        &self,
        left: &[f64],
        pl: Primitive,
        right: &[f64],
        pr: Primitive,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Face, Error> {
        poll(cancelled)?;
        let n = self.species_count();
        let m = n + 2 + TANGENTS;
        let energy = m - 1;
        let mut fl = [0.0; MAX_SPATIAL_COMPONENTS];
        let mut fr = [0.0; MAX_SPATIAL_COMPONENTS];
        physical(left, pl, &mut fl[..m], n, energy);
        physical(right, pr, &mut fr[..m], n, energy);
        let sl = (pl.velocity - pl.sound_speed).min(pr.velocity - pr.sound_speed);
        let sr = (pl.velocity + pl.sound_speed).max(pr.velocity + pr.sound_speed);
        let speed = sl.abs().max(sr.abs());
        if !positive(speed) || !fl[..m].iter().chain(&fr[..m]).all(|x| x.is_finite()) {
            return Err(Error::Invalid(
                "unrepresentable mixture flux or signal speed",
            ));
        }
        // Equal source-admitted traces have exactly the physical flux. A
        // redundant star-state reconstruction can round (s * energy) / s away
        // from energy and manufacture transport in an otherwise uniform flow.
        if left == right && pl == pr {
            out.copy_from_slice(&fl[..m]);
            return checked(out, speed, false);
        }
        // A stationary pressure contact can have arbitrary composition,
        // density and tangential shear, but transports only normal pressure.
        // Match the frozen-gas owner without treating nearby pressures as equal.
        // Moving unequal-caloric contacts retain their existing no-claim limit.
        if pl.velocity == 0.0 && pr.velocity == 0.0 && pl.pressure == pr.pressure {
            out.fill(0.0);
            out[n] = pl.pressure;
            return checked(out, speed, false);
        }
        if sl >= 0.0 {
            out.copy_from_slice(&fl[..m]);
            return Ok(Face {
                signal_speed: speed,
                used_hll_fallback: false,
            });
        }
        if sr <= 0.0 {
            out.copy_from_slice(&fr[..m]);
            return Ok(Face {
                signal_speed: speed,
                used_hll_fallback: false,
            });
        }
        let al = pl.density * (sl - pl.velocity);
        let ar = pr.density * (sr - pr.velocity);
        let contact = (pr.pressure - pl.pressure + al * pl.velocity - ar * pr.velocity) / (al - ar);
        let pstar = pl.pressure + al * (contact - pl.velocity);
        if contact.is_finite() && sl < contact && contact < sr && positive(pstar) {
            let mut ls = [0.0; MAX_SPATIAL_COMPONENTS];
            let mut rs = [0.0; MAX_SPATIAL_COMPONENTS];
            let star = |q: &[f64], p: Primitive, s: f64, target: &mut [f64]| {
                let factor = (s - p.velocity) / (s - contact);
                for i in 0..n {
                    target[i] = q[i] * factor;
                }
                target[n] = p.density * factor * contact;
                for k in n + 1..energy {
                    target[k] = q[k] * factor;
                }
                target[energy] = ((s - p.velocity) * q[energy] - p.pressure * p.velocity
                    + pstar * contact)
                    / (s - contact);
            };
            star(left, pl, sl, &mut ls[..m]);
            star(right, pr, sr, &mut rs[..m]);
            if self.admissible_face::<TANGENTS>(&ls[..m])
                && self.admissible_face::<TANGENTS>(&rs[..m])
            {
                let (q, qs, f, s) = if contact >= 0.0 {
                    (left, &ls[..m], &fl[..m], sl)
                } else {
                    (right, &rs[..m], &fr[..m], sr)
                };
                for k in 0..m {
                    out[k] = f[k] + s * (qs[k] - q[k]);
                }
                return checked(out, speed, false);
            }
        }
        let mut middle = [0.0; MAX_SPATIAL_COMPONENTS];
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
            for k in 0..m {
                middle[k] = (r * right[k] - l * left[k] + fl[k] - fr[k]) / (r - l);
            }
            if self.admissible_face::<TANGENTS>(&middle[..m]) {
                for k in 0..m {
                    out[k] = (r * fl[k] - l * fr[k] + l * r * (right[k] - left[k])) / (r - l);
                }
                return checked(out, l.abs().max(r.abs()), true);
            }
        }
        Err(Error::Invalid(
            "mixture HLL intermediate state outside source domain",
        ))
    }

    fn admissible_face<const TANGENTS: usize>(&self, q: &[f64]) -> bool {
        if TANGENTS == 0 {
            return self.admissible(q);
        }
        let n = self.species_count();
        if TANGENTS == 1 {
            self.admissible_normal_cell(NormalCell {
                partial_densities: &q[..n],
                normal_momentum: q[n],
                transverse_momentum: q[n + 1],
                energy: q[n + 2],
            })
        } else {
            self.admissible_spatial_cell(SpatialCell {
                partial_densities: &q[..n],
                normal_momentum: q[n],
                transverse_momentum: [q[n + 1], q[n + 2]],
                energy: q[n + 3],
            })
        }
    }
}
fn physical(q: &[f64], p: Primitive, out: &mut [f64], n: usize, energy: usize) {
    for i in 0..n {
        out[i] = q[i] * p.velocity;
    }
    out[n] = q[n] * p.velocity + p.pressure;
    for k in n + 1..energy {
        out[k] = q[k] * p.velocity;
    }
    out[energy] = (q[energy] + p.pressure) * p.velocity;
}
fn checked(out: &[f64], signal_speed: f64, used_hll_fallback: bool) -> Result<Face, Error> {
    if !out.iter().all(|x| x.is_finite()) {
        return Err(Error::Invalid("non-finite mixture numerical flux"));
    }
    Ok(Face {
        signal_speed,
        used_hll_fallback,
    })
}
