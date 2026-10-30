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

/// Conservative trace in a two-dimensional face's orthonormal frame.
///
/// Spatial owners rotate momentum into this frame. `energy` includes both
/// normal and transverse kinetic energy, with the same fixed caloric reference
/// as [`Conserved`]. The Riemann closure owns all thermochemical decisions.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NormalConserved {
    /// Mass density [kg/m³].
    pub density: f64,
    /// Momentum density along the oriented face normal [kg/(m² s)].
    pub normal_momentum: f64,
    /// Momentum density along its unit tangent [kg/(m² s)].
    pub transverse_momentum: f64,
    /// Referenced total energy, including both kinetic components [J/m³].
    pub energy: f64,
}

impl From<Conserved> for NormalConserved {
    fn from(q: Conserved) -> Self {
        Self {
            density: q.density,
            normal_momentum: q.momentum,
            transverse_momentum: 0.0,
            energy: q.energy,
        }
    }
}

impl NormalConserved {
    fn array(self) -> [f64; 4] {
        [
            self.density,
            self.normal_momentum,
            self.transverse_momentum,
            self.energy,
        ]
    }

    fn from_array(q: [f64; 4]) -> Self {
        Self {
            density: q[0],
            normal_momentum: q[1],
            transverse_momentum: q[2],
            energy: q[3],
        }
    }

    fn normal_state(self) -> Conserved {
        Conserved {
            density: self.density,
            momentum: self.normal_momentum,
            energy: self.energy
                - 0.5 * self.transverse_momentum * (self.transverse_momentum / self.density),
        }
    }
}

/// Reconstructed planar trace; all caloric properties remain source-computed.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NormalPrimitive {
    /// Mass density [kg/m³].
    pub density: f64,
    /// Velocity along the oriented face normal [m/s].
    pub normal_velocity: f64,
    /// Velocity along its unit tangent [m/s].
    pub transverse_velocity: f64,
    /// Absolute pressure [Pa].
    pub pressure: f64,
}

/// Shared HLLC/HLL flux in the same normal/tangent frame as its traces.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NormalFaceFlux {
    /// Mass flux [kg/(m² s)].
    pub mass: f64,
    /// Normal momentum flux, including pressure [Pa].
    pub normal_momentum: f64,
    /// Transverse momentum flux [Pa].
    pub transverse_momentum: f64,
    /// Referenced total-energy flux, including transverse kinetic energy [W/m²].
    pub energy: f64,
    /// Largest magnitude of the actual numerical face wave speeds [m/s].
    pub signal_speed: f64,
    /// Whether this face used the more dissipative HLL fallback.
    pub used_hll_fallback: bool,
}

impl From<NormalFaceFlux> for FaceFlux {
    fn from(f: NormalFaceFlux) -> Self {
        Self {
            mass: f.mass,
            momentum: f.normal_momentum,
            energy: f.energy,
            signal_speed: f.signal_speed,
            used_hll_fallback: f.used_hll_fallback,
        }
    }
}

impl FrozenGas {
    /// Build a planar trace from density, normal/transverse velocity and pressure.
    ///
    /// Refuses a kinetic-energy scale that cannot retain the admitted sensible
    /// energy precision in the full conservative energy representation.
    pub fn normal_conserved(
        &self,
        density: f64,
        normal_velocity: f64,
        transverse_velocity: f64,
        pressure: f64,
    ) -> Result<NormalConserved, Error> {
        self.normal_state_at(NormalPrimitive {
            density,
            normal_velocity,
            transverse_velocity,
            pressure,
        })
        .map(|(state, _)| state)
    }

    fn normal_state_at(
        &self,
        state: NormalPrimitive,
    ) -> Result<(NormalConserved, Primitive), Error> {
        let NormalPrimitive {
            density,
            normal_velocity,
            transverse_velocity,
            pressure,
        } = state;
        if !transverse_velocity.is_finite() {
            return Err(Error::Invalid("non-finite transverse velocity"));
        }
        let (normal, primitive) = self.state_at(density, normal_velocity, pressure)?;
        let q = NormalConserved {
            density,
            normal_momentum: normal.momentum,
            transverse_momentum: density * transverse_velocity,
            energy: normal.energy + density * (0.5 * transverse_velocity * transverse_velocity),
        };
        let reduced = q.normal_state();
        let recovered_velocity = reduced.momentum / reduced.density;
        let recovered_energy =
            reduced.energy / reduced.density - 0.5 * recovered_velocity * recovered_velocity;
        let expected_energy = self.sensible(primitive.temperature)?;
        let domain = self.domain();
        let tolerance = domain.absolute_energy_tolerance
            + domain.relative_energy_tolerance * expected_energy.abs();
        if !self.admissible_normal_state(q)
            || (recovered_energy - expected_energy).abs() > tolerance
        {
            return Err(Error::Invalid(
                "planar conservative state cannot retain the declared caloric precision",
            ));
        }
        Ok((q, primitive))
    }

    /// Decode a planar trace through the original caloric source.
    ///
    /// The returned `Primitive::velocity` is the normal component; the caller
    /// obtains transverse velocity from `transverse_momentum / density`.
    pub fn normal_primitive(
        &self,
        state: NormalConserved,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        poll(&mut cancelled)?;
        if !state.array().iter().all(|x| x.is_finite()) || !positive(state.density) {
            return Err(Error::Invalid("invalid planar conservative state"));
        }
        self.decode(state.normal_state(), None, &mut cancelled)
    }

    /// HLLC flux for a planar face, with full-energy HLL fallback admission.
    ///
    /// Transverse momentum and its kinetic energy use the same waves as mass,
    /// normal momentum and heat. In particular, HLL shear mixing contributes
    /// to the intermediate state's internal energy before source-domain checks;
    /// appending an upwind transverse flux to a one-dimensional HLL result would
    /// miss both that mixing and the conservative transverse flux.
    pub fn hllc_normal_flux(
        &self,
        left: NormalConserved,
        right: NormalConserved,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NormalFaceFlux, Error> {
        let pl = self.normal_primitive(left, &mut cancelled)?;
        let pr = self.normal_primitive(right, &mut cancelled)?;
        self.normal_interface(left, pl, right, pr, &mut cancelled)
    }

    /// Evaluate reconstructed traces without inverting their just-computed energy.
    ///
    /// Each trace passes the same source-domain and caloric-precision admission
    /// as `normal_conserved`. Pressure, temperature and sound speed are computed
    /// by this owner, then reused by the unchanged HLLC/HLL kernel. No caller
    /// supplies cached thermodynamics or an approximate equation of state.
    pub fn hllc_normal_flux_from_primitive(
        &self,
        left: NormalPrimitive,
        right: NormalPrimitive,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NormalFaceFlux, Error> {
        poll(&mut cancelled)?;
        let (left, pl) = self.normal_state_at(left)?;
        poll(&mut cancelled)?;
        let (right, pr) = self.normal_state_at(right)?;
        self.normal_interface(left, pl, right, pr, &mut cancelled)
    }

    /// Reflect a planar trace at a stationary, adiabatic slip wall.
    ///
    /// `outward_normal` is -1 or +1 relative to the supplied normal axis. The
    /// original thermally perfect shock/rarefaction wall closure supplies the
    /// pressure; no mass, transverse momentum or energy crosses the wall.
    pub fn stationary_wall_normal_flux(
        &self,
        state: NormalConserved,
        outward_normal: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NormalFaceFlux, Error> {
        let p = self.normal_primitive(state, &mut cancelled)?;
        self.normal_wall_interface(state, p, outward_normal, &mut cancelled)
    }

    /// Source-admitted reconstructed trace at a stationary adiabatic slip wall.
    /// Reuses source-computed thermodynamics without an energy inversion, while
    /// preserving the same thermally perfect wall shock/rarefaction solve.
    pub fn stationary_wall_normal_flux_from_primitive(
        &self,
        state: NormalPrimitive,
        outward_normal: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NormalFaceFlux, Error> {
        poll(&mut cancelled)?;
        let (state, p) = self.normal_state_at(state)?;
        self.normal_wall_interface(state, p, outward_normal, &mut cancelled)
    }

    fn normal_wall_interface(
        &self,
        state: NormalConserved,
        primitive: Primitive,
        outward_normal: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NormalFaceFlux, Error> {
        let wall = self.wall_flux(state.normal_state(), primitive, outward_normal, cancelled)?;
        Ok(NormalFaceFlux {
            mass: 0.0,
            normal_momentum: wall.momentum,
            transverse_momentum: 0.0,
            energy: 0.0,
            signal_speed: wall.signal_speed,
            used_hll_fallback: false,
        })
    }

    fn admissible_normal_state(&self, state: NormalConserved) -> bool {
        state.array().iter().all(|x| x.is_finite())
            && positive(state.density)
            && self.admissible_state(state.normal_state())
    }

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
        self.normal_interface(left.into(), pl, right.into(), pr, cancelled)
            .map(FaceFlux::from)
    }

    fn normal_interface(
        &self,
        left: NormalConserved,
        pl: Primitive,
        right: NormalConserved,
        pr: Primitive,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NormalFaceFlux, Error> {
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
            let star = |q: NormalConserved, p: Primitive, s: f64| NormalConserved {
                density: q.density * (s - p.velocity) / (s - contact),
                normal_momentum: q.density * (s - p.velocity) / (s - contact) * contact,
                transverse_momentum: q.transverse_momentum * (s - p.velocity) / (s - contact),
                energy: ((s - p.velocity) * q.energy - p.pressure * p.velocity + pstar * contact)
                    / (s - contact),
            };
            let ls = star(left, pl, sl);
            let rs = star(right, pr, sr);
            if self.admissible_normal_state(ls) && self.admissible_normal_state(rs) {
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
            let middle = NormalConserved::from_array(core::array::from_fn(|k| {
                (r * qr[k] - l * ql[k] + fl[k] - fr[k]) / (r - l)
            }));
            if self.admissible_normal_state(middle) {
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

fn physical_flux(q: NormalConserved, p: Primitive) -> [f64; 4] {
    [
        q.normal_momentum,
        q.normal_momentum * p.velocity + p.pressure,
        q.transverse_momentum * p.velocity,
        (q.energy + p.pressure) * p.velocity,
    ]
}
fn flux(f: [f64; 4], signal_speed: f64, used_hll_fallback: bool) -> Result<NormalFaceFlux, Error> {
    if !f.iter().all(|x| x.is_finite()) {
        return Err(Error::Invalid("non-finite numerical flux"));
    }
    Ok(NormalFaceFlux {
        mass: f[0],
        normal_momentum: f[1],
        transverse_momentum: f[2],
        energy: f[3],
        signal_speed,
        used_hll_fallback,
    })
}
