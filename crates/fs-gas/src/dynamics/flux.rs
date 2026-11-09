use super::{Conserved, Error, FrozenGas, Primitive, poll, positive};

mod spatial;
pub use spatial::{SpatialConserved, SpatialFaceFlux, SpatialPrimitive};

// The component count is fixed at monomorphization: [rho, rho*u_n,
// rho*u_t..., E]. Planar and volumetric faces use this same numerical kernel.
trait Trace<const N: usize>: Copy {
    fn array(self) -> [f64; N];
    fn from_array(values: [f64; N]) -> Self;
    fn normal_state(self) -> Conserved;
}

struct NumericalFlux<const N: usize> {
    values: [f64; N],
    signal_speed: f64,
    used_hll_fallback: bool,
}

/// One shared oriented face flux, per unit cross-sectional area.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FaceFlux {
    /// Mass flux [kg/(m² s)].
    pub mass: f64,
    /// Momentum flux [Pa].
    pub momentum: f64,
    /// Referenced total-energy flux [W/m²]. Add `e_ref * mass` for physical energy.
    pub energy: f64,
    /// Wave-speed bound [m/s], relative to the face for moving-wall APIs and
    /// in the laboratory frame for stationary-face APIs.
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

impl Trace<4> for NormalConserved {
    fn array(self) -> [f64; 4] {
        self.array()
    }
    fn from_array(values: [f64; 4]) -> Self {
        Self::from_array(values)
    }
    fn normal_state(self) -> Conserved {
        self.normal_state()
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

/// Shared face flux in the same normal/tangent frame as its traces.
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
    /// Wave-speed bound [m/s], relative to the face for moving-wall APIs and
    /// in the laboratory frame for stationary-face APIs.
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

impl From<NumericalFlux<4>> for NormalFaceFlux {
    fn from(f: NumericalFlux<4>) -> Self {
        Self {
            mass: f.values[0],
            normal_momentum: f.values[1],
            transverse_momentum: f.values[2],
            energy: f.values[3],
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
        self.admit_trace_precision(
            q,
            primitive,
            "planar conservative state cannot retain the declared caloric precision",
        )
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
        self.decode_trace(state, "invalid planar conservative state", &mut cancelled)
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

    /// Conservative planar flux through a face translating at `face_velocity`
    /// along the supplied normal axis [m/s].
    ///
    /// The original Riemann problem is sampled in the face frame, then mass,
    /// both momenta and full referenced energy are returned in laboratory
    /// coordinates. This is `F - face_velocity * U`, not a stationary-face
    /// flux with an independently upwinded geometric correction. Tangential
    /// momentum follows the same contact and HLL waves as all other storage.
    /// The wave-speed bound is relative to the face. Tangential face motion
    /// does not change the normal transport through an inviscid planar face.
    ///
    /// Zero face velocity retains the stationary API exactly. A spatial owner
    /// must still evolve cell volumes consistently and apply each shared face
    /// flux with opposite signs; this leaf supplies neither a moving mesh nor
    /// a discrete geometric-conservation or temporal-accuracy guarantee.
    pub fn moving_hllc_normal_flux(
        &self,
        left: NormalConserved,
        right: NormalConserved,
        face_velocity: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NormalFaceFlux, Error> {
        let pl = self.normal_primitive(left, &mut cancelled)?;
        let pr = self.normal_primitive(right, &mut cancelled)?;
        self.moving_trace_interface(left, pl, right, pr, face_velocity, &mut cancelled)
            .map(NormalFaceFlux::from)
    }

    /// Source-owned reconstructed-trace path to [`Self::moving_hllc_normal_flux`].
    /// Laboratory and face-frame energy must both retain the declared caloric
    /// precision. The original source computes every thermodynamic property.
    pub fn moving_hllc_normal_flux_from_primitive(
        &self,
        left: NormalPrimitive,
        right: NormalPrimitive,
        face_velocity: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NormalFaceFlux, Error> {
        poll(&mut cancelled)?;
        let (left, pl) = self.normal_state_at(left)?;
        poll(&mut cancelled)?;
        let (right, pr) = self.normal_state_at(right)?;
        self.moving_trace_interface(left, pl, right, pr, face_velocity, &mut cancelled)
            .map(NormalFaceFlux::from)
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

    /// Moving-face ALE flux at an adiabatic, impermeable slip wall.
    ///
    /// `wall_velocity` is along the supplied normal axis [m/s], and
    /// `outward_normal` selects the domain side (-1 or +1). The original
    /// caloric shock/rarefaction solve uses gas velocity relative to the wall.
    /// In this axis the flux is `(0, p_star, 0, p_star * wall_velocity)`;
    /// multiply by `outward_normal` for outward flux. The energy is mechanical
    /// pressure work in the laboratory frame, with no advected mass or heat.
    /// Tangential wall motion has no effect on this inviscid slip condition.
    ///
    /// `signal_speed` bounds waves RELATIVE TO THE MOVING FACE. A moving-volume
    /// integrator must supply its geometry, geometric conservation law and
    /// time quadrature; this flux must not be inserted into a stationary cell
    /// update as if its volume had remained fixed. Zero wall speed reproduces
    /// the stationary API. Unsupported caloric paths are refused unchanged.
    pub fn moving_wall_normal_flux(
        &self,
        state: NormalConserved,
        outward_normal: f64,
        wall_velocity: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NormalFaceFlux, Error> {
        let primitive = self.normal_primitive(state, &mut cancelled)?;
        self.trace_moving_wall_interface(
            state,
            primitive,
            outward_normal,
            wall_velocity,
            &mut cancelled,
        )
        .map(NormalFaceFlux::from)
    }

    /// Source-admitted reconstructed trace for [`Self::moving_wall_normal_flux`].
    /// Avoids caloric inversion while retaining full kinetic-energy admission.
    pub fn moving_wall_normal_flux_from_primitive(
        &self,
        state: NormalPrimitive,
        outward_normal: f64,
        wall_velocity: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NormalFaceFlux, Error> {
        poll(&mut cancelled)?;
        let (state, primitive) = self.normal_state_at(state)?;
        self.trace_moving_wall_interface(
            state,
            primitive,
            outward_normal,
            wall_velocity,
            &mut cancelled,
        )
        .map(NormalFaceFlux::from)
    }

    /// Pressure-only outlet retaining the outgoing tangential velocity and energy.
    /// See [`Self::pressure_outlet_flux`] for wave sampling and backflow policy.
    pub fn pressure_outlet_normal_flux(
        &self,
        state: NormalConserved,
        outward_normal: f64,
        backpressure: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NormalFaceFlux, Error> {
        let primitive = self.normal_primitive(state, &mut cancelled)?;
        self.trace_pressure_outlet(
            state,
            primitive,
            outward_normal,
            backpressure,
            &mut cancelled,
        )
        .map(NormalFaceFlux::from)
    }

    /// Source-admitted reconstructed trace for [`Self::pressure_outlet_normal_flux`].
    pub fn pressure_outlet_normal_flux_from_primitive(
        &self,
        state: NormalPrimitive,
        outward_normal: f64,
        backpressure: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NormalFaceFlux, Error> {
        poll(&mut cancelled)?;
        let (state, primitive) = self.normal_state_at(state)?;
        self.trace_pressure_outlet(
            state,
            primitive,
            outward_normal,
            backpressure,
            &mut cancelled,
        )
        .map(NormalFaceFlux::from)
    }

    /// Incoming acoustic characteristic retaining outgoing planar shear.
    ///
    /// Reference velocity is signed along the supplied normal axis. The
    /// original interior isentrope and caloric source determine the face
    /// pressure; this is the local one-dimensional smooth-wave condition of
    /// [`Self::characteristic_outlet_flux`], not an oblique-wave or shock
    /// nonreflection guarantee. Backflow and unsupported source paths refuse.
    pub fn characteristic_outlet_normal_flux(
        &self,
        state: NormalConserved,
        outward_normal: f64,
        reference_pressure: f64,
        reference_velocity: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NormalFaceFlux, Error> {
        let primitive = self.normal_primitive(state, &mut cancelled)?;
        self.trace_characteristic_outlet(
            state,
            primitive,
            outward_normal,
            reference_pressure,
            reference_velocity,
            &mut cancelled,
        )
        .map(NormalFaceFlux::from)
    }

    /// Source-admitted reconstructed path to [`Self::characteristic_outlet_normal_flux`].
    /// Retains complete kinetic-energy precision without a second caloric inversion.
    pub fn characteristic_outlet_normal_flux_from_primitive(
        &self,
        state: NormalPrimitive,
        outward_normal: f64,
        reference_pressure: f64,
        reference_velocity: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NormalFaceFlux, Error> {
        poll(&mut cancelled)?;
        let (state, primitive) = self.normal_state_at(state)?;
        self.trace_characteristic_outlet(
            state,
            primitive,
            outward_normal,
            reference_pressure,
            reference_velocity,
            &mut cancelled,
        )
        .map(NormalFaceFlux::from)
    }

    fn trace_pressure_outlet<const N: usize, Q: Trace<N>>(
        &self,
        state: Q,
        primitive: Primitive,
        outward_normal: f64,
        backpressure: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NumericalFlux<N>, Error> {
        let (face, signal) =
            self.pressure_outlet_state(primitive, outward_normal, backpressure, cancelled)?;
        self.trace_outgoing_flux(state, face, signal, cancelled)
    }

    fn trace_characteristic_outlet<const N: usize, Q: Trace<N>>(
        &self,
        state: Q,
        primitive: Primitive,
        outward_normal: f64,
        reference_pressure: f64,
        reference_velocity: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NumericalFlux<N>, Error> {
        let (face, signal) = self.characteristic_outlet_state(
            primitive,
            outward_normal,
            reference_pressure,
            reference_velocity,
            cancelled,
        )?;
        self.trace_outgoing_flux(state, face, signal, cancelled)
    }

    /// Shared physical transport after either pressure boundary has admitted
    /// its axial face state. Every outgoing tangent retains the same velocity.
    fn trace_outgoing_flux<const N: usize, Q: Trace<N>>(
        &self,
        state: Q,
        face: Primitive,
        signal: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NumericalFlux<N>, Error> {
        let normal = self.conserved(face.density, face.velocity, face.pressure)?;
        let old = state.array();
        let mut values = [0.0; N];
        values[0] = normal.density;
        values[1] = normal.momentum;
        values[N - 1] = normal.energy;
        // The outgoing contact transports each tangent without inventing an
        // exterior shear state. Its kinetic energy crosses the same face once.
        for k in 2..N - 1 {
            let velocity = old[k] / old[0];
            values[k] = face.density * velocity;
            values[N - 1] += 0.5 * values[k] * velocity;
        }
        let (state, face) = self.admit_trace_precision(
            Q::from_array(values),
            face,
            "pressure outlet cannot retain the declared caloric precision",
        )?;
        poll(cancelled)?;
        flux(physical_flux(state, face), signal, false)
    }

    fn normal_wall_interface(
        &self,
        state: NormalConserved,
        primitive: Primitive,
        outward_normal: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NormalFaceFlux, Error> {
        self.trace_wall_interface(state, primitive, outward_normal, cancelled)
            .map(NormalFaceFlux::from)
    }

    fn admissible_trace<const N: usize>(&self, state: impl Trace<N>) -> bool {
        state.array().iter().all(|x| x.is_finite())
            && positive(state.array()[0])
            && self.admissible_state(state.normal_state())
    }

    fn admit_trace_precision<const N: usize, Q: Trace<N>>(
        &self,
        state: Q,
        primitive: Primitive,
        reason: &'static str,
    ) -> Result<(Q, Primitive), Error> {
        let reduced = state.normal_state();
        let recovered_velocity = reduced.momentum / reduced.density;
        let recovered_energy =
            reduced.energy / reduced.density - 0.5 * recovered_velocity * recovered_velocity;
        let expected_energy = self.sensible(primitive.temperature)?;
        let domain = self.domain();
        let tolerance = domain.absolute_energy_tolerance
            + domain.relative_energy_tolerance * expected_energy.abs();
        if !self.admissible_trace(state) || (recovered_energy - expected_energy).abs() > tolerance {
            return Err(Error::Invalid(reason));
        }
        Ok((state, primitive))
    }

    fn decode_trace<const N: usize>(
        &self,
        state: impl Trace<N>,
        reason: &'static str,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        poll(cancelled)?;
        let values = state.array();
        if !values.iter().all(|x| x.is_finite()) || !positive(values[0]) {
            return Err(Error::Invalid(reason));
        }
        self.decode(state.normal_state(), None, cancelled)
    }

    fn trace_wall_interface<const N: usize>(
        &self,
        state: impl Trace<N>,
        primitive: Primitive,
        outward_normal: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NumericalFlux<N>, Error> {
        self.trace_moving_wall_interface(state, primitive, outward_normal, 0.0, cancelled)
    }

    fn trace_moving_wall_interface<const N: usize>(
        &self,
        state: impl Trace<N>,
        primitive: Primitive,
        outward_normal: f64,
        wall_velocity: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NumericalFlux<N>, Error> {
        let wall = self.moving_wall_flux(
            state.normal_state(),
            primitive,
            outward_normal,
            wall_velocity,
            cancelled,
        )?;
        let mut values = [0.0; N];
        values[1] = wall.momentum;
        values[N - 1] = wall.energy;
        flux(values, wall.signal_speed, false)
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

    /// F - w*U through a translating face, in laboratory conserved variables.
    /// Reconstruct source-owned caloric energy in the face frame rather than
    /// subtracting two potentially large laboratory kinetic energies.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn moving_interface(
        &self,
        left: Conserved,
        pl: Primitive,
        right: Conserved,
        pr: Primitive,
        velocity: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        self.moving_trace_interface(
            NormalConserved::from(left),
            pl,
            NormalConserved::from(right),
            pr,
            velocity,
            cancelled,
        )
        .map(|f| FaceFlux::from(NormalFaceFlux::from(f)))
    }

    // One Galilean transformation for axial, planar and spatial transport.
    // Solve in the face frame, then transform the numerical momentum and
    // energy FLUX, including pressure work, back to laboratory coordinates.
    // Springel (2010), section 3.4, arXiv:0901.4107.
    #[allow(clippy::too_many_arguments)]
    fn moving_trace_interface<const N: usize, Q: Trace<N>>(
        &self,
        left: Q,
        pl: Primitive,
        right: Q,
        pr: Primitive,
        velocity: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NumericalFlux<N>, Error> {
        if velocity == 0.0 {
            return self.trace_interface(left, pl, right, pr, cancelled);
        }
        if !velocity.is_finite() {
            return Err(Error::Invalid("nonfinite face velocity"));
        }
        let shift = |state: Q, primitive: Primitive| {
            let (normal, primitive) = self.state_at(
                primitive.density,
                primitive.velocity - velocity,
                primitive.pressure,
            )?;
            let mut values = state.array();
            values[1] = normal.momentum;
            values[N - 1] = normal.energy;
            let mut has_shear = false;
            for k in 2..N - 1 {
                // Preserve the original tangential momentum; only the normal
                // frame changed. Rebuild heat from the source rather than
                // subtracting large laboratory kinetic-energy terms.
                has_shear |= values[k] != 0.0;
                values[N - 1] += 0.5 * values[k] * (values[k] / values[0]);
            }
            let state = Q::from_array(values);
            if has_shear {
                self.admit_trace_precision(
                    state,
                    primitive,
                    "moving-face state cannot retain the declared caloric precision",
                )
            } else {
                // state_at already admitted the complete axial energy. Keep
                // the original one-dimensional source work and arithmetic.
                Ok((state, primitive))
            }
        };
        poll(cancelled)?;
        let (left, pl) = shift(left, pl)?;
        poll(cancelled)?;
        let (right, pr) = shift(right, pr)?;
        let relative = self.trace_interface(left, pl, right, pr, cancelled)?;
        let mut values = relative.values;
        values[1] = relative.values[1] + velocity * relative.values[0];
        values[N - 1] = relative.values[N - 1]
            + velocity * relative.values[1]
            + (0.5 * velocity * velocity) * relative.values[0];
        poll(cancelled)?;
        flux(values, relative.signal_speed, relative.used_hll_fallback)
    }

    fn normal_interface(
        &self,
        left: NormalConserved,
        pl: Primitive,
        right: NormalConserved,
        pr: Primitive,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NormalFaceFlux, Error> {
        self.trace_interface(left, pl, right, pr, cancelled)
            .map(NormalFaceFlux::from)
    }

    fn trace_interface<const N: usize, Q: Trace<N>>(
        &self,
        left: Q,
        pl: Primitive,
        right: Q,
        pr: Primitive,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NumericalFlux<N>, Error> {
        checked_trace_interface(
            left,
            pl,
            right,
            pr,
            |q, _| Ok(self.admissible_trace(q)),
            cancelled,
        )
    }
}

/// The original HLLC/HLL algebra with a constitutive owner's checked admission.
/// Every supplied trace/primitive pair must already describe one accepted EOS
/// state in one common energy reference. The callback resolves trial states;
/// it must distinguish cancellation from a recoverable star-state refusal.
#[cfg(feature = "material")]
pub(super) fn checked_euler_flux<E, C>(
    left: Conserved,
    pl: Primitive,
    right: Conserved,
    pr: Primitive,
    mut admit: impl FnMut(Conserved, &mut C) -> Result<bool, E>,
    cancelled: &mut C,
) -> Result<FaceFlux, E>
where
    E: From<Error>,
    C: FnMut() -> bool,
{
    checked_trace_interface(
        NormalConserved::from(left),
        pl,
        NormalConserved::from(right),
        pr,
        |q, cancelled| admit(q.normal_state(), cancelled),
        cancelled,
    )
    .map(|f| FaceFlux::from(NormalFaceFlux::from(f)))
}

fn checked_trace_interface<const N: usize, Q: Trace<N>, E, C>(
    left: Q,
    pl: Primitive,
    right: Q,
    pr: Primitive,
    mut admit: impl FnMut(Q, &mut C) -> Result<bool, E>,
    cancelled: &mut C,
) -> Result<NumericalFlux<N>, E>
where
    E: From<Error>,
    C: FnMut() -> bool,
{
    poll(cancelled).map_err(E::from)?;
    let fl = physical_flux(left, pl);
    let fr = physical_flux(right, pr);
    let sl = (pl.velocity - pl.sound_speed).min(pr.velocity - pr.sound_speed);
    let sr = (pl.velocity + pl.sound_speed).max(pr.velocity + pr.sound_speed);
    let speed = sl.abs().max(sr.abs());
    if !positive(speed) || !fl.iter().chain(&fr).all(|x| x.is_finite()) {
        return Err(Error::Invalid("unrepresentable Euler flux or signal speed").into());
    }
    let (ql, qr) = (left.array(), right.array());
    // Consistency is exact for equal traces. Forming a redundant HLLC
    // star state can round (rho*s)/s away from rho and manufacture flux
    // in a resting atmosphere or a uniform moving flow.
    if ql == qr && pl == pr {
        return flux(fl, speed, false);
    }
    // A stationary pressure contact (including a tangential shear jump)
    // transports no mass, tangential momentum or energy in inviscid Euler.
    if pl.velocity == 0.0 && pr.velocity == 0.0 && pl.pressure == pr.pressure {
        return flux(
            core::array::from_fn(|k| if k == 1 { pl.pressure } else { 0.0 }),
            speed,
            false,
        );
    }
    // An exact moving pressure contact transports the upwind trace,
    // including its tangential kinetic energy. Reconstructing unchanged
    // star traces can round heat outside the admitted source interval and
    // spuriously select HLL. Both original traces are already admitted.
    if pl.velocity == pr.velocity && pl.pressure == pr.pressure {
        return flux(if pl.velocity > 0.0 { fl } else { fr }, speed, false);
    }
    if sl >= 0.0 {
        return flux(fl, speed, false);
    }
    if sr <= 0.0 {
        return flux(fr, speed, false);
    }
    let al = ql[0] * (sl - pl.velocity);
    let ar = qr[0] * (sr - pr.velocity);
    let contact = (pr.pressure - pl.pressure + al * pl.velocity - ar * pr.velocity) / (al - ar);
    let pstar = pl.pressure + al * (contact - pl.velocity);
    if contact.is_finite() && sl < contact && contact < sr && positive(pstar) {
        let star = |q: Q, p: Primitive, s: f64| {
            let values = q.array();
            Q::from_array(core::array::from_fn(|k| match k {
                0 => values[0] * (s - p.velocity) / (s - contact),
                1 => values[0] * (s - p.velocity) / (s - contact) * contact,
                k if k == N - 1 => {
                    ((s - p.velocity) * values[k] - p.pressure * p.velocity + pstar * contact)
                        / (s - contact)
                }
                k => values[k] * (s - p.velocity) / (s - contact),
            }))
        };
        let ls = star(left, pl, sl);
        let rs = star(right, pr, sr);
        if admit(ls, cancelled)? && admit(rs, cancelled)? {
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
    for expansion in 0..=8 {
        poll(cancelled).map_err(E::from)?;
        let (l, r) = if expansion == 0 {
            (sl, sr)
        } else {
            let a = speed * (1_u32 << (expansion - 1)) as f64;
            (-a, a)
        };
        if !l.is_finite() || !r.is_finite() {
            return Err(Error::Invalid("HLL wave-speed overflow").into());
        }
        let middle = Q::from_array(core::array::from_fn(|k| {
            (r * qr[k] - l * ql[k] + fl[k] - fr[k]) / (r - l)
        }));
        if admit(middle, cancelled)? {
            return flux(
                core::array::from_fn(|k| {
                    (r * fl[k] - l * fr[k] + l * r * (qr[k] - ql[k])) / (r - l)
                }),
                l.abs().max(r.abs()),
                true,
            );
        }
    }
    Err(Error::Invalid("HLL intermediate state outside admitted source domain").into())
}

fn physical_flux<const N: usize>(q: impl Trace<N>, p: Primitive) -> [f64; N] {
    let q = q.array();
    core::array::from_fn(|k| match k {
        0 => q[1],
        1 => q[1] * p.velocity + p.pressure,
        k if k == N - 1 => (q[k] + p.pressure) * p.velocity,
        k => q[k] * p.velocity,
    })
}
fn flux<const N: usize, E: From<Error>>(
    f: [f64; N],
    signal_speed: f64,
    used_hll_fallback: bool,
) -> Result<NumericalFlux<N>, E> {
    if !f.iter().all(|x| x.is_finite()) {
        return Err(Error::Invalid("non-finite numerical flux").into());
    }
    Ok(NumericalFlux {
        values: f,
        signal_speed,
        used_hll_fallback,
    })
}
