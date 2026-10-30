//! Three-dimensional traces of the shared thermally perfect Riemann closure.

use super::{Conserved, Error, FrozenGas, NumericalFlux, Primitive, Trace, poll};

/// Conservative trace in a three-dimensional face's orthonormal frame.
///
/// The spatial owner supplies one normal and two orthogonal tangents. Total
/// energy includes all three kinetic components in the owner's fixed sensible
/// energy gauge; both shear components participate in HLL mixing admission.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SpatialConserved {
    /// Mass density [kg/m³].
    pub density: f64,
    /// Momentum density along the oriented face normal [kg/(m² s)].
    pub normal_momentum: f64,
    /// Momentum densities along the two unit tangents [kg/(m² s)].
    pub transverse_momentum: [f64; 2],
    /// Referenced total energy density, including all kinetic energy [J/m³].
    pub energy: f64,
}

impl Trace<5> for SpatialConserved {
    fn array(self) -> [f64; 5] {
        [
            self.density,
            self.normal_momentum,
            self.transverse_momentum[0],
            self.transverse_momentum[1],
            self.energy,
        ]
    }

    fn from_array(q: [f64; 5]) -> Self {
        Self {
            density: q[0],
            normal_momentum: q[1],
            transverse_momentum: [q[2], q[3]],
            energy: q[4],
        }
    }

    fn normal_state(self) -> Conserved {
        let [v, w] = self.transverse_momentum;
        let transverse_kinetic = 0.5 * v * (v / self.density) + 0.5 * w * (w / self.density);
        Conserved {
            density: self.density,
            momentum: self.normal_momentum,
            energy: self.energy - transverse_kinetic,
        }
    }
}

/// Reconstructed volumetric trace, without caller-supplied caloric properties.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SpatialPrimitive {
    /// Mass density [kg/m³].
    pub density: f64,
    /// Velocity along the oriented face normal [m/s].
    pub normal_velocity: f64,
    /// Velocities along the two unit tangents [m/s].
    pub transverse_velocity: [f64; 2],
    /// Absolute pressure [Pa].
    pub pressure: f64,
}

/// Shared three-dimensional face flux from the planar/volumetric HLLC kernel.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SpatialFaceFlux {
    /// Mass flux [kg/(m² s)].
    pub mass: f64,
    /// Normal momentum flux including pressure [Pa].
    pub normal_momentum: f64,
    /// Tangential momentum fluxes [Pa].
    pub transverse_momentum: [f64; 2],
    /// Referenced total-energy flux, including all kinetic energy [W/m²].
    pub energy: f64,
    /// Largest magnitude of the actual numerical wave speeds [m/s].
    pub signal_speed: f64,
    /// Whether the shared kernel used its admitted HLL fallback.
    pub used_hll_fallback: bool,
}

impl From<NumericalFlux<5>> for SpatialFaceFlux {
    fn from(f: NumericalFlux<5>) -> Self {
        Self {
            mass: f.values[0],
            normal_momentum: f.values[1],
            transverse_momentum: [f.values[2], f.values[3]],
            energy: f.values[4],
            signal_speed: f.signal_speed,
            used_hll_fallback: f.used_hll_fallback,
        }
    }
}

impl FrozenGas {
    /// Construct a three-dimensional trace using the original caloric source.
    /// Refuses a full kinetic-energy scale that erases admitted heat precision.
    pub fn spatial_conserved(
        &self,
        density: f64,
        normal_velocity: f64,
        transverse_velocity: [f64; 2],
        pressure: f64,
    ) -> Result<SpatialConserved, Error> {
        self.spatial_state_at(SpatialPrimitive {
            density,
            normal_velocity,
            transverse_velocity,
            pressure,
        })
        .map(|(state, _)| state)
    }

    fn spatial_state_at(
        &self,
        state: SpatialPrimitive,
    ) -> Result<(SpatialConserved, Primitive), Error> {
        let SpatialPrimitive {
            density,
            normal_velocity,
            transverse_velocity: [v, w],
            pressure,
        } = state;
        if !v.is_finite() || !w.is_finite() {
            return Err(Error::Invalid("non-finite transverse velocity"));
        }
        let (normal, primitive) = self.state_at(density, normal_velocity, pressure)?;
        let q = SpatialConserved {
            density,
            normal_momentum: normal.momentum,
            transverse_momentum: [density * v, density * w],
            energy: normal.energy + density * (0.5 * v * v + 0.5 * w * w),
        };
        self.admit_trace_precision(
            q,
            primitive,
            "spatial conservative state cannot retain the declared caloric precision",
        )
    }

    /// Decode full three-dimensional kinetic energy through the original source.
    /// `Primitive::velocity` is normal; each tangent velocity is momentum/rho.
    pub fn spatial_primitive(
        &self,
        state: SpatialConserved,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        self.decode_trace(state, "invalid spatial conservative state", &mut cancelled)
    }

    /// Three-dimensional HLLC/HLL flux sharing the existing planar kernel.
    /// Both shear momenta and full energy use the same Riemann waves; HLL star
    /// admission includes conversion of both transverse kinetic energies to heat.
    pub fn hllc_spatial_flux(
        &self,
        left: SpatialConserved,
        right: SpatialConserved,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialFaceFlux, Error> {
        let pl = self.spatial_primitive(left, &mut cancelled)?;
        let pr = self.spatial_primitive(right, &mut cancelled)?;
        self.trace_interface(left, pl, right, pr, &mut cancelled)
            .map(SpatialFaceFlux::from)
    }

    /// Source-computed primitive fast path; retains source/precision admission
    /// while avoiding inversion of just-reconstructed conserved energies.
    pub fn hllc_spatial_flux_from_primitive(
        &self,
        left: SpatialPrimitive,
        right: SpatialPrimitive,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialFaceFlux, Error> {
        poll(&mut cancelled)?;
        let (left, pl) = self.spatial_state_at(left)?;
        poll(&mut cancelled)?;
        let (right, pr) = self.spatial_state_at(right)?;
        self.trace_interface(left, pl, right, pr, &mut cancelled)
            .map(SpatialFaceFlux::from)
    }

    /// Exact thermally perfect shock/rarefaction at a stationary adiabatic slip
    /// wall. `outward_normal` is -1 or +1 along the supplied normal axis. There
    /// is no mass, tangential momentum or energy transfer through this wall.
    pub fn stationary_wall_spatial_flux(
        &self,
        state: SpatialConserved,
        outward_normal: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialFaceFlux, Error> {
        let primitive = self.spatial_primitive(state, &mut cancelled)?;
        self.trace_wall_interface(state, primitive, outward_normal, &mut cancelled)
            .map(SpatialFaceFlux::from)
    }

    /// Source-admitted primitive fast path to the same stationary wall closure.
    pub fn stationary_wall_spatial_flux_from_primitive(
        &self,
        state: SpatialPrimitive,
        outward_normal: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialFaceFlux, Error> {
        poll(&mut cancelled)?;
        let (state, primitive) = self.spatial_state_at(state)?;
        self.trace_wall_interface(state, primitive, outward_normal, &mut cancelled)
            .map(SpatialFaceFlux::from)
    }
}
