//! Acoustic outlet over the incident trace's frozen local composition.
//!
//! Only one boundary solve freezes composition. The next reconstruction and
//! every RK stage borrow their own current partial densities and original
//! caloric sources; no surrogate or cached fixed-mixture EOS is constructed.

use super::{LocalCaloric, characteristic_error};
use crate::dynamics::characteristic::acoustic::outlet;
use crate::multicomponent::{
    Cell, Error, FaceFlux, MAX_NORMAL_COMPONENTS, MAX_SPATIAL_COMPONENTS, MixtureGas, NormalCell,
    NormalPrimitive, Primitive, SpatialCell, SpatialPrimitive, poll,
};

impl MixtureGas {
    /// Incoming acoustic characteristic relative to a reference absolute
    /// pressure [Pa] and signed axial velocity [m/s].
    ///
    /// The incident trace supplies the outgoing characteristic, composition
    /// and entropy. The reference pressure is continued on that local
    /// isentrope using the original species Cp(T), so an outgoing entropy or
    /// composition contact does not acquire an artificial pressure reflection.
    /// Each species leaves with its current mass fraction; momentum and total
    /// energy use the same face state and original per-species energy datums.
    ///
    /// Both velocities are along increasing duct coordinate. The outward
    /// normal is -1 at a left end, +1 at a right end. Supersonic outflow retains
    /// the interior flux; an expansion crossing sonic speed uses its sonic
    /// state. Backflow refuses because this outlet supplies no incoming
    /// composition or thermodynamic state. Unsupported source paths refuse.
    ///
    /// This is a local one-dimensional smooth-wave closure, not a guarantee
    /// of nonreflection for shocks, composition/entropy-gradient interactions
    /// or transverse modes. Output has `species_count()+2` entries: species
    /// mass fluxes, axial momentum flux and referenced total-energy flux.
    /// Every refusal or cancellation leaves output unchanged. Successful
    /// evaluation allocates no heap workspace or replacement caloric owner.
    #[allow(clippy::too_many_arguments)]
    pub fn characteristic_outlet_flux(
        &self,
        interior: Cell<'_>,
        outward_normal: f64,
        reference_pressure: f64,
        reference_velocity: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let primitive = self.primitive(interior, &mut cancelled)?;
        self.characteristic_outlet_admitted_flux(
            interior.partial_densities,
            primitive,
            outward_normal,
            reference_pressure,
            reference_velocity,
            out,
            &mut cancelled,
        )
    }

    /// Planar incoming acoustic characteristic carrying the outgoing species,
    /// entropy and transverse velocity. Output order is species, normal
    /// momentum, transverse momentum and full referenced energy.
    /// Reference velocity is signed along the supplied normal axis; the
    /// smooth normal-wave and backflow policy is that of
    /// [`Self::characteristic_outlet_flux`]. Output is atomic on refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn characteristic_outlet_normal_flux(
        &self,
        interior: NormalCell<'_>,
        outward_normal: f64,
        reference_pressure: f64,
        reference_velocity: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let primitive = self.normal_primitive(interior, &mut cancelled)?;
        self.characteristic_outlet_components(
            interior.partial_densities,
            primitive,
            [interior.transverse_momentum / primitive.density],
            outward_normal,
            reference_pressure,
            reference_velocity,
            out,
            &mut cancelled,
        )
    }

    /// Source-admitted reconstruction path to [`Self::characteristic_outlet_normal_flux`].
    /// Retains full kinetic-energy precision without another caloric inversion.
    #[allow(clippy::too_many_arguments)]
    pub fn characteristic_outlet_normal_flux_from_primitive(
        &self,
        interior: NormalPrimitive<'_>,
        outward_normal: f64,
        reference_pressure: f64,
        reference_velocity: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(&mut cancelled)?;
        let mut state = [0.0; MAX_NORMAL_COMPONENTS];
        let primitive =
            self.normal_conserved_into(interior, &mut state[..self.species_count() + 3])?;
        self.characteristic_outlet_components(
            interior.partial_densities,
            primitive,
            [interior.transverse_velocity],
            outward_normal,
            reference_pressure,
            reference_velocity,
            out,
            &mut cancelled,
        )
    }

    /// Three-dimensional incoming acoustic characteristic with both outgoing
    /// shear velocities and full kinetic energy. Output order is species,
    /// normal momentum, two transverse momenta and referenced energy. The
    /// local smooth-wave condition does not claim nonreflection of oblique
    /// waves or shocks. All source failures leave output untouched.
    #[allow(clippy::too_many_arguments)]
    pub fn characteristic_outlet_spatial_flux(
        &self,
        interior: SpatialCell<'_>,
        outward_normal: f64,
        reference_pressure: f64,
        reference_velocity: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let primitive = self.spatial_primitive(interior, &mut cancelled)?;
        self.characteristic_outlet_components(
            interior.partial_densities,
            primitive,
            interior.transverse_momentum.map(|m| m / primitive.density),
            outward_normal,
            reference_pressure,
            reference_velocity,
            out,
            &mut cancelled,
        )
    }

    /// Source-admitted reconstruction path to [`Self::characteristic_outlet_spatial_flux`].
    /// Preserves both tangential-energy admission checks and atomic output.
    #[allow(clippy::too_many_arguments)]
    pub fn characteristic_outlet_spatial_flux_from_primitive(
        &self,
        interior: SpatialPrimitive<'_>,
        outward_normal: f64,
        reference_pressure: f64,
        reference_velocity: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(&mut cancelled)?;
        let mut state = [0.0; MAX_SPATIAL_COMPONENTS];
        let primitive =
            self.spatial_conserved_into(interior, &mut state[..self.species_count() + 4])?;
        self.characteristic_outlet_components(
            interior.partial_densities,
            primitive,
            interior.transverse_velocity,
            outward_normal,
            reference_pressure,
            reference_velocity,
            out,
            &mut cancelled,
        )
    }

    /// Reconstructed interior state and canonical partial densities from the
    /// same admitted RK stage; avoids a redundant total-energy inversion.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::multicomponent) fn characteristic_outlet_admitted_flux(
        &self,
        partial_densities: &[f64],
        inside: Primitive,
        outward_normal: f64,
        reference_pressure: f64,
        reference_velocity: f64,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        self.characteristic_outlet_components(
            partial_densities,
            inside,
            [],
            outward_normal,
            reference_pressure,
            reference_velocity,
            out,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn characteristic_outlet_components<const TANGENTS: usize>(
        &self,
        partial_densities: &[f64],
        inside: Primitive,
        transverse_velocity: [f64; TANGENTS],
        outward_normal: f64,
        reference_pressure: f64,
        reference_velocity: f64,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(cancelled)?;
        if TANGENTS > 2
            || out.len() != self.species_count() + TANGENTS + 2
            || !matches!(outward_normal, -1.0 | 1.0)
            || !transverse_velocity.iter().all(|v| v.is_finite())
        {
            return Err(Error::Invalid(
                "mixture characteristic outlet dimension or normal",
            ));
        }
        let law = LocalCaloric::new(self, partial_densities)?;
        let (face, signal) = outlet(
            &law,
            Primitive {
                velocity: outward_normal * inside.velocity,
                ..inside
            },
            reference_pressure,
            outward_normal * reference_velocity,
            cancelled,
        )
        .map_err(characteristic_error)?;
        self.outgoing_characteristic_flux(
            &law,
            face,
            transverse_velocity,
            outward_normal,
            signal,
            out,
            cancelled,
        )
    }
}
