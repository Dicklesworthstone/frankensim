//! Acoustic outlet over the incident trace's frozen local composition.
//!
//! Only one boundary solve freezes composition. The next reconstruction and
//! every RK stage borrow their own current partial densities and original
//! caloric sources; no surrogate or cached fixed-mixture EOS is constructed.

use super::{LocalCaloric, characteristic_error};
use crate::dynamics::characteristic::acoustic::outlet;
use crate::multicomponent::{Cell, Error, FaceFlux, MixtureGas, Primitive, poll};

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
        poll(cancelled)?;
        if out.len() != self.species_count() + 2 || !matches!(outward_normal, -1.0 | 1.0) {
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
        self.outgoing_characteristic_flux(&law, face, [], outward_normal, signal, out, cancelled)
    }
}
