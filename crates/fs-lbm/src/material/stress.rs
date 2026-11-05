//! Physical viscous loads from the native owner's accepted velocity gradient.
//!
//! This is a local constitutive observation on the original boundary-aware
//! stencil, not another collision, wall-force estimator or heating source.

use super::{MaterialLbmError, MaterialScaling};
use crate::kinematics::FlowKinematics3;
use fs_material::fluid_transport::{
    FluidTransport,
    rheology::ThermalCarreauYasuda,
    viscous::{DeformationRate, NewtonianViscosity, ThermalViscousResponse, ViscousResponse},
};

impl<F: FluidTransport + ?Sized> MaterialScaling<'_, F> {
    /// Resolve SI viscous stress, surface traction and physical dissipation
    /// from native LATTICE-UNIT kinematics and this retained material snapshot.
    ///
    /// Pass a sample from the owner's accepted `kinematics` query without
    /// calling `to_physical` first. `grad(u)_SI = grad(u)_lattice / dt`; dx
    /// cancels. Moving walls, sparse occupancy and the current half-force
    /// velocity remain the native gradient owner's responsibility. The same
    /// original transport snapshot is used on every call, with no source
    /// re-evaluation and no replacement constants. An evolving material must
    /// use a scaling bound to its current snapshot, or the explicit thermal
    /// rheology method below.
    ///
    /// `bulk_viscosity_pa_s` is independently supplied, including explicit
    /// zero. The response excludes thermodynamic pressure and modeled SGS
    /// stress. Its tensor order is xx/yy/zz/xy/yz/zx, with tensor shear.
    /// This continuum constitutive observation is distinct from native
    /// collision stress and link momentum exchange. It does not advance flow,
    /// apply a traction or insert heat, and adds no accuracy certification.
    pub fn viscous_stress_from_lattice_kinematics(
        &self,
        kinematics: FlowKinematics3,
        bulk_viscosity_pa_s: f64,
    ) -> Result<ViscousResponse, MaterialLbmError> {
        let rate = self.physical_deformation_rate(kinematics)?;
        Ok(
            NewtonianViscosity::from_transport(&self.transport, bulk_viscosity_pa_s)?
                .evaluate(rate)?,
        )
    }

    /// Evaluate an explicitly authored thermal shear law at the actual local
    /// temperature and resolved deviatoric strain rate, in physical SI units.
    ///
    /// Kinematics must still be in lattice units. The original temperature
    /// shift, shear support and differential viscosity remain in the returned
    /// material response. In particular this evaluates the law at the current
    /// rate instead of applying a viscosity cached at an earlier scalar rate.
    /// Isotropic dilation is handled only by the explicit bulk viscosity; it
    /// does not cause shear thinning. The same no-mutation/no-added-heat
    /// qualification as [`Self::viscous_stress_from_lattice_kinematics`] applies.
    pub fn thermal_viscous_stress_from_lattice_kinematics(
        &self,
        kinematics: FlowKinematics3,
        law: &ThermalCarreauYasuda,
        temperature_k: f64,
        bulk_viscosity_pa_s: f64,
    ) -> Result<ThermalViscousResponse, MaterialLbmError> {
        let rate = self.physical_deformation_rate(kinematics)?;
        Ok(law.viscous_response(temperature_k, rate, bulk_viscosity_pa_s)?)
    }

    fn physical_deformation_rate(
        &self,
        kinematics: FlowKinematics3,
    ) -> Result<DeformationRate, MaterialLbmError> {
        let lattice = kinematics.gradient();
        let gradient = lattice.map(|row| row.map(|x| x / self.timestep_s));
        if gradient
            .iter()
            .flatten()
            .zip(lattice.iter().flatten())
            .any(|(&si, &lu)| !si.is_finite() || (lu != 0.0 && si == 0.0))
        {
            return Err(MaterialLbmError::Invalid(
                "physical velocity gradient is unrepresentable",
            ));
        }
        Ok(DeformationRate::from_velocity_gradient(gradient)?)
    }
}
