//! A rigid-backed porous lining in the existing local wall-admittance path.
//!
//! The material owner supplies JCA(L) impedance with exp(+i omega t); the
//! duct uses exp(-i omega t). Preparation evaluates the original source once
//! per frequency and conjugates once. Local radius and surface slope remain
//! the duct geometry's responsibility, so every quadrature point consumes the
//! same specific impedance without repeating the constitutive evaluation.
//!
//! This is a planar, locally reacting, rigid-frame boundary approximation.
//! Tangential propagation inside the lining, cylindrical-layer curvature,
//! moving-frame Biot waves and mean flow are not resolved. Local reaction must
//! be justified for the application; a normal-incidence material spectrum is
//! not evidence of an arbitrary duct lining's angular response. See COMSOL's
//! local-versus-extended-reaction example:
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.models.aco.porous_absorber_time_domain/porous_absorber_time_domain.html>.

use super::{WallLaw, WallModel};
use crate::compose::{ComposeOptions, PreparedChain};
use crate::fluid::{self, FluidDuctError};
use crate::{ConeSlices, Duct, DuctError};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_material::porous::{PorousFluidLaw, PorousFluidResponse};
use fs_math::c64::C64;

/// One original porous material, coherent saturating fluid, thickness and
/// frequency, prepared for distributed local wall absorption.
///
/// The source law and fluid sample are retained. Its source frequency support
/// is enforced by the material owner, and the result cannot be reused at a
/// different frequency. Prepare a new value for each frequency in a sweep.
/// Geometry radii describe the remaining fluid core, after the lining is
/// installed; the adapter does not silently subtract a lining thickness.
#[derive(Debug, Clone)]
pub struct PorousLining<'a> {
    material: &'a PorousFluidLaw,
    response: PorousFluidResponse,
    thickness_m: f64,
    impedance: C64,
}

impl<'a> PorousLining<'a> {
    /// Prepare the original JCA(L) law once, retaining its material refusal.
    /// The reference viscosity remains the one on the original porous law;
    /// no room-air or material-name default is introduced.
    pub fn try_new<M: ThermoviscousFluid + ?Sized>(
        material: &'a PorousFluidLaw,
        fluid: &M,
        omega: f64,
        thickness_m: f64,
    ) -> Result<Self, MaterialError> {
        if !thickness_m.is_finite() || thickness_m <= 0.0 {
            return Err(MaterialError::Parameters {
                what: "porous lining thickness must be positive and finite".into(),
            });
        }
        let response = material.effective_fluid(fluid, omega)?;
        let impedance = response.surface_impedance_rigid(thickness_m)?.conj();
        if impedance == C64::ZERO {
            return Err(MaterialError::State {
                what: "porous lining requires a nonzero surface impedance".into(),
            });
        }
        Ok(Self {
            material,
            response,
            thickness_m,
            impedance,
        })
    }

    /// Original source, characterization viscosity and frequency support.
    pub const fn material(&self) -> &PorousFluidLaw {
        self.material
    }

    /// Exactly the coherent fluid and JCA(L) response used in preparation.
    pub const fn response(&self) -> &PorousFluidResponse {
        &self.response
    }

    /// Physical layer thickness [m], independent of the fluid-core radius.
    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }

    /// Cached specific impedance [Pa s/m] in the duct's exp(-i omega t)
    /// convention. Multiplication by surface area belongs to the duct owner.
    pub const fn surface_impedance_pa_s_m(&self) -> C64 {
        self.impedance
    }

    /// Distributed absorption in the existing cylindrical, conical and
    /// smooth-profile plane-wave chain, using the same fluid as the pores.
    ///
    /// This explicitly uses an inviscid core and the porous wall admittance.
    /// The rigid-wall no-slip/isothermal diffusion law is not also imposed on
    /// this permeable surface. For independently justified additional physics,
    /// callers may bind this value through the general WallModel::Law seam.
    /// The original cone/profile resolution, budgets and cancellation apply.
    pub fn duct_chain(
        &self,
        duct: &Duct,
        cone_slices: ConeSlices,
        options: ComposeOptions,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedChain, FluidDuctError> {
        fluid::duct_chain(
            duct,
            self.response.fluid(),
            self.response.angular_frequency_rad_s(),
            WallModel::Law(self),
            cone_slices,
            options,
            cancelled,
        )
    }
}

impl WallLaw for PorousLining<'_> {
    fn specific_impedance(&self, radius_m: f64, omega: f64) -> Result<C64, DuctError> {
        if !radius_m.is_finite() || radius_m <= 0.0 {
            return Err(DuctError::BadParameter {
                what: "porous lining needs a positive finite fluid-core radius",
            });
        }
        if omega != self.response.angular_frequency_rad_s() {
            return Err(DuctError::BadParameter {
                what: "porous lining was prepared at a different frequency",
            });
        }
        Ok(self.impedance)
    }
}

#[cfg(test)]
mod tests;
