//! Explicit physical volumetric viscosity for the native D3Q19 liquid.

use super::{FreeSurface3, MaterialFreeSurface3, MaterialFreeSurfaceScaling, MaterialLbmError};
use crate::d3q19::CollisionModel3;

#[derive(Clone, Copy, Debug)]
pub(super) struct BulkViscosity {
    viscosity_pa_s: f64,
    relaxation_rate: f64,
    higher_order_rate: f64,
}

impl MaterialFreeSurfaceScaling<'_> {
    /// Author a uniform physical bulk viscosity [Pa s], independent of the
    /// liquid source's shear viscosity. Every 3-D constructor retains this
    /// coefficient through local shear/temperature updates and uses the same
    /// coefficient for `zeta (div u)^2` when molecular heating is requested.
    ///
    /// This explicitly supplied scalar is not inferred from numerical damping,
    /// nor a source law varying with temperature, pressure or composition.
    /// It controls volumetric viscous stress, not reversible pressure work or
    /// compressible thermodynamics. Higher moments use the authored numerical
    /// rate; their damping is not deposited as heat. The native strict `(0,2)`
    /// window requires positive representable zeta; zero refuses without a
    /// viscosity floor. The 2-D and SGS owners refuse this collision family.
    pub fn with_bulk_viscosity(
        mut self,
        bulk_viscosity_pa_s: f64,
        higher_order_rate: f64,
    ) -> Result<Self, MaterialLbmError> {
        let tau = crate::material::physical_bulk_tau(
            self.material.density_kg_m3(),
            self.velocity_scale,
            self.cell_size_m,
            bulk_viscosity_pa_s,
        ).map_err(MaterialLbmError::Invalid)?;
        let relaxation_rate = 1.0 / tau;
        CollisionModel3::CentralMomentBulk {
            shear_rate: 1.0 / self.flow.tau,
            bulk_rate: relaxation_rate,
            higher_order_rate,
        }.validate().map_err(|_| MaterialLbmError::Invalid(
            "physical bulk collision rates are unrepresentable",
        ))?;
        self.bulk_viscosity = Some(BulkViscosity {
            viscosity_pa_s: bulk_viscosity_pa_s,
            relaxation_rate,
            higher_order_rate,
        });
        Ok(self)
    }

    /// Authored independent bulk viscosity [Pa s]. None retains `zeta=2mu/3`
    /// from the current molecular shear viscosity, including local updates.
    #[must_use]
    pub const fn bulk_viscosity_pa_s(&self) -> Option<f64> {
        match self.bulk_viscosity {
            Some(bulk) => Some(bulk.viscosity_pa_s),
            None => None,
        }
    }

    pub(super) fn apply_collision(&self, lattice: &mut FreeSurface3) -> Result<(), MaterialLbmError> {
        if let Some(bulk) = self.bulk_viscosity {
            lattice.set_central_moment_bulk(bulk.relaxation_rate, bulk.higher_order_rate)
                .map_err(|_| MaterialLbmError::Invalid("physical bulk collision rates are unrepresentable"))?;
        }
        Ok(())
    }
}

impl MaterialFreeSurface3<'_> {
    /// Bind independent physical bulk viscosity without advancing the native
    /// trajectory. Retains shear viscosity, source receipts, populations, mass
    /// and force frame. The new coefficient changes how retained trace stress
    /// is interpreted; use the scaling builder before initialization when an
    /// authored initial velocity gradient must be preserved.
    /// Invalid coefficients leave both physical and native bindings unchanged.
    pub fn set_bulk_viscosity(
        &mut self,
        bulk_viscosity_pa_s: f64,
        higher_order_rate: f64,
    ) -> Result<(), MaterialLbmError> {
        let scaling = self.scaling.with_bulk_viscosity(bulk_viscosity_pa_s, higher_order_rate)?;
        scaling.apply_collision(&mut self.lattice)?;
        self.scaling = scaling;
        Ok(())
    }
}
