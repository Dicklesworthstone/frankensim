//! Source-resolved molecular transport and local turbulent heat/stress share
//! the original native buoyant flow, thermal endpoint and publication.

use super::*;
use crate::lattice3::open_boundaries::sgs::{
    OpenThermalSgsConfig3, OpenThermalSgsError3, OpenThermalSgsWorkspace3,
};
use crate::sgs3::{SgsLattice3, SgsReport};

/// Original material/transport refusal or shared local SGS coefficient refusal.
#[derive(Clone, Debug)]
pub enum BoussinesqOpenMaterialSgsError3<E> {
    Material(BoussinesqOpenMaterialError3<E>),
    Sgs(OpenThermalSgsError3<E>),
}
impl<E> From<BoussinesqOpenMaterialError3<E>> for BoussinesqOpenMaterialSgsError3<E> {
    fn from(error: BoussinesqOpenMaterialError3<E>) -> Self {
        Self::Material(error)
    }
}
impl<E> From<BoussinesqOpenError3<E>> for BoussinesqOpenMaterialSgsError3<E> {
    fn from(error: BoussinesqOpenError3<E>) -> Self {
        Self::Material(error.into())
    }
}
impl<E> From<OpenThermalSgsError3<E>> for BoussinesqOpenMaterialSgsError3<E> {
    fn from(error: OpenThermalSgsError3<E>) -> Self {
        Self::Sgs(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for BoussinesqOpenMaterialSgsError3<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native buoyant material SGS: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for BoussinesqOpenMaterialSgsError3<E> {}

/// One accepted native endpoint with original material sources and SGS evidence.
#[derive(Debug)]
pub struct BoussinesqOpenMaterialSgsStep3<F> {
    /// Original heat/force/flow receipt and opening/actual-endpoint sources.
    /// Reported material relaxation values remain molecular, before SGS.
    pub material: BoussinesqOpenMaterialStep3<F>,
    /// Incoming resolved/modeled shear transfer; never an added heat source.
    pub sgs: SgsReport,
    pub maximum_eddy_diffusivity_lu: f64,
    pub maximum_diffusivity_to_molecular: f64,
}

impl<L: EquilibriumLattice3<Q> + SgsLattice3<Q>, const Q: usize> BoussinesqOpenGrid3<L, Q> {
    /// Advance with source-dependent molecular mu/k and explicit local SGS.
    ///
    /// Query each original material once at its accepted opening temperature;
    /// bounded reference rho/Cp and the original SI planner supply fresh native
    /// molecular rates. The existing Smagorinsky/WALE/Vreman preparation samples
    /// accepted physical velocities with their retained spatial half-force, then
    /// adds nu_t and alpha_t=nu_t/Pr_t to those fresh rates. Collision families,
    /// bulk/higher rates and TRT wall placement remain with their original owners.
    ///
    /// The original buoyancy stage transports heat once and advances momentum
    /// once. Actual final-temperature material queries must pass before its
    /// shared final cancellation gate and joint publication, including finite
    /// wall history, source/heat receipts and accepted force epochs.
    ///
    /// Reuse both existing workspaces. The material workspace always starts
    /// from freshly resolved molecular transport; an earlier augmented SGS
    /// field is never fed back into it. Both workspaces are scratch and may be
    /// overwritten on refusal. Only successful reports/arrays may be consumed.
    /// Source callbacks must be repeatable; their external effects are not
    /// rolled back. Only the two original material-report vectors allocate.
    ///
    /// Molecular and modeled rates are held over this interval. The existing
    /// affine buoyancy law remains authoritative. This adds no source-density
    /// force law, mechanical heating or higher-order coupled integration.
    #[allow(clippy::too_many_arguments)]
    pub fn step_serial_with_material_transport_and_sgs<F: FluidTransport>(
        &mut self,
        cx: &Cx<'_>,
        config: OpenThermalSgsConfig3,
        collision_templates: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        material_workspace: &mut BoussinesqOpenMaterialWorkspace3<L::Collision>,
        sgs_workspace: &mut OpenThermalSgsWorkspace3<L::Collision>,
        origin: [f64; 3],
        mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<BoussinesqOpenMaterialSgsStep3<F>, BoussinesqOpenMaterialSgsError3<L::Error>> {
        self.advance_material_sgs_with_poll(
            None,
            config,
            collision_templates,
            heat_per_step,
            material_workspace,
            sgs_workspace,
            origin,
            &mut resolve,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// The same source/SGS transaction with a prescribed next mechanical load.
    /// SGS uses the opening total force, including opening buoyancy. Successful
    /// publication alone swaps the caller's endpoint background with the old one;
    /// source, closure, native or late cancellation refusal preserves that input.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_serial_with_material_transport_and_sgs<F: FluidTransport>(
        &mut self,
        cx: &Cx<'_>,
        config: OpenThermalSgsConfig3,
        collision_templates: Option<&[L::Collision]>,
        background_at_end: &mut Vec<[f64; 3]>,
        heat_per_step: &[f64],
        material_workspace: &mut BoussinesqOpenMaterialWorkspace3<L::Collision>,
        sgs_workspace: &mut OpenThermalSgsWorkspace3<L::Collision>,
        origin: [f64; 3],
        mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<BoussinesqOpenMaterialSgsStep3<F>, BoussinesqOpenMaterialSgsError3<L::Error>> {
        self.advance_material_sgs_with_poll(
            Some(background_at_end),
            config,
            collision_templates,
            heat_per_step,
            material_workspace,
            sgs_workspace,
            origin,
            &mut resolve,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::lattice3::open_boundaries::thermal::buoyancy) fn advance_material_sgs_with_poll<
        F: FluidTransport,
    >(
        &mut self,
        background_at_end: Option<&mut Vec<[f64; 3]>>,
        config: OpenThermalSgsConfig3,
        collision_templates: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        material_workspace: &mut BoussinesqOpenMaterialWorkspace3<L::Collision>,
        sgs_workspace: &mut OpenThermalSgsWorkspace3<L::Collision>,
        origin: [f64; 3],
        resolve: &mut impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<BoussinesqOpenMaterialSgsStep3<F>, BoussinesqOpenMaterialSgsError3<L::Error>> {
        let applied = self.prepare_material(
            background_at_end.as_deref().map(Vec::as_slice),
            collision_templates,
            heat_per_step,
            material_workspace,
            resolve,
            poll,
        )?;
        let (sgs, maximum_eddy_diffusivity_lu, maximum_diffusivity_to_molecular) =
            self.flow.inner.prepare_open_thermal_sgs(
                config,
                Some(material_workspace.momentum_collisions()),
                material_workspace.thermal_relaxation(),
                sgs_workspace,
                Some(&self.flow.forces),
                &mut *poll,
            )?;
        let (transport, step) = self.stage_with_diffusion_poll(
            background_at_end.as_deref().map(Vec::as_slice),
            Some(sgs_workspace.momentum_collisions()),
            heat_per_step,
            Some(sgs_workspace.thermal_relaxation()),
            None,
            origin,
            poll,
        )?;
        let material = self.finish_material_with_poll(
            background_at_end,
            collision_templates,
            material_workspace,
            applied,
            transport,
            step,
            resolve,
            poll,
        )?;
        Ok(BoussinesqOpenMaterialSgsStep3 {
            material,
            sgs,
            maximum_eddy_diffusivity_lu,
            maximum_diffusivity_to_molecular,
        })
    }
}
