//! Original regional fluid sources in the native weighted DG owner.
//!
//! A material interface is an exact element face. Pressure and volume flow
//! use the original conservative flux; reflection follows the local physical
//! masses and endpoint impedances, without a separate interface integrator.

use super::*;
use std::sync::Arc;

/// One constant original material evaluation over a contiguous axial region.
/// The first region starts at zero; each later one starts at the previous end.
#[derive(Clone, Debug, PartialEq)]
pub struct DgFluidRegion<M> {
    /// Exclusive downstream coordinate [m], except at the final duct endpoint.
    pub end_m: f64,
    /// Complete original source, including its phase/query policy and receipts.
    pub material: M,
}

/// Immutable source field and the responses actually used by a regional duct.
///
/// This is deliberately a distinct medium from [`DgFluidMaterial`]: uniform
/// live-rebind methods cannot silently rescale a nonuniform operator. Regional
/// material changes currently require explicit reassembly and initialization;
/// no spatial history transfer or thermodynamic background evolution is implied.
#[derive(Clone, Debug)]
pub struct DgRegionalFluid<M> {
    regions: Arc<[DgFluidRegion<M>]>,
    acoustic: Arc<[FluidAcousticState]>,
    thermoviscous: Option<Arc<[ThermoviscousFluidState]>>,
    endpoint: usize,
}

impl<M> DgRegionalFluid<M> {
    /// Every original source and its exact spatial support.
    pub fn regions(&self) -> &[DgFluidRegion<M>] {
        &self.regions
    }
    /// Frozen acoustic responses in the same order as the original regions.
    pub fn acoustic_regions(&self) -> &[FluidAcousticState] {
        &self.acoustic
    }
    /// Frozen caloric/transport responses used by the actual diffusion arms.
    pub fn thermoviscous_regions(&self) -> Option<&[ThermoviscousFluidState]> {
        self.thermoviscous.as_deref()
    }
    fn index_at(&self, x_m: f64) -> usize {
        self.regions
            .partition_point(|region| region.end_m <= x_m)
            .min(self.regions.len() - 1)
    }
}

impl<M: Clone> sealed::Medium for DgRegionalFluid<M> {}
impl<M: Clone> DgMedium for DgRegionalFluid<M> {
    fn density_kg_m3(&self) -> f64 {
        self.acoustic[self.endpoint].density_kg_m3()
    }
    fn sound_speed_m_s(&self) -> f64 {
        self.acoustic[self.endpoint].sound_speed_m_s()
    }
    fn admit(&self) -> Result<(), DuctError> {
        Ok(())
    }
    fn has_spatial_acoustics(&self) -> bool {
        self.acoustic.iter().any(|state| *state != self.acoustic[0])
    }
}

impl DgDuct {
    /// Assemble contiguous original fluid regions with the same native DG
    /// masses, interface flux, ports, compliant walls and exact exponential.
    ///
    /// Each source is evaluated once, and its immutable acoustic response is
    /// integrated only over its own support. Interfaces split the original
    /// cylinder, cone or profile geometry at their exact coordinates; derived
    /// material coefficients are not averaged across interfaces. Matched ends
    /// use their own first/last region, while all original sources are retained.
    ///
    /// Regions describe a prescribed frozen background at equilibrium. They
    /// do not imply miscibility, mean flow, heat/species exchange, moving
    /// interfaces or surface tension. The caller must supply a physically
    /// compatible background. No fictitious gas state is constructed.
    ///
    /// # Errors
    /// Refuses empty/unordered/nonfinite regions, incomplete duct coverage,
    /// source/phase failures, insufficient elements for the combined geometry
    /// and material spans, and the original geometry/quadrature/DOF limits.
    pub fn from_material_regions<M: FluidAcoustics + Clone>(
        duct: &Duct,
        regions: &[DgFluidRegion<M>],
        mesh: DgMesh,
        inlet: DgBoundary,
        outlet: DgBoundary,
    ) -> Result<DgDuct<DgRegionalFluid<M>>, FluidDuctError> {
        if regions.is_empty() || regions.len() > MAX_DG_ELEMENTS {
            return Err(DuctError::BadParameter {
                what: "DG fluid regions must have between one and MAX_DG_ELEMENTS entries",
            }
            .into());
        }
        let length: f64 = duct.segments.iter().map(Segment::axial_length).sum();
        let mut previous = 0.0;
        for region in regions {
            if !region.end_m.is_finite() || region.end_m <= previous {
                return Err(DuctError::BadParameter {
                    what: "DG fluid region ends must be finite and strictly increasing from zero",
                }
                .into());
            }
            previous = region.end_m;
        }
        if !length.is_finite() || length <= 0.0 || previous != length {
            return Err(DuctError::BadParameter {
                what: "DG fluid regions must cover the exact complete duct length",
            }
            .into());
        }
        let acoustic: Vec<_> = regions
            .iter()
            .map(|region| region.material.acoustic_state())
            .collect::<Result<_, _>>()?;
        let medium = DgRegionalFluid {
            regions: regions.to_vec().into(),
            acoustic: acoustic.into(),
            thermoviscous: None,
            endpoint: 0,
        };
        let mut outlet_medium = medium.clone();
        outlet_medium.endpoint = regions.len() - 1;
        let stations: Vec<_> = regions[..regions.len() - 1]
            .iter()
            .map(|region| region.end_m)
            .collect();
        Ok(DgDuct::assemble_at_stations(
            duct,
            &medium,
            &outlet_medium,
            mesh,
            inlet,
            outlet,
            &stations,
            &mut |x| Ok::<_, DuctError>(medium.acoustic[medium.index_at(x)]),
        )?)
    }
}

impl<M: Clone> DgDuct<DgRegionalFluid<M>> {
    /// Complete original sources and their exact axial supports.
    pub fn material_regions(&self) -> &[DgFluidRegion<M>] {
        self.gas.regions()
    }
    /// Actual acoustic coefficients used in each region's physical mass.
    pub fn acoustic_regions(&self) -> &[FluidAcousticState] {
        self.gas.acoustic_regions()
    }
    /// Actual regional caloric/transport tuples used by attached diffusion.
    pub fn thermoviscous_regions(&self) -> Option<&[ThermoviscousFluidState]> {
        self.gas.thermoviscous_regions()
    }
}

impl<M: ThermoviscousFluid + Clone> DgDuct<DgRegionalFluid<M>> {
    /// Attach the original spatial Bessel kernel using each region's SAME
    /// acoustic, viscosity, conductivity and Cp/Cv evaluation. Every memory
    /// remains in the native pressure/flow state and exponential. The original
    /// local-cylinder, rigid isothermal wall approximation and band admission
    /// apply at each Gauss point. Cp == Cv omits only that region's thermal arms.
    pub fn with_spatial_viscothermal(
        mut self,
        resolution: DgWallResolution,
    ) -> Result<Self, FluidDuctError> {
        let states: Vec<_> = self
            .gas
            .regions
            .iter()
            .map(|region| region.material.thermoviscous_state())
            .collect::<Result<_, _>>()?;
        if states
            .iter()
            .zip(self.gas.acoustic.iter())
            .any(|(state, acoustic)| state.acoustic() != *acoustic)
        {
            return Err(DuctError::BadParameter {
                what: "DG regional thermoviscous and acoustic views must be the same original states",
            }.into());
        }
        let geometry = self.geometry.clone();
        let source = self.gas.clone();
        self.attach_spatial_viscothermal::<FluidDuctError, _>(&geometry, resolution, &mut |x| {
            Ok(states[source.index_at(x)])
        })?;
        let states: Arc<[ThermoviscousFluidState]> = states.into();
        self.gas.thermoviscous = Some(states.clone());
        self.outlet_gas.thermoviscous = Some(states);
        Ok(self)
    }
}
