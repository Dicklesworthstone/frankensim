//! Instantaneous impulses retaining every original Hex8 material point.
use super::*;

/// Complete motion candidate at unchanged thermal and quadrature states.
#[derive(Debug)]
pub struct FiniteThermoplasticHexMotionImpulse<
    'a,
    L: FiniteThermalJ2Material = AffineThermalFiniteJ2,
> {
    pub motion: FiniteThermoplasticHexMotion<'a, L>,
    /// Actual momentum change and applied/support work; no caloric deposition.
    pub impulse: SolidImpulseReport,
}

impl<'a, L: FiniteThermalJ2Material> FiniteThermoplasticHexMotion<'a, L> {
    /// Physical trilinear position/velocity at an original material point.
    /// Mean-dilatation constitutive coordinates never replace this geometry.
    pub fn point_kinematics(
        &self,
        point: HexImpulsePoint,
        cx: &Cx<'_>,
    ) -> Result<SolidPointKinematics, FiniteThermoplasticHexError> {
        crate::inertia::point::hex_kinematics(
            point,
            self.solid.hexahedra,
            self.solid.nodes_m,
            &self.solid.displacement_m,
            &self.velocity_mps,
            cx,
        )
        .map_err(inertia_error)
    }

    /// Directional mobility [1/kg] on the original supported Hex8 mass.
    pub fn point_directional_inverse_mass_per_kg(
        &self,
        held: &[usize],
        point: HexImpulsePoint,
        direction: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<SolidPointMobility, FiniteThermoplasticHexError> {
        crate::inertia::point::hex_directional_mobility(
            &self.mass,
            self.solid.hexahedra,
            held,
            point,
            direction,
            controls,
            cx,
        )
        .map_err(inertia_error)
    }

    /// Complete cross-point/direction mobility B M_ff^-1 B^T [1/kg].
    /// Original trilinear point rows retain their full shared-body coupling.
    pub fn point_mobility_matrix_per_kg(
        &self,
        held: &[usize],
        rows: &[(HexImpulsePoint, [f64; 3])],
        maximum_rows: usize,
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<SolidPointMobilityMatrix, FiniteThermoplasticHexError> {
        crate::inertia::point::hex_mobility_matrix(
            &self.mass,
            self.solid.hexahedra,
            held,
            rows,
            maximum_rows,
            controls,
            cx,
        )
        .map_err(inertia_error)
    }

    /// Apply a world impulse through the transpose of the physical point trace.
    /// No temperature, constitutive history, storage or material time changes.
    pub fn point_impulse_candidate(
        &self,
        held: &[usize],
        point: HexImpulsePoint,
        impulse_world: [f64; 3],
        origin: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticHexMotionImpulse<'a, L>, FiniteThermoplasticHexError> {
        let entries = crate::inertia::point::hex_distribution(
            point,
            self.solid.hexahedra,
            self.velocity_mps.len(),
            impulse_world,
            cx,
        )
        .map_err(inertia_error)?;
        self.impulse_candidate(held, &entries, origin, controls, cx)
    }

    /// Accumulate every material-point impulse before one native mass solve.
    /// Thermal, physical and assumed quadrature histories are all retained.
    pub fn point_impulses_candidate(
        &self,
        held: &[usize],
        impulses: &[(HexImpulsePoint, [f64; 3])],
        origin: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticHexMotionImpulse<'a, L>, FiniteThermoplasticHexError> {
        let entries = crate::inertia::point::hex_distributions(
            self.solid.hexahedra,
            self.velocity_mps.len(),
            impulses,
            controls,
            cx,
        )
        .map_err(inertia_error)?;
        self.impulse_candidate(held, &entries, origin, controls, cx)
    }

    /// Propose a zero-duration velocity jump through the original Hex8 mass.
    ///
    /// Impulses are full nodal world components [kg m/s]; repeated entries sum.
    /// Held DOFs retain their current velocity and their reaction impulse/work.
    /// Current physical positions define angular momentum about `origin_m`;
    /// the assumed mean-dilatation field never replaces that physical geometry.
    ///
    /// All original reference quadrature, density, formulation, source families
    /// and thermal reference survive. Temperature, physical/assumed deformation,
    /// each quadrature point's plastic/creep history, stress and free/internal
    /// energy remain exactly unchanged. No thermal/material interval is advanced
    /// and no heat is generated. The next native `advance` uses the new velocity
    /// with that complete accepted state. The caller owns the contact law and
    /// equal-and-opposite partner impulse.
    ///
    /// # Errors
    /// Invalid input/controls, native mass solve or momentum/work refusal, or
    /// cancellation. No part of the borrowed accepted motion changes.
    pub fn impulse_candidate(
        &self,
        held_velocity_dofs: &[usize],
        nodal_impulses_kg_m_s: &[(usize, f64)],
        origin_m: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticHexMotionImpulse<'a, L>, FiniteThermoplasticHexError> {
        checkpoint(cx)?;
        let old = &self.solid;
        let mut positions_m = Vec::with_capacity(old.nodes_m.len());
        for (node, (reference, displacement)) in
            old.nodes_m.iter().zip(&old.displacement_m).enumerate()
        {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            positions_m.push(core::array::from_fn(|axis| {
                reference[axis] + displacement[axis]
            }));
        }
        let (velocity_mps, impulse) = self
            .mass
            .impulse_candidate(
                &positions_m,
                &self.velocity_mps,
                held_velocity_dofs,
                nodal_impulses_kg_m_s,
                origin_m,
                controls,
                cx,
            )
            .map_err(inertia_error)?;
        let solid = FiniteThermoplasticHexSolid {
            nodes_m: old.nodes_m,
            hexahedra: old.hexahedra,
            elements: Arc::clone(&old.elements),
            families: Arc::clone(&old.families),
            formulation: old.formulation,
            budget: old.budget,
            temperatures_k: old.temperatures_k.clone(),
            displacement_m: old.displacement_m.clone(),
            physical_deformation: old.physical_deformation.clone(),
            assumed_deformation: old.assumed_deformation.clone(),
            states: old.states.clone(),
            first_piola_pa: old.first_piola_pa.clone(),
            free_energy_j: old.free_energy_j,
            internal_energy_j: old.internal_energy_j,
        };
        checkpoint(cx)?;
        Ok(FiniteThermoplasticHexMotionImpulse {
            motion: Self {
                solid,
                velocity_mps,
                reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                mass: Arc::clone(&self.mass),
            },
            impulse,
        })
    }
}
