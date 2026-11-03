//! Instantaneous impulses on the original thermal-plastic state and mass.
use super::*;

/// Complete inertial solid candidate at unchanged deformation and temperature.
#[derive(Debug)]
pub struct FiniteThermoplasticTetImpulse<'a, L: FiniteThermalJ2Material = AffineThermalFiniteJ2> {
    pub solid: FiniteThermoplasticTetSolid<'a, L>,
    /// Actual momentum change and applied/support work; no caloric deposition.
    pub impulse: SolidImpulseReport,
}

impl<'a, L: FiniteThermalJ2Material> FiniteThermoplasticTetSolid<'a, L> {
    /// Propose a zero-duration velocity jump using the already bound mass.
    ///
    /// Impulses are full nodal world components [kg m/s]; repeated entries sum.
    /// Held DOFs keep their accepted velocity, including moving supports. The
    /// shared mass owner accounts for their reaction impulse and signed work,
    /// and measures angular momentum at the current deformed positions.
    ///
    /// This requires `with_inertia` first. The candidate retains the original
    /// reference density, source families and thermal reference, temperature,
    /// deformation, all plastic/creep history, stress and free/internal storage.
    /// No material interval is advanced and no physical or numerical heat is
    /// generated. Subsequent `advance` starts from the retained thermal state
    /// and changed velocity. The caller owns any contact/restitution law and
    /// the opposite impulse on its interacting partner.
    ///
    /// # Errors
    /// Missing reference inertia, invalid input/controls, native solve or
    /// momentum/work refusal, or cancellation. The borrowed solid is unchanged.
    pub fn impulse_candidate(
        &self,
        held_velocity_dofs: &[usize],
        nodal_impulses_kg_m_s: &[(usize, f64)],
        origin_m: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticTetImpulse<'a, L>, FiniteThermoplasticTetError> {
        checkpoint(cx)?;
        let motion = self.motion.as_ref().ok_or_else(|| {
            invalid("finite thermal impulse requires already bound reference inertia")
        })?;
        let state = &self.state;
        let mut positions_m = Vec::with_capacity(state.nodes_m.len());
        for (node, (reference, displacement)) in
            state.nodes_m.iter().zip(&state.displacement_m).enumerate()
        {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            positions_m.push(core::array::from_fn(|axis| {
                reference[axis] + displacement[axis]
            }));
        }
        let (velocity_mps, impulse) = motion
            .mass
            .impulse_candidate(
                &positions_m,
                &motion.velocity_mps,
                held_velocity_dofs,
                nodal_impulses_kg_m_s,
                origin_m,
                controls,
                cx,
            )
            .map_err(inertia_error)?;
        let solid = Self {
            state: FinitePlasticTetState {
                nodes_m: state.nodes_m,
                tetrahedra: state.tetrahedra,
                elements: Arc::clone(&state.elements),
                budget: state.budget,
                displacement_m: state.displacement_m.clone(),
                element_deformation: state.element_deformation.clone(),
                element_states: state.element_states.clone(),
                element_first_piola_pa: state.element_first_piola_pa.clone(),
                stored_energy_j: state.stored_energy_j,
            },
            families: Arc::clone(&self.families),
            temperatures_k: self.temperatures_k.clone(),
            internal_energy_j: self.internal_energy_j,
            motion: Some(ThermalMotion {
                velocity_mps,
                reference_density_kg_m3: Arc::clone(&motion.reference_density_kg_m3),
                mass: Arc::clone(&motion.mass),
            }),
        };
        checkpoint(cx)?;
        Ok(FiniteThermoplasticTetImpulse { solid, impulse })
    }
}
