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
    fn impulse_motion(&self) -> Result<&ThermalMotion, FiniteThermoplasticTetError> {
        self.motion.as_ref().ok_or_else(|| {
            invalid("finite thermal impulse requires already bound reference inertia").into()
        })
    }

    /// Current physical material-point position and work-conjugate velocity.
    /// Requires previously bound inertia; quasistatic motion is not inferred.
    pub fn point_kinematics(
        &self,
        point: TetImpulsePoint,
        cx: &Cx<'_>,
    ) -> Result<SolidPointKinematics, FiniteThermoplasticTetError> {
        checkpoint(cx)?;
        let motion = self.impulse_motion()?;
        crate::inertia::point::kinematics(
            point,
            self.state.tetrahedra,
            self.state.nodes_m,
            &self.state.displacement_m,
            &motion.velocity_mps,
            cx,
        )
        .map_err(|error| FiniteThermoplasticTetError::from(inertia_error(error)))
    }

    /// Directional mobility [1/kg] from the retained supported reference mass.
    /// Temperature, thermal expansion and plastic volume never refresh inertia.
    pub fn point_directional_inverse_mass_per_kg(
        &self,
        held: &[usize],
        point: TetImpulsePoint,
        direction: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<SolidPointMobility, FiniteThermoplasticTetError> {
        checkpoint(cx)?;
        crate::inertia::point::directional_mobility(
            &self.impulse_motion()?.mass,
            self.state.tetrahedra,
            held,
            point,
            direction,
            controls,
            cx,
        )
        .map_err(|error| FiniteThermoplasticTetError::from(inertia_error(error)))
    }

    /// Complete cross-point/direction mobility B M_ff^-1 B^T [1/kg].
    /// Rows retain caller order and the same physical point trace as impulses.
    pub fn point_mobility_matrix_per_kg(
        &self,
        held: &[usize],
        rows: &[(TetImpulsePoint, [f64; 3])],
        maximum_rows: usize,
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<SolidPointMobilityMatrix, FiniteThermoplasticTetError> {
        checkpoint(cx)?;
        crate::inertia::point::mobility_matrix(
            &self.impulse_motion()?.mass,
            self.state.tetrahedra,
            held,
            rows,
            maximum_rows,
            controls,
            cx,
        )
        .map_err(|error| FiniteThermoplasticTetError::from(inertia_error(error)))
    }

    /// Apply one material-point impulse through the exact transpose trace.
    /// Thermal/plastic history and both stored energies remain unchanged.
    pub fn point_impulse_candidate(
        &self,
        held: &[usize],
        point: TetImpulsePoint,
        impulse_world: [f64; 3],
        origin: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticTetImpulse<'a, L>, FiniteThermoplasticTetError> {
        let entries = crate::inertia::point::distribution(
            point,
            self.state.tetrahedra,
            self.state.nodes_m.len(),
            impulse_world,
            cx,
        )
        .map_err(inertia_error)?;
        self.impulse_candidate(held, &entries, origin, controls, cx)
    }

    /// Apply all material-point impulses together on the original mass.
    /// Shared nodes accumulate before the one native velocity solve; this does
    /// not replace a simultaneous contact law with sequential pairwise events.
    pub fn point_impulses_candidate(
        &self,
        held: &[usize],
        impulses: &[(TetImpulsePoint, [f64; 3])],
        origin: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticTetImpulse<'a, L>, FiniteThermoplasticTetError> {
        let entries = crate::inertia::point::distributions(
            self.state.tetrahedra,
            self.state.nodes_m.len(),
            impulses,
            controls,
            cx,
        )
        .map_err(inertia_error)?;
        self.impulse_candidate(held, &entries, origin, controls, cx)
    }

    /// Propose a zero-duration velocity jump using the already bound mass.
    ///
    /// Impulses are full nodal world components [kg m/s]; repeated entries sum.
    /// Held DOFs keep their accepted velocity, including moving supports. The
    /// shared mass owner accounts for their reaction impulse and signed work,
    /// and measures angular momentum at the current deformed positions.
    ///
    /// This requires bound inertia first. The candidate retains the original
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
        let motion = self.impulse_motion()?;
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
                reference_nodal_density_kg_m3: motion
                    .reference_nodal_density_kg_m3
                    .as_ref()
                    .map(Arc::clone),
                mass: Arc::clone(&motion.mass),
            }),
        };
        checkpoint(cx)?;
        Ok(FiniteThermoplasticTetImpulse { solid, impulse })
    }
}
