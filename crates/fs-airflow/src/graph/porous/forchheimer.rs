//! One retained porous calibration supplies both coefficients to the native
//! frozen-density pressure graph. The graph's original nonlinear inversion,
//! conservation, fan intersection and tangent owners remain unchanged.

use super::*;

impl PassiveBranch {
    /// Bind the original Darcy–Forchheimer material to one passive graph edge.
    ///
    /// The same original skeleton, superficial area and coherent gas state give
    /// `A = mu sum(L/k)/area` and `B = rho sum(beta L)/area²` in
    /// `dp = A Q + B Q |Q|`. Neither coefficient is inferred from the other.
    /// The retained binding admits the actual signed flow and both face pressures
    /// under the declared frozen-density allowance. This does not activate a
    /// compressible graph or a spatial nonlinear Brinkman equation.
    ///
    /// Graph admission checks both retained coefficients together: extracting
    /// the linear component cannot turn a calibrated inertial filter into a
    /// Darcy-only edge, and altering the quadratic coefficient cannot silently
    /// detach it from its material. Independently authored mixed losses retain
    /// their existing API.
    ///
    /// # Errors
    /// Refuses invalid coefficient uncertainty/pressure allowance, unsupported
    /// material state, slip, unrepresentable coefficients, or a sealed path.
    #[allow(clippy::too_many_arguments)]
    pub fn from_darcy_forchheimer(
        name: impl Into<String>,
        from: usize,
        to: usize,
        barrier: DarcyGasBarrier,
        gas: GasState,
        area_m2: f64,
        maximum_relative_pressure_drop: f64,
        uncertainty_rel: f64,
        source: SourceProvenance,
        tolerance_basis: ToleranceBasis,
    ) -> Result<Self, GraphError> {
        if !(maximum_relative_pressure_drop.is_finite()
            && maximum_relative_pressure_drop > 0.0
            && maximum_relative_pressure_drop < 1.0)
        {
            return Err(GraphError::InvalidInput(
                "porous relative pressure allowance must be in (0, 1)",
            ));
        }
        let (viscous, inertial) = barrier
            .volume_resistance_coefficients(&gas, area_m2)
            .map_err(GraphError::Material)?
            .ok_or(GraphError::InvalidInput(
                "a sealed porous path must be omitted from the graph",
            ))?;
        let mut linear_loss = LinearLossElement::new(
            name,
            LinearLossResistance::new(viscous),
            uncertainty_rel,
            source,
            tolerance_basis,
        )?;
        linear_loss.material = Some(DarcyLossBinding {
            barrier,
            gas,
            area_m2,
            maximum_relative_pressure_drop,
        });
        linear_loss.quadratic_binding = Some(LossResistance::new(inertial));
        let mut loss = linear_loss.zero_quadratic_metadata();
        loss.resistance = LossResistance::new(inertial);
        if inertial == 0.0 {
            Ok(Self::Linear {
                from,
                to,
                loss: linear_loss,
            })
        } else {
            Ok(Self::LinearQuadratic {
                branch: GraphBranch { from, to, loss },
                linear_loss,
            })
        }
    }
}
