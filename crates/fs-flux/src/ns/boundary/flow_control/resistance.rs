//! Integral hydraulic resistance coupled to the native resolved port responses.

use super::{
    PreparedStokesFlowControl, StokesAssembly, StokesBoundaryError, StokesSolution, invalid,
    port_flow,
};
use fs_la::factor::lu;
use fs_math::dd::Dd;

/// One external reservoir and linear resistance attached to a prepared port.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StokesResistancePort {
    /// Reservoir pressure increment relative to the original port traction,
    /// divided by the SAME reference density [m²/s²]. Signed and finite.
    pub reservoir_pressure_m2_s2: f64,
    /// Nonnegative external resistance divided by reference density [1/s].
    /// The port law is p = p_reservoir + resistance * Q_outward, where Q is
    /// volume flow per unit depth [m²/s]. Zero prescribes reservoir pressure.
    pub resistance_s_inv: f64,
}

/// A resolved flow satisfying every external reservoir/resistance law.
pub struct StokesResistanceSolution {
    /// Accepted continuum field, pressure datum and discrete boundary reactions.
    pub flow: StokesSolution,
    /// Uniform pressure increments added to original port tractions [m²/s²].
    pub pressure_offsets_m2_s2: Vec<f64>,
    /// Signed outward flow per unit depth [m²/s], in prepared port order.
    pub outward_volume_rates_m2_s: Vec<f64>,
    /// Each external resistance's nonnegative R*Q² [m⁴/s³]. Multiply by the
    /// reference density for W/m. This is outside the resolved fluid domain;
    /// it is NOT added to its bulk viscous heat or wall-slip dissipation.
    pub resistance_dissipation_m4_s3: Vec<f64>,
    /// Largest |p - p_reservoir - R*Q| at the accepted field [m²/s²].
    pub maximum_pressure_residual_m2_s2: f64,
}

impl PreparedStokesFlowControl<'_, '_> {
    /// Couple resolved laminar flow to finite external hydraulic resistances.
    ///
    /// With signed outward Q, each port adds total traction -p*n with
    /// p = p_reservoir + R*Q. Inflow therefore drops pressure from the external
    /// reservoir towards the fluid; reversing flow reverses that drop while
    /// R*Q² remains nonnegative. The original tangential traction is retained.
    /// Reservoir pressures are INCREMENTS relative to the original traction.
    ///
    /// All ports are solved together: Q = Q0 + M*p uses the entire cached
    /// native response matrix, including off-diagonal interactions. Only the
    /// small system (I - diag(R)*M)*p = p_reservoir + diag(R)*Q0 is factored.
    /// No new full-flow solve, boundary/viscosity query or imposed inlet
    /// profile is needed. The cache's mesh, coefficients, body force, wall
    /// motion and any transient history remain fixed. Reprepare to change them.
    /// Each connected controlled component still needs the uncontrolled
    /// traction datum required by `prepare_flow_control`.
    ///
    /// This is a linear external pressure/flow law, not a fitted turbulent
    /// loss, compliance, fluid inertia or automatic heat transfer model.
    /// A frozen Oseen or rheological response remains frozen.
    ///
    /// # Errors
    /// Refuses incomplete/nonfinite inputs, negative resistance, singular
    /// reduced systems, unrepresentable outputs, failed native equations or
    /// any port-law defect above absolute_tolerance + native_tol*scale.
    /// The cache is immutable and no partial result is published.
    pub fn solve_resistance_ports(
        &self,
        loads: &[StokesResistancePort],
        absolute_pressure_tolerance_m2_s2: f64,
    ) -> Result<StokesResistanceSolution, StokesBoundaryError> {
        let count = self.ports.len();
        admit_resistance_loads(count, loads, absolute_pressure_tolerance_m2_s2)?;
        let mut matrix = vec![0.0; self.mobility.len()];
        let mut rhs = vec![0.0; count];
        let mut row_scales = vec![0.0_f64; count];
        for (i, load) in loads.iter().enumerate() {
            let resistance = Dd::from_f64(load.resistance_s_inv);
            rhs[i] = (Dd::from_f64(load.reservoir_pressure_m2_s2)
                + resistance * Dd::from_f64(self.baseline_flows[i]))
            .to_f64();
            for j in 0..count {
                let entry = (Dd::from_f64(if i == j { 1.0 } else { 0.0 })
                    - resistance * Dd::from_f64(self.mobility[i * count + j]))
                .to_f64();
                matrix[i * count + j] = entry;
                row_scales[i] = row_scales[i].max(entry.abs());
            }
        }
        if matrix.iter().chain(&rhs).any(|x| !x.is_finite())
            || row_scales.iter().any(|&scale| scale == 0.0)
        {
            return Err(invalid("unrepresentable external hydraulic response"));
        }
        let mut scaled = matrix.clone();
        let mut pressures = rhs.clone();
        for i in 0..count {
            pressures[i] /= row_scales[i];
            for j in 0..count {
                scaled[i * count + j] /= row_scales[i];
            }
        }
        let factor = lu(&scaled, count).map_err(|_| StokesBoundaryError::SingularSystem)?;
        factor.solve(&mut pressures);
        for _ in 0..2 {
            let mut correction = Vec::with_capacity(count);
            for i in 0..count {
                let mut residual = Dd::from_f64(rhs[i]);
                for (j, &pressure) in pressures.iter().enumerate() {
                    residual =
                        residual - Dd::from_f64(matrix[i * count + j]) * Dd::from_f64(pressure);
                }
                correction.push(residual.to_f64() / row_scales[i]);
            }
            factor.solve(&mut correction);
            for (pressure, delta) in pressures.iter_mut().zip(correction) {
                *pressure += delta;
            }
        }
        let flow = self.report_pressures(&pressures)?;
        let (flows, dissipation, maximum_pressure_residual_m2_s2) = checked_resistance_ports(
            &self.assembly,
            &self.ports,
            loads,
            absolute_pressure_tolerance_m2_s2,
            &pressures,
            &flow.field.x,
        )?;
        Ok(StokesResistanceSolution {
            flow,
            pressure_offsets_m2_s2: pressures,
            outward_volume_rates_m2_s: flows,
            resistance_dissipation_m4_s3: dissipation,
            maximum_pressure_residual_m2_s2,
        })
    }
}

pub(super) fn admit_resistance_loads(
    count: usize,
    loads: &[StokesResistancePort],
    absolute_pressure_tolerance_m2_s2: f64,
) -> Result<(), StokesBoundaryError> {
    if loads.len() != count
        || !absolute_pressure_tolerance_m2_s2.is_finite()
        || absolute_pressure_tolerance_m2_s2 < 0.0
        || loads.iter().any(|load| {
            !load.reservoir_pressure_m2_s2.is_finite()
                || !load.resistance_s_inv.is_finite()
                || load.resistance_s_inv < 0.0
        })
    {
        return Err(invalid(
            "finite complete reservoir pressures, nonnegative resistances and tolerance required",
        ));
    }
    Ok(())
}

// The frozen response and the reassembled nonlinear endpoint must use this
// SAME integrated-flow law check, including relaxed pressure unknowns.
pub(super) fn checked_resistance_ports(
    assembly: &StokesAssembly<'_, '_>,
    ports: &[impl AsRef<[usize]>],
    loads: &[StokesResistancePort],
    absolute_pressure_tolerance_m2_s2: f64,
    pressures: &[f64],
    velocity: &[f64],
) -> Result<(Vec<f64>, Vec<f64>, f64), StokesBoundaryError> {
    admit_resistance_loads(ports.len(), loads, absolute_pressure_tolerance_m2_s2)?;
    if pressures.len() != ports.len() || pressures.iter().any(|p| !p.is_finite()) {
        return Err(invalid("finite complete external port pressures required"));
    }
    let mut flows = Vec::with_capacity(ports.len());
    let mut dissipation = Vec::with_capacity(ports.len());
    let mut maximum_pressure_residual_m2_s2 = 0.0_f64;
    for (i, load) in loads.iter().enumerate() {
        let q = port_flow(assembly, ports[i].as_ref(), velocity);
        let drop = load.resistance_s_inv * q;
        let loss = drop * q;
        let residual = (Dd::from_f64(pressures[i])
            - Dd::from_f64(load.reservoir_pressure_m2_s2)
            - Dd::from_f64(load.resistance_s_inv) * Dd::from_f64(q))
        .to_f64()
        .abs();
        let scale = pressures[i]
            .abs()
            .max(load.reservoir_pressure_m2_s2.abs())
            .max(drop.abs());
        let tolerance = absolute_pressure_tolerance_m2_s2 + assembly.params.tol * scale;
        if ![q, drop, loss, residual, tolerance]
            .iter()
            .all(|x| x.is_finite())
        {
            return Err(invalid(
                "unrepresentable hydraulic port pressure, flow or power",
            ));
        }
        if residual > tolerance {
            return Err(StokesBoundaryError::ResistanceLawNotConverged {
                port: i,
                residual_m2_s2: residual,
                tolerance_m2_s2: tolerance,
            });
        }
        maximum_pressure_residual_m2_s2 = maximum_pressure_residual_m2_s2.max(residual);
        flows.push(q);
        dissipation.push(loss);
    }
    Ok((flows, dissipation, maximum_pressure_residual_m2_s2))
}
