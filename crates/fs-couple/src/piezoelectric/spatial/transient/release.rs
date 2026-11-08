//! Zero-time release of ideal voltage sources on the accepted spatial state.
//!
//! Removing a voltage constraint supplies no charge impulse or work. The
//! original electrode becomes charge-controlled at its actual conductor charge,
//! while mechanical motion, material memories and physical fields remain intact.
//! A remaining ideal voltage reference is required on every original component.

use super::{
    ElectrodeCondition, FixedBoundary, Result, TetPiezoelectricTransient, electrode_views, invalid,
    poll, prepare_boundary,
};
use fs_exec::Cx;
use std::sync::Arc;

impl<'model, 'mesh> TetPiezoelectricTransient<'model, 'mesh> {
    /// Disconnect selected ideal voltage sources at the accepted endpoint,
    /// retaining each electrode's actual physical conductor charge.
    ///
    /// This zero-time candidate preserves displacement, velocity, fields,
    /// Maxwell/Debye memories, heat history and all stored energy. It neither
    /// solves a new equilibrium nor resets a constitutive or circuit state.
    /// Subsequent ordinary steps use Charge conditions at the released
    /// electrodes; load-driven steps use Current drives (zero for isolation).
    ///
    /// Indices must be nonempty, unique, in range and currently voltage-driven.
    /// Every connected dielectric/electrode component must retain a voltage
    /// reference. Electrode nodes and mechanical support DOFs remain unchanged,
    /// including the actual displacement and velocity of moving supports.
    ///
    /// A changed reference rebases the retained relative potentials without
    /// subtracting large absolute gauges. A nonfinite shift or loss beyond one
    /// relative rounding unit of an original coordinate refuses. This is a
    /// representability check; subsequent steps retain their independent field
    /// and work audits. Every refusal leaves the accepted body untouched.
    pub fn release_electrodes_candidate(&self, electrodes: &[usize], cx: &Cx<'_>) -> Result<Self> {
        poll(cx)?;
        let count = self.state.electrode_conditions.len();
        if electrodes.is_empty() || electrodes.len() > count {
            return Err(invalid(
                "electrode release requires a nonempty bounded index list",
            ));
        }
        let mut conditions = self.state.electrode_conditions.clone();
        let mut seen = vec![false; count];
        for &id in electrodes {
            poll(cx)?;
            if id >= count || seen[id] {
                return Err(invalid(
                    "released electrode indices must be unique and in range",
                ));
            }
            if !matches!(conditions[id], ElectrodeCondition::Voltage(_)) {
                return Err(invalid(
                    "only an existing prescribed-voltage electrode can be released",
                ));
            }
            let charge = self.state.field.electrodes[id].charge_c;
            if !charge.is_finite() {
                return Err(invalid("released conductor charge is not representable"));
            }
            seen[id] = true;
            conditions[id] = ElectrodeCondition::Charge(charge);
        }

        let n = self.model.nodes_m.len();
        let mut supports = Vec::new();
        for dof in 0..3 * n {
            poll(cx)?;
            if self.boundary.map.index[dof].is_none() {
                supports.push((dof, self.state.field.displacement_m[dof / 3][dof % 3]));
            }
        }
        let views = electrode_views(&self.boundary, &conditions);
        let zero_force = vec![[0.0; 3]; n];
        // The native boundary owner retains complete electrode face traces and
        // their connectivity, and refuses a component with no voltage anchor.
        let map = prepare_boundary(
            self.model.nodes_m,
            self.model.tetrahedra,
            &supports,
            &views,
            &zero_force,
            &self.state.nodal_free_charge_c,
            cx,
        )?;
        let mut relative_voltage_v = Vec::with_capacity(n);
        for (node, &reference_electrode) in map.voltage_reference_electrode.iter().enumerate() {
            poll(cx)?;
            let reference_node = self.boundary.electrode_nodes[reference_electrode][0];
            let reference = self.state.relative_voltage_v[reference_node];
            let original = self.state.relative_voltage_v[node];
            let relative = original - reference;
            let reversed = relative + reference;
            // Ordinary subtraction roundoff is allowed. Erasing a weak retained
            // coordinate under a distant new reference must not silently alter
            // the next step's electric field while old storage stays unchanged.
            if !relative.is_finite()
                || !reversed.is_finite()
                || (reversed - original).abs() > f64::EPSILON * original.abs()
            {
                return Err(invalid(
                    "released electrode reference cannot retain the original relative potentials",
                ));
            }
            relative_voltage_v.push(relative);
        }
        let mut state = self.state.clone();
        state.electrode_conditions = conditions;
        state.relative_voltage_v = relative_voltage_v;
        let candidate = Self {
            model: self.model,
            boundary: Arc::new(FixedBoundary {
                map,
                electrode_nodes: self.boundary.electrode_nodes.clone(),
            }),
            state,
            intrinsic: self.intrinsic.clone(),
        };
        poll(cx)?;
        Ok(candidate)
    }
}
