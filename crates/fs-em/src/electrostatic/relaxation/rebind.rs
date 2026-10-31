//! Externally supplied material replacement at retained physical free charge.
//!
//! Held electrodes may exchange charge; their work is separate from the signed
//! change of the dielectric Hamiltonian. This operation advances neither time
//! nor Ohmic heat, and does not model material formation or an insertion path.

use super::*;

#[derive(Clone, Copy, Debug)]
pub struct DielectricMaterialRebindControls {
    pub linear: CheckedCgConfig,
    pub charge_tolerance_c: f64,
    pub energy_tolerance_j: f64,
    /// Bounds the absolute material work and electrode work separately [J].
    /// Zero admits a conductivity-only replacement with exactly zero work.
    pub maximum_absolute_work_j: f64,
}

/// A zero-duration replacement. Two-entry arrays are ordered [before, after].
#[derive(Clone, Debug)]
pub struct DielectricMaterialRebind {
    pub time_s: f64,
    pub stored_energy_j: [f64; 2],
    /// Signed work delivered by the held-voltage electrodes: sum(V * delta_q).
    pub supplied_work_j: f64,
    /// Signed work delivered by the external material-replacement owner.
    /// Negative work must be received by that owner, never discarded as heat.
    pub parameter_work_j: f64,
    /// Actual electrode charge transfer [C], not current over an invented dt.
    pub nodal_electrode_charge_c: Vec<Option<f64>>,
    /// Native recovered charge versus the conserved free-node charge [C].
    pub maximum_free_charge_residual_c: f64,
    pub electrode_charge_residual_c: f64,
    pub energy_residual_j: f64,
    pub electric_cells: Vec<DielectricEnergyCell>,
    /// Complete material-state identities in original triangle order.
    pub material_state_identities: Vec<[[u8; 32]; 2]>,
    pub linear: Option<SolveReport>,
}

impl<C: PlanarConductorMaterial> PlanarLeakyDielectric<C> {
    /// Replace dielectric/conductor data on the SAME mesh and accepted clock.
    ///
    /// Every free-node charge is retained bit-for-bit. The original scalar
    /// potential solver finds the new voltage from `K_epsilon_new V_new=q_old`
    /// at free nodes; original prescribed electrode potentials stay fixed.
    /// Electrode charges may change and their explicit work is admitted beside
    /// the independently calculated material work
    /// `P=-sum(volume * E_old.(epsilon_new-epsilon_old).E_new/2)`.
    /// Native cell storage must close `delta_U=W_electrodes+P` before publication.
    ///
    /// The new source cards, conductivity basis, principal coefficients and
    /// material orientations are data. Declared mesh-coordinate frames must
    /// agree. New dielectric sources must be validity-wide constant; normal
    /// later temperature/phase refresh remains pinned to these new sources.
    /// Conductivity-only substitutions preserve potentials/charges exactly and
    /// affect subsequent relaxation and local Joule heat through the native
    /// conductor owner. No state of `self` is modified, including on refusal.
    ///
    /// This is an externally controlled instantaneous constitutive replacement,
    /// not continuous thermal evolution, mass transport, geometry motion or a
    /// prediction of material formation energy. It models no intermediate
    /// insertion path, dielectric dispersion or transient radiation.
    pub fn rebind_material_candidate<D: PlanarConductorMaterial>(
        &self,
        conductors: Vec<D>,
        laws: Vec<PlanarDielectricLaw>,
        limits: WeightedAssemblyLimits,
        controls: DielectricMaterialRebindControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(PlanarLeakyDielectric<D>, DielectricMaterialRebind), Error> {
        poll(cancelled)?;
        if !controls.linear.tolerance.is_finite()
            || controls.linear.tolerance <= 0.0
            || controls.linear.tolerance >= 1.0
            || !controls.charge_tolerance_c.is_finite()
            || controls.charge_tolerance_c <= 0.0
            || !controls.energy_tolerance_j.is_finite()
            || controls.energy_tolerance_j <= 0.0
            || !controls.maximum_absolute_work_j.is_finite()
            || controls.maximum_absolute_work_j < 0.0
        {
            return Err(Error::Invalid(
                "positive finite replacement tolerances and nonnegative work allowance required",
            ));
        }
        let conduction = PlanarConduction::new(
            self.conduction.mesh().clone(),
            conductors,
            limits,
            cancelled,
        )?;
        let dielectric = dielectric(&conduction, laws, None, limits, cancelled)?;
        let frame0 = self
            .conduction
            .conductors()
            .iter()
            .find_map(|c| c.frame().map(|f| f.geometry_frame()))
            .or_else(|| {
                self.dielectric
                    .permittivity_laws()
                    .iter()
                    .find_map(|law| law.frame().map(|f| f.geometry_frame()))
            });
        let frame1 = conduction
            .conductors()
            .iter()
            .find_map(|c| c.frame().map(|f| f.geometry_frame()))
            .or_else(|| {
                dielectric
                    .permittivity_laws()
                    .iter()
                    .find_map(|law| law.frame().map(|f| f.geometry_frame()))
            });
        if matches!((frame0, frame1), (Some(a), Some(b)) if a != b) {
            return Err(Error::Invalid(
                "material replacement changed the mesh-coordinate frame",
            ));
        }
        let mut prescribed = vec![None; self.voltage_v.len()];
        for &v in &self.voltage_vertices {
            poll(cancelled)?;
            prescribed[v] = Some(self.voltage_v[v]);
        }
        let unchanged = self
            .dielectric
            .permittivity_coefficients()
            .iter()
            .zip(dielectric.permittivity_coefficients())
            .all(|(a, b)| {
                a.matrix_f_m()
                    .into_iter()
                    .flatten()
                    .zip(b.matrix_f_m().into_iter().flatten())
                    .all(|(x, y)| x.to_bits() == y.to_bits())
            });
        let (voltage_v, linear) = if unchanged {
            (self.voltage_v.clone(), None)
        } else {
            scalar_potential::solve(
                &dielectric.stiffness,
                &self.nodal_charge_c,
                &prescribed,
                controls.linear,
                cancelled,
            )?
        };
        let previous = self.dielectric.cell_response(&self.voltage_v, cancelled)?;
        let endpoint = dielectric.cell_response(&voltage_v, cancelled)?;
        let (electric_cells, parameter_work_j) = parameter_work_cells(
            &self.dielectric,
            &dielectric,
            &previous,
            &endpoint,
            cancelled,
        )?;
        let recovered_charge = apply(&dielectric.stiffness, &voltage_v, cancelled)?;
        let mut nodal_charge_c = self.nodal_charge_c.clone();
        let mut nodal_electrode_charge_c = vec![None; voltage_v.len()];
        let mut maximum_free_charge_residual_c: f64 = 0.0;
        let mut electrode_charge_residual_c = 0.0;
        let mut supplied_work_j = 0.0;
        let reference = self.voltage_v[self.voltage_vertices[0]];
        for (v, &charge) in recovered_charge.iter().enumerate() {
            poll(cancelled)?;
            let delta = finite(charge - self.nodal_charge_c[v])?;
            if unchanged && delta.abs() > controls.charge_tolerance_c {
                return Err(Error::Balance(
                    "unchanged dielectric no longer represents its retained nodal charge",
                ));
            }
            if prescribed[v].is_some() {
                // An unchanged Hodge has no electrode impulse, including after
                // a prior solve retained free charge within its declared gate.
                let delta = if unchanged { 0.0 } else { delta };
                if !unchanged {
                    nodal_charge_c[v] = charge;
                }
                nodal_electrode_charge_c[v] = Some(delta);
                electrode_charge_residual_c = finite(electrode_charge_residual_c + delta)?;
                supplied_work_j =
                    finite((voltage_v[v] - reference).mul_add(delta, supplied_work_j))?;
            } else {
                maximum_free_charge_residual_c = maximum_free_charge_residual_c.max(delta.abs());
            }
        }
        let energy_residual_j =
            finite(endpoint.energy_j - self.stored_energy_j - supplied_work_j - parameter_work_j)?;
        if maximum_free_charge_residual_c > controls.charge_tolerance_c
            || electrode_charge_residual_c.abs() > controls.charge_tolerance_c
            || energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(Error::Balance(
                "material replacement charge, electrode work and field storage do not close",
            ));
        }
        if supplied_work_j.abs() > controls.maximum_absolute_work_j
            || parameter_work_j.abs() > controls.maximum_absolute_work_j
        {
            return Err(Error::Invalid(
                "material replacement exceeds its external work allowance",
            ));
        }
        // The retained charge, not a reconstructed replacement, is the next
        // time step's canonical initial state. Its dual pairing still has to
        // agree with the original native cell storage independently.
        check_storage(
            &voltage_v,
            &nodal_charge_c,
            endpoint.energy_j,
            controls.charge_tolerance_c,
            controls.energy_tolerance_j,
        )?;
        let mut material_state_identities = Vec::with_capacity(electric_cells.len());
        for (old, new) in self
            .conduction
            .conductors()
            .iter()
            .zip(conduction.conductors())
        {
            poll(cancelled)?;
            material_state_identities
                .push([old.material().identity().0, new.material().identity().0]);
        }
        let report = DielectricMaterialRebind {
            time_s: self.time_s,
            stored_energy_j: [self.stored_energy_j, endpoint.energy_j],
            supplied_work_j,
            parameter_work_j,
            nodal_electrode_charge_c,
            maximum_free_charge_residual_c,
            electrode_charge_residual_c,
            energy_residual_j,
            electric_cells,
            material_state_identities,
            linear,
        };
        poll(cancelled)?;
        Ok((
            PlanarLeakyDielectric {
                dielectric: Arc::new(dielectric),
                thermal_laws: None,
                conduction: Arc::new(conduction),
                voltage_v,
                nodal_charge_c,
                cell_energy_j: endpoint.cell_energy_j,
                stored_energy_j: endpoint.energy_j,
                voltage_vertices: self.voltage_vertices.clone(),
                time_s: self.time_s,
            },
            report,
        ))
    }
}
