//! Instantaneous, externally supplied material replacement in the same circuit.
//! Flux and charge remain the canonical differential state. A changed magnetic
//! Hamiltonian has signed parameter work; it is neither Joule heat nor a native
//! source-port impulse. Continuous thermal evolution needs a caloric owner.

use super::*;

/// Independent work admission for an explicit material replacement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WindingMaterialRebindControls {
    /// Maximum absolute externally supplied/received parameter work [J].
    /// Zero admits an exactly zero-work change, including resistance-only swaps.
    pub maximum_absolute_work_j: f64,
    /// Agreement of the canonical work and native storage difference [J].
    pub energy_tolerance_j: f64,
}

/// Zero-duration replacement at the original accepted clock.
/// Every two-entry array is ordered [before, after]. No heat or time is advanced.
#[derive(Clone, Debug, PartialEq)]
pub struct WindingMaterialRebind {
    pub time_s: f64,
    /// Native magnetic/capacitor storage PLUS the declared zero-current field
    /// reference W0. Ordinary circuit steps report incremental storage instead.
    pub stored_energy_j: [f64; 2],
    /// Original zero-current field references [J], in the source model's
    /// ZeroFlux or FixedRemanence convention. Not material formation energies.
    pub bias_energy_j: [f64; 2],
    /// W = W0_after - W0_before + delta_lambda_after^T L_after^-1
    /// delta_lambda_after / 2 - delta_lambda_before^T L_before^-1
    /// delta_lambda_before / 2, at conserved TOTAL physical linkage.
    /// Positive work is supplied to the circuit by the replacement owner;
    /// negative work must be received by that owner, never discarded as heat.
    pub supplied_work_j: f64,
    /// Independently checked incremental-storage/work defect [J]. Evaluating
    /// this before adding W0 avoids masking electrical work by a large bias.
    pub energy_residual_j: f64,
    /// Winding order, including mutual inductance and declared end storage.
    pub currents_a: [Vec<f64>; 2],
    pub resistance_ohm: [Vec<f64>; 2],
    /// Original triangle order, including passive magnetic regions.
    pub material_state_identities: Vec<[[u8; 32]; 2]>,
    /// Original conducting-triangle order, at each actual electrical snapshot.
    pub conductor_state_identities: Vec<[[u8; 32]; 2]>,
}

fn same_fixed_model(
    before: &WindingCircuit<'_>,
    after: &WindingCircuit<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), OhmicWindingError> {
    let invalid = || {
        OhmicWindingError::Invalid(
            "material replacement requires the same mesh, boundaries, winding profiles, end model and circuit topology",
        )
    };
    let a = before.basis();
    let b = after.basis();
    let ma = a.problem().mesh();
    let mb = b.problem().mesh();
    if ma.metric() != mb.metric()
        || ma.vertices().len() != mb.vertices().len()
        || ma.faces().len() != mb.faces().len()
        || a.prescribed_vertices() != b.prescribed_vertices()
        || a.winding_count() != b.winding_count()
        || a.end_inductance_h() != b.end_inductance_h()
        || before.branch_indices() != after.branch_indices()
        || before.flux_coordinates() != after.flux_coordinates()
        || before.graph().node_count != after.graph().node_count
        || before.graph().branches.len() != after.graph().branches.len()
        || before.graph().transformers != after.graph().transformers
    {
        return Err(invalid());
    }
    for (x, y) in ma.vertices().iter().zip(mb.vertices()) {
        poll(cancelled)?;
        if x != y {
            return Err(invalid());
        }
    }
    for (x, y) in ma.faces().iter().zip(mb.faces()) {
        poll(cancelled)?;
        if x != y {
            return Err(invalid());
        }
    }
    for (x, y) in a.profiles_per_m2().iter().zip(b.profiles_per_m2()) {
        for (u, v) in x.iter().zip(y) {
            poll(cancelled)?;
            if u != v {
                return Err(invalid());
            }
        }
    }
    for (&(a0, a1, x), &(b0, b1, y)) in before.graph().branches.iter().zip(&after.graph().branches)
    {
        poll(cancelled)?;
        // Only field-derived inductance changes. Capacitors, graph resistors,
        // source/transformer conventions and every branch orientation stay fixed.
        if (a0, a1) != (b0, b1)
            || !(matches!((x, y), (Branch::Inductor { .. }, Branch::Inductor { .. })) || x == y)
        {
            return Err(invalid());
        }
    }
    Ok(())
}

impl OhmicWindingCircuit<'_, '_> {
    /// Substitute a newly source-resolved magnetic/electrical model at the SAME
    /// accepted instant. Reference geometry, magnetic boundaries, signed winding
    /// profiles, end inductances, graph topology and nonmagnetic storage must
    /// match. New material cards, selected electrical basis and state are data.
    ///
    /// The original native consistency solver recomputes only algebraic values;
    /// every capacitor charge is preserved bit-for-bit. Total physical flux
    /// linkage is conserved: a changed permanent-flux offset translates native
    /// delta_lambda by lambda0_before - lambda0_after. Unchanged offsets retain
    /// incremental linkage bit-for-bit. If the
    /// new currents cannot obey the original ideal constraints at those values,
    /// replacement refuses rather than inventing a source impulse. No step or
    /// heat is generated. Later evolution consumes the new actual field basis
    /// and conductor resistance through the original midpoint/heat owner.
    ///
    /// The explicit work is the magnetic Hamiltonian change at fixed TOTAL
    /// lambda, including full mutual coupling and the changed zero-current
    /// field reference W0. Actual endpoint fields must remain in their source
    /// recoil envelopes. This is an externally supplied material-change
    /// operation, not continuous heating, motion, hysteresis or a prediction
    /// of material formation energy. No path through intermediate magnet states
    /// or irreversible demagnetization is inferred.
    pub fn rebind_material_candidate<'c, 'm>(
        &self,
        magnetic: &'c WindingCircuit<'m>,
        sources: Vec<WindingConductor>,
        ports: &[f64],
        controls: WindingMaterialRebindControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(OhmicWindingCircuit<'c, 'm>, WindingMaterialRebind), OhmicWindingError> {
        poll(cancelled)?;
        if !controls.maximum_absolute_work_j.is_finite()
            || controls.maximum_absolute_work_j < 0.0
            || !positive(controls.energy_tolerance_j)
            || ports.len() != magnetic.dae().ports
            || ports.iter().any(|v| !v.is_finite())
        {
            return Err(OhmicWindingError::Invalid(
                "finite work bound, positive energy tolerance and complete finite ports required",
            ));
        }
        same_fixed_model(self.magnetic, magnetic, cancelled)?;
        // A new explicit source binding permits data-only material/basis swaps;
        // normal temperature refresh still pins its original source selection.
        let (conductors, conductor_windings, resistance_ohm) =
            bind(magnetic, sources, None, cancelled)?;
        let dae = lower(magnetic, &resistance_ohm)?;
        let old_dae = self.magnetic.dae();
        let nd = dae.system.differential_dim();
        if dae.system.state_dim() != self.state.len()
            || nd != old_dae.system.differential_dim()
            || dae.flux_index != old_dae.flux_index
            || dae.charge_index != old_dae.charge_index
            || dae.node_potential_index != old_dae.node_potential_index
            || conductor_windings != self.conductor_windings
        {
            return Err(OhmicWindingError::Invalid(
                "canonical circuit/material coordinates changed",
            ));
        }
        let mut seed = vec![0.0; self.state.len()];
        seed[..nd].copy_from_slice(&self.state[..nd]);
        let offsets: Vec<_> = self
            .magnetic
            .basis()
            .bias_flux_linkages_wb()
            .iter()
            .zip(magnetic.basis().bias_flux_linkages_wb())
            .map(|(before, after)| before - after)
            .collect();
        for (&coordinate, &offset) in magnetic.flux_coordinates().iter().zip(&offsets) {
            poll(cancelled)?;
            // Subtract offsets first: identical bias must not erase a small
            // live current by adding then subtracting a much larger lambda0.
            if offset != 0.0 {
                seed[coordinate] += offset;
            }
        }
        if seed.iter().chain(&offsets).any(|v| !v.is_finite()) {
            return Err(OhmicWindingError::Invalid(
                "total-flux coordinate transfer is unrepresentable",
            ));
        }
        let state = dae.consistent_initial_state(&seed, ports)?;
        if state[..nd]
            .iter()
            .zip(&seed[..nd])
            .any(|(a, b)| a.to_bits() != b.to_bits())
            || state.iter().any(|v| !v.is_finite())
        {
            return Err(OhmicWindingError::Balance(
                "material replacement changed conserved flux or charge",
            ));
        }
        self.magnetic.admit_state(&self.state, cancelled)?;
        magnetic.admit_state(&state, cancelled)?;
        let incremental_energy_j = [
            old_dae.system.hamiltonian(&self.state),
            dae.system.hamiltonian(&state),
        ];
        let bias_energy_j = [
            self.magnetic.basis().bias_energy_j(),
            magnetic.basis().bias_energy_j(),
        ];
        let stored_energy_j = [
            incremental_energy_j[0] + bias_energy_j[0],
            incremental_energy_j[1] + bias_energy_j[1],
        ];
        let before = self.magnetic.basis().inverse_inductance_per_h();
        let after = magnetic.basis().inverse_inductance_per_h();
        let n = magnetic.basis().winding_count();
        let mut incremental_work_j = 0.0;
        for (i, &row) in magnetic.flux_coordinates().iter().enumerate() {
            poll(cancelled)?;
            for (j, &col) in magnetic.flux_coordinates().iter().enumerate() {
                incremental_work_j += (0.5 * self.state[row])
                    * (after[i * n + j] - before[i * n + j])
                    * self.state[col];
                incremental_work_j +=
                    (self.state[row] + 0.5 * offsets[i]) * after[i * n + j] * offsets[j];
            }
        }
        let supplied_work_j = incremental_work_j + (bias_energy_j[1] - bias_energy_j[0]);
        let energy_residual_j =
            incremental_energy_j[1] - incremental_energy_j[0] - incremental_work_j;
        if stored_energy_j.iter().any(|v| !v.is_finite() || *v < 0.0)
            || incremental_energy_j
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0)
            || !incremental_work_j.is_finite()
            || !supplied_work_j.is_finite()
            || !energy_residual_j.is_finite()
            || energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(OhmicWindingError::Balance(
                "replacement work and native field/capacitor storage disagree",
            ));
        }
        if supplied_work_j.abs() > controls.maximum_absolute_work_j {
            return Err(OhmicWindingError::Invalid(
                "material replacement exceeds its external work allowance",
            ));
        }
        let currents_a = [self.currents_a()?, magnetic.currents_a(&state)?];
        let mut material_state_identities =
            Vec::with_capacity(magnetic.basis().problem().materials().len());
        for (a, b) in self
            .magnetic
            .basis()
            .problem()
            .materials()
            .iter()
            .zip(magnetic.basis().problem().materials())
        {
            poll(cancelled)?;
            material_state_identities.push([a.identity().0, b.identity().0]);
        }
        let mut conductor_state_identities = Vec::with_capacity(conductors.len());
        for (a, b) in self.conductors.iter().zip(&conductors) {
            poll(cancelled)?;
            conductor_state_identities.push([
                a.conductor.material().identity().0,
                b.conductor.material().identity().0,
            ]);
        }
        let report = WindingMaterialRebind {
            time_s: self.time_s,
            stored_energy_j,
            bias_energy_j,
            supplied_work_j,
            energy_residual_j,
            currents_a,
            resistance_ohm: [self.resistance_ohm.clone(), resistance_ohm.clone()],
            material_state_identities,
            conductor_state_identities,
        };
        poll(cancelled)?;
        Ok((
            OhmicWindingCircuit {
                magnetic,
                conductors,
                conductor_windings,
                resistance_ohm,
                state,
                time_s: self.time_s,
            },
            report,
        ))
    }
}
