//! Close a physical memoryless load in the original continuous generator.
use crate::linear_exact::positive_heat::factor_resistance;
use crate::{PhsError, PortHamiltonian, ResistanceWitness, require_symmetric};

impl PortHamiltonian {
    /// Add the passive law `u_selected = u_external - D*y_selected`.
    ///
    /// `D` is a complete row-major conductance in `selected` order. It adds
    /// `G_selected D G_selected^T` to the original resistance and dissipates
    /// `y_selected^T D y_selected`. All states, storage, port columns and their
    /// order remain unchanged; the retained ports accept the additive external
    /// input. Thus a dynamic load can be interconnected afterwards and the
    /// complete coupled system propagated in one exact flow.
    ///
    /// Preparation retains positive loss factors, including any original
    /// resistance witness. No load state, lag, timestep or frequency fit is
    /// introduced. This is algebraic admission of the declared physical law;
    /// its units, source validity and applicability belong to the caller.
    ///
    /// # Errors
    /// Duplicate/out-of-range ports, malformed/nonfinite or asymmetric data,
    /// a negative/unresolved resistance direction, overflow, or failure of
    /// componentwise factor reconstruction refuse without clipping.
    pub fn with_port_conductance(
        self,
        selected: &[usize],
        conductance: &[f64],
    ) -> Result<Self, PhsError> {
        self.admit_port_conductance(selected, conductance)?;
        if conductance.iter().all(|value| *value == 0.0) {
            return Ok(self);
        }
        let root = factor_resistance(conductance, selected.len())?;
        let witness =
            ResistanceWitness::try_from_source_rows(selected.len(), conductance.to_vec(), root)?;
        self.with_port_conductance_witness(selected, witness)
    }

    /// Same physical load, retaining its original positive source squares.
    ///
    /// Use this when a source assembles a low-rank `D`: its exact nullspace
    /// need not be rediscovered through rounded eigenvalues. The witness must
    /// describe precisely the selected-port matrix, not a state resistance.
    /// Original base loss and transformed load squares are combined through
    /// the existing checked positive QR owner, never by subtracting energies.
    ///
    /// # Errors
    /// The same port, arithmetic and reconstruction refusals as
    /// [`Self::with_port_conductance`].
    pub fn with_port_conductance_witness(
        mut self,
        selected: &[usize],
        witness: ResistanceWitness,
    ) -> Result<Self, PhsError> {
        self.admit_port_conductance(selected, &witness.matrix)?;
        if witness.matrix.iter().all(|value| *value == 0.0) {
            return Ok(self);
        }
        let (n, ports) = (self.n, selected.len());
        let square = n.checked_mul(n).ok_or(PhsError::Dimension {
            what: "loaded state resistance extent",
        })?;
        let count = ports.checked_mul(n).ok_or(PhsError::Dimension {
            what: "loaded port factor extent",
        })?;
        let rows_len = square.checked_add(count).ok_or(PhsError::Dimension {
            what: "combined source resistance factor extent",
        })?;
        if n == 0 {
            return Ok(self);
        }
        let base = match &self.resistance_witness {
            Some(original) if original.matrix == self.r => original.root.clone(),
            _ => factor_resistance(&self.r, n)?,
        };
        let mut rows = Vec::with_capacity(rows_len);
        rows.extend_from_slice(&base);
        rows.resize(rows_len, 0.0);
        let mut dg = vec![0.0; count];
        for a in 0..ports {
            for i in 0..n {
                for (b, &port) in selected.iter().enumerate() {
                    let g = self.g[i * self.m + port];
                    dg[a * n + i] += witness.matrix[a * ports + b] * g;
                    rows[square + a * n + i] += witness.root[a * ports + b] * g;
                }
            }
        }
        let mut resistance = self.r.clone();
        for i in 0..n {
            for j in 0..n {
                let mut addition = 0.0;
                for (a, &port) in selected.iter().enumerate() {
                    addition += self.g[i * self.m + port] * dg[a * n + j];
                }
                resistance[i * n + j] += addition;
            }
        }
        let combined = ResistanceWitness::try_from_source_rows(n, resistance.clone(), rows)?;
        self.r = resistance;
        self.resistance_witness = Some(combined);
        Ok(self)
    }

    fn admit_port_conductance(
        &self,
        selected: &[usize],
        conductance: &[f64],
    ) -> Result<(), PhsError> {
        if selected.len().checked_mul(selected.len()) != Some(conductance.len())
            || conductance.iter().any(|value| !value.is_finite())
        {
            return Err(PhsError::Dimension {
                what: "finite complete selected-port conductance required",
            });
        }
        for (i, &port) in selected.iter().enumerate() {
            if port >= self.m || selected[..i].contains(&port) {
                return Err(PhsError::BadPortPairing);
            }
        }
        require_symmetric(conductance, selected.len(), "port conductance D")
    }
}
