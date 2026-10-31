//! Algebraic elimination of one explicit same-clock terminal port group.
//! Propagation, rational memory and port work stay with their original owners.

use super::*;
use fs_la::prepared_lu::LuWorkspace;

const MAX_TERMINAL_CONDITION: f64 = 1e10;

#[derive(Debug, Clone)]
pub(super) enum TerminalSolve {
    // Preserve the original scalar division, rather than multiply by 1/Z.
    Scalar(f64),
    // Cold inverse columns come from the existing normalized prepared LU.
    Matrix { inverse: Vec<f64>, dimension: usize },
}

impl TerminalSolve {
    pub(super) fn new(matrix: &[f64], dimension: usize) -> Result<Self, AcousticRealizeError> {
        if dimension == 1 {
            let denominator = matrix[0];
            if !denominator.is_finite() || denominator <= 0.0 {
                return Err(invalid("DG terminal mean solve requires positive finite combined outlet impedance"));
            }
            return Ok(Self::Scalar(denominator));
        }
        let norm = matrix.chunks_exact(dimension)
            .map(|row| row.iter().map(|value| value.abs()).sum::<f64>()).fold(0.0, f64::max);
        if !norm.is_finite() || norm <= 0.0 || matrix.iter().any(|value| !value.is_finite()) {
            return Err(invalid("DG terminal group matrix must be finite and nonzero"));
        }
        let normalized: Vec<_> = matrix.iter().map(|value| value / norm).collect();
        let mut lu = LuWorkspace::new(dimension)
            .map_err(|_| invalid("DG terminal group solve storage refused"))?;
        lu.factor(&normalized).map_err(|_| invalid("DG terminal group matrix is singular or nonfinite"))?;
        let mut inverse = vec![0.0; dimension * dimension];
        let mut rhs = vec![0.0; dimension];
        let mut column = vec![0.0; dimension];
        for col in 0..dimension {
            rhs.fill(0.0);
            rhs[col] = 1.0;
            lu.solve_factored_into(&rhs, &mut column)
                .map_err(|_| invalid("DG terminal group inverse solve left the finite domain"))?;
            for row in 0..dimension { inverse[row * dimension + col] = column[row]; }
        }
        let inverse_norm = inverse.chunks_exact(dimension)
            .map(|row| row.iter().map(|value| value.abs()).sum::<f64>()).fold(0.0, f64::max);
        let normalized_norm = normalized.chunks_exact(dimension)
            .map(|row| row.iter().map(|value| value.abs()).sum::<f64>()).fold(0.0, f64::max);
        let condition = normalized_norm * inverse_norm;
        if !condition.is_finite() || condition > MAX_TERMINAL_CONDITION {
            return Err(invalid("DG terminal group infinity-norm condition exceeds 1e10"));
        }
        for value in &mut inverse { *value /= norm; }
        if inverse.iter().any(|value| !value.is_finite()) {
            return Err(invalid("DG terminal group inverse is not representable"));
        }
        Ok(Self::Matrix { inverse, dimension })
    }

    pub(super) fn apply(&self, rhs: &[f64]) -> Result<Vec<f64>, AcousticRealizeError> {
        let values = match self {
            Self::Scalar(denominator) => vec![rhs[0] / denominator],
            Self::Matrix { inverse, dimension } => inverse.chunks_exact(*dimension)
                .map(|row| row.iter().zip(rhs).map(|(a, b)| a * b).sum()).collect(),
        };
        if values.iter().any(|value| !value.is_finite()) {
            return Err(invalid("DG terminal group mean solve left the finite flow domain"));
        }
        Ok(values)
    }
}

/// Constant coupling matrices; no accepted physical history is stored here.
#[derive(Debug, Clone)]
pub(super) struct TerminalCoupling {
    pub(super) ports: Vec<PortId>,
    pub(super) inlet_index: usize,
    pub(super) indices: Vec<usize>,
    pub(super) outward_per_inlet: Vec<f64>,
    pub(super) driven_mouth: bool,
    pub(super) mouth_count: usize,
    impedance_pa_s_m3: f64,
    active: Vec<usize>,
    solve: Option<TerminalSolve>,
}

pub(super) struct TerminalRelation {
    pub(super) inlet: DgPortResponse,
    pub(super) outward_history_m3_s: Vec<f64>,
    pub(super) response: DgMultiportResponse,
}

impl TerminalCoupling {
    pub(super) fn new(
        response: &DgMultiportResponse,
        inlet: PortId,
        ports: Vec<PortId>,
        terminal: &PreparedMimoLoad,
    ) -> Result<Self, AcousticRealizeError> {
        Self::with_channels(response, inlet, ports, terminal, false, &[])
    }

    /// Channel zero is a driven exterior mouth with outward flow -Qin.
    /// Terminal channel order and original external DG port identities survive.
    pub(super) fn with_driven_mouth(
        response: &DgMultiportResponse, inlet: PortId, ports: Vec<PortId>,
        terminal: &PreparedMimoLoad, channels: &[DgTerminalChannelLoad],
    ) -> Result<Self, AcousticRealizeError> {
        Self::with_channels(response, inlet, ports, terminal, true, channels)
    }

    pub(super) fn with_channels(
        response: &DgMultiportResponse, inlet: PortId, ports: Vec<PortId>,
        terminal: &PreparedMimoLoad, driven_mouth: bool,
        terminal_channels: &[DgTerminalChannelLoad],
    ) -> Result<Self, AcousticRealizeError> {
        let count = response.ports.len();
        let n = ports.len();
        let offset = if driven_mouth {
            terminal.ports().checked_sub(n).filter(|count| *count > 0)
                .ok_or_else(|| invalid("joint terminal source has no leading mouth channels"))?
        } else { 0 };
        let channels = n + offset;
        if n == 0 || channels != terminal.ports() || count != n + 1 {
            return Err(invalid("DG terminal group must own every non-inlet port with exactly one load channel"));
        }
        let inlet_index = response.ports.iter().position(|id| *id == inlet)
            .ok_or_else(|| invalid("DG inlet missing from its original port basis"))?;
        let mut indices = Vec::with_capacity(n);
        for (i, id) in ports.iter().enumerate() {
            if *id == inlet || ports[..i].contains(id) {
                return Err(invalid("DG terminal ports must be distinct and exclude the source inlet"));
            }
            indices.push(response.ports.iter().position(|port| port == id)
                .ok_or_else(|| invalid("DG terminal group references an undeclared port"))?);
        }
        let mut series = vec![0.0; n];
        let mut open = vec![true; n];
        for channel in terminal_channels {
            let index = ports.iter().position(|port| *port == channel.input().port)
                .ok_or_else(|| invalid("DG terminal channel has no original exterior row"))?;
            match channel.pressure_load()? {
                Some(load) => series[index] = load.impedance_pa_s_m3,
                None => open[index] = false,
            }
        }
        let active: Vec<_> = (0..n).filter(|index| open[*index]).collect();
        let dimension = active.len();
        let z = &response.impedance_pa_s_m3;
        let direct = terminal.mean_direct_term();
        let mut matrix = vec![0.0; dimension * dimension];
        for (row, &original_row) in active.iter().enumerate() {
            for (col, &original_col) in active.iter().enumerate() {
                let value = z[count * indices[original_row] + indices[original_col]]
                    + direct[(original_row + offset) * channels + original_col + offset];
                matrix[row * dimension + col] = if row == col && series[original_row] != 0.0 {
                    value + series[original_row]
                } else { value };
            }
        }
        let solve = if dimension == 0 { None }
            else { Some(TerminalSolve::new(&matrix, dimension)?) };
        let rhs: Vec<_> = active.iter().map(|&i| {
            let mut value = z[count * indices[i] + inlet_index];
            if offset == 1 { value += direct[(i + 1) * channels]; }
            value
        }).collect();
        let mut outward_per_inlet = vec![0.0; n];
        if let Some(solve) = &solve {
            for (&index, value) in active.iter().zip(solve.apply(&rhs)?) {
                outward_per_inlet[index] = value;
            }
        }
        let mut impedance_pa_s_m3 = z[count * inlet_index + inlet_index];
        if offset == 1 { impedance_pa_s_m3 += direct[0]; }
        for (i, &column) in indices.iter().enumerate() {
            let mut transfer = z[count * inlet_index + column];
            if offset == 1 { transfer += direct[i + 1]; }
            impedance_pa_s_m3 -= transfer * outward_per_inlet[i];
        }
        if !impedance_pa_s_m3.is_finite() || impedance_pa_s_m3 <= 0.0 {
            return Err(invalid("DG terminal elimination left the finite positive inlet-load domain"));
        }
        Ok(Self { ports, inlet_index, indices, outward_per_inlet, driven_mouth,
            mouth_count: offset,
            impedance_pa_s_m3, active, solve })
    }

    pub(super) fn relation(
        &self, response: DgMultiportResponse, terminal: &PreparedMimoLoad,
        channels: &[DgTerminalChannelLoad],
    ) -> Result<TerminalRelation, AcousticRealizeError> {
        let count = response.ports.len();
        let offset = self.mouth_count;
        let mut load_history = vec![0.0; self.ports.len() + offset];
        terminal.mean_history_output_into(&mut load_history)
            .map_err(|_| invalid("DG terminal mean history is outside its admitted finite domain"))?;
        let mut rhs = Vec::with_capacity(self.active.len());
        for &index in &self.active {
            let mut value = response.history_pressure_pa[self.indices[index]] - load_history[index + offset];
            if let Some(channel) = channels.iter().find(|channel| channel.input().port == self.ports[index]) {
                let load = channel.pressure_load()?
                    .ok_or_else(|| invalid("closed DG channel appeared in the active flow solve"))?;
                value += load.history_pressure_pa;
            }
            rhs.push(value);
        }
        let mut outward_history_m3_s = vec![0.0; self.ports.len()];
        if let Some(solve) = &self.solve {
            for (&index, value) in self.active.iter().zip(solve.apply(&rhs)?) {
                outward_history_m3_s[index] = value;
            }
        }
        let mut history_pressure_pa = response.history_pressure_pa[self.inlet_index];
        if offset == 1 { history_pressure_pa -= load_history[0]; }
        let direct = terminal.mean_direct_term();
        for (i, &column) in self.indices.iter().enumerate() {
            let mut transfer = response.impedance_pa_s_m3[count * self.inlet_index + column];
            if offset == 1 { transfer += direct[i + 1]; }
            history_pressure_pa -= transfer * outward_history_m3_s[i];
        }
        if !history_pressure_pa.is_finite() {
            return Err(invalid("DG terminal elimination left the finite inlet-pressure domain"));
        }
        Ok(TerminalRelation {
            inlet: DgPortResponse { history_pressure_pa, impedance_pa_s_m3: self.impedance_pa_s_m3 },
            outward_history_m3_s, response,
        })
    }
}


/// Independent physical source coordinates after the SAME terminal solve.
pub(super) struct IndependentMouthsRelation {
    pub(super) efforts: DgIndependentMouthsResponse,
    pub(super) outward_history_m3_s: Vec<f64>,
    pub(super) outward_per_source: Vec<Vec<f64>>,
}

impl TerminalCoupling {
    pub(super) fn independent_mouths_relation(
        &self, response: DgMultiportResponse, terminal: &PreparedMimoLoad,
        channels: &[DgTerminalChannelLoad],
    ) -> Result<IndependentMouthsRelation, AcousticRealizeError> {
        if !self.driven_mouth {
            return Err(invalid("independent mouth response requires the original leading exterior channel"));
        }
        // Reuse the exact accepted history elimination, including branch
        // condensation, finite-gap memories and sealed-channel active rows.
        let original = self.relation(response, terminal, channels)?;
        let response = &original.response;
        let count = response.ports.len();
        let n = self.ports.len();
        let offset = self.mouth_count;
        let width = n + offset;
        let sources = offset + 1;
        let z = &response.impedance_pa_s_m3;
        let direct = terminal.mean_direct_term();
        let mut history = vec![0.0; width];
        terminal.mean_history_output_into(&mut history)
            .map_err(|_| invalid("independent exterior history left the finite domain"))?;
        let mut outward_per_source = vec![vec![0.0; sources]; n];
        for coordinate in 0..sources {
            let rhs: Vec<_> = self.active.iter().map(|&i| {
                if coordinate == 0 {
                    z[count * self.indices[i] + self.inlet_index]
                } else {
                    -direct[(i + offset) * width + coordinate - 1]
                }
            }).collect();
            if let Some(solve) = &self.solve {
                for (&index, value) in self.active.iter().zip(solve.apply(&rhs)?) {
                    outward_per_source[index][coordinate] = value;
                }
            }
        }
        let mut efforts = DgIndependentMouthsResponse {
            history_pressure_pa: vec![0.0; sources],
            impedance_pa_s_m3: vec![0.0; sources * sources],
        };
        efforts.history_pressure_pa[0] = response.history_pressure_pa[self.inlet_index];
        efforts.impedance_pa_s_m3[0] = z[count * self.inlet_index + self.inlet_index];
        for mouth in 0..offset {
            efforts.history_pressure_pa[mouth + 1] = history[mouth];
            for other in 0..offset {
                efforts.impedance_pa_s_m3[(mouth + 1) * sources + other + 1] =
                    direct[mouth * width + other];
            }
        }
        for (i, &column) in self.indices.iter().enumerate() {
            let bore_transfer = z[count * self.inlet_index + column];
            let q = original.outward_history_m3_s[i];
            efforts.history_pressure_pa[0] -= bore_transfer * q;
            for coordinate in 0..sources {
                let transfer = outward_per_source[i][coordinate];
                efforts.impedance_pa_s_m3[coordinate] -= bore_transfer * transfer;
            }
            for mouth in 0..offset {
                let mouth_transfer = direct[mouth * width + i + offset];
                efforts.history_pressure_pa[mouth + 1] += mouth_transfer * q;
                for coordinate in 0..sources {
                    efforts.impedance_pa_s_m3[(mouth + 1) * sources + coordinate] +=
                        mouth_transfer * outward_per_source[i][coordinate];
                }
            }
        }
        if !efforts.history_pressure_pa.iter()
            .chain(efforts.impedance_pa_s_m3.iter())
            .chain(outward_per_source.iter().flatten()).all(|value| value.is_finite())
        {
            return Err(invalid("independent mouth elimination left the finite effort domain"));
        }
        Ok(IndependentMouthsRelation { efforts,
            outward_history_m3_s: original.outward_history_m3_s, outward_per_source })
    }
}
