//! Eliminate only the algebraic junction flows; original DG states stay intact.
use super::*;
use crate::dg_port::terminal::TerminalSolve;

#[derive(Debug, Clone)]
pub(in crate::dg_port) struct BranchCoupling {
    pub(super) main_indices: Vec<usize>,
    pub(super) branch_inlets: Vec<usize>,
    pub(super) branch_outlets: Vec<Option<usize>>,
    external_indices: Vec<usize>,
    ports: Vec<PortId>,
    solve: TerminalSolve,
    // x = A^-1(H_main,j - H_branch,in) + transfer * q_external.
    transfer: Vec<f64>,
    // p_external = H_external + Z_base q_external + output * x.
    output: Vec<f64>,
    impedance: Vec<f64>,
}

impl BranchCoupling {
    pub(in crate::dg_port) fn new(
        main: &DgMultiportResponse,
        branches: &[DgBranchLoad],
    ) -> Result<Self, AcousticRealizeError> {
        let original_count = main.ports.len();
        let external_indices: Vec<_> = main
            .ports
            .iter()
            .enumerate()
            .filter_map(|(index, port)| {
                (!branches
                    .iter()
                    .any(|branch| branch.port == *port && branch.outlet.is_none()))
                .then_some(index)
            })
            .collect();
        let n = external_indices.len();
        let b = branches.len();
        let responses = responses(branches)?;
        let mut main_indices = Vec::with_capacity(b);
        let mut branch_inlets = Vec::with_capacity(b);
        let mut branch_outlets = Vec::with_capacity(b);
        for (branch, response) in branches.iter().zip(&responses) {
            main_indices.push(main.ports.iter().position(|id| *id == branch.port)
                .ok_or_else(|| invalid("DG branch junction is missing from the original bore basis"))?);
            branch_inlets.push(response.ports.iter().position(|id| *id == branch.inlet)
                .ok_or_else(|| invalid("DG branch omitted its original inlet"))?);
            branch_outlets.push(branch.outlet.map(|outlet|
                response.ports.iter().position(|id| *id == outlet)
                    .ok_or_else(|| invalid("DG branch omitted its original outlet"))).transpose()?);
        }
        let mut matrix = vec![0.0; b * b];
        let mut input = vec![0.0; b * n];
        let mut output = vec![0.0; n * b];
        let mut impedance = vec![0.0; n * n];
        for row in 0..b {
            let branch = &responses[row];
            let width = branch.ports.len();
            let inlet = branch_inlets[row];
            let outlet = branch_outlets[row];
            for col in 0..b {
                matrix[row * b + col] = main.impedance_pa_s_m3[
                    main_indices[row] * original_count + main_indices[col]];
                if row == col {
                    matrix[row * b + col] += branch.impedance_pa_s_m3[width * inlet + inlet];
                }
            }
            for col in 0..n {
                input[row * n + col] = match main_indices.iter().position(|&i| i == external_indices[col]) {
                    Some(index) if index == row => -branch.impedance_pa_s_m3[
                        width * inlet + outlet.ok_or_else(|| invalid("sealed DG branch acquired an exterior column"))?],
                    Some(_) => 0.0,
                    None => main.impedance_pa_s_m3[main_indices[row] * original_count + external_indices[col]],
                };
            }
        }
        for row in 0..n {
            if let Some(index) = main_indices.iter().position(|&i| i == external_indices[row]) {
                let branch = &responses[index];
                let width = branch.ports.len();
                let inlet = branch_inlets[index];
                let outlet = branch_outlets[index]
                    .ok_or_else(|| invalid("sealed DG branch acquired an exterior row"))?;
                output[row * b + index] = branch.impedance_pa_s_m3[width * outlet + inlet];
                impedance[row * n + row] = branch.impedance_pa_s_m3[width * outlet + outlet];
            } else {
                for col in 0..b {
                    output[row * b + col] =
                        -main.impedance_pa_s_m3[external_indices[row] * original_count + main_indices[col]];
                }
                for col in 0..n {
                    if !main_indices.contains(&external_indices[col]) {
                        impedance[row * n + col] = main.impedance_pa_s_m3[
                            external_indices[row] * original_count + external_indices[col]];
                    }
                }
            }
        }
        let solve = TerminalSolve::new(&matrix, b)?;
        let mut transfer = vec![0.0; b * n];
        for col in 0..n {
            let rhs: Vec<_> = (0..b).map(|row| input[row * n + col]).collect();
            let values = solve.apply(&rhs)?;
            for row in 0..b { transfer[row * n + col] = values[row]; }
        }
        for row in 0..n {
            for col in 0..n {
                for junction in 0..b {
                    impedance[row * n + col] +=
                        output[row * b + junction] * transfer[junction * n + col];
                }
            }
        }
        if impedance.iter().any(|z| !z.is_finite())
            || (0..n).any(|i| impedance[i * n + i] < 0.0)
        {
            return Err(invalid("DG branch condensation left the finite passive mean-port domain"));
        }
        Ok(Self {
            ports: external_indices.iter().map(|&index| main.ports[index]).collect(),
            main_indices, branch_inlets, branch_outlets, external_indices,
            solve, transfer, output, impedance,
        })
    }

    fn junction_history(
        &self, main: &DgMultiportResponse, branches: &[DgMultiportResponse],
    ) -> Result<Vec<f64>, AcousticRealizeError> {
        let rhs: Vec<_> = self.main_indices.iter().enumerate().map(|(i, &slot)|
            main.history_pressure_pa[slot] - branches[i].history_pressure_pa[self.branch_inlets[i]]
        ).collect();
        self.solve.apply(&rhs)
    }

    pub(in crate::dg_port) fn response(
        &self,
        main: DgMultiportResponse,
        branches: &[DgBranchLoad],
    ) -> Result<DgMultiportResponse, AcousticRealizeError> {
        let branches = responses(branches)?;
        let history = self.junction_history(&main, &branches)?;
        let b = branches.len();
        let mut pressure: Vec<_> = self
            .external_indices
            .iter()
            .map(|&original| {
                if let Some(index) = self.main_indices.iter().position(|&slot| slot == original) {
                    branches[index].history_pressure_pa[self.branch_outlets[index]
                        .expect("only unsealed branches have exterior rows")]
                } else {
                    main.history_pressure_pa[original]
                }
            })
            .collect();
        for (row, value) in pressure.iter_mut().enumerate() {
            for junction in 0..b {
                *value += self.output[row * b + junction] * history[junction];
            }
        }
        if pressure.iter().any(|value| !value.is_finite()) {
            return Err(invalid("DG branch history left the finite mean-pressure domain"));
        }
        Ok(DgMultiportResponse {
            ports: self.ports.clone(), history_pressure_pa: pressure,
            impedance_pa_s_m3: self.impedance.clone(),
        })
    }

    pub(in crate::dg_port) fn flows(
        &self, main: &DgMultiportResponse, branches: &[DgMultiportResponse],
        external: &[PortInput],
    ) -> Result<Vec<f64>, AcousticRealizeError> {
        let n = self.ports.len();
        let mut flows = self.junction_history(main, branches)?;
        for (junction, flow) in flows.iter_mut().enumerate() {
            for input in external {
                let slot = self.ports.iter().position(|id| *id == input.id)
                    .ok_or_else(|| invalid("DG branch received an undeclared external port"))?;
                *flow += self.transfer[junction * n + slot] * input.effort;
            }
        }
        if flows.iter().any(|value| !value.is_finite()) {
            return Err(invalid("DG branch junction solve left the finite flow domain"));
        }
        Ok(flows)
    }
}

pub(super) fn responses(
    branches: &[DgBranchLoad],
) -> Result<Vec<DgMultiportResponse>, AcousticRealizeError> {
    branches.iter().map(|branch| branch.propagator.multiport_response(&branch.state)
        .map_err(AcousticRealizeError::Duct)).collect()
}
