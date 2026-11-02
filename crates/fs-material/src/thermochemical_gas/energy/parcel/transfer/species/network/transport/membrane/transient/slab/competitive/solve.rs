//! Coupled nonlinear residual, structured correction and independent balances.
use super::*;

fn ratio(residual: f64, tolerance: f64) -> f64 {
    if residual == 0.0 {
        0.0
    } else {
        residual.abs() / tolerance
    }
}
fn check(residual: f64, tolerance: f64) -> Result<()> {
    if !residual.is_finite() || residual.abs() > tolerance {
        Err(MembraneDiffusionError::AmountBalance {
            residual_mol: residual,
            tolerance_mol: tolerance,
        }
        .into())
    } else {
        Ok(())
    }
}
impl CompetitiveMembraneDiffusionState {
    fn row(&self, species: usize, row: usize, boundary: Boundary<'_>) -> (f64, f64, Option<usize>) {
        let c = &self.channels[species];
        let n = self.capacities.len();
        let offset = boundary.offset();
        if let Boundary::Finite {
            pressures,
            capacities,
        } = boundary
        {
            if row == 0 {
                return (capacities[0], pressures[species][0], None);
            }
            if row == n + 1 {
                return (capacities[1], pressures[species][1], None);
            }
        }
        let i = row - offset;
        (
            c.mobile.model.cells[i].capacity_mol_m2_pa * c.mobile.model.area_m2,
            c.mobile.partial_pressures_pa[i],
            Some(i),
        )
    }
    fn evaluate(
        &self,
        boundary: Boundary<'_>,
        duration: f64,
        controls: CompetitiveMembraneSolve,
        work: &mut CompetitiveMembraneWorkspace,
        trial: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64)> {
        let ns = self.channels.len();
        let n = self.capacities.len();
        let rows = work.rows;
        let offset = boundary.offset();
        work.occupied.fill(0.0);
        work.adsorbed.fill(0.0);
        work.desorbed.fill(0.0);
        work.site_residual.fill(0.0);
        work.site_diagonal.fill(0.0);
        work.coupling_row.fill(0.0);
        work.coupling_column.fill(0.0);
        work.vacancy.fill(0.0);
        let mut site = CompetitiveSiteBuffer::default();
        let mut initial = [0.0; MAX_COMPETITIVE_SPECIES];
        let mut rates = [0.0; MAX_COMPETITIVE_SPECIES];
        let mut desorption = [0.0; MAX_COMPETITIVE_SPECIES];
        let mut alpha = [0.0; MAX_COMPETITIVE_SPECIES];
        for i in 0..n {
            poll(cancelled)?;
            if self.capacities[i] == 0.0 {
                continue;
            }
            for s in 0..ns {
                let c = &self.channels[s];
                let k = s * rows + offset + i;
                initial[s] = c.trapped[i];
                let (a, d) = c.traps[c.mobile.model.cells[i].layer]
                    .map(|t| {
                        (
                            t.affinity_pa_inverse() * t.desorption_rate_s_inverse(),
                            t.desorption_rate_s_inverse(),
                        )
                    })
                    .unwrap_or((0.0, 0.0));
                let q = if trial {
                    work.trial[k]
                } else {
                    work.current[k]
                };
                alpha[s] = a;
                desorption[s] = d;
                rates[s] = finite(a * q)?;
                if a > 0.0 && q > 0.0 && rates[s] == 0.0 {
                    return Err(invalid("competitive adsorption forcing underflow"));
                }
            }
            backward_euler_sites_into(
                &initial[..ns],
                self.capacities[i],
                &rates[..ns],
                &desorption[..ns],
                duration,
                &mut site,
                &mut *cancelled,
            )?;
            work.vacancy[offset + i] = site.vacancy_moles()[1];
            for s in 0..ns {
                let k = s * rows + offset + i;
                work.occupied[k] = site.final_moles()[s];
                work.adsorbed[k] = site.gross_adsorption_mol()[s];
                work.desorbed[k] = site.gross_desorption_mol()[s];
                work.site_residual[k] = site.kinetic_residuals_mol()[s];
                if alpha[s] > 0.0 {
                    work.site_diagonal[k] =
                        finite(site.adsorption_rate_derivative_mol_s(s, s)? * alpha[s])?;
                    if ns > 1 {
                        work.coupling_column[k] = finite(site.rate_response_mol_s(s)? * alpha[s])?;
                    }
                }
                if ns > 1 {
                    work.coupling_row[k] =
                        site.vacancy_coupling_fraction(s).expect("admitted species");
                }
            }
        }
        let mut maximum_ratio: f64 = 0.0;
        let mut maximum_residual: f64 = 0.0;
        for s in 0..ns {
            let c = &self.channels[s];
            let model = &c.mobile.model;
            let mut changes = Sum::default();
            let mut scales = Sum::default();
            for row in 0..rows {
                if row % 128 == 0 {
                    poll(cancelled)?;
                }
                let k = s * rows + row;
                let (capacity, old, cell) = self.row(s, row, boundary);
                let q = if trial {
                    work.trial[k]
                } else {
                    work.current[k]
                };
                let left = if row > 0 {
                    if trial {
                        work.trial[k - 1]
                    } else {
                        work.current[k - 1]
                    }
                } else {
                    boundary.pressures()[s][0]
                };
                let right = if row + 1 < rows {
                    if trial {
                        work.trial[k + 1]
                    } else {
                        work.current[k + 1]
                    }
                } else {
                    boundary.pressures()[s][1]
                };
                let l = if offset == 1 && row == 0 {
                    0.0
                } else {
                    finite((model.faces[row - offset].permeance * model.area_m2) * duration)?
                };
                let r = if offset == 1 && row + 1 == rows {
                    0.0
                } else {
                    finite((model.faces[row + 1 - offset].permeance * model.area_m2) * duration)?
                };
                let mobile = finite(capacity * (q - old))?;
                let trapped = cell.map(|i| work.occupied[k] - c.trapped[i]).unwrap_or(0.0);
                let incoming = finite(l * (left - q))?;
                let outgoing = finite(r * (q - right))?;
                let residual = finite((mobile + trapped) - (incoming - outgoing))?;
                // Reserve algebraic headroom for the independently sampled
                // original face transfers and grouped final species balances.
                let tolerance = 0.125
                    * controls.balance.tolerance(
                        mobile.abs() + trapped.abs() + incoming.abs() + outgoing.abs(),
                    )?;
                maximum_ratio = maximum_ratio.max(ratio(residual, tolerance));
                maximum_residual = maximum_residual.max(residual.abs());
                changes.add(mobile + trapped)?;
                scales.add(mobile.abs() + trapped.abs())?;
                let mass = finite(capacity + work.site_diagonal[k])?;
                let scale = mass.max(l).max(r);
                let (m, a, b) = (mass / scale, l / scale, r / scale);
                let diagonal = m + a + b;
                if !positive(m)
                    || diagonal <= a + b
                    || (l > 0.0 && a == 0.0)
                    || (r > 0.0 && b == 0.0)
                {
                    return Err(invalid(
                        "competitive Newton matrix lost positive storage or face coupling",
                    ));
                }
                work.diagonal[k] = diagonal;
                work.scale[k] = scale;
                work.rhs[k] = finite(-residual / scale)?;
                work.linear_allowance[k] = finite(
                    (controls.linear_relative_tolerance * residual.abs().max(tolerance)) / scale,
                )?;
                if row > 0 {
                    work.lower[k - 1] = -a;
                }
                if row + 1 < rows {
                    work.upper[k] = -b;
                }
            }
            let q0 = if trial {
                work.trial[s * rows]
            } else {
                work.current[s * rows]
            };
            let q1 = if trial {
                work.trial[(s + 1) * rows - 1]
            } else {
                work.current[(s + 1) * rows - 1]
            };
            let (incoming, outgoing) = if offset == 0 {
                (
                    finite(
                        ((model.faces[0].permeance * model.area_m2) * duration)
                            * (boundary.pressures()[s][0] - q0),
                    )?,
                    finite(
                        ((model.faces[n].permeance * model.area_m2) * duration)
                            * (q1 - boundary.pressures()[s][1]),
                    )?,
                )
            } else {
                (0.0, 0.0)
            };
            let residual = finite(changes.value - (incoming - outgoing))?;
            maximum_ratio = maximum_ratio.max(ratio(
                residual,
                0.125
                    * controls
                        .balance
                        .tolerance(scales.value + incoming.abs() + outgoing.abs())?,
            ));
            maximum_residual = maximum_residual.max(residual.abs());
        }
        Ok((maximum_ratio, maximum_residual))
    }

    pub(super) fn solve(
        &self,
        boundary: Boundary<'_>,
        duration: f64,
        controls: CompetitiveMembraneSolve,
        work: &mut CompetitiveMembraneWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompetitiveMembraneCandidate> {
        poll(&mut cancelled)?;
        controls.balance.validate()?;
        let ns = self.channels.len();
        let n = self.capacities.len();
        let rows = n + 2 * boundary.offset();
        let elapsed = self.elapsed_s() + duration;
        if !positive(duration)
            || !elapsed.is_finite()
            || elapsed <= self.elapsed_s()
            || boundary.pressures().len() != ns
            || boundary
                .pressures()
                .iter()
                .flatten()
                .any(|&p| !nonnegative(p))
            || work.rows != rows
            || work.species != ns
            || !(1..=256).contains(&controls.maximum_iterations)
            || !(1..=64).contains(&controls.maximum_backtracks)
            || !(1..=4096).contains(&controls.maximum_linear_sweeps)
            || !positive(controls.linear_relative_tolerance)
            || controls.linear_relative_tolerance > 0.5
        {
            return Err(invalid(
                "matched workspace/boundaries, positive interval and bounded competitive solve controls required",
            ));
        }
        for s in 0..ns {
            for row in 0..rows {
                if row % 128 == 0 {
                    poll(&mut cancelled)?;
                }
                work.current[s * rows + row] = self.row(s, row, boundary).1;
            }
        }
        let mut evaluations = 1;
        let mut iterations = 0;
        let mut sweeps = 0;
        let mut residual =
            self.evaluate(boundary, duration, controls, work, false, &mut cancelled)?;
        while residual.0 > 1.0 {
            poll(&mut cancelled)?;
            if iterations >= controls.maximum_iterations {
                return Err(CompetitiveMembraneError::NotConverged {
                    iterations,
                    residual_ratio: residual.0,
                });
            }
            sweeps += work.correct(controls.maximum_linear_sweeps, &mut cancelled)?;
            let mut fraction = 1.0;
            let mut accepted = None;
            for _ in 0..controls.maximum_backtracks {
                poll(&mut cancelled)?;
                let mut admissible = true;
                for k in 0..ns * rows {
                    if k % 128 == 0 {
                        poll(&mut cancelled)?;
                    }
                    work.trial[k] = work.current[k] + fraction * work.correction[k];
                    admissible &= nonnegative(work.trial[k]);
                }
                if admissible {
                    let next =
                        self.evaluate(boundary, duration, controls, work, true, &mut cancelled)?;
                    evaluations += 1;
                    if next.0 <= 1.0 || next.1 < residual.1 {
                        accepted = Some(next);
                        break;
                    }
                }
                fraction *= 0.5;
            }
            let Some(next) = accepted else {
                return Err(CompetitiveMembraneError::NotConverged {
                    iterations,
                    residual_ratio: residual.0,
                });
            };
            core::mem::swap(&mut work.current, &mut work.trial);
            residual = next;
            iterations += 1;
        }
        self.finish(
            boundary,
            duration,
            elapsed,
            controls,
            work,
            iterations,
            sweeps,
            evaluations,
            &mut cancelled,
        )
    }

    fn finish(
        &self,
        boundary: Boundary<'_>,
        duration: f64,
        elapsed: f64,
        controls: CompetitiveMembraneSolve,
        work: &CompetitiveMembraneWorkspace,
        iterations: usize,
        sweeps: usize,
        evaluations: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CompetitiveMembraneCandidate> {
        let n = self.capacities.len();
        let rows = work.rows;
        let offset = boundary.offset();
        let mut channels = Vec::with_capacity(self.channels.len());
        let mut reports = Vec::with_capacity(self.channels.len());
        for (s, old) in self.channels.iter().enumerate() {
            poll(cancelled)?;
            let start = s * rows + offset;
            let mobile = MembraneDiffusionState {
                model: old.mobile.model.clone(),
                partial_pressures_pa: work.current[start..start + n].to_vec(),
                elapsed_s: elapsed,
            };
            let trapped = work.occupied[start..start + n].to_vec();
            let pressures = if offset == 0 {
                boundary.pressures()[s]
            } else {
                [work.current[s * rows], work.current[(s + 1) * rows - 1]]
            };
            let transport = mobile.sample_fluxes(pressures, &mut *cancelled)?;
            let area = old.mobile.model.area_m2;
            let boundary_transfers = [
                finite((transport.face_molar_fluxes_mol_m2_s[0] * area) * duration)?,
                finite((transport.face_molar_fluxes_mol_m2_s[n] * area) * duration)?,
            ];
            let mut mobile_change = Sum::default();
            let mut trapped_change = Sum::default();
            let mut before_trapped = Sum::default();
            let mut after_trapped = Sum::default();
            let mut maximum: f64 = 0.0;
            for i in 0..n {
                if i % 128 == 0 {
                    poll(cancelled)?;
                }
                let dm = finite(
                    old.mobile.model.cells[i].capacity_mol_m2_pa
                        * area
                        * (mobile.partial_pressures_pa[i] - old.mobile.partial_pressures_pa[i]),
                )?;
                let dn = finite(trapped[i] - old.trapped[i])?;
                let incoming = finite((transport.face_molar_fluxes_mol_m2_s[i] * area) * duration)?;
                let outgoing =
                    finite((transport.face_molar_fluxes_mol_m2_s[i + 1] * area) * duration)?;
                let residual = finite((dm + dn) - (incoming - outgoing))?;
                let tolerance = controls
                    .balance
                    .tolerance(dm.abs() + dn.abs() + incoming.abs() + outgoing.abs())?;
                check(residual, tolerance)?;
                check(work.site_residual[start + i], tolerance)?;
                maximum = maximum
                    .max(residual.abs())
                    .max(work.site_residual[start + i].abs());
                mobile_change.add(dm)?;
                trapped_change.add(dn)?;
                before_trapped.add(old.trapped[i])?;
                after_trapped.add(trapped[i])?;
            }
            let membrane_residual = finite(
                (mobile_change.value + trapped_change.value)
                    - (boundary_transfers[0] - boundary_transfers[1]),
            )?;
            check(
                membrane_residual,
                controls.balance.tolerance(
                    mobile_change.value.abs()
                        + trapped_change.value.abs()
                        + boundary_transfers[0].abs()
                        + boundary_transfers[1].abs(),
                )?,
            )?;
            let finite_reservoirs = if let Boundary::Finite {
                pressures: before,
                capacities,
            } = boundary
            {
                let changes = [
                    finite(capacities[0] * (pressures[0] - before[s][0]))?,
                    finite(capacities[1] * (pressures[1] - before[s][1]))?,
                ];
                let residuals = [
                    finite(changes[0] + boundary_transfers[0])?,
                    finite(changes[1] - boundary_transfers[1])?,
                ];
                for i in 0..2 {
                    check(
                        residuals[i],
                        controls
                            .balance
                            .tolerance(changes[i].abs() + boundary_transfers[i].abs())?,
                    )?;
                }
                let mut closed = Sum::default();
                for v in [
                    changes[0],
                    mobile_change.value,
                    trapped_change.value,
                    changes[1],
                ] {
                    closed.add(v)?;
                }
                check(
                    closed.value,
                    controls.balance.tolerance(
                        changes[0].abs()
                            + changes[1].abs()
                            + mobile_change.value.abs()
                            + trapped_change.value.abs(),
                    )?,
                )?;
                Some(TrappingReservoirBalance {
                    pressures_before_pa: before[s],
                    pressures_after_pa: pressures,
                    amount_changes_mol: changes,
                    amount_residuals_mol: residuals,
                    closed_amount_residual_mol: closed.value,
                })
            } else {
                None
            };
            reports.push((
                old.species.clone(),
                TrappingMembraneStep {
                    duration_s: duration,
                    transport,
                    mobile_before_mol: old.mobile.stored_amount_mol()?,
                    mobile_after_mol: mobile.stored_amount_mol()?,
                    mobile_change_mol: mobile_change.value,
                    trapped_before_mol: before_trapped.value,
                    trapped_after_mol: after_trapped.value,
                    trapped_change_mol: trapped_change.value,
                    gross_adsorption_mol: work.adsorbed[start..start + n].to_vec(),
                    gross_desorption_mol: work.desorbed[start..start + n].to_vec(),
                    boundary_transfers_mol: boundary_transfers,
                    membrane_amount_residual_mol: membrane_residual,
                    maximum_cell_amount_residual_mol: maximum,
                    finite_reservoirs,
                    nonlinear_iterations: iterations,
                    residual_evaluations: evaluations,
                },
            ));
            channels.push(Channel {
                species: old.species.clone(),
                mobile,
                traps: old.traps.clone(),
                trapped,
            });
        }
        let mut vacancy_before = Vec::with_capacity(n);
        for i in 0..n {
            poll(cancelled)?;
            let mut sum = Sum::default();
            for c in &self.channels {
                sum.add(c.trapped[i])?;
            }
            vacancy_before.push(finite(self.capacities[i] - sum.value)?);
            let mut occupied = Sum::default();
            for c in &channels {
                occupied.add(c.trapped[i])?;
            }
            if occupied.value > self.capacities[i] {
                return Err(invalid(
                    "rounded final species inventories exceed the shared site capacity",
                ));
            }
        }
        let state = Self {
            channels,
            capacities: self.capacities.clone(),
        };
        let step = CompetitiveMembraneStep {
            species: reports,
            vacancy_before_mol: vacancy_before,
            vacancy_after_mol: work.vacancy[offset..offset + n].to_vec(),
            nonlinear_iterations: iterations,
            linear_sweeps: sweeps,
            residual_evaluations: evaluations,
        };
        poll(cancelled)?;
        Ok(CompetitiveMembraneCandidate { state, step })
    }
}

impl CompetitiveMembraneWorkspace {
    fn correct(
        &mut self,
        maximum_sweeps: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<usize> {
        let rows = self.rows;
        let ns = self.species;
        for s in 0..ns {
            let start = s * rows;
            self.factors[s].factor(
                &self.lower[start..start + rows - 1],
                &self.diagonal[start..start + rows],
                &self.upper[start..start + rows - 1],
                &mut *cancelled,
            )?;
        }
        self.correction.fill(0.0);
        self.cross.fill(0.0);
        self.cross_compensation.fill(0.0);
        let mut residual_ratio = f64::INFINITY;
        for sweep in 1..=maximum_sweeps {
            for s in 0..ns {
                poll(cancelled)?;
                let start = s * rows;
                for row in 0..rows {
                    if row % 128 == 0 {
                        poll(cancelled)?;
                    }
                    let k = start + row;
                    let other =
                        finite(self.cross[row] - self.coupling_column[k] * self.correction[k])?;
                    self.linear_rhs[row] =
                        finite(self.rhs[k] + (self.coupling_row[k] / self.scale[k]) * other)?;
                }
                self.factors[s].solve_factored_into(
                    &self.linear_rhs,
                    &mut self.next_correction,
                    &mut *cancelled,
                )?;
                for row in 0..rows {
                    let k = start + row;
                    let update = finite(
                        self.coupling_column[k] * (self.next_correction[row] - self.correction[k]),
                    )?;
                    let mut sum = Sum {
                        value: self.cross[row],
                        correction: self.cross_compensation[row],
                    };
                    sum.add(update)?;
                    self.cross[row] = sum.value;
                    self.cross_compensation[row] = sum.correction;
                    self.correction[k] = self.next_correction[row];
                }
            }
            // Independently reconstruct the original block Jacobian residual.
            // Prefix/suffix sums exclude each diagonal without subtracting a
            // dominant self term, still using O(cells*species) work.
            residual_ratio = 0.0;
            for row in 0..rows {
                if row % 128 == 0 {
                    poll(cancelled)?;
                }
                let mut prefix = [0.0; MAX_COMPETITIVE_SPECIES + 1];
                let mut suffix = [0.0; MAX_COMPETITIVE_SPECIES + 1];
                let mut sum = Sum::default();
                for s in 0..ns {
                    let k = s * rows + row;
                    sum.add(finite(self.coupling_column[k] * self.correction[k])?)?;
                    prefix[s + 1] = sum.value;
                }
                self.cross[row] = sum.value;
                self.cross_compensation[row] = sum.correction;
                let mut sum = Sum::default();
                for s in (0..ns).rev() {
                    let k = s * rows + row;
                    sum.add(finite(self.coupling_column[k] * self.correction[k])?)?;
                    suffix[s] = sum.value;
                }
                for s in 0..ns {
                    let k = s * rows + row;
                    let mut applied = finite(self.diagonal[k] * self.correction[k])?;
                    if row > 0 {
                        applied = finite(applied + self.lower[k - 1] * self.correction[k - 1])?;
                    }
                    if row + 1 < rows {
                        applied = finite(applied + self.upper[k] * self.correction[k + 1])?;
                    }
                    applied = finite(
                        applied
                            - (self.coupling_row[k] / self.scale[k]) * (prefix[s] + suffix[s + 1]),
                    )?;
                    residual_ratio = residual_ratio.max(ratio(
                        finite(applied - self.rhs[k])?,
                        self.linear_allowance[k],
                    ));
                }
            }
            if residual_ratio <= 1.0 {
                poll(cancelled)?;
                return Ok(sweep);
            }
        }
        Err(CompetitiveMembraneError::LinearWork {
            sweeps: maximum_sweeps,
            residual_ratio,
        })
    }
}
