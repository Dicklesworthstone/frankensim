//! Relaxed DC equilibrium on the original biased piezoelectric operators.
//!
//! Maxwell stresses vanish at equilibrium; each admitted Debye pole adds its
//! actual static permittivity. The shared harmonic stamp evaluates those limits
//! at zero frequency. Independent field recovery supplies the physical pole
//! charges and storage without creating an unrestricted transient history.

use super::super::super::{
    PiezoelectricElectrode, Result, TetPiezoelectricControls, TetPiezoelectricError,
    TetPiezoelectricModel, TetPiezoelectricSolution, dot3, equilibrated_solve, invalid, poll,
    prepare_boundary, validate_controls,
};
use super::super::relaxation;
use fs_exec::Cx;
use fs_material::piezoelectric::tensor::relaxation::RelaxingTensorPiezoelectric;
use fs_solid::linear3::{TetLinearElasticProblem, TetMaterialField};

impl TetPiezoelectricModel<'_> {
    /// One original-model relaxed solve shared by DC harmonic bias and transient
    /// initialization. The caller selects mechanical constraints; this operation
    /// does not create or transfer any constitutive history.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::piezoelectric::spatial) fn solve_relaxed_with_state(
        &self,
        intrinsic: &[RelaxingTensorPiezoelectric],
        static_supports: &[(usize, f64)],
        electrodes: &[PiezoelectricElectrode<'_>],
        nodal_forces_n: &[[f64; 3]],
        nodal_charges_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<(TetPiezoelectricSolution, Vec<f64>, f64)> {
        let model = self;
        poll(cx)?;
        validate_controls(controls)?;
        let n = model.nodes_m.len();
        let nu = 3 * n;
        if nodal_forces_n.len() != n
            || nodal_charges_c.len() != n
            || nodal_forces_n
                .iter()
                .flatten()
                .chain(nodal_charges_c)
                .any(|value| !value.is_finite())
        {
            return Err(invalid("complete finite static loads required"));
        }
        let mut boundary = prepare_boundary(
            model.nodes_m,
            model.tetrahedra,
            static_supports,
            electrodes,
            nodal_forces_n,
            nodal_charges_c,
            cx,
        )?;
        let fixed_dofs = static_supports
            .iter()
            .map(|&(dof, _)| dof)
            .collect::<Vec<_>>();
        let field = if model.elastic_materials.len() == 1 {
            TetMaterialField::Uniform(&model.elastic_materials[0])
        } else {
            TetMaterialField::PerElement(&model.elastic_materials)
        };
        TetLinearElasticProblem {
            nodes_m: model.nodes_m,
            tetrahedra: model.tetrahedra,
            materials: field,
            fixed_dofs: &fixed_dofs,
            budget: model.budget,
        }
        .validate_static_constraints(cx)?;
        model.stamp_thermal_bias(&mut boundary, cx)?;

        let (matrix, rhs) = {
            // Native preparation retains every original source/frame and requires
            // its band to admit DC. Finite upper frequency bounds remain valid.
            let coefficients = relaxation::prepare(model, intrinsic, 0.0, cx)?;
            let maximum_entries = model
                .tetrahedra
                .len()
                .checked_mul(16)
                .ok_or_else(|| invalid("static polarization entry budget overflows"))?;
            let mut entries = Vec::new();
            relaxation::stamp(
                model,
                &coefficients,
                &mut |row, column, value, _| {
                    if !value.re.is_finite() || value.im != 0.0 {
                        return Err(invalid("DC relaxation coefficient is not finite and real"));
                    }
                    if value.re != 0.0 {
                        // Only the 4-by-4 dielectric block survives at DC. Mesh
                        // admission already bounds the number of original cells.
                        if row < nu || column < nu || entries.len() >= maximum_entries {
                            return Err(invalid(
                                "DC relaxation exceeds its dielectric entry budget",
                            ));
                        }
                        entries.push((row, column, value.re));
                    }
                    Ok(())
                },
                cx,
            )?;
            model.reduced_system(&boundary, 0.0, &entries, cx)?
        };
        let (solution, iterations, relative_residual) =
            equilibrated_solve(matrix, rhs, controls, cx)?;
        let mut full = boundary.prescribed.clone();
        for (i, index) in boundary.index.iter().enumerate() {
            if i % 1024 == 0 {
                poll(cx)?;
            }
            if let Some(index) = index {
                full[i] = solution[*index];
            }
        }

        let electric = model
            .dielectric
            .electric_field_v_m(&full[nu..], &mut || cx.checkpoint().is_err())?;
        let mut pole_charge_c = vec![0.0; n];
        let mut memories = Vec::with_capacity(model.tetrahedra.len());
        let mut memory_energy_j = 0.0;
        for (cell, tet) in model.tetrahedra.iter().enumerate() {
            poll(cx)?;
            let material = &intrinsic[if intrinsic.len() == 1 { 0 } else { cell }];
            let memory = material.equilibrium_response(electric[cell])?;
            let volume = model.dielectric.volumes_m3()[cell];
            let gradients = &model.dielectric.gradients()[cell];
            for local in 0..4 {
                pole_charge_c[tet[local]] -=
                    volume * dot3(gradients[local], memory.polarization_c_m2);
            }
            memory_energy_j += volume * memory.memory_energy_j_m3;
            memories.push((memory.polarization_c_m2, memory.memory_energy_j_m3));
        }
        let mut effective_charges_c = Vec::with_capacity(n);
        let mut displacement_m = Vec::with_capacity(n);
        let mut voltage_v = Vec::with_capacity(n);
        let mut polarization_work_j = 0.0;
        for node in 0..n {
            poll(cx)?;
            // Native recovery then returns Q_base-(q_user-Q_pole), exactly the
            // physical boundary reaction Q_total-q_user, including electrodes.
            let effective = nodal_charges_c[node] - pole_charge_c[node];
            if !pole_charge_c[node].is_finite() || !effective.is_finite() {
                return Err(invalid("static polarization charge is not representable"));
            }
            effective_charges_c.push(effective);
            displacement_m.push([full[3 * node], full[3 * node + 1], full[3 * node + 2]]);
            voltage_v.push(full[nu + node] + boundary.voltage_reference[node]);
            // Remove the work of the fictitious NEGATIVE pole volume load used
            // only for native base recovery. Relative potentials preserve gauges.
            polarization_work_j += (0.5 * full[nu + node]) * pole_charge_c[node];
        }
        if !memory_energy_j.is_finite() || !polarization_work_j.is_finite() {
            return Err(invalid(
                "static polarization storage or work is not representable",
            ));
        }
        let mut state = model.recover::<false>(
            &full,
            displacement_m,
            voltage_v,
            &boundary,
            electrodes,
            nodal_forces_n,
            &effective_charges_c,
            controls,
            iterations,
            relative_residual,
            cx,
        )?;
        for ((cell, (polarization, energy_density)), volume) in state
            .cells
            .iter_mut()
            .zip(memories)
            .zip(model.dielectric.volumes_m3())
        {
            poll(cx)?;
            for (value, extra) in cell.response.displacement_c_m2.iter_mut().zip(polarization) {
                *value += extra;
            }
            cell.response.stored_energy_j_m3 += energy_density;
            cell.stored_energy_j += energy_density * volume;
            if !cell.response.stored_energy_j_m3.is_finite()
                || !cell.stored_energy_j.is_finite()
                || cell
                    .response
                    .displacement_c_m2
                    .iter()
                    .any(|value| !value.is_finite())
            {
                return Err(invalid("relaxed static cell response is not representable"));
            }
        }
        state.stored_energy_j += memory_energy_j;
        state.electrical_ramp_work_j += polarization_work_j;
        // Preserve the independently evaluated base defect when a large inventory
        // would round it away. Pole work comes from actual nodal charges, whereas
        // pole storage comes from the material's directed equilibrium states.
        state.energy_residual_j += memory_energy_j - polarization_work_j;
        if [
            state.stored_energy_j,
            state.electrical_ramp_work_j,
            state.energy_residual_j,
        ]
        .iter()
        .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "relaxed static energy account is not representable",
            ));
        }
        if state.energy_residual_j.abs() > controls.energy_tolerance_j {
            return Err(TetPiezoelectricError::Balance(
                "independent relaxed static storage and ramp work do not close",
            ));
        }
        poll(cx)?;
        Ok((state, full, memory_energy_j))
    }
}
