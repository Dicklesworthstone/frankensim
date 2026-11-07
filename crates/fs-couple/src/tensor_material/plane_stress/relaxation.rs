//! Exact coordinates for complete constrained material energy, not a fit.
use super::*;
use fs_math::dd::Dd;

/// Exact independent relaxation coordinates of one original plane-stress point.
/// The source still owns normal-strain elimination, its full 3D law and strain
/// limit. Every original history coordinate survives this invertible transform.
/// Its energy is `q^T K_eq q/2 + sum (c_j^T q-z_j)^2/2`, with
/// `z_j' = rate_j (c_j^T q-z_j)`. Rates can differ from the source branch rates:
/// enforcing zero normal traction couples the ORIGINAL complete histories.
#[derive(Clone, Debug)]
pub struct PlaneStressRelaxation {
    point: PlaneStressPoint,
    equilibrium: Vec<f64>,
    rows: Vec<Vec<f64>>,
    rates: Vec<f64>,
    original_from_memory: Vec<f64>,
    // A retained viscous normal strain also depends on displacement: n=y-L*q.
    // Empty for the original pure-Maxwell image, whose arithmetic is unchanged.
    original_from_displacement: Vec<f64>,
}

struct RelaxationChange {
    energy: Vec<f64>,
    rates: Vec<f64>,
    original_memory_scale: Vec<f64>,
    original_from_displacement: Vec<f64>,
}

impl PlaneStressPoint {
    /// Express this already-condensed energy in independent passive channels.
    /// This is an exact change of history coordinates, using the shared symmetric
    /// eigensolver. No branches are condensed separately, fitted, clipped or
    /// truncated. The complete energy reconstruction must pass at 1e-9 relative
    /// diagonal scale. Cancellation is polled before/after the bounded dense
    /// eigensolve; interruption within that low-level call is not promised.
    /// Kelvin viscosity remains in the complete body evolution; this pure
    /// Maxwell image cannot silently discard its normal-rate or dashpot terms.
    pub fn relaxation_image(
        self,
        max_scalar_terms: usize,
        gate: &CancelGate,
    ) -> Result<PlaneStressRelaxation> {
        poll(gate)?;
        if self.source.viscosity_matrix().iter().any(|v| *v != 0.0) {
            return Err(bad(
                "Kelvin plane-stress viscosity is not a pure Maxwell relaxation image",
            ));
        }
        self.relaxation_coordinates(max_scalar_terms, gate, None)
    }

    /// Express complete Kelvin-Maxwell material as passive relaxation channels
    /// PLUS its full, instantaneous Kelvin feedthrough. Nothing is fitted or
    /// discarded. The original pure-Maxwell entry and its refusal are unchanged.
    ///
    /// At most one retained viscous normal coordinate is supported, as in the
    /// five-strain shell. With L=D_nn^-1 D_nu, y=sqrt(V)*(n+L*q) removes the
    /// normal/velocity skew coupling. Its positive mobility is V*D_nn^-1; the
    /// existing symmetric relaxation transform then applies to [a,y]. The
    /// original factor-congruence owner supplies D_uu-D_un D_nn^-1 D_nu.
    /// Multiple coupled viscous normal histories require the complete body or
    /// descriptor and are explicitly refused by this coordinate image.
    pub fn viscoelastic_image(
        self,
        max_scalar_terms: usize,
        gate: &CancelGate,
    ) -> Result<PlaneStressRelaxation> {
        poll(gate)?;
        if self.source.viscosity_matrix().iter().all(|v| *v == 0.0) {
            return self.relaxation_image(max_scalar_terms, gate);
        }
        let Some(normal) = &self.normal_viscosity else {
            return self.relaxation_coordinates(max_scalar_terms, gate, None);
        };
        if normal.indices.len() != 1 {
            return Err(bad(
                "Kelvin-Maxwell image supports at most one retained viscous normal strain",
            ));
        }
        let n = self.coordinates;
        let r = self.memory_coordinates();
        let m = n + r;
        let last = m - 1;
        let extra = 32 * (m as u128).pow(2) + 8 * r as u128 * n as u128 + 128;
        if extra > max_scalar_terms as u128 {
            return Err(bad(
                "Kelvin-Maxwell coordinate change exceeds scalar-work budget",
            ));
        }
        let root = det::sqrt(self.source.volume_m3());
        let inverse_root = root.recip();
        let mobility = self.source.volume_m3() * normal.mobility[0];
        if !inverse_root.is_finite() || !mobility.is_finite() || mobility <= 0.0 {
            return Err(bad("normalized Kelvin normal mobility is unrepresentable"));
        }
        // Original [q,a,n] = [q,a,y/sqrt(V)-L*q]. This sparse congruence
        // uses the existing compensated scalar owner for its signed work sums.
        // Maxwell history-history entries are copied without new arithmetic.
        let mut energy = self.energy.clone();
        for i in 0..m {
            poll(gate)?;
            let left = if i < n { -normal.coupling[i] } else { 0.0 };
            for j in i..m {
                let right = if j < n { -normal.coupling[j] } else { 0.0 };
                let value = if i == last {
                    self.energy[last * m + last] / root / root
                } else if j == last {
                    (Dd::from_f64(self.energy[i * m + last])
                        + Dd::from_f64(left) * Dd::from_f64(self.energy[last * m + last]))
                    .to_f64()
                        / root
                } else if i < n || j < n {
                    (Dd::from_f64(self.energy[i * m + j])
                        + Dd::from_f64(left) * Dd::from_f64(self.energy[last * m + j])
                        + Dd::from_f64(right) * Dd::from_f64(self.energy[i * m + last])
                        + (Dd::from_f64(left) * Dd::from_f64(right))
                            * Dd::from_f64(self.energy[last * m + last]))
                    .to_f64()
                } else {
                    self.energy[i * m + j]
                };
                energy[i * m + j] = value;
                energy[j * m + i] = value;
            }
        }
        if !finite(&energy) {
            return Err(bad("Kelvin-Maxwell energy coordinate change overflow"));
        }
        let mut rates = self.rates.clone();
        rates.push(mobility);
        let mut original_memory_scale = vec![1.0; r];
        original_memory_scale[r - 1] = inverse_root;
        let mut original_from_displacement = vec![0.0; r * n];
        for i in 0..n {
            original_from_displacement[(r - 1) * n + i] = -normal.coupling[i];
        }
        self.relaxation_coordinates(
            max_scalar_terms - extra as usize,
            gate,
            Some(RelaxationChange {
                energy,
                rates,
                original_memory_scale,
                original_from_displacement,
            }),
        )
    }

    fn relaxation_coordinates(
        self,
        max_scalar_terms: usize,
        gate: &CancelGate,
        change: Option<RelaxationChange>,
    ) -> Result<PlaneStressRelaxation> {
        poll(gate)?;
        let n = self.coordinates;
        let r = self.memory_coordinates();
        let m = n + r;
        let energy = change.as_ref().map_or(&self.energy, |c| &c.energy);
        let internal_rates = change.as_ref().map_or(&self.rates, |c| &c.rates);
        let constrained = self.normal_map.len();
        let d = n + constrained;
        let work = 512 * (r as u128).pow(3)
            + 8 * (m as u128).pow(3)
            + (n as u128).pow(2) * (d as u128).pow(2);
        if work > max_scalar_terms as u128 {
            return Err(bad(
                "complete plane-stress relaxation transform exceeds scalar-work budget",
            ));
        }
        // The equilibrium normal map comes from C_eq alone. It is used only
        // for K_eq; all transient normal strains still use the complete owner.
        let mut t = vec![vec![0.0; d]; n];
        for i in 0..n {
            poll(gate)?;
            t[i][i] = 1.0;
            let all: Vec<_> = (0..constrained).collect();
            let rhs: Vec<_> = all
                .iter()
                .map(|&k| -self.source.equilibrium_stiffness()[(n + k) * d + i])
                .collect();
            let normal =
                normal_subset_solve(self.source.equilibrium_stiffness(), d, n, &all, &rhs)?;
            t[i][n..].copy_from_slice(&normal);
        }
        let mut equilibrium = vec![0.0; n * n];
        for i in 0..n {
            for j in i..n {
                let mut value = 0.0;
                for a in 0..d {
                    for b in 0..d {
                        value += t[i][a] * self.source.equilibrium_stiffness()[a * d + b] * t[j][b];
                    }
                }
                equilibrium[i * n + j] = value;
                equilibrium[j * n + i] = value;
            }
        }
        if !finite(&equilibrium) {
            return Err(bad("plane-stress equilibrium projection overflow"));
        }
        if r == 0 {
            return Ok(PlaneStressRelaxation {
                point: self,
                equilibrium,
                rows: Vec::new(),
                rates: Vec::new(),
                original_from_memory: Vec::new(),
                original_from_displacement: Vec::new(),
            });
        }
        // H=[q,a]^T[A B; B^T D][q,a]/2, a'=-R(B^T q+D a).
        // S=sqrt(R) D sqrt(R)=U Lambda U^T;
        // z=sqrt(Lambda) U^T R^-1/2 a,
        // c=-B sqrt(R) U Lambda^-1/2.
        // Therefore H=q^T K_eq q/2+|c^T q-z|^2/2 and z'=Lambda(c^T q-z).
        let scale = internal_rates.iter().copied().fold(0.0_f64, f64::max);
        let roots: Vec<_> = internal_rates
            .iter()
            .map(|v| det::sqrt(v / scale))
            .collect();
        if roots.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(bad(
                "plane-stress memory rates are unresolved under common scaling",
            ));
        }
        let matrix: Vec<_> = (0..r * r)
            .map(|k| roots[k / r] * energy[(n + k / r) * m + n + k % r] * roots[k % r])
            .collect();
        if !finite(&matrix) {
            return Err(bad("plane-stress rate matrix overflow"));
        }
        poll(gate)?;
        let (values, vectors) = fs_la::eigen::jacobi_eigh(&matrix, r);
        poll(gate)?;
        if values.iter().any(|v| !v.is_finite() || *v <= 0.0) || !finite(&vectors) {
            return Err(bad(
                "plane-stress memory has an unresolved positive relaxation rate",
            ));
        }
        for i in 0..r {
            for j in 0..r {
                let orthogonal: f64 = (0..r)
                    .map(|k| vectors[k * r + i] * vectors[k * r + j])
                    .sum();
                let rebuilt: f64 = (0..r)
                    .map(|k| vectors[i * r + k] * values[k] * vectors[j * r + k])
                    .sum();
                let diagonal = det::sqrt(matrix[i * r + i]) * det::sqrt(matrix[j * r + j]);
                if !diagonal.is_finite()
                    || diagonal <= 0.0
                    || (orthogonal - if i == j { 1.0 } else { 0.0 }).abs() > 1e-10
                    || (rebuilt - matrix[i * r + j]).abs() > 1e-10 * diagonal
                {
                    return Err(bad(
                        "plane-stress relaxation eigensolve failed reconstruction",
                    ));
                }
            }
        }
        let rates: Vec<_> = values.iter().map(|v| scale * v).collect();
        let mut original_from_memory = vec![0.0; r * r];
        let mut memory_from_original = vec![0.0; r * r];
        for i in 0..r {
            for j in 0..r {
                original_from_memory[i * r + j] =
                    roots[i] * vectors[i * r + j] / det::sqrt(values[j]);
                memory_from_original[j * r + i] =
                    det::sqrt(values[j]) * vectors[i * r + j] / roots[i];
            }
        }
        let rows: Vec<Vec<f64>> = (0..r)
            .map(|channel| {
                (0..n)
                    .map(|i| {
                        -(0..r)
                            .map(|j| energy[i * m + n + j] * original_from_memory[j * r + channel])
                            .sum::<f64>()
                    })
                    .collect()
            })
            .collect();
        if rates
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0 || !v.recip().is_finite())
            || !finite(&original_from_memory)
            || !finite(&memory_from_original)
            || rows.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(bad("plane-stress relaxation coordinates overflow"));
        }
        // Reconstruct the COMPLETE original Q, including normal-history cross
        // terms; the frequency/time comparison is not circular with this gate.
        for i in 0..m {
            poll(gate)?;
            for j in 0..m {
                let mut value = if i < n && j < n {
                    equilibrium[i * n + j]
                } else {
                    0.0
                };
                for channel in 0..r {
                    let a = if i < n {
                        rows[channel][i]
                    } else {
                        -memory_from_original[channel * r + i - n]
                    };
                    let b = if j < n {
                        rows[channel][j]
                    } else {
                        -memory_from_original[channel * r + j - n]
                    };
                    value += a * b;
                }
                let diagonal = det::sqrt(energy[i * m + i]) * det::sqrt(energy[j * m + j]);
                if !value.is_finite()
                    || !diagonal.is_finite()
                    || (value - energy[i * m + j]).abs() > 1e-9 * diagonal.max(f64::MIN_POSITIVE)
                {
                    return Err(bad(
                        "plane-stress channels do not reconstruct the complete original energy",
                    ));
                }
            }
        }
        let original_from_displacement = if let Some(change) = change {
            for i in 0..r {
                for j in 0..r {
                    let value = original_from_memory[i * r + j];
                    let scaled = value * change.original_memory_scale[i];
                    if !scaled.is_finite() || (value != 0.0 && scaled == 0.0) {
                        return Err(bad(
                            "original Kelvin-Maxwell history scale is unrepresentable",
                        ));
                    }
                    original_from_memory[i * r + j] = scaled;
                }
            }
            change.original_from_displacement
        } else {
            Vec::new()
        };
        poll(gate)?;
        Ok(PlaneStressRelaxation {
            point: self,
            equilibrium,
            rows,
            rates,
            original_from_memory,
            original_from_displacement,
        })
    }
}

impl PlaneStressRelaxation {
    /// Project one unit-volume five-strain shell image into physical displacement
    /// columns and quadrature volume, without another constitutive eigensolve.
    /// Prepare `self` with `new_shell` on the five identity engineering-strain
    /// columns at volume 1. Every rate and original-history transformation is
    /// reused exactly; only strain rows, physical-history scales and work are
    /// projected. Retained viscous normal strains keep their physical units.
    ///
    /// The underlying full-law point is prepared by the same `new_shell` owner,
    /// so source strain limits, physical stresses and normal reconstruction stay
    /// available in the caller's coordinates. The budget includes that complete
    /// preparation plus this congruence and a full original-energy comparison.
    pub fn project_shell(
        &self,
        columns: &[[f64; 5]],
        volume_m3: f64,
        budget: TensorProjectionBudget,
        gate: &CancelGate,
    ) -> Result<Self> {
        poll(gate)?;
        const UNIT_COLUMNS: [Voigt; 6] = [
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.5, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.5],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        ];
        if self.point.coordinates() != 5
            || self.point.normal_map.len() != 1
            || self.point.source.volume_m3() != 1.0
            || self.point.source.strain_columns() != UNIT_COLUMNS
        {
            return Err(bad(
                "shell image projection requires unit-volume identity five-strain coordinates",
            ));
        }
        let n = columns.len();
        let r = self.rates.len();
        let m = n
            .checked_add(r)
            .ok_or_else(|| bad("shell image projection extent"))?;
        let extra = 64 * (m as u128).pow(2) + 8 * r as u128 * n as u128 + 256;
        if extra > budget.max_scalar_terms as u128 {
            return Err(bad("shell image projection exceeds scalar-work budget"));
        }
        let point = PlaneStressPoint::new_shell(
            self.point.material(),
            columns,
            volume_m3,
            TensorProjectionBudget {
                max_scalar_terms: budget.max_scalar_terms - extra as usize,
                ..budget
            },
            gate,
        )?;
        let root = det::sqrt(volume_m3);
        let rows: Vec<Vec<f64>> = self
            .rows
            .iter()
            .map(|row| {
                columns
                    .iter()
                    .map(|column| root * (0..5).fold(0.0, |s, k| s + row[k] * column[k]))
                    .collect()
            })
            .collect();
        let mut equilibrium = vec![0.0; n * n];
        for i in 0..n {
            poll(gate)?;
            for j in i..n {
                let mut value = 0.0;
                for a in 0..5 {
                    for b in 0..5 {
                        value += (root * columns[i][a])
                            * self.equilibrium[5 * a + b]
                            * (root * columns[j][b]);
                    }
                }
                equilibrium[i * n + j] = value;
                equilibrium[j * n + i] = value;
            }
        }
        if !finite(&equilibrium) || rows.iter().flatten().any(|v| !v.is_finite()) {
            return Err(bad("shell image projection overflow"));
        }
        // Maxwell factor histories and image histories scale by sqrt(V), but
        // original normal strain does not. Check the full physical Q under
        // q_unit=S*q, a_unit=a_global/sqrt(V), n_unit=n_global.
        let original = 5 + r;
        let maxwell = self.point.rates.len();
        let physical_scale = |index: usize| if index < maxwell { 1.0 } else { root };
        for i in 0..m {
            poll(gate)?;
            for j in i..m {
                let expected = if j < n {
                    let mut value = 0.0;
                    for a in 0..5 {
                        for b in 0..5 {
                            value += (root * columns[i][a])
                                * self.point.energy[a * original + b]
                                * (root * columns[j][b]);
                        }
                    }
                    value
                } else if i < n {
                    (0..5).fold(0.0, |s, k| {
                        s + root * columns[i][k] * self.point.energy[k * original + 5 + j - n]
                    }) * physical_scale(j - n)
                } else {
                    physical_scale(i - n)
                        * self.point.energy[(5 + i - n) * original + 5 + j - n]
                        * physical_scale(j - n)
                };
                let scale = det::sqrt(point.energy[i * m + i]) * det::sqrt(point.energy[j * m + j]);
                if !expected.is_finite()
                    || !scale.is_finite()
                    || (expected - point.energy[i * m + j]).abs()
                        > 1e-9 * scale.max(f64::MIN_POSITIVE)
                {
                    return Err(bad(
                        "shell image does not reconstruct the original full-law point energy",
                    ));
                }
            }
        }
        let mut original_from_memory = self.original_from_memory.clone();
        for i in maxwell..r {
            for j in 0..r {
                let value = original_from_memory[i * r + j];
                let scaled = value / root;
                if !scaled.is_finite() || (value != 0.0 && scaled == 0.0) {
                    return Err(bad("projected physical normal history is unrepresentable"));
                }
                original_from_memory[i * r + j] = scaled;
            }
        }
        let original_from_displacement = if self.original_from_displacement.is_empty() {
            Vec::new()
        } else {
            (0..r * n)
                .map(|k| {
                    (0..5).fold(0.0, |s, j| {
                        s + self.original_from_displacement[(k / n) * 5 + j] * columns[k % n][j]
                    })
                })
                .collect()
        };
        if !finite(&original_from_displacement) {
            return Err(bad(
                "projected normal displacement history is unrepresentable",
            ));
        }
        poll(gate)?;
        Ok(Self {
            point,
            equilibrium,
            rows,
            rates: self.rates.clone(),
            original_from_memory,
            original_from_displacement,
        })
    }

    /// The original full-law point and its normal-traction constraint.
    pub fn point(&self) -> &PlaneStressPoint {
        &self.point
    }
    /// Equilibrium stiffness, with no material memory counted twice.
    pub fn equilibrium_stiffness(&self) -> &[f64] {
        &self.equilibrium
    }
    /// Complete instantaneous Kelvin matrix in the original point coordinates.
    /// Its power is v^T C v; every off-diagonal term remains. Physical displacement
    /// columns give [N s/m], mass-normalized modal columns give [1/s]. It is zero
    /// for a pure-Maxwell image. Use this PLUS the relaxation channels exactly
    /// once; the constrained normal Kelvin loss is already in those channels.
    pub fn feedthrough_viscosity(&self) -> &[f64] {
        &self.point.viscosity
    }
    /// Exact strain rows c_j; zero/dependent channels are retained too.
    pub fn coupling_rows(&self) -> &[Vec<f64>] {
        &self.rows
    }
    /// Exact constrained relaxation rates [1/s], one per retained coordinate.
    pub fn relaxation_rates_per_s(&self) -> &[f64] {
        &self.rates
    }
    /// Recover original history when it is independent of displacement.
    /// A retained Kelvin normal offset requires `original_memory_at` instead.
    pub fn original_memory(&self, memory: &[f64]) -> Result<Vec<f64>> {
        if self.original_from_displacement.iter().any(|v| *v != 0.0) {
            return Err(bad(
                "Kelvin normal history reconstruction requires the actual displacement",
            ));
        }
        self.memory_contribution(memory)
    }

    fn memory_contribution(&self, memory: &[f64]) -> Result<Vec<f64>> {
        let r = self.rates.len();
        if memory.len() != r || !finite(memory) {
            return Err(bad("complete finite plane-stress memory required"));
        }
        let original: Vec<_> = self
            .original_from_memory
            .chunks(r.max(1))
            .map(|row| row.iter().zip(memory).map(|(a, z)| a * z).sum())
            .collect();
        if !finite(&original) {
            return Err(bad("original plane-stress memory reconstruction overflow"));
        }
        Ok(original)
    }
    /// Recover original Maxwell factor histories followed by retained physical
    /// normal strain, including the exact n=y-L*q displacement contribution.
    pub fn original_memory_at(&self, displacement: &[f64], memory: &[f64]) -> Result<Vec<f64>> {
        let n = self.point.coordinates();
        if displacement.len() != n || !finite(displacement) {
            return Err(bad(
                "original material history requires matching finite displacement",
            ));
        }
        let mut original = self.memory_contribution(memory)?;
        if !self.original_from_displacement.is_empty() {
            for (i, value) in original.iter_mut().enumerate() {
                *value += (0..n).fold(0.0, |s, j| {
                    s + self.original_from_displacement[i * n + j] * displacement[j]
                });
            }
        }
        if !finite(&original) {
            return Err(bad("original displacement-dependent history overflow"));
        }
        Ok(original)
    }
    /// Reconstruct the full original strain and enforce its material domain.
    pub fn strain(&self, displacement: &[f64], memory: &[f64]) -> Result<Voigt> {
        self.point.strain(
            displacement,
            &self.original_memory_at(displacement, memory)?,
        )
    }
    /// Independently recompute normal traction from the original 3D law.
    pub fn normal_stress_pa(&self, displacement: &[f64], memory: &[f64]) -> Result<[f64; 3]> {
        self.point.normal_stress_pa(
            displacement,
            &self.original_memory_at(displacement, memory)?,
        )
    }

    /// Original physical normal/shear traction including the actual Kelvin
    /// rates. Shell sigma_yz and sigma_zx remain generally nonzero.
    pub fn normal_stress_pa_with_velocity(
        &self,
        displacement: &[f64],
        velocity: &[f64],
        memory: &[f64],
    ) -> Result<[f64; 3]> {
        self.point.normal_stress_pa_with_velocity(
            displacement,
            velocity,
            &self.original_memory_at(displacement, memory)?,
        )
    }

    /// Original thickness-normal traction. Shell transverse shear remains
    /// physical and is not included in the zero-traction constraint.
    pub fn thickness_normal_stress_pa(&self, displacement: &[f64], memory: &[f64]) -> Result<f64> {
        self.point.thickness_normal_stress_pa(
            displacement,
            &self.original_memory_at(displacement, memory)?,
        )
    }

    /// Thickness-normal traction from the same complete velocity-aware law.
    pub fn thickness_normal_stress_pa_with_velocity(
        &self,
        displacement: &[f64],
        velocity: &[f64],
        memory: &[f64],
    ) -> Result<f64> {
        Ok(self.normal_stress_pa_with_velocity(displacement, velocity, memory)?[0])
    }
}
