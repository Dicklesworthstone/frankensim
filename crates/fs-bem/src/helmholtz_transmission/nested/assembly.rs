//! One coupled two-interface operator for real or molecular outgoing waves.
use super::*;

pub(super) fn solve_correction(
    drive: &[C64],
    factor: Option<&LuComplex>,
    incoming: &[C64],
    gate: &CancelGate,
) -> Result<Vec<C64>, HelmholtzError> {
    let m = incoming.len();
    let mut correction = zero_matrix(m, "core-shell trace correction")?;
    if let Some(factor) = factor {
        for (row, value) in correction.iter_mut().enumerate() {
            poll(gate)?;
            for (col, &input) in incoming.iter().enumerate() {
                *value = *value + drive[row * m + col] * input;
            }
        }
        if correction.iter().any(|&value| !finite(value)) {
            return Err(bad("nonfinite core-shell contrast drive"));
        }
        factor
            .solve_with_cancel(&mut correction, || gate.is_requested())
            .map_err(lu_error)?;
    }
    if correction.iter().any(|&value| !finite(value)) {
        return Err(bad("nonfinite core-shell scattered solution"));
    }
    Ok(correction)
}

pub(super) struct NestedAssembly<'a, K> {
    pub surfaces: [&'a SpherePanels; 2],
    pub density: [f64; 3],
    pub waves: [K; 3],
    pub beta: [f64; 2],
    pub offsets: [Option<usize>; 2],
    pub dimension: usize,
}

impl<K: WaveNumber> NestedAssembly<'_, K> {
    pub fn assemble(
        &self,
        gate: &CancelGate,
        difference: impl Fn(K, K, [f64; 3], [f64; 3], [f64; 3], [f64; 3]) -> [C64; 4],
        self_contrast: impl Fn(K, K, f64) -> (C64, C64, C64),
    ) -> Result<(Vec<C64>, Vec<C64>), HelmholtzError> {
        let m = self.dimension;
        let mut operator = zero_matrix(m * m, "core-shell coupled matrix")?;
        let mut drive = zero_matrix(m * m, "core-shell contrast forcing")?;
        for boundary in 0..2 {
            let Some(offset) = self.offsets[boundary] else {
                continue;
            };
            let outer_density = self.density[boundary];
            let inner_density = self.density[boundary + 1];
            let ratio = inner_density / outer_density;
            let contrast = (inner_density - outer_density) / outer_density;
            if !ratio.is_finite() || ratio <= 0.0 || !contrast.is_finite() {
                return Err(bad("unrepresentable core-shell density contrast"));
            }
            let local = assemble_wave(
                self.surfaces[boundary],
                self.waves[boundary],
                self.waves[boundary + 1],
                self.beta[boundary],
                ratio,
                contrast,
                gate,
                &difference,
                &self_contrast,
            )?;
            let size = 2 * self.surfaces[boundary].areas().len();
            for row in 0..size {
                poll(gate)?;
                for col in 0..size {
                    let target = (offset + row) * m + offset + col;
                    operator[target] = local[row * size + col];
                    drive[target] = -local[row * size + col];
                }
            }
        }
        if self.offsets.iter().all(Option::is_some) {
            self.assemble_shell_coupling(&mut operator, &mut drive, gate, difference)?;
        }
        if operator.iter().chain(&drive).any(|&value| !finite(value)) {
            return Err(bad("nonfinite core-shell coupled operator or forcing"));
        }
        Ok((operator, drive))
    }

    fn assemble_shell_coupling(
        &self,
        operator: &mut [C64],
        drive: &mut [C64],
        gate: &CancelGate,
        kernel_difference: impl Fn(K, K, [f64; 3], [f64; 3], [f64; 3], [f64; 3]) -> [C64; 4],
    ) -> Result<(), HelmholtzError> {
        let m = self.dimension;
        let rho0 = self.density[0];
        let rho1 = self.density[1];
        let density_change = (rho1 - rho0) / rho0;
        for target in 0..2 {
            let source = 1 - target;
            let target_offset = self.offsets[target].expect("two active interfaces");
            let source_offset = self.offsets[source].expect("two active interfaces");
            let a = self.surfaces[target];
            let b = self.surfaces[source];
            let (na, nb) = (a.areas().len(), b.areas().len());
            let sign = if target == 0 { -1.0 } else { 1.0 };
            let beta = self.beta[target];
            for i in 0..na {
                poll(gate)?;
                for j in 0..nb {
                    if j % 64 == 0 {
                        poll(gate)?;
                    }
                    let (x, nx) = (a.centroids()[i], a.normals()[i]);
                    let (y, ny) = (b.centroids()[j], b.normals()[j]);
                    let (s, d, dp, n) = kernels(self.waves[1], x, nx, y, ny);
                    let (s0, _, dp0, _) = kernels(self.waves[0], x, nx, y, ny);
                    let [ds, dd, ddp, dn] = kernel_difference(self.waves[0], self.waves[1], x, nx, y, ny);
                    let area = b.areas()[j];
                    let transfer = [
                        d,
                        s.scale(-self.beta[1]),
                        n.scale(1.0 / beta),
                        dp.scale(-self.beta[1] / beta),
                    ];
                    // E_shell-E_exterior, with the TARGET derivative row's
                    // scale unchanged. Original differences cancel static
                    // singular terms and retain representable weak contrasts.
                    let difference = [
                        dd,
                        -s0.scale(self.beta[0] * density_change) - ds.scale(self.beta[1]),
                        dn.scale(1.0 / beta),
                        -dp0.scale((self.beta[0] / beta) * density_change)
                            - ddp.scale(self.beta[1] / beta),
                    ];
                    let row = target_offset + i;
                    let col = source_offset + j;
                    let entries = [
                        row * m + col,
                        row * m + col + nb,
                        (row + na) * m + col,
                        (row + na) * m + col + nb,
                    ];
                    for entry in 0..4 {
                        operator[entries[entry]] = transfer[entry].scale(sign * area);
                        drive[entries[entry]] = difference[entry].scale(-sign * area);
                    }
                }
            }
        }
        // The reference wave obeys exterior momentum. At the inner target its
        // derivative Green identity therefore returns (rho0/rho1)*w, not w.
        // Factor the missing density contribution instead of subtracting one.
        let core_offset = self.offsets[1].expect("two active interfaces");
        let n = self.surfaces[1].areas().len();
        let density_correction = (rho1 - rho0) / rho1;
        for i in 0..n {
            let diagonal = (core_offset + n + i) * m + core_offset + n + i;
            drive[diagonal] = drive[diagonal] - C64::from_re(density_correction);
        }
        Ok(())
    }
}
