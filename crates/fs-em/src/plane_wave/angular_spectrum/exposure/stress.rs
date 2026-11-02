//! Finite-window exterior Maxwell stress from the original vector fields.
//!
//! In the admitted homogeneous, lossless exterior,
//! T_ij = epsilon <E_i E_j> + mu <H_i H_j>
//!        - delta_ij (epsilon <E²> + mu <H²>)/2.
//! Products are of REAL fields, before any cycle average. Off-diagonal terms
//! retain polarization and spatial interference. A surface stress resultant is
//! not by itself a transient mechanical force: the field-momentum change in
//! the enclosed region must also be accounted for. No material momentum,
//! electrostriction or magnetostriction law is inferred here.

use super::*;

/// Additional admission for a stress query. Existing exposure preparation
/// already bounds its twelve scalar products per point/frequency pair, line
/// reconstruction and FFT work. This cap includes that original workspace,
/// six compensated dyadic sums per point and the complete stress output.
/// Borrowed owners, earlier outputs and allocator overhead are excluded.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumExposureStressLimits {
    pub maximum_workspace_bytes: usize,
}

/// Symmetric, real time/ensemble-averaged exterior Maxwell stress in pascals.
/// This is an electromagnetic stress, excluding mechanical medium stresses.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumExposureStressPoint {
    pub maxwell_stress_pa: [[f64; 3]; 3],
}

impl AngularSpectrumExposureStressPoint {
    /// Stress traction T*n [Pa] on a declared unit normal in the original
    /// Cartesian frame. Reversing the normal reverses the traction.
    pub fn traction_pa(&self, normal: [f64; 3]) -> Result<[f64; 3], AngularSpectrumError> {
        admit_normal(normal)?;
        let value = self.maxwell_stress_pa.map(|row| {
            row.into_iter()
                .zip(normal)
                .fold(Dd::ZERO, |s, (t, n)| s + Dd::from_f64(t) * Dd::from_f64(n))
                .to_f64()
        });
        finite_vector(value)
    }
}

/// Local stress and its periodic-cell surface quadrature. Spatial products
/// can need a finer grid than the underlying field; this reuses the original
/// grid without oversampling. It is not a continuous aperture integral,
/// finite-edge surface model or spatial-convergence certificate.
#[derive(Debug)]
pub struct AngularSpectrumExposureStressField {
    grid: AngularSpectrumGrid,
    points: Vec<AngularSpectrumExposureStressPoint>,
    // xx, yy, zz, xy, xz, yz. Retain low components when two boundary
    // resultants nearly cancel; conversion belongs to the public output.
    integrated_stress: [Dd; 6],
}

impl AngularSpectrumExposureStressField {
    pub fn grid(&self) -> AngularSpectrumGrid {
        self.grid
    }

    pub fn points(&self) -> &[AngularSpectrumExposureStressPoint] {
        &self.points
    }

    /// Original-grid quadrature of T_ij dA over the transverse periodic cell [N].
    pub fn integrated_stress_n(&self) -> [[f64; 3]; 3] {
        symmetric(self.integrated_stress.map(Dd::to_f64))
    }

    /// Contract the original transverse-cell tensor quadrature with a constant
    /// unit normal [N]. This does not rotate or resample the physical plane and
    /// does not automatically give a net force on a body.
    pub fn integrated_traction_n(
        &self,
        normal: [f64; 3],
    ) -> Result<[f64; 3], AngularSpectrumError> {
        admit_normal(normal)?;
        let [xx, yy, zz, xy, xz, yz] = self.integrated_stress;
        let tensor = [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]];
        finite_vector(tensor.map(|row| {
            row.into_iter()
                .zip(normal)
                .fold(Dd::ZERO, |s, (t, n)| s + t * Dd::from_f64(n))
                .to_f64()
        }))
    }
}

impl AngularSpectrumExposure<'_, '_, '_> {
    /// Conservative logical numerical storage for one complete stress plane.
    pub fn stress_workspace_bytes(&self) -> Result<usize, AngularSpectrumError> {
        stress_workspace(self.anchor.modes.len(), self.workspace_bytes())
    }

    /// Real-field Maxwell stress, including all carrier, beat, polarization
    /// and standing-wave terms. The two exterior half-spaces must retain the
    /// same original lossless material, as for the cycle-mean force owner.
    ///
    /// O(L*N + F*N*log N + F²*N), with twelve scalar temporal contractions
    /// per point/frequency pair. Original frequency fields are reconstructed
    /// once; no new constitutive or film solve is performed.
    pub fn stress_field(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureStressLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureStressField, AngularSpectrumError> {
        stress_field(
            self.anchor,
            core::iter::once((self, false)),
            self.stress_workspace_bytes()?,
            ExposureSide::OneSided(side),
            distance_from_face_m,
            window,
            limits,
            &mut cancelled,
        )
    }

    /// Stress of an incoming, outgoing or total field at either physical
    /// face. A total field combines both propagation directions BEFORE any
    /// dyadic product; separate incoming/outgoing stresses cannot replace it.
    pub fn illumination_stress_field(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureStressLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureStressField, AngularSpectrumError> {
        stress_field(
            self.anchor,
            core::iter::once((self, false)),
            self.stress_workspace_bytes()?,
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            window,
            limits,
            &mut cancelled,
        )
    }

    /// Grid-quadrature sum of exterior stress tractions on the two periodic-cell
    /// planes [N]: Q(T_right*ez - T_left*ez). Distances are outward
    /// from the left and right faces respectively; lateral periodic faces
    /// cancel. Only one field plane is retained at a time. Each plane uses the
    /// original per-plane work admission; total work is twice one plane query.
    ///
    /// This is NOT a net finite-window mechanical force. Even an empty
    /// lossless region can have a nonzero short-window resultant as field
    /// momentum changes. For a uniform shutter in vacuum, mechanical force
    /// requires subtracting (P_field,end-P_field,start)/duration; a shaped
    /// shutter instead requires its weighted field-momentum rate. A
    /// material-filled region requires its own consistent material/field
    /// momentum accounting. No such term
    /// is silently set to zero, and no material-force convention is inferred.
    pub fn boundary_stress_resultant_n(
        &self,
        distances_from_faces_m: [f64; 2],
        window: HarmonicWindow,
        limits: AngularSpectrumExposureStressLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<[f64; 3], AngularSpectrumError> {
        boundary_resultant(
            |side, distance| {
                self.illumination_stress_field(side, distance, window, limits, &mut cancelled)
            },
            distances_from_faces_m,
        )
    }
}

impl AngularSpectrumExposureEnsemble<'_, '_, '_, '_> {
    pub fn stress_workspace_bytes(&self) -> Result<usize, AngularSpectrumError> {
        stress_workspace(self.anchor.modes.len(), self.workspace_bytes())
    }

    /// Stress of the authored mean plus mutually uncorrelated proper factors.
    /// Proper factors retain difference-frequency beats; only their explicitly
    /// zero pseudocovariance removes sum-frequency terms. The mean retains
    /// both. No stationary, Gaussian or independent-face source is inferred.
    /// One factor's original frequency fields are retained at a time.
    pub fn stress_field(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureStressLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureStressField, AngularSpectrumError> {
        stress_field(
            self.anchor,
            self.mean
                .into_iter()
                .map(|bank| (bank, false))
                .chain(self.factors.iter().map(|bank| (bank, true))),
            self.stress_workspace_bytes()?,
            ExposureSide::OneSided(side),
            distance_from_face_m,
            window,
            limits,
            &mut cancelled,
        )
    }

    pub fn illumination_stress_field(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureStressLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureStressField, AngularSpectrumError> {
        stress_field(
            self.anchor,
            self.mean
                .into_iter()
                .map(|bank| (bank, false))
                .chain(self.factors.iter().map(|bank| (bank, true))),
            self.stress_workspace_bytes()?,
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            window,
            limits,
            &mut cancelled,
        )
    }

    /// The total-field surface-stress resultant, with the same finite-window
    /// momentum-storage distinction as the deterministic owner's method.
    pub fn boundary_stress_resultant_n(
        &self,
        distances_from_faces_m: [f64; 2],
        window: HarmonicWindow,
        limits: AngularSpectrumExposureStressLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<[f64; 3], AngularSpectrumError> {
        boundary_resultant(
            |side, distance| {
                self.illumination_stress_field(side, distance, window, limits, &mut cancelled)
            },
            distances_from_faces_m,
        )
    }
}

fn stress_workspace(points: usize, original: usize) -> Result<usize, AngularSpectrumError> {
    points
        .checked_mul(
            core::mem::size_of::<[Dd; 6]>()
                + core::mem::size_of::<AngularSpectrumExposureStressPoint>(),
        )
        .and_then(|v| v.checked_add(original))
        .ok_or_else(|| bad("exposure stress workspace overflow"))
}

fn stress_field<'e, 'b, 's, 'm>(
    anchor: &PlanarAngularSpectrum<'s, 'm>,
    banks: impl IntoIterator<Item = (&'e AngularSpectrumExposure<'b, 's, 'm>, bool)>,
    workspace: usize,
    side: ExposureSide,
    distance: f64,
    window: HarmonicWindow,
    limits: AngularSpectrumExposureStressLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AngularSpectrumExposureStressField, AngularSpectrumError>
where
    'b: 'e,
    's: 'b,
    'm: 's,
{
    poll(cancelled)?;
    if workspace > limits.maximum_workspace_bytes {
        return Err(bad("exposure stress exceeds its workspace budget"));
    }
    if !distance.is_finite() || distance < 0.0 {
        return Err(bad("stress field distance must be finite and nonnegative"));
    }
    anchor.source().admit_momentum_exterior()?;
    let medium = anchor.source().incident_material();
    let epsilon = linear_electric_permittivity(medium).map_err(ReflectionError::from)?;
    let mu = linear_magnetic_permeability(medium).map_err(ReflectionError::from)?;
    let n = anchor.modes.len();
    let mut sums = Vec::new();
    sums.try_reserve_exact(n)
        .map_err(|_| bad("exposure stress accumulation allocation failed"))?;
    for i in 0..n {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        sums.push([Dd::ZERO; 6]);
    }
    for (bank, proper) in banks {
        poll(cancelled)?;
        if proper {
            accumulate_stress_bank(
                bank,
                side,
                distance,
                epsilon,
                mu,
                |a, b| window.proper_product_kernel(a, b),
                proper_product,
                &mut sums,
                cancelled,
            )?;
        } else {
            accumulate_stress_bank(
                bank,
                side,
                distance,
                epsilon,
                mu,
                |a, b| window.product_kernel(a, b),
                product,
                &mut sums,
                cancelled,
            )?;
        }
    }
    let mut points = Vec::new();
    points
        .try_reserve_exact(n)
        .map_err(|_| bad("exposure stress result allocation failed"))?;
    let mut integrated = [Dd::ZERO; 6];
    for (i, dyad) in sums.into_iter().enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        let stress = stress_from_dyad(dyad);
        let value = stress.map(Dd::to_f64);
        if value.iter().any(|v| !v.is_finite()) {
            return Err(bad("exposure Maxwell stress is unrepresentable"));
        }
        for (sum, value) in integrated.iter_mut().zip(stress) {
            *sum = *sum + value;
        }
        points.push(AngularSpectrumExposureStressPoint {
            maxwell_stress_pa: symmetric(value),
        });
    }
    let grid = anchor.grid();
    let area = Dd::from_f64(grid.spacing_m[0]) * Dd::from_f64(grid.spacing_m[1]);
    let integrated_stress = integrated.map(|s| s * area);
    if integrated_stress.iter().any(|s| !s.to_f64().is_finite()) {
        return Err(bad("integrated exposure stress is unrepresentable"));
    }
    poll(cancelled)?;
    Ok(AngularSpectrumExposureStressField {
        grid,
        points,
        integrated_stress,
    })
}

fn accumulate_stress_bank<P: Copy>(
    bank: &AngularSpectrumExposure<'_, '_, '_>,
    side: ExposureSide,
    distance: f64,
    epsilon: f64,
    mu: f64,
    kernel: impl Fn(f64, f64) -> Result<P, fs_fft::harmonic_window::WindowError>,
    product: impl Fn(P, C64, C64) -> Result<Dd, MomentError>,
    sums: &mut [[Dd; 6]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AngularSpectrumError> {
    let f = bank.frequency_count();
    if f == 0 {
        return Ok(());
    }
    let mut kernels = Vec::new();
    kernels
        .try_reserve_exact(f * f)
        .map_err(|_| bad("exposure stress kernel allocation failed"))?;
    for &a in &bank.frequencies {
        for &b in &bank.frequencies {
            poll(cancelled)?;
            kernels.push(
                kernel(a, b)
                    .map_err(|_| bad("exposure stress temporal kernel is unrepresentable"))?,
            );
        }
    }
    let mut fields = Vec::new();
    fields
        .try_reserve_exact(f)
        .map_err(|_| bad("exposure stress field allocation failed"))?;
    for group in 0..f {
        fields.push(bank.frequency_field(group, side, distance, cancelled)?);
    }
    for (point, sum) in sums.iter_mut().enumerate() {
        poll(cancelled)?;
        let value = integrate_dyad(
            f,
            &kernels,
            |i| fields[i][point],
            epsilon,
            mu,
            &product,
            cancelled,
        )
        .map_err(|error| match error {
            MomentError::Cancelled => AngularSpectrumError::Cancelled,
            MomentError::Product(_) => bad("exposure stress product is unrepresentable"),
            MomentError::Unrepresentable => bad("exposure stress dyad is unrepresentable"),
        })?;
        for (sum, value) in sum.iter_mut().zip(value) {
            *sum = *sum + value;
        }
    }
    Ok(())
}

// Six independent components; symmetry is imposed only AFTER the complete
// ordered frequency-pair contraction, never on an individual cross term.
fn integrate_dyad<P: Copy>(
    frequencies: usize,
    kernels: &[P],
    field_at: impl Fn(usize) -> [C64; 6],
    epsilon: f64,
    mu: f64,
    product: impl Fn(P, C64, C64) -> Result<Dd, MomentError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Dd; 6], MomentError> {
    let epsilon = Dd::from_f64(epsilon);
    let mu = Dd::from_f64(mu);
    let mut sum = [Dd::ZERO; 6];
    for i in 0..frequencies {
        let a = field_at(i);
        for j in 0..frequencies {
            if j % 64 == 0 && cancelled() {
                return Err(MomentError::Cancelled);
            }
            let b = field_at(j);
            let kernel = kernels[i * frequencies + j];
            for (c, (u, v)) in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
                .into_iter()
                .enumerate()
            {
                sum[c] = sum[c]
                    + epsilon * product(kernel, a[u], b[v])?
                    + mu * product(kernel, a[u + 3], b[v + 3])?;
            }
        }
    }
    if sum.iter().any(|x| !x.to_f64().is_finite()) {
        return Err(MomentError::Unrepresentable);
    }
    Ok(sum)
}

fn stress_from_dyad(mut dyad: [Dd; 6]) -> [Dd; 6] {
    let half_trace = (dyad[0] + dyad[1] + dyad[2]) * Dd::from_f64(0.5);
    for value in &mut dyad[..3] {
        *value = *value - half_trace;
    }
    dyad
}

fn symmetric<T: Copy>([xx, yy, zz, xy, xz, yz]: [T; 6]) -> [[T; 3]; 3] {
    [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]
}

fn admit_normal(normal: [f64; 3]) -> Result<(), AngularSpectrumError> {
    if normal.iter().any(|n| !n.is_finite()) {
        return Err(bad("stress traction requires a finite unit normal"));
    }
    let norm2 = normal
        .into_iter()
        .fold(Dd::ZERO, |s, n| s + Dd::from_f64(n) * Dd::from_f64(n))
        .to_f64();
    if !norm2.is_finite() || (norm2 - 1.0).abs() > 1e-12 {
        return Err(bad("stress traction requires a finite unit normal"));
    }
    Ok(())
}

fn finite_vector(value: [f64; 3]) -> Result<[f64; 3], AngularSpectrumError> {
    if value.iter().any(|v| !v.is_finite()) {
        Err(bad("exposure stress traction is unrepresentable"))
    } else {
        Ok(value)
    }
}

fn boundary_resultant(
    mut field: impl FnMut(
        AngularSpectrumIlluminationSide,
        f64,
    ) -> Result<AngularSpectrumExposureStressField, AngularSpectrumError>,
    distances: [f64; 2],
) -> Result<[f64; 3], AngularSpectrumError> {
    if distances.iter().any(|d| !d.is_finite() || *d < 0.0) {
        return Err(bad(
            "stress boundary distances must be finite and nonnegative",
        ));
    }
    use AngularSpectrumIlluminationSide::{TotalLeft, TotalRight};
    let left = field(TotalLeft, distances[0])?.integrated_stress;
    let right = field(TotalRight, distances[1])?.integrated_stress;
    finite_vector([4, 5, 2].map(|c| (right[c] - left[c]).to_f64()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn evaluate(
        fields: &[[C64; 6]],
        omega: &[f64],
        window: HarmonicWindow,
        proper: bool,
    ) -> [[f64; 3]; 3] {
        let dyad = if proper {
            let kernels = omega
                .iter()
                .flat_map(|&a| {
                    omega
                        .iter()
                        .map(move |&b| window.proper_product_kernel(a, b).unwrap())
                })
                .collect::<Vec<_>>();
            integrate_dyad(
                omega.len(),
                &kernels,
                |i| fields[i],
                1.3,
                0.7,
                proper_product,
                &mut || false,
            )
            .unwrap()
        } else {
            let kernels = omega
                .iter()
                .flat_map(|&a| {
                    omega
                        .iter()
                        .map(move |&b| window.product_kernel(a, b).unwrap())
                })
                .collect::<Vec<_>>();
            integrate_dyad(
                omega.len(),
                &kernels,
                |i| fields[i],
                1.3,
                0.7,
                product,
                &mut || false,
            )
            .unwrap()
        };
        symmetric(stress_from_dyad(dyad).map(Dd::to_f64))
    }

    #[test]
    fn g1_full_stress_matches_real_time_and_authored_proper_ensemble() {
        // General complex Cartesian fields exercise all six dyadic entries.
        // The production owner supplies Maxwell-consistent source fields;
        // this independent contraction fixture makes no propagation claim.
        let fields = [
            [
                C64::new(0.8, 0.2),
                C64::new(-0.2, 0.5),
                C64::new(0.3, -0.4),
                C64::new(-0.5, 0.1),
                C64::new(0.6, 0.2),
                C64::new(-0.1, 0.3),
            ],
            [
                C64::new(-0.4, 0.7),
                C64::new(0.2, 0.3),
                C64::new(-0.6, 0.1),
                C64::new(0.3, -0.2),
                C64::new(0.4, -0.1),
                C64::new(0.5, 0.6),
            ],
        ];
        let omega = [2.0, 2.3];
        let window = HarmonicWindow::new(0.17, 0.31).unwrap();
        for proper in [false, true] {
            let phases = if proper {
                [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)]
            } else {
                [C64::ONE; 4]
            };
            let mut expected = [[0.0; 3]; 3];
            let intervals = 2048;
            for phase in phases {
                for k in 0..=intervals {
                    let t = window.start_s() + window.duration_s() * k as f64 / intervals as f64;
                    let mut f = [0.0; 6];
                    for (line, w) in fields.iter().zip(omega) {
                        let (s, c) = (w * t).sin_cos();
                        for j in 0..6 {
                            let a = line[j] * phase;
                            f[j] += a.re * c - a.im * s;
                        }
                    }
                    let weight = if k == 0 || k == intervals {
                        1.0
                    } else if k % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    let trace = (0..3)
                        .map(|j| 1.3 * f[j] * f[j] + 0.7 * f[j + 3] * f[j + 3])
                        .sum::<f64>();
                    for i in 0..3 {
                        for j in 0..3 {
                            let stress = 1.3 * f[i] * f[j] + 0.7 * f[i + 3] * f[j + 3]
                                - if i == j { 0.5 * trace } else { 0.0 };
                            expected[i][j] += weight * stress;
                        }
                    }
                }
            }
            let actual = evaluate(&fields, &omega, window, proper);
            for i in 0..3 {
                for j in 0..3 {
                    let reference = expected[i][j] / (12.0 * intervals as f64);
                    assert!(
                        (actual[i][j] - reference).abs() < 2e-13,
                        "{} != {}",
                        actual[i][j],
                        reference
                    );
                    assert_eq!(actual[i][j], actual[j][i]);
                }
            }
            assert!(actual[0][1].abs() > 0.01);
        }
    }

    #[test]
    fn g1_subcycle_sine_stress_retains_carrier_phase_and_properness() {
        let fields = [[
            C64::new(0.0, 1.0),
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
            C64::new(0.0, 1.0),
            C64::ZERO,
        ]];
        let window = HarmonicWindow::new(0.0, 1e-10).unwrap();
        let real = evaluate(&fields, &[10.0], window, false);
        let proper = evaluate(&fields, &[10.0], window, true);
        assert!((real[2][2] / (-1e-18 / 3.0) - 1.0).abs() < 2e-12);
        assert!((proper[2][2] + 0.5).abs() < 2e-15);
        let kernels = [window.product_kernel(10.0, 10.0).unwrap()];
        assert!(matches!(
            integrate_dyad(1, &kernels, |_| fields[0], 1.3, 0.7, product, &mut || true),
            Err(MomentError::Cancelled)
        ));
    }

    #[test]
    fn g3_stress_changes_covariantly_under_cartesian_axis_permutation() {
        let fields = [[
            C64::new(0.8, 0.2),
            C64::new(-0.4, 0.1),
            C64::new(0.3, 0.6),
            C64::new(-0.1, 0.7),
            C64::new(0.5, 0.3),
            C64::new(0.2, -0.6),
        ]];
        let window = HarmonicWindow::new(-0.13, 0.71).unwrap();
        let permutation = [2, 0, 1];
        let permuted = [std::array::from_fn(|i| {
            fields[0][permutation[i % 3] + 3 * (i / 3)]
        })];
        for proper in [false, true] {
            let original = evaluate(&fields, &[2.7], window, proper);
            let rotated = evaluate(&permuted, &[2.7], window, proper);
            for i in 0..3 {
                for j in 0..3 {
                    assert!(
                        (rotated[i][j] - original[permutation[i]][permutation[j]]).abs() < 2e-15
                    );
                }
            }
        }
    }
}
