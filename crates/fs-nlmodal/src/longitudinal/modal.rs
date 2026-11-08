//! Complete fixed-position string storage with retained axial inertia.
//!
//! Spatial strain products remain the original exact cosine identities. This
//! module only supplies the physical mass-normalized coordinate pullback and
//! the complete continuous/discrete derivatives needed by generic pHS solvers.

use super::{LongitudinalString, MAX_LONGITUDINAL_MODES, MAX_TRANSVERSE_MODES, invalid};
use crate::{KcStringParams, NlModalError, prestressed_beam_omega};
use fs_phs::{
    ConstitutiveDiscreteGradient, ContinuousStorageTangent, Storage, linear_exact::QuadraticPart,
};

const MAX_SLOPES: usize = 2 * (MAX_TRANSVERSE_MODES + 1);

/// Conservative moderate-slope string in complete mass-normalized coordinates.
///
/// State is interleaved `[q,p]`: all transverse modes in the first plane,
/// then the optional second plane, then the explicitly retained axial indices.
/// Every `p` is the physical mass-normalized velocity. Transverse positions
/// and axial displacements vanish at both ends; when EI is nonzero, the sine
/// basis also has zero bending moment (pinned, not clamped supports).
///
/// The original strain energy retains axial/transverse and cross-polarization
/// coupling, including sum/difference and second-harmonic generation. Omitted
/// axial strain channels remain statically equilibrated as in the source;
/// they are not simulated waves. This is neither a scalar Duffing surrogate
/// nor a geometrically exact/shear-deformable rod or an automatic band cutoff.
///
/// Source energy, force, continuous tangent and constitutive discrete gradient
/// share one law. Warm evaluation is read-only and allocation-free, with bounded
/// stack scratch and O(planes * (transverse_modes + 1) * (axial_modes + 1)) work.
#[derive(Clone, Debug)]
pub struct LongitudinalModalStorage {
    source: LongitudinalString,
    transverse_slope_scales: Vec<f64>,
    omegas: Vec<f64>,
    linear_stiffness: Vec<f64>,
}

impl LongitudinalModalStorage {
    /// Prepare the physical fixed-position modal pullback.
    ///
    /// One or two nonnegative bending rigidities EI [N m²] select the transverse
    /// planes. They share the original tension, mass and axial rigidity. For
    /// transverse mode n, `omega²=(T*kappa²+EI*kappa⁴)/mu`, `kappa=n*pi/L`;
    /// the retained axial frequencies remain `kappa*sqrt(EA/mu)`.
    /// Both directions use `phi=sqrt(2/(mu*L))*sin(kappa*x)`, so modal masses
    /// are exactly one and slope/force pullbacks are reciprocal.
    ///
    /// The source's limits apply: at most 512 transverse modes per plane and
    /// 32 explicit, strictly increasing axial indices in 1..=1024. Zero
    /// transverse modes is allowed only with a nonempty axial selection.
    /// No spatial quadrature or temporal aliasing assumption enters preparation;
    /// the caller must resolve the retained frequencies and generated harmonics.
    ///
    /// # Errors
    /// Invalid source material, EI, mode selection, empty system, or
    /// unrepresentable physical frequency, stiffness or normalization.
    pub fn new(
        params: KcStringParams,
        transverse_modes: usize,
        bending_rigidities_nm2: &[f64],
        axial_modes: &[usize],
    ) -> Result<Self, NlModalError> {
        if !(1..=2).contains(&bending_rigidities_nm2.len())
            || bending_rigidities_nm2
                .iter()
                .any(|ei| !ei.is_finite() || *ei < 0.0)
            || (transverse_modes == 0 && axial_modes.is_empty())
        {
            return Err(invalid(
                "longitudinal modal planes, bending rigidity or empty system",
            ));
        }
        let source = LongitudinalString::new(
            params,
            transverse_modes,
            bending_rigidities_nm2.len(),
            axial_modes,
        )?;
        let transverse_slope_scales: Vec<f64> = (1..=transverse_modes)
            .map(|n| {
                source.axial_displacement_scale * (n as f64 * core::f64::consts::PI / params.length)
            })
            .collect();
        let mut omegas =
            Vec::with_capacity(transverse_modes * bending_rigidities_nm2.len() + axial_modes.len());
        for &ei in bending_rigidities_nm2 {
            omegas.extend((1..=transverse_modes).map(|n| {
                prestressed_beam_omega(params.length, params.tension, params.lin_density, ei, n)
            }));
        }
        let mut linear_stiffness: Vec<f64> = omegas.iter().map(|omega| omega * omega).collect();
        omegas.extend_from_slice(source.axial_omegas());
        // Use the original strain owner's quadratic coefficient, not a second
        // independently rounded square of its reported axial frequency.
        linear_stiffness.extend(
            source
                .axial_slope_scales()
                .iter()
                .map(|d| source.strain_energy_coefficient() * d * d),
        );
        if transverse_slope_scales
            .iter()
            .chain(&omegas)
            .chain(&linear_stiffness)
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(invalid(
                "longitudinal modal frequency, stiffness or slope scale",
            ));
        }
        Ok(Self {
            source,
            transverse_slope_scales,
            omegas,
            linear_stiffness,
        })
    }

    /// Unchanged immutable strain source and retained axial-mode selection.
    pub fn source(&self) -> &LongitudinalString {
        &self.source
    }

    /// Number of interleaved displacement/velocity pairs, including axial modes.
    pub fn mode_count(&self) -> usize {
        self.omegas.len()
    }

    /// Number of retained transverse sine modes in each plane.
    pub fn transverse_modes_per_plane(&self) -> usize {
        self.source.transverse_modes
    }

    /// One or two transverse directions, coupled through the same axial strain.
    pub fn plane_count(&self) -> usize {
        self.source.planes
    }

    /// Physical small-amplitude angular frequencies, in the complete state order.
    pub fn omegas(&self) -> &[f64] {
        &self.omegas
    }

    /// One plane's mass-normalized slope coefficients, in increasing sine order.
    /// The same row maps every plane's q to the source's cosine slopes; chord
    /// coefficients are exactly zero because endpoint positions are fixed.
    pub fn transverse_slope_scales(&self) -> &[f64] {
        &self.transverse_slope_scales
    }

    fn transverse_count(&self) -> usize {
        self.source.transverse_modes * self.source.planes
    }

    fn valid_state(&self, state: &[f64]) -> bool {
        state.len() == 2 * self.mode_count() && state.iter().all(|x| x.is_finite())
    }

    fn project(
        &self,
        state: &[f64],
        slopes: &mut [f64; MAX_SLOPES],
        axial: &mut [f64; MAX_LONGITUDINAL_MODES],
    ) -> bool {
        if !self.valid_state(state) {
            return false;
        }
        let n = self.source.transverse_modes;
        for plane in 0..self.source.planes {
            slopes[plane * (n + 1)] = 0.0;
            for (j, &scale) in self.transverse_slope_scales.iter().enumerate() {
                slopes[plane * (n + 1) + j + 1] = scale * state[2 * (plane * n + j)];
            }
        }
        for (j, value) in axial
            .iter_mut()
            .enumerate()
            .take(self.source.axial_modes.len())
        {
            *value = state[2 * (self.transverse_count() + j)];
        }
        slopes[..self.source.slope_count()]
            .iter()
            .all(|x| x.is_finite())
    }

    fn pullback(
        &self,
        slope_effort: &[f64],
        axial_effort: &[f64],
        linear_coordinate: impl Fn(usize) -> f64,
        out: &mut [f64],
    ) -> bool {
        if out.len() != 2 * self.mode_count() {
            return false;
        }
        let n = self.source.transverse_modes;
        for plane in 0..self.source.planes {
            for (j, &scale) in self.transverse_slope_scales.iter().enumerate() {
                let mode = plane * n + j;
                out[2 * mode] = self.linear_stiffness[mode] * linear_coordinate(2 * mode)
                    + scale * slope_effort[plane * (n + 1) + j + 1];
            }
        }
        for (j, &effort) in axial_effort.iter().enumerate() {
            out[2 * (self.transverse_count() + j)] = effort;
        }
        for mode in 0..self.mode_count() {
            out[2 * mode + 1] = linear_coordinate(2 * mode + 1);
        }
        out.iter().all(|x| x.is_finite())
    }

    fn strain(&self, slopes: &[f64], axial: &[f64]) -> Result<Strain, NlModalError> {
        let count = self.source.axial_modes.len();
        let mut geometry = [0.0; MAX_LONGITUDINAL_MODES];
        let mean = self
            .source
            .geometric_strain(slopes, &mut geometry[..count])?;
        let mut result = Strain {
            mean,
            harmonics: [0.0; MAX_LONGITUDINAL_MODES],
        };
        self.source.storage_from_strain(
            mean,
            &geometry[..count],
            axial,
            &mut result.harmonics[..count],
        )?;
        Ok(result)
    }

    fn strain_direction(
        &self,
        slopes: &[f64],
        direction_slopes: &[f64],
        direction_axial: &[f64],
    ) -> Result<Strain, NlModalError> {
        let mut result = Strain {
            mean: 0.0,
            harmonics: [0.0; MAX_LONGITUDINAL_MODES],
        };
        result.mean = self.source.strain_rate(
            slopes,
            direction_slopes,
            direction_axial,
            &mut result.harmonics[..self.source.axial_modes.len()],
        )?;
        Ok(result)
    }

    // d(J^T W e) = J^T W de + (dJ)^T W e. Since geometric strain is
    // quadratic, dJ is the SAME original pullback evaluated at slope direction.
    fn strain_force_direction(
        &self,
        slopes: &[f64],
        direction_slopes: &[f64],
        strain: &Strain,
        direction: &Strain,
        slope_out: &mut [f64],
        axial_out: &mut [f64],
    ) -> Result<(), NlModalError> {
        let c = self.source.strain_energy_coefficient();
        let count = self.source.axial_modes.len();
        let mut weights = [0.0; MAX_LONGITUDINAL_MODES];
        for (j, weight) in weights.iter_mut().enumerate().take(count) {
            *weight = c * direction.harmonics[j];
        }
        self.source.geometric_strain_gradient(
            slopes,
            2.0 * c * direction.mean,
            &weights[..count],
            slope_out,
        )?;
        let mut curvature = [0.0; MAX_SLOPES];
        for (j, weight) in weights.iter_mut().enumerate().take(count) {
            *weight = c * strain.harmonics[j];
        }
        self.source.geometric_strain_gradient(
            direction_slopes,
            2.0 * c * strain.mean,
            &weights[..count],
            &mut curvature[..self.source.slope_count()],
        )?;
        for (value, correction) in slope_out.iter_mut().zip(curvature) {
            *value += correction;
        }
        for (j, value) in axial_out.iter_mut().enumerate() {
            *value = c * self.source.axial_slope_scales[j] * direction.harmonics[j];
        }
        if slope_out
            .iter()
            .chain(axial_out.iter())
            .any(|x| !x.is_finite())
        {
            Err(invalid("longitudinal modal strain tangent overflows"))
        } else {
            Ok(())
        }
    }
}

struct Strain {
    mean: f64,
    harmonics: [f64; MAX_LONGITUDINAL_MODES],
}

impl Storage for LongitudinalModalStorage {
    fn shared_clone(&self) -> Option<Box<dyn Storage>> {
        Some(Box::new(self.clone()))
    }

    fn hamiltonian(&self, state: &[f64]) -> f64 {
        let mut slopes = [0.0; MAX_SLOPES];
        let mut axial = [0.0; MAX_LONGITUDINAL_MODES];
        if !self.project(state, &mut slopes, &mut axial) {
            return f64::NAN;
        }
        let Ok(mut energy) = self.source.energy(
            &slopes[..self.source.slope_count()],
            &axial[..self.source.axial_modes.len()],
        ) else {
            return f64::NAN;
        };
        for mode in 0..self.mode_count() {
            energy += 0.5 * state[2 * mode + 1] * state[2 * mode + 1];
            if mode < self.transverse_count() {
                energy += 0.5 * self.linear_stiffness[mode] * state[2 * mode] * state[2 * mode];
            }
        }
        if energy.is_finite() { energy } else { f64::NAN }
    }

    fn gradient(&self, state: &[f64], out: &mut [f64]) {
        let mut slopes = [0.0; MAX_SLOPES];
        let mut axial = [0.0; MAX_LONGITUDINAL_MODES];
        let mut slope_effort = [0.0; MAX_SLOPES];
        let mut axial_effort = [0.0; MAX_LONGITUDINAL_MODES];
        let (s, a) = (self.source.slope_count(), self.source.axial_modes.len());
        if !self.project(state, &mut slopes, &mut axial) {
            out.fill(f64::NAN);
            return;
        }
        let Ok(strain) = self.strain(&slopes[..s], &axial[..a]) else {
            out.fill(f64::NAN);
            return;
        };
        // The one-state owner computes each geometric strain coefficient once,
        // using its contiguous dot products instead of evaluating two identical
        // endpoints of the discrete-gradient rule in every force callback.
        let c = self.source.strain_energy_coefficient();
        let mut weights = [0.0; MAX_LONGITUDINAL_MODES];
        for j in 0..a {
            weights[j] = c * strain.harmonics[j];
            axial_effort[j] = c * self.source.axial_slope_scales[j] * strain.harmonics[j];
        }
        if self
            .source
            .geometric_strain_gradient(
                &slopes[..s],
                2.0 * c * strain.mean,
                &weights[..a],
                &mut slope_effort[..s],
            )
            .is_err()
            || !self.pullback(&slope_effort[..s], &axial_effort[..a], |i| state[i], out)
        {
            out.fill(f64::NAN);
        }
    }

    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        Some(self)
    }
    fn discrete_gradient_rule(&self) -> Option<&dyn ConstitutiveDiscreteGradient> {
        Some(self)
    }

    fn quadratic_part(&self) -> Option<QuadraticPart> {
        let dimension = 2 * self.mode_count();
        let mut matrix = vec![0.0; dimension * dimension];
        for (mode, &k) in self.linear_stiffness.iter().enumerate() {
            matrix[(2 * mode) * dimension + 2 * mode] = k;
            matrix[(2 * mode + 1) * dimension + 2 * mode + 1] = 1.0;
        }
        Some(QuadraticPart {
            matrix,
            exact: self.transverse_count() == 0,
        })
    }
}

impl ContinuousStorageTangent for LongitudinalModalStorage {
    fn apply(&self, state: &[f64], direction: &[f64], out: &mut [f64]) -> bool {
        let mut slopes = [0.0; MAX_SLOPES];
        let mut axial = [0.0; MAX_LONGITUDINAL_MODES];
        let mut ds = [0.0; MAX_SLOPES];
        let mut du = [0.0; MAX_LONGITUDINAL_MODES];
        let (s, a) = (self.source.slope_count(), self.source.axial_modes.len());
        if !self.project(state, &mut slopes, &mut axial)
            || !self.project(direction, &mut ds, &mut du)
        {
            return false;
        }
        let Ok(strain) = self.strain(&slopes[..s], &axial[..a]) else {
            return false;
        };
        let Ok(delta) = self.strain_direction(&slopes[..s], &ds[..s], &du[..a]) else {
            return false;
        };
        let mut slope_out = [0.0; MAX_SLOPES];
        let mut axial_out = [0.0; MAX_LONGITUDINAL_MODES];
        self.strain_force_direction(
            &slopes[..s],
            &ds[..s],
            &strain,
            &delta,
            &mut slope_out[..s],
            &mut axial_out[..a],
        )
        .is_ok()
            && self.pullback(&slope_out[..s], &axial_out[..a], |i| direction[i], out)
    }
}

impl ConstitutiveDiscreteGradient for LongitudinalModalStorage {
    fn evaluate(&self, a: &[f64], b: &[f64], midpoint: &[f64], out: &mut [f64]) {
        let mut sa = [0.0; MAX_SLOPES];
        let mut sb = [0.0; MAX_SLOPES];
        let mut ua = [0.0; MAX_LONGITUDINAL_MODES];
        let mut ub = [0.0; MAX_LONGITUDINAL_MODES];
        let mut slope_out = [0.0; MAX_SLOPES];
        let mut axial_out = [0.0; MAX_LONGITUDINAL_MODES];
        let (s, m) = (self.source.slope_count(), self.source.axial_modes.len());
        if !self.valid_state(midpoint)
            || !self.project(a, &mut sa, &mut ua)
            || !self.project(b, &mut sb, &mut ub)
            || self
                .source
                .discrete_gradient(
                    &sa[..s],
                    &sb[..s],
                    &ua[..m],
                    &ub[..m],
                    &mut slope_out[..s],
                    &mut axial_out[..m],
                )
                .is_err()
            || !self.pullback(
                &slope_out[..s],
                &axial_out[..m],
                |i| f64::midpoint(a[i], b[i]),
                out,
            )
        {
            out.fill(f64::NAN);
        }
    }

    fn end_jacobian_action(
        &self,
        a: &[f64],
        b: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> bool {
        let mut sa = [0.0; MAX_SLOPES];
        let mut sb = [0.0; MAX_SLOPES];
        let mut ds = [0.0; MAX_SLOPES];
        let mut ua = [0.0; MAX_LONGITUDINAL_MODES];
        let mut ub = [0.0; MAX_LONGITUDINAL_MODES];
        let mut du = [0.0; MAX_LONGITUDINAL_MODES];
        let (s, m) = (self.source.slope_count(), self.source.axial_modes.len());
        if !self.project(a, &mut sa, &mut ua)
            || !self.project(b, &mut sb, &mut ub)
            || !self.project(direction, &mut ds, &mut du)
        {
            return false;
        }
        let Ok(mut mean) = self.strain(&sa[..s], &ua[..m]) else {
            return false;
        };
        let Ok(end) = self.strain(&sb[..s], &ub[..m]) else {
            return false;
        };
        let Ok(mut delta) = self.strain_direction(&sb[..s], &ds[..s], &du[..m]) else {
            return false;
        };
        // The exact source rule is J_geom((a+b)/2)^T W (e(a)+e(b))/2.
        // Its b derivative contains half the b strain rate, plus half the
        // direction pullback at the averaged strain. No midpoint Hessian
        // replaces this finite-endpoint Jacobian.
        mean.mean = f64::midpoint(mean.mean, end.mean);
        delta.mean *= 0.5;
        for j in 0..m {
            mean.harmonics[j] = f64::midpoint(mean.harmonics[j], end.harmonics[j]);
            delta.harmonics[j] *= 0.5;
        }
        for j in 0..s {
            sa[j] = f64::midpoint(sa[j], sb[j]);
            ds[j] *= 0.5;
        }
        let mut slope_out = [0.0; MAX_SLOPES];
        let mut axial_out = [0.0; MAX_LONGITUDINAL_MODES];
        self.strain_force_direction(
            &sa[..s],
            &ds[..s],
            &mean,
            &delta,
            &mut slope_out[..s],
            &mut axial_out[..m],
        )
        .is_ok()
            && self.pullback(
                &slope_out[..s],
                &axial_out[..m],
                |i| 0.5 * direction[i],
                out,
            )
    }
}
