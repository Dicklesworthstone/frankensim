//! Finite metal-sheet forming through the original three-dimensional return map.
//!
//! For `a = [F00,F01,F10,F11]`, solve `P[:,2] = 0` for the transverse
//! column in `F = [a0,a1,t0; a2,a3,t1; 0,0,t2]`. The thickness stretch
//! `t2` stays positive. All three tractions are relaxed: a retained plastic
//! metric can carry transverse shear even for an isotropic constitutive law.
//! Every local trial starts from the SAME accepted history. Only the returned
//! candidate is committed by the surrounding physical step.
//!
//! The tangent is the Schur complement of the original discrete `dP/dF`,
//! `Aaa - Aat Att^-1 Ata`. Stored energy, reduced incremental potential and
//! irreversible work retain their distinct native meanings. Thickness is a
//! solved kinematic variable, not an imposed total-volume constraint: plastic
//! flow is isochoric while elastic and thermal volume changes remain possible.
//!
//! This extends the existing plane-stress elimination to the multiplicative
//! return-map owner, following consistent algorithmic linearization as in
//! Simo and Taylor, doi:10.1002/nme.1620220310. It admits a locally stable
//! transverse branch, not sheet bending, wrinkling, necking, fracture, contact,
//! anisotropic yield or an empirical forming-limit prediction.

use super::combined::FiniteNortonJ2;
use super::creep::FiniteNortonCreep;
use super::rate::FiniteOverstressJ2;
use super::thermal::{FiniteJ2ThermalError, FrozenThermalFiniteJ2};
use super::{FiniteJ2, FiniteJ2Error, FiniteJ2Material, FiniteJ2Response, FiniteJ2State};
use crate::hyper::plane_stress::{PlaneStressSettings, PositiveBlock};

const IN_PLANE: [usize; 4] = [0, 1, 3, 4];
const TRANSVERSE: [usize; 3] = [2, 5, 8];

/// One converged transverse equilibrium at frozen accepted plastic history.
#[derive(Clone, Debug, PartialEq)]
pub struct FinitePlaneStressResponse {
    /// Full gradient, including the solved transverse shear and thickness.
    pub deformation_gradient: [f64; 9],
    /// Work-conjugate in-plane Piola components `[P00,P01,P10,P11]`, Pa.
    pub piola_pa: [f64; 4],
    /// Exact derivative of `piola_pa` with transverse equilibrium eliminated.
    pub tangent_pa: [[f64; 4]; 4],
    /// Original candidate history, full stresses/tangent and separate physical
    /// accounts, all per original reference volume. No history is reset.
    pub material: FiniteJ2Response,
    /// Transverse traction norm divided by a fixed compliance-based scale.
    pub relative_transverse_residual: f64,
    /// Accepted local Newton updates, excluding initialization and backtracks.
    pub iterations: usize,
}

impl FinitePlaneStressResponse {
    /// Positive current normal thickness per unit original thickness.
    pub fn thickness_stretch(&self) -> f64 {
        self.deformation_gradient[8]
    }

    /// Warm start for the next endpoint; this is geometry, not plastic history.
    pub fn transverse_column(&self) -> [f64; 3] {
        TRANSVERSE.map(|i| self.deformation_gradient[i])
    }
}

/// Complete fixed-temperature material behavior for retained plastic sheets.
///
/// Unlike an elastic/plastic backbone accessor, this interface retains the
/// selected flow law, physical interval, accepted-state admissibility and
/// independent stored-energy inspection. A viscous endpoint can legitimately
/// remain outside the equilibrium yield surface. Inspecting that endpoint
/// must neither project it nor advance its history through another interval.
///
/// The response potential is the integral used by the native return map; it
/// generally differs from stored energy plus physical inelastic loss. Spatial
/// consumers must retain both accounts. Evolving temperature additionally
/// requires the native two-temperature thermodynamic work account and is not
/// represented by this fixed-temperature interface.
pub trait FinitePlaneStressMaterial {
    fn initial_state(&self) -> FiniteJ2State;

    fn plane_stress_response(
        &self,
        in_plane: &[f64; 4],
        old: &FiniteJ2State,
        transverse_seed: [f64; 3],
        settings: PlaneStressSettings,
    ) -> Result<FinitePlaneStressResponse, FiniteJ2Error>;

    /// Fixed-history inspection, without constitutive evolution.
    fn stored_energy(&self, f: &[f64; 9], state: &FiniteJ2State) -> Result<f64, FiniteJ2Error>;

    /// Apply the actual selected law's accepted-state condition.
    fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2Error>;

    /// Physical time advanced by one accepted return; static laws have none.
    fn duration_s(&self) -> Option<f64> {
        None
    }
}

impl<L: FiniteJ2Material + ?Sized> FinitePlaneStressMaterial for L {
    fn initial_state(&self) -> FiniteJ2State {
        self.finite_j2().initial_state()
    }

    fn plane_stress_response(
        &self,
        in_plane: &[f64; 4],
        old: &FiniteJ2State,
        transverse_seed: [f64; 3],
        settings: PlaneStressSettings,
    ) -> Result<FinitePlaneStressResponse, FiniteJ2Error> {
        self.finite_j2()
            .plane_stress_response(in_plane, old, transverse_seed, settings)
    }

    fn stored_energy(&self, f: &[f64; 9], state: &FiniteJ2State) -> Result<f64, FiniteJ2Error> {
        self.finite_j2().stored_energy(f, state)
    }

    fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2Error> {
        self.finite_j2()
            .validate_committed(f, state, yield_tolerance_pa)
    }
}

impl FinitePlaneStressMaterial for FiniteOverstressJ2 {
    fn initial_state(&self) -> FiniteJ2State {
        Self::initial_state(self)
    }

    fn plane_stress_response(
        &self,
        in_plane: &[f64; 4],
        old: &FiniteJ2State,
        transverse_seed: [f64; 3],
        settings: PlaneStressSettings,
    ) -> Result<FinitePlaneStressResponse, FiniteJ2Error> {
        Self::plane_stress_response(self, in_plane, old, transverse_seed, settings)
    }

    fn stored_energy(&self, f: &[f64; 9], state: &FiniteJ2State) -> Result<f64, FiniteJ2Error> {
        Self::stored_energy(self, f, state)
    }

    fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2Error> {
        Self::validate_committed(self, f, state, yield_tolerance_pa)
    }

    fn duration_s(&self) -> Option<f64> {
        Some(Self::duration_s(self))
    }
}

impl FinitePlaneStressMaterial for FiniteNortonCreep {
    fn initial_state(&self) -> FiniteJ2State {
        Self::initial_state(self)
    }

    fn plane_stress_response(
        &self,
        in_plane: &[f64; 4],
        old: &FiniteJ2State,
        transverse_seed: [f64; 3],
        settings: PlaneStressSettings,
    ) -> Result<FinitePlaneStressResponse, FiniteJ2Error> {
        Self::plane_stress_response(self, in_plane, old, transverse_seed, settings)
    }

    fn stored_energy(&self, f: &[f64; 9], state: &FiniteJ2State) -> Result<f64, FiniteJ2Error> {
        Self::stored_energy(self, f, state)
    }

    fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2Error> {
        if !yield_tolerance_pa.is_finite() || yield_tolerance_pa < 0.0 {
            return Err(FiniteJ2Error::Parameters(
                "yield tolerance must be finite and nonnegative",
            ));
        }
        Self::validate_committed(self, f, state)
    }

    fn duration_s(&self) -> Option<f64> {
        Some(Self::duration_s(self))
    }
}

impl FinitePlaneStressMaterial for FiniteNortonJ2 {
    fn initial_state(&self) -> FiniteJ2State {
        Self::initial_state(self)
    }

    fn plane_stress_response(
        &self,
        in_plane: &[f64; 4],
        old: &FiniteJ2State,
        transverse_seed: [f64; 3],
        settings: PlaneStressSettings,
    ) -> Result<FinitePlaneStressResponse, FiniteJ2Error> {
        Self::plane_stress_response(self, in_plane, old, transverse_seed, settings)
    }

    fn stored_energy(&self, f: &[f64; 9], state: &FiniteJ2State) -> Result<f64, FiniteJ2Error> {
        Self::stored_energy(self, f, state)
    }

    fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2Error> {
        Self::validate_committed(self, f, state, yield_tolerance_pa)
    }

    fn duration_s(&self) -> Option<f64> {
        Some(Self::duration_s(self))
    }
}

fn numerical(what: &'static str) -> FiniteJ2Error {
    FiniteJ2Error::Numerical(what)
}

fn block(response: &FiniteJ2Response) -> Result<PositiveBlock, FiniteJ2Error> {
    PositiveBlock::new(TRANSVERSE.map(|i| TRANSVERSE.map(|j| response.tangent_pa[i][j]))).map_err(
        |_| {
            numerical(
                "plane-stress transverse tangent is not finite, symmetric and positive definite",
            )
        },
    )
}

fn residual(response: &FiniteJ2Response, scale: f64) -> f64 {
    TRANSVERSE.iter().fold(0.0_f64, |norm, &i| {
        norm.hypot(response.first_piola_pa[i] / scale)
    })
}

fn condense(
    f: [f64; 9],
    material: FiniteJ2Response,
    block: &PositiveBlock,
    relative_transverse_residual: f64,
    iterations: usize,
) -> Result<FinitePlaneStressResponse, FiniteJ2Error> {
    let a = &material.tangent_pa;
    let mut tangent = [[0.0; 4]; 4];
    for (j, &col) in IN_PLANE.iter().enumerate() {
        let extension = block
            .solve(TRANSVERSE.map(|k| a[k][col]))
            .map_err(|_| numerical("plane-stress tangent extension is not representable"))?;
        for (i, &row) in IN_PLANE.iter().enumerate() {
            tangent[i][j] = a[row][col]
                - (0..3)
                    .map(|k| a[row][TRANSVERSE[k]] * extension[k])
                    .sum::<f64>();
        }
    }
    if tangent.iter().flatten().any(|v| !v.is_finite()) {
        return Err(numerical(
            "condensed finite-plastic tangent is not representable",
        ));
    }
    Ok(FinitePlaneStressResponse {
        deformation_gradient: f,
        piola_pa: IN_PLANE.map(|i| material.first_piola_pa[i]),
        tangent_pa: tangent,
        material,
        relative_transverse_residual,
        iterations,
    })
}

// One reduction for all native finite laws. The callback captures the accepted
// history, endpoint temperature and actual physical interval immutably. Local
// Newton iterations therefore cannot advance creep or hardening more than once.
fn reduce<E: From<FiniteJ2Error>>(
    in_plane: &[f64; 4],
    transverse_seed: [f64; 3],
    natural_volume_ratio: f64,
    settings: PlaneStressSettings,
    evaluate: impl Fn(&[f64; 9]) -> Result<FiniteJ2Response, E>,
) -> Result<FinitePlaneStressResponse, E> {
    if !(settings.relative_tolerance.is_finite()
        && settings.relative_tolerance > 0.0
        && settings.relative_tolerance < 1.0)
        || settings.max_iterations == 0
        || settings.max_backtracks == 0
    {
        return Err(FiniteJ2Error::Parameters(
            "finite plane stress requires 0 < relative tolerance < 1 and positive iteration/backtracking budgets",
        ).into());
    }
    let determinant = in_plane[0] * in_plane[3] - in_plane[1] * in_plane[2];
    if in_plane
        .iter()
        .chain(&transverse_seed)
        .any(|v| !v.is_finite())
        || !determinant.is_finite()
        || determinant <= 0.0
        || transverse_seed[2] <= 0.0
        || !natural_volume_ratio.is_finite()
        || natural_volume_ratio <= 0.0
    {
        return Err(FiniteJ2Error::Parameters(
            "finite orientation-preserving in-plane gradient and positive finite thickness seed required",
        ).into());
    }
    let mut f = [
        in_plane[0],
        in_plane[1],
        transverse_seed[0],
        in_plane[2],
        in_plane[3],
        transverse_seed[1],
        0.0,
        0.0,
        transverse_seed[2],
    ];
    let initial = evaluate(&f).and_then(|value| {
        let block = block(&value).map_err(E::from)?;
        Ok((value, block))
    });
    let (mut value, mut positive) = match initial {
        Ok(candidate) => candidate,
        Err(initial_error) => {
            // A large tensile area change can make the old thickness unstable
            // or outside the elastic domain even when a stable root exists.
            // Try mechanical J=1, then bounded interpolants toward the supplied
            // geometry. This is only initialization; final J is never imposed.
            let target = natural_volume_ratio / determinant;
            if !target.is_finite() || target <= 0.0 {
                return Err(initial_error);
            }
            let mut seed = None;
            let mut fraction = 1.0;
            for _ in 0..settings.max_backtracks {
                let mut trial = f;
                trial[8] = (1.0 - fraction) * transverse_seed[2] + fraction * target;
                if let Ok(response) = evaluate(&trial)
                    && let Ok(stable) = block(&response)
                {
                    seed = Some((trial, response, stable));
                    break;
                }
                fraction *= 0.5;
            }
            let Some((gradient, response, stable)) = seed else {
                return Err(initial_error);
            };
            f = gradient;
            (response, stable)
        }
    };
    let gradient_norm = in_plane.iter().fold(0.0_f64, |n, v| n.hypot(*v));
    let scale = positive.traction_scale().map_err(|_| {
        E::from(numerical(
            "plane-stress traction scale is not representable",
        ))
    })? * gradient_norm.max(1.0);
    if !scale.is_finite() || scale <= 0.0 {
        return Err(numerical("plane-stress traction normalization is not representable").into());
    }
    for iteration in 0..=settings.max_iterations {
        let norm = residual(&value, scale);
        if norm <= settings.relative_tolerance {
            return condense(f, value, &positive, norm, iteration).map_err(E::from);
        }
        if iteration == settings.max_iterations || !norm.is_finite() {
            break;
        }
        let direction = positive
            .solve(TRANSVERSE.map(|i| -value.first_piola_pa[i]))
            .map_err(|_| {
                E::from(numerical(
                    "plane-stress Newton direction is not representable",
                ))
            })?;
        let mut step = 1.0;
        let mut accepted = None;
        for _ in 0..settings.max_backtracks {
            let mut trial = f;
            for i in 0..3 {
                trial[TRANSVERSE[i]] += step * direction[i];
            }
            if trial[8] > 0.0
                && let Ok(next) = evaluate(&trial)
                && residual(&next, scale) <= (1.0 - 1e-4 * step) * norm
                && let Ok(stable) = block(&next)
            {
                accepted = Some((trial, next, stable));
                break;
            }
            step *= 0.5;
        }
        let Some((gradient, response, stable)) = accepted else {
            return Err(numerical(
                "finite plane-stress line search exhausted its budget or material domain",
            )
            .into());
        };
        f = gradient;
        value = response;
        positive = stable;
    }
    Err(numerical("finite plane-stress Newton solve exhausted its iteration budget").into())
}

impl FiniteJ2 {
    /// Relax transverse shear and thickness at frozen accepted history using
    /// the original rate-independent return, storage and dissipation. The
    /// caller supplies a positive thickness seed (usually the previous column).
    /// Every local trial retains the original material bounds. No caller state
    /// changes on success or refusal; spatial callers poll between points.
    pub fn plane_stress_response(
        &self,
        in_plane: &[f64; 4],
        old: &FiniteJ2State,
        transverse_seed: [f64; 3],
        settings: PlaneStressSettings,
    ) -> Result<FinitePlaneStressResponse, FiniteJ2Error> {
        self.validate_history(old)?;
        reduce(in_plane, transverse_seed, 1.0, settings, |f| {
            self.response(f, old)
        })
    }
}

impl FiniteOverstressJ2 {
    /// Same traction-free reduction with the selected Perzyna duration/kinetics.
    /// Local Newton retries do not accumulate additional physical time or flow.
    pub fn plane_stress_response(
        &self,
        in_plane: &[f64; 4],
        old: &FiniteJ2State,
        transverse_seed: [f64; 3],
        settings: PlaneStressSettings,
    ) -> Result<FinitePlaneStressResponse, FiniteJ2Error> {
        self.validate_history(old)?;
        reduce(in_plane, transverse_seed, 1.0, settings, |f| {
            self.response(f, old)
        })
    }
}

impl FiniteNortonCreep {
    /// Threshold-free creep with its actual physical duration and native heat.
    pub fn plane_stress_response(
        &self,
        in_plane: &[f64; 4],
        old: &FiniteJ2State,
        transverse_seed: [f64; 3],
        settings: PlaneStressSettings,
    ) -> Result<FinitePlaneStressResponse, FiniteJ2Error> {
        self.validate_history(old)?;
        reduce(in_plane, transverse_seed, 1.0, settings, |f| {
            self.response(f, old)
        })
    }
}

impl FiniteNortonJ2 {
    /// Simultaneous creep/yielding retains both accumulations and both losses.
    pub fn plane_stress_response(
        &self,
        in_plane: &[f64; 4],
        old: &FiniteJ2State,
        transverse_seed: [f64; 3],
        settings: PlaneStressSettings,
    ) -> Result<FinitePlaneStressResponse, FiniteJ2Error> {
        self.validate_history(old)?;
        reduce(in_plane, transverse_seed, 1.0, settings, |f| {
            self.response(f, old)
        })
    }
}

impl FrozenThermalFiniteJ2 {
    /// Reduce the total-F response at the original frozen endpoint temperature.
    /// Thermal stretch and every selected creep/overstress mechanism survive.
    /// A kinetic inspection without a duration still returns `DurationRequired`.
    pub fn plane_stress_response(
        &self,
        in_plane: &[f64; 4],
        old: &FiniteJ2State,
        transverse_seed: [f64; 3],
        settings: PlaneStressSettings,
    ) -> Result<FinitePlaneStressResponse, FiniteJ2ThermalError> {
        self.validate_history(old)?;
        reduce(
            in_plane,
            transverse_seed,
            self.thermal_stretch().powi(3),
            settings,
            |f| self.response(f, old),
        )
    }
}
