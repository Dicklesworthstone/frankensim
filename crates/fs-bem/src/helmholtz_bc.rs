//! e^-iwt Robin data on the existing CBIE/HBIE pair A p = B q.
//! Dimensionless b=rho*c*Y, passive flow v=v_drive-Yp, hence
//! q=i*omega*rho*v_drive-i*k*b*p and (A+i*k*B*b)p=B*q_drive.
//! Applying b to only one combined row is not this boundary condition.

use super::*;
use fs_exec::CancelGate;

/// Finite-source diffraction and source-bound total receiver fields.
#[path = "helmholtz_incident.rs"]
pub mod incident;

/// Hook for a qualified nonlocal BL operator, e.g. a tangential surface
/// Laplacian. It maps pressure [Pa] to rho*c*v_absorbed [Pa]. Its weighted
/// quadratic form must be passive unless explicitly declared active. The
/// implementation owns its discretization/units/budget qualification.
pub trait BoundaryOperator: Sync {
    /// Boundary degrees of freedom.
    fn len(&self) -> usize;
    /// Identity must bind the operator coefficients/discretization, for
    /// checkpoint cross-wiring guards; this is not a scientific certificate.
    fn fingerprint(&self) -> u64;
    /// Fully overwrite output with the normalized admittance image.
    fn apply(&self, pressure: &[C64], output: &mut [C64]) -> Result<(), HelmholtzError>;
}

/// Explicit normalized surface admittance, evaluated at this frequency.
#[derive(Clone, Copy)]
pub enum Admittance<'a> {
    /// Exact legacy Neumann path.
    Rigid,
    /// Dimensionless b=rho*c*Y per boundary DOF. Active data may be explored
    /// but cannot be used as passive acceptance evidence.
    Local {
        beta: &'a [C64],
        declared_active: bool,
    },
    /// No local-diagonal approximation of the supplied nonlocal operator.
    NonLocal {
        operator: &'a dyn BoundaryOperator,
        declared_active: bool,
    },
}
impl Admittance<'_> {
    pub(crate) fn validate(self, n: usize) -> Result<(), HelmholtzError> {
        match self {
            Self::Rigid => Ok(()),
            Self::Local {
                beta,
                declared_active,
            } => {
                if beta.len() != n {
                    return Err(HelmholtzError::ShapeMismatch {
                        what: "one normalized beta per boundary DOF",
                    });
                }
                if beta
                    .iter()
                    .any(|v| !finite(*v) || (!declared_active && v.re < 0.0))
                {
                    return Err(bad("nonfinite or undeclared active admittance"));
                }
                Ok(())
            }
            Self::NonLocal { operator, .. } if operator.len() == n => Ok(()),
            _ => Err(HelmholtzError::ShapeMismatch {
                what: "nonlocal admittance dimension",
            }),
        }
    }
    pub(crate) fn is_zero(self) -> bool {
        matches!(self, Self::Rigid)
            || matches!(self,Self::Local {beta,..} if beta.iter().all(|v|*v==C64::ZERO))
    }
    /// Active declarations are carried through results, never hidden.
    pub fn declared_active(self) -> bool {
        match self {
            Self::Rigid => false,
            Self::Local {
                declared_active, ..
            }
            | Self::NonLocal {
                declared_active, ..
            } => declared_active,
        }
    }
    pub(crate) fn image(self, p: &[C64]) -> Result<Vec<C64>, HelmholtzError> {
        self.validate(p.len())?;
        let mut out = vec![C64::ZERO; p.len()];
        match self {
            Self::Rigid => {}
            Self::Local { beta, .. } => {
                for ((o, &b), &v) in out.iter_mut().zip(beta).zip(p) {
                    *o = b * v;
                }
            }
            Self::NonLocal { operator, .. } => operator.apply(p, &mut out)?,
        }
        if out.iter().any(|v| !finite(*v)) {
            return Err(bad("nonfinite admittance image"));
        }
        Ok(out)
    }
    pub(crate) fn stamp(self) -> u64 {
        let mut hash = 0xcbf29ce484222325u64;
        let mut word = |v: u64| {
            hash = (hash ^ v).wrapping_mul(0x100000001b3);
        };
        match self {
            Self::Rigid => word(0),
            Self::Local {
                beta,
                declared_active,
            } => {
                word(1);
                word(u64::from(declared_active));
                for b in beta {
                    word(b.re.to_bits());
                    word(b.im.to_bits());
                }
            }
            Self::NonLocal {
                operator,
                declared_active,
            } => {
                word(2);
                word(operator.fingerprint());
                word(u64::from(declared_active));
            }
        }
        hash
    }
}

/// A uniform, rigid-backed porous coating resolved from one material and fluid.
///
/// The JCA(L) constitutive law remains in `fs-material`. This adapter retains
/// its source card, frequency support and coherent saturating-fluid snapshot,
/// and converts its exp(+i omega t) specific impedance to this solver's
/// exp(-i omega t) normalized admittance: `beta = rho*c / conj(Zs)`.
/// Specific impedance is pressure / particle velocity; no panel-area factor
/// belongs in this conversion.
///
/// This is a locally reacting boundary. Lateral propagation within the layer,
/// frame motion, curvature corrections to its planar backing, and exterior
/// viscothermal attenuation are not resolved.
#[derive(Debug, Clone)]
pub struct PorousRigidBoundary<'a> {
    law: &'a fs_material::porous::PorousFluidLaw,
    response: fs_material::porous::PorousFluidResponse,
    thickness_m: f64,
    medium: Medium,
    k: f64,
    beta: C64,
}

impl<'a> PorousRigidBoundary<'a> {
    /// Evaluate the porous law once at the actual frequency and fluid state.
    ///
    /// Source support, coherent fluid data and layer thickness are admitted by
    /// the original material owner. No air default or absorption-only fit is
    /// substituted for the complex impedance.
    pub fn new<M: fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid + ?Sized>(
        law: &'a fs_material::porous::PorousFluidLaw,
        fluid: &M,
        omega: f64,
        thickness_m: f64,
    ) -> Result<Self, fs_material::MaterialError> {
        let response = law.effective_fluid(fluid, omega)?;
        let impedance = response.surface_impedance_rigid(thickness_m)?;
        let acoustic = response.fluid().acoustic();
        let medium = Medium {
            density: acoustic.density_kg_m3(),
            sound_speed: acoustic.sound_speed_m_s(),
        };
        let k = omega / medium.sound_speed;
        let beta = C64::from_re(acoustic.characteristic_impedance_pa_s_m()) / impedance.conj();
        if !k.is_finite() || k <= 0.0 || !finite(beta) || beta.re < 0.0 {
            return Err(fs_material::MaterialError::State {
                what: "porous BEM boundary needs finite passive admittance and positive wavenumber"
                    .into(),
            });
        }
        Ok(Self { law, response, thickness_m, medium, k, beta })
    }

    /// Original skeleton characterization, provenance and frequency support.
    pub const fn law(&self) -> &fs_material::porous::PorousFluidLaw {
        self.law
    }

    /// The single coherent material response used for both wall and exterior.
    pub const fn response(&self) -> &fs_material::porous::PorousFluidResponse {
        &self.response
    }

    /// Thickness measured along the locally planar backing normal [m].
    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }

    /// Exterior acoustic coefficients copied from the same saturating fluid.
    pub const fn medium(&self) -> Medium {
        self.medium
    }

    /// Exterior wavenumber at the material query frequency [rad/m].
    pub const fn wavenumber_rad_m(&self) -> f64 {
        self.k
    }

    /// Dimensionless exp(-i omega t) coefficient for `Admittance::Local`.
    ///
    /// This also supports mixed rigid/coated patches: put this value on coated
    /// DOFs and zero on rigid DOFs, retaining this medium and query frequency.
    pub const fn normalized_admittance(&self) -> C64 {
        self.beta
    }

    /// Complex pressure reflection of a plane wave at this locally reacting
    /// surface, `R = (cos(theta) - beta) / (cos(theta) + beta)`.
    ///
    /// The cosine is positive for incidence toward the wall. This preserves
    /// phase; absorption alone does not determine a reflected pressure wave.
    /// It is a planar diagnostic, not a replacement for the BEM field solve.
    /// See COMSOL Acoustics, Wall / Normal impedance:
    /// <https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_pressure.05.124.html>.
    pub fn plane_wave_reflection(&self, incidence_cosine: f64) -> Result<C64, HelmholtzError> {
        if !incidence_cosine.is_finite()
            || incidence_cosine <= 0.0
            || incidence_cosine > 1.0
        {
            return Err(bad("porous reflection needs incidence cosine in (0, 1]"));
        }
        let cosine = C64::from_re(incidence_cosine);
        let reflected = (cosine - self.beta) / (cosine + self.beta);
        if !finite(reflected) {
            return Err(bad("nonfinite porous pressure reflection"));
        }
        Ok(reflected)
    }

    /// Uniform coating on the existing dense closed-surface radiation owner.
    pub fn solve_radiation(
        &self,
        surface: &SpherePanels,
        drive_velocity: &[C64],
        formulation: Formulation,
    ) -> Result<BoundarySolution, HelmholtzError> {
        // Reuse the original mesh/wavelength/cap admission before allocation.
        admit_resolution(surface, self.k)?;
        let beta = vec![self.beta; surface.areas().len()];
        solve_radiation_bc(
            surface,
            self.k,
            self.medium,
            SurfaceBc {
                drive_velocity,
                admittance: Admittance::Local { beta: &beta, declared_active: false },
            },
            formulation,
        )
    }

    /// Uniform coating on the existing dense plane-wave scattering owner.
    /// Absorbed work is evaluated from total incident plus scattered pressure.
    pub fn solve_scattering(
        &self,
        surface: &SpherePanels,
        incident: PlaneWave,
        formulation: Formulation,
    ) -> Result<ScatteringSolution, HelmholtzError> {
        admit_resolution(surface, self.k)?;
        let beta = vec![self.beta; surface.areas().len()];
        solve_scattering_bc(
            surface,
            self.k,
            self.medium,
            Admittance::Local { beta: &beta, declared_active: false },
            incident,
            formulation,
        )
    }

    /// The same material boundary through the original bounded FMM/GMRES owner.
    pub fn solve_scattering_iterative(
        &self,
        surface: &SpherePanels,
        incident: PlaneWave,
        formulation: Formulation,
        options: &IterativeOptions,
        gate: &CancelGate,
    ) -> Result<
        (
            ScatteringSolution,
            fs_solver::krylov::SolveReport,
            ResourceReceipt,
        ),
        IterativeError,
    > {
        if gate.is_requested() {
            return Err(IterativeError::Cancelled { matvecs: 0 });
        }
        let beta = vec![self.beta; surface.areas().len()];
        solve_scattering_bc_iterative(
            surface,
            self.k,
            self.medium,
            Admittance::Local { beta: &beta, declared_active: false },
            incident,
            formulation,
            options,
            gate,
        )
    }
}

/// Drive velocity [m/s] and the frequency/space-dependent boundary law.
#[derive(Clone, Copy)]
pub struct SurfaceBc<'a> {
    pub drive_velocity: &'a [C64],
    pub admittance: Admittance<'a>,
}

/// Outgoing trace plus separate wall absorption and drive work, in watts.
pub struct BoundarySolution {
    pub radiation: RadiationSolution,
    pub absorbed_power_w: f64,
    pub driving_power_w: f64,
    pub declared_active: bool,
}

/// Explicit scheduling of the same closed Robin equation across the dense
/// ceiling. Open P1 screens are a different discretization and do not use it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClosedRoute {
    Dense,
    Iterative,
}
pub fn closed_route(panels: usize) -> ClosedRoute {
    if panels <= MAX_DENSE_PANELS {
        ClosedRoute::Dense
    } else {
        ClosedRoute::Iterative
    }
}
pub struct BoundaryRun {
    pub solution: BoundarySolution,
    pub iterative: Option<(fs_solver::krylov::SolveReport, ResourceReceipt)>,
}
/// Route before attempting dense allocation. Iterative options are explicit;
/// the .109 owner retains its own geometry, ppw, FMM and work refusals.
pub fn solve_radiation_bc_auto(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    bc: SurfaceBc<'_>,
    formulation: Formulation,
    options: &IterativeOptions,
    gate: &CancelGate,
) -> Result<BoundaryRun, IterativeError> {
    if gate.is_requested() {
        return Err(IterativeError::Cancelled { matvecs: 0 });
    }
    match closed_route(surface.areas().len()) {
        ClosedRoute::Dense => Ok(BoundaryRun {
            solution: solve_radiation_bc(surface, k, medium, bc, formulation)?,
            iterative: None,
        }),
        ClosedRoute::Iterative => {
            let (solution, report, receipt) =
                solve_radiation_bc_iterative(surface, k, medium, bc, formulation, options, gate)?;
            Ok(BoundaryRun {
                solution,
                iterative: Some((report, receipt)),
            })
        }
    }
}

/// Incident unit-direction plane wave p=A exp(ik d.x), in pascals.
#[derive(Clone, Copy)]
pub struct PlaneWave {
    pub direction: [f64; 3],
    pub amplitude_pa: C64,
}
/// Outgoing scattered trace and physical absorption of the TOTAL wall trace.
/// Incident power is not mislabeled as an outgoing radiation solution.
pub struct ScatteringSolution {
    pub scattered: RadiationSolution,
    pub incident_pressure_pa: Vec<C64>,
    pub absorbed_power_w: f64,
    pub declared_active: bool,
}
fn scattering_drive(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    law: Admittance<'_>,
    incident: PlaneWave,
) -> Result<(Vec<C64>, Vec<C64>), HelmholtzError> {
    let d = incident.direction;
    if d.iter().any(|v| !v.is_finite())
        || (d.iter().map(|v| v * v).sum::<f64>() - 1.0).abs() > 1e-10
        || !finite(incident.amplitude_pa)
        || !k.is_finite()
        || k <= 0.0
        || !medium.density.is_finite()
        || medium.density <= 0.0
        || !medium.sound_speed.is_finite()
        || medium.sound_speed <= 0.0
    {
        return Err(bad(
            "scattering requires finite plane wave, unit direction and positive k/medium",
        ));
    }
    let p: Vec<_> = surface
        .centroids()
        .iter()
        .map(|x| incident.amplitude_pa * expik(k * (0..3).map(|a| d[a] * x[a]).sum::<f64>()))
        .collect();
    let y = law.image(&p)?;
    let rho_c = medium.density * medium.sound_speed;
    let drive = surface
        .normals()
        .iter()
        .zip(&p)
        .zip(y)
        .map(|((normal, &p), y)| {
            (p.scale(-(0..3).map(|a| d[a] * normal[a]).sum::<f64>()) - y).scale(1.0 / rho_c)
        })
        .collect();
    Ok((p, drive))
}
fn finish_scattering(
    surface: &SpherePanels,
    solution: BoundarySolution,
    incident: Vec<C64>,
    law: Admittance<'_>,
) -> Result<ScatteringSolution, HelmholtzError> {
    let total: Vec<_> = solution
        .radiation
        .pressure
        .iter()
        .zip(&incident)
        .map(|(&a, &b)| a + b)
        .collect();
    let y = law.image(&total)?;
    let rho_c = solution.radiation.medium.density * solution.radiation.medium.sound_speed;
    let absorbed = total
        .iter()
        .zip(y)
        .zip(surface.areas())
        .map(|((&p, y), &area)| (p.conj() * y).re * area / (2.0 * rho_c))
        .sum::<f64>();
    if !absorbed.is_finite() || (!law.declared_active() && absorbed < -1e-10) {
        return Err(bad("scattering wall failed passive-work check"));
    }
    Ok(ScatteringSolution {
        scattered: solution.radiation,
        incident_pressure_pa: incident,
        absorbed_power_w: absorbed,
        declared_active: law.declared_active(),
    })
}
/// Scattering Robin data enter BOTH original BM rows through the same solve.
pub fn solve_scattering_bc(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    law: Admittance<'_>,
    incident: PlaneWave,
    formulation: Formulation,
) -> Result<ScatteringSolution, HelmholtzError> {
    let (p, drive) = scattering_drive(surface, k, medium, law, incident)?;
    let solution = solve_radiation_bc(
        surface,
        k,
        medium,
        SurfaceBc {
            drive_velocity: &drive,
            admittance: law,
        },
        formulation,
    )?;
    finish_scattering(surface, solution, p, law)
}
/// The same plane-wave problem through .109's bounded iterative owner.
pub fn solve_scattering_bc_iterative(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    law: Admittance<'_>,
    incident: PlaneWave,
    formulation: Formulation,
    options: &IterativeOptions,
    gate: &CancelGate,
) -> Result<
    (
        ScatteringSolution,
        fs_solver::krylov::SolveReport,
        ResourceReceipt,
    ),
    IterativeError,
> {
    let (p, drive) = scattering_drive(surface, k, medium, law, incident)?;
    let (solution, report, receipt) = solve_radiation_bc_iterative(
        surface,
        k,
        medium,
        SurfaceBc {
            drive_velocity: &drive,
            admittance: law,
        },
        formulation,
        options,
        gate,
    )?;
    Ok((
        finish_scattering(surface, solution, p, law)?,
        report,
        receipt,
    ))
}
fn finite(v: C64) -> bool {
    v.re.is_finite() && v.im.is_finite()
}
fn bad(what: &'static str) -> HelmholtzError {
    HelmholtzError::BadParameter { what }
}

pub(crate) fn validate_closed(surface: &SpherePanels) -> Result<(), HelmholtzError> {
    let triangles = surface
        .triangles()
        .ok_or_else(|| bad("Robin surfaces require retained closed triangles"))?;
    surface
        .validate()
        .map_err(|_| bad("invalid Robin closed geometry"))?;
    let _components = crate::near_field::components(triangles)?;
    Ok(())
}

pub(crate) fn finish_boundary(
    surface: &SpherePanels,
    solution: RadiationSolution,
    bc: SurfaceBc<'_>,
) -> Result<BoundarySolution, HelmholtzError> {
    let absorbed = bc.admittance.image(&solution.pressure)?;
    let scale = 0.5 / (solution.medium.density * solution.medium.sound_speed);
    let mut power = 0.0;
    let mut apparent = 0.0;
    let mut driving = 0.0;
    for i in 0..absorbed.len() {
        power += (solution.pressure[i] * absorbed[i].conj()).re * surface.areas()[i] * scale;
        apparent += solution.pressure[i].abs() * absorbed[i].abs() * surface.areas()[i] * scale;
        driving +=
            0.5 * (solution.pressure[i] * bc.drive_velocity[i].conj()).re * surface.areas()[i];
    }
    if !power.is_finite()
        || !driving.is_finite()
        || (!bc.admittance.declared_active() && power < -1e-10 * apparent.max(f64::MIN_POSITIVE))
    {
        return Err(bad(
            "nonlocal boundary failed its discrete passive-work check",
        ));
    }
    Ok(BoundarySolution {
        radiation: solution,
        absorbed_power_w: power,
        driving_power_w: driving,
        declared_active: bc.admittance.declared_active(),
    })
}

fn admit_bc(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    bc: SurfaceBc<'_>,
) -> Result<f64, HelmholtzError> {
    let n = surface.areas().len();
    bc.admittance.validate(n)?;
    if !k.is_finite()
        || k <= 0.0
        || !medium.density.is_finite()
        || medium.density <= 0.0
        || !medium.sound_speed.is_finite()
        || medium.sound_speed <= 0.0
        || bc.drive_velocity.len() != n
        || bc.drive_velocity.iter().any(|v| !finite(*v))
    {
        return Err(bad(
            "Robin boundary requires positive finite k/medium and one finite drive per panel",
        ));
    }
    // Dense cap/resolution admission precedes the quadratic assembly.
    admit_resolution(surface, k)
}

/// Dense Robin solve from the existing assembled CBIE+HBIE pair. beta=0 calls
/// solve_radiation directly, preserving every historical pressure/velocity bit.
pub fn solve_radiation_bc(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    bc: SurfaceBc<'_>,
    formulation: Formulation,
) -> Result<BoundarySolution, HelmholtzError> {
    Ok(solve_radiation_bc_batch(
        surface,
        k,
        medium,
        &[bc.drive_velocity],
        bc.admittance,
        formulation,
    )?
    .remove(0))
}

/// Several prescribed drives through one Robin assembly and factorization.
/// Each result retains its own physical velocity, wall work and driving work.
/// The zero-admittance route uses the original rigid batch without changing
/// its arithmetic. The existing dense panel, resolution and batch caps apply.
pub fn solve_radiation_bc_batch(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    drive_velocities: &[&[C64]],
    admittance: Admittance<'_>,
    formulation: Formulation,
) -> Result<Vec<BoundarySolution>, HelmholtzError> {
    let ppw = admit_operator(surface, k, medium, drive_velocities)?;
    for &drive_velocity in drive_velocities {
        admit_bc(
            surface,
            k,
            medium,
            SurfaceBc {
                drive_velocity,
                admittance,
            },
        )?;
    }
    validate_closed(surface)?;
    if admittance.is_zero() {
        return solve_radiation_batch(surface, k, medium, drive_velocities, formulation)?
            .into_iter()
            .zip(drive_velocities)
            .map(|(solution, &drive_velocity)| {
                finish_boundary(
                    surface,
                    solution,
                    SurfaceBc {
                        drive_velocity,
                        admittance,
                    },
                )
            })
            .collect();
    }
    let n = surface.areas().len();
    let (mut a, b) = assemble_dense_with_static(surface, k, alpha_for(formulation, k), None);
    let ik = C64::new(0.0, k);
    match admittance {
        Admittance::Local { beta, .. } => {
            for i in 0..n {
                for j in 0..n {
                    a[i * n + j] = a[i * n + j] + ik * b[i * n + j] * beta[j];
                }
            }
        }
        _ => {
            for j in 0..n {
                let mut e = vec![C64::ZERO; n];
                e[j] = C64::ONE;
                let col = admittance.image(&e)?;
                for i in 0..n {
                    for q in 0..n {
                        a[i * n + j] = a[i * n + j] + ik * b[i * n + q] * col[q];
                    }
                }
            }
        }
    }
    let lu = lu_complex(&a, n).map_err(|_| HelmholtzError::Singular)?;
    let condition_lower_bound = condition_lower_bound(&a, n, &lu);
    drive_velocities
        .iter()
        .map(|&drive_velocity| {
            let q: Vec<_> = drive_velocity
                .iter()
                .map(|&v| C64::new(0.0, k * medium.sound_speed * medium.density) * v)
                .collect();
            let mut pressure = vec![C64::ZERO; n];
            for i in 0..n {
                for j in 0..n {
                    pressure[i] = pressure[i] + b[i * n + j] * q[j];
                }
            }
            lu.solve(&mut pressure);
            if pressure.iter().any(|v| !finite(*v)) {
                return Err(bad("nonfinite Robin pressure"));
            }
            let absorbed = admittance.image(&pressure)?;
            let velocity: Vec<_> = drive_velocity
                .iter()
                .zip(absorbed)
                .map(|(&v, y)| v - y.scale(1.0 / (medium.density * medium.sound_speed)))
                .collect();
            let (radiated_power, radiated_power_roundoff_interval, power_diagnostics) =
                evaluate_radiated_power(&pressure, &velocity, surface.areas(), medium);
            finish_boundary(
                surface,
                RadiationSolution {
                    pressure,
                    velocity,
                    k,
                    medium,
                    surface_fingerprint: surface_fingerprint(surface),
                    panels_per_wavelength: ppw,
                    radiated_power,
                    radiated_power_roundoff_interval,
                    power_diagnostics,
                    condition_lower_bound,
                    dense_cap_utilization: n as f64 / MAX_DENSE_PANELS as f64,
                },
                SurfaceBc {
                    drive_velocity,
                    admittance,
                },
            )
        })
        .collect()
}

/// The .109 operator/GMRES path with exactly the same Robin substitution.
/// FMM order, work, near blocks, cancellation and residual gates remain owned
/// by IterativeRadiation; an unconverged trace never becomes a solution.
pub fn solve_radiation_bc_iterative(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    bc: SurfaceBc<'_>,
    formulation: Formulation,
    options: &IterativeOptions,
    gate: &CancelGate,
) -> Result<
    (
        BoundarySolution,
        fs_solver::krylov::SolveReport,
        ResourceReceipt,
    ),
    IterativeError,
> {
    let op = IterativeRadiation::prepare_bc(
        surface,
        k,
        medium,
        formulation,
        options,
        bc.admittance,
        gate,
    )?;
    let mut state = op.start(bc.drive_velocity, gate)?;
    op.advance(&mut state, options.max_iters.saturating_add(1), gate)?;
    let (solution, report, receipt) = op.finish(&state)?;
    Ok((finish_boundary(surface, solution, bc)?, report, receipt))
}
