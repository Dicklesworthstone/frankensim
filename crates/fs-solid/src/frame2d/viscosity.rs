//! Original material section viscosity on the native frame/Newmark trajectory.
//!
//! Integrate axial, curvature and explicitly bound transverse shear strain
//! rates in the same reference member. Shear-deformable members use the SAME
//! static Timoshenko interpolation as their elastic and inertial owners.
//! The full assembled C retains joint and directional cross terms;
//! heterogeneous members need not share one proportional damping ratio.
//! Material work is exported as dissipation, not silently deposited into heat.

use super::{Frame2d, Geometry, SectionModel, SolidError};
use fs_material::section::{
    ConcentricElasticSection, ConcentricKelvinVoigtSection, axial::ConcentricAxialViscosity,
    transverse::ConcentricTransverseViscosity,
};

/// Explicit axial AND bending viscosities, with optional transverse shear,
/// from one original material section.
/// Each original layer/source and its frequency applicability remain attached.
#[derive(Clone, Debug)]
pub struct FrameViscousSection {
    axial: ConcentricAxialViscosity,
    bending: ConcentricKelvinVoigtSection,
    transverse: Option<ConcentricTransverseViscosity>,
    band: (f64, f64),
}

/// One physical scalar dashpot and its work-conjugate nodal attachment.
/// Extension is `sum(weights[i]*u[dofs[i]])` in metres. Multiplying the
/// extension rate by `damping_n_s_m` gives its force in newtons; the SAME
/// weights return nodal forces and moments. No matrix factorization is needed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameViscousPort {
    /// Original member index, matching the retained `sections()` bindings.
    pub element: usize,
    /// Original node-major `[ux_i, uy_i, theta_i, ux_j, uy_j, theta_j]` DOFs.
    pub dofs: [usize; 6],
    /// Translation weights are dimensionless; rotation weights have metres.
    pub displacement_weights: [f64; 6],
    /// Nonnegative physical damping. Explicit source zero remains zero.
    pub damping_n_s_m: f64,
}

/// Immutable physical strain-rate operator with its original material sources.
/// Supports, modal projections and native dynamics all use the same full C.
#[derive(Clone, Debug)]
pub struct FrameMaterialDamping {
    ndof: usize,
    entries: Vec<f64>,
    ports: Vec<FrameViscousPort>,
    free_dofs: Vec<usize>,
    sections: Vec<(usize, FrameViscousSection)>,
    band: (f64, f64),
    // Frequency-only assembly has no timestep claim. Only the public transient
    // preparation fills this after the complete native spectrum is admitted.
    maximum_omega_rad_s: Option<f64>,
    covers_all_members: bool,
}

impl FrameMaterialDamping {
    /// Complete node-major coordinate count, including supported coordinates.
    pub const fn ndof(&self) -> usize {
        self.ndof
    }
    /// Full symmetric row-major matrix. C maps `[vx, vy, theta_dot]` to
    /// `[Fx, Fy, moment]`; rotational entries retain their work-conjugate units.
    pub fn entries(&self) -> &[f64] {
        &self.entries
    }
    /// Three exact positive-semidefinite factors per source-bound member:
    /// axial, symmetric bending/shear, antisymmetric bending. Shear and
    /// symmetric curvature share one displacement coordinate; their physical
    /// losses add in that factor. Explicit zero factors are retained, and
    /// `C = sum(damping_n_s_m * weights * weights^T)`.
    pub fn ports(&self) -> &[FrameViscousPort] {
        &self.ports
    }
    /// Original unconstrained full-coordinate indices; no condensation is used.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }
    /// Original member index and every original layered source binding.
    pub fn sections(&self) -> &[(usize, FrameViscousSection)] {
        &self.sections
    }
    /// Every original member, including members fixed at both ends, has its
    /// source-bound elastic and viscous section here. Thermal elastic owners
    /// need complete inventory, not only members visible in the free pencil.
    pub const fn covers_all_members(&self) -> bool {
        self.covers_all_members
    }
    /// Intersection of all participating viscosity/elastic/inertial source bands.
    pub const fn omega_band_rad_s(&self) -> (f64, f64) {
        self.band
    }

    /// Admit the actual receiver's clock. All represented frequencies, including
    /// the complete original spectrum before modal truncation, must fit its
    /// sampled Nyquist band, and that whole band must fit the source domain.
    /// This is not a time-discretization error estimate or anti-aliasing filter.
    pub fn admit_timestep(&self, dt: f64) -> Result<(), SolidError> {
        let maximum_omega_rad_s = self.maximum_omega_rad_s.ok_or_else(|| {
            invalid("transient frame damping requires complete natural-spectrum admission")
        })?;
        let nyquist = core::f64::consts::PI / dt;
        if !dt.is_finite()
            || dt <= 0.0
            || !nyquist.is_finite()
            || self.band.0 != 0.0
            || self.band.1 < nyquist
        {
            return Err(invalid(
                "frame viscosity must support DC through the actual sampled Nyquist band",
            ));
        }
        if maximum_omega_rad_s >= nyquist {
            return Err(invalid(
                "complete frame spectrum must lie below the sampled Nyquist limit",
            ));
        }
        Ok(())
    }
}

impl FrameViscousSection {
    /// Bind the two named solid viscosities without inferring either from the
    /// other, a fluid viscosity, loss angle, material name or fitted decay.
    /// Explicit zero is supported; every layer must supply both declarations.
    pub fn try_new(
        section: ConcentricElasticSection,
        axial_property: &str,
        bending_property: &str,
    ) -> Result<Self, SolidError> {
        let source = |error: fs_material::state_point::MaterialStatePointError| {
            SolidError::MaterialRefused {
                what: error.to_string(),
            }
        };
        let axial = section
            .clone()
            .with_axial_viscosity(axial_property)
            .map_err(source)?;
        let bending = section
            .with_kelvin_voigt_bending(bending_property)
            .map_err(source)?;
        let a = axial.omega_band_rad_s();
        let b = bending.omega_band_rad_s();
        let band = (a.0.max(b.0), a.1.min(b.1));
        if band.1 < band.0 {
            return Err(invalid(
                "axial and bending viscosity frequency domains do not intersect",
            ));
        }
        Ok(Self {
            axial,
            bending,
            transverse: None,
            band,
        })
    }

    /// Bind a shear-deformable Kelvin--Voigt member. The transverse law supplies
    /// its original shear modulus, correction and independently selected shear
    /// viscosity; axial and bending viscosities come from that SAME section.
    /// All three loss laws must have intersecting source frequency domains.
    /// Explicit zero shear viscosity still retains transverse shear elasticity.
    pub fn try_new_with_shear(
        transverse: ConcentricTransverseViscosity,
        axial_property: &str,
        bending_property: &str,
    ) -> Result<Self, SolidError> {
        let mut section = Self::try_new(
            transverse.shear().elastic().clone(),
            axial_property,
            bending_property,
        )?;
        let band = transverse.omega_band_rad_s();
        section.band = (section.band.0.max(band.0), section.band.1.min(band.1));
        if section.band.1 < section.band.0 {
            return Err(invalid(
                "axial, bending and transverse shear viscosity frequency domains do not intersect",
            ));
        }
        section.transverse = Some(transverse);
        Ok(section)
    }

    /// Same geometry, elastic stiffness and inertia used by every loss owner.
    pub const fn elastic(&self) -> &ConcentricElasticSection {
        self.axial.elastic()
    }
    /// Original layered axial viscosity, including its source selector.
    pub const fn axial(&self) -> &ConcentricAxialViscosity {
        &self.axial
    }
    /// Original layered bending viscosity, including its source selector.
    pub const fn bending(&self) -> &ConcentricKelvinVoigtSection {
        &self.bending
    }
    /// Original transverse shear stiffness/loss, when explicitly bound.
    pub const fn transverse(&self) -> Option<&ConcentricTransverseViscosity> {
        self.transverse.as_ref()
    }
    /// Intersection of all original coefficient bands [rad/s].
    pub const fn omega_band_rad_s(&self) -> (f64, f64) {
        self.band
    }
}

impl Frame2d {
    /// Original linear-reference longitudinal strain coordinates.
    ///
    /// Returns node-major DOFs, three work-conjugate displacement columns and
    /// reference length L [m]. Their dot products with displacement are x0
    /// (axial extension), x1 (symmetric bending) and x2 (antisymmetric bending),
    /// all in metres. Constant-section Euler-Bernoulli storage is exactly
    /// `EA*x0^2/(2*L) + EI*(3*x1^2+x2^2)/(2*L^3)`.
    /// The same columns map their scalar forces back to nodal forces/moments.
    ///
    /// This is geometry, not a material law. Elastic, viscous and constitutive
    /// memory consumers must supply their own admitted section integrals once.
    /// No source is evaluated, allocation made or accepted state changed.
    /// Corotational geometry needs current-configuration maps and is refused.
    /// Shear-deformable members also refuse: their section rotation and
    /// curvature are not these Euler–Bernoulli traces.
    pub fn reference_strain_traces(
        &self,
        element: usize,
    ) -> Result<([usize; 6], [[f64; 6]; 3], f64), SolidError> {
        let coordinates = self.reference_member_coordinates(element)?;
        if self.elements[element].transverse_shear.is_some() {
            return Err(invalid(
                "Euler-Bernoulli strain traces do not admit transverse shear",
            ));
        }
        Ok(coordinates)
    }

    // Geometric displacement coordinates, not a constitutive strain law.
    // Shear-deformable viscosity applies its own curvature/shear coefficients;
    // the public longitudinal traces above retain their original refusal.
    fn reference_member_coordinates(
        &self,
        element: usize,
    ) -> Result<([usize; 6], [[f64; 6]; 3], f64), SolidError> {
        if self.geometry != Geometry::Linear {
            return Err(invalid("strain traces require linear reference geometry"));
        }
        let element = self.elements.get(element)
            .ok_or_else(|| invalid("strain trace member does not exist"))?;
        let (l, c, s) = (element.l0, element.cos0, element.sin0);
        if !l.is_finite() || l <= 0.0 || !c.is_finite() || !s.is_finite() {
            return Err(invalid("strain trace geometry is not finite"));
        }
        let columns = [
            [-c, -s, 0.0, c, s, 0.0],
            [-2.0 * s, 2.0 * c, l, 2.0 * s, -2.0 * c, l],
            [0.0, 0.0, l, 0.0, 0.0, -l],
        ];
        let dofs = [
            3 * element.i, 3 * element.i + 1, 3 * element.i + 2,
            3 * element.j, 3 * element.j + 1, 3 * element.j + 2,
        ];
        Ok((dofs, columns, l))
    }

    /// Add an elastic member whose SAME sources also determine its physical
    /// strain-rate damping. Static equilibrium and mass are unchanged by eta;
    /// every later native Newmark call automatically includes these losses.
    /// The elastic coefficients retain `add_material_element`'s unrestricted
    /// frequency support, so static and ordinary undamped modal queries remain
    /// valid. Every explicit viscosity law retains its separate finite band.
    /// A section built with `FrameViscousSection::try_new_with_shear` also uses
    /// its original Timoshenko stiffness, consistent inertia and shear loss.
    ///
    /// This is frozen-temperature, small-strain linear geometry, with no
    /// Maxwell history, creep, thermal feedback or radiation model. All frame
    /// members must remain linear elastic during a material-damped history.
    /// Free coordinates must have positive reference inertia. Source-defined
    /// section spin also supplies rotational inertia in the lumped formulation;
    /// massless raw members can still fail this admission. Static massless
    /// condensation cannot represent a viscous first-order degree of freedom.
    ///
    /// The original source bands must cover DC through pi/dt and every native
    /// mode must be resolved below that Nyquist frequency. This admits the
    /// represented sampled band, not aliasing or temporal accuracy of a physical
    /// excitation. Refine the trajectory inside the declared source bandwidth.
    pub fn add_material_viscous_element(
        &mut self,
        i: usize,
        j: usize,
        np: usize,
        section: FrameViscousSection,
    ) -> Result<usize, SolidError> {
        if self.geometry != Geometry::Linear {
            return Err(invalid(
                "source frame viscosity requires fixed linear geometry",
            ));
        }
        let element = match &section.transverse {
            Some(transverse) => {
                self.add_material_shear_element(i, j, np, transverse.shear().clone())?
            }
            None => self.add_material_element(i, j, np, section.elastic().clone())?,
        };
        self.elements[element].viscous_section = Some(section);
        Ok(element)
    }

    /// Original explicit viscosity bindings; purely elastic members return None.
    pub fn material_viscous_section(&self, element: usize) -> Option<&FrameViscousSection> {
        self.elements.get(element)?.viscous_section.as_ref()
    }

    /// Prepare the same physical operator used by native Newmark. A receiving
    /// reduced/dynamic owner must call `admit_timestep` with its actual clock;
    /// it must preserve off-diagonal terms or explicitly bound their removal.
    /// Purely elastic frames return None without changing their existing model.
    pub fn material_damping_operator(&self) -> Result<Option<FrameMaterialDamping>, SolidError> {
        let Some(mut operator) = self.assemble_material_damping()? else {
            return Ok(None);
        };
        let mut maximum_omega_rad_s = 0.0_f64;
        for mode in self.modal_analysis()?.modes {
            let omega = mode.angular_frequency().ok_or_else(|| {
                invalid("complete material-damped frame spectrum must be positive and resolved")
            })?;
            maximum_omega_rad_s = maximum_omega_rad_s.max(omega);
        }
        operator.maximum_omega_rad_s = Some(maximum_omega_rad_s);
        Ok(Some(operator))
    }

    // One source-owned C assembly for both native transients and harmonic
    // queries. A single forcing frequency needs neither a Nyquist clock nor
    // a complete positive modal spectrum: real rigid motion remains admissible.
    pub(super) fn assemble_material_damping(&self) -> Result<Option<FrameMaterialDamping>, SolidError> {
        if self
            .elements
            .iter()
            .all(|element| element.viscous_section.is_none())
        {
            return Ok(None);
        }
        if self.geometry != Geometry::Linear
            || self.elements.iter().any(|element| {
                element.sections.iter().any(|section| {
                    !matches!(section,
                SectionModel::Elastic { ea, ei } if ea.is_finite() && *ea > 0.0
                    && ei.is_finite() && *ei > 0.0)
                })
            })
        {
            return Err(invalid(
                "every member in a material-damped frame must use fixed linear elasticity",
            ));
        }
        let mass = self.reference_mass()?;
        if !self.mass_partition(&mass).1.is_empty() {
            return Err(invalid(
                "material damping requires inertia on every free coordinate; use consistent member mass",
            ));
        }
        let n = self.ndof();
        let entries = n
            .checked_mul(n)
            .ok_or_else(|| invalid("frame damping size overflow"))?;
        let mut matrix = vec![0.0; entries];
        let mut ports = Vec::new();
        let mut sections = Vec::new();
        let mut band = (0.0_f64, f64::MAX);
        for (index, element) in self.elements.iter().enumerate() {
            let Some(section) = &element.viscous_section else {
                continue;
            };
            sections.push((index, section.clone()));
            band = (band.0.max(section.band.0), band.1.min(section.band.1));
            let (dofs, columns, l) = self.reference_member_coordinates(index)?;
            let ca = section.axial.viscous_rigidity_n_s() / l;
            let cb = section.bending.viscous_stiffness_n_m2_s() / l / l / l;
            if (section.axial.viscous_rigidity_n_s() > 0.0 && ca == 0.0)
                || (section.bending.viscous_stiffness_n_m2_s() > 0.0 && cb == 0.0)
                || !ca.is_finite()
                || !cb.is_finite()
            {
                return Err(invalid(
                    "material section viscosity is unrepresentable on the member geometry",
                ));
            }
            let symmetric = match (&element.transverse_shear, &section.transverse) {
                (Some(shear), Some(transverse)) => {
                    // For x1 = 2*(v_i-v_j)+L*(theta_i+theta_j), the SAME
                    // interpolation used by consistent_mass gives
                    // gamma_dot = -s*x1_dot/(2L) and
                    // kappa_dot = -(3*beta*(1-2*x/L)*x1_dot+x2_dot)/L².
                    // Exact integration leaves three PSD factors, including
                    // both independently sourced losses in the x1 factor.
                    let beta = shear.bending_fraction;
                    let fraction = shear.shear_fraction;
                    let bending = (cb * beta) * (3.0 * beta);
                    let shear_viscosity = transverse.viscous_rigidity_n_s();
                    let transverse = (shear_viscosity * fraction / l) * (0.25 * fraction);
                    if (cb > 0.0 && bending == 0.0)
                        || (shear_viscosity > 0.0 && transverse == 0.0)
                    {
                        return Err(invalid(
                            "Timoshenko section viscosity is unrepresentable on the member geometry",
                        ));
                    }
                    bending + transverse
                }
                (None, None) => 3.0 * cb,
                _ => return Err(invalid("frame shear stiffness and viscosity bindings differ")),
            };
            if !symmetric.is_finite() {
                return Err(invalid(
                    "material section viscosity is unrepresentable on the member geometry",
                ));
            }
            for (damping_n_s_m, displacement_weights) in
                [ca, symmetric, cb].into_iter().zip(columns)
            {
                ports.push(FrameViscousPort {
                    element: index,
                    dofs,
                    displacement_weights,
                    damping_n_s_m,
                });
            }
        }
        // These exact physical ports are also exported to reciprocal coupling.
        // Native forces and reduced ports cannot drift to different C operators.
        for port in &ports {
            let weights = port.displacement_weights;
            let dofs = port.dofs;
            for i in 0..6 {
                for j in i..6 {
                    let value = port.damping_n_s_m * weights[i] * weights[j];
                    let entry = &mut matrix[dofs[i] * n + dofs[j]];
                    *entry += value;
                    if !value.is_finite()
                        || !entry.is_finite()
                        || (port.damping_n_s_m > 0.0
                            && weights[i] != 0.0
                            && weights[j] != 0.0
                            && value == 0.0)
                    {
                        return Err(invalid("source-derived frame damping is unrepresentable"));
                    }
                    if i != j {
                        matrix[dofs[j] * n + dofs[i]] += value;
                    }
                }
            }
        }
        if band.1 < band.0 {
            return Err(invalid(
                "frame member viscosity source domains do not intersect",
            ));
        }
        Ok(Some(FrameMaterialDamping {
            ndof: n,
            entries: matrix,
            ports,
            free_dofs: (0..n).filter(|&dof| !self.fixed[dof]).collect(),
            sections,
            band,
            maximum_omega_rad_s: None,
            covers_all_members: self.elements.iter().all(|e| e.viscous_section.is_some()),
        }))
    }
}

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}

#[cfg(test)]
mod strain_tests;
